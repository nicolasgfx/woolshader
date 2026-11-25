// The knitting port on the device (gpu/knitting.h): GpuKnitNumerics implements KnitNumerics
// (knitting/numerics.h) by composing the device capabilities: gate (a) (validation/), the relaxation
// (solver/), the arc-length tables (geometry/), and the surfaces, placement, lattice, gate (d) and
// boundary checks (surface/). The solver's tuned constants live here. A failed device call is an
// exception, never a CPU stand-in. Its validations are GPU diagnostic verdicts.
#include "gpu/knitting.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "gpu/device.h"
#include "gpu/device/runtime.h"
#include "gpu/geometry/device_geometry.h"
#include "gpu/port/ncp_call_file.h"
#include "gpu/solver/ncp.h"
#include "gpu/surface/placement.h"
#include "gpu/surfaces.h"
#include "gpu/validation/checks.h"
#include "knitting/relax/settle.h"
#include "model/cancel.h"
#include "model/env.h"

namespace knit::gpu {

// Timings of the device numerics (wall time of each call including its uploads and downloads;
// summed over calls), and counts.
struct KnitNumericsStats {
    long long validateCalls = 0, storeCalls = 0, curvatureCalls = 0;
    double msUpload = 0, msCurvature = 0, msSampling = 0, msSignature = 0, msClearance = 0, msSummarize = 0;
    double msStore = 0, msCurvatureOnly = 0;
    // The relax calls.
    long long relaxCalls = 0, converged = 0, stopped = 0, unresolved = 0, infeasible = 0, pcgIterations = 0,
              lmSteps = 0, maxUnknowns = 0, stepBudgetHits = 0, maxLmStepsDone = 0;
    double msRelax = 0, msRelaxMax = 0, peakMB = 0;
    std::string firstUnresolved;
    std::vector<std::string> callLog;  // one line per relax call
    std::string batchLine() const;
    std::string relaxLine() const;
};

// The port on the device. One object per run (a worker): its placement statistics take no lock.
class GpuKnitNumerics final : public KnitNumerics {
public:
    explicit GpuKnitNumerics(const NcpNumericsOptions& options) : options_(options) {}
    const char* name() const override { return "ncp-gpu"; }
    ValidationReport validate(const CurveSet& cs, const PatternParams& pp, const ValidateOptions& opt) override;
    CurvatureResult curvature(const YarnPath& y, double radius, const Allowances& al, const ExecPolicy& exec) override;
    void relax(const RelaxRequest& request) override;
    void store(CurveSet& cs, const Controls& P, const ExecPolicy& exec) override;
    std::unique_ptr<MeshSurface> meshSurface(const TriMesh& mesh) override;
    std::unique_ptr<BaseSurface> pieceQuery(const Chart& chart) override;
    void classifyLattice(const UvAtlas& atlas, const LatticeChartQuery& q, const LatticeSolids& solids,
                         std::vector<LatticeNodeClass>& out) override;
    SurfaceCheck checkRenderedSurface(const CurveSet& cs, const PatternParams& ppR, const SurfaceParams& sp) override;
    BoundaryCheck checkBoundary(const CurveSet& cs, const SurfaceBinding& bind, const MeshSurface& surf, double Rs) override;
    std::string statsText() const;
    std::string placementText() const { return placement_.line(); }
    std::vector<std::string> relaxCallLog() const;

private:
    NcpNumericsOptions options_;
    KnitNumericsStats stats_;
    PlacementStats placement_;  // the surfaces made here write to it; one worker, no lock
    mutable std::mutex mutex_;  // the relax figures of stats_ and the capture count
    long long captured_ = 0;    // relax calls written to options_.captureDir
};

namespace {

constexpr double kAnchorWeight = 0.01;  // KnitProblemParams::wA
constexpr double kDelta = 0.3;          // KnitProblemParams::delta [mm]
// Above this many unknowns the coverage delta is kDeltaLarge: the candidate rows grow with the pad
// (sphere component, 480k unknowns: 9.2 GB peak at 0.3 mm, 6.3 GB at 0.15 mm, same time and quality).
constexpr long long kDeltaLargeAbove = 200000;
constexpr double kDeltaLarge = 0.15;
constexpr double kHardMarginContact = 0.03;    // [mm] hard contact = D + c + this (validation allowance 2 eps + our polyline sagitta)
constexpr double kHardMarginCurvature = 0.03;  // hard curvature = 1 / (r + c/2) / (1 + this) (sampled vs certified curvature)
constexpr int kCurvatureSamples = 9;           // KnitProblemParams::curvatureSamples

// The iteration limits are safety nets (reported as UNRESOLVED when hit): with the hybrid contact
// rows and the double-precision band the large patches converge (teddy body round 0: KKT 1e-4 in
// 17.6 s; earlier budgets of 15 / 25 only masked a stall, analysis/linacc2).
NcpSolveParams solveParams() {
    NcpSolveParams p;
    p.maxOuter = 40;
    // Early stops: the progress rule (1 %) kept every gate (a) verdict on
    // the measured patches and scenes (teddy 276 -> 123 s, sphere 861 -> 215 s); the infeasible rule
    // ends calls whose violation stopped shrinking at large rho (86 -> 16 s on a man head patch).
    p.progressTol = 1e-2;
    p.infeasibleWindow = 3;
    // Step budget per call: 2x the most any converged or stopped call
    // took on the calibrated teddy (450). A safety net against runaway calls: no measured call
    // reached it (teddy D 1.46: at most 563; its slow round-0 call took 282 steps).
    p.maxLmSteps = 900;
    return p;
}

// Settle's relax call (knit-perf 2026-11-07): the progress stop at 5 % and the
// contact rows feasible within 0.005 (rows are normalised by the move threshold, about 1.31 mm on the
// teddy, so 0.0066 mm; not covered by the 0.03 mm hard margin in the worst case);
// curvature, band and cap stay at 1e-3. On the Rustic teddy body / head / ear this removed 40-50 % of
// the steps (the augmented Lagrangian's endgame on the last contact violation) with the same
// clearance and signature verdicts; repair keeps solveParams().
NcpSolveParams settleSolveParams() {
    NcpSolveParams p = solveParams();
    p.progressTol = 0.05;
    p.tolViolContact = 0.005;
    return p;
}

// Whether a call ended above its feasibility allowances, measured as the device's stop test measures
// them (ncp_control.cu feasibilityViolation): contact against its own allowance, the other hard
// families against tolViol.
bool aboveAllowance(const KnitProblemScore& score, const NcpSolveParams& solve) {
    for (int f = 0; f < kRowFamilies; ++f) {
        const RowFamily family = static_cast<RowFamily>(f);
        if (!isHard(family)) continue;
        const double allowance = family == RowFamily::Contact && solve.tolViolContact > 0 ? solve.tolViolContact : solve.tolViol;
        if (score.maxViol[f] > allowance) return true;
    }
    return false;
}

// The cap of one solve [bytes]: the granted budget (solveCapBytes), lowered to the caller's cap when it
// sets one.
size_t solveCapFor(double memoryCapMB) {
    const size_t granted = solveCapBytes();
    if (memoryCapMB <= 0) return granted;
    return std::min(granted, static_cast<size_t>(memoryCapMB * 1048576.0));
}

// One relax call's line in the stage profile (ncpNumericsCalls).
std::string callLineOf(const RelaxCall& ci, long long freeN, const NcpStats& st, double ms, double capMB, const NcpSolveReport& rep) {
    char call[400];
    std::snprintf(call, sizeof call, "%s round %d stage %d: %lld free controls, %lld unknowns, %lld LM steps, %lld PCG iterations, %.0f ms, cap %.0f MiB: %s",
                  ci.phase && *ci.phase ? ci.phase : "relax", ci.round, ci.stage, freeN, st.unknowns, st.lmSteps, st.pcgIterations, ms, capMB, rep.status.c_str());
    // The solver's final maxima of the normalized hard rows, per family.
    const auto viol = [&](RowFamily f) { return rep.score.maxViol[static_cast<int>(f)]; };
    char fam[200];
    std::snprintf(fam, sizeof fam, "; final max violation: curvature %.3g, contact %.3g, band %.3g / %.3g, cap %.3g (rho %.3g)",
                  viol(RowFamily::Curvature), viol(RowFamily::Contact), viol(RowFamily::BandUpper), viol(RowFamily::BandLower),
                  viol(RowFamily::Cap), rep.rho);
    return std::string(call) + fam;
}

long long freeControlCount(const ControlMask& freeCtrl) {
    long long n = 0;
    for (const auto& y : freeCtrl)
        for (char f : y) n += f != 0;
    return n;
}

// The refusal of a solve its cap stopped (solveNcp's message `err`).
std::string refusalText(long long freeN, size_t capBytes, const std::string& err) {
    char b[200];
    std::snprintf(b, sizeof b, "ncp-gpu: this relaxation (%lld free controls) needs more GPU memory than its cap of %.0f MB; ", freeN, capBytes / 1048576.0);
    return b + err + ". Refused instead of overcommitting the GPU (it would page to system memory and run many times slower); a smaller part, a larger "
                     "gauge or a thinner yarn needs less.";
}

}  // namespace

ContactStorage contactStorageFromEnvironment() {
    const std::optional<std::string> value = envVar("KNIT_NCP_CONTACT_STORAGE");
    return value && *value == "live" ? ContactStorage::LiveRows : ContactStorage::CompactRows;
}

CurvatureRows curvatureRowsFromEnvironment() {
    const std::optional<std::string> value = envVar("KNIT_NCP_CURVATURE_ROWS");
    return value && *value == "span-max" ? CurvatureRows::SamplesAndSpanMaximum : CurvatureRows::Samples;
}

std::string ncpSolverIdentity() { return ncpSolverIdentity(contactStorageFromEnvironment(), curvatureRowsFromEnvironment()); }

std::string ncpSolverIdentity(ContactStorage storage, CurvatureRows curvature) {
    const NcpSolveParams p = solveParams(), s = settleSolveParams();
    char b[400];
    std::snprintf(b, sizeof b, "ncp-gpu v5 (GEOMETRY 8.21; device-controlled iteration, matrix-free PCG, compact contact rows; early stops: progress %g, infeasible window %d at rho %g, "
                  "step budget %lld; settle: one round, progress %g, contact tolViol %g; repair: one certified move)",
                  p.progressTol, p.infeasibleWindow, p.infeasibleRho, p.maxLmSteps, s.progressTol, s.tolViolContact);
    std::string id = b;
    switch (storage) {
    case ContactStorage::CompactRows: break;
    case ContactStorage::LiveRows: id += " + live contact rows r1 (experimental: near cutoff max(hardC, thr) + 2 delta, live list)"; break;
    }
    switch (curvature) {
    case CurvatureRows::Samples: break;
    case CurvatureRows::SamplesAndSpanMaximum: id += " + span-maximum curvature rows r1 (experimental: one row per span at its curvature maximum)"; break;
    }
    return id;
}

std::string KnitNumericsStats::batchLine() const {
    char b[512];
    std::snprintf(b, sizeof b,
                  "GPU knit numerics: validate %lld calls (upload %.0f, curvature %.0f, sampling %.0f, signature %.0f, clearance %.0f, "
                  "summarize %.0f ms); store %lld calls %.0f ms; curvature %lld calls %.0f ms",
                  validateCalls, msUpload, msCurvature, msSampling, msSignature, msClearance, msSummarize, storeCalls, msStore, curvatureCalls, msCurvatureOnly);
    return b;
}

std::string KnitNumericsStats::relaxLine() const {
    char b[768];
    std::snprintf(b, sizeof b,
                  "ncp-gpu: %lld relax calls %.0f ms (max %.0f), %lld converged, %lld stopped (progress), %lld unresolved (%lld infeasible, "
                  "%lld step budget %lld)%s%s; "
                  "%lld LM steps (max %lld in a converged or stopped call), %lld PCG iterations; max %lld unknowns, "
                  "device peak %.0f MB",
                  relaxCalls, msRelax, msRelaxMax, converged, stopped, unresolved, infeasible, stepBudgetHits, solveParams().maxLmSteps,
                  firstUnresolved.empty() ? "" : ", first: ", firstUnresolved.c_str(), lmSteps, maxLmStepsDone, pcgIterations, maxUnknowns, peakMB);
    return b;
}

ValidationReport GpuKnitNumerics::validate(const CurveSet& cs, const PatternParams& pp, const ValidateOptions& opt) {
    Reservation lease(Reservation::exclusiveBytes());  // runs alone on the device (runtime.h)
    ValidateStats& vs = validateStats();
    const auto t0 = Clock::now();
    ++stats_.validateCalls;
    ValidationReport r;
    r.allowances = Allowances::forExtent(cs.extent());
    r.gpuDiagnostic = true;
    cancelCheckpoint();  // after the device reservation
    auto t = Clock::now();
    const UploadedCurves d = uploadCurves(cs);
    stats_.msUpload += elapsedMs(t);
    cancelCheckpoint();
    t = Clock::now();
    r.curvature = checkCurvatureOnDevice(d, pp.radius(), r.allowances);
    const double msCurv = elapsedMs(t);
    stats_.msCurvature += msCurv;
    t = Clock::now();
    const DevicePolylines polylines = samplePolylinesOnDevice(d, pp.sampleTolerance, 512, r.polylines);
    const double msSamp = elapsedMs(t);
    stats_.msSampling += msSamp;
    cancelCheckpoint();
    t = Clock::now();
    if (cs.kind == CurveKind::FlatPatch)
        throw DeviceError("GPU validate: flat patches (the global projected crossing signature) have no device version");
    r.localSignature = checkSignatureOnDevice(d);
    // Crossing attribution (checks.h): the device port computes ownership-span attribution; the
    // correction re-evaluates the few checks that can change, on the CPU, with the same code as
    // checkLocalSignature.
    applyCrossingAttribution(cs, r.localSignature, opt.attributionWindow);
    const double msSig = elapsedMs(t);
    stats_.msSignature += msSig;
    cancelCheckpoint();
    t = Clock::now();
    r.clearance = checkClearanceOnDevice({d, polylines, r.polylines, r.curvature, pp, r.allowances});
    const double msClear = elapsedMs(t);
    stats_.msClearance += msClear;
    t = Clock::now();
    summarize(r, cs, pp, opt);
    stats_.msSummarize += elapsedMs(t);
    ++vs.calls;
    vs.msCurvature += msCurv;
    vs.msSampling += msSamp;
    vs.msSignature += msSig;
    vs.msClearance += msClear;
    vs.ms += elapsedMs(t0);
    return r;
}

void GpuKnitNumerics::relax(const RelaxRequest& request) {
    const auto t0 = Clock::now();
    // The problem (knit_problem.h): hard thresholds at what the next validation judges (the call's
    // validation parameters, RelaxMoveParams::call) plus certified margins; the move's own
    // thresholds as soft targets. Outside settle / repair (no call info) the move's thresholds are hard.
    KnitProblemParams kp;
    kp.wA = kAnchorWeight;
    const long long freeN = freeControlCount(request.freeCtrl);
    kp.delta = 3 * freeN > kDeltaLargeAbove ? kDeltaLarge : kDelta;
    kp.curvatureSamples = kCurvatureSamples;
    kp.contactModel = KnitProblemParams::ContactModel::Hybrid;
    const RelaxCall& ci = request.mp.call;
    if (ci.phase && *ci.phase) {
        const HardThresholds hard = hardThresholds(ci.validatePP, kHardMarginContact, kHardMarginCurvature);
        kp.hardContact = hard.contact;
        kp.hardKappa = hard.kappa;
    }
    NcpOptions options;
    options.contactStorage = options_.contactStorage;
    options.curvatureRows = options_.curvatureRows;
    const bool settleCall = ci.phase && std::string(ci.phase) == "settle";
    NcpSolveParams solve = settleCall ? settleSolveParams() : solveParams();
    if (options_.slowCall == SlowCall::UseStepBudget) solve.infeasibleWindow = 0;
    if (!options_.captureDir.empty()) {
        long long n = 0;
        {
            std::lock_guard<std::mutex> lock(mutex_);
            n = ++captured_;
        }
        char name[32];
        std::snprintf(name, sizeof name, "/call_%03lld.ncpcall", n);
        std::string why;
        if (!writeNcpCall(options_.captureDir + name, request, kp, solve, &why)) throw std::runtime_error(why);
    }
    NcpSolveReport rep;
    NcpStats st;
    std::string err;
    size_t capBytes = 0;
    {
        // The solve holds the whole device (peak memory grows with the problem): concurrent pieces
        // wait instead of overcommitting.
        Reservation lease(Reservation::exclusiveBytes());
        // The lease only orders the jobs; the cap stops the solve itself at the budget.
        capBytes = solveCapFor(options_.memoryCapMB);
        AllocationCap cap(capBytes);
        if (!solveNcp({request.cs, request.P, request.Prest, request.freeCtrl, request.surface, request.bind, request.mp, kp, solve, options}, rep, &st, &err)) {
            if (AllocationCap::capHit()) throw DeviceBudgetExceeded(refusalText(freeN, capBytes, err));
            throw std::runtime_error("ncp-gpu relax: " + err);
        }
    }
    const double ms = elapsedMs(t0);
    const std::string callLine = callLineOf(ci, freeN, st, ms, capBytes / 1048576.0, rep);
    std::lock_guard<std::mutex> lock(mutex_);
    KnitNumericsStats& t = stats_;
    ++t.relaxCalls;
    t.msRelax += ms;
    t.msRelaxMax = std::max(t.msRelaxMax, ms);
    if (rep.converged) ++t.converged;
    else if (rep.stopped) ++t.stopped;
    else ++t.unresolved;
    if (rep.budgetHit) ++t.stepBudgetHits;
    if (rep.converged || rep.stopped) t.maxLmStepsDone = std::max(t.maxLmStepsDone, st.lmSteps);
    if (!rep.converged && !rep.stopped && t.firstUnresolved.empty()) t.firstUnresolved = rep.status;
    if (aboveAllowance(rep.score, solve)) ++t.infeasible;
    t.pcgIterations += st.pcgIterations;
    t.lmSteps += st.lmSteps;
    t.maxUnknowns = std::max(t.maxUnknowns, st.unknowns);
    t.peakMB = std::max(t.peakMB, st.devPeakMB);
    // The solve's phases (NcpStats::line: refresh = broad phase, trial, swept test, pcg, ...).
    t.callLog.emplace_back(callLine + "\n         " + st.line() + "\n         " + st.memoryLine());
}

void GpuKnitNumerics::store(CurveSet& cs, const Controls& P, const ExecPolicy&) {
    Reservation lease(Reservation::exclusiveBytes());  // runs alone on the device (runtime.h)
    const auto t0 = Clock::now();
    ++stats_.storeCalls;
    for (size_t y = 0; y < cs.yarns.size(); ++y) cs.yarns[y].curve.ctrl = P[y];
    computeArcLengthTables(cs);
    stats_.msStore += elapsedMs(t0);
}

std::unique_ptr<MeshSurface> GpuKnitNumerics::meshSurface(const TriMesh& mesh) {
    std::string e;
    std::unique_ptr<MeshSurface> surface = makeMeshSurface(mesh, &placement_, &e);
    if (!surface) throw std::runtime_error("GPU mesh surface: " + e);
    return surface;
}

std::unique_ptr<BaseSurface> GpuKnitNumerics::pieceQuery(const Chart& chart) {
    std::string e;
    std::unique_ptr<BaseSurface> query = makePieceQuery(chart, &placement_, &e);
    if (!query) throw std::runtime_error("GPU placement: " + e);
    return query;
}

void GpuKnitNumerics::classifyLattice(const UvAtlas& atlas, const LatticeChartQuery& q, const LatticeSolids& solids,
                                      std::vector<LatticeNodeClass>& out) {
    std::string e;
    if (!gpu::classifyLattice({atlas, q, solids.meshes, solids.offset}, out, &e)) throw PartFailure("lattice node classification: " + e);
}

SurfaceCheck GpuKnitNumerics::checkRenderedSurface(const CurveSet& cs, const PatternParams& ppR, const SurfaceParams& sp) {
    SurfaceCheck check;
    std::string e;
    if (!gpu::checkSurface({cs, ppR, sp}, {}, check, nullptr, &e)) throw std::runtime_error("gate (d): " + e);
    return check;
}

BoundaryCheck GpuKnitNumerics::checkBoundary(const CurveSet& cs, const SurfaceBinding& bind, const MeshSurface& surf, double Rs) {
    BoundaryCheck check;
    std::string e;
    if (!gpu::checkBoundary({cs, bind, surf, Rs, 0.1, 1e-3}, check, &e)) throw std::runtime_error("GPU boundary check: " + e);
    return check;
}

CurvatureResult GpuKnitNumerics::curvature(const YarnPath& y, double radius, const Allowances& al, const ExecPolicy&) {
    Reservation lease(Reservation::exclusiveBytes());  // runs alone on the device (runtime.h)
    const auto t0 = Clock::now();
    ++stats_.curvatureCalls;
    CurveSet one;
    one.yarns.push_back(y);
    const std::vector<CurvatureResult> out = checkCurvatureOnDevice(uploadCurves(one), radius, al);
    stats_.msCurvatureOnly += elapsedMs(t0);
    return out.at(0);
}

std::unique_ptr<KnitNumerics> makeKnitNumerics(std::string* err, const NcpNumericsOptions& o) {
    if (!available(err)) return nullptr;
    return std::make_unique<GpuKnitNumerics>(o);
}

std::string GpuKnitNumerics::statsText() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return stats_.relaxLine() + "\n     " + stats_.batchLine();
}

std::vector<std::string> GpuKnitNumerics::relaxCallLog() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return stats_.callLog;
}

std::string knitNumericsStats(const KnitNumerics& n) {
    const auto* numerics = dynamic_cast<const GpuKnitNumerics*>(&n);
    return numerics ? numerics->statsText() : "not the device numerics";
}

std::string knitPlacementStats(const KnitNumerics& n) {
    const auto* numerics = dynamic_cast<const GpuKnitNumerics*>(&n);
    return numerics ? numerics->placementText() : "not the device numerics";
}

std::vector<std::string> ncpNumericsCalls(const KnitNumerics& n) {
    const auto* numerics = dynamic_cast<const GpuKnitNumerics*>(&n);
    return numerics ? numerics->relaxCallLog() : std::vector<std::string>{};
}

namespace {

// FNV-1a 64 of the controls' bits in yarn order: two solves gave the same result exactly when these agree.
unsigned long long controlsDigest(const Controls& P) {
    unsigned long long h = 1469598103934665603ull;
    for (const std::vector<Vec3>& y : P)
        for (const Vec3& p : y)
            for (double v : {p.x, p.y, p.z}) {
                unsigned char b[sizeof v];
                std::memcpy(b, &v, sizeof v);
                for (unsigned char c : b) h = (h ^ c) * 1099511628211ull;
            }
    return h;
}

bool writeControlsFile(const std::string& path, const Controls& P, std::string* err) {
    std::ofstream f(path, std::ios::binary);
    for (const std::vector<Vec3>& y : P)
        for (const Vec3& p : y)
            for (double v : {p.x, p.y, p.z}) f.write(reinterpret_cast<const char*>(&v), sizeof v);
    if (f) return true;
    if (err) *err = "replay: cannot write " + path;
    return false;
}

}  // namespace

bool replayNcpCall(const NcpReplayRequest& request, NcpReplayResult& result, std::string* err) {
    result = NcpReplayResult{};
    NcpCallRecord rec;
    if (!readNcpCall(request.callFile, rec, err)) return false;
    rec.move.call.phase = rec.phase.c_str();
    if (request.maxLmSteps > 0) rec.solve.maxLmSteps = request.maxLmSteps;
    const long long freeN = freeControlCount(rec.freeCtrl);
    const auto t0 = Clock::now();
    size_t capBytes = 0;
    try {
        std::unique_ptr<MeshSurface> surface;
        if (rec.hasMesh) {
            surface = makeMeshSurface(rec.mesh, nullptr, err);
            if (!surface) return false;
        }
        Controls P = rec.P;
        NcpOptions options;
        options.contactStorage = request.contactStorage;
        options.curvatureRows = request.curvatureRows;
        const NcpCall call{rec.curves, P, rec.Prest, rec.freeCtrl, surface.get(), nullptr, rec.move, rec.problem, rec.solve, options};
        Reservation lease(Reservation::exclusiveBytes());
        capBytes = solveCapFor(request.memoryCapMB);
        result.capMB = capBytes / 1048576.0;
        AllocationCap cap(capBytes);
        NcpStats st;
        NcpSolveReport rep;
        if (request.mode == NcpReplayMode::Frozen) {
            resetPeak();
            NcpDevice device(call);
            std::vector<double> p(size_t(3 * device.vars()));
            for (size_t i = 0; i < p.size(); ++i) p[i] = std::sin(0.37 * static_cast<double>(i) + 1.0);  // a fixed direction
            const std::vector<NcpDevice::FrozenPass> passes = device.frozenLinearizations(request.frozenRho, p);
            std::ofstream f(request.frozenOut, std::ios::binary);
            char line[300];
            for (const NcpDevice::FrozenPass& pass : passes) {
                const double head[4] = {pass.eval.merit, pass.eval.objective, static_cast<double>(pass.contactRows), static_cast<double>(pass.activeRows)};
                f.write(reinterpret_cast<const char*>(head), sizeof head);
                for (const std::vector<double>* v : {&pass.grad, &pass.hblock, &pass.product}) f.write(reinterpret_cast<const char*>(v->data()), v->size() * 8);
                std::snprintf(line, sizeof line, "frozen pass: merit %.17g, objective %.17g, contact rows %lld, active rows %lld, active entries %lld; ",
                              pass.eval.merit, pass.eval.objective, pass.contactRows, pass.activeRows, pass.activeEntries);
                result.statsLine += line;
            }
            if (!f) {
                if (err) *err = "replay: cannot write " + request.frozenOut;
                return false;
            }
            st = device.stats();
            st.devPeakMB = memStats().peak / 1048576.0;
            rep.status = "FROZEN: two linearizations at the start";
            result.ms = elapsedMs(t0);
            result.status = rep.status;
            result.memoryLine = st.memoryLine();
            return true;
        }
        if (request.mode == NcpReplayMode::CountOnly) {
            resetPeak();
            NcpDevice device(call, NcpSizing::CountOnly);
            st = device.stats();
            st.devPeakMB = memStats().peak / 1048576.0;
            rep.status = "COUNTED: the starting sizes only";
        } else {
            std::string why;
            if (!solveNcp(call, rep, &st, &why)) {
                if (!AllocationCap::capHit()) {
                    if (err) *err = "replay: " + why;
                    return false;
                }
                result.refused = true;
                result.status = refusalText(freeN, capBytes, why);
                result.ms = elapsedMs(t0);
                return true;
            }
            result.controlsDigest = controlsDigest(P);
            if (!request.controlsOut.empty() && !writeControlsFile(request.controlsOut, P, err)) return false;
        }
        result.ms = elapsedMs(t0);
        result.status = rep.status;
        result.callLine = callLineOf(rec.move.call, freeN, st, result.ms, result.capMB, rep);
        result.statsLine = st.line();
        result.memoryLine = st.memoryLine();
        return true;
    } catch (const DeviceError& e) {
        if (!AllocationCap::capHit()) return reportFailure(e, err);
        result.refused = true;
        result.status = refusalText(freeN, capBytes, e.what());
        result.ms = elapsedMs(t0);
        return true;
    }
}

bool ncpAvailable(std::string* why) { return available(why); }

ComputeCaps detectComputeCaps() {
    ComputeCaps r;
    r.device = available(&r.whyNoDevice);
    if (r.device) r.deviceName = deviceSummary();
    r.ncp = r.device && ncpAvailable(&r.whyNoNcp);
    if (!r.device) r.whyNoNcp = r.whyNoDevice;
    return r;
}

} // namespace knit::gpu
