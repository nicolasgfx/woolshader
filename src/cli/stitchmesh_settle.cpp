// The yarn-settle commands of stitchmesh_lab: the existing
// knitting solver applied to the embedded yarn of the first-yarn checkpoint.
//
//   yarn-settle-fixtures  yarn-fixtures' cases, each through the pipeline below; refused cases are reported as refused
//   yarn-settle           one frozen mesh (yarn-embed's input, read back, not rebuilt) through the same pipeline
//   ncp-replay            one relax call captured by --capture-calls, solved again or only sized (--mode count)
//
// Both settle commands take --gpu-budget-mib N (the solver's device cap: the smaller of N and the granted budget) and
// --capture-calls <dir> (every relax call's exact input as <dir>/call_<n>.ncpcall). A relax call the cap refuses is a
// refused stage: the verdict says so and the command exits 3.
//
// The pipeline, one directory per state so every source of change can be told apart:
//   raw/        embedYarn and checkYarn exactly as the first-yarn run (same settings)
//   converted/  the yarn as production curves (yarn_curves.h: uniform cubic B-splines through equal-arc points), re-embedded
//               on them and checked again; the conversion's deviation and length change measured
//   settled/    after the solve, re-embedded and checked again
// A case whose raw construction check fails is not settled (the solver is never asked to repair a crossing order), nor
// is one whose conversion changes the construction.
//
// The solve. KnitNumerics::relax (knitting/numerics.h) on the device: production's NCP problem (knit_problem.h: soft
// edge lengths at the converted control polygon, hard contact and curvature rows, a band around the reference, a
// shape anchor to each call's start, a displacement cap) with its own step-wise swept pass-through test. It is a
// constrained geometric post-pass, not a physical yarn simulation. The tubes start interpenetrating, so the target
// diameter is reached by radius continuation: from the start diameter D0 at which the converted curves pass the
// certified checks (the rule of certifiedRepair, knitting/relax/settle.h) in stages D_k = D0 (D / D0)^(k / K),
// K = ceil(log(D / D0) / log(growth)), plus `repeats` calls at D. Each call: hard contact D_k + c (+ the solver's margin),
// hard curvature 1 / (r_k + c / 2), soft targets as settle's, displacement cap displacementFraction (D_{k-1} + c). After
// each call the curves are validated at D_k (certified curvature and clearance, knitting/checks/checks.h; the local
// signature does not apply: there are no production stitch records) and the straight motion of the call is certified
// at D_{k-1} when both its end states pass there (certifyStraightMotion): a certified chain from a verified start proves
// that no centreline passed through another, except same-yarn pairs closer than (pi/2) delta in arc.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <unordered_map>
#include <string>
#include <vector>

#ifdef _WIN32
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <psapi.h>
#endif

#include "stitchmesh_lab_shared.h"
#include "gpu/knitting.h"
#include "knitting/checks/checks.h"
#include "knitting/layout/base_surface.h"
#include "knitting/numerics.h"
#include "knitting/relax/settle.h"
#include "knitting/stitch/stockinette.h"
#include "model/mesh_query.h"
#include "stitchmesh/typed_mesh.h"
#include "stitchmesh/virtual_shift.h"
#include "stitchmesh/yarn_curves.h"
#include "stitchmesh/yarn_embed.h"
#include "stitchmesh/yarn_settle.h"

namespace fs = std::filesystem;
using namespace knit;

namespace {

using Clock = std::chrono::steady_clock;

double msSince(Clock::time_point t0) { return std::chrono::duration<double, std::milli>(Clock::now() - t0).count(); }

void putf(std::ostream& o, const char* fmt, ...) {
    char b[4096];
    va_list ap;
    va_start(ap, fmt);
    std::vsnprintf(b, sizeof b, fmt, ap);
    va_end(ap);
    o << b;
}

std::map<std::string, std::string> parseArgs(int argc, char** argv) {
    std::map<std::string, std::string> kv;
    for (int i = 2; i + 1 < argc; i += 2) kv[argv[i]] = argv[i + 1];
    return kv;
}

double number(const std::map<std::string, std::string>& kv, const std::string& k, double d) {
    const auto it = kv.find(k);
    return it == kv.end() ? d : std::atof(it->second.c_str());
}

// The process's host memory peaks [MiB] (working set and committed private bytes); 0 where unknown.
struct HostPeak {
    double workingSetMiB = 0, privateMiB = 0;
};
HostPeak hostPeak() {
    HostPeak h;
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS c;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &c, sizeof c)) {
        h.workingSetMiB = c.PeakWorkingSetSize / 1048576.0;
        h.privateMiB = c.PeakPagefileUsage / 1048576.0;
    }
#endif
    return h;
}

SettleSettings settleFrom(const std::map<std::string, std::string>& kv) {
    SettleSettings s;
    s.spacing = number(kv, "--spacing", s.spacing);
    s.growth = number(kv, "--growth", s.growth);
    s.maxRepeats = static_cast<int>(number(kv, "--repeats", s.maxRepeats));   // target calls at most
    s.stallCalls = static_cast<int>(number(kv, "--stall-calls", s.stallCalls));
    s.budgetMinutes = number(kv, "--budget-minutes", s.budgetMinutes);
    s.displacementFraction = number(kv, "--cap-fraction", s.displacementFraction);
    s.band = number(kv, "--band", 1) != 0;
    s.convertOnly = number(kv, "--convert-only", 0) != 0;
    s.samplesPerSpan = static_cast<int>(number(kv, "--samples-per-span", s.samplesPerSpan));
    s.perturbUlps = static_cast<int>(number(kv, "--perturb-ulps", 0));
    s.pinnedEnd = static_cast<int>(number(kv, "--pinned-end", s.pinnedEnd));   // controls held at each yarn end
    // --settle-failed-construction 1: the diagnostic opt-in, a known bad crossing order settled anyway for research (the
    // default refuses it; a broken strand or port is refused either way)
    if (number(kv, "--settle-failed-construction", 0) != 0) s.failedConstruction = FailedConstruction::SettleAnyway;
    return s;
}

// Sampled geometry of chosen samples (the template calibration's interior comparison): for each masked sample, the nearest
// sample of any strand that is not within 3 D of it along its own strand, and the circumradius of the sample triplet.
struct SampleStats {
    int samples = 0;
    double minClearance = 1e300, pctBelowD = 0, pctBelowDc = 0, pctBelow08 = 0;
    double minRadius = 1e300, pctBelowR = 0, pctBelowRc = 0;   // radius below D / 2, below D / 2 + c / 2
    std::vector<double> radii;                                 // per masked interior sample (for the distribution)
};

SampleStats sampleStats(const std::vector<std::vector<Vec3>>& strands, const std::vector<std::vector<char>>& mask, double D, double c) {
    SampleStats r;
    struct S {
        int strand;
        double arc;
        Vec3 p;
    };
    std::vector<S> all;
    std::vector<std::vector<double>> arcs(strands.size());
    for (size_t k = 0; k < strands.size(); ++k) {
        arcs[k].assign(strands[k].size(), 0);
        for (size_t i = 1; i < strands[k].size(); ++i) arcs[k][i] = arcs[k][i - 1] + norm(strands[k][i] - strands[k][i - 1]);
        for (size_t i = 0; i < strands[k].size(); ++i) all.push_back({static_cast<int>(k), arcs[k][i], strands[k][i]});
    }
    auto cell = [&](double x) { return static_cast<long long>(std::floor(x / D)); };
    auto key = [](long long x, long long y, long long z) { return (x * 73856093LL) ^ (y * 19349663LL) ^ (z * 83492791LL); };
    std::unordered_map<long long, std::vector<int>> grid;
    for (size_t i = 0; i < all.size(); ++i) grid[key(cell(all[i].p.x), cell(all[i].p.y), cell(all[i].p.z))].push_back(static_cast<int>(i));
    int below = 0, belowDc = 0, below08 = 0, kinks = 0, kinksRc = 0, nr = 0;
    for (size_t k = 0; k < strands.size(); ++k)
        for (size_t i = 0; i < strands[k].size(); ++i) {
            if (!mask[k][i]) continue;
            const Vec3 p = strands[k][i];
            double best = 1e300;
            for (long long dx = -1; dx <= 1; ++dx)
                for (long long dy = -1; dy <= 1; ++dy)
                    for (long long dz = -1; dz <= 1; ++dz) {
                        const auto it = grid.find(key(cell(p.x) + dx, cell(p.y) + dy, cell(p.z) + dz));
                        if (it == grid.end()) continue;
                        for (int j : it->second) {
                            if (all[j].strand == static_cast<int>(k) && std::fabs(all[j].arc - arcs[k][i]) < 3 * D) continue;
                            best = std::min(best, norm(all[j].p - p));
                        }
                    }
            ++r.samples;
            r.minClearance = std::min(r.minClearance, best);
            below += best < D;
            belowDc += best < D + c;
            below08 += best < 0.8 * D;
            if (i == 0 || i + 1 >= strands[k].size()) continue;
            const Vec3 a = p - strands[k][i - 1], b = strands[k][i + 1] - p;
            const double la = norm(a), lb = norm(b), lc = norm(a + b), cr = norm(cross(a, b));
            if (la < 1e-12 || lb < 1e-12) continue;
            const double rad = cr > 1e-15 ? la * lb * lc / (2 * cr) : 1e300;
            r.radii.push_back(rad);
            r.minRadius = std::min(r.minRadius, rad);
            kinks += rad < 0.5 * D;
            kinksRc += rad < 0.5 * D + 0.5 * c;
            ++nr;
        }
    if (r.samples) {
        r.pctBelowD = 100.0 * below / r.samples;
        r.pctBelowDc = 100.0 * belowDc / r.samples;
        r.pctBelow08 = 100.0 * below08 / r.samples;
    }
    if (nr) {
        r.pctBelowR = 100.0 * kinks / nr;
        r.pctBelowRc = 100.0 * kinksRc / nr;
    }
    return r;
}

double polylineLength(const std::vector<Vec3>& p) {
    double L = 0;
    for (size_t i = 1; i < p.size(); ++i) L += norm(p[i] - p[i - 1]);
    return L;
}

// The interior of a stitch-mesh fixture: its faces' samples, and per face one repeat of yarn = its row pieces (the new loop's
// legs) plus the head their T ports continue into (the same loop's top, in the face above) - one loop, as production's
// stitch span from sinker to sinker is one loop.
struct InteriorStats {
    int repeats = 0;
    double lengthMean = 0, lengthMin = 1e300, lengthMax = 0;
    SampleStats samples;
};

InteriorStats interiorStats(const YarnEmbedding& e, const std::set<int>& interior, double D, double c) {
    InteriorStats r;
    std::map<YarnPort, int> headByStart;
    for (size_t i = 0; i < e.pieces.size(); ++i)
        if (e.pieces[i].kind == YarnPieceKind::Head) headByStart[e.pieces[i].start] = static_cast<int>(i);
    double sum = 0;
    for (int f : interior) {
        double L = 0;
        int row0 = -1;
        for (size_t i = 0; i < e.pieces.size(); ++i) {
            const YarnPiece& p = e.pieces[i];
            if (p.face != f || p.kind != YarnPieceKind::Row) continue;
            L += polylineLength(p.samples);
            if (p.index == 0) row0 = static_cast<int>(i);
        }
        if (row0 < 0) continue;
        const auto it = headByStart.find(e.pieces[row0].end);
        if (it == headByStart.end()) continue;
        L += polylineLength(e.pieces[it->second].samples);
        ++r.repeats;
        sum += L;
        r.lengthMin = std::min(r.lengthMin, L);
        r.lengthMax = std::max(r.lengthMax, L);
    }
    if (r.repeats) r.lengthMean = sum / r.repeats;
    std::vector<std::vector<Vec3>> strands;
    std::vector<std::vector<char>> mask;
    for (const YarnStrand& sd : e.strands) {
        strands.push_back(sd.samples);
        std::vector<char> m(sd.samples.size(), 0);
        for (size_t i = 0; i < sd.samples.size(); ++i) m[i] = interior.count(sd.sampleFace[i]) ? 1 : 0;
        mask.push_back(m);
    }
    r.samples = sampleStats(strands, mask, D, c);
    return r;
}

std::string interiorLine(const char* state, const InteriorStats& s) {
    char b[600];
    std::snprintf(b, sizeof b,
                  "interior %s: %d repeats, yarn per repeat %.3f mm (%.3f - %.3f); %d samples: nearest other yarn min %.3f mm, closer than D %.1f %%, than D + c %.1f %%, "
                  "than 0.8 D %.1f %%; bend radius min %.3f mm, below D/2 %.1f %%, below r + c/2 %.1f %%",
                  state, s.repeats, s.lengthMean, s.repeats ? s.lengthMin : 0.0, s.lengthMax, s.samples.samples, s.samples.minClearance < 1e299 ? s.samples.minClearance : -1.0,
                  s.samples.pctBelowD, s.samples.pctBelowDc, s.samples.pctBelow08, s.samples.minRadius < 1e299 ? s.samples.minRadius : -1.0, s.samples.pctBelowR, s.samples.pctBelowRc);
    return b;
}

void writeRadii(const std::string& path, const std::vector<double>& radii) {
    std::ofstream f(path, std::ios::binary);
    f << "radius_mm" << '\n';
    for (double r : radii) putf(f, "%.5f\n", std::min(r, 99.0));
}

struct CaseResult {
    std::string name, verdict;      // verdict: settled / refused (<why>) / failed (<why>)
    SettleResult s;                 // the settle (stitchmesh/yarn_settle.h)
    InteriorStats interiorRaw, interiorConverted, interiorSettled;   // with an interior (template calibration)
};

void writeStages(const std::string& path, const std::vector<SettleStage>& rows) {
    std::ofstream f(path, std::ios::binary);
    f << "stage,D_prev,D,cap_mm,ms,valid_at_D,regular,no_fold,rho_lo_mm,clearance,min_distance_mm,fail_regions,unresolved_regions,start_valid,certified,delta_mm,local_scale_mm,solver\n";
    for (const SettleStage& r : rows)
        putf(f, "%d,%.4f,%.4f,%.4f,%.0f,%d,%s,%s,%.4f,%s,%.4f,%d,%d,%d,%d,%.4f,%.4f,\"%s\"\n", r.k, r.Dprev, r.D, r.cap, r.ms, r.gate.valid() ? 1 : 0, toString(r.gate.regular),
             toString(r.gate.fold), r.gate.rhoLo, toString(r.gate.clearance), r.gate.minDistance, r.gate.failRegions, r.gate.unresolvedRegions, r.startValid ? 1 : 0,
             r.certified ? 1 : 0, r.delta, r.localScale, r.solver.c_str());
}

// One case: raw -> converted -> settled (settleYarn), written to `dir`. flipFace: the converted curves' first head of that face
// mirrored through the reference before the solve (a deliberate pass-through, a control of the checks).
CaseResult settleCase(const std::string& name, const TypedStitchMesh& s, const std::vector<std::vector<int>>& faces, const TriMesh& ref, double w, double h,
                      const YarnSettings& st, SettleSettings ss, KnitNumerics& num, const std::string& dir, std::ostream& log, int flipFace = -1,
                      const std::set<int>* interior = nullptr) {
    CaseResult R;
    R.name = name;
    fs::create_directories(dir);
    std::ofstream rep(dir + "/settle_report.txt", std::ios::binary);
    putf(rep, "%s\n", name.c_str());
    YarnEmbedding e;
    std::string err;
    if (!embedYarn(s, faces, ref, w, h, st, e, &err)) {
        R.verdict = "failed (no yarn: " + err + ")";
        putf(rep, "%s\n", R.verdict.c_str());
        return R;
    }
    ss.mirrorHeadOfFace = flipFace;
    SettleHooks hooks;
    hooks.lastSolverCall = [&num]() {
        const std::vector<std::string> calls = gpu::ncpNumericsCalls(num);
        return calls.empty() ? std::string() : calls.back();
    };
    R.s = settleYarn({e, s, faces, ref, st, name}, ss, num, hooks, rep, log);
    R.verdict = R.s.verdict;
    lab::writeYarnFiles(dir + "/raw", s, e, R.s.raw, st, name + ": raw (first-yarn embedding)");
    if (interior) {
        R.interiorRaw = interiorStats(e, *interior, st.diameter, ss.clearance);
        putf(rep, "%s\n", interiorLine("raw", R.interiorRaw).c_str());
        writeRadii(dir + "/raw/interior_radii.csv", R.interiorRaw.samples.radii);
    }
    if (R.s.end == SettleEnd::RawConstructionFailed || R.s.end == SettleEnd::YarnDisconnected || R.s.end == SettleEnd::ConversionRefused) return R;
    lab::writeYarnFiles(dir + "/converted", s, R.s.convertedEmbedding, R.s.converted, st, name + ": converted to uniform cubic B-splines, before the solve");
    if (interior) {
        R.interiorConverted = interiorStats(R.s.convertedEmbedding, *interior, st.diameter, ss.clearance);
        putf(rep, "%s\n", interiorLine("converted", R.interiorConverted).c_str());
        writeRadii(dir + "/converted/interior_radii.csv", R.interiorConverted.samples.radii);
    }
    if (R.s.end != SettleEnd::Solved) return R;
    writeStages(dir + "/settle_stages.csv", R.s.stages);
    {
        // The final gate's clearance regions and curvature, as the certified checks judged them (what fails
        // and why: distance, threshold, the check's reason).
        std::ofstream f(dir + "/final_gate.txt", std::ios::binary);
        const ValidationReport& v = R.s.gateFinal.report;
        putf(f, "%s\nclearance: %zu regions, %zu issues, allowance deltaD %.3g mm, min inter-yarn %.4f, min nonlocal self %.4f mm\n", R.s.gateFinal.line().c_str(),
             v.clearance.regions.size(), v.clearance.issues.size(), v.clearance.deltaD, v.clearance.minInterYarn, v.clearance.minSelfNonlocal);
        for (const ClearanceRegion& g : v.clearance.regions)
            putf(f, "  %s: yarn %d u[%.4f, %.4f] vs yarn %d u[%.4f, %.4f], worst distance %.5f vs threshold %.5f%s: %s\n", toString(g.status), g.yarnA, g.uA0, g.uA1,
                 g.yarnB, g.uB0, g.uB1, g.worstDistance, g.threshold, g.penetrating ? " (penetrating)" : "",
                 g.worstIssue >= 0 && g.worstIssue < static_cast<int>(v.clearance.issues.size()) ? v.clearance.issues[g.worstIssue].reason.c_str() : "");
        for (const CurvatureResult& c : v.curvature)
            putf(f, "curvature yarn %d: rho_lo %.5f (sampled minimum %.5f at u %.4f), regular %s, no fold %s\n", c.yarnId, c.rhoLo, c.rhoEst, c.uRhoEst, toString(c.regular),
                 toString(c.noFold));
    }
    lab::writeYarnFiles(dir + "/settled", s, R.s.settledEmbedding, R.s.settled, st, name + ": after the solve");
    if (interior) {
        R.interiorSettled = interiorStats(R.s.settledEmbedding, *interior, st.diameter, ss.clearance);
        putf(rep, "%s\n", interiorLine("settled", R.interiorSettled).c_str());
        writeRadii(dir + "/settled/interior_radii.csv", R.interiorSettled.samples.radii);
    }
    return R;
}

void summaryHeader(std::ostream& f) {
    f << "case,verdict,raw_status,raw_min_clearance_mm,raw_contact_pct,raw_severe,raw_min_bend_mm,raw_kinks,conv_spacing_mm,conv_max_dev_mm,conv_p99_dev_mm,conv_length_pct,"
         "conv_status,conv_min_clearance_mm,conv_contact_pct,conv_kinks,D0,start_verified,stages,certified_stages,chain_certified,solve_ms,final_valid,final_min_distance_mm,"
         "final_rho_lo_mm,final_fail_regions,settled_status,settled_faces_ok,settled_faces_own_plane_ok,faces,settled_interlocks,interlocks,faces_failing_both_frames,settled_min_clearance_mm,settled_contact_pct,"
         "settled_overlap,settled_severe,settled_min_bend_mm,settled_kinks,max_point_move_mm,p95_point_move_mm,length_change_pct,strain_max,mean_abs_offset_mm,"
         "refused_calls,host_peak_working_set_mib,host_peak_private_mib\n";
}

void summaryRow(std::ostream& f, const CaseResult& r) {
    const SettleResult& t = r.s;
    const YarnCheck &a = t.raw, &b = t.converted, &s = t.settled;
    auto pct = [](const YarnCheck& c) { return c.samples ? 100.0 * c.contactSamples / c.samples : 0.0; };
    auto mc = [](const YarnCheck& c) { return c.minClearance < 1e299 ? c.minClearance : -1.0; };
    auto mr = [](const YarnCheck& c) { return c.minRadius < 1e299 ? c.minRadius : -1.0; };
    const HostPeak host = hostPeak();
    putf(f, "%s,\"%s\",%s,%.4f,%.2f,%d,%.4f,%d,%.4f,%.4f,%.4f,%.4f,%s,%.4f,%.2f,%d,%.4f,%d,%d,%d,%d,%.0f,%d,%.4f,%.4f,%d,%s,%d,%d,%d,%d,%d,%d,%.4f,%.2f,%d,%d,%.4f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%.1f,%.1f\n",
         r.name.c_str(), r.verdict.c_str(), toString(a.status), mc(a), pct(a), a.severeSamples, mr(a), a.kinkSamples, t.conversion.spacing.empty() ? 0.0 : t.conversion.spacing[0],
         t.conversion.maxDeviation, t.conversion.p99Deviation, t.conversion.hermiteLength > 0 ? 100.0 * (t.conversion.curveLength / t.conversion.hermiteLength - 1) : 0.0,
         toString(b.status), mc(b), pct(b), b.kinkSamples, t.D0, t.startVerified ? 1 : 0, static_cast<int>(t.stages.size()), t.certifiedStages, t.chainCertified ? 1 : 0, t.solveMs, t.gateFinal.valid() ? 1 : 0,
         t.gateFinal.minDistance, t.gateFinal.rhoLo, t.gateFinal.failRegions, toString(s.status), s.facesCrossingsOk, s.facesOwnPlaneOk, s.facesChecked, s.interlocksVerified, s.interlocksExpected,
         t.letterSwapped, mc(s), pct(s), s.overlapSamples, s.severeSamples, mr(s), s.kinkSamples, t.drift.maxPointMove, t.drift.p95PointMove,
         t.drift.lengthRest > 0 ? 100.0 * (t.drift.length / t.drift.lengthRest - 1) : 0.0, t.drift.strainMax, t.drift.meanAbsDepth, t.refusedCalls, host.workingSetMiB,
         host.privateMiB);
}

// --contact-storage compact|live (ContactStorage); without it the settle's own (SettleSettings).
ContactStorage contactStorageFrom(const std::map<std::string, std::string>& kv) {
    if (!kv.count("--contact-storage")) return SettleSettings{}.contactStorage;
    return kv.at("--contact-storage") == "live" ? ContactStorage::LiveRows : ContactStorage::CompactRows;
}

// --curvature-rows samples|span-max (CurvatureRows); without it the settle's own (SettleSettings).
CurvatureRows curvatureRowsFrom(const std::map<std::string, std::string>& kv) {
    if (!kv.count("--curvature-rows")) return SettleSettings{}.curvatureRows;
    return kv.at("--curvature-rows") == "span-max" ? CurvatureRows::SamplesAndSpanMaximum : CurvatureRows::Samples;
}

// The solver identity of the commands' solver options.
std::string solverIdentityFrom(const std::map<std::string, std::string>& kv) { return gpu::ncpSolverIdentity(contactStorageFrom(kv), curvatureRowsFrom(kv)); }

// The device numerics with the commands' solver options: --gpu-budget-mib (the cap), --capture-calls (the call files),
// --contact-storage (the solver's contact storage), --curvature-rows (its curvature rows).
std::unique_ptr<KnitNumerics> deviceNumerics(const std::map<std::string, std::string>& kv) {
    std::string err;
    gpu::NcpNumericsOptions options;
    options.contactStorage = contactStorageFrom(kv);
    options.curvatureRows = curvatureRowsFrom(kv);
    // --slow-call stop|budget (SlowCall): whether a slowly converging call is ended early or runs to its step budget
    if (kv.count("--slow-call") && kv.at("--slow-call") == "budget") options.slowCall = SlowCall::UseStepBudget;
    options.memoryCapMB = number(kv, "--gpu-budget-mib", 0);
    if (kv.count("--capture-calls")) {
        options.captureDir = kv.at("--capture-calls");
        fs::create_directories(options.captureDir);
    }
    std::unique_ptr<KnitNumerics> num = gpu::makeKnitNumerics(&err, options);
    if (!num) std::fprintf(stderr, "yarn-settle: no device numerics: %s\n", err.c_str());
    return num;
}

void writeSettings(const std::string& path, const SettleSettings& ss, const YarnSettings& st, const std::string& solverIdentity) {
    std::ofstream f(path, std::ios::binary);
    putf(f, "templates %s; yarn diameter %.3f mm, clearance %.3f mm; conversion spacing %.4f mm, %d dense points per Hermite span; soft contact D + c + %.3f, soft bend factor %.2f, length "
            "stiffness %.2f, band %s %.2f mm, %d samples per span, %d iterations, refresh %d; cap %.2f (D_prev + c) per call, D0 min %.3f, growth %.2f, at most %d calls at the target (stop rule: gate valid and construction kept, or stalled), "
            "%d controls pinned per yarn end, budget %.0f min of solve time; a failed construction: %s; solver %s\n",
         st.templateSet.c_str(), st.diameter, ss.clearance, ss.spacing, ss.densePerSpan, ss.contactMargin, ss.bendFactor, ss.lengthStiffness, ss.band ? "on" : "off", ss.bandHalf, ss.samplesPerSpan,
         ss.iterations, ss.refresh, ss.displacementFraction, ss.D0Min, ss.growth, ss.maxRepeats, ss.pinnedEnd, ss.budgetMinutes, toString(ss.failedConstruction), solverIdentity.c_str());
}

int commandSettleFixtures(const std::map<std::string, std::string>& kv) {
    const std::string out = kv.count("--out") ? kv.at("--out") : "";
    if (out.empty()) return 1;
    fs::create_directories(out);
    const double w = 5.5, h = 3.4;
    const YarnSettings base = lab::yarnSettingsFrom(kv);
    const SettleSettings ss = settleFrom(kv);
    std::unique_ptr<KnitNumerics> num = deviceNumerics(kv);
    if (!num) return 1;
    writeSettings(out + "/settle_settings.txt", ss, base, solverIdentityFrom(kv));
    std::ofstream csv(out + "/settle_fixtures.csv", std::ios::binary);
    summaryHeader(csv);
    std::ofstream log(out + "/settle_log.txt", std::ios::binary);
    int refusedCalls = 0;
    std::vector<lab::YarnFixture> cases = lab::yarnFixtureCases();
    // the motion control: K_patch with one head mirrored through the reference after conversion, before the solve
    lab::YarnFixture motion = cases.front();
    motion.name = "neg_head_pushed_through";
    motion.purpose = "negative control: K_patch converted, then one head (strip 1, slot 4) mirrored through the surface: the face crossings (both frames) and the certificate of the mirroring must flag it";
    motion.negative = true;
    motion.flipHeadsOfSlot = 4;
    cases.push_back(motion);
    std::set<std::string> only;
    if (kv.count("--cases")) {
        std::string list = kv.at("--cases"), item;
        for (char ch : list + ",") {
            if (ch == ',') {
                if (!item.empty()) only.insert(item);
                item.clear();
            } else item += ch;
        }
    }
    for (const lab::YarnFixture& fx : cases) {
        if (!only.empty() && !only.count(fx.name)) continue;
        TypedStitchMesh m;
        std::vector<std::vector<int>> faces;
        TriMesh ref;
        RelaxProduct relax;
        std::string err;
        if (!lab::buildYarnFixture(fx, w, h, m, faces, ref, relax, &err)) {
            putf(log, "%s: not built: %s\n", fx.name.c_str(), err.c_str());
            continue;
        }
        YarnSettings st = base;
        int flipFace = -1;
        for (int f = 0; f < m.mesh.faces(); ++f) {
            if (!m.strip.empty() && m.strip[f] == 1 && m.slot[f] == fx.flipHeadsOfSlot) (fx.name == "neg_head_pushed_through" ? flipFace : st.flipHeadsOfFace) = f;
            if (!m.strip.empty() && m.strip[f] == 1 && m.slot[f] == fx.skipSlot) st.skipFace = f;
        }
        putf(log, "== %s: %s\n", fx.name.c_str(), fx.purpose.c_str());
        // the interior of the regular patch: strips two rows from either boundary
        std::set<int> interior;
        int maxStrip = 0;
        for (int f = 0; f < m.mesh.faces(); ++f) maxStrip = std::max(maxStrip, m.strip.empty() ? 0 : m.strip[f]);
        if (maxStrip >= 4)
            for (int f = 0; f < m.mesh.faces(); ++f)
                if (m.strip[f] >= 2 && m.strip[f] <= maxStrip - 2) interior.insert(f);
        const CaseResult r = settleCase(fx.name, m, faces, ref, w, h, st, ss, *num, out + "/" + fx.name, log, flipFace, interior.empty() ? nullptr : &interior);
        lab::writePolygonObj(out + "/" + fx.name + "/stitch_mesh.obj", m);
        lab::writeTriangleObj(out + "/" + fx.name + "/reference.obj", ref);
        summaryRow(csv, r);
        if (!interior.empty()) {
            std::ofstream ic(out + "/" + fx.name + "/interior.csv", std::ios::binary);
            ic << "state,repeats,length_mean_mm,length_min_mm,length_max_mm,samples,min_clearance_mm,pct_below_D,pct_below_Dc,pct_below_08D,min_radius_mm,pct_below_r,pct_below_rc" << '\n';
            for (const auto& [name, x] : {std::pair<const char*, const InteriorStats*>{"raw", &r.interiorRaw}, {"converted", &r.interiorConverted}, {"settled", &r.interiorSettled}})
                putf(ic, "%s,%d,%.4f,%.4f,%.4f,%d,%.4f,%.2f,%.2f,%.2f,%.4f,%.2f,%.2f\n", name, x->repeats, x->lengthMean, x->lengthMin, x->lengthMax, x->samples.samples,
                     x->samples.minClearance < 1e299 ? x->samples.minClearance : -1.0, x->samples.pctBelowD, x->samples.pctBelowDc, x->samples.pctBelow08,
                     x->samples.minRadius < 1e299 ? x->samples.minRadius : -1.0, x->samples.pctBelowR, x->samples.pctBelowRc);
        }
        putf(log, "   %s\n", r.verdict.c_str());
        refusedCalls += r.s.refusedCalls;
        log.flush();
        std::printf("%s: %s\n", fx.name.c_str(), r.verdict.c_str());
    }
    putf(log, "%s\n", gpu::knitNumericsStats(*num).c_str());
    for (const std::string& l : gpu::ncpNumericsCalls(*num)) putf(log, "  %s\n", l.c_str());
    return refusedCalls ? 3 : 0;
}

int commandSettle(const std::map<std::string, std::string>& kv) {
    const std::string caseDir = kv.count("--case") ? kv.at("--case") : "", which = kv.count("--which") ? kv.at("--which") : "after";
    const std::string out = kv.count("--out") ? kv.at("--out") : "";
    if (caseDir.empty() || out.empty()) return 1;
    const double w = number(kv, "--w", 5.5), h = number(kv, "--h", 3.4);
    TypedStitchMesh s;
    std::string saved;
    if (!lab::readTypedMeshFiles(caseDir + "/" + which, s, saved)) {
        std::fprintf(stderr, "yarn-settle: cannot read the typed mesh in %s/%s\n", caseDir.c_str(), which.c_str());
        return 1;
    }
    TriMesh ref;
    if (!lab::readReferenceObj(caseDir + "/reference.obj", ref)) {
        std::fprintf(stderr, "yarn-settle: no reference.obj in %s\n", caseDir.c_str());
        return 1;
    }
    char digest[32];
    std::snprintf(digest, sizeof digest, "%016llx", static_cast<unsigned long long>(typedDigest(s)));
    VirtualShift vs;
    std::string err;
    if (!buildVirtualShift(s, ShiftReading::PathSideDown, vs, &err)) {
        std::fprintf(stderr, "yarn-settle: %s\n", err.c_str());
        return 1;
    }
    const YarnSettings st = lab::yarnSettingsFrom(kv);
    const SettleSettings ss = settleFrom(kv);
    std::unique_ptr<KnitNumerics> num = deviceNumerics(kv);
    if (!num) return 1;
    fs::create_directories(out);
    writeSettings(out + "/settle_settings.txt", ss, st, solverIdentityFrom(kv));
    std::ofstream log(out + "/settle_log.txt", std::ios::binary);
    putf(log, "frozen mesh %s/%s, typed digest as read %s\n", caseDir.c_str(), which.c_str(), digest);
    const CaseResult r = settleCase("arm " + caseDir, s, vs.corners, ref, w, h, st, ss, *num, out, log);
    lab::writePolygonObj(out + "/stitch_mesh.obj", s);
    std::ofstream csv(out + "/settle_summary.csv", std::ios::binary);
    summaryHeader(csv);
    summaryRow(csv, r);
    putf(log, "%s\n%s\n", r.verdict.c_str(), gpu::knitNumericsStats(*num).c_str());
    for (const std::string& l : gpu::ncpNumericsCalls(*num)) putf(log, "  %s\n", l.c_str());
    const HostPeak host = hostPeak();
    putf(log, "host memory peak: working set %.1f MiB, private %.1f MiB\n", host.workingSetMiB, host.privateMiB);
    std::printf("yarn-settle: %s\n", r.verdict.c_str());
    return r.s.refusedCalls ? 3 : 0;
}

// One captured relax call again (gpu::replayNcpCall): --call <file> [--mode solve|count|frozen] [--gpu-budget-mib N] [--rho R] [--out <dir>].
// Writes <out>/replay.txt (status, call line, solver statistics, memory report, controls digest, host peak) and, when
// solved, <out>/controls.bin. Exit 0 solved or counted, 3 refused by the cap, 1 failed.
int commandReplay(const std::map<std::string, std::string>& kv) {
    gpu::NcpReplayRequest rq;
    rq.callFile = kv.count("--call") ? kv.at("--call") : "";
    if (rq.callFile.empty()) return 1;
    const std::string mode = kv.count("--mode") ? kv.at("--mode") : "solve";
    rq.mode = mode == "count" ? gpu::NcpReplayMode::CountOnly : mode == "frozen" ? gpu::NcpReplayMode::Frozen : gpu::NcpReplayMode::Solve;
    rq.frozenRho = number(kv, "--rho", 10);
    rq.maxLmSteps = static_cast<long long>(number(kv, "--max-lm-steps", 0));
    rq.memoryCapMB = number(kv, "--gpu-budget-mib", 0);
    rq.contactStorage = contactStorageFrom(kv);
    rq.curvatureRows = curvatureRowsFrom(kv);
    const std::string out = kv.count("--out") ? kv.at("--out") : "";
    if (!out.empty()) {
        fs::create_directories(out);
        if (rq.mode == gpu::NcpReplayMode::Solve) rq.controlsOut = out + "/controls.bin";
        if (rq.mode == gpu::NcpReplayMode::Frozen) rq.frozenOut = out + "/frozen.bin";
    }
    gpu::NcpReplayResult r;
    std::string err;
    if (!gpu::replayNcpCall(rq, r, &err)) {
        std::fprintf(stderr, "ncp-replay: %s\n", err.c_str());
        return 1;
    }
    const HostPeak host = hostPeak();
    char tail[200];
    std::snprintf(tail, sizeof tail, "controls digest %016llx; cap %.0f MiB; %.0f ms; host memory peak: working set %.1f MiB, private %.1f MiB",
                  r.controlsDigest, r.capMB, r.ms, host.workingSetMiB, host.privateMiB);
    const std::string text = "call " + rq.callFile + " (" + mode + ")\nstatus: " + (r.refused ? "REFUSED: " : "") + r.status + "\n" + r.callLine + "\n" + r.statsLine +
                             "\n" + r.memoryLine + "\n" + tail + "\n";
    std::printf("%s", text.c_str());
    if (!out.empty()) {
        std::ofstream f(out + "/replay.txt", std::ios::binary);
        f << text;
    }
    return r.refused ? 3 : 0;
}

// Production's stockinette (knitting/stitch/stockinette.h) on a flat n x n patch at the same gauge and diameter: the positive
// control of the template calibration. Its interior stitches (two from every side) measured as the stitch-mesh interior is,
// before and after production's settle pass (knitting/relax/settle.h, SettleParams defaults, no surface band: a flat patch).
int commandProductionPatch(const std::map<std::string, std::string>& kv) {
    const std::string out = kv.count("--out") ? kv.at("--out") : "";
    if (out.empty()) return 1;
    fs::create_directories(out);
    PatternParams pp;
    pp.wales = pp.courses = static_cast<int>(number(kv, "--n", 8));
    pp.waleSpacing = number(kv, "--w", 5.5);
    pp.courseSpacing = number(kv, "--h", 3.4);
    pp.yarnDiameter = number(kv, "--diameter", 1.16);
    pp.clearance = number(kv, "--c", 0.1);
    CurveSet cs = generateStockinette(pp);
    std::unique_ptr<KnitNumerics> num = deviceNumerics(kv);
    if (!num) return 1;
    std::ofstream rep(out + "/production_patch.txt", std::ios::binary);
    putf(rep, "production stockinette %d x %d, w %.2f h %.2f mm, D %.2f mm, c %.2f mm, template A %.2f D, head %.3f w, leg %.3f w (PatternParams defaults)\n", pp.wales, pp.courses,
         pp.waleSpacing, pp.courseSpacing, pp.yarnDiameter, pp.clearance, pp.depthAmplitude, pp.headHalfWidth, pp.legHalfWidth);
    const int lo = 2, hi = pp.wales - 3;
    auto measure = [&](const CurveSet& c, const char* state) {
        // samples: 8 per span; mask: inside an interior stitch's span; repeat: an interior stitch's arc length
        std::vector<std::vector<Vec3>> strands;
        std::vector<std::vector<char>> mask;
        InteriorStats st;
        double sum = 0;
        for (const YarnPath& y : c.yarns) {
            std::vector<Vec3> pts;
            std::vector<char> m;
            const int per = 8;
            for (int k = 0; k <= per * y.curve.numSpans(); ++k) {
                const double u = static_cast<double>(k) / per;
                pts.push_back(y.positionAtParam(u));
                bool in = false;
                for (const StitchSpan& sp : y.stitches)
                    if (sp.wale >= lo && sp.wale <= hi && sp.course >= lo && sp.course <= hi && u >= sp.uBegin && u < sp.uEnd) in = true;
                m.push_back(in ? 1 : 0);
            }
            for (const StitchSpan& sp : y.stitches)
                if (sp.wale >= lo && sp.wale <= hi && sp.course >= lo && sp.course <= hi) {
                    const double L = y.arcLengthAt(sp.uEnd) - y.arcLengthAt(sp.uBegin);
                    ++st.repeats;
                    sum += L;
                    st.lengthMin = std::min(st.lengthMin, L);
                    st.lengthMax = std::max(st.lengthMax, L);
                }
            strands.push_back(pts);
            mask.push_back(m);
        }
        if (st.repeats) st.lengthMean = sum / st.repeats;
        st.samples = sampleStats(strands, mask, pp.yarnDiameter, pp.clearance);
        const ValidationReport r = validate(c, pp, {.generatorChecks = false});
        putf(rep, "%s\n  production checks: %s", interiorLine(state, st).c_str(), formatReport(r, c, pp, false).c_str());
        writeRadii(out + "/" + state + "_interior_radii.csv", st.samples.radii);
        // centrelines and tubes (radius D / 2, 8 sides), as the stitch-mesh states are written
        YarnEmbedding view;
        for (const std::vector<Vec3>& p : strands) {
            YarnStrand sd;
            sd.samples = p;
            view.strands.push_back(sd);
        }
        lab::writeTriangleObj(out + "/" + state + "_tubes.obj", yarnTubes(view, 0.5 * pp.yarnDiameter, 8));
        std::ofstream f(out + "/" + state + "_centrelines.obj", std::ios::binary);
        int base = 1;
        for (const std::vector<Vec3>& p : strands) {
            for (const Vec3& q : p) putf(f, "v %.6f %.6f %.6f\n", q.x, q.y, q.z);
            f << "l";
            for (size_t i = 0; i < p.size(); ++i) putf(f, " %d", base + static_cast<int>(i));
            f << '\n';
            base += static_cast<int>(p.size());
        }
        // the template's controls of one interior stitch (for the control diagram)
        std::ofstream g(out + "/" + state + "_controls.csv", std::ios::binary);
        g << "yarn,control,x,y,z" << '\n';
        for (size_t yi = 0; yi < c.yarns.size(); ++yi)
            for (size_t k = 0; k < c.yarns[yi].curve.ctrl.size(); ++k)
                putf(g, "%zu,%zu,%.6f,%.6f,%.6f\n", yi, k, c.yarns[yi].curve.ctrl[k].x, c.yarns[yi].curve.ctrl[k].y, c.yarns[yi].curve.ctrl[k].z);
        return st;
    };
    const InteriorStats before = measure(cs, "placed");
    SettleParams sp;
    sp.exec.parallel = true;
    ValidationReport last;
    const auto t0 = Clock::now();
    rep.flush();
    std::fprintf(stderr, "production-patch: settling\n");
    // production settles surface knits (pieces); the device path does not take a flat patch's global crossing check, so the
    // patch is settled as a surface knit (the same curves and topology, the per-stitch local signature) and judged afterwards
    // with the flat checks again
    cs.kind = CurveKind::Surface;
    const SettleReport sr = settle({cs, pp, nullptr, nullptr, nullptr}, sp, *num, &last);
    cs.kind = CurveKind::FlatPatch;
    std::fprintf(stderr, "production-patch: settled\n");
    const double ms = msSince(t0);
    putf(rep, "settle: %d round(s), %d free stitches, max control move %.4f mm, %.0f ms; %s\n", sr.rounds, sr.freeStitches, sr.maxMove, ms,
         gpu::ncpNumericsCalls(*num).empty() ? "" : gpu::ncpNumericsCalls(*num).back().c_str());
    const InteriorStats after = measure(cs, "settled");
    std::ofstream csv(out + "/production_patch.csv", std::ios::binary);
    csv << "state,repeats,length_mean_mm,length_min_mm,length_max_mm,samples,min_clearance_mm,pct_below_D,pct_below_Dc,pct_below_08D,min_radius_mm,pct_below_r,pct_below_rc" << '\n';
    for (const auto& [name, x] : {std::pair<const char*, const InteriorStats*>{"placed", &before}, {"settled", &after}})
        putf(csv, "%s,%d,%.4f,%.4f,%.4f,%d,%.4f,%.2f,%.2f,%.2f,%.4f,%.2f,%.2f\n", name, x->repeats, x->lengthMean, x->lengthMin, x->lengthMax, x->samples.samples,
             x->samples.minClearance, x->samples.pctBelowD, x->samples.pctBelowDc, x->samples.pctBelow08, x->samples.minRadius, x->samples.pctBelowR, x->samples.pctBelowRc);
    std::printf("production-patch: placed %.3f mm per repeat, settled %.3f mm (%+.1f %%); settle %.0f ms\n", before.lengthMean, after.lengthMean,
                100.0 * (after.lengthMean / before.lengthMean - 1), ms);
    return 0;
}

} // namespace

int settleCommand(const std::string& cmd, int argc, char** argv) {
    const std::map<std::string, std::string> kv = parseArgs(argc, argv);
    if (cmd == "yarn-settle-fixtures") return commandSettleFixtures(kv);
    if (cmd == "yarn-settle") return commandSettle(kv);
    if (cmd == "ncp-replay") return commandReplay(kv);
    if (cmd == "production-patch") return commandProductionPatch(kv);
    return 2;
}
