#pragma once
// ncp-gpu: the knitting problem (knit_problem.h) solved with global sparse
// Levenberg-Marquardt steps, the whole nonlinear iteration resident on the device AND decided on the
// device. solveNcp (ncp_solver.cpp) uploads the call, builds one CUDA graph of the iteration
// (ncp_control.cu: conditional WHILE / IF nodes, one-thread controller kernels make every decision),
// launches it once and downloads the result. Between the upload and the download the host reads
// nothing: conjugate gradients, step acceptance, swept-step halving, multipliers, penalty, contact
// refresh, capacity checks and the stopping rules all run on the device.
//
// Converged = feasible (hard families within tolViol, no degenerate curvature sample) AND KKT
// stationarity <= tolStat AND complementarity <= tolViol at the accepted state; otherwise UNRESOLVED
// with the reason.
//
// Rows (double precision, value-consistent with KnitProblem::linearize, tests/ncp_reference.h):
//   contact rows      one per broad-phase candidate pair; hard row (hardC - d) / D and, with
//                     hard thresholds, its target companion max(0, thr - d) / D (same gradient);
//   general rows      edges, curvature samples (hard + target companion), band (two rows, the
//                     surface evaluated AT THE CURRENT CONTROL by continuation from the refresh
//                     state), cap, anchor (three rows); with CurvatureRows::SamplesAndSpanMaximum one more
//                     curvature row per span at its curvature maximum, after all others.
// Merit (augmented Lagrangian, PHR): objective rows (w / 2) r^2; hard rows (rho / 2) phi^2 with
// phi = max(0, g + lambda / rho).
// Linear solve: matrix-free preconditioned conjugate gradients (3x3 block Jacobi) on the damped
// Gauss-Newton system, the product summed per variable over its active rows (ncp_pcg.cu).
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "knitting/numerics.h"
#include "knitting/relax/knit_problem.h"

namespace knit::gpu {

struct NcpOptions {
    // ContactStorage (knit_problem.h): how the contact rows are held on the device.
    ContactStorage contactStorage = ContactStorage::CompactRows;
    // CurvatureRows (knit_problem.h): the fixed samples per span, or those and one row at the span's maximum.
    CurvatureRows curvatureRows = CurvatureRows::Samples;
    // Contact rows that outgrow their capacity in a refresh: the solve keeps its state, grows the
    // capacity within the memory cap and continues exactly where it stopped, at most this often; then
    // (or when the cap leaves no room) it ends UNRESOLVED ("capacity exceeded") or is refused by the cap.
    int maxRegrows = 8;
    // Tests: after every trial (accepted or rejected, and every halving), check on the device that
    // every pair of active sample segments closer than thr at the trial or at the midpoint of the
    // straight motion to it was within the broad phase's pad at the last refresh (a candidate with
    // rows). The count of pairs that were not lands in NcpStats::coverageViolations.
    bool auditCoverage = false;
};

// One relax call.
struct NcpCall {
    const CurveSet& curves;
    Controls& controls;  // the start, and the result
    const Controls& rest;
    const ControlMask& freeControls;
    const BaseSurface* surface;
    SurfaceBinding* binding;  // may be null; read at the start, written with the result
    const RelaxMoveParams& move;
    const KnitProblemParams& problem;
    const NcpSolveParams& solve;
    const NcpOptions& options;
};

// The bytes one group of the solve's device buffers holds at its largest (the memory report).
struct NcpMemoryFamily {
    std::string name;
    long long bytes = 0;
};

struct NcpStats {
    long long unknowns = 0;
    // Sizes for the memory report: the starting state's candidates and the contact rows with a term in
    // the merit there (multipliers 0), and the largest lists of the solve.
    long long initialCandidates = 0, initialActiveContactRows = 0, samples = 0, generalRows = 0;
    long long maxCandidates = 0, maxActiveRows = 0, maxActiveEntries = 0, registrations = 0;
    // The rows that must be resident: contact rows whose sample segments are closer than the near cutoff
    // (max(hardC, thr) + 2 delta) at a broad phase, rows that carried a positive multiplier from beyond
    // it, and live rows (a term in the merit or a positive multiplier) at a current-state evaluation.
    long long initialNearRows = 0, maxNearRows = 0, maxCarriedFar = 0, maxLiveRows = 0;
    long long capacityLiveRows = 0;  // ContactStorage::LiveRows: the live list's capacity (capacityContactRows: the rows')
    std::vector<NcpMemoryFamily> memory;  // per buffer group, at the end of the setup
    double setupPeakMB = 0;               // tracked process peak after the setup
    std::string memoryLine() const;
    long long lmSteps = 0, rejected = 0, rejectedMerit = 0, sweptRejected = 0, sweptUnresolved = 0, ownershipRejected = 0, ownershipPinned = 0,
              regainedTrials = 0, refreshes = 0, maxContactRows = 0, initialContactRows = 0, solves = 0;
    long long pcgIterations = 0, pcgCapped = 0, pcgNonSpdBlocks = 0, pcgBreakdowns = 0;  // breakdown: p.q <= 0 (the iterate is kept)
    long long capacityContactRows = 0;
    long long coverageChecks = 0, coverageViolations = 0;  // NcpOptions::auditCoverage
    double msSetup = 0, msGraph = 0, msIterate = 0, msDownload = 0;
    long long uploadBytes = 0, downloadBytes = 0;
    double devPeakMB = 0;
    long long regrows = 0, launches = 0;  // capacity growths, and graph launches (1 + regrows)
    double msRegrow = 0;
    std::string line() const;
};

// false with `err` when the device fails; an unconverged solve is a result (the report says why).
bool solveNcp(const NcpCall& call, NcpSolveReport& report, NcpStats* stats, std::string* err);
// Where one outer iteration left the solve.
struct OuterProgress {
    bool feasible;
    double objPrev, obj;  // the objective after the previous and the current outer iteration (objPrev < 0: none yet)
    double rho;
    double viol;          // the current maximum violation
    double violBack;      // the maximum violation infeasibleWindow outer iterations ago (< 0: not that many yet)
};
// The early-stop decision after an outer iteration (the controller applies the same rule).
NcpEarlyStop ncpEarlyStop(const NcpSolveParams& sp, const OuterProgress& progress);

// One evaluation's reductions.
struct NcpEvalSummary {
    double merit = 0;          // augmented-Lagrangian merit at the evaluated controls
    double objective = 0;      // objective rows only (edges, targets, anchor)
    double maxViol[knit::kRowFamilies] = {};  // as KnitProblemScore (Edge: max |r|; targets: shortfall)
    long long nViol[knit::kRowFamilies] = {};
    long long degenerate = 0, notOwned = 0;
    long long pinned = 0;        // bound controls pinned (pinLostOwnership): their band rows are off
    double pinnedBandViol = 0;   // the largest band violation among them (a constant of the call)
    double strainMax = 0, strainRms = 0;
    double maxConstraint() const;
};

enum class ControlSet { Current, Trial };
struct SweptResult {
    long long touching = 0;    // pairs whose segments come closer than eps
    long long unresolved = 0;  // pairs the conservative advancement could not clear within its budget
};

// The problem of one relax call resident on the device. solveNcp drives it through its graph; the
// methods below run single operations (each one synchronizes and reads its result), for the tests.
// What the constructor prepares: the whole solve, or only the starting state's sizes (the NcpStats
// counts, no contact buffers), which tells how large a call is without allocating its rows.
enum class NcpSizing { Solve, CountOnly };

class NcpDevice {
public:
    // Uploads one relax call (the inputs of knit::KnitProblem) and sizes the contact capacity from
    // the starting state. The binding is read here; refresh() advances the device copy as
    // KnitProblem does.
    explicit NcpDevice(const NcpCall& call, NcpSizing sizing = NcpSizing::Solve);
    ~NcpDevice();
    NcpDevice(const NcpDevice&) = delete;
    NcpDevice& operator=(const NcpDevice&) = delete;

    int vars() const;
    // Tests: lower the index limit (default 2^31 - 1) to exercise the explicit failure; set before
    // the constructor sizes anything (NcpDevice::State::checkIndexLimit).
    static void setIndexLimitForTests(long long limit);
    // Tests: shrink the contact-row capacity, so the next refresh overflows (the explicit stop).
    void setContactCapacityForTests(long long rows);
    // Tests: the capacity's headroom over the starting rows (factor, minimum) for the devices made
    // after the call; a negative factor restores the default.
    static void setContactHeadroomForTests(double factor, long long minimum);
    // Tests (ContactStorage::LiveRows): the live list's headroom over the starting live rows (factor,
    // minimum) for the devices made after the call; a negative factor restores the declared default.
    static void setLiveHeadroomForTests(double factor, long long minimum);

    // Broad phase at the current controls (pad thr + 2 delta), Xref := P, binding states continued
    // (advanced on Owned), contact rows rebuilt (multipliers carried over by row key).
    void refresh();
    double motionSinceRefresh();  // max over free controls |P - Xref|

    // NcpRows at the current controls (Current) or at the trial controls (Trial): values, gradients and
    // the augmented-Lagrangian terms at penalty rho; the summary of that evaluation.
    NcpEvalSummary linearize(ControlSet which, double rho);
    // Swept (continuous) contact test of the straight motion P -> Q for every candidate pair.
    SweptResult swept(double eps);
    // Tests: the trial controls given directly (all controls; fixed ones must equal P's).
    void setTrial(const knit::Controls& Q);
    // Accept the trial: P := Q and the trial rows become the current rows.
    void accept();

    // Controls and binding back to the host (end of a solve).
    void downloadResult(knit::Controls& P, knit::SurfaceBinding* bind);

    // Tests: every current row as KnitProblem::linearize would emit it (targets as separate rows,
    // degenerate / not-owned rows omitted), keys as ncp_reference.cpp rowKey. Downloads everything.
    struct DebugRow {
        int family = 0;
        uint64_t key = 0;
        double value = 0;
        int n = 0;
        int var[knit::kRowMaxEntries] = {};
        double g[knit::kRowMaxEntries][3] = {};
        bool hasGradient = true;  // false: a contact row without a term in the merit (no gradient evaluated)
    };
    // onlyVariable >= 0: the rows touching that variable only (small device buffers: full-size states);
    // kDebugGeneralRowsOnly: every general row and no contact row.
    static constexpr int kDebugGeneralRowsOnly = -2;
    std::vector<DebugRow> debugRows(int onlyVariable = -1);

    // Diagnostics (ncp-replay --mode frozen): linearizations held fixed, for comparing two layouts or a
    // reference at the same state. Pass 0: refresh, the current rows at rho, the assembly, the undamped
    // operator product H p. Pass 1: the multipliers updated (lambda := max(0, lambda + rho g)), a second
    // refresh at the same controls (carrying them by key), the same again. Downloads everything.
    struct FrozenPass {
        NcpEvalSummary eval;
        std::vector<double> grad, hblock, product;  // 3, 6, 3 per variable
        long long contactRows = 0, activeRows = 0, activeEntries = 0;
    };
    std::vector<FrozenPass> frozenLinearizations(double rho, const std::vector<double>& p);

    long long contactRows() const;
    // Diagnostics: the current contact rows' keys (ncp_reference.cpp rowKey), in row order.
    std::vector<uint64_t> contactKeys();
    // Tests (the omission control): rows [first, first + count) overwritten with row first - 1, so the
    // row count stays and only a key comparison can tell that those rows are gone.
    void overwriteContactRowsForTests(long long first, long long count);
    NcpStats& stats();

    struct State;
    State& state();

private:
    std::unique_ptr<State> state_;
};

} // namespace knit::gpu
