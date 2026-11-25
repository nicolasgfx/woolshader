#pragma once
// ncp-gpu: what the parts of NcpDevice share: rows, parameters, buffers, the control record (ncp.h).
// Double precision for the controls, rows, products and reductions (the solver's decisions compare
// merits to ~1e-12 relative; the RTX 40xx FP64 rate is ample for these per-row kernels); the broad
// phase uses float boxes with a guard added to the pad, and the surface continuation and projection
// are the double twins in surface.cuh.
//
// DEVICE-OWNED ITERATION: every size that changes during a solve (contact rows, candidates, active
// rows) lives in NcpCounts on the device, and every decision lives in NcpControl, written by
// one-thread controller kernels (ncp_control.cu). Kernels over a changing list run a fixed grid and
// loop up to the device count, so a captured launch serves every iteration. The host reads nothing
// between the solve's start and its end.
//
// DETERMINISM: no floating-point atomics. Sums are two-pass block reductions with a fixed grid;
// maxima of non-negative doubles are integer atomicMax on their bit patterns (order independent).
// Lists filled in arbitrary order are either summed per row in place, or sorted per variable before
// anything is summed over them (activeLists), so every sum adds the same numbers in the same order.
//
// CAPACITY: buffers are sized once per call from the starting state (host preparation) with
// headroom. A list that would exceed its capacity is not truncated: the refresh records the
// overflow and the iteration ends at the last accepted state with an explicit UNRESOLVED status.
//
// MEMORY: a contact row stores only what cannot be derived
// again: its segment pair, type, multiplier, values and the current terms (about 50 bytes). Its
// controls and key follow from the pair (contactControls, contactKey); its gradient is formed for
// the rows with a term in the merit only, in the active lists (per slot: the normal and the closest
// points' parameters; the coefficient per control from them, contactCoefficient). The broad phase
// writes the rows directly (no candidate list), and the multiplier carry-over sorts the previous
// keys in the active entries' memory, which is free between a refresh and the next assembly.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include <cuda_runtime.h>

#include "gpu/device/kernel.cuh"
#include "gpu/device/primitives.h"
#include "gpu/geometry/device_geometry.h"
#include "gpu/geometry/surface.cuh"
#include "gpu/solver/ncp.h"

namespace knit::gpu {

enum RowFamilyCode : int { FEdge = 0, FCurv = 1, FContact = 2, FBandU = 3, FBandL = 4, FCap = 5, FContactT = 6, FCurvT = 7, FAnchor = 8 };
constexpr int kContactRowEntries = knit::kRowMaxEntries;  // contact row entries
constexpr int kGeneralRowEntries = 4;                      // general row entries

// Kernels over a list whose length lives on the device: this many blocks, each thread striding.
constexpr int kLoopBlocks = 512;
__device__ __forceinline__ long long loopStart() { return static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x; }
__device__ __forceinline__ long long loopStride() { return static_cast<long long>(gridDim.x) * blockDim.x; }

struct NcpSamples {
    const int* base;     // global index of the span's first control
    const double* w;     // 4 per sample
    const int* yarn;
    const double* s;
    const unsigned char* active;
    const unsigned char* segOk;
    int n;
};

__device__ __forceinline__ d3 evalSample(const d3* P, const NcpSamples& S, int i) {
    const int b = S.base[i];
    const double* w = S.w + 4 * i;
    return P[b] * w[0] + P[b + 1] * w[1] + P[b + 2] * w[2] + P[b + 3] * w[3];
}
// Coefficient of control g in the sample's combination (0 outside its four controls).
__device__ __forceinline__ double sampleCoef(const NcpSamples& S, int i, int g) {
    const int k = g - S.base[i];
    return k >= 0 && k < 4 ? S.w[4 * i + k] : 0.0;
}

struct NcpParams {
    double D, thr, hardC, kappaMax, hardK, bandHalf, cap, wL, wT, wA, kDeg;
    bool targets;
    int curvQ;     // curvature samples per span (t = q / (curvQ - 1))
    int contactModel;  // KnitProblemParams::ContactModel: 0 segment, 1 point, 2 point-segment, 3 hybrid
    double parallelCos;
};

// What the contact values are computed from when a live row stores its distance (contactHardValue).
struct ContactValueScale {
    double hardC, thr, D;
    int targets;
};
__host__ __device__ inline ContactValueScale contactValueScale(const NcpParams& pr) { return {pr.hardC, pr.thr, pr.D, pr.targets ? 1 : 0}; }

// Per-evaluation reductions (atomics: maxima as bits, counts; sums by the fixed-grid reductions).
struct EvalReductions {
    double maxViol[knit::kRowFamilies];
    unsigned long long nViol[knit::kRowFamilies];
    unsigned long long degenerate, notOwned, pinnedBand;
    double strainMax, pinnedBandViol;
    double merit, objective, strain2;  // filled by the sum reductions
    unsigned long long edges;
    double contactMerit, contactObjective;  // a trial's contact rows (reduced, not stored)
    // The memory report: contact rows with a term in the merit or a positive multiplier at this
    // evaluation (the rows whose state cannot be derived from the positions alone).
    unsigned long long liveRows;
};

// Sizes that change during a solve. Written by the refresh and activeLists kernels, read by the
// kernels that loop over the lists.
struct NcpCounts {
    int contactRows;                  // rows of the last broad phase (unchanged by a refresh that overflowed)
    int previousRows;                 // the rows before it (multiplier carry-over)
    unsigned long long candidates;    // candidate pairs of the last broad phase (statistics)
    unsigned long long registrations; // grid registrations of the last broad phase
    int activeRows;                   // rows with a term in the merit (last activeLists)
    int activeEntries;                // their (variable, row) entries
    int contactEntryOverflow;         // a contact row touches more than kContactRowEntries free controls
    int overflow;                     // NcpCapacityKind of the first capacity exceeded (0: none)
    long long overflowNeeded, overflowCapacity;
    // The memory report of the last broad phase: its rows whose two sample segments are closer than the
    // near cutoff (nearCutoff), and its rows that carried a positive multiplier from farther away.
    unsigned long long nearRows, carriedFar;
    // ContactStorage::LiveRows: the live list's length (the contact list the terms, multipliers and
    // active lists run over), its length before the evaluation that renumbered it, and the multipliers
    // a refresh carries (the sorted carry-over table's length).
    int liveRows, previousLive, carried;
};
enum NcpCapacityKind : int { kCapNone = 0, kCapRegistrations, kCapCandidates, kCapContactRows, kCapActiveEntries, kCapLiveRows };

// Fixed capacities of one call (host preparation, never changed by the iteration). The active slots
// are one per row of the contact list (general and contact): the active lists cannot overflow.
struct NcpCapacity {
    long long registrations = 0, contactRows = 0, activeEntries = 0;
    long long liveRows = 0;  // ContactStorage::LiveRows: the live list (the contact list); CompactRows: 0
};

// The row buffers that depend on the controls (two sets: current and trial; the device index
// NcpControl::current says which set is current).
struct NcpRowValues {
    // general rows (one array per set)
    double* gVal;    // hard value / objective residual
    double* gTval;   // target companion value (curvature), 0 when inactive
    double* gAux;    // edge length
    d3* gGrad;       // kGeneralRowEntries per row
    unsigned char* gFlag;  // 1 degenerate / not owned (row off)
    double* wres;    // d merit / d value
    double* wgn;     // Gauss-Newton weight
    double* merit;   // merit contribution
    double* obj;     // objective contribution
    double* strain2; // edges' strain^2
    // contact rows: one storage for both sets (a trial reduces its contact rows without storing them,
    // and accepting it evaluates them again); their merit and objective are summed where needed
    // (contactMeritTerms), not stored
    double* cVal;    // CompactRows (LiveRows: null, see cDist)
    double* cTval;   // CompactRows (LiveRows: null, see cDist)
    double* cWres;
    double* cWgn;
    unsigned char* cGradOn;  // 1: the row had a term in the merit when evaluated (its gradient is used)
    // LiveRows: the distance d instead of cVal and cTval (8 bytes less per live row); both follow from it by
    // the evaluation's own expressions (contactHardValue, contactTargetValue), so they are the same doubles.
    double* cDist;
};
// The terms of row r (general rows first, then contact rows).
__device__ __forceinline__ double rowWres(const NcpRowValues& V, int nG, int r) { return r < nG ? V.wres[r] : V.cWres[r - nG]; }
__device__ __forceinline__ double rowWgn(const NcpRowValues& V, int nG, int r) { return r < nG ? V.wgn[r] : V.cWgn[r - nG]; }
struct RowValuePair {
    NcpRowValues set[2];
    const int* current;  // NcpControl::current
};
enum class RowSet { Current, Trial };
__device__ __forceinline__ const NcpRowValues& rowValues(const RowValuePair& v, RowSet which) {
    const int c = *v.current;
    return v.set[which == RowSet::Current ? c : 1 - c];
}

struct NcpRows {
    // general rows (fixed per call)
    int nG = 0;
    const int* gFam;
    const int* gIdx;     // control (edge: first control; curvature: span's first control; band/cap/anchor: the control)
    const int* gSub;     // curvature q; anchor axis
    const int* gN;
    const int* gVar;     // kGeneralRowEntries per row
    const double* gRest;
    double* gLam;
    // Two lists (ContactStorage). The row list: every contact row of the last refresh, its pair and type
    // (nRows of them, NcpCounts::contactRows). The contact list: the rows whose multiplier, values and
    // terms are held (nC of them; the terms, multipliers and active lists run over it). CompactRows:
    // both are the same list (cRow, rowSlot and rowLive null). LiveRows: the contact list holds the
    // live rows, c -> row cRow[c]; a row is live when rowLive[row], its contact index then rowSlot[row].
    const int* nRows;
    const int* nC;
    const int2* cPair;           // per row
    const unsigned char* cType;  // per row: 0 segment-segment, 1 point-point, 2 point i - segment j
    double* cLam;                // per contact
    const int* cRow;             // LiveRows: contact -> row
    const int* rowSlot;          // LiveRows: row -> contact (where rowLive)
    const unsigned char* rowLive;
};
// The row of contact c (the row list's index of the contact list's entry).
__device__ __forceinline__ int contactRowOf(const NcpRows& R, int c) { return R.cRow ? R.cRow[c] : c; }
// The multiplier of row `row` (0 for a row that is not live: LiveRows holds no multiplier for it).
__device__ __forceinline__ double rowMultiplier(const NcpRows& R, int row) {
    if (!R.rowLive) return R.cLam[row];
    return R.rowLive[row] ? R.cLam[R.rowSlot[row]] : 0.0;
}

// Hinge terms (augmented Lagrangian, PHR): the merit and its gradient are exact (max(0, .)^2).
// Acceptance uses the exact merit.
__device__ __forceinline__ void hardTerm(double g, double lam, double rho, double& wres, double& wgn, double& m) {
    const double phi = g + lam / rho;
    if (phi > 0) {
        wres += rho * phi;
        m += 0.5 * rho * phi * phi;
        wgn += rho;
    }
}

// A contact row's merit and objective from its stored values (the sums evalContactRowsKernel's terms
// gave: the hard term, then the target term when the distance is below the target threshold).
__device__ __forceinline__ void contactMeritTerms(double val, double tval, double lam, double rho, double wT, double& m, double& o) {
    double wres = 0, wgn = 0;
    m = 0;
    o = 0;
    hardTerm(val, lam, rho, wres, wgn, m);
    if (tval > 0) {
        m += 0.5 * wT * tval * tval;
        o = 0.5 * wT * tval * tval;
    }
}

// Contact c's hard value (hardC - d) / D and target value max(0, thr - d) / D (0 without targets), as
// evalContactRowsKernel computes them: stored (CompactRows) or from the stored distance (LiveRows).
__device__ __forceinline__ double contactHardValue(const NcpRowValues& V, const ContactValueScale& s, int c) {
    return V.cDist ? (s.hardC - V.cDist[c]) / s.D : V.cVal[c];
}
__device__ __forceinline__ double contactTargetValue(const NcpRowValues& V, const ContactValueScale& s, int c) {
    if (!V.cDist) return V.cTval[c];
    const double d = V.cDist[c];
    return s.targets && d < s.thr ? (s.thr - d) / s.D : 0.0;
}

// The key of contact row (i, j, type) (ncp_reference.cpp rowKey): the multipliers' identity across refreshes.
__device__ __forceinline__ unsigned long long contactKey(int i, int j, int type) {
    return (2ull << 60) ^ (static_cast<unsigned long long>(i) << 30) ^ static_cast<unsigned long long>(j) ^ (static_cast<unsigned long long>(type) << 57);
}

// The free controls of contact row (i, j) in entry order: side i, then side j, each from its first
// sample's span start to its second sample's span end, free ones only, each once. At most
// kContactRowEntries are written; `overflow` is set when the row touches more.
__device__ __forceinline__ int contactControls(const NcpSamples& S, const int* varOf, int i, int j, int* C, bool& overflow) {
    int n = 0;
    for (int side = 0; side < 2; ++side) {
        const int s0 = side == 0 ? i : j;
        const int g0 = S.base[s0], g1 = S.base[s0 + 1] + 3;
        for (int g = g0; g <= g1; ++g) {
            if (varOf[g] < 0) continue;
            bool dup = false;
            for (int e = 0; e < n; ++e) dup = dup || C[e] == g;
            if (dup) continue;
            if (n >= kContactRowEntries) {
                overflow = true;
                continue;
            }
            C[n++] = g;
        }
    }
    return n;
}

// d value / d x_g = n * this, for control g of contact row (i, j) with closest points at s and t:
// the value (hardC - d) / D falls with the distance along the normal n = (pa - pb) / d.
__device__ __forceinline__ double contactCoefficient(const NcpSamples& S, int i, int j, double s, double t, int g, double D) {
    const double ca = (1 - s) * sampleCoef(S, i, g) + s * sampleCoef(S, i + 1, g);
    const double cb = (1 - t) * sampleCoef(S, j, g) + t * sampleCoef(S, j + 1, g);
    return -(ca - cb) / D;
}

// The distance of contact row (i, j, type) between the sample points X: its closest points and their
// parameters (type 1: points i and j; type 2: point i and segment j; type 0: the two segments).
__device__ __forceinline__ double contactDistance(const d3* X, int i, int j, int ty, d3& pa, d3& pb, double& s, double& t) {
    const d3 a0 = X[i], a1 = X[i + 1], b0 = X[j], b1 = X[j + 1];
    s = 0;
    t = 0;
    const double d = ty == 1 ? norm(a0 - b0) : ty == 2 ? segmentDistance(a0, a0, b0, b1, s, t) : segmentDistance(a0, a1, b0, b1, s, t);
    pa = a0 + (a1 - a0) * s;
    pb = b0 + (b1 - b0) * t;
    return d;
}

// Per-evaluation tolerance: counts nViol (the score's allowance); maxima are exact.
constexpr double kViolationTolerance = 1e-3;

// The band's surface at every free control, evaluated once per evaluation (both band rows of a
// control read it): foot, unit normal, and whether it exists (bandSurface).
struct BandFeet {
    d3* foot;
    d3* nrm;
    unsigned char* ok;
};

// What kernels read, by value. What a kernel writes stays a parameter of its own.
struct NcpView {
    NcpParams params;
    NcpSamples samples;
    NcpRows rows;
    BandParams band;
    const DoubleProjectionState* bindState;
    const d3* P0;
    const int* varOf;
    const int* ctrlOfVar;
    const unsigned char* pin;
    int variables;
};
static_assert(sizeof(NcpView) < 4096, "CUDA limits a kernel's arguments to 4 KB");

// The rows with a term in the merit (the products and the assembly read only these), one slot each
// in row order; per variable its entries (slot << 4) | entry index, sorted, so every sum over a
// variable's rows adds in row order. A contact slot holds the row's gradient: the unit normal (0: the
// row had no term when evaluated, no gradient) and the closest points' parameters s, t.
struct ActiveLists {
    const int* start;  // variables + 1
    const int* entry;
    const int* rows;   // slot -> row, NcpCounts::activeRows of them
    const d3* normal;  // per slot (contact rows)
    const double* st;  // 2 per slot (contact rows)
};

// The gradient of variable a's entry (slot k, entry e). A contact entry's control is the variable's.
__device__ __forceinline__ d3 entryGradient(const NcpView& v, const NcpRowValues& V, const ActiveLists& L, int k, int e, int a) {
    const int r = L.rows[k];
    if (r < v.rows.nG) return V.gGrad[kGeneralRowEntries * r + e];
    const int2 pair = v.rows.cPair[contactRowOf(v.rows, r - v.rows.nG)];
    return L.normal[k] * contactCoefficient(v.samples, pair.x, pair.y, L.st[2 * k], L.st[2 * k + 1], v.ctrlOfVar[a], v.params.D);
}

// The decisions and the record of one solve (ncp_control.cu writes it; the host reads it at the end).
enum NcpStatus : int {
    kRunning = 0, kConverged, kStoppedProgress, kStalledFeasible, kInfeasible, kStepBudget, kIterationCap, kTimeLimit, kCapacity,
    kContactEntries, kCancelled
};
struct NcpControl {
    // parameters (host, once)
    double tolViol, tolViolContact, tolStat, tolObjRel, objFloor, rhoMax, progressTol, infeasibleRho;
    double delta, sweptEps, epsRelative, refreshAt, innerTol0, pcgTolerance;
    long long maxLmSteps, timeLimitNs, startNs;
    int maxOuter, stallWindow, infeasibleWindow, innerMax, maxAttempts, maxHalvings, pcgMaxIterations;
    // the iteration
    double rho, mu, nu, eps, firstGradient, previousViolation;
    double alpha, decrease, regained, predicted, gain, lastTrialRho;
    int current;  // which row value set is current
    int outer, inner, attempt, halvings, flatFeasible, stepTaken, accepted, pinnedNow;
    double objectives[8], violations[8];  // the last stallWindow + 1 / infeasibleWindow + 1 outer iterations
    int objectiveCount, violationCount;
    NcpEvalSummary currentEval, trialEval;
    double gradMax, diagMax, motion;
    // the step
    double gdx, dDd, alphaCov, dxMax, rDotX;
    unsigned long long swept[3];  // touching, unresolved, closer than eps at the start
    // conjugate gradients
    double pcgRz, pcgPq, pcgAlpha, pcgBeta, pcgRr, pcgBb;
    int pcgIteration, pcgBreakdown;
    // KKT
    double stationarity, complementarity;
    unsigned long long activeMultipliers;
    // the result
    int status;
    double statusValue[4];  // the numbers of the status message
    // statistics
    long long lmSteps, rejected, rejectedMerit, sweptRejected, sweptUnresolved, ownershipRejected, ownershipPinned, regainedTrials;
    long long refreshes, localSteps, backtracks, coverageLimited, passThroughRejected, pcgIterations, pcgCapped, nonSpdBlocks, solves;
    long long maxContactRows, pinnedCount, pcgBreakdowns;
    long long maxCandidates, maxActiveRows, maxActiveEntries;  // the memory report: the largest lists of the solve
    long long maxNearRows, maxLiveRows, maxCarriedFar;         // NcpCounts::nearRows / carriedFar, EvalReductions::liveRows
    unsigned long long coverageChecks, coverageViolations;  // NcpOptions::auditCoverage
    unsigned long long scratchCount;  // a count a kernel accumulates for the next controller
    // 1 while a launch continues a solve whose refresh outgrew the contact capacity (solveNcp grew it):
    // the start and outer-iteration resets are skipped, and the next inner iteration redoes that refresh's
    // contact rows (its binding continuation is already done). Cleared after that refresh.
    int resume;
};

// The early-stop rules (ncpEarlyStop), on the host and in the controller.
// Progress: the magnitude of the change, so an outer iteration that raised the objective by more than
// the tolerance (multipliers or rho trading objective for feasibility) continues.
__host__ __device__ inline NcpEarlyStop earlyStopRule(double progressTol, int infeasibleWindow, double infeasibleRho, bool feasible, double objPrev,
                                                      double obj, double rho, double viol, double violBack) {
    if (feasible && progressTol > 0 && objPrev > 0 && fabs(objPrev - obj) <= progressTol * objPrev) return NcpEarlyStop::Progress;
    if (!feasible && infeasibleWindow > 0 && rho >= infeasibleRho && violBack >= 0 && viol > 0.5 * violBack) return NcpEarlyStop::Infeasible;
    return NcpEarlyStop::None;
}

// --- buffers, by purpose
struct ControlBuffers {
    DeviceBuffer P, Q;              // current and trial controls
    DeviceBuffer P0, Xref;          // at the call's start, at the last broad phase
    DeviceBuffer varOf, ctrlOfVar;  // control -> variable (-1 fixed), variable -> control
    DeviceBuffer pin;               // per variable: 1 = held fixed for the rest of the call (pinLostOwnership)
    d3* current() const { return static_cast<d3*>(P.data()); }
    d3* trial() const { return static_cast<d3*>(Q.data()); }
    const d3* start() const { return static_cast<d3*>(P0.data()); }
    const d3* reference() const { return static_cast<d3*>(Xref.data()); }
};
struct SampleBuffers {
    DeviceBuffer base, w, yarn, s, active, segOk;
};
struct GeneralRowBuffers {
    DeviceBuffer fam, idx, sub, n, var, rest, lambda;
};
struct ContactRowBuffers {
    DeviceBuffer pair, type;   // per row
    DeviceBuffer lambda;       // per contact (CompactRows: per row)
    DeviceBuffer sampleRows, sampleRowOffset;  // per sample: its candidates' rows and their offsets (nS + 1)
};
// ContactStorage::LiveRows: the live list next to the row list. slot: per row (+1), the live flags'
// exclusive scan (the row's contact index where `live`); live: per row, 1 when the row has a contact
// index; row: per contact, its row.
struct LiveRowBuffers {
    DeviceBuffer slot, live, row;
};
// Which penalty an evaluation of the current contact rows uses: rho, or the rho its trial was evaluated
// with (an accepted trial, evaluated again at the new current state).
enum class EvaluationPenalty { CurrentRho, LastTrialRho };
struct RowValueSet {
    DeviceBuffer gVal, gTval, gAux, gGrad, gFlag, wres, wgn, merit, obj, strain2;  // general rows
    DeviceBuffer cVal, cTval, cWres, cWgn, cGradOn, cDist;                         // contact rows (set 0 only)
    NcpRowValues view{};
};
struct ActiveRowBuffers {  // the rows with a nonzero term in the merit (ActiveLists)
    DeviceBuffer slotOf;   // per row (+1): the exclusive scan of the active flags, the row's slot
    DeviceBuffer rows, count, start, cursor, entry, normal, st;
};
// The grid of the broad phase: its shape is fixed for the call (items outside it are clamped to
// its border cells, which keeps every query conservative).
struct FixedGrid {
    float cell = 1;
    long long lo[3] = {0, 0, 0}, n[3] = {1, 1, 1};
    int endBit = 64;  // sort keys (cell << 32 | item) on bits [0, endBit)
    DeviceBuffer boxLo, boxHi, counts, offsets, keys, sorted, cellCounts, cellStart, entries;
    GridView view;
};
// The multiplier carry-over's arrays (capacity n each): the previous rows' keys and multipliers and
// the second halves of their double buffer (NcpDevice::State::carryBuffers).
struct CarryBuffers {
    unsigned long long* keys;
    unsigned long long* keysAlt;
    double* lambda;
    double* lambdaAlt;
    long long n;
};
// LiveRows: a refresh's carried multipliers sorted by row key (the first NcpCounts::carried entries;
// the largest key after them), in the active entries' memory until the rows are rebuilt.
struct CarryTable {
    const unsigned long long* keys;
    const double* lambda;
};
struct LinearSystem {
    DeviceBuffer grad, gnorm, hdiag, hblock;  // per variable: gradient, its norm, H's diagonal, H's 3x3 diagonal block (6)
    DeviceBuffer rowT, rowS;                  // per row: the product's row scalar, an edge row's 3-vector
    DeviceBuffer x, r, z, p, q, precond;      // conjugate gradients; precond: the damped blocks' inverses (6 per variable)
    DeviceBuffer partials;                    // fixed-grid block partial sums
};
struct ReductionBuffers {
    DeviceBuffer evaluation, acceptEvaluation, blocks, regained, perVariable, perVariable2;
};

struct NcpDevice::State {
    NcpStats st;
    NcpParams pr{};
    double delta = 0.1, window = 0;
    int nCtrl = 0, N = 0, nG = 0, nSamples = 0;
    bool hasBind = false;
    long long capLimit = 0x7fffffffLL;  // index capacity, lowered by tests
    ContactStorage storage = ContactStorage::CompactRows;
    NcpCapacity capacity;
    cudaStream_t stream = nullptr;  // every launch of the call; the graph is captured from it

    ControlBuffers controls;
    SampleBuffers samples;
    NcpSamples S{};
    GeneralRowBuffers general;
    ContactRowBuffers contacts;
    LiveRowBuffers liveList;  // LiveRows only
    RowValueSet values[2];
    ActiveRowBuffers active;
    FixedGrid grid;
    LinearSystem system;
    ReductionBuffers reductions;
    DevicePrimitives prims;
    DeviceBuffer countsBuffer, controlBuffer, bandFootBuffer, bandNrmBuffer, bandOkBuffer;
    // The sample points at P and at Q (evalSample once per state): the contact kernels read these
    // instead of four controls and four weights per sample.
    DeviceBuffer samplesAtP, samplesAtQ;

    // The band.
    std::shared_ptr<DeviceMesh> mesh;
    BandParams bp{};
    DeviceBuffer bindState, bindResult;
    DeviceBuffer mPos, mVn, mFn;  // the mesh in double, next to the resident float mesh's grid and adjacency

    // Host copies.
    std::vector<int> ctrlBase;  // per yarn + 1
    knit::Controls Pshape;      // sizes per yarn (download)

    ~State();
    NcpCounts* counts() const { return static_cast<NcpCounts*>(countsBuffer.data()); }
    const d3* pointsAtP() const { return static_cast<const d3*>(samplesAtP.data()); }
    const d3* pointsAtQ() const { return static_cast<const d3*>(samplesAtQ.data()); }
    NcpControl* control() const { return static_cast<NcpControl*>(controlBuffer.data()); }
    BandFeet feet() const { return {static_cast<d3*>(bandFootBuffer.data()), static_cast<d3*>(bandNrmBuffer.data()), bandOkBuffer.chars()}; }
    NcpRows rows() const {
        NcpRows R;
        R.nG = nG;
        R.gFam = general.fam.ints();
        R.gIdx = general.idx.ints();
        R.gSub = general.sub.ints();
        R.gN = general.n.ints();
        R.gVar = general.var.ints();
        R.gRest = general.rest.doubles();
        R.gLam = general.lambda.doubles();
        R.nRows = &counts()->contactRows;
        R.cPair = static_cast<int2*>(contacts.pair.data());
        R.cType = contacts.type.chars();
        R.cLam = contacts.lambda.doubles();
        if (storage == ContactStorage::LiveRows) {
            R.nC = &counts()->liveRows;
            R.cRow = liveList.row.ints();
            R.rowSlot = liveList.slot.ints();
            R.rowLive = liveList.live.chars();
        } else {
            R.nC = &counts()->contactRows;
            R.cRow = nullptr;
            R.rowSlot = nullptr;
            R.rowLive = nullptr;
        }
        return R;
    }
    // The contact list's capacity: the live list (LiveRows) or the rows (CompactRows).
    long long contactListCapacity() const { return storage == ContactStorage::LiveRows ? capacity.liveRows : capacity.contactRows; }
    NcpView view() const {
        return {pr, S, rows(), bp, static_cast<DoubleProjectionState*>(bindState.data()), controls.start(), controls.varOf.ints(), controls.ctrlOfVar.ints(), controls.pin.chars(), N};
    }
    RowValuePair valuePair() const { return {{values[0].view, values[1].view}, &control()->current}; }
    ActiveLists activeLists() const {
        return {active.start.ints(), active.entry.ints(), active.rows.ints(), static_cast<const d3*>(active.normal.data()), active.st.doubles()};
    }

    // Totals that int-indexed arrays, kernels and library calls consume must stay below capLimit.
    // Checked at setup, before anything is allocated.
    void checkIndexLimit(long long total, const char* what) const {
        if (total >= 0 && total < capLimit) return;
        throw DeviceError(std::string("ncp-gpu: capacity exceeded: ") + what + " " + std::to_string(total) + " >= " + std::to_string(capLimit));
    }

    // ncp_setup.cu
    void uploadMesh(const knit::TriMesh& tm);
    void allocateCapacity(long long contactRows);
    // The capacity buffers' bytes at a capacity (plus the scratch growth they need), and now.
    long long capacityBytes(long long contactRows) const;
    long long capacityBytesNow() const;
    // The capacity for `rows` contact rows: the headroom where the memory cap allows it, else the most
    // that fits, but never below kCapacityFloor x rows (that allocation is then refused by the cap).
    long long chooseCapacity(long long rows) const;
    // After a refresh overflowed (kCapacity, contact rows): the rows kept on the host, every capacity
    // buffer reallocated for `needed` rows (chooseCapacity), the rows put back, the overflow cleared.
    void regrowContactRows(long long needed);
    // ncp_live.cu (ContactStorage::LiveRows): every buffer that follows the row and live capacities (no
    // regrow: an overflow ends the solve), its bytes, and the device sequences.
    void allocateLiveCapacity(long long rows, long long live);
    // Both row value sets: the general rows' values per set, the contact list's in set 0 (shared).
    void allocateValueSets();
    // The refresh's carry-over: the live rows' keys and positive multipliers, sorted, before the rows are
    // rebuilt (enqueueLiveCarryTable); afterwards each new row with a carried multiplier becomes live with
    // it (enqueueLiveCarryOver). Returns the sorted table's place, decided when enqueued.
    CarryTable enqueueLiveCarryTable();
    void enqueueLiveCarryOver(const CarryTable& table);
    // The current rows at P: every row's value decides whether it is live (a term in the merit or a
    // positive multiplier), the live list is renumbered in row order (multipliers move along), and the live
    // rows' values and terms are stored. Violation maxima and counts into `sum`.
    void enqueueLiveEvaluation(EvaluationPenalty penalty, EvalReductions* sum);
    // ncp_refresh.cu: the device sequences (enqueued on `stream`, no host reads)
    void enqueueSamplePoints(RowSet which);  // the sample points at P (Current) or Q (Trial)
    void enqueueMotion();                 // control()->motion := max |P - Xref| over the variables
    void enqueueRefresh();                // binding, Xref := P, broad phase, contact rows (multipliers carried)
    void enqueueGrid();                   // the broad phase's grid of the segments at P
    void enqueueContactRows(float pad);
    CarryBuffers carryBuffers() const;
    // The memory report: the starting state's candidates, contact rows, and rows with a term in the
    // merit at multipliers 0, into totals[0..2] (after enqueueGrid).
    void enqueueStartCounts(float pad, unsigned long long* totals);
    std::vector<NcpMemoryFamily> memoryFamilies() const;
    void enqueueActiveLists(RowSet which);
    // ncp_evaluate.cu
    void enqueueLinearize(RowSet which);  // rows at P (Current) or Q (Trial); summary into control()->currentEval / trialEval
    void enqueueAccept();
    void enqueueMultipliers();
    void enqueueTerms();                  // the current rows' terms at control()->rho; control()->currentEval.merit
    // ncp_step.cu
    void enqueueAssemble();               // gradient, diagonal blocks, gradMax, diagMax at the current rows
    void enqueueStepSummary();
    void enqueueTrial();                  // Q := proj(P + alpha x), alpha = control()->alpha
    void enqueueSwept();
    void enqueueCoverageAudit();         // tests: NcpOptions::auditCoverage
    void enqueuePinLost();
    void enqueueRegained();
    void enqueueKkt();
    // ncp_pcg.cu
    void enqueuePcgSetup(cudaGraphConditionalHandle loop);      // sets `loop` to 1 when there is anything to solve
    void enqueuePcgIteration(cudaGraphConditionalHandle loop);  // one iteration; sets `loop` to 0 when it stops
    void enqueuePcgFinish();
    void enqueueOperatorProduct();  // diagnostics: q := (H + mu diag(H) + eps I) p at the current rows
    void enqueueProductPasses();    // that product and the partials of p.q (the PCG iteration's first part)
};

} // namespace knit::gpu
