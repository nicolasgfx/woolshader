// ncp-gpu: the nonlinear iteration as one CUDA graph whose decisions are made on the device (ncp.h).
// An augmented Lagrangian outer loop over Levenberg-Marquardt steps:
//
//   refresh; startCheck (a capacity exceeded: no outer loop); evaluate the current rows
//   (a resumed launch, after solveNcp grew the contact capacity: only resume, then the loops below)
//   WHILE outer:                                     outer iteration (multipliers, penalty, stops)
//     WHILE inner:                                   Levenberg-Marquardt steps at fixed multipliers
//       motion; innerHead                            step budget, cancel; refresh needed?
//       IF refresh:  refresh; refreshCheck; IF ok: evaluate the current rows
//       IF step:     assemble; afterAssemble         inner stop on the gradient
//         WHILE attempt:                             damped attempts of one step
//           PCG setup; WHILE pcg: PCG iteration;     conjugate gradients to |r| <= tol |b|
//           PCG finish; step reductions; afterStep
//           WHILE halve: trial; swept test; halve    shorten while the motion touches a pair
//           afterSwept
//           IF evaluate: trial rows; regained merit; trialDecide
//             IF pin:    pin controls that lose their binding; IF pinned: evaluate the current rows
//             IF accept: accept the trial
//           attemptEnd                               damping update (Nielsen), next attempt?
//         innerEnd                                   next step?
//     IF outer tail: KKT; multipliers; outerDecide; terms at the new penalty
//     outerFinish                                    next outer iteration?
//
// A refresh whose contact rows outgrow their capacity ends the launch (kCapacity) with the previous rows
// intact; solveNcp grows the capacity and launches again with NcpControl::resume set: the outer and
// inner loops continue where they stopped (no resets) and the first inner iteration redoes that
// refresh's contact rows, so the solve takes the same path as with the larger capacity from the start.
//
// Each "Decide" is a one-thread controller kernel that reads NcpControl and the counts, writes the
// decision there and sets the conditions of the graph's WHILE / IF nodes (cudaGraphSetConditional).
// The host builds the graph from the same enqueue sequences the tests run one by one, launches it
// once and reads NcpControl when it has finished: no value read on the host decides anything in
// between.
#include <algorithm>
#include <cmath>
#include <cstdio>

#include "gpu/device/runtime.h"
#include "gpu/solver/ncp_device.cuh"
#include "model/cancel.h"

namespace knit::gpu {

namespace {

constexpr double kPcgTolerance = 0.05;   // PCG relative residual (inexact Levenberg-Marquardt steps; the merit test is exact)
constexpr int kPcgMaxIterations = 2000;  // PCG safety net (reported; the iterate is used)
constexpr double kSweptEps = 1e-3;       // swept test: a trial is shortened / rejected when a pair comes within kSweptEps * D
constexpr double kRefreshAt = 0.5;       // broad phase when the cumulative motion exceeds kRefreshAt * delta
constexpr int kInnerMax = 50;            // Levenberg-Marquardt steps per multiplier estimate (safety net, reported)
constexpr double kInnerTol0 = 0.1;       // inner stop: |grad| <= max(tolStat, g0 * kInnerTol0^(outer + 1))
constexpr int kMaxAttempts = 12;         // damped attempts of one Levenberg-Marquardt step
constexpr int kMaxHalvings = 6;          // a trial is halved at most this often while pairs touch
constexpr int kHistory = 8;              // NcpControl::objectives / violations

// The conditions of the graph's WHILE / IF nodes (each one created on the graph that holds its node).
struct Conditions {
    cudaGraphConditionalHandle outer, inner, refresh, afterRefresh, step, attempt, pcg, halve, evaluate, pin, relinearize, accept, outerTail;
};

__device__ __forceinline__ long long deviceNanoseconds() {
    unsigned long long t;
    asm volatile("mov.u64 %0, %%globaltimer;" : "=l"(t));
    return static_cast<long long>(t);
}

__device__ __forceinline__ bool isTerminal(int status) { return status != kRunning; }

// The maximum violation the feasibility test compares with tolViol: contact against its own allowance
// (rescaled to tolViol), the other hard families as they are (NcpSolveParams::tolViolContact).
__device__ double feasibilityViolation(const NcpControl* ctl, const NcpEvalSummary& s) {
    double v = 0;
    for (int f = 0; f < knit::kRowFamilies; ++f) {
        if (f == FEdge || f == FContactT || f == FCurvT || f == FAnchor) continue;
        const double x = f == FContact && ctl->tolViolContact > 0 ? s.maxViol[f] * ctl->tolViol / ctl->tolViolContact : s.maxViol[f];
        v = fmax(v, x);
    }
    return v;
}

__global__ void startKernel(NcpControl* ctl, Conditions c) {
    ctl->startNs = deviceNanoseconds();
    ctl->outer = 0;
    cudaGraphSetConditional(c.outer, ctl->maxOuter > 0 ? 1 : 0);
}

// A resumed launch: the outer loop continues (its counters and the clock are as the last launch left them).
__global__ void resumeKernel(Conditions c) { cudaGraphSetConditional(c.outer, 1); }

__global__ void outerBeginKernel(NcpControl* ctl, Conditions c) {
    if (!ctl->resume) {  // a resumed launch continues the outer iteration it stopped in
        ctl->nu = 2;
        ctl->inner = 0;
    }
    cudaGraphSetConditional(c.inner, 1);
}

// The head of a Levenberg-Marquardt step: cancel, the step budget, and whether the motion since the
// last broad phase needs a refresh.
__global__ void innerHeadKernel(NcpControl* ctl, Conditions c, const volatile int* abortFlag) {
    const bool cancelled = *abortFlag != 0;
    if (cancelled) ctl->status = kCancelled;
    const bool budgetSpent = ctl->maxLmSteps > 0 && ctl->lmSteps >= ctl->maxLmSteps;
    if (cancelled || budgetSpent) {
        cudaGraphSetConditional(c.inner, 0);
        cudaGraphSetConditional(c.refresh, 0);
        cudaGraphSetConditional(c.step, 0);
        return;
    }
    cudaGraphSetConditional(c.refresh, ctl->resume || ctl->motion > ctl->refreshAt * ctl->delta ? 1 : 0);
    cudaGraphSetConditional(c.step, 1);
}

// After a refresh: a capacity it exceeded ends the iteration at the last accepted state.
__global__ void refreshCheckKernel(NcpControl* ctl, const NcpCounts* counts, Conditions c) {
    ctl->resume = 0;
    const bool failed = counts->overflow != 0 || counts->contactEntryOverflow != 0;
    if (failed) {
        ctl->status = counts->contactEntryOverflow ? kContactEntries : kCapacity;
        ctl->statusValue[0] = counts->overflow;
        ctl->statusValue[1] = static_cast<double>(counts->overflowNeeded);
        ctl->statusValue[2] = static_cast<double>(counts->overflowCapacity);
        cudaGraphSetConditional(c.step, 0);
        cudaGraphSetConditional(c.inner, 0);
    }
    cudaGraphSetConditional(c.afterRefresh, failed ? 0 : 1);
}

// After the first refresh: a capacity it exceeded ends the solve before any step (nothing moved).
__global__ void startCheckKernel(NcpControl* ctl, const NcpCounts* counts, Conditions c) {
    if (counts->overflow == 0 && counts->contactEntryOverflow == 0) return;
    ctl->status = counts->contactEntryOverflow ? kContactEntries : kCapacity;
    ctl->statusValue[0] = counts->overflow;
    ctl->statusValue[1] = static_cast<double>(counts->overflowNeeded);
    ctl->statusValue[2] = static_cast<double>(counts->overflowCapacity);
    cudaGraphSetConditional(c.outer, 0);
}

// A live list that overflowed in an evaluation of the current state (ContactStorage::LiveRows) ends the
// iteration at the last accepted state, as an overflowing refresh does: the inner and attempt loops stop
// and the outer tail and finish see the terminal status. (CompactRows: a refresh's overflow already
// stopped the step, so this never fires there.)
__device__ bool stopOnLiveOverflow(NcpControl* ctl, const NcpCounts* counts, Conditions c) {
    if (counts->overflow != kCapLiveRows) return false;
    ctl->status = kCapacity;
    ctl->statusValue[0] = counts->overflow;
    ctl->statusValue[1] = static_cast<double>(counts->overflowNeeded);
    ctl->statusValue[2] = static_cast<double>(counts->overflowCapacity);
    cudaGraphSetConditional(c.attempt, 0);
    cudaGraphSetConditional(c.inner, 0);
    return true;
}

// After the assembly: the damping's regularization from the diagonal, and the inner stop on the gradient.
__global__ void afterAssembleKernel(NcpControl* ctl, const NcpCounts* counts, Conditions c) {
    if (stopOnLiveOverflow(ctl, counts, c)) return;
    ctl->eps = fmax(1e-12, ctl->epsRelative * ctl->diagMax);
    if (ctl->firstGradient < 0) ctl->firstGradient = ctl->gradMax;
    ctl->stationarity = ctl->gradMax;
    const double stop = fmax(ctl->tolStat, ctl->firstGradient * pow(ctl->innerTol0, ctl->outer + 1));
    if (ctl->gradMax <= stop) {
        ctl->stepTaken = 0;
        cudaGraphSetConditional(c.attempt, 0);
        cudaGraphSetConditional(c.inner, 0);
        return;
    }
    ctl->stepTaken = 1;
    ctl->attempt = 0;
    ctl->accepted = 0;
    cudaGraphSetConditional(c.attempt, 1);
}

// After the linear solve: the step counts, its first length is the coverage fraction.
__global__ void afterStepKernel(NcpControl* ctl, Conditions c) {
    ++ctl->lmSteps;
    ctl->alpha = fmin(1.0, ctl->alphaCov);
    if (ctl->alpha < 1) ++ctl->coverageLimited;
    ctl->halvings = 0;
    ctl->pinnedNow = 0;
    cudaGraphSetConditional(c.halve, 1);
}

// Shorten (halve) while the straight motion brings a pair within the swept test's eps.
__global__ void halveKernel(NcpControl* ctl, Conditions c) {
    const unsigned long long touching = ctl->swept[0] + ctl->swept[1];
    if (touching == 0 || ctl->halvings >= ctl->maxHalvings) {
        cudaGraphSetConditional(c.halve, 0);
        return;
    }
    ++ctl->halvings;
    ctl->alpha *= 0.5;
    cudaGraphSetConditional(c.halve, 1);
}

__global__ void afterSweptKernel(NcpControl* ctl, Conditions c) {
    ctl->sweptUnresolved += static_cast<long long>(ctl->swept[1]);
    const bool touching = ctl->swept[0] + ctl->swept[1] > 0;
    if (touching) {
        ++ctl->sweptRejected;
        ++ctl->passThroughRejected;
        ctl->accepted = 0;
    }
    cudaGraphSetConditional(c.evaluate, touching ? 0 : 1);
}

// The trial against the current state: ownership loss pins, a merit increase rejects, else accept.
// The acceptance compares both states on the rows both have (regained rows enter once accepted).
__global__ void trialDecideKernel(NcpControl* ctl, Conditions c) {
    const NcpEvalSummary& cur = ctl->currentEval;
    const NcpEvalSummary& tri = ctl->trialEval;
    if (ctl->regained > 0) ++ctl->regainedTrials;
    // Undamped model: dx.H.dx = dx.(H + mu D + eps I).dx - dDd = -g.dx - r.dx - dDd.
    const double a = ctl->alpha;
    ctl->predicted = -a * ctl->gdx + 0.5 * a * a * (ctl->gdx + ctl->dDd + ctl->rDotX);
    ctl->gain = ctl->predicted > 0 ? (cur.merit - (tri.merit - ctl->regained)) / ctl->predicted : -1;
    int pin = 0, accept = 0;
    if (tri.notOwned > cur.notOwned) {
        // Hold the controls that would lose their binding and retry at the same damping:
        // shortening the step (more damping) only creeps toward the ownership boundary.
        ++ctl->ownershipRejected;
        pin = 1;
    } else if (!(tri.merit - ctl->regained < cur.merit)) {
        ++ctl->rejectedMerit;
    } else {
        accept = 1;
        ctl->decrease = cur.merit - (tri.merit - ctl->regained);
    }
    ctl->accepted = accept;
    cudaGraphSetConditional(c.pin, pin);
    cudaGraphSetConditional(c.accept, accept);
}

// The pinned controls' band rows are now off: the current merit is evaluated again without them, so
// the retried trial is compared on the same rows.
__global__ void pinDecideKernel(NcpControl* ctl, Conditions c) {
    const long long pinned = static_cast<long long>(ctl->scratchCount);
    ctl->ownershipPinned += pinned;
    ctl->pinnedCount += pinned;
    ctl->pinnedNow = pinned > 0;
    cudaGraphSetConditional(c.relinearize, pinned > 0 ? 1 : 0);
}

// The damping after an attempt (Nielsen), and whether another attempt follows.
__global__ void attemptEndKernel(NcpControl* ctl, const NcpCounts* counts, Conditions c) {
    if (stopOnLiveOverflow(ctl, counts, c)) return;
    if (ctl->accepted) {
        const double g = fmin(ctl->gain, 1.0);
        const double t = 2 * g - 1;
        ctl->mu = fmax(ctl->mu * fmax(1.0 / 3, 1 - t * t * t), 1e-12);
        ctl->nu = 2;
        cudaGraphSetConditional(c.attempt, 0);
        return;
    }
    ++ctl->rejected;
    ++ctl->backtracks;
    if (!ctl->pinnedNow) {  // a retry without the pinned controls needs no more damping
        ctl->mu *= ctl->nu;
        ctl->nu *= 2;
    }
    ++ctl->attempt;
    cudaGraphSetConditional(c.attempt, ctl->attempt < ctl->maxAttempts ? 1 : 0);
}

// After a step: stop the inner loop on a rejected step, a negligible decrease or kInnerMax steps.
__global__ void innerEndKernel(NcpControl* ctl, const NcpCounts* counts, Conditions c) {
    if (stopOnLiveOverflow(ctl, counts, c)) return;
    if (!ctl->stepTaken) return;  // the gradient stop already ended the loop
    ++ctl->localSteps;
    const bool more = ctl->accepted && ctl->decrease > 1e-12 * fmax(ctl->currentEval.merit, 1e-300) && ctl->inner + 1 < ctl->innerMax;
    ++ctl->inner;
    cudaGraphSetConditional(c.inner, more ? 1 : 0);
}

__global__ void outerTailCheckKernel(const NcpControl* ctl, Conditions c) {
    const int s = ctl->status;
    cudaGraphSetConditional(c.outerTail, s == kCapacity || s == kContactEntries || s == kCancelled ? 0 : 1);
}

__device__ void pushHistory(double* values, int& count, int keep, double v) {
    if (count == keep) {
        for (int i = 1; i < count; ++i) values[i - 1] = values[i];
        --count;
    }
    values[count++] = v;
}

// The outer iteration's end (finishOuterIteration): convergence, the stalled
// and early stops, then the penalty for the next outer iteration. The KKT measures were taken at the
// accepted state before the multiplier update.
__global__ void outerDecideKernel(NcpControl* ctl) {
    const double stationarity = ctl->gradMax, complementarity = ctl->complementarity;
    ctl->stationarity = stationarity;
    const double violation = feasibilityViolation(ctl, ctl->currentEval);
    const bool feasible = violation <= ctl->tolViol && ctl->currentEval.degenerate == 0;
    pushHistory(ctl->objectives, ctl->objectiveCount, ctl->stallWindow + 1, ctl->currentEval.objective);
    const double front = ctl->objectives[0], back = ctl->objectives[ctl->objectiveCount - 1];
    const bool objectiveFlat = ctl->objectiveCount > ctl->stallWindow && fabs(back - front) <= ctl->tolObjRel * fmax(ctl->objFloor, fabs(front));
    if (feasible && stationarity <= ctl->tolStat && complementarity <= ctl->tolViol) {
        ctl->status = kConverged;
        ++ctl->outer;
        return;
    }
    ctl->flatFeasible = feasible && objectiveFlat ? ctl->flatFeasible + 1 : 0;
    if (ctl->flatFeasible > ctl->stallWindow) {
        ctl->status = kStalledFeasible;
        ctl->statusValue[0] = stationarity;
        ctl->statusValue[1] = complementarity;
        ++ctl->outer;
        return;
    }
    pushHistory(ctl->violations, ctl->violationCount, ctl->infeasibleWindow + 1, violation);
    const double previousObjective = ctl->objectiveCount >= 2 ? ctl->objectives[ctl->objectiveCount - 2] : -1;
    const double earlierViolation = ctl->infeasibleWindow > 0 && ctl->violationCount > ctl->infeasibleWindow ? ctl->violations[0] : -1;
    const NcpEarlyStop early = earlyStopRule(ctl->progressTol, ctl->infeasibleWindow, ctl->infeasibleRho, feasible, previousObjective, ctl->currentEval.objective,
                                             ctl->rho, violation, earlierViolation);
    if (early == NcpEarlyStop::Progress) {
        ctl->status = kStoppedProgress;
        ctl->statusValue[0] = fabs(previousObjective - ctl->currentEval.objective) / previousObjective;
        ctl->statusValue[1] = stationarity;
        ++ctl->outer;
        return;
    }
    if (early == NcpEarlyStop::Infeasible) {
        ctl->status = kInfeasible;
        ctl->statusValue[0] = violation;
        ctl->statusValue[1] = earlierViolation;
        ctl->statusValue[2] = ctl->rho;
        ++ctl->outer;
        return;
    }
    if (violation > ctl->tolViol && violation > 0.25 * ctl->previousViolation) ctl->rho = fmin(ctl->rho * 4, ctl->rhoMax);
    ctl->previousViolation = violation;
}

// After the terms at the new penalty: the time limit, the step budget, the iteration cap.
__global__ void outerFinishKernel(NcpControl* ctl, Conditions c) {
    if (isTerminal(ctl->status)) {
        cudaGraphSetConditional(c.outer, 0);
        return;
    }
    if (ctl->timeLimitNs > 0 && deviceNanoseconds() - ctl->startNs > ctl->timeLimitNs) {
        ctl->status = kTimeLimit;
        cudaGraphSetConditional(c.outer, 0);
        return;
    }
    ++ctl->outer;
    if (ctl->maxLmSteps > 0 && ctl->lmSteps >= ctl->maxLmSteps) {
        ctl->status = kStepBudget;
        cudaGraphSetConditional(c.outer, 0);
        return;
    }
    if (ctl->outer >= ctl->maxOuter) {
        ctl->status = kIterationCap;
        cudaGraphSetConditional(c.outer, 0);
        return;
    }
    cudaGraphSetConditional(c.outer, 1);
}

// --- building the graph
// Appends nodes to one graph in order: captured launch sequences and conditional nodes.
class GraphSequence {
public:
    GraphSequence(cudaGraph_t graph, cudaStream_t stream) : graph_(graph), stream_(stream) {}
    // A failure between beginCapture and endCapture (a check that throws) must not leave the stream
    // capturing: every later use of it, including the next solve's, would fail. The graph being
    // captured into belongs to SolveGraph, which destroys it.
    ~GraphSequence() {  // cleanup: never throws
        if (!capturing_) return;
        cudaGraph_t captured = nullptr;
        (void)cudaStreamEndCapture(stream_, &captured);
        cudaGetLastError();  // the invalidated capture's error is this solve's, already reported
    }
    GraphSequence(const GraphSequence&) = delete;
    GraphSequence& operator=(const GraphSequence&) = delete;
    cudaGraphConditionalHandle newCondition() {
        cudaGraphConditionalHandle h;
        check(cudaGraphConditionalHandleCreate(&h, graph_, 0, 0), "ncp-gpu graph: condition");
        return h;
    }
    // The launches enqueued on the stream between beginCapture and endCapture follow the last node.
    void beginCapture() {
        const cudaError_t begun = cudaStreamBeginCaptureToGraph(stream_, graph_, tail_.data(), nullptr, tail_.size(), cudaStreamCaptureModeThreadLocal);
        capturing_ = begun == cudaSuccess;  // before the check, which may throw
        check(begun, "ncp-gpu graph: capture");
    }
    void endCapture() {
        cudaStreamCaptureStatus status;
        const cudaGraphNode_t* deps = nullptr;
        size_t n = 0;
        check(cudaStreamGetCaptureInfo(stream_, &status, nullptr, nullptr, &deps, nullptr, &n), "ncp-gpu graph: capture info");
        tail_.assign(deps, deps + n);
        cudaGraph_t captured;
        capturing_ = false;
        check(cudaStreamEndCapture(stream_, &captured), "ncp-gpu graph: end capture");
    }
    // A WHILE or IF node after the last node; returns its body graph.
    cudaGraph_t addConditional(cudaGraphConditionalHandle h, cudaGraphConditionalNodeType type) {
        cudaGraphNodeParams p{};
        p.type = cudaGraphNodeTypeConditional;
        p.conditional.handle = h;
        p.conditional.type = type;
        p.conditional.size = 1;
        cudaGraphNode_t node;
        check(cudaGraphAddNode(&node, graph_, tail_.data(), nullptr, tail_.size(), &p), "ncp-gpu graph: conditional node");
        tail_.assign(1, node);
        return p.conditional.phGraph_out[0];
    }

private:
    cudaGraph_t graph_;
    cudaStream_t stream_;
    std::vector<cudaGraphNode_t> tail_;
    bool capturing_ = false;
};

// The first launch of a solve, or one that continues it after the contact capacity grew.
enum class LaunchKind { Start, Resume };

class SolveGraph {
public:
    SolveGraph(NcpDevice::State& m, const NcpOptions& options) : m_(m), options_(options) {
        check(cudaGraphCreate(&graph_, 0), "ncp-gpu graph");
    }
    ~SolveGraph() {  // cleanup: never throws
        if (exec_) cudaGraphExecDestroy(exec_);
        if (graph_) cudaGraphDestroy(graph_);
    }
    SolveGraph(const SolveGraph&) = delete;
    SolveGraph& operator=(const SolveGraph&) = delete;

    void build(LaunchKind kind);
    void launch() {
        check(cudaGraphInstantiate(&exec_, graph_, 0), "ncp-gpu graph: instantiate");
        check(cudaGraphLaunch(exec_, m_.stream), "ncp-gpu graph: launch");
    }

private:
    void buildOuter(cudaGraph_t body);
    void buildInner(cudaGraph_t body);
    void buildStep(cudaGraph_t body);
    void buildAttempt(cudaGraph_t body);
    void buildEvaluate(cudaGraph_t body);

    NcpDevice::State& m_;
    const NcpOptions& options_;
    cudaGraph_t graph_ = nullptr;
    cudaGraphExec_t exec_ = nullptr;
    Conditions c_{};
    // The thread's cancel flag as the device reads it: fetched before any capture, since its first
    // use on a thread allocates pinned host memory, which a capture forbids.
    const volatile int* abortFlag_ = abortFlagDevice();
};

void SolveGraph::build(LaunchKind kind) {
    GraphSequence top(graph_, m_.stream);
    c_.outer = top.newCondition();
    top.beginCapture();
    if (kind == LaunchKind::Start) {
        startKernel<<<1, 1, 0, m_.stream>>>(m_.control(), c_);
        m_.enqueueSamplePoints(RowSet::Current);
        m_.enqueueRefresh();
        startCheckKernel<<<1, 1, 0, m_.stream>>>(m_.control(), m_.counts(), c_);
        m_.enqueueLinearize(RowSet::Current);
    } else {
        resumeKernel<<<1, 1, 0, m_.stream>>>(c_);
    }
    top.endCapture();
    buildOuter(top.addConditional(c_.outer, cudaGraphCondTypeWhile));
}

void SolveGraph::buildOuter(cudaGraph_t body) {
    GraphSequence seq(body, m_.stream);
    c_.inner = seq.newCondition();
    c_.outerTail = seq.newCondition();
    seq.beginCapture();
    outerBeginKernel<<<1, 1, 0, m_.stream>>>(m_.control(), c_);
    seq.endCapture();
    buildInner(seq.addConditional(c_.inner, cudaGraphCondTypeWhile));
    seq.beginCapture();
    outerTailCheckKernel<<<1, 1, 0, m_.stream>>>(m_.control(), c_);
    seq.endCapture();
    GraphSequence tail(seq.addConditional(c_.outerTail, cudaGraphCondTypeIf), m_.stream);
    tail.beginCapture();
    m_.enqueueKkt();
    m_.enqueueMultipliers();
    outerDecideKernel<<<1, 1, 0, m_.stream>>>(m_.control());
    m_.enqueueTerms();
    tail.endCapture();
    seq.beginCapture();
    outerFinishKernel<<<1, 1, 0, m_.stream>>>(m_.control(), c_);
    seq.endCapture();
}

void SolveGraph::buildInner(cudaGraph_t body) {
    GraphSequence seq(body, m_.stream);
    c_.refresh = seq.newCondition();
    c_.step = seq.newCondition();
    seq.beginCapture();
    m_.enqueueMotion();
    innerHeadKernel<<<1, 1, 0, m_.stream>>>(m_.control(), c_, abortFlag_);
    seq.endCapture();
    {
        GraphSequence refresh(seq.addConditional(c_.refresh, cudaGraphCondTypeIf), m_.stream);
        c_.afterRefresh = refresh.newCondition();
        refresh.beginCapture();
        m_.enqueueRefresh();
        refreshCheckKernel<<<1, 1, 0, m_.stream>>>(m_.control(), m_.counts(), c_);
        refresh.endCapture();
        GraphSequence evaluate(refresh.addConditional(c_.afterRefresh, cudaGraphCondTypeIf), m_.stream);
        evaluate.beginCapture();
        m_.enqueueLinearize(RowSet::Current);
        evaluate.endCapture();
    }
    buildStep(seq.addConditional(c_.step, cudaGraphCondTypeIf));
}

void SolveGraph::buildStep(cudaGraph_t body) {
    GraphSequence seq(body, m_.stream);
    c_.attempt = seq.newCondition();
    seq.beginCapture();
    m_.enqueueAssemble();
    afterAssembleKernel<<<1, 1, 0, m_.stream>>>(m_.control(), m_.counts(), c_);
    seq.endCapture();
    buildAttempt(seq.addConditional(c_.attempt, cudaGraphCondTypeWhile));
    seq.beginCapture();
    innerEndKernel<<<1, 1, 0, m_.stream>>>(m_.control(), m_.counts(), c_);
    seq.endCapture();
}

void SolveGraph::buildAttempt(cudaGraph_t body) {
    GraphSequence seq(body, m_.stream);
    c_.pcg = seq.newCondition();
    c_.halve = seq.newCondition();
    c_.evaluate = seq.newCondition();
    seq.beginCapture();
    m_.enqueuePcgSetup(c_.pcg);
    seq.endCapture();
    GraphSequence pcg(seq.addConditional(c_.pcg, cudaGraphCondTypeWhile), m_.stream);
    pcg.beginCapture();
    m_.enqueuePcgIteration(c_.pcg);
    pcg.endCapture();
    seq.beginCapture();
    m_.enqueuePcgFinish();
    m_.enqueueStepSummary();
    afterStepKernel<<<1, 1, 0, m_.stream>>>(m_.control(), c_);
    seq.endCapture();
    GraphSequence halve(seq.addConditional(c_.halve, cudaGraphCondTypeWhile), m_.stream);
    halve.beginCapture();
    m_.enqueueTrial();
    m_.enqueueSwept();
    if (options_.auditCoverage) m_.enqueueCoverageAudit();
    halveKernel<<<1, 1, 0, m_.stream>>>(m_.control(), c_);
    halve.endCapture();
    seq.beginCapture();
    afterSweptKernel<<<1, 1, 0, m_.stream>>>(m_.control(), c_);
    seq.endCapture();
    buildEvaluate(seq.addConditional(c_.evaluate, cudaGraphCondTypeIf));
    seq.beginCapture();
    attemptEndKernel<<<1, 1, 0, m_.stream>>>(m_.control(), m_.counts(), c_);
    seq.endCapture();
}

void SolveGraph::buildEvaluate(cudaGraph_t body) {
    GraphSequence seq(body, m_.stream);
    c_.pin = seq.newCondition();
    c_.accept = seq.newCondition();
    seq.beginCapture();
    m_.enqueueLinearize(RowSet::Trial);
    m_.enqueueRegained();
    trialDecideKernel<<<1, 1, 0, m_.stream>>>(m_.control(), c_);
    seq.endCapture();
    {
        GraphSequence pin(seq.addConditional(c_.pin, cudaGraphCondTypeIf), m_.stream);
        c_.relinearize = pin.newCondition();
        pin.beginCapture();
        m_.enqueuePinLost();
        pinDecideKernel<<<1, 1, 0, m_.stream>>>(m_.control(), c_);
        pin.endCapture();
        GraphSequence again(pin.addConditional(c_.relinearize, cudaGraphCondTypeIf), m_.stream);
        again.beginCapture();
        m_.enqueueLinearize(RowSet::Current);
        again.endCapture();
    }
    GraphSequence accept(seq.addConditional(c_.accept, cudaGraphCondTypeIf), m_.stream);
    accept.beginCapture();
    m_.enqueueAccept();
    accept.endCapture();
}

// The solve's parameters and starting state, written once before the launch (host preparation).
NcpControl initialControl(const NcpCall& call) {
    const NcpSolveParams& sp = call.solve;
    if (sp.stallWindow + 1 > kHistory || sp.infeasibleWindow + 1 > kHistory) throw DeviceError("ncp-gpu: stall / infeasible window longer than the history");
    NcpControl c{};
    c.tolViol = sp.tolViol;
    c.tolViolContact = sp.tolViolContact;
    c.tolStat = sp.tolStat;
    c.tolObjRel = sp.tolObjRel;
    c.objFloor = sp.objFloor;
    c.rhoMax = sp.rhoMax;
    c.progressTol = sp.progressTol;
    c.infeasibleRho = sp.infeasibleRho;
    c.delta = call.problem.delta;
    c.sweptEps = kSweptEps * (call.move.thr > 0 ? call.move.thr : 1.0);
    // SPD with null directions of the Gauss-Newton matrix (sliding modes without anchor), resolvable
    // in double. Part of the step model (dDd).
    c.epsRelative = 1e-10;
    c.refreshAt = kRefreshAt;
    c.innerTol0 = kInnerTol0;
    c.pcgTolerance = kPcgTolerance;
    c.maxLmSteps = sp.maxLmSteps;
    c.timeLimitNs = static_cast<long long>(sp.timeLimitMs * 1e6);
    c.maxOuter = sp.maxOuter;
    c.stallWindow = sp.stallWindow;
    c.infeasibleWindow = sp.infeasibleWindow;
    c.innerMax = kInnerMax;
    c.maxAttempts = kMaxAttempts;
    c.maxHalvings = kMaxHalvings;
    c.pcgMaxIterations = kPcgMaxIterations;
    c.rho = sp.rho0;
    c.mu = 1e-4;
    c.nu = 2;
    c.eps = 1e-9;
    c.firstGradient = -1;
    c.previousViolation = 1e300;
    c.status = kRunning;
    return c;
}

const char* capacityName(int kind) {
    switch (kind) {
    case kCapRegistrations: return "grid registrations";
    case kCapCandidates: return "contact candidates";
    case kCapContactRows: return "contact rows";
    case kCapActiveEntries: return "active entries";
    case kCapLiveRows: return "live contact rows";
    default: return "?";
    }
}

std::string statusText(const NcpControl& c, const NcpSolveParams& sp) {
    char b[300];
    switch (c.status) {
    case kConverged: return "CONVERGED";
    case kStoppedProgress:
        std::snprintf(b, sizeof b, "STOPPED: feasible, objective changed by %.2g (<= %.2g) over the last outer iteration; not KKT (stationarity %.3g)", c.statusValue[0],
                      sp.progressTol, c.statusValue[1]);
        return b;
    case kStalledFeasible:
        std::snprintf(b, sizeof b, "UNRESOLVED: stalled feasible (stationarity %.3g > %.3g or complementarity %.3g > %.3g, objective flat)", c.statusValue[0], sp.tolStat,
                      c.statusValue[1], sp.tolViol);
        return b;
    case kInfeasible:
        std::snprintf(b, sizeof b, "UNRESOLVED: infeasible (violation %.3g, %.3g %d outer iterations earlier, rho %.3g)", c.statusValue[0], c.statusValue[1],
                      sp.infeasibleWindow, c.statusValue[2]);
        return b;
    case kStepBudget: return "UNRESOLVED: step budget (" + std::to_string(c.lmSteps) + " LM steps)";
    case kIterationCap: return "UNRESOLVED: iteration cap (" + std::to_string(sp.maxOuter) + " outer iterations)";
    case kTimeLimit: return "UNRESOLVED: time limit";
    case kCapacity:
        std::snprintf(b, sizeof b, "UNRESOLVED: capacity exceeded (%s: %.0f needed, %.0f allocated); the last accepted state is kept",
                      capacityName(static_cast<int>(c.statusValue[0])), c.statusValue[1], c.statusValue[2]);
        return b;
    default: return "UNRESOLVED: ended in state " + std::to_string(c.status);
    }
}

KnitProblemScore scoreOf(const NcpEvalSummary& s, double motion, double delta) {
    KnitProblemScore score;
    score.objective = s.objective;
    for (int f = 0; f < kRowFamilies; ++f) {
        score.maxViol[f] = s.maxViol[f];
        score.nViol[f] = s.nViol[f];
    }
    score.degenerate = s.degenerate;
    score.notOwned = s.notOwned;
    score.strainMax = s.strainMax;
    score.strainRms = s.strainRms;
    score.coverage = motion;
    score.covered = motion <= delta * (1 + 1e-12);
    return score;
}

}  // namespace

std::string NcpStats::line() const {
    char b[1000];
    std::snprintf(b, sizeof b,
                  "%lld unknowns; %lld LM steps (%lld rejected: %lld merit, %lld swept, %lld ownership; %lld swept unresolved; %lld controls pinned; "
                  "%lld trials regaining rows); %lld refreshes, contact rows %lld at the start, max %lld, capacity %lld; %lld PCG solves, %lld iterations "
                  "(%lld capped, %lld non-SPD blocks, %lld breakdowns); ms: setup %.0f, graph %.0f, iterate %.0f, regrow %.0f, download %.0f; transfers up %.1f MB, down %.2f MB; "
                  "device peak %.0f MB; %lld capacity regrows, %lld launches",
                  unknowns, lmSteps, rejected, rejectedMerit, sweptRejected, ownershipRejected, sweptUnresolved, ownershipPinned, regainedTrials, refreshes,
                  initialContactRows, maxContactRows, capacityContactRows, solves, pcgIterations, pcgCapped, pcgNonSpdBlocks, pcgBreakdowns, msSetup, msGraph, msIterate, msRegrow,
                  msDownload, uploadBytes / 1048576.0, downloadBytes / 1048576.0, devPeakMB, regrows, launches);
    return b;
}

// The memory report of one solve: bytes per buffer group (MiB, binary), the lists they hold and their
// capacities, and the tracked peaks.
std::string NcpStats::memoryLine() const {
    std::string s = "memory [MiB]:";
    long long total = 0;
    char b[400];
    for (const NcpMemoryFamily& f : memory) {
        std::snprintf(b, sizeof b, " %s %.1f;", f.name.c_str(), f.bytes / 1048576.0);
        s += b;
        total += f.bytes;
    }
    std::snprintf(b, sizeof b, " total %.1f; tracked peak after setup %.1f, solve %.1f", total / 1048576.0, setupPeakMB, devPeakMB);
    s += b;
    std::snprintf(b, sizeof b, "; samples %lld, general rows %lld, registrations %lld; start: candidates %lld, contact rows %lld, active contact rows %lld",
                  samples, generalRows, registrations, initialCandidates, initialContactRows, initialActiveContactRows);
    s += b;
    std::snprintf(b, sizeof b, "; max: candidates %lld, contact rows %lld, active rows %lld, active entries %lld; capacity: contact rows %lld",
                  maxCandidates, maxContactRows, maxActiveRows, maxActiveEntries, capacityContactRows);
    s += b;
    std::snprintf(b, sizeof b, "; resident: near rows start %lld, max %lld; live rows max %lld; carried from beyond the cutoff max %lld; live capacity %lld",
                  initialNearRows, maxNearRows, maxLiveRows, maxCarriedFar, capacityLiveRows);
    s += b;
    return s;
}

NcpEarlyStop ncpEarlyStop(const NcpSolveParams& sp, const OuterProgress& o) {
    return earlyStopRule(sp.progressTol, sp.infeasibleWindow, sp.infeasibleRho, o.feasible, o.objPrev, o.obj, o.rho, o.viol, o.violBack);
}

bool solveNcp(const NcpCall& call, NcpSolveReport& report, NcpStats* stats, std::string* err) try {
    const auto start = Clock::now();
    report = NcpSolveReport{};
    resetPeak();
    const MemStats before = memStats();
    NcpDevice device(call);
    NcpDevice::State& m = device.state();
    NcpStats& st = device.stats();
    st.unknowns = 3LL * m.N;
    st.msSetup = elapsedMs(start);
    if (m.N == 0) {
        report.converged = true;
        report.status = "CONVERGED";
    } else {
        const NcpControl initial = initialControl(call);
        check(cudaMemcpyAsync(m.control(), &initial, sizeof initial, cudaMemcpyHostToDevice, m.stream), "ncp-gpu control");
        // Launches until the solve ends; a refresh that outgrew the contact capacity (rows) grows it within
        // the memory cap and continues the same solve (NcpOptions::maxRegrows). The host decides only this.
        NcpControl c;
        LaunchKind kind = LaunchKind::Start;
        for (;;) {
            auto t = Clock::now();
            SolveGraph graph(m, call.options);
            graph.build(kind);
            st.msGraph += elapsedMs(t);
            t = Clock::now();
            graph.launch();
            waitForStream("ncp-gpu iteration", m.stream);
            st.msIterate += elapsedMs(t);
            ++st.launches;
            m.controlBuffer.download(&c, sizeof c);
            NcpCounts counts;
            m.countsBuffer.download(&counts, sizeof counts);
            // LiveRows never regrows (its capacities are declared per call; an overflow is a result).
            if (c.status != kCapacity || counts.overflow != kCapContactRows || st.regrows >= call.options.maxRegrows) break;
            if (call.options.contactStorage == ContactStorage::LiveRows) break;
            knit::cancelCheckpoint();
            t = Clock::now();
            m.regrowContactRows(counts.overflowNeeded);  // throws DeviceBudgetExceeded when the cap leaves no room
            st.msRegrow += elapsedMs(t);
            ++st.regrows;
            c.status = kRunning;
            for (double& v : c.statusValue) v = 0;
            c.resume = 1;
            --c.refreshes;  // the overflowed refresh is redone
            m.controlBuffer.upload(&c, sizeof c);
            kind = LaunchKind::Resume;
        }
        auto t = Clock::now();
        if (c.status == kCancelled) throw knit::Cancelled{};
        if (c.status == kContactEntries) throw DeviceError("ncp-gpu: a contact row touches more than 10 free controls");
        report.converged = c.status == kConverged;
        report.stopped = c.status == kStoppedProgress;
        report.budgetHit = c.status == kStepBudget;
        report.status = statusText(c, call.solve);
        report.outer = c.outer;
        report.localSteps = c.localSteps;
        report.backtracks = c.backtracks;
        report.coverageLimited = c.coverageLimited;
        report.passThroughRejected = c.passThroughRejected;
        report.refreshes = c.refreshes;
        report.stationarity = c.stationarity;
        report.rho = c.rho;
        report.score = scoreOf(c.currentEval, device.motionSinceRefresh(), call.problem.delta);
        device.downloadResult(call.controls, call.binding);
        st.lmSteps = c.lmSteps;
        st.rejected = c.rejected;
        st.rejectedMerit = c.rejectedMerit;
        st.sweptRejected = c.sweptRejected;
        st.sweptUnresolved = c.sweptUnresolved;
        st.ownershipRejected = c.ownershipRejected;
        st.ownershipPinned = c.ownershipPinned;
        st.regainedTrials = c.regainedTrials;
        st.refreshes = c.refreshes;
        st.maxContactRows = c.maxContactRows;
        st.maxCandidates = c.maxCandidates;
        st.maxActiveRows = c.maxActiveRows;
        st.maxActiveEntries = c.maxActiveEntries;
        st.maxNearRows = c.maxNearRows;
        st.maxLiveRows = c.maxLiveRows;
        st.maxCarriedFar = c.maxCarriedFar;
        st.solves = c.solves;
        st.pcgIterations = c.pcgIterations;
        st.pcgCapped = c.pcgCapped;
        st.pcgNonSpdBlocks = c.nonSpdBlocks;
        st.pcgBreakdowns = c.pcgBreakdowns;
        st.coverageChecks = static_cast<long long>(c.coverageChecks);
        st.coverageViolations = static_cast<long long>(c.coverageViolations);
        st.msDownload = elapsedMs(t);
    }
    report.ms = elapsedMs(start);
    if (stats) {
        const MemStats after = memStats();
        *stats = st;
        stats->uploadBytes = after.uploadedBytes - before.uploadedBytes;
        stats->downloadBytes = after.downloadedBytes - before.downloadedBytes;
        stats->devPeakMB = after.peak / 1048576.0;
    }
    return true;
} catch (const DeviceError& e) {
    return reportFailure(e, err);
}

}  // namespace knit::gpu
