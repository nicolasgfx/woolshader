#pragma once
// The written reference of the knitting problem (knitting/relax/knit_problem.h): its rows and its score on
// the CPU, in double precision. Only tests use it: the device's rows must match it
// (gpu_ncp_gpu_rows_match_cpu_*), and its gradients are checked by finite differences
// (ncp_jacobian_fd_*).
//
// Two operations are kept apart:
//   linearize      the model rows (values and gradients at the current controls, contact
//                  closest-point parameters s, t and normal n frozen) used ONLY to propose steps.
//                  CHANGE (ncp-gpu): the band row is evaluated at the current control on the
//                  actual surface (the binding continued from the refresh state, not advanced;
//                  without a binding the projection): z = +-|x - foot|, signed by the surface
//                  normal at the foot, with its exact gradient +-(x - foot) / |x - foot| (earlier:
//                  n . (x - foot) along the interpolated vertex normal of a tangent plane frozen at
//                  the refresh, whose gradient n is not the derivative of the measured value on a
//                  mesh: the solver stalled at stationarity ~1e-3 on the teddy body). scoreGeometry
//                  measures the same z;
//   scoreGeometry  the authoritative evaluation from the actual curves: true segment distances
//                  of every covered pair, the curvature samples, the band measured on the real
//                  surface (continuation for bound controls, without advancing their states),
//                  the cap, degenerate spans, edge strain. Every acceptance and every convergence
//                  claim uses it.
//
// Contact coverage (cumulative motion): a broad phase built at controls Xref with pad thr + 2
// delta finds every pair closer than thr as long as every control stays within delta of Xref
// (a B-spline sample is a convex combination of its controls, so it moves no farther than they
// do). Solvers check the CUMULATIVE displacement against Xref before evaluating or accepting
// a trial, and refresh (rebuild at the current controls) or shorten the trial otherwise.
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "knitting/relax/knit_problem.h"
#include "knitting/relax/settle.h"

namespace knit {

struct RowEntry {
    int var = -1;  // free-control index
    Vec3 g;        // d value / d x_var
};
struct Row {
    RowFamily family = RowFamily::Edge;
    double value = 0;   // Edge: the residual r; constraints: g (<= 0 satisfied)
    int n = 0;
    RowEntry e[kRowMaxEntries];
    uint64_t key = 0;   // identity across linearizations (multipliers): family, yarn / pair, index
    double aux = 0;     // Edge: the current length l (second-order terms of block solvers)
};

// One relax call (KnitNumerics::relax, numerics.h) as a problem. P holds every control (fixed ones
// never change); the binding's states are those at the call's start.
class KnitProblem {
public:
    KnitProblem(const CurveSet& cs, const Controls& P, const Controls& Prest, const ControlMask& freeCtrl,
                const BaseSurface* surface, const SurfaceBinding* bind, const RelaxMoveParams& mp, const KnitProblemParams& kp = {});

    int vars() const { return static_cast<int>(ctrlOf_.size()); }
    std::pair<int, int> ctrlOf(int v) const { return ctrlOf_[v]; }   // (yarn, control)
    int varOf(int y, int c) const { return varOf_[y][c]; }            // -1: fixed
    const RelaxMoveParams& moveParams() const { return mp_; }
    const KnitProblemParams& params() const { return kp_; }
    double D() const { return D_; }
    const Controls& start() const { return P0_; }                     // the cap origin

    // Broad phase at P: candidate segment pairs (legacy rules, pad thr + 2 delta); Xref := P.
    // With a binding, every free control's state continues from P (advanced only on Owned, as
    // the legacy refresh does) and its band plane is taken from the new foot; without one, the
    // plane comes from surface->project. bind may be null (no binding).
    void refresh(const Controls& P, SurfaceBinding* bind);
    // max over free controls of |P - Xref| (the cumulative motion since the last refresh).
    double motionSinceRefresh(const Controls& P) const;
    const Controls& refreshControls() const { return Xref_; }
    long long refreshes() const { return refreshes_; }
    size_t candidates() const { return cand_.size(); }
    // The contact rows of the last refresh (sample pairs, both directions of a point-segment pair) and
    // their distance types, and a row's distance at P: the reference key set of the contact-memory checks.
    const std::vector<std::pair<size_t, size_t>>& contactRowPairs() const { return cand_; }
    const std::vector<unsigned char>& contactRowTypes() const { return candType_; }
    double contactRowDistance(const Controls& P, size_t i, size_t j, int type) const {
        double s, t;
        return contactDistance(P, i, j, type, s, t);
    }

    // Model rows at P (gradients analytic). Only constraint rows whose value is above
    // `activeBelow` (e.g. -0.5) are emitted, edge rows always.
    void linearize(const Controls& P, std::vector<Row>& rows, double activeBelow = -1e300) const;
    // Rows touching variable v only (block solvers).
    void linearizeVar(const Controls& P, int v, std::vector<Row>& rows, double activeBelow = -1e300) const;

    // Authoritative evaluation (see the header comment). tol: the violation counted in nViol.
    KnitProblemScore scoreGeometry(const Controls& P, const SurfaceBinding* bind, double tol = 1e-9) const;

    // Pass-through motion test (stop point 1, class iii): true if moving variable v from `before`
    // to `after` (all other controls equal) makes a covered pair touching v closer than `near`
    // switch sides (the separation direction between the closest points reverses) or touch. A
    // crossing changes sides only through distance 0, so a step that keeps every near pair's
    // direction is not a pass-through of those pairs; the test is conservative (an end-point
    // closest feature can also turn the direction): it may reject steps, never accepts a reversal.
    // P holds `to` at v on return (it is temporarily set to `from`).
    bool passThrough(Controls& P, int v, const Vec3& from, const Vec3& to, double near) const;

private:
    struct Sample {
        int y, k;
        double t;  // parameter in the span
        double w[4];
        double s;
    };
    // The hard row and, with hard thresholds (KnitProblemParams), the target row (n = 0 when the
    // target is met: no objective contribution).
    // type: 0 segment-segment, 1 point-point, 2 point i - segment j (contact model, classified at refresh).
    void contactRows(const Controls& P, size_t i, size_t j, int type, Row& hard, Row* target) const;
    // The contact distance of row (i, j, type) with its closest parameters.
    double contactDistance(const Controls& P, size_t i, size_t j, int type, double& s, double& t) const;
    bool curvatureRows(const Controls& P, int y, int k, int q, Row& hard, Row* target, bool* degenerate) const;
    double hardC() const { return kp_.hardContact > 0 ? kp_.hardContact : mp_.thr; }
    double hardK() const { return kp_.hardKappa > 0 ? kp_.hardKappa : mp_.kappaMax; }
    bool targets() const { return kp_.hardContact > 0 || kp_.hardKappa > 0; }
    void anchorRows(const Controls& P, int v, std::vector<Row>& rows) const;
    // The band's signed distance at variable v's control and its gradient (false: not Owned / no surface).
    bool bandDistance(const Controls& P, int v, double& z, Vec3& u) const;
    Vec3 evalSample(const Controls& P, const Sample& s) const;
    void addSampleTerms(Row& r, const Sample& s, double f, const Vec3& n) const;

    const CurveSet& cs_;
    RelaxMoveParams mp_;
    KnitProblemParams kp_;
    double D_ = 1;
    Controls P0_, Xref_;
    ControlMask free_;
    const BaseSurface* surface_ = nullptr;
    std::vector<std::vector<int>> varOf_;
    std::vector<std::pair<int, int>> ctrlOf_;
    std::vector<std::vector<double>> rest_;
    std::vector<Sample> samples_;
    std::vector<char> active_;
    std::vector<std::pair<size_t, size_t>> cand_;  // contact rows (sample pairs)
    std::vector<unsigned char> candType_;          // their distance type (contactRows)
    std::vector<std::vector<int>> candOfVar_;   // candidate indices touching each variable
    const SurfaceBinding* bind_ = nullptr;      // the refresh's binding (band rows continue from its states)
    long long refreshes_ = 0;
};

} // namespace knit
