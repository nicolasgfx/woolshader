#pragma once
// The knitting relaxation as an explicit constrained problem. One
// KnitNumerics::relax call (numerics.h) defines it: the free controls are the variables, every other
// control stays where it is.
//
//   minimize   (w_L / 2) sum_edges r_e^2,   r_e = (|x_{c+1} - x_c| - rest_e) / D
//   subject to K  kappa(span, q) / kappaMax - 1 <= 0       q over the span's curvature samples
//              C  (thr - dist(segment a, segment b)) / D <= 0
//              B  |z| / bandHalf - 1 <= 0                  Owned controls only
//              M  |x_c - x_c^0|^2 / cap^2 - 1 <= 0         repair: cap > 0
//
// K covers every span with a free support control. A span with |C'| < kDegenerate at a sample is
// degenerate: reported as violated, never treated as satisfied.
//
// C covers the sample segments near each other. Two segments of one yarn count only when the arc gap
// between them reaches the move's window: neighbours share an end point, a distance no motion can raise.
//
// B measures z on the actual surface: z = +-|x - foot|, signed by the surface normal at the foot, the
// foot coming from the binding when there is one and from the projection otherwise.
//
// With hard thresholds (KnitProblemParams::hardContact, hardKappa) C and K are hard at what the next
// validation judges, and the move's thr and kappaMax become soft objective rows.
//
// The device solves this problem (gpu/solver/ncp.h). This header holds what the solver and its callers
// share. Its written CPU reference, class KnitProblem, is in tests/ncp_reference.h: the device's
// rows are tested against it.
#include <cstdint>
#include <cstdio>
#include <string>

#include "model/params.h"

namespace knit {

// Edge, ContactTarget, CurvatureTarget, Anchor: objective residuals; the others: hard rows (<= 0).
enum class RowFamily : uint8_t { Edge, Curvature, Contact, BandUpper, BandLower, Cap, ContactTarget, CurvatureTarget, Anchor };
constexpr int kRowFamilies = 9;
inline bool isHard(RowFamily f) {
    return f != RowFamily::Edge && f != RowFamily::ContactTarget && f != RowFamily::CurvatureTarget && f != RowFamily::Anchor;
}
const char* toString(RowFamily f);

constexpr int kRowMaxEntries = 10;  // contact: <= 5 controls per segment side

struct KnitProblemParams {
    double wL = 1.0;          // edge objective weight (normalized residuals)
    double delta = 0.1;       // [mm] contact coverage: max cumulative control motion between broad phases
    double kDegenerate = 1e-9;  // |C'| below this at a curvature sample: degenerate span
    // Hard thresholds: > 0 makes the C / K rows hard at these values (the ones gate (a) and the next
    // settle round judge: D + c and 1 / (r + c/2) of the validation parameters) and the move's own
    // thr / kappaMax soft targets (objective rows ContactTarget / CurvatureTarget, weight wT).
    // 0: the move's thr / kappaMax are hard.
    double hardContact = 0, hardKappa = 0;
    double wT = 1.0;
    // Shape preservation (Stitch Meshes' yarn relaxation): (wA / 2) |x - x0|^2 / D^2 per free control
    // (three Anchor rows). The edge objective leaves rigid slides of whole loops free; large slides
    // change the projected crossing pattern the local signature checks. 0: off.
    double wA = 0;
    // Curvature rows per span: at t = q / (curvatureSamples - 1), q = 0 .. curvatureSamples - 1.
    // 5 (t in quarters) missed the maximum of a yarn's first span by 13 % on the teddy body (the
    // certified curvature check then found a tight bend); more samples follow the certified bound.
    int curvatureSamples = 5;
    // The contact rows' distance. Segment: sample segment to sample segment. Its value is not smooth
    // in the controls: the segment point is bilinear in its parameter and the end points, and for
    // short, nearly parallel segments the closest parameters jump between the segment ends inside a
    // step (teddy body: a row's actual change matched the linear model up to 1/10 of a step and jumped
    // beyond; the solver crawled at gain ~0.45). Point: sample point to sample point (smooth, convex;
    // overestimates the polyline distance by up to ~L^2 / (2 d), 0.06 mm at the teddy's 0.3-1.3 mm
    // sample spacing: clearance FAIL). PointSegment: every sample point against every candidate segment
    // of the other strand, both directions (vertex-edge contact): the minimizer is unique, so each row
    // is smooth; it equals the segment distance unless both closest points are segment interiors
    // (crossings), where it overestimates by up to L^2 / (8 d) (measured: 0.019 mm on the teddy body ->
    // clearance UNRESOLVED). Hybrid: per candidate, classified at each broad phase by the angle of its
    // two segments: nearly parallel (|cos| > parallelCos) -> point-segment rows in both directions
    // (exact there: the minimum lies at an end point), else one segment row (exact; for crossing
    // segments the minimizer is unique and well conditioned). The candidate rules are the segments'
    // in every model.
    enum class ContactModel { Segment, Point, PointSegment, Hybrid };
    ContactModel contactModel = ContactModel::Segment;
    double parallelCos = 0.9;
    // > 0: a span gets max(samplesPerSpan, ceil(arc length / maxSampleSpacing)) samples (arc length on
    // the call's input curve), so no sample spacing exceeds about this [mm]. Bounds the point contacts'
    // overestimate (L^2 / (2 d)) and the polylines' sagitta where spans are long (cast-on turns, the
    // CLOSE tail: up to 1.4 mm per span on the teddy body). 0: samplesPerSpan everywhere.
    double maxSampleSpacing = 0;
};

// The problem evaluated at the actual curves: every acceptance and convergence claim uses this.
struct KnitProblemScore {
    double objective = 0;
    double maxViol[kRowFamilies] = {};   // max(0, g) per family (Edge: max |r|; targets: max shortfall)
    long long nViol[kRowFamilies] = {};  // g > tol
    long long degenerate = 0;     // curvature samples with |C'| < kDegenerate
    long long notOwned = 0;       // bound free controls whose continuation is not Owned
    double strainMax = 0, strainRms = 0;  // relative edge strain |l - rest| / rest
    double coverage = 0;          // cumulative motion since the last broad phase; > delta: pairs may be missing
    bool covered = true;
    double maxConstraint() const;  // over K, C, B (upper / lower), M
    std::string line() const;
};

// The hard thresholds of a relax call from the next validation's parameters vp:
// contact D_v + c_v + marginC, curvature 1 / (r_v + c_v / 2) / (1 + marginK).
struct HardThresholds {
    double contact = 0, kappa = 0;
};
inline HardThresholds hardThresholds(const PatternParams& vp, double marginC, double marginK) {
    return {vp.yarnDiameter + vp.clearance + marginC, 1.0 / (vp.radius() + 0.5 * vp.clearance) / (1 + marginK)};
}

// How the device solver holds its contact rows.
//   CompactRows  every broad-phase candidate is a row with its multiplier, values and terms (v5, the default).
//   LiveRows     experimental: a row is kept only when its two sample segments are within the near cutoff
//                max(hardC, thr) + 2 delta (or it carries a positive multiplier); a row stores its pair and
//                type only, and the multiplier, values and terms live in a smaller list of live rows (a
//                term in the merit or a positive multiplier). Same problem, same iteration: a row left
//                out has no term, no violation and no swept touch until the next refresh.
enum class ContactStorage { CompactRows, LiveRows };

// Which curvature rows bound a span.
//   Samples                the curvature at KnitProblemParams::curvatureSamples fixed parameters per span (v5, the
//                          default).
//   SamplesAndSpanMaximum  experimental: also one row at the span's curvature maximum, located again at every
//                          evaluation. The fixed samples let a bend sit between two of them (tube-144: radius
//                          0.724 mm at both samples, 0.613 mm between them, against the hard 0.649 mm).
enum class CurvatureRows { Samples, SamplesAndSpanMaximum };

// When a call whose violation shrinks slowly ends.
//   StopWhenSlow    the infeasible rule of NcpSolveParams: at rho >= infeasibleRho, a violation that did not halve over
//                   infeasibleWindow outer iterations ends the call (v5, the default: production's repair fails fast).
//   UseStepBudget   that rule is off; the call runs to convergence or to its step budget. For the stitch-mesh settle,
//                   whose calls at the target diameter keep converging slowly and were each cut at about half the budget.
enum class SlowCall { StopWhenSlow, UseStepBudget };

// The solver's stopping rules. Converged only when the score is feasible (hard families <= tolViol,
// no degenerate span, coverage held) and the stationarity and the objective change are small;
// otherwise UNRESOLVED with the reason.
struct NcpSolveParams {
    double tolViol = 1e-3;     // normalized violation allowance (K, C, B, M)
    // The contact family's own allowance (0: tolViol). Feasibility and the penalty increase measure
    // contact against it and the other hard families against tolViol.
    double tolViolContact = 0;
    double tolStat = 1e-4;     // max |augmented Lagrangian gradient| over the variables
    double tolObjRel = 1e-4;   // relative objective change over `stallWindow` outer iterations
    double objFloor = 1e-8;    // ... relative to max(objFloor, |objective|) (a zero objective is flat)
    int stallWindow = 5;
    double rho0 = 10, rhoMax = 1e8;
    double mu = 1e-6;          // Hessian floor (relative to its trace)
    int maxOuter = 400;        // safety net, reported when hit
    double timeLimitMs = 0;    // 0: none
    // Levenberg-Marquardt steps per call (0: none). Checked before each inner iteration; when spent,
    // the outer iteration's measures are taken and the call ends UNRESOLVED: step budget (unless it
    // converged or an early stop applies). Deterministic, unlike timeLimitMs.
    long long maxLmSteps = 0;
    // Early stops (0 = off). Progress: feasible and the objective changed by at most progressTol
    // (relative, either sign) over the last outer iteration -> STOPPED. Infeasible: not feasible,
    // rho >= infeasibleRho, and the violation did not halve over the last infeasibleWindow outer
    // iterations -> UNRESOLVED: infeasible.
    double progressTol = 0;
    int infeasibleWindow = 0;
    double infeasibleRho = 1e3;
};
enum class NcpEarlyStop { None, Progress, Infeasible };
struct NcpSolveReport {
    bool converged = false;
    bool stopped = false;      // STOPPED by the progress rule: feasible, not KKT
    bool budgetHit = false;    // UNRESOLVED: step budget (NcpSolveParams::maxLmSteps)
    std::string status;        // "CONVERGED", "STOPPED: <reason>" or "UNRESOLVED: <reason>"
    int outer = 0;
    long long localSteps = 0, backtracks = 0, coverageLimited = 0, passThroughRejected = 0;
    long long refreshes = 0;
    double stationarity = 0, rho = 0;
    KnitProblemScore score;
    double ms = 0;
};

} // namespace knit
