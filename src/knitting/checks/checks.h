#pragma once
// Geometric checks. Numerical contract:
//   PASS       predicate holds after subtracting every stated allowance
//   FAIL       predicate violated after adding every allowance
//   UNRESOLVED within the allowances, or a method precondition was not met
// Allowances are per quantity and scaled by L_ref (CurveSet::extent()).
//
// Two groups:
//   generator checks  apply to an unedited generated patch only (periodicity,
//                     identical interior signatures)
//   curve-set checks  apply to any CurveSet (regularity, curvature, crossings vs
//                     the expected signature, clearance)
#include <string>
#include <vector>

#include "knitting/checks/interval.h"
#include "model/parallel.h"
#include "model/params.h"
#include "model/status.h"
#include "model/yarn_path.h"

namespace knit {

struct Allowances {
    double Lref = 1;     // [mm]
    double etaX = 0;     // position allowance [mm]
    double rel = 1e-12;  // relative allowance for derived quantities
    static Allowances forExtent(double Lref) { return {Lref, 1e-12 * Lref, 1e-12}; }
};

// ---------------------------------------------------------------- curvature ----
struct CurvatureResult {
    int yarnId = 0;
    double speedLo = 0;      // certified lower bound of |C'| over the whole yarn
    double kappaHi = 0;      // certified upper bound of curvature
    double rhoLo = 0;        // = 1/kappaHi
    double rhoEst = 0;       // sampled minimum radius (estimate)
    double uRhoEst = 0;      // where the sampled minimum occurs
    double kappaWitnessLo = 0;  // best sampled curvature lower bound (for FAIL)
    double uWitness = 0;
    Status regular = Status::Unresolved;
    Status noFold = Status::Unresolved;  // rho >= r
    int leaves = 0;
    // Per span k (u in [k, k+1]): certified curvature upper bound over the span (max
    // over its leaves; 1e300 where the speed bound failed). Used for the local
    // same-yarn window of the clearance check (the bound over an intervening arc).
    std::vector<double> spanKappaHi;
};
// exec: spans in parallel, merged in span order (identical result).
CurvatureResult checkCurvature(const YarnPath& y, double radius, const Allowances& al, const ExecPolicy& exec = {});

// ---------------------------------------------------------------- crossings ----
struct Crossing {
    int yarnA = 0, yarnB = 0;           // yarnA <= yarnB; same yarn => uA < uB
    Interval boxA, boxB;                // Krawczyk root box (global parameters)
    double uA = 0, uB = 0;              // root estimate (inside the box)
    Vec3 pA, pB;                        // 3D points at the estimate (display)
    Interval depth;                     // z_A - z_B enclosure over the box [mm]
    int front = -1;                     // yarn id in front (+N side), -1 = unresolved
    int stitchA = -1, stitchB = -1;     // index into YarnPath::stitches, -1 = on boundary
    Status status = Status::Pass;       // Unresolved if sign or ownership unresolved
    std::string note;
};

struct UnresolvedItem {
    int yarnA = 0, yarnB = 0;
    double uA0 = 0, uA1 = 0, uB0 = 0, uB1 = 0;
    std::string reason;
};

struct CrossingOptions {
    int maxDepth = 40;
};

struct CrossingResult {
    std::vector<Crossing> crossings;
    std::vector<UnresolvedItem> unresolved;
    long pairsVisited = 0;
};
// Finds all projected crossings (front projection along N) between all yarns,
// including self-crossings of one yarn, on the actual spline geometry.
CrossingResult findCrossings(const CurveSet& cs, const Allowances& al, const CrossingOptions& opt = {});

// Comparison against expected_signature.h.
struct PairSignature {
    int child = -1, parent = -1;    // topology stitch ids (child is drawn through parent)
    int wale = 0, olderCourse = 0;  // child index in its turn, parent turn (flat: wale, older course)
    std::vector<int> crossingIdx;   // into CrossingResult::crossings, sorted along newer
    Status status = Status::Pass;
    std::string detail;
};
struct SignatureResult {
    std::vector<PairSignature> pairs;
    std::vector<int> unexpected;    // crossings not belonging to any intermeshed pair
    std::vector<int> overlap;       // DEC parent-parent crossings (allowed with the later parent in front)
    int overlapWrongOrder = 0;
    std::vector<int> layerOverlap;  // rib layers: yarn-neighbour stitches of different faces (allowed with the knit stitch in front)
    int layerWrongOrder = 0;
    Status status = Status::Pass;
    double minDepthGapLo = 0;       // min over crossings of |depth| lower bound [mm]
    Status depthGapVsClearance = Status::Pass;  // |depth| >= D + c (necessary condition)
    std::string detail;
};
SignatureResult checkSignature(const CurveSet& cs, const CrossingResult& cr, const PatternParams& pp);

// ---------------------------------------------------------------- clearance ----
struct ClearanceIssue {
    int yarnA = 0, yarnB = 0;
    double uA = 0, uB = 0;
    Vec3 pA, pB;
    double distance = 0;   // polyline distance [mm]
    double threshold = 0;  // D + c, or D in the tight-bend case
    Status status = Status::Fail;
    bool penetrating = false;
    std::string reason;
};
struct ClearanceRegion {
    int yarnA = 0, yarnB = 0;
    double uA0 = 0, uA1 = 0, uB0 = 0, uB1 = 0;
    double worstDistance = 0, threshold = 0;
    Status status = Status::Fail;
    bool penetrating = false;
    int worstIssue = -1;
};
struct ClearanceResult {
    std::vector<ClearanceIssue> issues;  // all FAIL / UNRESOLVED pairs
    std::vector<ClearanceRegion> regions;
    double minInterYarn = 1e300;         // min polyline distance between different yarns
    double minSelfNonlocal = 1e300;      // min polyline distance over nonlocal same-yarn pairs
    double deltaD = 0;                   // distance allowance used [mm]
    long segmentPairs = 0;
    long schurCertified = 0;             // same-yarn pairs certified by the chord bound
    Status status = Status::Pass;
};
struct ClearanceInput {
    const CurveSet& cs;
    const std::vector<Polyline>& polys;           // the yarns sampled (runSampling)
    const std::vector<CurvatureResult>& curvature;  // per yarn (runCurvature)
};
// exec: rows of candidate chunk pairs in parallel, merged in row order (identical issues, counters
// and regions).
ClearanceResult checkClearance(const ClearanceInput& in, const PatternParams& pp, const Allowances& al, const ExecPolicy& exec = {});

// ---------------------------------------------------------- local signature ----
// Projection-local signature for curve sets on a curved surface (CurveKind::Surface),
// where one global front projection is meaningless. For every stitch c drawn through
// parent(s) p:
//   * the child's own span and the parents' neighbourhood (from the stitch before the
//     first parent to the stitch after the last, along the parent yarn) are extracted
//     as EXACT sub-B-splines (the controls of their span ranges),
//   * transformed into the child's frame (origin, U, V, N; N = the surface normal),
//   * and the existing crossing search runs on this two-yarn mini set, so portions of
//     the SAME yarn (a spiral interlocking with its previous turn) become an ordinary
//     pair. The crossings are compared with the expected signature of the child's op
//     (checkSignature on the mini set): the child must cross exactly its parent(s),
//     with the expected count and depth order, and nothing else of the neighbourhood.
//     Crossings of the parent neighbourhood with itself are reported, not counted:
//     they say nothing about the child's interlock, and where the surface bends
//     within the neighbourhood the projection makes neighbouring heads overlap.
struct LocalPairResult {
    int child = -1, parent = -1;
    Status status = Status::Pass;
    int crossings = 0, unexpected = 0, overlap = 0;
    int terminal = 0;  // 1: a CLOSE drawstring pass (child = the gathered loop, parent = -1)
    int neighbourhoodSelf = 0;  // parent-row crossings among themselves (reported, not counted)
    std::string detail;
};
struct LocalSignatureResult {
    std::vector<LocalPairResult> pairs;
    Status status = Status::Pass;
    int pass = 0, fail = 0, unresolved = 0;
    int neighbourhoodSelf = 0;  // total of the per-pair diagnostic
    double ms = 0;
    std::string detail;
    // Crossing attribution (applyCrossingAttribution): how many checks it re-evaluated.
    int attributionReevaluated = 0;
};

// Crossing attribution. A crossing of a child's yarn with a stitch p of the
// parent row counts for the child whose ownership span [uBegin, uEnd) contains it if p is that
// child's topology parent (the original rule). A crossing inside ANOTHER stitch's span that is not
// that stitch's parent counts for the unique stitch of the yarn that has p as a parent and whose span
// lies within `window` B-spline spans of it (u in [uBegin - window, uEnd + window)); two such
// stitches: ambiguous (the check is UNRESOLVED); none: it stays with no pair (unexpected for the
// stitch whose span contains it). Crossings outside every stitch span (beyond a yarn's stitches) are
// not attributed, as before: the window moves crossings across a boundary between two owners only.
// kAttributionWindow: the largest excursion of a verified-correct crossing outside its span on the
// teddy ear, body and the man body is 1.76 spans (one stitch in the ear's CLOSE row, which fails for
// other crossings anyway); without it 1.10. 2 spans covers both; results are the same for windows
// 1.5 .. 4 on those fixtures. The window must stay below the shortest stitch span.
constexpr double kAttributionWindow = 2.0;
struct CrossingOwner {
    int owner = -1;          // stitch id (attributeCrossing) or CLOSE loop index (attributeCloseCrossing)
    bool ambiguous = false;
};
CrossingOwner attributeCrossing(const KnitTopology& topo, int yarn, double u, int parentStitch, double window);
// CLOSE: a drawstring crossing at strand parameter u with gathered head `headStitch`: the loop whose
// nominal strand range [first, second] contains u if that loop gathers this head; else, if u lies in
// another loop's range, the unique loop of this head whose range lies within `window` spans; outside
// every range (the strand beyond its first / last loop): not attributed, as before.
CrossingOwner attributeCloseCrossing(const TerminalRecord& tr, int headStitch, double u, double window);
// The local signature (computed with ownership spans only, CPU or device) corrected for the
// crossing attribution: re-evaluates every check that can change (failing children, their yarn
// neighbours, the first and last stitch of each yarn, failing CLOSE loops and their neighbouring
// loops) in the same frame, neighbourhood and expected sequences. A crossing claimed by two checks
// makes both UNRESOLVED. window <= 0: unchanged.
void applyCrossingAttribution(const CurveSet& cs, LocalSignatureResult& sig, double window);
// exec: stitches in parallel, accumulated in stitch order (identical result).
// window: crossing attribution (0: ownership spans only, the original rule).
LocalSignatureResult checkLocalSignature(const CurveSet& cs, const PatternParams& pp, const ExecPolicy& exec = {},
                                         double window = kAttributionWindow);

// --------------------------------------------------------------- generator ----
struct GeneratorChecks {
    double maxPeriodicityDeviation = 0;  // [mm] control points vs translated template
    Status periodicity = Status::Pass;
    Status interiorSignatures = Status::Pass;
    std::string detail;
};

// ------------------------------------------------------------------ report ----
struct StageTimes {
    double samplingMs = 0, curvatureMs = 0, crossingsMs = 0, clearanceMs = 0;
};
struct CheckLine {
    std::string name;
    Status status = Status::Pass;
    std::string detail;
};
struct ValidationReport {
    Status verdict = Status::Pass;
    std::vector<CheckLine> lines;
    std::vector<Polyline> polylines;
    std::vector<CurvatureResult> curvature;
    CrossingResult crossings;
    SignatureResult signature;
    ClearanceResult clearance;
    GeneratorChecks generator;
    LocalSignatureResult localSignature;  // CurveKind::Surface only
    StageTimes times;
    Allowances allowances;
    bool generatorChecksApplied = false;
    // Produced by the device (float32 evaluation of these checks): a diagnostic verdict, not the
    // certified CPU one.
    bool gpuDiagnostic = false;
};

struct ValidateOptions {
    bool generatorChecks = true;  // false for edited / relaxed curve sets
    ExecPolicy exec;              // curvature, local signature, clearance (identical results)
    double attributionWindow = kAttributionWindow;  // local signature crossing attribution (0: spans only)
};

// Time spent in validate() on the calling thread, by sub-check (thread-local; a part's settle, repair
// and gate (a) run on one thread). The scene pipelines read it before and after each stage. Profiling
// only; never affects results.
struct ValidateStats {
    long long calls = 0;
    double ms = 0, msCurvature = 0, msSampling = 0, msSignature = 0, msClearance = 0;
    // The calls and times added since `start`, an earlier reading of the same thread's stats.
    ValidateStats since(const ValidateStats& start) const {
        return {calls - start.calls, ms - start.ms, msCurvature - start.msCurvature, msSampling - start.msSampling,
                msSignature - start.msSignature, msClearance - start.msClearance};
    }
};
ValidateStats& validateStats();

// Runs every applicable check. `polylines` are produced at pp.sampleTolerance.
ValidationReport validate(const CurveSet& cs, const PatternParams& pp, const ValidateOptions& opt = {});

// The verdict and the report lines from the checks' results (the device's validate calls it too).
void summarize(ValidationReport& r, const CurveSet& cs, const PatternParams& pp, const ValidateOptions& opt = {});

std::string formatReport(const ValidationReport& r, const CurveSet& cs, const PatternParams& pp, bool verbose);
// Every result field of a report as text (doubles as exact hex floats; timings and the sample
// polylines' points excluded): equal text <=> equal results. Saved knits store its digest.
std::string validationSnapshot(const ValidationReport& r);

// ------------------------------------------------------------- flat swatches ----
// The flat reference verdict of a knit / purl pattern (PatternReference, model/params.h): validator on an n x n flat patch of `p` at pp's D, w, h, A, c, template and transition.
PatternReference flatPatternReference(const PatternParams& pp, StitchPattern p, int n = 8);

// ------------------------------------------------------------------------------
// Projection-independent topology criterion for a motion between two curve sets
// with identical control-point structure: the straight-line interpolation of all
// control points. Every curve point moves by at most delta = max control
// displacement (a B-spline point is a convex combination of its controls), so
// for any pair of points d(t) >= (d0 + d1)/2 - delta. If both end states pass the
// clearance and curvature checks, every nonlocal pair has d >= D + c at both ends
// (same-yarn pairs inside the Schur window by the chord bound), hence
//   delta < D + c   =>  no centreline passes through another during the motion;
//   delta <= c      =>  the tubes do not even touch during the motion.
// Same-yarn pairs closer than (pi/2) * delta in arc length are covered only by the
// end-state curvature bound (a loop that short would need curvature >= 4/delta).
// Unlike the projected crossing signature this does not depend on a projection;
// it certifies that no yarn passes through yarn, which is what preserves the
// interlocking (open row ends aside: the patch is anchored).
struct MotionCertificate {
    bool endStatesValid = false;   // clearance and curvature PASS at both ends
    double maxDisplacement = 0;    // delta [mm]
    bool noPassThrough = false;
    bool noTubeContact = false;
    double localScale = 0;         // (pi/2) delta [mm]
    std::string detail;
};
MotionCertificate certifyStraightMotion(const CurveSet& a, const ValidationReport& ra, const CurveSet& b,
                                        const ValidationReport& rb, const PatternParams& pp);

} // namespace knit
