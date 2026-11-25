// growYarnFibres: every candidate, from its anchor to an accepted fibre or the test that refused it.
// The tests run in the order of the GPU kernels (fibres.cu kernel A, then growFibre), cheapest first.
#include "fibres/fibre_generator.h"

#include <algorithm>
#include <chrono>
#include <cmath>

#include "fibres/fibre_islands.h"
#include "fibres/fibre_shape.h"
#include "fibres/strands_file.h"
#include "knitting/tube/plies.h"
#include "knitting/tube/yarn_surface.h"
#include "model/cancel.h"
#include "model/hash.h"

namespace knit {

namespace {

constexpr double kDeg = kPi / 180;
constexpr double kOriginOffset = 0.02;   // [mm] the exposure ray's origin above the surface
constexpr double kYarnEndMargin = 0.05;  // [mm] the directed curve and the root stay this far from a strand's ends
constexpr double kDeepContact = 0.2;     // a fibre may lie this share of the yarn radius inside a strand

// What every candidate of a run shares.
struct FibreRun {
    const CurveSet& cs;
    const FibreParams& fp;
    const FibreFrame* frame;
    SurfaceIndex index;
    uint64_t seed = 0;
    double r = 0;          // radius of the supporting strands
    double ownWindow = 0;  // own-yarn range handled with the exact surface
    bool curvedFabric = false;  // the frame gives the fabric's base surface, the ray's length and the back depth
    FibreRun(const CurveSet& cs_, const SurfaceModel& sm, const PatternParams& pp, const FibreParams& fp_, const FibreFrame* frame_)
        : cs(cs_), fp(fp_), frame(frame_), index(cs_, sm), seed(fibreSeed(fp_, frame_ ? frame_->pieceSeed : 0)), r(pp.radius()),
          ownWindow(3 * pp.yarnDiameter), curvedFabric(frame_ && frame_->fabric && frame_->exposureLength > 0) {}
    bool onPlies() const { return frame && frame->plies; }
    // How far q lies in front of the fabric's back; negative: behind it. Flat: above the lowest point
    // of the tubes. Curved: along the base surface's normal, plus the back depth.
    double frontOfBack(const Vec3& q) const {
        if (!curvedFabric) return q.z - index.zBottom();
        const FabricPoint f = frame->fabric(q);
        return dot(q - f.p, normalized(f.n)) + frame->backDepth;
    }
};

// One chunk of candidates of one strand.
struct ChunkRange {
    int strand = 0;
    int64_t begin = 0, end = 0;
};

// One candidate while it is tested.
struct Candidate {
    const YarnSurface* ys = nullptr;
    int strand = 0;                  // the strand's position in the surface model
    int i = 0;                       // index on the strand
    double s = 0, theta = 0;         // the anchor: arc length and ellipse parameter
    FibrePopulation population = FibrePopulation::Directed;
    YarnSurfacePoint sp;                 // the anchor's surface point
    IslandWeight island;             // what the fuzz is there
    Vec3 N;                          // the fabric's outward normal there
    SurfaceIndex::Exclude ownYarn;   // the own yarn near the anchor (handled with the exact surface)
    FibreSet* out = nullptr;         // the chunk's fibres and report
    uint8_t* trace = nullptr;        // FibreFrame::outcomes (null without it)
    bool random() const { return population == FibrePopulation::Random; }
    int yarn() const { return ys->yarnId; }
};

// The fibre of one candidate while it is shaped.
struct GrowingFibre {
    FibreDraw draw;
    int minSegments = 0, segments = 0;  // visible segments: the fewest allowed, the number chosen
    bool curled = false;
    bool crimped = false;
    double radius = 0;
    double buriedDepth = 0;   // of the root below the own surface
    double minTurnDeg = 0;    // smallest turn allowed at an interior vertex
    SmoothCurve curve;
    std::vector<Vec3> pts;    // [0] the root, [1] the curve's start, then the visible vertices
    bool rootOnYarn = true;   // the root's arc length stays inside the strand
    double firstChordElevation = 90;  // first visible chord against the surface [deg]
    double emergence = -1;    // the surface crossing on segment 1, as lerp(pts[0], pts[1], t); -1: none
    double minTurn = 1e300, maxTurn = 0;  // measured [deg]
    double deviation = 0;     // smooth curve to polyline [mm]
    bool truncated = false;
};

// A refusal goes to the chunk's report (FibreReport::addCandidate) and to the trace.
void record(Candidate& c, FibreOutcome o) {
    c.out->rep.addCandidate(c.population, o, false);
    if (c.trace) *c.trace = static_cast<uint8_t>(static_cast<int>(o) | (c.random() ? 0x80 : 0));
}

// An accepted fibre: counted as truncated if a cut shortened it.
void recordAccepted(Candidate& c, const GrowingFibre& g) {
    c.out->rep.addCandidate(c.population, FibreOutcome::Accepted, g.truncated);
    if (c.trace) *c.trace = static_cast<uint8_t>(static_cast<int>(FibreOutcome::Accepted) | (c.random() ? 0x80 : 0));
}

bool insideOwnYarn(const FibreRun& run, const Candidate& c, const Vec3& q) { return run.index.exactInside(q, 0.0, c.ownYarn, {}).inside; }

// --- the anchor

// The anchor on the strand and the population, drawn before any rejection (some depend on it).
void placeAnchor(const FibreRun& run, const CandidateSource& source, Candidate& c) {
    const size_t k = static_cast<size_t>(c.strand);
    const int y = source.yarnId[k], n = source.count[k];
    const bool explicitCounts = run.frame && !run.frame->candidates.empty();
    const double jittered = c.i + hash01(run.seed, y, c.i, 1);
    c.s = explicitCounts ? jittered * c.ys->length / n : jittered / run.fp.density;
    c.theta = 2 * kPi * hash01(run.seed, y, c.i, 2);
    c.sp = evalSurface(run.cs, *c.ys, c.s, c.theta);
    if (run.frame && !run.frame->islands.empty()) c.island = islandAt(run.frame->islands[y], c.sp.ring.u);
    c.population = drawPopulation(run.fp, run.seed, y, c.i, c.island);
    c.ownYarn = SurfaceIndex::Exclude{y, c.s - run.ownWindow, c.s + run.ownWindow};
}

bool facesOutward(const FibreRun& run, Candidate& c) {
    c.N = run.curvedFabric ? normalized(run.frame->fabric(c.sp.p).n) : Vec3{0, 0, 1};
    return (run.curvedFabric ? dot(c.sp.n, c.N) : c.sp.n.z) > run.fp.frontCos;
}

// Plies: only the side facing out of the yarn can grow a directed fibre along it.
bool onOuterPlySide(const FibreRun& run, const Candidate& c) {
    if (c.random() || !run.onPlies()) return true;
    const auto [py, ps0] = plyParent(*run.frame->plies, c.yarn(), c.sp.ring.u);
    const YarnSurface& parent = run.frame->parentSurface->yarns[py];
    const Ring R = ringAt(*run.frame->parentCurves, parent, std::clamp(ps0, 0.0, parent.length));
    Vec3 out = c.sp.p - R.c;
    out = out - R.T * dot(out, R.T);
    return norm(out) >= 1e-12 && dot(c.sp.n, normalized(out)) >= run.fp.plyOuterCos;
}

// The ray along the fabric normal: its origin just off the surface lies outside every tube, and it
// clears the other yarns (conservative enclosures). Accepted: exposed.
FibreOutcome anchorExposure(const FibreRun& run, const Candidate& c) {
    const Vec3 origin = c.sp.p + c.sp.n * kOriginOffset;
    if (run.index.pointInEnclosure(origin, 0.0, c.ownYarn) || insideOwnYarn(run, c, origin)) return FibreOutcome::OriginInside;
    const double rayLength = run.curvedFabric ? run.frame->exposureLength : run.index.zTop() - origin.z + 1.0;
    return run.index.segmentHitsEnclosure(origin, origin + c.N * rayLength, 0.0, c.ownYarn) ? FibreOutcome::Occluded : FibreOutcome::Accepted;
}

// Artistic thinning, optionally biased towards high centreline curvature.
bool keptByThinning(const FibreRun& run, const Candidate& c) {
    const BSpline& curve = run.cs.yarns[c.yarn()].curve;
    const Vec3 d1 = curve.d1(c.sp.ring.u), d2 = curve.d2(c.sp.ring.u);
    const double kappa = norm(cross(d1, d2)) / std::pow(std::max(1e-12, norm(d1)), 3);
    const double kappaRef = 1.0 / (4 * run.r);
    const double factor = (1 - run.fp.curvatureWeight) + run.fp.curvatureWeight * std::min(2.0, kappa / kappaRef);
    return hash01(run.seed, c.yarn(), c.i, 3) < std::clamp(run.fp.keepProbability * factor * (1 - double(c.island.worn)), 0.0, 1.0);
}

// --- the shape

// The radius and the number of visible segments: about one per segLenTarget of the drawn length,
// each at least kMinSegmentRadii radii. Random: at least one interior vertex with a real turn.
void chooseSize(const FibreRun& run, const Candidate& c, GrowingFibre& g) {
    const FibreParams& fp = run.fp;
    const FibreDraw& d = g.draw;
    g.minSegments = c.random() ? std::max(2, fp.minSegments) : fp.minSegments;
    const int maxSegments = c.random() ? std::max(g.minSegments, fp.maxSegments) : fp.maxSegments;
    g.curled = c.random() && d.curlRadius > 0;
    // A curled stray or a crimped fibre has a real shape whatever the turn at one vertex (a crimp's
    // inflections turn by nothing): no lower limit. Nor has a hugging fibre: along a straight leg of
    // yarn it is straight.
    g.crimped = !c.random() && fp.crimpAmplitude > 0;
    const bool hugging = !c.random() && fp.hugHeight > 0;
    g.minTurnDeg = g.curled || g.crimped || hugging ? 0.0 : c.random() ? std::max(1.0, fp.minTurnDeg) : fp.minTurnDeg;
    g.radius = 0.5 * (c.random() && fp.randomThickness > 0 ? fp.randomThickness : fp.thickness) * d.radiusFactor;
    const double wanted = g.curled ? std::max(kCurlChordsPerTurn * d.curlTurns, d.length / kCurlMaxChord) : d.length / fp.segLenTarget;
    const int segments = std::clamp(static_cast<int>(std::lround(wanted)), g.minSegments, maxSegments);
    g.segments = std::min(segments, static_cast<int>(std::floor(d.length / (kMinSegmentRadii * g.radius))));
    g.buriedDepth = fp.rootDepth + g.radius;  // the export's inset + the radius, so the round end cap is hidden too
}

// On a ply: the directed curve lives in the parent yarn's coordinates, from the anchor (its arc
// length refined so the anchor lies in the ring plane), at the anchor's own radius lambda0. False:
// the anchor faces into the yarn.
bool toParentYarn(const FibreRun& run, const Candidate& c, DirectedCurve& curve) {
    const YarnSurfacePoint& sp = c.sp;
    const auto [py, ps0] = plyParent(*run.frame->plies, c.yarn(), sp.ring.u);
    const YarnSurface& parent = run.frame->parentSurface->yarns[py];
    double sP = std::clamp(ps0, 0.0, parent.length);
    Ring R = ringAt(*run.frame->parentCurves, parent, sP);
    for (int it = 0; it < 4; ++it) {
        sP = std::clamp(sP + dot(sp.p - R.c, R.T), 0.0, parent.length);
        R = ringAt(*run.frame->parentCurves, parent, sP);
    }
    const Vec3 v = sp.p - R.c;
    const Vec3 ea = R.n1 * std::cos(R.phi) + R.n2 * std::sin(R.phi), eb = R.n2 * std::cos(R.phi) - R.n1 * std::sin(R.phi);
    const double x = dot(v, ea) / R.a, y = dot(v, eb) / R.b;
    curve.cs = run.frame->parentCurves;
    curve.ys = &parent;
    curve.sE = sP;
    curve.thetaE = R.phi + std::atan2(y, x);
    curve.lambda0 = std::hypot(x, y);
    curve.r = std::sqrt(R.a * R.b);
    const Vec3 out = v - R.T * dot(v, R.T);
    return norm(out) >= 1e-12 && dot(sp.n, normalized(out)) >= run.fp.plyOuterCos;
}

// The smooth curve. The own tube of SurfaceIndex is the undisplaced envelope, and both populations
// start on it. False: the anchor faces into the yarn.
bool setUpCurve(const FibreRun& run, const Candidate& c, GrowingFibre& g) {
    const FibreDraw& d = g.draw;
    SmoothCurve& curve = g.curve;
    curve.random = c.random();
    curve.bowRadial = d.bowRadial;
    if (c.random()) {
        curve.stray = randomCurve(d, c.sp.p, c.sp.n, c.sp.ring.T);
        curve.maxSpan = 10 * d.length + 1;  // no yarn-end bound: the root's arc is checked in placeRoot
        return true;
    }
    DirectedCurve& dc = curve.directed;
    dc.cs = &run.cs;
    dc.ys = c.ys;
    dc.sE = c.s;
    dc.thetaE = c.theta;
    dc.sgn = d.sgn;
    dc.p = std::tan(d.psi);
    dc.q = std::tan(d.lift);
    dc.bowTheta = d.bowTheta;
    dc.r = run.r;
    dc.crimpAmplitude = run.fp.crimpAmplitude;
    dc.crimpWavelength = run.fp.crimpWavelength;
    dc.crimpPhase = d.crimpPhase;
    dc.crimpPlane = d.crimpPlane;
    dc.hug = run.fp.hugHeight;
    if (run.onPlies() && !toParentYarn(run, c, dc)) return false;
    // Along the ply helix: d theta / ds = handedness 2 pi / pitch (plies.cpp), so tan(psi) = 2 pi r / pitch.
    if (run.fp.plyPitch > 0) {
        dc.p = std::tan(run.fp.handedness * std::atan(2 * kPi * dc.r / run.fp.plyPitch) + d.psi);
        dc.rampWinding = true;
    }
    curve.maxSpan = d.sgn > 0 ? dc.ys->length - kYarnEndMargin - dc.sE : dc.sE - kYarnEndMargin;
    return true;
}

// The root: a buried segment of its own, from a point behind and below the curve's start, in the
// plane of the inward normal and the first visible chord, tilted so that the turn at the start is
// rootTurnDeg (the tilt from the normal kept in [10, 70] deg), deep enough for buriedDepth. Its
// own-yarn coordinate must stay inside the strand: ringAt clamps outside.
void placeRoot(const FibreRun& run, const Candidate& c, GrowingFibre& g, const std::vector<Vec3>& visible) {
    const YarnSurfacePoint& sp = c.sp;
    const Vec3 start = g.curve.table[0] + sp.n * 1e-4;
    Vec3 chord = visible[0] - start;
    g.firstChordElevation = std::asin(std::clamp(dot(normalized(chord), sp.n), -1.0, 1.0)) / kDeg;
    chord = chord - sp.n * dot(chord, sp.n);
    chord = norm(chord) > 1e-12 ? normalized(chord) : normalized(sp.ring.T * double(g.draw.sgn));
    const double tilt = std::clamp(90.0 - g.firstChordElevation - run.fp.rootTurnDeg, 10.0, 70.0) * kDeg;
    const double rootLength = std::max(g.buriedDepth / std::cos(tilt), 4.5 * g.radius);
    g.pts.assign(1, start - (sp.n * std::cos(tilt) + chord * std::sin(tilt)) * rootLength);
    g.pts.push_back(start);
    g.pts.insert(g.pts.end(), visible.begin(), visible.end());
    const double rootArc = c.s - dot(chord, sp.ring.T) * rootLength * std::sin(tilt);
    g.rootOnYarn = rootArc >= kYarnEndMargin && rootArc <= c.ys->length - kYarnEndMargin;
}

// The surface crossing on segment 1 (root inside, curve start outside the own tube), by bisection.
void findEmergence(const FibreRun& run, const Candidate& c, GrowingFibre& g) {
    g.emergence = -1;
    const Vec3 &a = g.pts[0], &b = g.pts[1];
    if (!insideOwnYarn(run, c, a) || insideOwnYarn(run, c, b)) return;
    double lo = 0, hi = 1;
    for (int it = 0; it < 30; ++it) {
        const double mid = 0.5 * (lo + hi);
        if (insideOwnYarn(run, c, lerp(a, b, mid))) lo = mid;
        else hi = mid;
    }
    g.emergence = hi;
}

// The polyline from the emergence to the tip.
double visiblePolylineLength(const GrowingFibre& g) {
    double L = (1 - g.emergence) * norm(g.pts[1] - g.pts[0]);
    for (size_t k = 2; k < g.pts.size(); ++k) L += norm(g.pts[k] - g.pts[k - 1]);
    return L;
}

// The polyline whose visible length matches the drawn one: solve the span for a smooth target,
// sample, place the root, find the emergence, then correct the target once for the chord and
// emergence offsets. False: no span reaches the target; otherwise rootOnYarn and emergence tell how
// far the fit got.
bool fitPolyline(const FibreRun& run, const Candidate& c, GrowingFibre& g) {
    double target = g.draw.length;
    for (int it = 0; it < 2; ++it) {
        if (!g.curve.solveLength(target)) return false;
        placeRoot(run, c, g, g.curve.equalArcVertices(g.segments));
        if (!g.rootOnYarn) return true;
        findEmergence(run, c, g);
        if (g.emergence < 0) return true;
        if (it == 0) target += g.draw.length - visiblePolylineLength(g);
    }
    return true;
}

// --- the checks on the polyline

// The renderer's limits (strands_file.h): every segment of pts (g's polyline or its draped copy) at least
// kMinSegmentRadii radii, at most 200 radii at the thinnest radius the export applies (Fibre::vertexRadius).
bool meetsStrandRules(const FibreRun& run, const GrowingFibre& g, const std::vector<Vec3>& pts) {
    const double thinnest = std::max(run.fp.taper, std::pow(kTaperStep, -double(g.segments + 1)));
    for (size_t k = 1; k < pts.size(); ++k) {
        const double length = norm(pts[k] - pts[k - 1]);
        if (length < kMinSegmentRadii * g.radius || length > 200 * g.radius * thinnest) return false;
    }
    return true;
}

// The smallest and largest turn at the interior vertices of pts [deg], lowering minTurn and raising maxTurn.
void measureTurns(const std::vector<Vec3>& pts, double& minTurn, double& maxTurn) {
    for (size_t k = 1; k + 1 < pts.size(); ++k) {
        const Vec3 a = normalized(pts[k] - pts[k - 1]), b = normalized(pts[k + 1] - pts[k]);
        const double turn = std::acos(std::clamp(dot(a, b), -1.0, 1.0)) / kDeg;
        minTurn = std::min(minTurn, turn);
        maxTurn = std::max(maxTurn, turn);
    }
}

// Light contact: every segment midpoint and vertex after the curve's start may lie at most 20 % of
// the yarn radius inside any strand, its own included (the root segment starts inside it by
// construction). The push that lifts q out of the strand it sinks into deepest: along that strand's
// outward normal by the depth plus the fibre's radius; zero when q is not too deep.
Vec3 deepContactPush(const FibreRun& run, const Vec3& q, double fibreRadius) {
    const SurfaceIndex::Hit hit = run.index.exactInside(q, 0.0, {}, {});
    if (!hit.inside || hit.depth <= kDeepContact * run.r) return {0, 0, 0};
    return hit.normal * (hit.depth + fibreRadius);
}
bool isZero(const Vec3& v) { return v.x == 0 && v.y == 0 && v.z == 0; }

// The first segment of pts (from the curve's start) whose midpoint or end lies too deep; pts.size(): none.
size_t firstDeepSegment(const FibreRun& run, const std::vector<Vec3>& pts, double fibreRadius) {
    for (size_t k = 2; k < pts.size(); ++k)
        if (!isZero(deepContactPush(run, lerp(pts[k - 1], pts[k], 0.5), fibreRadius)) || !isZero(deepContactPush(run, pts[k], fibreRadius)))
            return k;
    return pts.size();
}

// The length of pts from the curve's start (the root segment, which the drape never moves, left out).
double lengthFromCurveStart(const std::vector<Vec3>& pts) {
    double L = 0;
    for (size_t k = 2; k < pts.size(); ++k) L += norm(pts[k] - pts[k - 1]);
    return L;
}

// A draped copy of a fibre's polyline and what its checks measure on it.
struct DrapedPolyline {
    std::vector<Vec3> pts;
    double minTurn = 1e300, maxTurn = 0;
    double deviation = 0;
};
constexpr int kDrapePushes = 3;  // per vertex; a sample still too deep after them ends the drape

// Drape (FibreParams::drape): a fibre that sinks into a strand rests on it. Walking the segments as the deep-contact
// cut does, a segment whose midpoint (or else its end) lies too deep lifts its end vertex k by that sample's push,
// and the next two vertices by 2/3 and 1/3 of it: the fibre lifts over the stitch and settles. The draped polyline
// keeps its point count and is never longer than the original (the tip segment is shortened by any excess), so the
// drawn length still bounds it. It must pass the strand rules, the turn limits and the deep-contact walk; then its
// turns and its deviation from the smooth curve are measured again (the exported geometry: the approximation check
// passed on the smooth fit, and the departure is the drape, at most r + radius per push). False: the drape failed,
// and the fibre is cut as without it.
bool drapeOverStrands(const FibreRun& run, const GrowingFibre& g, DrapedPolyline& d) {
    d.pts = g.pts;
    std::vector<Vec3>& P = d.pts;
    for (size_t k = 2; k < P.size(); ++k) {
        for (int push = 0;; ++push) {
            Vec3 lift = deepContactPush(run, lerp(P[k - 1], P[k], 0.5), g.radius);
            if (isZero(lift)) lift = deepContactPush(run, P[k], g.radius);
            if (isZero(lift)) break;
            if (push == kDrapePushes) return false;
            P[k] = P[k] + lift;
            if (k + 1 < P.size()) P[k + 1] = P[k + 1] + lift * (2.0 / 3.0);
            if (k + 2 < P.size()) P[k + 2] = P[k + 2] + lift * (1.0 / 3.0);
        }
    }
    const double excess = lengthFromCurveStart(P) - lengthFromCurveStart(g.pts);
    if (excess > 0) {
        const Vec3 tip = P.back() - P[P.size() - 2];
        const double tipLength = norm(tip) - excess;
        if (tipLength < kMinSegmentRadii * g.radius) return false;
        P.back() = P[P.size() - 2] + tip * (tipLength / norm(tip));
    }
    if (!meetsStrandRules(run, g, P)) return false;
    measureTurns(P, d.minTurn, d.maxTurn);
    if (d.maxTurn > run.fp.maxTurnDeg || d.minTurn < g.minTurnDeg) return false;
    if (firstDeepSegment(run, P, g.radius) < P.size()) return false;
    d.deviation = g.curve.deviationFrom(P);
    return true;
}

// The deep-contact rule: a directed fibre with FibreParams::drape first tries to rest on what it sinks into
// (drapeOverStrands); otherwise it is cut before the first segment that fails, and what remains must keep its
// minimum of visible segments. Accepted: the fibre stands, draped, cut or neither.
FibreOutcome cutAtDeepContact(const FibreRun& run, const Candidate& c, GrowingFibre& g) {
    const size_t keep = firstDeepSegment(run, g.pts, g.radius);
    if (keep == g.pts.size()) return FibreOutcome::Accepted;
    DrapedPolyline draped;
    if (run.fp.drape && !c.random() && drapeOverStrands(run, g, draped)) {
        g.pts = draped.pts;
        g.minTurn = draped.minTurn;
        g.maxTurn = draped.maxTurn;
        g.deviation = draped.deviation;
        return FibreOutcome::Accepted;
    }
    if (static_cast<int>(keep) - 2 < std::max(1, g.minSegments)) return keep >= 3 ? FibreOutcome::CollisionShort : FibreOutcome::Collision;
    g.pts.resize(keep);
    g.truncated = true;
    return FibreOutcome::Accepted;
}

// No vertex from the curve start on, and no segment midpoint, may lie behind the fabric's back. The
// fibre is cut before the first segment that fails; what remains must keep its minimum of visible
// segments. Accepted: the fibre stands, cut or not.
FibreOutcome cutBehindFabric(const FibreRun& run, GrowingFibre& g) {
    const auto behind = [&](const Vec3& q) { return run.frontOfBack(q) < 0; };
    size_t keep = behind(g.pts[1]) ? 2 : g.pts.size();
    for (size_t k = 2; k < g.pts.size() && keep == g.pts.size(); ++k)
        if (behind(lerp(g.pts[k - 1], g.pts[k], 0.5)) || behind(g.pts[k])) keep = k;
    if (keep == g.pts.size()) return FibreOutcome::Accepted;
    if (static_cast<int>(keep) - 2 < std::max(1, g.minSegments)) return FibreOutcome::BehindFabric;
    g.pts.resize(keep);
    g.truncated = true;
    return FibreOutcome::Accepted;
}

// The accepted fibre, with its draws and what the checks measured.
Fibre toFibre(const FibreRun& run, const Candidate& c, const GrowingFibre& g) {
    const FibreDraw& d = g.draw;
    const std::vector<Vec3>& P = g.pts;
    Fibre f;
    f.strand = c.strand;
    f.yarnId = c.yarn();
    f.parentYarn = run.onPlies() ? plyParent(*run.frame->plies, c.yarn(), c.sp.ring.u).first : c.yarn();
    f.candidate = c.i;
    f.population = c.population;
    f.s = c.s;
    f.theta = float(c.theta);
    f.anchor = toVec3f(lerp(P[0], P[1], g.emergence));
    f.normal = toVec3f(c.sp.n);
    const double emergeL = g.emergence * norm(P[1] - P[0]);
    f.emergeL = float(emergeL);
    f.radius0 = float(g.radius);
    f.radius1 = float(g.radius * run.fp.taper);
    f.targetLength = float(d.length);
    f.attach = float(emergeL + g.radius);
    f.attachArc = float(emergeL + run.r + g.radius);
    f.key = d.key;
    f.truncated = g.truncated;
    f.psi = float((c.random() ? d.azimuth : d.psi) / kDeg);
    f.lift = float((c.random() ? d.elev : d.lift) / kDeg);
    f.smoothLength = float(g.curve.length());
    f.emergeAngle = float(std::asin(std::clamp(dot(normalized(P[2] - P[1]), c.sp.n), -1.0, 1.0)) / kDeg);
    if (c.random()) {
        f.helixActual = float(std::acos(std::clamp(std::fabs(dot(g.curve.stray.d0, c.sp.ring.T)), 0.0, 1.0)) / kDeg);
    } else {
        const DirectedCurve& dc = g.curve.directed;
        const Vec3 T = normalized(yarnCoordinates(*dc.cs, *dc.ys, dc.sE, dc.thetaE, 0.0).dS);
        f.helixActual = float(std::acos(std::clamp(std::fabs(dot(normalized(dc.dX(0.0)), T)), 0.0, 1.0)) / kDeg);
    }
    f.maxDev = float(g.deviation);
    f.minTurn = float(g.minTurn);
    f.maxTurn = float(g.maxTurn);
    return f;
}

// --- one candidate, in the order it happens

void growCandidate(const FibreRun& run, Candidate& c) {
    if (!facesOutward(run, c)) return record(c, FibreOutcome::BackFacing);
    if (!onOuterPlySide(run, c)) return record(c, FibreOutcome::PlyInnerEarly);
    const FibreOutcome exposure = anchorExposure(run, c);
    if (exposure != FibreOutcome::Accepted) return record(c, exposure);
    if (!keptByThinning(run, c)) return record(c, FibreOutcome::Probability);

    GrowingFibre g;
    g.draw = drawFibre(run.fp, run.seed, c.yarn(), c.i, c.island);
    chooseSize(run, c, g);
    if (g.segments < std::max(1, g.minSegments)) return record(c, FibreOutcome::TooShortToBend);
    if (!setUpCurve(run, c, g)) return record(c, FibreOutcome::PlyInnerLate);
    if (g.buriedDepth >= 0.9 * run.r) return record(c, FibreOutcome::RootDepth);
    if (!fitPolyline(run, c, g)) return record(c, FibreOutcome::YarnEnd);
    if (!g.rootOnYarn) return record(c, FibreOutcome::RootRange);
    if (g.firstChordElevation < run.fp.minEmergeDeg) return record(c, FibreOutcome::NoEmergenceElev);
    if (g.emergence < 0) return record(c, insideOwnYarn(run, c, g.pts[0]) ? FibreOutcome::NoEmergence : FibreOutcome::RootNotInside);
    if (!meetsStrandRules(run, g, g.pts)) return record(c, FibreOutcome::RuleRejected);
    measureTurns(g.pts, g.minTurn, g.maxTurn);
    if (g.maxTurn > run.fp.maxTurnDeg) return record(c, FibreOutcome::TurnTooSharp);
    if (g.minTurn < g.minTurnDeg) return record(c, FibreOutcome::NoBend);
    g.deviation = g.curve.deviationFrom(g.pts);
    // A curled stray's limit follows its curl (kCurlChordsPerTurn: a sixteenth of a turn), a crimped
    // fibre's its crimp.
    const double allowedDev = g.curled    ? std::max(run.fp.maxCurveDev, 0.06 * g.draw.curlRadius)
                              : g.crimped ? std::max(run.fp.maxCurveDev, 0.06 * run.fp.crimpAmplitude)
                                          : run.fp.maxCurveDev;
    if (g.deviation > allowedDev)
        return record(c, FibreOutcome::ApproxError);
    const FibreOutcome contact = cutAtDeepContact(run, c, g);
    if (contact != FibreOutcome::Accepted) return record(c, contact);
    const FibreOutcome back = cutBehindFabric(run, g);
    if (back != FibreOutcome::Accepted) return record(c, back);
    recordAccepted(c, g);
    c.out->addFibre(toFibre(run, c, g), g.pts);
}

} // namespace

int64_t CandidateSource::total() const {
    int64_t t = 0;
    for (int n : count) t += std::max(0, n);
    return t;
}

float CandidateSource::key(int strand, int i) const { return float(hash01(seed, yarnId[static_cast<size_t>(strand)], i, 9)); }

uint64_t fibreSeed(const FibreParams& fp, uint64_t pieceSeed) { return pieceSeed ? splitmix64(fp.seed ^ splitmix64(pieceSeed)) : fp.seed; }

CandidateSource candidateSource(const SurfaceModel& sm, const FibreParams& fp, const FibreFrame* frame) {
    CandidateSource src;
    src.seed = fibreSeed(fp, frame ? frame->pieceSeed : 0);
    const bool explicitCounts = frame && !frame->candidates.empty();
    for (const YarnSurface& ys : sm.yarns) {
        src.yarnId.push_back(ys.yarnId);
        src.count.push_back(explicitCounts ? frame->candidates[ys.yarnId] : static_cast<int>(std::floor(fp.density * ys.length)));
    }
    return src;
}

// The candidates run in chunks, each into its own report and fibres (in parallel with a parallel
// FibreFrame::exec), merged in candidate order, so the result is the serial one bit for bit.
FibreSet growYarnFibres(const CurveSet& cs, const SurfaceModel& sm, const PatternParams& pp, const FibreParams& fp, const FibreFrame* frame) {
    const auto t0 = std::chrono::steady_clock::now();
    const FibreRun run(cs, sm, pp, fp, frame);
    const CandidateSource source = candidateSource(sm, fp, frame);
    const ExecPolicy policy = frame ? frame->exec : ExecPolicy{};
    const int64_t chunkSize = policy.chunk > 0 ? policy.chunk : 512;
    std::vector<ChunkRange> chunks;
    std::vector<int64_t> firstCandidate;
    int64_t candidates = 0;
    for (int s = 0; s < source.strands(); ++s) {
        const int n = std::max(0, source.count[static_cast<size_t>(s)]);
        firstCandidate.push_back(candidates);
        candidates += n;
        for (int64_t b = 0; b < n; b += chunkSize) chunks.push_back({s, b, std::min<int64_t>(b + chunkSize, n)});
    }
    std::vector<uint8_t>* const trace = frame ? frame->outcomes : nullptr;
    if (trace) trace->assign(static_cast<size_t>(candidates), static_cast<uint8_t>(FibreOutcome::NotProcessed));
    std::vector<FibreSet> results(chunks.size());
    parallelFor(static_cast<int64_t>(chunks.size()), 1, [&](int64_t begin, int64_t end, int64_t) {
        for (int64_t k = begin; k < end; ++k) {
            cancelCheckpoint();
            const ChunkRange& chunk = chunks[static_cast<size_t>(k)];
            FibreSet& result = results[static_cast<size_t>(k)];
            result.rep.pop.resize(2);
            for (int64_t t = chunk.begin; t < chunk.end; ++t) {
                Candidate c;
                c.ys = &sm.yarns[static_cast<size_t>(chunk.strand)];
                c.strand = chunk.strand;
                c.i = static_cast<int>(t);
                c.out = &result;
                if (trace) c.trace = &(*trace)[static_cast<size_t>(firstCandidate[static_cast<size_t>(chunk.strand)] + t)];
                placeAnchor(run, source, c);
                growCandidate(run, c);
            }
        }
    }, policy);
    FibreSet out;
    out.rep.pop.resize(2);
    for (const FibreSet& result : results) {
        out.rep.addCounters(result.rep);
        out.appendSet(result);
    }
    resumFibreLengths(out);
    out.rep.ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    return out;
}

} // namespace knit
