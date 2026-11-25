// Focused tests for the geometry library. Run one case: knit_tests <name>;
// all cases: knit_tests; list: knit_tests --list.
#include <algorithm>
#include <chrono>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <regex>
#include <stdexcept>
#include <set>
#include <string>
#include <vector>
#include <unordered_map>
#include <thread>

#include "cpu_numerics.h"
#include "knitting/layout/chart.h"
#include "knitting/knit/pieces.h"
#include "io/files.h"
#include "knitting/checks/checks.h"
#include "knitting/knit/fixtures.h"
#include "knitting/stitch/hand_tension.h"
#include "fibres/fibre_checks.h"
#include "model/choices.h"
#include "model/compute_caps.h"
#include "model/coordinates.h"
#include "model/hash.h"
#include "fibres/fibre.h"
#include "ncp_reference.h"
#include "ply_reference.h"
#include "fibres/fibre_generator.h"
#include "fibres/fibre_islands.h"
#include "fibres/fibre_shape.h"
#include "knitting/layout/stitch_layout.h"
#include "model/mesh_query.h"
#include "model/parallel.h"
#include "knitting/knit/piece_knit.h"
#include "knitting/tube/plies.h"
#include "knitting/tube/ply_profile.h"
#include "knitting/relax/settle.h"
#include "knitting/tube/yarn_surface.h"
#include "knitting/layout/uv_atlas.h"
#include "knitting/stitch/stitch_program.h"
#include "model/spatial_hash.h"
#include "knitting/stitch/stockinette.h"
#include "fibres/fibre_strands.h"
#include "knitting/knit/coverage.h"
#include "strand_fixtures.h"

#include "fibres/fibre_budget.h"

using namespace knit;

namespace {

struct Case { const char* name; void (*fn)(); };
std::vector<Case>& registry() { static std::vector<Case> r; return r; }
struct Reg { Reg(const char* n, void (*f)()) { registry().push_back({n, f}); } };
// Thrown by SKIP(): reported as skipped, never as passed.
struct Skipped { const char* why; };

#define TEST(name) static void name(); static Reg reg_##name(#name, name); static void name()
#define CHECK(cond) do { if (!(cond)) { char b_[512]; std::snprintf(b_, sizeof b_, "%s:%d: CHECK(%s)", __FILE__, __LINE__, #cond); throw std::runtime_error(b_); } } while (0)
#define CHECK_NEAR(a, b, tol) do { double a_ = (a), b_ = (b); if (!(std::fabs(a_ - b_) <= (tol))) { char m_[512]; std::snprintf(m_, sizeof m_, "%s:%d: |%s - %s| = |%.17g - %.17g| > %g", __FILE__, __LINE__, #a, #b, a_, b_, double(tol)); throw std::runtime_error(m_); } } while (0)
#define SKIP(why) throw Skipped{why}

// ---------------------------------------------------------------- helpers ----

// Uniform B-spline controls whose single span equals the cubic Bezier b.
std::vector<Vec3> bezierToBSpline(const Vec3 b[4]) {
    Vec3 P1 = b[1] * 2.0 - b[2], P2 = b[2] * 2.0 - b[1];
    Vec3 P0 = b[0] * 6.0 - P1 * 4.0 - P2, P3 = b[3] * 6.0 - P1 - P2 * 4.0;
    return {P0, P1, P2, P3};
}

PatternParams testParams(double D = 1.0, double c = 0.1, double eps = 0.005) {
    PatternParams p;
    p.yarnDiameter = D;
    p.clearance = c;
    p.sampleTolerance = eps;
    return p;
}

ClearanceResult clearanceOf(const CurveSet& cs, const PatternParams& pp) {
    Allowances al = Allowances::forExtent(cs.extent());
    std::vector<CurvatureResult> cv;
    std::vector<Polyline> pl;
    for (const YarnPath& y : cs.yarns) {
        cv.push_back(checkCurvature(y, pp.radius(), al));
        pl.push_back(samplePath(y, pp.sampleTolerance));
    }
    return checkClearance({cs, pl, cv}, pp, al);
}

CrossingResult crossingsOf(const CurveSet& cs) { return findCrossings(cs, Allowances::forExtent(cs.extent())); }

std::vector<Vec3> wavyCtrl() {
    return {{0, 0, 0}, {1, 2, 0.5}, {2, -1, 1}, {3.5, 0.5, -0.5}, {4, 3, 0}, {5, 1, 2}, {6, -2, 1}, {7, 0, 0}};
}

// The symmetric loop Bezier: self-intersects in projection at t1 + t2 = 1.
const Vec3 kLoop[4] = {{0, 0, 0}, {1.5, 1, 1.0 / 3}, {-0.5, 1, 2.0 / 3}, {1, 0, 1}};

// ---------------------------------------------------------------- B-spline ----

// Independent evaluation through the cubic B-spline basis functions.
Vec3 evalBasisReference(const BSpline& s, double u) {
    const SpanPosition at = s.locate(u);
    const int k = at.span;
    const double t = at.t;
    double t2 = t * t, t3 = t2 * t;
    double N0 = (1 - t) * (1 - t) * (1 - t) / 6.0;
    double N1 = (3 * t3 - 6 * t2 + 4) / 6.0;
    double N2 = (-3 * t3 + 3 * t2 + 3 * t + 1) / 6.0;
    double N3 = t3 / 6.0;
    return s.ctrl[k] * N0 + s.ctrl[k + 1] * N1 + s.ctrl[k + 2] * N2 + s.ctrl[k + 3] * N3;
}

TEST(bspline_reference) {
    BSpline s; s.ctrl = wavyCtrl();
    for (int i = 0; i <= 500; ++i) {
        double u = s.uMax() * i / 500.0;
        CHECK(norm(s.eval(u) - evalBasisReference(s, u)) < 1e-13);
    }
}

TEST(bspline_derivatives) {
    BSpline s; s.ctrl = wavyCtrl();
    const double h = 1e-5;
    for (int i = 1; i < 100; ++i) {
        double u = s.uMax() * i / 100.0;
        if (std::fabs(u - std::round(u)) < 2 * h) continue;
        Vec3 fd1 = (s.eval(u + h) - s.eval(u - h)) / (2 * h);
        Vec3 fd2 = (s.eval(u + h) - s.eval(u) * 2.0 + s.eval(u - h)) / (h * h);
        CHECK(norm(fd1 - s.d1(u)) < 1e-7);
        CHECK(norm(fd2 - s.d2(u)) < 1e-3);
    }
}

TEST(bspline_continuity) {
    // Regression test of the evaluator: C0/C1/C2 vector jumps at every knot.
    BSpline s; s.ctrl = wavyCtrl();
    for (int k = 0; k + 1 < s.numSpans(); ++k) {
        CubicBezier L = s.span(k), R = s.span(k + 1);
        CHECK(norm(L.eval(1) - R.eval(0)) < 1e-13);
        CHECK(norm(L.d1(1) - R.d1(0)) < 1e-12);
        CHECK(norm(L.d2(1) - R.d2(0)) < 1e-11);
    }
}

TEST(bezier_split) {
    CubicBezier b{{{0, 0, 0}, {1, 2, 3}, {4, -1, 2}, {5, 5, 0}}};
    const BezierHalves halves = b.split(0.3);
    const CubicBezier& L = halves.left;
    const CubicBezier& R = halves.right;
    CubicBezier S = b.sub(0.2, 0.7);
    QuadBezier H = b.hodograph().sub(0.2, 0.7);
    for (int i = 0; i <= 20; ++i) {
        double s = i / 20.0;
        CHECK(norm(L.eval(s) - b.eval(0.3 * s)) < 1e-13);
        CHECK(norm(R.eval(s) - b.eval(0.3 + 0.7 * s)) < 1e-13);
        CHECK(norm(S.eval(s) - b.eval(0.2 + 0.5 * s)) < 1e-13);
        CHECK(norm(H.eval(s) - b.d1(0.2 + 0.5 * s)) < 1e-12);
    }
}

// ------------------------------------------------------------- arc length ----

TEST(arclength_line) {
    YarnPath y = makeYarn(lineCtrl({0, 0, 0}, {3, 4, 0}, 7), 0);
    CHECK_NEAR(y.length(), 5.0, 1e-12);
    CHECK_NEAR(y.arcLengthAt(3.5), 2.5, 1e-12);
    CHECK(norm(y.start() - Vec3{0, 0, 0}) < 1e-13 && norm(y.end() - Vec3{3, 4, 0}) < 1e-13);
}

TEST(arclength_helix) {
    // A B-spline cannot represent a helix exactly, so the reference is a dense
    // chord sum of the same helix-like spline (2e6 chords; error ~ L*h^2*kappa^2/24).
    std::vector<Vec3> c;
    for (int i = -1; i <= 25; ++i) c.push_back({2 * std::cos(0.5 * i), 2 * std::sin(0.5 * i), 0.3 * i});
    YarnPath y = makeYarn(c, 0);
    const int N = 2'000'000;
    double ref = 0;
    Vec3 prev = y.curve.eval(0);
    for (int i = 1; i <= N; ++i) {
        Vec3 p = y.curve.eval(y.uEnd() * i / N);
        ref += norm(p - prev);
        prev = p;
    }
    CHECK_NEAR(y.length(), ref, 1e-8);
    CHECK(y.lengthErrorEstimate() < 1e-8);
}

TEST(arclength_roundtrip) {
    YarnPath y = makeYarn(wavyCtrl(), 0);
    for (int i = 0; i <= 200; ++i) {
        double u = y.uEnd() * i / 200.0;
        double s = y.arcLengthAt(u);
        CHECK_NEAR(y.paramAtArcLength(s), u, 1e-9);
    }
}

TEST(arclength_bounds) {
    YarnPath y = makeYarn(wavyCtrl(), 0);
    for (int i = 0; i < 50; ++i) {
        double u0 = std::fmod(0.37 * i, y.uEnd()), u1 = std::min(y.uEnd(), u0 + 0.05 + 0.11 * (i % 7));
        const ArcLengthBounds b = y.lengthBounds(u0, u1);
        double est = y.arcLengthAt(u1) - y.arcLengthAt(u0);
        CHECK(b.lo <= est + 1e-12 && est <= b.hi + 1e-12);
    }
}

TEST(sampling_bound) {
    YarnPath y = makeYarn(wavyCtrl(), 0);
    for (double eps : {0.1, 0.01, 0.001}) {
        Polyline pl = samplePath(y, eps);
        CHECK(pl.epsBound <= eps);
        double worst = 0;
        for (size_t k = 0; k < pl.segments(); ++k)
            for (int j = 1; j < 16; ++j) {
                double a = j / 16.0;
                Vec3 c = y.curve.eval(pl.u[k] + a * (pl.u[k + 1] - pl.u[k]));
                worst = std::max(worst, norm(c - lerp(pl.p[k], pl.p[k + 1], a)));
            }
        CHECK(worst <= pl.epsBound * (1 + 1e-9));
    }
}

// ---------------------------------------------------------------- extrema ----

TEST(extrema_speed_interior) {
    // Nearly cusped span: the speed minimum is interior, far below the end speeds.
    const Vec3 b[4] = {{0, 0, 0}, {1, 1, 0}, {0.1, 1, 0}, {1, 0, 0}};
    YarnPath y = makeYarn(bezierToBSpline(b), 0);
    double dense = 1e300;
    for (int i = 0; i <= 100000; ++i) dense = std::min(dense, norm(y.curve.d1(i / 100000.0)));
    double ends = std::min(norm(y.curve.d1(0)), norm(y.curve.d1(1)));
    CHECK(ends > 10 * dense);  // knot/vertex sampling would miss it
    CurvatureResult r = checkCurvature(y, 0.01, Allowances::forExtent(2));
    CHECK(r.speedLo <= dense);
    CHECK(r.speedLo > 0);
    CHECK(r.regular == Status::Pass);
}

TEST(extrema_curvature_interior) {
    const Vec3 b[4] = {{0, 0, 0}, {2, 1, 0}, {-1, 1, 0}, {1, 0, 0}};
    YarnPath y = makeYarn(bezierToBSpline(b), 0);
    double dense = 0;
    for (int i = 0; i <= 100000; ++i) {
        double t = i / 100000.0;
        Vec3 d1 = y.curve.d1(t), d2 = y.curve.d2(t);
        dense = std::max(dense, norm(cross(d1, d2)) / std::pow(norm(d1), 3));
    }
    auto kap = [&](double t) { Vec3 d1 = y.curve.d1(t), d2 = y.curve.d2(t); return norm(cross(d1, d2)) / std::pow(norm(d1), 3); };
    CHECK(std::max(kap(0), kap(1)) < 0.5 * dense);
    CurvatureResult r = checkCurvature(y, 1.0 / (2 * dense), Allowances::forExtent(2));
    CHECK(r.kappaHi >= dense);
    CHECK(r.kappaHi <= 1.05 * dense);
    // With r = 1/(2 kappa_max) the fold check passes; with r = 2/kappa_max it fails.
    CHECK(r.noFold == Status::Pass);
    CHECK(checkCurvature(y, 2.0 / dense, Allowances::forExtent(2)).noFold == Status::Fail);
}

// -------------------------------------------------------------- crossings ----

TEST(cross_two_close) {
    // Line y = 0 (z = 0) and parabola y = x^2 - a^2 (z = 1): crossings at x = +-a.
    const double a = 1e-3;
    const Vec3 bp[4] = {{-1, 1 - a * a, 1}, {-1.0 / 3, -1.0 / 3 - a * a, 1}, {1.0 / 3, -1.0 / 3 - a * a, 1}, {1, 1 - a * a, 1}};
    CurveSet cs = makeSet({makeYarn(lineCtrl({-1, 0, 0}, {1, 0, 0}, 1), 0), makeYarn(bezierToBSpline(bp), 1)});
    CrossingResult r = crossingsOf(cs);
    CHECK(r.unresolved.empty());
    CHECK(r.crossings.size() == 2);
    CHECK(!r.crossings[0].boxA.overlaps(r.crossings[1].boxA));
    for (const Crossing& x : r.crossings) {
        CHECK(x.front == 1);
        CHECK_NEAR(std::fabs(x.pA.x), a, 1e-9);
    }
}

TEST(cross_tangent) {
    const Vec3 bp[4] = {{-1, 1, 1}, {-1.0 / 3, -1.0 / 3, 1}, {1.0 / 3, -1.0 / 3, 1}, {1, 1, 1}};  // y = x^2
    CurveSet cs = makeSet({makeYarn(lineCtrl({-1, 0, 0}, {1, 0, 0}, 1), 0), makeYarn(bezierToBSpline(bp), 1)});
    CrossingResult r = crossingsOf(cs);
    CHECK(r.crossings.empty());
    CHECK(!r.unresolved.empty());
}

TEST(cross_at_boundary) {
    // Vertical line at a recursion subdivision boundary (u = 2.5 on a 4-span line).
    CurveSet cs = makeSet({makeYarn(lineCtrl({0, 0, 0}, {4, 0, 0}, 4), 0),
                           makeYarn(lineCtrl({2.5, -1, 1}, {2.5, 1, 1}, 2), 1)});
    CrossingResult r = crossingsOf(cs);
    CHECK(r.unresolved.empty());
    CHECK(r.crossings.size() == 1);
    CHECK_NEAR(r.crossings[0].uA, 2.5, 1e-9);
}

TEST(cross_collinear_overlap) {
    CurveSet cs = makeSet({makeYarn(lineCtrl({0, 0, 0}, {2, 0, 0}, 2), 0),
                           makeYarn(lineCtrl({1, 0, 1}, {3, 0, 1}, 2), 1)});
    CrossingResult r = crossingsOf(cs);
    CHECK(r.crossings.empty());
    CHECK(!r.unresolved.empty());
}

TEST(cross_single_span_loop) {
    CurveSet cs = makeSet({makeYarn(bezierToBSpline(kLoop), 0)});
    CrossingResult r = crossingsOf(cs);
    CHECK(r.unresolved.empty());
    CHECK(r.crossings.size() == 1);
    const Crossing& x = r.crossings[0];
    CHECK(x.yarnA == 0 && x.yarnB == 0);
    CHECK_NEAR(x.uA + x.uB, 1.0, 1e-9);
    CHECK(x.front == 0 && x.depth.hi < 0);  // later point is in front (z grows with t)
}

TEST(cross_adjacent_span) {
    // The same loop split at t = 0.5 into two C2-joined spans of a uniform B-spline:
    // the crossing is now between adjacent spans, away from their join.
    CubicBezier b{{kLoop[0], kLoop[1], kLoop[2], kLoop[3]}};
    const BezierHalves halves = b.split(0.5);
    const CubicBezier& L = halves.left;
    const CubicBezier& R = halves.right;
    std::vector<Vec3> c = bezierToBSpline(L.b);
    c.push_back(R.b[3] * 6.0 - c[2] - c[3] * 4.0);
    YarnPath y = makeYarn(c, 0);
    CHECK(norm(y.curve.span(1).eval(0.5) - b.eval(0.75)) < 1e-12);
    CrossingResult r = crossingsOf(makeSet({y}));
    CHECK(r.unresolved.empty());
    CHECK(r.crossings.size() == 1);
    CHECK(r.crossings[0].uA < 1.0 && r.crossings[0].uB > 1.0);
    CHECK(std::fabs(r.crossings[0].uA - 1.0) > 0.1);
}

TEST(cross_smooth_join) {
    CurveSet cs = makeSet({makeYarn(wavyCtrl(), 0)});
    CrossingResult r = crossingsOf(cs);
    CHECK(r.crossings.empty());
    CHECK(r.unresolved.empty());
}

TEST(cross_open_end_extension) {
    // Line A ends at x = 0; its polynomial extension would cross B at x = 1e-3.
    CurveSet cs = makeSet({makeYarn(lineCtrl({-1, 0, 0}, {0, 0, 0}, 1), 0),
                           makeYarn(lineCtrl({1e-3, -1, 1}, {1e-3, 1, 1}, 1), 1)});
    CrossingResult r = crossingsOf(cs);
    CHECK(r.crossings.empty());
    CHECK(r.unresolved.empty());
}

TEST(cross_root_box_spans_knot) {
    CurveSet cs = makeSet({makeYarn(lineCtrl({0, 0, 0}, {4, 0, 0}, 4), 0),
                           makeYarn(lineCtrl({2, -1, 1}, {2, 1, 1}, 2), 1)});
    CrossingResult r = crossingsOf(cs);
    CHECK(r.unresolved.empty());
    CHECK(r.crossings.size() == 1);
    CHECK(r.crossings[0].boxA.lo < 2.0 && r.crossings[0].boxA.hi > 2.0);  // box spans the knot u = 2
    CHECK(r.crossings[0].status == Status::Pass);
}

TEST(cross_depth_contraction) {
    // A: z = x (slope 1); B: vertical at x = 0.3 with z = 0.3 + 1e-9.
    // The first root box (radius >= 1e-8) cannot resolve a 1e-9 depth gap.
    CurveSet cs = makeSet({makeYarn(lineCtrl({0, 0, 0}, {1, 0, 1}, 1), 0),
                           makeYarn(lineCtrl({0.3, -1, 0.3 + 1e-9}, {0.3, 1, 0.3 + 1e-9}, 1), 1)});
    CrossingResult r = crossingsOf(cs);
    CHECK(r.crossings.size() == 1);
    const Crossing& x = r.crossings[0];
    if (x.status == Status::Pass) CHECK(x.front == 1);
    else CHECK(x.front == -1);
    CHECK(x.status == Status::Pass);  // contraction is expected to succeed here
}

// -------------------------------------------------------------- clearance ----

TEST(clear_parallel_states) {
    PatternParams pp = testParams();
    auto statusAt = [&](double d) {
        CurveSet cs = makeSet({makeYarn(lineCtrl({0, 0, 0}, {5, 0, 0}, 5), 0),
                               makeYarn(lineCtrl({0, d, 0}, {5, d, 0}, 5), 1)});
        return clearanceOf(cs, pp).status;
    };
    CHECK(statusAt(1.15) == Status::Pass);
    CHECK(statusAt(1.05) == Status::Fail);
    CHECK(statusAt(1.1) == Status::Unresolved);
}

TEST(clear_intersecting) {
    PatternParams pp = testParams();
    CurveSet cs = makeSet({makeYarn(lineCtrl({0, 0, 0}, {4, 0, 0}, 4), 0),
                           makeYarn(lineCtrl({2, -2, 0}, {2, 2, 0}, 4), 1)});
    ClearanceResult r = clearanceOf(cs, pp);
    CHECK(r.status == Status::Fail);
    CHECK(r.minInterYarn < 1e-12);
    bool pen = false;
    for (const ClearanceRegion& g : r.regions) pen = pen || g.penetrating;
    CHECK(pen);
}

std::vector<Vec3> stadium(double s) {
    std::vector<Vec3> c = {{-1, 4, 0},  {-1, 3, 0},    {-1, 2, 0},   {-1, 1, 0}, {-1, 0, 0}, {-0.9, -1.0, 0},
                           {0, -1.6, 0}, {0.9, -1.0, 0}, {1, 0, 0},    {1, 1, 0},  {1, 2, 0},  {1, 3, 0},
                           {1, 4, 0}};
    for (Vec3& v : c) v *= s;
    return c;
}

TEST(clear_open_end_self) {
    // A J-shaped yarn whose open end comes back next to its own first leg.
    std::vector<Vec3> c = {{0, 12, 0}, {0, 10, 0}, {0, 8, 0},   {0, 6, 0},   {0, 4, 0},   {0, 2, 0},
                           {1, 0, 0},  {3, -1, 0}, {5, 0, 0},   {6, 2, 0},   {6, 4, 0},   {5, 5.5, 0},
                           {3.5, 6, 0}, {2, 6, 0}, {1.5, 6, 0}, {0.8, 6, 0}, {0.0, 6, 0}};
    CurveSet cs = makeSet({makeYarn(c, 0)});
    CHECK(cs.yarns[0].end().x < 1.0);
    ClearanceResult r = clearanceOf(cs, testParams());
    CHECK(r.status == Status::Fail);
    bool endHit = false;
    for (const ClearanceIssue& is : r.issues)
        if (is.yarnA == 0 && is.yarnB == 0 && is.status == Status::Fail && std::max(is.uA, is.uB) > cs.yarns[0].uEnd() - 0.5)
            endHit = true;
    CHECK(endHit);
}

TEST(clear_hairpin_ok) {
    CurveSet cs = makeSet({makeYarn(stadium(1.0), 0)});
    PatternParams pp = testParams();
    CurvatureResult k = checkCurvature(cs.yarns[0], pp.radius(), Allowances::forExtent(cs.extent()));
    CHECK(k.rhoLo >= pp.radius() + 0.5 * pp.clearance);
    ClearanceResult r = clearanceOf(cs, pp);
    CHECK(r.status == Status::Pass);
}

TEST(clear_hairpin_tight) {
    // Choose D, c around the measured hairpin so that r <= rho < r + c/2 and the
    // legs (2 mm apart) satisfy D <= d < D + c: the result must be UNRESOLVED.
    CurveSet cs = makeSet({makeYarn(stadium(1.0), 0)});
    double rho = checkCurvature(cs.yarns[0], 0.1, Allowances::forExtent(cs.extent())).rhoLo;
    PatternParams pp = testParams(1.98 * rho, 0.0);
    pp.clearance = std::max(2.05 - pp.yarnDiameter, 2.2 * (rho - pp.radius())) + 0.05;
    CHECK(pp.radius() <= rho && rho < pp.radius() + 0.5 * pp.clearance);
    CHECK(pp.yarnDiameter <= 2.0 && 2.0 < pp.yarnDiameter + pp.clearance);
    ClearanceResult r = clearanceOf(cs, pp);
    CHECK(r.status == Status::Unresolved);
}

TEST(clear_distant_self) {
    // 1.25 turns of a circle (radius 3) in one plane: start and end overlap.
    std::vector<Vec3> c;
    for (int i = -1; i <= 22; ++i) c.push_back({3.2 * std::cos(0.4 * i), 3.2 * std::sin(0.4 * i), 0});
    CurveSet cs = makeSet({makeYarn(c, 0)});
    ClearanceResult r = clearanceOf(cs, testParams());
    CHECK(r.status == Status::Fail);
    bool pen = false;
    for (const ClearanceRegion& g : r.regions) pen = pen || g.penetrating;
    CHECK(pen);
}

TEST(clear_window_straddle) {
    // Coarse sampling makes segment pairs span both L and pi*rho; subdivision must
    // still certify the legitimate hairpin (never FAIL).
    CurveSet cs = makeSet({makeYarn(stadium(1.0), 0)});
    ClearanceResult r = clearanceOf(cs, testParams(1.0, 0.1, 0.2));
    CHECK(r.status == Status::Pass);
}

// ------------------------------------------------------------ stockinette ----

TEST(stockinette_signature) {
    PatternParams pp = findPreset("expanded")->params;
    CurveSet cs = generateStockinette(pp);
    CrossingResult cr = crossingsOf(cs);
    SignatureResult s = checkSignature(cs, cr, pp);
    CHECK(s.status == Status::Pass);
    CHECK(static_cast<int>(s.pairs.size()) == pp.wales * (pp.courses - 1));
    CHECK(s.unexpected.empty());
    CHECK(static_cast<int>(cr.crossings.size()) == 4 * pp.wales * (pp.courses - 1));
}

TEST(stockinette_periodicity) {
    PatternParams pp = findPreset("expanded")->params;
    CurveSet cs = generateStockinette(pp);
    ValidationReport r = validate(cs, pp);
    CHECK(r.generator.periodicity == Status::Pass);
    // Regression: dense evaluation of a translated interior stitch.
    const YarnPath& y = cs.yarns[2];  // forward course
    const StitchSpan& a = y.stitches[1];
    const StitchSpan& b = y.stitches[3];
    Vec3 shift{(b.wale - a.wale) * pp.waleSpacing, 0, 0};
    for (int i = 0; i <= 400; ++i) {
        double t = i / 400.0;
        double ua = a.uBegin + t * (a.uEnd - a.uBegin), ub = b.uBegin + t * (b.uEnd - b.uBegin);
        CHECK(norm(y.curve.eval(ub) - y.curve.eval(ua) - shift) < 1e-11);
        CHECK(norm(y.curve.d1(ub) - y.curve.d1(ua)) < 1e-11);
    }
    // Reversed course reproduces the same geometry one row up.
    const YarnPath& z = cs.yarns[3];
    for (int i = 0; i <= 100; ++i) {
        double u = y.uEnd() * i / 100.0;
        CHECK(norm(z.curve.eval(z.uEnd() - u) - y.curve.eval(u) - Vec3{0, pp.courseSpacing, 0}) < 1e-11);
    }
}

TEST(stockinette_scaling_invariance) {
    // Projected crossings are preserved by the axis scaling (x,y) -> (w x, h y), so
    // their parameters must not depend on the spacings.
    auto sorted = [](CrossingResult r) {
        std::sort(r.crossings.begin(), r.crossings.end(), [](const Crossing& a, const Crossing& b) {
            if (a.yarnA != b.yarnA) return a.yarnA < b.yarnA;
            if (a.yarnB != b.yarnB) return a.yarnB < b.yarnB;
            return a.uA < b.uA;
        });
        return r;
    };
    PatternParams base = findPreset("expanded")->params;
    CurveSet c0 = generateStockinette(base);
    CrossingResult r0 = sorted(crossingsOf(c0));
    for (auto [w, h] : {std::pair{6.0, 6.0}, {10.0, 5.0}, {5.5, 7.5}}) {
        PatternParams p = base;
        p.waleSpacing = w;
        p.courseSpacing = h;
        CurveSet c = generateStockinette(p);
        CrossingResult r = sorted(crossingsOf(c));
        CHECK(checkSignature(c, r, p).status == Status::Pass);
        CHECK(r.crossings.size() == r0.crossings.size());
        for (size_t k = 0; k < r.crossings.size(); ++k) {
            CHECK_NEAR(r.crossings[k].uA, r0.crossings[k].uA, 1e-9);
            CHECK_NEAR(r.crossings[k].uB, r0.crossings[k].uB, 1e-9);
        }
    }
}

TEST(presets_expected_verdicts) {
    for (const Preset& p : presets()) {
        CurveSet cs = generateStockinette(p.params);
        ValidationReport r = validate(cs, p.params);
        if (r.verdict != p.expected) throw std::runtime_error("preset " + p.name + " verdict " + toString(r.verdict));
        if (p.expected == Status::Fail) {
            bool anyFail = false;
            for (const ClearanceRegion& g : r.clearance.regions) anyFail = anyFail || g.status == Status::Fail;
            CHECK(anyFail);
        }
    }
}

// ------------------------------------------------------------ Experiment R ----

// --------------------------------------------------------- Milestones 3 - 5 ----

TEST(rmf_frames) {
    // A curve with an inflection: the Frenet normal flips there, the RMF must not.
    std::vector<Vec3> c;
    for (int i = -1; i <= 20; ++i) {
        double x = 0.5 * i;
        c.push_back({x, std::sin(0.6 * x), 0.2 * x});
    }
    CurveSet cs = makeSet({makeYarn(c, 0)});
    PatternParams pp = testParams();
    SurfaceModel sm = buildSurface(cs, pp, SurfaceParams{});
    const YarnSurface& ys = sm.yarns[0];
    double maxTurn = 0, maxTwist = 0;
    for (size_t k = 0; k < ys.rings.size(); ++k) {
        const Ring& R = ys.rings[k];
        CHECK(std::fabs(norm(R.n1) - 1) < 1e-12 && std::fabs(dot(R.n1, R.T)) < 1e-12);
        CHECK(std::fabs(dot(cross(R.T, R.n1), R.n2) - 1) < 1e-12);
        if (k == 0) continue;
        const Ring& P = ys.rings[k - 1];
        maxTurn = std::max(maxTurn, std::acos(std::clamp(dot(P.n1, R.n1), -1.0, 1.0)));
        maxTwist = std::max(maxTwist, std::fabs(dot(R.n1 - P.n1, P.n2)) / ys.ds);  // rotation about T
    }
    CHECK(maxTurn < 0.05);   // continuous: no flip at the inflection
    CHECK(maxTwist < 1e-2);  // rotation-minimizing to O(ds)
}

TEST(tube_uv_seamless) {
    PatternParams pp = findPreset("expanded")->params;
    CurveSet cs = generateStockinette(pp);
    SurfaceModel sm = buildSurface(cs, pp, SurfaceParams{});
    TubeParams tp;
    TubeMesh m = buildTubeMesh(cs, sm, tp);
    const int stride = tp.sides + 1;
    double lastS = -1;
    for (size_t v0 = m.yarnFirstVertex[0]; v0 < m.yarnFirstVertex[1]; v0 += stride) {
        CHECK(norm(m.pos[v0] - m.pos[v0 + tp.sides]) < 1e-12);    // seam closes
        CHECK(m.uv[2 * v0 + 1] == 0.0f && m.uv[2 * (v0 + tp.sides) + 1] == 1.0f);
        CHECK(m.s[v0] > lastS);                                    // u follows arc length across spans
        CHECK(std::fabs(m.uv[2 * v0] - m.s[v0] / tp.uvRepeat) < 1e-5);
        CHECK(std::fabs(m.t[v0] - m.s[v0] / cs.yarns[0].length()) < 1e-12);
        lastS = m.s[v0];
    }
    CHECK(std::fabs(lastS - cs.yarns[0].length()) < 1e-9);
}

TEST(compression_two_sided) {
    // Parallel yarns 0.8 mm apart, D = 1: both sides flattened to b = 0.4 so they touch.
    PatternParams pp = testParams();
    CurveSet cs = makeSet({makeYarn(lineCtrl({0, 0, 0}, {6, 0, 0}, 6), 0), makeYarn(lineCtrl({0, 0.8, 0}, {6, 0.8, 0}, 6), 1)});
    SurfaceModel sm = buildSurface(cs, pp, SurfaceParams{});
    const CompressionReport& rep = sm.compression;
    CHECK(rep.clamped == 0);
    CHECK_NEAR(rep.minAxisRatio, 0.8, 1e-9);
    CHECK_NEAR(rep.maxAreaChange, 0.2, 1e-9);  // the free axis keeps r: the section loses area
    for (int y = 0; y < 2; ++y) {
        Ring R = ringAt(cs, sm.yarns[y], 0.5 * cs.yarns[y].length());
        CHECK_NEAR(R.b, 0.4, 1e-9);
        CHECK_NEAR(R.a, 0.5, 1e-9);
        CHECK_NEAR(ringSupport(R, {0, y == 0 ? 1.0 : -1.0, 0}), 0.4, 1e-9);  // contact kept, no gap created
    }
}

TEST(compression_only_where_overlapping) {
    PatternParams pp = testParams();
    CurveSet cs = makeSet({makeYarn(lineCtrl({0, 0, 0}, {6, 0, 0}, 6), 0), makeYarn(lineCtrl({0, 1.05, 0}, {6, 1.05, 0}, 6), 1)});
    SurfaceModel sm = buildSurface(cs, pp, SurfaceParams{});
    CHECK(sm.compression.clamped == 0);
    CHECK(sm.compression.minAxisRatio == 1.0);
    CHECK(sm.compression.maxAreaChange == 0.0);  // still separated: no contact manufactured
}

struct FibreFixture {
    PatternParams pp = findPreset("tightest")->params;
    CurveSet cs = generateStockinette(pp);
    static FibreParams params() {
        FibreParams fp = pieceFibres();
        fp.density = 8;
        fp.rootDepth = 0.126;
        return fp;
    }
    FibreSet run(double eps, FibreParams fp = params(), const ExecPolicy& exec = {}) {
        PatternParams p = pp;
        p.sampleTolerance = eps;
        SurfaceModel sm = buildSurface(cs, p, SurfaceParams{});
        FibreFrame fr;  // only the execution policy: behaves as no frame
        fr.exec = exec;
        return growYarnFibres(cs, sm, p, fp, &fr);
    }
};

TEST(fibers_deterministic) {
    FibreFixture f;
    FibreSet a = f.run(0.005), b = f.run(0.002);  // different polyline sampling
    CHECK(a.fibres.size() == b.fibres.size() && !a.fibres.empty());
    for (size_t i = 0; i < a.fibres.size(); ++i) {
        CHECK(a.fibres[i].s == b.fibres[i].s && a.fibres[i].theta == b.fibres[i].theta);
        CHECK(norm(toVec3(a.fibres[i].anchor) - toVec3(b.fibres[i].anchor)) < 1e-9);
    }
    FibreParams other = FibreFixture::params();
    other.seed = 7;
    FibreSet c = f.run(0.005, other);
    CHECK(c.fibres.size() != a.fibres.size() || c.fibres[0].s != a.fibres[0].s);
}

TEST(fibre_islands_cover_what_they_ask) {
    // No rate, no island; a chance of 1 per stitch covers the yarn; a kind with no share never appears;
    // nothing passes the strength.
    FibreFixture f;
    FibreIslands islands{0, 8, 60, 1, 0, 1};
    const auto covered = [&]() {
        size_t weighted = 0, controls = 0;
        float largest = 0, matted = 0;
        for (const std::vector<IslandWeight>& strand : islandWeights({f.cs, f.cs.topo, islands, 18.7, 9}))
            for (const IslandWeight& w : strand) {
                ++controls;
                weighted += w.fuzzy > 0 || w.worn > 0;
                largest = std::max({largest, w.fuzzy, w.worn});
                matted = std::max(matted, w.matted);
            }
        CHECK(matted == 0 && largest <= 0.6f);
        return double(weighted) / double(controls);
    };
    CHECK(covered() == 0);
    islands.rate = 100.0 / 18.7;
    const double all = covered();
    islands.rate = 0.4;
    const double some = covered();
    std::printf("    covered: all %.3f, at 0.4 per cm2 %.3f\n", all, some);
    CHECK(all > 0.9 && some > 0.02 && some < 0.8 && some < all);
}

TEST(fibre_islands_zero_weight_unchanged) {
    // Where the weight is zero a fibre is the one the even fuzz has there, bit for bit: same candidates,
    // same points. Islands change fibres only inside themselves.
    FibreFixture f;
    const SurfaceModel sm = buildSurface(f.cs, f.pp, SurfaceParams{});
    const FibreParams even = FibreFixture::params();
    FibreParams patched = even;
    patched.islands = {0.6, 10, 100, 1, 1, 1};
    FibreFrame plainFrame, islandFrame;
    islandFrame.islands = islandWeights({f.cs, f.cs.topo, patched.islands, 18.7, 21});
    const FibreSet a = growYarnFibres(f.cs, sm, f.pp, even, &plainFrame);
    const FibreSet b = growYarnFibres(f.cs, sm, f.pp, patched, &islandFrame);
    std::map<std::pair<int, int>, const Fibre*> plain;
    for (const Fibre& fb : a.fibres) plain[{fb.yarnId, fb.candidate}] = &fb;
    int same = 0, inside = 0;
    for (const Fibre& fb : b.fibres) {
        const std::vector<IslandWeight>& strand = islandFrame.islands[fb.strand];
        const IslandWeight w = islandAt(strand, f.cs.yarns[fb.strand].paramAtArcLength(fb.s));
        if (w.fuzzy > 0 || w.matted > 0 || w.worn > 0) {
            ++inside;
            continue;
        }
        const auto it = plain.find({fb.yarnId, fb.candidate});
        CHECK(it != plain.end());
        if (it == plain.end()) continue;
        const std::span<const Vec3f> p = a.pointsOf(*it->second), q = b.pointsOf(fb);
        CHECK(p.size() == q.size() && std::memcmp(p.data(), q.data(), p.size_bytes()) == 0);
        ++same;
    }
    std::printf("    %d fibres outside the islands equal, %d inside\n", same, inside);
    CHECK(same > 100 && inside > 20);
}

TEST(hand_columns_slide_apart) {
    // A column moves together, up its wale, and differs from the columns beside it; the yarn stays whole
    // where it crosses from one column to the next.
    PatternParams pp = testParams();
    pp.wales = 12;
    pp.courses = 10;
    CurveSet cs = generateStockinette(pp);
    for (StitchRecord& st : cs.topo.stitches) st.lattice = {0, st.index, st.turn, st.turn};
    const Controls before = controlsOf(cs);
    HandRequest hand;
    hand.settings.columnShift = 1.0;
    hand.settings.columnRows = 6;
    hand.seed = 4;
    hand.wale = pp.waleSpacing;
    knitByHand(cs, hand, knit_test::cpuNumerics(), {});
    const Controls after = controlsOf(cs);
    const auto slideOf = [&](int wale, int row) {  // the stitch's mean motion along V
        const StitchRecord& st = cs.topo.stitches[row * pp.wales + wale];
        const int first = static_cast<int>(st.uBegin) + 1;
        double sum = 0;
        for (int j = 0; j < kControlsPerStitch; ++j) sum += dot(after[st.yarn][first + j] - before[st.yarn][first + j], st.V);
        return sum / kControlsPerStitch;
    };
    double alongColumn = 0, acrossColumns = 0, largestStep = 0;
    for (int r = 2; r < 7; ++r)
        for (int w = 2; w < 9; ++w) {
            alongColumn = std::max(alongColumn, std::abs(slideOf(w, r + 1) - slideOf(w, r)));
            acrossColumns = std::max(acrossColumns, std::abs(slideOf(w + 1, r) - slideOf(w, r)));
        }
    for (size_t y = 0; y < after.size(); ++y)
        for (size_t c = 2; c + 3 < after[y].size(); ++c)
            largestStep = std::max(largestStep, norm(after[y][c + 1] - after[y][c]) - norm(before[y][c + 1] - before[y][c]));
    std::printf("    a row up %.3f mm, a column across %.3f mm, longest edge +%.3f mm\n", alongColumn, acrossColumns, largestStep);
    CHECK(acrossColumns > 2 * alongColumn && acrossColumns > 0.3);
    CHECK(largestStep < 0.6);
}

TEST(fibers_front_only) {
    FibreFixture f;
    const FibreParams fp = FibreFixture::params();
    FibreSet a = f.run(0.005, fp);
    CHECK(a.rep.backFacing > 0);
    for (const Fibre& fb : a.fibres) CHECK(fb.normal.z > fp.frontCos);
}

TEST(fibers_occlusion) {
    // Yarn 1 lies directly in front (+N) of yarn 0: yarn 0's top is occluded.
    PatternParams pp = testParams();
    CurveSet cs = makeSet({makeYarn(lineCtrl({0, 0, 0}, {8, 0, 0}, 8), 0), makeYarn(lineCtrl({0, 0, 1.5}, {8, 0, 1.5}, 8), 1)});
    SurfaceModel sm = buildSurface(cs, pp, SurfaceParams{});
    FibreParams fp = FibreFixture::params();
    fp.density = 20;
    FibreSet fs = growYarnFibres(cs, sm, pp, fp);
    CHECK(fs.rep.occluded > 0);
    for (const Fibre& fb : fs.fibres)
        if (fb.yarnId == 0) CHECK(std::fabs(fb.anchor.y) > 0.4);  // on the sides only, never under yarn 1's tube
}

// --------------------------------------------- surface and fibre validation ----

TEST(surface_fold_detected) {
    CurveSet cs = makeSet({makeYarn(stadium(1.0), 0)});
    double rho = checkCurvature(cs.yarns[0], 0.1, Allowances::forExtent(cs.extent())).rhoLo;
    PatternParams fat = testParams(2.2 * rho, 0.0), thin = testParams(1.6 * rho, 0.0);
    SurfaceParams sp;
    sp.compression = false;  // round tubes: the compression flattens the fold away (compression_fold_limit)
    SurfaceCheck a = checkSurface(cs, buildSurface(cs, fat, sp));
    SurfaceCheck b = checkSurface(cs, buildSurface(cs, thin, sp));
    CHECK(a.fold == Status::Fail && a.minFoldRatio < 1 && a.foldRings > 0);
    CHECK(b.fold == Status::Pass && b.minFoldRatio > 1);
}

TEST(compression_fold_limit) {
    // surface_fold_detected's fat yarn (r = 1.1 rho on the arcs): the fit flattens the
    // cross-section along the principal normal to rho / 1.02 (above the clamp 0.6 r), so the
    // swept surface no longer folds; the round tube keeps the fold.
    CurveSet cs = makeSet({makeYarn(stadium(1.0), 0)});
    double rho = checkCurvature(cs.yarns[0], 0.1, Allowances::forExtent(cs.extent())).rhoLo;
    PatternParams fat = testParams(2.2 * rho, 0.0);
    SurfaceParams round;
    round.compression = false;
    CHECK(checkSurface(cs, buildSurface(cs, fat, round)).fold == Status::Fail);
    SurfaceCheck b = checkSurface(cs, buildSurface(cs, fat, SurfaceParams{}));
    CHECK(b.fold == Status::Pass && b.minFoldRatio > 1.0 && b.minFoldRatio < 1.03);
}

TEST(compression_parallel_066) {
    // Two parallel yarns at 0.66 D (D = 1): each flattens to b = 0.33 = 0.66 r (above the clamp
    // 0.6 r), the drawn surfaces touch and nothing penetrates; the free axis keeps r.
    PatternParams pp = testParams();
    CurveSet cs = makeSet({makeYarn(lineCtrl({0, 0, 0}, {6, 0, 0}, 6), 0), makeYarn(lineCtrl({0, 0.66, 0}, {6, 0.66, 0}, 6), 1)});
    SurfaceModel sm = buildSurface(cs, pp, SurfaceParams{});
    CHECK(sm.compression.clamped == 0);
    CHECK_NEAR(ringAt(cs, sm.yarns[0], 3.0).b, 0.33, 1e-9);
    CHECK_NEAR(ringAt(cs, sm.yarns[0], 3.0).a, 0.5, 1e-9);
    SurfaceCheck chk = checkSurface(cs, sm);
    CHECK(chk.penetratingPoints == 0 && chk.collisions == Status::Pass);
}

TEST(compression_third_yarn) {
    // Yarn 0 between yarn 1 (0.8 along y) and yarn 2 along z. The fit bounds the support towards
    // every contact within 2r, so at z = 0.9 yarn 0 flattens for yarn 2 too, and at z = 1.05 its
    // free axis stays out of yarn 2: nothing penetrates.
    PatternParams pp = testParams();
    for (double z : {1.05, 0.9}) {
        CurveSet cs = makeSet({makeYarn(lineCtrl({0, 0, 0}, {6, 0, 0}, 6), 0), makeYarn(lineCtrl({0, 0.8, 0}, {6, 0.8, 0}, 6), 1),
                               makeYarn(lineCtrl({0, 0, z}, {6, 0, z}, 6), 2)});
        const SurfaceModel sm = buildSurface(cs, pp, SurfaceParams{});
        const SurfaceCheck chk = checkSurface(cs, sm);
        CHECK(chk.penetratingPoints == 0 && chk.collisions == Status::Pass);
        CHECK(sm.compression.clamped == 0);
        // Yarn 0 keeps b = 0.4 towards yarn 1 and its support towards yarn 2 at most half the distance.
        const Ring R = ringAt(cs, sm.yarns[0], 3.0);
        CHECK_NEAR(R.b, 0.4, 1e-9);
        CHECK(ringSupport(R, {0, 0, 1}) <= 0.5 * z + 1e-9);
    }
}

TEST(surface_new_collision) {
    // Yarns 0 and 1 touch (0.8 apart), set by hand to ellipses that keep their area (b = 0.4,
    // a = r^2 / b = 0.625) and widen along z into yarn 2, 1.05 above yarn 0. No requirement
    // predicted that collision; the surface check must find it.
    PatternParams pp = testParams();
    CurveSet cs = makeSet({makeYarn(lineCtrl({0, 0, 0}, {6, 0, 0}, 6), 0), makeYarn(lineCtrl({0, 0.8, 0}, {6, 0.8, 0}, 6), 1),
                           makeYarn(lineCtrl({0, 0, 1.05}, {6, 0, 1.05}, 6), 2)});
    SurfaceParams round;
    round.compression = false;
    SurfaceModel sm = buildSurface(cs, pp, round);
    for (int y = 0; y < 2; ++y)
        for (Ring& R : sm.yarns[y].rings) {
            R.a = 0.625;
            R.b = 0.4;
            R.phi = std::atan2(R.n2.z, R.n1.z);  // the a-axis along z
        }
    SurfaceCheck chk = checkSurface(cs, sm);
    CHECK(chk.collisions == Status::Fail);
    CHECK(chk.penetratingPoints > 0);
    CHECK(chk.maxDepth > 0.05 && chk.maxDepth < 0.1);  // ~0.075 mm expected
}

TEST(presets_compressed_surface_valid) {
    // The compressed preset is INVALID for nominal round tubes (M1) but its M3
    // compressed surface passes the independent surface check.
    PatternParams pp = findPreset("compressed")->params;
    CurveSet cs = generateStockinette(pp);
    ValidationReport r = validate(cs, pp);
    CHECK(r.verdict == Status::Fail);
    SurfaceModel sm = buildSurface(cs, pp, SurfaceParams{});
    CHECK(sm.compression.minAxisRatio < 1);
    SurfaceCheck chk = checkSurface(cs, sm);
    CHECK(chk.fold == Status::Pass);
    CHECK(chk.collisions == Status::Pass);
}

// ------------------------------------------------ certificate, bump, study ----

TEST(motion_certificate) {
    PatternParams pp = findPreset("tightest")->params;
    CurveSet a = generateStockinette(pp);
    ValidationReport ra = validate(a, pp, {.generatorChecks = false});
    auto moved = [&](double dz) {
        CurveSet b = a;
        for (YarnPath& y : b.yarns) {
            for (Vec3& q : y.curve.ctrl) q.z += dz;  // rigid motion: valid end state
            y.finalize();
        }
        return b;
    };
    CurveSet b1 = moved(0.05), b2 = moved(2.0);
    MotionCertificate c1 = certifyStraightMotion(a, ra, b1, validate(b1, pp, {.generatorChecks = false}), pp);
    MotionCertificate c2 = certifyStraightMotion(a, ra, b2, validate(b2, pp, {.generatorChecks = false}), pp);
    CHECK(c1.endStatesValid && c1.noPassThrough && c1.noTubeContact);
    CHECK(c2.endStatesValid && !c2.noPassThrough);  // 2 mm > D + c: not certifiable this way
    CHECK_NEAR(c1.maxDisplacement, 0.05, 1e-12);
    // An invalid end state is never certified.
    PatternParams bad = findPreset("invalid")->params;
    CurveSet x = generateStockinette(bad);
    ValidationReport rx = validate(x, bad, {.generatorChecks = false});
    CHECK(!certifyStraightMotion(x, rx, x, rx, bad).noPassThrough);
}

// ------------------------------------------------------------ chart skew ----

namespace {

// An icosahedron split `levels` times onto the unit sphere, scaled and stretched along z: a closed genus-0 mesh that was
// never a UV-sphere grid.
TriMesh icosphere(int levels, double radius, double stretchZ) {
    TriMesh m;
    const double g = (1 + std::sqrt(5.0)) / 2;
    m.pos = {{-1, g, 0}, {1, g, 0}, {-1, -g, 0}, {1, -g, 0}, {0, -1, g}, {0, 1, g}, {0, -1, -g}, {0, 1, -g}, {g, 0, -1}, {g, 0, 1}, {-g, 0, -1}, {-g, 0, 1}};
    m.tri = {{0, 11, 5}, {0, 5, 1}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11}, {1, 5, 9}, {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
             {3, 9, 4}, {3, 4, 2}, {3, 2, 6}, {3, 6, 8}, {3, 8, 9}, {4, 9, 5}, {2, 4, 11}, {6, 2, 10}, {8, 6, 7}, {9, 8, 1}};
    for (int l = 0; l < levels; ++l) {
        std::map<std::pair<int, int>, int> middle;
        std::vector<std::array<int, 3>> next;
        for (const std::array<int, 3>& t : m.tri) {
            int mid[3];
            for (int k = 0; k < 3; ++k) {
                const int u = std::min(t[k], t[(k + 1) % 3]), v = std::max(t[k], t[(k + 1) % 3]);
                std::map<std::pair<int, int>, int>::iterator it = middle.find({u, v});
                if (it == middle.end()) {
                    it = middle.emplace(std::make_pair(u, v), static_cast<int>(m.pos.size())).first;
                    m.pos.push_back((m.pos[u] + m.pos[v]) * 0.5);
                }
                mid[k] = it->second;
            }
            next.push_back({t[0], mid[0], mid[2]});
            next.push_back({t[1], mid[1], mid[0]});
            next.push_back({t[2], mid[2], mid[1]});
            next.push_back({mid[0], mid[1], mid[2]});
        }
        m.tri = next;
    }
    for (Vec3& p : m.pos) {
        p = p * (radius / norm(p));
        p.z *= stretchZ;
    }
    m.computeNormals();
    return m;
}

} // namespace

// chartSkew (chart.h): a cylinder band whose columns run straight up is 0 deg off perpendicular; columns that wind
// (b shifted by s a, on a band as tall as it is round) are atan(s) off, the angle between grad a and grad b off 90.
TEST(chart_skew_measures_sheared_columns) {
    auto band = [](double shear) {
        const int rings = 32, columns = 128;
        const double radius = 10, height = 2 * kPi * radius;
        TriMesh m;
        for (int i = 0; i <= rings; ++i)
            for (int j = 0; j < columns; ++j) {
                const double a = 0.1 + 0.8 * i / rings, angle = 2 * kPi * (double(j) / columns + shear * a);
                m.pos.push_back({radius * std::cos(angle), radius * std::sin(angle), height * a});
            }
        for (int i = 0; i < rings; ++i)
            for (int j = 0; j < columns; ++j) {
                const int v00 = i * columns + j, v01 = i * columns + (j + 1) % columns, v10 = v00 + columns, v11 = v01 + columns;
                const double a0 = 0.1 + 0.8 * i / rings, a1 = 0.1 + 0.8 * (i + 1) / rings;
                const double b0 = double(j) / columns, b1 = double(j + 1) / columns;  // the seam's far side holds b = 1
                m.tri.push_back({v00, v01, v11});
                m.uv.push_back({Vec2{a0, b0}, Vec2{a0, b1}, Vec2{a1, b1}});
                m.tri.push_back({v00, v11, v10});
                m.uv.push_back({Vec2{a0, b0}, Vec2{a1, b1}, Vec2{a1, b0}});
            }
        return chartSkew(m);
    };
    const ChartSkew straight = band(0), sheared = band(0.25);
    const double expected = std::atan(0.25) * 180 / kPi;
    std::printf("    straight %.4f deg (95th %.4f), sheared %.4f deg (95th %.4f), expected %.4f\n", straight.mean, straight.p95, sheared.mean,
                sheared.p95, expected);
    CHECK(straight.mean < 0.01 && straight.p95 < 0.01);
    CHECK(std::fabs(sheared.mean - expected) < 0.3 && std::fabs(sheared.p95 - expected) < 0.3);
}

TEST(topology_flat_matches_spans) {
    // The generator's explicit topology equals the one synthesized from (wale,
    // course) spans, and u -> stitch lookups agree with YarnPath.
    PatternParams pp = findPreset("tightest")->params;
    CurveSet cs = generateStockinette(pp);
    CurveSet hand = cs;
    fillFlatTopology(hand);
    const KnitTopology& syn = hand.topo;
    CHECK(syn.stitches.size() == cs.topo.stitches.size());
    for (size_t k = 0; k < syn.stitches.size(); ++k) {
        const StitchRecord &a = syn.stitches[k], &b = cs.topo.stitches[k];
        CHECK(a.yarn == b.yarn && a.turn == b.turn && a.index == b.index && a.nParents == b.nParents);
        CHECK(a.uBegin == b.uBegin && a.uEnd == b.uEnd);
        if (a.nParents) CHECK(a.parents[0] == b.parents[0]);
    }
    for (size_t y = 0; y < cs.yarns.size(); ++y)
        for (double u = 0; u <= cs.yarns[y].uEnd(); u += 0.37) {
            int s = cs.yarns[y].stitchContaining(u, u);
            int t = cs.topo.stitchAt(static_cast<int>(y), u);
            CHECK((s < 0 && t < 0) || (s >= 0 && cs.yarns[y].stitches[s].id == t));
        }
    CHECK(cs.topo.isParent(0, pp.wales) && !cs.topo.isParent(pp.wales, 0));
}

TEST(hash_matches_bruteforce) {
    // The hash only enumerates candidates: a query visits every box that overlaps the query box (closed,
    // componentwise), each once, with the map and with the dense lookup. So querying each box inflated by
    // reach and applying the callers' exact box-distance filter gives the all-pairs loop's pairs, in its order.
    unsigned seed = 12345;
    auto rnd = [&]() { seed = seed * 1664525u + 1013904223u; return (seed >> 8) / double(1 << 24); };
    for (int trial = 0; trial < 6; ++trial) {
        const int n = 300;
        std::vector<Vec3> lo(n), hi(n);
        for (int i = 0; i < n; ++i) {
            Vec3 c{40 * rnd() - 20, 40 * rnd() - 20, (trial % 2 ? 40 : 2) * rnd()};
            Vec3 e{3 * rnd() * (trial + 1) / 3, 3 * rnd(), 0.5 * rnd()};
            lo[i] = c - e;
            hi[i] = c + e;
        }
        const double reach = 0.4 + trial;
        auto boxDist = [&](int a, int b) {
            double dx = std::max({0.0, lo[a].x - hi[b].x, lo[b].x - hi[a].x});
            double dy = std::max({0.0, lo[a].y - hi[b].y, lo[b].y - hi[a].y});
            double dz = std::max({0.0, lo[a].z - hi[b].z, lo[b].z - hi[a].z});
            return std::sqrt(dx * dx + dy * dy + dz * dz);
        };
        std::vector<std::pair<int, int>> brute, hashed;
        for (int i = 0; i < n; ++i)
            for (int j = i; j < n; ++j)
                if (boxDist(i, j) <= reach) brute.push_back({i, j});
        double diag = 0;
        for (int i = 0; i < n; ++i) diag += norm(hi[i] - lo[i]);
        SpatialHash3 hash(std::max(diag / n, 1e-9) + reach);
        for (int i = 0; i < n; ++i) hash.insert(lo[i], hi[i]);
        if (trial % 2 == 1) hash.finalize();  // the dense lookup must visit the same items
        const Vec3 r{reach, reach, reach};
        bool sameCandidates = true;
        for (int i = 0; i < n; ++i) {
            const Vec3 qLo = lo[i] - r, qHi = hi[i] + r;
            std::vector<int> visited;
            hash.query(qLo, qHi, [&](int j) { visited.push_back(j); });
            std::sort(visited.begin(), visited.end());
            std::vector<int> overlapping;  // by brute force; an id visited twice makes the lists differ
            for (int j = 0; j < n; ++j)
                if (!(lo[j].x > qHi.x || hi[j].x < qLo.x || lo[j].y > qHi.y || hi[j].y < qLo.y || lo[j].z > qHi.z ||
                      hi[j].z < qLo.z))
                    overlapping.push_back(j);
            if (visited != overlapping) sameCandidates = false;
            for (int j : visited)
                if (j >= i && boxDist(i, j) <= reach) hashed.push_back({i, j});
        }
        CHECK(sameCandidates);
        CHECK(brute.size() > static_cast<size_t>(n));
        CHECK(brute == hashed);
    }
}

TEST(schur_local_tight_bend) {
    // The stadium hairpin of clear_hairpin_tight with long, straight (gentle) legs.
    // Leg pairs far from the bend have gentle curvature at both ends, but their
    // intervening arc runs through the tight bend: the local window must use the
    // tight rho there (UNRESOLVED, as for the short hairpin), never the legs' own
    // rho (which would demand D + c and report a FAIL).
    std::vector<Vec3> c = stadium(1.0);
    for (int k = 5; k <= 24; ++k) {
        c.insert(c.begin(), Vec3{-1, double(k), 0});
        c.push_back({1, double(k), 0});
    }
    CurveSet cs = makeSet({makeYarn(c, 0)});
    CurvatureResult k = checkCurvature(cs.yarns[0], 0.1, Allowances::forExtent(cs.extent()));
    double rho = k.rhoLo;
    PatternParams pp = testParams(1.98 * rho, 0.0);
    pp.clearance = std::max(2.05 - pp.yarnDiameter, 2.2 * (rho - pp.radius())) + 0.05;
    CHECK(pp.radius() <= rho && rho < pp.radius() + 0.5 * pp.clearance);
    // Per-span bounds: the legs are far gentler than the bend.
    double legKappa = 0;
    for (int s = 0; s < 10; ++s) legKappa = std::max(legKappa, k.spanKappaHi[s]);
    CHECK(legKappa < 0.01 / rho);
    ClearanceResult r = clearanceOf(cs, pp);
    CHECK(r.status == Status::Unresolved);
    for (const ClearanceRegion& g : r.regions) CHECK(g.status != Status::Fail);
    // Issues are reported far up the legs too (not only next to the bend).
    double maxY = 0;
    for (const ClearanceIssue& is : r.issues) maxY = std::max(maxY, is.pA.y);
    CHECK(maxY > 15);
}

TEST(length_bounds_prefix) {
    // Prefix-sum length bounds on a long yarn: they enclose the integrated length and
    // stay within a tiny margin of the direct per-span sums (partial end spans exact).
    std::vector<Vec3> c;
    for (int i = -1; i <= 3001; ++i) c.push_back({0.7 * i, 2.0 * std::sin(0.3 * i), 1.5 * std::cos(0.17 * i)});
    YarnPath y = makeYarn(c, 0);
    auto direct = [&](double u0, double u1, double& lo, double& hi) {
        lo = hi = 0;
        for (int k = static_cast<int>(std::floor(u0)); k < y.curve.numSpans() && k < u1; ++k) {
            double t0 = std::max(0.0, u0 - k), t1 = std::min(1.0, u1 - k);
            if (t1 <= t0) continue;
            CubicBezier p = y.curve.span(k).sub(t0, t1);
            hi += p.polygonLength();
            Vec3 prev = p.b[0];
            for (int j = 1; j <= 8; ++j) { Vec3 q = p.eval(j / 8.0); lo += norm(q - prev); prev = q; }
        }
    };
    unsigned seed = 7;
    auto rnd = [&]() { seed = seed * 1664525u + 1013904223u; return (seed >> 8) / double(1 << 24); };
    const double U = y.uEnd();
    std::vector<std::pair<double, double>> cases = {{0, U}, {0.0, 1.0}, {2.0, 5.0}, {2.3, 5.0}, {2.3, 2.7}, {U - 1.5, U}};
    for (int i = 0; i < 200; ++i) {
        double a = rnd() * U, b = rnd() * U;
        cases.push_back({std::min(a, b), std::max(a, b)});
    }
    for (auto [u0, u1] : cases) {
        const ArcLengthBounds b = y.lengthBounds(u0, u1);
        const double lo = b.lo, hi = b.hi;
        double dlo, dhi;
        direct(u0, u1, dlo, dhi);
        double L = y.arcLengthAt(u1) - y.arcLengthAt(u0);
        CHECK(lo <= L + 1e-9 && L <= hi + 1e-9);
        CHECK(lo <= dlo && hi >= dhi);                      // the allowance only widens
        CHECK(dlo - lo <= 1e-9 * (1 + y.length()) && hi - dhi <= 1e-9 * (1 + y.length()));
    }
    // stitchContaining by binary search agrees with a scan over many stitches.
    YarnPath z = makeYarn(c, 0);
    z.stitches.clear();
    for (int k = 0; k < 300; ++k) z.stitches.push_back({k, 0, 10.0 * k, 10.0 * (k + 1), k});
    z.finalize();
    for (int i = 0; i < 2000; ++i) {
        double u0 = rnd() * z.uEnd(), u1 = std::min(z.uEnd(), u0 + 3 * rnd());
        int want = -1;
        for (size_t k = 0; k < z.stitches.size(); ++k) {
            const StitchSpan& st = z.stitches[k];
            bool last = st.uEnd >= z.uEnd();
            if (u0 >= st.uBegin && (u1 < st.uEnd || (last && u1 <= st.uEnd))) { want = int(k); break; }
        }
        CHECK(z.stitchContaining(u0, u1) == want);
    }
    CHECK(z.stitchContaining(z.uEnd(), z.uEnd()) == 299);
}

// ------------------------------------------------------- S1: spiral pieces ----

TEST(stitch_queue_balance) {
    // Programs from feasible schedules pass the graph checks; corruptions are caught.
    std::vector<std::vector<int>> cases = {{6, 9, 13, 19, 24, 24, 20, 14, 8}, {12, 12, 12, 12}, {10, 15, 22, 11, 6}};
    for (const auto& counts : cases) {
        StitchProgram p = buildProgram(counts);
        GraphCheck g = checkProgram(p);
        if (g.status != Status::Pass) throw std::runtime_error(g.detail);
        int total = 0;
        for (int n : counts) total += n;
        CHECK(g.stitches == total);
        for (size_t r = 1; r < counts.size(); ++r) {
            const TurnPlan& tp = p.turns[r];
            CHECK(tp.input == counts[r - 1] && tp.output == counts[r]);
            CHECK(tp.inc == 0 || tp.dec == 0);
        }
    }
    // A stitch consuming a loop twice / skipping one is detected.
    StitchProgram p = buildProgram({8, 8, 8});
    p.stitches[p.turnStart[2] + 3].parents[0] = p.stitches[p.turnStart[2] + 2].parents[0];
    CHECK(checkProgram(p).status == Status::Fail);
}

TEST(count_schedule_feasible) {
    // Sphere-like targets (fast growth from 6, then shrinking): every turn respects the
    // density limit, and the program built from the schedule balances.
    std::vector<double> t;
    for (int r = 0; r < 60; ++r) t.push_back(6 + 120 * std::sin(kPi * r / 59.0));
    std::vector<int> n = scheduleCounts(6, t);
    CHECK(n.size() == t.size() && n[0] == 6);
    for (size_t r = 1; r < n.size(); ++r) CHECK(std::abs(n[r] - n[r - 1]) <= n[r - 1] / 2);
    CHECK(*std::max_element(n.begin(), n.end()) > 100);  // the deficit is caught up
    GraphCheck g = checkProgram(buildProgram(n));
    if (g.status != Status::Pass) throw std::runtime_error(g.detail);
    CHECK(g.inc > 0 && g.dec > 0);
}

TEST(dec_turn_boundary) {
    // Halving a turn: every DEC takes two consecutive loops of the SAME input turn; the
    // last DEC ends exactly on the last loop of that turn.
    StitchProgram p = buildProgram({10, 5, 5});
    GraphCheck g = checkProgram(p);
    if (g.status != Status::Pass) throw std::runtime_error(g.detail);
    const int s0 = p.turnStart[1], s1 = p.turnStart[2];
    for (int s = s0; s < s1; ++s) {
        const StitchRecord& st = p.stitches[s];
        CHECK(st.op == StitchOp::Dec && st.nParents == 2 && st.parents[1] == st.parents[0] + 1);
        CHECK(st.parents[1] < s0);  // never a loop of the current turn
    }
    CHECK(p.stitches[s1 - 1].parents[1] == s0 - 1);
    // A DEC that mixes generations is rejected.
    StitchProgram bad = buildProgram({10, 5, 5});
    bad.stitches[s1 - 1].parents[1] = s1;  // a loop of turn 2's generation
    GraphCheck gb = checkProgram(bad);
    CHECK(gb.status == Status::Fail);
}

TEST(planar_mapping_equals_flat) {
    // On a plane with a regular grid layout (nodes at the template origins of the flat
    // generator), the surface placement path reproduces the flat controls exactly,
    // including the parent-frame blend of the lower controls.
    PatternParams pp = findPreset("tightest")->params;
    CurveSet flat = generateStockinette(pp);
    std::vector<int> counts(pp.courses, pp.wales);
    StitchProgram prog = buildProgram(counts);
    PieceLayout lay;
    for (const StitchRecord& st : prog.stitches) {
        LayoutNode nd;
        nd.p = templateToLocal({0, 0, 0}, st.index, st.turn, pp);
        nd.n = {0, 0, 1};
        nd.c = {1, 0, 0};
        nd.v = {0, 1, 0};
        nd.w = pp.waleSpacing;
        nd.h = pp.courseSpacing;
        lay.nodes.push_back(nd);
    }
    lay.before = lay.nodes.front();
    lay.after = lay.nodes.back();
    PlaneSurface plane(100, 100);
    const CurveSet placed = placeSpiral(plane, lay, prog, pp);
    const std::vector<Vec3>& ctrl = placed.yarns[0].curve.ctrl;
    std::vector<TemplateControl> tpl = stitchTemplate(pp);
    double maxDev = 0;
    for (const StitchRecord& st : prog.stitches)
        for (int j = 0; j < kControlsPerStitch; ++j) {
            const Vec3 ref = templateToLocal(tpl[j], st.index, st.turn, pp);
            maxDev = std::max(maxDev, norm(ctrl[1 + st.id * kControlsPerStitch + j] - ref));  // one padding control first
        }
    CHECK(maxDev < 1e-12);
}

TEST(local_signature_matches_flat) {
    // On the flat patch, the per-stitch check in each stitch's own frame agrees with
    // the global projected signature: every pair PASS, 4 crossings each.
    PatternParams pp = findPreset("tightest")->params;
    CurveSet cs = generateStockinette(pp);
    LocalSignatureResult ls = checkLocalSignature(cs, pp);
    CHECK(ls.status == Status::Pass);
    CHECK(static_cast<int>(ls.pairs.size()) == pp.wales * (pp.courses - 1));
    for (const LocalPairResult& lp : ls.pairs) CHECK(lp.crossings == 4 && lp.unexpected == 0);
    // A broken interlock is caught locally: the whole loop of stitch (wale 2, course
    // 2) moved in front of its parent no longer enters it from behind.
    CurveSet bad = cs;
    YarnPath& yb = bad.yarns[2];
    CHECK(!yb.reversed);
    for (int j = 0; j < kControlsPerStitch; ++j) yb.curve.ctrl[2 * kControlsPerStitch + j + 1].z = 1.5;
    yb.finalize();
    LocalSignatureResult lb = checkLocalSignature(bad, pp);
    CHECK(lb.status == Status::Fail);
    int failing = 0;
    for (const LocalPairResult& lp : lb.pairs) failing += lp.status != Status::Pass;
    CHECK(failing >= 1 && failing <= 3);  // the edited stitch, and at most its child / neighbours
}

struct CylinderPiece {
    PatternParams pp;
    GaugeParams g;
    StitchProgram prog;
    PieceLayout lay;
    CurveSet cs;
};

CylinderPiece cylinderPiece(int n, int turns) {
    CylinderPiece c;
    c.pp = findPreset("tightest")->params;
    c.pp.waleSpacing = c.g.w;
    c.pp.courseSpacing = c.g.h;
    const double R = n * c.g.w / (2 * kPi), H = (turns + 2) * c.g.h;
    ConeSurface cyl(R, R, H);
    c.prog = buildProgram(std::vector<int>(turns, n));
    c.lay = helicalLayout(cyl, c.prog, c.g, evenHelix(c.prog, c.g.h, H));
    c.cs = placeSpiral(cyl, c.lay, c.prog, c.pp);
    finalizeCurves(c.cs);
    return c;
}

TEST(spiral_cylinder_valid) {
    // One helical yarn, 40 stitches per turn, 8 turns: gate (a) on the construction.
    CylinderPiece c = cylinderPiece(40, 8);
    CHECK(c.cs.kind == CurveKind::Surface && c.cs.yarns.size() == 1);
    GraphCheck g = checkProgram(c.prog);
    CHECK(g.status == Status::Pass);
    ValidationReport r = validate(c.cs, c.pp);
    if (r.verdict != Status::Pass) throw std::runtime_error(formatReport(r, c.cs, c.pp, true));
    CHECK(r.localSignature.pass == 40 * 7);
}

TEST(spiral_turn_transition) {
    // Across the shift-path transition: the last stitch of turn r and the first of
    // turn r + 1 are yarn neighbours one stitch apart, the first stitch of turn r + 1
    // is drawn through the first stitch of turn r, and its local signature holds.
    CylinderPiece c = cylinderPiece(24, 5);
    LocalSignatureResult ls = checkLocalSignature(c.cs, c.pp);
    for (int r = 1; r < 5; ++r) {
        const int first = c.prog.turnStart[r], last = first - 1;
        CHECK(c.prog.stitches[first].parents[0] == c.prog.turnStart[r - 1]);
        double d = norm(c.lay.nodes[first].p - c.lay.nodes[last].p);
        CHECK(std::fabs(d - c.g.w) < 0.05 * c.g.w);
        bool found = false;
        for (const LocalPairResult& lp : ls.pairs)
            if (lp.child == first) { found = true; CHECK(lp.status == Status::Pass && lp.crossings == 4); }
        CHECK(found);
    }
}

TEST(same_yarn_interlock_signature) {
    // In a spiral every interlock is between portions of the same yarn: the global
    // check would have rejected each of those crossings; the topology-driven check
    // pairs them, finds 4 crossings with the knit order, and nothing unexpected.
    CylinderPiece c = cylinderPiece(30, 4);
    LocalSignatureResult ls = checkLocalSignature(c.cs, c.pp);
    CHECK(static_cast<int>(ls.pairs.size()) == 30 * 3);
    for (const LocalPairResult& lp : ls.pairs) {
        CHECK(c.cs.topo.stitches[lp.child].yarn == c.cs.topo.stitches[lp.parent].yarn);
        CHECK(lp.status == Status::Pass && lp.crossings == 4 && lp.unexpected == 0);
    }
}

TEST(flat_inc_signature) {
    // The INC construction ("k2 into the same loop") reproduces the signature derived
    // in expected_signature.h, checked by the independent GLOBAL projected search on a
    // flat patch (one yarn per row): B,F,F for slot 0, F,F,B for slot 1, everything
    // else plain, no unexpected crossing.
    ShapingDemo d;
    CHECK(makeShapingDemo("inc", {}, d));
    ValidationReport r = validate(d.curves, d.pp, {.generatorChecks = false});
    if (r.signature.status != Status::Pass) throw std::runtime_error(r.signature.detail);
    int incPairs = 0;
    for (const PairSignature& ps : r.signature.pairs)
        if (d.curves.topo.stitches[ps.child].op == StitchOp::Inc) {
            ++incPairs;
            CHECK(ps.crossingIdx.size() == 3);
        }
    CHECK(incPairs == 2);  // rows 8, 8, 9, 9, 9: one INC, two children
    CHECK(r.signature.unexpected.empty());
}

TEST(flat_dec_signature) {
    // k2tog: the child makes the plain 4 crossings with EACH parent; the two parent
    // loops overlap with the later one in front.
    ShapingDemo d;
    CHECK(makeShapingDemo("dec", {}, d));
    ValidationReport r = validate(d.curves, d.pp, {.generatorChecks = false});
    if (r.signature.status != Status::Pass) throw std::runtime_error(r.signature.detail);
    int decPairs = 0;
    for (const PairSignature& ps : r.signature.pairs)
        if (d.curves.topo.stitches[ps.child].op == StitchOp::Dec) {
            ++decPairs;
            CHECK(ps.crossingIdx.size() == 4);
        }
    CHECK(decPairs == 2);
    CHECK(!r.signature.overlap.empty() && r.signature.overlapWrongOrder == 0);
    // Swapping the layers is detected (the earlier parent in front).
    ShapingParams sw;
    sw.decDepth = -sw.decDepth;
    ShapingDemo e;
    makeShapingDemo("dec", sw, e);
    ValidationReport re = validate(e.curves, e.pp, {.generatorChecks = false});
    CHECK(re.signature.status == Status::Fail);
}

TEST(spiral_cone_shaping_signature) {
    // Cylinder -> cone -> cylinder on one helical yarn with INC and DEC turns: the
    // program balances and every interlock (K, INC, DEC, same yarn) has its expected
    // signature in its own frame. Curvature and clearance are NOT asserted here: the
    // construction leaves tight bends next to the shaping, which is the repair's job.
    ShapingDemo d;
    CHECK(makeShapingDemo("cone", {}, d));
    GraphCheck g = checkProgram(d.prog);
    if (g.status != Status::Pass) throw std::runtime_error(g.detail);
    CHECK(g.inc > 0 && g.dec > 0);
    LocalSignatureResult ls = checkLocalSignature(d.curves, d.pp);
    if (ls.status != Status::Pass) throw std::runtime_error(ls.detail);
    int dec = 0, inc = 0;
    for (const LocalPairResult& lp : ls.pairs) {
        const StitchOp op = d.curves.topo.stitches[lp.child].op;
        if (op == StitchOp::Dec) { ++dec; CHECK(lp.crossings == 8); }
        if (op == StitchOp::Inc) { ++inc; CHECK(lp.crossings == 3); }
    }
    CHECK(dec == g.dec && inc == 2 * g.inc);
}

// Crossing attribution (checks.h attributeCrossing): one yarn, a parent row P0..P3 on [0, 64), a gap
// [64, 84) without stitches, then Z [84, 100) (no parent in the row), C0 [100, 116) child of P0,
// C1 [116, 132) child of P1, C2 [132, 148) child of P2, G [148, 150) (short, no parent), C3
// [150, 166) child of P2.
KnitTopology attributionTopology() {
    KnitTopology t;
    auto add = [&](double u0, double u1, std::vector<int> parents) {
        StitchRecord s;
        s.id = static_cast<int>(t.stitches.size());
        s.uBegin = u0;
        s.uEnd = u1;
        for (int p : parents) s.parents[s.nParents++] = p;
        t.stitches.push_back(s);
    };
    for (int k = 0; k < 4; ++k) add(16.0 * k, 16.0 * k + 16, {});  // P0..P3: ids 0..3
    add(84, 100, {});    // 4: Z
    add(100, 116, {0});  // 5: C0
    add(116, 132, {1});  // 6: C1
    add(132, 148, {2});  // 7: C2
    add(148, 150, {});   // 8: G
    add(150, 166, {2});  // 9: C3
    linkChildren(t.stitches);
    t.buildIndex({200.0});
    return t;
}

TEST(attribution_crossing_rules) {
    const KnitTopology t = attributionTopology();
    const double w = 2;
    auto owner = [&](double u, int parent, double win = 2) { return attributeCrossing(t, 0, u, parent, win); };
    CHECK(owner(110, 0).owner == 5);                              // the own span (original rule)
    CHECK(owner(117, 0).owner == 5);                              // in C1's span, C1 not a child of P0: C0 within the window
    CHECK(owner(117.999, 0).owner == 5 && owner(118, 0).owner == -1);  // upper window boundary open: [uBegin - w, uEnd + w)
    CHECK(owner(98, 0).owner == 5 && owner(97.999, 0).owner == -1);    // lower window boundary closed
    CHECK(owner(117, 1).owner == 6);                              // own span first even inside the neighbour's window
    const CrossingOwner none = owner(110, 1);                     // in C0's span, parent P1: C1's window starts at 114
    CHECK(none.owner == -1 && !none.ambiguous);                   //   no stitch can take it (unexpected for C0)
    const CrossingOwner amb = owner(149, 2);                      // in G's span: C2 (to 150) and C3 (from 148) both reach it
    CHECK(amb.ambiguous && amb.owner == -1);
    CHECK(owner(70, 0).owner == -1 && owner(83, 0).owner == -1);  // outside every stitch span: not attributed
    CHECK(owner(117, 0, 0.0).owner == -1);                        // window 0: the original rule
    (void)w;
}

TEST(attribution_close_rules) {
    // CLOSE drawstring: loops gathering heads 20, 21, 22 on strand ranges [0, 5], [5, 10], [10, 15].
    TerminalRecord tr;
    tr.op = TerminalOp::Close;
    tr.stitches = {20, 21, 22};
    tr.segments = {{0, 5}, {5, 10}, {10, 15}};
    CHECK(attributeCloseCrossing(tr, 21, 7, 2).owner == 1);       // own range
    CHECK(attributeCloseCrossing(tr, 21, 11, 2).owner == 1);      // in loop 2's range, head 21: loop 1 within the window
    CHECK(attributeCloseCrossing(tr, 21, 12, 2).owner == 1);      // window boundary closed: [first - w, second + w]
    CHECK(attributeCloseCrossing(tr, 21, 12.001, 2).owner == -1);
    CHECK(attributeCloseCrossing(tr, 22, 15.5, 2).owner == -1);   // beyond the last range (strand tail): not attributed
    CHECK(attributeCloseCrossing(tr, 20, -0.5, 2).owner == -1);   // before the first range: not attributed
    CHECK(attributeCloseCrossing(tr, 21, 11, 0).owner == -1);     // window 0: the original rule
    // one head gathered twice (constructed): a crossing within reach of both loops is ambiguous
    TerminalRecord twice = tr;
    twice.stitches = {20, 21, 20};
    const CrossingOwner amb = attributeCloseCrossing(twice, 20, 7.5, 3);
    CHECK(amb.ambiguous && amb.owner == -1);
}

TEST(local_signature_attribution_moved_boundary) {
    // A valid spiral; the ownership boundary between two yarn neighbours moved by s spans (the
    // geometry unchanged). The template keeps its outer crossings ~0.86 span inside the boundary:
    // s = 1 moves one crossing of the later stitch into the earlier one's span. The original rule
    // then fails both (3 crossings; 1 unexpected), the attribution window (2 spans) gives it back;
    // s = 4 moves it 3.1 spans away: beyond the window, both rules fail. Nothing else changes.
    CylinderPiece c = cylinderPiece(24, 5);
    const LocalSignatureResult ref0 = checkLocalSignature(c.cs, c.pp, {}, 0.0), ref = checkLocalSignature(c.cs, c.pp);
    CHECK(ref0.status == Status::Pass && ref.status == Status::Pass);
    const int a = c.prog.turnStart[2] + 5;
    auto moved = [&](double s) {
        CurveSet m = c.cs;
        m.topo.stitches[a].uEnd += s;
        m.topo.stitches[a + 1].uBegin += s;
        m.topo.buildIndex({m.yarns[0].uEnd()});
        return m;
    };
    auto failing = [&](const LocalSignatureResult& r) {
        std::vector<int> f;
        for (const LocalPairResult& lp : r.pairs)
            if (lp.status != Status::Pass) f.push_back(lp.child);
        return f;
    };
    for (double s : {1.0, -1.0}) {
        const CurveSet m = moved(s);
        const LocalSignatureResult r0 = checkLocalSignature(m, c.pp, {}, 0.0), r = checkLocalSignature(m, c.pp);
        std::printf("    boundary moved by %+g spans: original rule %s (%zu failing), attribution window %s (%d re-evaluated)\n", s, toString(r0.status),
                    failing(r0).size(), toString(r.status), r.attributionReevaluated);
        CHECK(r0.status == Status::Fail);
        for (int ch : failing(r0)) CHECK(ch == a || ch == a + 1);
        CHECK(r.status == Status::Pass);
    }
    {
        const CurveSet m = moved(4.0);
        const LocalSignatureResult r0 = checkLocalSignature(m, c.pp, {}, 0.0), r = checkLocalSignature(m, c.pp);
        std::printf("    boundary moved by +4 spans: original rule %s, attribution window %s\n", toString(r0.status), toString(r.status));
        CHECK(r0.status == Status::Fail && r.status == Status::Fail);
        for (int ch : failing(r)) CHECK(ch == a || ch == a + 1);
    }
}

// ------------------------------------------------------------ S1: chart ----

#ifndef KNIT_TEST_DATA
#define KNIT_TEST_DATA "tests/data"
#endif

Chart fixtureChart() {
    TriMesh m;
    std::string err;
    if (!readObj(std::string(KNIT_TEST_DATA) + "/chart_fixture.obj", m, &err)) throw std::runtime_error(err);
    return Chart(m, 5);  // the fixture sphere has 5 rings
}

TEST(chart_seam_continuity) {
    // The fixture went through the same export as the scene pieces (face-
    // corner chart UVs, Catmull-Clark level 2): both sides of the seam must agree.
    Chart ch = fixtureChart();
    ChartCheck c = ch.check();
    if (c.seam != Status::Pass) throw std::runtime_error(c.detail);
    CHECK(c.seamSamples > 300 && c.domainGaps == 0);
    CHECK(c.seamVertexMismatches == 0);
    // Seam corners really are unwrapped (b reaches 1 on the far side), so the test
    // is not vacuous.
    double bmax = 0;
    for (const auto& u : ch.mesh().uv)
        for (const Vec2& q : u) bmax = std::max(bmax, q.y);
    CHECK(bmax == 1.0);
    // eval wraps b: just below 1 and just above 0 are neighbours on the surface.
    for (double a : {0.2, 0.5, 0.8}) {
        ChartPoint s0 = ch.eval(a, 1e-7), s1 = ch.eval(a, 1.0 - 1e-7);
        CHECK(s0.ok && s1.ok && s0.tri != s1.tri);
        CHECK(norm(s0.p - s1.p) < 1e-5);
    }
}

TEST(chart_orientation) {
    Chart ch = fixtureChart();
    ChartCheck c = ch.check();
    if (c.orientation != Status::Pass || c.monotoneA != Status::Pass) throw std::runtime_error(c.detail);
    CHECK(c.aPoles == 2);
    // eval -> closestPoint round trip recovers the point and its chart coordinates.
    for (int i = 0; i < 50; ++i) {
        double a = ch.aMin() + (ch.aMax() - ch.aMin()) * ((i * 37) % 50) / 49.0, b = ((i * 13) % 50) / 50.0;
        ChartPoint s = ch.eval(a, b);
        CHECK(s.ok);
        ChartPoint q = ch.closestPoint(s.p);
        CHECK(q.ok && q.distance < 1e-12);
        CHECK_NEAR(q.ab.x, a, 1e-9);
        double db = std::fabs(q.ab.y - b);
        CHECK(std::min(db, std::fabs(db - 1)) < 1e-9);
    }
    // Ring length is continuous in a and vanishes at the poles.
    CHECK(ch.ringLength(1e-6) < 1e-3 * ch.ringLength(0.5));
    CHECK(ch.ringLength(0.5) > 0);
}

TEST(capped_sphere_poles) {
    // Spiral from pole A (CAST ring of 6) to pole B (threaded drawstring, virtual row
    // endGap from the pole): the program balances, terminal records carry the right
    // stitches, the pieces's c x v points outwards everywhere, and the hole at each
    // pole is measured on the ACTUAL curve (certified lower bound), positive and
    // below the cast-ring / end-gap radius.
    FixturePiece f = fixturePiece(25.0);
    PatternParams pp;
    PieceKnitParams kp;
    PieceKnit pk;
    std::string err;
    if (!knitPiece({*f.surf, true}, {pp, kp}, knit_test::cpuNumerics(), pk, &err)) throw std::runtime_error(err);
    GraphCheck g = checkProgram(pk.prog);
    if (g.status != Status::Pass) throw std::runtime_error(g.detail);
    CHECK(pk.prog.castRing && pk.prog.turns.front().output == kp.castCount);
    CHECK(pk.curves.topo.terminals.size() == 2);
    CHECK(pk.curves.topo.terminals[0].op == TerminalOp::Cast &&
          static_cast<int>(pk.curves.topo.terminals[0].stitches.size()) == kp.castCount);
    CHECK(static_cast<int>(pk.curves.topo.terminals[1].stitches.size()) == pk.prog.turns.back().output);
    for (const LayoutNode& nd : pk.lay.nodes) CHECK(dot(cross(nd.c, nd.v), nd.n) > 0);
    const double holeA = poleHoleRadius(pk.curves, f.surf->poleA(), pp.radius());
    const double holeB = poleHoleRadius(pk.curves, f.surf->poleB(), pp.radius());
    CHECK(holeA > 0 && holeA < pk.startArc + kp.g.h);  // the ring's lower parts sit at the head base
    CHECK(holeB > 0 && holeB < kp.endGap + kp.g.h);
    // The certified bound is below every sampled distance.
    double sampled = 1e300;
    const YarnPath& y = pk.curves.yarns[0];
    for (double u = 0; u <= y.uEnd(); u += 0.05) sampled = std::min(sampled, norm(y.curve.eval(u) - f.surf->poleA()));
    CHECK(holeA + pp.radius() <= sampled + 1e-9);
}

// A gathered end mirrors the start, as a knitter closes a piece: the last turn decreases to the cast count, the turns fill the
// meridian with nothing left over, the drawstring's ring sits n (2 D + c) / 2 pi round the end pole, and the opening
// seen from outside at the end pole is no larger than at the start pole. (The retired construction left the virtual
// row 8.6 to 15 mm off the end pole, the remainder of the rows deciding where, and an opening of 6 to 10 mm.)
TEST(spiral_end_gathered) {
    FixturePiece f = fixturePiece(25.0);
    PatternParams pp;
    pp.yarnDiameter = 1.16;  // the Dense preset's knitted yarn
    pp.clearance = 0.1;
    PieceKnitParams kp;
    kp.g.w = 5.5;
    kp.g.h = 3.4;
    kp.end = EndClosure::Gathered;
    PieceKnit pk;
    std::string err;
    CHECK(knitPiece({*f.surf, true}, {pp, kp}, knit_test::cpuNumerics(), pk, &err));
    if (!err.empty()) std::printf("    %s\n", err.c_str());
    const Vec3 endPole = f.surf->poleB(), startPole = f.surf->poleA();
    double ringLow = 1e300, ringHigh = 0;
    for (const LayoutNode& v : pk.lay.virtualNodes) {
        ringLow = std::min(ringLow, norm(v.p - endPole));
        ringHigh = std::max(ringHigh, norm(v.p - endPole));
    }
    const double openStart = poleOpening(pk.curves, startPole, f.surf->project(startPole).n, pp.yarnDiameter);
    const double openEnd = poleOpening(pk.curves, endPole, f.surf->project(endPole).n, pp.yarnDiameter);
    std::printf("    %s\n    last turn %d loops; ring %.2f mm, virtual children %.2f..%.2f mm off the end pole; opening start %.2f, end %.2f mm\n",
                pk.detail.c_str(), pk.prog.turns.back().output, pk.closeRing, ringLow, ringHigh, openStart, openEnd);
    CHECK(pk.prog.turns.back().output == kp.castCount);
    CHECK(std::fabs(pk.closeRing - kp.castCount * (2 * pp.yarnDiameter + pp.clearance) / (2 * kPi)) < 1e-12);
    CHECK(ringLow > 0.8 * pk.closeRing && ringHigh < 1.25 * pk.closeRing);  // chords on the cap, near the ring's arc radius
    CHECK(std::fabs(pk.rowStep - kp.g.h) <= 0.5 * kp.g.h / std::round(f.surf->totalArc() / kp.g.h));
    CHECK(openEnd <= openStart);
    const ValidationReport r = validate(pk.curves, pp, ValidateOptions{.generatorChecks = true});
    std::printf("    gate a before settle: %s, local signature %d/%zu\n", toString(r.verdict), r.localSignature.pass, r.localSignature.pairs.size());
}

// A scene piece gathers its end only where the end pole shows: alone, its end is gathered; with another piece of
// its object covering the end pole, the end stays threaded (sewn on there, never seen).
TEST(scene_piece_gathers_only_a_visible_end) {
    FixturePiece f = fixturePiece(25.0);
    PatternParams pp;
    PieceKnitParams kp;
    TriMesh cover = icosphere(3, 12, 1);
    for (Vec3& p : cover.pos) p = p + f.chart->poleB();
    const PlacedPiece piece{&f.mesh, {0, 0, 0}, 0};
    const std::vector<PlacedPiece> alone = {piece};
    const std::vector<PlacedPiece> covered = {piece, {&cover, {0, 0, 0}, 0}};
    ScenePieceKnit open, hidden;
    std::string err;
    CHECK(knitScenePiece({piece, 5, alone}, {pp, kp}, knit_test::cpuNumerics(), open, &err));
    CHECK(knitScenePiece({piece, 5, covered}, {pp, kp}, knit_test::cpuNumerics(), hidden, &err));
    std::printf("    alone: %s\n    end pole covered: %s\n", open.knit.detail.c_str(), hidden.knit.detail.c_str());
    CHECK(!open.poleBInside && open.knit.closeRing > 0);
    CHECK(hidden.poleBInside && hidden.knit.startAtA && hidden.knit.closeRing == 0);
}

// A spiral piece's placement and hand before settle, bit for bit (pinned at 05960b1): the helix's nodes and virtual
// row, the placed and hand-moved controls, the hand's records and the yarn's thickness, for a gathered end (alone)
// and a threaded end (the end pole under a cover).
TEST(spiral_placement_pinned) {
    FixturePiece f = fixturePiece(25.0);
    const TriMesh cover = octahedronAround(f.chart->poleB(), 3.0);
    const PlacedPiece piece{&f.mesh, {0, 0, 0}, 0};
    const std::vector<PlacedPiece> alone = {piece}, covered = {piece, {&cover, {0, 0, 0}, 0}};
    PieceKnitParams kp;
    kp.hand = everyChannelHand();
    const uint64_t pinned[2] = {0xd59f680eea3f27b1ull, 0x36db43085c1e1fcdull};
    uint64_t digest[2] = {0, 0};
    int k = 0;
    for (const std::vector<PlacedPiece>* others : {&alone, &covered}) {
        ScenePieceKnit spk;
        std::string err;
        CHECK(knitScenePiece({piece, 5, *others}, {PatternParams{}, kp}, knit_test::cpuNumerics(), spk, &err));
        const std::vector<StitchTension>& tensions = spk.hand.tensions;
        int raised = 0, spots = 0;
        for (const StitchTension& t : tensions) {
            raised += t.raise > 0;
            spots += t.spot != 0;
        }
        digest[k] = fnv1a64(placementText(spk.knit, tensions));
        std::printf("    %s end: %zu stitches, %d raised loops, %d spots, close ring %.3f mm, digest %016llx\n", k == 0 ? "gathered" : "threaded",
                    spk.knit.prog.stitches.size(), raised, spots, spk.knit.closeRing, static_cast<unsigned long long>(digest[k]));
        CHECK(raised > 0 && spots > 0);  // the hand is really on
        ++k;
    }
    CHECK(digest[0] == pinned[0] && digest[1] == pinned[1]);
}

TEST(layout_distortion_bounds) {
    // The layout check flags cells outside the width/height bounds and inverted cells.
    FixturePiece f = fixturePiece(25.0);
    PatternParams pp;
    PieceKnitParams kp;
    PieceKnit pk;
    CHECK(knitPiece({*f.surf, true}, {pp, kp}, knit_test::cpuNumerics(), pk));
    LayoutCheck lc = checkLayout(pk.prog, pk.lay, kp.g, 0.0, 1e9);
    CHECK(lc.outOfBounds == 0 && lc.inverted == 0 && lc.cells == static_cast<int>(pk.prog.stitches.size()));
    CHECK(lc.minW > 0.5 && lc.maxW < 2.0);
    // Tight bounds around 1 must flag cells; a flipped frame must count as inverted.
    LayoutCheck tight = checkLayout(pk.prog, pk.lay, kp.g, 0.999, 1.001);
    CHECK(tight.outOfBounds > 0 && tight.status == Status::Fail);
    PieceLayout flipped = pk.lay;
    for (LayoutNode& nd : flipped.nodes) nd.c = -nd.c;
    CHECK(checkLayout(pk.prog, flipped, kp.g).inverted > 0);
}

TEST(settle_keeps_topology_on_sphere) {
    // Settle reduces construction defects on the capped sphere without breaking any
    // interlock that was correct; the result is judged by gate (a) again (settle
    // itself is not a certificate).
    FixturePiece f = fixturePiece(25.0);
    PatternParams pp;
    PieceKnitParams kp;
    pp.waleSpacing = kp.g.w;
    pp.courseSpacing = kp.g.h;
    PieceKnit pk;
    CHECK(knitPiece({*f.surf, true}, {pp, kp}, knit_test::cpuNumerics(), pk));
    ValidationReport before = validate(pk.curves, pp);
    ValidationReport after;
    SettleReport sr = settle({pk.curves, pp, f.surf.get()}, SettleParams{}, knit_test::cpuNumerics(), &after);
    // Never worse in topology: every stitch failing after settling already failed.
    std::set<int> failBefore;
    for (const LocalPairResult& lp : before.localSignature.pairs)
        if (lp.status != Status::Pass) failBefore.insert(lp.child);
    for (const LocalPairResult& lp : after.localSignature.pairs)
        if (lp.status != Status::Pass) CHECK(failBefore.count(lp.child) == 1);
    int fb = 0, fa = 0;
    for (const ClearanceRegion& g : before.clearance.regions) fb += g.status == Status::Fail;
    for (const ClearanceRegion& g : after.clearance.regions) fa += g.status == Status::Fail;
    CHECK(fa <= fb);
    CHECK(sr.defectsPerRound.empty() || sr.defectsPerRound.back() <= sr.defectsPerRound.front());
}

TEST(settle_is_one_relax_round) {
    // Settle relaxes once: with defects left after the relaxation (a no-op relax leaves all of them)
    // it neither repeats nor restores anything; one relax call, one round, the starting defect count.
    struct NoMove final : knit_test::ForwardingNumerics {
        NoMove() : ForwardingNumerics(knit_test::cpuNumerics()) {}
        void relax(const RelaxRequest&) override { ++calls; }
        int calls = 0;
    };
    FixturePiece f = fixturePiece(25.0);
    PatternParams pp;
    PieceKnitParams kp;
    pp.waleSpacing = kp.g.w;
    pp.courseSpacing = kp.g.h;
    pp.yarnDiameter = 1.5;  // thicker than the gauge holds: defects before and after
    PieceKnit pk;
    CHECK(knitPiece({*f.surf, true}, {pp, kp}, knit_test::cpuNumerics(), pk));
    NoMove nm;
    ValidationReport last;
    const SettleReport sr = settle({pk.curves, pp, f.surf.get()}, SettleParams{}, nm, &last);
    std::printf("    %d round(s), %zu defect counts, %d at the start; verdict after %s\n", sr.rounds, sr.defectsPerRound.size(),
                sr.defectsPerRound.empty() ? -1 : sr.defectsPerRound[0], toString(last.verdict));
    CHECK(sr.defectsPerRound.size() == 1 && sr.defectsPerRound[0] > 0);
    CHECK(sr.rounds == 1 && nm.calls == 1);
    CHECK(last.verdict != Status::Pass);
}

TEST(stitch_spread_measures_one_stray_control) {
    // A flat patch's stitches share one shape. One control moved along N by d among n stitches of 16
    // controls: its stitch is off by d (1 - 1/n), each other one by d / n at that control.
    PatternParams pp = testParams();
    pp.wales = 6;
    pp.courses = 7;
    CurveSet cs = generateStockinette(pp);
    const StitchSpread even = measureStitchSpread(cs);
    CHECK(even.stitches == 12 && even.shape < 1e-9 && even.height < 1e-9);  // rows 2 and 4: the rest are edges or run backwards
    const double d = 0.16, n = even.stitches;
    const StitchRecord& st = cs.topo.stitches[2 * pp.wales + 3];
    cs.yarns[st.yarn].curve.ctrl[static_cast<int>(st.uBegin) + 1 + 5].z += d;
    const StitchSpread stray = measureStitchSpread(cs);
    const double expected = d * std::sqrt((1 - 1 / n) / (n * kControlsPerStitch));
    std::printf("    height %.6f mm, expected %.6f mm\n", stray.height, expected);
    CHECK(std::abs(stray.height - expected) < 1e-9 && std::abs(stray.shape - expected) < 1e-9);
}

TEST(hand_tensions_are_pure) {
    // Two draws agree, and a stitch keeps its record when another size of fabric gives it another id.
    PatternParams pp = testParams();
    pp.wales = 6;
    pp.courses = 8;
    const CurveSet narrow = generateStockinette(pp);
    pp.wales = 7;
    const CurveSet wide = generateStockinette(pp);
    HandRequest hand;
    hand.settings.unevenness = 30;
    hand.seed = 7;
    const std::vector<StitchTension> a = drawStitchTensions(narrow.topo, hand), b = drawStitchTensions(narrow.topo, hand);
    const std::vector<StitchTension> c = drawStitchTensions(wide.topo, hand);
    const int turn = 3, index = 2;
    const StitchTension& inNarrow = a[turn * 6 + index];
    const StitchTension& inWide = c[turn * 7 + index];
    CHECK(inNarrow.tension != 0 && inNarrow.tension == inWide.tension && inNarrow.lean == inWide.lean);
    for (size_t t = 0; t < a.size(); ++t) CHECK(a[t].tension == b[t].tension && a[t].lean == b[t].lean);
    CHECK(a[0].tension == 0 && a[turn * 6].tension == 0);  // the first turn and a run's end stay as they are
}

TEST(hand_rate_per_cm2) {
    // A rate per cm2 is a chance per stitch: the count over a large fabric lies within 3 sigma of it.
    PatternParams pp = testParams();
    pp.wales = 60;
    pp.courses = 60;
    const CurveSet cs = generateStockinette(pp);
    HandRequest hand;
    hand.settings.loopRate = 1;
    hand.settings.loopHeight = 1;
    hand.seed = 3;
    const double chance = hand.settings.loopRate * hand.stitchArea / 100;
    int eligible = 0, raised = 0;
    const std::vector<StitchTension> records = drawStitchTensions(cs.topo, hand);
    for (const StitchRecord& st : cs.topo.stitches) {
        eligible += st.turn > 0 && st.turn + 1 < pp.courses && st.index > 0 && st.index + 1 < pp.wales;
        raised += records[st.id].raise != 0;
    }
    const double expected = eligible * chance, sigma = std::sqrt(eligible * chance * (1 - chance));
    std::printf("    %d raised of %d, expected %.0f +- %.0f\n", raised, eligible, expected, sigma);
    CHECK(std::abs(raised - expected) < 3 * sigma);
}

TEST(hand_raised_loop_is_loose) {
    // A raised loop holds spare yarn whatever its own draw says: a tight one is pulled back down.
    PatternParams pp = testParams();
    pp.wales = 30;
    pp.courses = 30;
    const CurveSet cs = generateStockinette(pp);
    HandRequest hand;
    hand.settings.unevenness = 30;
    hand.settings.loopRate = 2;
    hand.settings.loopHeight = 1;
    hand.seed = 11;
    int raised = 0;
    for (const StitchTension& t : drawStitchTensions(cs.topo, hand)) {
        if (t.raise == 0) continue;
        ++raised;
        CHECK(t.tension <= 0 && t.raise >= 0.5f && t.raise <= 1.0f);
    }
    CHECK(raised > 0);
}

TEST(hand_thickness_follows_tension_and_spots) {
    // A spot reaches the largest change at its middle and fades along the yarn; a tight stitch is thin;
    // the clamp holds; plies of one radius get none.
    PatternParams pp = testParams();
    pp.wales = 8;
    pp.courses = 5;
    CurveSet cs = generateStockinette(pp);
    std::vector<StitchTension> tensions(cs.topo.stitches.size());
    const StitchRecord& st = cs.topo.stitches[2 * pp.wales + 3];
    tensions[st.id].spot = 0.4f;
    HandRequest hand;
    hand.settings.spotLength = 10;
    hand.maxThickness = 1.3;
    setYarnThickness(cs, tensions, hand);
    const std::vector<float>& t = cs.yarns[st.yarn].thickness;
    const int middle = static_cast<int>(st.uBegin) + 1 + kControlsPerStitch / 2;
    CHECK(t.size() == cs.yarns[st.yarn].curve.ctrl.size());
    CHECK(std::abs(t[middle] - 1.3f) < 1e-6);                        // 1.4 clamped
    CHECK(t[middle + 2] > 1.0f && t[middle + 2] <= t[middle]);
    CHECK(std::abs(t[middle + 3 * kControlsPerStitch] - 1.0f) < 1e-3);  // far away
    CHECK(std::abs(cs.yarns[st.yarn].thicknessAt(middle - 1) - t[middle]) < 1e-6);  // control c sits at the knot u = c - 1
    tensions[st.id] = {};
    tensions[st.id].tension = 0.2f;
    setYarnThickness(cs, tensions, hand);
    CHECK(std::abs(cs.yarns[st.yarn].thickness[middle] - 0.9f) < 1e-6);
    CurveSet plied = generateStockinette(pp);
    hand.maxThickness = 0;
    setYarnThickness(plied, tensions, hand);
    CHECK(plied.yarns[st.yarn].thickness.empty());
}

TEST(hand_rest_has_no_jump) {
    // Every edge rests at its own length times its stitch's scale, the edges between two stitches too:
    // scaling a stitch about its centre would stretch those.
    PatternParams pp = testParams();
    pp.wales = 6;
    pp.courses = 4;
    const CurveSet cs = generateStockinette(pp);
    const Controls P = controlsOf(cs);
    std::vector<StitchTension> tensions(cs.topo.stitches.size());
    const StitchRecord& st = cs.topo.stitches[2 * pp.wales + 3];
    tensions[st.id].tension = 0.2f;
    const Controls rest = restOfTensions(P, cs.topo, tensions);
    const int first = static_cast<int>(st.uBegin) + 1;
    for (size_t y = 0; y < P.size(); ++y)
        for (int c = 0; c + 1 < static_cast<int>(P[y].size()); ++c) {
            const bool tight = static_cast<int>(y) == st.yarn && c >= first && c < first + kControlsPerStitch;
            CHECK(std::abs(norm(rest[y][c + 1] - rest[y][c]) - (tight ? 0.8f : 1.0f) * norm(P[y][c + 1] - P[y][c])) < 1e-6);
        }
}

TEST(hand_rest_counts_the_tension_once) {
    // A stitch of tension t asks for 1 - t of the yarn it held before the hand. knitByHand's shrink makes
    // it smaller already (fading to nothing at the sinkers); scaling that smaller stitch again would ask
    // for about (1 - t)^1.5 (measured 2026-11-15, off by 0.46 t (1 - t)).
    PatternParams pp = testParams();
    pp.wales = 12;
    pp.courses = 10;
    CurveSet cs = generateStockinette(pp);
    for (StitchRecord& st : cs.topo.stitches) st.lattice = {0, st.index, st.turn, st.turn};
    const Controls before = controlsOf(cs);
    HandRequest hand;
    hand.settings.unevenness = 30;
    hand.seed = 7;
    hand.wale = pp.waleSpacing;
    const HandResult r = knitByHand(cs, hand, knit_test::cpuNumerics(), {});
    const auto lengthOf = [](const std::vector<Vec3>& P, int first) {
        double sum = 0;
        for (int c = first; c < first + kControlsPerStitch; ++c) sum += norm(P[c + 1] - P[c]);
        return sum;
    };
    int strong = 0;
    double worst = 0;
    for (const StitchRecord& st : cs.topo.stitches) {
        const double t = r.tensions[st.id].tension;
        if (t == 0) continue;
        const int first = static_cast<int>(st.uBegin) + 1;
        const double asked = lengthOf(r.rest[st.yarn], first) / lengthOf(before[st.yarn], first);
        worst = std::max(worst, std::abs(asked - (1 - t)));
        strong += std::abs(t) > 0.1;
    }
    std::printf("    %d stitches with |t| > 0.1; worst |rest / before - (1 - t)| %.4f\n", strong, worst);
    CHECK(strong > 5 && worst < 0.01);
}

TEST(settle_hand_frees_the_whole_fabric) {
    // Settle frees only defects, and a hand frees every stitch.
    struct NoMove final : knit_test::ForwardingNumerics {
        NoMove() : ForwardingNumerics(knit_test::cpuNumerics()) {}
        void relax(const RelaxRequest&) override { ++calls; }
        int calls = 0;
    };
    FixturePiece f = fixturePiece(25.0);
    PatternParams pp;
    PieceKnitParams kp;
    pp.waleSpacing = kp.g.w;
    pp.courseSpacing = kp.g.h;
    pp.yarnDiameter = 0.8;
    PieceKnit pk;
    CHECK(knitPiece({*f.surf, true}, {pp, kp}, knit_test::cpuNumerics(), pk));
    const Controls rest = controlsOf(pk.curves);
    const int stitches = static_cast<int>(pk.curves.topo.stitches.size());
    PieceKnit plain = pk;
    NoMove idle, handed;
    SettleParams sp;
    const SettleReport none = settle({plain.curves, pp, f.surf.get()}, sp, idle);
    const SettleReport sr = settle({pk.curves, pp, f.surf.get(), nullptr, &rest}, sp, handed);
    std::printf("    defects %d, with a hand %d of %d stitches\n", none.defectsPerRound[0], sr.defectsPerRound[0], stitches);
    CHECK(none.defectsPerRound[0] < stitches && sr.defectsPerRound[0] == stitches && handed.calls == 1);
}

TEST(invalid_init_rejected) {
    // A mis-threaded construction (one loop of the spiral moved entirely in front of
    // its parent) fails gate (a) at every diameter: the repair finds no verified
    // start diameter and leaves the curves untouched.
    CylinderPiece c = cylinderPiece(24, 5);
    CurveSet bad = c.cs;
    const StitchRecord& st = bad.topo.stitches[2 * 24 + 7];
    for (int j = 0; j < kControlsPerStitch; ++j) {
        Vec3& q = bad.yarns[0].curve.ctrl[static_cast<size_t>(st.uBegin) + j + 1];
        q = q + st.N * 2.0;
    }
    bad.yarns[0].finalize();
    ValidationReport r = validate(bad, c.pp);
    CHECK(r.localSignature.status == Status::Fail);
    const std::vector<Vec3> before = bad.yarns[0].curve.ctrl;
    RepairReport rr = certifiedRepair({bad, c.pp, nullptr}, RepairParams{}, knit_test::cpuNumerics());
    CHECK(!rr.reachedTarget && rr.D0 == 0 && !rr.moveTried);
    CHECK(rr.detail.find("no verified start diameter") != std::string::npos);
    CHECK(bad.yarns[0].curve.ctrl == before);
}

TEST(repair_with_failing_signature_validates_nothing) {
    // A local signature FAIL fails at every diameter: given the current report,
    // repair ends at once with no start diameter, validates nothing and leaves the curves untouched.
    struct CountingValidate final : knit_test::ForwardingNumerics {
        CountingValidate() : ForwardingNumerics(knit_test::cpuNumerics()) {}
        ValidationReport validate(const CurveSet& cs, const PatternParams& pp, const ValidateOptions& opt) override {
            ++calls;
            return inner_.validate(cs, pp, opt);
        }
        int calls = 0;
    };
    CylinderPiece c = cylinderPiece(24, 5);
    CurveSet bad = c.cs;
    const StitchRecord& st = bad.topo.stitches[2 * 24 + 7];
    for (int j = 0; j < kControlsPerStitch; ++j) {
        Vec3& q = bad.yarns[0].curve.ctrl[static_cast<size_t>(st.uBegin) + j + 1];
        q = q + st.N * 2.0;
    }
    bad.yarns[0].finalize();
    const ValidationReport current = validate(bad, c.pp);
    CHECK(current.localSignature.fail > 0);
    const std::vector<Vec3> before = bad.yarns[0].curve.ctrl;
    CountingValidate counting;
    ValidationReport fin;
    const RepairReport rr = certifiedRepair({bad, c.pp, nullptr}, RepairParams{}, counting, &fin, &current);
    std::printf("    %s; %d validations\n", rr.detail.c_str(), counting.calls);
    CHECK(counting.calls == 0);
    CHECK(!rr.reachedTarget && !rr.moveTried);
    CHECK(rr.detail.find("no verified start diameter") != std::string::npos);
    CHECK(fin.localSignature.fail == current.localSignature.fail);
    CHECK(bad.yarns[0].curve.ctrl == before);
}

TEST(diameter_monotone_fixed_centreline) {
    // For a fixed centreline the set of diameters at which the curve set is valid is
    // an initial segment: on a grid, once a diameter fails, every larger one fails.
    CylinderPiece c = cylinderPiece(24, 5);
    bool failed = false;
    int passes = 0;
    for (double D = 0.3; D <= 1.6 + 1e-9; D += 0.1) {
        PatternParams q = c.pp;
        q.yarnDiameter = D;
        const bool ok = validate(c.cs, q).verdict == Status::Pass;
        if (failed) CHECK(!ok);
        if (!ok) failed = true;
        passes += ok;
    }
    CHECK(passes >= 7 && failed);  // valid at least up to D = 0.9 here, invalid by 1.6
}

TEST(repair_one_certified_move) {
    // Repair of a construction that is topologically right but too tight for the target: the cylinder
    // spiral at D = 1.3 mm (valid at smaller D). The start diameter read off the report validates; the
    // one move is kept only if it passes at the target and its straight motion is certified at D0
    // (no pass-through, delta < D0 + c); otherwise the curves are exactly the input's.
    CylinderPiece c = cylinderPiece(24, 5);
    PatternParams pp = c.pp;
    pp.yarnDiameter = 1.3;
    const ValidationReport atTarget = validate(c.cs, pp);
    CHECK(atTarget.verdict != Status::Pass);
    const std::vector<Vec3> before = c.cs.yarns[0].curve.ctrl;
    ValidationReport fin;
    RepairReport rr = certifiedRepair({c.cs, pp, nullptr}, RepairParams{}, knit_test::cpuNumerics(), &fin, &atTarget);
    std::printf("    %s\n", rr.detail.c_str());
    if (!(rr.D0 > 0)) throw std::runtime_error(rr.detail);
    CHECK(rr.D0 < pp.yarnDiameter);
    PatternParams atD0 = pp;
    atD0.yarnDiameter = rr.D0;
    CHECK(rr.moveTried);
    CHECK(rr.validations == 2);  // given the report at the target: D0 once, the moved state once, no search
    CHECK(rr.reachedTarget == (fin.verdict == Status::Pass));
    CHECK(rr.reachedTarget == rr.move.accepted);
    if (rr.move.accepted) {
        CHECK(rr.move.noPassThrough);
        CHECK(rr.move.delta < rr.D0 + pp.clearance);
        CHECK(rr.move.localScale >= 0.5 * kPi * rr.move.delta - 1e-9);  // the uncertified same-yarn scale is reported
        CHECK(validate(c.cs, pp).verdict == Status::Pass);
    } else {
        CHECK(c.cs.yarns[0].curve.ctrl == before);
        CHECK(validate(c.cs, atD0).verdict == Status::Pass);  // the start the motion would have been certified at
    }
}

TEST(repair_start_diameter_from_the_report) {
    // The start diameter comes from the report at the target, not from a search: below every distance
    // minus c and the allowance, below 2 rho_lo, and it validates (the cylinder at D = 1.3).
    CylinderPiece c = cylinderPiece(24, 5);
    PatternParams pp = c.pp;
    pp.yarnDiameter = 1.3;
    const ValidationReport atTarget = validate(c.cs, pp);
    const RepairReport rr = certifiedRepair({c.cs, pp, nullptr}, RepairParams{}, knit_test::cpuNumerics());
    CHECK(rr.D0 > 0);
    const double allowance = pp.clearance + atTarget.clearance.deltaD;
    CHECK(rr.D0 < atTarget.clearance.minInterYarn - allowance || atTarget.clearance.minInterYarn > 1e299);
    CHECK(rr.D0 < atTarget.clearance.minSelfNonlocal - allowance || atTarget.clearance.minSelfNonlocal > 1e299);
    for (const CurvatureResult& cr : atTarget.curvature) CHECK(rr.D0 < 2 * cr.rhoLo);
    std::printf("    D0 %.4f mm (inter-yarn %.4f, self %.4f, c %.3f, deltaD %.2g)\n", rr.D0, atTarget.clearance.minInterYarn,
                atTarget.clearance.minSelfNonlocal, pp.clearance, atTarget.clearance.deltaD);
}

TEST(final_surface_gate) {
    // Gate (d) on round tubes: the cylinder spiral passes gate (a) at D = 1, and the
    // sampled surface check of its round tubes passes too, at the default model
    // spacing and at half of it (the verdict is stable under finer sampling). A
    // thicker yarn on the same centreline (D = 1.6) makes the tubes overlap and the
    // surface check reports it.
    CylinderPiece c = cylinderPiece(24, 4);
    CHECK(validate(c.cs, c.pp).verdict == Status::Pass);
    for (double spacing : {0.025, 0.0125}) {
        SurfaceParams sp;
        sp.compression = false;
        sp.ringSpacing = spacing;
        SurfaceModel sm = buildSurface(c.cs, c.pp, sp);
        SurfaceCheck sc = checkSurface(c.cs, sm);
        CHECK(sc.fold == Status::Pass && sc.collisions == Status::Pass && sc.pointsTested > 100000);
    }
    PatternParams thick = c.pp;
    thick.yarnDiameter = 1.6;
    SurfaceParams sp;
    sp.compression = false;
    SurfaceModel sm = buildSurface(c.cs, thick, sp);
    SurfaceCheck sc = checkSurface(c.cs, sm);
    CHECK(sc.collisions == Status::Fail || sc.fold == Status::Fail);
}

// ---- S2: strand export contract ----

TEST(pptstrands_roundtrip) {
    const std::string d = scratchDir("roundtrip");
    StrandSetF s = tinyStrands();
    StrandCheck c;
    std::string err;
    CHECK(writeStrandsFile(d + "/a.pptstrands", s, StrandRules{}, &err, &c));
    CHECK(c.ok);
    CHECK(c.precisionWarnings == 0);
    CHECK(std::filesystem::file_size(d + "/a.pptstrands") == 40 + 16 * s.vertices() + 4 * (s.strands() + 1));
    StrandSetF t;
    CHECK(readStrandsFile(d + "/a.pptstrands", t, &err));
    CHECK(t.xyz == s.xyz && t.radius == s.radius && t.first == s.first);
    CHECK(t.strands() == 2 && t.vertices() == 9);
    // Trailing bytes are refused by the reader, like the renderer's.
    { std::ofstream f(d + "/a.pptstrands", std::ios::binary | std::ios::app); f.put('x'); }
    CHECK(!readStrandsFile(d + "/a.pptstrands", t, &err));
}

// ---- S2: fibres on surface pieces ----

void printFibreReport(const FibreReport& r) {
    std::printf("  %d candidates: %d back, %d origin inside, %d occluded, %d too short, %d yarn end, %d root, %d serialized, %d rules, %d no bend, %d too sharp, %d vertex inside, %d no emergence, "
                "%d approx, %d emergence, %d no clearance, %d collision, %d truncated too short, %d approaches, %d hidden, %d ply inner, %d behind -> "
                "%d accepted (%d truncated)\n",
                r.candidates, r.backFacing, r.originInside, r.occluded, r.tooShortToBend, r.yarnEnd, r.rootRejected, r.serializedRejected, r.ruleRejected,
                r.noBend, r.turnTooSharp, r.vertexInside, r.noEmergence, r.approxError, r.emergeRejected, r.noClearance, r.collisionRejected, r.truncatedTooShort, r.approaches,
                r.hiddenInPiece, r.plyInner, r.behindFabric, r.accepted, r.truncated);
}

// Capped-sphere spiral (settled), its surface model, its exported tubes and fibres with the piece preset.
struct SphereFibres {
    FixturePiece f;
    PatternParams pp;
    PieceKnit pk;
    SurfaceModel sm;
    FibreParams fp;
    FibreFrame fr;
    FibreSet fs;
    SceneMapping map;
    double eSer = 0;
    std::unique_ptr<ExportTubeIndex> tubes;
};
const SphereFibres& sphereFibres() {
    static std::unique_ptr<SphereFibres> c;
    if (c) return *c;
    c = std::make_unique<SphereFibres>();
    c->f = fixturePiece(25.0);
    PieceKnitParams kp;
    c->pp.waleSpacing = kp.g.w;
    c->pp.courseSpacing = kp.g.h;
    CHECK(knitPiece({*c->f.surf, true}, {c->pp, kp}, knit_test::cpuNumerics(), c->pk));
    settle({c->pk.curves, c->pp, c->f.surf.get()}, SettleParams{}, knit_test::cpuNumerics());
    SurfaceParams sp;
    sp.compression = false;
    c->sm = buildSurface(c->pk.curves, c->pp, sp);
    c->fp = pieceFibres();
    c->fp.density = 20.0;  // properties, not density, are tested here
    c->fp.randomFraction = 0;  // the directed S2b set (legacy fixture (a)); mixtures have their own tests
    c->fp.rootDepth = 0.126;
    c->map.sigma = 68.103;
    c->map.centre = {0.3, -1.2, 0.8};
    c->eSer = 1e-4;  // the float envelope as knit_check computes it, rounded up
    c->fp.serializedMargin = c->eSer;
    c->fp.exportChordFactor = kExportRing * kExportRing / 8;
    const SceneMapping map = c->map;
    c->tubes = std::make_unique<ExportTubeIndex>(c->pk.curves, c->pp, ExportTubeSpec{kExportModelRing, kExportRing, kExportSides, false,
                                                 [map](const Vec3& p) { return map.toExportF(p); }, c->eSer});
    const PieceSurface* surf = c->f.surf.get();
    c->fr.fabric = [surf](const Vec3& p) {
        const SurfaceSample q = surf->project(p);
        return FabricPoint{q.p, q.n};
    };
    c->fr.exposureLength = 2 * (c->pp.depthAmplitude + 1.0) + c->fp.lengthMax;
    c->fr.backDepth = c->pp.depthAmplitude * c->pp.yarnDiameter + c->pp.radius();
    c->fr.pieceSeed = 1;
    c->fs = growYarnFibres(c->pk.curves, c->sm, c->pp, c->fp, &c->fr);
    printFibreReport(c->fs.rep);
    return *c;
}

TEST(fibers_surface_all_sides) {
    // With the outward frame, fibres grow on every side of a closed piece (per unit
    // area, the lower and upper halves agree within 25 %), and every anchor faces
    // outwards (n . outward > frontCos). The flat patch's fixed +z would starve the
    // lower half.
    const SphereFibres& s = sphereFibres();
    CHECK(s.fs.rep.accepted > 1000);
    const Vec3 c = s.f.mesh.centroid();
    double areaUp = 0, areaDown = 0;
    for (size_t t = 0; t < s.f.mesh.tri.size(); ++t) {
        const auto& T = s.f.mesh.tri[t];
        const Vec3 m = (s.f.mesh.pos[T[0]] + s.f.mesh.pos[T[1]] + s.f.mesh.pos[T[2]]) * (1.0 / 3);
        const double a = 0.5 * norm(cross(s.f.mesh.pos[T[1]] - s.f.mesh.pos[T[0]], s.f.mesh.pos[T[2]] - s.f.mesh.pos[T[0]]));
        (m.z > c.z ? areaUp : areaDown) += a;
    }
    int up = 0, down = 0;
    for (const Fibre& f : s.fs.fibres) {
        (f.anchor.z > c.z ? up : down) += 1;
        CHECK(dot(toVec3(f.normal), normalized(s.fr.fabric(toVec3(f.anchor)).n)) > s.fp.frontCos);
    }
    const double ratio = (up / areaUp) / (down / areaDown);
    std::printf("  up %d, down %d fibres; density ratio %.3f\n", up, down, ratio);
    CHECK(ratio > 0.75 && ratio < 1.0 / 0.75);
    FibreSet flat = growYarnFibres(s.pk.curves, s.sm, s.pp, s.fp, nullptr);
    int flatDown = 0;
    for (const Fibre& f : flat.fibres) flatDown += f.anchor.z < c.z;
    CHECK(flatDown < down / 4);
}

TEST(fibers_segment_interior_surface) {
    // Endpoints outside every tube, interior through a DISTANT turn of the same yarn: the
    // reference check rejects it. At the attachment-window boundary: a fibre leaving the own tube
    // passes, one continuing into it past the window fails.
    const SphereFibres& s = sphereFibres();
    const YarnPath& y = s.pk.curves.yarns[0];
    const double r = s.pp.radius();
    const double sC = 0.5 * y.length();
    const double u = y.paramAtArcLength(sC);
    const Vec3 c = y.curve.eval(u), T = normalized(y.curve.d1(u));
    const Vec3 n = normalized(s.fr.fabric(c).n);
    const Vec3 e = normalized(n - T * dot(n, T));  // perpendicular to the yarn, outward
    const Vec3 a = c + e * (r + 0.12), b = c - e * (r + 0.12);
    SurfaceIndex idx(s.pk.curves, s.sm);
    CHECK(!idx.exactInside(a, 0, {}, {}).inside);
    FibreSet fake;
    Fibre f;
    f.yarnId = 0;
    f.s = sC - 40.0;
    f.radius0 = f.radius1 = 0.01f;
    f.attach = 0.15f;
    f.attachArc = 0.3f;
    fake.addFibre(f, std::vector<Vec3>{a, b});
    // Anchored ON this yarn at c + e r: outwards is fine, inwards fails past attachLength.
    Fibre out = f, in = f;
    out.s = in.s = sC;
    fake.addFibre(out, std::vector<Vec3>{c + e * r, c + e * (r + 0.3)});
    fake.addFibre(in, std::vector<Vec3>{c + e * r, c + e * (r - f.attach - 0.15)});  // continues past the window
    FibreReferenceResult rr = fibreReferenceCheck(*s.tubes, fake);
    CHECK(rr.fail == 2 && rr.ok == 1);
    CHECK(rr.failing.size() == 2 && rr.failing[0] == 0 && rr.failing[1] == 2);
}

TEST(fibers_root_attached_exported) {
    // Roots strictly inside the EXPORTED polygonal tube (8 sides, rings every 0.4 mm,
    // positions mapped to the scene and rounded through float32), tips outside it.
    const SphereFibres& s = sphereFibres();
    SceneMapping map;
    map.sigma = 68.103;
    map.centre = {0.3, -1.2, 0.8};
    const double margin = 1e-4 / map.sigma;
    AttachmentResult ar = rootAttachmentCheck(*s.tubes, s.fs, margin);
    std::printf("  %d/%d roots inside (min depth %.4f mm), %d tips outside\n", ar.rootInside, ar.fibres,
                ar.minRootDepth * map.sigma, ar.tipOutside);
    CHECK(ar.fibres == static_cast<int>(s.fs.fibres.size()));
    CHECK(ar.rootInside == ar.fibres);
    CHECK(ar.tipOutside == ar.fibres);
    // Without the buried root the polyline would start at the emergence E*, ON the
    // smooth surface, outside the inscribed polygon: not attached to the exported mesh.
    FibreSet bare = s.fs;
    for (const Fibre& f : bare.fibres) bare.points[f.firstPoint] = f.anchor;
    AttachmentResult br = rootAttachmentCheck(*s.tubes, bare, margin);
    CHECK(br.rootInside < br.fibres / 2);
}

// The analytic mean of fibreTargetLength: the truncated log-normal's, or the clamped exponential's.
double fibreTargetLengthMean(const FibreParams& fp) {
    const double m = fp.lengthMean;
    if (fp.lengthSpread > 0) {
        const double mu = std::log(m), sg = fp.lengthSpread;
        const double a = (std::log(fp.lengthMin) - mu) / sg, b = (std::log(fp.lengthMax) - mu) / sg;
        return std::exp(mu + 0.5 * sg * sg) * (normalCdf(b - sg) - normalCdf(a - sg)) / (normalCdf(b) - normalCdf(a));
    }
    return fp.lengthMin + m * (std::exp(-fp.lengthMin / m) - std::exp(-fp.lengthMax / m));
}

TEST(fibers_target_length_distribution) {
    // Before any collision test: clamped exponential, mean within 3 % of the analytic
    // value, the clamps hold, and the mass at lengthMin matches 1 - exp(-min / mean).
    FibreParams fp = pieceFibres();
    fp.lengthSpread = 0;  // the S2 exponential (kept for the flat and S2 parameter sets)
    fp.lengthMean = 1.0;
    fp.lengthMin = 0.3;
    fp.lengthMax = 2.0;
    const int n = 200000;
    double sum = 0, lo = 1e300, hi = 0;
    int atMin = 0;
    for (int i = 0; i < n; ++i) {
        const double L = fibreTargetLength(fp, 12345, 0, i);
        sum += L;
        lo = std::min(lo, L);
        hi = std::max(hi, L);
        atMin += L == fp.lengthMin;
    }
    CHECK(std::fabs(sum / n / fibreTargetLengthMean(fp) - 1) < 0.03);
    CHECK(lo >= fp.lengthMin && hi <= fp.lengthMax);
    CHECK(std::fabs(double(atMin) / n - (1 - std::exp(-fp.lengthMin / fp.lengthMean))) < 0.01);
}

TEST(fibers_seed_per_piece) {
    // Every spiral is yarn 0: without the piece seed all pieces would repeat one
    // sequence of anchors.
    const SphereFibres& s = sphereFibres();
    FibreParams fp = pieceFibres();
    CHECK(fibreSeed(fp, 1) != fibreSeed(fp, 2));
    CHECK(fibreSeed(fp, 0) == fp.seed);
    FibreFrame fr2 = s.fr;
    fr2.pieceSeed = 2;
    FibreSet other = growYarnFibres(s.pk.curves, s.sm, s.pp, s.fp, &fr2);
    CHECK(!other.fibres.empty() && !s.fs.fibres.empty());
    CHECK(other.fibres.front().s != s.fs.fibres.front().s);
}

TEST(inside_mesh_matches_winding) {
    // The accelerated pseudo-normal inside test agrees with the winding number on a
    // closed mesh (points within 1e-6 of the surface skipped).
    FixturePiece f = fixturePiece(25.0);
    MeshQuery mq(f.mesh);
    Vec3 lo, hi;
    f.mesh.bounds(lo, hi);
    const Vec3 ext = hi - lo;
    int checked = 0, inside = 0;
    for (int i = 0; i < 4000; ++i) {
        const double u = std::fmod(0.5 + i * 0.6180339887, 1.0), v = std::fmod(0.5 + i * 0.7548776662, 1.0),
                     w = std::fmod(0.5 + i * 0.5698402910, 1.0);
        const Vec3 q = lo - ext * 0.2 + Vec3{ext.x * 1.4 * u, ext.y * 1.4 * v, ext.z * 1.4 * w};
        const MeshQuery::Hit h = mq.closest(q);
        if (h.dist < 1e-6) continue;
        ++checked;
        const bool wn = windingNumber(f.mesh, q) > 0.5;
        inside += wn;
        CHECK(mq.inside(q) == wn);
    }
    CHECK(checked > 3900 && inside > 500);
}

// ---- S2b: hair flyaways in yarn coordinates ----

double ksUniform(std::vector<double> x, const std::function<double(double)>& F) {
    std::sort(x.begin(), x.end());
    double d = 0;
    const double n = double(x.size());
    for (size_t i = 0; i < x.size(); ++i) d = std::max({d, std::fabs(F(x[i]) - i / n), std::fabs(F(x[i]) - (i + 1) / n)});
    return d;
}
// E[z^2] of a standard normal clamped to +-c.
double clampedNormalVar(double c) {
    const double phi = std::exp(-0.5 * c * c) / std::sqrt(2 * kPi);
    return (2 * normalCdf(c) - 1) - 2 * c * phi + 2 * c * c * (1 - normalCdf(c));
}

// ---- S3: explicit plies ----

YarnPath arcYarn(double R, double degrees);

TEST(ply_interpolation) {
    // The interpolating B-spline passes through its points (u = i), with zero second
    // difference of the control polygon at both ends.
    std::vector<Vec3> P;
    for (int i = 0; i <= 40; ++i) P.push_back({0.3 * i, std::sin(0.4 * i), 0.2 * std::cos(0.7 * i)});
    const BSpline b = interpolateBSpline(P);
    CHECK(b.numSpans() == 40);
    for (int i = 0; i <= 40; ++i) CHECK(norm(b.eval(i) - P[i]) < 1e-12);
    for (const std::vector<Vec3>& Q : {std::vector<Vec3>{{0, 0, 0}, {1, 2, 3}}, std::vector<Vec3>{{0, 0, 0}, {1, 0, 0}, {2, 1, 0}}}) {
        const BSpline c = interpolateBSpline(Q);
        for (size_t i = 0; i < Q.size(); ++i) CHECK(norm(c.eval(double(i)) - Q[i]) < 1e-12);
    }
}

TEST(ply_containment) {
    // Capped sphere, 3 plies: every ply centre stays within e of the circle of radius R
    // around the yarn centreline (so the ply tube within r + e), e < c / 2; plies of one
    // yarn keep >= 2 r_p between their centres at equal s.
    const SphereFibres& s = sphereFibres();
    PlyGeometry g;
    const PlySet ps = buildPlies(s.pk.curves, s.sm, g, s.pp.radius());
    std::printf("  %zu plies, r_p %.4f, R %.4f; containment e %.2e mm; min centre distance %.4f (2 r_p %.4f); %.0f ms\n",
                ps.curves.yarns.size(), ps.layout.rp, ps.layout.R, ps.containment, ps.minCentreDistance, 2 * ps.layout.rp, ps.ms);
    CHECK(ps.curves.yarns.size() == 3 * s.pk.curves.yarns.size());
    CHECK_NEAR(ps.layout.rp0 + ps.layout.R, s.pp.radius(), 1e-12);
    CHECK_NEAR(2 * ps.layout.R * std::sin(kPi / 3), 2 * ps.layout.rp0, 1e-12);  // touching packing
    CHECK(ps.containment < 0.5 * s.pp.clearance);
    CHECK(ps.minCentreDistance >= 2 * ps.layout.rp);
    // Independent dense check: points on the ply tube surfaces lie within r + e of the yarn
    // centreline (brute force over a dense centreline polyline near each sample).
    const YarnPath& y = s.pk.curves.yarns[0];
    std::vector<Vec3> dense;
    for (int i = 0; i <= 200000; ++i) dense.push_back(y.curve.eval(y.uEnd() * i / 200000.0));
    SpatialHash3 h(1.0);
    auto lo3 = [](const Vec3& a, const Vec3& b) { return Vec3{std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z)}; };
    auto hi3 = [](const Vec3& a, const Vec3& b) { return Vec3{std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)}; };
    for (size_t i = 1; i < dense.size(); ++i) h.insert(lo3(dense[i - 1], dense[i]), hi3(dense[i - 1], dense[i]));
    double worst = -1e300;
    for (const YarnPath& p : ps.curves.yarns)
        for (double u = 5; u < p.uEnd() - 5; u += 7.31) {
            const Vec3 c = p.curve.eval(u), T = normalized(p.curve.d1(u));
            Vec3 e1 = normalized(cross(T, std::fabs(T.z) < 0.9 ? Vec3{0, 0, 1} : Vec3{1, 0, 0})), e2 = cross(T, e1);
            for (int k = 0; k < 8; ++k) {
                const Vec3 q = c + (e1 * std::cos(k * kPi / 4) + e2 * std::sin(k * kPi / 4)) * ps.layout.rp;
                double d = 1e300;
                const Vec3 w{1.5, 1.5, 1.5};
                h.query(q - w, q + w, [&](int i) {
                    d = std::min(d, closestSegmentPoints(q, q, dense[i], dense[i + 1]).distance);
                });
                worst = std::max(worst, d - s.pp.radius());
            }
        }
    std::printf("  ply surface points beyond the yarn envelope by at most %.2e mm\n", worst);
    CHECK(worst <= ps.containment + 1e-4);  // + the dense polyline's sagitta
}

TEST(ply_single_equals_tube) {
    // count = 1 is the plain yarn: the same curves, the same tube mesh.
    const SphereFibres& s = sphereFibres();
    PlyGeometry g;
    g.count = 1;
    const PlySet ps = buildPlies(s.pk.curves, s.sm, g, s.pp.radius());
    CHECK(ps.curves.yarns.size() == s.pk.curves.yarns.size());
    CHECK(ps.curves.yarns[0].curve.ctrl == s.pk.curves.yarns[0].curve.ctrl);
    TubeParams tp;
    const TubeMesh a = buildTubeMesh(s.pk.curves, s.sm, tp), b = buildTubeMesh(ps.curves, s.sm, tp);
    CHECK(a.pos.size() == b.pos.size() && a.tri == b.tri);
    for (size_t i = 0; i < a.pos.size(); ++i) CHECK(a.pos[i].x == b.pos[i].x && a.pos[i].y == b.pos[i].y && a.pos[i].z == b.pos[i].z);
}

TEST(ply_curvature_fold) {
    // Ply folds are decided by the certified curvature check (rho_lo >= r_p). On the
    // capped sphere every ply passes. A nearly untwisted ply on the inner side of a bend
    // tighter than r (centreline radius 0.9 r, pitch 50 mm, fill 1) has rho ~ 0.9 r - R <
    // r_p0 there and is flagged. (Twist does not add curvature at the inner side: the
    // normal curvature there is cos^2(phi) / (rho - R) + sin^2(phi) / R <= 1 / (rho - R).)
    const SphereFibres& s = sphereFibres();
    PlyGeometry g;
    const PlySet ps = buildPlies(s.pk.curves, s.sm, g, s.pp.radius());
    const Allowances al = Allowances::forExtent(ps.curves.extent());
    int fail = 0, unres = 0;
    double rhoLo = 1e300;
    for (const YarnPath& p : ps.curves.yarns) {
        const CurvatureResult cr = checkCurvature(p, ps.layout.rp, al);
        fail += cr.noFold == Status::Fail;
        unres += cr.noFold == Status::Unresolved;
        rhoLo = std::min(rhoLo, cr.rhoLo);
    }
    std::printf("  sphere plies: %d FAIL, %d UNRESOLVED, rho_lo %.4f vs r_p %.4f\n", fail, unres, rhoLo, ps.layout.rp);
    CHECK(fail == 0 && unres == 0);
    PatternParams pp;
    CurveSet cs = makeSet({arcYarn(pp.radius() * 0.9, 300)});
    SurfaceParams sp;
    sp.compression = false;
    const SurfaceModel sm = buildSurface(cs, pp, sp);
    PlyGeometry tight;
    tight.fill = 1.0;
    tight.pitch = 50.0;
    const PlySet pt = buildPlies(cs, sm, tight, pp.radius());
    int tf = 0;
    for (const YarnPath& p : pt.curves.yarns) {
        const CurvatureResult cr = checkCurvature(p, pt.layout.rp, Allowances::forExtent(cs.extent()));
        std::printf("  ply %s: %d spans, rho_est %.4f, rho_lo %.4f, noFold %s (r_p %.4f)\n", p.name.c_str(), p.curve.numSpans(),
                    cr.rhoEst, cr.rhoLo, toString(cr.noFold), pt.layout.rp);
        tf += cr.noFold == Status::Fail;
    }
    std::printf("  rho = 0.9 r bend, untwisted, fill 1: %d of %zu plies fold\n", tf, pt.curves.yarns.size());
    CHECK(tf > 0);
}

TEST(fibers_on_plies) {
    // Fibres on the 3 plies of the capped sphere: anchored on and buried in the plies, authored in
    // the yarn's coordinates, so they lie closer to the yarn axis than the plies do.
    const SphereFibres& s = sphereFibres();
    PlyGeometry g;
    const PlySet ps = buildPlies(s.pk.curves, s.sm, g, s.pp.radius());
    PatternParams pp = s.pp;
    pp.yarnDiameter = 2 * ps.layout.rp;
    SurfaceParams sp;
    sp.compression = false;
    const SurfaceModel sm = buildSurface(ps.curves, pp, sp);
    const double ring = kExportPlyRing;
    const int sides = kExportPlySides;
    double kappa = 0;
    for (const YarnPath& y : ps.curves.yarns)
        kappa = std::max(kappa, std::min(1e3, checkCurvature(y, pp.radius(), Allowances::forExtent(ps.curves.extent())).kappaHi));
    FibreParams fp = s.fp;
    fp.density = 20.0;
    fp.rootDepth = pp.radius() * (1 - std::cos(kPi / sides)) + ring * ring * kappa * (1 + pp.radius() * kappa) / 8 + 0.02;
    fp.exportChordFactor = ring * ring / 8;
    FibreFrame fr = s.fr;
    fr.parentCurves = &s.pk.curves;
    fr.parentSurface = &s.sm;
    fr.plies = &ps;
    const FibreSet fs = growYarnFibres(ps.curves, sm, pp, fp, &fr);
    printFibreReport(fs.rep);
    CHECK(fs.fibres.size() > 1000);
    std::vector<double> helix;  // a DIRECTED property (random fibres have no helix model)
    for (const Fibre& f : fs.fibres)
        if (f.population == FibrePopulation::Directed) helix.push_back(f.helixActual);
    std::nth_element(helix.begin(), helix.begin() + helix.size() / 2, helix.end());
    const double alpha = plyHelixAngle(g, s.pp.radius()) * 180 / kPi;
    std::printf("  %zu fibres; median angle to the yarn axis %.1f deg (ply helix %.1f)\n", fs.fibres.size(), helix[helix.size() / 2], alpha);
    CHECK(helix[helix.size() / 2] < alpha);
}

TEST(curvature_degenerate_hodograph) {
    // The first span of a leg's ply (an interpolating spline's natural end: control
    // points 0-2 exactly collinear, two equal hodograph controls). The certified speed
    // bound once divided by the zero area of the hodograph hull and left the span
    // without a curvature bound; it certifies.
    YarnPath y = makeYarn({{-36.851141902634268, -38.425005369105335, 2.3503366918914557},
                           {-36.751977075414217, -38.443851764246844, 2.2930168620698708},
                           {-36.652812248194167, -38.462698159388353, 2.2356970322482859},
                           {-36.547952612554511, -38.486701104314605, 2.196183903655625}}, 0);
    const CurvatureResult cr = checkCurvature(y, 0.247, Allowances::forExtent(109.37));
    std::printf("  speedLo %.6g kappaHi %.6g rhoLo %.6g regular %s noFold %s leaves %d\n", cr.speedLo, cr.kappaHi, cr.rhoLo,
                toString(cr.regular), toString(cr.noFold), cr.leaves);
    CHECK(cr.regular == Status::Pass && cr.noFold == Status::Pass);
    CHECK(cr.kappaHi < 2.0 && cr.spanKappaHi[0] < 2.0);  // the true maximum is 1.464 at u = 1
}

TEST(fiber_draw_orientation) {
    // Before any rejection: psi ~ N(psi0, sigma) clamped at 2.5 sigma, lift uniform,
    // direction 50/50, bow amplitudes in range.
    const FibreParams fp = pieceFibres();
    const int n = 100000;
    double sq = 0;
    int plus = 0;
    std::vector<double> lift;
    for (int i = 0; i < n; ++i) {
        const FibreDraw d = drawFibre(fp, 777, 0, i);
        const double dev = d.psi / (kPi / 180) - fp.helixDeg;
        CHECK(std::fabs(dev) <= 2.5 * fp.helixSigmaDeg + 1e-9);
        sq += dev * dev;
        plus += d.sgn > 0;
        lift.push_back(d.lift / (kPi / 180));
        CHECK(std::fabs(d.bowTheta) >= 0.5 * fp.bowTheta - 1e-12 && std::fabs(d.bowTheta) <= fp.bowTheta + 1e-12);
        CHECK(std::fabs(d.bowRadial) <= fp.bowRadial + 1e-12);
    }
    const double rms = std::sqrt(sq / n), expect = fp.helixSigmaDeg * std::sqrt(clampedNormalVar(2.5));
    const double ks = ksUniform(lift, [&](double x) { return (x - fp.liftMinDeg) / (fp.liftMaxDeg - fp.liftMinDeg); });
    std::printf("  rms(psi - psi0) %.3f deg (expected %.3f), lift KS %.4f, +s %.4f\n", rms, expect, ks, double(plus) / n);
    CHECK(std::fabs(rms / expect - 1) < 0.05);
    CHECK(ks <= 0.015);
    CHECK(std::fabs(double(plus) / n - 0.5) <= 0.015);
}

TEST(fiber_draw_thickness) {
    const FibreParams fp = pieceFibres();
    const int n = 100000;
    double sq = 0;
    for (int i = 0; i < n; ++i) {
        const double f = drawFibre(fp, 778, 0, i).radiusFactor;
        CHECK(f >= std::exp(-2 * fp.thicknessSigma) - 1e-12 && f <= std::exp(2 * fp.thicknessSigma) + 1e-12);
        sq += std::log(f) * std::log(f);
    }
    const double sd = std::sqrt(sq / n), expect = fp.thicknessSigma * std::sqrt(clampedNormalVar(2.0));
    std::printf("  sd(ln r) %.5f (expected %.5f)\n", sd, expect);
    CHECK(std::fabs(sd / expect - 1) < 0.05);
}

TEST(fibers_length_truncated_lognormal) {
    // Inverse-CDF sampling: KS distance to the analytic truncated log-normal, no mass
    // at the bounds (no clamping), mean within 1 % of the analytic mean.
    const FibreParams fp = pieceFibres();
    const int n = 100000;
    std::vector<double> L;
    int atBound = 0;
    double sum = 0;
    for (int i = 0; i < n; ++i) {
        const double l = fibreTargetLength(fp, 4242, 0, i);
        L.push_back(l);
        sum += l;
        atBound += l == fp.lengthMin || l == fp.lengthMax;
        CHECK(l >= fp.lengthMin && l <= fp.lengthMax);
    }
    const double mu = std::log(fp.lengthMean), sg = fp.lengthSpread;
    const double Fa = normalCdf((std::log(fp.lengthMin) - mu) / sg), Fb = normalCdf((std::log(fp.lengthMax) - mu) / sg);
    const double ks = ksUniform(L, [&](double x) { return (normalCdf((std::log(x) - mu) / sg) - Fa) / (Fb - Fa); });
    std::printf("  KS %.4f, at a bound %d, mean %.4f (analytic %.4f)\n", ks, atBound, sum / n, fibreTargetLengthMean(fp));
    CHECK(ks <= 0.01);
    CHECK(atBound == 0);
    CHECK(std::fabs(sum / n / fibreTargetLengthMean(fp) - 1) < 0.01);
    for (double p : {1e-6, 0.01, 0.3, 0.5, 0.9, 0.999}) CHECK_NEAR(normalCdf(normalQuantile(p)), p, 1e-12);
}

// Analytic derivative of the yarn-coordinate curve against a central difference, on a
// straight round yarn, a curved spiral and elliptical (compressed) rings.
TEST(fiber_curve_derivative) {
    auto check = [](const CurveSet& cs, const SurfaceModel& sm, double sE) {
        const YarnSurface& ys = sm.yarns[0];
        DirectedCurve fc;
        fc.cs = &cs;
        fc.ys = &ys;
        fc.sE = sE;
        fc.thetaE = 0.9;
        fc.S = 0.9;
        fc.p = std::tan(10 * kPi / 180);
        fc.q = std::tan(14 * kPi / 180);
        fc.bowTheta = 0.09;
        fc.bowLambda = 0.05;
        fc.r = 0.5;
        int tested = 0;
        for (int sgn : {1, -1}) {
            fc.sgn = sgn;
            for (double u = 0.05; u < 1; u += 0.1) {
                const double h = 1e-6;
                if (std::floor(fc.coords(u - h).s / ys.ds) != std::floor(fc.coords(u + h).s / ys.ds)) continue;  // a model-ring kink
                const Vec3 fd = (fc.X(u + h) - fc.X(u - h)) / (2 * h), an = fc.dX(u);
                CHECK(norm(fd - an) <= 1e-6 * std::max(1.0, norm(an)));
                ++tested;
            }
        }
        CHECK(tested >= 16);
    };
    PatternParams pp;
    {
        CurveSet cs = makeSet({makeYarn(lineCtrl({0, 0, 0}, {10, 0, 0}, 10), 0)});
        SurfaceParams sp;
        sp.compression = false;
        check(cs, buildSurface(cs, pp, sp), 4.0);
    }
    {
        const SphereFibres& s = sphereFibres();
        check(s.pk.curves, s.sm, 0.4 * s.pk.curves.yarns[0].length());
    }
    {
        PatternParams tp = testParams();
        CurveSet cs = makeSet({makeYarn(lineCtrl({0, 0, 0}, {6, 0, 0}, 6), 0), makeYarn(lineCtrl({0, 0.7, 0}, {6, 1.1, 0}, 6), 1)});
        SurfaceModel sm = buildSurface(cs, tp, SurfaceParams{});
        CHECK(sm.compression.minAxisRatio < 1);
        const Ring R = ringAt(cs, sm.yarns[0], 2.0);
        CHECK(R.b < R.a);  // elliptical here
        check(cs, sm, 2.0);
    }
}

TEST(fiber_curve_yarn_coords) {
    // Straight round yarn: X(lambda = 1) is the surface; the physical helix and lift
    // angles measured from dX match their closed forms (lambda and the bow included);
    // lambda increases strictly. psi and e are only the nominal slopes.
    CurveSet cs = makeSet({makeYarn(lineCtrl({0, 0, 0}, {10, 0, 0}, 10), 0)});
    PatternParams pp;
    SurfaceParams sp;
    sp.compression = false;
    const SurfaceModel sm = buildSurface(cs, pp, sp);
    const YarnSurface& ys = sm.yarns[0];
    const double r = pp.radius();
    double worst = 0;
    for (double s : {1.013, 3.3, 7.77})
        for (double th : {0.3, 2.0, 4.5})
            worst = std::max(worst, norm(yarnCoordinates(cs, ys, s, th, 1.0).X - evalSurface(cs, ys, s, th).p));
    std::printf("  yarn coordinates vs evalSurface: %.2e mm\n", worst);
    CHECK(worst < 1e-4);
    DirectedCurve fc;
    fc.cs = &cs;
    fc.ys = &ys;
    fc.sE = 4.0;
    fc.thetaE = 0.7;
    fc.S = 0.8;
    fc.p = std::tan(12 * kPi / 180);
    fc.q = std::tan(15 * kPi / 180);
    fc.r = r;
    const double deg = 180 / kPi;
    for (int bow = 0; bow < 2; ++bow)
        for (int sgn : {1, -1}) {
            fc.sgn = sgn;
            fc.bowTheta = bow ? 0.08 : 0;
            fc.bowLambda = bow ? 0.04 : 0;
            double prevLam = 0;
            for (double u = 0; u <= 1 + 1e-12; u += 0.05) {
                const YarnCoords y = fc.coords(u);
                const double s = y.s, lam = y.lambda, dth = y.dTheta, dlam = y.dLambda;
                CHECK(u == 0 || lam > prevLam);
                prevLam = lam;
                const Vec3 X = fc.X(u), d = fc.dX(u), c{s, 0, 0}, T{1, 0, 0};
                const Vec3 er = normalized(X - c), et = cross(T, er);
                CHECK_NEAR(norm(X - c), lam * r, 1e-6);
                const double helix = std::atan2(std::fabs(dot(d, et)), std::fabs(dot(d, T))) * deg;
                const double lift = std::atan2(dot(d, er), std::hypot(dot(d, T), dot(d, et))) * deg;
                const double helixC = std::atan(lam * r * std::fabs(dth) / fc.S) * deg;
                const double liftC = std::atan(r * dlam / std::hypot(fc.S, lam * r * dth)) * deg;
                CHECK(std::fabs(helix - helixC) < 0.1 && std::fabs(lift - liftC) < 0.1);
                if (!bow) CHECK(std::fabs(lift - std::atan(fc.q / std::sqrt(1 + lam * lam * fc.p * fc.p)) * deg) < 0.1);
            }
        }
    // Generated on this yarn: the visible polyline length matches the target (0.5 %).
    FibreParams fp = pieceFibres();
    fp.randomFraction = 0;  // the directed curve is tested here
    fp.density = 30;
    fp.rootDepth = 0.126;
    const FibreSet fs = growYarnFibres(cs, sm, pp, fp, nullptr);
    printFibreReport(fs.rep);
    int n = 0;
    double worstL = 0;
    for (const Fibre& f : fs.fibres) {
        if (f.truncated) continue;
        ++n;
        worstL = std::max(worstL, std::fabs(visibleLength(fs, f) / f.targetLength - 1));
    }
    std::printf("  %d untruncated fibres, visible length vs target within %.4f\n", n, worstL);
    CHECK(n > 50);
    CHECK(worstL <= 0.005);
}

TEST(fiber_polyline_fit) {
    // Every accepted fibre: curve-to-polyline deviation <= maxCurveDev and segments >= 4.2 r_i;
    // a shape forced beyond the allowance is rejected.
    const SphereFibres& s = sphereFibres();
    CHECK(s.fs.fibres.size() > 300);
    double maxDev = 0;
    for (const Fibre& f : s.fs.fibres) {
        maxDev = std::max(maxDev, double(f.maxDev));
        CHECK(f.maxDev <= s.fp.maxCurveDev);
        const std::vector<Vec3> pts = fibrePoints(s.fs, f);
        for (size_t k = 1; k < pts.size(); ++k) CHECK(norm(pts[k] - pts[k - 1]) >= 4.2 * f.radius0 - 1e-5);
    }
    std::printf("  max curve deviation %.4f mm over %zu fibres\n", maxDev, s.fs.fibres.size());
    FibreParams tight = s.fp;
    tight.density = 4;
    tight.maxCurveDev = 1e-5;
    FibreFrame fr = s.fr;
    const FibreSet t = growYarnFibres(s.pk.curves, s.sm, s.pp, tight, &fr);
    CHECK(t.rep.accepted == 0 && t.rep.approxError > 0);
}

TEST(fibers_bent_all_accepted) {
    // Every accepted fibre (truncated ones included): >= 2 segments, segment 1 with a visible part,
    // interior turns in [3, 35] deg; the strand rules hold; the exported floats ARE the vertices.
    const SphereFibres& s = sphereFibres();
    int truncated = 0;
    for (const Fibre& f : s.fs.fibres) {
        const std::vector<Vec3> pts = fibrePoints(s.fs, f);
        CHECK(pts.size() >= 4);  // root, A0 = the curve start, >= 2 visible chords
        truncated += f.truncated;
        CHECK(f.emergeL > 0 && f.emergeL < norm(pts[1] - pts[0]));
        for (size_t k = 1; k + 1 < pts.size(); ++k) {
            const double turn = std::acos(std::clamp(dot(normalized(pts[k] - pts[k - 1]), normalized(pts[k + 1] - pts[k])), -1.0, 1.0)) * 180 / kPi;
            CHECK(turn >= s.fp.minTurnDeg - 1e-2 && turn <= s.fp.maxTurnDeg + 1e-2);  // measured on float32 points
        }
    }
    const StrandSetF st = fibresToStrands(s.fs, s.map);
    CHECK(checkStrands(st).ok);
    CHECK(st.strands() == s.fs.fibres.size());
    size_t v = 0;
    for (const Vec3f& p : s.fs.points) {
        const Vec3 q = s.map.toExportF(toVec3(p));
        CHECK(st.xyz[3 * v] == float(q.x) && st.xyz[3 * v + 1] == float(q.y) && st.xyz[3 * v + 2] == float(q.z));
        ++v;
    }
    std::printf("  %zu fibres (%d truncated), %zu vertices\n", s.fs.fibres.size(), truncated, v);
}

TEST(fiber_emergence_actual) {
    // E* lies on segment 1 and the visible length is measured from it.
    const SphereFibres& s = sphereFibres();
    for (const Fibre& f : s.fs.fibres) {
        const std::vector<Vec3> pts = fibrePoints(s.fs, f);
        const Vec3 E = toVec3(f.anchor);
        CHECK(closestSegmentPoints(E, E, pts[0], pts[1]).distance < 1e-5);  // float32 points
        CHECK_NEAR(f.emergeL, norm(E - pts[0]), 1e-5);
        CHECK_NEAR(visibleLength(s.fs, f), fibreLength(s.fs, f) - norm(E - pts[0]), 1e-5);
    }
}

// Circular arc of centreline radius R (control polygon on a circle), in the xy plane.
YarnPath arcYarn(double R, double degrees) {
    const int per = 36;
    const double Rc = R * 3 / (2 + std::cos(2 * kPi / per));  // uniform cubic B-spline shrinkage
    std::vector<Vec3> c;
    const int n = static_cast<int>(std::lround(degrees / 360 * per));
    for (int i = 0; i <= n; ++i) c.push_back({Rc * std::cos(2 * kPi * i / per), Rc * std::sin(2 * kPi * i / per), 0});
    return makeYarn(c, 0);
}

TEST(fibers_yarn_end_rejected) {
    // Near the yarn ends the span bracket leaves [s_margin, length - s_margin]: those
    // candidates are rejected, never clamped (ringAt would clamp silently).
    const double L = 1.5;
    CurveSet cs = makeSet({makeYarn(lineCtrl({0, 0, 0}, {L, 0, 0}, 2), 0)});
    PatternParams pp;
    SurfaceParams sp;
    sp.compression = false;
    const SurfaceModel sm = buildSurface(cs, pp, sp);
    FibreParams fp = pieceFibres();
    fp.randomFraction = 0;  // the yarn-end span bound is a DIRECTED rule
    fp.density = 400;
    fp.rootDepth = 0.126;
    const FibreSet fs = growYarnFibres(cs, sm, pp, fp, nullptr);
    printFibreReport(fs.rep);
    std::printf("  %d candidates, %d yarn end, %d accepted\n", fs.rep.candidates, fs.rep.yarnEnd, fs.rep.accepted);
    CHECK(fs.rep.yarnEnd > 0);
    for (const Fibre& f : fs.fibres)
        for (const Vec3& p : fibrePoints(fs, f)) CHECK(p.x >= 0.05 - 1e-6 && p.x <= L - 0.05 + 1e-6);
}

TEST(strand_budget_thinning) {
    // Whole fibres only, total <= budget, and the same kept set for any piece order
    // (as with --jobs 1 and --jobs 4).
    std::vector<std::pair<float, int>> ks;
    StrandSetF a, b;
    std::vector<float> ka, kb;
    long long total = 0;
    for (int i = 0; i < 5000; ++i) {
        const float key = float((i * 0.6180339887) - std::floor(i * 0.6180339887));
        const int segs = 2 + (i % 2);
        ks.push_back({key, segs});
        total += segs;
        std::vector<float> p, r;
        for (int k = 0; k <= segs; ++k) {
            p.insert(p.end(), {float(i), float(k), 0.f});
            r.push_back(1e-3f);
        }
        (i % 3 ? a : b).add(p, r);
        (i % 3 ? ka : kb).push_back(key);
    }
    const long long budget = total / 3;
    const float q = budgetThreshold(ks, budget);
    std::vector<std::pair<float, int>> rev(ks.rbegin(), ks.rend());
    CHECK(budgetThreshold(rev, budget) == q);
    const StrandSetF ta = keepStrands(a, ka, q), tb = keepStrands(b, kb, q);
    const long long kept = (ta.vertices() - ta.strands()) + (tb.vertices() - tb.strands());
    std::printf("  budget %lld: kept %lld segments (q %.4f)\n", budget, kept, double(q));
    CHECK(kept <= budget && kept > budget - 10);
    for (size_t i = 0; i < ta.strands(); ++i) {  // whole strands: every kept strand is 3 or 4 vertices of one fibre
        const uint32_t n = ta.first[i + 1] - ta.first[i];
        CHECK(n == 3 || n == 4);
        for (uint32_t k = ta.first[i]; k < ta.first[i + 1]; ++k) CHECK(ta.xyz[3 * k] == ta.xyz[3 * ta.first[i]]);
    }
    CHECK(budgetThreshold(ks, 0) == 2.0f);
}

// ---- S4: mixed fibre populations ----

FibreSet populationOf(const FibreSet& fs, FibrePopulation p) {
    FibreSet o;
    for (const Fibre& f : fs.fibres)
        if (f.population == p) o.addFibre(f, fs.pointsOf(f));
    return o;
}

// World turn at every interior VISIBLE vertex (beyond A0 = pts[1]) [deg].
std::vector<double> visibleTurns(const std::vector<Vec3>& pts) {
    std::vector<double> t;
    for (size_t k = 2; k + 1 < pts.size(); ++k)
        t.push_back(std::acos(std::clamp(dot(normalized(pts[k] - pts[k - 1]), normalized(pts[k + 1] - pts[k])), -1.0, 1.0)) * 180 / kPi);
    return t;
}

TEST(fibers_mix_endpoints) {
    // 0 % -> only Directed, 100 % -> only Random; the population is drawn before any
    // rejection, so the per-population counters split every report counter exactly.
    const StraightSet s = straightSet({0});
    for (double frac : {0.0, 1.0, 0.5}) {
        const FibreSet fs = growYarnFibres(s.cs, s.sm, s.pp, mixParams(frac, 30), nullptr);
        const FibreReport& r = fs.rep;
        CHECK(r.pop.size() == 2);
        for (int FibreReport::*m : fibreReportCounters()) CHECK(r.pop[0].*m + r.pop[1].*m == r.*m);
        int nD = 0, nR = 0;
        for (const Fibre& f : fs.fibres) (f.population == FibrePopulation::Random ? nR : nD) += 1;
        std::printf("  random %.1f: %d candidates (%d / %d), accepted directed %d, random %d\n", frac, r.candidates,
                    r.pop[0].candidates, r.pop[1].candidates, nD, nR);
        CHECK(nD == r.pop[0].accepted && nR == r.pop[1].accepted);
        if (frac == 0.0) CHECK(r.pop[1].candidates == 0 && nR == 0 && nD > 100);
        if (frac == 1.0) CHECK(r.pop[0].candidates == 0 && nD == 0 && nR > 100);
        if (frac == 0.5) CHECK(nD > 50 && nR > 50);
    }
}

TEST(fibers_mix_invariance) {
    // Changing the fraction never moves an anchor or changes a length, radius or key, and
    // the Random set only grows: draw level (all candidates) and generator level (every
    // candidate accepted in both runs).
    const double f1 = 0.3, f2 = 0.7;
    FibreParams a = pieceFibres(), b = pieceFibres();
    a.randomFraction = f1;
    b.randomFraction = f2;
    const int n = 100000;
    int r1 = 0, r2 = 0;
    for (int i = 0; i < n; ++i) {
        const FibreDraw da = drawFibre(a, 99, 0, i), db = drawFibre(b, 99, 0, i);
        CHECK(da.length == db.length && da.radiusFactor == db.radiusFactor && da.key == db.key && da.psi == db.psi &&
              da.lift == db.lift && da.azimuth == db.azimuth && da.elev == db.elev);
        const bool ra = da.population == FibrePopulation::Random, rb = db.population == FibrePopulation::Random;
        CHECK(!ra || rb);
        r1 += ra;
        r2 += rb;
    }
    std::printf("  random share %.4f (f1 %.1f), %.4f (f2 %.1f)\n", double(r1) / n, f1, double(r2) / n, f2);
    CHECK(std::fabs(double(r1) / n - f1) < 0.006 && std::fabs(double(r2) / n - f2) < 0.006);
    const StraightSet s = straightSet({0});
    const FibreSet g1 = growYarnFibres(s.cs, s.sm, s.pp, mixParams(f1, 30), nullptr);
    const FibreSet g2 = growYarnFibres(s.cs, s.sm, s.pp, mixParams(f2, 30), nullptr);
    CHECK(g1.rep.candidates == g2.rep.candidates);
    CHECK(g1.rep.pop[1].candidates <= g2.rep.pop[1].candidates);
    std::map<double, const Fibre*> by2;
    for (const Fibre& f : g2.fibres) by2[f.s] = &f;
    int common = 0;
    for (const Fibre& f : g1.fibres) {
        const auto it = by2.find(f.s);
        if (it == by2.end()) continue;
        ++common;
        const Fibre& h = *it->second;
        CHECK(h.theta == f.theta && h.targetLength == f.targetLength && h.radius0 == f.radius0 && h.key == f.key);
        if (f.population == FibrePopulation::Random) CHECK(h.population == FibrePopulation::Random);
    }
    std::printf("  %d candidates accepted in both runs\n", common);
    CHECK(common > 50);
}

// Curled strays (RandomCurve): off, nothing moves; on, the helix stays within its caps, the curve's
// derivative is right, the anchor and start tangent are those of the plain curve, and long fibres
// come out long, with their segments.
TEST(fibers_curled_strays) {
    FibreParams fp = pieceFibres();
    fp.randomFraction = 1;
    fp.lengthMin = 2;
    fp.lengthMean = 6;
    fp.lengthMax = 45;
    fp.lengthSpread = 1;
    fp.thickness = 0.05;
    fp.maxSegments = 48;
    FibreParams plain = fp;
    fp.randomCurlRadius = 1.5;
    fp.randomCurlLength = 9;
    const Vec3 A{1, 2, 3}, nrm = normalized(Vec3{0.3, -0.2, 1.0}), T = normalized(Vec3{1, 0.4, 0.1});
    double worstD = 0;
    for (int i = 0; i < 2000; ++i) {
        const FibreDraw d = drawFibre(fp, 4242, 0, i), o = drawFibre(plain, 4242, 0, i);
        CHECK(o.curlRadius == 0 && d.length == o.length && d.azimuth == o.azimuth && d.elev == o.elev && d.bendAmp == o.bendAmp);
        CHECK(d.curlRadius > 0 && d.curlRadius <= fp.randomCurlRadius && d.curlRadius <= 0.15 * d.length + 1e-12);
        CHECK(d.curlTurns > 0 && d.curlTurns <= 48.0 / kCurlChordsPerTurn);
        CHECK(2 * kPi * d.curlRadius * d.curlTurns <= 0.9 * d.length + 1e-9);
        RandomCurve c = randomCurve(d, A, nrm, T), b = randomCurve(o, A, nrm, T);
        c.S = b.S = 0.6 * d.length;
        CHECK(norm(b.X(0.5) - (b.A + (b.d0 * 0.5 + b.w * (b.a * 0.25)) * b.S)) == 0);  // the plain curve, unchanged
        CHECK(norm(c.X(0.0) - A) == 0 && norm(c.dX(0.0) - c.d0 * c.S) < 1e-12);
        for (double u = 0.03; u < 1; u += 0.07) {
            const double h = 1e-6;
            worstD = std::max(worstD, norm((c.X(u + h) - c.X(u - h)) / (2 * h) - c.dX(u)) / d.length);
        }
    }
    std::printf("  dX to central differences, relative to the length: %.1e\n", worstD);
    CHECK(worstD < 1e-6);
    const StraightSet s = straightSet({0});
    fp.density = 30;
    const FibreSet fs = growYarnFibres(s.cs, s.sm, s.pp, fp, nullptr);
    printFibreReport(fs.rep);
    double longest = 0;
    size_t most = 0;
    for (const Fibre& f : fs.fibres) {
        longest = std::max(longest, visibleLength(fs, f));
        most = std::max(most, static_cast<size_t>(f.pointCount));
        CHECK(f.pointCount <= 50 && visibleLength(fs, f) <= 1.1 * f.targetLength);
    }
    std::printf("  %zu curled strays, longest %.1f mm, at most %zu points\n", fs.fibres.size(), longest, most);
    CHECK(fs.fibres.size() > 100 && longest > 20 && most > 20);
    CHECK(fs.rep.pop[1].noBend == 0);  // a curled stray has no lower limit on its turns
}

// fibres.random.spread_deg: below 180 a stray starts within +-spread of one of the two directions
// along the yarn (seen from above the surface), both ways equally often; only the azimuth changes,
// and at 180 the draw and the params text are the original ones.
TEST(fiber_stray_spread) {
    FibreParams any = pieceFibres();
    any.randomFraction = 1;
    FibreParams spread = any;
    spread.randomSpreadDeg = 30;
    CHECK(fibreParamsText(any).find("spread") == std::string::npos);
    CHECK(fibreParamsText(spread) != fibreParamsText(any));
    const Vec3 nrm = normalized(Vec3{0.3, -0.2, 1.0}), T = normalized(Vec3{1, 0.4, 0.1});
    const Vec3 t1 = normalized(T - nrm * dot(T, nrm)), t2 = cross(nrm, t1);
    const int n = 20000;
    int forward = 0;
    double worst = 0, worstAny = 0;
    for (int i = 0; i < n; ++i) {
        const FibreDraw a = drawFibre(any, 4242, 0, i), d = drawFibre(spread, 4242, 0, i);
        CHECK(d.length == a.length && d.elev == a.elev && d.bendAmp == a.bendAmp && d.radiusFactor == a.radiusFactor);
        const RandomCurve c = randomCurve(d, {1, 2, 3}, nrm, T), ca = randomCurve(a, {1, 2, 3}, nrm, T);
        auto offYarn = [&](const Vec3& d0) {  // degrees between the start, seen from above, and the yarn line
            const double az = std::fabs(std::atan2(dot(d0, t2), dot(d0, t1))) * 180 / kPi;
            return std::min(az, 180 - az);
        };
        worst = std::max(worst, offYarn(c.d0));
        worstAny = std::max(worstAny, offYarn(ca.d0));
        if (dot(c.d0, t1) > 0) ++forward;
    }
    std::printf("  spread 30: largest angle to the yarn %.3f deg (spread 180: %.1f), forward share %.4f\n", worst, worstAny, double(forward) / n);
    CHECK(worst <= 30 + 1e-9 && worstAny > 85);
    CHECK(std::fabs(double(forward) / n - 0.5) < 0.02);
}

TEST(fiber_random_draws) {
    // Before any rejection: the start-tangent azimuth is uniform over the 8 sectors of
    // the local tangent plane, the elevation within its bounds, X'(0) is exactly the
    // start tangent d0 (the bow's derivative vanishes at the root), and lengths follow the
    // truncated log-normal (channel 6, shared with the directed fibres).
    FibreParams fp = pieceFibres();
    fp.randomFraction = 1;
    const int n = 100000;
    int sector[8] = {0};
    std::vector<double> L;
    const Vec3 nrm = normalized(Vec3{0.3, -0.2, 1.0}), T = normalized(Vec3{1, 0.4, 0.1});
    const Vec3 t1 = normalized(T - nrm * dot(T, nrm)), t2 = cross(nrm, t1);
    double worstTan = 0, worstElev = 0, worstAz = 0, worstW = 0;
    for (int i = 0; i < n; ++i) {
        const FibreDraw d = drawFibre(fp, 4242, 0, i);
        CHECK(d.population == FibrePopulation::Random);
        CHECK(d.azimuth >= 0 && d.azimuth < 2 * kPi);
        CHECK(d.elev >= fp.randomElevMinDeg * kPi / 180 - 1e-12 && d.elev <= fp.randomElevMaxDeg * kPi / 180 + 1e-12);
        CHECK(std::fabs(d.bendAmp) >= 0.5 * fp.randomBend - 1e-12 && std::fabs(d.bendAmp) <= fp.randomBend + 1e-12);
        ++sector[std::min(7, static_cast<int>(d.azimuth / (kPi / 4)))];
        L.push_back(d.length);
        if (i % 10) continue;
        RandomCurve c = randomCurve(d, {1, 2, 3}, nrm, T);
        c.S = 0.9;
        worstTan = std::max(worstTan, norm(normalized(c.dX(0.0)) - c.d0));
        worstElev = std::max(worstElev, std::fabs(std::asin(dot(c.d0, nrm)) - d.elev));
        double az = std::atan2(dot(c.d0, t2), dot(c.d0, t1));
        if (az < 0) az += 2 * kPi;
        worstAz = std::max(worstAz, std::fabs(std::remainder(az - d.azimuth, 2 * kPi)));
        worstW = std::max({worstW, std::fabs(dot(c.w, c.d0)), std::fabs(norm(c.w) - 1)});
        CHECK(norm(c.X(0.0) - Vec3{1, 2, 3}) < 1e-15);
    }
    std::printf("  sectors:");
    for (int k = 0; k < 8; ++k) std::printf(" %.4f", double(sector[k]) / n);
    std::printf("\n  |dX(0)/|dX(0)| - d0| <= %.1e, elevation %.1e, azimuth %.1e, w %.1e\n", worstTan, worstElev, worstAz, worstW);
    for (int k = 0; k < 8; ++k) CHECK(std::fabs(double(sector[k]) / n - 0.125) < 0.005);
    CHECK(worstTan < 1e-12 && worstElev < 1e-12 && worstAz < 1e-9 && worstW < 1e-12);
    const double mu = std::log(fp.lengthMean), sg = fp.lengthSpread;
    const double Fa = normalCdf((std::log(fp.lengthMin) - mu) / sg), Fb = normalCdf((std::log(fp.lengthMax) - mu) / sg);
    const double ks = ksUniform(L, [&](double x) { return (normalCdf((std::log(x) - mu) / sg) - Fa) / (Fb - Fa); });
    std::printf("  length KS %.4f\n", ks);
    CHECK(ks <= 0.01);
    // Generated on a straight yarn without quantization: the first visible chord rises at
    // least elevMin - atan(randomBend * 0.55) = 10.6 deg, so the shape alone never fails
    // the 5 deg emergence test.
    const StraightSet s = straightSet({0});
    const FibreSet fs = growYarnFibres(s.cs, s.sm, s.pp, mixParams(1.0, 30), nullptr);
    printFibreReport(fs.rep);
    double minEmerge = 90;
    for (const Fibre& f : fs.fibres) minEmerge = std::min(minEmerge, double(f.emergeAngle));
    std::printf("  %zu random fibres; min emergence angle %.2f deg; no-emergence rejections %d\n", fs.fibres.size(), minEmerge,
                fs.rep.pop[1].noEmergence);
    CHECK(fs.fibres.size() > 100);
    CHECK(fs.rep.pop[1].noEmergence == 0);
    CHECK(minEmerge >= fp.minEmergeDeg);
}

// Plies (or the single tube, count 1) of the capped sphere with a mixture, with the shared setup,
// the export's tube polygon and density per mm of the parent yarn.
struct PlyMix {
    PlySet ps;
    PatternParams pp;
    SurfaceModel sm;
    double ring = 0;
    int sides = 0;
    FibreParams fp;
    FibreSet fs;
};
std::unique_ptr<PlyMix> plyMix(int count, double frac, double density, const ExecPolicy& exec = {}) {
    const SphereFibres& s = sphereFibres();
    auto m = std::make_unique<PlyMix>();
    PlyGeometry g;
    g.count = count;
    m->ps = buildPlies(s.pk.curves, s.sm, g, s.pp.radius());
    m->pp = s.pp;
    if (count > 1) m->pp.yarnDiameter = 2 * m->ps.layout.rp;
    SurfaceParams sp;
    sp.compression = false;
    m->sm = buildSurface(m->ps.curves, m->pp, sp);
    m->ring = count > 1 ? kExportPlyRing : kExportRing;
    m->sides = count > 1 ? kExportPlySides : kExportSides;
    FibreSetup set;
    set.tubeRadius = m->pp.radius();
    set.sides = m->sides;
    set.ring = m->ring;
    const Allowances al = Allowances::forExtent(m->ps.curves.extent());
    for (const YarnPath& y : m->ps.curves.yarns)
        set.kappaMax = std::max(set.kappaMax, std::min(1e3, checkCurvature(y, m->pp.radius(), al).kappaHi));
    if (count > 1) {
        set.plies = &m->ps;
        set.parentCurves = &s.pk.curves;
        set.parentSurface = &s.sm;
    }
    m->fp = s.fp;
    m->fp.randomFraction = frac;
    m->fp.density = density;
    m->fp.densityPerParent = true;
    m->fp.exportChordFactor = m->ring * m->ring / 8;
    FibreFrame fr = s.fr;
    setupPieceFibres(m->fp, fr, s.pp, set);
    fr.exec = exec;
    m->fs = growYarnFibres(m->ps.curves, m->sm, m->pp, m->fp, &fr);
    return m;
}

TEST(fibers_mix_on_plies) {
    // Single tube and 2, 3, 6 plies at 50 %: per population a NON-TRIVIAL accepted count.
    for (int count : {1, 2, 3, 6}) {
        const std::unique_ptr<PlyMix> m = plyMix(count, 0.5, 5.0);
        std::printf(" %d plies (r_p %.3f, root depth %.3f mm):\n", count, m->pp.radius(), m->fp.rootDepth);
        printFibreReport(m->fs.rep);
        for (FibrePopulation p : {FibrePopulation::Directed, FibrePopulation::Random}) {
            const FibreSet part = populationOf(m->fs, p);
            std::printf("   %-8s %5zu accepted of %6d candidates\n", toString(p), part.fibres.size(),
                        m->fs.rep.pop[static_cast<int>(p)].candidates);
            // About half the smallest measured count (6 plies, directed: 497 at 5 / mm);
            // never lowered without recording why.
            const int nMin = 240;
            CHECK(static_cast<int>(part.fibres.size()) >= nMin);
        }
    }
}

TEST(fibers_random_shape) {
    // Accepted random fibres on the 3-ply sphere: every start-tangent sector of the local
    // tangent plane keeps >= 1/32 of them (a quarter of the uniform share), visible
    // lengths vary, and every one bends at an interior visible vertex within the turn limits.
    const std::unique_ptr<PlyMix> m = plyMix(3, 0.5, 5.0);
    const FibreSet r = populationOf(m->fs, FibrePopulation::Random);
    int sector[8] = {0};
    std::vector<double> vis;
    for (const Fibre& f : r.fibres) {
        double az = f.psi;  // Random: azimuth [deg]
        ++sector[std::clamp(static_cast<int>(az / 45.0), 0, 7)];
        vis.push_back(visibleLength(r, f));
        CHECK(f.pointCount >= 4);
        const std::vector<double> t = visibleTurns(fibrePoints(r, f));
        CHECK(!t.empty());
        for (double a : t) CHECK(a >= std::max(1.0, m->fp.minTurnDeg) - 1e-2 && a <= m->fp.maxTurnDeg + 1e-2);  // float32 points
    }
    std::sort(vis.begin(), vis.end());
    const double iqr = vis[3 * vis.size() / 4] - vis[vis.size() / 4];
    std::printf("  %zu random fibres; sectors:", r.fibres.size());
    for (int k = 0; k < 8; ++k) std::printf(" %.3f", double(sector[k]) / r.fibres.size());
    std::printf("; visible length IQR %.3f mm (median %.3f)\n", iqr, vis[vis.size() / 2]);
    CHECK(r.fibres.size() > 200);
    for (int k = 0; k < 8; ++k) CHECK(sector[k] * 32 >= static_cast<int>(r.fibres.size()));
    CHECK(iqr > 0.15);
}

TEST(fibers_random_truncation) {
    // Two parallel straight yarns 0.3 mm apart: random fibres growing sideways hit the neighbour
    // and are cut. Some truncated random fibres must actually occur, and each keeps >= 2 visible
    // segments with an interior turn within the limits.
    const StraightSet s = straightSet({0, 1.3});
    FibreParams fp = mixParams(1.0, 150);
    const FibreSet fs = growYarnFibres(s.cs, s.sm, s.pp, fp);
    printFibreReport(fs.rep);
    int truncated = 0;
    for (const Fibre& f : fs.fibres) {
        if (!f.truncated) continue;
        ++truncated;
        CHECK(f.population == FibrePopulation::Random);
        CHECK(f.pointCount >= 4);
        const std::vector<double> t = visibleTurns(fibrePoints(fs, f));
        CHECK(!t.empty());
        for (double a : t) CHECK(a >= std::max(1.0, fp.minTurnDeg) - 1e-2 && a <= fp.maxTurnDeg + 1e-2);  // float32 points
    }
    std::printf("  %zu random fibres, %d truncated\n", fs.fibres.size(), truncated);
    CHECK(truncated >= 1);
}

TEST(fibres_cut_before_deep_contact) {
    // The light rule, measured without the generator's own index: on straight round yarns the depth
    // of a point inside a strand is r minus its distance to that yarn's axis. After the curve start
    // no vertex or segment midpoint of any fibre lies deeper than 20 % of r in any strand, and the
    // cut really happens.
    const StraightSet s = straightSet({0, 1.3});
    const FibreParams fp = mixParams(0.5, 150);
    const FibreSet fs = growYarnFibres(s.cs, s.sm, s.pp, fp);
    const double r = s.pp.radius();
    double deepest = 0;
    int truncated = 0;
    for (const Fibre& f : fs.fibres) {
        truncated += f.truncated;
        const std::vector<Vec3> pts = fibrePoints(fs, f);
        for (size_t k = 2; k < pts.size(); ++k)
            for (const Vec3& q : {lerp(pts[k - 1], pts[k], 0.5), pts[k]})
                for (double y : {0.0, 1.3}) deepest = std::max(deepest, r - std::hypot(q.y - y, q.z));
    }
    std::printf("  %zu fibres, %d truncated, deepest %.4f mm (limit %.4f)\n", fs.fibres.size(), truncated, deepest, 0.2 * r);
    CHECK(truncated >= 1);
    CHECK(deepest <= 0.2 * r + 1e-3);  // the capsule-axis chord of the index
}

TEST(fibers_density_per_parent) {
    // Density per mm of the parent yarn: every yarn's plies share floor(density L_y)
    // candidates exactly, for any pitch and ply count; per mm of each ply (S3) the total
    // grows with the twist.
    const StraightSet s = straightSet({0, 3}, {12, 17});
    const double density = 2.0;
    int expect = 0;
    for (const YarnPath& y : s.cs.yarns) expect += static_cast<int>(std::floor(density * y.length()));
    for (double pitch : {1.5, 2.5, 5.0})
        for (int count : {2, 3, 6}) {
            PlyGeometry g;
            g.count = count;
            g.pitch = pitch;
            const PlySet ps = buildPlies(s.cs, s.sm, g, s.pp.radius());
            const std::vector<int> budget = plyCandidateBudget(ps, s.cs, density);
            std::vector<int> perYarn(s.cs.yarns.size(), 0);
            for (size_t j = 0; j < budget.size(); ++j) perYarn[ps.yarnOf[j]] += budget[j];
            for (size_t y = 0; y < perYarn.size(); ++y) CHECK(perYarn[y] == static_cast<int>(std::floor(density * s.cs.yarns[y].length())));
            PatternParams pp = s.pp;
            pp.yarnDiameter = 2 * ps.layout.rp;
            SurfaceParams sp;
            sp.compression = false;
            const SurfaceModel sm = buildSurface(ps.curves, pp, sp);
            FibreParams fp = mixParams(0.5, density);
            FibreFrame fr;
            FibreSetup set;
            set.tubeRadius = pp.radius();
            set.plies = &ps;
            set.parentCurves = &s.cs;
            set.parentSurface = &s.sm;
            setupPieceFibres(fp, fr, s.pp, set);
            fp.rootDepth = 0.05;
            const FibreSet a = growYarnFibres(ps.curves, sm, pp, fp, &fr);
            fr.candidates.clear();  // per mm of each ply
            const FibreSet b = growYarnFibres(ps.curves, sm, pp, fp, &fr);
            int perPly = 0;
            for (const YarnSurface& ys : sm.yarns) perPly += static_cast<int>(std::floor(density * ys.length));
            std::printf("  pitch %.1f, %d plies: %d candidates per yarn length (expected %d), %d per ply length\n", pitch, count,
                        a.rep.candidates, expect, b.rep.candidates);
            CHECK(a.rep.candidates == expect);
            CHECK(b.rep.candidates == perPly);
            for (const Fibre& f : a.fibres) CHECK(f.s >= 0 && f.s <= sm.yarns[f.yarnId].length);
        }
}

TEST(fibers_parent_yarn_on_plies) {
    // On plies every fibre carries its PARENT yarn (viewer colour and selection index the yarn
    // array with it), in both populations.
    const StraightSet s = straightSet({0, 3});
    PlyGeometry g;
    const PlySet ps = buildPlies(s.cs, s.sm, g, s.pp.radius());
    PatternParams pp = s.pp;
    pp.yarnDiameter = 2 * ps.layout.rp;
    SurfaceParams sp;
    sp.compression = false;
    const SurfaceModel sm = buildSurface(ps.curves, pp, sp);
    FibreParams fp = mixParams(0.5, 20);
    FibreFrame fr;
    FibreSetup set;
    set.tubeRadius = pp.radius();
    set.plies = &ps;
    set.parentCurves = &s.cs;
    set.parentSurface = &s.sm;
    setupPieceFibres(fp, fr, s.pp, set);
    const FibreSet fs = growYarnFibres(ps.curves, sm, pp, fp, &fr);
    int pops[2] = {0, 0}, parents[2] = {0, 0};
    for (const Fibre& f : fs.fibres) {
        CHECK(f.parentYarn == ps.yarnOf[f.yarnId]);
        ++pops[static_cast<int>(f.population)];
        ++parents[f.parentYarn];
    }
    std::printf("  %zu fibres (directed %d, random %d) on yarns %d / %d\n", fs.fibres.size(), pops[0], pops[1], parents[0], parents[1]);
    CHECK(parents[0] > 10 && parents[1] > 10);
    CHECK(pops[0] > 10 && pops[1] > 10);
}

// ---- Complete snapshots (fibre_snapshot.h): the reference for the performance work ----
// Four configurations through the current generator, compared field by field (every Fibre and
// FibreReport field, hex floats; timing excluded). Persisted as a digest of the whole snapshot plus
// its first 300 lines, scoped to the recorded build configuration.

std::vector<std::pair<std::string, FibreSet>> snapshotSets(const ExecPolicy& exec = {}) {
    std::vector<std::pair<std::string, FibreSet>> v;
    const SphereFibres& s = sphereFibres();
    FibreFrame sphereFrame = s.fr;
    sphereFrame.exec = exec;
    v.push_back({"sphere_tube", growYarnFibres(s.pk.curves, s.sm, s.pp, s.fp, &sphereFrame)});
    {
        FibreFixture f;
        FibreParams fp = FibreFixture::params();
        fp.randomFraction = 0;
        v.push_back({"patch_along", f.run(0.005, fp, exec)});
    }
    v.push_back({"sphere_plies3_mix", std::move(plyMix(3, 0.5, 5.0, exec)->fs)});
    FibreFrame execOnly;  // only the execution policy: behaves as no frame
    execOnly.exec = exec;
    const StraightSet straight = straightSet({0, 1.3});
    v.push_back({"straight_random", growYarnFibres(straight.cs, straight.sm, straight.pp, mixParams(1.0, 60), &execOnly)});
    return v;
}

std::string snapshotPath(const std::string& name) { return std::string(KNIT_TEST_DATA) + "/fibers_snapshot_" + name + ".txt"; }

// Persisted form: digest and line count of the whole snapshot, then its first 300 lines.
std::string snapshotFixtureText(const std::string& snap) {
    size_t lines = 0, cut = 0;
    for (size_t i = 0; i < snap.size(); ++i)
        if (snap[i] == '\n' && ++lines == 300) cut = i + 1;
    char b[96];
    std::snprintf(b, sizeof b, "digest %016llx lines %zu\n", static_cast<unsigned long long>(fnv1a64(snap)), lines);
    return b + (cut ? snap.substr(0, cut) : snap);
}

// Not registered with CTest: writes missing snapshot fixtures (create-only).
TEST(fibers_write_snapshots) {
    for (const auto& [name, fs] : snapshotSets()) {
        if (std::filesystem::exists(snapshotPath(name))) {
            std::printf("  %s: exists, not rewritten\n", name.c_str());
            continue;
        }
        std::ofstream f(snapshotPath(name), std::ios::binary);
        f << "# build: " << buildConfiguration() << "\n" << snapshotFixtureText(fibreSnapshot(fs));
        std::printf("  %s: %zu fibres written\n", name.c_str(), fs.fibres.size());
    }
}

TEST(fibers_snapshots) {
    // Every configuration reproduces its recorded snapshot exactly (digest of the whole
    // snapshot and the first 300 lines). Scope: the recorded build configuration; under
    // another configuration a mismatch is UNRESOLVED, not a failure.
    for (const auto& [name, fs] : snapshotSets()) {
        std::ifstream f(snapshotPath(name), std::ios::binary);
        CHECK(f.good());
        std::string header;
        std::getline(f, header);
        const std::string want((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
        const std::string got = snapshotFixtureText(fibreSnapshot(fs));
        if (got == want) {
            std::printf("  %s: %zu fibres identical\n", name.c_str(), fs.fibres.size());
            continue;
        }
        size_t i = 0;
        while (i < got.size() && i < want.size() && got[i] == want[i]) ++i;
        std::printf("  %s differs at byte %zu:\n    want %.160s\n    got  %.160s\n", name.c_str(), i,
                    want.c_str() + std::min(i, want.size()), got.c_str() + std::min(i, got.size()));
        if (header != "# build: " + buildConfiguration()) {
            std::printf("  UNRESOLVED: recorded '%s', this build '%s'\n", header.c_str(), buildConfiguration().c_str());
            continue;
        }
        CHECK(got == want);
    }
}

// ---- SpatialHash3: canonical-cell visits (fibre perf A2) ----

// The previous implementation (a per-query visit stamp), kept as the reference order.
class StampHashRef {
public:
    explicit StampHashRef(double cell) : cell_(cell) {}
    void insert(const Vec3& lo, const Vec3& hi) {
        const int id = static_cast<int>(lo_.size());
        lo_.push_back(lo);
        hi_.push_back(hi);
        stamp_.push_back(0);
        forCells(lo, hi, [&](std::int64_t key) { cells_[key].push_back(id); });
    }
    std::vector<int> query(const Vec3& lo, const Vec3& hi) {
        std::vector<int> v;
        ++epoch_;
        forCells(lo, hi, [&](std::int64_t key) {
            auto it = cells_.find(key);
            if (it == cells_.end()) return;
            for (int id : it->second) {
                if (stamp_[id] == epoch_) continue;
                stamp_[id] = epoch_;
                const Vec3 &a = lo_[id], &b = hi_[id];
                if (a.x > hi.x || b.x < lo.x || a.y > hi.y || b.y < lo.y || a.z > hi.z || b.z < lo.z) continue;
                v.push_back(id);
            }
        });
        return v;
    }

private:
    double cell_;
    std::unordered_map<std::int64_t, std::vector<int>> cells_;
    std::vector<Vec3> lo_, hi_;
    std::vector<std::uint64_t> stamp_;
    std::uint64_t epoch_ = 0;
    static std::int64_t key(std::int64_t i, std::int64_t j, std::int64_t k) {
        const std::int64_t o = 1 << 20;
        return ((i + o) & 0x1FFFFF) | (((j + o) & 0x1FFFFF) << 21) | (((k + o) & 0x1FFFFF) << 42);
    }
    template <class F>
    void forCells(const Vec3& lo, const Vec3& hi, F&& f) const {
        const auto c = [&](double v) { return static_cast<std::int64_t>(std::floor(v / cell_)); };
        for (std::int64_t i = c(lo.x); i <= c(hi.x); ++i)
            for (std::int64_t j = c(lo.y); j <= c(hi.y); ++j)
                for (std::int64_t k = c(lo.z); k <= c(hi.z); ++k) f(key(i, j, k));
    }
};

// Random boxes: half of them snapped to cell boundaries (faces, edges, corners), around
// negative and positive coordinates, from point-sized to several cells.
struct HashCase {
    std::vector<std::pair<Vec3, Vec3>> items, queries;
};
HashCase hashCase(double cell, uint64_t seed) {
    HashCase c;
    uint64_t st = seed;
    auto u = [&] {
        st = st * 6364136223846793005ull + 1442695040888963407ull;
        return double(st >> 11) * (1.0 / 9007199254740992.0);
    };
    auto box = [&] {
        Vec3 lo{(u() - 0.5) * 12 * cell, (u() - 0.5) * 12 * cell, (u() - 0.5) * 12 * cell};
        Vec3 ext{u() * 2.5 * cell, u() * 2.5 * cell, u() * 2.5 * cell};
        const double kind = u();
        if (kind < 0.15) ext = Vec3{0, 0, 0};  // a point
        if (kind > 0.5) {                        // snap the corners to cell boundaries
            auto snap = [&](double v) { return std::round(v / cell) * cell; };
            lo = {snap(lo.x), snap(lo.y), snap(lo.z)};
            if (kind > 0.75) ext = {snap(ext.x), snap(ext.y), snap(ext.z)};
        }
        return std::make_pair(lo, lo + ext);
    };
    for (int i = 0; i < 400; ++i) c.items.push_back(box());
    for (int i = 0; i < 400; ++i) c.queries.push_back(box());
    return c;
}

TEST(topology_stitch_near) {
    // Two stitches on yarn 0 with a gap between them; yarn 1 has none.
    KnitTopology topo;
    StitchRecord a, b;
    a.id = 0;
    a.uBegin = 2;
    a.uEnd = 3;
    b.id = 1;
    b.uBegin = 5;
    b.uEnd = 6;
    topo.stitches = {a, b};
    topo.buildIndex({8, 8});
    CHECK(topo.stitchNear(0, 0) == 0);    // before the first stitch: the first
    CHECK(topo.stitchNear(0, 2) == 0);
    CHECK(topo.stitchNear(0, 4) == 0);    // in the gap: the last one started
    CHECK(topo.stitchAt(0, 4) == -1);
    CHECK(topo.stitchNear(0, 5) == 1);
    CHECK(topo.stitchNear(0, 7) == 1);    // on the tail
    CHECK(topo.stitchNear(1, 4) == -1);   // a yarn without stitches
    CHECK(topo.stitchNear(2, 4) == -1);   // no such yarn
}

TEST(hash_canonical_cell_order) {
    // Every query visits exactly the items, in exactly the order, of the stamp version.
    int queries = 0;
    long long visits = 0;
    for (double cell : {1.0, 0.37, 2.5})
        for (uint64_t seed = 1; seed <= 20; ++seed) {
            const HashCase hc = hashCase(cell, seed);
            SpatialHash3 h(cell);
            StampHashRef ref(cell);
            for (const auto& [lo, hi] : hc.items) {
                h.insert(lo, hi);
                ref.insert(lo, hi);
            }
            // Every lookup mode: map (occupied range clipped), dense after finalize(), and a copy
            // of the finalized hash.
            auto compare = [&](const SpatialHash3& hh) {
                for (const auto& [lo, hi] : hc.queries) {
                    std::vector<int> got;
                    hh.query(lo, hi, [&](int id) { got.push_back(id); });
                    CHECK(got == ref.query(lo, hi));
                    ++queries;
                    visits += static_cast<long long>(got.size());
                }
            };
            compare(h);
            h.finalize();
            compare(h);
            const SpatialHash3 copy = h;
            compare(copy);
        }
    std::printf("  %d queries, %lld visits: identical to the stamp order (map, dense, copy)\n",
                queries, visits);
    CHECK(visits > 10000);
}

TEST(hash_concurrent_queries) {
    // 8 threads querying one hash at once reproduce the serial visit sequences.
    const HashCase hc = hashCase(0.5, 77);
    SpatialHash3 h(0.5);
    for (const auto& [lo, hi] : hc.items) h.insert(lo, hi);
    std::vector<std::vector<int>> serial;
    for (const auto& [lo, hi] : hc.queries) {
        serial.emplace_back();
        h.query(lo, hi, [&](int id) { serial.back().push_back(id); });
    }
    std::vector<int> mismatches(8, 0);
    std::vector<std::thread> pool;
    for (int t = 0; t < 8; ++t)
        pool.emplace_back([&, t] {
            for (int rep = 0; rep < 50; ++rep)
                for (size_t q = t; q < hc.queries.size(); q += 3) {
                    std::vector<int> got;
                    h.query(hc.queries[q].first, hc.queries[q].second, [&](int id) { got.push_back(id); });
                    mismatches[t] += got != serial[q];
                }
        });
    for (std::thread& th : pool) th.join();
    int total = 0;
    for (int m : mismatches) total += m;
    std::printf("  8 threads x 50 rounds: %d mismatching queries\n", total);
    CHECK(total == 0);
}

TEST(mesh_closest_skip_revisits) {
    // MeshQuery::closest without re-evaluating earlier rounds' triangles returns exactly the
    // distance of a loop over every triangle, and its hit is the closest point of that hit's
    // triangle: random points
    // in and around the piece (near, inside, far outside the bounding box) and points exactly
    // on vertices and edge midpoints (ties between triangles sharing them).
    FixturePiece f = fixturePiece(25.0);
    MeshQuery mq(f.mesh);
    Vec3 lo, hi;
    f.mesh.bounds(lo, hi);
    const Vec3 ext = hi - lo;
    std::vector<Vec3> pts;
    for (int i = 0; i < 6000; ++i) {
        const double a = std::fmod(0.5 + i * 0.6180339887, 1.0), b = std::fmod(0.5 + i * 0.7548776662, 1.0),
                     c = std::fmod(0.5 + i * 0.5698402910, 1.0);
        const double spread = i % 3 == 0 ? 3.0 : 1.4;
        pts.push_back(lo + ext * 0.5 - ext * (0.5 * spread) + Vec3{ext.x * spread * a, ext.y * spread * b, ext.z * spread * c});
    }
    for (size_t v = 0; v < f.mesh.pos.size(); v += 7) pts.push_back(f.mesh.pos[v]);
    for (size_t t = 0; t < f.mesh.tri.size(); t += 11)
        pts.push_back((f.mesh.pos[f.mesh.tri[t][0]] + f.mesh.pos[f.mesh.tri[t][1]]) * 0.5);
    int insideCount = 0;
    for (const Vec3& q : pts) {
        const MeshQuery::Hit a = mq.closest(q);
        double nearest = 1e300;
        for (const auto& t : f.mesh.tri)
            nearest =
                std::min(nearest, norm(q - closestOnTriangle(q, f.mesh.pos[t[0]], f.mesh.pos[t[1]], f.mesh.pos[t[2]]).p));
        CHECK(a.dist == nearest);
        const auto& t = f.mesh.tri[a.tri];
        const Vec3 onHit = closestOnTriangle(q, f.mesh.pos[t[0]], f.mesh.pos[t[1]], f.mesh.pos[t[2]]).p;
        CHECK(onHit.x == a.p.x && onHit.y == a.p.y && onHit.z == a.p.z);
        insideCount += a.inside;
    }
    std::printf("  %zu points (%d inside): identical hits\n", pts.size(), insideCount);
    CHECK(insideCount > 500);
}

// ---- Executor (fibre perf, Stage B) ----

// Threads currently inside a chunk (counted once per thread, nested chunks included).
struct ThreadGauge {
    std::atomic<int> now{0}, peak{0};
    void enter(int& depth) {
        if (depth++ == 0) {
            const int v = ++now;
            int p = peak.load();
            while (v > p && !peak.compare_exchange_weak(p, v)) {}
        }
    }
    void leave(int& depth) {
        if (--depth == 0) --now;
    }
};
thread_local int tGaugeDepth = 0;

void spin(int us) {
    const auto t = std::chrono::steady_clock::now();
    while (std::chrono::steady_clock::now() - t < std::chrono::microseconds(us)) {}
}

TEST(executor_basics) {
    // Zero work, uneven chunks (every index exactly once, chunk indices in range), serial
    // policy in chunk order on the calling thread.
    ExecPolicy par;
    par.parallel = true;
    int calls = 0;
    parallelFor(0, 8, [&](int64_t, int64_t, int64_t) { ++calls; }, par);
    CHECK(calls == 0);
    for (int64_t n : {1, 7, 1000, 1001}) {
        std::vector<std::atomic<int>> hit(static_cast<size_t>(n));
        parallelFor(n, 64, [&](int64_t b, int64_t e, int64_t k) {
            CHECK(b == k * 64 && e == std::min<int64_t>(n, b + 64));
            for (int64_t i = b; i < e; ++i) ++hit[static_cast<size_t>(i)];
        }, par);
        for (auto& h : hit) CHECK(h.load() == 1);
    }
    std::vector<int64_t> order;
    const auto self = std::this_thread::get_id();
    parallelFor(10, 3, [&](int64_t, int64_t, int64_t k) {
        CHECK(std::this_thread::get_id() == self);
        order.push_back(k);
    }, ExecPolicy{});
    CHECK((order == std::vector<int64_t>{0, 1, 2, 3}));
}

TEST(executor_budget_and_nesting) {
    // Several external submitters, each with nested ranges: never more than H threads
    // inside chunks at once; one range alone reaches all H (the one-piece tail); the per-
    // range cap holds.
    const int H = threadBudget();
    ExecPolicy par;
    par.parallel = true;
    ThreadGauge g;
    std::atomic<long long> work{0};
    auto inner = [&](int64_t b, int64_t e, int64_t) {
        g.enter(tGaugeDepth);
        for (int64_t i = b; i < e; ++i) work += 1;
        spin(200);
        g.leave(tGaugeDepth);
    };
    auto outer = [&](int64_t, int64_t, int64_t) {
        g.enter(tGaugeDepth);
        parallelFor(64, 4, inner, par);  // nested
        g.leave(tGaugeDepth);
    };
    std::vector<std::thread> submitters;
    for (int s = 0; s < 4; ++s) submitters.emplace_back([&] { parallelFor(12, 1, outer, par); });
    for (std::thread& t : submitters) t.join();
    std::printf("  budget H = %d; 4 submitters x nested ranges: peak %d threads in chunks\n", H, g.peak.load());
    CHECK(work.load() == 4LL * 12 * 64);
    CHECK(g.peak.load() <= H);
    // One range alone can use all H slots: H chunks that each wait until all H have
    // arrived (a barrier) complete only if H run at once (10 s timeout instead of a hang).
    std::atomic<int> arrived{0};
    std::atomic<bool> timedOut{false};
    parallelFor(H, 1, [&](int64_t, int64_t, int64_t) {
        ++arrived;
        const auto t = std::chrono::steady_clock::now();
        while (arrived.load() < H)
            if (std::chrono::steady_clock::now() - t > std::chrono::seconds(10)) {
                timedOut = true;
                break;
            }
    }, par);
    std::printf("  one range: %d of %d chunks ran at once%s\n", arrived.load(), H, timedOut ? " (TIMEOUT)" : "");
    CHECK(!timedOut.load());
    // The per-range cap.
    ThreadGauge capped;
    ExecPolicy cap = par;
    cap.maxConcurrency = 3;
    parallelFor(48, 1, [&](int64_t, int64_t, int64_t) {
        capped.enter(tGaugeDepth);
        spin(1000);
        capped.leave(tGaugeDepth);
    }, cap);
    CHECK(capped.peak.load() <= 3);
}

TEST(executor_exceptions) {
    // The lowest throwing chunk's exception is rethrown on the caller and the executor stays
    // usable. Cancellation of chunks not yet started is checked deterministically with a
    // concurrency cap of 1 (chunks start strictly in order: a throw in chunk 5 leaves exactly
    // 6 run); without a cap, how many started before the throw depends on scheduling.
    ExecPolicy par;
    par.parallel = true;
    std::atomic<int> ran{0};
    std::string what;
    try {
        parallelFor(200, 1, [&](int64_t, int64_t, int64_t k) {
            ++ran;
            spin(100);
            if (k == 5 || k == 9) throw std::runtime_error("chunk " + std::to_string(k));
        }, par);
    } catch (const std::runtime_error& e) {
        what = e.what();
    }
    std::printf("  uncapped: rethrown '%s'; %d of 200 chunks ran\n", what.c_str(), ran.load());
    CHECK(what == "chunk 5");  // chunk 5 always starts before chunk 9
    ExecPolicy one = par;
    one.maxConcurrency = 1;
    std::atomic<int> ranOne{0};
    what.clear();
    try {
        parallelFor(200, 1, [&](int64_t, int64_t, int64_t k) {
            ++ranOne;
            if (k == 5) throw std::runtime_error("chunk 5");
        }, one);
    } catch (const std::runtime_error& e) {
        what = e.what();
    }
    std::printf("  cap 1: rethrown '%s'; %d of 200 chunks ran\n", what.c_str(), ranOne.load());
    CHECK(what == "chunk 5" && ranOne.load() == 6);
    std::atomic<int> after{0};
    parallelFor(100, 1, [&](int64_t, int64_t, int64_t) { ++after; }, par);
    CHECK(after.load() == 100);
}

TEST(fibers_parallel_identical) {
    // All snapshot configurations generated in parallel (at most 1, 2, 8 and all H chunks
    // at once, small chunks): the complete snapshot strings equal the serial ones exactly.
    const auto serial = snapshotSets();
    std::vector<std::string> want;
    for (const auto& [name, fs] : serial) want.push_back(fibreSnapshot(fs));
    for (int conc : {1, 2, 8, 0}) {
        ExecPolicy par;
        par.parallel = true;
        par.maxConcurrency = conc;
        par.chunk = 97;
        const auto got = snapshotSets(par);
        CHECK(got.size() == serial.size());
        for (size_t k = 0; k < got.size(); ++k) {
            const bool same = fibreSnapshot(got[k].second) == want[k];
            if (!same) std::printf("  %s differs at concurrency %d\n", got[k].first.c_str(), conc);
            CHECK(same);
        }
        std::printf("  concurrency %s: %zu configurations identical\n", conc ? std::to_string(conc).c_str() : "H", got.size());
    }
}

TEST(fibers_serial_policy_single_thread) {
    // The default (and a serial) policy calls every callback from the calling thread only,
    // in generation and in each check.
    const SphereFibres& s = sphereFibres();
    std::mutex m;
    std::set<std::thread::id> threads;
    auto note = [&] {
        std::lock_guard<std::mutex> lk(m);
        threads.insert(std::this_thread::get_id());
    };
    FibreFrame fr = s.fr;
    const auto fabric = fr.fabric;
    fr.fabric = [&](const Vec3& p) {
        note();
        return fabric(p);
    };
    FibreParams fp = s.fp;
    fp.density = 2;
    const FibreSet fs = growYarnFibres(s.pk.curves, s.sm, s.pp, fp, &fr);
    CHECK(!fs.fibres.empty());
    const SceneMapping map = s.map;
    const ExportTubeIndex tubes(s.pk.curves, s.pp, ExportTubeSpec{kExportModelRing, kExportRing, kExportSides, false,
                                [&](const Vec3& p) {
                                    note();
                                    return map.toExportF(p);
                                },
                                s.eSer});
    rootAttachmentCheck(tubes, fs, 1e-4 / map.sigma);
    std::printf("  callbacks ran on %zu thread(s)\n", threads.size());
    CHECK(threads.size() == 1 && *threads.begin() == std::this_thread::get_id());
}

TEST(checks_parallel_identical) {
    // The reference check and the attachment check in parallel at 1, 2, 8 and H chunks at once:
    // complete result snapshots (counters, index lists, minima / maxima) equal the serial ones.
    const std::unique_ptr<PlyMix> m = plyMix(3, 0.5, 5.0);
    const SphereFibres& s = sphereFibres();
    const SceneMapping map = s.map;
    auto toExp = [map](const Vec3& p) { return map.toExportF(p); };
    const double modelRing = kExportModelRing;
    auto run = [&](const ExecPolicy& exec) {
        const ExportTubeIndex tubes(m->ps.curves, m->pp, ExportTubeSpec{modelRing, m->ring, m->sides, false, toExp, s.eSer}, exec);
        std::string o = checkSnapshot(fibreReferenceCheck(tubes, m->fs, 1e-4, exec)) + "\n";
        return o + checkSnapshot(rootAttachmentCheck(tubes, m->fs, 1e-4 / map.sigma, exec));
    };
    const std::string serial = run(ExecPolicy{});
    for (int conc : {1, 2, 8, 0}) {
        ExecPolicy par;
        par.parallel = true;
        par.maxConcurrency = conc;
        par.chunk = 37;
        CHECK(run(par) == serial);
    }
    std::printf("  %zu fibres: identical at 1, 2, 8 and H\n", m->fs.fibres.size());
}

// ---- recompute 2c: parallel validation inside a part (identical results) ----

// Every result field of a validation report (doubles as exact hex floats; timings excluded).
std::string reportSnapshot(const ValidationReport& r) {
    std::string o;
    char b[512];
    auto put = [&](const char* fmt, auto... a) {
        std::snprintf(b, sizeof b, fmt, a...);
        o += b;
    };
    put("verdict %s\n", toString(r.verdict));
    for (const CheckLine& l : r.lines) put("line %s %s %s\n", l.name.c_str(), toString(l.status), l.detail.c_str());
    for (const CurvatureResult& c : r.curvature) {
        put("curv %d %a %a %a %a %a %a %a %s %s %d\n", c.yarnId, c.speedLo, c.kappaHi, c.rhoLo, c.rhoEst, c.uRhoEst,
            c.kappaWitnessLo, c.uWitness, toString(c.regular), toString(c.noFold), c.leaves);
        for (double k : c.spanKappaHi) put(" %a", k);
        o += "\n";
    }
    for (const Polyline& p : r.polylines) put("poly %zu %a\n", p.p.size(), p.epsBound);
    const ClearanceResult& c = r.clearance;
    put("clr %s %a %a %a %ld %ld %zu %zu\n", toString(c.status), c.minInterYarn, c.minSelfNonlocal, c.deltaD, c.segmentPairs,
        c.schurCertified, c.issues.size(), c.regions.size());
    for (const ClearanceIssue& is : c.issues)
        put("is %d %d %a %a %a %a %a %a %a %a %a %a %s %d %s\n", is.yarnA, is.yarnB, is.uA, is.uB, is.pA.x, is.pA.y, is.pA.z, is.pB.x,
            is.pB.y, is.pB.z, is.distance, is.threshold, toString(is.status), int(is.penetrating), is.reason.c_str());
    for (const ClearanceRegion& g : c.regions)
        put("rg %d %d %a %a %a %a %a %a %s %d %d\n", g.yarnA, g.yarnB, g.uA0, g.uA1, g.uB0, g.uB1, g.worstDistance, g.threshold,
            toString(g.status), int(g.penetrating), g.worstIssue);
    const LocalSignatureResult& s = r.localSignature;
    put("sig %s %d %d %d %d %s\n", toString(s.status), s.pass, s.fail, s.unresolved, s.neighbourhoodSelf, s.detail.c_str());
    for (const LocalPairResult& lp : s.pairs)
        put("lp %d %d %s %d %d %d %d %d %s\n", lp.child, lp.parent, toString(lp.status), lp.crossings, lp.unexpected, lp.overlap,
            lp.terminal, lp.neighbourhoodSelf, lp.detail.c_str());
    return o;
}

std::vector<ExecPolicy> parallelPolicies(int64_t chunk) {
    std::vector<ExecPolicy> v;
    for (int conc : {1, 2, 8, 0}) {
        ExecPolicy par;
        par.parallel = true;
        par.maxConcurrency = conc;
        par.chunk = chunk;
        v.push_back(par);
    }
    return v;
}

TEST(validate_parallel_identical) {
    // validate() with curvature (per span), local signature (per stitch) and clearance (per
    // row of candidate chunk pairs) in parallel at 1, 2, 8 and H chunks at once, small chunks:
    // the complete report equals the serial one, on inputs with clearance FAIL and UNRESOLVED
    // issues, a local-signature FAIL and certified curvature over many spans.
    FixturePiece f = fixturePiece(25.0);
    PatternParams pp;
    PieceKnitParams kp;
    pp.waleSpacing = kp.g.w;
    pp.courseSpacing = kp.g.h;
    PieceKnit pk;
    CHECK(knitPiece({*f.surf, true}, {pp, kp}, knit_test::cpuNumerics(), pk));
    PatternParams thick = pp;
    thick.yarnDiameter = 1.5;
    CylinderPiece cyl = cylinderPiece(24, 5);
    CurveSet bad = cyl.cs;  // one loop moved in front of its parent (invalid_init_rejected)
    const StitchRecord& st = bad.topo.stitches[2 * 24 + 7];
    for (int j = 0; j < kControlsPerStitch; ++j) {
        Vec3& q = bad.yarns[0].curve.ctrl[static_cast<size_t>(st.uBegin) + j + 1];
        q = q + st.N * 2.0;
    }
    bad.yarns[0].finalize();
    struct Case { const char* name; const CurveSet* cs; PatternParams pp; ValidateOptions opt; };
    const std::vector<Case> cases{{"sphere (knitted, not settled)", &pk.curves, pp, {}},
                                  {"sphere at D 1.5", &pk.curves, thick, {false, false}},
                                  {"cylinder, mis-threaded loop", &bad, cyl.pp, {}},
                                  {"cylinder at D 1.3", &cyl.cs, [&] { PatternParams q = cyl.pp; q.yarnDiameter = 1.3; return q; }(), {}}};
    int fails = 0, unres = 0, sigFail = 0;
    for (const Case& c : cases) {
        const ValidationReport serial = validate(*c.cs, c.pp, c.opt);
        const std::string want = reportSnapshot(serial);
        for (const ClearanceIssue& is : serial.clearance.issues) {
            fails += is.status == Status::Fail;
            unres += is.status == Status::Unresolved;
        }
        sigFail += serial.localSignature.fail;
        for (const ExecPolicy& par : parallelPolicies(7)) {
            ValidateOptions o = c.opt;
            o.exec = par;
            const bool same = reportSnapshot(validate(*c.cs, c.pp, o)) == want;
            if (!same) std::printf("  %s differs at concurrency %d\n", c.name, par.maxConcurrency);
            CHECK(same);
        }
        std::printf("  %s: %s, %zu clearance issues, %zu regions, %d signature FAIL: identical at 1, 2, 8 and H\n", c.name,
                    toString(serial.verdict), serial.clearance.issues.size(), serial.clearance.regions.size(),
                    serial.localSignature.fail);
    }
    CHECK(fails > 0 && unres > 0 && sigFail > 0);  // the cases exercise every list
}

TEST(tube_mesh_parallel_identical) {
    // buildTubeMesh with rings in parallel: every array equal to the serial mesh (plain; with the
    // uv wrap, which duplicates rings).
    const SphereFibres& s = sphereFibres();
    const CurveSet& cs = s.pk.curves;
    auto same = [](const TubeMesh& a, const TubeMesh& b) {
        auto eq = [](const auto& x, const auto& y) {
            return x.size() == y.size() && (x.empty() || std::memcmp(x.data(), y.data(), x.size() * sizeof(x[0])) == 0);
        };
        return eq(a.pos, b.pos) && eq(a.nrm, b.nrm) && eq(a.uv, b.uv) && eq(a.s, b.s) && eq(a.t, b.t) && eq(a.theta, b.theta) &&
               eq(a.yarn, b.yarn) && eq(a.axisRatio, b.axisRatio) && eq(a.tri, b.tri) && eq(a.yarnFirstVertex, b.yarnFirstVertex);
    };
    int copies = 0;
    for (int variant = 0; variant < 2; ++variant) {
        TubeParams tp;
        tp.ringSpacing = 0.2;
        tp.sides = 8;
        if (variant == 1) tp.uvWrap = 3;  // u wraps every 7.5 mm of yarn
        const TubeMesh want = buildTubeMesh(cs, s.sm, tp);
        if (variant == 1) copies = static_cast<int>(want.pos.size()) - static_cast<int>(buildTubeMesh(cs, s.sm, TubeParams{0.2, 8}).pos.size());
        for (const ExecPolicy& par : parallelPolicies(5)) CHECK(same(buildTubeMesh(cs, s.sm, tp, par), want));
    }
    CHECK(copies > 0);  // the wrap variant duplicated rings
    std::printf("  %d duplicated wrap vertices; plain and wrap meshes identical at 1, 2, 8 and H\n", copies);
}

TEST(plies_parallel_identical) {
    // buildPlies with samples and plies in parallel: control points, arc-length tables,
    // containment, minimum centre distance and the ply maps equal the serial set.
    const SphereFibres& s = sphereFibres();
    PlyGeometry g;
    g.count = 3;
    const PlySet want = buildPlies(s.pk.curves, s.sm, g, s.pp.radius());
    for (const ExecPolicy& par : parallelPolicies(7)) {
        const PlySet got = buildPlies(s.pk.curves, s.sm, g, s.pp.radius(), par);
        CHECK(got.containment == want.containment && got.minCentreDistance == want.minCentreDistance);
        CHECK(got.yarnOf == want.yarnOf && got.plyOf == want.plyOf && got.parentScale == want.parentScale &&
              got.parentOffset == want.parentOffset);
        CHECK(got.curves.yarns.size() == want.curves.yarns.size());
        for (size_t y = 0; y < want.curves.yarns.size(); ++y) {
            const YarnPath &a = want.curves.yarns[y], &b = got.curves.yarns[y];
            CHECK(a.id == b.id && a.name == b.name && a.curve.ctrl == b.curve.ctrl && a.length() == b.length());
            const YarnPath::Tables ta = a.tables(), tb = b.tables();
            CHECK(ta.tu == tb.tu && ta.ts == tb.ts && ta.preLo == tb.preLo && ta.preHi == tb.preHi);
        }
    }
    std::printf("  %zu plies, containment %.3e mm: identical at 1, 2, 8 and H\n", want.curves.yarns.size(), want.containment);
}

TEST(sampling_parallel_identical) {
    // samplePath with spans and arc lengths in parallel: every array of the polyline equals the
    // serial one (points, parameters, arc lengths, segment bounds, prefix sums, error bounds).
    const SphereFibres& s = sphereFibres();
    for (const YarnPath& y : s.pk.curves.yarns) {
        const Polyline want = samplePath(y, 0.005);
        for (const ExecPolicy& par : parallelPolicies(3)) {
            const Polyline got = samplePath(y, 0.005, 512, par);
            CHECK(got.p == want.p && got.u == want.u && got.s == want.s && got.segLo == want.segLo && got.segHi == want.segHi &&
                  got.cumLo == want.cumLo && got.cumHi == want.cumHi && got.segEps == want.segEps && got.epsBound == want.epsBound);
        }
        std::printf("  %zu segments: identical at 1, 2, 8 and H\n", want.segments());
    }
}

TEST(finalize_parallel_identical) {
    // YarnPath::finalize with spans in parallel: the arc-length table, the length, its error
    // estimate and the bound prefix sums equal the serial ones bit for bit.
    const SphereFibres& s = sphereFibres();
    for (const YarnPath& y0 : s.pk.curves.yarns) {
        YarnPath a = y0, b = y0;
        a.finalize();
        for (const ExecPolicy& par : parallelPolicies(7)) {
            b.finalize(1e-9, par);
            const YarnPath::Tables ta = a.tables(), tb = b.tables();
            CHECK(ta.tu == tb.tu && ta.ts == tb.ts && ta.preLo == tb.preLo && ta.preHi == tb.preHi);
            CHECK(std::memcmp(&ta.length, &tb.length, 8) == 0 && std::memcmp(&ta.lengthErr, &tb.lengthErr, 8) == 0);
            CHECK(ta.stitchesSorted == tb.stitchesSorted && ta.sortedCount == tb.sortedCount);
        }
        std::printf("  %d spans, %zu table nodes: identical at 1, 2, 8 and H\n", a.curve.numSpans(), a.tables().tu.size());
    }
}

TEST(gate_d_parallel_identical) {
    // Gate (d): buildSurface (round and compressed) and checkSurface in parallel at 1, 2, 8 and
    // H chunks at once: every ring field, the compression report and the check result (counters,
    // maxima, first worst point) equal the serial ones, on an unsettled sphere with collisions.
    FixturePiece f = fixturePiece(25.0);
    PatternParams pp;
    PieceKnitParams kp;
    pp.waleSpacing = kp.g.w;
    pp.courseSpacing = kp.g.h;
    PieceKnit pk;
    CHECK(knitPiece({*f.surf, true}, {pp, kp}, knit_test::cpuNumerics(), pk));
    PatternParams thick = pp;
    thick.yarnDiameter = 1.4;  // more contacts: collisions and compression clamps
    auto snap = [&](const ExecPolicy& exec, bool compression) {
        SurfaceParams sp;
        sp.compression = compression;
        sp.ringSpacing = 0.05;
        const SurfaceModel sm = buildSurface(pk.curves, thick, sp, exec);
        const SurfaceCheck sc = checkSurface(pk.curves, sm, {}, exec);
        std::string o;
        char b[400];
        for (const YarnSurface& ys : sm.yarns) {
            std::snprintf(b, sizeof b, "yarn %d %a %a %a %zu\n", ys.yarnId, ys.radius, ys.length, ys.ds, ys.rings.size());
            o += b;
            for (const Ring& R : ys.rings) {
                std::snprintf(b, sizeof b, "%a %a %a %a %a %a %a %a %a %a %a %a %a %a %a %a %a %a\n", R.s, R.u, R.c.x, R.c.y, R.c.z,
                              R.T.x, R.T.y, R.T.z, R.n1.x, R.n1.y, R.n1.z, R.n2.x, R.n2.y, R.n2.z, R.a, R.b, R.phi, 0.0);
                o += b;
            }
        }
        const CompressionReport& c = sm.compression;
        std::snprintf(b, sizeof b, "comp %d %d %a %a\n", c.clamped, c.dropped, c.minAxisRatio, c.maxAreaChange);
        o += b;
        std::snprintf(b, sizeof b, "check %a %d %d %a %s %ld %ld %a %d %d %a %a %a %a %a %s %s\n", sc.minFoldRatio, sc.foldRings,
                      sc.worstFoldYarn, sc.worstFoldS, toString(sc.fold), sc.pointsTested, sc.penetratingPoints,
                      sc.maxDepth, sc.worstYarnA, sc.worstYarnB, sc.worstSA, sc.worstSB, sc.worstPoint.x,
                      sc.worstPoint.y, sc.worstPoint.z, toString(sc.collisions), sc.resolution.c_str());
        o += b;
        return std::make_pair(o, sc.penetratingPoints);
    };
    for (bool compression : {false, true}) {
        const auto [want, pen] = snap(ExecPolicy{}, compression);
        for (const ExecPolicy& par : parallelPolicies(7)) CHECK(snap(par, compression).first == want);
        std::printf("  %s: %ld penetrating points; identical at 1, 2, 8 and H\n", compression ? "compressed" : "round", pen);
        if (!compression) CHECK(pen > 0);  // the check's worst-point merge is exercised
    }
}

TEST(settle_repair_parallel_identical) {
    // Settle and certified repair with parallel inner validations: control points (bitwise),
    // the reports and the final validation equal the serial run's.
    FixturePiece f = fixturePiece(25.0);
    PatternParams pp;
    PieceKnitParams kp;
    pp.waleSpacing = kp.g.w;
    pp.courseSpacing = kp.g.h;
    PieceKnit pk;
    CHECK(knitPiece({*f.surf, true}, {pp, kp}, knit_test::cpuNumerics(), pk));
    auto settleRun = [&](const ExecPolicy& exec, ValidationReport& last) {
        CurveSet cs = pk.curves;
        SettleParams sp;
        sp.exec = exec;
        const SettleReport sr = settle({cs, pp, f.surf.get()}, sp, knit_test::cpuNumerics(), &last);
        std::string o;
        char b[160];
        std::snprintf(b, sizeof b, "settle %d %d %a:", sr.rounds, sr.freeStitches, sr.maxMove);
        o += b;
        for (int d : sr.defectsPerRound) o += " " + std::to_string(d);
        for (const YarnPath& y : cs.yarns)
            for (const Vec3& q : y.curve.ctrl) {
                std::snprintf(b, sizeof b, " %a %a %a", q.x, q.y, q.z);
                o += b;
            }
        return o;
    };
    ValidationReport lastSerial;
    const std::string wantSettle = settleRun(ExecPolicy{}, lastSerial);
    CylinderPiece cyl = cylinderPiece(24, 5);
    PatternParams target = cyl.pp;
    target.yarnDiameter = 1.3;
    auto repairRun = [&](const ExecPolicy& exec, ValidationReport& fin) {
        CurveSet cs = cyl.cs;
        RepairParams rp;
        rp.settle.exec = exec;
        const RepairReport rr = certifiedRepair({cs, target, nullptr}, rp, knit_test::cpuNumerics(), &fin);
        std::string o;
        char b[256];
        std::snprintf(b, sizeof b, "repair %a %a %d %d %d %s\n", rr.D0, rr.targetD, int(rr.reachedTarget), int(rr.moveTried), rr.validations,
                      rr.detail.c_str());
        o += b;
        const RepairMove& m = rr.move;
        std::snprintf(b, sizeof b, "move %d %a %d %d %a %a %d %s\n", int(m.accepted), m.delta, int(m.noPassThrough), int(m.noTubeContact),
                      m.localScale, m.midRhoLo, m.freeStitches, m.why.c_str());
        o += b;
        for (const Vec3& q : cs.yarns[0].curve.ctrl) {
            std::snprintf(b, sizeof b, " %a %a %a", q.x, q.y, q.z);
            o += b;
        }
        return o;
    };
    ValidationReport finSerial;
    const std::string wantRepair = repairRun(ExecPolicy{}, finSerial);
    for (int conc : {2, 0}) {
        ExecPolicy par;
        par.parallel = true;
        par.maxConcurrency = conc;
        par.chunk = 5;
        ValidationReport last, fin;
        CHECK(settleRun(par, last) == wantSettle);
        CHECK(reportSnapshot(last) == reportSnapshot(lastSerial));
        CHECK(repairRun(par, fin) == wantRepair);
        CHECK(reportSnapshot(fin) == reportSnapshot(finSerial));
    }
    std::printf("  settle (3 rounds, gate a %s) and repair (%s): identical at 2 and H\n", toString(lastSerial.verdict),
                toString(finSerial.verdict));
}

// ---- recompute step 1: the knit cache (knit_cache.h) ----

std::string fileBytes(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    return std::string((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
}

// ---- S5: plies as shading maps ----

// Model (1) against a brute-force radial march: the first point (from outside) within r_p
// of any circle centre, refined by bisection.
TEST(ply_profile_circles) {
    // The plan's numbers: r = 1, n = 3, fill 0.95.
    {
        const PlyLayout l = plyLayout(PlyGeometry{}, 1.0);
        CHECK_NEAR(l.rp * l.rp - l.R * l.R * std::sin(kPi / 3) * std::sin(kPi / 3), -0.0210, 5e-4);
        CHECK_NEAR(2 * (60 - std::asin(l.rp / l.R) * 180 / kPi), 9.28, 0.01);
    }
    const double r = 0.56;
    for (int n : {1, 2, 3, 4})
        for (double fill : {0.7, 0.85, 0.95, 1.0}) {
            PlyGeometry g;
            g.count = n;
            g.fill = fill;
            const PlyLayout l = plyLayout(g, r);
            auto inside = [&](double rho, double phi) {
                for (int j = 0; j < n; ++j) {
                    const double a = 2 * kPi * j / n;
                    if (std::hypot(rho * std::cos(phi) - l.R * std::cos(a), rho * std::sin(phi) - l.R * std::sin(a)) <= l.rp) return true;
                }
                return false;
            };
            double worst = 0;
            int hits = 0, misses = 0;
            for (int k = 0; k < 3600; ++k) {
                const double phi = 2 * kPi * (k + 0.37) / 3600;
                const RadialHit h = circleLayoutHit(l, n, phi);
                double hi = 2 * r, lo = -1;
                for (int i = 1; i <= 4000; ++i) {
                    const double rho = 2 * r * (1 - i / 4000.0);
                    if (inside(rho, phi)) { lo = rho; break; }
                    hi = rho;
                }
                if (lo < 0) {
                    // A miss of the march: the exact intersection, if any, is a sliver thinner
                    // than the march step (within 1e-3 of the edge).
                    if (h.hit) CHECK(std::fabs(l.rp * l.rp - l.R * l.R * std::pow(std::sin(phi - 2 * kPi / n * std::round(phi / (2 * kPi / n))), 2)) < 1e-3 * r);
                    ++misses;
                    continue;
                }
                for (int i = 0; i < 80; ++i) (inside(0.5 * (lo + hi), phi) ? lo : hi) = 0.5 * (lo + hi);
                CHECK(h.hit);
                worst = std::max(worst, std::fabs(h.rho - lo));
                ++hits;
            }
            std::printf("  n %d fill %.2f: %d hits, %d misses, max |rho - brute| %.2e mm\n", n, fill, hits, misses, worst);
            CHECK(worst < 1e-9);
            if (n == 1) CHECK(misses == 0);
        }
}

TEST(ply_profile_bridge) {
    const double r = 0.56;
    for (int n : {1, 2, 3, 4})
        for (double fill : {0.7, 0.85, 0.95, 1.0})
            for (PlyProfileModel model : {PlyProfileModel::Circles, PlyProfileModel::Swept}) {
                PlyGeometry g;
                g.count = n;
                g.fill = fill;
                PlyProfileParams pp;
                pp.model = model;
                pp.aoDirections = 32;
                const PlyShadingProfile p = buildPlyShadingProfile(g, r, pp);
                if (n == 1) {
                    for (double phi : {0.0, 1.0, 4.0}) {
                        double s;
                        CHECK(p.eval(phi, &s) == r && s == 0);
                    }
                    continue;
                }
                // Bound on |slope|: the flank up to aL has tilt <= blendTiltDeg; the Hermite piece
                // (Fritsch-Carlson alpha <= 3) has |slope| <= max(|m0|, 1.5 (p0 - floor) / Lh).
                const double tiltBound = std::tan(pp.blendTiltDeg * kPi / 180);
                const double hermBound = std::max(std::fabs(p.m0), 1.5 * (p.p0 - p.floorRho) / p.Lh);
                double maxFd = 0, minRho = 1e300;
                for (int k = 0; k <= 20000; ++k) {
                    const double phi = -kPi + 2 * kPi * (k + 0.5) / 20001;
                    double s;
                    const double rho = p.eval(phi, &s);
                    CHECK(std::isfinite(rho) && std::isfinite(s) && rho > 0);
                    CHECK(std::fabs(s) <= std::max(tiltBound * rho, hermBound) * (1 + 1e-9));
                    minRho = std::min(minRho, rho);
                    // Periodicity 2 pi / n and evenness.
                    double s2, s3;
                    CHECK_NEAR(p.eval(phi + 2 * kPi / n, &s2), rho, 1e-12);
                    CHECK_NEAR(s2, s, 1e-9 * (1 + std::fabs(s)));
                    CHECK_NEAR(p.eval(-phi, &s3), rho, 1e-12);
                    // (At phi = +-pi / n, rounding may put -phi on the other side of the
                    // junction, where the slope is ~1e-9: an absolute floor.)
                    CHECK_NEAR(s3, -s, 1e-7 + 1e-9 * std::fabs(s));
                    // The slope is the derivative of this function.
                    const double d = 1e-6;
                    const double fd = (p.eval(phi + d) - p.eval(phi - d)) / (2 * d);
                    maxFd = std::max(maxFd, std::fabs(fd - s));
                    // On the flank the profile is the model itself (all plies).
                    if (!inBridge(p, phi)) {
                        const RadialHit h = model == PlyProfileModel::Circles ? circleLayoutHit(p.layout, n, phi) : sweptPlyUnion(p.layout, g, phi);
                        CHECK(h.hit);
                        CHECK_NEAR(h.rho, rho, 1e-12);
                    }
                }
                CHECK(maxFd < 1e-4);
                CHECK(minRho >= p.floorRho - 1e-12);
                CHECK(p.floorRho >= p.gapFloor - 1e-9 * r);
                // C1 at every junction: aL, aL + Lh, pi / n (and their mirrors), 0.
                for (double a : {p.aL, p.aL + p.Lh, p.half, 2 * p.half - p.aL - p.Lh, 2 * p.half - p.aL, 0.0}) {
                    const double e = 1e-9;
                    double sm, sp;
                    const double vm = p.eval(a - e, &sm), vp = p.eval(a + e, &sp);
                    CHECK(std::fabs(vp - vm) < 1e-7);
                    CHECK(std::fabs(sp - sm) < 1e-4 * (1 + std::fabs(p.m0)));
                }
                std::printf("  n %d fill %.2f %s: %s, bridge from %.2f deg (half %.1f), Hermite %.2f deg, floor %.4f (d_gap %.4f), rho [%.4f, %.4f] mm, max |fd - slope| %.1e\n",
                            n, fill, model == PlyProfileModel::Circles ? "circles" : "swept  ", p.gap ? "gap   " : "crease", p.aL * 180 / kPi,
                            p.half * 180 / kPi, p.Lh * 180 / kPi, p.floorRho, p.gapFloor, p.rhoMin, p.rhoMax, maxFd);
            }
}

// TubeParams::uvWrap: the same surface (positions, normals, triangles in the same order), u
// differing by integer multiples of uvWrap, no triangle across a wrap, u bounded.
TEST(tube_uv_wrap) {
    const StrandDemo d = makeStrandDemo(60.0, 1.12);
    SurfaceParams sp;
    sp.compression = false;
    sp.ringSpacing = 0.2;
    const SurfaceModel sm = buildSurface(d.curves, d.pp, sp);
    TubeParams tp;
    tp.sides = kExportSides;
    tp.ringSpacing = kExportRing;
    tp.uvRepeat = 2.5;
    const TubeMesh a = buildTubeMesh(d.curves, sm, tp);
    tp.uvWrap = 4;
    const TubeMesh b = buildTubeMesh(d.curves, sm, tp);
    for (size_t i = 0; i < a.pos.size(); ++i) CHECK(a.uv[2 * i] == float(a.s[i] / tp.uvRepeat));
    const size_t stride = tp.sides + 1, extra = b.pos.size() - a.pos.size();
    CHECK(extra % stride == 0);
    const double L = sm.yarns[0].length, du = (L / std::ceil(L / tp.ringSpacing)) / tp.uvRepeat;
    // A copy where a ring's u reaches 4 (not at the end ring): rings j < N with u_j >= 4 k.
    CHECK(extra / stride == static_cast<size_t>(std::floor((L / tp.uvRepeat - du) / tp.uvWrap)));
    CHECK(a.tri.size() == b.tri.size());
    double uMax = 0;
    for (size_t t = 0; t + 2 < b.tri.size(); t += 3) {
        double lo = 1e300, hi = -1e300;
        for (int k = 0; k < 3; ++k) {
            const uint32_t ia = a.tri[t + k], ib = b.tri[t + k];
            CHECK(a.pos[ia] == b.pos[ib] && a.nrm[ia] == b.nrm[ib] && a.uv[2 * ia + 1] == b.uv[2 * ib + 1]);
            const double shift = double(a.uv[2 * ia]) - b.uv[2 * ib];
            CHECK(std::fabs(shift - tp.uvWrap * std::round(shift / tp.uvWrap)) < 1e-4);
            lo = std::min(lo, double(b.uv[2 * ib]));
            hi = std::max(hi, double(b.uv[2 * ib]));
        }
        CHECK(hi - lo <= du + 1e-5);
        uMax = std::max(uMax, hi);
    }
    std::printf("  %zu copied rings; u in [0, %.4f] (unwrapped up to %.4f)\n", extra / stride, uMax, L / tp.uvRepeat);
    CHECK(uMax <= tp.uvWrap + du + 1e-5);
}

// TubeParams::uvPhaseJitter, uvPitchJitter (strand.ply_*_jitter): no jitter is the identity, bit for
// bit; with jitter each yarn's u is (s / uvRepeat) / pitchFactor + phase with both in range, the
// surface unchanged.
TEST(tube_uv_jitter) {
    TubeParams none;
    const YarnUv id = yarnUvOf(none, 7);
    CHECK(id.pitchFactor == 1 && id.phase == 0);
    TubeParams tp;
    tp.uvPhaseJitter = 1;
    tp.uvPitchJitter = 0.1;
    double lo = 1e300, hi = -1e300;
    for (int y = 0; y < 200; ++y) {
        const YarnUv a = yarnUvOf(tp, y), b = yarnUvOf(tp, y);
        CHECK(a.pitchFactor == b.pitchFactor && a.phase == b.phase);
        CHECK(a.pitchFactor >= 0.9 && a.pitchFactor <= 1.1 && a.phase >= 0 && a.phase < 1);
        lo = std::min(lo, a.pitchFactor);
        hi = std::max(hi, a.pitchFactor);
    }
    CHECK(lo < 0.92 && hi > 1.08);  // the range is used
    const StrandDemo d = makeStrandDemo(60.0, 1.12);
    SurfaceParams sp;
    sp.compression = false;
    sp.ringSpacing = 0.2;
    const SurfaceModel sm = buildSurface(d.curves, d.pp, sp);
    tp.sides = kExportSides;
    tp.ringSpacing = kExportRing;
    tp.uvRepeat = 2.5;
    TubeParams plain = tp;
    plain.uvPhaseJitter = plain.uvPitchJitter = 0;
    const TubeMesh a = buildTubeMesh(d.curves, sm, plain), b = buildTubeMesh(d.curves, sm, tp);
    CHECK(a.pos.size() == b.pos.size());
    for (size_t i = 0; i < a.pos.size() && i < b.pos.size(); ++i) {
        CHECK(a.pos[i] == b.pos[i] && a.uv[2 * i + 1] == b.uv[2 * i + 1]);
        const YarnUv y = yarnUvOf(tp, b.yarn[i]);
        CHECK(b.uv[2 * i] == float(b.s[i] / tp.uvRepeat / y.pitchFactor + y.phase));
    }
}

// ------------------------------------------------ explicit knitting problem ----
// Plan rev. 2 Phase 1 (knit_problem.h): the problem's rows, their Jacobians, the contact
// coverage rule, model vs authoritative scoring, and the CPU reference solver.

// The cylinder spiral with a free block of stitches, the cylinder as the band surface.
struct NcpFixture {
    CylinderPiece c;
    std::unique_ptr<ConeSurface> surf;
    ControlMask freeCtrl;
    RelaxMoveParams mp;
    Controls P;
};

NcpFixture ncpFixture(double D, double thrExtra, double cap, int firstStitch, int stitches, double jitter) {
    NcpFixture f;
    f.c = cylinderPiece(24, 5);
    const double R = 24 * f.c.g.w / (2 * kPi), H = (5 + 2) * f.c.g.h;
    f.surf = std::make_unique<ConeSurface>(R, R, H);
    const PatternParams pp = [&] { PatternParams q = f.c.pp; q.yarnDiameter = D; return q; }();
    const double Dc = pp.yarnDiameter + pp.clearance;
    f.mp.thr = Dc + 0.05 + thrExtra;
    f.mp.kappaMax = 1.0 / (1.15 * (pp.radius() + 0.5 * pp.clearance));
    f.mp.window = 0.5 * kPi * Dc * 1.3;
    f.mp.maxDisplacement = cap;
    const YarnPath& y = f.c.cs.yarns[0];
    f.freeCtrl.assign(1, std::vector<char>(y.curve.ctrl.size(), 0));
    for (int t = firstStitch; t < firstStitch + stitches; ++t) {
        const StitchRecord& st = f.c.cs.topo.stitches[t];
        for (int c = static_cast<int>(st.uBegin) + 1; c <= static_cast<int>(st.uEnd) && c + 1 < static_cast<int>(y.curve.ctrl.size()); ++c)
            f.freeCtrl[0][c] = 1;
    }
    f.P = {y.curve.ctrl};
    // Deterministic jitter of the free controls (generic curvature, nonzero cap distances).
    for (size_t c = 0; c < f.P[0].size(); ++c)
        if (f.freeCtrl[0][c]) {
            const double a = 0.7 * c, b = 1.3 * c;
            f.P[0][c] += Vec3{std::sin(a), std::cos(b), std::sin(a + b)} * jitter;
        }
    return f;
}

// Central differences of every row of `fam` touching a variable, away from nonsmooth points
// (both one-sided differences agree): |fd - analytic| <= 1e-6 max(1, |analytic|).
int ncpFdCheck(const KnitProblem& pr, const Controls& P, RowFamily fam, int maxRows) {
    std::vector<Row> rows;
    pr.linearize(P, rows);
    const double h = 1e-6;
    int tested = 0, rowsSeen = 0;
    auto valueOf = [&](const Controls& Q, int v, uint64_t key, double* out) {
        std::vector<Row> rs;
        pr.linearizeVar(Q, v, rs);
        for (const Row& r : rs)
            if (r.key == key) { *out = r.value; return true; }
        return false;
    };
    for (const Row& r : rows) {
        if (r.family != fam || r.n == 0) continue;
        if (++rowsSeen > maxRows) break;
        for (int e = 0; e < r.n; ++e) {
            const auto [y, c] = pr.ctrlOf(r.e[e].var);
            for (int a = 0; a < 3; ++a) {
                Controls Pp = P, Pm = P;
                double* pp = a == 0 ? &Pp[y][c].x : a == 1 ? &Pp[y][c].y : &Pp[y][c].z;
                double* pm = a == 0 ? &Pm[y][c].x : a == 1 ? &Pm[y][c].y : &Pm[y][c].z;
                *pp += h;
                *pm -= h;
                double vp, vm;
                if (!valueOf(Pp, r.e[e].var, r.key, &vp) || !valueOf(Pm, r.e[e].var, r.key, &vm)) continue;
                const double fwd = (vp - r.value) / h, bwd = (r.value - vm) / h;
                if (std::fabs(fwd - bwd) > 1e-4 * (1 + std::fabs(fwd))) continue;  // a kink: not a smooth point
                const double fd = (vp - vm) / (2 * h);
                const double an = a == 0 ? r.e[e].g.x : a == 1 ? r.e[e].g.y : r.e[e].g.z;
                if (std::fabs(fd - an) > 1e-6 * std::max(1.0, std::fabs(an))) {
                    char b[256];
                    std::snprintf(b, sizeof b, "%s row: d/dx%d of control %d: fd %.12g analytic %.12g", toString(fam), a, c, fd, an);
                    throw std::runtime_error(b);
                }
                ++tested;
            }
        }
    }
    return tested;
}

void ncpJacobianCase(RowFamily fam) {
    NcpFixture f = ncpFixture(0.8, 1.0, 0.5, 30, 6, 0.05);
    KnitProblem pr(f.c.cs, f.P, f.P, f.freeCtrl, f.surf.get(), nullptr, f.mp);
    pr.refresh(f.P, nullptr);
    // Edges and caps need a state away from the start / rest configuration.
    Controls Q = f.P;
    for (int v = 0; v < pr.vars(); ++v) {
        const auto [y, c] = pr.ctrlOf(v);
        Q[y][c] += Vec3{0.03 * std::cos(v), 0.02 * std::sin(2.0 * v), 0.025} ;
    }
    const int n = ncpFdCheck(pr, Q, fam, 300);
    std::printf("    %s: %d derivative components checked\n", toString(fam), n);
    CHECK(n >= 30);
}

TEST(ncp_jacobian_fd_edge) { ncpJacobianCase(RowFamily::Edge); }
TEST(ncp_jacobian_fd_curv) { ncpJacobianCase(RowFamily::Curvature); }
TEST(ncp_jacobian_fd_contact) { ncpJacobianCase(RowFamily::Contact); }
TEST(ncp_jacobian_fd_band) {
    ncpJacobianCase(RowFamily::BandUpper);
    ncpJacobianCase(RowFamily::BandLower);
}
TEST(ncp_jacobian_fd_cap) { ncpJacobianCase(RowFamily::Cap); }

TEST(ncp_contact_sign_support) {
    // A contact row's gradient points towards closing the pair: a step against it separates the
    // segments (g decreases); its entries are exactly the free controls the row depends on.
    NcpFixture f = ncpFixture(0.8, 1.0, 0, 30, 6, 0.0);
    KnitProblem pr(f.c.cs, f.P, f.P, f.freeCtrl, nullptr, nullptr, f.mp);
    pr.refresh(f.P, nullptr);
    std::vector<Row> rows;
    pr.linearize(f.P, rows);
    int checked = 0;
    for (const Row& r : rows) {
        if (r.family != RowFamily::Contact || r.n == 0 || checked >= 50) continue;
        ++checked;
        Controls Q = f.P;
        double gNorm2 = 0;
        for (int e = 0; e < r.n; ++e) gNorm2 += dot(r.e[e].g, r.e[e].g);
        if (gNorm2 < 1e-20) { --checked; continue; }  // closest points on fixed controls only: no first-order dependence
        const double step = 1e-6 / std::sqrt(gNorm2);  // a step of length 1e-6: first order dominates
        for (int e = 0; e < r.n; ++e) {
            const auto [y, c] = pr.ctrlOf(r.e[e].var);
            CHECK(pr.varOf(y, c) == r.e[e].var);  // entries are free controls
            Q[y][c] -= r.e[e].g * step;
        }
        std::vector<Row> after;
        pr.linearizeVar(Q, r.e[0].var, after);
        bool found = false;
        for (const Row& a : after)
            if (a.key == r.key) { found = true; CHECK(a.value < r.value); }
        CHECK(found);
        // Moving a free control that is not an entry leaves the row's value unchanged.
        const auto [y0, c0] = pr.ctrlOf(r.e[0].var);
        for (int dc = -8; dc <= 8; ++dc) {
            const int c = c0 + dc;
            if (c < 0 || c >= static_cast<int>(f.P[y0].size())) continue;
            const int v = pr.varOf(y0, c);
            bool entry = v < 0;
            for (int e = 0; e < r.n; ++e) entry = entry || r.e[e].var == v;
            if (entry) continue;
            Controls M = f.P;
            M[y0][c] += Vec3{0.3, -0.2, 0.1};
            std::vector<Row> rm;
            pr.linearizeVar(M, r.e[0].var, rm);
            for (const Row& a : rm)
                if (a.key == r.key) CHECK(a.value == r.value);
        }
    }
    CHECK(checked > 0);
}

TEST(ncp_fixed_eliminated) {
    NcpFixture f = ncpFixture(0.8, 1.0, 0.5, 30, 6, 0.0);
    KnitProblem pr(f.c.cs, f.P, f.P, f.freeCtrl, f.surf.get(), nullptr, f.mp);
    pr.refresh(f.P, nullptr);
    int nFree = 0;
    for (char m : f.freeCtrl[0]) nFree += m != 0;
    CHECK(pr.vars() == nFree);
    std::vector<Row> rows;
    pr.linearize(f.P, rows);
    for (const Row& r : rows)
        for (int e = 0; e < r.n; ++e) {
            CHECK(r.e[e].var >= 0 && r.e[e].var < pr.vars());
            const auto [y, c] = pr.ctrlOf(r.e[e].var);
            CHECK(f.freeCtrl[y][c]);
        }
    // A fixed control next to the free block enters the rows' values, not their variables.
    int cFixed = -1;
    for (int c = 1; c + 1 < static_cast<int>(f.freeCtrl[0].size()); ++c)
        if (!f.freeCtrl[0][c] && f.freeCtrl[0][c + 1]) { cFixed = c; break; }
    CHECK(cFixed >= 0);
    Controls Q = f.P;
    Q[0][cFixed] += Vec3{0.2, 0.1, -0.1};
    std::vector<Row> before, after;
    pr.linearizeVar(f.P, pr.varOf(0, cFixed + 1), before);
    pr.linearizeVar(Q, pr.varOf(0, cFixed + 1), after);
    bool changed = false;
    for (size_t i = 0; i < before.size() && i < after.size(); ++i) changed = changed || before[i].value != after[i].value;
    CHECK(changed);
}

TEST(ncp_curv_degenerate_flagged) {
    // A span whose controls coincide has no tangent: counted degenerate, never a satisfied row.
    std::vector<Vec3> ctrl = lineCtrl({0, 0, 0}, {10, 0, 0}, 10);
    for (int c = 5; c <= 8; ++c) ctrl[c] = ctrl[5];
    CurveSet cs = makeSet({makeYarn(ctrl, 0)});
    ControlMask freeCtrl(1, std::vector<char>(ctrl.size(), 0));
    freeCtrl[0][6] = 1;
    RelaxMoveParams mp;
    mp.thr = 1.0;
    mp.kappaMax = 1.0;
    mp.window = 100;
    Controls P{ctrl};
    KnitProblem pr(cs, P, P, freeCtrl, nullptr, nullptr, mp);
    pr.refresh(P, nullptr);
    const KnitProblemScore sc = pr.scoreGeometry(P, nullptr);
    CHECK(sc.degenerate > 0);
    std::vector<Row> rows;
    pr.linearize(P, rows);
    for (const Row& r : rows)
        if (r.family == RowFamily::Curvature) CHECK(std::isfinite(r.value));
}

TEST(ncp_band_owned_only) {
    // A bound control whose continuation is not Owned (here: past the mesh's open boundary) gets
    // no band row and keeps its state; the others get both band rows.
    TriMesh m;
    const int N = 10;
    for (int j = 0; j <= N; ++j)
        for (int i = 0; i <= N; ++i) m.pos.push_back({2.0 * i, 2.0 * j, 0});
    for (int j = 0; j < N; ++j)
        for (int i = 0; i < N; ++i) {
            const int a = j * (N + 1) + i;
            m.tri.push_back({a, a + 1, a + N + 2});
            m.tri.push_back({a, a + N + 2, a + N + 1});
        }
    m.computeNormals();
    MeshSurface ms(m);
    std::vector<Vec3> ctrl = lineCtrl({1, 10, 0.5}, {19, 10, 0.5}, 9);
    CurveSet cs = makeSet({makeYarn(ctrl, 0)});
    ControlMask freeCtrl(1, std::vector<char>(ctrl.size(), 0));
    for (size_t c = 1; c + 1 < ctrl.size(); ++c) freeCtrl[0][c] = 1;
    SurfaceBinding bind;
    bind.surface = &ms;
    bind.state.assign(1, std::vector<ProjState>(ctrl.size()));
    bind.result.assign(1, std::vector<ProjResult>(ctrl.size(), ProjResult::Owned));
    for (size_t c = 0; c < ctrl.size(); ++c) {
        const Vec3 q = ctrl[c];
        const int i = std::clamp(static_cast<int>(q.x / 2), 0, N - 1), j = std::clamp(static_cast<int>(q.y / 2), 0, N - 1);
        const int t = 2 * (j * N + i) + 1;  // a triangle of the cell; the walk finds the foot
        const double b[3] = {1.0 / 3, 1.0 / 3, 1.0 / 3};
        bind.state[0][c] = ms.stateAt(t, b, q);
    }
    RelaxMoveParams mp;
    mp.thr = 0.5;
    mp.kappaMax = 10;
    mp.window = 100;
    mp.bandHalf = 0.2;
    Controls P{ctrl};
    const int cOut = 9;  // moved past the mesh's edge at x = 20
    P[0][cOut] = {26, 10, 0.5};
    const ProjState stateBefore = bind.state[0][cOut];
    KnitProblem pr(cs, P, P, freeCtrl, &ms, &bind, mp);
    pr.refresh(P, &bind);
    CHECK(bind.result[0][cOut] != ProjResult::Owned);
    CHECK(bind.state[0][cOut] == stateBefore);
    std::vector<Row> rows;
    pr.linearize(P, rows);
    std::set<int> banded;
    for (const Row& r : rows)
        if (r.family == RowFamily::BandUpper || r.family == RowFamily::BandLower) banded.insert(r.e[0].var);
    CHECK(banded.count(pr.varOf(0, cOut)) == 0);
    CHECK(static_cast<int>(banded.size()) == pr.vars() - 1);
    const KnitProblemScore sc = pr.scoreGeometry(P, &bind);
    CHECK(sc.notOwned == 1);
    CHECK(sc.maxViol[static_cast<int>(RowFamily::BandUpper)] > 0);  // z = 0.5 > bandHalf 0.2 on the others
}

TEST(ncp_refresh_cumulative_motion) {
    // thr 1, delta 0.1: two parallel yarns 1.21 apart are not candidates (pad thr + 2 delta =
    // 1.2). Two moves of 0.09 per side bring them to 0.85 < thr: after the first the coverage
    // holds, after the second it does not and the score says so (pairs may be missing) instead
    // of reporting a feasible state; a refresh restores the pairs. The reference solver, which
    // checks cumulative motion, ends with no pair below thr (brute force over all segments).
    auto line = [](double y) { return lineCtrl({0, y, 0}, {12, y, 0}, 12); };
    CurveSet cs = makeSet({makeYarn(line(0), 0), makeYarn(line(1.21), 1)});
    ControlMask freeCtrl(2, std::vector<char>(cs.yarns[0].curve.ctrl.size(), 0));
    for (int y = 0; y < 2; ++y)
        for (size_t c = 1; c + 1 < freeCtrl[y].size(); ++c) freeCtrl[y][c] = 1;
    RelaxMoveParams mp;
    mp.thr = 1.0;
    mp.kappaMax = 10;
    mp.window = 1e9;  // same-yarn pairs excluded
    KnitProblemParams kp;
    kp.delta = 0.1;
    Controls P{cs.yarns[0].curve.ctrl, cs.yarns[1].curve.ctrl};
    KnitProblem pr(cs, P, P, freeCtrl, nullptr, nullptr, mp, kp);
    pr.refresh(P, nullptr);
    CHECK(pr.candidates() == 0);
    auto moveTogether = [&](Controls& Q, double d) {
        for (size_t c = 0; c < Q[0].size(); ++c) {
            if (freeCtrl[0][c]) Q[0][c].y += d;
            if (freeCtrl[1][c]) Q[1][c].y -= d;
        }
    };
    moveTogether(P, 0.09);
    KnitProblemScore s1 = pr.scoreGeometry(P, nullptr);
    CHECK(s1.covered && s1.maxConstraint() == 0);
    moveTogether(P, 0.09);
    KnitProblemScore s2 = pr.scoreGeometry(P, nullptr);
    CHECK(!s2.covered);                      // flagged, not silently feasible
    pr.refresh(P, nullptr);
    CHECK(pr.candidates() > 0);
    KnitProblemScore s3 = pr.scoreGeometry(P, nullptr);
    CHECK(s3.covered && s3.maxViol[static_cast<int>(RowFamily::Contact)] > 0);

    // Inside delta of the refresh the covered pairs hold the closest one: the score's contact
    // violation equals a brute force over every pair of sample segments of the two yarns.
    moveTogether(P, 0.04);
    const KnitProblemScore s4 = pr.scoreGeometry(P, nullptr);
    CHECK(s4.covered);
    CurveSet out = cs;
    for (int y = 0; y < 2; ++y) { out.yarns[y].curve.ctrl = P[y]; out.yarns[y].finalize(); }
    double minD = 1e300;
    const int sps = mp.samplesPerSpan;
    auto pts = [&](const YarnPath& yp) {
        std::vector<Vec3> v;
        for (int k = 0; k < yp.curve.numSpans(); ++k)
            for (int j = 0; j < sps; ++j) v.push_back(yp.curve.eval(k + double(j) / sps));
        return v;
    };
    const std::vector<Vec3> a = pts(out.yarns[0]), b = pts(out.yarns[1]);
    for (size_t i = 0; i + 1 < a.size(); ++i)
        for (size_t j = 0; j + 1 < b.size(); ++j)
            minD = std::min(minD, closestSegmentPoints(a[i], a[i + 1], b[j], b[j + 1]).distance);
    std::printf("    brute force %.6f, score %.6f\n", minD, mp.thr - s4.maxViol[static_cast<int>(RowFamily::Contact)] * pr.D());
    CHECK_NEAR(s4.maxViol[static_cast<int>(RowFamily::Contact)], (mp.thr - minD) / pr.D(), 1e-9);
}

TEST(ncp_model_vs_score_divergence) {
    // The linear model of a contact (s, t, n frozen) can call a moved pair separated while the
    // actual segments come closer than thr: the authoritative score sees the violation, and the
    // rows re-linearized at the moved controls (what acceptance uses) carry it.
    auto line = [](Vec3 a, Vec3 b) { return lineCtrl(a, b, 3); };
    CurveSet cs = makeSet({makeYarn(line({-6, 0, 0}, {6, 0, 0}), 0), makeYarn(line({0, -6, 1.2}, {0, 6, 1.2}), 1)});
    ControlMask freeCtrl{std::vector<char>(cs.yarns[0].curve.ctrl.size(), 0), std::vector<char>(cs.yarns[1].curve.ctrl.size(), 0)};
    for (size_t c = 1; c + 1 < freeCtrl[1].size(); ++c) freeCtrl[1][c] = 1;
    RelaxMoveParams mp;
    mp.thr = 1.0;
    mp.kappaMax = 100;
    mp.window = 1e9;
    KnitProblemParams kp;
    kp.delta = 10;  // coverage is not the subject here
    Controls P{cs.yarns[0].curve.ctrl, cs.yarns[1].curve.ctrl};
    KnitProblem pr(cs, P, P, freeCtrl, nullptr, nullptr, mp, kp);
    pr.refresh(P, nullptr);
    std::vector<Row> rows;
    pr.linearize(P, rows);
    bool found = false;
    uint64_t seed = 12345;
    auto rnd = [&] { seed = seed * 6364136223846793005ull + 1442695040888963407ull; return ((seed >> 11) * (1.0 / 9007199254740992.0)) * 2 - 1; };
    for (int trial = 0; trial < 2000 && !found; ++trial) {
        Controls Q = P;
        for (int v = 0; v < pr.vars(); ++v) {
            const auto [y, c] = pr.ctrlOf(v);
            Q[y][c] += Vec3{rnd(), rnd(), rnd()} * 1.5;
        }
        std::vector<Row> moved;
        pr.linearize(Q, moved);
        for (const Row& r : rows) {
            if (r.family != RowFamily::Contact) continue;
            double pred = r.value;
            for (int e = 0; e < r.n; ++e) {
                const auto [y, c] = pr.ctrlOf(r.e[e].var);
                pred += dot(r.e[e].g, Q[y][c] - P[y][c]);
            }
            for (const Row& m : moved)
                if (m.key == r.key && pred <= 0 && m.value > 0.05) {
                    found = true;
                    CHECK(pr.scoreGeometry(Q, nullptr).maxViol[static_cast<int>(RowFamily::Contact)] >= m.value);
                }
        }
    }
    CHECK(found);
}

TEST(ncp_same_yarn_arc_gap) {
    // Stop point 1, class (ii): in long spans (here 4 mm, ~1 mm per sample segment, as in the
    // CLOSE tail) adjacent segments of one yarn share an endpoint (distance 0). The legacy rule
    // (first samples' arc distance < window) keeps them as candidates; the problem's rule (arc gap
    // between the segments < window) excludes them, so no row is unsatisfiable by construction.
    std::vector<Vec3> ctrl = lineCtrl({0, 0, 0}, {40, 0, 0}, 10);
    CurveSet cs = makeSet({makeYarn(ctrl, 0)});
    ControlMask freeCtrl(1, std::vector<char>(ctrl.size(), 0));
    for (size_t c = 1; c + 1 < ctrl.size(); ++c) freeCtrl[0][c] = 1;
    RelaxMoveParams mp;
    mp.thr = 0.4;     // the measured CLOSE case: segment ~1 mm > window 0.69 > thr 0.39 (repair at D 0.24)
    mp.kappaMax = 1.0;
    mp.window = 0.7;
    Controls P{ctrl};
    KnitProblem pr(cs, P, P, freeCtrl, nullptr, nullptr, mp);
    pr.refresh(P, nullptr);
    CHECK(pr.candidates() == 0);  // a straight yarn: no nonlocal pair, and adjacent ones are excluded
    const KnitProblemScore sc = pr.scoreGeometry(P, nullptr);
    CHECK(sc.maxViol[static_cast<int>(RowFamily::Contact)] == 0);
}

TEST(ncp_hard_target_split) {
    // With hard thresholds the C / K rows are hard at hardContact / hardKappa and the move's
    // thr / kappaMax become objective rows that exist only while short of the target, with the
    // same gradients (finite differences on the target families).
    NcpFixture f = ncpFixture(0.8, 1.0, 0.5, 30, 6, 0.05);
    KnitProblemParams kp;
    kp.hardContact = 0.9;                // below the move's thr (~1.95): targets active
    kp.hardKappa = 2 * f.mp.kappaMax;    // looser than kappaMax
    f.mp.kappaMax *= 0.05;               // a target the spiral's spans exceed: curvature-target rows active
    KnitProblem pr(f.c.cs, f.P, f.P, f.freeCtrl, f.surf.get(), nullptr, f.mp, kp);
    pr.refresh(f.P, nullptr);
    std::vector<Row> rows;
    pr.linearize(f.P, rows);
    int hardC = 0, tgtC = 0;
    for (const Row& r : rows) {
        if (r.family == RowFamily::ContactTarget) {
            ++tgtC;
            CHECK(r.value > 0);  // only while short of thr
        }
        if (r.family == RowFamily::Contact) ++hardC;
    }
    CHECK(hardC > 0 && tgtC > 0 && tgtC <= hardC);
    for (const Row& r : rows)
        if (r.family == RowFamily::Contact) {
            for (const Row& t : rows)
                if (t.family == RowFamily::ContactTarget && (t.key ^ (uint64_t(RowFamily::ContactTarget) << 60)) ==
                                                                (r.key ^ (uint64_t(RowFamily::Contact) << 60))) {
                    CHECK_NEAR(t.value - r.value, (f.mp.thr - kp.hardContact) / f.mp.thr, 1e-12);  // both normalized by thr
                }
        }
    Controls Q = f.P;
    for (int v = 0; v < pr.vars(); ++v) {
        const auto [y, c] = pr.ctrlOf(v);
        Q[y][c] += Vec3{0.03 * std::cos(v), 0.02 * std::sin(2.0 * v), 0.025};
    }
    CHECK(ncpFdCheck(pr, Q, RowFamily::ContactTarget, 300) >= 30);
    const int nK = ncpFdCheck(pr, Q, RowFamily::CurvatureTarget, 300);
    std::printf("    contact target and %d curvature-target derivative components checked\n", nK);
    CHECK(nK >= 30);
    const KnitProblemScore sc = pr.scoreGeometry(Q, nullptr);
    CHECK(sc.maxViol[static_cast<int>(RowFamily::ContactTarget)] > 0);
    CHECK(sc.maxConstraint() <= std::max({sc.maxViol[1], sc.maxViol[2], sc.maxViol[3], sc.maxViol[4], sc.maxViol[5]}));
}

TEST(ncp_pass_through_detected) {
    // Yarn B crosses above yarn A (0.6 mm). Moving B's controls through A to 0.6 mm below is a
    // pass-through (flagged); moving them up or sideways is not.
    CurveSet cs = makeSet({makeYarn(lineCtrl({-6, 0, 0}, {6, 0, 0}, 6), 0), makeYarn(lineCtrl({0, -6, 0.6}, {0, 6, 0.6}, 6), 1)});
    ControlMask freeCtrl{std::vector<char>(cs.yarns[0].curve.ctrl.size(), 0), std::vector<char>(cs.yarns[1].curve.ctrl.size(), 0)};
    for (size_t c = 1; c + 1 < freeCtrl[1].size(); ++c) freeCtrl[1][c] = 1;
    RelaxMoveParams mp;
    mp.thr = 1.0;
    mp.kappaMax = 100;
    mp.window = 1e9;
    Controls P{cs.yarns[0].curve.ctrl, cs.yarns[1].curve.ctrl};
    KnitProblem pr(cs, P, P, freeCtrl, nullptr, nullptr, mp);
    pr.refresh(P, nullptr);
    const int c = static_cast<int>(P[1].size()) / 2;  // the control above the crossing
    const int v = pr.varOf(1, c);
    CHECK(v >= 0);
    const Vec3 x0 = P[1][c];
    const Vec3 through = x0 + Vec3{0, 0, -2.4}, away = x0 + Vec3{0, 0, 0.5}, along = x0 + Vec3{0, 0.4, 0};
    CHECK(pr.passThrough(P, v, x0, through, 2.0));   // through A
    CHECK(P[1][c] == through);                       // P holds `to`
    P[1][c] = x0;
    CHECK(!pr.passThrough(P, v, x0, away, 2.0));
    CHECK(!pr.passThrough(P, v, x0, along, 2.0));
    P[1][c] = x0;
}

} // namespace

int main(int argc, char** argv) {
    if (argc > 1 && std::strcmp(argv[1], "--list") == 0) {
        for (const Case& c : registry()) std::printf("%s\n", c.name);
        return 0;
    }
    int failed = 0, ran = 0, skipped = 0;
    for (const Case& c : registry()) {
        if (argc > 1 && std::strcmp(argv[1], c.name) != 0) continue;
        ++ran;
        try {
            c.fn();
            std::printf("[ OK ] %s\n", c.name);
        } catch (const Skipped& s) {
            ++skipped;
            std::printf("[SKIP] %s: %s\n", c.name, s.why);
        } catch (const std::exception& e) {
            ++failed;
            std::printf("[FAIL] %s: %s\n", c.name, e.what());
        }
    }
    if (ran == 0) { std::printf("no test named %s\n", argc > 1 ? argv[1] : "?"); return 2; }
    if (skipped) std::printf("%d skipped\n", skipped);
    if (failed) return 1;
    return skipped == ran ? 77 : 0;  // 77: every case run was skipped (CTest SKIP_RETURN_CODE)
}

TEST(fibers_root_bridge_serialized) {
    // On the float32 export of the capped sphere: a grazing fibre is caught and a clear one passes;
    // the searched slab range follows the root's reach (more than +-3 slabs at a 0.1 mm ring
    // spacing).
    const SphereFibres& s = sphereFibres();
    const double margin = 1e-4 / s.map.sigma;
    AttachmentResult ar = rootAttachmentCheck(*s.tubes, s.fs, margin);
    std::printf("  %d fibres: buried %d, single exit %d, clear %d; %d slabs at most\n", ar.fibres, ar.buried, ar.singleExit, ar.clear,
                ar.maxSlabs);
    CHECK(ar.fibres == static_cast<int>(s.fs.fibres.size()));
    const ExportTubeIndex fineTubes(s.pk.curves, s.pp, ExportTubeSpec{kExportModelRing, 0.1, kExportSides, false,
                                    [&](const Vec3& p) { return s.map.toExportF(p); }, s.eSer});
    AttachmentResult fine = rootAttachmentCheck(fineTubes, s.fs, margin);
    CHECK(fine.maxSlabs > 7);
    // Grazing: a polyline passing within r_i of a distant turn of the yarn.
    const YarnPath& y = s.pk.curves.yarns[0];
    const double sC = 0.5 * y.length();
    const double u = y.paramAtArcLength(sC);
    const Vec3 c = y.curve.eval(u), T = normalized(y.curve.d1(u));
    const Vec3 n = normalized(s.fr.fabric(c).n);
    const Vec3 e = normalized(n - T * dot(n, T));
    const double r = s.pp.radius(), ri = 0.035;
    Fibre g;
    g.yarnId = 0;
    g.s = sC - 40.0;
    g.radius0 = float(ri);
    g.radius1 = float(0.5 * ri);
    g.attach = 0.1f;
    g.attachArc = 0.5;
    g.emergeL = 0;
    // r cos(pi / 8) is the exported octagon's inscribed radius: at r_i / 2 beyond it, the
    // point is within r_i / 2 of a facet whatever its angle (the smooth tube is no reference).
    const double rIn = r * std::cos(kPi / 8);
    FibreSet grazing, clear;
    grazing.addFibre(g, std::vector<Vec3>{c + e * (rIn + 0.5 * ri) - T * 0.3, c + e * (rIn + 0.5 * ri) + T * 0.3});
    CHECK(!s.tubes->check(grazing, grazing.fibres[0], 1e-4).clear);
    clear.addFibre(g, std::vector<Vec3>{c + e * (r + 3 * ri + 0.05) - T * 0.3, c + e * (r + 3 * ri + 0.05) + T * 0.3});
    CHECK(s.tubes->check(clear, clear.fibres[0], 1e-4).clear);
}

TEST(export_tube_accepts_matches_check) {
    // Every fibre grown on the 3-ply sphere (the light rules leave some that the serialized check
    // refuses): accepts() == check().ok().
    const SphereFibres& s = sphereFibres();
    PlyGeometry g;
    const PlySet ps = buildPlies(s.pk.curves, s.sm, g, s.pp.radius());
    PatternParams pp = s.pp;
    pp.yarnDiameter = 2 * ps.layout.rp;
    SurfaceParams sp;
    sp.compression = false;
    const SurfaceModel sm = buildSurface(ps.curves, pp, sp);
    const SceneMapping map = s.map;
    const ExportTubeIndex tubes(ps.curves, pp, ExportTubeSpec{kExportModelRing, kExportPlyRing, kExportPlySides, false,
                                [map](const Vec3& p) { return map.toExportF(p); }, s.eSer});
    FibreSetup set;
    set.tubeRadius = pp.radius();
    set.sides = kExportPlySides;
    set.ring = kExportPlyRing;
    set.plies = &ps;
    set.parentCurves = &s.pk.curves;
    set.parentSurface = &s.sm;
    FibreParams fp = s.fp;
    fp.randomFraction = 0.5;
    fp.density = 10;
    FibreFrame fr = s.fr;
    setupPieceFibres(fp, fr, s.pp, set);
    int checked = 0, rejected = 0, mismatches = 0;
    const FibreSet fs = growYarnFibres(ps.curves, sm, pp, fp, &fr);
    for (const Fibre& f : fs.fibres) {
        const ExportTubeIndex::TubeVerdict a = tubes.check(fs, f, 1e-4);
        ++checked;
        rejected += !a.ok();
        mismatches += tubes.accepts(fs, f, 1e-4) != a.ok();
    }
    std::printf("  %d fibres checked (%d rejected): %d mismatches\n", checked, rejected, mismatches);
    CHECK(checked > 1000 && rejected > 10);
    CHECK(mismatches == 0);
}
