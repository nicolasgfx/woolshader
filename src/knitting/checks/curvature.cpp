// Speed and curvature bounds (curve-set check).
//
// Per piece [t0,t1] of a span:
//   C'  restricted hodograph (quadratic Bernstein, values are d/du)
//   C'' linear: its Bernstein coefficients are its end values
//   C' x C''  degree-3 Bernstein with coefficients b_k (product formula)
//   speed_lo = dist(0, hull{h0,h1,h2}) - eta1
//   kappa_hi = (max|b_k| + err) / speed_lo^3
// No root solving: everything is a convex-hull bound plus an explicit allowance.
#include <algorithm>
#include <cmath>

#include "knitting/checks/checks.h"

namespace knit {

namespace {

double originSegmentDistance(const Vec3& a, const Vec3& b) {
    const Vec3 ab = b - a;
    const double l2 = dot(ab, ab);
    const double t = l2 > 0 ? std::clamp(-dot(a, ab) / l2, 0.0, 1.0) : 0.0;
    return norm(a + ab * t);
}

// Distance from the origin to triangle (a,b,c) (Ericson, Real-Time Collision Detection).
// A degenerate (collinear) triangle is its longest edge's segment: Ericson's
// barycentric branch would divide by its zero area (e.g. a span whose first control
// points are collinear, as at a natural end, has two equal hodograph controls).
double originTriangleDistance(const Vec3& a, const Vec3& b, const Vec3& c) {
    const Vec3 p{0, 0, 0};
    Vec3 ab = b - a, ac = c - a, ap = p - a;
    const Vec3 n = cross(ab, ac);
    if (dot(n, n) <= 1e-24 * dot(ab, ab) * dot(ac, ac) + 1e-300)
        return std::min({originSegmentDistance(a, b), originSegmentDistance(a, c), originSegmentDistance(b, c)});
    double d1 = dot(ab, ap), d2 = dot(ac, ap);
    if (d1 <= 0 && d2 <= 0) return norm(a);
    Vec3 bp = p - b;
    double d3 = dot(ab, bp), d4 = dot(ac, bp);
    if (d3 >= 0 && d4 <= d3) return norm(b);
    double vc = d1 * d4 - d3 * d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) return norm(a + ab * (d1 / (d1 - d3)));
    Vec3 cp = p - c;
    double d5 = dot(ab, cp), d6 = dot(ac, cp);
    if (d6 >= 0 && d5 <= d6) return norm(c);
    double vb = d5 * d2 - d1 * d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) return norm(a + ac * (d2 / (d2 - d6)));
    double va = d3 * d6 - d5 * d4;
    if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0)
        return norm(b + (c - b) * ((d4 - d3) / ((d4 - d3) + (d5 - d6))));
    double denom = 1.0 / (va + vb + vc);
    double v = vb * denom, w = vc * denom;
    return norm(a + ab * v + ac * w);
}

struct PieceBound {
    double speedLo, kappaHi;
};

PieceBound boundPiece(const CubicBezier& span, double t0, double t1, double floorAbs) {
    QuadBezier h = span.hodograph().sub(t0, t1);
    Vec3 s0 = lerp(span.d2At0(), span.d2At1(), t0);
    Vec3 s1 = lerp(span.d2At0(), span.d2At1(), t1);
    double hMax = std::max({norm(h.h[0]), norm(h.h[1]), norm(h.h[2])});
    double sMax = std::max(norm(s0), norm(s1));
    double eta1 = 1e-12 * hMax + floorAbs, eta2 = 1e-12 * sMax + floorAbs;
    Vec3 b[4] = {cross(h.h[0], s0), (cross(h.h[0], s1) + cross(h.h[1], s0) * 2.0) / 3.0,
                 (cross(h.h[1], s1) * 2.0 + cross(h.h[2], s0)) / 3.0, cross(h.h[2], s1)};
    double numHi = 0;
    for (const Vec3& v : b) numHi = std::max(numHi, norm(v));
    numHi += eta1 * sMax + eta2 * hMax + eta1 * eta2;
    double speedLo = originTriangleDistance(h.h[0], h.h[1], h.h[2]) - eta1;
    double kappaHi = speedLo > 0 ? numHi / (speedLo * speedLo * speedLo) : 1e300;
    return {speedLo, kappaHi};
}

// One span's accumulator: the running values start as for the whole yarn.
struct SpanBound {
    CurvatureResult r;
    double spanKappa = 0;  // the span's certified bound
    double sampledSpeedMin = 1e300;
};

// A piece [t0, t1] of a span on the explicit stack.
struct SpanPiece {
    double t0, t1;
    int depth;
};

// Sampled curvature on a piece of span k (estimate + lower-bound witness); returns its maximum.
double sampleCurvature(const CubicBezier& sp, int k, const SpanPiece& piece, double floorAbs, SpanBound& acc) {
    const double t0 = piece.t0, t1 = piece.t1;
    CurvatureResult& r = acc.r;
    double kSamp = 0;
    for (int j = 0; j <= 4; ++j) {
        double t = t0 + (t1 - t0) * j / 4.0;
        Vec3 d1 = sp.d1(t), d2 = sp.d2(t);
        double spd = norm(d1);
        acc.sampledSpeedMin = std::min(acc.sampledSpeedMin, spd);
        if (spd <= 0) continue;
        double num = norm(cross(d1, d2));
        double kap = num / (spd * spd * spd);
        kSamp = std::max(kSamp, kap);
        double eta1 = 1e-12 * spd + floorAbs, eta2 = 1e-12 * norm(d2) + floorAbs;
        double errNum = eta1 * norm(d2) + eta2 * spd + eta1 * eta2;
        double kLo = std::max(0.0, num - errNum) / std::pow(spd + eta1, 3);
        if (kLo > r.kappaWitnessLo) { r.kappaWitnessLo = kLo; r.uWitness = k + t; }
        if (kap > 0 && 1.0 / kap < r.rhoEst) { r.rhoEst = 1.0 / kap; r.uRhoEst = k + t; }
    }
    return kSamp;
}

// Span k: pieces halved until the bound is within 1 % of the sampled curvature, or far below 1/r.
SpanBound boundSpan(const YarnPath& y, int k, double floorAbs, double kappaFloor) {
    const int maxDepth = 18;
    SpanBound acc;
    acc.r.speedLo = 1e300;
    acc.r.rhoEst = 1e300;
    const CubicBezier sp = y.curve.span(k);
    std::vector<SpanPiece> stack{{0.0, 1.0, 0}};
    while (!stack.empty()) {
        const SpanPiece it = stack.back();
        stack.pop_back();
        PieceBound pb = boundPiece(sp, it.t0, it.t1, floorAbs);
        const double kSamp = sampleCurvature(sp, k, it, floorAbs, acc);
        bool tight = pb.speedLo > 0 && pb.kappaHi <= std::max(1.01 * kSamp, kappaFloor);
        if (tight || it.depth >= maxDepth) {
            acc.r.speedLo = std::min(acc.r.speedLo, pb.speedLo);
            acc.r.kappaHi = std::max(acc.r.kappaHi, pb.kappaHi);
            acc.spanKappa = std::max(acc.spanKappa, pb.speedLo > 0 ? pb.kappaHi : 1e300);
            ++acc.r.leaves;
        } else {
            double m = 0.5 * (it.t0 + it.t1);
            stack.push_back({it.t0, m, it.depth + 1});
            stack.push_back({m, it.t1, it.depth + 1});
        }
    }
    return acc;
}

} // namespace

CurvatureResult checkCurvature(const YarnPath& y, double radius, const Allowances& al, const ExecPolicy& exec) {
    CurvatureResult r;
    r.yarnId = y.id;
    r.speedLo = 1e300;
    r.rhoEst = 1e300;
    const double floorAbs = 1e-15 * al.Lref;
    const double kappaFloor = 0.05 / radius;  // no need to tighten bounds far below 1/r
    double sampledSpeedMin = 1e300;
    const int spans = y.curve.numSpans();
    r.spanKappaHi.assign(spans, 0.0);
    // Merged in span order with the same strict comparisons, so extrema and their first locations
    // are those of one pass over the spans.
    std::vector<SpanBound> acc(static_cast<size_t>(spans));
    parallelFor(spans, exec.chunk > 0 ? exec.chunk : 64, [&](int64_t kBegin, int64_t kEnd, int64_t) {
        for (int k = static_cast<int>(kBegin); k < static_cast<int>(kEnd); ++k) acc[k] = boundSpan(y, k, floorAbs, kappaFloor);
    }, exec);
    for (int k = 0; k < spans; ++k) {
        const SpanBound& a = acc[k];
        r.spanKappaHi[k] = a.spanKappa;
        r.speedLo = std::min(r.speedLo, a.r.speedLo);
        r.kappaHi = std::max(r.kappaHi, a.r.kappaHi);
        r.leaves += a.r.leaves;
        sampledSpeedMin = std::min(sampledSpeedMin, a.sampledSpeedMin);
        if (a.r.kappaWitnessLo > r.kappaWitnessLo) { r.kappaWitnessLo = a.r.kappaWitnessLo; r.uWitness = a.r.uWitness; }
        if (a.r.rhoEst < r.rhoEst) { r.rhoEst = a.r.rhoEst; r.uRhoEst = a.r.uRhoEst; }
    }
    r.rhoLo = r.kappaHi > 0 ? 1.0 / r.kappaHi : 1e300;
    // Regularity.
    if (r.speedLo > 0) r.regular = Status::Pass;
    else if (sampledSpeedMin <= 1e-12 * al.Lref) r.regular = Status::Fail;
    else r.regular = Status::Unresolved;
    // Local fold: tube of radius r folds into itself where curvature radius < r.
    if (r.regular == Status::Pass && r.rhoLo >= radius) r.noFold = Status::Pass;
    else if (r.kappaWitnessLo > 1.0 / radius) r.noFold = Status::Fail;
    else r.noFold = Status::Unresolved;
    return r;
}

} // namespace knit
