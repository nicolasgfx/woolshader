#include "model/bspline.h"

#include <cassert>
#include <cmath>

namespace knit {

Vec3 QuadBezier::eval(double t) const {
    Vec3 a = lerp(h[0], h[1], t), b = lerp(h[1], h[2], t);
    return lerp(a, b, t);
}

QuadBezier QuadBezier::sub(double t0, double t1) const {
    // Blossom of a quadratic: control i of the restriction is B(t0^(2-i), t1^i).
    auto blossom = [&](double s, double t) {
        Vec3 a = lerp(h[0], h[1], s), b = lerp(h[1], h[2], s);
        return lerp(a, b, t);
    };
    return {{blossom(t0, t0), blossom(t0, t1), blossom(t1, t1)}};
}

Vec3 CubicBezier::eval(double t) const {
    Vec3 a = lerp(b[0], b[1], t), c = lerp(b[1], b[2], t), d = lerp(b[2], b[3], t);
    Vec3 e = lerp(a, c, t), f = lerp(c, d, t);
    return lerp(e, f, t);
}

Vec3 CubicBezier::d1(double t) const { return hodograph().eval(t); }

Vec3 CubicBezier::d2(double t) const { return lerp(d2At0(), d2At1(), t); }

QuadBezier CubicBezier::hodograph() const {
    return {{(b[1] - b[0]) * 3.0, (b[2] - b[1]) * 3.0, (b[3] - b[2]) * 3.0}};
}

BezierHalves CubicBezier::split(double t) const {
    Vec3 a = lerp(b[0], b[1], t), c = lerp(b[1], b[2], t), d = lerp(b[2], b[3], t);
    Vec3 e = lerp(a, c, t), f = lerp(c, d, t);
    Vec3 g = lerp(e, f, t);
    return {CubicBezier{{b[0], a, e, g}}, CubicBezier{{g, f, d, b[3]}}};
}

CubicBezier CubicBezier::sub(double t0, double t1) const {
    // Blossom of a cubic evaluated at (t0^(3-i), t1^i).
    auto blossom = [&](double s, double t, double v) {
        Vec3 a = lerp(b[0], b[1], s), c = lerp(b[1], b[2], s), d = lerp(b[2], b[3], s);
        Vec3 e = lerp(a, c, t), f = lerp(c, d, t);
        return lerp(e, f, v);
    };
    return {{blossom(t0, t0, t0), blossom(t0, t0, t1), blossom(t0, t1, t1), blossom(t1, t1, t1)}};
}

CubicBezier BSpline::span(int k) const {
    assert(k >= 0 && k < numSpans());
    const Vec3& P0 = ctrl[k];
    const Vec3& P1 = ctrl[k + 1];
    const Vec3& P2 = ctrl[k + 2];
    const Vec3& P3 = ctrl[k + 3];
    return {{(P0 + P1 * 4.0 + P2) / 6.0, (P1 * 2.0 + P2) / 3.0, (P1 + P2 * 2.0) / 3.0,
             (P1 + P2 * 4.0 + P3) / 6.0}};
}

SpanPosition BSpline::locate(double u) const {
    int n = numSpans();
    if (u <= 0.0) return {0, 0.0};
    if (u >= n) return {n - 1, 1.0};
    int k = static_cast<int>(std::floor(u));
    if (k >= n) k = n - 1;
    return {k, u - k};
}

Vec3 BSpline::eval(double u) const { const SpanPosition at = locate(u); return span(at.span).eval(at.t); }
Vec3 BSpline::d1(double u) const { const SpanPosition at = locate(u); return span(at.span).d1(at.t); }
Vec3 BSpline::d2(double u) const { const SpanPosition at = locate(u); return span(at.span).d2(at.t); }

} // namespace knit
