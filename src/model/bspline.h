#pragma once
// Uniform cubic B-spline (C2) and per-span cubic Bezier utilities.
//
// Control points Q[0..m-1] define m-3 spans. Span k uses Q[k..k+3] and covers the
// global parameter interval u in [k, k+1]. The knot at u = k is centred on Q[k+1].
// Every span converts exactly to a cubic Bezier, which gives convex-hull
// enclosures for positions and (via the hodograph) for derivatives.
#include <vector>

#include "model/vec.h"

namespace knit {

struct BezierHalves;  // defined after CubicBezier

struct QuadBezier {
    Vec3 h[3];
    Vec3 eval(double t) const;
    QuadBezier sub(double t0, double t1) const;  // restriction to [t0,t1], no rescaling of values
};

struct CubicBezier {
    Vec3 b[4];
    Vec3 eval(double t) const;
    Vec3 d1(double t) const;  // derivative with respect to the local parameter
    Vec3 d2(double t) const;
    BezierHalves split(double t) const;  // de Casteljau at t; sub(t0, t1) differs by ulps, so both stay
    CubicBezier sub(double t0, double t1) const;
    // Hodograph (derivative w.r.t. the local parameter) as a quadratic Bezier.
    QuadBezier hodograph() const;
    // Second derivative is linear: values at t=0 and t=1.
    Vec3 d2At0() const { return (b[2] - b[1] * 2.0 + b[0]) * 6.0; }
    Vec3 d2At1() const { return (b[3] - b[2] * 2.0 + b[1]) * 6.0; }
    double polygonLength() const { return norm(b[1] - b[0]) + norm(b[2] - b[1]) + norm(b[3] - b[2]); }
};

// The two halves of a cubic Bezier split at t: [0, t] and [t, 1].
struct BezierHalves {
    CubicBezier left, right;
};

// Where u lies on a B-spline: span k, and the local parameter on it.
struct SpanPosition {
    int span = 0;
    double t = 0;  // in [0, 1]
};

class BSpline {
public:
    std::vector<Vec3> ctrl;

    int numSpans() const { return ctrl.size() >= 4 ? static_cast<int>(ctrl.size()) - 3 : 0; }
    double uMax() const { return static_cast<double>(numSpans()); }
    CubicBezier span(int k) const;
    // Locate the span of u (clamped to the physical domain [0, uMax]).
    // Values outside the domain are never extrapolated.
    SpanPosition locate(double u) const;
    Vec3 eval(double u) const;
    Vec3 d1(double u) const;  // span length is 1, so d/du == d/dt
    Vec3 d2(double u) const;
};

} // namespace knit
