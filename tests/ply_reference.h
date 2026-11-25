#pragma once
// Reference values for the ply tests: closed forms and a texel sampler that production does not
// need, against which the tests check what production builds (the shading profile, the maps).
#include <algorithm>
#include <cmath>
#include <vector>

#include "knitting/tube/plies.h"
#include "knitting/tube/ply_profile.h"

namespace knit {

// Helix angle of the ply centres to the yarn axis on a straight yarn: atan(2 pi R / pitch).
inline double plyHelixAngle(const PlyGeometry& g, double r) { return std::atan(2 * kPi * plyLayout(g, r).R / g.pitch); }

// The circle of ply 0 (centre at angle 0, distance R, radius r_p) against the radial ray at angle a.
inline RadialHit plyCircleHit(const PlyLayout& l, double a) {
    RadialHit h;
    const double sa = std::sin(a), ca = std::cos(a), d = l.rp * l.rp - l.R * l.R * sa * sa;
    if (d < 0) return h;
    const double sq = std::sqrt(d);
    h.rho = l.R * ca + sq;
    if (!(h.rho > 0)) return h;
    h.hit = true;
    h.slope = sq > 0 ? -l.R * sa - l.R * l.R * sa * ca / sq : (sa > 0 ? -1e300 : 1e300);
    return h;
}

// The n circles of the layout (R, r_p) centred at phi = 2 pi j / n: the outermost hit.
inline RadialHit circleLayoutHit(const PlyLayout& l, int n, double phi) {
    RadialHit best;
    for (int j = 0; j < std::max(1, n); ++j) {
        const RadialHit h = plyCircleHit(l, phi - 2 * kPi * j / std::max(1, n));
        if (h.hit && (!best.hit || h.rho > best.rho)) best = h;
    }
    return best;
}

// The swept plies, all n of them: the outermost hit.
inline RadialHit sweptPlyUnion(const PlyLayout& l, const PlyGeometry& g, double phi) {
    const double w = g.handedness * 2 * kPi / g.pitch;
    RadialHit best;
    for (int j = 0; j < std::max(1, g.count); ++j) {
        const RadialHit h = sweptPlyHit(l, w, phi - 2 * kPi * j / std::max(1, g.count));
        if (h.hit && (!best.hit || h.rho > best.rho)) best = h;
    }
    return best;
}

// Whether phi lies on the bridge between two plies rather than on a flank (the fold of phi to
// [0, pi / n] as in PlyShadingProfile::eval).
inline bool inBridge(const PlyShadingProfile& p, double phi) {
    if (p.g.count <= 1) return false;
    double a = std::fmod(phi, 2 * p.half);
    if (a < 0) a += 2 * p.half;
    if (a > p.half) a = 2 * p.half - a;
    return a > p.aL;
}

// Bilinear sample of a grid of PlyMaps with wrap addressing (the texel centres of PlyMaps).
inline double sampleWrap(const std::vector<float>& g, int width, int height, double u, double v, int channels = 1, int channel = 0) {
    const double fx = u * width - 0.5, fy = (1 - v) * height - 0.5;
    const double x0 = std::floor(fx), y0 = std::floor(fy), ax = fx - x0, ay = fy - y0;
    auto at = [&](double x, double y) {
        const int xi = static_cast<int>(((static_cast<long long>(x) % width) + width) % width);
        const int yi = static_cast<int>(((static_cast<long long>(y) % height) + height) % height);
        return double(g[(size_t(yi) * width + xi) * channels + channel]);
    };
    return (at(x0, y0) * (1 - ax) + at(x0 + 1, y0) * ax) * (1 - ay) + (at(x0, y0 + 1) * (1 - ax) + at(x0 + 1, y0 + 1) * ax) * ay;
}

} // namespace knit
