// The plies as shading maps on the single yarn tube (ply_profile.h).
#include "knitting/tube/ply_profile.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include "model/choices.h"

namespace knit {

namespace {

constexpr double kNoHit = -std::numeric_limits<double>::infinity();

double wrapPi(double a) {
    a = std::fmod(a + kPi, 2 * kPi);
    if (a < 0) a += 2 * kPi;
    return a - kPi;
}

// The swept ply 0 at the radial ray alpha, as a function of the ball's position t along the ply.
struct SweptRay {
    const PlyLayout& l;
    double w, alpha;
};

enum class SweptObjective {
    Discriminant,  // r_p^2 - t^2 - R^2 sin^2(alpha - w t): >= 0 iff that ball's disc meets the ray
    OuterRadius,   // R cos(alpha - w t) + sqrt(discriminant), kNoHit where the root is not real
};

double sweptValue(SweptObjective objective, const SweptRay& ray, double t) {
    const PlyLayout& l = ray.l;
    const double b = ray.alpha - ray.w * t, sb = std::sin(b), d = l.rp * l.rp - t * t - l.R * l.R * sb * sb;
    switch (objective) {
    case SweptObjective::Discriminant: return d;
    case SweptObjective::OuterRadius: return d < 0 ? kNoHit : l.R * std::cos(b) + std::sqrt(d);
    }
    return kNoHit;
}

struct SweptMaximum {
    double t = 0, value = kNoHit;
};

// Maximum of a unimodal objective on [-r_p, r_p]: a scan of `scan` points, then golden section on
// the bracket around the best sample (the value may be kNoHit outside a sub-interval).
SweptMaximum maximizeSwept(SweptObjective objective, const SweptRay& ray, int scan) {
    const double lo = -ray.l.rp, hi = ray.l.rp;
    int kb = -1;
    double fb = kNoHit;
    for (int k = 0; k <= scan; ++k) {
        const double t = lo + (hi - lo) * k / scan, v = sweptValue(objective, ray, t);
        if (v > fb) { fb = v; kb = k; }
    }
    if (kb < 0) return {0.5 * (lo + hi), kNoHit};
    double a = lo + (hi - lo) * std::max(0, kb - 1) / scan, b = lo + (hi - lo) * std::min(scan, kb + 1) / scan;
    const double gr = 0.5 * (std::sqrt(5.0) - 1);
    double c = b - gr * (b - a), d = a + gr * (b - a), fc = sweptValue(objective, ray, c), fd = sweptValue(objective, ray, d);
    for (int it = 0; it < 90 && b - a > 1e-15 * (1 + std::fabs(a)); ++it) {
        if (fc >= fd) { b = d; d = c; fd = fc; c = b - gr * (b - a); fc = sweptValue(objective, ray, c); }
        else { a = c; c = d; fc = fd; d = a + gr * (b - a); fd = sweptValue(objective, ray, d); }
    }
    double t = fc >= fd ? c : d, v = std::max(fc, fd);
    if (fb > v) { t = lo + (hi - lo) * kb / scan; v = fb; }
    return {t, v};
}

// max_t of the swept discriminant: >= 0 iff the radial ray at alpha meets ply 0's union of balls in
// the parent plane.
double sweptMaxDisc(const PlyLayout& l, double w, double alpha) {
    return maximizeSwept(SweptObjective::Discriminant, {l, w, alpha}, 64).value;
}

RadialHit circleHit(const PlyLayout& l, double a) {
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

} // namespace

RadialHit sweptPlyHit(const PlyLayout& l, double w, double alpha) {
    alpha = wrapPi(alpha);
    const SweptRay ray{l, w, alpha};
    const SweptMaximum m = maximizeSwept(SweptObjective::OuterRadius, ray, 96);
    double t = m.t, best = m.value;
    RadialHit h;
    if (!(best > 0)) return h;
    // The golden section pins t only to ~sqrt(eps) (the maximum is flat), and the slope below
    // depends on t to first order: polish with Newton on dg/dt = 0 (dg/dt analytic, its
    // derivative by central difference), keeping a step only if it does not lower g.
    auto dg = [&](double tt) {
        const double b = alpha - w * tt, sb = std::sin(b), d = l.rp * l.rp - tt * tt - l.R * l.R * sb * sb;
        return d > 0 ? l.R * w * sb + (-tt + l.R * l.R * w * sb * std::cos(b)) / std::sqrt(d) : 0.0;
    };
    for (int it = 0; it < 4; ++it) {
        const double e = 1e-6 * l.rp, g1 = dg(t), g2 = (dg(t + e) - dg(t - e)) / (2 * e);
        if (!(g2 < 0)) break;
        const double tn = std::clamp(t - g1 / g2, -l.rp, l.rp), vn = sweptValue(SweptObjective::OuterRadius, ray, tn);
        if (!(vn >= best - 1e-15 * best)) break;
        t = tn;
        best = std::max(best, vn);
        if (std::fabs(g1 / g2) < 1e-15 * l.rp) break;
    }
    h.hit = true;
    h.rho = best;
    // Envelope theorem: d rho / d alpha = dg/dalpha at the maximizing t.
    const double b = alpha - w * t, sb = std::sin(b), cb = std::cos(b);
    const double d = l.rp * l.rp - t * t - l.R * l.R * sb * sb;
    h.slope = d > 0 ? -l.R * sb - l.R * l.R * sb * cb / std::sqrt(d) : (sb > 0 ? -1e300 : 1e300);
    return h;
}

namespace {

RadialHit flankHit(const PlyShadingProfile& p, double a) {
    return p.p.model == PlyProfileModel::Circles ? circleHit(p.layout, a) : sweptPlyHit(p.layout, p.w, a);
}

double tiltDeg(const RadialHit& h) { return h.hit ? std::atan(std::fabs(h.slope) / h.rho) * 180 / kPi : 90.0; }

} // namespace

double PlyShadingProfile::phiAt(double s, double theta) const { return theta - w * s - g.phase; }

double PlyShadingProfile::eval(double phi, double* slope) const {
    if (g.count <= 1) {
        if (slope) *slope = 0;
        return r;
    }
    double a = std::fmod(phi, 2 * half);
    if (a < 0) a += 2 * half;
    double sign = 1;
    if (a > half) { a = 2 * half - a; sign = -1; }  // mirror: the neighbouring ply's flank
    double rho, sl;
    if (a <= aL) {
        const RadialHit h = flankHit(*this, a);
        rho = h.rho;
        sl = h.slope;
    } else if (a < aL + Lh) {
        const double x = (a - aL) / Lh, x2 = x * x, x3 = x2 * x;
        rho = (2 * x3 - 3 * x2 + 1) * p0 + (x3 - 2 * x2 + x) * Lh * m0 + (-2 * x3 + 3 * x2) * floorRho;
        sl = ((6 * x2 - 6 * x) * p0 + (3 * x2 - 4 * x + 1) * Lh * m0 + (-6 * x2 + 6 * x) * floorRho) / Lh;
    } else {
        rho = floorRho;
        sl = 0;
    }
    if (slope) *slope = sign * sl;
    return rho;
}

double PlyShadingProfile::surface(double phi, double* slope) const {
    if (g.count <= 1) return eval(phi, slope);
    double a = std::fmod(phi, 2 * half);
    if (a < 0) a += 2 * half;
    double sign = 1;
    if (a > half) { a = 2 * half - a; sign = -1; }
    const RadialHit h = a <= edge ? flankHit(*this, std::min(a, half * (1 - 1e-9))) : RadialHit{};
    if (!h.hit) return eval(phi, slope);
    if (slope) *slope = sign * h.slope;
    return h.rho;
}

double PlyShadingProfile::fastSurface(double phi) const {
    if (g.count <= 1 || surfaceTable.size() < 2) return r;
    double a = std::fmod(phi, 2 * half);
    if (a < 0) a += 2 * half;
    if (a > half) a = 2 * half - a;
    const double x = a / half * (surfaceTable.size() - 1);
    const size_t i = std::min(surfaceTable.size() - 2, static_cast<size_t>(x));
    const double f = x - i;
    return surfaceTable[i] * (1 - f) + surfaceTable[i + 1] * f;
}

double PlyShadingProfile::ao(double phi) const {
    if (g.count <= 1 || aoTable.size() < 2) return 1.0;
    double a = std::fmod(phi, 2 * half);
    if (a < 0) a += 2 * half;
    if (a > half) a = 2 * half - a;
    const double x = a / half * (aoTable.size() - 1);
    const size_t i = std::min(aoTable.size() - 2, static_cast<size_t>(x));
    const double f = x - i;
    return aoTable[i] * (1 - f) + aoTable[i + 1] * f;
}

double heightFieldAO(const PlyShadingProfile& prof, double phi, int rays) {
    if (prof.g.count <= 1) return 1.0;
    double sl;
    const double rho = prof.surface(phi, &sl);
    // The point at s = 0, theta = phi + phase, in a frame with the yarn axis along z.
    const double th = phi + prof.g.phase, ct = std::cos(th), st = std::sin(th);
    const Vec3 er{ct, st, 0}, et{-st, ct, 0}, ez{0, 0, 1};
    const Vec3 P = er * rho;
    const Vec3 N = normalized(er - et * (sl / rho) - ez * (-prof.w * sl));
    const Vec3 a1 = normalized(cross(N, std::fabs(N.z) < 0.9 ? ez : er)), a2 = cross(N, a1);
    const double step = prof.r / 300, tol = 1e-9 * prof.r, rOut = prof.rhoMax + tol, tMax = 2 * prof.g.pitch + 4 * prof.r;
    // Stratified cosine-weighted directions: a sqrt(rays) x sqrt(rays) grid on the unit disk
    // (Malley), cell centres.
    const int m = std::max(1, static_cast<int>(std::lround(std::sqrt(double(rays)))));
    int vis = 0;
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < m; ++j) {
            const double u1 = (i + 0.5) / m, u2 = (j + 0.5) / m, rr = std::sqrt(u1), ph = 2 * kPi * u2;
            const Vec3 d = a1 * (rr * std::cos(ph)) + a2 * (rr * std::sin(ph)) + N * std::sqrt(1 - u1);
            bool occluded = false;
            for (double t = step; t < tMax; t += step) {
                const Vec3 q = P + d * t;
                const double L = std::hypot(q.x, q.y);
                if (L > rOut && q.x * d.x + q.y * d.y >= 0) break;
                if (L < prof.fastSurface(std::atan2(q.y, q.x) - prof.w * q.z - prof.g.phase) - tol) { occluded = true; break; }
            }
            if (!occluded) ++vis;
        }
    return double(vis) / (m * m);
}

PlyShadingProfile buildPlyShadingProfile(const PlyGeometry& g, double r, const PlyProfileParams& pp) {
    PlyShadingProfile p;
    p.g = g;
    p.p = pp;
    p.r = r;
    p.layout = plyLayout(g, r);
    p.w = g.handedness * 2 * kPi / g.pitch;
    if (g.count <= 1) {
        p.half = kPi;
        p.edge = p.aL = kPi;
        p.rhoMin = p.rhoMax = p.p0 = p.floorRho = p.gapFloor = r;
        return p;
    }
    const int n = g.count;
    p.half = kPi / n;
    // End of ply 0's hit region in [0, pi / n].
    auto hits = [&](double a) {
        if (pp.model == PlyProfileModel::Circles) return circleHit(p.layout, a).hit;
        return sweptMaxDisc(p.layout, p.w, a) >= 0;
    };
    if (hits(p.half)) {
        p.edge = p.half;
    } else {
        double lo = 0, hi = p.half;
        for (int it = 0; it < 100; ++it) (hits(0.5 * (lo + hi)) ? lo : hi) = 0.5 * (lo + hi);
        p.edge = lo;
        p.gap = true;
    }
    // Bridge start: where the flank tilts blendTiltDeg, at most pi / n (1 - minBlend).
    double aL = p.half * (1 - pp.minBlend);
    const double hiA = std::min(aL, p.edge * (1 - 1e-12));
    if (tiltDeg(flankHit(p, hiA)) >= pp.blendTiltDeg) {
        double lo = 0, hi = hiA;
        for (int it = 0; it < 100; ++it) (tiltDeg(flankHit(p, 0.5 * (lo + hi))) < pp.blendTiltDeg ? lo : hi) = 0.5 * (lo + hi);
        aL = lo;
    } else {
        aL = hiA;
    }
    p.aL = aL;
    const RadialHit h0 = flankHit(p, aL);
    p.p0 = h0.rho;
    p.m0 = std::min(0.0, h0.slope);
    const double b = p.half - aL;
    p.gapFloor = pp.gapFloor < 0 ? p.layout.R * std::cos(p.half) : pp.gapFloor;
    p.floorRho = std::max(p.gapFloor, p.p0 - std::fabs(p.m0) * b / 2);
    p.floorRho = std::min(p.floorRho, p.p0 - 1e-9 * r);  // a floor above the flank start would need Lh = 0
    p.Lh = std::fabs(p.m0) > 0 ? std::min(b, 3 * (p.p0 - p.floorRho) / std::fabs(p.m0)) : b;
    // Tables (half period).
    const int N = std::max(2, pp.tableSamples);
    p.table.resize(N);
    for (int i = 0; i < N; ++i) p.table[i] = p.eval(p.half * i / (N - 1));
    p.rhoMin = *std::min_element(p.table.begin(), p.table.end());
    p.rhoMax = *std::max_element(p.table.begin(), p.table.end());
    p.surfaceTable.resize(N);
    for (int i = 0; i < N; ++i) p.surfaceTable[i] = p.surface(p.half * i / (N - 1));
    const int M = 257;
    p.aoTable.resize(M);
    for (int i = 0; i < M; ++i)
        p.aoTable[i] = heightFieldAO(p, p.half * i / (M - 1), pp.aoDirections);
    return p;
}

PlyMapSample plyMapSample(const PlyShadingProfile& prof, double uvRepeat, double u, double v) {
    PlyMapSample m;
    const double phi = prof.phiAt(u * uvRepeat, 2 * kPi * v);
    m.rho = prof.eval(phi, &m.slope);
    const double dTheta = m.slope, dS = -prof.w * m.slope;
    m.nTangent = normalized(Vec3{-dS, -dTheta / m.rho, 1.0});
    m.H = prof.rhoMax > prof.rhoMin ? (m.rho - prof.rhoMin) / (prof.rhoMax - prof.rhoMin) : 0.0;
    m.ao = prof.ao(phi);
    return m;
}

PlyMaps makePlyMaps(const PlyShadingProfile& prof, double uvRepeat, int width, int height) {
    PlyMaps m;
    m.width = width;
    m.height = height;
    m.H.resize(size_t(width) * height);
    m.ao.resize(m.H.size());
    m.normal.resize(3 * m.H.size());
    for (int y = 0; y < height; ++y)
        for (int x = 0; x < width; ++x) {
            const PlyMapSample s = plyMapSample(prof, uvRepeat, (x + 0.5) / width, 1 - (y + 0.5) / height);
            const size_t i = size_t(y) * width + x;
            m.H[i] = float(s.H);
            m.ao[i] = float(s.ao);
            m.normal[3 * i] = float(s.nTangent.x);
            m.normal[3 * i + 1] = float(s.nTangent.y);
            m.normal[3 * i + 2] = float(s.nTangent.z);
        }
    return m;
}

} // namespace knit
