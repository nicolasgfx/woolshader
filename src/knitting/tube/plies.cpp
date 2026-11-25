// Explicit twisted plies inside the yarn tube (plies.h).
#include "knitting/tube/plies.h"

#include <algorithm>
#include <chrono>
#include <cmath>

namespace knit {

PlyLayout plyLayout(const PlyGeometry& g, double r) {
    PlyLayout l;
    if (g.count <= 1) {
        l.rp0 = l.rp = r;
        l.R = 0;
        return l;
    }
    const double sn = std::sin(kPi / g.count);
    l.rp0 = r * sn / (1 + sn);
    l.R = r - l.rp0;
    l.rp = g.fill * l.rp0;
    return l;
}

BSpline interpolateBSpline(const std::vector<Vec3>& P) {
    // Unknowns Q[1..m+1]; Q[1] = P[0], Q[m+1] = P[m] (from the end conditions
    // Q[0] = 2 Q[1] - Q[2], Q[m+2] = 2 Q[m+1] - Q[m]); interior rows
    // (Q[i] + 4 Q[i+1] + Q[i+2]) / 6 = P[i], i = 1..m-1: tridiagonal in Q[2..m].
    const int m = static_cast<int>(P.size()) - 1;
    BSpline b;
    b.ctrl.assign(m + 3, Vec3{0, 0, 0});
    std::vector<Vec3>& Q = b.ctrl;
    Q[1] = P[0];
    Q[m + 1] = P[m];
    if (m >= 2) {
        // Rows i = 1..m-1 for x_i = Q[i+1]: x_{i-1} + 4 x_i + x_{i+1} = 6 P[i], with the
        // known x_0 = Q[1] and x_m = Q[m+1] moved to the right-hand side (Thomas).
        const int n = m - 1;
        std::vector<double> c(n);
        std::vector<Vec3> d(n);
        for (int i = 0; i < n; ++i) {
            Vec3 rhs = P[i + 1] * 6.0;
            if (i == 0) rhs -= Q[1];
            if (i == n - 1) rhs -= Q[m + 1];
            const double denom = 4.0 - (i ? c[i - 1] : 0.0);
            c[i] = 1.0 / denom;
            d[i] = (rhs - (i ? d[i - 1] : Vec3{0, 0, 0})) / denom;
        }
        for (int i = n - 1; i >= 0; --i) Q[i + 2] = d[i] - (i + 1 < n ? Q[i + 3] * c[i] : Vec3{0, 0, 0});
    }
    Q[0] = Q[1] * 2.0 - Q[2];
    Q[m + 2] = Q[m + 1] * 2.0 - Q[m];
    return b;
}

std::pair<int, double> plyParent(const PlySet& ps, int ply, double u) {
    return {ps.yarnOf[ply], ps.parentOffset[ply] + u * ps.parentScale[ply]};
}

namespace {

// One yarn's plies: the samples run at s0 + i ds, i = 0..m. The plies start and end half a yarn radius
// inside the yarn's ends (open ends, hidden at the pole holes or inside other pieces): at the very
// end the ring frame and the curvature certificate can degenerate.
struct YarnPlies {
    const CurveSet& cs;
    const YarnSurface& ys;
    const YarnPath& yp;
    const PlyGeometry& g;
    double s0, ds;
    int m;
};

YarnPlies yarnPlies(const CurveSet& cs, const YarnSurface& ys, const PlyGeometry& g, double r) {
    const double s0 = std::min(0.5 * r, 0.25 * ys.length), L = ys.length - 2 * s0;
    const int m = std::max(1, static_cast<int>(std::ceil(L / std::min(g.pitch / g.samplesPerTwist, g.maxSampleStep))));
    return {cs, ys, cs.yarns[ys.yarnId], g, s0, L / m, m};
}

struct PlyCentres {
    std::vector<std::vector<Vec3>> points;  // per ply, per sample
    double minDistance = 1e300;             // between plies at equal s
};

// The ply centres at every sample, in parallel (each sample its own; the minimum is exact in any
// order). lam: R / r, how far out in the ring ellipse the centres sit.
PlyCentres samplePlyCentres(const YarnPlies& y, double lam, const ExecPolicy& exec) {
    const PlyGeometry& g = y.g;
    PlyCentres out;
    out.points.assign(g.count, std::vector<Vec3>(y.m + 1));
    const int64_t chunk = exec.chunk > 0 ? exec.chunk : 256;
    std::vector<double> minPart(static_cast<size_t>((int64_t(y.m) + 1 + chunk - 1) / chunk), 1e300);
    parallelFor(y.m + 1, chunk, [&](int64_t b, int64_t e, int64_t ci) {
        for (int i = static_cast<int>(b); i < static_cast<int>(e); ++i) {
            const double s = y.s0 + i * y.ds;
            for (int j = 0; j < g.count; ++j) {
                const double th = g.phase + g.handedness * 2 * kPi * s / g.pitch + 2 * kPi * j / g.count;
                out.points[j][i] = yarnCoordinates(y.cs, y.ys, s, th, lam).X;
            }
            for (int j = 0; j < g.count; ++j)
                for (int k = j + 1; k < g.count; ++k) minPart[ci] = std::min(minPart[ci], norm(out.points[j][i] - out.points[k][i]));
        }
    }, exec);
    for (double v : minPart) out.minDistance = std::min(out.minDistance, v);
    return out;
}

// Distance of a point to the yarn centreline near arc length s (Newton on the spline parameter from
// the model ring's u).
double centrelineDistance(const YarnPlies& y, const Vec3& q, double s) {
    const BSpline& C = y.yp.curve;
    const int N = static_cast<int>(y.ys.rings.size()) - 1;
    const int k = std::clamp(static_cast<int>(std::lround(s / y.ys.ds)), 0, N);
    double u = y.ys.rings[k].u;
    for (int it = 0; it < 8; ++it) {
        const Vec3 dv = C.eval(u) - q, d1 = C.d1(u), d2 = C.d2(u);
        const double gq = dot(dv, d1), hq = dot(d1, d1) + dot(dv, d2);
        const double step = hq > 0 ? -gq / hq : -gq / std::max(1e-12, dot(d1, d1));
        u = std::clamp(u + step, 0.0, C.uMax());
        if (std::fabs(step) < 1e-13) break;
    }
    return norm(q - C.eval(u));
}

// Ply j's curve through its centres, and its containment, sampled at the nodes and between them
// (u = i + 1/2): the largest |distance to the yarn centreline - R|.
double fitPly(const YarnPlies& y, const std::vector<Vec3>& centres, double R, YarnPath& ply) {
    ply.curve = interpolateBSpline(centres);
    ply.finalize();
    double containment = 0;
    for (int i = 0; i <= 2 * y.m; ++i) {
        const double u = 0.5 * i, s = y.s0 + u * y.ds;
        containment = std::max(containment, std::fabs(centrelineDistance(y, ply.curve.eval(u), s) - R));
    }
    return containment;
}

// Plies of one yarn: fitted in parallel, each its own curve and containment, appended in ply order.
void appendYarnPlies(const YarnPlies& y, double lam, const ExecPolicy& exec, PlySet& out) {
    const PlyCentres centres = samplePlyCentres(y, lam, exec);
    out.minCentreDistance = std::min(out.minCentreDistance, centres.minDistance);
    const int count = y.g.count;
    std::vector<YarnPath> plyPaths(static_cast<size_t>(count));
    std::vector<double> containment(static_cast<size_t>(count), 0.0);
    const int firstId = static_cast<int>(out.curves.yarns.size());
    const double R = out.layout.R;
    parallelFor(count, 1, [&](int64_t b, int64_t e, int64_t) {
        for (int j = static_cast<int>(b); j < static_cast<int>(e); ++j) {
            YarnPath& ply = plyPaths[j];
            ply.id = firstId + j;
            ply.name = y.yp.name + "_ply" + std::to_string(j);
            ply.course = y.yp.course;
            containment[j] = fitPly(y, centres.points[j], R, ply);
        }
    }, exec);
    for (int j = 0; j < count; ++j) {
        out.containment = std::max(out.containment, containment[j]);
        out.curves.yarns.push_back(std::move(plyPaths[j]));
        out.yarnOf.push_back(y.ys.yarnId);
        out.plyOf.push_back(j);
        out.parentScale.push_back(y.ds);
        out.parentOffset.push_back(y.s0);
    }
}

} // namespace

PlySet buildPlies(const CurveSet& cs, const SurfaceModel& sm, const PlyGeometry& g, double r, const ExecPolicy& exec) {
    const auto t0 = std::chrono::steady_clock::now();
    PlySet out;
    out.layout = plyLayout(g, r);
    out.curves.kind = CurveKind::Surface;
    out.curves.waleSpacing = cs.waleSpacing;
    out.curves.courseSpacing = cs.courseSpacing;
    out.minCentreDistance = 1e300;
    if (g.count <= 1) {
        out.curves.yarns = cs.yarns;
        for (size_t y = 0; y < cs.yarns.size(); ++y) {
            out.yarnOf.push_back(static_cast<int>(y));
            out.plyOf.push_back(0);
            out.parentScale.push_back(0);  // count 1: the plies are the yarns (no parent map)
            out.parentOffset.push_back(0);
        }
        out.minCentreDistance = 0;
        return out;
    }
    const double lam = out.layout.R / r;
    for (const YarnSurface& ys : sm.yarns) appendYarnPlies(yarnPlies(cs, ys, g, r), lam, exec, out);
    out.ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    return out;
}

std::vector<int> plyCandidateBudget(const PlySet& ps, const CurveSet& parent, double density) {
    const size_t np = ps.curves.yarns.size();
    std::vector<int> n(np, 0);
    std::vector<std::vector<int>> plies(parent.yarns.size());
    for (size_t j = 0; j < np; ++j) plies[ps.yarnOf[j]].push_back(static_cast<int>(j));
    for (size_t y = 0; y < parent.yarns.size(); ++y) {
        const std::vector<int>& J = plies[y];
        if (J.empty()) continue;
        const int N = static_cast<int>(std::floor(density * parent.yarns[y].length()));
        double total = 0;
        for (int j : J) total += ps.curves.yarns[j].length();
        if (N <= 0 || total <= 0) continue;
        // Largest remainder: floors first, then one more for the largest fractional parts.
        std::vector<std::pair<double, int>> frac;
        int used = 0;
        for (int j : J) {
            const double q = N * ps.curves.yarns[j].length() / total;
            n[j] = static_cast<int>(std::floor(q));
            used += n[j];
            frac.push_back({q - n[j], j});
        }
        std::sort(frac.begin(), frac.end(), [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
            return a.first != b.first ? a.first > b.first : a.second < b.second;
        });
        for (int k = 0; k < N - used; ++k) ++n[frac[k % frac.size()].second];
        // (Rounding of q can push a floor up by one: take it back from the smallest parts.)
        for (int k = static_cast<int>(frac.size()) - 1; used > N && k >= 0; --k)
            if (n[frac[k].second] > 0) { --n[frac[k].second]; --used; }
    }
    return n;
}

} // namespace knit
