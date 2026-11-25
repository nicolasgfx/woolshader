// Yarn tubes: rings along each yarn, display meshes, the surface index, gate (d), the final surface (yarn_surface.h).
#include "knitting/tube/yarn_surface.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>

#include "knitting/tube/compression_fit.h"
#include "model/mesh.h"

namespace knit {

namespace {

using Clock = std::chrono::steady_clock;

// Wrap an axis angle difference into (-pi/2, pi/2] (an ellipse axis is a line).
double wrapAxis(double x) {
    while (x > 0.5 * kPi) x -= kPi;
    while (x <= -0.5 * kPi) x += kPi;
    return x;
}

void axes(const Ring& r, Vec3& ea, Vec3& eb) {
    double c = std::cos(r.phi), s = std::sin(r.phi);
    ea = r.n1 * c + r.n2 * s;
    eb = r.n2 * c - r.n1 * s;
}

} // namespace

double ringSupport(const Ring& r, const Vec3& dir) {
    Vec3 ea, eb;
    axes(r, ea, eb);
    double x = r.a * dot(dir, ea), y = r.b * dot(dir, eb);
    return std::sqrt(x * x + y * y);
}

Ring ringAt(const CurveSet& cs, const YarnSurface& ys, double s) {
    const YarnPath& yp = cs.yarns[ys.yarnId];
    const int N = static_cast<int>(ys.rings.size()) - 1;
    s = std::clamp(s, 0.0, ys.length);
    int k = std::min(N - 1, std::max(0, static_cast<int>(std::floor(s / ys.ds))));
    double f = std::clamp((s - ys.rings[k].s) / ys.ds, 0.0, 1.0);
    const Ring &A = ys.rings[k], &B = ys.rings[k + 1];
    Ring r;
    r.s = s;
    r.u = yp.paramAtArcLength(s);
    r.c = yp.curve.eval(r.u);
    r.T = normalized(yp.curve.d1(r.u));
    Vec3 n = A.n1 * (1 - f) + B.n1 * f;
    r.n1 = normalized(n - r.T * dot(n, r.T));
    r.n2 = cross(r.T, r.n1);
    r.a = A.a * (1 - f) + B.a * f;
    r.b = A.b * (1 - f) + B.b * f;
    r.phi = A.phi + f * wrapAxis(B.phi - A.phi);
    return r;
}

YarnSurfacePoint evalSurface(const CurveSet& cs, const YarnSurface& ys, double s, double theta) {
    YarnSurfacePoint sp;
    sp.ring = ringAt(cs, ys, s);
    const Ring& r = sp.ring;
    Vec3 ea, eb;
    axes(r, ea, eb);
    double ph = theta - r.phi;
    double c = std::cos(ph), sn = std::sin(ph);
    sp.p = r.c + ea * (r.a * c) + eb * (r.b * sn);
    sp.n = normalized(ea * (c / r.a) + eb * (sn / r.b));
    return sp;
}

YarnCoordPoint yarnCoordinates(const CurveSet& cs, const YarnSurface& ys, double s, double theta, double lambda) {
    const YarnPath& yp = cs.yarns[ys.yarnId];
    const int N = static_cast<int>(ys.rings.size()) - 1;
    s = std::clamp(s, 0.0, ys.length);
    const int k = std::min(N - 1, std::max(0, static_cast<int>(std::floor(s / ys.ds))));
    const double f = std::clamp((s - ys.rings[k].s) / ys.ds, 0.0, 1.0);
    const Ring &A = ys.rings[k], &B = ys.rings[k + 1];
    // Centre and tangent on the spline. The curve parameter is interpolated linearly
    // between the model rings' exact u (no arc-length inversion per call; c(s) stays
    // within the interpolation error of the exact arc-length point, see the tests), so
    // dc/ds = d1 du/ds and dT/ds = (d2 - T (T.d2)) du/ds / |d1|.
    const double dus = (B.u - A.u) / ys.ds, u = A.u + f * (B.u - A.u);
    const Vec3 d1 = yp.curve.d1(u), d2 = yp.curve.d2(u);
    const Vec3 c = yp.curve.eval(u);
    const double sp = norm(d1);
    const Vec3 T = d1 / sp, dT = (d2 - T * dot(d2, T)) * (dus / sp), dc = d1 * dus;
    // n1 = normalize(m - T (m.T)) with m interpolated linearly between the rings.
    const Vec3 m = A.n1 * (1 - f) + B.n1 * f, dm = (B.n1 - A.n1) / ys.ds;
    const double mT = dot(m, T);
    const Vec3 w = m - T * mT, dw = dm - dT * mT - T * (dot(dm, T) + dot(m, dT));
    const double wl = norm(w);
    const Vec3 n1 = w / wl, dn1 = (dw - n1 * dot(n1, dw)) / wl;
    const Vec3 n2 = cross(T, n1), dn2 = cross(dT, n1) + cross(T, dn1);
    const double a = A.a * (1 - f) + B.a * f, b = A.b * (1 - f) + B.b * f;
    const double da = (B.a - A.a) / ys.ds, db = (B.b - A.b) / ys.ds;
    const double dphi = wrapAxis(B.phi - A.phi) / ys.ds, phi = A.phi + f * wrapAxis(B.phi - A.phi);
    const double cp = std::cos(phi), sph = std::sin(phi);
    const Vec3 ea = n1 * cp + n2 * sph, eb = n2 * cp - n1 * sph;
    const Vec3 dea = dn1 * cp + dn2 * sph + eb * dphi, deb = dn2 * cp - dn1 * sph - ea * dphi;
    const double ph = theta - phi, co = std::cos(ph), si = std::sin(ph);
    YarnCoordPoint out;
    out.dLambda = ea * (a * co) + eb * (b * si);
    out.X = c + out.dLambda * lambda;
    out.dTheta = (eb * (b * co) - ea * (a * si)) * lambda;
    // d/ds of cos(theta - phi) is sin(ph) dphi, of sin(theta - phi) is -cos(ph) dphi.
    out.dS = dc + (dea * (a * co) + ea * (da * co + a * si * dphi) + deb * (b * si) + eb * (db * si - b * co * dphi)) * lambda;
    return out;
}

Vec3 yarnPosition(const CurveSet& cs, const YarnSurface& ys, double s, double theta, double lambda) {
    // yarnCoordinates' expressions for X, in its order.
    const YarnPath& yp = cs.yarns[ys.yarnId];
    const int N = static_cast<int>(ys.rings.size()) - 1;
    s = std::clamp(s, 0.0, ys.length);
    const int k = std::min(N - 1, std::max(0, static_cast<int>(std::floor(s / ys.ds))));
    const double f = std::clamp((s - ys.rings[k].s) / ys.ds, 0.0, 1.0);
    const Ring &A = ys.rings[k], &B = ys.rings[k + 1];
    const double u = A.u + f * (B.u - A.u);
    const Vec3 d1 = yp.curve.d1(u);
    const Vec3 c = yp.curve.eval(u);
    const double sp = norm(d1);
    const Vec3 T = d1 / sp;
    const Vec3 m = A.n1 * (1 - f) + B.n1 * f;
    const double mT = dot(m, T);
    const Vec3 w = m - T * mT;
    const double wl = norm(w);
    const Vec3 n1 = w / wl;
    const Vec3 n2 = cross(T, n1);
    const double a = A.a * (1 - f) + B.a * f, b = A.b * (1 - f) + B.b * f;
    const double phi = A.phi + f * wrapAxis(B.phi - A.phi);
    const double cp = std::cos(phi), sph = std::sin(phi);
    const Vec3 ea = n1 * cp + n2 * sph, eb = n2 * cp - n1 * sph;
    const double ph = theta - phi, co = std::cos(ph), si = std::sin(ph);
    const Vec3 dLambda = ea * (a * co) + eb * (b * si);
    return c + dLambda * lambda;
}

namespace {

// One yarn's rings at the fixed spacing: position and tangent, round at radius r. Rings are
// independent here, so they run in parallel.
YarnSurface sampleRings(const YarnPath& yp, int yarnId, double r, double ringSpacing, const ExecPolicy& exec) {
    YarnSurface ys;
    ys.yarnId = yarnId;
    ys.radius = r;
    ys.length = yp.length();
    const int N = std::max(1, static_cast<int>(std::ceil(ys.length / ringSpacing)));
    ys.ds = ys.length / N;
    ys.rings.resize(N + 1);
    parallelFor(N + 1, exec.chunk > 0 ? exec.chunk : 512, [&](int64_t kBegin, int64_t kEnd, int64_t) {
        for (int64_t k = kBegin; k < kEnd; ++k) {
            Ring& R = ys.rings[k];
            R.s = k * ys.ds;
            R.u = yp.paramAtArcLength(R.s);
            R.c = yp.curve.eval(R.u);
            R.T = normalized(yp.curve.d1(R.u));
            R.a = R.b = r * yp.thicknessAt(R.u);
        }
    }, exec);
    return ys;
}

// The rotation-minimizing frame along the rings (double reflection, Wang et al. 2008). Each ring's
// frame depends on the previous one, so this runs in ring order.
void transportFrames(YarnSurface& ys) {
    for (size_t k = 0; k < ys.rings.size(); ++k) {
        Ring& R = ys.rings[k];
        if (k == 0) {
            // Start the frame from the fabric normal (front) projected off the tangent.
            Vec3 n = Vec3{0, 0, 1} - R.T * R.T.z;
            if (norm(n) < 1e-6) n = Vec3{1, 0, 0} - R.T * R.T.x;
            R.n1 = normalized(n);
        } else {
            const Ring& P = ys.rings[k - 1];
            Vec3 v1 = R.c - P.c;
            double c1 = dot(v1, v1);
            Vec3 rL = P.n1, tL = P.T;
            if (c1 > 0) {
                rL = P.n1 - v1 * (2 / c1 * dot(v1, P.n1));
                tL = P.T - v1 * (2 / c1 * dot(v1, P.T));
            }
            Vec3 v2 = R.T - tL;
            double c2 = dot(v2, v2);
            R.n1 = c2 > 0 ? rL - v2 * (2 / c2 * dot(v2, rL)) : rL;
            R.n1 = normalized(R.n1 - R.T * dot(R.n1, R.T));
        }
        R.n2 = cross(R.T, R.n1);
    }
}

// What one ring's cross-section must respect: up to kMaxReq contact directions and the fold
// requirement (compression_fit.h).
struct RingRequirements {
    cfit::RequirementSet set;
    double dropped = 1e300;  // the smallest h a full set gave up (diagnostics)
    void add(double h, double psi) {
        double d = 1e300;
        cfit::insert(set, h, psi, d);
        dropped = std::min(dropped, d);
    }
};

// The fold requirement: the support along the principal normal stays below rho.
void addFoldRequirement(const BSpline& C, const Ring& R, double r, RingRequirements& rq) {
    const Vec3 d1 = C.d1(R.u), d2 = C.d2(R.u);
    const double spd = norm(d1);
    if (spd <= 0) return;
    const Vec3 T = d1 / spd, kv = (d2 - T * dot(d2, T)) / (spd * spd);
    const double kappa = norm(kv);
    if (kappa <= 1e-12) return;
    const double h = 1.0 / (kappa * (1 + cfit::kFoldMargin));
    if (h < r) rq.add(h, std::atan2(dot(kv, R.n2), dot(kv, R.n1)));
}

// The spline parameter closest to p, by Newton from u.
double closestSplineParam(const BSpline& C, double u, const Vec3& p) {
    for (int it = 0; it < 8; ++it) {
        Vec3 dv = C.eval(u) - p, d1 = C.d1(u), d2 = C.d2(u);
        double g = dot(dv, d1), h = dot(d1, d1) + dot(dv, d2);
        double step = h > 0 ? -g / h : -g / std::max(1e-12, dot(d1, d1));
        u = std::clamp(u + step, 0.0, C.uMax());
        if (std::fabs(step) < 1e-12) break;
    }
    return u;
}

// What the contact search of one model needs. window: the own-yarn arc window, half a turn at the
// maximum extent; portions of the same yarn closer than that along the yarn are not contacts.
struct ContactSearch {
    const CurveSet& cs;
    const SurfaceModel& sm;
    const SurfaceIndex& index;
    double window;
    const ExecPolicy& exec;
};

// Per ring, the requirement from the actual centreline distance to every nearby tube portion (other
// yarns, and portions of the same yarn beyond the window). Both sides of a contact compute their own
// requirement h = d / 2, so the flattened supports meet. Each ring writes only its own requirements.
void collectContactRequirements(const ContactSearch& search, int y, std::vector<RingRequirements>& req) {
    const YarnSurface& ys = search.sm.yarns[y];
    const YarnPath& yp = search.cs.yarns[y];
    const double W = search.window;
    req.assign(ys.rings.size(), RingRequirements{});
    const ExecPolicy& exec = search.exec;
    parallelFor(static_cast<int64_t>(ys.rings.size()), exec.chunk > 0 ? exec.chunk : 256, [&](int64_t kBegin, int64_t kEnd, int64_t) {
        std::vector<SurfaceIndex::NearbyAxis> nearby;
        for (size_t k = static_cast<size_t>(kBegin); k < static_cast<size_t>(kEnd); ++k) {
            const Ring& R = ys.rings[k];
            const double r = R.a, rad = 2 * r;  // round here: this ring's own radius
            // The own window in the spline parameter (as the device's requirement search): Newton
            // started on a distant portion of the own yarn can converge back into the window, onto
            // a point near the ring (a spurious requirement far below the lower bound).
            const double uOwnLo = yp.paramAtArcLength(std::max(0.0, R.s - W));
            const double uOwnHi = yp.paramAtArcLength(std::min(ys.length, R.s + W));
            addFoldRequirement(yp.curve, R, r, req[k]);
            search.index.nearbyAxes(R.c, rad, {y, R.s - W, R.s + W}, nearby);
            for (const SurfaceIndex::NearbyAxis& axis : nearby) {
                // The capsule axis is a chord of the centreline (sag ~ 2 um); refine the closest
                // point on the actual spline (Newton from the nearest model ring).
                const YarnSurface& ys2 = search.sm.yarns[axis.yarn];
                const BSpline& C2 = search.cs.yarns[axis.yarn].curve;
                const int k2 = std::clamp(static_cast<int>(std::round(axis.s / ys2.ds)), 0, static_cast<int>(ys2.rings.size()) - 1);
                const double u = closestSplineParam(C2, ys2.rings[k2].u, R.c);
                const Vec3 q = C2.eval(u);
                const double d = norm(q - R.c);
                if (d >= rad || d <= 1e-12 || (axis.yarn == y && u >= uOwnLo && u <= uOwnHi)) continue;
                Vec3 dir = q - R.c;
                dir = dir - R.T * dot(dir, R.T);
                if (norm(dir) < 1e-12) continue;
                req[k].add(0.5 * d, std::atan2(dot(dir, R.n2), dot(dir, R.n1)));
            }
        }
    }, exec);
}

// A max-filter of the flattening deficit (r - h) along s, so blending can only flatten more, never
// leave an overlap.
void smoothRequirements(const YarnSurface& ys, double sigma, std::vector<RingRequirements>& req) {
    const int half = static_cast<int>(std::ceil(3 * sigma / ys.ds));
    std::vector<RingRequirements> out = req;
    for (int k = 0; k < int(ys.rings.size()); ++k) {
        const double r = ys.rings[k].a;
        for (int j = std::max(0, k - half); j <= std::min(int(ys.rings.size()) - 1, k + half); ++j) {
            const RingRequirements& rj = req[j];
            if (rj.set.n == 0) continue;
            double w = std::exp(-0.5 * std::pow((j - k) * ys.ds / sigma, 2));
            for (int i = 0; i < rj.set.n; ++i) {
                const double h = cfit::smoothed(rj.set.e[i].h, w, r);
                if (h >= r) continue;
                out[k].add(h, rj.set.e[i].psi);
            }
        }
    }
    req = std::move(out);
}

// Each ring's ellipse from its requirements; lowerRatio is the smallest allowed semi-axis over the
// ring's round radius.
void fitSections(const std::vector<std::vector<RingRequirements>>& req, double lowerRatio, SurfaceModel& sm) {
    CompressionReport& rep = sm.compression;
    for (size_t y = 0; y < sm.yarns.size(); ++y) {
        for (size_t k = 0; k < sm.yarns[y].rings.size(); ++k) {
            const RingRequirements& rq = req[y][k];
            Ring& R = sm.yarns[y].rings[k];
            const double r = R.a;
            const cfit::Section f = cfit::fit(rq.set, r, lowerRatio * r);
            if (rq.dropped < r) ++rep.dropped;
            if (!f.flattened) continue;
            if (f.clamped) ++rep.clamped;
            R.b = f.b;
            R.phi = f.phi;  // a-axis perpendicular to the primary contact
            R.a = f.a;
        }
    }
}

} // namespace

SurfaceModel buildSurface(const CurveSet& cs, const PatternParams& pp, const SurfaceParams& sp, const ExecPolicy& exec) {
    auto t0 = Clock::now();
    SurfaceModel sm;
    const double r = pp.radius();
    for (size_t y = 0; y < cs.yarns.size(); ++y) {
        YarnSurface ys = sampleRings(cs.yarns[y], static_cast<int>(y), r, sp.ringSpacing, exec);
        transportFrames(ys);
        sm.yarns.push_back(std::move(ys));
    }
    if (sp.compression) {
        const SurfaceIndex index(cs, sm);  // circular rings at this point
        const ContactSearch search{cs, sm, index, kPi * (1 + sp.maxExpansion) * r, exec};
        std::vector<std::vector<RingRequirements>> req(sm.yarns.size());
        for (size_t y = 0; y < sm.yarns.size(); ++y) {
            collectContactRequirements(search, static_cast<int>(y), req[y]);
            if (sp.smoothing > 0) smoothRequirements(sm.yarns[y], sp.smoothing, req[y]);
        }
        fitSections(req, 1 - sp.maxCompression, sm);
    }
    CompressionReport& rep = sm.compression;
    for (const YarnSurface& ys : sm.yarns)
        for (const Ring& R : ys.rings) {
            rep.minAxisRatio = std::min(rep.minAxisRatio, R.b / r);
            rep.maxAreaChange = std::max(rep.maxAreaChange, std::fabs(R.a * R.b / (r * r) - 1));
        }
    sm.ms = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
    return sm;
}

namespace {

// The tube's u at arc length s before the uvWrap offset (yarnUvOf).
double tubeU(const YarnSurface& ys, const TubeParams& tp, double s) {
    const YarnUv y = yarnUvOf(tp, ys.yarnId);
    return s / tp.uvRepeat / y.pitchFactor + y.phase;
}

// One tube ring's place in the mesh: at arc length s, u = s / uvRepeat - uOff0 (uvWrap), emitted once
// or, where u wraps, twice; `first` is its first vertex.
struct TubeRing {
    double s = 0, uOff0 = 0;
    int copies = 1;
    uint32_t first = 0;
};

// The rings in order, from vertex `first` on.
std::vector<TubeRing> planTubeRings(const YarnSurface& ys, const TubeParams& tp, int N, uint32_t stride, uint32_t first) {
    std::vector<TubeRing> plan(static_cast<size_t>(N) + 1);
    double uOff = 0;
    uint32_t at = first;
    for (int j = 0; j <= N; ++j) {
        const double sj = ys.length * j / N;
        const bool wrapHere = tp.uvWrap > 0 && j > 0 && j < N && tubeU(ys, tp, sj) - uOff >= tp.uvWrap;
        plan[j] = {sj, uOff, wrapHere ? 2 : 1, at};
        if (wrapHere) uOff += tp.uvWrap;
        at += plan[j].copies * stride;
    }
    return plan;
}

void resizeTubeVertices(TubeMesh& m, size_t total) {
    m.pos.resize(total);
    m.nrm.resize(total);
    m.uv.resize(2 * total);
    m.s.resize(total);
    m.t.resize(total);
    m.theta.resize(total);
    m.yarn.resize(total);
    m.axisRatio.resize(total);
}

// The vertices of one ring: every copy, seam duplicated so that v runs 0..1.
void writeTubeRing(const CurveSet& cs, const YarnSurface& ys, const TubeRing& plan, const TubeParams& tp, TubeMesh& m) {
    const int sides = std::max(3, tp.sides);
    const uint32_t stride = sides + 1;
    for (int copy = 0; copy < plan.copies; ++copy) {
        const double uOff = copy ? plan.uOff0 + tp.uvWrap : plan.uOff0;
        size_t v = plan.first + static_cast<size_t>(copy) * stride;
        const double s = plan.s;
        Ring R = ringAt(cs, ys, s);
        Vec3 ea, eb;
        axes(R, ea, eb);
        for (int i = 0; i <= sides; ++i) {
            double theta = 2 * kPi * i / sides;
            double ph = theta - R.phi, c = std::cos(ph), sn = std::sin(ph);
            m.pos[v] = R.c + ea * (R.a * c) + eb * (R.b * sn);
            m.nrm[v] = normalized(ea * (c / R.a) + eb * (sn / R.b));
            m.uv[2 * v] = float(tubeU(ys, tp, s) - uOff);
            m.uv[2 * v + 1] = float(double(i) / sides);
            m.s[v] = s;
            m.t[v] = ys.length > 0 ? s / ys.length : 0;
            m.theta[v] = theta;
            m.yarn[v] = ys.yarnId;
            m.axisRatio[v] = float(R.b / ys.radius);
            ++v;
        }
    }
}

} // namespace

namespace {

// A uniform in [0, 1) from a yarn's index and a channel (splitmix64).
double yarnUniform(int yarnIndex, uint64_t channel) {
    uint64_t z = static_cast<uint64_t>(static_cast<uint32_t>(yarnIndex)) * 0x9E3779B97F4A7C15ull + channel * 0xBF58476D1CE4E5B9ull + 1;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    z ^= z >> 31;
    return static_cast<double>(z >> 11) * (1.0 / 9007199254740992.0);
}

} // namespace

YarnUv yarnUvOf(const TubeParams& tp, int yarnIndex) {
    YarnUv y;
    if (tp.uvPitchJitter > 0) y.pitchFactor = 1 + tp.uvPitchJitter * (2 * yarnUniform(yarnIndex, 1) - 1);
    if (tp.uvPhaseJitter > 0) y.phase = tp.uvPhaseJitter * yarnUniform(yarnIndex, 2);
    return y;
}

TubeMesh buildTubeMesh(const CurveSet& cs, const SurfaceModel& sm, const TubeParams& tp, const ExecPolicy& exec) {
    auto t0 = Clock::now();
    TubeMesh m;
    const int sides = std::max(3, tp.sides);
    const uint32_t stride = sides + 1;
    for (const YarnSurface& ys : sm.yarns) {
        m.yarnFirstVertex.push_back(static_cast<uint32_t>(m.pos.size()));
        const int N = std::max(1, static_cast<int>(std::ceil(ys.length / tp.ringSpacing)));
        // The ring plan in order, then the rings' vertices in parallel at those offsets, then the
        // strips in order.
        const std::vector<TubeRing> plan = planTubeRings(ys, tp, N, stride, static_cast<uint32_t>(m.pos.size()));
        resizeTubeVertices(m, plan.back().first + static_cast<size_t>(plan.back().copies) * stride);
        parallelFor(N + 1, exec.chunk > 0 ? exec.chunk : 64, [&](int64_t jBegin, int64_t jEnd, int64_t) {
            for (int j = static_cast<int>(jBegin); j < static_cast<int>(jEnd); ++j) writeTubeRing(cs, ys, plan[j], tp, m);
        }, exec);
        for (int j = 1; j <= N; ++j) {
            // The strip from the previous ring's last copy ends at this ring's first copy (same u offset).
            const uint32_t prev = plan[j - 1].first + (plan[j - 1].copies - 1) * stride, ringFirst = plan[j].first;
            for (int i = 0; i < sides; ++i) {
                uint32_t a = prev + i, b = a + 1, c = ringFirst + i, d = c + 1;
                m.tri.insert(m.tri.end(), {a, c, b, b, c, d});
            }
        }
    }
    m.yarnFirstVertex.push_back(static_cast<uint32_t>(m.pos.size()));
    m.ms = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
    return m;
}

// ------------------------------------------------------------------ index ----

SurfaceIndex::SurfaceIndex(const CurveSet& cs, const SurfaceModel& sm) : cs_(cs), sm_(sm) {
    zTop_ = -1e300;
    zBottom_ = 1e300;
    for (const YarnSurface& ys : sm.yarns) {
        const int N = static_cast<int>(ys.rings.size()) - 1;
        const int step = std::max(1, static_cast<int>(std::round(0.1 / ys.ds)));  // ~0.1 mm capsules
        for (int k = 0; k < N; k += step) {
            int k1 = std::min(N, k + step);
            const Vec3 a = ys.rings[k].c, b = ys.rings[k1].c;
            double R = 0, dev = 0;
            for (int j = k; j <= k1; ++j) {
                R = std::max({R, ys.rings[j].a, ys.rings[j].b});
                dev = std::max(dev, closestSegmentPoints(ys.rings[j].c, ys.rings[j].c, a, b).distance);
            }
            // Deviation of the centreline from the capsule axis (measured at the model
            // rings) plus a margin for the curve between rings.
            R += dev + 1e-3;
            Capsule c{a, b, R, ys.rings[k].s, ys.rings[k1].s, ys.yarnId};
            caps_.push_back(c);
            zTop_ = std::max({zTop_, a.z + R, b.z + R});
            zBottom_ = std::min({zBottom_, a.z - R, b.z - R});
        }
    }
    double maxR = 0;
    for (const Capsule& c : caps_) maxR = std::max(maxR, c.R);
    hash_ = SpatialHash3(1.0 + maxR);
    for (const Capsule& c : caps_) {
        Vec3 lo{std::min(c.a.x, c.b.x) - c.R, std::min(c.a.y, c.b.y) - c.R, std::min(c.a.z, c.b.z) - c.R};
        Vec3 hi{std::max(c.a.x, c.b.x) + c.R, std::max(c.a.y, c.b.y) + c.R, std::max(c.a.z, c.b.z) + c.R};
        hash_.insert(lo, hi);
    }
    hash_.finalize();
}

bool SurfaceIndex::pointInEnclosure(const Vec3& q, double inflate, const Exclude& ex) const {
    bool hit = false;
    hash_.query(q - Vec3{inflate, inflate, inflate}, q + Vec3{inflate, inflate, inflate}, [&](int i) {
        if (hit) return;
        const Capsule& c = caps_[i];
        const ClosestSegmentPoints closestPair = closestSegmentPoints(q, q, c.a, c.b);
        hit = !ex.covers(c.yarn, c.s0 + closestPair.t * (c.s1 - c.s0)) && closestPair.distance < c.R + inflate;
    });
    return hit;
}

bool SurfaceIndex::segmentHitsEnclosure(const Vec3& a, const Vec3& b, double inflate, const Exclude& ex) const {
    bool hit = false;
    const Vec3 lo = Vec3{std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z)} - Vec3{inflate, inflate, inflate};
    const Vec3 hi = Vec3{std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)} + Vec3{inflate, inflate, inflate};
    hash_.query(lo, hi, [&](int i) {
        if (hit) return;
        const Capsule& c = caps_[i];
        const ClosestSegmentPoints closestPair = closestSegmentPoints(a, b, c.a, c.b);
        hit = !ex.covers(c.yarn, c.s0 + closestPair.t * (c.s1 - c.s0)) && closestPair.distance < c.R + inflate;
    });
    return hit;
}

void SurfaceIndex::nearbyAxes(const Vec3& q, double radius, const Exclude& ex, std::vector<NearbyAxis>& out) const {
    out.clear();
    hash_.query(q - Vec3{radius, radius, radius}, q + Vec3{radius, radius, radius}, [&](int i) {
        const Capsule& c = caps_[i];
        const ClosestSegmentPoints closestPair = closestSegmentPoints(q, q, c.a, c.b);
        const double sc = c.s0 + closestPair.t * (c.s1 - c.s0);
        if (closestPair.distance < radius && !ex.covers(c.yarn, sc)) out.push_back({c.yarn, sc});
    });
}

bool SurfaceIndex::ringInside(const Capsule& c, const Vec3& q, double inflate, Hit& hit) const {
    const YarnSurface& ys = sm_.yarns[c.yarn];
    double sEst = c.s0 + closestSegmentPoints(q, q, c.a, c.b).t * (c.s1 - c.s0);
    const int N = static_cast<int>(ys.rings.size()) - 1;
    int k0 = std::clamp(static_cast<int>(std::round(sEst / ys.ds)), 0, N);
    int best = -1;
    double bd = 1e300;
    for (int k = std::max(0, k0 - 6); k <= std::min(N, k0 + 6); ++k) {
        double d = std::fabs(dot(q - ys.rings[k].c, ys.rings[k].T));
        if (d < bd) { bd = d; best = k; }
    }
    if (best < 0 || bd > ys.ds) return false;  // not within the tube's extent (e.g. past an open end)
    const Ring& R = ys.rings[best];
    Vec3 v = q - R.c;
    v = v - R.T * dot(v, R.T);
    Vec3 ea, eb;
    axes(R, ea, eb);
    double A = R.a + inflate, B = R.b + inflate;
    double x = dot(v, ea) / A, y = dot(v, eb) / B;
    double qq = x * x + y * y;
    if (qq >= 1.0) return false;
    // The gradient of x^2 + y^2 in the ring plane; on the axis itself any direction is out: ea.
    const Vec3 grad = ea * (x / A) + eb * (y / B);
    hit = {true, c.yarn, R.s, (1 - std::sqrt(qq)) * std::min(A, B), norm(grad) > 1e-12 ? normalized(grad) : ea};
    return true;
}

SurfaceIndex::Hit SurfaceIndex::exactInside(const Vec3& q, double inflate, const Exclude& only, const Exclude& ex) const {
    Hit h;
    hash_.query(q - Vec3{inflate, inflate, inflate}, q + Vec3{inflate, inflate, inflate}, [&](int i) {
        const Capsule& c = caps_[i];
        if (only.yarn >= 0 && c.yarn != only.yarn) return;
        if (closestSegmentPoints(q, q, c.a, c.b).distance >= c.R + inflate) return;
        Hit ring;
        if (!ringInside(c, q, inflate, ring)) return;
        if (only.yarn >= 0 && !only.covers(c.yarn, ring.s)) return;
        if (ex.covers(c.yarn, ring.s)) return;
        if (!h.inside || ring.depth > h.depth) h = ring;
    });
    return h;
}

// ------------------------------------------------------------ surface check ----

namespace {

// Local folds: the sweep of a cross-section is locally injective where support(kappa_hat) * kappa < 1.
void checkFolds(const CurveSet& cs, const SurfaceModel& sm, SurfaceCheck& out) {
    for (const YarnSurface& ys : sm.yarns) {
        const BSpline& C = cs.yarns[ys.yarnId].curve;
        for (const Ring& R : ys.rings) {
            Vec3 d1 = C.d1(R.u), d2 = C.d2(R.u);
            double sp = norm(d1);
            if (sp <= 0) continue;
            Vec3 T = d1 / sp;
            Vec3 kv = (d2 - T * dot(d2, T)) / (sp * sp);
            double kappa = norm(kv);
            if (kappa < 1e-12) continue;
            double h = ringSupport(R, kv / kappa);
            double ratio = 1.0 / (kappa * h);
            if (ratio < 1.0) ++out.foldRings;
            if (ratio < out.minFoldRatio) { out.minFoldRatio = ratio; out.worstFoldYarn = ys.yarnId; out.worstFoldS = R.s; }
        }
    }
    out.fold = out.minFoldRatio < 1.0 ? Status::Fail : Status::Pass;
}

// The surface points of one ring, `angles` around, against every tube outside the local window own.
void checkRingPoints(const SurfaceIndex& index, const Ring& R, const SurfaceIndex::Exclude& own, const SurfaceCheckParams& p, SurfaceCheck& out) {
    Vec3 ea, eb;
    axes(R, ea, eb);
    for (int i = 0; i < p.angles; ++i) {
        double theta = 2 * kPi * i / p.angles, ph = theta - R.phi;
        double cosPh = std::cos(ph), sinPh = std::sin(ph);
        const Vec3 q = R.c + ea * (R.a * cosPh) + eb * (R.b * sinPh);
        ++out.pointsTested;
        SurfaceIndex::Hit h = index.exactInside(q, 0.0, {}, own);
        if (!h.inside || h.depth <= p.depthTolerance) continue;
        ++out.penetratingPoints;
        if (h.depth > out.maxDepth) {
            out.maxDepth = h.depth;
            out.worstYarnA = own.yarn; out.worstYarnB = h.yarn;
            out.worstSA = R.s; out.worstSB = h.s;
            out.worstPoint = q;
        }
    }
}

// One yarn's sampled rings in parallel: each chunk accumulates its own counters and worst point,
// merged in ring order (strict comparisons: the first worst point, as in one pass).
void checkYarnCollisions(const SurfaceIndex& index, const YarnSurface& ys, const SurfaceCheckParams& p, const ExecPolicy& exec, SurfaceCheck& out) {
    double ext = 0;
    for (const Ring& R : ys.rings) ext = std::max({ext, R.a, R.b});
    const double W = kPi * ext;  // same-yarn local window (half a turn at the extent)
    const size_t stride = static_cast<size_t>(std::max(1, p.ringStride));
    const int64_t nSampled = static_cast<int64_t>((ys.rings.size() + stride - 1) / stride);
    const int64_t chunk = exec.chunk > 0 ? exec.chunk : 64;
    std::vector<SurfaceCheck> parts(static_cast<size_t>((nSampled + chunk - 1) / chunk));
    parallelFor(nSampled, chunk, [&](int64_t iBegin, int64_t iEnd, int64_t pi) {
        for (int64_t si = iBegin; si < iEnd; ++si) {
            const Ring& R = ys.rings[static_cast<size_t>(si) * stride];
            checkRingPoints(index, R, {ys.yarnId, R.s - W, R.s + W}, p, parts[pi]);
        }
    }, exec);
    for (const SurfaceCheck& a : parts) {
        out.pointsTested += a.pointsTested;
        out.penetratingPoints += a.penetratingPoints;
        if (a.maxDepth > out.maxDepth) {
            out.maxDepth = a.maxDepth;
            out.worstYarnA = a.worstYarnA; out.worstYarnB = a.worstYarnB;
            out.worstSA = a.worstSA; out.worstSB = a.worstSB;
            out.worstPoint = a.worstPoint;
        }
    }
}

void checkCollisions(const CurveSet& cs, const SurfaceModel& sm, const SurfaceCheckParams& p, const ExecPolicy& exec, SurfaceCheck& out) {
    const SurfaceIndex index(cs, sm);
    for (const YarnSurface& ys : sm.yarns) checkYarnCollisions(index, ys, p, exec, out);
    out.collisions = out.penetratingPoints > 0 ? Status::Fail : Status::Pass;
}

} // namespace

SurfaceCheck checkSurface(const CurveSet& cs, const SurfaceModel& sm, const SurfaceCheckParams& p, const ExecPolicy& exec) {
    auto t0 = Clock::now();
    SurfaceCheck out;
    checkFolds(cs, sm, out);
    checkCollisions(cs, sm, p, exec, out);
    char buf[160];
    std::snprintf(buf, sizeof buf, "sampled: every %d model rings (%.3f mm), %d angles, depth tolerance %.0e mm",
                  p.ringStride, sm.yarns.empty() ? 0.0 : p.ringStride * sm.yarns[0].ds, p.angles, p.depthTolerance);
    out.resolution = buf;
    out.ms = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
    return out;
}

// ------------------------------------------------------- final surface spec ----

double FinalSurfaceSpec::supportRadius() const {
    if (plies > 1) return 0.5 * renderedDiameter();
    const double r = 0.5 * renderedDiameter();
    return compression ? (1 + SurfaceParams{}.maxExpansion) * r : r;
}

PatternParams FinalSurfaceSpec::rendered(const PatternParams& centreline) const {
    PatternParams p = centreline;
    p.yarnDiameter = renderedDiameter();
    return p;
}

bool FinalSurfaceSpec::valid(std::string* why) const {
    auto no = [&](const char* w) {
        if (why) *why = w;
        return false;
    };
    if (inflate < 0) return no("inflation e must be >= 0");
    if (inflate > 0 && plies > 1 && !compressedYarn) return no("contact flattening (inflation) with explicit plies is not supported: the plies are round");
    if (inflate > 0 && !compression) return no("inflation needs compression (inflated round tubes would overlap)");
    if (!(maxCompression >= kMinMaxCompression && maxCompression <= kMaxMaxCompression))
        return no("strand.max_flattening must be 0.05 to 0.6");
    return true;
}

} // namespace knit
