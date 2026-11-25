// GPU backend: explicit plies (plies.cpp buildPlies) on the device (gpu/surfaces.h buildPlies).
//   samples   C_j(s) = yarnCoordinates(s, theta_j(s), R / r) per (yarn, sample, ply) in parallel;
//             theta_j is reduced modulo 2 pi in double before the float evaluation (s / pitch
//             reaches 2e4 turns on a body yarn: float would lose the angle);
//   splines   the interpolating B-spline per ply (a tridiagonal Thomas recurrence, sequential:
//             on the host, reported as remaining CPU work), then the arc-length tables on the
//             device (computeArcLengthTables);
//   checks    the minimum centre distance at equal s and the containment (distance of the fitted
//             ply centre to the yarn centreline, Newton from the model ring) per sample, device.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>

#include "gpu/geometry/surface.cuh"
#include "gpu/geometry/device_geometry.h"
#include "gpu/device.h"
#include "gpu/surfaces.h"

namespace knit::gpu {

namespace {

struct PlyYarn {
    int yarn;        // parent yarn (model index)
    int sampleOff;   // first sample (of m + 1)
    int m;
    double s0, ds;
};

struct PlySampling {
    CurveSetView cs;
    SurfaceView sm;
    const PlyYarn* yarns;
    int nYarns, nSamples, count, hand;
    double phase, pitch;
    float lam;
};
// The ply yarn that owns sample g: the last whose first sample is <= g (an empty item shares its successor's first
// index, and the last of equal entries owns the slots).
__device__ __forceinline__ int plyYarnOfSample(const PlySampling& in, int g) {
    int lo = 0, hi = in.nYarns - 1;
    while (lo < hi) {
        const int mid = (lo + hi + 1) >> 1;
        if (in.yarns[mid].sampleOff <= g) lo = mid; else hi = mid - 1;
    }
    return lo;
}
__global__ void plySampleKernel(PlySampling in, float* P, float* minDist) {
    const PlyYarn* py = in.yarns;
    const int nSamples = in.nSamples, count = in.count;
    const int g = threadItem();
    if (g >= nSamples) return;
    const int y = plyYarnOfSample(in, g);
    const PlyYarn Y = py[y];
    const int i = g - Y.sampleOff;
    const double s = Y.s0 + i * Y.ds;
    const YarnSurfaceView ys = devYarnSurf(in.cs, in.sm, Y.yarn);
    f3 pts[8];
    float md = 3.0e38f;
    for (int j = 0; j < count; ++j) {
        double th = in.phase + in.hand * 2 * kPiD * s / in.pitch + 2 * kPiD * j / count;
        th = fmod(th, 2 * kPiD);
        const f3 X = yarnCoordinates(ys, s, static_cast<float>(th), in.lam).X;
        if (j < 8) pts[j] = X;
        P[3 * (static_cast<long long>(j) * nSamples + g)] = X.x;
        P[3 * (static_cast<long long>(j) * nSamples + g) + 1] = X.y;
        P[3 * (static_cast<long long>(j) * nSamples + g) + 2] = X.z;
    }
    for (int j = 0; j < count && j < 8; ++j)
        for (int k = j + 1; k < count && k < 8; ++k) md = fminf_(md, norm(pts[j] - pts[k]));
    minDist[g] = md;
}

// Containment per (ply, sample u = i / 2): |distance(ply centre, yarn centreline) - R|.
struct PlyCurve {
    int plyCurve;    // index into the ply curve set
    int parent;      // parent yarn (model index)
    int sampleOff, n;  // samples (2 m + 1)
    double s0, ds;
};
struct ContainmentInput {
    CurveSetView plies, parents;
    SurfaceView sm;
    const PlyCurve* curves;
    int nCurves, nSamples;
    float R;
};
// The ply curve that owns sample g: the last whose first sample is <= g (an empty item shares its successor's first
// index, and the last of equal entries owns the slots).
__device__ __forceinline__ int plyCurveOfSample(const ContainmentInput& in, int g) {
    int lo = 0, hi = in.nCurves - 1;
    while (lo < hi) {
        const int mid = (lo + hi + 1) >> 1;
        if (in.curves[mid].sampleOff <= g) lo = mid; else hi = mid - 1;
    }
    return lo;
}
__global__ void containmentKernel(ContainmentInput in, float* dev) {
    const PlyCurve* pc = in.curves;
    const SurfaceView& sm = in.sm;
    const int g = threadItem();
    if (g >= in.nSamples) return;
    const int p = plyCurveOfSample(in, g);
    const PlyCurve C = pc[p];
    const int i = g - C.sampleOff;
    const double u = 0.5 * i, s = C.s0 + u * C.ds;
    const CurveView ply = devCurve(in.plies, C.plyCurve), yarn = devCurve(in.parents, C.parent);
    const f3 q = curveEval(ply, u);
    const YarnSurfaceRecord& ys = sm.yarns[C.parent];
    const int N = ys.nRings - 1;
    int k = static_cast<int>(llround(s / ys.ds));
    k = k < 0 ? 0 : (k > N ? N : k);
    double uu = sm.rings[ys.ringOff + k].u;
    for (int it = 0; it < 8; ++it) {
        const f3 dv = curveEval(yarn, uu) - q, d1 = curveD1(yarn, uu), d2 = curveD2(yarn, uu);
        const float gq = dot(dv, d1), hq = dot(d1, d1) + dot(dv, d2);
        const float step = hq > 0 ? -gq / hq : -gq / fmaxf_(1e-12f, dot(d1, d1));
        uu = clampd(uu + step, 0.0, static_cast<double>(yarn.nSpans));
        if (fabsf(step) < 1e-7f) break;  // CPU 1e-13 (float: 1e-7 of a span)
    }
    dev[g] = fabsf(norm(q - curveEval(yarn, uu)) - in.R);
}

void identityPlies(const CurveSet& cs, PlySet& out) {
    out.curves.yarns = cs.yarns;
    for (size_t y = 0; y < cs.yarns.size(); ++y) {
        out.yarnOf.push_back(static_cast<int>(y));
        out.plyOf.push_back(0);
        out.parentScale.push_back(0);
        out.parentOffset.push_back(0);
    }
    out.minCentreDistance = 0;
}

// The sample points of every ply; out.minCentreDistance from them.
// The parent yarns on the device.
struct PlyParents {
    const DeviceCurves& curves;
    const DeviceSurface& surface;
};
// The plies' samples: per parent yarn its first sample (PlyYarn::sampleOff), and the total.
struct PlySamples {
    const std::vector<PlyYarn>& yarns;
    int total;
};

std::vector<float> samplePlies(const PlyParents& parents, const PlySamples& sampled, const PlyGeometry& g, double lam, PlySet& out) {
    const std::vector<PlyYarn>& yarns = sampled.yarns;
    const int samples = sampled.total;
    DeviceBuffer dpy, P, minD;
    dpy.upload(yarns.data(), yarns.size() * sizeof(yarns[0]));
    P.alloc(12ull * g.count * samples + 1);
    minD.alloc(4ull * samples + 1);
    plySampleKernel<<<blocksFor(samples), kThreadsPerBlock>>>({parents.curves.view, parents.surface.view, static_cast<PlyYarn*>(dpy.data()), static_cast<int>(yarns.size()), samples, g.count, g.handedness, g.phase, g.pitch,
                                                               static_cast<float>(lam)},
                                                              P.floats(), minD.floats());
    std::vector<float> points, minima;
    checkLaunch("GPU ply samples");
    P.download(points, 3ull * g.count * samples);
    minD.download(minima, samples);
    for (float v : minima) out.minCentreDistance = std::min(out.minCentreDistance, double(v));
    return points;
}

// The interpolating splines (host: tridiagonal recurrence per ply).
CurveSet fitPlySplines(const CurveSet& cs, const PlySamples& sampled, const std::vector<float>& points, int count, PlySet& out) {
    CurveSet plies;
    for (const PlyYarn& Y : sampled.yarns) {
        const YarnPath& parent = cs.yarns[Y.yarn];
        const int firstId = static_cast<int>(plies.yarns.size());
        for (int j = 0; j < count; ++j) {
            std::vector<Vec3> pts(Y.m + 1);
            for (int i = 0; i <= Y.m; ++i) {
                const size_t at = 3 * (static_cast<size_t>(j) * sampled.total + Y.sampleOff + i);
                pts[i] = {points[at], points[at + 1], points[at + 2]};
            }
            YarnPath ply;
            ply.id = firstId + j;
            ply.name = parent.name + "_ply" + std::to_string(j);
            ply.course = parent.course;
            ply.curve = interpolateBSpline(pts);
            plies.yarns.push_back(std::move(ply));
            out.yarnOf.push_back(Y.yarn);
            out.plyOf.push_back(j);
            out.parentScale.push_back(Y.ds);
            out.parentOffset.push_back(Y.s0);
        }
    }
    return plies;
}

double worstContainment(const CurveSet& plies, const PlyParents& parents, const std::vector<PlyYarn>& yarns, int count, double R) {
    DeviceCurves dp;
    dp.upload(flattenCurves(plies, ArcBounds::Omitted));
    std::vector<PlyCurve> pc;
    int samples = 0;
    for (size_t p = 0; p < plies.yarns.size(); ++p) {
        const PlyYarn& Y = yarns[p / count];
        pc.push_back({static_cast<int>(p), Y.yarn, samples, 2 * Y.m + 1, Y.s0, Y.ds});
        samples += 2 * Y.m + 1;
    }
    DeviceBuffer dpc, dev;
    dpc.upload(pc.data(), pc.size() * sizeof(pc[0]));
    dev.alloc(4ull * samples + 1);
    containmentKernel<<<blocksFor(samples), kThreadsPerBlock>>>({dp.view, parents.curves.view, parents.surface.view, static_cast<PlyCurve*>(dpc.data()), static_cast<int>(pc.size()), samples, static_cast<float>(R)},
                                                                dev.floats());
    std::vector<float> deviations;
    checkLaunch("GPU ply containment");
    dev.download(deviations, samples);
    double worst = 0;
    for (float v : deviations) worst = std::max(worst, double(v));
    return worst;
}

} // namespace

bool buildPlies(const PlyRequest& request, PlySet& out, PliesGpuStats* st, std::string* err) try {
    const CurveSet& cs = request.curves;
    const SurfaceModel& sm = request.surface;
    const PlyGeometry& g = request.geometry;
    const double r = request.radius;
    if (!available(err)) return false;
    Reservation lease(Reservation::exclusiveBytes());  // runs alone on the device (runtime.h)
    const auto t0 = Clock::now();
    out = PlySet{};
    out.layout = plyLayout(g, r);
    out.curves.kind = CurveKind::Surface;
    out.curves.waleSpacing = cs.waleSpacing;
    out.curves.courseSpacing = cs.courseSpacing;
    out.minCentreDistance = 1e300;
    if (g.count <= 1) {
        identityPlies(cs, out);
        return true;
    }
    if (g.count > 8) throw DeviceError("GPU plies: at most 8 plies");
    DeviceCurves curves;
    DeviceSurface surface;
    curves.upload(flattenCurves(cs, ArcBounds::Omitted));
    surface.upload(flattenSurface(sm));
    std::vector<PlyYarn> yarns;
    int samples = 0;
    for (const YarnSurface& ys : sm.yarns) {
        const double s0 = std::min(0.5 * r, 0.25 * ys.length), L = ys.length - 2 * s0;
        const int m = std::max(1, static_cast<int>(std::ceil(L / std::min(g.pitch / g.samplesPerTwist, g.maxSampleStep))));
        yarns.push_back({ys.yarnId, samples, m, s0, L / m});
        samples += m + 1;
    }
    const PlyParents parents{curves, surface};
    const PlySamples sampled{yarns, samples};
    const std::vector<float> points = samplePlies(parents, sampled, g, out.layout.R / r, out);
    const double msSamples = elapsedMs(t0);
    const auto t1 = Clock::now();
    CurveSet plies = fitPlySplines(cs, sampled, points, g.count, out);
    const double msSplines = elapsedMs(t1);
    const auto t2 = Clock::now();
    computeArcLengthTables(plies);
    const double msFinalize = elapsedMs(t2);
    const auto t3 = Clock::now();
    out.containment = std::max(out.containment, worstContainment(plies, parents, yarns, g.count, out.layout.R));
    out.curves.yarns = std::move(plies.yarns);
    out.ms = elapsedMs(t0);
    if (st) {
        st->msSamples += msSamples;
        st->msSplinesHost += msSplines;
        st->msFinalize += msFinalize;
        st->msContainment += elapsedMs(t3);
        st->plies += static_cast<int>(out.curves.yarns.size());
    }
    return true;
} catch (const DeviceError& e) {
    return reportFailure(e, err);
}

std::string PliesGpuStats::line() const {
    char b[220];
    std::snprintf(b, sizeof b, "GPU plies: %d plies; samples %.0f ms, splines (host Thomas) %.0f ms, tables %.0f ms, containment %.0f ms", plies,
                  msSamples, msSplinesHost, msFinalize, msContainment);
    return b;
}

} // namespace knit::gpu
