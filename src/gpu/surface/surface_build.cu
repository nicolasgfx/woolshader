// GPU backend: the surface model (buildSurface) and gate (d) (checkSurface) on the device
// (gpu/surfaces.h buildSurface / checkSurface), single precision (math.cuh).
//
// Rings: per ring in parallel (arc-length parameter by the device inversion, centre, tangent).
// Frames: the rotation-minimizing double reflection (Wang et al. 2008), a sequential recurrence
// along each yarn, is evaluated IN ORDER within chunks of kChunk rings; the frame at each chunk
// start comes from the composed linear maps of the double reflections before it (H2 H1 per step,
// independent of the frame), followed by the CPU's projection and normalization. In exact
// arithmetic this is the CPU recurrence; only the rounding differs.
// Compression (SurfaceParams::compression): the per-ring requirements from the capsule axes
// within 2 r (outside the own half-turn window), refined by Newton on the neighbour's spline,
// plus the fold requirement, the smoothing max-filter, the axes (compression_fit.h, the twin of
// the CPU's). The capsule index (SurfaceIndex) is built on the device (device_geometry.h).
#include <chrono>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <cstring>

#include "gpu/geometry/device_geometry.h"
#include "gpu/geometry/surface.cuh"
#include "knitting/tube/compression_fit.h"
#include "gpu/device.h"
#include "gpu/surfaces.h"

namespace knit::gpu {

namespace {

constexpr int kChunk = 256;

struct Mat3 {
    float m[9];
};
__device__ __forceinline__ f3 mul(const Mat3& A, f3 v) {
    return {A.m[0] * v.x + A.m[1] * v.y + A.m[2] * v.z, A.m[3] * v.x + A.m[4] * v.y + A.m[5] * v.z, A.m[6] * v.x + A.m[7] * v.y + A.m[8] * v.z};
}
__device__ __forceinline__ Mat3 mul(const Mat3& A, const Mat3& B) {
    Mat3 C;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) C.m[3 * i + j] = A.m[3 * i] * B.m[j] + A.m[3 * i + 1] * B.m[3 + j] + A.m[3 * i + 2] * B.m[6 + j];
    return C;
}
__device__ __forceinline__ Mat3 eye() { return Mat3{{1, 0, 0, 0, 1, 0, 0, 0, 1}}; }
// I - (2 / c) v v^T (identity for c <= 0).
__device__ __forceinline__ Mat3 householder(f3 v, float c) {
    if (!(c > 0)) return eye();
    const float k = 2 / c;
    return Mat3{{1 - k * v.x * v.x, -k * v.x * v.y, -k * v.x * v.z, -k * v.y * v.x, 1 - k * v.y * v.y, -k * v.y * v.z, -k * v.z * v.x,
               -k * v.z * v.y, 1 - k * v.z * v.z}};
}

// --- the model: rings and their rotation-minimizing frames
struct RingYarn {
    int ringOff, nRings;  // nRings = N + 1
    double ds, length;
    int chunkOff, nChunks;
};

struct RingYarns {
    const RingYarn* items;
    int count;
};

__device__ __forceinline__ int yarnOfRing(const RingYarns& ys, long long g) {
    int lo = 0, hi = ys.count - 1;
    while (lo < hi) {
        const int mid = (lo + hi + 1) >> 1;
        if (ys.items[mid].ringOff <= g) lo = mid; else hi = mid - 1;
    }
    return lo;
}

__device__ __forceinline__ int yarnOfChunk(const RingYarns& ys, int chunk) {
    int lo = 0, hi = ys.count - 1;
    while (lo < hi) {
        const int mid = (lo + hi + 1) >> 1;
        if (ys.items[mid].chunkOff <= chunk) lo = mid; else hi = mid - 1;
    }
    return lo;
}

// Ring centres and tangents (yarn_surface.cpp buildSurface, first loop).
__global__ void ringKernel(CurveSetView cs, RingYarns ys, long long n, float r, RingRecord* rings) {
    const long long g = threadItem64();
    if (g >= n) return;
    const int y = yarnOfRing(ys, g);
    const RingYarn Y = ys.items[y];
    const int k = static_cast<int>(g - Y.ringOff);
    const CurveView c = devCurve(cs, y);
    RingRecord R{};
    R.s = k * Y.ds;
    R.u = paramAtArcLength(c, R.s);
    const f3 p = curveEval(c, R.u), T = normalized(curveD1(c, R.u));
    R.c[0] = p.x; R.c[1] = p.y; R.c[2] = p.z;
    R.T[0] = T.x; R.T[1] = T.y; R.T[2] = T.z;
    R.a = R.b = r * thicknessAt(cs, y, R.u);
    R.phi = 0;
    rings[g] = R;
}

// The step map of ring k (k >= 1 within its yarn): n_k = proj(H2 H1 n_{k-1}).
__device__ __forceinline__ Mat3 stepMap(const RingRecord& P, const RingRecord& R) {
    const f3 v1 = fromF(R.c) - fromF(P.c);
    const float c1 = dot(v1, v1);
    const Mat3 H1 = householder(v1, c1);
    const f3 tL = mul(H1, fromF(P.T));
    const f3 v2 = fromF(R.T) - tL;
    return mul(householder(v2, dot(v2, v2)), H1);
}

// One step exactly as the CPU (double reflection, then project and normalize).
__device__ __forceinline__ f3 rmfStep(const RingRecord& P, f3 n1P, const RingRecord& R) {
    const f3 v1 = fromF(R.c) - fromF(P.c);
    const float c1 = dot(v1, v1);
    f3 rL = n1P, tL = fromF(P.T);
    if (c1 > 0) {
        rL = n1P - v1 * (2 / c1 * dot(v1, n1P));
        tL = fromF(P.T) - v1 * (2 / c1 * dot(v1, fromF(P.T)));
    }
    const f3 v2 = fromF(R.T) - tL;
    const float c2 = dot(v2, v2);
    f3 n = c2 > 0 ? rL - v2 * (2 / c2 * dot(v2, rL)) : rL;
    const f3 T = fromF(R.T);
    return normalized(n - T * dot(n, T));
}

// Per chunk: the composed map of its steps (rings start+1 .. end).
__global__ void chunkMapKernel(RingYarns ys, int nChunks, const RingRecord* rings, Mat3* maps) {
    const int j = threadItem();
    if (j >= nChunks) return;
    const RingYarn Y = ys.items[yarnOfChunk(ys, j)];
    const int c = j - Y.chunkOff;
    const int k0 = c * kChunk, k1 = min(Y.nRings - 1, k0 + kChunk);
    Mat3 A = eye();
    for (int k = k0 + 1; k <= k1; ++k) A = mul(stepMap(rings[Y.ringOff + k - 1], rings[Y.ringOff + k]), A);
    maps[j] = A;
}

// Per yarn, sequentially over its chunks: the frame at every chunk start.
__global__ void chunkStartKernel(RingYarns ys, const RingRecord* rings, const Mat3* maps, f3* starts) {
    const int y = threadItem();
    if (y >= ys.count) return;
    const RingYarn Y = ys.items[y];
    const RingRecord& R0 = rings[Y.ringOff];
    const f3 T0 = fromF(R0.T);
    f3 n = f3{0, 0, 1} - T0 * T0.z;
    if (norm(n) < 1e-6f) n = f3{1, 0, 0} - T0 * T0.x;
    n = normalized(n);
    for (int c = 0; c < Y.nChunks; ++c) {
        starts[Y.chunkOff + c] = n;
        const int k1 = min(Y.nRings - 1, c * kChunk + kChunk);
        const f3 T = fromF(rings[Y.ringOff + k1].T);
        f3 m = mul(maps[Y.chunkOff + c], n);
        n = normalized(m - T * dot(m, T));
    }
}

// Per chunk: the CPU recurrence from its start frame.
__global__ void chunkPropagateKernel(RingYarns ys, int nChunks, RingRecord* rings, const f3* starts) {
    const int j = threadItem();
    if (j >= nChunks) return;
    const RingYarn Y = ys.items[yarnOfChunk(ys, j)];
    const int c = j - Y.chunkOff;
    const int k0 = c * kChunk, k1 = min(Y.nRings - 1, k0 + kChunk);
    f3 n = starts[j];
    if (c == 0) {
        RingRecord& R = rings[Y.ringOff];
        R.n1[0] = n.x; R.n1[1] = n.y; R.n1[2] = n.z;
        const f3 n2 = cross(fromF(R.T), n);
        R.n2[0] = n2.x; R.n2[1] = n2.y; R.n2[2] = n2.z;
    }
    for (int k = k0 + 1; k <= k1; ++k) {
        RingRecord& R = rings[Y.ringOff + k];
        n = rmfStep(rings[Y.ringOff + k - 1], n, R);
        R.n1[0] = n.x; R.n1[1] = n.y; R.n1[2] = n.z;
        const f3 n2 = cross(fromF(R.T), n);
        R.n2[0] = n2.x; R.n2[1] = n2.y; R.n2[2] = n2.z;
    }
}

// --- the capsule index (built by the compression, and again by the checks)
// SurfaceIndex capsules: every `step` rings, radius max(a, b) + chord
// deviation + 1e-3 mm; their boxes (padded by 1e-5 relative for float rounding).
struct CapsuleYarn {
    int ringOff, nRings, step, capOff, nCaps, yarnId;
};

struct CapsuleInput {
    const CapsuleYarn* ys;
    int nY, nCaps;
    const RingRecord* rings;
};

// The yarn that owns capsule g: the last whose first capsule is <= g (an empty item shares its successor's first
// index, and the last of equal entries owns the slots).
__device__ __forceinline__ int yarnOfCapsule(const CapsuleInput& in, int g) {
    int a = 0, b = in.nY - 1;
    while (a < b) {
        const int mid = (a + b + 1) >> 1;
        if (in.ys[mid].capOff <= g) a = mid; else b = mid - 1;
    }
    return a;
}

__global__ void capsuleKernel(CapsuleInput in, CapsuleRecord* caps, float* lo, float* hi) {
    const CapsuleYarn* ys = in.ys;
    const int g = threadItem();
    if (g >= in.nCaps) return;
    const int y = yarnOfCapsule(in, g);
    const CapsuleYarn Y = ys[y];
    const int N = Y.nRings - 1;
    const int k = (g - Y.capOff) * Y.step, k1 = min(N, k + Y.step);
    const RingRecord* R = in.rings + Y.ringOff;
    const f3 A = fromF(R[k].c), B = fromF(R[k1].c);
    float rad = 0, dev = 0;
    for (int j = k; j <= k1; ++j) {
        rad = fmaxf_(rad, fmaxf_(R[j].a, R[j].b));
        float t;
        dev = fmaxf_(dev, pointSegment(fromF(R[j].c), A, B, t));
    }
    rad += dev + 1e-3f;
    CapsuleRecord c{};
    c.a[0] = A.x; c.a[1] = A.y; c.a[2] = A.z;
    c.b[0] = B.x; c.b[1] = B.y; c.b[2] = B.z;
    c.R = rad;
    c.yarn = Y.yarnId;
    c.s0 = R[k].s;
    c.s1 = R[k1].s;
    caps[g] = c;
    const f3 l = vmin(A, B) - f3{rad, rad, rad}, h = vmax(A, B) + f3{rad, rad, rad};
    const float e = 1e-5f * (1 + fmaxf_(fmaxf_(fabsf(l.x), fabsf(l.y)), fabsf(l.z)));
    lo[3 * g] = l.x - e; lo[3 * g + 1] = l.y - e; lo[3 * g + 2] = l.z - e;
    hi[3 * g] = h.x + e; hi[3 * g + 1] = h.y + e; hi[3 * g + 2] = h.z + e;
}

// --- the compression
// Compression requirements per ring (yarn_surface.cpp buildSurface, compression_fit.h): half distances h
// to nearby centrelines and their directions psi in the ring frame: up to kMaxReq directions and
// the fold requirement.
using CompressionRequirement = knit::cfit::RequirementSetF;

struct CompressionInput {
    CurveSetView cs;
    SurfaceIndexQuery I;
    RingYarns ys;
    long long n;
    float W;
};

__global__ void compressReqKernel(CompressionInput in, CompressionRequirement* req) {
    const SurfaceIndexQuery& I = in.I;
    const CurveSetView& cs = in.cs;
    const float W = in.W;
    const long long g = threadItem64();
    if (g >= in.n) return;
    const int y = yarnOfRing(in.ys, g);
    const DeviceRing R = loadRing(I.sm.rings[g]);
    CompressionRequirement best;
    const Exclude own{y, R.s - W, R.s + W};
    const float r = R.a, rad = 2 * r;  // round here: this ring's own radius
    // The own window in the spline parameter: Newton started on a distant portion of the own
    // yarn can converge back onto the query point (in float the converged residual is ~1e-5 mm),
    // so the window is tested on the converged parameter (the CPU does the same).
    const CurveView C1 = devCurve(cs, y);
    const double uOwnLo = paramAtArcLength(C1, R.s - W), uOwnHi = paramAtArcLength(C1, R.s + W);
    // The fold requirement: support along the principal normal below rho.
    {
        const f3 d1 = curveD1(C1, R.u), d2 = curveD2(C1, R.u);
        const float spd = norm(d1);
        if (spd > 0) {
            const f3 T = d1 / spd, kv = (d2 - T * dot(d2, T)) / (spd * spd);
            const float kappa = norm(kv);
            if (kappa > 1e-12f) {
                const float h = 1.0f / (kappa * (1 + static_cast<float>(knit::cfit::kFoldMargin)));
                if (h < r) knit::cfit::insert(best, h, atan2f(dot(kv, R.n2), dot(kv, R.n1)));
            }
        }
    }
    for (GridCursor cursor(I.idx.grid, R.c - f3{rad, rad, rad}, R.c + f3{rad, rad, rad}); cursor.next();) {
        const int i = cursor.item();
        const CapsuleRecord& c = I.idx.caps[i];
        float t;
        const float dAxis = pointSegment(R.c, fromF(c.a), fromF(c.b), t);
        const double sc = c.s0 + t * (c.s1 - c.s0);
        if (!(dAxis < rad) || own.covers(c.yarn, sc)) continue;
        if (dAxis >= rad + 0.01f) continue;
        const YarnSurfaceRecord& ys2 = I.sm.yarns[c.yarn];
        const CurveView C2 = devCurve(cs, c.yarn);
        int k2 = static_cast<int>(llround(sc / ys2.ds));
        k2 = k2 < 0 ? 0 : (k2 > ys2.nRings - 1 ? ys2.nRings - 1 : k2);
        double u = I.sm.rings[ys2.ringOff + k2].u;
        for (int it = 0; it < 8; ++it) {
            const f3 dv = curveEval(C2, u) - R.c, d1 = curveD1(C2, u), d2 = curveD2(C2, u);
            const float gq = dot(dv, d1), hq = dot(d1, d1) + dot(dv, d2);
            const float step = hq > 0 ? -gq / hq : -gq / fmaxf_(1e-12f, dot(d1, d1));
            u = clampd(u + step, 0.0, static_cast<double>(C2.nSpans));
            if (fabsf(step) < 1e-7f) break;  // CPU 1e-12 (float: 1e-7 of a span)
        }
        const f3 q = curveEval(C2, u);
        const float d = norm(q - R.c);
        if (d >= rad || d <= kTinyLen || (c.yarn == y && u >= uOwnLo && u <= uOwnHi)) continue;
        const float h = 0.5f * d;
        f3 dir = q - R.c;
        dir = dir - R.T * dot(dir, R.T);
        if (norm(dir) < kTinyLen) continue;
        knit::cfit::insert(best, h, atan2f(dot(dir, R.n2), dot(dir, R.n1)));
    }
    req[g] = best;
}

// The smoothing max-filter of the flattening deficit along s, then the axes (compression_fit.h).
struct CompressionFitInput {
    RingYarns ys;
    long long n;
    const CompressionRequirement* req;
    float sigma, lowerRatio;  // lowerRatio: the smallest semi-axis over the ring's round radius
};

__global__ void compressApplyKernel(CompressionFitInput in, RingRecord* rings, unsigned long long* clamped) {
    const CompressionRequirement* req = in.req;
    const float sigma = in.sigma;
    const long long g = threadItem64();
    if (g >= in.n) return;
    const float r = rings[g].a;  // this ring's own, before the fit
    const RingYarn Y = in.ys.items[yarnOfRing(in.ys, g)];
    const int k = static_cast<int>(g - Y.ringOff);
    CompressionRequirement out = req[g];
    if (sigma > 0) {
        const int half = static_cast<int>(ceil(3 * sigma / Y.ds));
        const int j0 = max(0, k - half), j1 = min(Y.nRings - 1, k + half);
        for (int j = j0; j <= j1; ++j) {
            const CompressionRequirement& q = req[Y.ringOff + j];
            if (q.n == 0) continue;
            const float x = static_cast<float>((j - k) * Y.ds / sigma);
            const float w = expf(-0.5f * x * x);
            for (int i = 0; i < q.n; ++i) {
                const float h = knit::cfit::smoothed(q.e[i].h, w, r);
                if (h >= r) continue;
                knit::cfit::insert(out, h, q.e[i].psi);
            }
        }
    }
    const knit::cfit::SectionF f = knit::cfit::fit(out, r, in.lowerRatio * r);
    if (!f.flattened) return;
    RingRecord& R = rings[g];
    if (f.clamped) atomicAdd(clamped, 1ull);
    R.b = f.b;
    R.phi = f.phi;
    R.a = f.a;
}

// --- the checks: folds, collisions
// checkSurface, folds: per ring rho / support ratio (3e38: no curvature).
__global__ void foldKernel(CurveSetView cs, RingYarns ys, long long n, const RingRecord* rings, float* ratio) {
    const long long g = threadItem64();
    if (g >= n) return;
    const int y = yarnOfRing(ys, g);
    const CurveView c = devCurve(cs, y);
    const DeviceRing R = loadRing(rings[g]);
    const f3 d1 = curveD1(c, R.u), d2 = curveD2(c, R.u);
    const float sp = norm(d1);
    ratio[g] = 3.0e38f;
    if (sp <= 0) return;
    const f3 T = d1 / sp;
    const f3 kv = (d2 - T * dot(d2, T)) / (sp * sp);
    const float kappa = norm(kv);
    if (kappa < 1e-12f) return;
    f3 ea, eb;
    ringAxes(R, ea, eb);
    const f3 dir = kv / kappa;
    const float x = R.a * dot(dir, ea), yy = R.b * dot(dir, eb);
    const float h = sqrtf(x * x + yy * yy);
    ratio[g] = 1.0f / (kappa * h);
}

// Reductions in the CPU's order (first minimum / maximum in sample order): packed keys
// (value bits << 32 | index) with atomicMin / atomicMax; values are non-negative floats, so
// their bits order as the values.
__global__ void foldReduceKernel(const float* ratio, long long n, unsigned long long* minKey, unsigned long long* folds) {
    const long long i = threadItem64();
    if (i >= n) return;
    const float q = ratio[i];
    if (q >= 3.0e38f) return;
    if (q < 1.0f) atomicAdd(folds, 1ull);
    atomicMin(minKey, (static_cast<unsigned long long>(__float_as_uint(q)) << 32) | static_cast<unsigned long long>(i));
}

// checkSurface, collisions: sample i = (sampled ring, angle). out: depth (-1 if not penetrating
// beyond the tolerance) and the hit (yarn, s).
struct SampleYarn {
    int ringOff, nRings, sampleOff, nSampledRings;
    float W;
};

struct CollisionSampling {
    const SampleYarn* ys;
    int nY;
    long long n;
    int stride, angles;
    float tol;
};

struct CollisionOutput {
    float* depth;
    int* hitYarn;
    double* hitS;
    float* points;
};

// The yarn that owns sampled ring sr: the last whose first sampled ring is <= sr (an empty item shares its successor's first
// index, and the last of equal entries owns the slots).
__device__ __forceinline__ int yarnOfCollisionSample(const CollisionSampling& s, long long sr) {
    int lo = 0, hi = s.nY - 1;
    while (lo < hi) {
        const int mid = (lo + hi + 1) >> 1;
        if (s.ys[mid].sampleOff <= sr) lo = mid; else hi = mid - 1;
    }
    return lo;
}

__global__ void collisionKernel(SurfaceIndexQuery I, CollisionSampling s, CollisionOutput out) {
    const SampleYarn* ys = s.ys;
    const int angles = s.angles;
    const long long i = threadItem64();
    if (i >= s.n) return;
    const long long sr = i / angles;
    const int a = static_cast<int>(i % angles);
    const int y = yarnOfCollisionSample(s, sr);
    const SampleYarn Y = ys[y];
    const int k = static_cast<int>(sr - Y.sampleOff) * s.stride;
    const DeviceRing R = loadRing(I.sm.rings[Y.ringOff + k]);
    f3 ea, eb;
    ringAxes(R, ea, eb);
    const float theta = 2 * kPiF * a / angles;
    float sn, cs;
    sincosf(theta - R.phi, &sn, &cs);
    const f3 q = R.c + ea * (R.a * cs) + eb * (R.b * sn);
    out.points[3 * i] = q.x; out.points[3 * i + 1] = q.y; out.points[3 * i + 2] = q.z;
    const Exclude own{y, R.s - Y.W, R.s + Y.W};
    const Hit h = exactInside(I, q, 0.0f, Exclude{}, own);
    out.depth[i] = h.inside && h.depth > s.tol ? h.depth : -1.0f;
    out.hitYarn[i] = h.yarn;
    out.hitS[i] = h.s;
}

__global__ void collisionReduceKernel(const float* depth, long long n, unsigned long long* maxKey, unsigned long long* count) {
    const long long i = threadItem64();
    if (i >= n) return;
    const float d = depth[i];
    if (d < 0) return;
    atomicAdd(count, 1ull);
    atomicMax(maxKey, (static_cast<unsigned long long>(__float_as_uint(d)) << 32) | (0xFFFFFFFFull - static_cast<unsigned long long>(i)));
}

// --- the tube
// buildTubeMesh: one thread per tube ring copy (its sides + 1 vertices), ringAt on the model.
struct TubeRing {
    int yarn;          // index into the model yarns
    double s;          // arc length
    double uOff;       // uv u offset of this copy
    unsigned first;    // first vertex
    double pitchFactor, phase;  // the yarn's ply-map variation (yarnUvOf)
};

struct TubeInput {
    CurveSetView cs;
    SurfaceView sm;
    const TubeRing* rings;
    int n, sides;
    float uvRepeat, r;
};

struct TubeOutput {
    float* pos;
    float* nrm;
    float* uv;
    double* s;
    double* t;
    double* theta;
    int* yarn;
    float* axis;
};

__global__ void tubeKernel(TubeInput in, TubeOutput out) {
    const int sides = in.sides;
    const int j = threadItem();
    if (j >= in.n) return;
    const TubeRing T = in.rings[j];
    const YarnSurfaceView ys = devYarnSurf(in.cs, in.sm, T.yarn);
    const DeviceRing R = ringAt(ys, T.s);
    f3 ea, eb;
    ringAxes(R, ea, eb);
    for (int i = 0; i <= sides; ++i) {
        const unsigned v = T.first + i;
        const float theta = 2 * kPiF * i / sides;
        float sn, c;
        sincosf(theta - R.phi, &sn, &c);
        const f3 nn = normalized(ea * (c / R.a) + eb * (sn / R.b));
        const f3 p = R.c + ea * (R.a * c) + eb * (R.b * sn);
        out.pos[3 * v] = p.x; out.pos[3 * v + 1] = p.y; out.pos[3 * v + 2] = p.z;
        out.nrm[3 * v] = nn.x; out.nrm[3 * v + 1] = nn.y; out.nrm[3 * v + 2] = nn.z;
        out.uv[2 * v] = static_cast<float>(T.s / in.uvRepeat / T.pitchFactor + T.phase - T.uOff);
        out.uv[2 * v + 1] = static_cast<float>(double(i) / sides);
        out.s[v] = T.s;
        out.t[v] = ys.length > 0 ? T.s / ys.length : 0;
        out.theta[v] = 2 * kPiD * i / sides;
        out.yarn[v] = T.yarn;
        out.axis[v] = R.b / in.r;
    }
}

// --- the host side, in the same order

// A surface model on the device, built by the constructor's stages in order.
class DeviceSurfaceModel {
public:
    explicit DeviceSurfaceModel(const SurfaceRequest& request);
    SurfaceModel download(const PatternParams& pp) const;

    DeviceCurves curves;
    DeviceSurface surf;
    std::vector<YarnSurfaceRecord> yarns;  // host copy
    long long nRings = 0;
    float r = 0;
    unsigned long long clamped = 0;  // compression: rings clamped at maxCompression

    RingYarns ringYarns() const { return {static_cast<RingYarn*>(ringYarnsDev_.data()), static_cast<int>(ringYarns_.size())}; }
    RingRecord* rings() const { return static_cast<RingRecord*>(surf.rings.data()); }

private:
    void planRings(const CurveSet& cs, const SurfaceParams& sp);
    void computeRings();
    void propagateFrames();
    void compressAxes(const SurfaceParams& sp);

    std::vector<RingYarn> ringYarns_;
    DeviceBuffer ringYarnsDev_;
    int chunks_ = 0;
};

// The capsule index of a surface model (SurfaceIndex).
struct DeviceCapsuleIndex {
    DeviceBuffer caps;
    DeviceGrid grid;

    explicit DeviceCapsuleIndex(const DeviceSurfaceModel& d) {
        std::vector<CapsuleYarn> cy;
        int nCaps = 0;
        for (size_t y = 0; y < d.yarns.size(); ++y) {
            const int N = d.yarns[y].nRings - 1;
            const int step = std::max(1, static_cast<int>(std::lround(0.1 / d.yarns[y].ds)));
            CapsuleYarn c{d.yarns[y].ringOff, d.yarns[y].nRings, step, nCaps, (N + step - 1) / step, static_cast<int>(y)};
            nCaps += c.nCaps;
            cy.push_back(c);
        }
        DeviceBuffer cyDev;
        cyDev.upload(cy.data(), cy.size() * sizeof(cy[0]));
        caps.alloc(sizeof(CapsuleRecord) * std::max(nCaps, 1));
        grid.boxLo.alloc(12ull * std::max(nCaps, 1));
        grid.boxHi.alloc(12ull * std::max(nCaps, 1));
        if (nCaps > 0) {
            capsuleKernel<<<blocksFor(nCaps), kThreadsPerBlock>>>({static_cast<CapsuleYarn*>(cyDev.data()), static_cast<int>(cy.size()), nCaps, d.rings()}, static_cast<CapsuleRecord*>(caps.data()),
                                                                  grid.boxLo.floats(), grid.boxHi.floats());
            checkLaunch("GPU surface capsules");
        }
        // The CPU cell: 1 + max capsule radius (the model's extent + small deviation + 1e-3).
        float ext = d.r;
        for (const YarnSurfaceRecord& f : d.yarns) ext = std::max(ext, f.ext);
        grid.build(nCaps, 1.0f + ext + 0.01f);
    }

    SurfaceIndexQuery query(const DeviceSurfaceModel& d) const {
        SurfaceIndexQuery q;
        q.idx.caps = static_cast<CapsuleRecord*>(caps.data());
        q.idx.grid = grid.view;
        q.sm = d.surf.view;
        return q;
    }
};

DeviceSurfaceModel::DeviceSurfaceModel(const SurfaceRequest& request) : r(float(request.pattern.radius())) {
    curves.upload(flattenCurves(request.curves, ArcBounds::Omitted));
    planRings(request.curves, request.surface);
    if (nRings == 0) return;
    computeRings();
    propagateFrames();
    if (request.surface.compression) compressAxes(request.surface);
}

void DeviceSurfaceModel::planRings(const CurveSet& cs, const SurfaceParams& sp) {
    long long total = 0;
    for (size_t y = 0; y < cs.yarns.size(); ++y) {
        const double L = cs.yarns[y].length();
        const int N = std::max(1, static_cast<int>(std::ceil(L / sp.ringSpacing)));
        RingYarn ry{};
        ry.ringOff = static_cast<int>(total);
        ry.nRings = N + 1;
        ry.ds = L / N;
        ry.length = L;
        ry.chunkOff = chunks_;
        ry.nChunks = (N + kChunk - 1) / kChunk;
        chunks_ += ry.nChunks;
        total += N + 1;
        ringYarns_.push_back(ry);
        YarnSurfaceRecord f{};
        f.ringOff = ry.ringOff;
        f.nRings = ry.nRings;
        f.ds = ry.ds;
        f.length = L;
        f.radius = r;
        f.ext = r * static_cast<float>(cs.yarns[y].largestThickness());  // the fit never grows a ring (compression_fit.h)
        yarns.push_back(f);
    }
    if (total >= (1ll << 31)) throw DeviceError("GPU surface: too many rings");
    nRings = total;
    ringYarnsDev_.upload(ringYarns_.data(), ringYarns_.size() * sizeof(ringYarns_[0]));
    surf.rings.alloc(sizeof(RingRecord) * std::max(1ll, total));
    surf.yarns.upload(yarns.data(), yarns.size() * sizeof(yarns[0]));
    surf.view.yarns = static_cast<YarnSurfaceRecord*>(surf.yarns.data());
    surf.view.nYarns = static_cast<int>(yarns.size());
    surf.view.rings = rings();
}

void DeviceSurfaceModel::computeRings() {
    ringKernel<<<blocksFor(nRings), kThreadsPerBlock>>>(curves.view, ringYarns(), nRings, r, rings());
    checkLaunch("GPU surface rings");
}

void DeviceSurfaceModel::propagateFrames() {
    const RingYarns ys = ringYarns();
    DeviceBuffer maps, starts;
    maps.alloc(sizeof(Mat3) * std::max(chunks_, 1));
    starts.alloc(sizeof(f3) * std::max(chunks_, 1));
    chunkMapKernel<<<blocksFor(chunks_), kThreadsPerBlock>>>(ys, chunks_, rings(), static_cast<Mat3*>(maps.data()));
    chunkStartKernel<<<blocksFor(ys.count, 64), 64>>>(ys, rings(), static_cast<Mat3*>(maps.data()), static_cast<f3*>(starts.data()));
    chunkPropagateKernel<<<blocksFor(chunks_), kThreadsPerBlock>>>(ys, chunks_, rings(), static_cast<f3*>(starts.data()));
    checkLaunch("GPU surface frames");
    check(cudaDeviceSynchronize(), "GPU surface frames (sync)");
}

// The requirement on the circular model (its capsule index), then the axes.
void DeviceSurfaceModel::compressAxes(const SurfaceParams& sp) {
    const DeviceCapsuleIndex index(*this);
    DeviceBuffer requirements, clampedCount;
    requirements.alloc(sizeof(CompressionRequirement) * nRings);
    clampedCount.alloc(8);
    check(cudaMemset(clampedCount.data(), 0, 8), "GPU surface");
    const float W = static_cast<float>(knit::kPi * (1 + sp.maxExpansion) * r);
    CompressionRequirement* req = static_cast<CompressionRequirement*>(requirements.data());
    compressReqKernel<<<blocksFor(nRings), kThreadsPerBlock>>>({curves.view, index.query(*this), ringYarns(), nRings, W}, req);
    compressApplyKernel<<<blocksFor(nRings), kThreadsPerBlock>>>({ringYarns(), nRings, req, static_cast<float>(sp.smoothing), static_cast<float>(1 - sp.maxCompression)},
                                                                 rings(), clampedCount.ulongs());
    checkLaunch("GPU surface compression");
    clampedCount.download(&clamped, 8);
}

SurfaceModel DeviceSurfaceModel::download(const PatternParams& pp) const {
    std::vector<RingRecord> records(static_cast<size_t>(nRings));
    surf.rings.download(records.data(), records.size() * sizeof(records[0]));
    SurfaceModel out;
    for (size_t y = 0; y < yarns.size(); ++y) {
        YarnSurface ys;
        ys.yarnId = static_cast<int>(y);
        ys.radius = pp.radius();
        ys.length = yarns[y].length;
        ys.ds = yarns[y].ds;
        ys.rings.resize(yarns[y].nRings);
        for (int k = 0; k < yarns[y].nRings; ++k) {
            const RingRecord& f = records[static_cast<size_t>(yarns[y].ringOff) + k];
            Ring& R = ys.rings[k];
            R.s = f.s;
            R.u = f.u;
            R.c = {f.c[0], f.c[1], f.c[2]};
            R.T = {f.T[0], f.T[1], f.T[2]};
            R.n1 = {f.n1[0], f.n1[1], f.n1[2]};
            R.n2 = {f.n2[0], f.n2[1], f.n2[2]};
            R.a = f.a;
            R.b = f.b;
            R.phi = f.phi;
        }
        out.yarns.push_back(std::move(ys));
    }
    // Compression report (yarn_surface.cpp).
    out.compression.clamped = static_cast<int>(clamped);
    for (const YarnSurface& ys : out.yarns)
        for (const Ring& R : ys.rings) {
            out.compression.minAxisRatio = std::min(out.compression.minAxisRatio, R.b / pp.radius());
            out.compression.maxAreaChange = std::max(out.compression.maxAreaChange, std::fabs(R.a * R.b / (pp.radius() * pp.radius()) - 1));
        }
    return out;
}

// Gate (d): folds. Ring order = (yarn, ring): the first minimum, as the CPU.
void checkFolds(const DeviceSurfaceModel& d, SurfaceCheck& out) {
    const int Y = static_cast<int>(d.yarns.size());
    DeviceBuffer ratio, keys;
    ratio.alloc(4ull * std::max(1ll, d.nRings));
    keys.alloc(16);
    const unsigned long long start[2] = {~0ull, 0ull};
    check(cudaMemcpy(keys.data(), start, sizeof start, cudaMemcpyHostToDevice), "GPU surface");
    if (d.nRings > 0) {
        foldKernel<<<blocksFor(d.nRings), kThreadsPerBlock>>>(d.curves.view, d.ringYarns(), d.nRings, d.rings(), ratio.floats());
        foldReduceKernel<<<blocksFor(d.nRings), kThreadsPerBlock>>>(ratio.floats(), d.nRings, keys.ulongs(), keys.ulongs() + 1);
    }
    unsigned long long result[2];
    checkLaunch("GPU gate d folds");
    keys.download(result, sizeof result);
    out.foldRings = static_cast<int>(result[1]);
    if (result[0] != ~0ull) {
        const long long g = static_cast<long long>(result[0] & 0xFFFFFFFFull);
        const unsigned bits = static_cast<unsigned>(result[0] >> 32);
        float q;
        std::memcpy(&q, &bits, 4);
        int y = 0;
        while (y + 1 < Y && d.yarns[y + 1].ringOff <= g) ++y;
        out.minFoldRatio = q;
        out.worstFoldYarn = y;
        out.worstFoldS = (g - d.yarns[y].ringOff) * d.yarns[y].ds;
    }
    out.fold = out.minFoldRatio < 1.0 ? Status::Fail : Status::Pass;
}

// Gate (d): sampled collisions against the capsule index.
void checkCollisions(const DeviceSurfaceModel& d, const SurfaceIndexQuery& index, const SurfaceCheckParams& p, SurfaceCheck& out) {
    const int Y = static_cast<int>(d.yarns.size());
    std::vector<SampleYarn> sy;
    long long sampled = 0;
    const int stride = std::max(1, p.ringStride);
    for (int y = 0; y < Y; ++y) {
        SampleYarn s{};
        s.ringOff = d.yarns[y].ringOff;
        s.nRings = d.yarns[y].nRings;
        s.sampleOff = static_cast<int>(sampled);
        s.nSampledRings = (s.nRings + stride - 1) / stride;
        s.W = float(kPi * d.yarns[y].ext);  // half a turn at the model's extent max(a, b)
        sampled += s.nSampledRings;
        sy.push_back(s);
    }
    const long long n = sampled * p.angles;
    DeviceBuffer syDev, depth, hitYarn, hitS, q, keys;
    syDev.upload(sy.data(), sy.size() * sizeof(sy[0]));
    depth.alloc(4ull * std::max(n, 1ll));
    hitYarn.alloc(4ull * std::max(n, 1ll));
    hitS.alloc(8ull * std::max(n, 1ll));
    q.alloc(12ull * std::max(n, 1ll));
    keys.alloc(16);
    check(cudaMemset(keys.data(), 0, 16), "GPU surface");
    if (n > 0) {
        collisionKernel<<<blocksFor(n), kThreadsPerBlock>>>(index, {static_cast<SampleYarn*>(syDev.data()), Y, n, stride, p.angles, float(p.depthTolerance)},
                                                            {depth.floats(), hitYarn.ints(), hitS.doubles(), q.floats()});
        collisionReduceKernel<<<blocksFor(n), kThreadsPerBlock>>>(depth.floats(), n, keys.ulongs(), keys.ulongs() + 1);
    }
    checkLaunch("GPU gate d collisions");
    unsigned long long result[2];
    keys.download(result, sizeof result);
    out.pointsTested = static_cast<long>(n);
    out.penetratingPoints = static_cast<long>(result[1]);
    if (result[1] > 0) {
        const unsigned bits = static_cast<unsigned>(result[0] >> 32);
        float dmax;
        std::memcpy(&dmax, &bits, 4);
        out.maxDepth = dmax;
        const long long worst = static_cast<long long>(0xFFFFFFFFull - (result[0] & 0xFFFFFFFFull));
        int hy = 0;
        double hs = 0;
        float hq[3];
        check(cudaMemcpy(&hy, hitYarn.ints() + worst, 4, cudaMemcpyDeviceToHost), "GPU surface");
        check(cudaMemcpy(&hs, hitS.doubles() + worst, 8, cudaMemcpyDeviceToHost), "GPU surface");
        check(cudaMemcpy(hq, q.floats() + 3 * worst, 12, cudaMemcpyDeviceToHost), "GPU surface");
        const long long sr = worst / p.angles;
        int y = 0;
        while (y + 1 < Y && sy[y + 1].sampleOff <= sr) ++y;
        const int k = static_cast<int>(sr - sy[y].sampleOff) * stride;
        out.worstYarnA = y;
        out.worstYarnB = hy;
        out.worstSA = k * d.yarns[y].ds;
        out.worstSB = hs;
        out.worstPoint = {hq[0], hq[1], hq[2]};
    }
    out.collisions = out.penetratingPoints > 0 ? Status::Fail : Status::Pass;
}

// The tube's ring plan (yarn_surface.cpp buildTubeMesh: uv wrap copies), in order on the host.
struct TubePlanRing {
    int copies;
    uint32_t first;
};
struct TubePlan {
    std::vector<TubeRing> rings;
    std::vector<std::vector<TubePlanRing>> yarns;
    uint32_t vertices = 0;
};

TubePlan planTube(const std::vector<YarnSurfaceRecord>& yarns, const TubeParams& tp, uint32_t stride, TubeMesh& m) {
    TubePlan plan;
    plan.yarns.resize(yarns.size());
    uint32_t at = 0;
    for (size_t y = 0; y < yarns.size(); ++y) {
        m.yarnFirstVertex.push_back(at);
        const double L = yarns[y].length;
        const int N = std::max(1, static_cast<int>(std::ceil(L / tp.ringSpacing)));
        double uOff = 0;
        const YarnUv yuv = yarnUvOf(tp, static_cast<int>(y));
        plan.yarns[y].resize(static_cast<size_t>(N) + 1);
        for (int j = 0; j <= N; ++j) {
            const double sj = L * j / N;
            const bool wrapHere = tp.uvWrap > 0 && j > 0 && j < N && sj / tp.uvRepeat / yuv.pitchFactor + yuv.phase - uOff >= tp.uvWrap;
            const int copies = wrapHere ? 2 : 1;
            plan.yarns[y][j] = {copies, at};
            for (int c = 0; c < copies; ++c)
                plan.rings.push_back({static_cast<int>(y), sj, c ? uOff + tp.uvWrap : uOff, at + c * stride, yuv.pitchFactor, yuv.phase});
            if (wrapHere) uOff += tp.uvWrap;
            at += copies * stride;
        }
    }
    m.yarnFirstVertex.push_back(at);
    plan.vertices = at;
    return plan;
}

// The vertices of the plan's rings on the device, into the mesh.
void evaluateTube(const DeviceSurfaceModel& d, const TubePlan& plan, const TubeParams& tp, int sides, TubeMesh& m) {
    const size_t nv = plan.vertices;
    DeviceBuffer rings, pos, nrm, uv, sArr, tArr, thArr, yArr, axis;
    rings.upload(plan.rings.data(), plan.rings.size() * sizeof(plan.rings[0]));
    pos.alloc(12 * nv + 1);
    nrm.alloc(12 * nv + 1);
    uv.alloc(8 * nv + 1);
    sArr.alloc(8 * nv + 1);
    tArr.alloc(8 * nv + 1);
    thArr.alloc(8 * nv + 1);
    yArr.alloc(4 * nv + 1);
    axis.alloc(4 * nv + 1);
    const int n = static_cast<int>(plan.rings.size());
    if (n > 0)
        tubeKernel<<<blocksFor(n), kThreadsPerBlock>>>({d.curves.view, d.surf.view, static_cast<TubeRing*>(rings.data()), n, sides, float(tp.uvRepeat), d.r},
                                                       {pos.floats(), nrm.floats(), uv.floats(), sArr.doubles(), tArr.doubles(), thArr.doubles(), yArr.ints(), axis.floats()});
    std::vector<float> hp, hn;
    checkLaunch("GPU tube mesh");
    pos.download(hp, 3 * nv);
    nrm.download(hn, 3 * nv);
    uv.download(m.uv, 2 * nv);
    sArr.download(m.s, nv);
    tArr.download(m.t, nv);
    thArr.download(m.theta, nv);
    yArr.download(m.yarn, nv);
    axis.download(m.axisRatio, nv);
    m.pos.resize(nv);
    m.nrm.resize(nv);
    for (size_t v = 0; v < nv; ++v) {
        m.pos[v] = {hp[3 * v], hp[3 * v + 1], hp[3 * v + 2]};
        m.nrm[v] = {hn[3 * v], hn[3 * v + 1], hn[3 * v + 2]};
    }
}

} // namespace

std::string SurfaceGpuStats::line() const {
    char b[256];
    std::snprintf(b, sizeof b, "GPU surface: %lld rings; build %.0f ms, index %.0f ms, check %.0f ms, download %.0f ms", rings, msBuild,
                  msIndex, msCheck, msDownload);
    return b;
}

bool buildSurface(const SurfaceRequest& request, SurfaceModel& out, SurfaceGpuStats* st, std::string* err) try {
    if (!available(err)) return false;
    Reservation lease(Reservation::exclusiveBytes());  // runs alone on the device (runtime.h)
    const auto t0 = Clock::now();
    const DeviceSurfaceModel d(request);
    const double msBuild = elapsedMs(t0);
    const auto t1 = Clock::now();
    out = d.download(request.pattern);
    out.ms = elapsedMs(t0);
    if (st) {
        st->rings += d.nRings;
        st->msBuild += msBuild;
        st->msDownload += elapsedMs(t1);
    }
    return true;
} catch (const DeviceError& e) {
    return reportFailure(e, err);
}

bool checkSurface(const SurfaceRequest& request, const SurfaceCheckParams& p, SurfaceCheck& out, SurfaceGpuStats* st, std::string* err) try {
    if (!available(err)) return false;
    Reservation lease(Reservation::exclusiveBytes());  // runs alone on the device (runtime.h)
    const auto t0 = Clock::now();
    out = SurfaceCheck{};
    const DeviceSurfaceModel d(request);
    const double msBuild = elapsedMs(t0);
    const auto t1 = Clock::now();
    const DeviceCapsuleIndex index(d);
    const double msIndex = elapsedMs(t1);
    const auto t2 = Clock::now();
    checkFolds(d, out);
    out.compressionClamped = request.surface.compression ? static_cast<int>(d.clamped) : 0;
    checkCollisions(d, index.query(d), p, out);
    char b[200];
    std::snprintf(b, sizeof b, "sampled: every %d model rings (%.3f mm), %d angles, depth tolerance %.0e mm [GPU float32]", p.ringStride,
                  d.yarns.empty() ? 0.0 : p.ringStride * d.yarns[0].ds, p.angles, p.depthTolerance);
    out.resolution = b;
    out.ms = elapsedMs(t0);
    if (st) {
        st->rings += d.nRings;
        st->msBuild += msBuild;
        st->msIndex += msIndex;
        st->msCheck += elapsedMs(t2);
    }
    return true;
} catch (const DeviceError& e) {
    return reportFailure(e, err);
}

bool buildTubeMesh(const SurfaceRequest& request, const TubeParams& tp, TubeMesh& m, SurfaceGpuStats* st, std::string* err) try {
    if (!available(err)) return false;
    Reservation lease(Reservation::exclusiveBytes());  // runs alone on the device (runtime.h)
    const auto t0 = Clock::now();
    const DeviceSurfaceModel d(request);
    const double msBuild = elapsedMs(t0);
    const auto t1 = Clock::now();
    m = TubeMesh{};
    const int sides = std::max(3, tp.sides);
    const uint32_t stride = sides + 1;
    const TubePlan plan = planTube(d.yarns, tp, stride, m);
    evaluateTube(d, plan, tp, sides, m);
    // Strips (the CPU's order).
    for (const std::vector<TubePlanRing>& rings : plan.yarns)
        for (size_t j = 1; j < rings.size(); ++j) {
            const uint32_t prev = rings[j - 1].first + (rings[j - 1].copies - 1) * stride, ringFirst = rings[j].first;
            for (int i = 0; i < sides; ++i) {
                const uint32_t a = prev + i, b = a + 1, c = ringFirst + i, e = c + 1;
                m.tri.insert(m.tri.end(), {a, c, b, b, c, e});
            }
        }
    m.ms = elapsedMs(t0);
    if (st) {
        st->rings += d.nRings;
        st->msBuild += msBuild;
        st->msDownload += elapsedMs(t1);
    }
    return true;
} catch (const DeviceError& e) {
    return reportFailure(e, err);
}

} // namespace knit::gpu
