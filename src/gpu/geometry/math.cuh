#pragma once
// GPU backend: device helpers shared by the kernels (single precision geometry).
//
// Precision policy: 3-D geometry, frames, distances and local integrals
// are float. One-dimensional yarn coordinates that grow with the yarn - the arc length s and
// the spline parameter u (a body yarn is ~52 m long: float would quantize s to 4 um, half the
// fibre candidate spacing) - and arc-length table entries are double scalars; u is split into
// (span, float local t) before any geometric evaluation. Integer hashing is exact.
//
// Guards: the CPU's 1e-300 / 1e-12 / 1e-15 guards become kTiny* constants that are finite in
// float (listed with each use).
#include <cstdint>
#include <cuda_runtime.h>

#include "gpu/device/kernel.cuh"
#include "gpu/geometry/device_geometry.h"

namespace knit::gpu {

// --- vectors
struct f3 {
    float x, y, z;
};
__host__ __device__ __forceinline__ f3 operator+(f3 a, f3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
__host__ __device__ __forceinline__ f3 operator-(f3 a, f3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
__host__ __device__ __forceinline__ f3 operator-(f3 a) { return {-a.x, -a.y, -a.z}; }
__host__ __device__ __forceinline__ f3 operator*(f3 a, float s) { return {a.x * s, a.y * s, a.z * s}; }
__host__ __device__ __forceinline__ f3 operator*(float s, f3 a) { return {a.x * s, a.y * s, a.z * s}; }
__host__ __device__ __forceinline__ f3 operator/(f3 a, float s) { return {a.x / s, a.y / s, a.z / s}; }
__host__ __device__ __forceinline__ f3& operator+=(f3& a, f3 b) { a.x += b.x; a.y += b.y; a.z += b.z; return a; }
__host__ __device__ __forceinline__ f3& operator-=(f3& a, f3 b) { a.x -= b.x; a.y -= b.y; a.z -= b.z; return a; }
__host__ __device__ __forceinline__ float dot(f3 a, f3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
__host__ __device__ __forceinline__ f3 cross(f3 a, f3 b) { return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }
__device__ __forceinline__ float norm(f3 a) { return sqrtf(dot(a, a)); }
__device__ __forceinline__ f3 normalized(f3 a) {
    const float n = norm(a);
    return n > 0 ? a / n : f3{0, 0, 0};
}
__host__ __device__ __forceinline__ f3 lerp(f3 a, f3 b, float t) { return a + (b - a) * t; }
__host__ __device__ __forceinline__ f3 fromF(const float* p) { return {p[0], p[1], p[2]}; }
__host__ __device__ __forceinline__ float clampf(float x, float lo, float hi) { return x < lo ? lo : (x > hi ? hi : x); }
__host__ __device__ __forceinline__ double clampd(double x, double lo, double hi) { return x < lo ? lo : (x > hi ? hi : x); }
__host__ __device__ __forceinline__ float fminf_(float a, float b) { return b < a ? b : a; }  // std::min semantics
__host__ __device__ __forceinline__ float fmaxf_(float a, float b) { return a < b ? b : a; }  // std::max semantics
__host__ __device__ __forceinline__ f3 vmin(f3 a, f3 b) { return {fminf_(a.x, b.x), fminf_(a.y, b.y), fminf_(a.z, b.z)}; }
__host__ __device__ __forceinline__ f3 vmax(f3 a, f3 b) { return {fmaxf_(a.x, b.x), fmaxf_(a.y, b.y), fmaxf_(a.z, b.z)}; }

constexpr float kPiF = 3.14159265358979323846f;
constexpr double kPiD = 3.14159265358979323846;
constexpr float kDegF = kPiF / 180.0f;
// Squared-length guard of segmentDistance (CPU: 1e-300 mm^2): a segment shorter than 1e-12 mm.
constexpr float kTinySq = 1e-24f;
// Normalization guards (CPU: 1e-12 mm or 1e-15 mm): float keeps 1e-12.
constexpr float kTinyLen = 1e-12f;

// --- B-spline
// Uniform cubic B-spline span k (controls k..k+3) at local t, as the CPU's CubicBezier
// (same Bezier control points, de Casteljau order).
struct Bezier {
    f3 b[4];
};
__device__ __forceinline__ Bezier spanBez(const float* ctrl, int k) {
    const f3 P0 = fromF(ctrl + 3 * k), P1 = fromF(ctrl + 3 * (k + 1)), P2 = fromF(ctrl + 3 * (k + 2)), P3 = fromF(ctrl + 3 * (k + 3));
    Bezier z;
    z.b[0] = (P0 + P1 * 4.0f + P2) / 6.0f;
    z.b[1] = (P1 * 2.0f + P2) / 3.0f;
    z.b[2] = (P1 + P2 * 2.0f) / 3.0f;
    z.b[3] = (P1 + P2 * 4.0f + P3) / 6.0f;
    return z;
}
__device__ __forceinline__ f3 bezEval(const Bezier& z, float t) {
    const f3 a = lerp(z.b[0], z.b[1], t), c = lerp(z.b[1], z.b[2], t), d = lerp(z.b[2], z.b[3], t);
    const f3 e = lerp(a, c, t), f = lerp(c, d, t);
    return lerp(e, f, t);
}
__device__ __forceinline__ f3 bezD1(const Bezier& z, float t) {
    const f3 h0 = (z.b[1] - z.b[0]) * 3.0f, h1 = (z.b[2] - z.b[1]) * 3.0f, h2 = (z.b[3] - z.b[2]) * 3.0f;
    return lerp(lerp(h0, h1, t), lerp(h1, h2, t), t);
}
__device__ __forceinline__ f3 bezD2(const Bezier& z, float t) {
    const f3 a = (z.b[2] - z.b[1] * 2.0f + z.b[0]) * 6.0f, b = (z.b[3] - z.b[2] * 2.0f + z.b[1]) * 6.0f;
    return lerp(a, b, t);
}
// BSpline::locate: clamped to [0, nSpans] (never extrapolated).
__device__ __forceinline__ int locate(double u, int nSpans, float& t) {
    if (u <= 0.0) { t = 0.0f; return 0; }
    if (u >= nSpans) { t = 1.0f; return nSpans - 1; }
    int k = static_cast<int>(floor(u));
    if (k >= nSpans) k = nSpans - 1;
    t = static_cast<float>(u - k);
    return k;
}

// One curve of a curve set on the device.
struct CurveView {
    const float* ctrl;   // 3 per control
    int nSpans;
    const double* tu;    // arc-length table nodes (u, s)
    const double* ts;
    int nTab;
    double length;
};
__device__ __forceinline__ CurveView devCurve(const CurveSetView& c, int y) {
    const CurveRecord& fc = c.curves[y];
    return {c.ctrl + 3 * fc.ctrlOff, fc.nCtrl - 3, c.tu + fc.tabOff, c.ts + fc.tabOff, fc.nTab, fc.length};
}
// YarnPath::thicknessAt.
__device__ __forceinline__ float thicknessAt(const CurveSetView& cs, int y, double u) {
    if (!cs.thickness) return 1;
    const CurveRecord& fc = cs.curves[y];
    const double c = clampd(u + 1, 0.0, static_cast<double>(fc.nCtrl - 1));
    const int k = static_cast<int>(c);
    const float f = static_cast<float>(c - k);
    const float* t = cs.thickness + fc.ctrlOff;
    return k + 1 < fc.nCtrl ? t[k] * (1 - f) + t[k + 1] * f : t[k];
}
__device__ __forceinline__ f3 curveEval(const CurveView& c, double u) {
    float t;
    const int k = locate(u, c.nSpans, t);
    return bezEval(spanBez(c.ctrl, k), t);
}
__device__ __forceinline__ f3 curveD1(const CurveView& c, double u) {
    float t;
    const int k = locate(u, c.nSpans, t);
    return bezD1(spanBez(c.ctrl, k), t);
}
__device__ __forceinline__ f3 curveD2(const CurveView& c, double u) {
    float t;
    const int k = locate(u, c.nSpans, t);
    return bezD2(spanBez(c.ctrl, k), t);
}

// QUADPACK G7/K15 of the speed over [a, b] inside ONE span k (local parameters, float).
__device__ __forceinline__ float gk15Local(const Bezier& z, float a, float b) {
    const float xgk[8] = {0.991455371120812639f, 0.949107912342758525f, 0.864864423359769073f, 0.741531185599394440f,
                          0.586087235467691130f, 0.405845151377397167f, 0.207784955007898468f, 0.0f};
    const float wgk[8] = {0.022935322010529225f, 0.063092092629978553f, 0.104790010322250184f, 0.140653259715525919f,
                          0.169004726639267903f, 0.190350578064785410f, 0.204432940075298892f, 0.209482141084727828f};
    const float h = 0.5f * (b - a), m = 0.5f * (a + b);
    float K = norm(bezD1(z, m)) * wgk[7];
    for (int j = 0; j < 7; ++j) K += wgk[j] * (norm(bezD1(z, m - h * xgk[j])) + norm(bezD1(z, m + h * xgk[j])));
    return K * h;
}

// YarnPath::paramAtArcLength: the table interval by binary search (double), then a
// safeguarded Newton step on the local integral. The table leaves lie inside one span, so
// the integral is a float evaluation in local parameters. Tolerance: |f| <= 1e-6 mm or 30
// iterations (CPU: 1e-12 L or 60 iterations; float cannot reach the CPU's tolerance).
__device__ inline double paramAtArcLength(const CurveView& c, double s) {
    if (s <= 0) return 0.0;
    if (s >= c.length) return static_cast<double>(c.nSpans);
    // upper_bound(ts, s) - 1
    int lo = 0, hi = c.nTab;  // first index with ts > s in [lo, hi)
    while (lo < hi) {
        const int mid = (lo + hi) >> 1;
        if (c.ts[mid] > s) hi = mid; else lo = mid + 1;
    }
    int i = lo - 1;
    if (i > c.nTab - 2) i = c.nTab - 2;
    if (i < 0) i = 0;
    const double a = c.tu[i], b = c.tu[i + 1];
    int k = static_cast<int>(floor(a));
    if (k > c.nSpans - 1) k = c.nSpans - 1;
    const Bezier z = spanBez(c.ctrl, k);
    const float ta = static_cast<float>(a - k), tb = static_cast<float>(b - k);
    const double ds = c.ts[i + 1] - c.ts[i];
    const float target = static_cast<float>(s - c.ts[i]);  // small: inside one table leaf
    float lo2 = ta, hi2 = tb;
    float t = ta + (tb - ta) * static_cast<float>((s - c.ts[i]) / (ds > 1e-300 ? ds : 1e-300));
    for (int it = 0; it < 30; ++it) {
        const float f = (t > ta ? gk15Local(z, ta, t) : 0.0f) - target;
        if (fabsf(f) <= 1e-6f) break;
        if (f > 0) hi2 = t; else lo2 = t;
        const float sp = norm(bezD1(z, t));
        float tn = sp > 0 ? t - f / sp : 0.5f * (lo2 + hi2);
        if (!(tn > lo2 && tn < hi2)) tn = 0.5f * (lo2 + hi2);
        if (tn == t) break;
        t = tn;
    }
    return k + static_cast<double>(t);
}

// --- distances
// Ericson 5.1.9 (model/mesh.cpp closestSegmentPoints), float, kTinySq for 1e-300.
__device__ __forceinline__ float segmentDistance(f3 p1, f3 q1, f3 p2, f3 q2, float& s, float& t) {
    const f3 d1 = q1 - p1, d2 = q2 - p2, r = p1 - p2;
    const float a = dot(d1, d1), e = dot(d2, d2), f = dot(d2, r);
    if (a <= kTinySq && e <= kTinySq) { s = t = 0; return norm(p1 - p2); }
    if (a <= kTinySq) {
        s = 0;
        t = clampf(f / e, 0.0f, 1.0f);
    } else {
        const float c = dot(d1, r);
        if (e <= kTinySq) {
            t = 0;
            s = clampf(-c / a, 0.0f, 1.0f);
        } else {
            const float b = dot(d1, d2), denom = a * e - b * b;
            s = denom > 0 ? clampf((b * f - c * e) / denom, 0.0f, 1.0f) : 0.0f;
            t = (b * s + f) / e;
            if (t < 0) { t = 0; s = clampf(-c / a, 0.0f, 1.0f); }
            else if (t > 1) { t = 1; s = clampf((b - c) / a, 0.0f, 1.0f); }
        }
    }
    return norm((p1 + d1 * s) - (p2 + d2 * t));
}
// Point to segment (segmentDistance(q, q, a, b)).
__device__ __forceinline__ float pointSegment(f3 q, f3 a, f3 b, float& t) {
    float s;
    return segmentDistance(q, q, a, b, s, t);
}

// Ericson 5.1.5 (mesh.cpp closestOnTriangle).
__device__ __forceinline__ f3 closestOnTriangle(f3 p, f3 a, f3 b, f3 c, float bary[3]) {
    const f3 ab = b - a, ac = c - a, ap = p - a;
    const float d1 = dot(ab, ap), d2 = dot(ac, ap);
    if (d1 <= 0 && d2 <= 0) { bary[0] = 1; bary[1] = bary[2] = 0; return a; }
    const f3 bp = p - b;
    const float d3 = dot(ab, bp), d4 = dot(ac, bp);
    if (d3 >= 0 && d4 <= d3) { bary[1] = 1; bary[0] = bary[2] = 0; return b; }
    const float vc = d1 * d4 - d3 * d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) {
        const float v = d1 / (d1 - d3);
        bary[0] = 1 - v; bary[1] = v; bary[2] = 0;
        return a + ab * v;
    }
    const f3 cp = p - c;
    const float d5 = dot(ab, cp), d6 = dot(ac, cp);
    if (d6 >= 0 && d5 <= d6) { bary[2] = 1; bary[0] = bary[1] = 0; return c; }
    const float vb = d5 * d2 - d1 * d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) {
        const float w = d2 / (d2 - d6);
        bary[0] = 1 - w; bary[1] = 0; bary[2] = w;
        return a + ac * w;
    }
    const float va = d3 * d6 - d5 * d4;
    if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
        const float w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        bary[0] = 0; bary[1] = 1 - w; bary[2] = w;
        return b + (c - b) * w;
    }
    const float denom = 1.0f / (va + vb + vc);
    const float v = vb * denom, w = vc * denom;
    bary[0] = 1 - v - w; bary[1] = v; bary[2] = w;
    return a + ab * v + ac * w;
}

// --- hashing
__host__ __device__ __forceinline__ uint64_t splitmix(uint64_t x) {
    x += 0x9e3779b97f4a7c15ull;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
    return x ^ (x >> 31);
}
// model/hash.h hash01: exact (53-bit uniform as a double).
__host__ __device__ __forceinline__ double hash01(uint64_t seed, int yarn, int bin, int channel) {
    uint64_t h = splitmix(seed);
    h = splitmix(h ^ (uint64_t(uint32_t(yarn)) << 32 | uint32_t(bin)));
    h = splitmix(h ^ uint64_t(channel));
    return double(h >> 11) * (1.0 / 9007199254740992.0);
}

// --- grid
// The items of a uniform grid (GridView) whose boxes overlap [lo, hi]: each item once, in the first
// cell it shares with the query, as SpatialHash3::query does.
//     for (GridCursor cursor(grid, lo, hi); cursor.next();) {
//         const int id = cursor.item();
//         ...                            // continue: the next item; break: stop early
//     }
class GridCursor {
public:
    __host__ __device__ GridCursor(const GridView& g, f3 lo, f3 hi) : g_(g), lo_(lo), hi_(hi) {
        if (g.nItems == 0) {
            done_ = true;
            return;
        }
        // The CPU divides (floor(v / cell)); v * (1 / cell) differs in the last ulp only at cell
        // boundaries, where both sides are conservative (the box overlap test decides).
        const float inv = 1.0f / g.cell;
        qi_ = cellOf(lo.x, inv);
        qj_ = cellOf(lo.y, inv);
        qk_ = cellOf(lo.z, inv);
        i_ = i0_ = qi_ > g.lo[0] ? qi_ : g.lo[0];
        j_ = j0_ = qj_ > g.lo[1] ? qj_ : g.lo[1];
        k0_ = qk_ > g.lo[2] ? qk_ : g.lo[2];
        k_ = k0_ - 1;  // next() starts by stepping into the first cell
        i1_ = cellOf(hi.x, inv);
        j1_ = cellOf(hi.y, inv);
        k1_ = cellOf(hi.z, inv);
        if (i1_ > g.hi[0]) i1_ = g.hi[0];
        if (j1_ > g.hi[1]) j1_ = g.hi[1];
        if (k1_ > g.hi[2]) k1_ = g.hi[2];
        done_ = i0_ > i1_ || j0_ > j1_ || k0_ > k1_;
    }

    // Moves to the next item; false when there is none.
    __host__ __device__ __forceinline__ bool next() {
        if (done_) return false;
        for (;;) {
            while (at_ < end_) {
                const unsigned entry = g_.entries[at_++];
                if (((entry >> 29) | atQueryLow_) != 7u) continue;
                item_ = static_cast<int>(entry & ((1u << 29) - 1));
                const f3 a = fromF(g_.boxLo + 3 * item_), b = fromF(g_.boxHi + 3 * item_);
                if (a.x > hi_.x || b.x < lo_.x || a.y > hi_.y || b.y < lo_.y || a.z > hi_.z || b.z < lo_.z) continue;
                return true;
            }
            if (!nextCell()) {
                done_ = true;
                return false;
            }
        }
    }
    __host__ __device__ __forceinline__ int item() const { return item_; }

private:
    __host__ __device__ static long long cellOf(float v, float inv) { return static_cast<long long>(floorf(v * inv)); }
    __host__ __device__ __forceinline__ bool nextCell() {
        if (++k_ > k1_) {
            k_ = k0_;
            if (++j_ > j1_) {
                j_ = j0_;
                if (++i_ > i1_) return false;
            }
        }
        const long long nj = g_.hi[1] - g_.lo[1] + 1, nk = g_.hi[2] - g_.lo[2] + 1;
        const long long c = ((i_ - g_.lo[0]) * nj + (j_ - g_.lo[1])) * nk + (k_ - g_.lo[2]);
        at_ = g_.cellStart[c];
        end_ = g_.cellStart[c + 1];
        atQueryLow_ = (i_ == qi_ ? 1u : 0u) | (j_ == qj_ ? 2u : 0u) | (k_ == qk_ ? 4u : 0u);
        return true;
    }

    const GridView& g_;
    f3 lo_, hi_;
    long long qi_ = 0, qj_ = 0, qk_ = 0;  // the cell of the query's low corner
    long long i_ = 0, j_ = 0, k_ = 0, i0_ = 0, j0_ = 0, k0_ = 0, i1_ = -1, j1_ = -1, k1_ = -1;
    unsigned at_ = 0, end_ = 0, atQueryLow_ = 0;
    int item_ = -1;
    bool done_ = false;
};

// --- local spans
// A span's Bezier controls are formed from the B-spline controls RELATIVE to P1 (P_i - P1 is exact
// in float for nearby controls, Sterbenz), so derivatives carry only the rounding of the stored
// float controls, not that of the absolute coordinates (a teddy body has coordinates up to ~100 mm).
// Positions add the origin back (one rounding).
struct LocalSpan {
    Bezier z;     // local Bezier controls (relative to `o`)
    f3 o;      // origin (the span's control P1)
};

__device__ __forceinline__ LocalSpan localSpan(const float* ctrl, int k) {
    const f3 P0 = fromF(ctrl + 3 * k), P1 = fromF(ctrl + 3 * (k + 1)), P2 = fromF(ctrl + 3 * (k + 2)), P3 = fromF(ctrl + 3 * (k + 3));
    const f3 q0 = P0 - P1, q2 = P2 - P1, q3 = P3 - P1;
    LocalSpan s;
    s.o = P1;
    s.z.b[0] = (q0 + q2) / 6.0f;
    s.z.b[1] = q2 / 3.0f;
    s.z.b[2] = (q2 * 2.0f) / 3.0f;
    s.z.b[3] = (q2 * 4.0f + q3) / 6.0f;
    return s;
}

// CubicBezier::sub (blossom at (t0^(3-i), t1^i)).
__device__ __forceinline__ Bezier bezSub(const Bezier& z, float t0, float t1) {
    auto blossom = [&](float s, float t, float v) {
        const f3 a = lerp(z.b[0], z.b[1], s), c = lerp(z.b[1], z.b[2], s), d = lerp(z.b[2], z.b[3], s);
        const f3 e = lerp(a, c, t), f = lerp(c, d, t);
        return lerp(e, f, v);
    };
    Bezier r;
    r.b[0] = blossom(t0, t0, t0);
    r.b[1] = blossom(t0, t0, t1);
    r.b[2] = blossom(t0, t1, t1);
    r.b[3] = blossom(t1, t1, t1);
    return r;
}

__device__ __forceinline__ float polygonLength(const Bezier& z) {
    return norm(z.b[1] - z.b[0]) + norm(z.b[2] - z.b[1]) + norm(z.b[3] - z.b[2]);
}

// YarnPath::spanBounds (splits 8): chord sum (lo) and control-polygon length (hi) of [t0, t1].
__device__ __forceinline__ void spanBounds(const float* ctrl, int k, float t0, float t1, float& lo, float& hi) {
    const Bezier piece = bezSub(localSpan(ctrl, k).z, t0, t1);
    hi = polygonLength(piece);
    lo = 0;
    f3 prev = piece.b[0];
    for (int j = 1; j <= 8; ++j) {
        const f3 p = bezEval(piece, j / 8.0f);
        lo += norm(p - prev);
        prev = p;
    }
}

// Float allowance of a float spanBounds pair (CPU: 1e-12 relative): 4e-6 relative (a few
// dozen float roundings of the eight chords / three polygon legs).
constexpr double kBoundRel = 4e-6;

// YarnPath::lengthBounds(u0, u1), 8 chords per partial span (the validation's only use): the
// partial end spans in float, whole spans from the double prefix sums when present.
__device__ inline ArcLengthBounds lengthBounds(const CurveView& c, const double* preLo, const double* preHi, double u0, double u1) {
    ArcLengthBounds out;
    if (u1 <= u0) return out;
    u0 = u0 > 0 ? u0 : 0;
    u1 = u1 < c.nSpans ? u1 : c.nSpans;
    const int n = c.nSpans;
    int k0 = static_cast<int>(floor(u0));
    if (k0 > n - 1) k0 = n - 1;
    if (!preLo) {
        for (int k = k0; k < n && k < u1; ++k) {
            const double t0 = u0 - k > 0 ? u0 - k : 0.0, t1 = u1 - k < 1 ? u1 - k : 1.0;
            if (t1 <= t0) continue;
            float l, h;
            spanBounds(c.ctrl, k, static_cast<float>(t0), static_cast<float>(t1), l, h);
            out.lo += l;
            out.hi += h;
        }
        out.lo *= (1.0 - kBoundRel);
        out.hi *= (1.0 + kBoundRel);
        return out;
    }
    int k1 = static_cast<int>(ceil(u1)) - 1;
    k1 = k1 < k0 ? k0 : (k1 > n - 1 ? n - 1 : k1);
    float l, h;
    if (k0 == k1) {
        spanBounds(c.ctrl, k0, static_cast<float>(u0 - k0), static_cast<float>(u1 - k0), l, h);
        out.lo = l * (1.0 - kBoundRel);
        out.hi = h * (1.0 + kBoundRel);
        return out;
    }
    spanBounds(c.ctrl, k0, static_cast<float>(u0 - k0), 1.0f, l, h);
    double pl = l, ph = h;
    spanBounds(c.ctrl, k1, 0.0f, static_cast<float>(u1 - k1), l, h);
    pl += l;
    ph += h;
    out.lo = (pl + (preLo[k1] - preLo[k0 + 1])) * (1.0 - kBoundRel) - 1e-13 * preLo[k1];
    out.hi = (ph + (preHi[k1] - preHi[k0 + 1])) * (1.0 + kBoundRel) + 1e-13 * preHi[k1];
    out.lo = out.lo > 0 ? out.lo : 0;
    return out;
}

// YarnPath::arcLengthAt: table node by upper_bound, plus the local G7/K15 over [tu[i], u]
// (inside one span: the table's leaves never cross a knot) in float.
__device__ inline double arcLengthAt(const CurveView& c, double u) {
    const double uEnd = c.nSpans;
    u = u < 0 ? 0 : (u > uEnd ? uEnd : u);
    int lo = 0, hi = c.nTab;  // first index with tu > u
    while (lo < hi) {
        const int mid = (lo + hi) >> 1;
        if (c.tu[mid] > u) hi = mid; else lo = mid + 1;
    }
    const int i = lo == 0 ? 0 : lo - 1;
    if (i >= c.nTab - 1) return c.length;
    if (!(u > c.tu[i])) return c.ts[i];
    int k = static_cast<int>(floor(c.tu[i]));
    if (k > c.nSpans - 1) k = c.nSpans - 1;
    const Bezier z = localSpan(c.ctrl, k).z;
    return c.ts[i] + gk15Local(z, static_cast<float>(c.tu[i] - k), static_cast<float>(u - k));
}

// Curve position (absolute, float) at global parameter u.
__device__ __forceinline__ f3 evalAt(const CurveView& c, double u) {
    float t;
    const int k = locate(u, c.nSpans, t);
    const LocalSpan s = localSpan(c.ctrl, k);
    return bezEval(s.z, t) + s.o;
}

// --- double precision: the solver's row values
struct d3 {
    double x, y, z;
};
__host__ __device__ __forceinline__ d3 operator+(d3 a, d3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
__host__ __device__ __forceinline__ d3 operator-(d3 a, d3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
__host__ __device__ __forceinline__ d3 operator*(d3 a, double s) { return {a.x * s, a.y * s, a.z * s}; }
__host__ __device__ __forceinline__ d3 operator/(d3 a, double s) { return {a.x / s, a.y / s, a.z / s}; }
__host__ __device__ __forceinline__ d3& operator+=(d3& a, d3 b) { a.x += b.x; a.y += b.y; a.z += b.z; return a; }
__host__ __device__ __forceinline__ double dot(d3 a, d3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
__host__ __device__ __forceinline__ d3 cross(d3 a, d3 b) { return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }
__device__ __forceinline__ double norm(d3 a) { return sqrt(dot(a, a)); }
__device__ __forceinline__ double comp(d3 a, int i) { return i == 0 ? a.x : i == 1 ? a.y : a.z; }
__device__ __forceinline__ f3 tof(d3 a) { return {static_cast<float>(a.x), static_cast<float>(a.y), static_cast<float>(a.z)}; }

// model/mesh.cpp closestSegmentPoints (Ericson 5.1.9) in double, the CPU's 1e-300 guard.
__device__ inline double segmentDistance(d3 p1, d3 q1, d3 p2, d3 q2, double& s, double& t) {
    const d3 d1 = q1 - p1, d2 = q2 - p2, r = p1 - p2;
    const double a = dot(d1, d1), e = dot(d2, d2), f = dot(d2, r);
    if (a <= 1e-300 && e <= 1e-300) {
        s = t = 0;
        return norm(p1 - p2);
    }
    if (a <= 1e-300) {
        s = 0;
        t = clampd(f / e, 0.0, 1.0);
    } else {
        const double c = dot(d1, r);
        if (e <= 1e-300) {
            t = 0;
            s = clampd(-c / a, 0.0, 1.0);
        } else {
            const double b = dot(d1, d2), denom = a * e - b * b;
            s = denom > 0 ? clampd((b * f - c * e) / denom, 0.0, 1.0) : 0.0;
            t = (b * s + f) / e;
            if (t < 0) {
                t = 0;
                s = clampd(-c / a, 0.0, 1.0);
            } else if (t > 1) {
                t = 1;
                s = clampd((b - c) / a, 0.0, 1.0);
            }
        }
    }
    return norm((p1 + d1 * s) - (p2 + d2 * t));
}

// Uniform cubic B-spline basis and derivatives (tests/ncp_reference.cpp basis()).
__device__ __forceinline__ void bsplineBasis(double t, double w[4], double d1[4], double d2[4]) {
    const double s = 1 - t;
    w[0] = s * s * s / 6;
    w[1] = (3 * t * t * t - 6 * t * t + 4) / 6;
    w[2] = (-3 * t * t * t + 3 * t * t + 3 * t + 1) / 6;
    w[3] = t * t * t / 6;
    d1[0] = -s * s / 2;
    d1[1] = (3 * t * t - 4 * t) / 2;
    d1[2] = (-3 * t * t + 2 * t + 1) / 2;
    d1[3] = t * t / 2;
    d2[0] = s;
    d2[1] = 3 * t - 2;
    d2[2] = -3 * t + 1;
    d2[3] = t;
}

} // namespace knit::gpu
