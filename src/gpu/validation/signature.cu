// Gate (a): checkLocalSignature on the device (checks.h).
//
// The CPU algorithm (src/knitting/checks/signature.cpp checkLocalSignature, crossings.cpp findCrossings,
// signature.cpp checkSignature) in single precision, in three kernels:
//   0. per item (a stitch with parents, or one CLOSE drawstring pass): the mini set's controls,
//      extracted as the CPU does and transformed into the frame stitch's (origin, U, V, N) in
//      double from the host's double controls, stored as float; its extent L_ref;
//   1. per (item, top-level span pair) in the CPU's loop order: the Self / Adjacent / Distinct
//      recursion with an explicit stack in the CPU's depth-first order, hull rejection, cones,
//      Newton and the Krawczyk test (float intervals with directed rounding), depth contraction;
//      every certified root becomes a CANDIDATE keyed (pair index, order inside the pair);
//   2. per item: the candidates in key order (= the CPU's order of verifyAndRecord) through the
//      CPU's deduplication, then checkSignature on the mini set and the LocalPairResult
//      accounting (or the CLOSE pass rule).
// The deduplication only depends on earlier crossings, the verification does not, so splitting
// verification (1) from deduplication (2) gives the CPU's sequence.
//
// Precision (CPU value -> GPU value; L_ref = the mini set's extent, ~10-30 mm). LocalCurves controls
// are the double controls transformed into the frame in double, then stored as float. Curve
// parameters stay double; geometry (enclosures, hodographs, cones, Jacobians) is float. The
// residual F = A(a) - B(b) of Newton and of the Krawczyk test is evaluated in DOUBLE from the
// float controls (evalAtD): the test then certifies roots of the float-control curve, and
// near-tangent crossings certify as on the CPU (with float F they failed at every radius and the
// recursion exploded to the budget on the ears).
//   etaX (float enclosures)    1e-12 L_ref        -> 2e-7 L_ref
//   F pad (double evaluation)  2 etaX             -> 2e-11 L_ref
//   Krawczyk derivative pad    1e-12 |J| + 1e-15 L -> 1e-5 |J| + 1e-6 L
//   angular hull min norm      1e-9 L_ref         -> 1e-6 L_ref;  angle tolerance 1e-9 -> 2e-6 rad
//   Newton |F| (loop / final)  1e-13 / 1e-12 L    -> unchanged (F in double); float iterations first
//                                                    until |F| <= 4e-6 L_ref; det guard 1e-300 -> 1e-30
//   Krawczyk radii             1e-8 .. 1e-4       -> 1e-7, 1e-6, 1e-5, 1e-4, 3e-4, 1e-3 (parameter units)
//   deduplication pad          1e-12              -> 1e-6
//   self piece min width       1e-9               -> 1e-6
//   max depth                  40                 -> 24 (float local span parameters: 2^-24)
//   uncertified transversal    subdivide to depth -> UNRESOLVED below width 1e-4 (as non-transversal):
//                                                    a near-tangent crossing the float pads cannot certify
//   work budget                4e6 per findCrossings -> 2.5e5 per top-level span pair (and a 96-deep stack)
// Budget, stack, candidate or crossing overflow -> the item is UNRESOLVED and its detail says
// "GPU work budget exhausted" (never a silent PASS). The candidate buffer grows and the
// crossing kernel reruns if it overflows. These are GPU diagnostic verdicts.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>


#include "gpu/geometry/math.cuh"
#include "gpu/device/primitives.h"
#include "gpu/device/runtime.h"
#include "gpu/validation/checks.h"

namespace knit::gpu {

namespace {


constexpr float kEtaRel = 2e-7f;
constexpr float kDerRel = 1e-5f, kDerAbs = 1e-6f;
constexpr float kMinNormRel = 1e-6f, kAngTol = 2e-6f;
constexpr double kNewtonIn = 1e-13, kNewtonOut = 1e-12;  // F in double (evalAtD): the CPU's tolerances
constexpr double kEtaF = 1e-11;  // F pad (double evaluation of the float-control curve) [L_ref]
constexpr double kSelfMinWidth = 1e-6, kMinNonTransversalWidth = 1e-4, kDedupPad = 1e-6;
constexpr int kMaxDepth = 24;
constexpr int kBudget = 250000;
constexpr int kStack = 96;
constexpr int kMaxCross = 48;    // accepted crossings per item
constexpr int kSeqBits = 4;      // candidates per top-level pair: < 16

// Item flags.
constexpr unsigned kFlagUnres = 1, kFlagBudget = 2;

// --- float intervals
struct FloatInterval {
    float lo, hi;
};
__device__ __forceinline__ FloatInterval fi(float v) { return {v, v}; }
__device__ __forceinline__ FloatInterval iadd(FloatInterval a, FloatInterval b) { return {__fadd_rd(a.lo, b.lo), __fadd_ru(a.hi, b.hi)}; }
__device__ __forceinline__ FloatInterval isub(FloatInterval a, FloatInterval b) { return {__fsub_rd(a.lo, b.hi), __fsub_ru(a.hi, b.lo)}; }
__device__ __forceinline__ FloatInterval ineg(FloatInterval a) { return {-a.hi, -a.lo}; }
__device__ __forceinline__ FloatInterval imul(FloatInterval a, FloatInterval b) {
    const float l = fminf(fminf(__fmul_rd(a.lo, b.lo), __fmul_rd(a.lo, b.hi)), fminf(__fmul_rd(a.hi, b.lo), __fmul_rd(a.hi, b.hi)));
    const float h = fmaxf(fmaxf(__fmul_ru(a.lo, b.lo), __fmul_ru(a.lo, b.hi)), fmaxf(__fmul_ru(a.hi, b.lo), __fmul_ru(a.hi, b.hi)));
    return {l, h};
}
__device__ __forceinline__ FloatInterval ipad(FloatInterval a, float e) { return {__fsub_rd(a.lo, e), __fadd_ru(a.hi, e)}; }
__device__ __forceinline__ float imag(FloatInterval a) { return fmaxf(fabsf(a.lo), fabsf(a.hi)); }
__device__ __forceinline__ bool strictlyInside(FloatInterval a, FloatInterval o) { return a.lo > o.lo && a.hi < o.hi; }
__device__ __forceinline__ bool containsZero(FloatInterval a) { return a.lo <= 0.0f && 0.0f <= a.hi; }
__device__ __forceinline__ FloatInterval iintersect(FloatInterval a, FloatInterval b) { return {fmaxf(a.lo, b.lo), fminf(a.hi, b.hi)}; }

// Double parameter intervals (the CPU's Interval on curve parameters).
struct DoubleInterval {
    double lo, hi;
};
__device__ __forceinline__ DoubleInterval dintersect(DoubleInterval a, DoubleInterval b) { return {a.lo > b.lo ? a.lo : b.lo, a.hi < b.hi ? a.hi : b.hi}; }
__device__ __forceinline__ bool doverlaps(DoubleInterval a, DoubleInterval b) { return a.lo <= b.hi && b.lo <= a.hi; }
__device__ __forceinline__ DoubleInterval dhull(DoubleInterval a, DoubleInterval b) { return {a.lo < b.lo ? a.lo : b.lo, a.hi > b.hi ? a.hi : b.hi}; }
// Offsets of a global parameter interval from the centre c, rounded outward to float.
__device__ __forceinline__ FloatInterval offsetOf(DoubleInterval x, double c) { return {__double2float_rd(x.lo - c), __double2float_ru(x.hi - c)}; }

// --- the mini set
struct SignatureItem {
    int type;                       // 0: stitch pair, 1: CLOSE pass
    int ctrlA, nCtrlA, ctrlB, nCtrlB;  // into the mini control buffer (float, 3 per control)
    int srcA, srcB;                 // first source control in the double control buffer
    double frame[12];               // origin, U, V, N
    long long pairBase;             // first global top-level pair
    int spanOff, nSpans;            // yarn-1 stitch spans (SignatureSpan), sorted by uBegin
    // pair: the child (mini stitch 0)
    int op, incSlot, face, nParents, parentMini[2];
    // CLOSE: the strand's offset s0 and nominal range, the gathered loop's mini id, k parity
    double s0, segLo, segHi;
    int headMini, even;
};
struct SignatureSpan {
    double uBegin, uEnd;
    int id;
};
struct CrossingCandidate {
    int ya, yb, status, front, stA, stB;
    double uA, uB;
    DoubleInterval boxA, boxB;
    FloatInterval depth;
};
struct SignatureOutcome {
    int status, crossings, unexpected, overlap, neighbourhoodSelf, wrongOrder, budget;
    int pairStatus[2], pairCode[2], pairFound[2];
    int closeCode;                  // 0 ok, 1 count, 2 order
};

struct LocalCurves {
    const float* ctrl[2];
    int nSpans[2];
    float lref, etaX;
    const SignatureItem* it;
    const SignatureSpan* spans;
};

__device__ __forceinline__ Bezier spanOf(const LocalCurves& m, int y, int k) { return spanBez(m.ctrl[y], k); }
__device__ __forceinline__ f3 evalAt(const LocalCurves& m, int y, double u) {
    float t;
    const int k = locate(u, m.nSpans[y], t);
    return bezEval(spanOf(m, y, k), t);
}
// The curve point in double from the (float) mini controls: F = A(a) - B(b) of Newton and of the
// Krawczyk test is the exact (to double rounding) value of the float-control curve, so the test
// certifies roots of the curve the GPU validates; enclosures and derivatives stay float.
__device__ d3 evalAtD(const LocalCurves& m, int y, double u) {
    const int n = m.nSpans[y];
    int k;
    double t;
    if (u <= 0.0) { k = 0; t = 0.0; }
    else if (u >= n) { k = n - 1; t = 1.0; }
    else { k = static_cast<int>(floor(u)); if (k >= n) k = n - 1; t = u - k; }
    const float* c = m.ctrl[y] + 3 * k;
    d3 P[4], b[4];
    for (int q = 0; q < 4; ++q) P[q] = {double(c[3 * q]), double(c[3 * q + 1]), double(c[3 * q + 2])};
    auto mix = [](d3 a, d3 b2, d3 c2, double wa, double wb, double wc, double d) {
        return d3{(a.x * wa + b2.x * wb + c2.x * wc) / d, (a.y * wa + b2.y * wb + c2.y * wc) / d, (a.z * wa + b2.z * wb + c2.z * wc) / d};
    };
    b[0] = mix(P[0], P[1], P[2], 1, 4, 1, 6);
    b[1] = mix(P[1], P[2], P[2], 2, 1, 0, 3);
    b[2] = mix(P[1], P[2], P[2], 1, 2, 0, 3);
    b[3] = mix(P[1], P[2], P[3], 1, 4, 1, 6);
    for (int r = 3; r > 0; --r)
        for (int q = 0; q < r; ++q) b[q] = {b[q].x + (b[q + 1].x - b[q].x) * t, b[q].y + (b[q + 1].y - b[q].y) * t, b[q].z + (b[q + 1].z - b[q].z) * t};
    return b[0];
}
__device__ __forceinline__ f3 d1At(const LocalCurves& m, int y, double u) {
    float t;
    const int k = locate(u, m.nSpans[y], t);
    return bezD1(spanOf(m, y, k), t);
}

// crossings.cpp pieceRange: the span pieces of [u0, u1], each with its local sub-interval.
struct SpanPieces {
    int first, last;
    double u0, u1;
};
__device__ SpanPieces spanPieces(const LocalCurves& m, int y, double u0, double u1) {
    const int n = m.nSpans[y];
    u0 = u0 > 0.0 ? u0 : 0.0;
    u1 = u1 < double(n) ? u1 : double(n);
    int k0 = static_cast<int>(floor(u0));
    k0 = k0 < n - 1 ? k0 : n - 1;
    int k1 = static_cast<int>(floor(u1));
    k1 = k1 < n - 1 ? k1 : n - 1;
    if (u1 == floor(u1) && k1 > k0) k1 = static_cast<int>(u1) - 1;
    return {k0, k1, u0, u1};
}
__device__ __forceinline__ void pieceAt(const SpanPieces& p, int k, float& t0, float& t1) {
    t0 = static_cast<float>(p.u0 - k > 0.0 ? p.u0 - k : 0.0);
    t1 = static_cast<float>(p.u1 - k < 1.0 ? p.u1 - k : 1.0);
}
__device__ __forceinline__ f3 blossom3(const Bezier& z, float s, float t, float v) {
    const f3 a = lerp(z.b[0], z.b[1], s), c = lerp(z.b[1], z.b[2], s), d = lerp(z.b[2], z.b[3], s);
    const f3 e = lerp(a, c, t), f = lerp(c, d, t);
    return lerp(e, f, v);
}
__device__ __forceinline__ f3 blossom2(f3 h0, f3 h1, f3 h2, float s, float t) {
    const f3 a = lerp(h0, h1, s), b = lerp(h1, h2, s);
    return lerp(a, b, t);
}

struct Box3F {
    FloatInterval x, y, z;
};
__device__ void growBox(Box3F& b, bool& first, f3 v) {
    if (first) {
        b = {fi(v.x), fi(v.y), fi(v.z)};
        first = false;
        return;
    }
    b.x = {fminf(b.x.lo, v.x), fmaxf(b.x.hi, v.x)};
    b.y = {fminf(b.y.lo, v.y), fmaxf(b.y.hi, v.y)};
    b.z = {fminf(b.z.lo, v.z), fmaxf(b.z.hi, v.z)};
}
__device__ Box3F enclosePosition(const LocalCurves& m, int y, double u0, double u1) {
    Box3F b{};
    bool first = true;
    const SpanPieces pieces = spanPieces(m, y, u0, u1);
    for (int k = pieces.first; k <= pieces.last; ++k) {
        float t0, t1;
        pieceAt(pieces, k, t0, t1);
        const Bezier z = spanOf(m, y, k);
        growBox(b, first, blossom3(z, t0, t0, t0));
        growBox(b, first, blossom3(z, t0, t0, t1));
        growBox(b, first, blossom3(z, t0, t1, t1));
        growBox(b, first, blossom3(z, t1, t1, t1));
    }
    return b;
}
__device__ __forceinline__ void hodo(const Bezier& z, f3& h0, f3& h1, f3& h2) {
    h0 = (z.b[1] - z.b[0]) * 3.0f;
    h1 = (z.b[2] - z.b[1]) * 3.0f;
    h2 = (z.b[3] - z.b[2]) * 3.0f;
}
__device__ Box3F encloseDerivative(const LocalCurves& m, int y, double u0, double u1) {
    Box3F b{};
    bool first = true;
    const SpanPieces pieces = spanPieces(m, y, u0, u1);
    for (int k = pieces.first; k <= pieces.last; ++k) {
        float t0, t1;
        pieceAt(pieces, k, t0, t1);
        f3 h0, h1, h2;
        hodo(spanOf(m, y, k), h0, h1, h2);
        growBox(b, first, blossom2(h0, h1, h2, t0, t0));
        growBox(b, first, blossom2(h0, h1, h2, t0, t1));
        growBox(b, first, blossom2(h0, h1, h2, t1, t1));
    }
    return b;
}

struct Box2 {
    FloatInterval x, y;
};
__device__ __forceinline__ Box2 projBox(const LocalCurves& m, int y, double a0, double a1) {
    const Box3F b = enclosePosition(m, y, a0, a1);
    return {ipad(b.x, m.etaX), ipad(b.y, m.etaX)};
}
__device__ __forceinline__ bool overlaps2(const Box2& A, const Box2& B) {
    return A.x.lo <= B.x.hi && B.x.lo <= A.x.hi && A.y.lo <= B.y.hi && B.y.lo <= A.y.hi;
}

// Angular hull of projected vectors (crossings.cpp angularHull).
constexpr int kHullMax = 24;
__device__ bool angularHull(const f3* v, int n, float minNorm, float& lo, float& hi) {
    float ang[kHullMax];
    for (int i = 0; i < n; ++i) {
        if (hypotf(v[i].x, v[i].y) <= minNorm) return false;
        ang[i] = atan2f(v[i].y, v[i].x);
    }
    for (int a = 1; a < n; ++a) {  // insertion sort
        const float x = ang[a];
        int b = a - 1;
        while (b >= 0 && ang[b] > x) { ang[b + 1] = ang[b]; --b; }
        ang[b + 1] = x;
    }
    float bestGap = ang[0] + 2 * kPiF - ang[n - 1];
    int bestIdx = 0;
    for (int i = 1; i < n; ++i) {
        const float g = ang[i] - ang[i - 1];
        if (g > bestGap) { bestGap = g; bestIdx = i; }
    }
    if (bestGap <= kPiF + kAngTol) return false;
    lo = ang[bestIdx] - kAngTol;
    hi = (bestIdx == 0 ? ang[n - 1] : ang[bestIdx - 1] + 2 * kPiF) + kAngTol;
    return true;
}
// Hodograph controls of [a0, a1] appended to v (false: too many pieces: not certified).
__device__ bool hodoControls(const LocalCurves& m, int y, double a0, double a1, f3* v, int& n) {
    bool ok = true;
    const SpanPieces pieces = spanPieces(m, y, a0, a1);
    for (int k = pieces.first; k <= pieces.last; ++k) {
        float t0, t1;
        pieceAt(pieces, k, t0, t1);
        if (n + 3 > kHullMax) { ok = false; continue; }
        f3 h0, h1, h2;
        hodo(spanOf(m, y, k), h0, h1, h2);
        v[n++] = blossom2(h0, h1, h2, t0, t0);
        v[n++] = blossom2(h0, h1, h2, t0, t1);
        v[n++] = blossom2(h0, h1, h2, t1, t1);
    }
    return ok;
}
__device__ bool cone(const LocalCurves& m, int y, double a0, double a1, float& lo, float& hi) {
    f3 v[kHullMax];
    int n = 0;
    if (!hodoControls(m, y, a0, a1, v, n)) return false;
    return angularHull(v, n, kMinNormRel * m.lref, lo, hi);
}
__device__ bool unionCone(const LocalCurves& m, int y, double a0, double a1, double b0, double b1) {
    f3 v[kHullMax];
    int n = 0;
    if (!hodoControls(m, y, a0, a1, v, n) || !hodoControls(m, y, b0, b1, v, n)) return false;
    float lo, hi;
    return angularHull(v, n, kMinNormRel * m.lref, lo, hi);
}
__device__ bool disjointModPi(float a0, float a1, float b0, float b1) {
    for (int k = -4; k <= 4; ++k) {
        const float s = k * kPiF;
        if (a0 <= b1 + s && b0 + s <= a1) return false;
    }
    return true;
}

// Krawczyk operator for F(a, b) = A(a) - B(b) (projected), on the box (a + dXa) x (b + dXb);
// the result as offsets from (a, b).
// The two curves of the interval test (yarn 0 or 1 of the mini set), and a point of their parameter plane.
struct CurvePair {
    int ya, yb;
};
struct PairPoint {
    double uA, uB;
};
struct KrawczykResult {
    bool ok;
    FloatInterval Ka, Kb;
};
__device__ KrawczykResult krawczyk(const LocalCurves& m, CurvePair pair, PairPoint at, FloatInterval dXa, FloatInterval dXb) {
    KrawczykResult out{false, {0, 0}, {0, 0}};
    const d3 fa = evalAtD(m, pair.ya, at.uA), fb = evalAtD(m, pair.yb, at.uB);
    const double fx = fa.x - fb.x, fy = fa.y - fb.y, e = 2 * kEtaF * m.lref;
    const FloatInterval F0 = {__double2float_rd(fx - e), __double2float_ru(fx + e)}, F1 = {__double2float_rd(fy - e), __double2float_ru(fy + e)};
    const Box3F dA = encloseDerivative(m, pair.ya, at.uA + dXa.lo, at.uA + dXa.hi), dB = encloseDerivative(m, pair.yb, at.uB + dXb.lo, at.uB + dXb.hi);
    auto padD = [&](FloatInterval v) { return ipad(v, kDerRel * imag(v) + kDerAbs * m.lref); };
    const FloatInterval J00 = padD(dA.x), J10 = padD(dA.y), J01 = ineg(padD(dB.x)), J11 = ineg(padD(dB.y));
    const f3 pa = d1At(m, pair.ya, at.uA), pb = d1At(m, pair.yb, at.uB);
    const float m00 = pa.x, m01 = -pb.x, m10 = pa.y, m11 = -pb.y;
    const float det = m00 * m11 - m01 * m10;
    if (!(fabsf(det) > 0)) return out;
    const FloatInterval Y00 = fi(m11 / det), Y01 = fi(-m01 / det), Y10 = fi(-m10 / det), Y11 = fi(m00 / det);
    const FloatInterval YF0 = iadd(imul(Y00, F0), imul(Y01, F1)), YF1 = iadd(imul(Y10, F0), imul(Y11, F1));
    const FloatInterval M00 = isub(fi(1.0f), iadd(imul(Y00, J00), imul(Y01, J10)));
    const FloatInterval M01 = ineg(iadd(imul(Y00, J01), imul(Y01, J11)));
    const FloatInterval M10 = ineg(iadd(imul(Y10, J00), imul(Y11, J10)));
    const FloatInterval M11 = isub(fi(1.0f), iadd(imul(Y10, J01), imul(Y11, J11)));
    out.Ka = iadd(ineg(YF0), iadd(imul(M00, dXa), imul(M01, dXb)));
    out.Kb = iadd(ineg(YF1), iadd(imul(M10, dXa), imul(M11, dXb)));
    out.ok = strictlyInside(out.Ka, dXa) && strictlyInside(out.Kb, dXb);
    return out;
}

// crossings.cpp newton: float iterations first (cheap), then the CPU's tolerance with F in
// double (at most 50 iterations in all, as the CPU).
struct NewtonResult {
    bool converged;
    PairPoint at;  // where the iteration ended (inside the boxes)
};
__device__ NewtonResult newton(const LocalCurves& m, CurvePair pair, PairPoint start, DoubleInterval Ia, DoubleInterval Ib) {
    double a = start.uA, b = start.uB;
    int it = 0;
    for (; it < 50; ++it) {
        const f3 F = evalAt(m, pair.ya, a) - evalAt(m, pair.yb, b);
        if (hypotf(F.x, F.y) <= 4e-6f * m.lref) break;  // float resolution reached: refine in double
        const f3 da = d1At(m, pair.ya, a), db = d1At(m, pair.yb, b);
        const float j00 = da.x, j01 = -db.x, j10 = da.y, j11 = -db.y;
        const float det = j00 * j11 - j01 * j10;
        if (!(fabsf(det) > 1e-30f)) return {false, {a, b}};
        a = clampd(a - (j11 * F.x - j01 * F.y) / det, Ia.lo, Ia.hi);
        b = clampd(b - (-j10 * F.x + j00 * F.y) / det, Ib.lo, Ib.hi);
    }
    for (; it < 50; ++it) {
        const d3 A = evalAtD(m, pair.ya, a), B = evalAtD(m, pair.yb, b);
        const double Fx = A.x - B.x, Fy = A.y - B.y;
        if (hypot(Fx, Fy) <= kNewtonIn * m.lref) return {true, {a, b}};
        const f3 da = d1At(m, pair.ya, a), db = d1At(m, pair.yb, b);
        const float j00 = da.x, j01 = -db.x, j10 = da.y, j11 = -db.y;
        const float det = j00 * j11 - j01 * j10;
        if (!(fabsf(det) > 1e-30f)) return {false, {a, b}};
        const double sa = -(j11 * Fx - j01 * Fy) / det;
        const double sb = -(-j10 * Fx + j00 * Fy) / det;
        a = clampd(a + sa, Ia.lo, Ia.hi);
        b = clampd(b + sb, Ib.lo, Ib.hi);
    }
    const d3 A = evalAtD(m, pair.ya, a), B = evalAtD(m, pair.yb, b);
    return {hypot(A.x - B.x, A.y - B.y) <= kNewtonOut * m.lref, {a, b}};
}

__device__ FloatInterval depthOver(const LocalCurves& m, int ya, int yb, DoubleInterval Xa, DoubleInterval Xb) {
    const Box3F A = enclosePosition(m, ya, Xa.lo, Xa.hi), B = enclosePosition(m, yb, Xb.lo, Xb.hi);
    return ipad(isub(A.z, B.z), 2 * m.etaX);
}

// YarnPath::stitchContaining on the mini yarn -> the mini topology id (-1: none).
__device__ int stitchContaining(const LocalCurves& m, int y, double u0, double u1) {
    if (y == 0) {
        const double uE = m.nSpans[0];
        return u0 >= 0.0 && u1 <= uE ? 0 : -1;  // one span [0, uMax], which is the last
    }
    const double uE = m.nSpans[1];
    for (int s = 0; s < m.it->nSpans; ++s) {
        const SignatureSpan& st = m.spans[m.it->spanOff + s];
        const bool last = st.uEnd >= uE;
        if (u0 >= st.uBegin && (u1 < st.uEnd || (last && u1 <= st.uEnd))) return st.id;
    }
    return -1;
}

// verifyAndRecord without the deduplication: true and the candidate if a root is certified.
__device__ bool verify(const LocalCurves& m, CurvePair pair, PairPoint at, CrossingCandidate& x) {
    const double ua = m.nSpans[pair.ya], ub = m.nSpans[pair.yb];
    const double radii[6] = {1e-7, 1e-6, 1e-5, 1e-4, 3e-4, 1e-3};
    for (int ri = 0; ri < 6; ++ri) {
        const double r = radii[ri];
        const DoubleInterval Xa = dintersect({at.uA - r, at.uA + r}, {0.0, ua}), Xb = dintersect({at.uB - r, at.uB + r}, {0.0, ub});
        if (pair.ya == pair.yb && doverlaps(Xa, Xb)) return false;  // diagonal: not at.uA crossing
        const KrawczykResult k = krawczyk(m, pair, at, offsetOf(Xa, at.uA), offsetOf(Xb, at.uB));
        if (!k.ok) continue;
        FloatInterval Ka = k.Ka, Kb = k.Kb;
        auto glob = [](FloatInterval o, double c) { return DoubleInterval{c + double(o.lo), c + double(o.hi)}; };
        FloatInterval depth = depthOver(m, pair.ya, pair.yb, glob(Ka, at.uA), glob(Kb, at.uB));
        for (int it = 0; it < 40 && containsZero(depth); ++it) {
            const KrawczykResult k2 = krawczyk(m, pair, at, Ka, Kb);
            const FloatInterval na = iintersect(k2.Ka, Ka), nb = iintersect(k2.Kb, Kb);
            if (na.hi - na.lo >= Ka.hi - Ka.lo && nb.hi - nb.lo >= Kb.hi - Kb.lo) break;
            Ka = na;
            Kb = nb;
            depth = depthOver(m, pair.ya, pair.yb, glob(Ka, at.uA), glob(Kb, at.uB));
        }
        x.ya = pair.ya;
        x.yb = pair.yb;
        x.uA = at.uA;
        x.uB = at.uB;
        x.boxA = glob(Ka, at.uA);
        x.boxB = glob(Kb, at.uB);
        x.depth = depth;
        x.status = 0;
        if (depth.lo > 0) x.front = pair.ya;
        else if (depth.hi < 0) x.front = pair.yb;
        else { x.front = -1; x.status = 2; }
        x.stA = stitchContaining(m, pair.ya, x.boxA.lo, x.boxA.hi);
        x.stB = stitchContaining(m, pair.yb, x.boxB.lo, x.boxB.hi);
        if (x.stA < 0 || x.stB < 0) x.status = 2;
        return true;
    }
    return false;
}

struct Work {
    double a0, a1, b0, b1;
    signed char ya, yb, kind, depth;  // kind 0 Self, 1 Adjacent, 2 Distinct
};

// The items of a check, as the kernels read them.
struct SignatureView {
    const SignatureItem* items;
    int count;
    const float* ctrl;
    const float* lref;
    const SignatureSpan* spans;
};
struct CrossingOutput {
    unsigned* flags;
    CrossingCandidate* candidates;
    unsigned long long* keys;
    unsigned* count;
    unsigned capacity;
};
struct SortedCandidates {
    const CrossingCandidate* candidates;
    const unsigned long long* keys;
    unsigned count;
};

__device__ __forceinline__ int itemOfPair(const SignatureItem* items, int nItems, long long p) {
    int lo = 0, hi = nItems - 1;
    while (lo < hi) {
        const int mid = (lo + hi + 1) >> 1;
        if (items[mid].pairBase <= p) lo = mid; else hi = mid - 1;
    }
    return lo;
}

__device__ LocalCurves miniOf(const SignatureView& view, int i) {
    const SignatureItem* items = view.items;
    const float* ctrl = view.ctrl;
    const float* lref = view.lref;
    const SignatureSpan* spans = view.spans;
    const SignatureItem& it = items[i];
    LocalCurves m;
    m.ctrl[0] = ctrl + 3 * it.ctrlA;
    m.ctrl[1] = ctrl + 3 * it.ctrlB;
    m.nSpans[0] = it.nCtrlA - 3;
    m.nSpans[1] = it.nCtrlB - 3;
    m.lref = lref[i];
    m.etaX = kEtaRel * m.lref;
    m.it = &it;
    m.spans = spans;
    return m;
}

// --- kernel 0
__global__ void transformKernel(const SignatureItem* items, int nItems, const double* src, float* ctrl, float* lref) {
    const int i = threadItem();
    if (i >= nItems) return;
    const SignatureItem& it = items[i];
    const double* F = it.frame;
    double ext = 0;
    for (int y = 0; y < 2; ++y) {
        const int n = y ? it.nCtrlB : it.nCtrlA, s = y ? it.srcB : it.srcA, dst = y ? it.ctrlB : it.ctrlA;
        for (int c = 0; c < n; ++c) {
            const double dx = src[3 * (s + c)] - F[0], dy = src[3 * (s + c) + 1] - F[1], dz = src[3 * (s + c) + 2] - F[2];
            const double X = dx * F[3] + dy * F[4] + dz * F[5], Y = dx * F[6] + dy * F[7] + dz * F[8], Z = dx * F[9] + dy * F[10] + dz * F[11];
            ctrl[3 * (dst + c)] = static_cast<float>(X);
            ctrl[3 * (dst + c) + 1] = static_cast<float>(Y);
            ctrl[3 * (dst + c) + 2] = static_cast<float>(Z);
            ext = fmax(ext, fmax(fabs(X), fmax(fabs(Y), fabs(Z))));
        }
    }
    lref[i] = static_cast<float>(fmax(1.0, 2.0 * ext));  // CurveSet::extent
}

// --- kernel 1
__global__ void findCrossingsKernel(SignatureView view, long long first, long long n, CrossingOutput out, const volatile int* abortFlag) {
    const SignatureItem* items = view.items;
    const long long t = threadItem64();
    if (t >= n) return;
    const long long p = first + t;
    const int i = itemOfPair(items, view.count, p);
    const LocalCurves m = miniOf(view, i);
    // Decode the top-level pair (crossings.cpp loop order: ya, yb >= ya, ka, kb).
    long long q = p - items[i].pairBase;
    const int na = m.nSpans[0], nb = m.nSpans[1];
    int ya = 0, yb = 0, ka = 0, kb = 0;
    const long long s00 = static_cast<long long>(na) * (na + 1) / 2, s01 = static_cast<long long>(na) * nb;
    if (q < s00) {
        ya = yb = 0;
        for (ka = 0; q >= na - ka; ++ka) q -= na - ka;
        kb = ka + static_cast<int>(q);
    } else if (q < s00 + s01) {
        q -= s00;
        ya = 0;
        yb = 1;
        ka = static_cast<int>(q / nb);
        kb = static_cast<int>(q % nb);
    } else {
        q -= s00 + s01;
        ya = yb = 1;
        for (ka = 0; q >= nb - ka; ++ka) q -= nb - ka;
        kb = ka + static_cast<int>(q);
    }
    Work stack[kStack];
    int top = 0;
    if (ya == yb && kb == ka) stack[top++] = {double(ka), double(ka + 1), double(ka), double(ka + 1), (signed char)ya, (signed char)ya, 0, 0};
    else if (ya == yb && kb == ka + 1)
        stack[top++] = {double(ka), double(ka + 1), double(kb), double(kb + 1), (signed char)ya, (signed char)yb, 1, 0};
    else {
        if (!overlaps2(projBox(m, ya, ka, ka + 1), projBox(m, yb, kb, kb + 1))) return;
        stack[top++] = {double(ka), double(ka + 1), double(kb), double(kb + 1), (signed char)ya, (signed char)yb, 2, 0};
    }
    int budget = kBudget, seq = 0;
    unsigned fl = 0;
    auto push = [&](const Work& w) {
        if (top < kStack) stack[top++] = w;
        else fl |= kFlagUnres | kFlagBudget;
    };
    while (top > 0) {
        const Work w = stack[--top];
        if (--budget < 0) {
            fl |= kFlagUnres | kFlagBudget;
            break;
        }
        if (abortRequested(abortFlag, static_cast<unsigned>(budget))) break;  // cancelled: the result is discarded
        const int depth = w.depth;
        if (w.kind == 0) {
            float lo, hi;
            if (cone(m, w.ya, w.a0, w.a1, lo, hi)) continue;
            if (depth >= kMaxDepth || w.a1 - w.a0 < kSelfMinWidth) { fl |= kFlagUnres; continue; }
            const double mid = 0.5 * (w.a0 + w.a1);
            const signed char d1 = static_cast<signed char>(depth + 1);
            // Pushed in reverse: processed L-L, R-R, L-R.
            push({w.a0, mid, mid, w.a1, w.ya, w.ya, 1, d1});
            push({mid, w.a1, mid, w.a1, w.ya, w.ya, 0, d1});
            push({w.a0, mid, w.a0, mid, w.ya, w.ya, 0, d1});
            continue;
        }
        if (w.kind == 1) {
            if (unionCone(m, w.ya, w.a0, w.a1, w.b0, w.b1)) continue;
            if (depth >= kMaxDepth || w.a1 - w.a0 < kSelfMinWidth) { fl |= kFlagUnres; continue; }
            const double ma = 0.5 * (w.a0 + w.a1), mb = 0.5 * (w.b0 + w.b1);
            const signed char d1 = static_cast<signed char>(depth + 1);
            // Processed: (A2, B1) Adjacent, (A1, B1), (A1, B2), (A2, B2) Distinct.
            push({ma, w.a1, mb, w.b1, w.ya, w.yb, 2, d1});
            push({w.a0, ma, mb, w.b1, w.ya, w.yb, 2, d1});
            push({w.a0, ma, w.b0, mb, w.ya, w.yb, 2, d1});
            push({ma, w.a1, w.b0, mb, w.ya, w.yb, 1, d1});
            continue;
        }
        if (!overlaps2(projBox(m, w.ya, w.a0, w.a1), projBox(m, w.yb, w.b0, w.b1))) continue;
        float la, ha, lb, hb;
        const bool transversal = cone(m, w.ya, w.a0, w.a1, la, ha) && cone(m, w.yb, w.b0, w.b1, lb, hb) && disjointModPi(la, ha, lb, hb);
        if (transversal) {
            const double qa = 0.25 * (w.a1 - w.a0), qb = 0.25 * (w.b1 - w.b0);
            const DoubleInterval Ea = dintersect({w.a0 - qa, w.a1 + qa}, {0.0, double(m.nSpans[w.ya])});
            const DoubleInterval Eb = dintersect({w.b0 - qb, w.b1 + qb}, {0.0, double(m.nSpans[w.yb])});
            const bool sameYarnOverlap = w.ya == w.yb && doverlaps(Ea, Eb);
            float l2a, h2a, l2b, h2b;
            const bool expT = !sameYarnOverlap && cone(m, w.ya, Ea.lo, Ea.hi, l2a, h2a) && cone(m, w.yb, Eb.lo, Eb.hi, l2b, h2b) &&
                              disjointModPi(l2a, h2a, l2b, h2b);
            if (expT) {
                const CurvePair pair{w.ya, w.yb};
                const NewtonResult root = newton(m, pair, {0.5 * (w.a0 + w.a1), 0.5 * (w.b0 + w.b1)}, Ea, Eb);
                if (root.converged) {
                    CrossingCandidate x;
                    if (verify(m, pair, root.at, x)) {
                        if (seq >= (1 << kSeqBits)) {
                            fl |= kFlagUnres | kFlagBudget;
                        } else {
                            const unsigned at = atomicAdd(out.count, 1u);
                            if (at < out.capacity) {
                                out.candidates[at] = x;
                                out.keys[at] = (static_cast<unsigned long long>(p) << kSeqBits) | static_cast<unsigned>(seq);
                            }
                            ++seq;
                        }
                        continue;
                    }
                }
            }
            // A transversal box whose root Newton / Krawczyk cannot certify: verify's radii do
            // not depend on the box, so subdividing below kMinNonTransversalWidth cannot help.
            // (A near-tangent crossing, angle below the float derivative pad: the recursion
            // otherwise ran to depth 24 at ~0.25 ms per step, 60 s per launch on teddy piece 5.)
            if (depth >= kMaxDepth || fmax(w.a1 - w.a0, w.b1 - w.b0) < kMinNonTransversalWidth) { fl |= kFlagUnres; continue; }
        } else if (depth >= kMaxDepth || fmax(w.a1 - w.a0, w.b1 - w.b0) < kMinNonTransversalWidth) {
            fl |= kFlagUnres;
            continue;
        }
        const double ma = 0.5 * (w.a0 + w.a1), mb = 0.5 * (w.b0 + w.b1);
        const signed char d1 = static_cast<signed char>(depth + 1);
        // Processed: (A1, B1), (A1, B2), (A2, B1), (A2, B2).
        push({ma, w.a1, mb, w.b1, w.ya, w.yb, 2, d1});
        push({ma, w.a1, w.b0, mb, w.ya, w.yb, 2, d1});
        push({w.a0, ma, mb, w.b1, w.ya, w.yb, 2, d1});
        push({w.a0, ma, w.b0, mb, w.ya, w.yb, 2, d1});
    }
    if (fl) atomicOr(&out.flags[i], fl);
}

// --- kernel 2
__device__ int combineSt(int a, int b) { return a == 1 || b == 1 ? 1 : (a == 2 || b == 2 ? 2 : 0); }

__global__ void evaluateSignatureKernel(SignatureView view, const unsigned* flags, SortedCandidates sorted, SignatureOutcome* outs) {
    const SignatureItem* items = view.items;
    const int nItems = view.count;
    const int i = threadItem();
    if (i >= nItems) return;
    const SignatureItem& it = items[i];
    const LocalCurves m = miniOf(view, i);
    SignatureOutcome o{};
    unsigned fl = flags[i];
    // This item's candidates: keys in [pairBase, next pairBase) << seq bits.
    const unsigned long long kLo = static_cast<unsigned long long>(it.pairBase) << kSeqBits;
    unsigned lo = 0, hi = sorted.count;
    while (lo < hi) {
        const unsigned mid = (lo + hi) >> 1;
        if (sorted.keys[mid] < kLo) lo = mid + 1; else hi = mid;
    }
    const unsigned long long kHi = i + 1 < nItems ? static_cast<unsigned long long>(items[i + 1].pairBase) << kSeqBits : ~0ull;
    CrossingCandidate X[kMaxCross];
    int nx = 0;
    for (unsigned c = lo; c < sorted.count && sorted.keys[c] < kHi; ++c) {
        CrossingCandidate x = sorted.candidates[c];
        bool drop = false;
        for (int k = 0; k < nx && !drop; ++k) {
            CrossingCandidate& e = X[k];
            if (e.ya != x.ya || e.yb != x.yb) continue;
            if (!doverlaps(e.boxA, x.boxA) || !doverlaps(e.boxB, x.boxB)) continue;
            DoubleInterval Ha = dhull(e.boxA, x.boxA), Hb = dhull(e.boxB, x.boxB);
            const double pa = 0.5 * (Ha.hi - Ha.lo) + kDedupPad, pb = 0.5 * (Hb.hi - Hb.lo) + kDedupPad;
            Ha = dintersect({Ha.lo - pa, Ha.hi + pa}, {0.0, double(m.nSpans[x.ya])});
            Hb = dintersect({Hb.lo - pb, Hb.hi + pb}, {0.0, double(m.nSpans[x.yb])});
            if (krawczyk(m, {x.ya, x.yb}, {e.uA, e.uB}, offsetOf(Ha, e.uA), offsetOf(Hb, e.uB)).ok) { drop = true; break; }
            e.status = 2;  // ambiguous duplicate
            x.status = 2;
        }
        if (drop) continue;
        if (nx == kMaxCross) { fl |= kFlagUnres | kFlagBudget; continue; }
        X[nx++] = x;
    }
    o.budget = (fl & kFlagBudget) ? 1 : 0;
    const bool anyUnres = (fl & kFlagUnres) != 0;
    if (it.type == 1) {
        // CLOSE pass: twice through the gathered loop, B,F (even) or F,B (odd), nothing else.
        int st = anyUnres ? 2 : 0;
        double hu[kMaxCross];
        int hf[kMaxCross], nh = 0;
        for (int k = 0; k < nx; ++k) {
            const CrossingCandidate& x = X[k];
            if (x.status != 0) { st = combineSt(st, 2); continue; }
            if (x.ya != x.yb) {
                const bool strandA = x.ya == 0;
                const double us = (strandA ? x.uA : x.uB) + it.s0;
                if (us < it.segLo || us > it.segHi) continue;
                const int other = stitchContaining(m, 1, strandA ? x.uB : x.uA, strandA ? x.uB : x.uA);
                const bool strandFront = strandA ? x.depth.lo > 0 : x.depth.hi < 0;
                if (other == it.headMini) {
                    const double u = strandA ? x.uA : x.uB;
                    int a = nh - 1;  // sorted insert
                    while (a >= 0 && hu[a] > u) { hu[a + 1] = hu[a]; hf[a + 1] = hf[a]; --a; }
                    hu[a + 1] = u;
                    hf[a + 1] = strandFront;
                    ++nh;
                } else {
                    ++o.unexpected;
                }
            } else {
                ++o.unexpected;
            }
        }
        o.crossings = nh;
        const bool even = it.even != 0;
        if (o.unexpected > 0) st = 1;
        if (nh != 2) { st = 1; o.closeCode = 1; }
        else if ((hf[0] != 0) == even || (hf[1] != 0) != even) { st = 1; o.closeCode = 2; }
        o.status = st;
        outs[i] = o;
        return;
    }
    // checkSignature on the mini set: pairs (child 0, parent) in parent order.
    int pairX[2][kMaxCross], pairN[2] = {0, 0};
    int unexp[kMaxCross], nUnexp = 0;
    for (int k = 0; k < nx; ++k) {
        const CrossingCandidate& x = X[k];
        const int ta = x.stA, tb = x.stB;
        int pi = -1;
        if (ta >= 0 && tb >= 0)
            for (int q = 0; q < it.nParents; ++q)
                if ((ta == 0 && tb == it.parentMini[q]) || (tb == 0 && ta == it.parentMini[q])) { pi = q; break; }
        if (pi < 0) {
            bool decChild = false;
            if (ta >= 0 && tb >= 0 && ta != tb && it.op == 2 && it.nParents == 2) {
                const bool taParent = ta == it.parentMini[0] || ta == it.parentMini[1];
                const bool tbParent = tb == it.parentMini[0] || tb == it.parentMini[1];
                decChild = taParent && tbParent;
            }
            if (decChild && x.status == 0) {
                const bool aIsLater = it.parentMini[1] == ta;
                const bool aFront = x.depth.lo > 0;
                if (aFront != aIsLater) ++o.wrongOrder;
                ++o.overlap;
                continue;
            }
            unexp[nUnexp++] = k;
            continue;
        }
        pairX[pi][pairN[pi]++] = k;
    }
    int st = 0;
    for (int q = 0; q < it.nParents; ++q) {
        int ps = 0, code = 0;
        const int n = pairN[q];
        o.crossings += n;
        // Expected sequence (expected_signature.h): 1 = Front, 0 = Behind.
        int exp[4] = {0, 1, 1, 0}, nExp = 4;
        if (it.face == 1 && it.op != 0) {
            ps = 2;
            code = 1;  // purl shaping: no derived expectation
        } else {
            if (it.face == 1) { exp[0] = 1; exp[1] = 0; exp[2] = 0; exp[3] = 1; }
            if (it.op == 1) {
                nExp = 3;
                if (it.incSlot == 0) { exp[0] = 0; exp[1] = 1; exp[2] = 1; }
                else { exp[0] = 1; exp[1] = 1; exp[2] = 0; }
            }
            if (n != nExp) {
                ps = 1;
                code = 2;
            } else {
                int byN[4], byO[4];
                for (int k = 0; k < n; ++k) byN[k] = byO[k] = pairX[q][k];
                auto newerA = [&](int k) { return X[k].stA == 0; };
                auto newerU = [&](int k) { return newerA(k) ? X[k].uA : X[k].uB; };
                auto olderU = [&](int k) { return newerA(k) ? X[k].uB : X[k].uA; };
                for (int a = 1; a < n; ++a)
                    for (int b = a; b > 0 && newerU(byN[b]) < newerU(byN[b - 1]); --b) { const int t = byN[b]; byN[b] = byN[b - 1]; byN[b - 1] = t; }
                for (int a = 1; a < n; ++a)
                    for (int b = a; b > 0 && olderU(byO[b]) < olderU(byO[b - 1]); --b) { const int t = byO[b]; byO[b] = byO[b - 1]; byO[b - 1] = t; }
                auto depthOf = [&](int k, int& s) {
                    const CrossingCandidate& x = X[k];
                    if (x.front < 0 || x.status != 0) { s = combineSt(s, 2); return 0; }
                    const bool aFront = x.depth.lo > 0;
                    return aFront == newerA(k) ? 1 : 0;
                };
                for (int k = 0; k < nExp; ++k) {
                    int s1 = 0, s2 = 0;
                    const int d1 = depthOf(byN[k], s1), d2 = depthOf(byO[k], s2);
                    ps = combineSt(ps, combineSt(s1, s2));
                    if (s1 == 0 && d1 != exp[k]) { ps = 1; code = 3; }
                    if (s2 == 0 && d2 != exp[k]) { ps = 1; code = 4; }
                }
                if (ps == 2 && code == 0) code = 5;
            }
        }
        o.pairStatus[q] = ps;
        o.pairCode[q] = code;
        o.pairFound[q] = n;
        st = combineSt(st, ps);
    }
    if (o.wrongOrder) st = combineSt(st, 1);
    if (anyUnres) st = combineSt(st, 2);
    for (int u = 0; u < nUnexp; ++u) {
        const CrossingCandidate& x = X[unexp[u]];
        if (x.ya != 0 && x.yb != 0) { ++o.neighbourhoodSelf; continue; }
        ++o.unexpected;
        st = combineSt(st, x.status == 0 ? 1 : 2);
    }
    o.status = st;
    outs[i] = o;
}

__global__ void gatherKernel(const CrossingCandidate* in, const int* order, unsigned n, CrossingCandidate* out) {
    const unsigned q = threadItem();
    if (q < n) out[q] = in[order[q]];
}

Status toStatus(int s) { return s == 1 ? Status::Fail : s == 2 ? Status::Unresolved : Status::Pass; }

// One item per child stitch, then one per CLOSE segment: each is a small two-curve problem in the
// stitch's own frame.
class SignatureItems {
public:
    SignatureItems(const CurveSet& cs, const KnitTopology& topology);
    std::vector<double> src;  // the controls of every yarn, in double
    std::vector<SignatureItem> items;
    std::vector<SignatureSpan> spans;
    std::vector<int> itemStitch;  // pair: child id; CLOSE: head id
    size_t pairItems = 0;
    int ctrlTotal = 0;
    long long pairTotal = 0;

private:
    struct Neighbourhood {
        int lo, hi;
        double u0, u1;
    };
    Neighbourhood neighbourhoodOf(int lo, int hi, const StitchRecord& first, const StitchRecord& last) const;
    void addChild(const StitchRecord& child);
    void addCloseSegments(const TerminalRecord& terminal);
    void addItem(SignatureItem& it, int yarnA, double uA0, double uA1, int yarnB, const Neighbourhood& n);

    const CurveSet& cs_;
    const KnitTopology& topology_;
    std::vector<int> srcOffset_;
};

void copyFrame(const StitchRecord& stitch, double* frame) {
    const Vec3* v[4] = {&stitch.origin, &stitch.U, &stitch.V, &stitch.N};
    for (int q = 0; q < 4; ++q) {
        frame[3 * q] = v[q]->x;
        frame[3 * q + 1] = v[q]->y;
        frame[3 * q + 2] = v[q]->z;
    }
}

SignatureItems::SignatureItems(const CurveSet& cs, const KnitTopology& topology) : cs_(cs), topology_(topology) {
    for (const YarnPath& y : cs.yarns) {
        srcOffset_.push_back(static_cast<int>(src.size() / 3));
        for (const Vec3& p : y.curve.ctrl) src.insert(src.end(), {p.x, p.y, p.z});
    }
    for (const StitchRecord& child : topology.stitches)
        if (child.nParents != 0) addChild(child);
    pairItems = items.size();
    for (const TerminalRecord& terminal : topology.terminals)
        if (terminal.op == TerminalOp::Close && terminal.segments.size() == terminal.stitches.size()) addCloseSegments(terminal);
}

// One stitch on either side of [lo, hi] on the yarn of `first` / `last`.
SignatureItems::Neighbourhood SignatureItems::neighbourhoodOf(int lo, int hi, const StitchRecord& first, const StitchRecord& last) const {
    const std::vector<StitchRecord>& stitches = topology_.stitches;
    if (lo > 0 && stitches[lo - 1].yarn == first.yarn) --lo;
    if (hi + 1 < static_cast<int>(stitches.size()) && stitches[hi + 1].yarn == last.yarn) ++hi;
    Neighbourhood n{lo, hi, 1e300, -1e300};
    for (int q = lo; q <= hi; ++q) {
        n.u0 = std::min(n.u0, stitches[q].uBegin);
        n.u1 = std::max(n.u1, stitches[q].uEnd);
    }
    return n;
}

void SignatureItems::addChild(const StitchRecord& child) {
    const std::vector<StitchRecord>& stitches = topology_.stitches;
    const StitchRecord& p0 = stitches[child.parents[0]];
    // The last parent whose neighbourhood is extracted (signature.cpp lastParentInNeighbourhood):
    // parents on both ends of a cut ring (the last before the first) or on two yarns: the first
    // parent's neighbourhood only (otherwise lo > hi: no span, an invalid control offset).
    const StitchRecord& last = stitches[child.parents[child.nParents - 1]];
    const StitchRecord& p1 = last.id < p0.id || last.yarn != p0.yarn ? p0 : last;
    const Neighbourhood n = neighbourhoodOf(p0.id, p1.id, p0, p1);
    SignatureItem it{};
    it.type = 0;
    copyFrame(child, it.frame);
    it.op = child.op == StitchOp::K ? 0 : child.op == StitchOp::Inc ? 1 : 2;
    it.incSlot = child.incSlot;
    it.face = child.face == StitchFace::Purl ? 1 : 0;
    // The mini child's parents in neighbourhood order (signature.cpp: q outer, k inner).
    it.nParents = 0;
    for (int q = n.lo; q <= n.hi; ++q)
        for (int k = 0; k < child.nParents; ++k)
            if (child.parents[k] == q && it.nParents < 2) it.parentMini[it.nParents++] = 1 + (q - n.lo);
    addItem(it, child.yarn, child.uBegin, child.uEnd, p0.yarn, n);
    items.push_back(it);
    itemStitch.push_back(child.id);
}

void SignatureItems::addCloseSegments(const TerminalRecord& terminal) {
    for (size_t k = 0; k < terminal.stitches.size(); ++k) {
        const StitchRecord& head = topology_.stitches[terminal.stitches[k]];
        const Neighbourhood n = neighbourhoodOf(head.id, head.id, head, head);
        const double s0 = std::max(0.0, std::floor(terminal.segments[k].first));
        const double s1 = std::min(cs_.yarns[terminal.yarn].uEnd(), std::ceil(terminal.segments[k].second));
        SignatureItem it{};
        it.type = 1;
        copyFrame(head, it.frame);
        it.s0 = s0;
        it.segLo = terminal.segments[k].first;
        it.segHi = terminal.segments[k].second;
        it.headMini = 1 + (head.id - n.lo);
        it.even = k % 2 == 0;
        addItem(it, terminal.yarn, s0, s1, head.yarn, n);
        items.push_back(it);
        itemStitch.push_back(head.id);
    }
}

void SignatureItems::addItem(SignatureItem& it, int yarnA, double uA0, double uA1, int yarnB, const Neighbourhood& n) {
    const int ka0 = static_cast<int>(std::floor(uA0)), ka1 = static_cast<int>(std::ceil(uA1)) - 1;
    const int kb0 = static_cast<int>(std::floor(n.u0)), kb1 = static_cast<int>(std::ceil(n.u1)) - 1;
    it.nCtrlA = ka1 - ka0 + 4;
    it.nCtrlB = kb1 - kb0 + 4;
    it.srcA = srcOffset_[yarnA] + ka0;
    it.srcB = srcOffset_[yarnB] + kb0;
    it.ctrlA = ctrlTotal;
    it.ctrlB = ctrlTotal + it.nCtrlA;
    ctrlTotal += it.nCtrlA + it.nCtrlB;
    const double bOff = std::floor(n.u0);
    it.spanOff = static_cast<int>(spans.size());
    for (int q = n.lo; q <= n.hi; ++q) {
        const StitchRecord& s = topology_.stitches[q];
        spans.push_back({s.uBegin - bOff, s.uEnd - bOff, 1 + (q - n.lo)});
    }
    std::sort(spans.begin() + it.spanOff, spans.end(), [](const SignatureSpan& a, const SignatureSpan& b) { return a.uBegin < b.uBegin; });
    it.nSpans = n.hi - n.lo + 1;
    const long long na = it.nCtrlA - 3, nb = it.nCtrlB - 3;
    it.pairBase = pairTotal;
    pairTotal += na * (na + 1) / 2 + na * nb + nb * (nb + 1) / 2;
}

// What the crossing kernel writes, with its capacity.
struct CrossingBuffers {
    DeviceBuffer flags, candidates, keys, count;
};

// Batches of top-level pairs, paced below the display watchdog (runtime.h LaunchPacer). Reruns with
// a larger buffer when the candidates overflow: none is dropped.
unsigned findCrossings(const SignatureView& view, long long pairTotal, CrossingBuffers& b, LaunchPacer& pace) {
    unsigned capacity = std::max<unsigned>(1u << 16, static_cast<unsigned>(view.count) * 16u);
    for (;;) {
        b.candidates.alloc(capacity * sizeof(CrossingCandidate));
        b.keys.alloc(capacity * 8ull);
        check(cudaMemset(b.flags.data(), 0, view.count * 4), "GPU local signature");
        check(cudaMemset(b.count.data(), 0, 4), "GPU local signature");
        const CrossingOutput out{b.flags.uints(), static_cast<CrossingCandidate*>(b.candidates.data()), b.keys.ulongs(), b.count.uints(), capacity};
        for (long long first = 0, n = 0; first < pairTotal; first += n) {
            n = std::min(pace.next(), pairTotal - first);
            const auto started = Clock::now();
            findCrossingsKernel<<<blocksFor(n), kThreadsPerBlock>>>(view, first, n, out, abortFlagDevice());
            checkLaunch("local signature: crossings");
            // Waits with the abort flag: a degenerate pair can take seconds in one batch.
            waitForKernels("local signature: crossings");
            pace.measured(n, elapsedMs(started));
        }
        unsigned found = 0;
        check(cudaMemcpy(&found, b.count.data(), 4, cudaMemcpyDeviceToHost), "GPU local signature");
        if (found <= capacity) return found;
        capacity = found + found / 2;
    }
}

std::vector<SignatureOutcome> evaluateItems(const SignatureItems& s) {
    const int nItems = static_cast<int>(s.items.size());
    std::vector<SignatureOutcome> outs(nItems);
    if (nItems == 0) return outs;
    DeviceBuffer dItems, dSpans, dSrc, dCtrl, dLref, dOuts;
    CrossingBuffers crossings;
    dItems.upload(s.items.data(), s.items.size() * sizeof(s.items[0]));
    dSpans.upload(s.spans.data(), s.spans.size() * sizeof(s.spans[0]));
    dSrc.upload(s.src);
    dCtrl.alloc(static_cast<size_t>(s.ctrlTotal) * 12);
    dLref.alloc(nItems * 4);
    crossings.flags.alloc(nItems * 4);
    crossings.count.alloc(4);
    dOuts.alloc(nItems * sizeof(SignatureOutcome));
    transformKernel<<<blocksFor(nItems), kThreadsPerBlock>>>(static_cast<SignatureItem*>(dItems.data()), nItems, dSrc.doubles(), dCtrl.floats(), dLref.floats());
    checkLaunch("local signature: transform");
    const SignatureView view{static_cast<SignatureItem*>(dItems.data()), nItems, dCtrl.floats(), dLref.floats(), static_cast<SignatureSpan*>(dSpans.data())};
    LaunchPacer pace(1ll << 13, 256, 1ll << 21);
    const unsigned nCand = findCrossings(view, s.pairTotal, crossings, pace);
    // Candidates in key order.
    DeviceBuffer dCand2, dKeys2, dIdx, dIdx2;
    dCand2.alloc(std::max(1u, nCand) * sizeof(CrossingCandidate));
    dKeys2.alloc(std::max(1u, nCand) * 8ull);
    std::vector<int> ids(nCand);
    for (unsigned q = 0; q < nCand; ++q) ids[q] = static_cast<int>(q);
    dIdx.upload(ids);
    dIdx2.alloc(std::max(1u, nCand) * 4ull);
    DevicePrimitives prims;
    prims.sortPairs({crossings.keys.ulongs(), dKeys2.ulongs(), nCand}, dIdx.ints(), dIdx2.ints());
    if (nCand > 0)
        gatherKernel<<<blocksFor(nCand), kThreadsPerBlock>>>(static_cast<CrossingCandidate*>(crossings.candidates.data()), dIdx2.ints(), nCand,
                                                             static_cast<CrossingCandidate*>(dCand2.data()));
    checkLaunch("local signature: gather");
    evaluateSignatureKernel<<<blocksFor(nItems), kThreadsPerBlock>>>(view, crossings.flags.uints(), {static_cast<CrossingCandidate*>(dCand2.data()), dKeys2.ulongs(), nCand},
                                                                     static_cast<SignatureOutcome*>(dOuts.data()));
    checkLaunch("local signature: evaluation");
    check(cudaDeviceSynchronize(), "GPU local signature");
    dOuts.download(outs.data(), outs.size() * sizeof(outs[0]));
    return outs;
}

void appendDetail(LocalPairResult& result, const std::string& text) {
    result.detail += (result.detail.empty() ? "" : "; ") + text;
}

std::string expectedOrder(const StitchRecord& child) {
    if (child.face == StitchFace::Purl) return "F,B,B,F (purl)";
    if (child.op == StitchOp::Inc) return child.incSlot == 0 ? "B,F,F (INC slot 0)" : "F,F,B (INC slot 1)";
    return "B,F,F,B";
}

// The detail text of one parent of a child, as checkLocalSignature words it.
std::string parentDetail(const SignatureOutcome& o, int q, const StitchRecord& child) {
    switch (o.pairCode[q]) {
    case 1: return std::string("purl ") + toString(child.op) + ": no derived expectation";
    case 2: return "expected " + std::to_string(child.op == StitchOp::Inc ? 3 : 4) + " crossings (" + toString(child.op) + "), found " + std::to_string(o.pairFound[q]);
    case 3: return "order along newer yarn does not match " + expectedOrder(child);
    case 4: return "order along older yarn does not match " + expectedOrder(child);
    case 5: return "a crossing sign or ownership is unresolved";
    default: return "";
    }
}

LocalPairResult childResult(const SignatureOutcome& o, const SignatureItem& item, const StitchRecord& child) {
    LocalPairResult r;
    r.child = child.id;
    r.parent = child.parents[0];
    r.overlap = o.overlap;
    r.neighbourhoodSelf = o.neighbourhoodSelf;
    for (int q = 0; q < item.nParents; ++q) {
        const std::string text = parentDetail(o, q, child);
        if (!text.empty()) appendDetail(r, text);
    }
    if (o.wrongOrder) appendDetail(r, "DEC layers in the wrong order");
    return r;
}

LocalPairResult closeResult(const SignatureOutcome& o, const SignatureItem& item, const StitchRecord& head) {
    LocalPairResult r;
    r.child = head.id;
    r.parent = -1;
    r.terminal = 1;
    if (o.closeCode == 1) r.detail = "CLOSE: expected 2 crossings with the gathered loop, found " + std::to_string(o.crossings);
    else if (o.closeCode == 2) r.detail = std::string("CLOSE: order along the strand does not match ") + (item.even ? "B,F" : "F,B");
    return r;
}

} // namespace

LocalSignatureResult checkSignatureOnDevice(const UploadedCurves& d) {
    const auto start = Clock::now();
    LocalSignatureResult res;
    const KnitTopology& topology = d.host->topo;
    const SignatureItems s(*d.host, topology);
    const std::vector<SignatureOutcome> outs = evaluateItems(s);
    // The result as checkLocalSignature builds it.
    int closeChecks = 0, budgetItems = 0;
    for (size_t i = 0; i < outs.size(); ++i) {
        const SignatureOutcome& o = outs[i];
        const StitchRecord& stitch = topology.stitches[s.itemStitch[i]];
        const bool isChild = i < s.pairItems;
        LocalPairResult lp = isChild ? childResult(o, s.items[i], stitch) : closeResult(o, s.items[i], stitch);
        lp.status = toStatus(o.status);
        lp.crossings = o.crossings;
        lp.unexpected = o.unexpected;
        budgetItems += o.budget;
        if (isChild) res.neighbourhoodSelf += lp.neighbourhoodSelf;
        else ++closeChecks;
        if (lp.unexpected) appendDetail(lp, std::to_string(lp.unexpected) + " unexpected crossings");
        if (o.budget) appendDetail(lp, "GPU work budget exhausted");
        res.pass += lp.status == Status::Pass;
        res.fail += lp.status == Status::Fail;
        res.unresolved += lp.status == Status::Unresolved;
        res.status = combine(res.status, lp.status);
        res.pairs.push_back(std::move(lp));
    }
    res.ms = elapsedMs(start);
    char buf[240];
    std::snprintf(buf, sizeof buf, "%d/%zu stitches match in their own frame (%d FAIL, %d UNRESOLVED)%s", res.pass, res.pairs.size(), res.fail,
                  res.unresolved, closeChecks ? (", incl. " + std::to_string(closeChecks) + " CLOSE loops").c_str() : "");
    res.detail = buf;
    if (budgetItems) res.detail += "; GPU work budget exhausted in " + std::to_string(budgetItems) + " items";
    return res;
}

} // namespace knit::gpu
