#pragma once
// GPU fibres: the candidate body of fibre_generator.cpp growYarnFibres, transcribed in
// source order. Kernel A: placement, population, front, hidden, ply side, exposure and
// thinning. Kernel B: draws, curve, length solve, polyline, emergence and every acceptance
// test, except the serialized root check (ExportTubeIndex::accepts), omitted by request.
#include "gpu/geometry/surface.cuh"
#include "gpu/fibres/fibres.h"

namespace knit::gpu {

constexpr float kCurlChordsPerTurnF = 16;  // fibre_shape.h kCurlChordsPerTurn
constexpr float kCurlMaxChordF = 2;        // fibre_shape.h kCurlMaxChord
constexpr float kWindRampLengthF = 0.8f;   // fibre_shape.h kWindRampLength
constexpr int kStraySideChannel = 30;      // fibre_shape.cpp kStraySideChannel

struct FibreContext {
    FibreViews d;
    FibreParamsFlat p;
    SurfaceIndexQuery index;
};

// --- islands
// fibre_islands.cpp islandAt.
__device__ __forceinline__ IslandWeight islandAt(const FibreViews& d, int yarn, double u) {
    if (!d.islands) return {};
    const CurveRecord& fc = d.cs.curves[yarn];
    const double c = clampd(u + 1, 0.0, static_cast<double>(fc.nCtrl - 1));
    const int k = static_cast<int>(c), next = min(k + 1, fc.nCtrl - 1);
    const float f = static_cast<float>(c - k);
    const IslandWeight& a = d.islands[fc.ctrlOff + k];
    const IslandWeight& b = d.islands[fc.ctrlOff + next];
    return {a.fuzzy * (1 - f) + b.fuzzy * f, a.matted * (1 - f) + b.matted * f, a.worn * (1 - f) + b.worn * f};
}

// --- draws
__device__ __forceinline__ float normalCdfF(float x) { return 0.5f * erfcf(-x / sqrtf(2.0f)); }
// fibre_shape.cpp normalQuantile (Acklam + two Newton steps) in float; p clamped to
// [FLT_MIN, 1 - 2^-24] (CPU: [1e-300, 1 - 1e-16], not finite in float).
__device__ inline float normalQuantileF(float p) {
    const float a[] = {-3.969683028665376e+01f, 2.209460984245205e+02f, -2.759285104469687e+02f, 1.383577518672690e+02f,
                       -3.066479806614716e+01f, 2.506628277459239e+00f};
    const float b[] = {-5.447609879822406e+01f, 1.615858368580409e+02f, -1.556989798598866e+02f, 6.680131188771972e+01f,
                       -1.328068155288572e+01f};
    const float c[] = {-7.784894002430293e-03f, -3.223964580411365e-01f, -2.400758277161838e+00f, -2.549732539343734e+00f,
                       4.374664141464968e+00f, 2.938163982698783e+00f};
    const float d[] = {7.784695709041462e-03f, 3.224671290700398e-01f, 2.445134137142996e+00f, 3.754408661907416e+00f};
    p = clampf(p, 1.17549435e-38f, 1.0f - 5.9604645e-8f);
    float x;
    if (p < 0.02425f) {
        const float q = sqrtf(-2 * logf(p));
        x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
    } else if (p > 1 - 0.02425f) {
        const float q = sqrtf(-2 * log1pf(-p));
        x = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
    } else {
        const float q = p - 0.5f, r = q * q;
        x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
            (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
    }
    for (int it = 0; it < 2; ++it) {
        const float e = normalCdfF(x) - p, pdf = expf(-0.5f * x * x) / sqrtf(2 * kPiF);
        if (pdf > 0) x -= e / pdf;
    }
    return x;
}
__device__ __forceinline__ float clampedNormalF(uint64_t seed, int yarn, int i, int ch0, int ch1, float c) {
    const float u1 = static_cast<float>(1.0 - hash01(seed, yarn, i, ch0)), u2 = static_cast<float>(hash01(seed, yarn, i, ch1));
    return clampf(sqrtf(-2 * logf(u1)) * cosf(2 * kPiF * u2), -c, c);
}

struct LenParams {
    float lengthMean, lengthMin, lengthMax, lengthSpread, length, thicknessSigma;
};
__device__ inline float targetLengthF(const LenParams& lp, uint64_t seed, int yarn, int i) {
    if (lp.lengthMean <= 0) return lp.length;
    const float u = static_cast<float>(hash01(seed, yarn, i, 6));
    if (lp.lengthSpread > 0) {
        const float mu = logf(lp.lengthMean), sg = lp.lengthSpread;
        const float Fa = normalCdfF((logf(lp.lengthMin) - mu) / sg), Fb = normalCdfF((logf(lp.lengthMax) - mu) / sg);
        const float L = expf(mu + sg * normalQuantileF(Fa + u * (Fb - Fa)));
        return clampf(L, lp.lengthMin, lp.lengthMax);
    }
    return clampf(-lp.lengthMean * log1pf(-u), lp.lengthMin, lp.lengthMax);
}

struct Draw {
    int sgn;
    float psi, lift, bowTheta, bowRadial, length, radiusFactor;
    float crimpPhase, crimpPlane;
    double key;
    bool random;
    float azimuth, elev, bendAngle, bendAmp;
    float curlRadius, curlTurns, curlPhase, curlDrift;
    int curlSense;
};
// fibre_shape.cpp drawPopulation.
__device__ __forceinline__ bool drawsRandom(const FibreParamsFlat& p, int yarn, int i, const IslandWeight& island) {
    const double share = static_cast<double>(island.fuzzy) * kFuzzyStrays + island.matted;
    return hash01(p.seed, yarn, i, 17) < p.randomFractionD + (1 - p.randomFractionD) * share;
}

// fibre_shape.cpp drawFibre: the same hash channels; the transforms in float, the key exact.
__device__ inline Draw drawFibre(const FibreParamsFlat& p, int yarn, int i, const IslandWeight& island) {
    const uint64_t seed = p.seed;
    Draw d;
    d.sgn = hash01(seed, yarn, i, 7) < 0.5 ? 1 : -1;
    d.psi = p.handedness * (p.helixDeg + p.helixSigmaDeg * clampedNormalF(seed, yarn, i, 8, 10, 2.5f)) * kDegF;
    d.lift = (p.liftMinDeg + (p.liftMaxDeg - p.liftMinDeg) * static_cast<float>(hash01(seed, yarn, i, 11))) * kDegF;
    d.bowTheta = (hash01(seed, yarn, i, 12) < 0.5 ? -1.0f : 1.0f) * (0.5f + 0.5f * static_cast<float>(hash01(seed, yarn, i, 13))) * p.bowTheta;
    d.bowRadial = (2 * static_cast<float>(hash01(seed, yarn, i, 14)) - 1) * p.bowRadial;
    d.random = drawsRandom(p, yarn, i, island);
    d.crimpPhase = d.crimpPlane = 0;
    if (!d.random && p.crimpAmplitude > 0) {
        d.crimpPhase = 2 * kPiF * static_cast<float>(hash01(seed, yarn, i, 28));
        d.crimpPlane = (2 * static_cast<float>(hash01(seed, yarn, i, 29)) - 1) * 0.25f * kPiF;
    }
    // populationParams: the Random population's overrides (fp.length is shared).
    LenParams lp{p.lengthMean, p.lengthMin, p.lengthMax, p.lengthSpread, p.length, p.thicknessSigma};
    if (d.random && p.hasRandomOverrides)
        lp = {p.rLengthMean, p.rLengthMin, p.rLengthMax, p.rLengthSpread, p.length, p.rThicknessSigma};
    d.length = targetLengthF(lp, seed, yarn, i) * (1 + kFuzzyLength * island.fuzzy) * (1 - kMattedLength * island.matted);
    d.radiusFactor = expf(lp.thicknessSigma * clampedNormalF(seed, yarn, i, 15, 16, 2.0f));
    d.key = hash01(seed, yarn, i, 9);
    // fibre_shape.cpp drawFibre: every direction at spread 180, else within +-spread along the yarn.
    if (p.randomSpreadDeg >= 180) {
        d.azimuth = 2 * kPiF * static_cast<float>(hash01(seed, yarn, i, 18));
    } else {
        const float side = hash01(seed, yarn, i, kStraySideChannel) < 0.5 ? 0.0f : kPiF;
        d.azimuth = side + p.randomSpreadDeg * kDegF * (2 * static_cast<float>(hash01(seed, yarn, i, 18)) - 1);
    }
    d.elev = (p.randomElevMinDeg + (p.randomElevMaxDeg - p.randomElevMinDeg) * static_cast<float>(hash01(seed, yarn, i, 19))) * kDegF * (1 - kMattedFlat * island.matted);
    d.bendAngle = 2 * kPiF * static_cast<float>(hash01(seed, yarn, i, 20));
    d.bendAmp = (hash01(seed, yarn, i, 21) < 0.5 ? -1.0f : 1.0f) * (0.5f + 0.5f * static_cast<float>(hash01(seed, yarn, i, 22))) * p.randomBend;
    d.curlRadius = d.curlTurns = d.curlPhase = d.curlDrift = 0;
    d.curlSense = 1;
    if (d.random && p.rCurlRadius > 0) {
        const float turnLength = p.rCurlLength * (0.7f + 0.8f * static_cast<float>(hash01(seed, yarn, i, 23)));
        d.curlTurns = fminf_(d.length / turnLength, static_cast<float>(p.maxSegments) / kCurlChordsPerTurnF);
        d.curlRadius = fminf_(fminf_(p.rCurlRadius * (0.5f + 0.5f * static_cast<float>(hash01(seed, yarn, i, 24))),
                                     0.9f * d.length / (2 * kPiF * d.curlTurns)),
                              0.15f * d.length);
        d.curlPhase = 2 * kPiF * static_cast<float>(hash01(seed, yarn, i, 25));
        d.curlDrift = 2 * kPiF * static_cast<float>(hash01(seed, yarn, i, 26));
        d.curlSense = hash01(seed, yarn, i, 27) < 0.5 ? -1 : 1;
    }
    return d;
}

// --- curves
// fibre_shape.h RandomCurve, the curled stray included.
struct RCurve {
    f3 A, d0, w, e1, e2;
    float a, S, rho, turns, phase, drift, rampEnd;
    int sense;
    __device__ f3 X(float u) const {
        const f3 axis = A + (d0 * u + w * (a * u * u)) * S;
        if (rho <= 0) return axis;
        const float t = 2 * kPiF * turns * u, r = fminf_(1.0f, u / rampEnd);
        const float radius = rho * r * r * (3 - 2 * r) * (1 + 0.3f * cosf(0.37f * t + drift));
        const float th = phase + sense * t + 0.7f * sinf(0.6f * t + drift);
        return axis + (e1 * cosf(th) + e2 * sinf(th)) * radius;
    }
};
__device__ inline RCurve randomCurve(const Draw& d, f3 A, f3 n, f3 T) {
    f3 t1 = T - n * dot(T, n);
    t1 = norm(t1) > kTinyLen ? normalized(t1) : normalized(cross(n, fabsf(n.x) < 0.9f ? f3{1, 0, 0} : f3{0, 1, 0}));
    const f3 t2 = cross(n, t1);
    RCurve c;
    c.A = A;
    float sa, ca, se, ce, sb, cb;
    sincosf(d.azimuth, &sa, &ca);
    sincosf(d.elev, &se, &ce);
    c.d0 = normalized((t1 * ca + t2 * sa) * ce + n * se);
    f3 e1 = n - c.d0 * dot(n, c.d0);
    e1 = norm(e1) > 1e-9f ? normalized(e1) : t1;
    const f3 e2 = cross(c.d0, e1);
    sincosf(d.bendAngle, &sb, &cb);
    c.w = e1 * cb + e2 * sb;
    c.a = d.bendAmp;
    c.S = 0;
    c.e1 = e1;
    c.e2 = e2;
    c.rho = d.curlRadius;
    c.turns = d.curlTurns;
    c.phase = d.curlPhase;
    c.drift = d.curlDrift;
    c.sense = d.curlSense;
    c.rampEnd = c.rho > 0 ? clampf(0.5f / c.turns, 0.05f, 0.5f) : 1.0f;
    return c;
}

// A point of the directed curve in its strand's coordinates (arc length, angle, radial level) and their
// derivatives along u.
struct CurveCoordinates {
    double s;
    float th, lam, ds, dth, dlam;
};
struct FCurve {
    YarnSurfaceView ys;
    double sE;
    float thetaE;
    int sgn;
    float S, p, q, bowTheta, bowLambda, r, lambda0;
    float crimpAmplitude, crimpWavelength, crimpPhase, crimpPlane;  // fibre_shape.h DirectedCurve
    float hug;  // fibre_shape.h DirectedCurve::hug; 0: the linear rise
    bool rampWinding;  // DirectedCurve::rampWinding: no winding at the anchor
    __device__ CurveCoordinates coords(float u) const {
        CurveCoordinates k;
        if (u < 0) {
            k.ds = sgn * S;
            k.dth = (rampWinding && S > 0 ? 0.0f : sgn * (S * p / r)) + bowTheta * kPiF;
            k.dlam = S * q / r + bowLambda * kPiF;
            k.s = sE + static_cast<double>(k.ds) * u;
            k.th = thetaE + k.dth * u;
            k.lam = lambda0 + k.dlam * u;
            return k;
        }
        float sn, cs;
        sincosf(kPiF * u, &sn, &cs);
        k.s = sE + static_cast<double>(sgn * S) * u;
        k.th = thetaE + sgn * (S * p / r) * u + bowTheta * sn;
        k.lam = lambda0 + (S * q / r) * u + bowLambda * sn;
        k.ds = sgn * S;
        k.dth = sgn * (S * p / r) + bowTheta * kPiF * cs;
        k.dlam = S * q / r + bowLambda * kPiF * cs;
        if (hug > 0) {
            const float level = expf(-(S * q / r) / hug * u);
            k.lam = lambda0 + hug * (1 - level) + bowLambda * sn;
            k.dlam = (S * q / r) * level + bowLambda * kPiF * cs;
        }
        if (rampWinding && S > 0) {
            const float a = S / kWindRampLengthF, e = expf(-a * u);
            const float wound = u - 2 * (1 - e) / a + (1 - e * e) / (2 * a);
            k.th = thetaE + sgn * (S * p / r) * wound + bowTheta * sn;
            k.dth = sgn * (S * p / r) * (1 - e) * (1 - e) + bowTheta * kPiF * cs;
        }
        if (crimpAmplitude > 0 && S > 0) {
            const float cycles = S / crimpWavelength, ramp = clampf(0.5f / cycles, 0.05f, 0.5f);
            const float t = fminf_(1.0f, u / ramp), f = t * t * (3 - 2 * t), dF = u < ramp ? 6 * t * (1 - t) / ramp : 0.0f;
            const float w = 2 * kPiF * cycles;
            float wsn, wcs, psn, pcs;
            sincosf(w * u + crimpPhase, &wsn, &wcs);
            sincosf(crimpPlane, &psn, &pcs);
            const float offset = crimpAmplitude * f * wsn / r, dOffset = crimpAmplitude * (dF * wsn + f * w * wcs) / r;
            k.th += offset * pcs;
            k.lam += offset * psn;
            k.dth += dOffset * pcs;
            k.dlam += dOffset * psn;
        }
        return k;
    }
    __device__ YarnCoordinates point(double s, float th, float lam) const { return yarnCoordinates(ys, s, th, lam); }
    __device__ f3 X(float u) const {
        const CurveCoordinates k = coords(u);
        return point(k.s, k.th, k.lam).X;
    }
    __device__ f3 dX(float u) const {
        const CurveCoordinates k = coords(u);
        const YarnCoordinates y = point(k.s, k.th, k.lam);
        return y.dS * k.ds + y.dTheta * k.dth + y.dLambda * k.dlam;
    }
};

// --- frame tests
__device__ __forceinline__ bool hiddenTest(const FibreContext& c, f3 p) {
    const f3 q = p * c.p.hidScale + fromF(c.p.hidOffset);
    for (int m = 0; m < c.p.nHidden; ++m)
        if (meshInside(c.d.hiddenMeshes[m], q)) return true;
    return false;
}
__device__ __forceinline__ f3 outwardAt(const FibreContext& c, f3 p) {
    return c.p.surf ? normalized(meshOutwardNormal(c.d.outward, p)) : f3{0, 0, 1};
}
// fibre_generator.cpp FibreRun::frontOfBack: how far p lies in front of the fabric's back. Asked along a fibre:
// `tri` is the previous sample's closest triangle of the outward mesh (-1: none yet), where the next search starts
// (meshHeightFrom: the same height as a fresh closest-point search, without scanning a cube twice the sample's height wide).
__device__ __forceinline__ float frontOfBack(const FibreContext& c, f3 p, int& tri) {
    return c.p.surf ? meshHeightFrom(c.d.outward, p, tri) + c.p.backDepth : p.z - c.p.zBottom;
}
// The fabric-back test along one fibre on a surface scene. A sample far in front of the back needs no search: its
// height changes by at most about kFabricBackSlope per mm it moves (a signed distance changes by 1; the interpolated
// normals of a curved mesh add a little), so a sample whose distance from the last searched one is below that one's
// clearance / kFabricBackSlope is taken as in front. Long strays rise far off the fabric, and every search there scans
// a ball about their height wide: skipping them is what makes long strays fast (GPU only; the CPU twin searches every
// sample, so the two may differ where the slope bound fails: a fold, a sharp edge).
constexpr float kFabricBackSlope = 2.0f;
struct FabricBackWalk {
    int tri = -1;           // the last searched sample's closest triangle: the next search starts there
    bool searched = false;  // a sample was searched
    f3 at{0, 0, 0};         // that sample
    float clearance = 0;    // its frontOfBack (>= 0: it was in front)
};
// True when p lies behind the fabric's back (frontOfBack(p) < 0), searched or bounded as above.
__device__ __forceinline__ bool behindFabric(const FibreContext& c, FabricBackWalk& w, f3 p) {
    if (!c.p.surf) return p.z - c.p.zBottom < 0;
    if (w.searched && w.clearance - kFabricBackSlope * norm(p - w.at) > 0) return false;
    const float front = frontOfBack(c, p, w.tri);
    w.searched = true;
    w.at = p;
    w.clearance = front;
    return front < 0;
}
// fibre_generator.cpp cutBehindFabric: the vertex count to keep (P.n: all; 2: none visible).
__device__ inline int keptInFront(const FibreContext& c, const f3* v, int n) {
    FabricBackWalk w;
    if (behindFabric(c, w, v[1])) return 2;
    for (int k = 2; k < n; ++k)
        if (behindFabric(c, w, lerp(v[k - 1], v[k], 0.5f)) || behindFabric(c, w, v[k])) return k;
    return n;
}
// The export's rounding in float: library -> scene (float) -> library.
__device__ __forceinline__ f3 quant(const FibreContext& c, f3 p) {
    if (!c.p.quantize) return p;
    const f3 cen = fromF(c.p.qCentre);
    const f3 q = p * c.p.qInvSigma + cen;
    return (q - cen) * c.p.qSigma;
}

// The exposure test (0 exposed, 1 origin inside, 2 occluded, 3 occluded by the
// own yarn), with the own-yarn march pruned to the own capsules' boxes along the ray
// (conservative: boxes widened by 1e-3 mm; more than 32 boxes: no pruning).
__device__ inline int exposure(const FibreContext& c, const DeviceSurfacePoint& sp, int yarn, double sAt, f3 N, float ext) {
    const float delta = 0.02f;
    const f3 p0 = sp.p + sp.n * delta;
    const Exclude ownLocal{yarn, sAt - c.p.ownWindow, sAt + c.p.ownWindow};
    if (pointInEnclosure(c.index, p0, 0.0f, ownLocal) || exactInside(c.index, p0, 0.0f, ownLocal, Exclude{}).inside) return 1;
    const float rayLen = c.p.surf ? c.p.exposureLength : c.p.zTop - p0.z + 1.0f;
    if (segmentHitsEnclosure(c.index, p0, p0 + N * rayLen, 0.0f, ownLocal)) return 2;
    if (c.p.light) return 0;  // light checks: no march along the own yarn (self-occlusion)
    const float ownSpan = fminf_(rayLen, 2 * c.p.ownWindow + 2 * ext);
    const float march = fminf_(0.01f, 0.25f * c.p.r);
    const float d0 = 2 * delta;
    // Own-yarn capsule boxes overlapping the marched segment: parameter intervals.
    constexpr int kMaxIv = 32;
    float iv0[kMaxIv], iv1[kMaxIv];
    int nIv = 0;
    bool prune = true;
    {
        const f3 a = p0 + N * d0, b = p0 + N * ownSpan;
        const float m = 1e-3f;
        for (GridCursor cursor(c.index.idx.grid, vmin(a, b) - f3{m, m, m}, vmax(a, b) + f3{m, m, m}); cursor.next();) {
            const int i = cursor.item();
            if (c.index.idx.caps[i].yarn != yarn) continue;
            const f3 lo = fromF(c.index.idx.grid.boxLo + 3 * i) - f3{m, m, m}, hi = fromF(c.index.idx.grid.boxHi + 3 * i) + f3{m, m, m};
            float t0 = d0, t1 = ownSpan;
            const float P[3] = {p0.x, p0.y, p0.z}, D[3] = {N.x, N.y, N.z}, L[3] = {lo.x, lo.y, lo.z}, H[3] = {hi.x, hi.y, hi.z};
            for (int k = 0; k < 3 && t0 <= t1; ++k) {
                if (fabsf(D[k]) < 1e-6f) {
                    if (P[k] < L[k] || P[k] > H[k]) t0 = t1 + 1;
                    continue;
                }
                float u0 = (L[k] - P[k]) / D[k], u1 = (H[k] - P[k]) / D[k];
                if (u0 > u1) { const float tmp = u0; u0 = u1; u1 = tmp; }
                t0 = fmaxf_(t0, u0 - 1e-3f);
                t1 = fminf_(t1, u1 + 1e-3f);
            }
            if (t0 > t1) continue;
            if (nIv == kMaxIv) { prune = false; break; }
            iv0[nIv] = t0;
            iv1[nIv] = t1;
            ++nIv;
        }
    }
    if (prune) {
        // Sort by start (insertion sort) and merge.
        for (int a = 1; a < nIv; ++a) {
            const float x0 = iv0[a], x1 = iv1[a];
            int b = a - 1;
            while (b >= 0 && iv0[b] > x0) { iv0[b + 1] = iv0[b]; iv1[b + 1] = iv1[b]; --b; }
            iv0[b + 1] = x0;
            iv1[b + 1] = x1;
        }
        int m = 0;
        for (int k = 0; k < nIv; ++k) {
            if (m && iv0[k] <= iv1[m - 1]) iv1[m - 1] = fmaxf_(iv1[m - 1], iv1[k]);
            else { iv0[m] = iv0[k]; iv1[m] = iv1[k]; ++m; }
        }
        nIv = m;
    }
    int iv = 0;
    const int nSteps = static_cast<int>(floorf((ownSpan - d0) / march));
    for (int j = 0; j <= nSteps; ++j) {
        const float d = d0 + j * march;
        if (d > ownSpan) break;
        if (prune) {
            while (iv < nIv && iv1[iv] < d) ++iv;
            if (iv == nIv) break;
            if (d < iv0[iv]) continue;
        }
        if (exactInside(c.index, p0 + N * d, 0.0f, ownLocal, Exclude{}).inside) return 3;
    }
    return 0;
}

} // namespace knit::gpu
