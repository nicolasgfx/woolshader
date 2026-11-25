// GPU fibres: kernels, the batched runner and generateFibres (fibres.h, gpu/fibre_generator.h). The candidate body with
// the light rules is fibre_generator.cpp's, in the same order; see fibre_candidate.cuh.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <functional>
#include <memory>

#include <cub/device/device_scan.cuh>
#include <cub/device/device_select.cuh>

#include "model/mesh_query.h"
#include "fibres/fibre_generator.h"
#include "fibres/fibre_shape.h"
#include "gpu/device/runtime.h"
#include "gpu/fibres/fibre_candidate.cuh"
#include "gpu/fibres/fibres.h"
#include "gpu/geometry/device_geometry.h"
#include "gpu/device.h"
#include "gpu/fibre_generator.h"
#include "knitting/tube/plies.h"

namespace knit::gpu {

namespace {

__device__ __forceinline__ int strandOf(const FibreStrand* st, int nStrands, long long cand) {
    int lo = 0, hi = nStrands - 1;
    while (lo < hi) {  // last strand with first <= cand
        const int mid = (lo + hi + 1) >> 1;
        if (st[mid].first <= cand) lo = mid; else hi = mid - 1;
    }
    return lo;
}

// --- kernel A
struct StrandList {
    const FibreStrand* items;
    int count;
};
struct PlaceOutput {
    uint8_t* code;
    uint8_t* flags;
    FibreSeed* seeds;
    uint8_t* survive;
};

__global__ void placeCandidatesKernel(FibreContext c, StrandList strandList, long long first, long long n, PlaceOutput out) {
    const FibreStrand* strands = strandList.items;
    uint8_t* code = out.code;
    uint8_t* flags = out.flags;
    FibreSeed* seeds = out.seeds;
    uint8_t* survive = out.survive;
    const long long t = threadItem64();
    if (t >= n) return;
    const long long cand = first + t;
    const int si = strandOf(strands, strandList.count, cand);
    const FibreStrand st = strands[si];
    const int i = static_cast<int>(cand - st.first);
    const int yarn = st.yarn;
    const FibreParamsFlat& p = c.p;
    survive[t] = 0;
    // Budget order: the key (Fibre::key, fibre_candidate.cuh) decides membership before any geometry.
    const float key = static_cast<float>(hash01(p.seed, yarn, i, 9));
    if (key < p.keyLo || key >= p.keyHi) {
        flags[t] = 0;
        code[t] = kOutOfWindowCode;
        return;
    }
    const YarnSurfaceView ys = devYarnSurf(c.d.cs, c.d.sm, yarn);
    const double h1 = hash01(p.seed, yarn, i, 1);
    // Candidate placement (fibre_generator.cpp): budget mode (n per strand) or density mode.
    const double sA = p.budget ? (i + h1) * ys.length / st.n : (i + h1) / p.density;
    const float theta = 2 * kPiF * static_cast<float>(hash01(p.seed, yarn, i, 2));
    const DeviceSurfacePoint sp = evalSurface(ys, sA, theta);
    const IslandWeight island = islandAt(c.d, yarn, sp.ring.u);
    const bool random = drawsRandom(p, yarn, i, island);
    flags[t] = random ? 1 : 0;
    const f3 N = outwardAt(c, sp.p);
    if ((p.surf ? dot(sp.n, N) : sp.n.z) <= p.frontCos) { code[t] = outcomeCode(FibreOutcome::BackFacing); return; }
    if (p.hidden && hiddenTest(c, sp.p)) { code[t] = outcomeCode(FibreOutcome::HiddenAnchor); return; }
    int parentYarn = yarn;
    if (p.parents) {
        const int py = c.d.plyYarnOf[yarn];
        const double ps0 = c.d.plyOffset[yarn] + sp.ring.u * c.d.plyScale[yarn];
        parentYarn = py;
        if (!random) {
            const YarnSurfaceView pys = devYarnSurf(c.d.parentCs, c.d.parentSm, py);
            const DeviceRing R = ringAt(pys, clampd(ps0, 0.0, pys.length));
            f3 vr = sp.p - R.c;
            vr = vr - R.T * dot(vr, R.T);
            if (norm(vr) < kTinyLen || dot(sp.n, normalized(vr)) < p.plyOuterCos) { code[t] = outcomeCode(FibreOutcome::PlyInnerEarly); return; }
        }
    }
    const int ex = exposure(c, sp, yarn, sA, N, ys.ext);
    if (ex == 1) { code[t] = outcomeCode(FibreOutcome::OriginInside); return; }
    if (ex == 2) { code[t] = outcomeCode(FibreOutcome::Occluded); return; }
    if (ex == 3) { code[t] = outcomeCode(FibreOutcome::OccludedSelf); return; }
    // Artistic thinning (curvature-biased keep probability).
    const f3 d1 = curveD1(ys.curve, sp.ring.u), d2 = curveD2(ys.curve, sp.ring.u);
    const float l1 = fmaxf_(kTinyLen, norm(d1));
    const float kappa = norm(cross(d1, d2)) / (l1 * l1 * l1);
    const float factor = (1 - p.curvatureWeight) + p.curvatureWeight * fminf_(2.0f, kappa / p.kappaRef);
    if (hash01(p.seed, yarn, i, 3) >= clampf(p.keepProbability * factor * (1 - island.worn), 0.0f, 1.0f)) { code[t] = outcomeCode(FibreOutcome::Probability); return; }
    code[t] = outcomeCode(FibreOutcome::NotProcessed);
    FibreSeed sd;
    sd.cand = cand;
    sd.strand = si;
    sd.i = i;
    sd.s = sA;
    sd.theta = theta;
    sd.parentYarn = parentYarn;
    seeds[t] = sd;
    survive[t] = 1;
}

// --- kernel B: one candidate, the device twin of fibre_generator.cpp growCandidate, stage by stage in its order
constexpr int K = 32;
constexpr float kYarnEndMargin = 0.05f;  // a root or a curve this close to the strand's end is refused

struct Poly {
    f3 v[kMaxFibreVerts];
    int n;
};

// One anchor of kernel A, read once (fibre_generator.cpp Candidate).
struct DeviceCandidate {
    FibreSeed seed;
    int yarn;               // the supporting strand
    YarnSurfaceView ys;
    DeviceSurfacePoint sp;  // the anchor on the surface
    Exclude ownOnly;        // the own strand around the anchor: the seed's s, not the curve's (they differ on plies)
};

// The curve, its span limit and its samples (fibre_shape.h SmoothCurve).
struct DeviceSmoothCurve {
    RCurve random{};        // a random stray's curve
    FCurve directed;        // a directed fibre's curve
    double spanMax = 0;     // the longest span the strand allows
    float spanGuess = -1;   // the last solved span: the second length-fit pass starts there
    f3 samples[K + 1];      // the curve at u = k / K
    float sampleArc[K + 1]; // their arc length
};

// What the stages decide (fibre_generator.cpp GrowingFibre).
struct GrowingFibre {
    Draw draw;  // draw.random: the population
    bool curled = false, crimped = false, hugging = false;
    int minSegments = 0, maxSegments = 0, segments = 0;
    float minTurnDeg = 0;
    float radius = 0;           // the root radius
    float sampleStep = 0;       // the clearance tests' step
    float clearanceRadius = 0;  // the radius the full checks keep clear
    DeviceSmoothCurve curve;
    float buriedDepth = 0;
    Poly pts;                   // [0] the buried root, [1] the curve start, then the visible vertices
    bool rootOnYarn = true;
    float firstChordElevation = 90;
    float emergence = -1;       // where the root segment leaves the own strand, as a fraction of it; -1: never
    float minTurn = 3.0e38f, maxTurn = 0, deviation = 0;
    bool truncated = false;
};

// Where an accepted fibre goes: its record and its point slots in the batch.
struct FibreSlot {
    knit::Fibre& record;
    knit::Vec3f* points;
};
// What the light and the full rules write differently: the normal at the emergence (o.normal and the
// emergence angle) and the attachment window.
struct FibreAttachment {
    f3 normal;
    float attach, attachArc;
};

__device__ __forceinline__ bool insideOwn(const FibreContext& c, const Exclude& ownOnly, f3 q) {
    return exactInside(c.index, q, 0.0f, ownOnly, Exclude{}).inside;
}

__device__ __forceinline__ DeviceCandidate candidateOf(const FibreContext& c, const FibreSeed& seed, const FibreStrand* strands) {
    DeviceCandidate cand;
    cand.seed = seed;
    cand.yarn = strands[seed.strand].yarn;
    cand.ys = devYarnSurf(c.d.cs, c.d.sm, cand.yarn);
    cand.sp = evalSurface(cand.ys, seed.s, seed.theta);
    cand.ownOnly = Exclude{cand.yarn, seed.s - c.p.ownWindow, seed.s + c.p.ownWindow};
    return cand;
}

// fibre_generator.cpp chooseSize: the segment count from the drawn length, the radius and the steps.
__device__ __forceinline__ void chooseSize(const FibreContext& c, GrowingFibre& g) {
    const FibreParamsFlat& p = c.p;
    const Draw& dr = g.draw;
    const bool random = dr.random;
    g.minSegments = random ? max(2, p.minSegments) : p.minSegments;
    g.maxSegments = random ? max(g.minSegments, p.maxSegments) : p.maxSegments;
    g.curled = random && dr.curlRadius > 0;
    g.crimped = !random && p.crimpAmplitude > 0;
    g.hugging = !random && p.hugHeight > 0;
    g.minTurnDeg = g.curled || g.crimped || g.hugging ? 0.0f : random ? fmaxf_(1.0f, p.minTurnDeg) : p.minTurnDeg;
    const float rf = 0.5f * p.thickness;
    g.radius = (random && p.rThickness > 0 ? 0.5f * p.rThickness : rf) * dr.radiusFactor;
    g.sampleStep = fmaxf_(p.sampleStepFactor * g.radius, 0.005f);
    g.clearanceRadius = g.radius + p.serializedMargin;
    int nv = static_cast<int>(lroundf(g.curled ? fmaxf_(kCurlChordsPerTurnF * dr.curlTurns, dr.length / kCurlMaxChordF) : dr.length / p.segLenTarget));
    nv = nv < g.minSegments ? g.minSegments : (nv > g.maxSegments ? g.maxSegments : nv);
    g.segments = min(nv, static_cast<int>(floorf(dr.length / (4.2f * g.radius))));
}

// fibre_generator.cpp setUpCurve: the random stray's curve, or the directed curve on its strand (on plies moved
// onto the parent yarn, and along the ply helix); then the buried depth and the span limit. False: on plies, the
// anchor faces into the parent yarn at its ring (PlyInnerLate).
__device__ __forceinline__ bool setUpCurve(const FibreContext& c, const DeviceCandidate& cand, GrowingFibre& g) {
    const FibreParamsFlat& p = c.p;
    const Draw& dr = g.draw;
    const DeviceSurfacePoint& sp = cand.sp;
    const bool random = dr.random;
    if (random) g.curve.random = randomCurve(dr, sp.p, sp.n, sp.ring.T);
    FCurve& fc = g.curve.directed;
    fc.ys = cand.ys;
    fc.sE = cand.seed.s;
    fc.thetaE = cand.seed.theta;
    fc.sgn = dr.sgn;
    fc.S = 0;
    fc.p = tanf(dr.psi);
    fc.q = tanf(dr.lift);
    fc.bowTheta = dr.bowTheta;
    fc.bowLambda = 0;
    fc.crimpAmplitude = random ? 0.0f : p.crimpAmplitude;
    fc.crimpWavelength = p.crimpWavelength;
    fc.crimpPhase = dr.crimpPhase;
    fc.crimpPlane = dr.crimpPlane;
    fc.r = p.r;
    fc.lambda0 = 1;
    fc.hug = p.hugHeight;
    fc.rampWinding = false;
    if (!random && p.parents) {
        const int py = c.d.plyYarnOf[cand.yarn];
        const double ps0 = c.d.plyOffset[cand.yarn] + sp.ring.u * c.d.plyScale[cand.yarn];
        const YarnSurfaceView pys = devYarnSurf(c.d.parentCs, c.d.parentSm, py);
        double sP = clampd(ps0, 0.0, pys.length);
        DeviceRing R = ringAt(pys, sP);
        for (int it = 0; it < 4; ++it) {
            sP = clampd(sP + dot(sp.p - R.c, R.T), 0.0, pys.length);
            R = ringAt(pys, sP);
        }
        const f3 v = sp.p - R.c;
        f3 ea, eb;
        ringAxes(R, ea, eb);
        const float x = dot(v, ea) / R.a, y = dot(v, eb) / R.b;
        fc.ys = pys;
        fc.sE = sP;
        fc.thetaE = R.phi + atan2f(y, x);
        fc.lambda0 = hypotf(x, y);
        fc.r = sqrtf(R.a * R.b);
        const f3 vr = v - R.T * dot(v, R.T);
        if (norm(vr) < kTinyLen || dot(sp.n, normalized(vr)) < p.plyOuterCos) return false;
    }
    // fibre_generator.cpp setUpCurve: along the ply helix.
    if (!random && p.plyPitch > 0) {
        fc.p = tanf(p.handedness * atanf(2 * kPiF * fc.r / p.plyPitch) + dr.psi);
        fc.rampWinding = true;
    }
    g.buriedDepth = p.rootDepth + g.radius;
    g.curve.spanMax = random ? 10.0 * dr.length + 1 : dr.sgn > 0 ? fc.ys.length - kYarnEndMargin - fc.sE : fc.sE - kYarnEndMargin;
    return true;
}

__device__ __forceinline__ f3 curvePoint(const GrowingFibre& g, float u) {
    return g.draw.random ? g.curve.random.X(u) : g.curve.directed.X(u);
}

// fibre_shape.cpp SmoothCurve::setSpan: the curve at span S, sampled at every (K / chords)-th u = k / K, and the
// samples' arc length (the last one: the curve's length).
__device__ __forceinline__ float sampleCurve(GrowingFibre& g, float S, int chords) {
    DeviceSmoothCurve& sc = g.curve;
    if (g.draw.random) {
        sc.random.S = S;
    } else {
        FCurve& fc = sc.directed;
        fc.S = S;
        const float cap = 0.9f * S * fc.q / (kPiF * fc.r);
        fc.bowLambda = fc.hug > 0 ? 0.0f : copysignf(fminf_(fabsf(g.draw.bowRadial / fc.r), cap), g.draw.bowRadial);
    }
    const int step = K / chords;
    sc.sampleArc[0] = 0;
    for (int k = 0; k <= K; k += step) {
        sc.samples[k] = curvePoint(g, static_cast<float>(k) / K);
        if (k) sc.sampleArc[k] = sc.sampleArc[k - step] + norm(sc.samples[k] - sc.samples[k - step]);
    }
    return sc.sampleArc[K];
}

// fibre_shape.cpp solveLength: the span whose curve has the target length, by secants on 16 chords, then
// sampled on all K. Loop tolerance 4e-6 of the target (CPU 1e-6: a float 16-chord sum carries ~2e-6 relative
// error); the acceptance tolerance 1e-4 is unchanged. False: no span within the strand's limit.
__device__ __forceinline__ bool solveLength(GrowingFibre& g, float target) {
    const double sMax = g.curve.spanMax;
    float& Sguess = g.curve.spanGuess;
    const float sMaxF = static_cast<float>(sMax);
    if (sMax <= 0) return false;
    float S0 = fminf_(Sguess > 0 ? Sguess : target, sMaxF), l0 = sampleCurve(g, S0, 16) - target;
    if (l0 < 0 && S0 >= sMaxF) return false;
    float S1 = fminf_(sMaxF, S0 * (target / fmaxf_(1e-9f, l0 + target))), l1 = sampleCurve(g, S1, 16) - target;
    for (int it = 0; it < 12 && fabsf(l1) > 4e-6f * target; ++it) {
        float S2 = l1 != l0 ? S1 - l1 * (S1 - S0) / (l1 - l0) : S1 * target / (l1 + target);
        S2 = clampf(S2, 1e-6f, sMaxF);
        S0 = S1;
        l0 = l1;
        S1 = S2;
        l1 = sampleCurve(g, S1, 16) - target;
        if (l1 < 0 && S1 >= sMaxF) return false;
    }
    if (fabsf(l1) > 1e-4f * target) return false;
    Sguess = S1;
    sampleCurve(g, S1, K);
    return true;
}

// fibre_shape.cpp equalArcVertices: the visible vertices at equal arc length along the samples (the last one at
// the curve's end), rounded as the export rounds them.
__device__ __forceinline__ void equalArcVertices(const FibreContext& c, GrowingFibre& g) {
    const DeviceSmoothCurve& sc = g.curve;
    Poly& P = g.pts;
    P.n = 2;
    for (int k = 1; k <= g.segments; ++k) {
        const float target = sc.sampleArc[K] * k / g.segments;
        int j = 1;
        while (j < K && sc.sampleArc[j] < target) ++j;
        const float tt = (target - sc.sampleArc[j - 1]) / fmaxf_(1e-15f, sc.sampleArc[j] - sc.sampleArc[j - 1]);
        const float uk = (j - 1 + clampf(tt, 0.0f, 1.0f)) / K;
        P.v[P.n++] = quant(c, k == g.segments ? sc.samples[K] : curvePoint(g, uk));
    }
}

// fibre_generator.cpp placeRoot: the first chord's elevation, the buried root tilted against it, and whether the
// root's arc length stays on the strand (by the seed's s).
__device__ __forceinline__ void placeRoot(const FibreContext& c, const DeviceCandidate& cand, GrowingFibre& g) {
    const FibreParamsFlat& p = c.p;
    const DeviceSurfacePoint& sp = cand.sp;
    Poly& P = g.pts;
    const f3 A0 = g.curve.samples[0] + sp.n * 1e-4f;
    f3 tf = P.v[2] - A0;
    const float elev = asinf(clampf(dot(normalized(tf), sp.n), -1.0f, 1.0f)) / kDegF;
    g.firstChordElevation = elev;
    tf = tf - sp.n * dot(tf, sp.n);
    tf = norm(tf) > kTinyLen ? normalized(tf) : normalized(sp.ring.T * static_cast<float>(g.draw.sgn));
    const float tilt = clampf(90.0f - elev - p.rootTurnDeg, 10.0f, 70.0f) * kDegF;
    float sTilt, cTilt;
    sincosf(tilt, &sTilt, &cTilt);
    const float rootLen = fmaxf_(g.buriedDepth / cTilt, 4.5f * g.radius);
    const f3 R = A0 - (sp.n * cTilt + tf * sTilt) * rootLen;
    P.v[0] = quant(c, R);
    P.v[1] = quant(c, A0);
    const double sR = cand.seed.s - static_cast<double>(dot(tf, sp.ring.T) * rootLen * sTilt);
    g.rootOnYarn = sR >= kYarnEndMargin && sR <= cand.ys.length - kYarnEndMargin;
}

// fibre_generator.cpp findEmergence: where the root segment leaves the own strand (bisection); -1 when the root
// is not inside it or the curve start is not outside it.
__device__ __forceinline__ void findEmergence(const FibreContext& c, const DeviceCandidate& cand, GrowingFibre& g) {
    const Poly& P = g.pts;
    g.emergence = -1;
    if (!insideOwn(c, cand.ownOnly, P.v[0]) || insideOwn(c, cand.ownOnly, P.v[1])) return;
    float lo = 0, hi = 1;
    for (int it = 0; it < 30; ++it) {
        const float m = 0.5f * (lo + hi);
        if (insideOwn(c, cand.ownOnly, lerp(P.v[0], P.v[1], m))) lo = m; else hi = m;
    }
    g.emergence = hi;
}

// fibre_generator.cpp visiblePolylineLength: from the emergence to the tip.
__device__ __forceinline__ float visiblePolylineLength(const GrowingFibre& g) {
    const Poly& P = g.pts;
    float L = (1 - g.emergence) * norm(P.v[1] - P.v[0]);
    for (int k = 2; k < P.n; ++k) L += norm(P.v[k] - P.v[k - 1]);
    return L;
}

// fibre_generator.cpp fitPolyline: the span whose visible polyline has the drawn length, in two passes (the second
// adds what the buried part took). False: no span fits on the strand (YarnEnd). It stops early when the root
// leaves the strand or never emerges; growFibre judges both.
__device__ __forceinline__ bool fitPolyline(const FibreContext& c, const DeviceCandidate& cand, GrowingFibre& g) {
    float Ls = g.draw.length;
    for (int it = 0; it < 2; ++it) {
        if (!solveLength(g, Ls)) return false;
        equalArcVertices(c, g);
        placeRoot(c, cand, g);
        if (!g.rootOnYarn) return true;
        findEmergence(c, cand, g);
        if (g.emergence < 0) return true;
        if (it < 1) Ls += g.draw.length - visiblePolylineLength(g);
    }
    return true;
}

// fibre_generator.cpp meetsStrandRules on the rounded polyline P (g's or its draped copy): every segment at least
// 4.2 r, and at most the renderer's limit (200 r) at the thinnest exported radius. Full checks assume the steepest
// taper 1.45^-(nv+1); light checks take the taper the export applies (Fibre::vertexRadius: tau = max(taper,
// 1.45^-(nv+1))), which lets long fibres keep their segments.
__device__ __forceinline__ bool meetsStrandRules(const FibreParamsFlat& p, const GrowingFibre& g, const Poly& P) {
    bool rules = true;
    const float tauMin = powf(1.45f, -static_cast<float>(g.segments + 1));
    const float lenMax = 200 * g.radius * (p.light || g.curled ? fmaxf_(p.taper, tauMin) : tauMin);
    for (int k = 1; k < P.n; ++k) {
        const float len = norm(P.v[k] - P.v[k - 1]);
        rules = rules && len >= 4.2f * g.radius && len <= lenMax;
    }
    return rules;
}

// fibre_generator.cpp measureTurns: the smallest and largest turn at the interior vertices of P [deg], lowering
// minTurn and raising maxTurn.
__device__ __forceinline__ void measureTurns(const Poly& P, float& minTurn, float& maxTurn) {
    for (int k = 1; k + 1 < P.n; ++k) {
        const float cc = clampf(dot(normalized(P.v[k] - P.v[k - 1]), normalized(P.v[k + 1] - P.v[k])), -1.0f, 1.0f);
        const float turn = acosf(cc) / kDegF;
        minTurn = fminf_(minTurn, turn);
        maxTurn = fmaxf_(maxTurn, turn);
    }
}

// Full checks, not curled: no vertex after the root inside any strand (no CPU twin; light checks cut instead).
__device__ __forceinline__ bool anyVertexInside(const FibreContext& c, const GrowingFibre& g) {
    for (int k = 1; k < g.pts.n; ++k)
        if (exactInside(c.index, g.pts.v[k], 0.0f, Exclude{}, Exclude{}).inside) return true;
    return false;
}

// fibre_shape.cpp deviationFrom: the largest distance of the curve's samples from the polyline P.
__device__ __forceinline__ float deviationFrom(const GrowingFibre& g, const Poly& P) {
    float dev = 0;
    for (int k = 0; k <= K; ++k) {
        float dk = 3.0e38f;
        for (int j = 1; j < P.n; ++j) {
            float sa, ta;
            dk = fminf_(dk, segmentDistance(g.curve.samples[k], g.curve.samples[k], P.v[j - 1], P.v[j], sa, ta));
        }
        dev = fmaxf_(dev, dk);
    }
    return dev;
}

// --- the tails: the light rules (fibre_generator.cpp) and the full checks (no CPU twin)

// fibre_generator.cpp deepContactPush: every vertex after the anchor and every segment midpoint may lie at most
// 0.2 r inside any strand, its own included (the root segment, which starts inside the own strand by
// construction, is not tested). The push that lifts q out of the strand it sinks into deepest: along that strand's
// outward normal by the depth plus the fibre's radius; zero when q is not too deep.
__device__ __forceinline__ f3 deepContactPush(const FibreContext& c, f3 q, float fibreRadius) {
    const Hit h = exactInside(c.index, q, 0.0f, Exclude{}, Exclude{});
    if (!h.inside || h.depth <= 0.2f * c.p.r) return f3{0, 0, 0};
    return h.normal * (h.depth + fibreRadius);
}
__device__ __forceinline__ bool isZero(f3 v) { return v.x == 0 && v.y == 0 && v.z == 0; }

// fibre_generator.cpp firstDeepSegment: the first segment of P whose midpoint or end lies too deep; P.n: none.
__device__ __forceinline__ int firstDeepSegment(const FibreContext& c, const Poly& P, float fibreRadius) {
    for (int k = 2; k < P.n; ++k)
        if (!isZero(deepContactPush(c, lerp(P.v[k - 1], P.v[k], 0.5f), fibreRadius)) || !isZero(deepContactPush(c, P.v[k], fibreRadius))) return k;
    return P.n;
}

// fibre_generator.cpp lengthFromCurveStart: the length of P from the curve's start.
__device__ __forceinline__ float lengthFromCurveStart(const Poly& P) {
    float L = 0;
    for (int k = 2; k < P.n; ++k) L += norm(P.v[k] - P.v[k - 1]);
    return L;
}

// fibre_generator.cpp DrapedPolyline: a draped copy of a fibre's polyline and what its checks measure on it.
struct DrapedPolyline {
    Poly pts;
    float minTurn = 3.0e38f, maxTurn = 0, deviation = 0;
};
constexpr int kDrapePushes = 3;  // per vertex; a sample still too deep after them ends the drape

// fibre_generator.cpp drapeOverStrands, step for step: a segment whose midpoint (or else its end) lies too deep lifts
// its end vertex by that sample's push and the next two by 2/3 and 1/3 of it; the tip segment gives back any excess
// length; the copy must pass the strand rules, the turn limits (only with interior vertices, as in growFibre) and the
// deep-contact walk; its turns and deviation are measured again. False: the drape failed.
__device__ bool drapeOverStrands(const FibreContext& c, const GrowingFibre& g, DrapedPolyline& d) {
    d.pts = g.pts;
    Poly& P = d.pts;
    for (int k = 2; k < P.n; ++k) {
        for (int push = 0;; ++push) {
            f3 lift = deepContactPush(c, lerp(P.v[k - 1], P.v[k], 0.5f), g.radius);
            if (isZero(lift)) lift = deepContactPush(c, P.v[k], g.radius);
            if (isZero(lift)) break;
            if (push == kDrapePushes) return false;
            P.v[k] = P.v[k] + lift;
            if (k + 1 < P.n) P.v[k + 1] = P.v[k + 1] + lift * (2.0f / 3.0f);
            if (k + 2 < P.n) P.v[k + 2] = P.v[k + 2] + lift * (1.0f / 3.0f);
        }
    }
    const float excess = lengthFromCurveStart(P) - lengthFromCurveStart(g.pts);
    if (excess > 0) {
        const f3 tip = P.v[P.n - 1] - P.v[P.n - 2];
        const float tipLength = norm(tip) - excess;
        if (tipLength < 4.2f * g.radius) return false;
        P.v[P.n - 1] = P.v[P.n - 2] + tip * (tipLength / norm(tip));
    }
    if (!meetsStrandRules(c.p, g, P)) return false;
    measureTurns(P, d.minTurn, d.maxTurn);
    if (P.n > 2 && (d.maxTurn > c.p.maxTurnDeg || d.minTurn < g.minTurnDeg)) return false;
    if (firstDeepSegment(c, P, g.radius) < P.n) return false;
    d.deviation = deviationFrom(g, P);
    return true;
}

// fibre_generator.cpp cutAtDeepContact: a directed fibre with the drape (light checks) first tries to rest on what
// it sinks into; otherwise it is cut before the first segment that fails; shorter than its minimum: rejected.
__device__ __forceinline__ FibreOutcome cutAtDeepContact(const FibreContext& c, GrowingFibre& g) {
    Poly& P = g.pts;
    const int keep = firstDeepSegment(c, P, g.radius);
    if (keep == P.n) return FibreOutcome::Accepted;
    if (c.p.drape && !g.draw.random) {
        DrapedPolyline draped;
        if (drapeOverStrands(c, g, draped)) {
            g.pts = draped.pts;
            g.minTurn = draped.minTurn;
            g.maxTurn = draped.maxTurn;
            g.deviation = draped.deviation;
            return FibreOutcome::Accepted;
        }
    }
    if (keep - 2 < max(1, g.minSegments)) return keep >= 3 ? FibreOutcome::CollisionShort : FibreOutcome::Collision;
    P.n = keep;
    g.truncated = true;
    return FibreOutcome::Accepted;
}

// fibre_generator.cpp cutBehindFabric: the fibre is cut where it passes behind the fabric's back; shorter than
// its minimum: rejected.
__device__ __forceinline__ FibreOutcome cutBehindFabric(const FibreContext& c, GrowingFibre& g) {
    Poly& P = g.pts;
    const int front = keptInFront(c, P.v, P.n);
    if (front < P.n) {
        if (front - 2 < max(1, g.minSegments)) return FibreOutcome::BehindFabric;
        P.n = front;
        g.truncated = true;
    }
    return FibreOutcome::Accepted;
}

// fibre_generator.cpp toFibre: the accepted record and its points.
__device__ __forceinline__ void writeFibre(const FibreContext& c, const DeviceCandidate& cand, const GrowingFibre& g, const FibreAttachment& at,
                                           FibreSlot& slot) {
    const FibreParamsFlat& p = c.p;
    const Poly& P = g.pts;
    const Draw& dr = g.draw;
    const bool random = dr.random;
    const f3 E = lerp(P.v[0], P.v[1], g.emergence);
    const float eL = g.emergence * norm(P.v[1] - P.v[0]);
    knit::Fibre& o = slot.record;
    o.strand = cand.seed.strand;
    o.candidate = cand.seed.i;
    o.yarnId = cand.yarn;
    o.parentYarn = cand.seed.parentYarn;
    o.pointCount = P.n;
    o.population = random ? knit::FibrePopulation::Random : knit::FibrePopulation::Directed;
    o.truncated = g.truncated;
    o.s = cand.seed.s;
    o.key = dr.key;
    o.theta = cand.seed.theta;
    o.anchor = knit::Vec3f{E.x, E.y, E.z};
    o.normal = knit::Vec3f{at.normal.x, at.normal.y, at.normal.z};
    for (int k = 0; k < P.n; ++k) slot.points[k] = knit::Vec3f{P.v[k].x, P.v[k].y, P.v[k].z};
    o.targetLength = dr.length;
    o.radius0 = g.radius;
    o.radius1 = g.radius * p.taper;
    o.emergeL = eL;
    o.attach = at.attach;
    o.attachArc = at.attachArc;
    o.psi = (random ? dr.azimuth : dr.psi) / kDegF;
    o.lift = (random ? dr.elev : dr.lift) / kDegF;
    o.smoothLength = g.curve.sampleArc[K];
    o.emergeAngle = asinf(clampf(dot(normalized(P.v[2] - P.v[1]), at.normal), -1.0f, 1.0f)) / kDegF;
    if (random) {
        o.helixActual = acosf(clampf(fabsf(dot(g.curve.random.d0, cand.sp.ring.T)), 0.0f, 1.0f)) / kDegF;
    } else {
        const FCurve& fc = g.curve.directed;
        const f3 d0 = normalized(fc.dX(0.0f));
        const f3 T0 = normalized(fc.point(fc.sE, fc.thetaE, 0.0f).dS);
        o.helixActual = acosf(clampf(fabsf(dot(d0, T0)), 0.0f, 1.0f)) / kDegF;
    }
    o.maxDev = g.deviation;
    o.minTurn = P.n > 2 ? g.minTurn : 0;
    o.maxTurn = P.n > 2 ? g.maxTurn : 0;
}

// Light checks (FibreParams::lightChecks, the CPU's rules): no clearance sampling, approaches, attachment window
// or exposure march at the emergence point.
__device__ FibreOutcome finishLight(const FibreContext& c, const DeviceCandidate& cand, GrowingFibre& g, FibreSlot& slot) {
    const FibreOutcome contact = cutAtDeepContact(c, g);
    if (contact != FibreOutcome::Accepted) return contact;
    const FibreOutcome back = cutBehindFabric(c, g);
    if (back != FibreOutcome::Accepted) return back;
    if (c.p.hidden && hiddenTest(c, g.pts.v[g.pts.n - 1])) return FibreOutcome::HiddenTip;
    const float eL = g.emergence * norm(g.pts.v[1] - g.pts.v[0]);
    writeFibre(c, cand, g, {cand.sp.n, eL + g.radius, eL + c.p.r + g.radius}, slot);
    return FibreOutcome::Accepted;
}

// --- the full checks' state: each made once by finishFull and read by the stages after it

// The polyline's cumulative arc length and the emergence's arc length, measured once before any cut. A cut only
// lowers pts.n and never removes v0 or v1, so cum[0 .. pts.n - 1] and emergeL stay valid after it.
struct PolylineArc {
    float cum[kMaxFibreVerts];
    float emergeL;
};
// E*: where the root segment leaves the own strand, and the surface under it.
struct EmergencePoint {
    f3 at;
    double s;
    DeviceSurfacePoint surface;
};
// Where the fibre first clears its own strand, and the strand arc it attaches to.
struct AttachmentWindow {
    bool found;
    int yarn;
    double sAnchor;
    float length, arc;  // fibre arc length from E*, strand arc around sAnchor
};
struct FibreSegment {
    f3 a, b;
    float l0, l1;  // their arc length from E*
};
enum class SegmentSampling { Uniform, Lipschitz };

// The full rules' collision test of one segment (Lipschitz: every sample's ball grown by half a step, and the
// attachment end sampled exactly).
__device__ bool fibreSegmentHits(const FibreContext& c, const FibreSegment& seg, const AttachmentWindow& window, float radius, SegmentSampling sampling) {
    const FibreParamsFlat& p = c.p;
    const f3 a = seg.a, b = seg.b;
    const float L0 = seg.l0, L1 = seg.l1, attachLength = window.length;
    const Exclude ownLocal{window.yarn, window.sAnchor - p.ownWindow, window.sAnchor + p.ownWindow};
    if (segmentHitsEnclosure(c.index, a, b, radius, ownLocal)) return true;
    const float len = norm(b - a);
    const float step = fmaxf_(p.sampleStepFactor * radius, 0.005f);
    const int n = max(1, static_cast<int>(ceilf(len / step)));
    const Exclude attach{window.yarn, window.sAnchor - window.arc, window.sAnchor + window.arc};
    const Exclude none{};
    if (sampling == SegmentSampling::Lipschitz) {
        const float h = len / n, grown = radius + 0.5f * h;
        const bool extra = L1 > L0 && attachLength > L0 && attachLength < L1;
        for (int k = 0; k <= n + (extra ? 1 : 0); ++k) {
            const float f = k <= n ? static_cast<float>(k) / n : (attachLength - L0) / (L1 - L0);
            const float l = L0 + f * (L1 - L0);
            const Exclude& ex = l + 0.5f * h <= attachLength ? attach : none;
            if (exactInside(c.index, lerp(a, b, f), grown, ownLocal, ex).inside) return true;
        }
        return false;
    }
    for (int k = 0; k <= n; ++k) {
        const float f = static_cast<float>(k) / n;
        const float l = L0 + f * (L1 - L0);
        const Exclude& ex = l <= attachLength ? attach : none;
        if (exactInside(c.index, lerp(a, b, f), radius, ownLocal, ex).inside) return true;
    }
    return false;
}

// E* by Newton on the ring plane: stop at |ds| < 1e-5 mm (CPU 1e-10; below float resolution of the
// coordinates), 8 iterations; and the surface point under it.
__device__ __forceinline__ EmergencePoint emergencePoint(const FibreContext& c, const DeviceCandidate& cand, const GrowingFibre& g) {
    const YarnSurfaceView& ys = cand.ys;
    EmergencePoint e;
    e.at = lerp(g.pts.v[0], g.pts.v[1], g.emergence);
    double sStar = cand.seed.s;
    for (int it = 0; it < 8; ++it) {
        const DeviceRing R = ringAt(ys, sStar);
        const float dsN = dot(e.at - R.c, R.T);
        sStar = clampd(sStar + dsN, 0.0, ys.length);
        if (fabsf(dsN) < 1e-5f) break;
    }
    const DeviceRing RE = ringAt(ys, sStar);
    f3 ea, eb;
    ringAxes(RE, ea, eb);
    const float thStar = RE.phi + atan2f(dot(e.at - RE.c, eb) / RE.b, dot(e.at - RE.c, ea) / RE.a);
    e.s = sStar;
    e.surface = evalSurface(ys, sStar, thStar);
    return e;
}

// The anchor's tests again at E*: front, hidden, exposure. False: EmergeRejected.
__device__ __forceinline__ bool emergenceVisible(const FibreContext& c, const DeviceCandidate& cand, const EmergencePoint& e) {
    const FibreParamsFlat& p = c.p;
    const f3 NE = outwardAt(c, e.surface.p);
    return !((p.surf ? dot(e.surface.n, NE) : e.surface.n.z) <= p.frontCos || (p.hidden && hiddenTest(c, e.at)) ||
             exposure(c, e.surface, cand.yarn, e.s, NE, cand.ys.ext) != 0);
}

__device__ __forceinline__ PolylineArc polylineArc(const Poly& P, float emergence) {
    PolylineArc arc;
    arc.cum[0] = 0;
    for (int k = 1; k < P.n; ++k) arc.cum[k] = arc.cum[k - 1] + norm(P.v[k] - P.v[k - 1]);
    arc.emergeL = emergence * arc.cum[1];
    return arc;
}

// The point at arc length l along the polyline (its current vertices).
__device__ __forceinline__ f3 pointAtArc(const Poly& P, const PolylineArc& arc, float l) {
    int k = 1;
    while (k + 1 < P.n && arc.cum[k] < l) ++k;
    return lerp(P.v[k - 1], P.v[k], clampf((l - arc.cum[k - 1]) / fmaxf_(1e-15f, arc.cum[k] - arc.cum[k - 1]), 0.0f, 1.0f));
}

// The attachment window: from E*, the first point (steps of sampleStep, at most attachMax) whose ball of the
// clearance radius, grown by half a step and by the own strand's bulge between export chords, leaves the own
// strand. Not found: NoClearance.
__device__ __forceinline__ AttachmentWindow attachmentWindow(const FibreContext& c, const DeviceCandidate& cand, const GrowingFibre& g, const PolylineArc& arc) {
    const FibreParamsFlat& p = c.p;
    const YarnSurfaceView& ys = cand.ys;
    const double s = cand.seed.s;
    const float r = p.r, h = g.sampleStep, eL = arc.emergeL;
    float kOwn = 0;
    {
        const int N = ys.nRings - 1;
        const int k0 = max(0, static_cast<int>((s - p.attachMax - r) / ys.ds));
        const int k1 = min(N, static_cast<int>((s + p.attachMax + r) / ys.ds));
        for (int kk = k0; kk <= k1; kk += 2) {
            const double uk = ys.rings[kk].u;
            const f3 d1 = curveD1(ys.curve, uk), d2 = curveD2(ys.curve, uk);
            const float l1 = fmaxf_(kTinyLen, norm(d1));
            kOwn = fmaxf_(kOwn, norm(cross(d1, d2)) / (l1 * l1 * l1));
        }
    }
    const float ownBulge = p.exportChordFactor * kOwn;
    float lc = -1;
    for (int j = 0;; ++j) {
        const float l = j * h;
        if (!(eL + l <= arc.cum[g.pts.n - 1] && l <= p.attachMax)) break;
        if (!exactInside(c.index, pointAtArc(g.pts, arc, eL + l), g.clearanceRadius + 0.5f * h + ownBulge, cand.ownOnly, Exclude{}).inside) { lc = l; break; }
    }
    AttachmentWindow w;
    w.found = lc >= 0;
    w.yarn = cand.yarn;
    w.sAnchor = s;
    w.length = lc + 0.5f * h;
    w.arc = fmaxf_(w.length, eL) + r + h;
    return w;
}

// Collisions, segment by segment (arc length from E*): cut before the first segment that hits; shorter than its
// minimum: rejected.
__device__ __forceinline__ FibreOutcome cutAtCollision(const FibreContext& c, GrowingFibre& g, const PolylineArc& arc, const AttachmentWindow& w) {
    // FibreParams::lipschitzSampling; a random stray always samples so.
    const SegmentSampling sampling = c.p.lipschitz != 0 || g.draw.random ? SegmentSampling::Lipschitz : SegmentSampling::Uniform;
    Poly& P = g.pts;
    int keep = P.n;
    for (int k = 1; k < P.n; ++k)
        if (fibreSegmentHits(c, {P.v[k - 1], P.v[k], arc.cum[k - 1] - arc.emergeL, arc.cum[k] - arc.emergeL}, w, g.clearanceRadius, sampling)) {
            keep = k;
            break;
        }
    if (keep < P.n) {
        if (keep - 2 < max(1, g.minSegments)) return keep >= 3 ? FibreOutcome::CollisionShort : FibreOutcome::Collision;
        P.n = keep;
        g.truncated = true;
    }
    return FibreOutcome::Accepted;
}

// Directed fibres: no approach to another strand along the visible length (the distance to the nearest other
// axis may not fall by more than a step below its running maximum).
__device__ __forceinline__ bool approachesAnotherStrand(const FibreContext& c, const DeviceCandidate& cand, const GrowingFibre& g, const PolylineArc& arc) {
    const float r = c.p.r, h = g.sampleStep;
    float runMax = -3.0e38f;
    const float clEnd = arc.cum[g.pts.n - 1];
    for (int j = 0;; ++j) {
        const float l = arc.emergeL + j * h;
        if (!(l <= clEnd + 1e-6f)) break;
        const f3 q = pointAtArc(g.pts, arc, fminf_(l, clEnd));
        float d = 3.0e38f;
        for (float R = r + 0.3f; d == 3.0e38f && R < 50; R *= 2) d = nearestAxisDistance(c.index, q, R, cand.yarn);
        d -= r;
        if (d < runMax - h) return true;
        runMax = fmaxf_(runMax, d);
    }
    return false;
}

// Full checks (fibres.checks = full): the anchor's tests again at E*, the attachment window, sampled collisions,
// the fabric's back, approaches, the hidden tip. The serialized root check is not applied (GPU experiment).
__device__ __forceinline__ FibreOutcome finishFull(const FibreContext& c, const DeviceCandidate& cand, GrowingFibre& g, FibreSlot& slot) {
    const EmergencePoint e = emergencePoint(c, cand, g);
    if (!emergenceVisible(c, cand, e)) return FibreOutcome::EmergeRejected;
    const PolylineArc arc = polylineArc(g.pts, g.emergence);
    const AttachmentWindow window = attachmentWindow(c, cand, g, arc);
    if (!window.found) return FibreOutcome::NoClearance;
    const FibreOutcome hit = cutAtCollision(c, g, arc, window);
    if (hit != FibreOutcome::Accepted) return hit;
    const FibreOutcome back = cutBehindFabric(c, g);
    if (back != FibreOutcome::Accepted) return back;
    if (!g.draw.random && approachesAnotherStrand(c, cand, g, arc)) return FibreOutcome::Approaches;
    if (c.p.hidden && hiddenTest(c, g.pts.v[g.pts.n - 1])) return FibreOutcome::HiddenTip;
    writeFibre(c, cand, g, {e.surface.n, window.length, window.arc}, slot);
    return FibreOutcome::Accepted;
}

// One candidate, in the order of fibre_generator.cpp growCandidate. An accepted fibre's record and polyline go
// to its slot.
__device__ FibreOutcome growFibre(const FibreContext& c, const FibreSeed& seed, const FibreStrand* strands, FibreSlot& slot) {
    const FibreParamsFlat& p = c.p;
    const DeviceCandidate cand = candidateOf(c, seed, strands);
    GrowingFibre g;
    g.draw = drawFibre(p, cand.yarn, seed.i, islandAt(c.d, cand.yarn, cand.sp.ring.u));
    chooseSize(c, g);
    if (g.segments < max(1, g.minSegments)) return FibreOutcome::TooShortToBend;
    if (!setUpCurve(c, cand, g)) return FibreOutcome::PlyInnerLate;
    if (g.buriedDepth >= 0.9f * p.r) return FibreOutcome::RootDepth;
    if (!fitPolyline(c, cand, g)) return FibreOutcome::YarnEnd;
    if (!g.rootOnYarn) return FibreOutcome::RootRange;
    if (g.firstChordElevation < p.minEmergeDeg) return FibreOutcome::NoEmergenceElev;
    if (g.emergence < 0) return !insideOwn(c, cand.ownOnly, g.pts.v[0]) ? FibreOutcome::RootNotInside : FibreOutcome::NoEmergence;
    if (!meetsStrandRules(p, g, g.pts)) return FibreOutcome::RuleRejected;
    measureTurns(g.pts, g.minTurn, g.maxTurn);
    if (g.pts.n > 2 && g.maxTurn > p.maxTurnDeg) return FibreOutcome::TurnTooSharp;
    if (g.pts.n > 2 && g.minTurn < g.minTurnDeg) return FibreOutcome::NoBend;
    if (!p.light && !g.curled && anyVertexInside(c, g)) return FibreOutcome::VertexInside;  // light: the deep-contact cut replaces it
    g.deviation = deviationFrom(g, g.pts);
    // A curled stray's limit follows its curl, a crimped fibre's its crimp (inline, as in growCandidate).
    const float allowedDev = g.curled    ? fmaxf_(p.maxCurveDev, 0.06f * g.draw.curlRadius)
                             : g.crimped ? fmaxf_(p.maxCurveDev, 0.06f * p.crimpAmplitude)
                                         : p.maxCurveDev;
    if (g.deviation > allowedDev) return FibreOutcome::ApproxError;
    return p.light ? finishLight(c, cand, g, slot) : finishFull(c, cand, g, slot);
}

struct GrowOutput {
    int stride;  // point slots per fibre
    uint8_t* code;
    uint8_t* truncated;
    knit::Fibre* fibres;
    knit::Vec3f* slotPoints;
    int* pointCount;
    uint8_t* accepted;
};

__global__ void growFibresKernel(FibreContext c, const FibreStrand* strands, const FibreSeed* seeds, long long n, GrowOutput out) {
    const long long t = threadItem64();
    if (t >= n) return;
    FibreSlot slot{out.fibres[t], out.slotPoints + t * out.stride};
    const FibreOutcome r = growFibre(c, seeds[t], strands, slot);
    const bool accepted = r == FibreOutcome::Accepted;
    out.code[t] = outcomeCode(r);
    out.truncated[t] = accepted && out.fibres[t].truncated ? 1 : 0;  // the record's flag: set only on acceptance
    out.accepted[t] = accepted ? 1 : 0;
    out.pointCount[t] = accepted ? out.fibres[t].pointCount : 0;
}

struct PackInput {
    long long count;
    int stride;
    const int* pointCount;
    const int* offset;  // the running sum of the point counts
    const knit::Vec3f* slotPoints;
};

// Every accepted fibre's points from its slot to their packed place, and its record told where.
__global__ void packPointsKernel(PackInput in, knit::Fibre* fibres, knit::Vec3f* packed) {
    const long long t = threadItem64();
    if (t >= in.count || in.pointCount[t] == 0) return;
    fibres[t].firstPoint = static_cast<uint32_t>(in.offset[t]);
    for (int k = 0; k < in.pointCount[t]; ++k) packed[in.offset[t] + k] = in.slotPoints[t * in.stride + k];
}

// A stream with two blocking events, for the device time of one section. The host thread sleeps in
// cudaEventSynchronize instead of spinning a core for the whole kernel (default events busy-wait).
class TimedStream {
public:
    TimedStream() {
        check(cudaStreamCreateWithFlags(&stream_, cudaStreamNonBlocking), "GPU fibres: stream");
        check(cudaEventCreateWithFlags(&start_, cudaEventBlockingSync), "GPU fibres: event");
        check(cudaEventCreateWithFlags(&end_, cudaEventBlockingSync), "GPU fibres: event");
    }
    ~TimedStream() {  // cleanup in a destructor: never throws
        cudaEventDestroy(end_);
        cudaEventDestroy(start_);
        cudaStreamDestroy(stream_);
    }
    TimedStream(const TimedStream&) = delete;
    TimedStream& operator=(const TimedStream&) = delete;
    cudaStream_t stream() const { return stream_; }
    void begin() { check(cudaEventRecord(start_, stream_), "GPU fibres: event"); }
    // Waits for the section; its device time in ms.
    double finish() {
        check(cudaEventRecord(end_, stream_), "GPU fibres: event");
        check(cudaEventSynchronize(end_), "GPU fibres");
        float ms = 0;
        check(cudaEventElapsedTime(&ms, start_, end_), "GPU fibres: event");
        return ms;
    }

private:
    cudaStream_t stream_ = nullptr;
    cudaEvent_t start_ = nullptr, end_ = nullptr;
};

} // namespace

FibreFrameGpu frameFrom(const FibreFrame& fr) {
    FibreFrameGpu g;
    g.exposureLength = fr.exposureLength;
    g.backDepth = fr.backDepth;
    g.pieceSeed = fr.pieceSeed;
    g.plies = fr.plies;
    g.parentCurves = fr.parentCurves;
    g.parentSurface = fr.parentSurface;
    g.candidates = fr.candidates;
    g.islands = &fr.islands;
    return g;
}

std::string FibreRunReport::line() const {
    char b[640];
    std::snprintf(b, sizeof b,
                  "GPU fibres: %lld candidates -> %lld exposed -> %lld accepted (%lld without the serialized root check); "
                  "%.0f ms total: host prep %.0f, admission wait %.0f, upload %.0f (+resident %.0f), kernel A %.0f, select %.0f, "
                  "kernel B %.0f, compact %.0f, download %.0f, merge %.0f; %d + %d batches (slowest launch %.0f ms); input %.0f MB, resident %.0f MB, "
                  "reserved %.0f MB, device peak %.0f MB",
                  candidates, exposed, accepted, rootCheckSkipped, totalMs, prepMs, waitMs, uploadMs, residentMs, kernelAMs, selectMs,
                  kernelBMs, compactMs, downloadMs, mergeMs, batchesA, batchesB, slowestLaunchMs, inputMB, residentMB, reservedMB, peakMB);
    return b;
}

// One run's input on the host.
struct FibreRunInput {
    CurveArrays cs;
    SurfaceArrays sm;
    IndexArrays index;
    bool parents = false;
    CurveArrays parentCs;
    SurfaceArrays parentSm;
    std::vector<int> plyYarnOf;
    std::vector<double> plyOffset, plyScale;
    std::vector<IslandWeight> islands;     // per control, in the order of cs; empty: none
    const MeshArrays* outward = nullptr;   // required with p.surf
    std::vector<const MeshQuery*> hidden;  // resident (shared by the pieces' runs)
    FibreParamsFlat p{};
    std::vector<FibreStrand> strands;
    long long candidates = 0;
    long long batchA = 1ll << 21;          // candidates per kernel-A batch
    long long batchB = 1ll << 18;          // survivors per kernel-B batch
    std::function<void()> checkpoint;      // between batches (may throw knit::Cancelled)
    std::function<void(double)> progress;  // after each batch: covered fraction of the candidates
};

struct GpuFibreInputs {
    FibreRunInput in;
    MeshArrays outward;
};

namespace {

struct FibreRunOutput {
    std::vector<uint8_t> code;        // per candidate (outcomeCode)
    std::vector<uint8_t> flags;       // per candidate: bit0 random population, bit1 truncated
    std::vector<knit::Fibre> fibres;  // accepted, in candidate order
    std::vector<knit::Vec3f> points;  // their packed points (Fibre::firstPoint)
    long long exposed = 0;
};

struct BatchSizes {
    long long a, b;  // the buffers hold this many candidates, survivors
    int stride;      // points per fibre: root, curve start, segments
};

BatchSizes batchSizes(const FibreRunInput& in) {
    const long long a = std::max<long long>(1, std::min(in.batchA, std::max<long long>(1, in.candidates)));
    return {a, std::max<long long>(1, std::min(in.batchB, a)), 2 + std::max(in.p.maxSegments, std::max(2, in.p.minSegments))};
}

// A run's working set: its inputs plus the batch buffers.
size_t inputBytes(const FibreRunInput& in, size_t hiddenMeshes) {
    return in.cs.bytes() + in.sm.bytes() + in.index.bytes() + in.parentCs.bytes() + in.parentSm.bytes() + in.islands.size() * sizeof(IslandWeight) + in.plyYarnOf.size() * 4 +
           (in.plyOffset.size() + in.plyScale.size()) * 8 + (in.outward ? in.outward->bytes() : 0) + in.strands.size() * sizeof(FibreStrand) + hiddenMeshes * sizeof(MeshView);
}
size_t batchBytes(const BatchSizes& s) {
    return static_cast<size_t>(s.a) * (2 + 1 + 2 * sizeof(FibreSeed)) + static_cast<size_t>(s.b) * (3 + 2 * sizeof(int) + 2 * sizeof(knit::Fibre) + 2 * s.stride * sizeof(knit::Vec3f)) +
           (64u << 20);  // + cub temp
}

// The run's inputs on the device, and the context the kernels read.
struct DeviceFibreInputs {
    DeviceCurves curves, parentCurves;
    DeviceSurface surface, parentSurface;
    DeviceBuffer caps, strands, plyYarnOf, plyOffset, plyScale, islands, hiddenViews;
    DeviceGrid capsuleGrid;
    DeviceMesh outward;
    FibreContext context{};

    DeviceFibreInputs(const FibreRunInput& in, const std::vector<MeshView>& hidden) {
        curves.upload(in.cs);
        surface.upload(in.sm);
        caps.upload(in.index.caps.data(), in.index.caps.size() * sizeof(CapsuleRecord));
        capsuleGrid.upload(in.index.boxes);
        strands.upload(in.strands.data(), in.strands.size() * sizeof(FibreStrand));
        hiddenViews.upload(hidden.data(), hidden.size() * sizeof(MeshView));
        if (in.parents) {
            parentCurves.upload(in.parentCs);
            parentSurface.upload(in.parentSm);
            plyYarnOf.upload(in.plyYarnOf);
            plyOffset.upload(in.plyOffset);
            plyScale.upload(in.plyScale);
        }
        if (in.outward) outward.upload(*in.outward);
        FibreViews& d = context.d;
        context.p = in.p;
        d.cs = curves.view;
        d.sm = surface.view;
        d.index.caps = static_cast<CapsuleRecord*>(caps.data());
        d.index.grid = capsuleGrid.view;
        d.parentCs = parentCurves.view;
        d.parentSm = parentSurface.view;
        d.plyYarnOf = plyYarnOf.ints();
        d.plyOffset = plyOffset.doubles();
        d.plyScale = plyScale.doubles();
        islands.upload(in.islands.data(), in.islands.size() * sizeof(IslandWeight));
        d.islands = in.islands.empty() ? nullptr : static_cast<IslandWeight*>(islands.data());
        if (in.outward) d.outward = outward.view;
        d.hiddenMeshes = static_cast<MeshView*>(hiddenViews.data());
        context.p.nHidden = static_cast<int>(hidden.size());
        context.index.idx = d.index;
        context.index.sm = d.sm;
        check(cudaDeviceSynchronize(), "GPU fibres");
    }
};

// The two-stage generator in paced batches: kernel A places the candidates and keeps the exposed
// ones, kernel B grows those into fibres. The buffers hold one batch of each.
class FibreBatches {
public:
    FibreBatches(const FibreRunInput& in, const DeviceFibreInputs& device, const BatchSizes& sizes, FibreRunOutput& out, FibreRunReport& report);
    void run();

private:
    long long placeCandidates(long long first, long long n);
    void growBatch(long long firstSurvivor, long long m);

    const FibreRunInput& in_;
    const FibreContext& context_;
    const FibreStrand* strands_;
    const int stride_;
    FibreRunOutput& out_;
    FibreRunReport& report_;
    DeviceBuffer code_, flags_, seeds_, survivors_, survive_, survivorCount_;
    DeviceBuffer codeB_, truncatedB_, fibresB_, accepted_, slotPoints_, packed_, pointCount_, offset_, acceptedFlags_, acceptedCount_, temp_;
    size_t tempBytes_ = 0;
    TimedStream timing_;
    std::vector<FibreSeed> hostSurvivors_;
    std::vector<uint8_t> hostCodeB_, hostTruncatedB_;
    std::vector<knit::Fibre> hostAccepted_;
    // Batch sizes paced below the display watchdog (runtime.h LaunchPacer); the buffers hold more.
    LaunchPacer paceA_, paceB_;
};

FibreBatches::FibreBatches(const FibreRunInput& in, const DeviceFibreInputs& device, const BatchSizes& sizes, FibreRunOutput& out, FibreRunReport& report)
    : in_(in),
      context_(device.context),
      strands_(static_cast<FibreStrand*>(device.strands.data())),
      stride_(sizes.stride),
      out_(out),
      report_(report),
      paceA_(std::min(sizes.a, 1ll << 18), 1024, sizes.a),
      paceB_(std::min(sizes.b, 1ll << 15), 256, sizes.b) {
    const long long bA = sizes.a, bB = sizes.b;
    code_.alloc(bA);
    flags_.alloc(bA);
    seeds_.alloc(bA * sizeof(FibreSeed));
    survivors_.alloc(bA * sizeof(FibreSeed));
    survive_.alloc(bA);
    survivorCount_.alloc(sizeof(long long));
    codeB_.alloc(bB);
    truncatedB_.alloc(bB);
    fibresB_.alloc(bB * sizeof(knit::Fibre));
    accepted_.alloc(bB * sizeof(knit::Fibre));
    slotPoints_.alloc(bB * stride_ * sizeof(knit::Vec3f));
    packed_.alloc(bB * stride_ * sizeof(knit::Vec3f));
    pointCount_.alloc(bB * sizeof(int));
    offset_.alloc(bB * sizeof(int));
    acceptedFlags_.alloc(bB);
    acceptedCount_.alloc(sizeof(long long));
    size_t bytes = 0;
    check(cub::DeviceSelect::Flagged(nullptr, tempBytes_, static_cast<FibreSeed*>(seeds_.data()), survive_.chars(), static_cast<FibreSeed*>(survivors_.data()), survivorCount_.longs(), bA),
          "GPU fibres");
    check(cub::DeviceSelect::Flagged(nullptr, bytes, static_cast<knit::Fibre*>(fibresB_.data()), acceptedFlags_.chars(), static_cast<knit::Fibre*>(accepted_.data()), acceptedCount_.longs(), bB),
          "GPU fibres");
    tempBytes_ = std::max(tempBytes_, bytes);
    check(cub::DeviceScan::ExclusiveSum(nullptr, bytes, pointCount_.ints(), offset_.ints(), bB), "GPU fibres");
    tempBytes_ = std::max(tempBytes_, bytes);
    temp_.alloc(tempBytes_);
}

void FibreBatches::run() {
    out_.code.assign(static_cast<size_t>(in_.candidates), outcomeCode(FibreOutcome::NotProcessed));
    out_.flags.assign(static_cast<size_t>(in_.candidates), 0);
    out_.fibres.clear();
    for (long long first = 0, n = 0; first < in_.candidates; first += n) {
        if (in_.checkpoint) in_.checkpoint();
        n = std::min(paceA_.next(), in_.candidates - first);
        const long long survivors = placeCandidates(first, n);
        // Progress over the fixed candidate range: kernel A's batch [first, first + n) is covered once
        // its survivors are (kernel B below fills the batch's share in proportion).
        if (survivors == 0 && in_.progress) in_.progress(double(first + n) / double(in_.candidates));
        for (long long grown = 0, m = 0; grown < survivors; grown += m) {
            if (in_.checkpoint) in_.checkpoint();
            m = std::min(paceB_.next(), survivors - grown);
            growBatch(grown, m);
            if (in_.progress) in_.progress((double(first) + double(n) * double(grown + m) / double(survivors)) / double(in_.candidates));
        }
    }
    // Every candidate must have a final code (no lost work).
    for (size_t q = 0; q < out_.code.size(); ++q)
        if (out_.code[q] == outcomeCode(FibreOutcome::NotProcessed)) throw DeviceError("GPU fibres: candidate " + std::to_string(q) + " was not processed (internal error)");
    report_.slowestLaunchMs = std::max(paceA_.slowestMs(), paceB_.slowestMs());
}

// Kernel A over candidates [first, first + n), its survivors compacted; returns how many.
long long FibreBatches::placeCandidates(long long first, long long n) {
    const cudaStream_t stream = timing_.stream();
    ++report_.batchesA;
    timing_.begin();
    placeCandidatesKernel<<<blocksFor(n), kThreadsPerBlock, 0, stream>>>(
        context_, {strands_, static_cast<int>(in_.strands.size())}, first, n, {code_.chars(), flags_.chars(), static_cast<FibreSeed*>(seeds_.data()), survive_.chars()});
    checkLaunch("GPU fibres: kernel A");
    const double ms = timing_.finish();
    report_.kernelAMs += ms;
    paceA_.measured(n, ms);
    timing_.begin();
    size_t bytes = tempBytes_;
    check(cub::DeviceSelect::Flagged(temp_.data(), bytes, static_cast<FibreSeed*>(seeds_.data()), survive_.chars(), static_cast<FibreSeed*>(survivors_.data()), survivorCount_.longs(), n, stream),
          "GPU fibres: select");
    report_.selectMs += timing_.finish();
    const auto started = Clock::now();
    long long survivors = 0;
    survivorCount_.download(&survivors, sizeof survivors);
    code_.download(out_.code.data() + first, static_cast<size_t>(n));
    flags_.download(out_.flags.data() + first, static_cast<size_t>(n));
    hostSurvivors_.resize(static_cast<size_t>(survivors));
    survivors_.download(hostSurvivors_.data(), static_cast<size_t>(survivors) * sizeof(FibreSeed));
    report_.downloadMs += elapsedMs(started);
    out_.exposed += survivors;
    return survivors;
}

// Kernel B over survivors [firstSurvivor, firstSurvivor + m): the accepted fibres and their packed
// points on the host, in candidate order.
void FibreBatches::growBatch(long long firstSurvivor, long long m) {
    const cudaStream_t stream = timing_.stream();
    knit::Fibre* fibres = static_cast<knit::Fibre*>(fibresB_.data());
    knit::Vec3f* slots = static_cast<knit::Vec3f*>(slotPoints_.data());
    ++report_.batchesB;
    timing_.begin();
    growFibresKernel<<<blocksFor(m), kThreadsPerBlock, 0, stream>>>(context_, strands_, static_cast<FibreSeed*>(survivors_.data()) + firstSurvivor, m,
                                                                    {stride_, codeB_.chars(), truncatedB_.chars(), fibres, slots, pointCount_.ints(), acceptedFlags_.chars()});
    checkLaunch("GPU fibres: kernel B");
    const double ms = timing_.finish();
    report_.kernelBMs += ms;
    paceB_.measured(m, ms);
    timing_.begin();
    size_t bytes = tempBytes_;
    check(cub::DeviceScan::ExclusiveSum(temp_.data(), bytes, pointCount_.ints(), offset_.ints(), m, stream), "GPU fibres: offsets");
    packPointsKernel<<<blocksFor(m), kThreadsPerBlock, 0, stream>>>({m, stride_, pointCount_.ints(), offset_.ints(), slots}, fibres, static_cast<knit::Vec3f*>(packed_.data()));
    checkLaunch("GPU fibres: pack");
    bytes = tempBytes_;
    check(cub::DeviceSelect::Flagged(temp_.data(), bytes, fibres, acceptedFlags_.chars(), static_cast<knit::Fibre*>(accepted_.data()), acceptedCount_.longs(), m, stream), "GPU fibres: compact");
    report_.compactMs += timing_.finish();
    const auto started = Clock::now();
    long long nAccepted = 0;
    acceptedCount_.download(&nAccepted, sizeof nAccepted);
    hostCodeB_.resize(static_cast<size_t>(m));
    hostTruncatedB_.resize(static_cast<size_t>(m));
    codeB_.download(hostCodeB_.data(), static_cast<size_t>(m));
    truncatedB_.download(hostTruncatedB_.data(), static_cast<size_t>(m));
    hostAccepted_.resize(static_cast<size_t>(nAccepted));
    accepted_.download(hostAccepted_.data(), static_cast<size_t>(nAccepted) * sizeof(knit::Fibre));
    size_t nPoints = 0;
    for (const knit::Fibre& f : hostAccepted_) nPoints += static_cast<size_t>(f.pointCount);
    const size_t base = out_.points.size();
    out_.points.resize(base + nPoints);
    packed_.download(out_.points.data() + base, nPoints * sizeof(knit::Vec3f));
    for (long long q = 0; q < m; ++q) {
        const long long candidate = hostSurvivors_[static_cast<size_t>(firstSurvivor + q)].cand;
        out_.code[static_cast<size_t>(candidate)] = hostCodeB_[static_cast<size_t>(q)];
        if (hostTruncatedB_[static_cast<size_t>(q)]) out_.flags[static_cast<size_t>(candidate)] |= 2;
    }
    for (knit::Fibre& f : hostAccepted_) f.firstPoint += static_cast<uint32_t>(base);
    out_.fibres.insert(out_.fibres.end(), hostAccepted_.begin(), hostAccepted_.end());
    report_.downloadMs += elapsedMs(started);
}

void runFibres(const FibreRunInput& in, FibreRunOutput& out, FibreRunReport& report) {
    // Resident meshes (hidden tests).
    const auto residentStart = Clock::now();
    std::vector<std::shared_ptr<DeviceMesh>> keepAlive;
    std::vector<MeshView> hiddenViews;
    size_t residentBytes = 0;
    for (const MeshQuery* query : in.hidden) {
        keepAlive.push_back(residentMesh(*query));
        hiddenViews.push_back(keepAlive.back()->view);
        residentBytes += keepAlive.back()->bytes;
    }
    report.residentMs = elapsedMs(residentStart);
    report.residentMB = residentBytes / 1048576.0;
    const BatchSizes sizes = batchSizes(in);
    const size_t input = inputBytes(in, hiddenViews.size()), reserved = input + batchBytes(sizes);
    report.inputMB = input / 1048576.0;
    report.reservedMB = reserved / 1048576.0;
    Reservation lease(reserved);
    report.waitMs = lease.waitedMs();
    const auto uploadStart = Clock::now();
    const DeviceFibreInputs device(in, hiddenViews);
    report.uploadMs = elapsedMs(uploadStart);
    FibreBatches(in, device, sizes, out, report).run();
    report.peakMB = memStats().peak / 1048576.0;
}

void checkRequest(const CurveSet& cs, const SurfaceModel& sm, const FibreParams& fp, const FibreFrameGpu& frame, bool withOutcomes) {
    if (2 + std::max(fp.maxSegments, std::max(2, fp.minSegments)) > kMaxFibreVerts)
        throw DeviceError("GPU fibres: maxSegments " + std::to_string(fp.maxSegments) + " exceeds the device polyline bound (" + std::to_string(kMaxFibreVerts - 2) + ")");
    if (frame.plies && (!frame.parentCurves || !frame.parentSurface)) throw DeviceError("GPU fibres: plies need their parent curves and surface");
    if (withOutcomes && (frame.keyLo > 0.0f || frame.keyHi <= 1.0f)) throw DeviceError("GPU fibres: the candidate outcome trace needs the whole key range (no budget window)");
    for (size_t k = 0; k < sm.yarns.size(); ++k)
        if (sm.yarns[k].yarnId != static_cast<int>(k) || sm.yarns.size() != cs.yarns.size())
            throw DeviceError("GPU fibres: the surface model's yarns must be the curve set's, in order");
}

// The curves, surface model, capsules and outward mesh on the host.
std::shared_ptr<GpuFibreInputs> flattenInputs(const CurveSet& cs, const SurfaceModel& sm, const FibreFrameGpu& frame) {
    auto inputs = std::make_shared<GpuFibreInputs>();
    FibreRunInput& f = inputs->in;
    f.cs = flattenCurves(cs, ArcBounds::Omitted);
    f.sm = flattenSurface(sm);
    f.index = flattenIndex(sm);
    if (frame.islands)
        for (const std::vector<IslandWeight>& strand : *frame.islands) f.islands.insert(f.islands.end(), strand.begin(), strand.end());
    if (frame.outward) {
        inputs->outward = flattenMesh(*frame.outward);
        f.outward = &inputs->outward;
    }
    if (frame.plies) {
        f.parents = true;
        f.parentCs = flattenCurves(*frame.parentCurves, ArcBounds::Omitted);
        f.parentSm = flattenSurface(*frame.parentSurface);
        f.plyYarnOf = frame.plies->yarnOf;
        f.plyOffset = frame.plies->parentOffset;
        f.plyScale = frame.plies->parentScale;
    }
    return inputs;
}

FibreParamsFlat flatParameters(const FibreParams& fp, const PatternParams& pp, const FibreFrameGpu& frame, const FibreRunInput& in) {
    FibreParamsFlat p{};
    p.density = fp.density;
    p.budget = frame.candidates.empty() ? 0 : 1;
    p.frontCos = float(fp.frontCos);
    p.keepProbability = float(fp.keepProbability);
    p.curvatureWeight = float(fp.curvatureWeight);
    p.length = float(fp.length);
    p.lengthMean = float(fp.lengthMean);
    p.lengthMin = float(fp.lengthMin);
    p.lengthMax = float(fp.lengthMax);
    p.lengthSpread = float(fp.lengthSpread);
    p.thickness = float(fp.thickness);
    p.thicknessSigma = float(fp.thicknessSigma);
    p.taper = float(fp.taper);
    p.rootDepth = float(fp.rootDepth);
    p.helixDeg = float(fp.helixDeg);
    p.helixSigmaDeg = float(fp.helixSigmaDeg);
    p.liftMinDeg = float(fp.liftMinDeg);
    p.liftMaxDeg = float(fp.liftMaxDeg);
    p.bowTheta = float(fp.bowTheta);
    p.bowRadial = float(fp.bowRadial);
    p.hugHeight = float(fp.hugHeight);
    p.plyPitch = float(fp.plyPitch);
    p.segLenTarget = float(fp.segLenTarget);
    p.maxCurveDev = float(fp.maxCurveDev);
    p.minTurnDeg = float(fp.minTurnDeg);
    p.maxTurnDeg = float(fp.maxTurnDeg);
    p.attachMax = float(fp.attachMax);
    p.sampleStepFactor = float(fp.sampleStepFactor);
    p.serializedMargin = float(fp.serializedMargin);
    p.exportChordFactor = float(fp.exportChordFactor);
    p.rootTurnDeg = float(fp.rootTurnDeg);
    p.minEmergeDeg = float(fp.minEmergeDeg);
    p.plyOuterCos = float(fp.plyOuterCos);
    p.randomFraction = float(fp.randomFraction);
    p.randomFractionD = fp.randomFraction;
    p.randomElevMinDeg = float(fp.randomElevMinDeg);
    p.randomElevMaxDeg = float(fp.randomElevMaxDeg);
    p.randomBend = float(fp.randomBend);
    p.randomSpreadDeg = float(fp.randomSpreadDeg);
    p.rCurlRadius = float(fp.randomCurlRadius);
    p.rCurlLength = float(fp.randomCurlLength);
    p.crimpAmplitude = float(fp.crimpAmplitude);
    p.crimpWavelength = float(fp.crimpWavelength);
    const FibreParams rv = populationParams(fp, FibrePopulation::Random);
    p.rLengthMean = float(rv.lengthMean);
    p.rLengthMin = float(rv.lengthMin);
    p.rLengthMax = float(rv.lengthMax);
    p.rLengthSpread = float(rv.lengthSpread);
    p.rThickness = float(fp.randomThickness);
    p.rThicknessSigma = float(rv.thicknessSigma);
    p.handedness = fp.handedness;
    p.minSegments = fp.minSegments;
    p.maxSegments = fp.maxSegments;
    p.hasRandomOverrides = hasRandomOverrides(fp) ? 1 : 0;
    p.lipschitz = fp.lipschitzSampling ? 1 : 0;
    p.light = fp.lightChecks ? 1 : 0;
    p.drape = fp.drape ? 1 : 0;
    p.seed = fibreSeed(fp, frame.pieceSeed);
    p.keyLo = frame.keyLo;
    p.keyHi = frame.keyHi;
    p.r = float(pp.radius());
    p.ownWindow = float(3 * pp.yarnDiameter);
    p.kappaRef = float(1.0 / (4 * pp.radius()));
    p.surf = frame.outward && frame.exposureLength > 0 ? 1 : 0;
    p.exposureLength = float(frame.exposureLength);
    p.backDepth = float(frame.backDepth);
    p.zTop = in.index.zTop;
    p.zBottom = in.index.zBottom;
    p.hidden = frame.hidden.empty() ? 0 : 1;
    p.hidScale = float(frame.hiddenScale);
    p.hidOffset[0] = float(frame.hiddenOffset.x);
    p.hidOffset[1] = float(frame.hiddenOffset.y);
    p.hidOffset[2] = float(frame.hiddenOffset.z);
    p.quantize = frame.quantize ? 1 : 0;
    p.qSigma = float(frame.quantSigma);
    p.qInvSigma = float(1.0 / frame.quantSigma);
    p.qCentre[0] = float(frame.quantCentre.x);
    p.qCentre[1] = float(frame.quantCentre.y);
    p.qCentre[2] = float(frame.quantCentre.z);
    p.parents = in.parents ? 1 : 0;
    return p;
}

// Strands and their candidate counts (floor(density L) or the budget); returns the total.
long long listStrands(const SurfaceModel& sm, const FibreParams& fp, const FibreFrameGpu& frame, std::vector<FibreStrand>& strands) {
    strands.clear();
    long long total = 0;
    for (const YarnSurface& ys : sm.yarns) {
        FibreStrand st;
        st.first = total;
        st.yarn = ys.yarnId;
        st.n = !frame.candidates.empty() ? frame.candidates[ys.yarnId] : static_cast<int>(std::floor(fp.density * ys.length));
        st.n = std::max(0, st.n);
        total += st.n;
        strands.push_back(st);
    }
    return total;
}

// Fibres in candidate order (= the CPU's strand / candidate order), the report from the per-candidate
// outcomes.
void mergeFibres(FibreRunOutput& run, FibreSet& out, std::vector<uint8_t>* outcomes, FibreRunReport& report) {
    out = FibreSet{};
    out.rep.pop.resize(2);
    out.fibres = std::move(run.fibres);
    out.points = std::move(run.points);
    FibreReport& R = out.rep;
    long long inWindow = 0;
    for (size_t q = 0; q < run.code.size(); ++q) {
        const int code = run.code[q];
        if (code == kOutOfWindowCode) continue;
        ++inWindow;
        const FibreOutcome outcome = static_cast<FibreOutcome>(code);
        if (outcome >= FibreOutcome::NotProcessed && outcome != FibreOutcome::SerializedRejected) throw DeviceError("GPU fibres: unexpected outcome code " + std::to_string(code));
        R.addCandidate((run.flags[q] & 1) ? FibrePopulation::Random : FibrePopulation::Directed, outcome, (run.flags[q] & 2) != 0);
    }
    resumFibreLengths(out);
    if (outcomes) {  // FibreFrame::outcomes: bit 7 = the Random population
        *outcomes = std::move(run.code);
        for (size_t q = 0; q < outcomes->size(); ++q) (*outcomes)[q] |= (run.flags[q] & 1) << 7;
    }
    report.candidates = inWindow;
    report.exposed = run.exposed;
    report.accepted = static_cast<long long>(out.fibres.size());
    report.rootCheckSkipped = report.accepted;
}

} // namespace

bool generateFibres(const FibreRequest& request, FibreSet& out, FibreRunReport* reportOut, std::vector<uint8_t>* outcomes, std::string* err) try {
    const CurveSet& cs = request.curves;
    const SurfaceModel& sm = request.surface;
    const FibreParams& fp = request.fibres;
    const FibreFrameGpu& frame = request.frame;
    const auto start = Clock::now();
    FibreRunReport report;
    if (!available(err)) return false;
    checkRequest(cs, sm, fp, frame, outcomes != nullptr);
    // Built once per FibreFrameGpu::inputs holder and reused by its later calls (same cs / sm / meshes).
    std::shared_ptr<GpuFibreInputs> inputs = frame.inputs ? *frame.inputs : nullptr;
    if (!inputs) inputs = flattenInputs(cs, sm, frame);
    if (frame.inputs) *frame.inputs = inputs;
    FibreRunInput& in = inputs->in;
    in.hidden = frame.hidden;
    in.p = flatParameters(fp, request.pattern, frame, in);
    in.candidates = listStrands(sm, fp, frame, in.strands);
    in.checkpoint = frame.checkpoint;
    in.progress = frame.progress;
    report.prepMs = elapsedMs(start);
    FibreRunOutput run;
    if (in.candidates > 0) runFibres(in, run, report);
    const auto mergeStart = Clock::now();
    mergeFibres(run, out, outcomes, report);
    report.mergeMs = elapsedMs(mergeStart);
    report.totalMs = elapsedMs(start);
    out.rep.ms = report.totalMs;
    if (reportOut) *reportOut = report;
    return true;
} catch (const DeviceError& e) {
    return reportFailure(e, err);
}

namespace {

// One point and starting triangle per thread: both searches' results (compareMeshClosest).
__global__ void compareMeshClosestKernel(MeshView m, const float* points, const int* startTris, int n, MeshClosestPair* out) {
    const int i = threadItem();
    if (i >= n) return;
    const f3 q = fromF(points + 3 * i);
    const MeshHit hits[2] = {meshClosest(m, q), meshClosestFrom(m, q, startTris[i])};
    MeshClosestPair& o = out[i];
    for (int k = 0; k < 2; ++k) {
        o.tri[k] = hits[k].tri;
        o.dist[k] = hits[k].dist;
        for (int c = 0; c < 3; ++c) o.bary[k][c] = hits[k].bary[c];
        o.p[k] = knit::Vec3f{hits[k].p.x, hits[k].p.y, hits[k].p.z};
    }
}

} // namespace

bool compareMeshClosest(const TriMesh& mesh, const std::vector<Vec3>& points, const std::vector<int>& startTris, std::vector<MeshClosestPair>& out,
                        std::string* err) try {
    if (!available(err)) return false;
    if (points.size() != startTris.size()) throw DeviceError("compareMeshClosest: one starting triangle per point");
    const int n = static_cast<int>(points.size());
    DeviceMesh m;
    m.upload(flattenMesh(mesh));
    std::vector<float> p;
    appendFloat3(points, p);
    DeviceBuffer pointsDev, trisDev, outDev;
    pointsDev.upload(p);
    trisDev.upload(startTris);
    outDev.alloc(static_cast<size_t>(n) * sizeof(MeshClosestPair));
    compareMeshClosestKernel<<<blocksFor(n), kThreadsPerBlock>>>(m.view, pointsDev.floats(), trisDev.ints(), n, static_cast<MeshClosestPair*>(outDev.data()));
    check(cudaGetLastError(), "compareMeshClosest");
    check(cudaDeviceSynchronize(), "compareMeshClosest");
    out.resize(static_cast<size_t>(n));
    outDev.download(out.data(), out.size() * sizeof(MeshClosestPair));
    return true;
} catch (const DeviceError& e) {
    return reportFailure(e, err);
}

} // namespace knit::gpu
