// Gate (a): checkCurvature (src/knitting/checks/curvature.cpp) per span on the device
// (checks.h). One thread per global span runs the CPU's piece
// subdivision (explicit stack, same push / pop order, maxDepth 18) and writes the span's
// accumulator; the host merges the spans in span order with the CPU's strict comparisons.
//
// Float allowances (the CPU's are double rounding allowances):
//   kRel       1e-6  relative eta of |C'| and |C''| (CPU 1e-12)
//   floorAbs   1e-8 Lref (CPU 1e-15 Lref)
//   degenerate hodograph triangle: |n|^2 <= 1e-10 |ab|^2 |ac|^2 + 1e-30 (CPU 1e-24, 1e-300)
//   regularity FAIL: sampled speed <= 1e-6 Lref (CPU 1e-12 Lref)
// kappaHi is formed in double from the float numerator and speed bound (no float overflow).
#include <cmath>
#include <cstdio>

#include "gpu/validation/checks.h"
#include "gpu/geometry/math.cuh"

namespace knit::gpu {

namespace {


constexpr float kRel = 1e-6f;
constexpr int kMaxDepth = 18;

struct SpanCurvature {
    double speedLo, kappaHi, spanKappa, kappaWitnessLo, uWitness, rhoEst, uRhoEst, sampledSpeedMin;
    int leaves;
    int pad;
};

__device__ float originSegmentDistance(f3 a, f3 b) {
    const f3 ab = b - a;
    const float l2 = dot(ab, ab);
    const float t = l2 > 0 ? clampf(-dot(a, ab) / l2, 0.0f, 1.0f) : 0.0f;
    return norm(a + ab * t);
}

__device__ float originTriangleDistance(f3 a, f3 b, f3 c) {
    const f3 ab = b - a, ac = c - a, ap = f3{0, 0, 0} - a;
    const f3 n = cross(ab, ac);
    if (dot(n, n) <= 1e-10f * dot(ab, ab) * dot(ac, ac) + 1e-30f)
        return fminf_(fminf_(originSegmentDistance(a, b), originSegmentDistance(a, c)), originSegmentDistance(b, c));
    const float d1 = dot(ab, ap), d2 = dot(ac, ap);
    if (d1 <= 0 && d2 <= 0) return norm(a);
    const f3 bp = f3{0, 0, 0} - b;
    const float d3 = dot(ab, bp), d4 = dot(ac, bp);
    if (d3 >= 0 && d4 <= d3) return norm(b);
    const float vc = d1 * d4 - d3 * d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) return norm(a + ab * (d1 / (d1 - d3)));
    const f3 cp = f3{0, 0, 0} - c;
    const float d5 = dot(ab, cp), d6 = dot(ac, cp);
    if (d6 >= 0 && d5 <= d6) return norm(c);
    const float vb = d5 * d2 - d1 * d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) return norm(a + ac * (d2 / (d2 - d6)));
    const float va = d3 * d6 - d5 * d4;
    if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) return norm(b + (c - b) * ((d4 - d3) / ((d4 - d3) + (d5 - d6))));
    const float denom = 1.0f / (va + vb + vc);
    const float v = vb * denom, w = vc * denom;
    return norm(a + ab * v + ac * w);
}

// QuadBezier::sub of the hodograph (h0, h1, h2) at (t0, t1).
__device__ __forceinline__ void quadSub(const f3 h[3], float t0, float t1, f3 o[3]) {
    auto blossom = [&](float s, float t) { return lerp(lerp(h[0], h[1], s), lerp(h[1], h[2], s), t); };
    o[0] = blossom(t0, t0);
    o[1] = blossom(t0, t1);
    o[2] = blossom(t1, t1);
}

__global__ void curvKernel(SpanList spans, float floorAbs, float kappaFloorScale, SpanCurvature* out, const volatile int* abortFlag) {
    const CurveSetView& cs = spans.curves;
    const int g = threadItem();
    if (g >= spans.count) return;
    const int y = spans.yarn[g];
    const CurveRecord fc = cs.curves[y];
    const int k = g - fc.spanOff;
    const float* ctrl = cs.ctrl + 3 * fc.ctrlOff;
    const Bezier sp = localSpan(ctrl, k).z;
    const f3 hod[3] = {(sp.b[1] - sp.b[0]) * 3.0f, (sp.b[2] - sp.b[1]) * 3.0f, (sp.b[3] - sp.b[2]) * 3.0f};
    const f3 d2At0 = (sp.b[2] - sp.b[1] * 2.0f + sp.b[0]) * 6.0f, d2At1 = (sp.b[3] - sp.b[2] * 2.0f + sp.b[1]) * 6.0f;
    SpanCurvature r;
    r.speedLo = 1e300;
    r.kappaHi = 0;
    r.spanKappa = 0;
    r.kappaWitnessLo = 0;
    r.uWitness = 0;
    r.rhoEst = 1e300;
    r.uRhoEst = 0;
    r.sampledSpeedMin = 1e300;
    r.leaves = 0;
    r.pad = 0;
    const float kappaFloor = kappaFloorScale;  // 0.05 / radius
    struct Item {
        float t0, t1;
        int depth;
    };
    Item stack[kMaxDepth + 2];
    int top = 0;
    stack[top++] = {0.0f, 1.0f, 0};
    unsigned steps = 0;
    while (top > 0) {
        if (abortRequested(abortFlag, ++steps)) break;  // cancelled: the result is discarded
        const Item it = stack[--top];
        // boundPiece
        f3 h[3];
        quadSub(hod, it.t0, it.t1, h);
        const f3 s0 = lerp(d2At0, d2At1, it.t0), s1 = lerp(d2At0, d2At1, it.t1);
        const float hMax = fmaxf_(fmaxf_(norm(h[0]), norm(h[1])), norm(h[2]));
        const float sMax = fmaxf_(norm(s0), norm(s1));
        const float eta1 = kRel * hMax + floorAbs, eta2 = kRel * sMax + floorAbs;
        const f3 b[4] = {cross(h[0], s0), (cross(h[0], s1) + cross(h[1], s0) * 2.0f) / 3.0f, (cross(h[1], s1) * 2.0f + cross(h[2], s0)) / 3.0f,
                         cross(h[2], s1)};
        float numHi = 0;
        for (int q = 0; q < 4; ++q) numHi = fmaxf_(numHi, norm(b[q]));
        numHi += eta1 * sMax + eta2 * hMax + eta1 * eta2;
        const float speedLo = originTriangleDistance(h[0], h[1], h[2]) - eta1;
        const double kappaHi = speedLo > 0 ? double(numHi) / (double(speedLo) * speedLo * speedLo) : 1e300;
        // Sampled curvature (estimate + lower-bound witness).
        float kSamp = 0;
        for (int j = 0; j <= 4; ++j) {
            const float t = it.t0 + (it.t1 - it.t0) * j / 4.0f;
            const f3 d1 = bezD1(sp, t), d2 = bezD2(sp, t);
            const float spd = norm(d1);
            if (spd < r.sampledSpeedMin) r.sampledSpeedMin = spd;
            if (spd <= 0) continue;
            const float num = norm(cross(d1, d2));
            const float kap = num / (spd * spd * spd);
            kSamp = fmaxf_(kSamp, kap);
            const float e1 = kRel * spd + floorAbs, e2 = kRel * norm(d2) + floorAbs;
            const float errNum = e1 * norm(d2) + e2 * spd + e1 * e2;
            const float sp1 = spd + e1;
            const double kLo = double(fmaxf_(0.0f, num - errNum)) / (double(sp1) * sp1 * sp1);
            if (kLo > r.kappaWitnessLo) { r.kappaWitnessLo = kLo; r.uWitness = k + double(t); }
            if (kap > 0 && 1.0 / double(kap) < r.rhoEst) { r.rhoEst = 1.0 / double(kap); r.uRhoEst = k + double(t); }
        }
        const bool tight = speedLo > 0 && kappaHi <= fmaxf_(1.01f * kSamp, kappaFloor);
        if (tight || it.depth >= kMaxDepth) {
            if (double(speedLo) < r.speedLo) r.speedLo = speedLo;
            if (kappaHi > r.kappaHi) r.kappaHi = kappaHi;
            const double sk = speedLo > 0 ? kappaHi : 1e300;
            if (sk > r.spanKappa) r.spanKappa = sk;
            ++r.leaves;
        } else {
            const float m = 0.5f * (it.t0 + it.t1);
            stack[top++] = {it.t0, m, it.depth + 1};
            stack[top++] = {m, it.t1, it.depth + 1};
        }
    }
    out[g] = r;
}

} // namespace

std::vector<CurvatureResult> checkCurvatureOnDevice(const UploadedCurves& d, double radius, const Allowances& al) {
    std::vector<CurvatureResult> out;
    const CurveSet& cs = *d.host;
    out.clear();
    const int nSpans = d.nSpans;
    std::vector<SpanCurvature> acc;
    if (nSpans > 0) {
        DeviceBuffer sy, res;
        sy.upload(yarnOfSpans(d));
        res.alloc(sizeof(SpanCurvature) * nSpans);
        curvKernel<<<blocksFor(nSpans), kThreadsPerBlock>>>({d.dev.view, sy.ints(), nSpans}, float(1e-8 * al.Lref), float(0.05 / radius), static_cast<SpanCurvature*>(res.data()), abortFlagDevice());
        // Waits with the abort flag: degenerate spans subdivide to the full depth (seconds in one
        // launch on a widely failing knit).
        checkLaunch("curvature kernel");
        waitForKernels("curvature kernel");
        acc.resize(nSpans);
        res.download(acc.data(), acc.size() * sizeof(acc[0]));
    }
    const double zeroSpeed = 1e-6 * al.Lref;
    for (size_t y = 0; y < cs.yarns.size(); ++y) {
        CurvatureResult r;
        r.yarnId = cs.yarns[y].id;
        r.speedLo = 1e300;
        r.rhoEst = 1e300;
        double sampledSpeedMin = 1e300;
        const int k0 = d.dev.spanOffset[y], k1 = d.dev.spanOffset[y + 1];
        r.spanKappaHi.assign(static_cast<size_t>(k1 - k0), 0.0);
        for (int g = k0; g < k1; ++g) {
            const SpanCurvature& a = acc[g];
            r.spanKappaHi[g - k0] = a.spanKappa;
            r.speedLo = std::min(r.speedLo, a.speedLo);
            r.kappaHi = std::max(r.kappaHi, a.kappaHi);
            r.leaves += a.leaves;
            sampledSpeedMin = std::min(sampledSpeedMin, a.sampledSpeedMin);
            if (a.kappaWitnessLo > r.kappaWitnessLo) { r.kappaWitnessLo = a.kappaWitnessLo; r.uWitness = a.uWitness; }
            if (a.rhoEst < r.rhoEst) { r.rhoEst = a.rhoEst; r.uRhoEst = a.uRhoEst; }
            if (!std::isfinite(a.speedLo) || !std::isfinite(a.kappaHi)) {
                throw DeviceError("curvature: non-finite span bound (yarn " + std::to_string(y) + ", span " + std::to_string(g - k0) + ")");
            }
        }
        r.rhoLo = r.kappaHi > 0 ? 1.0 / r.kappaHi : 1e300;
        if (r.speedLo > 0) r.regular = Status::Pass;
        else if (sampledSpeedMin <= zeroSpeed) r.regular = Status::Fail;
        else r.regular = Status::Unresolved;
        if (r.regular == Status::Pass && r.rhoLo >= radius) r.noFold = Status::Pass;
        else if (r.kappaWitnessLo > 1.0 / radius) r.noFold = Status::Fail;
        else r.noFold = Status::Unresolved;
        out.push_back(std::move(r));
    }
    return out;
}


} // namespace knit::gpu
