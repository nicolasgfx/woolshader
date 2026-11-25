// The arc-length tables: YarnPath::finalize (src/model/yarn_path.cpp) on the device
// (device_geometry.h computeArcLengthTables), for every yarn of a curve set at once.
//
// Per global span (one thread): the CPU's adaptive G7/K15 over [k, k + 1] (an explicit stack,
// left leaves first, as the recursion) and the whole-span bounds (spanBounds, splits 8). Two
// passes: leaf counts, then the leaves at their offsets. The host sums the leaves in span and
// leaf order (double), builds the Kahan prefix sums of the span bounds and the stitch flags,
// and installs the tables with YarnPath::setTables.
//
// Float tolerance (the CPU's tolS / n = 1e-9 / n is below float resolution, every span would
// subdivide to the depth limit): a leaf is accepted when err <= tolS / n (unchanged, also
// split by half per level) OR err <= 2e-6 |integral of the leaf| (~32 float roundings of the
// fifteen-point sum), and the depth limit is 16 (CPU 30: below 2^-16 the local parameter
// interval approaches float resolution). Leaf integrals and errors are summed in double.
#include <cmath>

#include "gpu/geometry/device_geometry.h"
#include "gpu/geometry/math.cuh"

namespace knit::gpu {

namespace {


constexpr int kMaxDepth = 16;
constexpr float kRelTol = 2e-6f;

// G7/K15 over [a, b] in local parameters: value and |K - G| h.
__device__ float gk15Err(const Bezier& z, float a, float b, float& err) {
    const float xgk[8] = {0.991455371120812639f, 0.949107912342758525f, 0.864864423359769073f, 0.741531185599394440f,
                          0.586087235467691130f, 0.405845151377397167f, 0.207784955007898468f, 0.0f};
    const float wgk[8] = {0.022935322010529225f, 0.063092092629978553f, 0.104790010322250184f, 0.140653259715525919f,
                          0.169004726639267903f, 0.190350578064785410f, 0.204432940075298892f, 0.209482141084727828f};
    const float wg[4] = {0.129484966168869693f, 0.279705391489276668f, 0.381830050505118945f, 0.417959183673469388f};
    const float h = 0.5f * (b - a), m = 0.5f * (a + b);
    const float fc = norm(bezD1(z, m));
    float K = fc * wgk[7], G = fc * wg[3];
    for (int j = 0; j < 7; ++j) {
        const float f1 = norm(bezD1(z, m - h * xgk[j])), f2 = norm(bezD1(z, m + h * xgk[j]));
        K += wgk[j] * (f1 + f2);
        if (j % 2 == 1) G += wg[j / 2] * (f1 + f2);
    }
    err = fabsf((K - G) * h);
    return K * h;
}

// Runs the adaptive rule of one span; with `out` null only counts the leaves.
__device__ int adaptiveSpan(const Bezier& z, int k, float tol0, double* nodes, double* vals, double* errs) {
    struct Item {
        float a, b, tol;
        int depth;
    };
    Item stack[kMaxDepth + 2];
    int n = 0, leaves = 0;
    stack[n++] = {0.0f, 1.0f, tol0, 0};
    while (n > 0) {
        const Item it = stack[--n];
        float err;
        const float v = gk15Err(z, it.a, it.b, err);
        if (err <= it.tol || err <= kRelTol * fabsf(v) || it.depth >= kMaxDepth) {
            if (nodes) {
                nodes[leaves] = k + double(it.b);
                vals[leaves] = v;
                errs[leaves] = err;
            }
            ++leaves;
            continue;
        }
        const float m = 0.5f * (it.a + it.b);
        stack[n++] = {m, it.b, 0.5f * it.tol, it.depth + 1};
        stack[n++] = {it.a, m, 0.5f * it.tol, it.depth + 1};
    }
    return leaves;
}

__global__ void countKernel(SpanList spans, const float* tol, int* count, float* lo, float* hi) {
    const CurveSetView& cs = spans.curves;
    const int g = threadItem();
    if (g >= spans.count) return;
    const int y = spans.yarn[g];
    const CurveRecord fc = cs.curves[y];
    const int k = g - fc.spanOff;
    const float* ctrl = cs.ctrl + 3 * fc.ctrlOff;
    count[g] = adaptiveSpan(localSpan(ctrl, k).z, k, tol[y], nullptr, nullptr, nullptr);
    spanBounds(ctrl, k, 0.0f, 1.0f, lo[g], hi[g]);
}

struct LeafOutput {
    double* nodes;
    double* vals;
    double* errs;
};

__global__ void leafKernel(SpanList spans, const float* tol, const long long* first, LeafOutput out) {
    const CurveSetView& cs = spans.curves;
    const int g = threadItem();
    if (g >= spans.count) return;
    const int y = spans.yarn[g];
    const CurveRecord fc = cs.curves[y];
    const int k = g - fc.spanOff;
    const long long f = first[g];
    adaptiveSpan(localSpan(cs.ctrl + 3 * fc.ctrlOff, k).z, k, tol[y], out.nodes + f, out.vals + f, out.errs + f);
}

// The leaves of every span and the spans' length bounds, on the host.
struct SpanLeaves {
    std::vector<float> lo, hi;
    std::vector<double> nodes, vals, errs;
    std::vector<long long> first;  // per span + 1
};

SpanLeaves integrateSpans(const CurveArrays& curves, const std::vector<int>& spanYarn, const std::vector<float>& tol) {
    const int spans = static_cast<int>(spanYarn.size());
    SpanLeaves out;
    out.first.assign(static_cast<size_t>(spans) + 1, 0);
    if (spans == 0) return out;
    std::vector<int> count;
    DeviceCurves dc;
    DeviceBuffer sy, tolB, cnt, loB, hiB, firstB, nodeB, valB, errB;
    dc.upload(curves);
    sy.upload(spanYarn);
    tolB.upload(tol);
    cnt.alloc(4 * size_t(spans));
    loB.alloc(4 * size_t(spans));
    hiB.alloc(4 * size_t(spans));
    const SpanList list{dc.view, sy.ints(), spans};
    countKernel<<<blocksFor(spans), kThreadsPerBlock>>>(list, tolB.floats(), cnt.ints(), loB.floats(), hiB.floats());
    checkLaunch("finalize count kernel");
    cnt.download(count, spans);
    loB.download(out.lo, spans);
    hiB.download(out.hi, spans);
    for (int g = 0; g < spans; ++g) out.first[g + 1] = out.first[g] + count[g];
    const size_t leaves = static_cast<size_t>(out.first[spans]);
    firstB.upload(out.first);
    nodeB.alloc(8 * leaves);
    valB.alloc(8 * leaves);
    errB.alloc(8 * leaves);
    leafKernel<<<blocksFor(spans), kThreadsPerBlock>>>(list, tolB.floats(), firstB.longs(), {nodeB.doubles(), valB.doubles(), errB.doubles()});
    checkLaunch("finalize leaf kernel");
    nodeB.download(out.nodes, leaves);
    valB.download(out.vals, leaves);
    errB.download(out.errs, leaves);
    return out;
}

// One yarn's tables: the sums in span and leaf order (as one serial pass of finalize), the prefix
// sums of the span bounds (Kahan) and the stitch flags. Its spans start at global span g0.
YarnPath::Tables tablesOf(const YarnPath& yarn, int y, int g0, const SpanLeaves& leaves) {
    const int n = yarn.curve.numSpans();
    YarnPath::Tables t;
    t.tu.assign(1, 0.0);
    t.ts.assign(1, 0.0);
    double s = 0, errSum = 0;
    for (int g = g0; g < g0 + n; ++g)
        for (long long q = leaves.first[g]; q < leaves.first[g + 1]; ++q) {
            if (!std::isfinite(leaves.vals[q])) throw DeviceError("finalize: non-finite arc length (yarn " + std::to_string(y) + ", span " + std::to_string(g - g0) + ")");
            s += leaves.vals[q];
            t.tu.push_back(leaves.nodes[q]);
            t.ts.push_back(s);
            errSum += leaves.errs[q];
        }
    t.length = s;
    t.lengthErr = errSum;
    t.preLo.assign(n + 1, 0.0);
    t.preHi.assign(n + 1, 0.0);
    double cLo = 0, cHi = 0;
    for (int k = 0; k < n; ++k) {
        const double l = leaves.lo[g0 + k], h = leaves.hi[g0 + k];
        const double yl = l - cLo, tl = t.preLo[k] + yl;
        cLo = (tl - t.preLo[k]) - yl;
        t.preLo[k + 1] = tl;
        const double yh = h - cHi, th = t.preHi[k] + yh;
        cHi = (th - t.preHi[k]) - yh;
        t.preHi[k + 1] = th;
    }
    t.stitchesSorted = true;
    for (size_t i = 0; i < yarn.stitches.size(); ++i) {
        if (yarn.stitches[i].uBegin >= yarn.stitches[i].uEnd) t.stitchesSorted = false;
        if (i > 0 && yarn.stitches[i].uBegin < yarn.stitches[i - 1].uEnd) t.stitchesSorted = false;
    }
    t.sortedCount = yarn.stitches.size();
    return t;
}

} // namespace

void computeArcLengthTables(CurveSet& cs) {
    const int Y = static_cast<int>(cs.yarns.size());
    CurveArrays curves;
    std::vector<int> spanYarn;
    std::vector<float> tol(Y);
    for (int y = 0; y < Y; ++y) {
        const YarnPath& yarn = cs.yarns[y];
        CurveRecord c;
        c.ctrlOff = static_cast<int>(curves.ctrl.size() / 3);
        c.nCtrl = static_cast<int>(yarn.curve.ctrl.size());
        c.spanOff = static_cast<int>(spanYarn.size());
        c.preOff = -1;
        appendFloat3(yarn.curve.ctrl, curves.ctrl);
        const int n = yarn.curve.numSpans();
        spanYarn.insert(spanYarn.end(), n, y);
        tol[y] = float(1e-9 / std::max(1, n));  // YarnPath::finalize(1e-9): tolS / n per span
        curves.curves.push_back(c);
    }
    const SpanLeaves leaves = integrateSpans(curves, spanYarn, tol);
    int g0 = 0;
    for (int y = 0; y < Y; ++y) {
        cs.yarns[y].setTables(tablesOf(cs.yarns[y], y, g0, leaves));
        g0 += cs.yarns[y].curve.numSpans();
    }
}

} // namespace knit::gpu
