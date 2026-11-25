// GPU backend: the placement's surface queries on the device (placement.h makePieceQuery):
// PieceSurface::at (chart evaluation, pole caps by the closest point to the cap chord) and
// PieceSurface::project (Chart::closestPoint: closest point, interpolated vertex normals), in
// batches (BaseSurface::atBatch / projectBatch). Single precision (math.cuh).
#include <chrono>
#include <cmath>
#include <stdexcept>

#include "gpu/geometry/surface.cuh"
#include "gpu/geometry/device_geometry.h"
#include "knitting/knit/piece_knit.h"
#include "gpu/device.h"
#include "gpu/surface/placement.h"

namespace knit::gpu {

namespace {

struct ChartView {
    const float* uv = nullptr;          // 6 per triangle (a, b per corner)
    const unsigned* cellStart = nullptr;  // na * nb + 1
    const int* cellTris = nullptr;
    int na = 0, nb = 0;
    float b0 = 0, b1 = 1;
    float aMin = 0, aMax = 1;
    float poleA[3] = {0, 0, 0}, poleB[3] = {0, 0, 0};
    MeshView mesh;                      // positions, triangles, vertex normals, closest-point grid
};

struct ChartHit {
    bool ok = false;
    f3 p{0, 0, 0}, n{0, 0, 0};
};

// Chart::inTriangle (tolerance: CPU 1e-12, float 1e-6 of the unit (a, b) square).
__device__ __forceinline__ bool inTri(const float* u, float a, float b, float w[3], float tol) {
    const float e1x = u[2] - u[0], e1y = u[3] - u[1], e2x = u[4] - u[0], e2y = u[5] - u[1];
    const float qx = a - u[0], qy = b - u[1];
    const float det = e1x * e2y - e1y * e2x;
    if (det == 0) return false;
    w[1] = (qx * e2y - qy * e2x) / det;
    w[2] = (e1x * qy - e1y * qx) / det;
    w[0] = 1 - w[1] - w[2];
    return w[0] >= -tol && w[1] >= -tol && w[2] >= -tol;
}

__device__ ChartHit chartEvalRaw(const ChartView& c, float a, float b) {
    ChartHit h;
    const float sb = c.nb / (c.b1 - c.b0);
    int i = static_cast<int>(floorf(a * c.na)), j = static_cast<int>(floorf((b - c.b0) * sb));
    i = i < 0 ? 0 : (i > c.na - 1 ? c.na - 1 : i);
    j = j < 0 ? 0 : (j > c.nb - 1 ? c.nb - 1 : j);
    const float tol = 1e-6f;
    if (b < c.b0 - tol || b > c.b1 + tol || a < -tol || a > 1 + tol) return h;
    int best = -1;
    float bestMin = -3.0e38f, wb[3] = {0, 0, 0};
    const unsigned cell = static_cast<unsigned>(i) * c.nb + j;
    for (unsigned q = c.cellStart[cell]; q < c.cellStart[cell + 1]; ++q) {
        const int t = c.cellTris[q];
        float w[3];
        if (!inTri(c.uv + 6 * t, a, b, w, tol)) continue;
        const float m = fminf_(fminf_(w[0], w[1]), w[2]);
        if (m > bestMin) {
            bestMin = m;
            best = t;
            wb[0] = w[0]; wb[1] = w[1]; wb[2] = w[2];
        }
    }
    if (best < 0) return h;
    const int* tv = c.mesh.tri + 3 * best;
    h.ok = true;
    h.p = fromF(c.mesh.pos + 3 * tv[0]) * wb[0] + fromF(c.mesh.pos + 3 * tv[1]) * wb[1] + fromF(c.mesh.pos + 3 * tv[2]) * wb[2];
    h.n = normalized(fromF(c.mesh.vnormal + 3 * tv[0]) * wb[0] + fromF(c.mesh.vnormal + 3 * tv[1]) * wb[1] +
                     fromF(c.mesh.vnormal + 3 * tv[2]) * wb[2]);
    return h;
}

__device__ ChartHit chartEval(const ChartView& c, float a, float b) {
    const float bw = b - floorf(b);
    ChartHit h = chartEvalRaw(c, a, bw);
    if (!h.ok) h = chartEvalRaw(c, a, bw + 1.0f);
    return h;
}

// Chart::closestPoint: the closest point and the interpolated vertex normal (no face-normal
// fallback: the chart's convention).
__device__ ChartHit chartClosest(const ChartView& c, f3 q) {
    ChartHit h;
    const MeshHit m = meshClosest(c.mesh, q);
    if (m.tri < 0) return h;
    const int* tv = c.mesh.tri + 3 * m.tri;
    h.ok = true;
    h.p = m.p;
    h.n = normalized(fromF(c.mesh.vnormal + 3 * tv[0]) * m.bary[0] + fromF(c.mesh.vnormal + 3 * tv[1]) * m.bary[1] +
                     fromF(c.mesh.vnormal + 3 * tv[2]) * m.bary[2]);
    return h;
}

// PieceSurface::at.
__device__ ChartHit pieceAt(const ChartView& c, float a, float b) {
    if (a >= c.aMin && a <= c.aMax) return chartEval(c, a, b);
    const bool capA = a < c.aMin;
    const f3 pole = fromF(capA ? c.poleA : c.poleB);
    const f3 edge = chartEval(c, capA ? c.aMin : c.aMax, b).p;
    const float t = capA ? fmaxf_(0.0f, a) / c.aMin : fmaxf_(0.0f, 1.0f - a) / (1.0f - c.aMax);
    return chartClosest(c, pole + (edge - pole) * t);
}

__global__ void atKernel(ChartView c, const float* ab, int n, float* out) {
    const int t = threadItem();
    if (t >= n) return;
    const ChartHit h = pieceAt(c, ab[2 * t], ab[2 * t + 1]);
    float* o = out + 7 * t;
    o[0] = h.p.x; o[1] = h.p.y; o[2] = h.p.z;
    o[3] = h.n.x; o[4] = h.n.y; o[5] = h.n.z;
    o[6] = h.ok ? 1.0f : 0.0f;
}

__global__ void projectKernel(ChartView c, const float* q, int n, float* out) {
    const int t = threadItem();
    if (t >= n) return;
    const ChartHit h = chartClosest(c, fromF(q + 3 * t));
    float* o = out + 7 * t;
    o[0] = h.p.x; o[1] = h.p.y; o[2] = h.p.z;
    o[3] = h.n.x; o[4] = h.n.y; o[5] = h.n.z;
    o[6] = h.ok ? 1.0f : 0.0f;
}

} // namespace

namespace {

// A PieceSurface's queries on the device. Geometry uploaded once (the object's lifetime).
class GpuPieceQuery final : public BaseSurface {
public:
    GpuPieceQuery(const Chart& ch, PlacementStats* st) : ch_(ch), st_(st) {
        const TriMesh& m = ch.mesh();
        mesh_.upload(flattenMesh(m));
        std::vector<float> uv;
        uv.reserve(m.tri.size() * 6);
        for (const auto& u : m.uv)
            for (int k = 0; k < 3; ++k) uv.insert(uv.end(), {float(u[k].x), float(u[k].y)});
        const Chart::GridView gv = ch.gridView();
        std::vector<unsigned> start{0};
        std::vector<int> tris;
        for (const auto& cell : *gv.cells) {
            tris.insert(tris.end(), cell.begin(), cell.end());
            start.push_back(static_cast<unsigned>(tris.size()));
        }
        uv_.upload(uv);
        cellStart_.upload(start);
        cellTris_.upload(tris);
        c_.uv = uv_.floats();
        c_.cellStart = cellStart_.uints();
        c_.cellTris = cellTris_.ints();
        c_.na = gv.na;
        c_.nb = gv.nb;
        c_.b0 = float(gv.b0);
        c_.b1 = float(gv.b1);
        c_.aMin = float(ch.aMin());
        c_.aMax = float(ch.aMax());
        const Vec3 pa = ch.poleA(), pb = ch.poleB();
        c_.poleA[0] = float(pa.x); c_.poleA[1] = float(pa.y); c_.poleA[2] = float(pa.z);
        c_.poleB[0] = float(pb.x); c_.poleB[1] = float(pb.y); c_.poleB[2] = float(pb.z);
        c_.mesh = mesh_.view;
    }
    SurfaceSample at(double a, double b) const override {
        std::vector<SurfaceSample> o;
        atBatch({{a, b}}, o);
        return o[0];
    }
    SurfaceSample project(const Vec3& q) const override {
        std::vector<SurfaceSample> o;
        projectBatch({q}, o);
        return o[0];
    }
    double aMin() const override { return ch_.aMin(); }
    double aMax() const override { return ch_.aMax(); }
    void atBatch(const std::vector<Vec2>& ab, std::vector<SurfaceSample>& out) const override {
        std::vector<float> in;
        in.reserve(ab.size() * 2);
        for (const Vec2& v : ab) in.insert(in.end(), {float(v.x), float(v.y)});
        runBatch(true, in, ab.size(), out);
    }
    void projectBatch(const std::vector<Vec3>& q, std::vector<SurfaceSample>& out) const override {
        std::vector<float> in;
        in.reserve(q.size() * 3);
        for (const Vec3& v : q) in.insert(in.end(), {float(v.x), float(v.y), float(v.z)});
        runBatch(false, in, q.size(), out);
    }

private:
    void runBatch(bool at, const std::vector<float>& in, size_t n, std::vector<SurfaceSample>& out) const {
        const auto t0 = Clock::now();
        out.assign(n, SurfaceSample{});
        if (n == 0) return;
        Reservation lease(n * 40 + (1 << 20));  // the batch's buffers (the chart is resident)
        DeviceBuffer din, dout;
        din.upload(in.data(), in.size() * 4);
        dout.alloc(n * 7 * 4);
        if (at) atKernel<<<blocksFor(static_cast<long long>(n)), kThreadsPerBlock>>>(c_, din.floats(), static_cast<int>(n), dout.floats());
        else projectKernel<<<blocksFor(static_cast<long long>(n)), kThreadsPerBlock>>>(c_, din.floats(), static_cast<int>(n), dout.floats());
        checkLaunch("GPU placement kernel");
        std::vector<float> h;
        dout.download(h, n * 7);
        long long misses = 0;
        for (size_t i = 0; i < n; ++i) {
            const float* o = h.data() + 7 * i;
            out[i].p = {o[0], o[1], o[2]};
            out[i].n = {o[3], o[4], o[5]};
            misses += o[6] == 0.0f;
        }
        if (st_) {
            ++st_->batches;
            st_->queries += static_cast<long long>(n);
            st_->misses += misses;
            st_->ms += elapsedMs(t0);
        }
    }

    const Chart& ch_;
    PlacementStats* st_;
    DeviceMesh mesh_;
    DeviceBuffer uv_, cellStart_, cellTris_;
    ChartView c_;
};

} // namespace

std::unique_ptr<BaseSurface> makePieceQuery(const Chart& c, PlacementStats* st, std::string* err) {
    if (!available(err)) return nullptr;
    try {
        return std::make_unique<GpuPieceQuery>(c, st);
    } catch (const DeviceError& e) {
        reportFailure(e, err);
        return nullptr;
    }
}

std::string PlacementStats::line() const {
    char b[200];
    std::snprintf(b, sizeof b, "GPU placement queries: %lld in %lld batches, %.0f ms (%lld outside the chart)", queries, batches, ms, misses);
    return b;
}

} // namespace knit::gpu
