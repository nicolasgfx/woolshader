// GPU backend: the tiled pipeline's mesh surface on the device (placement.h makeMeshSurface):
// MeshSurface::project (closest point, MeshSurface::normalAt) and continueFrom (the ownership-
// preserving continuation, surface.cuh) in batches, and the boundary checks (tests/cpu_numerics.cpp:
// B1 certified distance of the centrelines to the open-boundary edges with the Lipschitz interval
// refinement, B2 sampled continuations) on the device. Single precision.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>

#include "gpu/geometry/device_geometry.h"
#include "knitting/layout/uv_atlas.h"
#include "knitting/layout/chart_lattice.h"
#include "gpu/device.h"
#include "gpu/surface/placement.h"
#include "gpu/geometry/surface.cuh"

namespace knit::gpu {

namespace {

__global__ void meshProjectKernel(MeshView m, const float* q, int n, float* out) {
    const int i = threadItem();
    if (i >= n) return;
    const MeshHit h = meshClosest(m, fromF(q + 3 * i));
    float* o = out + 7 * i;
    if (h.tri < 0) {
        o[6] = 0;
        return;
    }
    const f3 nn = meshNormalAt(m, h.tri, h.bary);
    o[0] = h.p.x; o[1] = h.p.y; o[2] = h.p.z;
    o[3] = nn.x; o[4] = nn.y; o[5] = nn.z;
    o[6] = 1;
}

struct ContinuationInput {
    MeshView mesh;
    AdjacencyView adjacency;
    DeviceContinuationParams params;
};

struct ContinuationOutput {
    ProjectionState* next;
    int* results;
};

__global__ void continueKernel(ContinuationInput in, const ProjectionState* prev, const float* q, int n, ContinuationOutput out) {
    const int i = threadItem();
    if (i >= n) return;
    ProjectionState nx;
    out.results[i] = continueFrom(in.mesh, in.adjacency, prev[i], fromF(q + 3 * i), nx, in.params);
    out.next[i] = nx;
}

ProjectionState toDev(const knit::ProjState& s) {
    ProjectionState d{};
    d.tri = s.tri;
    for (int k = 0; k < 3; ++k) d.b[k] = static_cast<float>(s.b[k]);
    d.foot = {float(s.foot.x), float(s.foot.y), float(s.foot.z)};
    d.n = {float(s.n.x), float(s.n.y), float(s.n.z)};
    d.query = {float(s.query.x), float(s.query.y), float(s.query.z)};
    return d;
}
knit::ProjState fromDev(const ProjectionState& d) {
    knit::ProjState s;
    s.tri = d.tri;
    for (int k = 0; k < 3; ++k) s.b[k] = d.b[k];
    s.foot = {d.foot.x, d.foot.y, d.foot.z};
    s.n = {d.n.x, d.n.y, d.n.z};
    s.query = {d.query.x, d.query.y, d.query.z};
    return s;
}
DeviceContinuationParams contParams(const knit::ContinuationParams& cp) {
    DeviceContinuationParams c{};
    c.lambdaW = static_cast<float>(cp.lambdaW);
    c.minNormalCos = static_cast<float>(cp.minNormalCos);
    c.band = cp.band > 1e30 ? INFINITY : static_cast<float>(cp.band);
    c.maxSteps = cp.maxSteps;
    return c;
}

// --- boundary check
struct BoundarySegment {
    float a[3], b[3];
};
// Distance from x to the open-boundary edges within `reach` (reach if none).
__device__ __forceinline__ float boundaryDist(const GridView& g, const BoundarySegment* segs, f3 x, float reach) {
    float d = reach;
    for (GridCursor cursor(g, x - f3{reach, reach, reach}, x + f3{reach, reach, reach}); cursor.next();) {
        const int id = cursor.item();
        const f3 a = fromF(segs[id].a), ab = fromF(segs[id].b) - a;
        const float L2 = dot(ab, ab);
        const float t = L2 > 0 ? clampf(dot(x - a, ab) / L2, 0.0f, 1.0f) : 0.0f;
        d = fminf_(d, norm(x - (a + ab * t)));
    }
    return d;
}

struct BoundaryYarn {
    int yarn;          // curve index
    int sampleOff, n;  // samples k = 0 .. n
    double L;
};

struct BoundarySamples {
    const BoundaryYarn* yarns;
    int nYarns, count;
    double* us;
    float* xs;
};
struct BoundaryEdges {
    GridView grid;
    const BoundarySegment* segments;
    float reach, Rs, minStep;
};
struct BoundaryCounters {
    unsigned long long* minKey;
    int* fail;
    int* unresolved;
    int* stackFull;
};
struct BindingStates {
    const ProjectionState* states;
    const int* offset;
    const int* count;
};

__device__ __forceinline__ int yarnOfSample(const BoundarySamples& s, int sample) {
    int lo = 0, hi = s.nYarns - 1;
    while (lo < hi) {
        const int mid = (lo + hi + 1) >> 1;
        if (s.yarns[mid].sampleOff <= sample) lo = mid; else hi = mid - 1;
    }
    return lo;
}

__global__ void boundarySampleKernel(CurveSetView cs, BoundarySamples in) {
    const int g = threadItem();
    if (g >= in.count) return;
    const int y = yarnOfSample(in, g);
    const BoundaryYarn Y = in.yarns[y];
    const int k = g - Y.sampleOff;
    const CurveView c = devCurve(cs, Y.yarn);
    const double u = paramAtArcLength(c, Y.L * k / Y.n);
    in.us[g] = u;
    const f3 x = curveEval(c, u);
    in.xs[3 * g] = x.x; in.xs[3 * g + 1] = x.y; in.xs[3 * g + 2] = x.z;
}

// B1 per sample (FAIL samples, minimum), and per interval [k, k + 1] the refinement.
__global__ void boundaryDistanceKernel(CurveSetView cs, BoundarySamples in, BoundaryEdges edges, BoundaryCounters out) {
    const GridView& g = edges.grid;
    const BoundarySegment* segs = edges.segments;
    const float reach = edges.reach, Rs = edges.Rs;
    const double* us = in.us;
    const int i = threadItem();
    if (i >= in.count) return;
    const int y = yarnOfSample(in, i);
    const BoundaryYarn Y = in.yarns[y];
    const int k = i - Y.sampleOff;
    const f3 x0 = fromF(in.xs + 3 * i);
    const float d0 = boundaryDist(g, segs, x0, reach);
    atomicMin(out.minKey, static_cast<unsigned long long>(__float_as_uint(d0)));
    if (d0 < Rs) atomicAdd(out.fail, 1);
    if (k >= Y.n) return;  // the last sample has no interval
    const CurveView c = devCurve(cs, Y.yarn);
    const CurveRecord& fc = cs.curves[Y.yarn];
    const double* preLo = cs.preLo && fc.preOff >= 0 ? cs.preLo + fc.preOff : nullptr;
    const double* preHi = cs.preHi && fc.preOff >= 0 ? cs.preHi + fc.preOff : nullptr;
    constexpr int kStack = 48;
    double s0[kStack], s1[kStack];
    int top = 0;
    s0[top] = us[i];
    s1[top] = us[i + 1];
    ++top;
    while (top > 0) {
        --top;
        const double u0 = s0[top], u1 = s1[top];
        const float da = boundaryDist(g, segs, curveEval(c, u0), reach), db = boundaryDist(g, segs, curveEval(c, u1), reach);
        if (da < Rs || db < Rs) continue;
        const double hi = lengthBounds(c, preLo, preHi, u0, u1).hi;
        if (0.5f * (da + db - static_cast<float>(hi)) >= Rs) continue;
        if (hi < edges.minStep) {
            atomicAdd(out.unresolved, 1);
            continue;
        }
        const double um = 0.5 * (u0 + u1);
        const float dm = boundaryDist(g, segs, curveEval(c, um), reach);
        atomicMin(out.minKey, static_cast<unsigned long long>(__float_as_uint(dm)));
        if (dm < Rs) {
            atomicAdd(out.fail, 1);
            continue;
        }
        if (top + 2 > kStack) {  // never silent: counted as UNRESOLVED
            atomicAdd(out.stackFull, 1);
            atomicAdd(out.unresolved, 1);
            continue;
        }
        s0[top] = u0; s1[top] = um; ++top;
        s0[top] = um; s1[top] = u1; ++top;
    }
}

// B2 per sample: the continuation from the nearest control's state.
__global__ void boundaryContinuationKernel(ContinuationInput mesh, BoundarySamples in, BindingStates bound, int* counts) {
    const int i = threadItem();
    if (i >= in.count) return;
    const int y = yarnOfSample(in, i);
    if (bound.count[y] <= 0) return;
    int c = static_cast<int>(llround(in.us[i])) + 1;
    c = c < 0 ? 0 : (c > bound.count[y] - 1 ? bound.count[y] - 1 : c);
    ProjectionState next;
    const int r = continueFrom(mesh.mesh, mesh.adjacency, bound.states[bound.offset[y] + c], fromF(in.xs + 3 * i), next, mesh.params);
    atomicAdd(counts + r, 1);
}

} // namespace

namespace {

class GpuMeshSurface final : public MeshSurface {
public:
    GpuMeshSurface(const TriMesh& m, PlacementStats* st) : MeshSurface(m), st_(st) {
        mesh_ = residentMesh(m, MeshWalk::Adjacency);
    }
    void projectBatch(const std::vector<Vec3>& q, std::vector<SurfaceSample>& out) const override {
        const auto t0 = Clock::now();
        out.assign(q.size(), SurfaceSample{});
        if (q.empty()) return;
        std::vector<float> in;
        in.reserve(q.size() * 3);
        for (const Vec3& v : q) in.insert(in.end(), {float(v.x), float(v.y), float(v.z)});
        Reservation lease(q.size() * 40 + (1 << 20));  // the batch's buffers (the mesh is resident)
        DeviceBuffer din, dout;
        din.upload(in.data(), in.size() * 4);
        dout.alloc(q.size() * 28);
        const int n = static_cast<int>(q.size());
        meshProjectKernel<<<blocksFor(n), kThreadsPerBlock>>>(mesh_->view, din.floats(), n, dout.floats());
        std::vector<float> h;
        checkLaunch("GPU mesh project");
        dout.download(h, q.size() * 7);
        for (size_t i = 0; i < q.size(); ++i) {
            out[i].p = {h[7 * i], h[7 * i + 1], h[7 * i + 2]};
            out[i].n = {h[7 * i + 3], h[7 * i + 4], h[7 * i + 5]};
        }
        account(q.size(), t0);
    }
    void continueBatch(const std::vector<ProjState>& prev, const std::vector<Vec3>& q, const ContinuationParams& cp, std::vector<ProjState>& next,
                       std::vector<ProjResult>& res) const override {
        const auto t0 = Clock::now();
        next.assign(prev.size(), ProjState{});
        res.assign(prev.size(), ProjResult::Budget);
        if (prev.empty()) return;
        std::vector<ProjectionState> hp(prev.size());
        std::vector<float> hq;
        hq.reserve(q.size() * 3);
        for (size_t i = 0; i < prev.size(); ++i) {
            hp[i] = toDev(prev[i]);
            hq.insert(hq.end(), {float(q[i].x), float(q[i].y), float(q[i].z)});
        }
        Reservation lease(prev.size() * (2 * sizeof(ProjectionState) + 16) + (1 << 20));
        DeviceBuffer dp, dq, dn, dr;
        dp.upload(hp.data(), hp.size() * sizeof(hp[0]));
        dq.upload(hq);
        dn.alloc(sizeof(ProjectionState) * prev.size());
        dr.alloc(4 * prev.size());
        const int n = static_cast<int>(prev.size());
        continueKernel<<<blocksFor(n), kThreadsPerBlock>>>({mesh_->view, mesh_->adjacency, contParams(cp)}, static_cast<ProjectionState*>(dp.data()), dq.floats(), n,
                                                           {static_cast<ProjectionState*>(dn.data()), dr.ints()});
        std::vector<ProjectionState> hn;
        std::vector<int> hr;
        checkLaunch("GPU continuation");
        hn.resize(prev.size());
        dn.download(hn.data(), hn.size() * sizeof(hn[0]));
        dr.download(hr, prev.size());
        for (size_t i = 0; i < prev.size(); ++i) {
            res[i] = static_cast<ProjResult>(hr[i]);
            next[i] = res[i] == ProjResult::Owned ? fromDev(hn[i]) : prev[i];
        }
        account(prev.size(), t0);
    }

private:
    void account(size_t n, Clock::time_point t0) const {
        if (!st_) return;
        ++st_->batches;
        st_->queries += static_cast<long long>(n);
        st_->ms += elapsedMs(t0);
    }
    PlacementStats* st_;
    std::shared_ptr<DeviceMesh> mesh_;
};

} // namespace

std::unique_ptr<MeshSurface> makeMeshSurface(const TriMesh& m, PlacementStats* st, std::string* err) {
    if (!available(err)) return nullptr;
    try {
        return std::make_unique<GpuMeshSurface>(m, st);
    } catch (const DeviceError& e) {
        reportFailure(e, err);
        return nullptr;
    }
}

namespace {

// B2: the sampled continuations from the binding's states. Per result the unresolved samples go to
// `counts`; returns how many samples were tested.
long boundaryContinuations(const BoundaryRequest& request, const std::vector<BoundaryYarn>& yarns, const BoundarySamples& samples, int* counts) {
    const CurveSet& cs = request.curves;
    const SurfaceBinding& bind = request.binding;
    std::vector<ProjectionState> states;
    std::vector<int> offset, count;
    long tested = 0;
    for (size_t y = 0; y < cs.yarns.size(); ++y) {
        offset.push_back(static_cast<int>(states.size()));
        const int c = y < bind.state.size() ? static_cast<int>(bind.state[y].size()) : 0;
        count.push_back(c);
        if (c > 0) {
            tested += yarns[y].n + 1;
            for (const ProjState& s : bind.state[y]) states.push_back(toDev(s));
        }
    }
    std::shared_ptr<DeviceMesh> mesh = residentMesh(request.surface.mesh(), MeshWalk::Adjacency);
    DeviceBuffer dStates, dOffset, dCount;
    dStates.upload(states.data(), states.size() * sizeof(states[0]));
    dOffset.upload(offset);
    dCount.upload(count);
    boundaryContinuationKernel<<<blocksFor(samples.count), kThreadsPerBlock>>>({mesh->view, mesh->adjacency, contParams(bind.cp)}, samples,
                                                                               {static_cast<ProjectionState*>(dStates.data()), dOffset.ints(), dCount.ints()}, counts);
    checkLaunch("GPU boundary B2");
    return tested;
}

// The open-boundary edges, as segments and boxes.
std::vector<BoundarySegment> boundarySegments(const MeshSurface& surf, BoxArrays& boxes) {
    const TriMesh& m = surf.mesh();
    std::vector<BoundarySegment> segs;
    for (size_t t = 0; t < m.tri.size(); ++t)
        for (int e = 0; e < 3; ++e)
            if (surf.isBoundaryEdge(static_cast<int>(t), e)) {
                const Vec3 a = m.pos[m.tri[t][e]], b = m.pos[m.tri[t][(e + 1) % 3]];
                segs.push_back({{float(a.x), float(a.y), float(a.z)}, {float(b.x), float(b.y), float(b.z)}});
                const float pad = 1e-5f * (1 + float(std::max({std::fabs(a.x), std::fabs(a.y), std::fabs(a.z)})));
                boxes.lo.insert(boxes.lo.end(), {float(std::min(a.x, b.x)) - pad, float(std::min(a.y, b.y)) - pad, float(std::min(a.z, b.z)) - pad});
                boxes.hi.insert(boxes.hi.end(), {float(std::max(a.x, b.x)) + pad, float(std::max(a.y, b.y)) + pad, float(std::max(a.z, b.z)) + pad});
            }
    return segs;
}

void writeVerdicts(int* counts, double Rs, double step, size_t edges, BoundaryCheck& out) {
    out.b1Fail = counts[0];
    out.b1Unresolved = counts[1];
    for (int k = 0; k < 5; ++k) out.b2Unresolved[k] = k == 0 ? 0 : counts[3 + k];
    out.b1 = out.b1Fail ? Status::Fail : out.b1Unresolved ? Status::Unresolved : Status::Pass;
    const int overhang = out.b2Unresolved[static_cast<int>(ProjResult::OpenBoundary)];
    int other = 0;
    for (int k = 1; k < 5; ++k) other += out.b2Unresolved[k];
    out.b2 = overhang ? Status::Fail : other ? Status::Unresolved : Status::Pass;
    char b[400];
    std::snprintf(b, sizeof b,
                  "B1 (GPU float32 diagnostic) %s: %zu open-boundary edges, min centreline distance %.4f mm vs R_s %.4f (%d FAIL, %d UNRESOLVED%s); "
                  "B2 (sampled, step %.3f mm) %s: %ld samples, %d overhang, %d budget, %d normal-turn, %d depth-band",
                  toString(out.b1), edges, out.minDistance, Rs, out.b1Fail, out.b1Unresolved,
                  counts[2] ? (", " + std::to_string(counts[2]) + " refinement stacks full").c_str() : "", step, toString(out.b2),
                  out.b2Samples, overhang, out.b2Unresolved[1], out.b2Unresolved[2], out.b2Unresolved[4]);
    out.detail = b;
}

} // namespace

bool checkBoundary(const BoundaryRequest& request, BoundaryCheck& out, std::string* err) try {
    const CurveSet& cs = request.curves;
    const SurfaceBinding& bind = request.binding;
    const MeshSurface& surf = request.surface;
    const double Rs = request.Rs, step = request.step;
    out = BoundaryCheck{};
    if (!available(err)) return false;
    Reservation lease(Reservation::exclusiveBytes());  // runs alone on the device (runtime.h)
    BoxArrays boxes;
    const std::vector<BoundarySegment> segs = boundarySegments(surf, boxes);
    const double reach = Rs + 4 * step;
    boxes.cell = float(std::max(reach, 1e-6));
    DeviceGrid grid;
    DeviceBuffer dSegs, dYarns, us, xs, reductions;
    grid.upload(boxes);
    dSegs.upload(segs.data(), segs.size() * sizeof(segs[0]));
    DeviceCurves curves;
    curves.upload(flattenCurves(cs, ArcBounds::Included));
    std::vector<BoundaryYarn> yarns;
    int ns = 0;
    for (size_t y = 0; y < cs.yarns.size(); ++y) {
        const double L = cs.yarns[y].length();
        const int n = std::max(1, static_cast<int>(std::ceil(L / step)));
        yarns.push_back({static_cast<int>(y), ns, n, L});
        ns += n + 1;
    }
    dYarns.upload(yarns.data(), yarns.size() * sizeof(yarns[0]));
    us.alloc(8ull * ns + 8);
    xs.alloc(12ull * ns + 12);
    const BoundarySamples samples{static_cast<BoundaryYarn*>(dYarns.data()), static_cast<int>(yarns.size()), ns, us.doubles(), xs.floats()};
    boundarySampleKernel<<<blocksFor(ns), kThreadsPerBlock>>>(curves.view, samples);
    // On the device: the distance key, then fail, unresolved, stackFull and the five B2 results.
    reductions.alloc(64);
    const unsigned long long start[2] = {~0ull, 0ull};
    check(cudaMemcpy(reductions.data(), start, sizeof start, cudaMemcpyHostToDevice), "GPU boundary check");
    int* counters = reinterpret_cast<int*>(reductions.ulongs() + 1);
    check(cudaMemset(counters, 0, 32), "GPU boundary check");
    boundaryDistanceKernel<<<blocksFor(ns), kThreadsPerBlock>>>(curves.view, samples, {grid.view, static_cast<BoundarySegment*>(dSegs.data()), float(reach), float(Rs), float(request.minStep)},
                                                                {reductions.ulongs(), counters, counters + 1, counters + 2});
    checkLaunch("GPU boundary B1");
    if (bind.surface) out.b2Samples = boundaryContinuations(request, yarns, samples, counters + 3);
    unsigned long long distanceKey = 0;
    int counts[8];
    reductions.download(&distanceKey, 8);
    reductions.download(counts, sizeof counts, 8);
    if (distanceKey != ~0ull) {
        const unsigned bits = static_cast<unsigned>(distanceKey);
        float d;
        std::memcpy(&d, &bits, 4);
        out.minDistance = d;
    }
    writeVerdicts(counts, Rs, step, segs.size(), out);
    return true;
} catch (const DeviceError& e) {
    return reportFailure(e, err);
}

} // namespace knit::gpu
