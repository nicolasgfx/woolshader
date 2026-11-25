// GPU backend: the tiled lattice's node classification (chart_lattice.cpp tileMesh) on the device:
// per candidate node (one thread) the chart location in uv' (UvAtlas::locate), the inside test
// against the other eligible components (MeshQuery::inside), the footprint test (clear: the
// minimum s2 under the rectangle, iterated, then the chart's boundary edges against it) and, for
// kept nodes, the two run-end footprint tests. Single precision; tolerances below.
#include <algorithm>
#include <cmath>
#include <memory>
#include <mutex>
#include <unordered_map>

#include "model/mesh_query.h"
#include "gpu/geometry/surface.cuh"
#include "gpu/geometry/device_geometry.h"
#include "knitting/layout/uv_atlas.h"
#include "knitting/layout/chart_lattice.h"
#include "gpu/device.h"
#include "gpu/surface/placement.h"

namespace knit::gpu {

namespace {

struct LatticeChartView {
    float glo[2];
    float cell;
    int nx, ny;
    const unsigned* cellStart;  // ny * nx + 1
    const int* cellTris;        // chart-local triangle indices q
    const int* tris;            // q -> mesh triangle
    const float* uvp;           // 6 per q (uv' corners)
    const unsigned char* usable;  // per q
    const float* s2;            // per q (usable only)
    const float* bnd;           // 4 per boundary edge (uv' a, b)
    int nBnd;
    const int* meshTri;         // 3 per mesh triangle
    const float* meshPos;       // 3 per vertex
};

struct LatticeNodeQuery {
    float c2[2], v2[2];
    float xmin, ymin, w, h;
    int ni, nr;
    float margin, Ymin, Ymax, Xpad;
    int footprint;
    int nSolids;
    float offset[3];
};

__device__ __forceinline__ float cross2(float ax, float ay, float bx, float by) { return ax * by - ay * bx; }

// chart_lattice.cpp segmentMeetsRect (Liang-Barsky).
struct UvRect {
    float x0, x1, y0, y1;
};
__device__ bool segmentMeetsRect(float2 a, float2 b, const UvRect& rect) {
    float t0 = 0, t1 = 1;
    const float dx = b.x - a.x, dy = b.y - a.y;
    auto clip = [&](float p, float q) {
        if (p == 0) return q >= 0;
        const float r = q / p;
        if (p < 0) {
            if (r > t1) return false;
            t0 = fmaxf_(t0, r);
        } else {
            if (r < t0) return false;
            t1 = fminf_(t1, r);
        }
        return true;
    };
    return clip(-dx, a.x - rect.x0) && clip(dx, rect.x1 - a.x) && clip(-dy, a.y - rect.y0) && clip(dy, rect.y1 - a.y) && t0 <= t1;
}

// UvAtlas::minS2In (3e38: none).
__device__ float minS2In(const LatticeChartView& c, float lox, float loy, float hix, float hiy) {
    int x0 = static_cast<int>(floorf((lox - c.glo[0]) / c.cell)), x1 = static_cast<int>(floorf((hix - c.glo[0]) / c.cell));
    int y0 = static_cast<int>(floorf((loy - c.glo[1]) / c.cell)), y1 = static_cast<int>(floorf((hiy - c.glo[1]) / c.cell));
    x0 = max(0, min(c.nx - 1, x0)); x1 = max(0, min(c.nx - 1, x1));
    y0 = max(0, min(c.ny - 1, y0)); y1 = max(0, min(c.ny - 1, y1));
    float s = 3.0e38f;
    for (int y = y0; y <= y1; ++y)
        for (int x = x0; x <= x1; ++x) {
            const int cell = y * c.nx + x;
            for (unsigned at = c.cellStart[cell]; at < c.cellStart[cell + 1]; ++at) {
                const int q = c.cellTris[at];
                if (c.usable[q]) s = fminf_(s, c.s2[q]);
            }
        }
    return s;
}

// tileMesh's clear(x, y, exL, exR).
__device__ bool clearAt(const LatticeChartView& c, const LatticeNodeQuery& Q, float x, float y, float exL, float exR) {
    if (!Q.footprint) return true;
    float mU = Q.margin;
    const float h = Q.h;
    for (int it = 0; it < 2; ++it) {
        const float x0 = x - exL - mU, x1 = x + exR + mU, y0 = y + Q.Ymin * h - mU, y1 = y + Q.Ymax * h + mU;
        float lox = 3.0e38f, loy = 3.0e38f, hix = -3.0e38f, hiy = -3.0e38f;
        for (float cx : {x0, x1})
            for (float cy : {y0, y1}) {
                const float qx = Q.c2[0] * cx + Q.v2[0] * cy, qy = Q.c2[1] * cx + Q.v2[1] * cy;
                lox = fminf_(lox, qx); loy = fminf_(loy, qy);
                hix = fmaxf_(hix, qx); hiy = fmaxf_(hiy, qy);
            }
        const float s2 = minS2In(c, lox, loy, hix, hiy);
        const float next = Q.margin / fmaxf_(1e-3f, fminf_(1.0f, s2));
        if (!(next > mU)) break;
        mU = next;
    }
    const float x0 = x - exL - mU, x1 = x + exR + mU, y0 = y + Q.Ymin * h - mU, y1 = y + Q.Ymax * h + mU;
    for (int e = 0; e < c.nBnd; ++e) {
        const float* b = c.bnd + 4 * e;
        const float ax = b[0] * Q.c2[0] + b[1] * Q.c2[1], ay = b[0] * Q.v2[0] + b[1] * Q.v2[1];
        const float bx = b[2] * Q.c2[0] + b[3] * Q.c2[1], by = b[2] * Q.v2[0] + b[3] * Q.v2[1];
        if (segmentMeetsRect(make_float2(ax, ay), make_float2(bx, by), {x0, x1, y0, y1})) return false;
    }
    return true;
}

struct NodeOut {
    int code, tri;
    float bary[3], p[3];
    signed char clearStart, clearEnd;
};

__global__ void latticeKernel(LatticeChartView c, LatticeNodeQuery Q, const MeshView* solids, NodeOut* out) {
    const int g = threadItem();
    if (g >= Q.ni * Q.nr) return;
    const int r = g / Q.ni, i = g % Q.ni;
    NodeOut o{};
    o.code = 1;
    o.tri = -1;
    o.clearStart = o.clearEnd = -1;
    const float x = Q.xmin + i * Q.w, y = Q.ymin + r * Q.h;
    const float px = Q.c2[0] * x + Q.v2[0] * y, py = Q.c2[1] * x + Q.v2[1] * y;
    // UvAtlas::locate (barycentric tolerance: CPU -1e-12, float -1e-6).
    const int gx = static_cast<int>(floorf((px - c.glo[0]) / c.cell)), gy = static_cast<int>(floorf((py - c.glo[1]) / c.cell));
    int hits[8], nHits = 0;
    float bary[3] = {0, 0, 0};
    if (gx >= 0 && gy >= 0 && gx < c.nx && gy < c.ny) {
        const int cell = gy * c.nx + gx;
        for (unsigned at = c.cellStart[cell]; at < c.cellStart[cell + 1]; ++at) {
            const int q = c.cellTris[at];
            if (!c.usable[q]) continue;
            const float* P = c.uvp + 6 * q;
            const float A = cross2(P[2] - P[0], P[3] - P[1], P[4] - P[0], P[5] - P[1]);
            if (A == 0) continue;
            const float b1 = cross2(px - P[0], py - P[1], P[4] - P[0], P[5] - P[1]) / A;
            const float b2 = cross2(P[2] - P[0], P[3] - P[1], px - P[0], py - P[1]) / A;
            const float b0 = 1 - b1 - b2;
            const float eps = -1e-6f;
            if (b0 >= eps && b1 >= eps && b2 >= eps) {
                if (nHits == 0) {
                    bary[0] = fmaxf_(0.0f, b0);
                    bary[1] = fmaxf_(0.0f, b1);
                    bary[2] = fmaxf_(0.0f, b2);
                    const float s = bary[0] + bary[1] + bary[2];
                    for (int k = 0; k < 3; ++k) bary[k] /= s;
                }
                if (nHits < 8) hits[nHits] = c.tris[q];
                ++nHits;
            }
        }
    }
    if (nHits > 1 && nHits <= 8) {  // one location on a shared edge / vertex
        bool allAdjacent = true;
        for (int h = 1; h < nHits && allAdjacent; ++h) {
            bool shares = false;
            for (int u = 0; u < 3; ++u)
                for (int v = 0; v < 3; ++v) shares = shares || c.meshTri[3 * hits[0] + u] == c.meshTri[3 * hits[h] + v];
            allAdjacent = shares;
        }
        if (allAdjacent) nHits = 1;
    }
    if (nHits == 0) { out[g] = o; return; }
    if (nHits > 1) { o.code = 2; out[g] = o; return; }
    o.tri = hits[0];
    const int* tv = c.meshTri + 3 * o.tri;
    const f3 p = fromF(c.meshPos + 3 * tv[0]) * bary[0] + fromF(c.meshPos + 3 * tv[1]) * bary[1] + fromF(c.meshPos + 3 * tv[2]) * bary[2];
    for (int k = 0; k < 3; ++k) o.bary[k] = bary[k];
    o.p[0] = p.x; o.p[1] = p.y; o.p[2] = p.z;
    const f3 q = p + fromF(Q.offset);
    for (int s = 0; s < Q.nSolids; ++s)
        if (meshInside(solids[s], q)) { o.code = 3; out[g] = o; return; }
    if (!clearAt(c, Q, x, y, 0.5f * Q.w, 0.5f * Q.w)) { o.code = 4; out[g] = o; return; }
    o.code = 0;
    o.clearStart = clearAt(c, Q, x, y, Q.Xpad * Q.w, 0.5f * Q.w) ? 1 : 0;
    o.clearEnd = clearAt(c, Q, x, y, 0.5f * Q.w, Q.Xpad * Q.w) ? 1 : 0;
    out[g] = o;
}

// The chart of a lattice on the device.
struct DeviceLatticeChart {
    DeviceBuffer start, cellTris, tris, uvp, usable, s2, boundary, meshTri, meshPos;
    LatticeChartView view{};

    DeviceLatticeChart(const UvAtlas& A, const LatticeChartQuery& lq) {
        const UvChart& chart = A.charts()[lq.chart];
        const UvAtlas::GridView grid = A.gridView(lq.chart);
        const TriMesh& m = A.mesh();
        std::vector<unsigned> hostStart{0};
        std::vector<int> hostCellTris;
        for (const auto& cell : *grid.cells) {
            hostCellTris.insert(hostCellTris.end(), cell.begin(), cell.end());
            hostStart.push_back(static_cast<unsigned>(hostCellTris.size()));
        }
        std::vector<float> hostUvp, hostS2, hostBoundary;
        std::vector<unsigned char> hostUsable;
        for (int t : chart.tris) {
            const std::array<Vec2, 3> P = A.uvp(t);
            for (const Vec2& v : P) hostUvp.insert(hostUvp.end(), {float(v.x), float(v.y)});
            hostUsable.push_back(A.usable(t) ? 1 : 0);
            hostS2.push_back(A.usable(t) ? float(A.s2()[t]) : 0.0f);
        }
        for (const auto& e : A.boundary(lq.chart)) hostBoundary.insert(hostBoundary.end(), {float(e[0].x), float(e[0].y), float(e[1].x), float(e[1].y)});
        std::vector<int> hostMeshTri;
        for (const auto& t : m.tri) hostMeshTri.insert(hostMeshTri.end(), {t[0], t[1], t[2]});
        start.upload(hostStart);
        cellTris.upload(hostCellTris);
        tris.upload(std::vector<int>(chart.tris.begin(), chart.tris.end()));
        uvp.upload(hostUvp);
        usable.upload(hostUsable);
        s2.upload(hostS2);
        boundary.upload(hostBoundary);
        meshTri.upload(hostMeshTri);
        std::vector<float> hostMeshPos;
        appendFloat3(m.pos, hostMeshPos);
        meshPos.upload(hostMeshPos);
        view.glo[0] = float(grid.lo.x);
        view.glo[1] = float(grid.lo.y);
        view.cell = float(grid.cell);
        view.nx = grid.nx;
        view.ny = grid.ny;
        view.cellStart = start.uints();
        view.cellTris = cellTris.ints();
        view.tris = tris.ints();
        view.uvp = uvp.floats();
        view.usable = usable.chars();
        view.s2 = s2.floats();
        view.bnd = boundary.floats();
        view.nBnd = static_cast<int>(hostBoundary.size() / 4);
        view.meshTri = meshTri.ints();
        view.meshPos = meshPos.floats();
    }
};

LatticeNodeQuery latticeQuery(const LatticeChartQuery& lq, int nSolids, const Vec3& offset) {
    LatticeNodeQuery Q{};
    Q.c2[0] = float(lq.c2.x); Q.c2[1] = float(lq.c2.y);
    Q.v2[0] = float(lq.v2.x); Q.v2[1] = float(lq.v2.y);
    Q.xmin = float(lq.xmin);
    Q.ymin = float(lq.ymin);
    Q.w = float(lq.w);
    Q.h = float(lq.h);
    Q.ni = lq.ni;
    Q.nr = lq.nr;
    Q.margin = float(lq.margin);
    Q.Ymin = float(lq.Ymin);
    Q.Ymax = float(lq.Ymax);
    Q.Xpad = float(lq.Xpad);
    Q.footprint = lq.footprintCheck ? 1 : 0;
    Q.nSolids = nSolids;
    Q.offset[0] = float(offset.x); Q.offset[1] = float(offset.y); Q.offset[2] = float(offset.z);
    return Q;
}

} // namespace

bool classifyLattice(const LatticeRequest& request, std::vector<LatticeNodeClass>& out, std::string* err) try {
    const LatticeChartQuery& lq = request.query;
    if (!available(err)) return false;
    Reservation lease(Reservation::exclusiveBytes());  // runs alone on the device (runtime.h)
    const DeviceLatticeChart chart(request.atlas, lq);
    std::vector<std::shared_ptr<DeviceMesh>> keep;
    std::vector<MeshView> solidViews;
    for (const MeshQuery* mq : request.solids) {
        keep.push_back(residentMesh(*mq));
        solidViews.push_back(keep.back()->view);
    }
    const int n = lq.ni * lq.nr;
    DeviceBuffer dSolids, dOut;
    dSolids.upload(solidViews.data(), solidViews.size() * sizeof(solidViews[0]));
    dOut.alloc(sizeof(NodeOut) * std::max(n, 1));
    if (n > 0)
        latticeKernel<<<blocksFor(n), kThreadsPerBlock>>>(chart.view, latticeQuery(lq, static_cast<int>(solidViews.size()), request.offset), static_cast<MeshView*>(dSolids.data()),
                                                          static_cast<NodeOut*>(dOut.data()));
    std::vector<NodeOut> nodes(n);
    checkLaunch("GPU lattice");
    dOut.download(nodes.data(), nodes.size() * sizeof(nodes[0]));
    out.assign(n, LatticeNodeClass{});
    for (int g = 0; g < n; ++g) {
        LatticeNodeClass& o = out[g];
        o.code = nodes[g].code;
        o.tri = nodes[g].tri;
        for (int k = 0; k < 3; ++k) o.bary[k] = nodes[g].bary[k];
        o.p = {nodes[g].p[0], nodes[g].p[1], nodes[g].p[2]};
        o.clearStart = nodes[g].clearStart;
        o.clearEnd = nodes[g].clearEnd;
    }
    return true;
} catch (const DeviceError& e) {
    return reportFailure(e, err);
}

} // namespace knit::gpu
