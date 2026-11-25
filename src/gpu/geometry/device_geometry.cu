// GPU backend: from knit types to the device (device_geometry.h). Top to bottom: flatten on the host,
// build a grid on the device, upload, resident meshes.
#include "gpu/geometry/device_geometry.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <map>
#include <vector>

#include <cuda_runtime.h>

#include "gpu/device/kernel.cuh"
#include "gpu/geometry/math.cuh"
#include "model/hash.h"
#include "model/mesh_query.h"
#include "model/parallel.h"

namespace knit::gpu {

// --- flatten

CurveArrays flattenCurves(const CurveSet& cs, ArcBounds bounds) {
    CurveArrays out;
    size_t controls = 0;
    for (const YarnPath& y : cs.yarns) controls += y.curve.ctrl.size();
    out.ctrl.reserve(3 * controls);
    int spans = 0;
    for (const YarnPath& y : cs.yarns) {
        CurveRecord c;
        c.ctrlOff = static_cast<int>(out.ctrl.size() / 3);
        c.nCtrl = static_cast<int>(y.curve.ctrl.size());
        appendFloat3(y.curve.ctrl, out.ctrl);
        out.thickness.insert(out.thickness.end(), y.thickness.begin(), y.thickness.end());
        const YarnPath::Tables& t = y.tables();
        c.tabOff = static_cast<int>(out.tu.size());
        c.nTab = static_cast<int>(t.tu.size());
        out.tu.insert(out.tu.end(), t.tu.begin(), t.tu.end());
        out.ts.insert(out.ts.end(), t.ts.begin(), t.ts.end());
        c.length = y.length();
        c.spanOff = spans;
        spans += y.curve.numSpans();
        c.preOff = -1;
        if (bounds == ArcBounds::Included && static_cast<int>(t.preLo.size()) == y.curve.numSpans() + 1) {
            c.preOff = static_cast<int>(out.preLo.size());
            out.preLo.insert(out.preLo.end(), t.preLo.begin(), t.preLo.end());
            out.preHi.insert(out.preHi.end(), t.preHi.begin(), t.preHi.end());
        }
        out.curves.push_back(c);
    }
    return out;
}

SurfaceArrays flattenSurface(const SurfaceModel& sm) {
    // Per yarn in parallel into its own slots (offsets first): the arrays are those of a serial fill.
    SurfaceArrays out;
    const size_t ny = sm.yarns.size();
    out.yarns.resize(ny);
    size_t total = 0;
    for (size_t k = 0; k < ny; ++k) {
        out.yarns[k].ringOff = static_cast<int>(total);
        total += sm.yarns[k].rings.size();
    }
    out.rings.resize(total);
    parallelFor(static_cast<int64_t>(ny), 1, [&](int64_t b0, int64_t e0, int64_t) {
        for (int64_t k = b0; k < e0; ++k) {
            const YarnSurface& ys = sm.yarns[size_t(k)];
            YarnSurfaceRecord& f = out.yarns[size_t(k)];
            f.nRings = static_cast<int>(ys.rings.size());
            f.ds = ys.ds;
            f.length = ys.length;
            f.radius = float(ys.radius);
            double ext = 0;
            RingRecord* dst = out.rings.data() + f.ringOff;
            for (const Ring& R : ys.rings) {
                ext = std::max({ext, R.a, R.b});
                RingRecord r{};
                r.s = R.s;
                r.u = R.u;
                const Vec3* v[4] = {&R.c, &R.T, &R.n1, &R.n2};
                float* d[4] = {r.c, r.T, r.n1, r.n2};
                for (int q = 0; q < 4; ++q) {
                    d[q][0] = float(v[q]->x);
                    d[q][1] = float(v[q]->y);
                    d[q][2] = float(v[q]->z);
                }
                r.a = float(R.a);
                r.b = float(R.b);
                r.phi = float(R.phi);
                *dst++ = r;
            }
            f.ext = float(ext);
        }
    }, ExecPolicy{true});
    return out;
}

namespace {
// Boxes rounded outward to float (a conservative grid filter).
void putBox(std::vector<float>& lo, std::vector<float>& hi, const Vec3& a, const Vec3& b) {
    auto down = [](double v) {
        float f = float(v);
        return double(f) > v ? std::nextafter(f, -INFINITY) : f;
    };
    auto up = [](double v) {
        float f = float(v);
        return double(f) < v ? std::nextafter(f, INFINITY) : f;
    };
    lo.insert(lo.end(), {down(a.x), down(a.y), down(a.z)});
    hi.insert(hi.end(), {up(b.x), up(b.y), up(b.z)});
}
} // namespace

// SurfaceIndex's capsules (yarn_surface.cpp) and their boxes.
IndexArrays flattenIndex(const SurfaceModel& sm) {
    // Per yarn in parallel into its own slots (capsule offsets first); zTop and maxR are maxima, zBottom a minimum, so
    // the result is that of the serial loop.
    IndexArrays out;
    const size_t ny = sm.yarns.size();
    std::vector<size_t> first(ny + 1, 0);
    auto stepOf = [](const YarnSurface& ys) { return std::max(1, static_cast<int>(std::round(0.1 / ys.ds))); };
    for (size_t y = 0; y < ny; ++y) {
        const YarnSurface& ys = sm.yarns[y];
        const int N = static_cast<int>(ys.rings.size()) - 1, step = stepOf(ys);
        first[y + 1] = first[y] + (N > 0 ? size_t((N + step - 1) / step) : 0);
    }
    const size_t nc = first[ny];
    out.caps.resize(nc);
    out.boxes.lo.resize(3 * nc);
    out.boxes.hi.resize(3 * nc);
    std::vector<double> zTopY(ny, -1e300), zBottomY(ny, 1e300), maxRY(ny, 0);
    parallelFor(static_cast<int64_t>(ny), 1, [&](int64_t b0, int64_t e0, int64_t) {
        for (int64_t y = b0; y < e0; ++y) {
            const YarnSurface& ys = sm.yarns[size_t(y)];
            const int N = static_cast<int>(ys.rings.size()) - 1, step = stepOf(ys);
            size_t q = first[size_t(y)];
            for (int k = 0; k < N; k += step, ++q) {
                const int k1 = std::min(N, k + step);
                const Vec3 a = ys.rings[k].c, b = ys.rings[k1].c;
                double R = 0, dev = 0;
                for (int j = k; j <= k1; ++j) {
                    R = std::max({R, ys.rings[j].a, ys.rings[j].b});
                    dev = std::max(dev, closestSegmentPoints(ys.rings[j].c, ys.rings[j].c, a, b).distance);
                }
                R += dev + 1e-3;
                CapsuleRecord c{};
                c.a[0] = float(a.x); c.a[1] = float(a.y); c.a[2] = float(a.z);
                c.b[0] = float(b.x); c.b[1] = float(b.y); c.b[2] = float(b.z);
                c.R = float(R);
                c.yarn = ys.yarnId;
                c.s0 = ys.rings[k].s;
                c.s1 = ys.rings[k1].s;
                out.caps[q] = c;
                zTopY[size_t(y)] = std::max({zTopY[size_t(y)], a.z + R, b.z + R});
                zBottomY[size_t(y)] = std::min({zBottomY[size_t(y)], a.z - R, b.z - R});
                maxRY[size_t(y)] = std::max(maxRY[size_t(y)], R);
                std::vector<float> l1, h1;
                putBox(l1, h1, Vec3{std::min(a.x, b.x) - R, std::min(a.y, b.y) - R, std::min(a.z, b.z) - R},
                       Vec3{std::max(a.x, b.x) + R, std::max(a.y, b.y) + R, std::max(a.z, b.z) + R});
                for (int d = 0; d < 3; ++d) {
                    out.boxes.lo[3 * q + d] = l1[size_t(d)];
                    out.boxes.hi[3 * q + d] = h1[size_t(d)];
                }
            }
        }
    }, ExecPolicy{true});
    double zTop = -1e300, zBottom = 1e300, maxR = 0;
    for (size_t y = 0; y < ny; ++y) {
        zTop = std::max(zTop, zTopY[y]);
        zBottom = std::min(zBottom, zBottomY[y]);
        maxR = std::max(maxR, maxRY[y]);
    }
    out.zTop = float(zTop);
    out.zBottom = float(zBottom);
    out.boxes.cell = float(1.0 + maxR);
    return out;
}

namespace {
// The triangles' boxes and the mesh's bounds.
void addBoxes(const TriMesh& m, MeshArrays& out) {
    double edges = 0;
    for (const auto& t : m.tri) {
        const Vec3 &a = m.pos[t[0]], &b = m.pos[t[1]], &c = m.pos[t[2]];
        for (int k = 0; k < 3; ++k) edges += norm(m.pos[t[(k + 1) % 3]] - m.pos[t[k]]);
        putBox(out.boxes.lo, out.boxes.hi, Vec3{std::min({a.x, b.x, c.x}), std::min({a.y, b.y, c.y}), std::min({a.z, b.z, c.z})},
               Vec3{std::max({a.x, b.x, c.x}), std::max({a.y, b.y, c.y}), std::max({a.z, b.z, c.z})});
    }
    const double mean = m.tri.empty() ? 1.0 : edges / (3.0 * m.tri.size());
    out.boxes.cell = float(2 * mean);
    Vec3 blo, bhi;
    m.bounds(blo, bhi);
    out.lo[0] = float(blo.x); out.lo[1] = float(blo.y); out.lo[2] = float(blo.z);
    out.hi[0] = float(bhi.x); out.hi[1] = float(bhi.y); out.hi[2] = float(bhi.z);
}
} // namespace

MeshArrays flattenMesh(const TriMesh& m) {
    MeshArrays out;
    appendFloat3(m.pos, out.pos);
    for (const auto& t : m.tri) out.tri.insert(out.tri.end(), {t[0], t[1], t[2]});
    TriMesh copy;
    const std::vector<Vec3>* vn = &m.vnormal;
    if (m.vnormal.size() != m.pos.size()) {
        copy = m;
        copy.computeNormals();
        vn = &copy.vnormal;
    }
    appendFloat3(*vn, out.vnormal);
    for (size_t f = 0; f < m.tri.size(); ++f) {
        const Vec3 n = m.faceNormal(static_cast<int>(f));
        out.faceN.insert(out.faceN.end(), {float(n.x), float(n.y), float(n.z)});
    }
    addBoxes(m, out);
    return out;
}

MeshArrays flattenMesh(const MeshQuery& mq) {
    const TriMesh& m = mq.mesh();
    MeshArrays out;
    appendFloat3(m.pos, out.pos);
    for (const auto& t : m.tri) out.tri.insert(out.tri.end(), {t[0], t[1], t[2]});
    appendFloat3(mq.faceNormals(), out.faceN);
    appendFloat3(mq.cornerEdgeNormals(), out.edgeN);
    appendFloat3(mq.vertexNormals(), out.pvertN);
    addBoxes(m, out);
    return out;
}

// --- memory

size_t CurveArrays::bytes() const {
    return curves.size() * sizeof(CurveRecord) + (ctrl.size() + thickness.size()) * 4 + (tu.size() + ts.size() + preLo.size() + preHi.size()) * 8;
}
size_t SurfaceArrays::bytes() const { return yarns.size() * sizeof(YarnSurfaceRecord) + rings.size() * sizeof(RingRecord); }
size_t IndexArrays::bytes() const { return caps.size() * sizeof(CapsuleRecord) + boxes.gridBytes(); }
size_t MeshArrays::bytes() const {
    return (pos.size() + vnormal.size() + faceN.size() + pvertN.size() + edgeN.size() + tri.size() + nbr.size() + vtStart.size() + vtList.size()) * 4 +
           boxes.gridBytes();
}

// --- the grid

namespace {

// The float expression of GridCursor, so item and query cells agree bit for bit.
__host__ __device__ __forceinline__ long long cellOf(float v, float inv) { return static_cast<long long>(floorf(v * inv)); }

struct GridBoxes {
    const float* lo;
    const float* hi;
    int count;
    __host__ __device__ bool leftOut(int i) const { return lo[3 * i] > hi[3 * i]; }
};

struct GridShape {
    float inv;              // 1 / cell
    long long lo[3], n[3];  // the first cell and the number of cells per axis
};

__global__ void boxCellRangeKernel(GridBoxes b, float inv, long long* lowest, long long* highest) {
    const int i = threadItem();
    if (i >= b.count || b.leftOut(i)) return;
    for (int a = 0; a < 3; ++a) {
        atomicMin(&lowest[a], cellOf(b.lo[3 * i + a], inv));
        atomicMax(&highest[a], cellOf(b.hi[3 * i + a], inv));
    }
}

// The cells of item i, relative to the grid's first cell.
__device__ __forceinline__ void itemCells(const GridBoxes& b, int i, const GridShape& g, long long c0[3], long long c1[3]) {
    for (int a = 0; a < 3; ++a) {
        c0[a] = cellOf(b.lo[3 * i + a], g.inv) - g.lo[a];
        c1[a] = cellOf(b.hi[3 * i + a], g.inv) - g.lo[a];
    }
}

__global__ void countRegistrationsKernel(GridBoxes b, GridShape g, unsigned long long* count) {
    const int i = threadItem();
    if (i >= b.count) return;
    if (b.leftOut(i)) {
        count[i] = 0;
        return;
    }
    long long c0[3], c1[3];
    itemCells(b, i, g, c0, c1);
    count[i] = static_cast<unsigned long long>((c1[0] - c0[0] + 1) * (c1[1] - c0[1] + 1) * (c1[2] - c0[2] + 1));
}

// One key per cell of the item: (cell << 32) | item, so sorting orders the items within a cell.
__global__ void emitRegistrationsKernel(GridBoxes b, GridShape g, const unsigned long long* offsets, unsigned long long* keys, unsigned* cellCounts) {
    const int i = threadItem();
    if (i >= b.count || b.leftOut(i)) return;
    long long c0[3], c1[3];
    itemCells(b, i, g, c0, c1);
    unsigned long long at = offsets[i];
    for (long long x = c0[0]; x <= c1[0]; ++x)
        for (long long y = c0[1]; y <= c1[1]; ++y)
            for (long long z = c0[2]; z <= c1[2]; ++z) {
                const unsigned cell = static_cast<unsigned>((x * g.n[1] + y) * g.n[2] + z);
                keys[at++] = (static_cast<unsigned long long>(cell) << 32) | static_cast<unsigned>(i);
                atomicAdd(&cellCounts[cell], 1u);
            }
}

// A sorted key becomes an entry: the item with its "lowest cell" flags.
__global__ void gridEntriesKernel(GridBoxes b, GridShape g, const unsigned long long* keys, unsigned registrations, unsigned* entries) {
    const unsigned p = threadItem();
    if (p >= registrations) return;
    const long long cell = static_cast<long long>(keys[p] >> 32);
    const int i = static_cast<int>(keys[p] & 0xffffffffu);
    const long long z = cell % g.n[2], y = (cell / g.n[2]) % g.n[1], x = cell / (g.n[2] * g.n[1]);
    long long c0[3], c1[3];
    itemCells(b, i, g, c0, c1);
    const unsigned flags = (x == c0[0] ? 1u : 0u) | (y == c0[1] ? 2u : 0u) | (z == c0[2] ? 4u : 0u);
    entries[p] = static_cast<unsigned>(i) | flags << 29;
}

// The occupied cell range. The cell grows until the range has at most limits.maxCells cells; false
// when every item is left out.
bool fitCells(const GridBoxes& boxes, const GridLimits& limits, DeviceBuffer& range, float& cell, GridShape& g) {
    range.reserve(6 * sizeof(long long));
    for (;;) {
        g.inv = 1.0f / cell;
        long long ends[6] = {LLONG_MAX, LLONG_MAX, LLONG_MAX, LLONG_MIN, LLONG_MIN, LLONG_MIN};
        check(cudaMemcpy(range.data(), ends, sizeof ends, cudaMemcpyHostToDevice), "grid range");
        boxCellRangeKernel<<<blocksFor(boxes.count), kThreadsPerBlock>>>(boxes, g.inv, range.longs(), range.longs() + 3);
        range.download(ends, sizeof ends);
        if (ends[0] > ends[3]) return false;
        long long cells = 1;
        for (int a = 0; a < 3; ++a) {
            g.lo[a] = ends[a];
            g.n[a] = ends[3 + a] - ends[a] + 1;
            cells *= g.n[a];
        }
        if (cells <= limits.maxCells) return true;
        cell *= std::cbrt(static_cast<float>(cells) / static_cast<float>(limits.maxCells)) * 1.01f;
    }
}

} // namespace

size_t BoxArrays::gridBytes() const {
    const float inv = 1.0f / cell;
    long long lowest[3] = {LLONG_MAX, LLONG_MAX, LLONG_MAX}, highest[3] = {LLONG_MIN, LLONG_MIN, LLONG_MIN};
    size_t registrations = 0;
    for (int i = 0; i < count(); ++i) {
        if (lo[3 * i] > hi[3 * i]) continue;
        size_t cellsOfItem = 1;
        for (int a = 0; a < 3; ++a) {
            const long long c0 = cellOf(lo[3 * i + a], inv), c1 = cellOf(hi[3 * i + a], inv);
            lowest[a] = std::min(lowest[a], c0);
            highest[a] = std::max(highest[a], c1);
            cellsOfItem *= static_cast<size_t>(c1 - c0 + 1);
        }
        registrations += cellsOfItem;
    }
    long long cells = 0;
    if (registrations) cells = std::min((highest[0] - lowest[0] + 1) * (highest[1] - lowest[1] + 1) * (highest[2] - lowest[2] + 1), GridLimits{}.maxCells);
    // Per item: two boxes, count, offset. Per cell: start, count. Per registration: entry, two keys,
    // the sort's scratch. Buffers keep a quarter in reserve (DeviceBuffer::reserve).
    const size_t bytes = size_t(count() + 1) * (24 + 16) + size_t(cells + 1) * 8 + registrations * (4 + 16 + 8);
    return bytes + bytes / 4;
}

void DeviceGrid::upload(const BoxArrays& boxes) {
    boxLo.upload(boxes.lo);
    boxHi.upload(boxes.hi);
    build(boxes.count(), boxes.cell);
}

void DeviceGrid::build(int items, float cell) {
    DevicePrimitives prims;
    GridScratch scratch;
    build(items, cell, prims, scratch);
}

void DeviceGrid::build(int items, float cell, DevicePrimitives& prims, GridScratch& s, const GridLimits& limits) {
    if (items >= (1 << 29)) throw DeviceError("device grid: too many items");
    if (!(cell > 0)) cell = 1;
    view = GridView{};
    view.cell = cell;
    const GridBoxes boxes{boxLo.floats(), boxHi.floats(), items};
    GridShape g{};
    if (items == 0 || !fitCells(boxes, limits, s.range, cell, g)) return;

    // How many (cell, item) registrations? The total is checked in 64 bits before it is allocated.
    const long long nCells = g.n[0] * g.n[1] * g.n[2];
    s.counts.reserve(size_t(items + 1) * 8);
    s.offsets.reserve(size_t(items + 1) * 8);
    countRegistrationsKernel<<<blocksFor(items), kThreadsPerBlock>>>(boxes, g, s.counts.ulongs());
    const unsigned long long total = prims.countsToOffsets(s.counts.ulongs(), s.offsets.ulongs(), items);
    if (total >= static_cast<unsigned long long>(limits.maxEntries))
        throw DeviceError("device grid: capacity exceeded: cell registrations " + std::to_string(total) + " >= " + std::to_string(limits.maxEntries));
    const unsigned registrations = static_cast<unsigned>(total);

    // Emit, sort, and the sorted registrations are the buckets.
    s.keys.reserve(size_t(registrations) * 8 + 8);
    s.sorted.reserve(size_t(registrations) * 8 + 8);
    s.cellCounts.reserve(size_t(nCells + 1) * 4);
    entries.reserve(size_t(registrations) * 4 + 4);
    cellStart.reserve(size_t(nCells + 1) * 4);
    check(cudaMemset(s.cellCounts.data(), 0, size_t(nCells + 1) * 4), "grid counts");
    emitRegistrationsKernel<<<blocksFor(items), kThreadsPerBlock>>>(boxes, g, s.offsets.ulongs(), s.keys.ulongs(), s.cellCounts.uints());
    int endBit = 32;
    while (endBit < 64 && (1ull << (endBit - 32)) <= static_cast<unsigned long long>(nCells)) ++endBit;
    prims.sortKeys({s.keys.ulongs(), s.sorted.ulongs(), static_cast<long long>(registrations), endBit});
    prims.countsToOffsets(s.cellCounts.uints(), cellStart.uints(), nCells);
    if (registrations) gridEntriesKernel<<<blocksFor(registrations), kThreadsPerBlock>>>(boxes, g, s.sorted.ulongs(), registrations, entries.uints());
    checkLaunch("device grid");

    view.cell = cell;
    for (int a = 0; a < 3; ++a) {
        view.lo[a] = g.lo[a];
        view.hi[a] = g.lo[a] + g.n[a] - 1;
    }
    view.cellStart = cellStart.uints();
    view.entries = entries.uints();
    view.boxLo = boxes.lo;
    view.boxHi = boxes.hi;
    view.nItems = items;
}

std::vector<int> itemsOverlapping(const GridView& hostGrid, const float lo[3], const float hi[3]) {
    std::vector<int> items;
    for (GridCursor cursor(hostGrid, fromF(lo), fromF(hi)); cursor.next();) items.push_back(cursor.item());
    return items;
}

// --- upload

void DeviceMesh::upload(const MeshArrays& m) {
    pos.upload(m.pos);
    tri.upload(m.tri);
    vnormal.upload(m.vnormal);
    faceN.upload(m.faceN);
    pvertN.upload(m.pvertN);
    edgeN.upload(m.edgeN);
    nbr.upload(m.nbr);
    vtStart.upload(m.vtStart);
    vtList.upload(m.vtList);
    grid.upload(m.boxes);
    view.pos = pos.floats();
    view.tri = tri.ints();
    view.vnormal = vnormal.floats();
    view.faceN = faceN.floats();
    view.pvertN = pvertN.floats();
    view.edgeN = edgeN.floats();
    view.nTri = static_cast<int>(m.tri.size() / 3);
    view.nPos = static_cast<int>(m.pos.size() / 3);
    for (int a = 0; a < 3; ++a) {
        view.lo[a] = m.lo[a];
        view.hi[a] = m.hi[a];
    }
    view.grid = grid.view;
    adjacency = {nbr.ints(), vtStart.ints(), vtList.ints()};
    bytes = pos.bytes() + tri.bytes() + vnormal.bytes() + faceN.bytes() + pvertN.bytes() + edgeN.bytes() + nbr.bytes() + vtStart.bytes() + vtList.bytes() +
            grid.bytes();
}

void DeviceCurves::upload(const CurveArrays& c) {
    curves.upload(c.curves.data(), c.curves.size() * sizeof(CurveRecord));
    ctrl.upload(c.ctrl);
    tu.upload(c.tu);
    ts.upload(c.ts);
    preLo.upload(c.preLo);
    preHi.upload(c.preHi);
    thickness.upload(c.thickness);
    view.curves = static_cast<CurveRecord*>(curves.data());
    view.nCurves = static_cast<int>(c.curves.size());
    view.ctrl = ctrl.floats();
    view.tu = tu.doubles();
    view.ts = ts.doubles();
    view.preLo = c.preLo.empty() ? nullptr : preLo.doubles();
    view.preHi = c.preHi.empty() ? nullptr : preHi.doubles();
    view.thickness = c.thickness.empty() ? nullptr : thickness.floats();
    spanOffset.clear();
    int total = 0;
    for (const CurveRecord& f : c.curves) {
        spanOffset.push_back(total);
        total += f.nCtrl >= 4 ? f.nCtrl - 3 : 0;
    }
    spanOffset.push_back(total);
}

void DeviceSurface::upload(const SurfaceArrays& s) {
    yarns.upload(s.yarns.data(), s.yarns.size() * sizeof(YarnSurfaceRecord));
    rings.upload(s.rings.data(), s.rings.size() * sizeof(RingRecord));
    view.yarns = static_cast<YarnSurfaceRecord*>(yarns.data());
    view.nYarns = static_cast<int>(s.yarns.size());
    view.rings = static_cast<RingRecord*>(rings.data());
}

// --- resident meshes

namespace {

enum class MeshUse : uint64_t { InsideQueries = 2, ClosestPoints = 4, Continuation = 5 };

uint64_t residentKey(const TriMesh& m, MeshUse use) {
    uint64_t h = fnv1a64Bytes(&use, sizeof use);
    h = fnv1a64Bytes(m.pos.data(), m.pos.size() * sizeof(Vec3), h);
    return fnv1a64Bytes(m.tri.data(), m.tri.size() * sizeof(m.tri[0]), h);
}

// MeshSurface's adjacency (base_surface.cpp): neighbour across edge k (corners k, k + 1) where exactly
// two triangles share the edge; triangles around each vertex in triangle order.
void addAdjacency(const TriMesh& m, MeshArrays& out) {
    const int T = static_cast<int>(m.tri.size());
    out.nbr.assign(3 * static_cast<size_t>(T), -1);
    std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> edges;
    std::vector<std::vector<int>> around(m.pos.size());
    for (int t = 0; t < T; ++t)
        for (int k = 0; k < 3; ++k) {
            const int a = m.tri[t][k], b = m.tri[t][(k + 1) % 3];
            edges[{std::min(a, b), std::max(a, b)}].push_back({t, k});
            around[a].push_back(t);
        }
    for (const auto& [edge, sides] : edges)
        if (sides.size() == 2) {
            out.nbr[3 * sides[0].first + sides[0].second] = sides[1].first;
            out.nbr[3 * sides[1].first + sides[1].second] = sides[0].first;
        }
    out.vtStart.assign(m.pos.size() + 1, 0);
    for (size_t v = 0; v < around.size(); ++v) {
        out.vtStart[v] = static_cast<int>(out.vtList.size());
        out.vtList.insert(out.vtList.end(), around[v].begin(), around[v].end());
    }
    out.vtStart[around.size()] = static_cast<int>(out.vtList.size());
}

std::shared_ptr<DeviceMesh> makeResident(uint64_t key, const MeshArrays& arrays) {
    Reservation lease(arrays.bytes());  // the grid's kernels and scratch
    auto mesh = std::make_shared<DeviceMesh>();
    mesh->upload(arrays);
    residentPut(key, mesh);
    return mesh;
}

} // namespace

std::shared_ptr<DeviceMesh> residentMesh(const TriMesh& m, MeshWalk walk) {
    const uint64_t key = residentKey(m, walk == MeshWalk::Adjacency ? MeshUse::Continuation : MeshUse::ClosestPoints);
    if (std::shared_ptr<Resident> found = residentGet(key)) return std::static_pointer_cast<DeviceMesh>(found);
    MeshArrays arrays = flattenMesh(m);
    if (walk == MeshWalk::Adjacency) addAdjacency(m, arrays);
    return makeResident(key, arrays);
}

std::shared_ptr<DeviceMesh> residentMesh(const MeshQuery& mq) {
    const uint64_t key = residentKey(mq.mesh(), MeshUse::InsideQueries);
    if (std::shared_ptr<Resident> found = residentGet(key)) return std::static_pointer_cast<DeviceMesh>(found);
    return makeResident(key, flattenMesh(mq));
}

std::vector<int> yarnOfSpans(const UploadedCurves& d) {
    std::vector<int> yarn(static_cast<size_t>(d.nSpans));
    for (size_t y = 0; y + 1 < d.dev.spanOffset.size(); ++y)
        for (int k = d.dev.spanOffset[y]; k < d.dev.spanOffset[y + 1]; ++k) yarn[k] = static_cast<int>(y);
    return yarn;
}

UploadedCurves uploadCurves(const CurveSet& cs) {
    UploadedCurves out;
    out.host = &cs;
    out.dev.upload(flattenCurves(cs, ArcBounds::Included));
    out.nSpans = out.dev.spanOffset.empty() ? 0 : out.dev.spanOffset.back();
    return out;
}

} // namespace knit::gpu
