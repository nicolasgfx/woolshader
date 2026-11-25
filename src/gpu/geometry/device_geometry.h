#pragma once
// GPU backend: yarns and meshes on the device, in three forms. A view is a plain struct of device
// pointers, passed by value to the kernels. Arrays are the host's copy of exactly those bytes: a run
// sizes its reservation from them before it uploads. A Device* class owns the buffers behind a view.
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "gpu/device/primitives.h"
#include "gpu/device/runtime.h"
#include "knitting/tube/yarn_surface.h"
#include "model/mesh_query.h"

namespace knit::gpu {

// --- views
// Uniform grid over item boxes (SpatialHash3's visit rule): cells over the occupied range
// [lo, hi] (inclusive, per axis), CSR buckets cellStart[c] .. cellStart[c + 1] into
// `entries`; an entry is the item id (low 29 bits) with its "lowest cell" flags (i 1, j 2,
// k 4) in bits 29..31. Cells are floor(float(v) * float(1 / cell)) on host and device alike.
struct GridView {
    float cell = 1;
    long long lo[3] = {0, 0, 0}, hi[3] = {-1, -1, -1};
    const unsigned* cellStart = nullptr;  // nCells + 1
    const unsigned* entries = nullptr;
    const float* boxLo = nullptr;  // 3 per item
    const float* boxHi = nullptr;
    int nItems = 0;
};

// A curve set: control points (float), arc-length tables (double; see math.cuh).
struct CurveRecord {
    int ctrlOff = 0, nCtrl = 0;
    int tabOff = 0, nTab = 0;
    double length = 0;
    int preOff = 0;     // lengthBounds prefix sums (nCtrl - 3 + 1 entries), -1: none
    int spanOff = 0;    // global span index of span 0
};
struct CurveSetView {
    const CurveRecord* curves = nullptr;
    int nCurves = 0;
    const float* ctrl = nullptr;
    const double* tu = nullptr;
    const double* ts = nullptr;
    const double* preLo = nullptr;  // YarnPath's compensated prefix sums of the span bounds
    const double* preHi = nullptr;
    const float* thickness = nullptr;  // per control, at ctrlOff (YarnPath::thickness); null: 1
};

// Surface model rings (knit::Ring), one array for all yarns.
struct RingRecord {
    double s, u;
    float c[3], T[3], n1[3], n2[3];
    float a, b, phi;
    float pad;
};
struct YarnSurfaceRecord {
    int ringOff = 0, nRings = 0;
    double ds = 0, length = 0;
    float radius = 0, ext = 0;  // nominal r; max over rings of max(a, b)
};
struct SurfaceView {
    const YarnSurfaceRecord* yarns = nullptr;
    int nYarns = 0;
    const RingRecord* rings = nullptr;
};

// SurfaceIndex: capsules around the centrelines (knit::SurfaceIndex::Capsule) and their grid.
struct CapsuleRecord {
    float a[3], b[3];
    float R;
    int yarn;
    double s0, s1;
};
struct SurfaceIndexView {
    const CapsuleRecord* caps = nullptr;
    GridView grid;
};

// A triangle mesh with its closest-point grid. vnormal: area-weighted vertex normals (the
// chart / MeshSurface shading normals); faceN, pvertN, edgeN: MeshQuery's pseudo-normals
// (inside test); edgeN holds 3 per triangle (the edge opposite corner k).
struct MeshView {
    const float* pos = nullptr;    // 3 per vertex
    const int* tri = nullptr;      // 3 per triangle
    const float* vnormal = nullptr;
    const float* faceN = nullptr;
    const float* pvertN = nullptr;
    const float* edgeN = nullptr;
    int nTri = 0, nPos = 0;
    float lo[3] = {0, 0, 0}, hi[3] = {0, 0, 0};
    GridView grid;
};

// --- host arrays
// The boxes of a grid's items, 3 floats per item. The grid itself is built on the device (DeviceGrid).
struct BoxArrays {
    std::vector<float> lo, hi;
    float cell = 1;
    int count() const { return static_cast<int>(lo.size() / 3); }
    size_t gridBytes() const;  // the grid and what building it takes
};

struct CurveArrays {
    std::vector<CurveRecord> curves;
    std::vector<float> ctrl;
    std::vector<double> tu, ts;
    std::vector<double> preLo, preHi;  // empty: not flattened (the fibres do not need them)
    std::vector<float> thickness;      // empty: no yarn has one
    size_t bytes() const;
};
struct SurfaceArrays {
    std::vector<YarnSurfaceRecord> yarns;
    std::vector<RingRecord> rings;
    size_t bytes() const;
};
struct IndexArrays {
    std::vector<CapsuleRecord> caps;
    BoxArrays boxes;
    float zTop = 0, zBottom = 0;  // the highest and lowest points of the capsules
    size_t bytes() const;
};
struct MeshArrays {
    std::vector<float> pos, vnormal, faceN, pvertN, edgeN;
    std::vector<int> tri;
    std::vector<int> nbr, vtStart, vtList;  // AdjacencyView; empty without MeshWalk::Adjacency
    float lo[3] = {0, 0, 0}, hi[3] = {0, 0, 0};
    BoxArrays boxes;
    size_t bytes() const;
};

// --- owners
// What building a grid needs for a moment. A caller that rebuilds often keeps one, so nothing is
// reallocated.
struct GridScratch {
    DeviceBuffer range, counts, offsets, keys, sorted, cellCounts;
};
struct GridLimits {
    long long maxCells = 1ll << 24;     // the cell grows until the occupied range has at most this many
    long long maxEntries = 0x7fffffff;  // more cell registrations fail: arrays and kernels index with int
};

// A uniform grid of boxes, built on the device. Every bucket lists its items in increasing id order,
// so a rerun gives the same grid. The caller fills boxLo and boxHi (3 floats per item), or passes
// BoxArrays to upload(); an item with lo.x > hi.x stays out of the grid.
struct DeviceGrid {
    DeviceBuffer boxLo, boxHi, cellStart, entries;
    GridView view;
    size_t bytes() const { return cellStart.bytes() + entries.bytes() + boxLo.bytes() + boxHi.bytes(); }
    void upload(const BoxArrays& boxes);
    void build(int items, float cell);
    void build(int items, float cell, DevicePrimitives& prims, GridScratch& scratch, const GridLimits& limits = {});
};

// MeshSurface's neighbours, in MeshSurface's order.
struct AdjacencyView {
    const int* nbr;      // 3 per triangle: neighbour across edge k (corners k, k + 1), -1 boundary
    const int* vtStart;  // per vertex + 1: CSR into vtList
    const int* vtList;   // triangles around a vertex (MeshSurface::vertexTris_ order)
};

struct DeviceMesh : Resident {
    DeviceBuffer pos, tri, vnormal, faceN, pvertN, edgeN, nbr, vtStart, vtList;
    DeviceGrid grid;
    MeshView view;
    AdjacencyView adjacency{};
    void upload(const MeshArrays& m);
};

// A curve set on the device: controls (float), arc-length tables (double), the lengthBounds
// prefix sums (double, may be empty) and the global span numbering (span k of yarn y is global
// span spanOffset[y] + k).
struct DeviceCurves {
    DeviceBuffer curves, ctrl, tu, ts, preLo, preHi, thickness;
    CurveSetView view;
    std::vector<int> spanOffset;  // host copy, per yarn plus the end
    size_t bytes() const { return curves.bytes() + ctrl.bytes() + tu.bytes() + ts.bytes() + preLo.bytes() + preHi.bytes() + thickness.bytes(); }
    void upload(const CurveArrays& c);
};

struct DeviceSurface {
    DeviceBuffer yarns, rings;
    SurfaceView view;
    size_t bytes() const { return yarns.bytes() + rings.bytes(); }
    void upload(const SurfaceArrays& s);
};

// --- flatten: knit types to host arrays
enum class ArcBounds { Omitted, Included };  // YarnPath's lengthBounds prefix sums (CurveRecord::preOff)
CurveArrays flattenCurves(const CurveSet& cs, ArcBounds bounds);
SurfaceArrays flattenSurface(const SurfaceModel& sm);
IndexArrays flattenIndex(const SurfaceModel& sm);  // SurfaceIndex's capsules
MeshArrays flattenMesh(const TriMesh& m);          // shading and face normals: closest-point queries
MeshArrays flattenMesh(const MeshQuery& mq);       // MeshQuery's pseudo-normals: inside queries

// Appends `in` as float triples, the layout of the host arrays.
inline void appendFloat3(const std::vector<Vec3>& in, std::vector<float>& out) {
    size_t at = out.size();
    out.resize(at + 3 * in.size());
    for (const Vec3& p : in) {
        out[at++] = float(p.x);
        out[at++] = float(p.y);
        out[at++] = float(p.z);
    }
}

// --- a curve set with its host source (shared by the validation checks), its global spans, its tables
// One curve set on the device (controls, arc-length tables, lengthBounds prefix sums), with its host
// source. Built once per validation and shared by its checks.
struct UploadedCurves {
    DeviceCurves dev;
    const CurveSet* host = nullptr;
    int nSpans = 0;
};
UploadedCurves uploadCurves(const CurveSet& cs);

// The global spans of a curve set on the device: span g belongs to yarn `yarn[g]`.
struct SpanList {
    CurveSetView curves;
    const int* yarn;
    int count;
};
std::vector<int> yarnOfSpans(const UploadedCurves& d);

// Every yarn's arc-length tables (YarnPath::setTables) from its current controls (arc_length.cu).
void computeArcLengthTables(CurveSet& cs);

// --- resident meshes: shared by the pieces of a scene, on the device once per content and use, until
// releaseResident().
enum class MeshWalk { None, Adjacency };  // Adjacency: MeshSurface's continuation walks the triangles
std::shared_ptr<DeviceMesh> residentMesh(const TriMesh& m, MeshWalk walk);
std::shared_ptr<DeviceMesh> residentMesh(const MeshQuery& mq);

// The items a GridCursor visits for the box [lo, hi], walked on the host (tests).
std::vector<int> itemsOverlapping(const GridView& hostGrid, const float lo[3], const float hi[3]);

} // namespace knit::gpu
