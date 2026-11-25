#pragma once
// The placement's surface queries on the device, private to the backend (surface/ implements them):
// the knitting port (port/knit_numerics.cpp) hands them to knitting through KnitNumerics.
#include <memory>
#include <string>
#include <vector>

#include "knitting/layout/base_surface.h"
#include "knitting/layout/chart.h"
#include "knitting/layout/chart_lattice.h"
#include "knitting/layout/uv_atlas.h"
#include "model/mesh_query.h"

namespace knit::gpu {

// Placement (knitPiece via PieceSurface::placementQuery): the piece surface's queries (at: chart
// evaluation and pole caps; project: closest point, interpolated vertex normals) in batches on
// the device; geometry uploaded once per object.
struct PlacementStats {
    long long queries = 0, batches = 0, misses = 0;  // misses: (a, b) outside every chart triangle
    double ms = 0;
    std::string line() const;
};
// From the piece's chart (PieceSurface semantics; KnitNumerics::pieceQuery).
std::unique_ptr<BaseSurface> makePieceQuery(const Chart& c, PlacementStats* st, std::string* err);
// The tiled pipeline's mesh surface (KnitNumerics::meshSurface): a MeshSurface whose batched
// projections and continuations run on the device (the settle's band uses the device mesh too).
std::unique_ptr<MeshSurface> makeMeshSurface(const TriMesh& m, PlacementStats* st, std::string* err);
// The boundary checks (BoundaryCheck, base_surface.h) on the device: B1 (float32 diagnostic of the certified distance
// check) and B2 (sampled continuations).
struct BoundaryRequest {
    const CurveSet& curves;
    const SurfaceBinding& binding;
    const MeshSurface& surface;
    double Rs, step, minStep;
};
bool checkBoundary(const BoundaryRequest& request, BoundaryCheck& out, std::string* err);
// The tiled lattice's node classification (KnitNumerics::classifyLattice): chart location, the inside
// test against `solids` (points shifted by `offset` into their frame), the footprint tests.
struct LatticeRequest {
    const UvAtlas& atlas;
    const LatticeChartQuery& query;
    const std::vector<const MeshQuery*>& solids;
    Vec3 offset;
};
bool classifyLattice(const LatticeRequest& request, std::vector<LatticeNodeClass>& out, std::string* err);

} // namespace knit::gpu
