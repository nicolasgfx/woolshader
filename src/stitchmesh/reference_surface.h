#pragma once
// The imported triangle surface as the stitch-mesh layout sees it: the one truth for projection, dimensions, component
// ownership and colour. Every coarse vertex is bound to a triangle and barycentric weights of this surface
// (VertexBinding, as knitting's ProjState holds tri and b), so defect attribution and later stages read the same surface.
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "model/mesh.h"
#include "model/mesh_check.h"
#include "model/mesh_query.h"

namespace knit {

struct ReferenceSurface {
    const TriMesh* mesh = nullptr;        // must outlive the surface (query keeps a reference)
    MeshCheck check;                      // per component: manifold, orientation, closed, self-intersection
    int degenerateTriangles = 0;          // zero-area triangles
    bool closed = false;                  // every component closed; an open sheet is allowed and reported
    std::unique_ptr<MeshQuery> query;     // closest point
    std::vector<Vec3> faceNormal;         // per triangle, unit, outward
    double area = 0, meanEdge = 0;
    uint64_t digest = 0;                  // meshDigest
};

// False, with the reason, when a component is non-manifold or not consistently oriented: nothing can be remeshed on it.
bool prepareReferenceSurface(const TriMesh& m, ReferenceSurface& out, std::string* why);

struct VertexBinding {
    int tri = -1;                  // -1: nothing found
    double bary[3] = {0, 0, 0};
    double distance = 0;           // from the point to its foot [mm]
};
std::vector<VertexBinding> bindToSurface(const ReferenceSurface& s, const std::vector<Vec3>& points);

} // namespace knit
