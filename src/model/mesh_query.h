#pragma once
// The search over a mesh (mesh.h): the closest point to a query and whether it lies inside.
#include <vector>

#include "model/mesh.h"
#include "model/spatial_hash.h"

namespace knit {

// Closest-point and inside queries on a closed triangle mesh: a 3D hash of the triangle boxes, the closest
// point by an expanding search, and the sign from the angle-weighted pseudo-normal of the closest feature
// (face, edge or vertex; Baerentzen & Aanaes 2005), exact for closed, consistently oriented manifold meshes.
// Replaces the O(triangles) winding number for many queries. The mesh must outlive the query object, which
// does not change after its constructor: queries are stateless, so concurrent const use is safe (spatial_hash.h).
class MeshQuery {
public:
    explicit MeshQuery(const TriMesh& m, double cell = 0);  // cell 0: 2 x the mean edge length
    struct Hit {
        int tri = -1;        // -1: nothing within maxDist
        Vec3 p;
        double dist = 1e300;
        bool inside = false;  // q is on the inner side of the closest feature's pseudo-normal
    };
    // Expanding search (query box doubles until the best hit is within it). Triangles that
    // cannot improve the strict minimum are not evaluated: those already evaluated in an
    // earlier round (box overlapping the previous round's box), and those whose box is
    // farther than the current best by more than a rounding margin.
    Hit closest(const Vec3& q, double maxDist = 1e300) const;
    bool inside(const Vec3& q) const;  // false outside the bounding box
    // Whether the segment from a to b crosses a triangle at a point more than `skip` from a (Moeller-Trumbore,
    // "Fast, Minimum Storage Ray/Triangle Intersection", 1997, on the triangles whose boxes overlap the segment's
    // box). Floating-point geometry, not an exact predicate: a segment parallel to a triangle (|det| < 1e-18)
    // never crosses it, and a hit on a shared edge counts once per triangle. skip > 0 lets a segment start on the
    // surface, as a fibre's first segment does.
    bool segmentCrosses(const Vec3& a, const Vec3& b, double skip = 0) const;
    // Read-only views for the GPU port (gpu/geometry/device_geometry.cu flattens the mesh and its normals).
    const TriMesh& mesh() const { return m_; }
    const std::vector<Vec3>& faceNormals() const { return faceN_; }
    const std::vector<Vec3>& vertexNormals() const { return vertN_; }
    // Edge pseudo-normals (the sum of the edge's face normals), 3 per triangle: the edge opposite
    // corner k. The layout of MeshArrays::edgeN.
    const std::vector<Vec3>& cornerEdgeNormals() const { return cornerEdgeN_; }

private:
    const TriMesh& m_;
    SpatialHash3 hash_;
    std::vector<Vec3> faceN_, vertN_, cornerEdgeN_;
    Vec3 lo_, hi_;
};

} // namespace knit
