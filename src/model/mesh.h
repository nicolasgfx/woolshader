#pragma once
// Triangle mesh with per-corner chart coordinates, and the closest-point primitives that queries and checks
// are built on. The search over a mesh is mesh_query.h.
//
// Geometry is shared per vertex (pos, area-weighted vertex normals); the chart
// coordinates (a, b) are stored per triangle CORNER, so a seam is an explicit
// discontinuity: the same vertex may carry b on one side and b + 1 on the other
// (as the scene pieces are exported).
#include <array>
#include <cstdint>
#include <vector>

#include "model/vec.h"

namespace knit {

struct TriMesh {
    std::vector<Vec3> pos;
    std::vector<std::array<int, 3>> tri;   // vertex indices, counter-clockwise seen from outside
    std::vector<std::array<Vec2, 3>> uv;   // per corner (a, b); empty if the file had none
    std::vector<Vec3> vnormal;             // area-weighted, normalized

    void computeNormals();
    // p -> s * (p - origin); normals are unchanged (s > 0).
    void transform(double s, const Vec3& origin);
    double area() const;
    Vec3 centroid() const;  // area-weighted
    Vec3 faceNormal(int t) const;  // unit
    void bounds(Vec3& lo, Vec3& hi) const;
};

// The closest point on triangle (a, b, c) to p (Ericson 5.1.5), and its weights on the corners.
struct TrianglePoint {
    Vec3 p;
    double bary[3] = {0, 0, 0};  // exact 0 and 1 where p is on an edge or a corner
};
TrianglePoint closestOnTriangle(const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c);

// The closest pair of points of two segments, p0 + s (p1 - p0) and q0 + t (q1 - q0), and their distance
// (Ericson 5.1.9; exact for segments, handles degenerate ones).
struct ClosestSegmentPoints {
    double distance = 0;
    double s = 0;  // on p0 p1, in [0, 1]
    double t = 0;  // on q0 q1, in [0, 1]
};
ClosestSegmentPoints closestSegmentPoints(const Vec3& p0, const Vec3& p1, const Vec3& q0, const Vec3& q1);

// FNV-1a 64 of the positions and triangles, bit for bit: a mesh's identity in stored knits.
uint64_t meshDigest(const TriMesh& m);

} // namespace knit
