#pragma once
// Validity of source meshes before they are knitted or used for inside tests. Per connected component (triangles sharing vertices) of a WELDED mesh:
// manifoldness, orientation consistency, closedness, signed volume and self-
// intersection. A component is inside-test ELIGIBLE only if it is closed, consistently
// oriented, has positive volume and shows no self-intersection (no FAIL and no
// UNRESOLVED triangle pair); otherwise inside-based culling and fibre hiding must not
// use it (the reason is reported).
#include <string>
#include <vector>

#include "model/status.h"
#include "model/mesh.h"

namespace knit {

// A mesh after welding, and where each input vertex went.
struct WeldedMesh {
    TriMesh mesh;
    std::vector<int> remap;    // remap[v]: the welded index of input vertex v
    int droppedTriangles = 0;  // degenerate after welding (two equal corners)
};

// Welds vertices closer than eps (exact pairwise distance, 3D hash; the first vertex of
// a cluster wins, chains are not merged transitively beyond eps of that vertex).
// Triangles keep their per-corner uv; degenerate triangles are dropped and counted.
WeldedMesh weldVertices(const TriMesh& in, double eps);

struct ComponentCheck {
    std::vector<int> triangles;
    Vec3 lo, hi;              // bounds
    double area = 0;
    int vertices = 0, edges = 0;
    int boundaryEdges = 0, nonManifoldEdges = 0, inconsistentEdges = 0;  // inconsistent: two uses in the same direction
    int euler = 0;            // V - E + F
    bool closed = false, consistent = false;
    double signedVolume = 0;  // closed components only (else 0)
    int uvCharts = 0;         // components under UV-continuous edges (0 without uv)
    Status selfIntersection = Status::Unresolved;  // PASS: no pair intersects; FAIL; UNRESOLVED: near-degenerate pairs
    int intersectingPairs = 0, unresolvedPairs = 0;
    bool insideEligible = false;
    std::string reason;       // why not eligible
};

struct MeshCheck {
    std::vector<ComponentCheck> components;  // largest (triangle count) first
    std::string summary() const;
};

// Skip leaves every component's selfIntersection UNRESOLVED, hence not eligible for inside tests.
enum class SelfIntersectionTest { Run, Skip };

MeshCheck checkMesh(const TriMesh& m, SelfIntersectionTest test = SelfIntersectionTest::Run);

// Triangle-triangle intersection of non-adjacent triangles: PASS (disjoint), FAIL
// (they intersect), UNRESOLVED (within tol of touching: not decided).
Status triangleIntersection(const Vec3 a[3], const Vec3 b[3], double tol);

} // namespace knit
