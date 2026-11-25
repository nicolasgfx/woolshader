#pragma once
// Paths on a triangle mesh between two surface points: what joins a marking's control points (surface_marking.h).
//
// The path between a and b is the section of the surface with the plane through a and b that contains the mean of
// their triangles' normals (a "straight" line drawn on the surface), walked triangle by triangle from a until it
// enters b's triangle. A closed section can be walked either way round; the shorter walk wins. Every point between
// a and b lies on a mesh edge, so a closed path splits triangles cleanly (stitchmesh/surface_trim.h). When neither
// walk reaches b (the plane's section through a does not pass b, an open edge on the way), the path falls back to the
// shortest path along mesh edges (Dijkstra), which is marked: it is drawn but cannot trim a surface.
#include <array>
#include <string>
#include <vector>

#include "model/mesh.h"
#include "model/surface_marking.h"

namespace knit {

// across[t][k]: the triangle across edge (corner k, corner k + 1) of triangle t; -1 on an open or non-manifold edge.
struct TriangleNeighbours {
    std::vector<std::array<int, 3>> across;
};
TriangleNeighbours triangleNeighbours(const TriMesh& m);

// One point of a path. A point between control points lies on edge (edgeFrom, edgeTo) at pos = (1 - s) from + s to.
struct TracePoint {
    Vec3 pos;
    int tri = -1;                    // the triangle the path runs through after this point (the last point: its own)
    int edgeFrom = -1, edgeTo = -1;  // -1: a control point inside a triangle, or a vertex of an edge path
    double s = 0;
};

// A marking's whole path: its control points joined in order (a closed curve also joins the last to the first).
struct MarkingTrace {
    std::vector<TracePoint> points;
    std::vector<int> controlAt;      // index in points of each control point
    bool edgePathUsed = false;       // some segment fell back to mesh edges (drawn, but cannot trim)
    std::vector<int> failedSegments; // segment k (control point k to k + 1) was not traced at all
};

Vec3 surfacePosition(const TriMesh& m, const SurfacePoint& p);
// The surface point of m closest to p (brute force over the triangles: for scripts, tests and single clicks).
SurfacePoint closestSurfacePoint(const TriMesh& m, const Vec3& p);

// The path from a to b (both on m), appended to out without its first point when out is not empty (so segments
// chain). Returns false when even the edge path fails (different connected pieces); usedEdgePath is set when the plane
// walk failed and the edge path was taken.
[[nodiscard]] bool traceSurfacePath(const TriMesh& m, const TriangleNeighbours& nb, const SurfacePoint& a, const SurfacePoint& b,
                                    std::vector<TracePoint>& out, bool& usedEdgePath);

MarkingTrace traceMarking(const TriMesh& m, const TriangleNeighbours& nb, const std::vector<SurfacePoint>& points, bool closed);

// Whether every control point lies on one and the same border loop of m: on a border edge (one weight exactly 0, the edge
// used by one triangle only, so a non-manifold edge is no border) or on a border vertex (one weight exactly 1). loop
// receives that loop's vertex ids, walked the way the control points go round it. A closed Start or Finish whose points
// all lie on a border loop means that loop itself: an open tube's real opening (marked_build.h).
bool pointsOnOneBorderLoop(const TriMesh& m, const std::vector<SurfacePoint>& points, std::vector<int>& loop);
// The loop as a closed trace on m: its vertices in order with each control point inserted on its edge (controlAt in the
// points' order; edgeFrom -1), ending where it began. points: the control points pointsOnOneBorderLoop accepted.
MarkingTrace borderLoopTrace(const TriMesh& m, const std::vector<int>& loop, const std::vector<SurfacePoint>& points);

} // namespace knit
