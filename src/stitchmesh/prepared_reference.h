#pragma once
// The prepared reference of a marked piece and its correspondence with the file:
// a piece's triangles are midpoint-refined so the row field is sampled finely
// enough, without moving the surface; every refined triangle knows the piece triangle it lies in and its corners'
// barycentric coordinates there, so a marking's file point maps to the refined mesh and back exactly. Clipping below a
// level of the row field gives the knit region of a finish-level request, its border on that level.
#include <array>
#include <string>
#include <vector>

#include "model/mesh.h"
#include "model/surface_marking.h"

namespace knit {

struct RefinedMesh {
    TriMesh mesh;
    std::vector<int> parent;                                    // per refined triangle: the input triangle it lies in
    std::vector<std::array<std::array<double, 3>, 3>> cornerBary;  // per refined triangle: its corners in the parent
    int levels = 0;
};
// How the piece becomes the reference the mesher lays out on and relaxes on.
//   Facets     the piece as imported, midpoint-refined: the surface does not move (the integration's reference)
//   LoopLimit  Loop subdivision of the piece's triangles, the last level's vertices moved to their limit positions: a
//              triangle mesh inscribed in the Loop limit surface (an approximation of that smooth surface, finer with
//              every level). The split is the same as Facets', so parent and corner coordinates are identical: a point
//              of the piece maps to the same refined triangle and weights under either preparation (the subdivision's
//              own parametrization); only the positions differ.
enum class ReferencePreparation : char { Facets, LoopLimit };
const char* toString(ReferencePreparation p);   // "facets", "loop_limit"
bool parsePreparation(const std::string& text, ReferencePreparation& out);

// levels times: every triangle into four at its edge midpoints (on its plane: the surface does not change).
RefinedMesh refineMidpoint(const TriMesh& in, int levels);
// levels Loop steps, then the limit positions (see ReferencePreparation); levels 0 gives the limit of the piece's own vertices.
RefinedMesh refineLoop(const TriMesh& in, int levels);
RefinedMesh prepareReference(const TriMesh& in, int levels, ReferencePreparation p);
// A point of one prepared reference of a piece as a point of another of the same piece (any preparation, any levels),
// through the piece's parametrization: the same triangle and weights when both have the same split.
SurfacePoint mapBetween(const RefinedMesh& from, const RefinedMesh& to, const SurfacePoint& p);
// The smallest number of levels that brings the mean edge length to at most maxEdge (0 if it already is), capped at 4.
int refinementLevels(const TriMesh& in, double maxEdge);
// A point of the input mesh as a point of the refined one (the refined triangle containing it) and back.
SurfacePoint toRefined(const RefinedMesh& r, const SurfacePoint& input);
SurfacePoint fromRefined(const RefinedMesh& r, const SurfacePoint& refined);

// The mesh below a level of a vertex field: triangles split along the level (the edge points extractLevel takes), so the
// border is the level's contour. parent: per output triangle, its triangle in m.
TriMesh clipBelowLevel(const TriMesh& m, const std::vector<double>& t, double level, std::vector<int>& parent);
// The border loops of a triangle mesh, vertex ids in the order the triangles walk them.
std::vector<std::vector<int>> borderLoops(const TriMesh& m);

} // namespace knit
