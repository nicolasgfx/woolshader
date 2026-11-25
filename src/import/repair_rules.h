#pragma once
// The rules the repairs apply (prepare_import.cpp) and the analysis counts by (import_analysis.cpp), written once so the two
// cannot drift: coincident vertices, zero-area, duplicate and unreferenced elements, and the mesh measures they rest on (the
// bounding box, a triangle's area, the edge list). Both work on the mesh as the file has it, in file units or library mm.
#include <array>
#include <cstdint>
#include <utility>
#include <vector>

#include "model/vec.h"

namespace knit {

constexpr double kWeldRelative = 1e-6;       // coincident: within this share of the bounding box diagonal
constexpr double kZeroAreaRelative = 1e-12;  // zero area: at most this share of the diagonal squared

struct BoundingBox {
    Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};  // as they start when there are no points
    double diagonal = 0;                                       // 0 when there are no points
};

BoundingBox boundingBoxOf(const std::vector<Vec3>& pos);

double triangleArea(const std::vector<Vec3>& pos, const std::array<int, 3>& t);

// A repeated corner, or an area of at most kZeroAreaRelative * diagonal^2.
bool isZeroArea(const std::vector<Vec3>& pos, const std::array<int, 3>& t, double diagonal);

// The corners in ascending order: two triangles with the same three corners, in any order, are duplicates.
std::array<int, 3> sortedCorners(const std::array<int, 3>& t);

// The vertices no triangle uses.
long long unreferencedVertices(const std::vector<std::array<int, 3>>& tri, size_t vertexCount);

// Coincident vertices (within eps, exact distance): remap[v] = the first vertex of v's cluster (an INPUT index; chains are not
// merged beyond eps of that vertex). Not model's weldVertices: that one returns compacted indices and breaks ties in another
// order, and import identities depend on this one.
std::vector<int> coincidentVertices(const std::vector<Vec3>& pos, double eps);

// An undirected edge as one number: the lower vertex in the high 32 bits, the higher in the low 32.
uint64_t edgeKey(int a, int b);

// Every triangle corner's edge (k: corners k, k + 1) as (edgeKey, 3 t + k), sorted by key: the uses of one edge are adjacent.
std::vector<std::pair<uint64_t, int>> sortedEdges(const std::vector<std::array<int, 3>>& tri);

} // namespace knit
