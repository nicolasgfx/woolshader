#pragma once
// Trimming a piece along closed curves (the termination rule): a start curve and a finish or
// boundary curve split a closed piece into caps and the band between them, and the band is what Method A knits (its
// cast-on and bind-off loops are the two curves). The curves are marking paths (model/surface_path.h) whose points
// between control points lie on mesh edges; each crossed triangle is split into three along the chord between its
// entry and exit crossing, so the new boundary follows the drawn curve, not the mesh's zig-zag edges. A control point
// inside a triangle is dropped from the cut: the chord between the triangle's entry and exit replaces it (an error
// below one triangle's size).
#include <string>
#include <vector>

#include "model/mesh.h"
#include "model/surface_path.h"

namespace knit {

struct TrimmedSurface {
    TriMesh mesh;                          // the piece with every crossed triangle split
    std::vector<int> parentTri;            // per triangle: the piece triangle it came from
    std::vector<int> region;               // per triangle: its region (connected without crossing a curve)
    int regions = 0;
    std::vector<std::vector<int>> loops;   // per curve: its vertices in the curve's order
    std::vector<std::vector<char>> touches;  // touches[region][curve]: the region borders on the curve
};

// curves: closed marking paths on `piece`. False when a curve is not usable for trimming (it ran along mesh edges, a
// triangle is crossed twice, two curves cross the same edge), with a message for the user.
[[nodiscard]] bool splitAlongCurves(const TriMesh& piece, const std::vector<MarkingTrace>& curves, TrimmedSurface& out, std::string* err);

// The triangles of one region as their own mesh (vertices renumbered), with each triangle's piece triangle and the
// loops renumbered into it (-1 entries dropped: a loop vertex not in the region).
struct RegionMesh {
    TriMesh mesh;
    std::vector<int> parentTri;
    std::vector<std::vector<int>> loops;
    std::vector<int> local;   // per trimmed-surface vertex: its index in the region (-1: not in it)
};
RegionMesh regionMesh(const TrimmedSurface& s, int region);

} // namespace knit
