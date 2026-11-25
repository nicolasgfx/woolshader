#pragma once
// The filling of a knitted part: a surface behind the yarn that shows through its gaps, written by
// the export (filling.enabled) for rendering only. It is the part's source surface moved by an
// offset along -N (vertex normals), fillingOffset (stage_params.h), so it lies behind the yarn:
//   * closed, consistently oriented, positive volume: moved inward, "stuffing";
//   * open (boundary edges), consistently oriented, manifold: moved to the side opposite the
//     knitted +N side, "backing";
//   * otherwise none, with the reason.
// The moved surface is reported, never repaired: triangles whose normal flipped against the source,
// its triangle-pair self-intersection (model/mesh_check.h) and, for stuffing, the sign of its volume
// (negative: turned inside out). It is not a check of the knit.
#include <string>

#include "model/mesh.h"
#include "model/mesh_check.h"
#include "model/status.h"

namespace knit {

struct Filling {
    enum class Kind { None, Stuffing, Backing } kind = Kind::None;
    std::string reason;  // why there is none
    TriMesh mesh;
    double offset = 0;   // [mm] along -N
    int inverted = 0;    // triangles whose normal flipped against the source's
    Status selfIntersection = Status::Unresolved;
    int intersectingPairs = 0, unresolvedPairs = 0;
    double signedVolume = 0;  // stuffing only
    // One line for the export log.
    std::string summary() const;
};

// The filling of `mesh`, one component whose check is `component`, at `offset` [mm].
Filling makeFilling(const TriMesh& mesh, const ComponentCheck& component, double offset);

} // namespace knit
