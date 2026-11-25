#pragma once
// A built stitch mesh as a Wavefront OBJ file (the Mark step's Export Stitch Mesh), and read back (tests and checks).
//
// The file keeps what the Stitch Mesh view shows: one `f` line per stitch polygon (quads, the five-sided increase and
// decrease faces and the pole-start triangles stay polygons, never split), the vertices exactly as built (`%.17g`, the
// imported file's units and axes, no transform and no display lift), and one group per row and stitch kind,
// `g row_<NNN>_<kind>` with the kind named as stitchmesh_lab names it (K, INC_kp, DEC_d12k, special_quad, unresolved,
// pole_start), so rows and shaping can be read from the file itself. `#` lines carry the caller's comments (source,
// piece, unit, the quality headline). Like the authors' exportResult (Wu et al., stitchMeshing), a mesh with findings is
// written as it is.
#include <string>
#include <vector>

#include "model/stitch_mesh_data.h"

namespace knit {

// The file name extension of a stitch mesh file (without the dot): the save dialog's filter and default extension.
constexpr const char* kStitchMeshExtension = "obj";

// The group name of a row and stitch kind: "row_007_INC_kp" (a face without a row: "row_none_K").
std::string stitchGroupName(int row, StitchKind kind);

// Writes `mesh` to `path` with `comments` as `#` lines and `objectName` as its `o` line. False with a message when the mesh
// cannot be drawn (stitchMeshUse) or the file cannot be written.
[[nodiscard]] bool writeStitchMeshObj(const std::string& path, const StitchMeshData& mesh, const std::vector<std::string>& comments, const std::string& objectName,
                                      std::string* err);

// Reads a polygon OBJ (v, f and g lines; `f` corners as v, v/vt, v//vn or v/vt/vn, 1-based or negative). Row and kind come
// from row groups named by stitchGroupName (row -1 and Regular otherwise); findings and deviations are not in the file.
[[nodiscard]] bool readPolygonObj(const std::string& path, StitchMeshData& out, std::string* err);

} // namespace knit
