#pragma once
// What the user chooses for a file (ImportOptions) and what a preparation makes of it (PreparedMesh): the repairs, the unit
// and up-axis contract and the identity of the result. What the file is before any choice is in import_analysis.h.
//
// Units: OBJ and PLY carry no unit. The caller states the size of one file unit in mm (the
// Import dialog asks; a suggestion from the bounding box is only a suggestion) and the file's
// up axis (Y: OBJ convention; Z). Library space is Z-up mm, the bounding box centred over the
// origin (x, y) with its lowest point at z = 0, as the GLB scenes (glb_load.h).
//
// The colour uv is kept per corner, separately from the knitting uv: a rebuilt surface knits on its
// own charts, while the file's colour is still sampled where the file put it.
#include <array>
#include <cstdint>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "io/mesh_read.h"
#include "model/mesh.h"
#include "model/coordinates.h"

namespace knit {

struct ImportOptions {
    double unitMm = 1;                           // mm per file unit
    UpAxis up = UpAxis::Y;
    bool weld = true, removeDegenerate = true, removeDuplicateFaces = true, removeUnreferenced = true, fixOrientation = true;
    // Wrapped: one closed surface at stitch scale with new charts and baked colours (import/wrap_surface.h); every knitted
    // import of the dialog. KeepFile: the file's own mesh and uv, never chosen in the dialog: an import whose knitting is
    // skipped (fibres grow on it), and scripted runs that knit the file's own uv (import.surface = file, the settings'
    // default, which this default matches).
    enum class Surface { KeepFile, Wrapped } surface = Surface::KeepFile;
    double wrapStitchMm = 5.5;                   // the stitch width a wrapped surface is prepared for
    std::set<std::string> plainNodes;            // ImportNode::name: kept as the file has them, not knitted
    bool operator==(const ImportOptions&) const = default;
};
std::string importOptionsText(const ImportOptions& o);   // stable "name = value" lines (identity, settings)

// import.plain_nodes: the names with | between them.
std::string plainNodesText(const std::set<std::string>& nodes);
std::set<std::string> parsePlainNodes(const std::string& text);

// The prepared mesh: library mm, Z-up, repaired; knit uv (the file's or the rebuilt surface's charts) and
// colour uv per corner; per-triangle material; per-vertex colour.
struct PreparedMesh {
    TriMesh mesh;                                // pos, tri, uv = KNITTING uv (per corner)
    std::vector<std::array<Vec2, 3>> colourUv;   // per corner (empty: none)
    std::vector<Vec3> vcolor;                    // per vertex (empty: none)
    std::vector<int> triMaterial;
    std::vector<ImportMaterial> materials;
    std::vector<std::string> textures;           // per material: copied texture file name in the import dir ("" none)
    bool chartsGenerated = false;
    std::vector<std::string> changes;            // what the repairs changed, one line each
    double scale = 1;                            // mm per file unit
    Vec3 offset;                                 // library = scale * zUp(file) + offset
    UpAxis up = UpAxis::Y;
};

// What one stored entry holds. plain has no triangles when every node is knitted.
struct PreparedImport {
    PreparedMesh knitted, plain;
};

// One preparation, as the import session hands it over: the file, the choices, the limits it was read under (its textures
// are read under them again) and the identity of the result, which names the stored entry.
struct ImportRequest {
    std::shared_ptr<const ImportedMesh> file;
    ImportOptions options;
    ImportLimits limits;
    uint64_t identity = 0;
};

// The file's vertices in library space: mm, Z up, centred over the origin with the lowest point at z = 0.
struct LibraryPlacement {
    std::vector<Vec3> pos;
    Vec3 offset;           // library = scale * zUp(file) + offset
    double diagonal = 0;   // of the bounding box [mm]
};
LibraryPlacement toLibrarySpace(const ImportedMesh& in, const ImportOptions& o);

bool prepareImport(const ImportedMesh& in, const ImportOptions& o, PreparedMesh& out, std::string* err);

// Library point -> file coordinates (OBJ export in the source's units and axes).
Vec3 toFileCoordinates(const PreparedMesh& p, const Vec3& lib);

// The import's identity: importer version, a frozen token, options, the source file's bytes and
// every dependency's bytes (MTL files, images) in reading order. Any change -> another identity.
uint64_t importIdentity(const ImportedMesh& m, const ImportOptions& o);

// Per-triangle base colour (linear RGB) sampled at the triangle's interior points (texture,
// vertex colours or Kd; default where none): what a part's "one colour per part" averages.
std::vector<Vec3> triangleColours(const PreparedMesh& p, const std::vector<ColourImage>& textures, const Vec3& fallback);

} // namespace knit
