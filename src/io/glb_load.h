#pragma once
// S6: GLB scenes into library space.
//
// Units contract: ONE scene-wide similarity T is applied to every mesh (knitted, plain,
// socket, shells): glTF world (Y-up, node transforms applied) -> library Z-up by
// (x, y, z) -> (x, -z, y), then a uniform scale s [mm per glTF unit] chosen so that the
// registry size holds in mm (GlbSelection::sizeRef over ALL included meshes), then a
// translation placing the bbox centre over the origin (x, y) and the lowest point at
// z = 0 (the sphere scene's nodes put it ~20 m away). Library units are mm; the exporter
// maps mm -> metres with
// SceneMapping{sigma = 1000, centre = component centroid / 1000} and back to Y-up, so a
// 300 mm figure exports 0.300 m tall.
//
// Welding is per source (node, primitive): coincident vertices of different sources are
// never merged. Corner uvs are kept per triangle. Plain meshes keep their source corner
// normals (hard edges); knitted meshes use welded area-weighted normals (frames only).
// Mirrored nodes (det < 0) get their winding flipped. Triangles only.
#include <array>
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "io/glb_selection.h"
#include "io/image_io.h"
#include "model/mesh.h"
#include "model/mesh_check.h"

namespace knit {


// glTF KHR_texture_transform of a texture reference (the identity without the extension).
struct TextureTransform {
    Vec2 offset{0, 0};
    double rotation = 0;  // [rad]
    Vec2 scale{1, 1};
};
// The texture coordinate a mesh's uv (glTF's, v down) reads: scale, then rotation, then offset (the
// extension's matrix T R S).
Vec2 transformTextureUv(const TextureTransform& t, const Vec2& uv);

struct SourcePrimitive {
    std::string node, mesh, material;
    int primitive = 0;
    // The source node (glTF node index; an import's one object: 0; -1 unknown): the object a part
    // is grouped under in the export. Names are not identities (empty, or one mesh on several nodes).
    int nodeIndex = -1;
    enum class Class { Knit, Plain, Skip } cls = Class::Knit;
    Vec3 baseColor{0.8, 0.8, 0.8};
    std::string baseColorImage;    // uri or name ("" if none)
    TriMesh welded;                // library mm, Z-up; per-corner uv (empty if none)
    std::vector<std::array<Vec3, 3>> cornerNormals;  // library Z-up, per welded triangle (plain export)
    bool mirrored = false;
    int droppedTriangles = 0;      // degenerate after welding
    MeshCheck check;               // per component (checkMesh)
    // Base colour for display and the per-part yarn colour (imports, prepare_import.h; GLB: the
    // decoded baseColorTexture). colourUv per welded triangle, OBJ convention (v up), with the
    // texture's KHR_texture_transform applied (empty: welded.uv places the colour); triColor: the base colour per welded triangle (empty: baseColor everywhere);
    // vertexColor per welded vertex (empty: none).
    std::shared_ptr<const ColourImage> albedo;
    std::vector<std::array<Vec2, 3>> colourUv;
    std::vector<Vec3> triColor;
    std::vector<Vec3> vertexColor;
};

struct SourceScene {
    std::string path;
    uint64_t fileHash = 0;         // FNV-1a 64 of the file bytes
    double scale = 0;              // s [mm per glTF unit]
    Vec3 offset;                   // library = s * zUp(glTF world) + offset [mm]
    Vec3 lo, hi;                   // bounds of all included meshes, library mm
    std::vector<SourcePrimitive> prims;
    std::vector<std::string> log;  // ignored extensions, classification, statistics
    // An imported mesh (import_store.h): library = scale * zUp(file) + offset with the file's up axis.
    bool imported = false, importYUp = true;
};

struct SourceOptions {
    bool selfIntersection = true;  // run checkMesh's triangle-pair test (cost ~ triangles)
    double weldRel = 1e-6;         // weld eps relative to the primitive's bbox diagonal
    bool requireCharts = true;     // imports: refuse a mesh without knitting charts (false: fur)
};

bool loadGlb(const std::string& path, const GlbSelection& selection, SourceScene& out, std::string* err, const SourceOptions& o = {});

} // namespace knit
