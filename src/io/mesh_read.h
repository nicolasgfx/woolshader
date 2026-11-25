#pragma once
// Reading foreign meshes (OBJ with its MTL, PLY) with bounds, for the viewer's Import step.
//
// Bounds (ImportLimits): the file size is checked before the file is
// opened; counts are checked before anything is allocated for them (PLY header counts, binary
// payload size against the file size) or while reading (OBJ, per statement). Knitting cost
// follows the surface area divided by the stitch area, not the triangle count: detail below a
// stitch (5.5 x 3.4 mm) is invisible in the knit, so 2 M triangles is far more than any
// knittable object needs.
//
// Colour: only the base colour is carried (MTL Kd and map_Kd, PLY vertex colours and a
// "comment TextureFile"); every other material property is listed as ignored. What happens to a mesh
// after reading is import/import_analysis.h (what it is) and import/prepare_import.h (what is made of it).
#include <array>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "io/image_io.h"
#include "model/vec.h"

namespace knit {

struct ImportLimits {
    uint64_t maxFileBytes = 200ull << 20;   // 200 MiB of mesh data: OBJ and PLY hold no textures
    int maxNodes = 200;                     // components (ImportNode): ImportSession refuses above it, after finishNodes
    size_t maxTriangles = 2000000;          // after triangulation
    size_t maxVertices = 2000000;
    int maxPolygon = 64;                    // corners of one face
    int maxTextureSide = 8192;
    uint64_t maxTextureBytes = 512ull << 20;  // decoded RGBA, all textures together
};

struct ImportMaterial {
    std::string name;
    Vec3 kd{0.8, 0.8, 0.8};
    bool hasKd = false;
    std::string textureRef;   // as written (map_Kd, TextureFile)
    std::string texture;      // resolved absolute path ("" if none or missing)
    bool textureMissing = false;
};

// One entry of the scene graph: a name the file gives (OBJ o or g), or a loose part of a file without names
// (finishNodes, import_analysis.h). The Import step calls it a component.
struct ImportNode {
    std::string name;
    long long triangles = 0;  // finishNodes
    double area = 0;          // finishNodes. File units^2: the stitch estimate of a selection
};

// What a file holds, in file units and file axes.
struct ImportedMesh {
    std::string path, format;                    // format: "OBJ", "PLY ascii", "PLY binary little endian", ...
    uint64_t fileBytes = 0;
    std::vector<Vec3> pos;
    std::vector<Vec3> vcolor;                    // per vertex in [0, 1] (empty: none)
    std::vector<std::array<int, 3>> tri;
    std::vector<std::array<Vec2, 3>> uv;         // per corner (empty: none)
    std::vector<int> triMaterial;                // per triangle, index into materials (-1: none)
    std::vector<int> quadFirstTriangle;          // the first of the two triangles of each face written with four corners
    std::vector<ImportMaterial> materials;
    std::vector<ImportNode> nodes;               // a name counts once a face uses it
    std::vector<int> nodeOfTriangle;             // per triangle, index into nodes
    int faces = 0, polygons = 0, maxPolygon = 0; // faces as written; those with more than 3 corners
    int facesWithoutUv = 0;                      // OBJ: faces without vt while others have them
    std::vector<std::string> ignored;            // statements / properties read but not used
    std::vector<std::string> dependencies;       // MTL files and images, resolved, in reading order
    std::vector<std::string> missingDependencies;
};

using ImportProgress = std::function<void(const std::string& stage, double fraction)>;

// Reads an OBJ or PLY file (by extension). false: a FATAL problem, with a message for the user
// (what is wrong, where, and what to do). Cancellable (cancel.h checkpoints).
bool readMeshFile(const std::string& path, const ImportLimits& lim, ImportedMesh& out, std::string* err,
                  const ImportProgress& progress = {});
// Decodes the material textures, within the texture limits. Missing or undecodable images are
// reported (warnings) and their materials fall back to Kd.
std::vector<ColourImage> loadImportTextures(ImportedMesh& m, const ImportLimits& lim, std::vector<std::string>* warnings);

std::string countText(uint64_t n);  // 1,234,567: the counts of messages to the user

} // namespace knit
