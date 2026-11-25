#pragma once
// BuiltScene: what the viewer draws and shows after a job, as plain CPU arrays.
//
// Built from the knit and fibre results by the jobs in this folder (display_scene.h). The viewer's
// main thread is the only place that uploads it to GL, so nothing here touches GL.
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "jobs/results.h"
#include "model/params.h"
#include "model/status.h"
#include "model/vec.h"
#include "settings/scene_registry.h"
#include "settings/scene_settings.h"

namespace knit {

struct DisplayMesh {
    std::string name;
    std::vector<float> pos, nrm;     // scene mm
    std::vector<uint32_t> idx;
    std::vector<float> uv;           // yarn with a texture profile: 2 per vertex (u = s / pitch, v = theta / 2 pi)
    float color[3] = {0.7f, 0.7f, 0.7f};
    enum class Kind { Base, Yarn, Plain = 3 } kind = Kind::Base;  // the number is stored in scene.bin
    // Yarn of a tiled part: its object colour (material.yarn_colour_source = object; the viewer
    // picks it or the slot colour from the CURRENT settings).
    bool hasPartColor = false;
    float partColor[3] = {0.7f, 0.7f, 0.7f};
    // Base / plain meshes with the object's base colour (imports, GLB textures): per-vertex
    // colour (empty: `color`) and a texture with per-vertex uv (the viewer uploads `albedo`).
    std::vector<float> vcol;         // 3 per vertex
    std::vector<float> albedoUv;     // 2 per vertex, image space (0, 0: the texture's first row, left)
    int albedo = -1;                 // index into BuiltScene::albedo
};

// A spatially compact run of fibre lines (arrangeFibreLines): the viewer skips a chunk outside
// the view, and draws a prefix of it while the camera moves (its fibres are in random order, so
// any prefix is an even thinning).
struct FibreChunk {
    uint32_t firstVertex = 0, vertexCount = 0;  // even: whole segments
    Vec3 lo, hi;                                // bounds [scene mm]
};

// Fibre lines: GL_LINES pairs, scene mm; per vertex the population (0 along the yarn, 1 stray, 2 fur) and
// its part's object colour, so the viewer colours them from the CURRENT materials.
struct DisplayLines {
    std::vector<float> pos;
    std::vector<uint8_t> population;
    std::vector<float> partColor;           // 3 per vertex: the part's object colour (linear)
    std::vector<float> ownColor;            // 3 per vertex: the fibre's stitch or root colour, else the part's
    std::vector<uint32_t> fibreStarts;      // while collecting: each fibre's first vertex (emptied by arrangeFibreLines)
    std::vector<FibreChunk> chunks;         // empty: one run, never culled (lines not arranged)
};

// Reorders the collected fibre lines into chunks (FibreChunk); every line is kept.
void arrangeFibreLines(DisplayLines& lines);

// A decoded base-colour texture (sRGB RGBA, rows top first) for the viewer.
struct DisplayTexture {
    std::string name;
    int width = 0, height = 0;
    std::vector<uint8_t> rgba;
};
struct PartRow {
    std::string name;
    bool knitted = false;
    Status verdict = Status::Unresolved;
    int stitches = 0;
    double yarnD = 0;
    std::string note;
    std::vector<std::string> lines;  // the part's report lines
    Vec3 lo, hi;                     // display bounds (scene mm; focus)
};

// The texture profile preview: the AO map of the export's cached map set, shown
// on the yarn tubes' uvs (the viewer multiplies the yarn colour by it).
struct TexturePreview {
    int width = 0, height = 0;
    std::vector<uint8_t> rgb;        // rows top first (the PNG); empty: the map set could not be made
    double pitch = 0;                // [mm] the map's u period (the ply pitch); 0: none
};

// Whether the scene's yarn tubes carry their baked cavity maps (a yarn mesh with an albedo texture): the
// Display popover's Yarn maps row and the Knit step's Bake button ask; such a scene is stored without them.
struct BuiltScene;
bool carriesYarnMaps(const BuiltScene& b);

struct BuiltScene {
    enum class Kind { Knit, Fibres, Export };
    Kind kind = Kind::Knit;
    int revision = 0;
    bool partial = false;  // published while the job runs: the parts finished so far
    bool cached = false;   // read from the persistent cache (scene_cache.h), not computed
    SceneSettings settings;
    std::vector<DisplayMesh> meshes;
    std::vector<std::shared_ptr<const DisplayTexture>> albedo;  // base-colour textures of the meshes
    DisplayLines fibres;
    std::vector<PartRow> parts;
    std::shared_ptr<const TexturePreview> texture;  // null: no texture profile (or explicit plies)
    std::string report, exportLog;
    // The display's stated deviation from the strand tessellation (strand.tube_*): above
    // the display budget (kDisplayTriangles, display_scene.h) the tubes take a coarser ring spacing (the export keeps the setting).
    std::string displayNote;
    Vec3 lo{-100, -100, 0}, hi{100, 100, 200};
    // [mm] the curves the yarn tubes follow, over the knitted parts: the plies when the knit has them (their
    // helix included), else the yarn. The Export step forecasts the tube triangles from it.
    double yarnLength = 0;
    double ms = 0;
    // The job's result: Knit -> knit (null while partial); Fibres -> fibreResult (and the knit
    // it used); Export -> both used, and the published folder.
    std::shared_ptr<const KnitResult> knit;
    std::shared_ptr<const FibreResult> fibreResult;
    std::string exportedTo;
};

} // namespace knit
