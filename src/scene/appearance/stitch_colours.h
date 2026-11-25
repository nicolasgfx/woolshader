#pragma once
// Texture colourwork (material.yarn_colour_source = texture): the colour of every stitch of a knitted
// part, the object's base colour texture averaged over the stitch's cell.
//
// A stitch's cell is origin + a w U + b h V, a in [-1/2, 1/2], b in [0, 1], in its layout frame
// (StitchRecord; w, h the knit's wale and course spacing). kCellSamples x kCellSamples points of the
// cell are projected to the part's surface and read bilinearly, and the mean is taken in linear
// light. Sampling on the surface, not in a chart, also serves stitches past a chart's edge (overlapped
// seams) and gap patches. The renderer has one colour per material, so the
// scene's colours are reduced to at most kYarnPaletteColours yarn colours; the viewer shows the same
// ones. Computed from the knit: the knit, its cache and its keys are unaffected.
//
// sceneStitchColours is read by the knit display, the tube rebuild, the fibre lines and the export
// materials. The colours are made by the first reader and kept on the Scene (Scene::stitchColours):
// always after knitScene, restoreKnit or loadBaseParts returned, so the parts are complete and nothing
// invalidates them. Lazy because colourwork is an Export setting (switching it on must not knit again);
// a lock, because two builders in one process can reach one scene (jobs/job_steps.cpp). fur_stage.cpp samples a hair's root with
// partSurfaceColour. Nothing here fails: a stitch that finds no surface takes its part's mean colour and
// is counted, and a knit without a texture has an empty palette. The texture must have pixels.
#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "io/image_io.h"
#include "model/vec.h"

namespace knit {

struct SourcePrimitive;
struct SceneSettings;

struct Scene;
struct ScenePart;
struct ColourImage;
struct CurveSet;
struct PlySet;

constexpr int kCellSamples = 4;          // per side of a stitch's cell
constexpr int kYarnPaletteColours = 64;  // per scene

struct StitchColours {
    std::vector<uint32_t> key;  // per stitch (id order): 0xRRGGBB, sRGB bytes; empty: the part has no colours
};

struct SceneStitchColours {
    std::vector<StitchColours> parts;  // one per Scene::parts entry
    std::vector<uint32_t> palette;     // the distinct colours, sorted; empty: no knitted part has a texture
    int stitchesOffSurface = 0;        // stitches with no cell sample near the surface: their part's mean colour
};

// The area-weighted mean base colour of a primitive's triangles (null: all of them); false if it has
// no colour per triangle.
bool meanTriangleColour(const SourcePrimitive& g, const std::vector<int>* tris, Vec3& out);
// The yarn colour of a knitted part: material.yarn_white, or with material.yarn_colour_source =
// object the area-weighted mean base colour of the part's source triangles (SourcePrimitive::triColor;
// the primitive's base colour where there is none). One colour per part: no colourwork.
Vec3 partYarnColour(const Scene& t, const ScenePart& r, const SceneSettings& s);
// sRGB bytes 0xRRGGBB -> linear colour (the material colour space).
Vec3 stitchColourLinear(uint32_t key);
// Linear colour -> sRGB bytes 0xRRGGBB.
uint32_t stitchColourKey(const Vec3& linear);
// At most maxColours colours for the weighted colours `stitchesOfColour` (key -> stitches): median cut
// (Heckbert 1982), always splitting the box with the widest channel, so a small distinct colour keeps
// its own entry; an entry is the stitch-weighted linear mean of its box. Every key maps to its entry;
// with at most maxColours colours each maps to itself.
std::map<uint32_t, uint32_t> reducePalette(const std::map<uint32_t, long long>& stitchesOfColour, int maxColours);
// Whether part r's surface has a base colour texture and a uv on its triangles.
bool partSurfaceHasTexture(const Scene& t, const ScenePart& r);
// The texture's colour (linear, read bilinearly, as a stitch reads it) at barycentric weights `bary`
// of triangle `tri` of part r's mesh. Needs partSurfaceHasTexture.
Vec3 partSurfaceColour(const Scene& t, const ScenePart& r, int tri, const double bary[3]);
// The colours of every knitted part's stitches. A part without a texture takes its object colour
// (partYarnColour); with no texture on any knitted part the result is empty.
const SceneStitchColours& sceneStitchColours(const Scene& t);

// The stitch on a rendered strand at arc length `s`: the stitch owning it, or before the yarn's first
// / after its last stitch (the run's padding) that stitch; -1: the yarn has no stitch. A strand is a yarn of cs, or with `plies` (separate ply tubes) a ply of
// that set, mapped to its parent yarn.
class StitchAlongYarn {
public:
    explicit StitchAlongYarn(const CurveSet& cs, const PlySet* plies = nullptr);
    int at(int strand, double s) const;

private:
    const CurveSet& cs_;
    const PlySet* plies_;
};

} // namespace knit
