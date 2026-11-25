#pragma once
// How the jobs build a BuiltScene (built_scene.h) from the knit and fibre results: private to this
// folder. Functions return false with the reason when a part cannot be drawn; a cancelled job throws
// Cancelled at a checkpoint.
#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include "fibres/fibre.h"
#include "io/glb_load.h"
#include "jobs/built_scene.h"
#include "scene/appearance/baked_maps.h"
#include "settings/scene_settings.h"
#include "scene/fibres/scene_fibres.h"
#include "scene/scene.h"
#include "scene/appearance/rendered_part.h"
#include "scene/appearance/stitch_colours.h"

namespace knit {

struct RenderedPart;

// At most this many fibre segments are drawn, the smallest keys first (displayThreshold); the export writes all.
constexpr long long kDisplaySegments = 2000000;
// The display tubes' triangle budget (yarn tubes of one scene; display only, never the export).
constexpr long long kDisplayTriangles = 20000000;

struct ColourImage;  // io/image_io.h
// A display copy of a base-colour texture, at most 2048 texels on a side (box filter).
std::shared_ptr<const DisplayTexture> makeDisplayTexture(const ColourImage& t);

// The colour of a mesh that stands in for a part still being knitted.
inline const Vec3 kDim{0.55, 0.55, 0.55};

// How a source surface is drawn: its mesh kind, a tint on its colours, and how far it is moved along -N.
struct SourceLook {
    DisplayMesh::Kind kind;
    Vec3 tint;
    double inset;  // [mm]
};
// Plain primitives and unknitted parts, as the source looks.
constexpr SourceLook kPlainLook{DisplayMesh::Kind::Plain, {1, 1, 1}, 0};
// A knitted part's own surface (the viewer's "Mesh" view), as the source looks.
constexpr SourceLook kSourceLook{DisplayMesh::Kind::Base, {1, 1, 1}, 0};
// A stand-in for a part still being knitted: dimmed, `inset` behind the fabric (fillingOffset).
inline SourceLook standInLook(double inset) { return {DisplayMesh::Kind::Base, kDim, inset}; }

// The partial scene once a scene's parts are known: its plain primitives (fixed), one stand-in per part
// (base), and the textures those meshes index. Made afresh on every call, so a second part list (gap
// patches) indexes its own textures.
struct KnitStandIns {
    std::vector<DisplayMesh> fixed, base;
    std::vector<std::shared_ptr<const DisplayTexture>> albedo;
};
KnitStandIns knitStandIns(const Scene& scene, double inset);

// The object's surface as the source looks: per corner (colour uv, vertex colours), the base colour texture when
// there is one, drawn as `look` says. triangles: the primitive's triangles to draw (null: all). Each texture goes
// into out.albedo once; textureSlots remembers where.
void addSourceMesh(BuiltScene& out, std::map<const ColourImage*, int>& textureSlots, const std::string& name, const SourcePrimitive& g,
                   const std::vector<int>* triangles, const SourceLook& look);

// The texture profile's preview at yarn diameter partD: the AO map of the cached map set the export uses (generated
// on a miss), with its u period in TexturePreview::pitch. Null: no texture profile applies, so the caller keeps its
// texture. Empty (no pixels, pitch 0): the map set failed, and the shown texture is replaced by nothing (the viewer
// treats an empty image as no texture).
std::shared_ptr<const TexturePreview> texturePreview(const SceneSettings& s, double partD);

// How this result draws a knit's tubes, one value for the Knit job and the strand mesh rebuild: the ring spacing
// (the strand tessellation's, or coarser where the scene's tubes would exceed kDisplayTriangles, with a note saying
// so; the export keeps the setting) and the texture preview for the parts of yarn diameter textureD.
struct KnitTubeLook {
    double ring = 0;                                // [mm]
    std::string displayNote;                        // "": the ring is the setting's
    double textureD = 0;                            // [mm] parts of this yarn diameter carry the texture uvs
    std::shared_ptr<const TexturePreview> texture;  // null: no texture profile
};
KnitTubeLook knitTubeLook(const Scene& t, const SceneSettings& s, double textureD);

// The yarn tube of one knitted part, pushed onto out.meshes as exactly one mesh: its strands (the plies when it has
// them) at the look's ring spacing, at the part's centre, with its object colour and its stitch colours per vertex.
// out.yarnLength grows by the strands' lengths yarn by yarn, in this order: returning the length instead would change
// the order of the floating-point additions. The callers keep their own cancel checkpoints and progress. False, with
// the reason, when the part cannot be rendered.
[[nodiscard]] bool addKnittedPartTube(BuiltScene& out, const Scene& t, const ScenePart& part, const SceneSettings& s,
                                      const KnitTubeLook& look, std::string* err);
// A knitted part's tubes as the settings render the yarn: the plain tube, or for a yarn of fibres rendered as
// maps the per-yarn tubes under their cached cavity maps (below); the plain tube and a display note while the
// cache holds no maps yet (the Knit step's Bake the yarn maps, a render preview or an export bakes them).
[[nodiscard]] bool addKnittedPartTubes(BuiltScene& out, const Scene& t, const ScenePart& part, const SceneSettings& s,
                                       const KnitTubeLook& look, std::string* err);
// A yarn of fibres rendered as maps: the part's tube as one mesh per yarn (`<part>#yNNN`), each under its cavity
// map from the baked entry as the viewer's albedo (the yarn colour multiplies it), the atlas uv the export
// uses; the display's tessellation is the look's (the maps fit any tessellation of the same curves).
[[nodiscard]] bool addKnittedPartMappedTubes(BuiltScene& out, const Scene& t, const ScenePart& part, const SceneSettings& s,
                                             const KnitTubeLook& look, const RenderedPart& rendered, const BakedMapEntry& entry,
                                             std::string* err);
// Puts rebuilt tubes (one mesh per part, or per yarn under the maps) where the scene's yarn meshes were, keeping
// each part's colour choices; the tubes' textures are appended. Every part with yarn meshes must get tubes again.
[[nodiscard]] bool replaceYarnMeshes(BuiltScene& out, BuiltScene& tubes, std::string* err);
// A scene restored from the scene cache under a yarn of fibres rendered as maps: its tubes come up under the
// cached maps (the stored scene never holds the maps' pixels, withoutYarnMaps). Nothing to do otherwise.
[[nodiscard]] bool showCachedYarnMaps(BuiltScene& c, const Scene& t, const SceneSettings& s, std::string* err);
// The scene as the scene cache stores it: the yarn tubes without their maps (the bake cache holds the pixels;
// the other meshes keep their textures and indices).
BuiltScene withoutYarnMaps(const BuiltScene& b);
// A knitted part's coarse tube for a partial scene (rings at least 0.8 mm apart, at most 6 sides, no colours; the
// result gets the strand tessellation). nullopt when the part cannot be rendered.
std::optional<DisplayMesh> coarsePartTube(const ScenePart& part, const SceneSettings& s);
// The part's object colour (material.yarn_colour_source = object), whatever the current colour source: the viewer
// picks it or the slot colour from the current settings.
Vec3 objectColourOf(const Scene& t, const ScenePart& part, const SceneSettings& s);

// One fibre as GL_LINES pairs from its emergence point (a polyline's buried root segment is skipped), moved by
// `offset`, with its population (0 along the yarn, 1 stray, 2 fur) and its colours per vertex; records its first
// vertex for arrangeFibreLines.
void addFibreLine(DisplayLines& out, std::span<const Vec3f> points, const Vec3& offset, uint8_t population, const Vec3& partColour,
                  const Vec3& ownColour);

// An axis-aligned box [scene mm], empty until something is added (lo above hi).
struct DisplayBounds {
    Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    bool empty() const { return lo.x > hi.x; }
};
void addToBounds(DisplayBounds& box, const DisplayMesh& d);
// The scene's box from all its meshes; an empty box leaves b.lo, b.hi as they were.
void updateSceneBounds(BuiltScene& b);

// The display's key threshold: at most kDisplaySegments segments drawn, the smallest keys first.
float displayThreshold(const SceneFibres& tf);

// The yarn diameter most of the knitted rows have (rows with stitches and a diameter): the texture preview's D.
double mostCommonD(const std::vector<PartRow>& parts);

// The display rows and tubes of a tiled knit (shared by a computed knit and a restored one).
bool knitDisplay(BuiltScene& out, const SceneEntry& e, const Scene& t, const SceneSettings& s, std::string* err);

} // namespace knit
