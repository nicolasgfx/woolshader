#pragma once
// From reviewed results to an ExportScene. It computes no knit, no fibres and no colour: those are
// given. What it computes is the form exports want: tube meshes in metres.
#include <cstdint>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <vector>

#include "model/coordinates.h"
#include "export/export.h"
#include "export/export_scene.h"
#include "fibres/yarn_bake.h"
#include "scene/appearance/ply_maps.h"
#include "scene/appearance/stitch_colours.h"
#include "scene/appearance/yarn_tube.h"
#include "scene/scene.h"

namespace knit {

// A part's yarn tube (scene/appearance/yarn_tube.h) in the exported scene: Y-up metres, counter-clockwise
// seen along the normals.
std::optional<ExportTube> buildExportTube(const CurveSet& cs, const PatternParams& pp, const SceneMapping& map,
                                          const TubeStyle& style, std::string* err);
// The tube under the baked fibre maps (fibres/yarn_bake.h): `style`'s compression at kBakeRing and
// kBakeSides, u running 0..1 along each atlas strip of the yarn and v naming the strip (relabelAtlasUv).
std::optional<ExportTube> buildBakedExportTube(const CurveSet& cs, const PatternParams& pp, const SceneMapping& map,
                                               const TubeStyle& style, const BakeParams& bake, std::string* err);

// --ply-mode texture|geometry|none and the texture options, shared by --strand-demo and --pieces:
//   texture   (default) the single validated yarn tube, plies as maps;
//   geometry  explicit ply tubes (the reference);
//   none      the plain yarn tube, no maps.
//   --ply-shading normal|normal-flipg|bump   the frame map (normal-flipg: convention fixture)
//   --ply-cavity on|off, --ply-cavity-k k, --ply-map-d mm, --ply-map-size px
struct PlyModeOptions {
    enum class Mode { Texture, Geometry, None } mode = Mode::Texture;
    enum class Shading { Normal, NormalFlipGreen, Bump } shading = Shading::Normal;
    bool cavity = false;  // the AO map is not attached: renderer schema 12 has no AO (RDR-320)
    PlyMapSpec spec;  // g is set by the caller, sigma by setupPlyTextures
    bool textured() const { return mode == Mode::Texture && spec.g.count > 1; }
};

// Copies the maps a material can use from the cached set `r` into <packageDir>/textures (height, depth,
// both normal maps, and the cavity map with o.cavity; never the AO map, which only the cache and the
// viewer's preview hold), writes ply_maps.json for `sigma` and gives their shading. Empty with the reason.
std::optional<PlyShading> attachPlySet(const PlyMapSetResult& r, const PlyModeOptions& o, const std::string& packageDir, double sigma,
                                       std::string* err);
// Texture mode (o.spec.g set): the cached set of o.spec (ensurePlyMapSet), attached (attachPlySet), and its
// shading; the default shading unless o.textured(). The strand demo's way in (pieces_cli.cpp). Empty with
// the reason.
std::optional<PlyShading> setupPlyTextures(const PlyModeOptions& o, const std::string& packageDir, double sigma, std::string* err);

// The material id of a name: letters, digits and '_' (every other character becomes '_').
std::string materialIdOf(std::string id);
// The yarn material of a stitch colour (texture colourwork): "yarn_c<rrggbb>".
std::string stitchColourId(uint32_t key);

// Every material of the scene, and which one each part takes. collectScene calls addSceneMaterials
// first and asks for a material id whenever it adds a mesh or strands: it decides nothing about colour.
// The order of the entries is part of the output. Nothing here can fail.
class PartMaterials {
public:
    PartMaterials(const Scene& t, const SceneSettings& s, ExportScene& scene, std::string& log);

    // The yarn, filling, plain and fur entries of the scene, and the scene's fibre look.
    void addSceneMaterials();

    // A plain surface in the source's texture, or in its part's mean colour.
    std::string texturedPlainMaterial(const SourcePrimitive& g, const std::string& texturePath);
    std::string meanColourMaterial(const std::string& name, const Vec3& colour);
    // The fur of one palette colour, shared by the parts.
    std::string furColourMaterial(uint32_t colour);

    std::string yarnOf(int part) const;

    // Which material each part takes, as the collect stages read it.
    std::map<int, StitchColours> stitchColours;  // texture colourwork: part index -> its stitch colours
    std::map<int, uint32_t> dominant;            // part index -> its most frequent stitch colour
    std::map<uint32_t, int> colourIndex;         // palette colour -> index (sorted keys)
    std::vector<uint32_t> colourKeys;
    std::map<std::string, std::string> plain;    // source material -> material id
    std::map<int, std::string> fur;              // fur: part index -> its fur material id

private:
    void addObjectYarnMaterials();
    void addStitchColourMaterials();
    void addPlainMaterials();
    void addFurMaterials();
    std::string addPlainColour(const std::string& id, const Vec3& colour);

    const Scene& t_;
    const SceneSettings& s_;
    ExportScene& scene_;
    std::string& log_;
    std::map<int, std::string> yarn_;                 // part index -> its yarn material id
    std::set<uint32_t> furColours_;                   // fur in texture colours: the colours with a material so far
    std::map<std::string, std::string> texMaterial_;  // "<plain id>|<texture path>" -> material id
};

// The scene graph object of source primitive `prim`: its node (key "node<index>", or "prim<index>"
// without one), named after the node, else the mesh, else the index.
Placement exportObjectOf(const SourceScene& source, int prim, const std::string& group = "");

// The part the close-up cameras look at (closeup_<part>, mid, medium_<part>, far): `requested` if
// >= 0, else the knitted part of the largest mesh area (the first on a tie); -1 if none is knitted.
int closeupPartOf(const std::vector<ScenePart>& parts, int requested);

// Fills `scene` stage by stage, in the order the output lists them, sorts it and reads every strand
// file back. False with the reason.
[[nodiscard]] bool collectScene(const ExportSource& source, const ExportOptions& options, ExportScene& scene,
                                std::string& log, std::string* err);

} // namespace knit
