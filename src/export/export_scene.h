#pragma once
// ExportScene: the scene on its way out, the one thing every format reads. Records live here; geometry
// goes to a file as soon as it is added, because a whole scene does not fit in memory, and a format
// reads the parts back one at a time. Everything is in Y-up metres. See ARCHITECTURE.md.
#include <array>
#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "fibres/strands_file.h"
#include "io/image_io.h"
#include "model/mesh.h"
#include "model/coordinates.h"
#include "settings/material_spec.h"

namespace knit {

Vec3 yUp(const Vec3& p);

// A mesh as a part file holds it: x y z nx ny nz, and u v when uv is not empty.
struct SceneMesh {
    std::vector<Vec3> pos, nrm;
    std::vector<float> uv;
    std::vector<uint32_t> tri;
};

// A part's mesh read back: `stride` floats per vertex (6, or 8 with uv).
struct StoredMesh {
    std::vector<float> vertices;
    int stride = 6;
    std::vector<uint32_t> tri;
};

// The one place the .ply layout is written and read (binary little endian).
[[nodiscard]] bool writePartMesh(const std::string& path, const SceneMesh& mesh);
std::optional<StoredMesh> readPartMesh(const std::string& path, std::string* err);

// A yarn tube in scene units, wound for the renderer. yarn and s: per vertex.
struct ExportTube {
    SceneMesh mesh;
    std::vector<int> yarn;
    std::vector<double> s;
};

// A tube split by colour: every triangle names its group, written as part "<name><suffix>" with its
// material, in increasing group.
struct TubeGroup {
    std::string suffix, material;
};
struct TubeGroups {
    std::vector<int> groupOfTriangle;  // index into groups, per triangle of the tube
    std::vector<TubeGroup> groups;
};

// One file under meshes/ or strands/.
struct PartFile {
    std::string name, path, material;
};
// A part under the group of its source object (`object`: a stable key, `objectName` the group's
// display name), in a child group `group` of that object ("fibres_oriented", "fibres_random", "fur";
// empty: directly under the object). An empty object leaves the part at the root.
struct Placement {
    std::string object, objectName, group;
};
struct ExportCamera {
    std::string name;
    Vec3 eye, lookAt;
    double fov = 30;
};
struct Box {
    Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    void grow(const Vec3& p);
    void grow(const Box& b);
};

struct BaseMaterial {
    std::string id;
    MaterialSpec spec;
};
// The table a scene starts from: its first two entries (wool_white, wool_brown) are yarn colours.
std::vector<BaseMaterial> defaultMaterials();
MaterialSpec defaultRandomFibre();

// How fibre parts take their colour from their yarn.
struct FibreLook {
    bool hairMaterials = false;         // fibre parts are hair BSDFs (`hair`), not their yarn's diffuse + sheen
    MaterialSpec hair;
    MaterialSpec random = defaultRandomFibre();
    std::optional<MaterialSpec> orientedOverride;  // the oriented fibres' own material, in place of their yarn's
    bool randomInheritsYarn = false;    // random fibres take their yarn's colour as it stands
    bool hairKeepsYarnColour = false;   // a dense coat that keeps its yarn's colour
};

// A resolved material, ready for a format. PlyMappedYarn: a yarn colour under the ply maps. BakedYarn:
// one yarn's tube under its baked fibre maps (ExportScene::bakedMaps).
enum class MaterialKind { Surface, Hair, PlyMappedYarn, BakedYarn };
struct ExportMaterial {
    std::string id;
    MaterialSpec spec;
    MaterialKind kind = MaterialKind::Surface;
    std::string colourTexture;  // its base colour texture ("": none)
};

// What a package needs to shade its yarn tubes with a ply map set (export/collect_scene.h attaches one):
// the renderer uses one frame map per scope, and the tubes' u must repeat with the maps. Paths are relative
// to the package; empty: not attached. Default: no maps.
struct PlyShading {
    bool on = false;
    enum class Frame { Normal, Bump } frame = Frame::Normal;
    std::string normal, height, cavity;
    double bumpScale = 0;    // [scene units] height amplitude of H in [0, 1]
    double normalScale = 1;
    double uvRepeat = 2.5;   // [mm] of yarn per u
    int uvWrap = 0;
};

// What a package needs to shade one baked yarn tube (fibres/yarn_bake.h, a pilot): its maps, relative to
// the package, and the renderer's constants. The highlight's x axis follows the measured fibre direction:
// the renderer turns it by anisoRotation half-turns from the tube's u axis towards w x u, which on the
// tube is -theta, so a fibre angle of +psi from +s towards +theta is -psi / 180 here.
// Two materials (yarn_body.map_material), both with the highlight from the roughness and specular maps: fabric, the
// yarn's own fabric (its sheen in front of the highlight, its translucency) under the cavity and normal maps; shiny,
// a plastic with the highlight alone (opaque: the renderer has translucency on fabric only).
struct BakedMaps {
    bool fabric = true;           // false: shiny (plastic)
    std::string cavity, normal, specular, roughness, height;  // relative to the package (baked_maps.h attachBakedYarn)
    double specularScale = 0.35;  // the renderer's `specular`: F0 / (2 F0(1.5)), F0 of the yarn's fibre ior
    double roughnessAlong = 0.5;  // GGX roughness along the yarn over across it (the map scales both)
    double anisoRotation = 0;     // [half-turns]
};

// An import's own coordinates. Unset on any other scene: its files stay in Y-up metres, numbers untouched.
struct ImportFrame {
    double scale = 1;   // library = scale * zUp(file) + offset
    Vec3 offset;
    UpAxis up = UpAxis::Z;
};

class ExportScene {
public:
    explicit ExportScene(std::string folder);
    const std::string& folder() const { return folder_; }

    // What the scene says: set while it is collected, read by the formats.
    std::string projectName = "knitted_bears";
    std::vector<BaseMaterial> baseMaterials = defaultMaterials();  // the first yarnEntries are yarn colours
    int yarnEntries = 2;
    FibreLook fibreLook;
    std::map<std::string, std::string> albedoTextures;             // material id -> its colour texture
    PlyShading plyShading;
    std::map<std::string, BakedMaps> bakedMaps;                    // baked yarn material id -> its maps
    // The renderer's tip_scale on every wool fibre (part_fibre_shapes): the wool look's taper over the
    // whole fibre, or 1 where the export already carries a late taper (fibres.tail_taper), so the tip
    // is not narrowed twice.
    double fibreTipScale = 0.2;
    std::vector<ExportCamera> extraCameras;                        // after the default front camera
    double skyIntensity = 1.0;
    bool groundPlane = true;
    double groundAlbedo = 0.5;
    std::optional<ImportFrame> importFrame;

    // Parts: each writes its file at once.
    bool addMeshPart(const std::string& name, const std::string& material, const SceneMesh& mesh, std::string* err);
    bool addMesh(const std::string& name, const std::string& material, const TriMesh& sceneMesh, std::string* err);
    // cornerUv, one per triangle, OBJ convention (v up): the vertices split where their corners' uvs
    // differ (texture seams).
    bool addTexturedMesh(const std::string& name, const std::string& material, const TriMesh& sceneMesh,
                         const std::vector<std::array<Vec2, 3>>& cornerUv, std::string* err);
    // The part names written, in group order. A triangle without a group is refused.
    std::optional<std::vector<std::string>> addTubeGroups(const std::string& name, const ExportTube& tube,
                                                          const TubeGroups& groups, std::string* err);
    // An empty set is left out (the renderer rejects empty strand parts): true, nothing added. A set
    // that breaks the strand rules is refused.
    bool addStrands(const std::string& name, const std::string& material, const StrandSetF& strands, std::string* err);
    // The rules addStrands writes with (taperRules of the late taper; StrandRules{} until set).
    void setStrandRules(const StrandRules& rules) { strandRules_ = rules; }
    // A base colour texture as kColourTexturePrefix<name>.png: its path, for albedoTextures.
    std::optional<std::string> addColourTexture(const std::string& name, const ColourImage& image, std::string* err);
    // Placement of the part `name`, added before or after this call. A name that never becomes a
    // part is ignored (a strand part emptied by the budget).
    void place(const std::string& name, const Placement& at);
    // Fibre parts' materials, one per population and yarn colour, so the renderer edits each on its own.
    std::string orientedFibreMaterial(const std::string& yarnId);
    std::string randomFibreMaterial(const std::string& yarnId);
    // The yarn body's materials (fibres/yarn_body.h, a prototype), base entries in the yarn's colour, so
    // no wool tip taper is written for them. Fibres: hair with the yarn's fibre_* parameters, its colour
    // read as the dense mass the renderer's hair colour means. Core proxy: the yarn's colour, plain
    // diffuse (opaque), or the yarn's fabric with kBodyCoreTranslucency (translucent: the renderer has
    // translucency on fabric only).
    std::string bodyFibreMaterial(const std::string& yarnId);
    std::string bodyCoreMaterial(const std::string& yarnId, bool translucent);
    // One baked yarn's material `id`: the yarn's spec (its colour is the albedo under the cavity map) and the
    // maps the cache attached under the package (scene/appearance/baked_maps.h), kind BakedYarn.
    void addBakedYarn(const std::string& id, const MaterialSpec& yarnSpec, const BakedMaps& maps);
    // The last step of collecting: the same files whichever order the parts came in.
    void sortParts();

    // What the formats read.
    const std::vector<PartFile>& meshParts() const { return meshes_; }
    const std::vector<PartFile>& strandParts() const { return strands_; }
    const std::map<std::string, Placement>& placements() const { return placements_; }
    const std::map<std::string, Box>& partBounds() const { return partBounds_; }
    const Box& bounds() const { return bounds_; }
    size_t triangleCount() const { return triangles_; }
    size_t strandSegmentCount() const { return strandSegments_; }  // over every strand part
    std::string pathOf(const PartFile& part) const { return folder_ + "/" + part.path; }
    // A wool fibre part's material, never fur: the only parts that take the wool look's thin ends.
    bool isWoolFibreMaterial(const std::string& id) const { return fibreMaterials_.count(id) > 0; }
    // The resolved spec of a registered fibre material id; its model says hair (Fibre) or diffuse +
    // sheen. Oriented: orientedOverride, else hair in the yarn's colour (hairMaterials: the mean of a
    // dense mass and one fibre of it), else the yarn entry. Random: its own hair model, else hair in
    // its colour (hairMaterials), else `random`; with randomInheritsYarn, in the yarn's colour as it
    // stands (hair: at least the renderer's floor). Refused if the id is not registered or a yarn it
    // needs is not a base material.
    std::optional<MaterialSpec> fibreMaterialSpec(const std::string& id, std::string* err) const;
    // The materials the parts reference, resolved, then the ground's: the base table, the ply-mapped
    // yarns, the fibre parts' materials. An entry left out can still colour a variant or a fibre. Call
    // it after the last part is added.
    std::optional<std::vector<ExportMaterial>> materials(std::string* err) const;

    static std::string yarnMaterial(const std::string& entry) { return entry + "_yarn"; }
    // The strand files of a strand table (first: per strand, plus the end; vertices): file k holds the
    // strands [cuts[k], cuts[k + 1]), whole strands of at most maxSegments segments together (a longer
    // strand alone). cuts = {0, N} when everything fits one file.
    static std::vector<size_t> strandFileCuts(const std::vector<uint64_t>& first, long long maxSegments);
    static constexpr const char* kOrientedGroup = "fibres_oriented";
    static constexpr const char* kRandomGroup = "fibres_random";
    static constexpr const char* kBodyGroup = "yarn_body";
    // The most segments one yarn-body strand file holds: the renderer builds one acceleration structure per
    // file and its memory peak follows the largest (the swatch's 35.7 M segments fit a 12 GB card only in
    // files of 2.5 M, 2026-11-07).
    static constexpr long long kBodyFileSegments = 2'500'000;
    static constexpr double kBodyCoreTranslucency = 0.3;     // as the soft wool presets' yarn
    static constexpr const char* kColourTexturePrefix = "textures/colour_";
    static constexpr const char* kGround = "_ground";         // the ground's mesh and material id
    static constexpr const char* kGroundObject = "~ground";   // sorts after source objects: its group last
    static constexpr double kGroundHalfSize = 3.0;            // ground half-size, in diagonals of the bounds
    static constexpr double kRandomFibreGrey = 0.55;          // the random fibres' default colour

private:
    struct FibreMaterial {
        std::string yarn;
        bool random = false;
    };

    // The one way a mesh enters the scene: its file, then the bookkeeping.
    bool addMeshFile(const std::string& name, const std::string& material, const Box& box, const SceneMesh& mesh,
                     std::string* err);

    std::string folder_;
    std::vector<PartFile> meshes_, strands_;
    std::map<std::string, Placement> placements_;
    std::map<std::string, Box> partBounds_;
    std::map<std::string, FibreMaterial> fibreMaterials_;
    Box bounds_;
    size_t triangles_ = 0, strandSegments_ = 0;
    StrandRules strandRules_;
};

// The default front camera, in front of the parts (Blender's front view looks along +Y; in Y-up that
// is -Z). Refused if there is nothing to frame or a camera is not finite: the renderer's JSON reader
// rejects inf and nan.
struct SceneFrame {
    Vec3 centre, eye;
    double radius = 0, fov = 30;
};
std::optional<SceneFrame> frameScene(const ExportScene& scene, std::string* err);

// A JSON string body: quotes, backslashes and control characters escaped.
std::string jsonEscaped(const std::string& text);

// --- the formats. Each reads an ExportScene; none sees how it was made.

// scene.json in the scene's folder, the renderer's scene format (format_ppt.cpp).
[[nodiscard]] bool writeSceneJson(const ExportScene& scene, std::string* err);
// A material's colour as scene.json rounds it, which the OBJ's Kd has always followed.
Vec3 colourAsWritten(const ExportMaterial& material);

// Where a format writes its own files.
struct FormatTarget {
    std::string folder;
    std::string baseName;  // the stem of the format's files, never empty
};

// <baseName>.obj and <baseName>.mtl: every part, fibres as polylines (format_obj.cpp).
[[nodiscard]] bool writeObjFile(const ExportScene& scene, const FormatTarget& target, std::string& log,
                                std::string* err);

// Blender and its conversion script are both here.
bool blenderConversionAvailable();
// <baseName>.blend from a package folder, not a scene: Blender reads the files. Throws Cancelled when
// the job is cancelled (Blender is stopped first). Blender's own message becomes *err on a failure.
[[nodiscard]] bool writeBlendFile(const std::string& packageFolder, const FormatTarget& target, std::string& log,
                                  std::string* err);

} // namespace knit
