#include "settings/scene_registry.h"

#include <algorithm>
#include <filesystem>
#include <map>
#include <memory>
#include <mutex>

#include "io/cache_paths.h"
#include "io/files.h"

namespace knit {

namespace {
// Sets key in a defaultSettings list: its value replaced where the list has it, else appended.
void setListedSetting(std::vector<std::pair<std::string, std::string>>& list, const std::string& key, const std::string& value) {
    for (std::pair<std::string, std::string>& kv : list)
        if (kv.first == key) {
            kv.second = value;
            return;
        }
    list.push_back({key, value});
}
} // namespace

const std::vector<SceneEntry>& sceneRegistry() {
    static const std::vector<SceneEntry> r = [] {
        std::vector<SceneEntry> v;
        {
            // The two teddy bears (user, 2026-10-28: every object goes through the same workflow).
            // "Teddy bears" by hectopod, https://sketchfab.com/3d-models/teddy-bears-e84b12b4ac20402aaf4d40f2219cd0e2,
            // CC-BY-4.0 (the credit below; the model is not distributed with the source).
            // Knitted, they are the spiral pieces (scenes/teddy-bears/pieces, made by the scene
            // preparation tools): the UV charts of bears.glb knit with holes at their
            // borders (0 of 18 PASS against 14 of 18 as pieces). With the Knit
            // step skipped they are bears.glb (converted from
            // source/bears.fbx: posed, the fabric subdivided, the yarn colours as the fabrics' base
            // colours). Size: 3.8146 scene units high x 68.1 mm per unit (the calibrated gauge, 65 stitches
            // of 5.5 mm around the body), the same as the pieces.
            SceneEntry e;
            e.id = "teddy";
            e.label = "Teddy bears (26 cm)";
            e.projectName = "knitted_bears";
            e.kind = SceneEntry::Kind::Glb;
            e.path = "scenes/teddy-bears/bears.glb";
            e.glb.sizeRef = GlbSelection::SizeRef::Height;
            e.glb.sizeMm = 259.8;
            // The mesh selection, as an import's: the fabrics knitted, the eyes, noses, mouth threads and the bow kept as
            // the file has them (knitted, the bears are the pieces below; the kept nodes are placed with them).
            e.glb.knitAllMaterials = true;
            e.glb.plainNodes = {"white_bear_eyes", "brown_bear_eyes", "white_bear_nose", "brown_bear_nose", "white_bear_mouth", "brown_bear_mouth", "bow"};
            e.defaultYarnColourSource = "object";
            // Natural islands; density 135, the highest at which neither bear reaches the Reduced
            // export's 7 M segments (owner, 2026-11-05: Handmade Dense knit, Soft wool, 6.88 M and
            // 6.68 M; at 140 bear.002 is cut). Fuller coats: one bear with a larger budget, by hand.
            e.defaultFibreDensity = 135;
            e.defaultFibreIslands = "natural";
            e.piecesManifest = "scenes/teddy-bears/pieces/manifest.txt";
            e.gaugePiece = "piece_bear_002_00.obj";  // the body of the rest-pose bear
            e.gaugeStitches = 65;
            e.piecesW = 5.5;
            e.piecesH = 3.4;
            // The Dense preset: drawn 1.86, knitted compressed at 1.16, the thickest round yarn that
            // passes gates (a) and (d) at 65 stitches around the body (the Thin preset).
            e.piecesD = 1.86;
            e.piecesDFallback = {1.80, 1.73, 1.67, 1.60};
            e.pieceColourMaterials = {{"wool_white", "fabric_white"}, {"wool_brown", "fabric_brown"}};
            e.credit = {"Teddy bears", "hectopod", "https://sketchfab.com/3d-models/teddy-bears-e84b12b4ac20402aaf4d40f2219cd0e2"};
            v.push_back(e);
        }
        {
            // The brown bear alone, a test model for both workflows (owner, 2026-11-19: model only, no workflow
            // defaults). The white bear and the bow are left out of the file (glb.skipMeshes: their glTF meshes);
            // the brown bear's eyes, nose and mouth stay plain meshes (no fur or yarn on them; they block the fur).
            // It is the file's full height. Knitted, it is the brown bear's spiral pieces: manifest_brown.txt, the
            // brown bear's lines of manifest.txt; its body piece is the same mesh as the white
            // bear's gauge piece (area 10.032129 against 10.032136 scene units), so 65 stitches give the pair's scale.
            SceneEntry e;
            e.id = "teddy-brown";
            e.label = "Teddy, one bear (26 cm)";
            e.projectName = "teddy_brown";
            e.kind = SceneEntry::Kind::Glb;
            e.path = "scenes/teddy-bears/bears.glb";
            e.glb.sizeRef = GlbSelection::SizeRef::Height;
            e.glb.sizeMm = 259.8;
            e.glb.knitAllMaterials = true;
            e.glb.plainNodes = {"brown_bear_eyes", "brown_bear_nose", "brown_bear_mouth"};
            e.glb.skipMeshes = {"Sphere.001", "bow_bear.002", "Sphere.003", "Sphere.004", "Sphere.005"};  // the white bear and the bow
            e.piecesManifest = "scenes/teddy-bears/pieces/manifest_brown.txt";
            e.gaugePiece = "piece_bear_009_00.obj";
            e.gaugeStitches = 65;
            e.piecesW = 5.5;
            e.piecesH = 3.4;
            e.piecesD = 1.86;
            e.piecesDFallback = {1.80, 1.73, 1.67, 1.60};
            e.pieceColourMaterials = {{"wool_brown", "fabric_brown"}};
            e.credit = {"Teddy bears", "hectopod", "https://sketchfab.com/3d-models/teddy-bears-e84b12b4ac20402aaf4d40f2219cd0e2"};
            v.push_back(e);
        }
        {
            // The wool swatch on a needle (owner, 2026-11-11: an ordinary preset, knitted in the app like every
            // other). scenes/wool-swatch/wool_swatch.glb (the author's own model, scenes/wool-swatch/LICENSE.md) holds
            // the open-edge plane of the swatch (standing, wales up) and its needle as a plain part (clear of the yarn
            // and inside the plane's box, so the knit sees the plane's own frame). Knitted on the plane's own uv chart
            // with the settings below: one yarn, 528 stitches in 24 runs (verdict FAIL: gate (a) and B1). The app adds
            // no loops around the needle, no working thread and no filler; the rendered swatch scene adds those.
            SceneEntry e;
            e.id = "wool-swatch";
            e.label = "Wool swatch on a needle (12 cm)";
            e.projectName = "wool_swatch";
            e.kind = SceneEntry::Kind::Glb;
            e.path = "scenes/wool-swatch/wool_swatch.glb";
            e.glb.sizeRef = GlbSelection::SizeRef::Height;
            e.glb.sizeMm = 104.993824;  // the height of plane and needle (export_swatch_glb.py prints it): 1 mm stays 1 mm
            // The mesh selection, as an import's: the swatch knitted, the needle kept as the file has it.
            e.glb.knitAllMaterials = true;
            e.glb.plainNodes = {"needle"};
            // Each part its object's colour: the GLB's yarn is the white scene's stored albedo (0.45, 0.44, 0.42).
            e.defaultYarnColourSource = "object";
            e.defaultFibres = true;
            // The swatch scene's knit settings with its cream fibres and materials, left out: the yarn body (owner,
            // 2026-11-11: the app's wool fibres), material.yarn_white.color (the colour is the GLB's),
            // filling.enabled, strand.diameter_fallback (the profile's) and everything scene-level (props, filler,
            // cameras, light). strand.stitch_size comes first (applySceneSettings). The blue swatch
            // below is this entry with the yarn made of fibres.
            e.defaultSettings = {
                {"strand.stitch_size", "custom"},
                {"strand.thickness_preset", "custom"},
                {"strand.diameter", "1.93"},
                {"strand.wale", "5.5"},
                {"strand.course", "3.95"},
                {"strand.tightness", "standard"},
                {"strand.compress_to_fit", "true"},
                {"strand.flatten", "true"},
                {"strand.flatten_inflate", "0.35"},
                {"strand.texture_profile", "2-swept"},
                {"strand.plies", "1"},
                {"strand.normal_strength", "0.6"},
                {"knit.coverage", "overlap"},
                {"knit.gap_patches", "false"},
                {"knit.hand.seed", "1"},
                {"knit.hand.unevenness", "20"},
                {"knit.hand.column_shift", "0.8"},
                {"knit.hand.column_rows", "5"},
                {"knit.hand.loop_rate", "0.6"},
                {"knit.hand.loop_height", "1"},
                {"knit.hand.wave_height", "0.5"},
                {"knit.hand.spot_rate", "0.8"},
                {"knit.hand.spot_change", "25"},
                {"knit.hand.spot_length", "16"},
                // No islands, as hR2 (the new objects' Worn would otherwise grow another coat).
                {"fibres.islands.rate", "0"},
                {"fibres.islands.size", "12"},
                {"fibres.islands.strength", "0"},
                {"fibres.islands.fuzzy", "1"},
                {"fibres.islands.matted", "1"},
                {"fibres.islands.worn", "1"},
                {"fibres.density", "350"},
                {"fibres.random_fraction", "0.03"},
                {"fibres.max_segments", "48"},
                {"fibres.segment_length", "0.45"},
                {"fibres.directed.length_min", "0.6"},
                {"fibres.directed.length_median", "1.5"},
                {"fibres.directed.length_max", "3.0"},
                {"fibres.directed.length_spread", "0.35"},
                {"fibres.directed.thickness", "0.036"},
                {"fibres.directed.thickness_sigma", "0.04"},
                {"fibres.directed.curvature", "1.0"},
                {"fibres.directed.max_angle_deg", "30"},
                {"fibres.directed.hug_height", "0.15"},
                {"fibres.random.length_min", "0.8"},
                {"fibres.random.length_median", "1.6"},
                {"fibres.random.length_max", "3.5"},
                {"fibres.random.length_spread", "0.5"},
                {"fibres.random.thickness", "0.034"},
                {"fibres.random.thickness_sigma", "0.08"},
                {"fibres.random.curvature", "0.4"},
                {"fibres.random.max_angle_deg", "25"},
                {"fibres.random.spread_deg", "30"},
                {"fibres.random.curl_radius", "0.5"},
                {"fibres.random.curl_length", "6"},
                {"fibres.tail_taper", "0.2"},
                {"fibres.tail_taper_min_diameters", "8"},
                {"fibres.tail_taper_tip", "0.3"},
                {"fibres.crimp_amplitude", "0"},
                {"fibres.crimp_wavelength", "1.5"},
                {"material.yarn_white.fibre_longitudinal_roughness", "0.6"},
                {"material.yarn_white.fibre_azimuthal_roughness", "0.6"},
                {"material.yarn_white.fibre_ior", "1.4"},
                {"material.yarn_white.translucency", "0.3"},
                {"material.fibre_directed.model", "fibre"},
                {"material.fibre_directed.fibre_longitudinal_roughness", "0.45"},
                {"material.fibre_directed.fibre_azimuthal_roughness", "0.6"},
                {"material.fibre_directed.fibre_ior", "1.4"},
                {"material.fibre_random.fibre_longitudinal_roughness", "0.6"},
                {"material.fibre_random.fibre_azimuthal_roughness", "0.6"},
                {"material.fibre_random.fibre_ior", "1.4"},
                {"material.fibre_hair_colour", "yarn"},
            };
            e.credit = {"Wool swatch on a needle", "Nicolas Menzel", "https://github.com/nicolasgfx/woolshader", "CC BY-NC 4.0",
                        "https://creativecommons.org/licenses/by-nc/4.0/"};
            v.push_back(e);
        }
        {
            // The blue swatch (owner, 2026-11-11: a second swatch preset whose yarn is made of fibres by default):
            // the wool swatch's GLB, knit and credit, with the blue scene's look: 400 fibres of 36 um per yarn across 3 plies
            // with an opaque core (yarn_body.*), the thick strands flattened by at most 15 % so the fabric closes
            // (strand.max_flattening, formerly KNIT_EXP_MAX_COMPRESSION), no filling, the fuzz's own angles and
            // taper, and the blue dye fitted to a reference photograph (0.0007, 0.145, 0.585 before the fit). Its body writes about 36 M segments.
            SceneEntry e = v.back();
            e.id = "wool-swatch-blue";
            e.label = "Blue wool swatch on a needle (12 cm)";
            e.projectName = "wool_swatch_blue";
            e.defaultYarnColourSource = "single";
            const std::vector<std::pair<std::string, std::string>> blue = {
                {"strand.texture_profile", "none"},  // the fibres are the plies (yarn_look.h: Fibre strands)
                {"strand.max_flattening", "0.15"},
                {"yarn_body.enabled", "true"},
                {"yarn_body.plies", "3"},
                {"yarn_body.fibres", "400"},
                {"yarn_body.fibre_diameter_um", "36"},
                {"yarn_body.shell_depth", "0.6"},
                {"yarn_body.irregularity", "1"},
                {"yarn_body.core", "opaque"},
                {"filling.enabled", "false"},
                {"fibres.directed.max_angle_deg", "15"},
                {"fibres.tail_taper", "0.3"},
                {"fibres.tail_taper_tip", "0.2"},
                {"material.fibre_directed.fibre_longitudinal_roughness", "0.6"},
                {"material.yarn_white.color", "0.0005, 0.1523, 0.4479"},
            };
            for (const std::pair<std::string, std::string>& kv : blue) setListedSetting(e.defaultSettings, kv.first, kv.second);
            v.push_back(e);
        }
        return v;
    }();
    return r;
}

namespace {
// The entry of knit_check --pieces (manifestScene).
std::unique_ptr<SceneEntry>& manifestEntry() {
    static std::unique_ptr<SceneEntry> e;
    return e;
}

const SceneEntry* sceneOfManifest(const std::string& manifestPath) {
    for (const SceneEntry& e : sceneRegistry()) {
        if (!e.knitsAsPieces()) continue;
        const std::string registered = findProjectFile(e.piecesManifest);
        std::error_code ec;
        if (!registered.empty() && std::filesystem::equivalent(registered, manifestPath, ec)) return &e;
    }
    return nullptr;
}
} // namespace

const SceneEntry* findScene(const std::string& id) {
    for (const SceneEntry& e : sceneRegistry())
        if (e.id == id) return &e;
    return manifestEntry() && manifestEntry()->id == id ? manifestEntry().get() : nullptr;
}

namespace {
std::vector<std::string>& exampleRestriction() {
    static std::vector<std::string> ids;
    return ids;
}
} // namespace

bool exampleModelPresent(const SceneEntry& e) {
    const std::vector<std::string>& only = exampleRestriction();
    if (!only.empty() && std::find(only.begin(), only.end(), e.id) == only.end()) return false;
    if (findProjectFile(e.path).empty()) return false;
    return e.piecesManifest.empty() || !findProjectFile(e.piecesManifest).empty();
}

void restrictExamplesForTests(const std::vector<std::string>& ids) { exampleRestriction() = ids; }

const SceneEntry* manifestScene(const std::string& manifestPath, const std::string& gaugeFile, double gaugeStitches, std::string* err) {
    const SceneEntry* registered = sceneOfManifest(manifestPath);
    if (registered && gaugeFile.empty()) return registered;
    if (gaugeFile.empty() || !(gaugeStitches > 0)) {
        if (err) *err = manifestPath + ": a manifest of no registered scene needs --gauge <piece file> <stitches>";
        return nullptr;
    }
    auto e = registered ? std::make_unique<SceneEntry>(*registered) : std::make_unique<SceneEntry>();
    const std::string stem = std::filesystem::path(manifestPath).parent_path().filename().string();
    e->id = "pieces";
    e->label = stem + " (pieces)";
    e->projectName = stem;
    if (!registered) {  // the spiral pieces' stitch and Dense yarn: no GLB, no colour sources of its own
        e->piecesW = 5.5;
        e->piecesH = 3.4;
        e->piecesD = 1.86;
        e->piecesDFallback = {1.80, 1.73, 1.67, 1.60};
        e->defaultFibrePreset = "standard";
    }
    e->piecesManifest = manifestPath;
    e->gaugePiece = gaugeFile;
    e->gaugeStitches = gaugeStitches;
    manifestEntry() = std::move(e);
    return manifestEntry().get();
}

const SceneEntry* sceneEntryFor(const SceneSettings& s) {
    if (s.scene != "import") return findScene(s.scene);
    uint64_t id = 0;
    if (!parseIdentityHex(s.importIdentity, id)) return nullptr;
    static std::mutex m;
    static std::map<std::string, std::unique_ptr<SceneEntry>> entries;
    std::lock_guard<std::mutex> lk(m);
    std::unique_ptr<SceneEntry>& e = entries[s.importIdentity + "|" + s.importSource];
    if (!e) {
        e = std::make_unique<SceneEntry>();
        const std::filesystem::path src(s.importSource);
        const std::string stem = src.stem().string().empty() ? std::string("mesh") : src.stem().string();
        e->id = "import";
        e->label = stem + " (imported)";
        e->projectName = stem;
        e->kind = SceneEntry::Kind::Mesh;
        std::error_code ec;
        e->path = std::filesystem::absolute(importCacheDir(importCacheRoot(), id), ec).string();
        e->glb.knitAllMaterials = true;
    }
    return e.get();
}

std::string creditLine(const ModelCredit& c) { return "\"" + c.title + "\" by " + c.author + ", " + c.licence; }

std::string creditText(const ModelCredit& c, Workflow w) {
    const char* changed = w == Workflow::Fur ? "a fur of curled fibres was grown on the model" : "the model was converted to knitted yarn and fibres";
    return "Made from \"" + c.title + "\" by " + c.author + "\n" + c.source + "\nLicensed under " + c.licence + ": " + c.licenceUrl +
           "\nChanged: " + changed + ".\n";
}

} // namespace knit
