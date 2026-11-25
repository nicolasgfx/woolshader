#pragma once
// The selectable scenes. A C++ table: typed and compile-checked (the tree
// has no JSON reader). Paths are repo-relative and resolved against the working
// directory, then the configure-time source tree (io/files.h sourceTreeDir).
#include <string>
#include <vector>

#include "io/glb_selection.h"
#include "model/vec.h"
#include "settings/scene_settings.h"

namespace knit {

// The model a preset is made from: what its exports must name.
struct ModelCredit {
    std::string title, author, source;
    std::string licence = "CC BY 4.0", licenceUrl = "https://creativecommons.org/licenses/by/4.0/";
};

struct SceneEntry {
    std::string id, label, projectName;
    // Glb: a registered GLB scene; Mesh: an imported OBJ / PLY (the prepared import directory,
    // import_store.h). Both go through the one workflow (UV-chart tiling); they differ only in where the
    // mesh comes from.
    enum class Kind { Glb, Mesh } kind = Kind::Glb;
    std::string path;  // GLB file, or the prepared import directory

    GlbSelection glb;  // what a GLB scene takes from its file; imports: every material knitted
    // No entry carries a workflow (owner, 2026-11-19): every model opens with neither chosen. The defaults
    // below are the Knit workflow's look; they change nothing under Fur.

    // The yarn colour source the scene opens with ("" = the settings default, "single" and "object"). A scene
    // with a base colour texture opens with "texture": the texture controls the yarn and the fibres. The teddy
    // has none (source/bears.fbx holds normal maps only) and takes "object": each bear its fabric's colour.
    std::string defaultYarnColourSource;
    // The wool fibre preset the scene opens with (presets/fibres id; "" = the first preset) and
    // whether its fibres are on (in the export) from the start.
    std::string defaultFibrePreset;
    bool defaultFibres = false;
    // The fibre density (candidates per mm; 0 = the settings default) and the fibre islands preset
    // (presets/fibre_islands id; "" = the new objects' one) the scene opens with in the viewer.
    double defaultFibreDensity = 0;
    std::string defaultFibreIslands;

    // ---- knitted as spiral pieces (the teddy; user, 2026-10-28): the same workflow and controls;
    // the Knit step knits these pieces as spirals (the spiral layout, sewn_pieces.h) instead of
    // tiling `path`'s UV charts, whose unsewn chart borders leave holes. The
    // fibres and the export are every scene's. Under Fur the scene is `path`'s mesh.
    std::string piecesManifest;  // the pieces manifest ("" = not knitted as pieces)
    std::string gaugePiece;      // the piece whose longest ring sets the scale
    // The calibration: gaugeStitches stitches of width piecesW around the gauge
    // piece set the size, which the bears keep at every stitch width (loadSpiralSource);
    // the scene opens with yarn piecesD at piecesW x piecesH and the fallback below it
    // (applySceneDefaults; its Dense preset, thicknessPresets).
    double gaugeStitches = 0, piecesW = 0, piecesH = 0, piecesD = 0;
    std::vector<double> piecesDFallback;
    // The pieces' materials (manifest material id -> the GLB material whose base colour they take
    // with yarn colour = object). The parts the selection keeps (glb.plainNodes) come from the GLB, placed
    // with the pieces.
    std::vector<std::pair<std::string, std::string>> pieceColourMaterials;
    bool knitsAsPieces() const { return !piecesManifest.empty(); }

    // Settings the scene opens with, as settings keys and values (one row each of settingFields(),
    // scene_settings.cpp): applied last wherever the scene is opened (applySceneSettings: the command line, the
    // viewer's start, Import > Examples), so all three give them. Never a key the profile owns (Draft / Final:
    // settle, repair, gates, fallback, attempts); the scene_settings_defaults_intended test checks both.
    std::vector<std::pair<std::string, std::string>> defaultSettings;

    // A preset: a model in scenes/ (exampleModelPresent says whether this checkout has it). Its exports name the model (CREDITS.txt). No credit: an
    // import, or a manifest of no registered scene.
    ModelCredit credit;
    bool isPreset() const { return !credit.title.empty(); }
};

const std::vector<SceneEntry>& sceneRegistry();
const SceneEntry* findScene(const std::string& id);
// Whether an example's model is in this checkout: its file and, for spiral pieces, its manifest resolve
// (findProjectFile). The published source ships only the wool swatch's model; the other examples' models keep
// their own licences and are not distributed, so the Import step offers an example only when this is true.
bool exampleModelPresent(const SceneEntry& e);
// Tests: only these example ids count as present, whatever is on disk (empty = no restriction), so a test's
// text is the same in every checkout.
void restrictExamplesForTests(const std::vector<std::string>& ids);
// knit_check --pieces: a pieces manifest as a scene of spiral pieces. A registered manifest is its
// scene; with another gauge, or for any other manifest (the gauge is then required), the entry "pieces"
// holds it for the process (findScene finds it), with the spiral pieces' stitch and Dense yarn. A
// registered manifest's entry keeps that scene's colours and its credit (so it stays a preset). Null
// with the reason.
const SceneEntry* manifestScene(const std::string& manifestPath, const std::string& gaugeFile, double gaugeStitches, std::string* err);
// The entry of the settings' scene: a registry entry, or for scene = "import" the entry of the
// prepared import named by import.identity (import_store.h; id "import", every material knitted). Null:
// unknown scene, or an import without a valid identity. Entries stay valid for the process.
const SceneEntry* sceneEntryFor(const SceneSettings& s);
// The credit in one line, and as the text of the file every export of the preset carries (its "Changed:" line
// says what the workflow made of the model).
std::string creditLine(const ModelCredit& c);
std::string creditText(const ModelCredit& c, Workflow w);

// The scene is the teddy's sewn pieces (SceneEntry::knitsAsPieces, scene/source/sewn_pieces.h); under Fur it is
// its mesh.
inline bool knitsSewnPieces(const SceneEntry& e, const SceneSettings& s) { return knitWorkflow(s) && e.knitsAsPieces(); }

} // namespace knit
