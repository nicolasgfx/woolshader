#pragma once
// The viewer's state and everything that touches it, with no ImGui and no GL, so all of it runs in
// the tests without a window. SceneUi is plain data; workflow.cpp holds the functions that change it
// (they take SceneUi&), presentation.cpp the ones that only read it (they take const SceneUi&).
// The panels draw what the readings say and call one action per click. src/viewer/ARCHITECTURE.md
// is the way into the folder.
#include <cstdint>
#include <filesystem>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <vector>

#include "model/coordinates.h"
#include "model/params.h"
#include "viewer/camera.h"
#include "export/export.h"
#include "jobs/fur_preview.h"
#include "jobs/import_session.h"
#include "jobs/previews.h"
#include "jobs/render_preview.h"
#include "jobs/scene_builder.h"
#include "settings/yarn_look.h"
#include "settings/scene_registry.h"
#include "viewer/marking_flow.h"
#include "viewer/viewer_prefs.h"

namespace knit {
struct FilePreset;
enum class PresetGroup;
}

namespace viewer {

// The workflow's steps in their order; Mark (surface markings for the stitch mesher, marking_flow.h) is optional.
enum class Stage { Import, Mark, Knit, Fibres, Export };
constexpr int kStages = 5;
enum class ViewMode { Mesh, Knit, Fibres };
// The Fibres step's viewport: the small reference sample (fibre_preview.h), the whole model, or under Fur the live
// disc preview on the model (jobs/fur_preview.h).
enum class FibreView { Reference, Model, Discs };
// The viewport's shading: Studio (per-pixel key light, hemisphere and ambient occlusion,
// mesh_lighting.h) or Flat (the fixed-function headlight; also the fallback without shaders).
enum class Lighting { Studio, Flat };

// UTF-8 in the Latin-1 font atlas (it has no ellipsis: "..." is written out).
const char* const kTimes = "\xC3\x97";
const char* const kDot = " \xC2\xB7 ";
// The fix, the stage word and the footer's reason while no workflow is chosen (SceneUi::workflowChosen).
const char* const kChooseWorkflow = "Choose Knit or Fur";

enum class UnitChoice { Unset, Mm, Cm, M, Inch, Other };  // the Import step's file unit (Set height is ImportFlow::sizeByHeight)

// The last attempt of a stage's job (what the stage card says besides its result).
struct Attempt {
    knit::JobOutcome outcome = knit::JobOutcome::None;
    std::string message;
};

// The one job the builder runs.
struct RunningJob {
    bool active = false, cancelRequested = false;
    bool meshOnly = false;         // the strand mesh rebuild: only the yarn tubes again (StrandMeshRequest)
    bool bakeMaps = false;         // that rebuild bakes the yarn maps first (startYarnMapBake): "Baking yarn maps"
    knit::BuiltScene::Kind kind = knit::BuiltScene::Kind::Knit;
    std::string stoppingNote;      // what a stopping job is finishing
    knit::JobProgress progress;  // one snapshot per frame (pollSceneUi)
};

// The other workflow's results, parked while this one is shown: switching the workflow never discards or
// recomputes. They come back only on the same object and prepared source (`source`); after a source change
// they are dropped, never shown or exported as current.
struct ParkedResults {
    knit::Workflow workflow = knit::Workflow::Knit;
    std::string source;  // parkedSource: the scene and the prepared import they were made on
    // An import's prepared source of that workflow (Knit: the rebuilt surface, Fur: the file's own mesh), selected
    // again on switching back when the import cache still holds it.
    std::string importIdentity, importSurface;
    std::shared_ptr<const knit::BuiltScene> knit, fibres;
    Attempt knitAttempt, fibreAttempt;
};

// Each completed and immutable, except `partial`: the running knit's scene so far, never a result.
struct Results {
    std::shared_ptr<const knit::BuiltScene> knit, partial, fibres, exported;
    ParkedResults parked;
    std::shared_ptr<const knit::BuiltScene> lastSeen;  // the builder's latest complete, already integrated
    Attempt knitAttempt, fibreAttempt, exportAttempt;
};

// The Import step's open file and the import it prepares.
struct ImportFlow {
    knit::ImportSession session;
    int generation = -1;
    bool fileOpen = false;          // stage 2: a file is read or waits for its choices
    bool revalidating = false;      // re-preparing an import named in the settings
    // Opened by the Knit step's Review Import: the prepared surface replaces the import named in the
    // settings and every Knit choice stays (importPrepared); any other open or a close ends it.
    bool keepKnitChoices = false;
    // Switched to Fur on a rebuilt import (prepareFileMesh): once the file is read, its own mesh is prepared with the
    // stored options; every setting stays.
    bool preparingFileMesh = false;
    // A scripted import asked to keep the file's surface under Knit too (--keep-file: the guided route knits from it).
    bool keepFileSurface = false;
    // The Knit workflow's rebuild of this file failed and the user chose to keep its own surface instead
    // (keepFileAfterFailedRebuild): the Mark step and the guided stitch mesh need no rebuilt surface. Once per opened file.
    bool keptAfterFailedRebuild = false;
    // The import kept the file's mesh: knitting needs a rebuilt surface first (pollSceneUi).
    bool lacksCharts = false;
    knit::ImportOptions options;    // the choices of stage 2; the surface follows the workflow (chooseWorkflow)
    UnitChoice unit = UnitChoice::Unset;  // the file unit; Set height never changes it
    double customUnitMm = 1, targetHeightMm = 300;
    bool unitTouched = false;
    bool sizeByHeight = false;      // "Set height": options.unitMm follows targetHeightMm along the up axis
    std::string focusedNode;        // ImportNode::name highlighted in the preview ("" none); never a setting
    bool revealFocus = false;       // the list scrolls the focused row into view (set by a pick in the preview)
    std::shared_ptr<const knit::BuiltScene> preview;  // file coordinates
    std::string note;               // e.g. the source file is missing
    std::vector<std::string> changes;
    // A preset's components: its fixed mesh selection, listed once per preset (presetNodes).
    std::string presetNodesOf;
    std::vector<knit::GlbNode> presetNodes;
};

// Where the export writes: the folder (default: <scene>_knit next to the imported file, or in
// <repo>/out/exports; a chosen one is kept until the scene changes) and the name of its files.
struct ExportTarget {
    std::string dir, name, format = knit::defaultExportFormat();
    bool dirChosen = false, nameEdited = false;
};

// What the viewport shows and how. None of these is a setting.
struct ViewChoices {
    ViewMode mode = ViewMode::Knit;
    FibreView fibres = FibreView::Reference;
    Lighting lighting = Lighting::Studio;
    bool wire = false;
    bool plyTexture = false;        // opt-in: high-frequency noise
    bool yarnMaps = true;           // the baked cavity maps on the tubes of a yarn of fibres rendered as maps
    bool unknittedParts = true;
    bool coverage = false;          // fur: the coverage overlay on the model
    bool fibrePopulations = false;  // wool fibres on the model yellow and blue (fibre_lines.h)
    knit::FurDiscColour discColour = knit::FurDiscColour::Plain;  // the Fur step's disc preview colouring
    float occlusion = 1.0f;         // Studio: the ambient occlusion's strength (0 none, 1 full)
    bool compare = false;           // Compare with source; kept while the comparison cannot be drawn
    // The shown camera's preset as of the last frame (presetOf): the window writes it after the camera
    // moved, the View button reads it. None: the camera was turned freely (Orbit).
    std::optional<CameraView> cameraPreset = CameraView::Front;
    bool studioAvailable = false;   // the viewer compiled the Studio shaders (studioNote: why not)
    std::string studioNote;
};

// Fur: the areas without fur, tinted on the model.
struct CoverageOverlay {
    std::shared_ptr<const knit::BuiltScene> scene;
    std::string key, error;
    double share = 1;
};

// What follows the settings without a job.
struct Previews {
    knit::FibrePreviewWorker fibres;
    knit::FlatReferences flat;
    std::shared_ptr<const knit::FibrePreview> shown;
    knit::FurPreviewWorker discs;                       // the Fur step's disc preview (FibreView::Discs)
    std::shared_ptr<const knit::FurPreview> discShown;  // its last complete preview
    knit::RenderPreviewWorker render;  // the live render preview's renderer process (Render preview window)
};

// The Render preview window: what it shows and what it asked for. The renderer (renderer.cpp) uploads
// the newest frame into `texture`; the interface draws it; updateRenderPreview builds the request.
struct RenderPreviewView {
    bool paused = false;
    int scalePercent = 50;            // the render's size against the window's image area
    float imageW = 0, imageH = 0;     // the window's image area [px], from the last frame
    unsigned texture = 0;             // GL texture id (0: none yet)
    int textureW = 0, textureH = 0;
    uint64_t uploadedFrame = 0;       // RenderPreviewState::frameNumber in the texture
    knit::RenderPreviewState state;   // copied every frame
    std::string hint;                 // why the preview shows an older state ("" none)
    std::string busy;                 // what the preview is waiting for or doing, over the image ("" none)
    // The package key waits until the edits pause (a slider drag would export on every frame).
    std::string settingsText;
    double settingsChangedAt = 0;
    uint64_t packageKey = 0;
    std::shared_ptr<const knit::SceneSettings> packageSettings;  // made when settingsText changes
};

// The fur mask image as it is on disk (checkMaskImage): hashed again only when its path, size or write
// time changes.
struct MaskOnDisk {
    std::string path, hash;
    std::filesystem::file_time_type writeTime{};
    std::uintmax_t size = 0;
};

// The Ply style sketches of the last plies and fill drawn: plyCrossSection takes 0.3 to 0.6 s for both
// styles (measured 2026-11-06), so they are kept, not made per frame; one plies / fill at a time.
struct PlySketch {
    int plies = 0;
    double fill = -1;
    std::vector<knit::Vec2> swept, circles;
};

// What the Fibres step's Yarn fibres section scrolls to and opens on its next frame, once: the section's
// heading, or one of its detail nodes (a link from the Knit step or the Export step's fix).
enum class Reveal { None, Section, Structure, Appearance };

// What the panels keep that is no setting.
struct PanelState {
    bool customTessellation = false;  // Yarn tessellation shows its custom sliders
    bool reportOpen = false, logOpen = false, previewOpen = false, aboutOpen = false;
    // The Fibres step's population tab to select next frame (-1 none; the tab only chooses what is
    // edited: never a setting, a key or a visibility).
    int fibreTabRequest = -1;
    int aboutTabRequest = -1;  // the About window's tab to select next frame (-1 none; --do about:<n>)
    knit::YarnLookMemory yarnMemory;  // the yarn look's mode not shown (restored on switching back)
    Reveal revealYarnFibres = Reveal::None;  // the Fibres step's Yarn fibres section scrolls there next frame
    MaskOnDisk maskOnDisk;            // the Fibres step's mask image check
    PlySketch plySketch;              // the Knit step's Ply style sketches
};

// Sizes of the last frame: drawSceneUi writes them, viewportRect and the footers read them.
struct PanelLayout {
    float panelWidth = 0, statusHeight = 0, stripHeight = 0, toolbarHeight = 0;
    float viewStatusHeight = 0;  // the line along the viewport's bottom edge: the object, the view's state
    float footerHeight[kStages] = {};  // by Stage
};

// What the interface asks of the window this frame. main.cpp answers after drawSceneUi
// (answerPlatformRequest) and gives a dialog's answer to the named action; a test answers it itself.
enum class PlatformRequest { None, OpenMesh, OpenMaskImage, ChooseExportFolder, ShowExportFolder, SaveMarkings, OpenMarkings, SaveStitchMesh };
struct PlatformServices {
    bool dialogs = false;      // native file and folder dialogs (else the controls say how instead)
    bool fileManager = false;  // "Open folder" after an export
    PlatformRequest request = PlatformRequest::None;
};

struct SceneUi {
    knit::SceneSettings s;  // being edited
    // The draft's workflow: unchosen, or s.workflow (knit | fur). A newly opened model is unchosen: the Import step
    // asks, the later steps say "Choose Knit or Fur". A settings file always carries a choice.
    bool workflowChosen = false;
    std::unique_ptr<knit::SceneBuilder> builder;
    std::string cacheRoot;  // the persistent scene cache ("" none)
    std::string loadKnit;   // --load-knit (scripts): every Knit job takes its spiral pieces from this knit cache directory
    Stage stage = Stage::Knit;
    Results results;
    RunningJob job;
    ImportFlow import;
    MarkingFlow mark;
    ExportTarget exportTo;
    ViewChoices view;
    CoverageOverlay coverage;
    Previews previews;
    RenderPreviewView renderPreview;
    PanelState panels;
    PanelLayout layout;
    CameraRequests camera;  // applied by the window after the frame
    PlatformServices platform;
    ViewerPrefs* prefs = nullptr;
};

// What a scripted import decides at confirmation (--units; the unit itself: presetImportUnits).
struct ImportChoices {
    std::set<std::string> plainNodes;  // the components left unticked
    bool keepFile = false;             // the file's surface kept as it is in either workflow (no rebuild; the guided route knits from it)
};

// Export detail: the renderer's GPU memory is mostly the yarn tubes' triangles and the fibre
// segments. Each level sets the yarn tessellation level (presets/strand_mesh) and the segment budget.
struct ExportDetail {
    const char* label;
    const char* meshLevel;  // a yarn tessellation level id
    long long budget;       // fibres.segment_budget
    const char* note;
};
extern const ExportDetail kExportDetails[3];  // Full, Reduced, Light

// --- What changes SceneUi (workflow.cpp)

// Integrates finished jobs and imports (call once per frame, before the interface).
void pollSceneUi(SceneUi& ui);
// The Render preview's request from the current results, settings and camera (every frame; the worker
// acts on what changed), and its state copied into ui.renderPreview.
void updateRenderPreview(SceneUi& ui, const Camera& camera, bool cameraMoving, double now);

// The window's start when nothing names an object: the Import step and no object (an empty scene
// id). The settings stay as they are, and the first object opened starts from them.
void startWithoutObject(SceneUi& ui);
void openExample(SceneUi& ui, const std::string& sceneId);
void openMeshFile(SceneUi& ui, const std::string& path);
void changeImport(SceneUi& ui);            // "Change the import ...": the choices stay, units are suggested again
void prepareKnittingCharts(SceneUi& ui);   // "Rebuild the surface ...": units and axis stay
// The Knit step's "Review Import": the same file opened to prepare its surface for the current stitch
// width (as prepareKnittingCharts); once prepared, the Knit step's settings are kept as they are.
void reviewImport(SceneUi& ui);
// Fur on an import whose surface was rebuilt: the same file prepared as its own mesh (KeepFile), every setting kept.
void prepareFileMesh(SceneUi& ui);
// A scripted import, in two steps on purpose: the unit is set BEFORE the first poll, which then
// skips the suggestion (integrateImport); the rest is confirmed once the file is analysed.
void presetImportUnits(SceneUi& ui, double unitMm, knit::UpAxis up);
void confirmImport(SceneUi& ui, const ImportChoices& choices);
void closeImport(SceneUi& ui);
// The Import step's choices that carry a rule. The plain ones (engine, repairs) are bound to
// ui.import.options directly, as settings are. A size that cannot be used sets options.unitMm to 0,
// which disables the import (importFooterView says why).
void chooseImportUnit(SceneUi& ui, UnitChoice unit);      // "File units": mm, cm, m, in: that unit; Custom: starts from the current size
void setCustomImportUnit(SceneUi& ui, double mmPerUnit);  // "mm per file unit": 0 or less is no size
void chooseSizeByHeight(SceneUi& ui, bool byHeight);      // "Set height" starts from the current height; "File units" brings the unit back
void setImportHeight(SceneUi& ui, double heightMm);       // "Target height": the unit follows the file's height along the up axis
void chooseImportUp(SceneUi& ui, knit::UpAxis up);        // and a size set by height follows the new axis
void toggleNode(SceneUi& ui, const std::string& name);    // the row's checkbox: knitted or kept as a mesh
void selectAllNodes(SceneUi& ui);                         // "Knit all"
// A preset's components and its fixed selection (the Import step's rows), read once per preset; empty when
// its file is missing.
const std::vector<knit::GlbNode>& presetNodes(SceneUi& ui, const knit::SceneEntry& e);
void selectNoNodes(SceneUi& ui);                          // "Keep all as mesh"
void focusNode(SceneUi& ui, const std::string& name);     // the row's name: highlighted in the preview, the tick unchanged
void pickNode(SceneUi& ui, const Camera& camera, const ViewPixel& pixel);  // focuses the component under the pixel and shows its row
void startImport(SceneUi& ui);                            // the footer's action: prepares importRequest(ui.import)
// After a failed rebuild in the Knit workflow (an open tube has no inside to wrap), the user may keep the file's own
// surface instead: the Mark step and the guided stitch mesh need no rebuilt surface. Never automatic: the failure stays
// on screen until the user chooses this, another setting, or Cancel.
bool canKeepFileAfterFailedRebuild(const SceneUi& ui);
void keepFileAfterFailedRebuild(SceneUi& ui);
// After loading a settings file: an import named in it is re-validated (its source re-hashed).
void adoptSettings(SceneUi& ui);

void startKnit(SceneUi& ui);
void startStrandMesh(SceneUi& ui);         // rebuilds the shown knit's yarn tubes (no knitting)
void startYarnMapBake(SceneUi& ui);        // bakes the yarn maps the cache lacks and shows them on the tubes (no knitting)
// Opens the Fibres step at its Yarn fibres section (a link in the Knit step, the Export step's fix).
void revealYarnFibres(SceneUi& ui, Reveal what);
void startFibres(SceneUi& ui);
void startExport(SceneUi& ui);
void cancelJob(SceneUi& ui);
void knitAtFinalQuality(SceneUi& ui);  // the Final profile, then the knit (the Knit footer; an export that needs a Final knit)
void knitAtDraftQuality(SceneUi& ui);  // the Draft profile, then the knit (the Knit footer's "Knit draft")
void runExportFix(SceneUi& ui, const knit::ExportEligibility& el);  // what a blocked export's button does
// A completed fibre job becomes the shown fibres; wool fibres that were generated are included in
// the export (Export > Contents > Include fibres; unticking it afterwards is kept until the next result).
void takeFibreResult(SceneUi& ui, std::shared_ptr<const knit::BuiltScene> done);
// Knit or Fur for the object (the Import step's choice, or a switch later): the settings and last results of
// both workflows are kept (the other's results parked, ParkedResults); nothing is computed but a missing prepared
// source. An open file's import follows it: Fur grows on the file's own mesh, Knit rebuilds the surface. A prepared
// import switches to the workflow's own source: the parked one when the cache holds it, else Fur prepares the file's
// mesh (Knit says to rebuild the surface, as before). The component choices stay (plain in both). Refused while a
// job runs, so no result lands in the other workflow.
void chooseWorkflow(SceneUi& ui, knit::Workflow w);
bool applyExportDetail(SceneUi& ui, const ExportDetail& d);
void chooseExportFolder(SceneUi& ui, const std::string& dir);  // the folder dialog's answer: kept until the scene changes

void loadMaskImage(SceneUi& ui, const std::string& path);  // the image dialog's answer: the mask and its hash
void reloadMaskImage(SceneUi& ui);                          // "Reload": accepts the image as it is now
// Checks the loaded mask image against the file on disk now: true when it differs (edited, replaced,
// deleted, unreadable). It reads the file's size and write time on every call and hashes the file only
// when they changed (ui.panels.maskOnDisk). The Fibres step calls it where it draws the warning, so the
// disk is read only while the warning can show.
bool checkMaskImage(SceneUi& ui);

// --- What reads it (presentation.cpp)

enum class StageStatus { NotStarted, Running, Cancelling, Cancelled, Failed, NeedsUpdate, ViewOnly, Blocked, Ready, Optional, Diagnostic };

struct StageState {
    StageStatus status = StageStatus::NotStarted;
    std::string word, detail;
};

// What is current.
const knit::SceneEntry* entryOf(const SceneUi& ui);
const knit::KnitResult* knitResult(const SceneUi& ui);
const knit::FibreResult* fibreResult(const SceneUi& ui);
bool furShown(const SceneUi& ui);  // the Fur workflow is chosen
bool knitCurrent(const SceneUi& ui);
bool fibresInputReady(const SceneUi& ui);
bool runningStage(const SceneUi& ui, knit::BuiltScene::Kind k);
bool knitFooterShown(const SceneUi& ui);
bool sameTube(const knit::SceneSettings& a, const knit::SceneSettings& b);
bool shownTubesDiffer(const SceneUi& ui);  // the shown, current knit's tubes are not at the (valid) tessellation in the settings
std::string waitReason(const SceneUi& ui);

// The strip's steps: Import -> Mark -> Knit -> Fibres -> Export, or Import -> Fur -> Export under Fur (Mark and Knit
// left out, the Fibres step's place holds the Fur step); each with its label ("1  Knit", "1  Fur").
bool stageInStrip(const SceneUi& ui, Stage s);
const char* stageLabel(const SceneUi& ui, Stage s);

// The stages' words, as the strip and the footers show them.
StageState importState(const SceneUi& ui);
StageState markState(const SceneUi& ui);
StageState knitState(const SceneUi& ui);
StageState fibreOutputState(const SceneUi& ui);
StageState fibreState(const SceneUi& ui);
StageState exportState(const SceneUi& ui);
knit::ExportEligibility eligibility(const SceneUi& ui);

// --- The Knit step's readings (presentation.cpp)

// The stitch width an imported surface was rebuilt for, against the current one ("Prepared surface").
struct PreparedSurface {
    bool applies = false;   // an import whose surface was rebuilt (wrapped) for a stitch width
    double preparedMm = 0;  // import.wrap_stitch_mm
    bool tooWide = false;   // the stitch is more than 1.3 x that width: its yarns may tangle on the surface
};
PreparedSurface preparedSurfaceCheck(const knit::SceneSettings& s);

// What one footer action would knit with: the settings with that quality's profile applied to a copy,
// so its notes and words are about the job it submits, not about the profile the settings hold now.
enum class KnitQuality { Draft, Final };
struct KnitActionView {
    std::vector<std::string> notes;  // what cannot run on this machine at that quality (knit::computeNotes)
    std::string tooltip;             // what the action does; Final adds its steps and diameter fallback
};
KnitActionView knitActionView(const SceneUi& ui, KnitQuality quality, const knit::ComputeCaps& caps);

// The generated fibres were rooted in other yarn tubes than the tessellation in the settings (a level
// changed, applied or not): they must be generated again. False without fibres.
bool fibresNeedNewTubes(const SceneUi& ui);

// The Knit step's labels and the one-line summaries of its collapsed sections.
const char* patternLabel(knit::StitchPattern p);                                        // "Rib 1 x 1"
std::string gaugeChoiceLabel(const knit::SceneSettings& s, const knit::SceneEntry& e);  // "Regular", "Dense", "Custom"
std::string yarnGeometrySummary(const knit::SceneSettings& s);  // "0.80 D · 0.10 mm clearance · contact flattening"
std::string characterSummary(const knit::SceneSettings& s);     // "20 % uneven · 0.30 raised loops / cm2 · ..."
std::string stitchLayoutSummary(const knit::SceneSettings& s, const knit::SceneEntry& e);  // "Overlapped seams · 0 deg", "Sewn pieces"

// Fit check 2, the flat swatch (jobs/previews.h FlatReferences, the references for the current settings, null
// while they are computed): a pattern's verdict at the current stitch as a state and one line. `outcome` is
// "fails (clearance)" or "is unresolved (...)" when it does not pass, else "". Information, never a lock: the bare
// template fails at every gauge preset (thickness_presets_flat_matrix) and knits.
struct FlatSwatchView {
    StageStatus status = StageStatus::Running;
    std::string text, outcome;
};
FlatSwatchView flatSwatchOf(const SceneUi& ui, const std::vector<knit::PatternReference>* refs, knit::StitchPattern p);

// The yarn made of fibres (yarn_body.*): what it costs the renderer, beside the fuzz, from the shown knit's
// yarn length and the generated fibres. One owner of the arithmetic for the Knit step's summary line
// (yarnFibresSummary) and the Fibres step's cost line (yarnFibresCostLine). Arithmetic only per frame.
struct YarnFibresCost {
    bool knitted = false;           // a shown knit with a yarn length: without it nothing is known
    bool maps = false;              // rendered as maps: texels and tube triangles instead of yarn segments
    double yarnMm = 0;              // the knit's yarn length
    double yarnSegments = 0;        // the yarn's fibre segments (strands)
    double fuzzSegments = 0;        // the generated fuzz's exported segments, 0 when none count
    double texels = 0, triangles = 0;  // maps: the atlases' texels, the tubes' triangles
    double texturesGib = 0;         // maps: the atlases' share of the memory
    double gib = 0;                 // the renderer's GPU memory for all of it
    bool yarnDominates = false;     // the yarn's fibres (or maps) take most of the memory
    std::string fuzzWord;           // "4.1 M segments" | "not generated yet" | "not exported"
};
YarnFibresCost yarnFibresCost(const SceneUi& ui);
std::string yarnFibresSummary(const SceneUi& ui);             // Knit step: "Smooth · 3 plies · 400 fibres · ~34.3 M yarn segments + 4.1 M fuzz · ~7.5 GiB while rendering"
std::string yarnFibresCostLine(const YarnFibresCost& c);      // Fibres step: the full line
std::string yarnFibresMemoryAdvice(const YarnFibresCost& c);  // above the cap: what takes most of it; "" below
// The Fibres step's node that holds a refused yarn body control (the Export fix and the Knit step's link open it).
Reveal revealOfBodyField(knit::YarnBodyField f);
std::string meshDetailSummary(const knit::SceneSettings& s, bool customShown);  // the level's label, or "Custom (5 sides, 0.40 mm)"

// A running job's progress as the interface shows it (one rule for the footers and the status bar):
// the counters are called parts only on the knitting steps ("knitting", "knitting <part>", "knitted
// <part>"); another step with a total shows a bar without unit words; without a total the bar is
// indeterminate (fraction < 0). Never a time estimate.
struct JobProgressView {
    std::string title;      // "Knitting ...", "Meshing ...", "Stopping ..."
    std::string count;      // "7 of 18 parts" or ""
    std::string qualifier;  // what the job is doing, or that a partial build is shown
    std::string status;     // the status bar's line: "Knitting piece_x · Elapsed 6 min 40 s"
    float fraction = -1;    // [0, 1], or < 0: indeterminate
    bool stopping = false;  // a cancel was requested
};
JobProgressView progressView(const RunningJob& job, bool partialShown);

// Numbers as text.
std::string fmtCount(long long n);
std::string fmtSeconds(double s);
// The stray share (fibres.random_fraction, 0..1) as the Fibres step shows it (percent) and back:
// the conversion happens only at the interface (the setting keeps the fraction).
double sharePercent(double fraction);
double shareFraction(double percent);
// The candidate workload line: density [candidates / mm of yarn] x yarn length [mm], e.g.
// "Approx. 16.2M candidates · 129.5 m yarn" (candidates, not the fibres that remain).
std::string candidateEstimate(double densityPerMm, double yarnLengthMm);

// What the export does with the wool fibres.
enum class FibreExport { Excluded, Omitted, Included, Generating, Update, Generate };
FibreExport fibreExport(const SceneUi& ui);
bool fibresGeneratedExcluded(const SceneUi& ui);
std::string fibreFooterExportNote(const SceneUi& ui);  // "Export: included, generate the fibres first" (the fuzz on a yarn of fibres)

// The Export step's texts.
bool exportDetailSelected(const SceneUi& ui, const ExportDetail& d);
const knit::FilePreset* meshLevelOf(const ExportDetail& d);
// The forecast of one level: the tube triangles, the fibre segments as the export writes them and the
// renderer's GPU memory (knit::rendererMemoryGiB); the warning above knit::kRendererMemoryCapGiB ("": none).
std::string exportDetailForecast(const SceneUi& ui, const ExportDetail& d);
std::string exportDetailMemoryWarning(const SceneUi& ui, const ExportDetail& d);
std::string exportDetailNote(const SceneUi& ui);
std::string fibreExportWord(const SceneUi& ui);
std::string exportReadyLine(const SceneUi& ui);
std::string exportOutputText(const SceneUi& ui);
std::string exportDoneText(const SceneUi& ui);
std::string exportExistingNote(const SceneUi& ui);
std::string exportTargetProblem(const SceneUi& ui);  // "" usable, else why not
std::string exportFolderStart(const SceneUi& ui);    // where the folder dialog opens: the destination, or its parent while it does not exist
std::string formatSummary(const SceneUi& ui, const knit::ExportFormat& f);
std::string formatShort(const knit::ExportFormat* f);
std::string defaultExportName(const SceneUi& ui);
bool exportNamesFiles(const SceneUi& ui);            // the format's files take the name (the package's names are fixed)
std::string exportHoldsText(const SceneUi& ui);      // the files that take the name, "" if none do
bool exportsMesh(const SceneUi& ui);
const char* blockedHeading(const std::string& fix);

// The Import step's open file.
bool importFooterShown(const SceneUi& ui);
// The one rule for a component kept as a mesh: ticked off, in both workflows (no yarn and no fur on it; it still
// blocks the fur).
bool nodeKeptAsMesh(const ImportFlow& im, const std::string& name);
knit::ImportOptions importRequest(const ImportFlow& im);  // what Import prepares: the choices
double tickedNodeArea(const SceneUi& ui);  // file units^2: what the stitch estimate counts; 0 with nothing ticked
int knittedNodeCount(const SceneUi& ui);
double importStitchCount(const SceneUi& ui);  // the knitted area over the stitch area at the chosen size (an upper bound when rebuilt)
// The object's size after the unit and up axis, as the import places it [mm]. valid: the unit is usable.
struct ImportSize {
    double width = 0, depth = 0, height = 0;
    bool valid = false;
};
ImportSize importSize(const SceneUi& ui);
std::string sizeProblem(const SceneUi& ui);       // "" usable, else what to correct
std::string suggestedUnitName(double unitMm);     // "Millimetres", "Centimetres", "Metres", "Inches"
std::string importSourceLine(const SceneUi& ui);  // "OBJ · 53.5 MB · 2 components"
std::string importColourSummary(const SceneUi& ui);  // "Base colour · 8 textures"
// The footer: what the import will make, how, the action's label and why it is disabled ("" enabled).
struct ImportFooterView {
    std::string summary, status, action, whyDisabled;
};
ImportFooterView importFooterView(const SceneUi& ui);

} // namespace viewer
