// What changes SceneUi (workflow.h): opening an object, submitting a job to the builder, switching
// the scene, and pollSceneUi, which integrates finished jobs and imports once per frame. What only
// reads it is in presentation.cpp. No ImGui and no GL.
#include "model/coordinates.h"
#include "io/child_process.h"
#include "io/files.h"
#include "model/env.h"
#include "settings/scene_registry.h"
#include "settings/settings_profiles.h"
#include "viewer/view_content.h"
#include "viewer/workflow.h"

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <string>
#include <utility>

#include "settings/file_presets.h"
#include "export/export.h"
#include "settings/thickness_presets.h"
#include "settings/settings_mapping.h"

namespace fs = std::filesystem;

namespace viewer {

namespace {

// A job of this kind starts (meshOnly: the strand mesh rebuild): its stage's attempt is reset and the
// running state set; returns the settings snapshot its request runs on.
std::shared_ptr<const knit::SceneSettings> jobStarting(SceneUi& ui, knit::BuiltScene::Kind kind, bool meshOnly) {
    ui.job.kind = kind;
    ui.job.meshOnly = meshOnly;
    ui.job.bakeMaps = false;
    ui.job.active = true;
    ui.job.cancelRequested = false;
    (kind == knit::BuiltScene::Kind::Knit ? ui.results.knitAttempt : kind == knit::BuiltScene::Kind::Fibres ? ui.results.fibreAttempt : ui.results.exportAttempt) = Attempt{};
    return std::make_shared<const knit::SceneSettings>(ui.s);
}

std::string defaultExportParent(const SceneUi& ui) {
    if (ui.s.scene == "import" && !ui.s.importSource.empty()) return fs::path(ui.s.importSource).parent_path().make_preferred().string();
    // The examples: <repo>/out/exports (git-ignored), created when it is the destination.
    std::error_code ec;
    const std::string tree = knit::sourceTreeDir();
    const fs::path d = tree.empty() ? fs::absolute("out/exports", ec).make_preferred() : (fs::path(tree) / "out" / "exports").make_preferred();
    fs::create_directories(d, ec);
    return d.string();
}

void resetDestination(SceneUi& ui) {
    if (!ui.exportTo.dirChosen) ui.exportTo.dir = (fs::path(defaultExportParent(ui)) / defaultExportName(ui)).make_preferred().string();
    if (!ui.exportTo.nameEdited) ui.exportTo.name = defaultExportName(ui);
}

void resetResults(SceneUi& ui) {
    ui.results.knit.reset();
    ui.results.partial.reset();
    ui.results.fibres.reset();
    ui.results.parked = ParkedResults{};
    ui.results.exported.reset();
    ui.coverage.scene.reset();
    ui.results.knitAttempt = ui.results.fibreAttempt = ui.results.exportAttempt = Attempt{};
    ui.camera.frame = true;
}

// The viewer's settings for a scene: the scene defaults on this machine's compute defaults, the
// user's explicit compute choices, and the current knitting / fibre / material choices kept.
void switchScene(SceneUi& ui, const knit::SceneSettings& next, const knit::SceneEntry& e) {
    knit::SceneSettings s = next;
    const knit::SceneEntry* prev = knit::sceneEntryFor(ui.s);
    // The previous scene's own settings (SceneEntry::defaultSettings) stay with it: a key that still holds the value
    // it listed takes this scene's start value (viewerDefaults); a key the user changed is kept, like every choice.
    if (prev && !prev->defaultSettings.empty()) {
        const knit::SceneSettings start = knit::viewerDefaults(e);
        for (const auto& [key, value] : prev->defaultSettings) {
            knit::SceneSettings listed = ui.s;
            std::string err;
            if (knit::setSetting(listed, key, value, &err) && knit::sameSetting(ui.s, listed, key)) knit::copySetting(s, start, key);
        }
    }
    // The thickness preset is kept by id (a scene without it, such as a block scene, takes its
    // first); a Custom choice is kept.
    std::string presetId = "thick-dense";
    if (prev) {
        const knit::ThicknessPreset* p = knit::matchingPreset(s, *prev);
        presetId = p ? p->id : std::string();
    }
    const std::vector<knit::ThicknessPreset> presets = knit::thicknessPresets(e);
    if (!presetId.empty()) {
        const knit::ThicknessPreset* use = &presets.front();
        for (const knit::ThicknessPreset& q : presets)
            if (q.id == presetId) use = &q;
        knit::applyThicknessPreset(s, *use);
    }
    s.angleDeg = 0;
    // The seams (knit.coverage): spiral pieces have no charts to continue (only "component"); a
    // uv-chart scene keeps the user's choice unless it comes from them, then it takes its default
    // (overlapped seams).
    if (e.knitsAsPieces() || !prev || prev->knitsAsPieces()) {
        knit::SceneSettings def;
        knit::applySceneDefaults(def, e);
        s.coverage = def.coverage;
    }
    // Spiral pieces (the teddy) knit stockinette in one colour per bear or one colour.
    if (e.knitsAsPieces()) {
        s.pattern = "stockinette";
        if (s.yarnColourSource == "texture") s.yarnColourSource = "object";
    }
    if (knit::isPreviewProfile(ui.s)) knit::applyPreviewProfile(s);
    else knit::applyFullProfile(s, e);
    if (ui.s.yarnColourSource.empty()) s.yarnColourSource = "single";
    // The scene's own wool preset and fibres.
    if (e.defaultFibres) s.fibres = true;
    if (const knit::FilePreset* p = knit::presetWithId(knit::PresetGroup::Fibres, e.defaultFibrePreset)) knit::applyPreset(knit::PresetGroup::Fibres, s, *p);
    knit::applySceneFibreDefaults(s, e);
    knit::applySceneSettings(s, e);  // the scene's own settings, last (as viewerDefaults and the command line)
    knit::applyStitchSize(s);  // yarn mode: w and h follow the preset set above
    ui.s = s;
    resetResults(ui);
    ui.exportTo.dirChosen = false;  // the export writes straight into a chosen folder: another scene's would overwrite it
    resetDestination(ui);
}

// A knitted object is rebuilt; fur grows on the file's own mesh (no knitting charts), so the file's mesh is kept.
void applySurface(SceneUi& ui) {
    const bool fur = ui.workflowChosen && knit::furWorkflow(ui.s);
    ui.import.options.surface = fur || ui.import.keepFileSurface ? knit::ImportOptions::Surface::KeepFile : knit::ImportOptions::Surface::Wrapped;
}

// The object and prepared source a workflow's results were made on: parked results come back only on the same.
std::string parkedSource(const knit::SceneSettings& s) { return s.scene + "|" + s.importIdentity + "|" + s.importSurface; }

// The size of one file unit as the "File units" choice says [mm]; 0 for a Custom value that is no size.
double fileUnitMm(const ImportFlow& im) {
    switch (im.unit) {
    case UnitChoice::Unset: return im.options.unitMm;
    case UnitChoice::Mm: return 1;
    case UnitChoice::Cm: return 10;
    case UnitChoice::M: return 1000;
    case UnitChoice::Inch: return 25.4;
    case UnitChoice::Other: return im.customUnitMm > 0 && std::isfinite(im.customUnitMm) ? im.customUnitMm : 0;
    }
    return 0;
}

// While the size is set by height, the unit follows the file's extent along the up axis. A height of 0 or
// less, or a file that is flat along that axis, gives no size (unit 0: Import disabled, sizeProblem says why).
void unitFromHeight(SceneUi& ui) {
    ImportFlow& im = ui.import;
    const std::shared_ptr<const knit::ImportAnalysis> an = im.session.snapshot().analysis;
    if (!im.sizeByHeight || !an) return;
    const knit::Vec3 extent = an->hi - an->lo;
    const double height = im.options.up == knit::UpAxis::Y ? extent.y : extent.z;
    const bool usable = height > 0 && im.targetHeightMm > 0 && std::isfinite(im.targetHeightMm);
    im.options.unitMm = usable ? im.targetHeightMm / height : 0;
}

} // namespace

void startWithoutObject(SceneUi& ui) {
    ui.s.scene.clear();
    ui.stage = Stage::Import;
}

void openExample(SceneUi& ui, const std::string& sceneId) {
    const knit::SceneEntry* e = knit::findScene(sceneId);
    if (!e) return;
    knit::SceneSettings s = ui.s;
    s.scene = e->id;
    s.yarnColourSource = e->defaultYarnColourSource.empty() ? "single" : e->defaultYarnColourSource;
    s.pileColourSource = "object";  // each part its model's stored colour
    switchScene(ui, s, *e);
    ui.import.preview.reset();
    ui.import.note.clear();
    ui.import.changes.clear();
    // No model carries a workflow (owner, 2026-11-19): the Import step asks for Knit or Fur.
    ui.workflowChosen = false;
    ui.stage = Stage::Import;
    ui.view.mode = ViewMode::Mesh;
    // With the Fur switch off (knit::furEnabled) there is nothing to ask: the example opens on its knit.
    if (!knit::furEnabled()) chooseWorkflow(ui, knit::Workflow::Knit);
}

void chooseWorkflow(SceneUi& ui, knit::Workflow w) {
    if (ui.job.active) return;  // a running job's result belongs to the workflow it started in
    const bool firstChoice = !ui.workflowChosen;
    const bool switching = ui.workflowChosen && knit::workflowOf(ui.s) != w;
    if (switching) {
        // This workflow's results are parked; the other's come back when they were made on the same source.
        ParkedResults now{knit::workflowOf(ui.s), parkedSource(ui.s), ui.s.importIdentity, ui.s.importSurface, ui.results.knit,
                          ui.results.fibres, ui.results.knitAttempt, ui.results.fibreAttempt};
        const ParkedResults back = ui.results.parked;
        // An import goes back to that workflow's own prepared source while the cache holds it.
        if (back.workflow == w && ui.s.scene == "import" && !back.importIdentity.empty()) {
            knit::SceneSettings other = ui.s;
            other.importIdentity = back.importIdentity;
            other.importSurface = back.importSurface;
            if (knit::sceneEntryFor(other)) {
                ui.s.importIdentity = back.importIdentity;
                ui.s.importSurface = back.importSurface;
            }
        }
        const bool restore = back.workflow == w && back.source == parkedSource(ui.s);
        ui.results.knit = restore ? back.knit : nullptr;
        ui.results.fibres = restore ? back.fibres : nullptr;
        ui.results.knitAttempt = restore ? back.knitAttempt : Attempt{};
        ui.results.fibreAttempt = restore ? back.fibreAttempt : Attempt{};
        ui.results.parked = std::move(now);
        ui.results.partial.reset();
        ui.results.exported.reset();
        ui.coverage.scene.reset();
        ui.coverage.key.clear();
    }
    ui.s.workflow = knit::toString(w);
    ui.workflowChosen = true;
    if (ui.import.fileOpen) applySurface(ui);  // the open file's import follows the workflow
    resetDestination(ui);                      // the default folder says it ("_knit" / "_fur")
    if (ui.import.fileOpen || !entryOf(ui)) return;  // the Import step goes on: Import prepares the file
    // Fur on a rebuilt import: the file's own mesh is prepared (in the background; every setting stays).
    if (w == knit::Workflow::Fur && ui.s.scene == "import" && ui.s.importSurface == "wrapped") prepareFileMesh(ui);
    if (w == knit::Workflow::Fur) {
        ui.stage = Stage::Fibres;  // the Fur step: nothing is knitted
        ui.view.mode = ViewMode::Fibres;
        ui.view.fibres = FibreView::Discs;  // the live disc preview first
        return;
    }
    ui.stage = Stage::Knit;
    ui.view.mode = ViewMode::Knit;
    if (firstChoice && !ui.results.knit) startKnit(ui);  // a new object opens from the scene cache when it holds these settings
}

void openMeshFile(SceneUi& ui, const std::string& path) {
    ui.import.session.open(path);
    ui.import.fileOpen = true;
    ui.import.keptAfterFailedRebuild = false;
    ui.import.revalidating = false;
    ui.import.keepKnitChoices = false;  // another file, or the same one opened anew: a scene switch
    ui.import.unit = UnitChoice::Unset;
    ui.import.unitTouched = false;
    ui.import.customUnitMm = 1;
    ui.import.sizeByHeight = false;
    ui.import.focusedNode.clear();
    ui.import.revealFocus = false;
    ui.import.options = knit::ImportOptions{};
    ui.workflowChosen = false;  // a newly opened model: Knit or Fur is asked (changeImport keeps the choice)
    if (!knit::furEnabled()) {  // the Fur switch is off: Knit is the only workflow, nothing is asked
        ui.s.workflow = knit::toString(knit::Workflow::Knit);
        ui.workflowChosen = true;
    }
    applySurface(ui);          // every knitted mesh is rebuilt; fur grows on the file's mesh
    ui.import.options.wrapStitchMm = ui.s.wale;
    ui.stage = Stage::Import;
    if (ui.prefs) {
        ui.prefs->lastImportDir = fs::path(path).parent_path().string();
        saveViewerPrefs(*ui.prefs);
    }
}

void changeImport(SceneUi& ui) {
    const knit::ImportOptions stored = knit::importOptionsFrom(ui.s);
    const bool chosen = ui.workflowChosen;
    openMeshFile(ui, ui.s.importSource);
    ui.workflowChosen = chosen;  // the same object: its workflow stays, and the surface follows it
    applySurface(ui);
    ui.import.options.plainNodes = stored.plainNodes;  // the choices stay; units and axis are suggested again
}

void prepareKnittingCharts(SceneUi& ui) {
    const bool chosen = ui.workflowChosen;
    openMeshFile(ui, ui.s.importSource);
    ui.workflowChosen = chosen;
    ui.import.options = knit::importOptionsFrom(ui.s);
    ui.import.options.surface = knit::ImportOptions::Surface::Wrapped;
    ui.import.options.wrapStitchMm = ui.s.wale;
    ui.import.unit = UnitChoice::Other;
    ui.import.customUnitMm = ui.s.importUnitMm;
    ui.import.unitTouched = true;
}

void reviewImport(SceneUi& ui) {
    prepareKnittingCharts(ui);
    ui.import.keepKnitChoices = true;
}

void prepareFileMesh(SceneUi& ui) {
    std::error_code ec;
    if (ui.s.importSource.empty() || !fs::exists(ui.s.importSource, ec)) {
        ui.import.note = "The source file " + ui.s.importSource + " is missing: the fur grows on the rebuilt surface until it is imported again.";
        return;
    }
    knit::ImportOptions o = knit::importOptionsFrom(ui.s);
    o.surface = knit::ImportOptions::Surface::KeepFile;
    ui.import.session.open(ui.s.importSource);  // prepared once read (importAnalysed)
    ui.import.options = o;
    ui.import.revalidating = false;
    ui.import.preparingFileMesh = true;  // the same object: importPrepared keeps every setting
}

void presetImportUnits(SceneUi& ui, double unitMm, knit::UpAxis up) {
    ui.import.options.unitMm = unitMm;
    ui.import.options.up = up;
    ui.import.unit = UnitChoice::Other;
    ui.import.customUnitMm = unitMm;
    ui.import.unitTouched = true;
}

void confirmImport(SceneUi& ui, const ImportChoices& choices) {
    if (choices.keepFile) ui.import.keepFileSurface = true;
    knit::ImportOptions& o = ui.import.options;
    applySurface(ui);
    o.wrapStitchMm = ui.s.wale;
    o.plainNodes = choices.plainNodes;
    startImport(ui);
}

void closeImport(SceneUi& ui) {
    ui.import.session.reset();
    ui.import.fileOpen = false;
    ui.import.keepKnitChoices = false;
    ui.import.focusedNode.clear();
    ui.import.revealFocus = false;
}

void startImport(SceneUi& ui) { ui.import.session.prepare(importRequest(ui.import)); }

void chooseImportUnit(SceneUi& ui, UnitChoice unit) {
    ImportFlow& im = ui.import;
    if (unit == UnitChoice::Unset) return;
    // Custom starts from the size on screen, so choosing it changes nothing until a value is typed.
    if (unit == UnitChoice::Other && im.options.unitMm > 0) im.customUnitMm = im.options.unitMm;
    im.unit = unit;
    im.unitTouched = true;
    im.sizeByHeight = false;
    im.options.unitMm = fileUnitMm(im);
}

void setCustomImportUnit(SceneUi& ui, double mmPerUnit) {
    ui.import.customUnitMm = mmPerUnit;
    if (!ui.import.sizeByHeight) ui.import.options.unitMm = fileUnitMm(ui.import);  // 0 for a value that is no size
}

// Set height keeps the size on screen (its height is the starting value) and never touches the file unit
// or the Custom value, so File units brings back exactly what was there.
void chooseSizeByHeight(SceneUi& ui, bool byHeight) {
    ImportFlow& im = ui.import;
    if (byHeight == im.sizeByHeight) return;
    im.unitTouched = true;
    if (!byHeight) {
        im.sizeByHeight = false;
        im.options.unitMm = fileUnitMm(im);
        return;
    }
    const ImportSize size = importSize(ui);
    if (size.valid && size.height > 0) im.targetHeightMm = size.height;
    im.sizeByHeight = true;
    unitFromHeight(ui);
}

void setImportHeight(SceneUi& ui, double heightMm) {
    ui.import.targetHeightMm = heightMm;
    unitFromHeight(ui);
}

void chooseImportUp(SceneUi& ui, knit::UpAxis up) {
    ui.import.options.up = up;
    unitFromHeight(ui);
}

void toggleNode(SceneUi& ui, const std::string& name) {
    std::set<std::string>& plain = ui.import.options.plainNodes;
    if (plain.erase(name) == 0) plain.insert(name);
}

void focusNode(SceneUi& ui, const std::string& name) { ui.import.focusedNode = name; }

// A click on the object only points at a component (its row comes into view); the tick is the checkbox's.
// The background clears the focus.
void pickNode(SceneUi& ui, const Camera& camera, const ViewPixel& pixel) {
    const knit::ImportSession::Snapshot file = ui.import.session.snapshot();
    if (!ui.import.fileOpen || !file.mesh || file.phase == knit::ImportSession::Phase::Preparing) return;
    const std::optional<int> node = nodeUnderPixel(ui, camera, pixel);
    ui.import.focusedNode = node ? file.mesh->nodes[*node].name : std::string();
    ui.import.revealFocus = node.has_value();
}

void selectAllNodes(SceneUi& ui) { ui.import.options.plainNodes.clear(); }

const std::vector<knit::GlbNode>& presetNodes(SceneUi& ui, const knit::SceneEntry& e) {
    if (ui.import.presetNodesOf != e.id) {
        ui.import.presetNodesOf = e.id;
        std::string err;
        const std::string path = knit::findProjectFile(e.path);
        if (path.empty() || !knit::listGlbNodes(path, e.glb, ui.import.presetNodes, &err)) ui.import.presetNodes.clear();
    }
    return ui.import.presetNodes;
}

void selectNoNodes(SceneUi& ui) {
    const std::shared_ptr<const knit::ImportedMesh> mesh = ui.import.session.snapshot().mesh;
    if (!mesh) return;
    for (const knit::ImportNode& node : mesh->nodes) ui.import.options.plainNodes.insert(node.name);
}

void adoptSettings(SceneUi& ui) {
    resetDestination(ui);
    if (ui.s.scene != "import") return;
    // The import named in the settings: its source is read again; the same identity uses the
    // stored entry, another one (the file, its MTL or images changed) is prepared again with the
    // stored options; a missing source falls back to the stored copy, labelled.
    std::error_code ec;
    if (!fs::exists(ui.s.importSource, ec)) {
        if (knit::sceneEntryFor(ui.s)) ui.import.note = "The source file " + ui.s.importSource + " is missing: using the prepared copy.";
        else ui.import.note = "The source file " + ui.s.importSource + " and its prepared copy are missing: import it again.";
        return;
    }
    ui.import.session.open(ui.s.importSource);
    ui.import.revalidating = true;
    ui.import.options = knit::importOptionsFrom(ui.s);
}

void startKnit(SceneUi& ui) {
    if (!entryOf(ui) || ui.job.active || !ui.workflowChosen || !knit::knitWorkflow(ui.s)) return;
    knit::KnitRequest r;
    const knit::KnitResult* shown = knitResult(ui);
    r.knitAgain = shown && !shown->computational;  // "Knit again" on a view-only result
    r.settings = jobStarting(ui, knit::BuiltScene::Kind::Knit, false);
    r.cacheRoot = ui.cacheRoot;
    r.loadKnitCache = ui.loadKnit;
    ui.builder->startKnit(std::move(r));
    ui.view.mode = ViewMode::Knit;
}

void startFibres(SceneUi& ui) {
    if (ui.job.active || !fibresInputReady(ui)) return;
    knit::FibreRequest r;
    if (knit::knitWorkflow(ui.s)) r.knit = ui.results.knit->knit;
    r.settings = jobStarting(ui, knit::BuiltScene::Kind::Fibres, false);
    ui.builder->startFibres(std::move(r));
    // "Apply to model": the preview target stays (the reference patch keeps following the edits);
    // the model view shows the fibres when it is chosen.
    ui.view.mode = ViewMode::Fibres;
}

void startExport(SceneUi& ui) {
    if (ui.job.active || !exportTargetProblem(ui).empty()) return;
    const knit::ExportFormat* f = knit::findExportFormat(ui.exportTo.format);
    knit::ExportRequest r;
    const bool fur = knit::furWorkflow(ui.s);
    r.knit = !fur && ui.results.knit ? ui.results.knit->knit : nullptr;
    r.includeFibres = (ui.s.fibres || fur) && f && f->fibres;
    // Fur: the result also carries the base meshes the export writes.
    r.fibres = (r.includeFibres || fur) && fibreResult(ui) ? ui.results.fibres->fibreResult : nullptr;
    r.dir = ui.exportTo.dir;
    r.name = ui.exportTo.name;
    r.format = ui.exportTo.format;
    r.settings = jobStarting(ui, knit::BuiltScene::Kind::Export, false);
    ui.builder->startExport(std::move(r));
}

// Rebuilds the shown knit's yarn tubes at the current tessellation (no knitting).
void startStrandMesh(SceneUi& ui) {
    knit::StrandMeshRequest r;
    r.shown = ui.results.knit;
    r.settings = jobStarting(ui, knit::BuiltScene::Kind::Knit, true);
    r.cacheRoot = ui.cacheRoot;
    ui.builder->startStrandMesh(std::move(r));
}

void startYarnMapBake(SceneUi& ui) {
    if (!ui.results.knit || ui.job.active) return;  // the button appears with a shown knit; a --do step may come early
    knit::StrandMeshRequest r;
    r.shown = ui.results.knit;
    r.settings = jobStarting(ui, knit::BuiltScene::Kind::Knit, true);
    ui.job.bakeMaps = true;
    r.cacheRoot = ui.cacheRoot;
    r.bakeMaps = true;
    ui.builder->startStrandMesh(std::move(r));
}

void revealYarnFibres(SceneUi& ui, Reveal what) {
    ui.stage = Stage::Fibres;
    ui.panels.revealYarnFibres = what;
}

void cancelJob(SceneUi& ui) {
    if (!ui.job.active) return;
    ui.builder->cancel();
    ui.job.cancelRequested = true;
    ui.job.stoppingNote = ui.builder->snapshot().progress.stage;
}

void knitAtFinalQuality(SceneUi& ui) {
    const knit::SceneEntry* e = entryOf(ui);
    if (!e) return;
    if (!ui.workflowChosen || knit::furWorkflow(ui.s)) chooseWorkflow(ui, knit::Workflow::Knit);  // asked for a Final knit
    knit::applyFullProfile(ui.s, *e);
    ui.stage = Stage::Knit;
    startKnit(ui);
}

void knitAtDraftQuality(SceneUi& ui) {
    if (!entryOf(ui)) return;
    knit::applyPreviewProfile(ui.s);
    ui.stage = Stage::Knit;
    startKnit(ui);
}

void runExportFix(SceneUi& ui, const knit::ExportEligibility& el) {
    if (el.fix == "Knit at Final quality" && entryOf(ui)) {
        knitAtFinalQuality(ui);
    } else if (el.fix == "Generate fibres" || el.fix == "Update fibres" || el.fix == "Grow fur" || el.fix == "Update fur") {
        ui.stage = Stage::Fibres;
        startFibres(ui);
    } else if (el.fix == kChooseWorkflow) {
        ui.stage = Stage::Import;  // the workflow is chosen there; nothing to compute
    } else if (el.fix == "Correct the yarn's fibres") {
        // The settings are changed in the Fibres step (Yarn fibres); nothing to compute. A refused body control
        // opens its node; a refused map texel (or an accepted body) the section.
        const knit::YarnBodyCheck body = knit::checkYarnBody(ui.s);
        revealYarnFibres(ui, body.accepted ? Reveal::Section : revealOfBodyField(body.field));
    } else {
        ui.stage = Stage::Knit;
        startKnit(ui);
    }
}

// Sets the level's yarn tessellation and segment budget; rebuilds the shown tubes as the Knit step's
// Yarn tessellation does (no knitting). True if anything changed.
bool applyExportDetail(SceneUi& ui, const ExportDetail& d) {
    const knit::FilePreset* level = meshLevelOf(d);
    if (!level || exportDetailSelected(ui, d)) return false;
    knit::applyPreset(knit::PresetGroup::StrandMesh, ui.s, *level);
    ui.s.segmentBudget = d.budget;
    ui.panels.customTessellation = false;
    if (!ui.job.active && shownTubesDiffer(ui)) startStrandMesh(ui);
    return true;
}

void chooseExportFolder(SceneUi& ui, const std::string& dir) {
    ui.exportTo.dir = fs::path(dir).make_preferred().string();
    ui.exportTo.dirChosen = true;
}

void loadMaskImage(SceneUi& ui, const std::string& path) {
    ui.s.pile.maskPath = path;
    ui.s.pile.maskHash = knit::fileHashHex(path);
}

void reloadMaskImage(SceneUi& ui) { ui.s.pile.maskHash = knit::fileHashHex(ui.s.pile.maskPath); }

bool checkMaskImage(SceneUi& ui) {
    const std::string& path = ui.s.pile.maskPath;
    if (path.empty()) return false;
    MaskOnDisk& disk = ui.panels.maskOnDisk;
    std::error_code ec;
    const fs::file_time_type time = fs::last_write_time(path, ec);
    const std::uintmax_t size = ec ? 0 : fs::file_size(path, ec);
    if (ec) return knit::fileHashHex(path) != ui.s.pile.maskHash;  // deleted or unreadable: "" (fileHashHex); the check keeps what it had
    if (path != disk.path || time != disk.writeTime || size != disk.size) {
        disk.path = path;
        disk.writeTime = time;
        disk.size = size;
        disk.hash = knit::fileHashHex(path);
    }
    return disk.hash != ui.s.pile.maskHash;
}

void takeFibreResult(SceneUi& ui, std::shared_ptr<const knit::BuiltScene> done) {
    const bool wool = done && done->fibreResult && done->fibreResult->fibres > 0 && knit::knitWorkflow(done->settings);
    ui.results.fibres = std::move(done);
    if (wool) ui.s.fibres = true;
}

namespace {

// The running job's partial scene (a knit shows its parts as they finish) and its end: the stage's
// attempt records the outcome; a cancelled or failed knit keeps the last completed result.
void integrateJobEnd(SceneUi& ui, const knit::JobSnapshot& job) {
    if (ui.job.active && ui.job.kind == knit::BuiltScene::Kind::Knit) {
        const std::shared_ptr<const knit::BuiltScene>& latest = job.latest;
        if (latest && latest->partial && latest != ui.results.partial) ui.results.partial = latest;
    }
    if (ui.job.active && !job.busy) {
        ui.job.active = false;
        const knit::JobOutcome oc = job.lastOutcome;
        Attempt& a = ui.job.kind == knit::BuiltScene::Kind::Knit ? ui.results.knitAttempt
                     : ui.job.kind == knit::BuiltScene::Kind::Fibres ? ui.results.fibreAttempt
                                                                       : ui.results.exportAttempt;
        a.outcome = oc;
        a.message = oc == knit::JobOutcome::Failed ? job.lastError : std::string();
        if (oc == knit::JobOutcome::Cancelled) {
            char b[200];
            if (ui.job.kind == knit::BuiltScene::Kind::Knit && ui.results.knit)
                std::snprintf(b, sizeof b, "Cancelled - showing the previous result #%d.", ui.results.knit->knit ? ui.results.knit->knit->id : 0);
            else if (ui.job.kind == knit::BuiltScene::Kind::Knit)
                std::snprintf(b, sizeof b, "Cancelled - the parts knitted so far stay visible as an incomplete preview (not usable).");
            else std::snprintf(b, sizeof b, "Cancelled.");
            a.message = b;
        }
        // A cancelled or failed knit keeps the last completed result on screen; without one, the
        // partial stays as an incomplete preview (never a result: it has no KnitResult).
        if (oc != knit::JobOutcome::Done && ui.results.knit) ui.results.partial.reset();
        ui.job.cancelRequested = false;
    }
}

// A result the builder completed: it becomes its stage's shown result.
void integrateCompletedResult(SceneUi& ui, const knit::JobSnapshot& job) {
    if (const std::shared_ptr<const knit::BuiltScene>& done = job.latestComplete; done && done != ui.results.lastSeen) {
        ui.results.lastSeen = done;
        Attempt* a = nullptr;
        switch (done->kind) {
        case knit::BuiltScene::Kind::Knit: {
            const bool first = !ui.results.knit;
            ui.results.knit = done;
            ui.results.partial.reset();
            a = &ui.results.knitAttempt;
            if (first) ui.camera.frame = true;
            break;
        }
        case knit::BuiltScene::Kind::Fibres:
            takeFibreResult(ui, done);
            a = &ui.results.fibreAttempt;
            break;
        case knit::BuiltScene::Kind::Export:
            ui.results.exported = done;
            a = &ui.results.exportAttempt;
            break;
        }
        if (a) a->outcome = knit::JobOutcome::Done;
    }
}

// One progress snapshot per frame: the Knit footer, the Fibres panel and the status bar show the
// same numbers.
void takeProgressSnapshot(SceneUi& ui, const knit::JobSnapshot& job) {
    ui.job.progress = ui.job.active ? job.progress : knit::JobProgress{};
    if (ui.job.active) {
        ui.job.cancelRequested = ui.job.cancelRequested || job.cancelRequested;
        if (ui.job.cancelRequested) ui.job.stoppingNote = ui.job.progress.stage;
    }
}

// A file analysed: the import named in the settings is re-validated; a new one gets its suggestions.
void importAnalysed(SceneUi& ui) {
    ImportFlow& im = ui.import;
    const knit::ImportSession::Snapshot file = im.session.snapshot();
    const std::shared_ptr<const knit::ImportAnalysis>& an = file.analysis;
    im.preview = file.preview;
    if (im.preparingFileMesh) {
        if (an && !an->fatal()) im.session.prepare(im.options);
        else {
            im.preparingFileMesh = false;
            im.note = "The file's own mesh cannot be prepared (the file has a problem): the fur grows on the rebuilt surface.";
        }
        return;
    }
    if (im.revalidating) {
        // The same identity keeps the stored entry.
        if (im.session.storedImportMatches(im.options, ui.s.importIdentity)) {
            im.revalidating = false;
        } else if (an && !an->fatal()) {
            im.note = "The source file or its material / textures changed since it was imported: it is prepared again with the same "
                      "options, and the knit needs an update.";
            im.session.prepare(im.options);
        } else {
            im.fileOpen = true;  // something the user must decide
            im.revalidating = false;
        }
        return;
    }
    ui.camera.frame = true;  // the object is on show, and clickable, while its choices are made
    if (im.unit != UnitChoice::Unset || !an) return;
    const double u = an->suggestedUnitMm;
    im.unit = u == 1 ? UnitChoice::Mm : u == 10 ? UnitChoice::Cm : u == 1000 ? UnitChoice::M : UnitChoice::Inch;
    im.options.unitMm = u;
    im.targetHeightMm = 300;
}

// A file prepared: it becomes the scene.
void importPrepared(SceneUi& ui) {
    ImportFlow& im = ui.import;
    const knit::ImportSession::Snapshot file = im.session.snapshot();
    knit::SceneSettings next = im.session.applyTo(ui.s);
    const knit::SceneEntry* e = knit::sceneEntryFor(next);
    if (!e) return;
    // A re-validation or the Knit step's Review Import: the same object, prepared again. applyTo wrote only
    // the import's own keys, so every Knit choice (colour, row direction, gauge, character) stays.
    const bool sameScene = (im.revalidating || im.keepKnitChoices || im.preparingFileMesh) && next.importSource == ui.s.importSource;
    if (sameScene) {
        ui.s = next;  // the knit and fibres stay, marked "needs update" (the asset changed)
        if (im.keepKnitChoices) {
            ui.stage = Stage::Knit;
            ui.view.mode = ViewMode::Mesh;  // the surface just prepared
        }
    } else {
        // A rebuilt surface carries its baked colours per stitch, even from a file with flat colours.
        const std::shared_ptr<const knit::ImportAnalysis>& an = file.analysis;
        const bool rebuilt = next.importSurface == "wrapped" && an && an->colourSource != "none";
        next.yarnColourSource = rebuilt ? "texture" : an ? knit::defaultYarnColourSource(*an) : "single";
        next.pileColourSource = next.yarnColourSource == "texture" ? "texture" : "object";
        switchScene(ui, next, *e);
        ui.workflowChosen = true;  // Import is enabled only after the choice (importFooterView)
        resetDestination(ui);
        ui.stage = knit::furWorkflow(ui.s) ? Stage::Fibres : Stage::Knit;
        ui.view.mode = ViewMode::Mesh;
    }
    im.changes = file.changes;
    im.fileOpen = false;
    im.revalidating = false;
    im.keepKnitChoices = false;
    im.preparingFileMesh = false;
    ui.camera.frame = true;
}

// The import session's new phase.
void integrateImport(SceneUi& ui) {
    const knit::ImportSession::Snapshot file = ui.import.session.snapshot();
    if (file.generation == ui.import.generation) return;
    ui.import.generation = file.generation;
    if (file.phase == knit::ImportSession::Phase::Analysed) {
        importAnalysed(ui);
    } else if (file.phase == knit::ImportSession::Phase::Prepared) {
        importPrepared(ui);
    } else if (file.phase == knit::ImportSession::Phase::Failed && (ui.import.revalidating || ui.import.preparingFileMesh)) {
        ui.import.note = "The source could not be read again (" + file.error + "): using the prepared copy.";
        ui.import.revalidating = false;
        ui.import.preparingFileMesh = false;
    }
}

// An import that kept the file's mesh (prepared for Fur) is rebuilt before it knits: the viewer knits only
// rebuilt surfaces. Scripted runs of knit_check may still knit the file's own uv.
void refreshImportCharts(SceneUi& ui) { ui.import.lacksCharts = ui.s.scene == "import" && ui.s.importSurface == "file"; }

// What follows the settings without a job: the flat-swatch references, the Fibres step's reference
// sample and the fur coverage overlay.
void updatePreviews(SceneUi& ui) {
    ui.previews.flat.request(ui.s);
    // The reference preview while the Fibres step shows it, after every edit, also while the model's
    // fibres are generated (the sample is small: about a second).
    if (ui.stage == Stage::Fibres && ui.view.fibres == FibreView::Reference && entryOf(ui)) ui.previews.fibres.request(ui.s);
    ui.previews.shown = ui.previews.fibres.latest();
    // The disc preview while the Fur step shows it, after every edit (coalesced: fur_preview.h).
    if (ui.stage == Stage::Fibres && ui.view.fibres == FibreView::Discs && furShown(ui) && entryOf(ui))
        ui.previews.discs.request(ui.s, ui.view.discColour);
    ui.previews.discShown = ui.previews.discs.latest();
    // The coverage overlay of the fur result (rebuilt when the result or the coverage changes).
    if (furShown(ui) && ui.view.coverage && ui.results.fibres && fibreResult(ui)) {
        const std::string key = std::to_string(reinterpret_cast<uintptr_t>(ui.results.fibres.get())) + "|" + ui.s.pile.coverage + "|" +
                                ui.s.pile.maskPath + "|" + ui.s.pile.maskHash + "|" + ui.s.pile.excludedParts;
        if (key != ui.coverage.key) {
            ui.coverage.key = key;
            ui.coverage.error.clear();
            std::string err;
            ui.coverage.scene = knit::furCoverageOverlay(*ui.results.fibres, ui.s, &ui.coverage.share, &err);
            if (!ui.coverage.scene) ui.coverage.error = err;
        }
    }
}

} // namespace

bool canKeepFileAfterFailedRebuild(const SceneUi& ui) {
    const knit::ImportSession::Snapshot file = ui.import.session.snapshot();
    return file.phase == knit::ImportSession::Phase::Failed && file.analysis && ui.import.fileOpen && !ui.import.keptAfterFailedRebuild &&
           !ui.import.revalidating && !ui.import.preparingFileMesh && ui.import.options.surface == knit::ImportOptions::Surface::Wrapped;
}

void keepFileAfterFailedRebuild(SceneUi& ui) {
    if (!canKeepFileAfterFailedRebuild(ui)) return;
    const knit::ImportSession::Snapshot file = ui.import.session.snapshot();
    ui.import.keptAfterFailedRebuild = true;
    ui.import.note = "The surface could not be rebuilt for knitting (" + file.error + "). You kept the file's own surface: the Mark step and "
                     "the guided stitch mesh work on it, surface tiling does not.";
    ui.import.options.surface = knit::ImportOptions::Surface::KeepFile;
    startImport(ui);
}

// Once per frame, in this order: the job's end, its completed result, its progress (all three from
// one snapshot of the builder, so they agree), the import, the previews.
// --- the live render preview

namespace {

// The ppt.exe the preview runs: KNIT_PPT_EXE, else the build-time default.
std::string previewRenderer() {
    if (const std::optional<std::string> chosen = knit::envVar("KNIT_PPT_EXE")) return *chosen;
#ifdef KNIT_DEFAULT_PPT_EXE
    return KNIT_DEFAULT_PPT_EXE;
#else
    return "";
#endif
}

uint64_t hashText(const std::string& text) {
    uint64_t h = 1469598103934665603ull;
    for (unsigned char c : text) h = (h ^ c) * 1099511628211ull;
    return h == 0 ? 1 : h;
}

// The settings a preview package is made with: the current ones (materials, export choices), with the
// knitting settings of the knit when they changed since (the preview shows that knit). Should the copy not
// give the knit's geometry, the knit's settings as a whole.
knit::SceneSettings previewSettings(const knit::SceneSettings& current, const knit::KnitResult& knitted) {
    knit::SceneSettings s = current;
    for (const knit::SettingField& f : knit::settingFields())
        if (f.stage == knit::Stage::Load || f.stage == knit::Stage::Knit || f.stage == knit::Stage::Surface)
            knit::setSetting(s, f.key, knit::settingValue(knitted.settings, f), nullptr);
    return knit::knitOutdated(knitted, s) ? knitted.settings : s;
}

// What the preview waits for or does, said over its image; "" when the image is current.
std::string previewBusy(const SceneUi& ui, const knit::RenderPreviewState& s, bool editsPending) {
    if (ui.job.active) return progressView(ui.job, false).title + " - the preview follows when it ends";
    if (s.phase == "exporting") return "Preparing the scene for the renderer...";
    if (s.phase == "loading") return "Loading the scene into the renderer...";
    if (s.phase == "starting") return "Starting the renderer...";
    if (s.phase == "waiting" || s.phase == "released" || s.phase == "suspended") return "Paused for a GPU job...";
    if (s.phase == "refused") return "Not enough GPU memory for the preview just now - it tries again";
    if (editsPending) return "Updating...";
    return "";
}

// A size the renderer gets: a multiple of 8, at least 64.
int previewSide(float pixels, int percent) {
    const int v = static_cast<int>(pixels * percent / 100.0f) / 8 * 8;
    return std::max(64, v);
}

} // namespace

void updateRenderPreview(SceneUi& ui, const Camera& camera, bool moving, double now) {
    RenderPreviewView& v = ui.renderPreview;
    knit::RenderPreviewRequest r;
    r.enabled = ui.panels.previewOpen;
    r.paused = v.paused;
    r.ppt = previewRenderer();
    const std::filesystem::path root = std::filesystem::temp_directory_path() / ("knit-preview-" + std::to_string(knit::currentProcessId()));
    r.logPath = (root / "ppt.log").string();
    r.packageRoot = (root / "packages").string();
    r.textureCache = (std::filesystem::temp_directory_path() / "knit-preview-textures").string();
    r.width = previewSide(v.imageW > 0 ? v.imageW : 512, v.scalePercent);
    r.height = previewSide(v.imageH > 0 ? v.imageH : 320, v.scalePercent);
    v.hint.clear();
    const knit::SceneEntry* e = entryOf(ui);
    if (e) {
        // Whatever results there are, drafts included (a preview is no deliverable: the Export step's
        // Final-quality rule does not apply). A knit made with other knitting settings is rendered with the
        // settings it was made with, and the window says so.
        const bool fur = knit::furWorkflow(ui.s);
        const std::shared_ptr<const knit::KnitResult> knitResult = !fur && ui.results.knit ? ui.results.knit->knit : nullptr;
        const std::shared_ptr<const knit::FibreResult> fibres = fibreResult(ui) ? ui.results.fibres->fibreResult : nullptr;
        const bool includeFibres = (ui.s.fibres || fur) && fibres;
        const bool haveResults = !fur ? knitResult != nullptr : (fibres && fibres->base);
        if (!ui.workflowChosen) {
            v.hint = "choose Knit or Fur first";
        } else if (!haveResults) {
            v.hint = !fur ? "knit first: the preview renders the knitted result" : "grow the fur first";
        } else {
            const bool outdated = knitResult && knit::knitOutdated(*knitResult, ui.s);
            if (outdated) v.hint = "the knitting settings changed: the preview shows the knit as it was made";
            const knit::SceneSettings settings = outdated ? previewSettings(ui.s, *knitResult) : ui.s;
            const std::string text = knit::writeSettings(settings) + "|" + std::to_string(reinterpret_cast<uintptr_t>(knitResult.get())) + "|" +
                                     std::to_string(reinterpret_cast<uintptr_t>(fibres.get())) + (includeFibres ? "|f" : "|-");
            if (text != v.settingsText) {
                v.settingsText = text;
                v.settingsChangedAt = now;
                v.packageSettings = std::make_shared<const knit::SceneSettings>(settings);
            }
            if (now - v.settingsChangedAt > 0.5) v.packageKey = hashText(text);
        }
        r.package.settings = v.packageSettings;
        r.package.knit = knitResult;
        r.package.fibres = (includeFibres || fur) ? fibres : nullptr;
        r.package.includeFibres = includeFibres;
        r.packageKey = v.packageSettings ? v.packageKey : 0;
    }
    // The viewer's camera (library mm, Z up) in the package's frame (metres, Y up): the export places
    // every part at its centre, as the display does, so one mapping serves the whole scene.
    const knit::Vec3 eye = camera.eye(), target = camera.target, up = camera.up;
    const double eyeM[3] = {eye.x / 1000, eye.z / 1000, -eye.y / 1000};
    const double targetM[3] = {target.x / 1000, target.z / 1000, -target.y / 1000};
    const double upM[3] = {up.x, up.z, -up.y};
    for (int i = 0; i < 3; ++i) {
        r.eye[i] = eyeM[i];
        r.target[i] = targetM[i];
        r.up[i] = upM[i];
    }
    r.fovDeg = camera.fovDeg;
    r.interacting = moving;
    char key[256];
    std::snprintf(key, sizeof key, "%.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.4g", eyeM[0], eyeM[1], eyeM[2], targetM[0], targetM[1],
                  targetM[2], upM[0], upM[1], upM[2], r.fovDeg);
    r.cameraKey = hashText(key) ^ (moving ? 1u : 0u);
    const bool editsPending = e && now - v.settingsChangedAt <= 0.5;
    r.hold = ui.job.active || editsPending;
    ui.previews.render.request(r);
    v.state = ui.previews.render.state();
    v.busy = r.enabled ? previewBusy(ui, v.state, editsPending) : "";
}

void pollSceneUi(SceneUi& ui) {
    const knit::JobSnapshot job = ui.builder->snapshot();
    integrateJobEnd(ui, job);
    integrateCompletedResult(ui, job);
    takeProgressSnapshot(ui, job);
    integrateImport(ui);
    refreshImportCharts(ui);
    updatePreviews(ui);
    integrateMarking(ui);
}

} // namespace viewer
