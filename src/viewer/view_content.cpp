// What the viewport shows (view_content.h): the window area it is drawn in, the meshes and fibres of the
// view on show with their caption, the toolbar's stages and Display rows, the comparison with the source,
// the scenes whose buffers may stay on the GPU, the counts overlay, and picking: the ray of a click taken
// back through the import preview's transform (previewPlacement and rayInFile are each other's inverse).
// Every function takes const SceneUi&. No ImGui and no GL.
#include "viewer/view_content.h"
#include "jobs/built_scene.h"

#include <algorithm>
#include <array>
#include <cstdio>
#include <filesystem>
#include <string>

#include "jobs/results.h"
#include "settings/scene_registry.h"
#include "settings/settings_profiles.h"

namespace fs = std::filesystem;
using knit::Vec3;

namespace viewer {

Vec3 toLibrary(const PreviewTransform& p, const Vec3& v) { return (p.yUp ? Vec3{v.x, -v.z, v.y} : v) * p.scale + p.offset; }

// --- what the viewport shows

ViewportRect viewportRect(const SceneUi& ui, float width, float height) {
    ViewportRect r;
    r.x = ui.layout.panelWidth;
    r.y = ui.layout.stripHeight + ui.layout.toolbarHeight;
    r.w = std::max(1.0f, width - r.x);
    r.h = std::max(1.0f, height - r.y - ui.layout.statusHeight - ui.layout.viewStatusHeight);
    return r;
}

namespace {

struct MeshKinds {
    bool base, yarn, plain;
};

// The scene's meshes of these kinds; the scene gives the view its bounds and its albedo textures.
void addMeshes(ViewContent& c, const knit::BuiltScene& b, MeshKinds kinds) {
    for (const knit::DisplayMesh& d : b.meshes) {
        const bool isBase = d.kind == knit::DisplayMesh::Kind::Base, isYarn = d.kind == knit::DisplayMesh::Kind::Yarn;
        if ((isBase && kinds.base) || (isYarn && kinds.yarn) || (!isBase && !isYarn && kinds.plain)) c.meshes.push_back(&d);
    }
    c.meshOwner = &b;
    c.bounds = {b.lo, b.hi};
}

bool holdsKind(const knit::BuiltScene& b, knit::DisplayMesh::Kind kind) {
    for (const knit::DisplayMesh& d : b.meshes)
        if (d.kind == kind) return true;
    return false;
}



// The import preview in library space: the chosen unit and up axis, centred, floor at z = 0.
struct PreviewPlacement {
    PreviewTransform transform;
    Bounds bounds;
};

// A point of the file in library space (the transform's forward direction; rayInFile is its inverse).

// A file-coordinate box placed with the unit and up axis: centred, floor at z = 0.
PreviewPlacement placeFileBox(const Vec3& fileLo, const Vec3& fileHi, double scale, bool yUp) {
    PreviewTransform p;
    p.scale = scale;
    p.yUp = yUp;
    const Vec3 a = toLibrary(p, fileLo), b = toLibrary(p, fileHi);
    const Vec3 lo{std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z)}, hi{std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)};
    p.offset = {-(lo.x + hi.x) * 0.5, -(lo.y + hi.y) * 0.5, -lo.z};
    return PreviewPlacement{p, {lo + p.offset, hi + p.offset}};
}

std::optional<PreviewPlacement> previewPlacement(const SceneUi& ui) {
    if (!ui.import.preview) return std::nullopt;
    const double scale = ui.import.fileOpen ? ui.import.options.unitMm : ui.s.importUnitMm;
    const bool yUp = ui.import.fileOpen ? ui.import.options.up == knit::UpAxis::Y : ui.s.importUp != "z";
    return placeFileBox(ui.import.preview->lo, ui.import.preview->hi, scale, yUp);
}

// The Mark step's pieces, placed as the import is (the imported unit and up axis).
std::optional<PreviewPlacement> markPlacement(const SceneUi& ui) {
    const knit::MarkingSession::Snapshot snap = ui.mark.session.snapshot();
    if (!snap.preview) return std::nullopt;
    return placeFileBox(snap.preview->lo, snap.preview->hi, ui.s.importUnitMm, ui.s.importUp != "z");
}

// The Mark step: every piece (the isolated one opaque, the others as chosen), the last build of the isolated
// piece in the Rows and Stitch mesh views, and the markings over them. The bounds are the isolated piece's, so
// framing goes to it.
void markContent(ViewContent& c, const SceneUi& ui) {
    c.viewKey = "mark";
    const knit::MarkingSession::Snapshot snap = ui.mark.session.snapshot();
    const std::optional<PreviewPlacement> placed = markPlacement(ui);
    if (!placed || !snap.pieces) {
        c.caption = snap.phase == knit::MarkingSession::Phase::Failed ? "Cannot read the file: " + snap.error
                    : ui.s.scene == "import"                       ? "Reading the file ..."
                                                                   : "The Mark step works on an imported file";
        return;
    }
    c.preview = placed->transform;
    c.bounds = placed->bounds;
    c.overlay = &ui.mark.overlay;
    const int isolated = ui.mark.piece;
    const knit::MarkedPieceView* built = shownBuild(ui);
    const bool stitchView = ui.mark.view == MarkView::StitchMesh;
    // the Stitch Mesh view draws the stitch faces only when the build has a mesh to draw; without one the piece stays
    // as imported, and the caption says why there is no mesh (never an earlier build's faces)
    const bool stitchFaces = stitchView && built && built->scene && built->use.drawable && built->scene->meshes.size() > 2;
    const bool rowsShown = ui.mark.view == MarkView::Rows && built && built->scene;
    for (size_t i = 0; i < snap.preview->meshes.size(); ++i) {
        const bool own = int(i) == isolated;
        if (!own && isolated >= 0 && ui.mark.around == NeighbourLook::Hidden) continue;
        if (own && (rowsShown || stitchFaces)) continue;  // the build's own meshes stand in for it
        c.meshes.push_back(&snap.preview->meshes[i]);
        c.meshAlpha.push_back(!own && isolated >= 0 && ui.mark.around == NeighbourLook::Transparent ? 0.22f : 1.0f);
    }
    if (rowsShown || (stitchFaces && ui.mark.stitchSource)) {
        c.meshes.push_back(&built->scene->meshes[0]);
        c.meshAlpha.push_back(1.0f);
    }
    if (stitchFaces) {
        c.meshes.push_back(&built->scene->meshes[ui.mark.stitchHighlight ? 2 : 1]);
        c.meshAlpha.push_back(1.0f);
    }
    c.meshOwner = snap.preview.get();
    if (isolated < 0) {
        c.caption = "Click a piece to isolate it";
        return;
    }
    Vec3 lo, hi;
    snap.pieces->pieces[isolated].mesh.bounds(lo, hi);
    const Vec3 a = toLibrary(placed->transform, lo), b = toLibrary(placed->transform, hi);
    c.bounds = {{std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z)}, {std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)}};
    c.caption = ui.mark.pieceLabels[isolated];
    if (stitchView) {
        c.caption += " - " + (built ? knit::stitchMeshHeadline(built->use, built->quality) : std::string("No stitch mesh for this piece yet: press Build stitch mesh"));
        if (built && built->use.drawable && !built->quality.relaxed) c.caption += " (unrelaxed layout)";
        if (built && shownBuildOutdated(ui)) c.caption += " (outdated)";
        if (snap.phase == knit::MarkingSession::Phase::Building) c.caption += " (building ...)";
    } else if (built) {
        c.caption += built->use.drawable ? " - stitch mesh built" : " - not built (see the findings)";
    }
}

bool addImportPreview(ViewContent& c, const SceneUi& ui, const char* caption) {
    const std::optional<PreviewPlacement> placed = previewPlacement(ui);
    if (!placed) return false;
    c.preview = placed->transform;
    c.bounds = placed->bounds;
    for (const knit::DisplayMesh& d : ui.import.preview->meshes) c.meshes.push_back(&d);
    c.meshOwner = ui.import.preview.get();
    c.caption = caption;  // the file's name and size are in the status line (objectLine)
    return true;
}

// The caption's note when the display holds fewer fibres than the result.
std::string decimationNote(const knit::FibreResult& f, bool fur) {
    if (f.displayedFibres >= f.fibres) return "";
    return " - displaying " + fmtCount(f.displayedFibres) + " of " + fmtCount(f.fibres) + (fur ? " hairs" : " fibres") + "; the export includes all " +
           fmtCount(f.fibres) + " (the display looks thinner than the real density)";
}

void importFileContent(ViewContent& c, const SceneUi& ui) {
    if (addImportPreview(c, ui, "Import preview")) return;
    c.caption = ui.import.session.snapshot().phase == knit::ImportSession::Phase::Failed ? "Cannot import this file (see the Import step)" : "Reading the file ...";
}

// The Fibres step's reference sample (its own camera) is on show.
bool referenceShown(const SceneUi& ui) { return ui.stage == Stage::Fibres && ui.view.fibres == FibreView::Reference && entryOf(ui); }

bool discsShown(const SceneUi& ui) { return ui.stage == Stage::Fibres && ui.view.fibres == FibreView::Discs && furShown(ui) && entryOf(ui); }

// The disc preview: the object's meshes and one disc per loop fibre, one needle per fine fibre (jobs/fur_preview.h).
void discContent(ViewContent& c, const SceneUi& ui) {
    const std::shared_ptr<const knit::FurPreview>& p = ui.previews.discShown;
    const std::string error = ui.previews.discs.error();
    if (!p) {
        c.caption = error.empty() ? "Building the disc preview ..." : "Discs: " + error;
        return;
    }
    if (p->body) {
        addMeshes(c, *p->body, {true, false, true});
        c.bounds = {p->body->lo, p->body->hi};
        c.meshOwner = p->body.get();
    }
    c.furDiscs = p;
    char b[200];
    std::snprintf(b, sizeof b, "Discs: %s loops, %s fine fibres%s", fmtCount(p->loops).c_str(), fmtCount(p->fibres).c_str(),
                  p->keyBelow < 1 ? " (a quarter shown while editing)" : "");
    c.caption = b;
    if (ui.previews.discs.busy()) c.caption += " - updating ...";
    if (!error.empty()) c.caption += " - the latest settings cannot grow fur: " + error;
}

void referenceContent(ViewContent& c, const SceneUi& ui) {
    c.viewKey = furShown(ui) ? "reference-fur" : "reference-wool";
    const knit::FibrePreview* p = ui.previews.shown.get();
    const std::string error = ui.previews.fibres.error();
    if (p && p->workflow == knit::workflowOf(ui.s)) {
        for (const knit::DisplayMesh& d : p->meshes) c.meshes.push_back(&d);
        c.bounds = {p->lo, p->hi};
        const bool current = p->key == knit::fibrePreviewKey(ui.s);
        c.caption = (current ? std::string() : std::string("Updating ... (showing the previous settings) - ")) + p->note;
    } else {
        c.caption = error.empty() ? "Building the reference ..." : "Reference: " + error;
    }
    if (!error.empty() && p) c.caption += " - the latest settings failed: " + error;
}

// The fur result's scene, when there is one under Fur.
const knit::BuiltScene* hairScene(const SceneUi& ui) { return ui.results.fibres && fibreResult(ui) ? ui.results.fibres.get() : nullptr; }

// Fur on the model: the base meshes (or the coverage overlay) and the fur's fibres.
void hairContent(ViewContent& c, const SceneUi& ui) {
    const knit::FibreResult* result = fibreResult(ui);
    const bool growing = runningStage(ui, knit::BuiltScene::Kind::Fibres);
    const knit::BuiltScene* scene = hairScene(ui);
    if (scene && ui.view.coverage && ui.coverage.scene) scene = ui.coverage.scene.get();
    if (!scene) {
        if (ui.s.scene == "import" && addImportPreview(c, ui, "Imported mesh")) c.caption += growing ? " - growing fur ..." : " - no fur yet";
        else c.caption = growing ? "Growing fur ..." : "No fur yet: Grow fur grows it on the mesh";
        return;
    }
    addMeshes(c, *scene, {true, true, true});
    if (!ui.view.coverage) c.fibres = &ui.results.fibres->fibres;  // the coverage shows on the bare surface
    c.caption = "Fur #" + std::to_string(result->id) + ": " + fmtCount(result->fibres) + " fibres";
    if (growing) c.caption = "Updating ... (fur #" + std::to_string(result->id) + " shown)";
    else if (knit::fibresOutdated(*result, ui.s)) c.caption += " - fur settings changed since: grow it again";
    c.caption += decimationNote(*result, true);
    if (!ui.view.coverage) return;
    if (!ui.coverage.error.empty()) {
        c.caption += " - coverage: " + ui.coverage.error;
    } else {
        char text[160];
        std::snprintf(text, sizeof text, " - coverage shown (hairs hidden): red areas grow no hair (%.0f %% of the texture kept)", 100 * ui.coverage.share);
        c.caption += text;
    }
}

// The object as imported: the result's source surfaces (the knit's, or the hair's), or the import preview.
// `mode`: the stage asked for; the source stands in for a knit not made yet.
void sourceMeshContent(ViewContent& c, const SceneUi& ui, const knit::BuiltScene* resultScene, ViewMode mode) {
    const bool importScene = ui.s.scene == "import";
    const bool hasBase = resultScene && holdsKind(*resultScene, knit::DisplayMesh::Kind::Base);
    if (hasBase && !(ui.results.partial && importScene && ui.import.preview)) {
        addMeshes(c, *resultScene, {true, false, true});
        c.caption = importScene ? "Imported mesh" : "Source mesh";
    } else if (importScene && addImportPreview(c, ui, "Imported mesh")) {
    } else if (resultScene && !furShown(ui)) {
        addMeshes(c, *resultScene, {true, true, true});
        c.caption = "Source mesh";
    } else {
        c.caption = entryOf(ui) ? "Nothing to show yet: knit to build the scene" : "No object: open a mesh (Import)";
    }
    if (mode != ViewMode::Mesh && !furShown(ui) && !resultScene && entryOf(ui)) c.caption += " - not knitted yet";
}

// The knit's caption: which result, its verdicts, whether it is still current.
std::string knitCaption(const SceneUi& ui, const knit::KnitResult& k) {
    char text[256];
    std::snprintf(text, sizeof text, "Knitted result #%d - %s - %d PASS", k.id, k.final ? "Final" : knit::isPreviewProfile(k.settings) ? "Draft" : "Custom profile", k.pass);
    std::string caption = text;
    if (k.fail) caption += ", " + std::to_string(k.fail) + " FAIL";
    if (k.unresolved) caption += ", " + std::to_string(k.unresolved) + " UNRESOLVED";
    if (k.plain) caption += ", " + std::to_string(k.plain) + " plain";
    if (!knitCurrent(ui)) caption += k.computational ? " - settings changed since" : " - view only";
    return caption;
}

// The knitted yarn, and in the Fibres view its fibres. A partial knit keeps the unfinished parts'
// (dimmed) source surfaces, so the object stays whole.
void knitContent(ViewContent& c, const SceneUi& ui, const knit::BuiltScene& knitScene, ViewMode mode) {
    const bool partial = ui.results.partial != nullptr;
    addMeshes(c, knitScene, {partial, true, ui.view.unknittedParts});
    const knit::KnitResult* k = ui.results.knit ? ui.results.knit->knit.get() : nullptr;
    // The progress itself is in the Knit footer (one place).
    if (partial) c.caption = ui.job.active ? "Partial build (knitting, not a result)" : "Incomplete knit (cancelled, not a result)";
    else if (k) c.caption = knitCaption(ui, *k);
    if (mode != ViewMode::Fibres) return;
    const knit::FibreResult* f = shownFibres(ui);
    if (!f) {
        c.caption += " - no fibres for this knit";
        return;
    }
    c.fibres = &ui.results.fibres->fibres;
    c.caption += " + fibres #" + std::to_string(f->id) + " (" + fmtCount(f->fibres) + ")";
    if (runningStage(ui, knit::BuiltScene::Kind::Fibres)) c.caption += " - updating ...";
    c.caption += decimationNote(*f, false);
}

// The knit on show: the running knit's partial scene, else the last result.
const knit::BuiltScene* knitSceneShown(const SceneUi& ui) { return ui.results.partial ? ui.results.partial.get() : ui.results.knit.get(); }

// The chosen stage as the model view shows it: under Fur there is no knit, so a Knit chosen before shows the fur
// (the Fibres stage).
ViewMode shownMode(const SceneUi& ui) { return furShown(ui) && ui.view.mode == ViewMode::Knit ? ViewMode::Fibres : ui.view.mode; }

} // namespace

const knit::FibreResult* shownFibres(const SceneUi& ui) {
    if (furShown(ui) || ui.results.partial) return nullptr;
    const knit::FibreResult* f = fibreResult(ui);
    const knit::KnitResult* k = ui.results.knit ? ui.results.knit->knit.get() : nullptr;
    return f && k && f->knitId == k->id ? f : nullptr;
}

ViewContent stageContent(const SceneUi& ui, ViewMode mode) {
    ViewContent c;
    if (furShown(ui)) {  // the Fur workflow: the mesh, or the mesh with its fur
        if (mode == ViewMode::Fibres) hairContent(c, ui);
        else sourceMeshContent(c, ui, hairScene(ui), mode);
        return c;
    }
    const knit::BuiltScene* knitScene = knitSceneShown(ui);
    if (mode == ViewMode::Mesh || !knitScene) sourceMeshContent(c, ui, knitScene, mode);
    else knitContent(c, ui, *knitScene, mode);
    return c;
}

ViewContent viewContent(const SceneUi& ui) {
    if (ui.import.fileOpen) {
        ViewContent c;
        importFileContent(c, ui);
        return c;
    }
    if (ui.stage == Stage::Mark) {
        ViewContent c;
        markContent(c, ui);
        return c;
    }
    if (referenceShown(ui)) {
        ViewContent c;
        referenceContent(c, ui);
        return c;
    }
    if (discsShown(ui)) {
        ViewContent c;
        discContent(c, ui);
        return c;
    }
    return stageContent(ui, shownMode(ui));
}

// --- the toolbar's stages and Display rows

ViewMode markedViewMode(const SceneUi& ui) { return ui.import.fileOpen || ui.stage == Stage::Mark ? ViewMode::Mesh : shownMode(ui); }

namespace {

ViewModeState unavailable(const std::string& reason) { return {false, reason}; }

ViewModeState fibreViewState(const SceneUi& ui) {
    const bool computing = runningStage(ui, knit::BuiltScene::Kind::Fibres);
    if (furShown(ui)) {
        if (hairScene(ui)) return {true, ""};
        return unavailable(computing ? "Growing fur ..." : "No fur yet: Grow fur (the Fur step).");
    }
    if (shownFibres(ui)) return {true, ""};
    if (ui.results.partial) {
        const bool older = fibreResult(ui) != nullptr;
        if (ui.job.active) return unavailable(older ? "Knitting again: the fibres belong to the previous knit." : "Knitting ...");
        return unavailable(older ? "Incomplete knit: the fibres belong to the previous knit." : "Incomplete knit: knit again (step 1).");
    }
    if (!ui.results.knit) return unavailable("Not knitted yet: Knit (step 1), then Apply to model (step 2).");
    return unavailable(computing ? "Computing fibres ..." : "No fibres for this knit yet: Apply to model (step 2).");
}

// "" when the source and the chosen stage can be drawn side by side, else why not.
std::string compareProblem(const SceneUi& ui) {
    if (ui.import.fileOpen) return "The open file shows its preview until it is imported (Import step).";
    if (ui.stage == Stage::Mark) return "The Mark step shows the pieces and their markings.";
    if (referenceShown(ui)) return "The reference patch has no source: choose Model.";
    if (shownMode(ui) == ViewMode::Mesh) return "Choose Knit or Fibres to compare them with the source.";
    const ViewModeState stage = viewModeState(ui, shownMode(ui));
    if (!stage.available) return stage.reason;
    if (ui.results.partial && !furShown(ui)) return "A partial build is not a result.";
    const ViewContent source = stageContent(ui, ViewMode::Mesh);
    bool sourceSurface = source.preview.has_value();
    for (const knit::DisplayMesh* d : source.meshes) sourceSurface |= d->kind == knit::DisplayMesh::Kind::Base;
    if (!sourceSurface) return "This object has no source mesh to compare with.";
    return "";
}

} // namespace

ViewModeState viewModeState(const SceneUi& ui, ViewMode mode) {
    if (ui.import.fileOpen) return unavailable("The open file shows its preview until it is imported (Import step).");
    if (referenceShown(ui)) return unavailable("The reference patch has no stages: choose Model.");
    switch (mode) {
    case ViewMode::Mesh:
        if (!stageContent(ui, ViewMode::Mesh).meshes.empty()) return {true, ""};
        return unavailable(entryOf(ui) ? "Nothing to show yet: knit to build the scene." : "No object: open a mesh (Import).");
    case ViewMode::Knit:
        if (furShown(ui)) return unavailable("The Fur workflow knits nothing (step 1).");
        if (knitSceneShown(ui)) return {true, ""};
        return unavailable(runningStage(ui, knit::BuiltScene::Kind::Knit) ? "Knitting ..." : "Not knitted yet: Knit (step 1).");
    case ViewMode::Fibres: return fibreViewState(ui);
    }
    return unavailable("");
}

DisplayRows displayRows(const SceneUi& ui) {
    DisplayRows r;
    r.occlusion = ui.view.lighting == Lighting::Studio;
    if (ui.import.fileOpen) return r;
    const bool reference = referenceShown(ui);
    const knit::BuiltScene* knitScene = knitSceneShown(ui);
    // Judged on the scene, not on the view: hiding the unknitted parts must not hide their checkbox.
    const ViewMode mode = shownMode(ui);
    const bool knitShown = !furShown(ui) && !reference && mode != ViewMode::Mesh && knitScene;
    r.plyTexture = knitShown && holdsKind(*knitScene, knit::DisplayMesh::Kind::Yarn);
    if (ui.s.plies > 1) r.plyTextureWhyNot = "Separate ply tubes are real geometry: no twisted-ply maps.";
    if (ui.s.yarnBody.enabled) r.plyTextureWhyNot = "The yarn is made of fibres: no twisted-ply maps.";
    r.yarnMaps = knitShown && knit::carriesYarnMaps(*knitScene);
    r.unknittedParts = knitShown && holdsKind(*knitScene, knit::DisplayMesh::Kind::Plain);
    r.fibrePopulations = mode == ViewMode::Fibres && !furShown(ui) && !reference;
    r.compare = mode != ViewMode::Mesh && !reference;
    r.compareWhyNot = compareProblem(ui);
    return r;
}

ViewportFrame viewportFrame(const SceneUi& ui) {
    ViewportFrame f;
    f.shown = viewContent(ui);
    f.bounds = f.shown.bounds;
    if (!ui.view.compare || !compareProblem(ui).empty()) return f;  // off, or suspended: the flag stays
    f.source = stageContent(ui, ViewMode::Mesh);
    const Bounds& a = f.shown.bounds;
    const Bounds& b = f.source->bounds;
    f.bounds = {{std::min(a.lo.x, b.lo.x), std::min(a.lo.y, b.lo.y), std::min(a.lo.z, b.lo.z)},
                {std::max(a.hi.x, b.hi.x), std::max(a.hi.y, b.hi.y), std::max(a.hi.z, b.hi.z)}};
    return f;
}

ObjectLine objectLine(const SceneUi& ui) {
    ObjectLine o;
    if (ui.import.fileOpen || ui.s.scene == "import") {
        const std::string path = ui.import.session.snapshot().path;
        o.name = fs::path(path.empty() ? ui.s.importSource : path).filename().string();
    } else if (const knit::SceneEntry* e = entryOf(ui)) {
        o.name = e->label;
    }
    std::optional<Bounds> box;
    if (ui.import.fileOpen) {
        if (const std::optional<PreviewPlacement> placed = previewPlacement(ui)) box = placed->bounds;
    } else if (const ViewContent source = stageContent(ui, ViewMode::Mesh); !source.meshes.empty()) {
        box = source.bounds;
    }
    if (box) {
        char text[96];
        std::snprintf(text, sizeof text, "%.0f x %.0f x %.0f mm", box->hi.x - box->lo.x, box->hi.y - box->lo.y, box->hi.z - box->lo.z);
        o.size = text;
    }
    return o;
}

PopoverRect popoverRect(const ViewportRect& viewport, float anchorRight, float preferredWidth, float margin) {
    PopoverRect p;
    p.width = std::max(1.0f, std::min(preferredWidth, viewport.w - 2 * margin));
    const float left = viewport.x + margin, right = std::max(left, viewport.x + viewport.w - margin - p.width);
    p.x = std::clamp(anchorRight - p.width, left, right);
    p.y = viewport.y + margin;
    p.maxHeight = std::max(1.0f, viewport.h - 2 * margin);
    return p;
}

// --- the GPU's scenes

// Every scene viewContent can take its meshes from: their buffers stay uploaded across view changes.
MeshOwners meshOwners(const SceneUi& ui) {
    MeshOwners o;
    for (const std::shared_ptr<const knit::BuiltScene>& s : {ui.results.knit, ui.results.partial, ui.results.fibres, ui.import.preview, ui.coverage.scene})
        if (s) o.scenes.push_back(s);
    const knit::MarkingSession::Snapshot mark = ui.mark.session.snapshot();
    if (mark.preview) o.scenes.push_back(mark.preview);
    // every piece's last build, not only the newest: switching back to a piece draws its mesh again
    for (const MarkingFlow::ArrivedBuild& b : ui.mark.builds)
        if (b.view->scene) o.scenes.push_back(b.view->scene);
    o.preview = ui.previews.shown;
    if (ui.previews.discShown && ui.previews.discShown->body) o.scenes.push_back(ui.previews.discShown->body);
    return o;
}

namespace {

// The ray in the file's coordinates: the inverse of the preview's transform (loadCamera in renderer.cpp).
Ray rayInFile(const Ray& ray, const PreviewTransform& p) {
    const Vec3 o = (ray.origin - p.offset) * (1.0 / p.scale), d = ray.direction * (1.0 / p.scale);
    if (!p.yUp) return {o, d};
    return {{o.x, o.z, -o.y}, {d.x, d.z, -d.y}};
}

// Where the ray meets a triangle: how far along it, and the weights (u on b, v on c) of the point.
struct RayTriangleHit {
    double t = 0, u = 0, v = 0;
};

// How far along the ray it meets the triangle, from either side (Moller, Trumbore 1997).
std::optional<RayTriangleHit> rayMeetsTriangle(const Ray& ray, const Vec3& a, const Vec3& b, const Vec3& c) {
    const Vec3 e1 = b - a, e2 = c - a, p = cross(ray.direction, e2);
    const double det = dot(e1, p);
    if (det == 0) return std::nullopt;
    const Vec3 s = ray.origin - a, q = cross(s, e1);
    const double u = dot(s, p) / det, v = dot(ray.direction, q) / det, t = dot(e2, q) / det;
    if (u < 0 || v < 0 || u + v > 1 || t <= 0) return std::nullopt;
    return RayTriangleHit{t, u, v};
}

} // namespace

// Every triangle, nearest hit: one click's work, so no index to build or keep.
std::optional<int> nodeUnderPixel(const SceneUi& ui, const Camera& camera, const ViewPixel& pixel) {
    const std::shared_ptr<const knit::ImportedMesh> mesh = ui.import.session.snapshot().mesh;
    const ViewContent shown = viewContent(ui);
    if (!mesh || !shown.preview) return std::nullopt;
    const Ray ray = rayInFile(viewRay(camera, pixel), *shown.preview);
    std::optional<int> node;
    double nearest = 1e300;
    for (size_t t = 0; t < mesh->tri.size(); ++t) {
        const std::array<int, 3>& corners = mesh->tri[t];
        const std::optional<RayTriangleHit> hit = rayMeetsTriangle(ray, mesh->pos[corners[0]], mesh->pos[corners[1]], mesh->pos[corners[2]]);
        if (!hit || hit->t >= nearest) continue;
        nearest = hit->t;
        node = mesh->nodeOfTriangle[t];
    }
    return node;
}

// Every triangle of the pieces on show, nearest hit (one click's work, as nodeUnderPixel).
std::optional<SurfaceHit> surfaceUnderPixel(const SceneUi& ui, const Camera& camera, const ViewPixel& pixel) {
    const knit::MarkingSession::Snapshot snap = ui.mark.session.snapshot();
    const ViewContent shown = viewContent(ui);
    if (!snap.pieces || !shown.preview) return std::nullopt;
    const Ray ray = rayInFile(viewRay(camera, pixel), *shown.preview);
    std::optional<SurfaceHit> best;
    double nearest = 1e300;
    for (size_t pi = 0; pi < snap.pieces->pieces.size(); ++pi) {
        if (ui.mark.piece >= 0 && int(pi) != ui.mark.piece && ui.mark.around == NeighbourLook::Hidden) continue;
        const knit::SurfacePiece& piece = snap.pieces->pieces[pi];
        for (size_t t = 0; t < piece.mesh.tri.size(); ++t) {
            const std::array<int, 3>& c = piece.mesh.tri[t];
            const std::optional<RayTriangleHit> hit = rayMeetsTriangle(ray, piece.mesh.pos[c[0]], piece.mesh.pos[c[1]], piece.mesh.pos[c[2]]);
            if (!hit || hit->t >= nearest) continue;
            nearest = hit->t;
            SurfaceHit h;
            h.piece = int(pi);
            h.point.tri = piece.fileTriangle[t];
            h.point.bary[0] = 1 - hit->u - hit->v;
            h.point.bary[1] = hit->u;
            h.point.bary[2] = hit->v;
            h.position = knit::surfacePosition(piece.mesh, knit::SurfacePoint{int(t), {h.point.bary[0], h.point.bary[1], h.point.bary[2]}});
            best = h;
        }
    }
    return best;
}

bool previewNodeIsPlain(const SceneUi& ui, const knit::DisplayMesh& d) {
    return ui.import.fileOpen && nodeKeptAsMesh(ui.import, d.name);  // a preview mesh is named after its node
}

bool previewNodeIsFocused(const SceneUi& ui, const knit::DisplayMesh& d) {
    return ui.import.fileOpen && !ui.import.focusedNode.empty() && d.name == ui.import.focusedNode;
}

std::string viewCounts(const SceneUi& ui, const ViewContent& c) {
    long long tri[4] = {0, 0, 0, 0};  // by DisplayMesh::Kind: base, yarn, none, plain
    for (const knit::DisplayMesh* m : c.meshes) tri[static_cast<int>(m->kind)] += static_cast<long long>(m->idx.size() / 3);
    const long long total = tri[0] + tri[1] + tri[2] + tri[3];
    const bool reference = c.viewKey.rfind("reference", 0) == 0;
    std::string t = (reference ? "Reference patch: " : "Shown: ") + fmtCount(total) + " triangles";
    const char* kinds[4] = {"mesh", "yarn", "", "plain"};
    std::string split;
    for (int k = 0; k < 4; ++k)
        if (tri[k] && tri[k] != total) split += (split.empty() ? "" : ", ") + std::string(kinds[k]) + " " + fmtCount(tri[k]);
    if (!split.empty()) t += " (" + split + ")";
    if (!reference && c.meshOwner && !c.meshOwner->displayNote.empty()) t += " - coarser than the export (Knit step, Yarn tessellation)";
    if (reference) return t;
    if (const knit::FibreResult* f = fibreResult(ui)) {
        const char* what = f->workflow == knit::Workflow::Fur ? "fur fibres" : "fibres";
        t += "; " + fmtCount(f->fibres) + " " + what + " generated (" + fmtCount(f->segments) + " segments), ";
        t += c.fibres ? fmtCount(f->displayedFibres) + " shown" : std::string("not shown in this view");
    }
    return t;
}

} // namespace viewer
