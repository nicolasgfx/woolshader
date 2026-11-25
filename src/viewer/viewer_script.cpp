// The command line, the --do steps and when a shot is due (viewer_script.h).
#include "settings/scene_registry.h"
#include "settings/scene_settings.h"
#include "settings/settings_profiles.h"
#include "viewer/viewer_script.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>


namespace viewer {

std::string retiredOption(int argc, char** argv) {
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--lab" || a == "--pieces" || a == "--demo" || a == "--preset" || a == "--compaction") return a;
    }
    return "";
}

ViewerOptions parseViewerOptions(int argc, char** argv) {
    ViewerOptions o;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&] { return i + 1 < argc ? std::string(argv[++i]) : std::string(); };
        if (a == "--scene") o.scene = next();
        else if (a == "--shot") o.shotPath = next();
        else if (a == "--view") o.startView = cameraViewNamed(next()).value_or(CameraView::Front);
        else if (a == "--settings") o.settingsFile = next();
        else if (a == "--load-knit") o.loadKnit = next();
        else if (a == "--import") o.importFile = next();
        else if (a == "--units") o.units = next();
        else if (a == "--up") o.up = next();
        else if (a == "--plain-nodes") o.importChoices.plainNodes = knit::parsePlainNodes(next());
        else if (a == "--keep-file") o.importChoices.keepFile = true;  // the file's own surface, no rebuild (the guided route)
        else if (a == "--keep-file-if-rebuild-fails") o.keepFileIfRebuildFails = true;  // the user's choice after a failed rebuild
        else if (a == "--stage") o.stage = next();
        else if (a == "--do") o.steps = next();
        else if (a == "--size") {
            const std::string v = next();
            o.width = std::atoi(v.c_str());
            o.height = v.find('x') != std::string::npos ? std::atoi(v.c_str() + v.find('x') + 1) : 0;
        }
        else if (a == "--ui-scale") o.uiScale = static_cast<float>(std::atof(next().c_str()));
        else if (a == "--zoom") o.zoom = std::atof(next().c_str());
        else if (a == "--plain") o.plain = next() == "off" ? 0 : 1;
        else if (a == "--lighting") o.lighting = next();
        else if (a == "--occlusion") o.occlusion = std::atof(next().c_str());
        else if (a == "--bench") o.benchFrames = std::atoi(next().c_str());
        else if (a == "--target-frac") {
            const std::string v = next();
            sscanf_s(v.c_str(), "%lf,%lf,%lf", &o.targetFraction[0], &o.targetFraction[1], &o.targetFraction[2]);
        }
    }
    return o;
}

double parseUnits(const std::string& unit) {
    if (unit == "mm") return 1;
    if (unit == "cm") return 10;
    if (unit == "m") return 1000;
    if (unit == "in" || unit == "inch") return 25.4;
    return std::atof(unit.c_str());
}

Stage stageNamed(const std::string& name) {
    if (name == "import") return Stage::Import;
    if (name == "mark") return Stage::Mark;
    if (name == "fibres") return Stage::Fibres;
    if (name == "export") return Stage::Export;
    return Stage::Knit;
}

namespace {

void writeTextFile(const std::string& path, const std::string& text) {
    std::ofstream f(path, std::ios::binary);
    f << text;
    if (!f) std::fprintf(stderr, "%s: not written\n", path.c_str());
}

} // namespace

void applyScriptStep(SceneUi& ui, const std::string& st) {
    const knit::SceneEntry* e = knit::sceneEntryFor(ui.s);
    if (st == "final" && e) knit::applyFullProfile(ui.s, *e);
    else if (st == "draft") knit::applyPreviewProfile(ui.s);
    else if (st == "knit") startKnit(ui);
    else if (st == "bake-maps") startYarnMapBake(ui);  // the Knit step's Bake the yarn maps
    else if (st == "fibres") {
        ui.s.fibres = true;
        startFibres(ui);
    } else if (st == "export" || knit::findExportFormat(st)) {
        ui.exportTo.format = st == "export" ? std::string(knit::kPackageFormat) : st;
        startExport(ui);
    } else if (st == "mesh") ui.view.mode = ViewMode::Mesh;
    else if (st == "cancel") cancelJob(ui);
    else if (st == "knitview") ui.view.mode = ViewMode::Knit;
    else if (st == "fibreview") ui.view.mode = ViewMode::Fibres;
    else if (st == "reference") ui.view.fibres = FibreView::Reference;  // the Fibres step's views
    else if (st == "model") ui.view.fibres = FibreView::Model;
    // The object's workflow (the Import step's Knit / Fur choice).
    else if (st == "choose-fur" && !knit::furEnabled()) std::fprintf(stderr, "choose-fur: the Fur workflow is switched off (KNIT_FUR=1)\n");
    else if (st == "choose-fur") chooseWorkflow(ui, knit::Workflow::Fur);
    else if (st == "choose-knit") chooseWorkflow(ui, knit::Workflow::Knit);
    else if (st == "coverage") ui.view.coverage = true;
    else if (st == "populations") ui.view.fibrePopulations = true;
    else if (st == "wire") ui.view.wire = true;
    else if (st == "preview") ui.panels.previewOpen = true;  // the Render preview window
    else if (st == "appearance") revealYarnFibres(ui, Reveal::Appearance);  // the Fibres step, Fibre appearance open
    else if (st == "compare") ui.view.compare = true;  // Compare with source (suspended where it cannot be drawn)
    else if (st == "about") ui.panels.aboutOpen = true;  // the About window
    else if (st.rfind("about:", 0) == 0) {  // the About window on its tab n (0 Wool Shader ... 4 Blender and privacy)
        ui.panels.aboutOpen = true;
        ui.panels.aboutTabRequest = std::atoi(st.c_str() + 6);
    }
    else if (st.rfind("stage:", 0) == 0) ui.stage = stageNamed(st.substr(6));
    // The Mark step: a markings file, the isolated piece, the other pieces' look, the preview, a build.
    else if (st.rfind("mark-load:", 0) == 0) loadMarkingsFrom(ui, st.substr(10));
    else if (st.rfind("mark-piece:", 0) == 0) isolatePiece(ui, std::atoi(st.c_str() + 11));
    else if (st.rfind("mark-component:", 0) == 0) setKnitComponent(ui, std::atoi(st.c_str() + 15) - 1);   // numbered from 1, as the panel shows
    else if (st == "mark-around:hidden") setNeighbourLook(ui, NeighbourLook::Hidden);
    else if (st == "mark-around:shown") setNeighbourLook(ui, NeighbourLook::Shown);
    else if (st == "mark-view:rows") setMarkView(ui, MarkView::Rows);
    else if (st == "mark-view:stitch") setMarkView(ui, MarkView::StitchMesh);
    else if (st == "mark-view:markings") setMarkView(ui, MarkView::Markings);
    else if (st == "mark-build") buildMarkedPiece(ui);
    else if (st.rfind("mark-export-mesh:", 0) == 0) {  // mark-export-mesh:<path>: Export Stitch Mesh to that OBJ file
        if (!exportStitchMesh(ui, st.substr(17))) std::fprintf(stderr, "%s: %s\n", st.c_str(), ui.mark.note.c_str());
    }
    else if (st == "mark-stitch-source:on") ui.mark.stitchSource = true;
    else if (st == "mark-stitch-source:off") ui.mark.stitchSource = false;
    else if (st == "mark-highlight-findings:on") ui.mark.stitchHighlight = true;
    else if (st == "mark-highlight-findings:off") ui.mark.stitchHighlight = false;
    else if (st == "mark-knit-piece") knitThisPiece(ui);  // the footer's Knit this piece (guided)
    else if (st == "print-report") {  // the shown knit's report (the Report window's text) on stdout
        if (ui.results.knit) std::fputs(ui.results.knit->report.c_str(), stdout);
    }
    else if (st.rfind("export-dir:", 0) == 0) {  // export-dir:<folder>: the Export step's folder, as if chosen
        ui.exportTo.dir = st.substr(11);
        ui.exportTo.dirChosen = true;
    }
    else if (st.rfind("settings-save:", 0) == 0) {  // settings-save:<path>: the settings as a file that --settings reads back
        std::ofstream f(st.substr(14), std::ios::binary);
        f << knit::writeSettings(ui.s);
        if (!f) std::fprintf(stderr, "%s: not written\n", st.c_str());
    }
    else if (st.rfind("mark-add:", 0) == 0) {  // mark-add:<meaning>[:<shape>] (start, finish, row_start, boundary, cut, attachment; point, path, closed)
        const std::string rest = st.substr(9);
        const size_t colon = rest.find(':');
        knit::MarkingMeaning meaning = knit::MarkingMeaning::Start;
        knit::MarkingShape shape = knit::MarkingShape::Point;
        if (!knit::parseMeaning(rest.substr(0, colon), meaning)) std::fprintf(stderr, "%s: unknown meaning\n", st.c_str());
        else {
            addMarking(ui, meaning);
            if (colon != std::string::npos && knit::parseShape(rest.substr(colon + 1), shape)) setMarkingShape(ui, ui.mark.activeMarking, shape);
        }
    } else if (st.rfind("mark-point:", 0) == 0) {  // mark-point:<x>/<y>/<z> (file coordinates): a click on the isolated piece there
        double x = 0, y = 0, z = 0;
        if (std::sscanf(st.c_str() + 11, "%lf/%lf/%lf", &x, &y, &z) != 3 || !placeMarkingPoint(ui, knit::Vec3{x, y, z}))
            std::fprintf(stderr, "%s: not placed (an active marking and an isolated piece are needed)\n", st.c_str());
    }
    else if (st.rfind("mark-accept:", 0) == 0) acceptProposal(ui, std::atoi(st.c_str() + 12));
    else if (st.rfind("mark-reject:", 0) == 0) rejectProposal(ui, std::atoi(st.c_str() + 12));
    else if (st.rfind("mark-save:", 0) == 0) saveMarkingsAs(ui, st.substr(10));
    else if (st.rfind("mark-preparation:", 0) == 0) setMarkingPreparation(ui, st.substr(17));   // facets | loop_limit
    else if (st.rfind("mark-move:", 0) == 0) {  // mark-move:<id>:<point>:<x>/<y>/<z> (file coordinates; '/' because ',' separates steps)
        int id = 0, point = 0;
        double x = 0, y = 0, z = 0;
        if (std::sscanf(st.c_str() + 10, "%d:%d:%lf/%lf/%lf", &id, &point, &x, &y, &z) != 5 || !moveMarkingPoint(ui, id, point, knit::Vec3{x, y, z}))
            std::fprintf(stderr, "%s: no such marking, point or piece\n", st.c_str());
    } else if (st.rfind("mark-report:", 0) == 0) {
        std::FILE* f = std::fopen(st.substr(12).c_str(), "wb");
        if (f) {
            const std::string text = markReviewText(ui);
            std::fwrite(text.data(), 1, text.size(), f);
            std::fclose(f);
        }
    }
    // Assisted stockinette: suggestions (cards numbered from 1), previews, the quick check, the latency record, bursts.
    else if (st == "mark-propose") suggestMarkings(ui);
    else if (st.rfind("mark-refine:", 0) == 0) refineSuggestion(ui, std::atoi(st.c_str() + 12));
    else if (st == "mark-cancel") cancelSuggestionWork(ui);
    else if (st == "mark-preview:off") previewSuggestion(ui, 0);
    else if (st.rfind("mark-preview:", 0) == 0) previewSuggestion(ui, std::atoi(st.c_str() + 13));
    else if (st.rfind("mark-accept-setup:", 0) == 0) {  // mark-accept-setup:<card>[:replace]
        if (!acceptSuggestion(ui, std::atoi(st.c_str() + 18), st.find(":replace") != std::string::npos))
            std::fprintf(stderr, "%s: not accepted (outdated, or your own markings would be replaced: add :replace)\n", st.c_str());
    }
    else if (st.rfind("mark-reject-setup:", 0) == 0) rejectSuggestion(ui, std::atoi(st.c_str() + 18));
    else if (st.rfind("mark-quick-report:", 0) == 0) writeTextFile(st.substr(18), quickReviewText(ui));
    else if (st.rfind("mark-inventory:", 0) == 0) writeTextFile(st.substr(15), markInventoryText(ui));   // the pieces as read, and the size
    else if (st == "mark-latency-reset") {  // a new latency record (each background condition its own)
        ui.mark.latencyMs.clear();
        ui.mark.latencyBackground.clear();
        ui.mark.coalesced = 0;
        ui.mark.coalescedBackground = 0;
        ui.mark.pendingStamps.clear();
    }
    else if (st.rfind("mark-latency:", 0) == 0) writeTextFile(st.substr(13), latencyText(ui));
    else if (st == "mark-refresh:on") setLayoutRefresh(ui, true);
    else if (st == "mark-refresh:off") setLayoutRefresh(ui, false);
    else if (st == "mark-relax:on") ui.mark.relax = true;
    else if (st == "mark-relax:off") ui.mark.relax = false;
    else if (st == "mark-highlight:problem") {  // the first problem with a place: the quick check's if current, else the layout's
        const std::shared_ptr<const knit::MarkingSession::QuickView> q = currentQuick(ui);
        const knit::MarkedPieceView* b = shownBuild(ui);
        const std::vector<knit::MarkingFinding>* f = q ? &q->review.findings : b ? &b->review.findings : nullptr;
        for (size_t i = 0; f && i < f->size(); ++i)
            if ((*f)[i].level == knit::FindingLevel::Problem && ((*f)[i].firstTracePoint >= 0 || !(*f)[i].at.empty())) {
                highlightFinding(ui, static_cast<int>(i), q != nullptr);
                break;
            }
    }
    else if (st.rfind("mark-highlight:", 0) == 0) {  // mark-highlight:<finding index>[:quick] (the shown layout's findings, or the quick check's)
        highlightFinding(ui, std::atoi(st.c_str() + 15), st.find(":quick") != std::string::npos);
    }
    else if (st.rfind("mark-burst:", 0) == 0) {  // mark-burst:<id>:<point>:<x>/<y>/<z>;<x>/<y>/<z>;... one move per frame (a drag)
        int id = 0, point = 0, used = 0;
        if (std::sscanf(st.c_str() + 11, "%d:%d:%n", &id, &point, &used) < 2) std::fprintf(stderr, "%s: mark-burst:<id>:<point>:<positions>\n", st.c_str());
        std::vector<knit::Vec3> positions;
        std::string rest = st.substr(11 + used);
        while (!rest.empty()) {
            const size_t semi = rest.find(';');
            double x = 0, y = 0, z = 0;
            if (std::sscanf(rest.substr(0, semi).c_str(), "%lf/%lf/%lf", &x, &y, &z) == 3) positions.push_back(knit::Vec3{x, y, z});
            rest = semi == std::string::npos ? std::string() : rest.substr(semi + 1);
        }
        queueBurst(ui, id, point, positions);
    }
    else if (st.rfind("set:", 0) == 0 && st.find('=') != std::string::npos) {  // set:<settings key>=<value>
        const size_t eq = st.find('=');
        std::string err;
        if (!knit::setSetting(ui.s, st.substr(4, eq - 4), st.substr(eq + 1), &err)) std::fprintf(stderr, "%s: %s\n", st.c_str(), err.c_str());
    }
}

// With the Render preview open, a shot waits until it shows a settled image (16 samples) or said why not.
static bool renderPreviewSettled(const SceneUi& ui) {
    if (!ui.panels.previewOpen) return true;
    const knit::RenderPreviewState& p = ui.renderPreview.state;
    if (p.phase == "refused" || p.phase == "error") return true;
    if (!ui.renderPreview.hint.empty() && p.frameNumber == 0 && (p.phase == "idle" || p.phase == "off")) return true;  // nothing to render
    return (p.phase == "rendering" || p.phase == "converged") && p.spp >= 16 && ui.renderPreview.uploadedFrame > 0;
}

ScriptedRun::ScriptedRun(const ViewerOptions& options) : options_(options) {
    for (size_t p = 0; p < options_.steps.size();) {
        const size_t q = options_.steps.find(',', p);
        todo_.push_back(options_.steps.substr(p, q == std::string::npos ? std::string::npos : q - p));
        p = q == std::string::npos ? options_.steps.size() : q + 1;
    }
    importConfirmed_ = options_.importFile.empty() || options_.units.empty();
}

bool ScriptedRun::beforeFrame(SceneUi& ui, double now) {
    using Phase = knit::ImportSession::Phase;
    const bool scriptedImport = !options_.importFile.empty() && !options_.units.empty();
    knit::ImportSession::Snapshot file = ui.import.session.snapshot();
    if (!importConfirmed_ && file.phase == Phase::Analysed) {
        confirmImport(ui, options_.importChoices);
        importConfirmed_ = true;
        file = ui.import.session.snapshot();
    }
    if (scriptedImport && file.phase == Phase::Failed && options_.keepFileIfRebuildFails && canKeepFileAfterFailedRebuild(ui)) {
        std::fprintf(stderr, "rebuild failed (%s): pressing \"Keep the file's own surface\"\n", file.error.c_str());
        keepFileAfterFailedRebuild(ui);
        file = ui.import.session.snapshot();
    }
    if (scriptedImport && file.phase == Phase::Failed) {
        std::fprintf(stderr, "import failed: %s\n", file.error.c_str());
        return false;
    }
    const knit::MarkingSession::Snapshot markSnap = ui.mark.session.snapshot();
    const knit::MarkingSession::Phase marking = markSnap.phase;
    const bool deepBusy = marking == knit::MarkingSession::Phase::Reading || marking == knit::MarkingSession::Phase::Building ||
                          marking == knit::MarkingSession::Phase::Proposing || marking == knit::MarkingSession::Phase::Refining;
    const bool markingBusy = deepBusy || markSnap.quickBusy || markSnap.quickPreparing || !ui.mark.burst.empty() || ui.mark.layoutDueAt > 0;
    // a burst of point moves measures the quick check under background work: it starts as soon as a build, suggestion,
    // refinement or knit runs (or at idle when nothing runs); each edit records whether background work ran when it was
    // made, and mark-latency reports those edits apart (a suggestion on the arm takes a fraction of a second)
    if (deepBusy || ui.job.active) {
        if (backgroundSince_ == 0) backgroundSince_ = now;
    } else {
        backgroundSince_ = 0;
    }
    if (!todo_.empty() && todo_.front().rfind("mark-burst:", 0) == 0 && backgroundSince_ > 0 && ui.mark.burst.empty()) {
        std::fprintf(stderr, "mark-burst under background work: Mark step's deep work %s, job %s (%s), %.1f s into it\n", deepBusy ? "running" : "idle",
                     ui.job.active ? "running" : "idle", ui.job.active ? "knit" : "-", now - backgroundSince_);
        applyScriptStep(ui, todo_.front());
        todo_.erase(todo_.begin());
    }
    idle_ = !ui.job.active && !ui.builder->snapshot().busy && file.phase != Phase::Reading && file.phase != Phase::Preparing && importConfirmed_ &&
            !ui.import.revalidating && !ui.previews.fibres.busy() && !ui.previews.discs.busy() && !markingBusy && renderPreviewSettled(ui);
    // "running": the shot is taken while the job runs; "cancel": cancels the running job; both after
    // it ran 3 s (right after a job step).
    if (ui.job.active && runningSince_ == 0) runningSince_ = now;
    if (!ui.job.active) runningSince_ = 0;
    if (!todo_.empty() && (todo_.front() == "cancel" || todo_.front() == "running") && ui.job.active && now - runningSince_ > 3 && ui.mark.burst.empty()) {
        if (todo_.front() == "cancel") cancelJob(ui);
        else shotNow_ = true;
        todo_.erase(todo_.begin());
    }
    if (idle_ && !todo_.empty() && (!ui.import.fileOpen || options_.importFile.empty())) {
        const std::string step = todo_.front();
        if (step.rfind("mark-burst:", 0) == 0) std::fprintf(stderr, "mark-burst at idle\n");
        todo_.erase(todo_.begin());
        applyScriptStep(ui, step);
    }
    return true;
}

void ScriptedRun::frameShot(Camera& camera, const Bounds& content) const {
    if (options_.shotPath.empty() || options_.zoom <= 0) return;
    frameBox(camera, content);
    const Vec3 size = content.hi - content.lo;
    const double* f = options_.targetFraction;
    camera.target = {content.lo.x + size.x * f[0], content.lo.y + size.y * f[1], content.lo.z + size.z * f[2]};
    camera.dist /= options_.zoom;
}

bool ScriptedRun::shotDue() {
    if (options_.shotPath.empty()) return false;
    quietFrames_ = idle_ && todo_.empty() ? quietFrames_ + 1 : 0;
    return quietFrames_ > 6 || shotNow_;
}

int ScriptedRun::reportShot(const SceneUi& ui, bool written) const {
    const Results& r = ui.results;
    const bool failed = r.knitAttempt.outcome == knit::JobOutcome::Failed || r.fibreAttempt.outcome == knit::JobOutcome::Failed ||
                        r.exportAttempt.outcome == knit::JobOutcome::Failed;
    for (const Attempt* a : {&r.knitAttempt, &r.fibreAttempt, &r.exportAttempt})
        if (!a->message.empty()) std::printf("  %s\n", a->message.c_str());
    if (r.knit && r.knit->knit && !r.knit->knit->frozenKnitDir.empty()) std::printf("knit: loaded from %s\n", r.knit->knit->frozenKnitDir.c_str());
    if (r.exported) std::printf("  exported to %s\n", r.exported->exportedTo.c_str());
    if (r.fibres && r.fibres->fibreResult) std::printf("fibre log (%.1f s):\n%s", r.fibres->ms / 1000, r.fibres->fibreResult->log.c_str());
    return written && !failed ? 0 : 1;
}

// --- the command line's object

// The settings: this machine's compute defaults, then a settings file. False: refused, said on stderr.
bool loadSettings(SceneUi& ui, const ViewerOptions& o) {
    const knit::SceneEntry* e = knit::findScene(o.scene.empty() ? knit::SceneSettings{}.scene : o.scene);
    if (!e) {
        std::fprintf(stderr, "unknown scene %s\n", o.scene.c_str());
        return false;
    }
    // Without a scene or a settings file the viewer starts with no object: plain defaults, not the default
    // scene's listed settings (they belong to that example).
    const bool noObject = o.scene.empty() && o.settingsFile.empty();
    ui.s = knit::viewerDefaults(noObject ? knit::SceneEntry{} : *e);
    if (noObject) startWithoutObject(ui);
    // A scene named on the command line is chosen as the command-line tool chooses it: the settings' workflow
    // (knit unless a settings file says fur), so scripts run as before. Only the interactive open asks.
    // With the Fur switch off Knit is the only workflow: nothing is ever asked.
    ui.workflowChosen = !knit::furEnabled() || !o.scene.empty() || !o.settingsFile.empty();
    if (o.settingsFile.empty()) return true;
    std::ifstream f(o.settingsFile);
    std::stringstream ss;
    ss << f.rdbuf();
    std::string err;
    std::vector<std::string> retired;
    if (!knit::readSettings(ss.str(), ui.s, &err, &retired)) {
        std::fprintf(stderr, "%s: %s\n", o.settingsFile.c_str(), err.c_str());
        return false;
    }
    for (const std::string& k : retired) std::fprintf(stderr, "%s: %s is no longer used (ignored)\n", o.settingsFile.c_str(), k.c_str());
    return true;
}

// Opens what the command line names, or what the settings name.
void openObject(SceneUi& ui, const ViewerOptions& o) {
    adoptSettings(ui);
    if (!o.importFile.empty()) {
        const bool chosen = ui.workflowChosen;  // the command line's settings carry the workflow (loadSettings)
        openMeshFile(ui, o.importFile);
        ui.workflowChosen = chosen;
        // Scripted: the options are given, no dialog decision needed.
        if (!o.units.empty()) presetImportUnits(ui, parseUnits(o.units), o.up == "z" ? knit::UpAxis::Z : knit::UpAxis::Y);
    } else if (!ui.s.scene.empty() && ui.s.scene != "import") {
        if (knit::furWorkflow(ui.s)) {  // Fur: nothing to knit; it opens on the disc preview (--do fibres grows the fur)
            ui.stage = Stage::Fibres;
            ui.view.fibres = FibreView::Discs;
        } else {
            startKnit(ui);  // the scene opens from the cache when it holds these settings
        }
    }
    if (!o.stage.empty()) ui.stage = stageNamed(o.stage);
}

} // namespace viewer
