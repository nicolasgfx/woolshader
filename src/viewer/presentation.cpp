// What reads SceneUi (workflow.h): what each step's state is, what a running job says, the export's
// texts, and what the viewport shows with its captions. Every function here takes const SceneUi&,
// so it can be tested without a window, and the panels only draw what it says. No ImGui and no GL.
#include "model/coordinates.h"
#include "jobs/results.h"
#include "settings/scene_registry.h"
#include "settings/settings_profiles.h"
#include "viewer/workflow.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <string>

#include "export/export.h"
#include "settings/file_presets.h"
#include "settings/settings_mapping.h"
#include "settings/thickness_presets.h"

namespace fs = std::filesystem;
using knit::Vec3;

namespace viewer {

// --- what is current

const knit::SceneEntry* entryOf(const SceneUi& ui) { return knit::sceneEntryFor(ui.s); }

const knit::KnitResult* knitResult(const SceneUi& ui) { return ui.results.knit ? ui.results.knit->knit.get() : nullptr; }

bool furShown(const SceneUi& ui) { return ui.workflowChosen && knit::furWorkflow(ui.s); }

const knit::FibreResult* fibreResult(const SceneUi& ui) {
    if (!ui.results.fibres || !ui.results.fibres->fibreResult) return nullptr;
    const knit::FibreResult* f = ui.results.fibres->fibreResult.get();
    if (f->workflow != knit::workflowOf(ui.s)) return nullptr;
    if (f->workflow == knit::Workflow::Fur) return f;  // fur needs no knit
    return knitResult(ui) ? f : nullptr;
}

bool knitCurrent(const SceneUi& ui) {
    const knit::KnitResult* k = knitResult(ui);
    return k && k->computational && !knit::knitOutdated(*k, ui.s);
}

// The Fibres (Fur) step's input is ready: wool needs the current knit; fur only the object and the choice.
bool fibresInputReady(const SceneUi& ui) {
    if (!entryOf(ui) || !ui.workflowChosen) return false;
    return knit::furWorkflow(ui.s) || knitCurrent(ui);
}

bool runningStage(const SceneUi& ui, knit::BuiltScene::Kind k) { return ui.job.active && ui.job.kind == k; }

bool sameTube(const knit::SceneSettings& a, const knit::SceneSettings& b) {
    return a.tubeSides == b.tubeSides && a.tubeRing == b.tubeRing && a.plyTubeSides == b.plyTubeSides && a.plyTubeRing == b.plyTubeRing;
}

// The shown knit's yarn tubes differ from the tessellation in the settings and could be rebuilt now (no
// knitting): the knit is current and the tessellation valid. The Knit step's Yarn tessellation and the
// Export step's detail levels both rebuild on it.
bool shownTubesDiffer(const SceneUi& ui) {
    const bool differ = ui.results.knit && !ui.results.partial && !sameTube(ui.results.knit->settings, ui.s);
    return differ && knitResult(ui) && knitCurrent(ui) && knit::tubeSettingsValid(ui.s, nullptr);
}

// The Knit footer below the settings (the stage's one place for its job): shown when the step knits.
bool knitFooterShown(const SceneUi& ui) {
    return ui.workflowChosen && knit::knitWorkflow(ui.s) && (!ui.import.lacksCharts || knit::knitsGuidedPiece(ui.s));
}

std::string waitReason(const SceneUi& ui) {
    if (!ui.job.active) return "";
    switch (ui.job.kind) {
    case knit::BuiltScene::Kind::Knit: return "Knitting is running: wait for it or cancel it (status bar).";
    case knit::BuiltScene::Kind::Fibres:
        return knit::furWorkflow(ui.s) ? "The fur is growing: wait or cancel it (status bar)." : "Fibres are being generated: wait or cancel them (status bar).";
    case knit::BuiltScene::Kind::Export: return "An export is running: wait for it or cancel it (status bar).";
    }
    return "";
}

// --- the stages' words

bool stageInStrip(const SceneUi& ui, Stage s) { return !(furShown(ui) && (s == Stage::Mark || s == Stage::Knit)); }

const char* stageLabel(const SceneUi& ui, Stage s) {
    const bool fur = furShown(ui);
    switch (s) {
    case Stage::Import: return "Import";
    case Stage::Mark: return "Mark";
    case Stage::Knit: return "1  Knit";
    case Stage::Fibres: return fur ? "1  Fur" : "2  Fibres";
    case Stage::Export: return fur ? "2  Export" : "3  Export";
    }
    return "";
}

StageState importState(const SceneUi& ui) {
    StageState r;
    const knit::ImportSession::Snapshot file = ui.import.session.snapshot();
    if (file.phase == knit::ImportSession::Phase::Reading || file.phase == knit::ImportSession::Phase::Preparing) return {StageStatus::Running, "Working", file.stage};
    if (ui.import.fileOpen) return {StageStatus::NeedsUpdate, "Confirm", "a file is open in the Import step"};
    if (entryOf(ui) && !ui.workflowChosen) return {StageStatus::NeedsUpdate, "Choose", "choose Knit or Fur for the object"};
    if (const knit::SceneEntry* e = entryOf(ui)) {
        r.status = StageStatus::Ready;
        r.word = "Ready";
        r.detail = e->label;
        return r;
    }
    if (ui.s.scene == "import") return {StageStatus::Failed, "Missing", "the prepared import is missing: import the file again"};
    return {StageStatus::NotStarted, "No object", ""};
}

// The Mark step is optional: it tells the stitch mesher where to start, finish and begin each row (marking_flow.h).
StageState markState(const SceneUi& ui) {
    if (ui.s.scene != "import") return {StageStatus::Optional, "Imports only", "surface markings are drawn on an imported file's pieces"};
    const knit::MarkingSession::Snapshot snap = ui.mark.session.snapshot();
    switch (snap.phase) {
    case knit::MarkingSession::Phase::Reading: return {StageStatus::Running, "Reading", ""};
    case knit::MarkingSession::Phase::Building: return {StageStatus::Running, "Building", ""};
    case knit::MarkingSession::Phase::Proposing: return {StageStatus::Running, "Suggesting", ""};
    case knit::MarkingSession::Phase::Refining: return {StageStatus::Running, "Refining", ""};
    default: break;
    }
    if (ui.mark.stale) return {StageStatus::NeedsUpdate, "Stale", "the markings were made on another mesh"};
    const size_t n = ui.mark.set.markings.size();
    if (n == 0) return {StageStatus::Optional, "Optional", "markings for the stitch mesher: start, finish, row-start path, ..."};
    return {StageStatus::Ready, std::to_string(n) + (n == 1 ? " marking" : " markings"), ""};
}

StageState knitState(const SceneUi& ui) {
    if (!ui.workflowChosen) return {StageStatus::Blocked, kChooseWorkflow, "the Import step asks for Knit or Fur"};
    if (knit::furWorkflow(ui.s)) return {StageStatus::Optional, "Not used", "the Fur workflow grows fur on the mesh: nothing is knitted"};
    if (runningStage(ui, knit::BuiltScene::Kind::Knit)) {
        if (ui.job.cancelRequested) return {StageStatus::Cancelling, "Stopping", ui.job.stoppingNote};
        return {StageStatus::Running, ui.job.bakeMaps ? "Baking" : ui.job.meshOnly ? "Meshing" : "Knitting", ""};
    }
    const knit::KnitResult* k = knitResult(ui);
    if (!k) {
        if (ui.results.knitAttempt.outcome == knit::JobOutcome::Failed) return {StageStatus::Failed, "Failed", ui.results.knitAttempt.message};
        if (ui.results.knitAttempt.outcome == knit::JobOutcome::Cancelled) return {StageStatus::Cancelled, "Cancelled", ui.results.knitAttempt.message};
        return {StageStatus::NotStarted, "Not knitted", ""};
    }
    if (!k->computational) return {StageStatus::ViewOnly, "View only", k->viewOnlyReason};
    if (knit::knitOutdated(*k, ui.s)) return {StageStatus::NeedsUpdate, "Outdated", "knitting settings changed since this result: update the knit"};
    if (k->final) return {StageStatus::Ready, "Final", ""};
    // A Draft (or custom) knit is a result to look at, not to export: said where it is shown.
    return {StageStatus::Diagnostic, knit::isPreviewProfile(k->settings) ? "Draft" : "Custom", "not exportable: the export needs a Final knit"};
}

// The generated fibres / fur themselves (the Fibres or Fur footer), whatever the export includes:
// running, not possible here, outdated (knit or settings), current, failed, cancelled, waiting for
// the knit, not generated.
StageState fibreOutputState(const SceneUi& ui) {
    if (runningStage(ui, knit::BuiltScene::Kind::Fibres)) return ui.job.cancelRequested ? StageState{StageStatus::Cancelling, "Stopping", ui.job.stoppingNote}
                                                                                       : StageState{StageStatus::Running, furShown(ui) ? "Growing" : "Generating", ""};
    if (!ui.workflowChosen) return {StageStatus::Blocked, kChooseWorkflow, "the Import step asks for Knit or Fur"};
    if (knit::furWorkflow(ui.s)) {
        const knit::FibreResult* f = fibreResult(ui);
        if (f && knit::fibresOutdated(*f, ui.s)) return {StageStatus::NeedsUpdate, "Outdated", "fur settings changed since: grow it again"};
        if (f) return {StageStatus::Ready, "Current", ""};
        if (ui.results.fibreAttempt.outcome == knit::JobOutcome::Failed) return {StageStatus::Failed, "Failed", ui.results.fibreAttempt.message};
        if (ui.results.fibreAttempt.outcome == knit::JobOutcome::Cancelled) return {StageStatus::Cancelled, "Cancelled", ""};
        return {StageStatus::NotStarted, "Not grown", ""};
    }
    const knit::KnitResult* k = knitResult(ui);
    const knit::FibreResult* f = fibreResult(ui);
    if (f && k && f->knitId != k->id) return {StageStatus::NeedsUpdate, "Outdated", "the knit changed since these fibres: generate them again"};
    if (f && knit::fibresOutdated(*f, ui.s)) return {StageStatus::NeedsUpdate, "Outdated", "fibre settings changed since: generate them again"};
    if (f) return {StageStatus::Ready, "Current", ""};
    if (ui.results.fibreAttempt.outcome == knit::JobOutcome::Failed) return {StageStatus::Failed, "Failed", ui.results.fibreAttempt.message};
    if (ui.results.fibreAttempt.outcome == knit::JobOutcome::Cancelled) return {StageStatus::Cancelled, "Cancelled", ""};
    if (!knitCurrent(ui)) return {StageStatus::Blocked, "Knit first", ""};
    return {StageStatus::NotStarted, "Not generated", ""};
}

// The stage bar's Fibres word: the output's state, except that wool fibres the export does not
// include read "Off" (generated or not), or "Yarn fibres" when the yarn itself is made of fibres (they
// need no job: the export and the render preview build them).
StageState fibreState(const SceneUi& ui) {
    StageState o = fibreOutputState(ui);
    if (!ui.workflowChosen || knit::furWorkflow(ui.s)) {
        if (o.status == StageStatus::Ready) o.word = "Ready";
        return o;
    }
    const bool noFuzz = !ui.s.fibres && (o.status == StageStatus::Ready || o.status == StageStatus::NotStarted ||
                                         (o.status == StageStatus::Blocked && o.word == "Knit first"));
    if (noFuzz && ui.s.yarnBody.enabled) return {StageStatus::Optional, "Yarn fibres", "the yarn's own fibres; no fuzz in the export"};
    if (o.status == StageStatus::Ready) return ui.s.fibres ? StageState{StageStatus::Ready, "Ready", ""} : StageState{StageStatus::Optional, "Off", "fibres generated; not in the export"};
    if (noFuzz) return {StageStatus::Optional, "Off", "not in the export"};
    return o;
}

knit::ExportEligibility eligibility(const SceneUi& ui) {
    const knit::ExportFormat* f = knit::findExportFormat(ui.exportTo.format);
    if (!ui.workflowChosen) {
        knit::ExportEligibility e;
        e.reasons.push_back("No workflow is chosen for the object yet.");
        e.fix = kChooseWorkflow;
        return e;
    }
    const bool fur = knit::furWorkflow(ui.s);
    const bool fibres = (ui.s.fibres || fur) && f && f->fibres;
    // Fur: the result carries the base meshes, so the export always needs it.
    return knit::exportEligibility(ui.s, knitResult(ui), fibres || fur ? fibreResult(ui) : nullptr,
                                   fibres ? knit::FibreChoice::Included : knit::FibreChoice::LeftOut);
}

StageState exportState(const SceneUi& ui) {
    if (runningStage(ui, knit::BuiltScene::Kind::Export)) return ui.job.cancelRequested ? StageState{StageStatus::Cancelling, "Stopping", ui.job.stoppingNote}
                                                                                       : StageState{StageStatus::Running, "Exporting", ""};
    const knit::ExportEligibility el = eligibility(ui);
    if (el.state == knit::ExportEligibility::State::Blocked) {
        // The stage bar's word names what is missing (the panel lists every reason).
        const std::string word = el.fix == "Knit at Final quality" ? "Needs Final knit"
                                 : el.fix == "Generate fibres" || el.fix == "Update fibres" ? "Needs fibres"
                                 : el.fix == "Grow fur" || el.fix == "Update fur"     ? "Needs fur"
                                 : el.fix == kChooseWorkflow                          ? kChooseWorkflow
                                 : el.fix == "Correct the yarn's fibres"              ? "Fibres refused"
                                                                                      : "Not ready";
        return {StageStatus::Blocked, word, el.reasons.empty() ? "" : el.reasons.front()};
    }
    const bool same = furShown(ui) ? ui.results.exported && ui.results.exported->fibreResult.get() == fibreResult(ui)
                                     : ui.results.exported && ui.results.exported->knit.get() == knitResult(ui);
    if (same && ui.results.exportAttempt.outcome == knit::JobOutcome::Done) return {StageStatus::Ready, "Exported", ui.results.exported->exportedTo};
    if (el.state == knit::ExportEligibility::State::Diagnostic) return {StageStatus::Diagnostic, "Diagnostic", "parts failed validation"};
    return {StageStatus::Ready, "Ready", ""};
}

// --- the Knit step's readings

PreparedSurface preparedSurfaceCheck(const knit::SceneSettings& s) {
    PreparedSurface p;
    p.applies = s.scene == "import" && s.importSurface == "wrapped" && s.importWrapStitchMm > 0;
    if (!p.applies) return p;
    p.preparedMm = s.importWrapStitchMm;
    // The rebuilt surface's detail follows the stitch it was prepared for: much wider stitches span
    // features the surface kept, and their yarns may tangle there.
    p.tooWide = s.wale > 1.3 * s.importWrapStitchMm;
    return p;
}

namespace {

// The Final profile's steps as they run for these settings (a Final copy): settle, repair on spiral
// pieces, one diameter or the fallback list, and the surface check.
std::string finalStepsText(const knit::SceneSettings& s, const knit::SceneEntry& e) {
    std::string t = std::string("This knit: settle, ") + (knit::knitsSewnPieces(e, s) ? "repair, " : "") +
                    (s.diameterAttempts == 1 ? "one diameter" : "diameter fallback") + " and the surface check run automatically.";
    if (!s.diameterFallback.empty()) {
        t += "\nDiameter fallback:";
        for (double d : s.diameterFallback) {
            char b[16];
            std::snprintf(b, sizeof b, " %.2f", d);
            t += b;
        }
        t += " mm";
        if (s.diameterAttempts > 0) t += " (" + std::to_string(s.diameterAttempts) + " attempt(s) with the explicit solver)";
        t += ".";
    }
    return t;
}

} // namespace

KnitActionView knitActionView(const SceneUi& ui, KnitQuality quality, const knit::ComputeCaps& caps) {
    KnitActionView v;
    const knit::SceneEntry* e = entryOf(ui);
    if (!e) return v;
    knit::SceneSettings copy = ui.s;
    switch (quality) {
    case KnitQuality::Draft:
        knit::applyPreviewProfile(copy);
        v.tooltip = "Knits each part once and checks it: seconds to minutes. The checks are provisional: parts that would pass "
                    "after settling or a smaller yarn show FAIL here. Not exportable.";
        break;
    case KnitQuality::Final:
        knit::applyFullProfile(copy, *e);
        v.tooltip = "Adds settling (relaxing the yarn around defects), repair on the teddy's spiral pieces, a smaller yarn where a "
                    "part still fails (when the solver allows several tries) and the rendered-surface check. Minutes to hours on "
                    "large objects. Required for export.\n\n" +
                    finalStepsText(copy, *e);
        break;
    }
    v.notes = knit::computeNotes(copy, caps);
    return v;
}

bool fibresNeedNewTubes(const SceneUi& ui) {
    const knit::FibreResult* f = fibreResult(ui);
    return f && !furShown(ui) && !sameTube(f->settings, ui.s);
}

const char* patternLabel(knit::StitchPattern p) {
    switch (p) {
    case knit::StitchPattern::Stockinette: return "Stockinette";
    case knit::StitchPattern::ReverseStockinette: return "Reverse stockinette";
    case knit::StitchPattern::Garter: return "Garter";
    case knit::StitchPattern::Rib1x1: return "Rib 1 x 1";
    case knit::StitchPattern::Rib2x2: return "Rib 2 x 2";
    }
    return "Stockinette";
}

std::string gaugeChoiceLabel(const knit::SceneSettings& s, const knit::SceneEntry& e) {
    const std::string current = knit::currentGaugeChoice(s, e);
    for (const knit::GaugeChoice& c : knit::gaugeChoices(e))
        if (c.id == current) return c.label;
    return "Custom";
}

std::string yarnGeometrySummary(const knit::SceneSettings& s) {
    char b[96];
    std::snprintf(b, sizeof b, "%.2f D%s%.2f mm clearance", s.depthAmplitude, kDot, s.clearance);
    std::string t = b;
    if (s.flatten && s.plies <= 1) t += std::string(kDot) + "contact flattening";
    return t;
}

std::string characterSummary(const knit::SceneSettings& s) {
    const knit::HandKnit& h = s.hand;
    if (h.unevenness == 0 && h.columnShift == 0 && h.loopRate == 0 && h.spotRate == 0 && h.waveHeight == 0) return "Every stitch alike";
    char b[160];
    std::snprintf(b, sizeof b, "%.0f %% uneven%s%.2f raised loops / cm2%s%.2f spots / cm2", h.unevenness, kDot, h.loopRate, kDot, h.spotRate);
    return b;
}

std::string stitchLayoutSummary(const knit::SceneSettings& s, const knit::SceneEntry& e) {
    if (knit::knitsSewnPieces(e, s)) return "Sewn pieces";
    char b[32];
    std::snprintf(b, sizeof b, "%.0f deg", s.angleDeg);
    return std::string(s.coverage == "overlap" ? "Overlapped seams" : "Open seams") + kDot + b;
}

// --- the flat swatch

namespace {

// A flat swatch that did not pass, in words: what failed, or what the check could not decide.
std::string swatchOutcome(const knit::PatternReference& r) {
    if (r.verdict == knit::Status::Fail) return "fails (" + knit::failCategoryText(r.fail) + ")";
    return "is unresolved (" + knit::failCategoryText(r.unresolved) + ")";
}

// The flat swatch's verdict on a pattern (null while it is computed).
const knit::PatternReference* swatchOf(const std::vector<knit::PatternReference>* refs, knit::StitchPattern p) {
    if (!refs) return nullptr;
    for (const knit::PatternReference& r : *refs)
        if (r.pattern == p) return &r;
    return nullptr;
}

} // namespace

FlatSwatchView flatSwatchOf(const SceneUi& ui, const std::vector<knit::PatternReference>* refs, knit::StitchPattern p) {
    const knit::SceneSettings& s = ui.s;
    FlatSwatchView v;
    if (!(s.wale > 0 && s.course > 0)) {
        v.status = StageStatus::Blocked;
        v.text = "Flat swatch not checked: the stitch width and row height must be positive.";
        return v;
    }
    const knit::PatternReference* r = swatchOf(refs, p);
    if (!r) {
        v.status = StageStatus::Running;
        v.text = "Checking the pattern on a flat swatch ...";
        return v;
    }
    if (r->verdict == knit::Status::Pass) {
        v.status = StageStatus::Ready;
        v.text = "Flat swatch fits at these values.";
        return v;
    }
    // The teddy's calibrated gauge fails flat by design: its curved pieces are validated by their knit.
    const knit::SceneEntry* e = entryOf(ui);
    v.status = StageStatus::Optional;
    v.outcome = swatchOutcome(*r);
    v.text = "Flat swatch " + v.outcome + (e && knit::knitsSewnPieces(*e, s) ? "; the teddy's pieces are checked by their knit." : "; the parts are checked after knitting.");
    return v;
}

std::string meshDetailSummary(const knit::SceneSettings& s, bool customShown) {
    const knit::FilePreset* level = knit::matchingFilePreset(knit::PresetGroup::StrandMesh, s);
    if (level && !customShown) return level->label;
    char b[64];
    std::snprintf(b, sizeof b, "Custom (%d sides, %.2f mm)", s.tubeSides, s.tubeRing);
    return b;
}

// --- the yarn made of fibres

// Its fibre segments on the knit's yarn or, rendered as maps, its map texels and tube triangles, beside the
// fuzz's exported segments, and the renderer's GPU memory for all of it (with the knit's tubes gone).
YarnFibresCost yarnFibresCost(const SceneUi& ui) {
    YarnFibresCost c;
    const knit::SceneSettings& s = ui.s;
    c.yarnMm = ui.results.knit ? ui.results.knit->yarnLength : 0;
    c.knitted = c.yarnMm > 0;
    c.maps = knit::yarnBodyRenderOf(s) == knit::YarnBodyRender::Maps;
    const knit::FibreResult* fr = fibreResult(ui);
    long long fuzz = 0;
    if (s.fibres && fr)
        for (const knit::ObjectSegments& o : fr->objectSegments) fuzz += o.exported;
    c.fuzzSegments = double(fuzz);
    char w[64];
    if (!s.fibres) std::snprintf(w, sizeof w, "not exported");
    else if (!fr) std::snprintf(w, sizeof w, "not generated yet");
    else std::snprintf(w, sizeof w, "%.1f M segments", double(fuzz) / 1e6);
    c.fuzzWord = w;
    if (!c.knitted) return c;
    const double base = knit::rendererMemoryGiB(0, 0, 0);
    if (c.maps) {
        const long long texels = knit::yarnMapTexels(s, c.yarnMm), triangles = knit::yarnMapTriangles(c.yarnMm);
        c.texels = double(texels);
        c.triangles = double(triangles);
        c.gib = knit::rendererMemoryGiB(fuzz, triangles, texels);
        c.texturesGib = knit::rendererMemoryGiB(0, 0, texels) - base;
        c.yarnDominates = c.texturesGib >= knit::rendererMemoryGiB(fuzz, triangles, 0) - base;
    } else {
        const long long body = knit::yarnBodySegments(s, c.yarnMm);
        c.yarnSegments = double(body);
        c.gib = knit::rendererMemoryGiB(body + fuzz, 0);
        c.yarnDominates = body >= fuzz;
    }
    return c;
}

// The Knit step's one line under the Construction tiles: the look, the plies and fibres, the yarn's cost and
// the fuzz's beside it (the generated fuzz counts as in the Fibres step's line; fuzz not generated yet is said).
std::string yarnFibresSummary(const SceneUi& ui) {
    const knit::SceneSettings& s = ui.s;
    const YarnFibresCost c = yarnFibresCost(ui);
    const knit::FilePreset* look = knit::matchingFilePreset(knit::PresetGroup::YarnFibres, s);
    char b[200];
    std::snprintf(b, sizeof b, "%s%s%d %s%s%d fibres", look ? look->label.c_str() : "Custom look", kDot, s.yarnBody.plies,
                  s.yarnBody.plies == 1 ? "ply" : "plies", kDot, s.yarnBody.fibres);
    std::string t = b;
    if (!c.knitted) return t + kDot + "knit to see the cost";
    if (c.maps) std::snprintf(b, sizeof b, "%smaps %.0f um%s~%.0f M texels", kDot, s.yarnBody.mapTexelUm, kDot, c.texels / 1e6);
    else std::snprintf(b, sizeof b, "%s~%.1f M yarn segments", kDot, c.yarnSegments / 1e6);
    t += b;
    if (s.fibres) {
        if (fibreResult(ui)) std::snprintf(b, sizeof b, " + %.1f M fuzz", c.fuzzSegments / 1e6);
        else std::snprintf(b, sizeof b, " + fuzz not generated yet");
        t += b;
    }
    std::snprintf(b, sizeof b, "%s~%.1f GiB while rendering", kDot, c.gib);
    return t + b;
}

// The Fibres step's cost line (the same numbers in full).
std::string yarnFibresCostLine(const YarnFibresCost& c) {
    if (!c.knitted) return "Knit to see what the yarn's fibres cost and the renderer's GPU memory.";
    char line[300];
    if (c.maps)
        std::snprintf(line, sizeof line,
                      "Yarn maps: ~%.0f M texels (~%.1f GiB of textures) and ~%.1f M tube triangles on %.1f m of yarn. Fuzz: %s. GPU memory "
                      "while rendering: ~%.1f GiB. The maps are baked once per knit and look (minutes), then cached.",
                      c.texels / 1e6, c.texturesGib, c.triangles / 1e6, c.yarnMm / 1000, c.fuzzWord.c_str(), c.gib);
    else
        std::snprintf(line, sizeof line, "Yarn fibres: ~%.1f M segments on %.1f m of yarn. Fuzz: %s. GPU memory while rendering: ~%.1f GiB.",
                      c.yarnSegments / 1e6, c.yarnMm / 1000, c.fuzzWord.c_str(), c.gib);
    return line;
}

// Above the renderer's memory cap: what takes most of it and what to lower (the controls are in the Fibres step).
std::string yarnFibresMemoryAdvice(const YarnFibresCost& c) {
    if (!c.knitted || c.gib <= knit::kRendererMemoryCapGiB) return "";
    char advice[260];
    std::snprintf(advice, sizeof advice, "About %.1f GiB: above the %.0f GiB a 12 GB card can give the renderer. %s", c.gib, knit::kRendererMemoryCapGiB,
                  !c.yarnDominates ? "The fuzz takes most of it: lower the fibre density (Fibres step, Coverage) or choose a lighter Export level."
                  : c.maps         ? "The yarn maps take most of it: choose a coarser Map texel (Fibres step, Yarn fibres)."
                                   : "The yarn's fibres take most of it: lower Fibres per yarn or Interior fibres (Fibres step, Yarn fibres).");
    return advice;
}

Reveal revealOfBodyField(knit::YarnBodyField f) {
    switch (f) {
    case knit::YarnBodyField::ShellDepth:
    case knit::YarnBodyField::FibreAngle:
    case knit::YarnBodyField::InteriorFibres:
    case knit::YarnBodyField::Core:
    case knit::YarnBodyField::Irregularity: return Reveal::Structure;
    case knit::YarnBodyField::None:
    case knit::YarnBodyField::Plies:
    case knit::YarnBodyField::Fibres:
    case knit::YarnBodyField::FibreDiameter:
    case knit::YarnBodyField::MapTexel: return Reveal::Section;
    }
    return Reveal::Section;
}

// --- a running job's progress

namespace {

std::string capitalised(std::string s) {
    if (!s.empty()) s[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(s[0])));
    return s;
}

} // namespace

JobProgressView progressView(const RunningJob& job, bool partialShown) {
    using Kind = knit::BuiltScene::Kind;
    JobProgressView v;
    const knit::JobProgress& p = job.progress;
    const std::string& st = p.stage;
    auto starts = [&](const char* w) { return st.rfind(w, 0) == 0; };
    // The knitting steps count parts (scene_builder.cpp: "knitting", "knitting <part>", "knitted
    // <part>"); the other steps' counters are not parts, and most have none.
    const bool parts = job.kind == Kind::Knit && !job.meshOnly && (st == "knitting" || starts("knitting ") || starts("knitted "));
    if (p.total > 0) v.fraction = std::clamp(float(p.done) / float(p.total), 0.0f, 1.0f);
    if (parts && p.total > 0) v.count = std::to_string(p.done) + " of " + std::to_string(p.total) + (p.total == 1 ? " part" : " parts");
    const char* what = job.kind == Kind::Knit ? (job.bakeMaps ? "Baking yarn maps" : job.meshOnly ? "Meshing" : "Knitting")
                       : job.kind == Kind::Fibres ? "Fibres" : "Exporting";
    const std::string step = st.empty() ? std::string(what) : capitalised(st);
    v.stopping = job.cancelRequested;
    v.title = job.cancelRequested ? "Stopping ..." : std::string(what) + " ...";
    if (job.cancelRequested) v.qualifier = "Cancel requested: stopping after " + (st.empty() ? std::string("the current step") : st);
    else if (parts) v.qualifier = partialShown ? std::string("Partial build shown") + kDot + "not a result yet" : std::string();
    else v.qualifier = step;
    v.status = (job.cancelRequested ? "Cancel requested" + std::string(kDot) + "stopping after: " + (st.empty() ? std::string(what) : st) : step) + kDot +
               "Elapsed " + fmtSeconds(p.seconds);
    return v;
}

// --- numbers as text

std::string fmtCount(long long n) {
    std::string s = std::to_string(n), out;
    for (size_t i = 0; i < s.size(); ++i) {
        if (i && (s.size() - i) % 3 == 0) out += ',';
        out += s[i];
    }
    return out;
}

std::string fmtSeconds(double s) {
    char b[48];
    if (s < 60) std::snprintf(b, sizeof b, "%.0f s", s);
    else if (s < 3600) std::snprintf(b, sizeof b, "%d min %02d s", int(s) / 60, int(s) % 60);
    else std::snprintf(b, sizeof b, "%d h %02d min", int(s) / 3600, int(s) % 3600 / 60);
    return b;
}

double sharePercent(double fraction) { return fraction * 100.0; }

double shareFraction(double percent) { return percent / 100.0; }

std::string candidateEstimate(double densityPerMm, double yarnLengthMm) {
    char b[96];
    std::snprintf(b, sizeof b, "Approx. %.1fM candidates%s%.1f m yarn", densityPerMm * yarnLengthMm / 1e6, kDot, yarnLengthMm / 1000);
    return b;
}

// --- the export's texts

// The detail levels (ExportDetail), from the finest to the lightest; Reduced is the default.
const ExportDetail kExportDetails[] = {
    {"Full", "medium", 14000000, "8 sides, a ring every 0.4 mm, at most 14 M fibre segments per object."},
    {"Reduced", "reduced", 7000000,
     "6 sides, a ring every 0.6 mm, at most 7 M fibre segments per object (the default): about half the GPU memory of Full; the tube's facets "
     "stay within the 0.1 mm clearance."},
    {"Light", "coarse", 5000000,
     "6 sides, a ring every 0.8 mm, at most 5 M fibre segments per object: the tube's facets may reach about 0.13 mm into neighbouring yarns "
     "(beyond the 0.1 mm clearance)."},
};

namespace {

std::string millions(double n) {
    char b[32];
    std::snprintf(b, sizeof b, n >= 1e6 ? "%.1f M" : "%.0f k", n >= 1e6 ? n / 1e6 : n / 1e3);
    return b;
}

} // namespace

const knit::FilePreset* meshLevelOf(const ExportDetail& d) {
    for (const knit::FilePreset& p : knit::presets(knit::PresetGroup::StrandMesh))
        if (p.id == d.meshLevel) return &p;
    return nullptr;  // its yarn tessellation level file is missing (listed in step 1)
}

// The format's short name (the action, "Omitted by ...").
std::string formatShort(const knit::ExportFormat* f) { return f ? f->shortName : std::string("export"); }

// What the export does with the wool fibres (the Fibres footer, Export > Contents): the request
// (s.fibres) against the format and the generated fibres' state. It never changes the request.
FibreExport fibreExport(const SceneUi& ui) {
    const knit::ExportFormat* fmt = knit::findExportFormat(ui.exportTo.format);
    const StageState st = fibreOutputState(ui);
    if (!ui.s.fibres) return FibreExport::Excluded;
    if (fmt && !fmt->fibres) return FibreExport::Omitted;
    if (st.status == StageStatus::Ready) return FibreExport::Included;
    if (runningStage(ui, knit::BuiltScene::Kind::Fibres)) return FibreExport::Generating;
    return st.status == StageStatus::NeedsUpdate ? FibreExport::Update : FibreExport::Generate;
}

// Wool fibres generated for the current knit and settings that the export leaves out (generating
// does not include them): said in the warning colour.
bool fibresGeneratedExcluded(const SceneUi& ui) { return !furShown(ui) && !ui.s.fibres && fibreOutputState(ui).status == StageStatus::Ready; }

// The Fibres footer's export line (read only: the switch is Export > Contents > Include fibres / Include fuzz).
std::string fibreFooterExportNote(const SceneUi& ui) {
    const char* fibres = ui.s.yarnBody.enabled ? "fuzz" : "fibres";  // on a yarn of fibres the generated ones are the fuzz
    switch (fibreExport(ui)) {
    case FibreExport::Excluded: return "Export: not included (Export step, Contents)";
    case FibreExport::Omitted: return std::string("Export: not written, the export format cannot carry ") + fibres;
    case FibreExport::Included: return "Export: included";
    case FibreExport::Generating: return "Export: included once generated";
    case FibreExport::Update: return std::string("Export: included, update the ") + fibres + " first";
    case FibreExport::Generate: return std::string("Export: included, generate the ") + fibres + " first";
    }
    return "";
}

bool exportDetailSelected(const SceneUi& ui, const ExportDetail& d) {
    const knit::FilePreset* level = meshLevelOf(d);
    return level && knit::matchingFilePreset(knit::PresetGroup::StrandMesh, ui.s) == level && ui.s.segmentBudget == d.budget;
}

namespace {

// What an export at one level would hold. A count below zero is not known yet.
struct LevelForecast {
    double triangles = -1;  // the yarn tubes
    double segments = -1;   // the fibre segments as the export writes them (0 without fibres), the yarn's fibres included
    double yarnFibres = 0;  // of them, the yarn made of fibres' (yarn_body; 0 for a tube)
    double textureTexels = 0;  // the yarn maps' texels, every map (a yarn of fibres rendered as maps)
    bool estimate = false;  // the segments scaled from another budget's
    std::string fibres;     // the fibre part of the forecast's text
};

// The late taper the result's counts were made with is not the current one.
bool taperChanged(const knit::FibreResult& fr, const knit::SceneSettings& s) {
    return fr.settings.tailTaper != s.tailTaper || fr.settings.tailTaperMinDiameters != s.tailTaperMinDiameters ||
           fr.settings.tailTaperTip != s.tailTaperTip;
}

// The tube triangles: 2 sides per ring spacing along every tubed curve. BuiltScene::yarnLength sums the
// curves the tubes follow: the plies when the knit has them, their helix included. The fibre
// segments, per object (the budget's unit): the result's exported count at its own budget; at another
// budget min(grown, budget) at the result's ratio of exported to kept (the late taper's vertices).
LevelForecast levelForecast(const SceneUi& ui, const ExportDetail& d, const knit::FilePreset& level) {
    LevelForecast f;
    const knit::SceneSettings& v = level.values;
    const double tubed = ui.results.knit ? ui.results.knit->yarnLength : 0;  // [mm]
    if (tubed > 0) f.triangles = 2.0 * (ui.s.plies > 1 ? v.plyTubeSides / v.plyTubeRing : v.tubeSides / v.tubeRing) * tubed;
    const knit::FibreResult* fr = fibreResult(ui);
    if (!ui.s.fibres) {
        f.segments = 0;
        f.fibres = ui.s.yarnBody.enabled ? "no fuzz" : "no fibres";
    } else if (!fr) {
        f.fibres = "<= " + millions(double(d.budget)) + " fibre segments per object";
    } else if (taperChanged(*fr, ui.s)) {
        f.fibres = "fibre segments: export to count (the late taper changed)";
    } else {
        f.segments = 0;
        for (const knit::ObjectSegments& o : fr->objectSegments) {
            if (d.budget == fr->settings.segmentBudget || o.kept == 0) {
                f.segments += double(o.exported);
                continue;
            }
            const double kept = d.budget > 0 ? std::min(double(o.grown), double(d.budget)) : double(o.grown);
            f.segments += kept * double(o.exported) / double(o.kept);
            f.estimate = true;
        }
        f.fibres = (f.estimate ? "~" : "") + millions(f.segments) + " fibre segments";
    }
    // Made of fibres: the body replaces the tubes (its cores' few triangles are not counted); as maps, its tubes
    // keep their own tessellation (not the level's) and its maps' texels count.
    if (ui.s.yarnBody.enabled && tubed > 0) {
        if (knit::yarnBodyRenderOf(ui.s) == knit::YarnBodyRender::Maps) {
            f.textureTexels = double(knit::yarnMapTexels(ui.s, tubed));
            f.triangles = double(knit::yarnMapTriangles(tubed));
            f.fibres = "~" + millions(f.textureTexels) + " yarn-map texels + " + (ui.s.fibres ? "fuzz: " + f.fibres : f.fibres);
        } else {
            f.yarnFibres = double(knit::yarnBodySegments(ui.s, tubed));
            f.triangles = 0;
            if (f.segments >= 0) f.segments += f.yarnFibres;
            f.fibres = "~" + millions(f.yarnFibres) + " yarn-fibre segments + " + (ui.s.fibres ? "fuzz: " + f.fibres : f.fibres);
        }
    }
    return f;
}

double memoryOf(const LevelForecast& f) {
    if (f.triangles < 0 || f.segments < 0) return -1;
    return knit::rendererMemoryGiB(std::llround(f.segments), std::llround(f.triangles), std::llround(f.textureTexels));
}

} // namespace

std::string exportDetailForecast(const SceneUi& ui, const ExportDetail& d) {
    const knit::FilePreset* level = meshLevelOf(d);
    if (!level) return "level file missing";
    const LevelForecast f = levelForecast(ui, d, *level);
    std::string out = (f.triangles >= 0 ? "~" + millions(f.triangles) + " triangles" : std::string("triangles: knit first")) + ", " + f.fibres;
    if (const double gib = memoryOf(f); gib >= 0) {
        char b[64];
        std::snprintf(b, sizeof b, ", ~%.1f GiB GPU while rendering", gib);
        out += b;
    }
    return out;
}

std::string exportDetailMemoryWarning(const SceneUi& ui, const ExportDetail& d) {
    const knit::FilePreset* level = meshLevelOf(d);
    if (!level) return "";
    const LevelForecast f = levelForecast(ui, d, *level);
    const double gib = memoryOf(f);
    if (gib <= knit::kRendererMemoryCapGiB) return "";
    // The remedy follows what takes the memory: the yarn's fibres (only Fibres per yarn lowers them) or the
    // tubes and fuzz (the level and the fibre density).
    const bool yarnFibresDominate = f.yarnFibres > 0 && f.yarnFibres >= f.segments - f.yarnFibres;
    const double texturesGiB = knit::rendererMemoryGiB(0, 0, std::llround(f.textureTexels)) - knit::rendererMemoryGiB(0, 0, 0);
    const bool mapsDominate = f.textureTexels > 0 && texturesGiB >= gib - texturesGiB - knit::rendererMemoryGiB(0, 0, 0);
    char b[400];
    std::snprintf(b, sizeof b, "About %.1f GiB of GPU memory while the renderer renders, more while it builds: above the %.0f GiB the editor "
                  "and the renderer can share on a 12 GB card. %s", gib, knit::kRendererMemoryCapGiB,
                  yarnFibresDominate ? "The yarn's fibres take most of it: lower Fibres per yarn (Fibres step, Yarn fibres)."
                  : mapsDominate     ? "The yarn maps take most of it: choose a coarser Map texel (Fibres step, Yarn fibres)."
                                     : "Choose a lighter level or a lower fibre density.");
    return b;
}

// What the chosen level sets, or what makes the settings Custom (no level sets this tessellation
// and segment budget together).
std::string exportDetailNote(const SceneUi& ui) {
    const knit::FilePreset* current = knit::matchingFilePreset(knit::PresetGroup::StrandMesh, ui.s);
    const std::string budget = millions(double(ui.s.segmentBudget));
    for (const ExportDetail& d : kExportDetails)
        if (meshLevelOf(d) && exportDetailSelected(ui, d))
            return "Sets the Knit step's Yarn tessellation to " + meshLevelOf(d)->label + " and at most " + millions(double(d.budget)) +
                   " fibre segments.";
    for (const ExportDetail& d : kExportDetails)
        if (current && meshLevelOf(d) == current)
            return std::string("Custom: the tessellation is ") + d.label + "'s (" + current->label + "), but the fibre segment budget is " + budget +
                   " instead of " + millions(double(d.budget)) + " (Fibres step, Generation details).";
    return "Custom: the Knit step's Yarn tessellation is " + (current ? current->label : std::string("custom")) + ", which no level uses; fibre "
           "segment budget " + budget + ".";
}

// What Export > Contents says about the wool fibres (the effective output, fibreExport).
std::string fibreExportWord(const SceneUi& ui) {
    switch (fibreExport(ui)) {
    case FibreExport::Excluded: return fibresGeneratedExcluded(ui) ? "Excluded: fibres are generated, tick to include" : "Excluded";
    case FibreExport::Omitted: return "Omitted by " + formatShort(knit::findExportFormat(ui.exportTo.format));
    case FibreExport::Included: return "Included";
    case FibreExport::Generating: return "Requested, being generated";
    case FibreExport::Update: return "Requested: update them (Fibres)";
    case FibreExport::Generate: return "Requested: generate them (Fibres)";
    }
    return "";
}

// The ready line: the format and what goes with the knit (or the mesh).
std::string exportReadyLine(const SceneUi& ui) {
    const knit::ExportFormat* fmt = knit::findExportFormat(ui.exportTo.format);
    std::string with;
    if (furShown(ui)) {
        with = fmt && fmt->fibres ? "with the fur" : "fur omitted by " + formatShort(fmt);
    } else {
        // A yarn made of fibres: they come with the yarn; the switch is the fuzz's.
        const std::string fuzz = ui.s.yarnBody.enabled ? "fuzz" : "fibres";
        switch (fibreExport(ui)) {
        case FibreExport::Excluded: with = "no " + fuzz + (fibresGeneratedExcluded(ui) ? " (generated, not included)" : ""); break;
        case FibreExport::Omitted: with = fuzz + " omitted by " + formatShort(fmt); break;
        case FibreExport::Included: with = fuzz + " included"; break;
        default: with = fuzz + " requested"; break;
        }
        if (ui.s.yarnBody.enabled)
            with = (knit::yarnBodyRenderOf(ui.s) == knit::YarnBodyRender::Maps ? "yarn of fibres as maps" : "yarn made of fibres") + std::string(kDot) + with;
    }
    return (fmt ? fmt->name : std::string("?")) + kDot + with;
}

// Destination: the folder the export job writes into, or why it cannot.
std::string exportOutputText(const SceneUi& ui) {
    const std::string problem = exportTargetProblem(ui);
    return problem.empty() ? ui.exportTo.dir : "Cannot export there: " + problem + ".";
}

// The last successful export: the folder that job published (its own result, not the current
// destination), "" if there is none.
std::string exportDoneText(const SceneUi& ui) {
    return ui.results.exported && ui.results.exportAttempt.outcome == knit::JobOutcome::Done ? "Exported to " + ui.results.exported->exportedTo : std::string();
}

// Fur: the export writes the mesh itself.
bool exportsMesh(const SceneUi& ui) { return furShown(ui); }

// The format's contents in one line; Fur: the mesh and its fur (no yarn or ply maps).
std::string formatSummary(const SceneUi& ui, const knit::ExportFormat& f) { return exportsMesh(ui) ? f.meshSummary : f.summary; }

// Destination: the warning on an output folder with files in it (publish.h, IntoFolder), "" if
// it is missing or empty.
std::string exportExistingNote(const SceneUi& ui) {
    std::error_code ec;
    if (!exportTargetProblem(ui).empty() || !fs::exists(ui.exportTo.dir, ec) || fs::is_empty(ui.exportTo.dir, ec)) return "";
    return "This folder is not empty. The export writes into it: files of the same name are replaced, all other files stay.";
}

// The footer's heading of a blocked export, named after its fix (wording only: the eligibility and
// the action read el.fix itself).
const char* blockedHeading(const std::string& fix) {
    if (fix == "Knit at Final quality") return "Final knit required";
    if (fix == "Generate fibres") return "Fibres required";
    if (fix == "Update fibres") return "Fibres outdated";
    if (fix == "Grow fur") return "Fur required";
    if (fix == "Update fur") return "Fur outdated";
    if (fix == kChooseWorkflow) return "Knit or Fur not chosen";
    if (fix == "Update the knit") return "Knit outdated";
    if (fix == "Correct the yarn's fibres") return "Yarn fibres refused";
    return "Export not available";
}

std::string defaultExportName(const SceneUi& ui) {
    const std::string suffix = furShown(ui) ? "_fur" : "_knit";
    if (ui.s.scene == "import" && !ui.s.importSource.empty()) return fs::path(ui.s.importSource).stem().string() + suffix;
    return ui.s.scene + suffix;
}

bool exportNamesFiles(const SceneUi& ui) {
    const knit::ExportFormat* f = knit::findExportFormat(ui.exportTo.format);
    return f && !f->extensions.empty();
}

std::string exportHoldsText(const SceneUi& ui) {
    if (!exportNamesFiles(ui)) return "";
    const std::vector<std::string>& extensions = knit::findExportFormat(ui.exportTo.format)->extensions;
    std::string files;
    for (size_t i = 0; i < extensions.size(); ++i) files += (i ? " and " : "") + ui.exportTo.name + extensions[i];
    return "Holds " + files + ", with report.json.";
}

std::string exportFolderStart(const SceneUi& ui) {
    std::error_code ec;
    const fs::path dir(ui.exportTo.dir);  // the default folder may not exist yet
    return (fs::is_directory(dir, ec) ? dir : dir.parent_path()).string();
}

std::string exportTargetProblem(const SceneUi& ui) {
    std::string why;
    if (exportNamesFiles(ui) && !knit::validateExportName(ui.exportTo.name, &why)) return why;
    // The folder itself, or where the default one will be created.
    std::error_code ec;
    const fs::path dir(ui.exportTo.dir);
    if (ui.exportTo.dir.empty() || !fs::is_directory(fs::exists(dir, ec) ? dir : dir.parent_path(), ec))
        return "the destination folder does not exist: choose another one";
    return "";
}

// --- the Import step's open file

bool importFooterShown(const SceneUi& ui) {
    const knit::ImportSession::Snapshot file = ui.import.session.snapshot();
    return ui.import.fileOpen && file.mesh && file.analysis;
}

bool nodeKeptAsMesh(const ImportFlow& im, const std::string& name) { return im.options.plainNodes.count(name) > 0; }

knit::ImportOptions importRequest(const ImportFlow& im) {
    return im.options;
}

double tickedNodeArea(const SceneUi& ui) {
    const std::shared_ptr<const knit::ImportedMesh> mesh = ui.import.session.snapshot().mesh;
    double area = 0;
    if (!mesh) return area;
    for (const knit::ImportNode& node : mesh->nodes)
        if (!nodeKeptAsMesh(ui.import, node.name)) area += node.area;
    return area;
}

int knittedNodeCount(const SceneUi& ui) {
    const std::shared_ptr<const knit::ImportedMesh> mesh = ui.import.session.snapshot().mesh;
    int count = 0;
    if (!mesh) return count;
    for (const knit::ImportNode& node : mesh->nodes) count += nodeKeptAsMesh(ui.import, node.name) ? 0 : 1;
    return count;
}

// The file's area counts its interiors too, which a rebuilt surface removes: an upper bound.
double importStitchCount(const SceneUi& ui) {
    const double u = ui.import.options.unitMm;
    return tickedNodeArea(ui) * u * u / (ui.s.wale * ui.s.course);
}

// The same mapping as the import itself (prepareImport) and its preview (addImportPreview): x stays the
// width; with Y up the file's z is the depth and y the height, with Z up y is the depth.
ImportSize importSize(const SceneUi& ui) {
    const std::shared_ptr<const knit::ImportAnalysis> an = ui.import.session.snapshot().analysis;
    ImportSize size;
    if (!an) return size;
    const Vec3 ext = an->hi - an->lo;
    const double u = ui.import.options.unitMm;
    const bool yUp = ui.import.options.up == knit::UpAxis::Y;
    size.width = ext.x * u;
    size.depth = (yUp ? ext.z : ext.y) * u;
    size.height = (yUp ? ext.y : ext.z) * u;
    size.valid = u > 0 && std::isfinite(u);
    return size;
}

std::string sizeProblem(const SceneUi& ui) {
    const double u = ui.import.options.unitMm;
    if (u > 0 && std::isfinite(u)) return "";
    if (!ui.import.sizeByHeight) return "Enter a positive size for one file unit.";
    const std::shared_ptr<const knit::ImportAnalysis> an = ui.import.session.snapshot().analysis;
    const Vec3 ext = an ? an->hi - an->lo : Vec3{};
    const double fileHeight = ui.import.options.up == knit::UpAxis::Y ? ext.y : ext.z;
    if (an && !(fileHeight > 0)) return "The object is flat along the up axis: set its size with File units.";
    return "Enter a positive height.";
}

// The unit ImportAnalysis suggests is one of these four (suggestUnit in import_analysis.cpp).
std::string suggestedUnitName(double unitMm) {
    if (unitMm == 1) return "Millimetres";
    if (unitMm == 10) return "Centimetres";
    if (unitMm == 1000) return "Metres";
    return "Inches";
}

std::string importSourceLine(const SceneUi& ui) {
    const std::shared_ptr<const knit::ImportedMesh> mesh = ui.import.session.snapshot().mesh;
    if (!mesh) return "";
    char b[160];
    std::snprintf(b, sizeof b, "%s%s%.1f MB%s%zu component%s", mesh->format.c_str(), kDot, mesh->fileBytes / 1048576.0, kDot, mesh->nodes.size(),
                  mesh->nodes.size() == 1 ? "" : "s");
    return b;
}

std::string importColourSummary(const SceneUi& ui) {
    const knit::ImportSession::Snapshot file = ui.import.session.snapshot();
    if (!file.mesh || !file.analysis) return "";
    const std::string& source = file.analysis->colourSource;
    const std::string head = std::string("Base colour") + kDot;
    if (source.rfind("texture", 0) == 0) {
        int textures = 0;  // the materials colourSource names, one per material with a texture
        for (const knit::ImportMaterial& mat : file.mesh->materials) textures += mat.texture.empty() ? 0 : 1;
        return head + std::to_string(textures) + (textures == 1 ? " texture" : " textures");
    }
    if (source == "none") return head + "none in the file (the default yarn colour)";
    return head + source;
}

ImportFooterView importFooterView(const SceneUi& ui) {
    const knit::ImportSession::Snapshot file = ui.import.session.snapshot();
    const ImportSize size = importSize(ui);
    char high[64] = "";
    if (size.valid) std::snprintf(high, sizeof high, "%s%.1f mm high", kDot, size.height);
    ImportFooterView v;
    const int knitted = knittedNodeCount(ui), total = file.mesh ? static_cast<int>(file.mesh->nodes.size()) : 0;
    const bool fur = ui.workflowChosen && knit::furWorkflow(ui.s);
    if (fur) {
        v.summary = std::to_string(knitted) + " with fur" + kDot + std::to_string(total - knitted) + " kept as mesh" + high;
        v.status = "The file's own mesh";
        v.action = "Import for fur";
    } else {
        v.summary = std::to_string(knitted) + " to knit" + kDot + std::to_string(total - knitted) + " kept as mesh" + high;
        v.status = "Closed surface";
        v.action = "Rebuild & import";
    }
    if (file.analysis && file.analysis->fatal()) v.whyDisabled = "The file has a problem that cannot be imported (above).";
    else if (!size.valid) v.whyDisabled = sizeProblem(ui);
    else if (!ui.workflowChosen) v.whyDisabled = std::string(kChooseWorkflow) + " for the object (above).";
    else if (!(tickedNodeArea(ui) > 0)) v.whyDisabled = fur ? "Nothing to grow fur on: tick a component." : "Nothing to knit: tick a component.";
    if (!v.whyDisabled.empty()) v.status = v.whyDisabled;
    return v;
}

} // namespace viewer
