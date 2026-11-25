// The Fibres step (interface.h): the settings panel (wool fibres on the yarn, or hair on the mesh)
// and the footer with the result's state and its one action.
#include "viewer/interface.h"
#include "jobs/scene_builder.h"
#include "settings/scene_registry.h"
#include "settings/settings_mapping.h"

#include <algorithm>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

#include "settings/file_presets.h"

namespace viewer {

namespace {

// The mask file's hash as it is on disk now, for the warning drawn every frame: the file is read again only
// when its time stamp or size changes (fileHashHex reads all of it).

// A strand material's look as label / value rows (inside a beginRows table): the shading model and
// its parameters (the sheen, or the hair BSDF's roughnesses).
void materialAppearance(knit::MaterialSpec& m) {
    int model = m.model == knit::MaterialSpec::Model::Fibre ? 1 : 0;
    const char* models[] = {"Matte with sheen", "Hair (fibre) shading"};
    if (comboRow("Shading", &model, models, 2, "How the renderer shades the fibres: a diffuse surface with a sheen layer, or a hair BSDF (light travels along the fibre)."))
        m.model = model ? knit::MaterialSpec::Model::Fibre : knit::MaterialSpec::Model::DiffuseSheen;
    if (m.model == knit::MaterialSpec::Model::Fibre) {
        dragRow("Roughness along", &m.fibreLongRoughness, {0.01, 1, "%.2f"});
        dragRow("Roughness around", &m.fibreAzimRoughness, {0.01, 1, "%.2f"});
        dragRow("IOR", &m.fibreIor, {1.3, 1.8, "%.2f"}, "The fibres' refractive index: lower gives a softer, less glassy shine.");
    } else {
        dragRow("Sheen", &m.sheen, {0, 1, "%.2f"}, "The soft highlight of fuzzy surfaces at grazing angles.");
    }
}

// A group of related number fields under one title ("Length · mm": Minimum / Typical / Maximum):
// side by side under muted headers when the panel is wide enough, else as label / value rows.
// The unit is in the title. True when a value changed.
struct NumberField {
    const char* label;  // the header / row label
    const char* id;     // the field's id ("##..."), unique in its population
    double* v;
    double lo, hi;
    const char* fmt;
};

bool fieldGroup(const char* title, const char* helpText, const NumberField* f, int n) {
    ImGui::TextUnformatted(title);
    if (helpText) help(helpText);
    bool changed = false;
    if (ImGui::GetContentRegionAvail().x >= ImGui::GetFontSize() * 5.5f * n) {
        if (ImGui::BeginTable(title, n, ImGuiTableFlags_SizingStretchSame)) {
            ImGui::TableNextRow();
            for (int i = 0; i < n; ++i) {
                ImGui::TableNextColumn();
                ImGui::TextDisabled("%s", f[i].label);
            }
            ImGui::TableNextRow();
            for (int i = 0; i < n; ++i) {
                ImGui::TableNextColumn();
                ImGui::SetNextItemWidth(-FLT_MIN);
                changed |= dragField(f[i].id, f[i].v, {f[i].lo, f[i].hi, f[i].fmt});
            }
            ImGui::EndTable();
        }
    } else if (beginRows(title)) {
        for (int i = 0; i < n; ++i) {
            row(f[i].label);
            changed |= dragField(f[i].id, f[i].v, {f[i].lo, f[i].hi, f[i].fmt});
        }
        ImGui::EndTable();
    }
    return changed;
}

// One fibre population's editor (a tab of Fibre shape & appearance): its shape, then its
// appearance. Directed = the yarn fuzz (bow, lift), else the stray fibres (bend, start angle).
void populationEditor(const char* id, const char* explanation, knit::SceneSettings::Population& p, bool directed, knit::SceneSettings& s) {
    ImGui::PushID(id);
    mutedWrapped(explanation);
    const double share = directed ? 1 - s.randomFraction : s.randomFraction;
    if (share <= 0) wrapped(kWarn, "Share 0 % (Coverage): none of these are generated; their values are kept.");
    const NumberField lengths[] = {{"Minimum", "##Shortest", &p.lengthMin, 0.1, 4, "%.2f"},
                                   {"Typical", "##Typical", &p.lengthMedian, 0.1, directed ? 8.0 : 20.0, "%.2f"},
                                   {"Maximum", "##Longest", &p.lengthMax, 0.1, directed ? 10.0 : 60.0, "%.2f"}};
    fieldGroup("Length \xC2\xB7 mm", "Lengths follow a log-normal distribution around the typical length, cut at the minimum and the maximum.",
               lengths, 3);
    if (beginRows("##spreadRows")) {
        dragRow("Length spread", &p.lengthSpread, {0, directed ? 0.8 : 1.5, "%.2f"}, "The log-normal spread (sigma) of the lengths around the typical length.");
        ImGui::EndTable();
    }
    double tmin = p.thicknessMin(), tmax = p.thicknessMax();
    const NumberField thickness[] = {{"Minimum", "##Thinnest", &tmin, 0.01, 0.2, "%.3f"}, {"Maximum", "##Thickest", &tmax, 0.01, 0.2, "%.3f"}};
    if (fieldGroup("Thickness \xC2\xB7 mm", "The fibre diameter: a log-normal distribution whose +-2 sigma limits are the minimum and the maximum.",
                   thickness, 2))
        p.setThicknessRange(std::min(tmin, tmax), std::max(tmin, tmax));
    if (beginRows("##shapeRows")) {
        row(directed ? "Bow" : "Bend", directed ? "How much each fibre bows away from a straight line." : "How much each stray fibre curls along its length.");
        dragField(directed ? "##Bow" : "##Bend", &p.curvature, {0, directed ? 3 : 0.6, "%.2f"});
        row(directed ? "Max lift angle" : "Max angle",
            directed ? "The largest angle at which a fibre leaves the yarn (nominal: bow and the yarn's curvature change the measured angle)."
                     : "The largest angle between a stray fibre's start and the surface.");
        dragField(directed ? "##Max lift angle" : "##Max angle", &p.maxAngleDeg, {directed ? 12.0 : 20.0, directed ? 45.0 : 80.0, "%.0f deg"});
        if (directed) {
            row("Hug height", "How far above the yarn the fuzz stays, in yarn radii: the fibres lie along the yarn instead of standing off. "
                              "0: off (the lift angle alone).");
            dragField("##Hug height", &s.directedHugHeight, {0, 1, "%.2f"});
            row("Lie on the stitches", "A fibre that reaches a neighbouring stitch rests on it and lifts over it, instead of being cut "
                                       "where it would sink in (long fuzz, as mohair). Light fibre checks only.");
            ImGui::Checkbox("##Lie on the stitches", &s.directedDrape);
        }
        if (!directed) {
            row("Direction spread", "Seen from above the yarn, how far a stray fibre's start may turn away from the yarn's direction. "
                                    "180: any direction; lower: the strays follow the yarn.");
            dragField("##Direction spread", &s.randomSpreadDeg, {5, 180, "%.0f deg"});
            row("Curl radius \xC2\xB7 mm", "How far a stray fibre winds around its line; the largest radius, each fibre draws its own. 0: no curl.");
            dragField("##Curl radius", &p.curlRadius, {0, 4, "%.2f"});
            row("Curl length \xC2\xB7 mm", "The fibre length of one turn of the curl. Very long fibres stretch it: they get at most a sixteenth of "
                                            "Segments per fibre in turns.");
            dragField("##Curl length", &p.curlLength, {2, 40, "%.1f"});
        }
        ImGui::EndTable();
    }
    if (ImGui::TreeNode("Appearance")) {
        bool& matchYarn = directed ? s.fibreInheritsYarn : s.fibreRandomInheritsYarn;
        ImGui::Checkbox("Match yarn colour", &matchYarn);
        tooltip(directed ? "The yarn fuzz takes the colour of the yarn it grows on. Off: the colour below."
                         : "Each stray fibre shows the colour of the yarn it grows on. Off: the colour below; with hair "
                           "shading that is the colour of a dense mass, so one fibre is nearly clear.");
        knit::MaterialSpec& m = directed ? s.fibreDirected : s.fibreRandom;
        if (beginRows("##appearanceRows")) {
            if (!matchYarn) colourRow("Colour", m.color);
            materialAppearance(m);
            ImGui::EndTable();
        }
        ImGui::TreePop();
    }
    ImGui::PopID();
}

// The wool presets (sets of values within the wool fibres).
void fibreSetupSection(SceneUi& ui) {
    knit::SceneSettings& s = ui.s;
    ImGui::SeparatorText("Fibre setup");
    const knit::FilePreset* match = nullptr;
    if (beginRows("##setupRows")) {
        match = presetRow(s, knit::PresetGroup::Fibres, "Wool preset",
                          "Your own values (editing a value below makes the wool preset Custom; Direction spread, Max angle and "
                          "Length spread of the stray fibres belong to the stray preset).");
        presetRow(s, knit::PresetGroup::Strays, "Stray preset",
                  "Your own values (Direction spread, Max angle and Length spread on the Stray fibres tab).");
        presetRow(s, knit::PresetGroup::Density, "Density preset", "Your own density (Coverage, below).");
        ImGui::EndTable();
    }
    s.fibrePreset = match ? match->id : "custom";
    for (const std::string& err : knit::presetErrors(knit::PresetGroup::Fibres)) wrapped(kWarn, "Preset not loaded: " + err);
    for (const std::string& err : knit::presetErrors(knit::PresetGroup::Strays)) wrapped(kWarn, "Stray preset not loaded: " + err);
    for (const std::string& err : knit::presetErrors(knit::PresetGroup::Density)) wrapped(kWarn, "Density preset not loaded: " + err);
}

// Wool: how many candidates (per mm of yarn) and which share of them are stray fibres.
void coverageSection(SceneUi& ui) {
    knit::SceneSettings& s = ui.s;
    ImGui::SeparatorText("Coverage");
    if (beginRows("##coverageRows")) {
        dragRow("Density", &s.fibreDensity, {1, 1000, "%.0f / mm"}, "Candidate fibres per millimetre of yarn; some are rejected (inside the object, colliding), so fewer remain.");
        row("Stray fibres", "The share of candidates that are stray fibres; the rest are the yarn fuzz.");
        double percent = sharePercent(s.randomFraction);
        if (dragField("##Stray fibres", &percent, {0, 100, "%.0f%%"})) s.randomFraction = shareFraction(percent);
        ImGui::EndTable();
    }
    if (ui.results.knit && ui.results.knit->yarnLength > 0) {
        ImGui::TextDisabled("%s", candidateEstimate(s.fibreDensity, ui.results.knit->yarnLength).c_str());
        help("Candidates on the knitted yarn, before any are rejected: not the number of fibres generated or exported.");
    }
}

// Wool: patches where the fuzz differs from the even fuzz around them.
void islandSection(SceneUi& ui) {
    knit::SceneSettings& s = ui.s;
    ImGui::SeparatorText("Fibre islands");
    if (beginRows("##islandRows")) {
        presetRow(s, knit::PresetGroup::FibreIslands, "Island preset", "Your own values (editing any value below makes the settings Custom).");
        dragRow("Islands", &s.islands.rate, {0, 2, "%.2f / cm2"}, "Patches of different fibres, per square centimetre.");
        dragRow("Island size", &s.islands.size, {2, 40, "%.0f mm"}, "How wide an island is.");
        dragRow("Island strength", &s.islands.strength, {0, 100, "%.0f%%"}, "How far an island's centre departs from the fuzz around it.");
        dragRow("Fuzzy share", &s.islands.fuzzy, {0, 5, "%.1f"}, "Islands with more and longer stray fibres.");
        dragRow("Matted share", &s.islands.matted, {0, 5, "%.1f"}, "Islands whose fibres lie flat on the yarn.");
        dragRow("Worn share", &s.islands.worn, {0, 5, "%.1f"}, "Islands with fewer fibres, where the yarn shows.");
        ImGui::EndTable();
    }
    for (const std::string& err : knit::presetErrors(knit::PresetGroup::FibreIslands)) wrapped(kWarn, "Island preset not loaded: " + err);
}

// Wool: the two populations' shape and appearance, one editor tab each. The tab only chooses what
// is edited: both populations are generated (by their share), whichever tab is open.
void fibreShapeSection(SceneUi& ui) {
    knit::SceneSettings& s = ui.s;
    ImGui::SeparatorText("Fibre shape & appearance");
    if (!ImGui::BeginTabBar("##populations")) return;
    const int request = std::exchange(ui.panels.fibreTabRequest, -1);
    char fuzz[64], stray[64];
    std::snprintf(fuzz, sizeof fuzz, "Yarn fuzz %.0f%%###fuzz", sharePercent(1 - s.randomFraction));
    std::snprintf(stray, sizeof stray, "Stray fibres %.0f%%###stray", sharePercent(s.randomFraction));
    // (A tab reports its label before it is added: testLabel names it afterwards.)
    const ImGuiID fuzzId = ImGui::GetID(fuzz), strayId = ImGui::GetID(stray);
    const bool fuzzOpen = ImGui::BeginTabItem(fuzz, nullptr, request == 0 ? ImGuiTabItemFlags_SetSelected : 0);
    testLabel(fuzzId, fuzz);
    if (fuzzOpen) {
        populationEditor("dir", "Short fibre ends that follow the yarn.", s.directed, true, s);
        ImGui::EndTabItem();
    }
    const bool strayOpen = ImGui::BeginTabItem(stray, nullptr, request == 1 ? ImGuiTabItemFlags_SetSelected : 0);
    testLabel(strayId, stray);
    if (strayOpen) {
        populationEditor("rnd", "Loose fibre ends that protrude from the yarn.", s.random, false, s);
        ImGui::EndTabItem();
    }
    ImGui::EndTabBar();
}

// Wool: the generation-wide settings (both populations), behind their own node; an unavailable GPU
// is said outside it too.
void generationDetailsSection(SceneUi& ui) {
    knit::SceneSettings& s = ui.s;
    const knit::ComputeCaps& caps = knit::computeCaps();
    ImGui::Spacing();
    if (!caps.device) wrapped(kWarn, "No usable GPU (" + caps.whyNoDevice + "): fibres cannot be generated here.");
    if (!ImGui::TreeNode("Generation details")) return;
    if (beginRows("##generationRows")) {
        row("Segment budget", "The most fibre segments the export may hold per object (each bear, an imported model; whole fibres beyond it are dropped, never shortened).");
        long long budget = s.segmentBudget;
        if (ImGui::InputScalar("##Segment budget", ImGuiDataType_S64, &budget)) s.segmentBudget = std::max(0LL, budget);
        dragRow("More at curved yarn", &s.curvatureBias, {0, 1, "%.2f"});
        dragRow("Segments per fibre", &s.fibreMaxSegments, 2, 48, "The most straight pieces a fibre is made of (long fibres need more to stay smooth; 48 is the GPU's limit).");
        row("Segment length \xC2\xB7 mm", "The length of one straight piece. Longer pieces spend the segment budget on more fibres, which "
                                            "cover the yarn better. Curled stray fibres follow their curl instead.");
        dragField("##Segment length", &s.fibreSegmentLength, {0.1, 4, "%.2f"});
        int gen = s.fibreGeneration == "full" ? 1 : 0;
        const char* gens[] = {"Budget order (only what the budget keeps)", "Every candidate"};
        if (comboRow("Generation", &gen, gens, 2, "Budget order evaluates only the candidates the segment budget keeps (the same fibres, faster); on the CPU and the GPU."))
            s.fibreGeneration = gen ? "full" : "budget";
        ImGui::EndTable();
    }
    ImGui::TreePop();
}

} // namespace

// The Fibres step: the heading, for a yarn made of fibres its Yarn fibres section first (panel_yarn_fibres.cpp),
// the wool presets, then the wool fibres' settings; its job and result are in the Fibres footer (fibreFooter).
// Under Fur the same place holds the Fur step (panel_fur.cpp).
void fibrePanel(SceneUi& ui, const knit::SceneEntry&) {
    if (furShown(ui)) {
        furPanel(ui);
        return;
    }
    ImGui::PushID("fibres");
    const bool yarnFibres = ui.s.yarnBody.enabled;
    stepHeading("Fibre settings", std::string(yarnFibres ? "Yarn fibres and the fuzz on them" : "Wool on knitted yarn") + kDot + "applies to the whole object");
    if (!ui.workflowChosen) wrapped(kWarn, std::string(kChooseWorkflow) + " for the object in the Import step.");
    if (yarnFibres) yarnFibresSection(ui);
    fibreSetupSection(ui);
    coverageSection(ui);
    islandSection(ui);
    fibreShapeSection(ui);
    generationDetailsSection(ui);
    ImGui::PopID();
}

namespace {

// The generated fibres or fur, once nothing runs for them.
void fibreStatus(SceneUi& ui, const StageState& st) {
    const bool fur = furShown(ui);
    const knit::FibreResult* f = fibreResult(ui);
    std::string text = st.detail;
    if (f) {
        char b[200];
        if (fur) std::snprintf(b, sizeof b, "Fur #%d: %s fibres", f->id, fmtCount(f->fibres).c_str());
        else std::snprintf(b, sizeof b, "#%d of knit #%d: %s fibres", f->id, f->knitId, fmtCount(f->fibres).c_str());
        text = b + (st.detail.empty() ? std::string() : " - " + st.detail);
    }
    ImGui::TextDisabled(fur ? "Fur" : ui.s.yarnBody.enabled ? "Fuzz" : "Fibres");
    ImGui::SameLine();
    stateCard(st, text);
    // (Once: without a result the attempt's message is already the state's detail.)
    attemptNote(ui.results.fibreAttempt, st.detail);
    if (!ui.workflowChosen) mutedWrapped(std::string(kChooseWorkflow) + " in the Import step.");
    else if (!fibresInputReady(ui)) mutedWrapped("Fibres grow on a knitted result that matches the current knitting settings.");
    else if (!fur) mutedWrapped("Reuses the current knit.");
    if (f && st.status != StageStatus::NeedsUpdate && ui.view.fibres == FibreView::Reference)
        mutedWrapped(std::string("Applied to the model: choose Model (Preview target) to see ") + (fur ? "the fur." : "the fibres."));
}

// The export: wool fibres only when included (the switch is Export > Contents > Include fibres; this
// line only says what the export does); the fur always goes with its mesh.
void fibreExportLine(const SceneUi& ui) {
    if (furShown(ui)) {
        mutedWrapped("The export writes the mesh with this fur.");
        return;
    }
    const FibreExport fx = fibreExport(ui);
    const bool warn = fx == FibreExport::Update || fx == FibreExport::Generate || fibresGeneratedExcluded(ui);
    wrapped(warn ? kWarn : muted(), fibreFooterExportNote(ui));
    tooltip(ui.s.yarnBody.enabled ? "Whether the export writes the fuzz is chosen in the Export step (Contents, Include fuzz). Generating does not include it."
                                  : "Whether the export writes the fibres is chosen in the Export step (Contents, Include fibres). Generating does not include them.");
}

// The one action: Cancel, what is missing, or applying the settings to the model. On a yarn made of fibres
// the action is the fuzz's by name (Generate fuzz / Update fuzz): the yarn's own fibres need no action.
void fibreAction(SceneUi& ui, const StageState& st) {
    const bool fur = furShown(ui);
    const bool fuzz = !fur && ui.s.yarnBody.enabled;
    const knit::FibreResult* f = fibreResult(ui);
    if (runningStage(ui, knit::BuiltScene::Kind::Fibres)) {
        if (cancelButton(ui.job.cancelRequested ? "Cancel requested ..." : (fur ? "Cancel fur" : "Cancel fibres"), ui.job.cancelRequested)) cancelJob(ui);
    } else if (!ui.workflowChosen) {
        if (primaryButton(kChooseWorkflow, true)) ui.stage = Stage::Import;
    } else if (!fibresInputReady(ui)) {
        if (primaryButton(knitResult(ui) ? "Update the knit first" : "Knit first", !ui.job.active, waitReason(ui).c_str())) {
            ui.stage = Stage::Knit;
            startKnit(ui);
        }
    } else {
        const bool update = f && st.status == StageStatus::NeedsUpdate;
        // One action from either preview target: the reference sample follows the edits by itself;
        // this computes the whole model with the same settings (the view stays where it is).
        const char* label = fur    ? (!f ? "Grow fur" : update ? "Update fur" : "Grow fur again")
                            : fuzz ? (!f ? "Generate fuzz" : update ? "Update fuzz" : "Generate fuzz again")
                                   : !f || update ? "Apply to model" : "Apply again";
        if (primaryButton(label, !ui.job.active, waitReason(ui).c_str())) startFibres(ui);
        tooltip(fur  ? "Grows the fur on the whole mesh with these settings; the view stays. The export writes exactly the fur grown "
                       "here."
                : fuzz ? "Generates the fuzz on the whole knitted model with these settings (the current knit is reused); the view "
                         "stays. The export writes exactly the fuzz generated here. The yarn's own fibres (Yarn fibres) need no action."
                       : "Generates the fibres on the whole knitted model with these settings (the current knit is reused); the view "
                         "stays. The export writes exactly the fibres generated here.");
    }
}

} // namespace

// The Fibres footer (both profiles; the step's one place for its job): the generated fibres' state,
// whether the export includes them, and the one action. The status block scrolls past its share of
// the panel; the export and action rows keep their room.
void fibreFooter(SceneUi& ui) {
    ImGui::PushID("fibreFooter");
    const StageState st = fibreOutputState(ui);
    const float actionH = ImGui::GetFrameHeight() * 1.35f + ImGui::GetStyle().ItemSpacing.y;
    beginFooterStatus("##fibreStatus", actionH + ImGui::GetFrameHeightWithSpacing() * 2);
    if (runningStage(ui, knit::BuiltScene::Kind::Fibres)) {
        // Work done (FibreProgress, jobs/job_steps.h: a share of the work, not a time estimate) and the elapsed time.
        JobProgressView v = progressView(ui.job, false);
        if (!v.stopping) v.title = furShown(ui) ? "Growing fur ..." : "Generating fibres ...";
        char work[48] = "";
        if (ui.job.progress.total > 0) std::snprintf(work, sizeof work, "%d %% of the work", int(v.fraction * 100));
        jobProgressBlock(v, work);
        mutedWrapped("Elapsed " + fmtSeconds(ui.job.progress.seconds) + kDot + "the previous result stays visible until the new one is complete.");
    } else {
        fibreStatus(ui, st);
    }
    ImGui::EndChild();
    fibreExportLine(ui);
    fibreAction(ui, st);  // always drawn last, with its room
    ImGui::PopID();
}

} // namespace viewer
