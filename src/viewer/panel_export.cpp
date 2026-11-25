// The Export step (interface.h): the settings panel (format, contents, tessellation, destination) and
// the footer with the export's readiness and its one action.
#include "viewer/interface.h"
#include "settings/scene_registry.h"

#include <algorithm>
#include <cstdio>
#include <string>

#include "export/export.h"
#include "settings/yarn_look.h"
#include "settings/file_presets.h"

namespace viewer {

namespace {

// The Export step's heading.
void exportContext(const knit::SceneEntry& e) {
    stepHeading("Export settings", e.label + kDot + "format, contents and destination");
    if (e.isPreset()) mutedWrapped("Made from " + knit::creditLine(e.credit) + ". Exports name it in CREDITS.txt.");
}

// File format: the format, its contents in one line and in full (a disclosure), and what it
// cannot carry (always visible).
void formatSection(SceneUi& ui) {
    const knit::SceneSettings& s = ui.s;
    ImGui::SeparatorText("File format");
    const knit::ExportFormat* fmt = knit::findExportFormat(ui.exportTo.format);
    ImGui::SetNextItemWidth(-FLT_MIN);
    if (beginList("##Format", fmt ? fmt->name.c_str() : "?")) {
        for (const knit::ExportFormat& f : knit::exportFormats()) {
            if (ImGui::Selectable(f.name.c_str(), fmt == &f, f.available ? 0 : ImGuiSelectableFlags_Disabled)) ui.exportTo.format = f.id;
            if (!f.available) tooltip(f.unavailable.c_str());
        }
        ImGui::EndCombo();
    }
    fmt = knit::findExportFormat(ui.exportTo.format);
    if (!fmt) return;
    mutedWrapped(formatSummary(ui, *fmt));
    if (ImGui::TreeNode(fmt->contentsHeading.c_str())) {
        mutedWrapped(fmt->contents);
        if (exportsMesh(ui)) mutedWrapped("With knitting skipped: the mesh, with its base colour texture, instead of the yarn tubes.");
        ImGui::TreePop();
    }
    if (!fmt->omits.empty()) wrapped(kWarn, "Not included: " + fmt->omits);
    if (furShown(ui) && !fmt->fibres)
        wrapped(kWarn, "The fur is not included: this format holds the base meshes only (the fur needs the photon_path_tracer package).");
    if (const std::string& coordinates = s.scene == "import" ? fmt->importCoordinates : fmt->coordinates; !coordinates.empty())
        mutedWrapped(coordinates);
}

// Contents: what the export writes besides the knitted yarn. Knit skipped: no fibre or filling choice.
// The colour source of the surfaces that are not knitted does not depend on the Knit step.
void contentsSection(SceneUi& ui) {
    knit::SceneSettings& s = ui.s;
    ImGui::SeparatorText("Contents");
    const knit::ExportFormat* fmt = knit::findExportFormat(ui.exportTo.format);
    if (exportsMesh(ui)) {
        // Fur: the mesh and its fur; nothing knitted, so no fibre choice.
        mutedWrapped(fmt && fmt->fibres ? "The mesh and the fur grown in the Fur step." : "The mesh only.");
        mutedWrapped("The Fur workflow knits nothing: no yarn tubes or ply maps are written.");
    }
    if (!beginRows("##contentsRows")) return;
    if (!furShown(ui)) {
        // A yarn made of fibres: its fibres come with the yarn (Knit step, Construction); the switch is the
        // fuzz's, the fibres grown in step 2 on top of them.
        const bool fibreYarn = s.yarnBody.enabled;
        if (fibreYarn) {
            const bool maps = knit::yarnBodyRenderOf(s) == knit::YarnBodyRender::Maps;
            row("Yarn fibres", maps ? "The yarn is made of fibres rendered as maps (Knit step, Construction: Fibre maps; set up in the Fibres "
                                      "step, Yarn fibres): every export writes one tube per yarn under the maps baked from its fibres, in "
                                      "place of the fibres."
                                    : "The yarn is made of fibre strands (Knit step, Construction; set up in the Fibres step, Yarn fibres): "
                                      "every export writes them in place of the yarn tubes.");
            ImGui::AlignTextToFramePadding();
            mutedWrapped(maps ? "Baked into maps (Construction)" : "Included through Construction");
        }
        ImGui::TableNextRow();
        ImGui::TableNextColumn();
        ImGui::BeginDisabled(!fmt || !fmt->fibres);
        ImGui::Checkbox(fibreYarn ? "Include fuzz" : "Include fibres", &s.fibres);
        ImGui::EndDisabled();
        if (fmt && !fmt->fibres) tooltip("This format cannot carry fibres.");
        else tooltip(fibreYarn ? "The fuzz: the fibres generated in step 2 for this knit (they must be current), on top of the yarn's fibres."
                               : "The fibres generated in step 2 for this knit (they must be current).");
        ImGui::TableNextColumn();
        ImGui::AlignTextToFramePadding();
        const FibreExport fx = fibreExport(ui);
        wrapped(fx == FibreExport::Update || fx == FibreExport::Generate || fibresGeneratedExcluded(ui) ? kWarn : muted(), fibreExportWord(ui));
        ImGui::TableNextRow();
        ImGui::TableNextColumn();
        ImGui::Checkbox("Filling", &s.filling);
        tooltip("A surface behind the knit that shows through its gaps: stuffing inside closed parts, a backing behind open sheets.");
        ImGui::TableNextColumn();
        ImGui::BeginDisabled(!s.filling);
        colorEdit("Colour", s.stuffing.color);
        ImGui::EndDisabled();
        tooltip(s.filling ? "The filling's colour." : "The filling's colour (no filling is written).");
    }
    const bool fur = furShown(ui);
    {
        // The surfaces that are not knitted: plain parts, parts too small to knit, the fur's body.
        const char* sources[] = {"texture", "average", "material"};
        const char* labels[] = {"Texture", "Average colour", "Material colour"};
        const char* tips[] = {"The object's own texture, as in the Mesh view. Parts without a texture get their average colour.",
                              "One colour per part: the average of its texture or vertex colours.",
                              "The material's base colour only, without its texture (often white)."};
        const char* name = fur ? "Body colour" : "Unknitted parts";
        row(name, fur ? "The colour of the surface under the fur: it shows between the fibres."
                      : "The colour of the parts that are not knitted (plain parts, parts too small for the gauge).");
        int cur = -1;
        for (int k = 0; k < 3; ++k)
            if (s.plainColourSource == sources[k]) cur = k;
        const std::string id = std::string("##") + name;
        if (beginList(id.c_str(), cur >= 0 ? labels[cur] : s.plainColourSource.c_str())) {
            for (int k = 0; k < 3; ++k) {
                if (ImGui::Selectable(labels[k], cur == k)) s.plainColourSource = sources[k];
                tooltip(tips[k]);
            }
            ImGui::EndCombo();
        }
    }
    ImGui::EndTable();
    mutedWrapped(fur ? "Colour source of the surface under the fur."
                     : "Unknitted parts are always written, also when hidden in the view; this chooses their colour.");
}

// A level in the list: the Knit step's tessellation level it sets and its forecast.
std::string exportDetailRow(const SceneUi& ui, const ExportDetail& d) {
    const knit::FilePreset* level = meshLevelOf(d);
    return (level ? level->label + " tessellation" + kDot : std::string()) + exportDetailForecast(ui, d);
}

// Tessellation: the export's polygon detail (a level sets the yarn tessellation, shared with the
// Knit step, and the fibre segment budget). Not the knit's Draft / Final quality.
void meshDetailSection(SceneUi& ui) {
    ImGui::SeparatorText("Tessellation");
    ImGui::PushID("exportDetail");
    const ExportDetail* chosen = nullptr;
    float labelW = 0;
    for (const ExportDetail& d : kExportDetails) {
        if (!meshLevelOf(d)) continue;
        if (exportDetailSelected(ui, d)) chosen = &d;
        labelW = std::max(labelW, ImGui::CalcTextSize(d.label).x);
    }
    if (beginRows("##detailRows")) {
        row("Level", "Each level sets the Knit step's Yarn tessellation and the most fibre segments the export holds per object. "
                     "Triangles: the yarn tubes (about 2 x sides x yarn length / ring spacing). Fibre segments: as the export "
                     "writes them, the thin ends' extra vertices included. GPU: the renderer's memory for the package at 1600 x 1000 "
                     "(measured on photon_path_tracer). "
                     "Changing the level rebuilds the yarn tubes (no knitting) and the fibres must be generated again: their roots sit in "
                     "the tube.");
        // The list is as wide as a level and its forecast (at least the field).
        const ImGuiStyle& sty = ImGui::GetStyle();
        float listW = ImGui::CalcItemWidth();
        for (const ExportDetail& d : kExportDetails)
            if (meshLevelOf(d))
                listW = std::max(listW, labelW + sty.ItemSpacing.x * 2 + ImGui::CalcTextSize(exportDetailRow(ui, d).c_str()).x + sty.WindowPadding.x * 2);
        ImGui::SetNextWindowSizeConstraints(ImVec2(listW, 0), ImVec2(FLT_MAX, FLT_MAX));
        if (beginList("##Level", chosen ? chosen->label : "Custom")) {
            for (const ExportDetail& d : kExportDetails) {
                if (!meshLevelOf(d)) continue;
                if (ImGui::Selectable(d.label, chosen == &d)) applyExportDetail(ui, d);
                tooltip(d.note);
                ImGui::SameLine(labelW + sty.ItemSpacing.x * 2 + sty.WindowPadding.x);
                ImGui::TextDisabled("%s", exportDetailRow(ui, d).c_str());
            }
            ImGui::EndCombo();
        }
        ImGui::EndTable();
    }
    if (chosen) mutedWrapped("Estimate: " + exportDetailForecast(ui, *chosen));
    if (chosen && !exportDetailMemoryWarning(ui, *chosen).empty()) wrapped(kWarn, exportDetailMemoryWarning(ui, *chosen));
    mutedWrapped(exportDetailNote(ui));
    mutedWrapped("Changes tessellation, not knit quality.");
    if (ui.s.fibres || fibreResult(ui)) mutedWrapped("Changing the level means generating the fibres again.");
    ImGui::PopID();
}

// Destination: the folder the export writes into, the name of its files for a format that takes
// one, and the warning on a folder with files in it, with the rule (publish.h) in a disclosure.
void destinationSection(SceneUi& ui) {
    ImGui::SeparatorText("Destination");
    const ImGuiStyle& sty = ImGui::GetStyle();
    ImGui::AlignTextToFramePadding();
    ImGui::TextUnformatted("Output folder");
    const char* choose = "Choose folder ...";
    ImGui::SameLine(std::max(ImGui::GetCursorPosX(), ImGui::GetContentRegionMax().x - ImGui::CalcTextSize(choose).x - sty.FramePadding.x * 2));
    if (ImGui::Button(choose) && ui.platform.dialogs) ui.platform.request = PlatformRequest::ChooseExportFolder;  // then chooseExportFolder
    if (!ui.platform.dialogs) tooltip("The folder dialog is not available here.");
    ImGui::TextWrapped("%s", ui.exportTo.dir.c_str());
    tooltip(ui.exportTo.dir.c_str());
    if (exportNamesFiles(ui) && beginRows("##destinationRows")) {
        row("File name");
        char buf[256];
        std::snprintf(buf, sizeof buf, "%s", ui.exportTo.name.c_str());
        if (ImGui::InputText("##File name", buf, sizeof buf)) {
            ui.exportTo.name = buf;
            ui.exportTo.nameEdited = true;
        }
        tooltip("The name of the exported files, without an extension.");
        ImGui::EndTable();
    }
    if (exportTargetProblem(ui).empty()) {
        if (const std::string holds = exportHoldsText(ui); !holds.empty()) mutedWrapped(holds);
        if (const std::string existing = exportExistingNote(ui); !existing.empty()) wrapped(kWarn, existing);
    } else {
        wrapped(kBad, exportOutputText(ui));
    }
    if (ImGui::TreeNode("Existing files")) {
        mutedWrapped("The export writes directly into this folder, and creates it if it is missing. Files of the same name are replaced; nothing "
                     "else in the folder is touched. A cancelled or failed export changes nothing.");
        ImGui::TreePop();
    }
}

} // namespace

void exportPanel(SceneUi& ui, const knit::SceneEntry& e) {
    ImGui::PushID("export");
    exportContext(e);
    formatSection(ui);
    contentsSection(ui);
    if (!furShown(ui)) meshDetailSection(ui);
    destinationSection(ui);
    ImGui::PopID();
}

namespace {

// A running export: its progress and that nothing changes until it is done.
void exportRunningStatus(const SceneUi& ui) {
    const knit::JobProgress& p = ui.job.progress;
    const JobProgressView v = progressView(ui.job, false);
    char count[48] = "";
    if (v.fraction >= 0) std::snprintf(count, sizeof count, "%d / %d", p.done, p.total);
    jobProgressBlock(v, count);
    mutedWrapped("Elapsed " + fmtSeconds(p.seconds) + kDot + "a cancelled or failed export leaves the earlier one as it was.");
}

// No export running: its readiness with every reason, the last attempt, the folder it published.
void exportStatus(SceneUi& ui, const knit::ExportEligibility& el, const std::string& problem) {
    const bool blocked = el.state == knit::ExportEligibility::State::Blocked;
    const bool diagnostic = el.state == knit::ExportEligibility::State::Diagnostic;
    if (ui.job.active) {
        // Another job (a knit started from here, a mesh rebuild, fibres): the export waits for it.
        const JobProgressView v = progressView(ui.job, ui.results.partial != nullptr);
        glyph(ui.job.cancelRequested ? StageStatus::Cancelling : StageStatus::Running);
        ImGui::SameLine();
        ImGui::TextUnformatted(v.title.c_str());
        mutedWrapped("The export is checked again when it finishes; its progress and Cancel are in the status bar.");
    } else if (blocked) {
        glyph(StageStatus::Blocked);
        ImGui::SameLine();
        wrapped(kWarn, blockedHeading(el.fix));
    } else {
        glyph(diagnostic ? StageStatus::Diagnostic : StageStatus::Ready);
        ImGui::SameLine();
        ImGui::TextUnformatted(diagnostic ? "Diagnostic export" : "Ready to export");
        mutedWrapped(exportReadyLine(ui));
    }
    for (const std::string& r : el.reasons) mutedWrapped(r);
    if (diagnostic) mutedWrapped("report.json marks the export as diagnostic and lists the failed parts.");
    if ((blocked || diagnostic) && ui.results.knit && !ui.results.partial && !ui.results.knit->parts.empty() && ImGui::TextLink("Parts and checks (Report)"))
        ui.panels.reportOpen = true;
    if (!blocked && !problem.empty()) wrapped(kBad, "Cannot export there: " + problem + ". Choose a valid destination first.");
    if (!blocked && !exportExistingNote(ui).empty()) wrapped(kWarn, "The output folder is not empty (see Destination).");
    attemptNote(ui.results.exportAttempt, "");
    if (const std::string done = exportDoneText(ui); !done.empty()) {
        glyph(StageStatus::Ready);
        ImGui::SameLine();
        ImGui::TextWrapped("%s", done.c_str());
        if (ui.platform.fileManager && ImGui::Button("Open folder")) ui.platform.request = PlatformRequest::ShowExportFolder;
    }
}

// The one action: Cancel, waiting, the fix of a blocked export, or Export.
void exportAction(SceneUi& ui, const knit::ExportEligibility& el, bool targetOk) {
    if (runningStage(ui, knit::BuiltScene::Kind::Export)) {
        if (cancelButton(ui.job.cancelRequested ? "Cancel requested ..." : "Cancel export", ui.job.cancelRequested)) cancelJob(ui);
    } else if (ui.job.active) {
        primaryButton("Waiting for the running job", false, waitReason(ui).c_str());
    } else if (el.state == knit::ExportEligibility::State::Blocked) {
        if (primaryButton(el.fix.c_str(), true)) runExportFix(ui, el);
    } else {
        const std::string label = "Export " + formatShort(knit::findExportFormat(ui.exportTo.format));
        if (primaryButton(label.c_str(), !ui.job.active && targetOk, !targetOk ? "Choose a valid destination first." : waitReason(ui).c_str()))
            startExport(ui);
    }
}

} // namespace

// The Export footer (the step's one place for its state and action): the status block, then the
// action row, always last and with its room.
void exportFooter(SceneUi& ui) {
    ImGui::PushID("exportFooter");
    const knit::ExportEligibility el = eligibility(ui);
    const std::string problem = exportTargetProblem(ui);
    const float actionH = ImGui::GetFrameHeight() * 1.35f + ImGui::GetStyle().ItemSpacing.y;
    beginFooterStatus("##exportStatus", actionH);
    if (runningStage(ui, knit::BuiltScene::Kind::Export)) exportRunningStatus(ui);
    else exportStatus(ui, el, problem);
    ImGui::EndChild();
    exportAction(ui, el, problem.empty());
    ImGui::PopID();
}

} // namespace viewer
