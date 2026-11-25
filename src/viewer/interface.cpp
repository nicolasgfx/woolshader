// The scene viewer's frame (interface.h): the workflow strip along the top, the shown step's panel
// and footer, the view toolbar with its View and Display popovers, the labels in the viewport and its
// status line, the status bar and the Report and Log windows. The steps' panels are panel_*.cpp, the
// widgets widgets.cpp. No GL.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

#include "imgui.h"
#include "viewer/interface.h"
#include "viewer/view_content.h"
#include "jobs/scene_builder.h"

namespace viewer {

namespace {

// --- the workflow strip

// The workflow along the top of the window: the steps, each with its state glyph and word (stageInStrip: Import ->
// [Mark] -> Knit -> Fibres -> Export, or Import -> Fur -> Export; the optional Mark step says so). Wraps on a narrow window.
void workflowStrip(SceneUi& ui, float width) {
    ImGui::SetNextWindowPos(ImVec2(0, 0), ImGuiCond_Always);
    ImGui::SetNextWindowSizeConstraints(ImVec2(width, 0), ImVec2(width, FLT_MAX));
    ImGui::Begin("##workflowStrip", nullptr, ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings |
                                                 ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoBringToFrontOnFocus);
    const StageState states[kStages] = {importState(ui), markState(ui), knitState(ui), fibreState(ui), exportState(ui)};
    const ImGuiStyle& sty = ImGui::GetStyle();
    const float h = ImGui::GetFrameHeight() * 1.15f, gap = sty.ItemSpacing.x * 0.5f, glyphW = ImGui::GetFontSize() * 0.85f;
    const float right = ImGui::GetWindowPos().x + ImGui::GetWindowContentRegionMax().x;
    ImGui::AlignTextToFramePadding();
    ImGui::TextUnformatted("Knit viewer");
    for (int i = 0; i < kStages; ++i) {
        if (!stageInStrip(ui, static_cast<Stage>(i))) continue;
        ImGui::PushID(i);
        const char* name = stageLabel(ui, static_cast<Stage>(i));
        const float textW = ImGui::CalcTextSize(name).x, wordW = ImGui::CalcTextSize(states[i].word.c_str()).x;
        const float w = sty.FramePadding.x * 2 + textW + gap * 2 + glyphW + 2 + wordW;
        if (ImGui::GetItemRectMax().x + sty.ItemSpacing.x * (i == 0 ? 3.0f : 1.0f) + w <= right) ImGui::SameLine(0, sty.ItemSpacing.x * (i == 0 ? 3.0f : 1.0f));
        const bool sel = static_cast<int>(ui.stage) == i;
        const ImVec2 p0 = ImGui::GetCursorScreenPos();
        ImGui::BeginDisabled(ui.import.fileOpen && i != 0);  // an open file waits for its choices
        if (ImGui::InvisibleButton("##stage", ImVec2(w, h))) ui.stage = static_cast<Stage>(i);
        ImGui::EndDisabled();
        const bool hov = ImGui::IsItemHovered();
        ImDrawList* dl = ImGui::GetWindowDrawList();
        if (sel || hov) dl->AddRectFilled(p0, {p0.x + w, p0.y + h}, ImGui::GetColorU32(sel ? ImGuiCol_Header : ImGuiCol_HeaderHovered), sty.FrameRounding);
        if (sel) dl->AddRectFilled({p0.x, p0.y + h - std::max(2.0f, uiPx(2))}, {p0.x + w, p0.y + h}, ImGui::GetColorU32(kAccent));
        const float ty = p0.y + (h - ImGui::GetTextLineHeight()) * 0.5f;
        ImGui::SetCursorScreenPos({p0.x + sty.FramePadding.x, ty});
        ImGui::TextUnformatted(name);
        ImGui::SetCursorScreenPos({p0.x + sty.FramePadding.x + textW + gap * 2, ty});
        glyph(states[i].status, glyphW);
        ImGui::SameLine(0, 2);
        ImGui::PushStyleColor(ImGuiCol_Text, muted());
        ImGui::TextUnformatted(states[i].word.c_str());
        ImGui::PopStyleColor();
        ImGui::SetCursorScreenPos(p0);
        ImGui::Dummy({w, h});  // the step's rect again is the last item (the next step's flow test)
        if (hov && !states[i].detail.empty()) tooltip(states[i].detail.c_str());
        ImGui::PopID();
    }
    ui.layout.stripHeight = ImGui::GetWindowHeight();
    ImGui::End();
}

// --- the status bar

// Text cut to a width (with "..."), for one-line bars.
std::string fitText(const std::string& s, float width) {
    if (ImGui::CalcTextSize(s.c_str()).x <= width) return s;
    std::string t = s;
    while (!t.empty() && ImGui::CalcTextSize((t + "...").c_str()).x > width) t.pop_back();
    return t + "...";
}

// The stage panel that shows a job's progress and Cancel itself (the Knit, Fibres and Export
// footers): the status bar then names the step and the time only.
bool jobPanelShown(const SceneUi& ui) {
    switch (ui.job.kind) {
    case knit::BuiltScene::Kind::Knit: return ui.stage == Stage::Knit && entryOf(ui) && knitFooterShown(ui);
    case knit::BuiltScene::Kind::Fibres: return ui.stage == Stage::Fibres && entryOf(ui);
    case knit::BuiltScene::Kind::Export: return ui.stage == Stage::Export && entryOf(ui);
    }
    return false;
}

// The job's line on the left: its state, its step and time; a bar and Cancel only when no panel
// shows them. `right`: the width kept for the tools.
void jobStatus(SceneUi& ui, float width, float right) {
    const JobProgressView v = progressView(ui.job, ui.results.partial != nullptr);
    const ImGuiStyle& sty = ImGui::GetStyle();
    const bool own = jobPanelShown(ui);
    const bool bar = !own;
    const bool cancel = !own && !ui.job.cancelRequested;
    const float barW = uiPx(160), cancelW = ImGui::CalcTextSize("Cancel").x + sty.FramePadding.x * 2;
    glyph(ui.job.cancelRequested ? StageStatus::Cancelling : StageStatus::Running);
    ImGui::SameLine();
    const float textW = width - ImGui::GetCursorPosX() - right - (bar ? barW + sty.ItemSpacing.x : 0) - (cancel ? cancelW + sty.ItemSpacing.x : 0);
    ImGui::TextUnformatted(fitText(v.status, std::max(uiPx(80), textW)).c_str());
    if (bar) {
        ImGui::SameLine();
        char overlay[48] = "";
        if (!v.count.empty()) std::snprintf(overlay, sizeof overlay, "%s", v.count.c_str());
        else if (v.fraction >= 0 && ui.job.kind == knit::BuiltScene::Kind::Fibres) std::snprintf(overlay, sizeof overlay, "%d %%", int(v.fraction * 100));
        else if (v.fraction >= 0) std::snprintf(overlay, sizeof overlay, "%d / %d", ui.job.progress.done, ui.job.progress.total);
        ImGui::ProgressBar(v.fraction >= 0 ? v.fraction : -1.0f * float(ImGui::GetTime()), ImVec2(barW, 0), overlay);
    }
    if (cancel) {
        ImGui::SameLine();
        if (ImGui::Button("Cancel")) cancelJob(ui);
    }
}

// Idle: the last job that did not finish, on one line up to the tools (the whole message on hover
// and in the Log).
void idleStatus(const SceneUi& ui, float width, float right) {
    const Attempt* last = nullptr;
    for (const Attempt* a : {&ui.results.knitAttempt, &ui.results.fibreAttempt, &ui.results.exportAttempt})
        if (a->outcome == knit::JobOutcome::Cancelled || a->outcome == knit::JobOutcome::Failed) last = a;
    if (!last) {
        glyph(StageStatus::Ready);
        ImGui::SameLine();
        ImGui::TextUnformatted("Idle");
        return;
    }
    glyph(last->outcome == knit::JobOutcome::Cancelled ? StageStatus::Cancelled : StageStatus::Failed);
    ImGui::SameLine();
    const std::string m = last->outcome == knit::JobOutcome::Cancelled ? last->message : "Failed: " + last->message;
    ImGui::TextUnformatted(fitText(m, std::max(uiPx(80), width - ImGui::GetCursorPosX() - right)).c_str());
    tooltip(m.c_str());
}

// The status bar's buttons on the right, after the compute device.
const char* const kStatusTools[] = {"Report", "Log", "Preview", "Interface", "About"};

// The thin divider before About: its width, and the extra room it takes beyond the usual spacing.
float aboutDividerWidth() { return std::max(1.0f, uiPx(1)); }
float aboutExtraRoom() { return ImGui::GetStyle().ItemSpacing.x + aboutDividerWidth(); }

// About is not a tool, so it does not look like one: a divider sets it apart, and it is an outline in the
// accent colour with light-blue text and no fill, faintly tinted while its window is open. Hover and press
// fill it a little more, as the filled buttons brighten.
bool aboutButton(bool open) {
    const ImGuiStyle& sty = ImGui::GetStyle();
    ImGui::SameLine(0, sty.ItemSpacing.x);
    const ImVec2 p = ImGui::GetCursorScreenPos();
    const float h = ImGui::GetTextLineHeight();
    ImGui::GetWindowDrawList()->AddRectFilled(ImVec2(p.x, p.y + h * 0.15f), ImVec2(p.x + aboutDividerWidth(), p.y + h * 0.85f),
                                              ImGui::GetColorU32(ImGuiCol_Separator));
    ImGui::Dummy(ImVec2(aboutDividerWidth(), h));
    ImGui::SameLine(0, sty.ItemSpacing.x);
    const ImVec4 accent = kAccent;
    ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(accent.x, accent.y, accent.z, open ? 0.22f : 0.0f));
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(accent.x, accent.y, accent.z, 0.30f));
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(accent.x, accent.y, accent.z, 0.42f));
    ImGui::PushStyleColor(ImGuiCol_Border, ImVec4(accent.x, accent.y, accent.z, 0.75f));
    ImGui::PushStyleColor(ImGuiCol_Text, kInfo);
    ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, aboutDividerWidth());
    const bool hit = ImGui::SmallButton("About");
    ImGui::PopStyleVar();
    ImGui::PopStyleColor(5);
    return hit;
}

// The right side: the compute device, Report, Log, Preview, Interface (its size) and About.
void statusTools(SceneUi& ui, const std::string& device) {
    const knit::ComputeCaps& caps = knit::computeCaps();
    ImGui::TextDisabled("%s", device.c_str());
    tooltip(caps.device ? (caps.deviceName + (caps.ncp ? std::string() : "; the knitting solver is unavailable: " + caps.whyNoNcp)).c_str()
                        : ("No GPU computing: " + caps.whyNoDevice + ". Knitting needs a usable GPU; without one only the CPU fibre placements run.").c_str());
    ImGui::SameLine();
    if (ImGui::SmallButton("Report")) ui.panels.reportOpen = !ui.panels.reportOpen;
    ImGui::SameLine();
    if (ImGui::SmallButton("Log")) ui.panels.logOpen = !ui.panels.logOpen;
    ImGui::SameLine();
    if (ImGui::SmallButton("Preview")) ui.panels.previewOpen = !ui.panels.previewOpen;
    tooltip("A live path-traced render of the scene in a window of its own (the photon path tracer)");
    ImGui::SameLine();
    if (ImGui::SmallButton("Interface")) ImGui::OpenPopup("##display");
    if (ImGui::BeginPopup("##display")) {
        ImGui::TextUnformatted("Interface size");
        const float before = uiScaleUser();
        uiScaleSlider();
        mutedWrapped("1.00 follows the operating system's display scaling; larger or smaller multiplies it. Kept for the next session.");
        if (ui.prefs && uiScaleUser() != before) {
            ui.prefs->uiScale = uiScaleUser();
            saveViewerPrefs(*ui.prefs);
        }
        ImGui::EndPopup();
    }
    if (aboutButton(ui.panels.aboutOpen)) ui.panels.aboutOpen = !ui.panels.aboutOpen;
    tooltip("The author, the licence, the third-party software and the example models' credits");
}

void statusBar(SceneUi& ui, float width, float height) {
    const float h = ImGui::GetFrameHeight() + ImGui::GetStyle().WindowPadding.y * 2;
    ui.layout.statusHeight = h;
    ImGui::SetNextWindowPos(ImVec2(0, height - h));
    ImGui::SetNextWindowSize(ImVec2(width, h));
    ImGui::Begin("##status", nullptr, ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoScrollbar);
    const knit::ComputeCaps& caps = knit::computeCaps();
    const std::string device = caps.device ? "GPU: " + caps.deviceName.substr(0, caps.deviceName.find(" (")) : "CPU only";
    const ImGuiStyle& sty = ImGui::GetStyle();
    float right = ImGui::CalcTextSize(device.c_str()).x + aboutExtraRoom() + uiPx(8);
    for (const char* tool : kStatusTools) right += ImGui::CalcTextSize(tool).x + sty.FramePadding.x * 2 + sty.ItemSpacing.x;
    if (ui.job.active) jobStatus(ui, width, right);
    else idleStatus(ui, width, right);
    ImGui::SameLine(std::max(ImGui::GetCursorPosX() + uiPx(8), width - right));
    statusTools(ui, device);
    ImGui::End();
}

// --- the view toolbar: the stage on the left; View, Fit and Display on the right

float buttonWidth(const char* label) { return ImGui::CalcTextSize(label).x + ImGui::GetStyle().FramePadding.x * 2; }

// The toolbar's groups in rows: a group that does not fit goes to the next one.
struct ToolbarFlow {
    float right;  // the window's right edge
    void next(float groupWidth, float gap) const {
        if (ImGui::GetItemRectMax().x + gap + groupWidth <= right) ImGui::SameLine(0, gap);
    }
};

// The Fibres step's preview target: the small reference sample (recomputed after each edit, in
// about a second) or the model with its last generated result. Choosing one starts no job; the step
// always opens on the reference sample.
void previewTargetGroup(SceneUi& ui, float right) {
    const ImGuiStyle& sty = ImGui::GetStyle();
    const bool fur = furShown(ui);
    const char* reference = fur ? "Reference sphere" : "Reference patch";
    int t = ui.view.fibres == FibreView::Reference ? 0 : ui.view.fibres == FibreView::Model ? 1 : 2;
    ImGui::BeginGroup();
    // The title only where the whole group fits (a narrow viewport keeps the buttons).
    const float buttonsW = (fur ? buttonWidth("Discs") + 2 : 0) + buttonWidth(reference) + buttonWidth("Model") + 2;
    if (ImGui::CalcTextSize("Preview target").x + sty.ItemSpacing.x + buttonsW <= right - ImGui::GetCursorScreenPos().x) {
        ImGui::AlignTextToFramePadding();
        ImGui::TextDisabled("Preview target");
        ImGui::SameLine();
    }
    if (fur) {
        if (choiceButton("Discs", t == 2)) t = 2;
        tooltip("Preview target: one disc per loop fibre and one needle per fine fibre on the model, following every edit (no settle: "
                "gravity and the exact bend show after Grow fur).");
        ImGui::SameLine(0, 2);
    }
    if (choiceButton(reference, t == 0)) t = 0;
    tooltip(fur ? "Preview target: a small sample that follows every edit in about a second (the whole model is computed by Grow fur "
                  "on the model)."
                : "Preview target: a small sample that follows every edit in about a second (the whole model is computed by Apply to "
                  "model).");
    ImGui::SameLine(0, 2);
    if (choiceButton("Model", t == 1)) t = 1;
    tooltip("Preview target: the whole object with its last generated result. Each view keeps its own camera; Fit refits it.");
    ImGui::EndGroup();
    ui.view.fibres = t == 0 ? FibreView::Reference : t == 1 ? FibreView::Model : FibreView::Discs;
}

float stageSegmentWidth() { return buttonWidth("Source") + buttonWidth("Knit") + buttonWidth("Fibres") + 4; }

// Source, Knit, Fibres: one segment. A stage with nothing to draw is disabled and says why; the marked one
// is the stage on show (markedViewMode). Only a viewport narrower than the segment itself breaks it into
// rows, so that every stage stays reachable.
void stageSegment(SceneUi& ui, const ToolbarFlow& flow) {
    const char* names[3] = {"Source", "Knit", "Fibres"};
    const char* about[3] = {"The object as imported, with its colour.", "The knitted yarn.",
                            !furShown(ui) ? "The knitted yarn with its fibres." : "The mesh with its fur."};
    const ViewMode marked = markedViewMode(ui);
    ImGui::BeginGroup();
    for (int i = 0; i < 3; ++i) {
        const ViewMode mode = static_cast<ViewMode>(i);
        const ViewModeState state = viewModeState(ui, mode);
        if (i > 0) flow.next(buttonWidth(names[i]), 2);
        ImGui::BeginDisabled(!state.available);
        if (choiceButton(names[i], mode == marked)) ui.view.mode = mode;
        ImGui::EndDisabled();
        tooltip(state.available ? about[i] : state.reason.c_str());
    }
    ImGui::EndGroup();
}

// The Mark step's views in place of Source, Knit, Fibres (which show nothing different there): Markings, Rows and the
// Stitch Mesh view, the same choice as setMarkView. Stitch Mesh is open as soon as the file is read: before a build it
// says so; Rows needs a build.
void markViewSegment(SceneUi& ui, const ToolbarFlow& flow) {
    const char* names[3] = {"Markings", "Rows", "Stitch Mesh"};
    const MarkView views[3] = {MarkView::Markings, MarkView::Rows, MarkView::StitchMesh};
    const bool read = ui.mark.session.snapshot().pieces != nullptr;
    const bool built = shownBuild(ui) != nullptr;
    const bool available[3] = {read, read && built, read};
    const char* about[3] = {"The markings on the piece.", "The row contours of the isolated piece's last build.",
                            "The generated stitch polygons on the piece: regular stitches blue, increases green, decreases orange; the findings of "
                            "the checks next to them."};
    const char* why = !read ? "Import a file first (Import)." : "Build the isolated piece first (Build stitch mesh).";
    ImGui::BeginGroup();
    for (int i = 0; i < 3; ++i) {
        if (i > 0) flow.next(buttonWidth(names[i]), 2);
        ImGui::BeginDisabled(!available[i]);
        if (choiceButton(names[i], ui.mark.view == views[i])) setMarkView(ui, views[i]);
        ImGui::EndDisabled();
        tooltip(available[i] ? about[i] : why);
    }
    ImGui::EndGroup();
}

const char* cameraViewLabel(CameraView view) {
    switch (view) {
    case CameraView::ThreeQuarter: return "3/4 view";
    case CameraView::Front: return "Front";
    case CameraView::Side: return "Side";
    case CameraView::Top: return "Top";
    case CameraView::Back: return "Back";
    }
    return "";
}

// A popover of the toolbar, placed inside the viewport (popoverRect); Escape closes it as a click outside does.
bool beginToolbarPopover(const char* id, const ViewportRect& viewport, float anchorRight) {
    const PopoverRect p = popoverRect(viewport, anchorRight, uiPx(246), uiPx(5));
    ImGui::SetNextWindowPos(ImVec2(p.x, p.y), ImGuiCond_Always);
    ImGui::SetNextWindowSizeConstraints(ImVec2(p.width, 0), ImVec2(p.width, p.maxHeight));
    if (!ImGui::BeginPopup(id)) return false;
    if (ImGui::IsKeyPressed(ImGuiKey_Escape)) ImGui::CloseCurrentPopup();
    return true;
}

// The camera presets; choosing one closes the popover.
void viewPopover(SceneUi& ui, const ViewportRect& viewport, float anchorRight) {
    if (!beginToolbarPopover("##viewMenu", viewport, anchorRight)) return;
    ImGui::TextDisabled("Camera view");
    const CameraView views[5] = {CameraView::ThreeQuarter, CameraView::Front, CameraView::Side, CameraView::Top, CameraView::Back};
    const float cell = (ImGui::GetContentRegionAvail().x - ImGui::GetStyle().ItemSpacing.x) * 0.5f;
    for (int i = 0; i < 5; ++i) {
        if (i % 2 == 1) ImGui::SameLine();
        if (choiceButton(cameraViewLabel(views[i]), ui.view.cameraPreset == views[i], cell)) {
            ui.camera.view = views[i];
            ImGui::CloseCurrentPopup();
        }
    }
    ImGui::Separator();
    if (ImGui::Button("Upright", ImVec2(-FLT_MIN, 0))) {
        ui.camera.upright = true;
        ImGui::CloseCurrentPopup();
    }
    tooltip("Level the view: look horizontally at the object, facing the same way.");
    ImGui::EndPopup();
}

// Lighting and what is drawn; the popover stays open while these change. A row that does not apply to the
// stage on show is hidden and keeps its value (displayRows).
void displayPopover(SceneUi& ui, const ViewportRect& viewport, float anchorRight) {
    if (!beginToolbarPopover("##displayMenu", viewport, anchorRight)) return;
    const DisplayRows rows = displayRows(ui);
    ImGui::TextDisabled("Lighting");
    ImGui::BeginDisabled(!ui.view.studioAvailable);
    if (choiceButton("Studio", ui.view.lighting == Lighting::Studio)) ui.view.lighting = Lighting::Studio;
    ImGui::EndDisabled();
    tooltip(ui.view.studioAvailable ? "Light from the upper left with ambient occlusion: the grooves between the yarns darken, so the stitches "
                                      "read in depth."
                                    : ("Studio light is unavailable on this graphics driver: " + ui.view.studioNote).c_str());
    ImGui::SameLine(0, 2);
    if (choiceButton("Flat", ui.view.lighting == Lighting::Flat)) ui.view.lighting = Lighting::Flat;
    tooltip("A plain headlight from the camera, without occlusion.");
    if (rows.occlusion) {
        ImGui::BeginDisabled(!ui.view.studioAvailable || ui.view.wire);
        ImGui::SetNextItemWidth(ImGui::GetContentRegionAvail().x - ImGui::CalcTextSize("Occlusion").x - ImGui::GetStyle().ItemInnerSpacing.x);
        ImGui::SliderFloat("Occlusion", &ui.view.occlusion, 0.0f, 1.0f, "%.2f");
        ImGui::EndDisabled();
        tooltip(ui.view.wire ? "Wireframe draws with the Flat light."
                             : "How strongly the ambient occlusion darkens the grooves between yarns and the contact between parts.");
    }
    ImGui::Separator();
    ImGui::Checkbox("Wireframe", &ui.view.wire);
    tooltip("The triangles' edges (drawn with the Flat light).");
    if (rows.plyTexture) {
        ImGui::BeginDisabled(!rows.plyTextureWhyNot.empty());
        ImGui::Checkbox("Ply texture", &ui.view.plyTexture);
        ImGui::EndDisabled();
        tooltip(!rows.plyTextureWhyNot.empty() ? rows.plyTextureWhyNot.c_str()
                                               : "Shows the twisted-ply maps' ambient occlusion on the yarn (the renderer also uses the normal and "
                                                 "height maps). Hiding it here does not change the ply texture (Knit step, Appearance).");
    }
    if (rows.yarnMaps) {
        ImGui::Checkbox("Yarn maps", &ui.view.yarnMaps);
        tooltip("The baked cavity maps on the yarn tubes (the renderer also uses their normal, specular and roughness maps). Display only.");
    }
    if (rows.unknittedParts) {
        ImGui::Checkbox("Unknitted parts", &ui.view.unknittedParts);
        tooltip("Parts that are not knitted (too small for the stitch size, or not knittable) as plain surfaces.");
    }
    if (rows.fibrePopulations) {
        ImGui::Checkbox("Fibre populations", &ui.view.fibrePopulations);
        tooltip("Tells the fibres apart on the model: yellow along the yarn, blue stray, as on the reference patch. Display only: "
                "the fibres' colours and the export do not change.");
    }
    if (rows.compare) {
        // Ticked while it cannot be drawn: suspended, and it can still be unticked.
        ImGui::BeginDisabled(!ui.view.compare && !rows.compareWhyNot.empty());
        ImGui::Checkbox("Compare with source", &ui.view.compare);
        ImGui::EndDisabled();
        tooltip(rows.compareWhyNot.empty() ? "Source in the left half, this stage in the right half: two halves of one view, one camera."
                                           : rows.compareWhyNot.c_str());
        if (ui.view.compare && !rows.compareWhyNot.empty()) mutedWrapped("Suspended: " + rows.compareWhyNot);
    }
    ImGui::EndPopup();
}

// View, Fit and Display: one group at the row's right end, or the next row's start when it does not fit;
// only a viewport narrower than the group itself breaks it. The popovers come after the group, placed
// by their buttons' right edges.
void rightGroup(SceneUi& ui, const ToolbarFlow& flow, const ViewportRect& viewport) {
    const ImGuiStyle& sty = ImGui::GetStyle();
    const char* view = ui.view.cameraPreset ? cameraViewLabel(*ui.view.cameraPreset) : "Orbit";
    const float width = menuButtonWidth(view) + buttonWidth("Fit") + menuButtonWidth("Display") + sty.ItemSpacing.x * 2;
    const float gap = ImGui::GetFontSize() * 1.5f;
    if (ImGui::GetItemRectMax().x + gap + width <= flow.right) ImGui::SameLine(flow.right - width - ImGui::GetWindowPos().x);
    ImGui::BeginGroup();
    if (menuButton(view, "viewMenu", ImGui::IsPopupOpen("##viewMenu"))) ImGui::OpenPopup("##viewMenu");
    tooltip("The camera's view: a preset, or Orbit once the view was turned by dragging.");
    const float viewRight = ImGui::GetItemRectMax().x;
    flow.next(buttonWidth("Fit"), sty.ItemSpacing.x);
    ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 1.0f);
    if (ImGui::Button("Fit")) ui.camera.frame = true;
    ImGui::PopStyleVar();
    tooltip("Frame the visible model (F). Drag: left turns around the vertical and tilts (Shift: one of the two only), right / middle "
            "pans; the wheel zooms.");
    flow.next(menuButtonWidth("Display"), sty.ItemSpacing.x);
    if (menuButton("Display", "displayMenu", ImGui::IsPopupOpen("##displayMenu"))) ImGui::OpenPopup("##displayMenu");
    tooltip("Lighting, wireframe, the yarn's and the parts' display, and the comparison with the source.");
    const float displayRight = ImGui::GetItemRectMax().x;
    ImGui::EndGroup();
    viewPopover(ui, viewport, viewRight);
    displayPopover(ui, viewport, displayRight);
}

// F frames the visible model, as Fit does, unless a field is being typed in.
void fitShortcut(SceneUi& ui) {
    const ImGuiIO& io = ImGui::GetIO();
    if (!io.WantTextInput && !io.KeyCtrl && !io.KeyAlt && !io.KeySuper && ImGui::IsKeyPressed(ImGuiKey_F, false)) ui.camera.frame = true;
}

// The Mark step's keys, unless a field is being typed in: Ctrl+Z undo, Ctrl+Y (or Ctrl+Shift+Z) redo, Delete the
// selected point.
void markShortcuts(SceneUi& ui) {
    const ImGuiIO& io = ImGui::GetIO();
    if (ui.stage != Stage::Mark || io.WantTextInput) return;
    if (io.KeyCtrl && ImGui::IsKeyPressed(ImGuiKey_Z, false)) {
        if (io.KeyShift) redoMarking(ui);
        else undoMarking(ui);
    }
    if (io.KeyCtrl && ImGui::IsKeyPressed(ImGuiKey_Y, false)) redoMarking(ui);
    if (ImGui::IsKeyPressed(ImGuiKey_Delete, false)) deleteSelectedPoint(ui);
}

void viewToolbar(SceneUi& ui, float width, float height) {
    const float x = ui.layout.panelWidth, w = std::max(1.0f, width - x);
    ImGui::SetNextWindowPos(ImVec2(x, ui.layout.stripHeight), ImGuiCond_Always);
    ImGui::SetNextWindowSizeConstraints(ImVec2(w, 0), ImVec2(w, FLT_MAX));
    ImGui::Begin("##viewToolbar", nullptr, ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings |
                                               ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoFocusOnAppearing);
    const ToolbarFlow flow{ImGui::GetWindowPos().x + ImGui::GetWindowContentRegionMax().x};
    const float gap = ImGui::GetFontSize() * 1.5f;
    if (ui.stage != Stage::Fibres) ui.view.fibres = FibreView::Reference;
    if (ui.stage == Stage::Fibres && entryOf(ui)) {
        previewTargetGroup(ui, flow.right);
        flow.next(stageSegmentWidth(), gap);
    }
    if (ui.stage == Stage::Mark) markViewSegment(ui, flow);
    else stageSegment(ui, flow);
    rightGroup(ui, flow, viewportRect(ui, width, height));  // the viewport of the last frame's layout
    fitShortcut(ui);
    markShortcuts(ui);
    ui.layout.toolbarHeight = ImGui::GetWindowHeight();
    ImGui::End();
}

// --- the viewport's status line, labels, Report and Log

// The line along the viewport's bottom edge: the object's name and size on the left, the view's state on the
// right (cut to its room; the whole on hover). A partial build is said in the warning colour.
void viewStatusLine(SceneUi& ui, float width, float height) {
    const float x = ui.layout.panelWidth, w = std::max(1.0f, width - x), h = ui.layout.viewStatusHeight;
    ImGui::SetNextWindowPos(ImVec2(x, height - ui.layout.statusHeight - h));
    ImGui::SetNextWindowSize(ImVec2(w, h));
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(uiPx(10), uiPx(3)));
    ImGui::Begin("##viewStatus", nullptr, ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings |
                                              ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoNav);
    ImGui::PopStyleVar();
    const ObjectLine object = objectLine(ui);
    const ViewContent c = viewContent(ui);
    if (!object.name.empty()) {
        ImGui::TextUnformatted(object.name.c_str());
        ImGui::SameLine(0, 0);
    }
    if (!object.size.empty()) {
        ImGui::PushStyleColor(ImGuiCol_Text, muted());
        ImGui::TextUnformatted(((object.name.empty() ? std::string() : std::string(kDot)) + object.size).c_str());
        ImGui::PopStyleColor();
        ImGui::SameLine(0, 0);
    }
    const float right = ImGui::GetWindowContentRegionMax().x, room = right - ImGui::GetCursorPosX() - uiPx(24);
    if (!c.caption.empty() && room > uiPx(40)) {
        const std::string shown = fitText(c.caption, room);
        ImGui::SameLine(right - ImGui::CalcTextSize(shown.c_str()).x);
        const bool partial = ui.results.partial && c.meshOwner == ui.results.partial.get();
        ImGui::PushStyleColor(ImGuiCol_Text, partial ? kWarn : muted());
        ImGui::TextUnformatted(shown.c_str());
        ImGui::PopStyleColor();
        if (shown != c.caption) tooltip(c.caption.c_str());
    } else {
        ImGui::NewLine();
    }
    ImGui::End();
}

// Passive labels in the viewport (they take no input: the camera drags through them): the triangle counts,
// bottom left, and while comparing the two halves' names over a divider at the viewport's middle (the
// renderer splits there too).
void viewportLabels(SceneUi& ui, float width, float height) {
    const ViewportRect r = viewportRect(ui, width, height);
    const float pad = uiPx(10), wrap = std::max(uiPx(160), r.w - 4 * pad);
    const ImGuiWindowFlags passive = ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoSavedSettings |
                                     ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoInputs | ImGuiWindowFlags_NoNav;
    const ViewportFrame frame = viewportFrame(ui);
    if (frame.source) {
        const float middle = r.x + std::floor(r.w * 0.5f);
        ImGui::GetBackgroundDrawList()->AddLine(ImVec2(middle, r.y), ImVec2(middle, r.y + r.h), ImGui::GetColorU32(muted()), std::max(1.0f, uiPx(1)));
        const char* stage = markedViewMode(ui) == ViewMode::Fibres ? "Fibres" : "Knit";
        const char* names[2] = {"Source", stage};
        const char* ids[2] = {"##compareSource", "##compareStage"};
        for (int i = 0; i < 2; ++i) {
            ImGui::SetNextWindowPos(ImVec2(r.x + r.w * (i == 0 ? 0.25f : 0.75f), r.y + pad), ImGuiCond_Always, ImVec2(0.5f, 0));
            ImGui::SetNextWindowBgAlpha(0.7f);
            ImGui::Begin(ids[i], nullptr, passive);
            ImGui::PushStyleColor(ImGuiCol_Text, muted());
            ImGui::TextUnformatted(names[i]);
            ImGui::PopStyleColor();
            ImGui::End();
        }
    }
    const ViewContent& c = frame.shown;
    const float line = ImGui::GetTextLineHeightWithSpacing() + ImGui::GetStyle().WindowPadding.y * 2;
    if (r.h < line * 4 + 3 * pad || (c.meshes.empty() && !entryOf(ui))) return;
    ImGui::SetNextWindowPos(ImVec2(r.x + pad, r.y + r.h - pad), ImGuiCond_Always, ImVec2(0, 1));
    ImGui::SetNextWindowSizeConstraints(ImVec2(0, 0), ImVec2(std::max(1.0f, r.w - 2 * pad), r.h * 0.5f - 2 * pad));
    ImGui::SetNextWindowBgAlpha(0.5f);
    ImGui::Begin("##viewCounts", nullptr, passive);
    ImGui::PushTextWrapPos(wrap);
    ImGui::PushStyleColor(ImGuiCol_Text, muted());
    ImGui::TextUnformatted(viewCounts(ui, c).c_str());
    ImGui::PopStyleColor();
    ImGui::PopTextWrapPos();
    ImGui::End();
}

// The preview image's "updating" cover: the image dimmed, a turning arc and what the preview waits for.
void previewBusyCover(const ImVec2& a, const ImVec2& b, const std::string& text) {
    ImDrawList* draw = ImGui::GetWindowDrawList();
    draw->AddRectFilled(a, b, IM_COL32(0, 0, 0, 120));
    const float radius = uiPx(9);
    const float wrap = std::max(1.0f, b.x - a.x - 4 * radius);
    const ImVec2 textSize = ImGui::CalcTextSize(text.c_str(), nullptr, false, wrap);
    const ImVec2 middle((a.x + b.x) * 0.5f, (a.y + b.y) * 0.5f);
    const float turn = static_cast<float>(ImGui::GetTime()) * 4;
    draw->PathArcTo(ImVec2(middle.x, middle.y - textSize.y * 0.5f - radius * 1.5f), radius, turn, turn + 4.2f, 24);
    draw->PathStroke(IM_COL32(255, 255, 255, 220), 0, uiPx(2.5f));
    draw->AddText(nullptr, 0, ImVec2(middle.x - textSize.x * 0.5f, middle.y - textSize.y * 0.5f + radius * 0.5f), IM_COL32(255, 255, 255, 235),
                  text.c_str(), nullptr, wrap);
}

// The Render preview window: the renderer's newest frame, fitted to the window, and one status line.
void renderPreviewWindow(SceneUi& ui, float width, float height) {
    // First opened: in the viewport's top right corner, over the 3D view and never under a panel.
    const ViewportRect view = viewportRect(ui, width, height);
    const float pad = uiPx(10);
    ImGui::SetNextWindowPos(ImVec2(view.x + view.w - pad, view.y + pad), ImGuiCond_FirstUseEver, ImVec2(1, 0));
    ImGui::SetNextWindowSize(ImVec2(std::min(uiPx(540), view.w - 2 * pad), std::min(uiPx(400), view.h - 2 * pad)), ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Render preview", &ui.panels.previewOpen)) {
        // Wholly inside the main window: a place kept from a larger screen would leave only the title in view.
        const ImVec2 display = ImGui::GetIO().DisplaySize, windowAt = ImGui::GetWindowPos(), extent = ImGui::GetWindowSize();
        const ImVec2 inside(std::clamp(windowAt.x, 0.0f, std::max(0.0f, display.x - extent.x)), std::clamp(windowAt.y, 0.0f, std::max(0.0f, display.y - extent.y)));
        if ((inside.x != windowAt.x || inside.y != windowAt.y) && !ImGui::IsMouseDragging(ImGuiMouseButton_Left)) ImGui::SetWindowPos(inside);
        RenderPreviewView& v = ui.renderPreview;
        const knit::RenderPreviewState& s = v.state;
        const float footer = ImGui::GetFrameHeightWithSpacing() * 2;
        const ImVec2 area(std::max(1.0f, ImGui::GetContentRegionAvail().x), std::max(1.0f, ImGui::GetContentRegionAvail().y - footer));
        v.imageW = area.x;
        v.imageH = area.y;
        if (v.texture && v.textureW > 0 && v.textureH > 0) {
            const float k = std::min(area.x / v.textureW, area.y / v.textureH);
            const ImVec2 size(v.textureW * k, v.textureH * k);
            ImGui::SetCursorPosX(ImGui::GetCursorPosX() + (area.x - size.x) * 0.5f);
            ImGui::Image(static_cast<ImTextureID>(static_cast<intptr_t>(v.texture)), size);
            if (!v.busy.empty()) previewBusyCover(ImGui::GetItemRectMin(), ImGui::GetItemRectMax(), v.busy);
            ImGui::Dummy(ImVec2(1, std::max(0.0f, area.y - size.y)));
        } else {
            // Nothing rendered yet: what the preview waits for, in the middle of the image area.
            const std::string why = !v.busy.empty() ? v.busy : !v.hint.empty() ? v.hint : s.phase == "off" ? "starting the renderer" : s.message;
            const ImVec2 at = ImGui::GetCursorPos();
            ImGui::PushTextWrapPos(at.x + area.x);
            ImGui::SetCursorPos(ImVec2(at.x, at.y + area.y * 0.45f));
            ImGui::TextDisabled("%s", why.c_str());
            ImGui::PopTextWrapPos();
            ImGui::SetCursorPos(at);
            ImGui::Dummy(area);
        }
        std::string line = s.phase;
        if (s.phase == "rendering" || s.phase == "converged") {
            char b[96];
            std::snprintf(b, sizeof b, " - %d spp, %.1f ms per sample", s.spp, s.msPerSpp);
            line += b;
        }
        if (!s.message.empty() && s.phase != "rendering") line += ": " + s.message;
        ImGui::TextUnformatted(fitText(line, ImGui::GetContentRegionAvail().x).c_str());
        if (!s.message.empty()) tooltip(s.message.c_str());
        ImGui::Checkbox("Pause", &v.paused);
        ImGui::SameLine();
        ImGui::SetNextItemWidth(uiPx(90));
        const char* sizes[] = {"25 %", "50 %", "100 %"};
        const int percents[] = {25, 50, 100};
        int pick = v.scalePercent == 25 ? 0 : v.scalePercent == 100 ? 2 : 1;
        if (ImGui::Combo("Resolution", &pick, sizes, 3)) v.scalePercent = percents[pick];
        if (!v.hint.empty() && v.texture) {
            ImGui::SameLine();
            ImGui::TextDisabled("(note)");
            tooltip(v.hint.c_str());
        }
    }
    ImGui::End();
}

void popOuts(SceneUi& ui, float width, float height) {
    if (ui.panels.previewOpen) renderPreviewWindow(ui, width, height);
    if (ui.panels.aboutOpen) aboutWindow(ui);
    if (ui.panels.reportOpen) {
        ImGui::SetNextWindowSize(ImVec2(uiPx(620), uiPx(420)), ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Report", &ui.panels.reportOpen)) {
            if (ui.results.knit) {
                if (ui.results.knit->knit) ImGui::Text("Knitted result #%d", ui.results.knit->knit->id);
                partsTable(ui, *ui.results.knit);
                if (ImGui::TreeNode("Validation report")) {
                    ImGui::TextUnformatted(ui.results.knit->report.c_str());
                    ImGui::TreePop();
                }
            } else {
                ImGui::TextDisabled("No knitted result yet.");
            }
        }
        ImGui::End();
    }
    if (ui.panels.logOpen) {
        ImGui::SetNextWindowSize(ImVec2(uiPx(620), uiPx(360)), ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Log", &ui.panels.logOpen)) {
            for (const Attempt* a : {&ui.results.knitAttempt, &ui.results.fibreAttempt, &ui.results.exportAttempt}) attemptNote(*a, "");
            if (ui.results.fibres && ui.results.fibres->fibreResult) ImGui::TextUnformatted(ui.results.fibres->fibreResult->log.c_str());
            if (ui.results.exported) ImGui::TextUnformatted(ui.results.exported->exportLog.c_str());
        }
        ImGui::End();
    }
}

// --- the step on show

// The step's panel, or what is missing for it.
void stepPanel(SceneUi& ui, const knit::SceneEntry* e) {
    if (ui.stage == Stage::Import) {
        importPanel(ui);
        return;
    }
    if (ui.stage == Stage::Mark) {  // it says itself what it needs (an imported file)
        markPanel(ui);
        return;
    }
    if (!e) {
        mutedWrapped("Import an object first (Import).");
        return;
    }
    switch (ui.stage) {
    case Stage::Mark: break;
    case Stage::Knit: knitPanel(ui, *e); break;
    case Stage::Fibres: fibrePanel(ui, *e); break;
    case Stage::Export: exportPanel(ui, *e); break;
    case Stage::Import: break;
    }
}

// Whether the shown step keeps a footer below its scrolling settings.
bool stepHasFooter(const SceneUi& ui, const knit::SceneEntry* e) {
    if (ui.stage == Stage::Import) return importFooterShown(ui);  // there may be no object yet
    if (ui.stage == Stage::Mark) return markFooterShown(ui);
    return e && (ui.stage == Stage::Fibres || ui.stage == Stage::Export || (ui.stage == Stage::Knit && knitFooterShown(ui)));
}

void stepFooter(SceneUi& ui, const knit::SceneEntry* e) {
    switch (ui.stage) {
    case Stage::Mark: markFooter(ui); break;
    case Stage::Knit: knitFooter(ui, *e); break;
    case Stage::Fibres: fibreFooter(ui); break;
    case Stage::Export: exportFooter(ui); break;
    case Stage::Import: importFooter(ui); break;
    }
}

} // namespace

void drawSceneUi(SceneUi& ui, float width, float height) {
    // The layout: the workflow strip along the top; the panel left, between the strip and the
    // status bar; the view toolbar above the viewport and its status line below it (viewportRect).
    ui.layout.statusHeight = ImGui::GetFrameHeight() + ImGui::GetStyle().WindowPadding.y * 2;
    ui.layout.viewStatusHeight = ImGui::GetTextLineHeight() + uiPx(3) * 2;
    workflowStrip(ui, width);
    const float top = ui.layout.stripHeight, panelH = std::max(ImGui::GetFrameHeight() * 4, height - top - ui.layout.statusHeight);
    ImGui::SetNextWindowPos(ImVec2(0, top), ImGuiCond_Always);
    // The width scales with the interface (a width kept in imgui.ini at another scale is widened);
    // at most 60 % of the window.
    ImGui::SetNextWindowSize(ImVec2(uiPx(480), panelH), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSizeConstraints(ImVec2(std::min(uiPx(440), width * 0.6f), panelH), ImVec2(std::max(uiPx(440), width * 0.6f), panelH));
    // "Workflow": its own imgui.ini entry (the width the user drags to is kept between sessions).
    // Never a scrollbar of its own: the stage panel scrolls. A footer one frame taller than the space
    // left for it would otherwise show one for that frame, narrowing (and re-wrapping) everything.
    ImGui::Begin("Workflow", nullptr, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoTitleBar |
                                          ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);
    ui.layout.panelWidth = ImGui::GetWindowWidth();
    const knit::SceneEntry* e = entryOf(ui);
    // A footer keeps the height it took last frame below the scrolling settings.
    const bool footer = stepHasFooter(ui, e);
    float& footerHeight = ui.layout.footerHeight[static_cast<int>(ui.stage)];
    if (footer && footerHeight <= 0) footerHeight = ImGui::GetWindowHeight() * 0.4f;  // first frame: its budget
    ImGui::BeginChild("##stagePanel", ImVec2(0, footer ? -footerHeight : 0.0f), ImGuiChildFlags_None);
    ImGui::PushItemWidth(std::max(ImGui::GetContentRegionAvail().x * 0.45f, ImGui::GetContentRegionAvail().x - ImGui::GetFontSize() * 10.5f));
    stepPanel(ui, e);
    ImGui::PopItemWidth();
    ImGui::EndChild();
    if (footer) {
        const float y0 = ImGui::GetCursorPosY();
        stepFooter(ui, e);
        footerHeight = ImGui::GetCursorPosY() - y0;
    }
    ImGui::End();
    viewToolbar(ui, width, height);
    viewportLabels(ui, width, height);
    viewStatusLine(ui, width, height);
    statusBar(ui, width, height);
    popOuts(ui, width, height);
}

// The frame readout in the viewport's top right corner (F3).
void frameStatsLabel(const SceneUi& ui, float width, float height, double frameMs, long long segments) {
    const ViewportRect r = viewportRect(ui, width, height);
    char b[160];
    std::snprintf(b, sizeof b, "%.1f ms / frame, %lld fibre segments drawn", frameMs, segments);
    const ImVec2 size = ImGui::CalcTextSize(b);
    ImGui::GetForegroundDrawList()->AddText(ImVec2(r.x + r.w - size.x - 8, r.y + 8), IM_COL32(255, 255, 160, 255), b);
}

} // namespace viewer
