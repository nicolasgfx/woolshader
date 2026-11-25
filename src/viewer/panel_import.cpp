// The Import step (interface.h), in two stages. Stage 1: the object so far, Open mesh, the examples, the method.
// Stage 2: the open file on the live 3D view, in five sections (what to make, its components, the real-world size,
// the surface preparation, colour and compatibility), and a footer with the import's summary and its action.
#include "model/coordinates.h"
#include "viewer/interface.h"
#include "settings/scene_registry.h"

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <string>


namespace fs = std::filesystem;

namespace viewer {

namespace {

constexpr size_t kNodeRowsShown = 8;  // the component list's height: it scrolls, the sections below stay in reach

// The x where a widget of this width ends flush with the right edge (or right here, when there is no room).
float rightAligned(float width) { return std::max(ImGui::GetCursorPosX(), ImGui::GetContentRegionMax().x - width); }

// --- what to make

// Import > What to make: knitted (the surface is rebuilt), or hair on the file's own mesh. Only that mesh is
// repaired, so the repairs and the warnings about it show under Fur.
void importRepairs(SceneUi& ui, const knit::ImportAnalysis& an) {
    for (const knit::ImportIssue& is : an.issues) {
        ImGui::PushID(is.id.c_str());
        if (is.kind == knit::ImportIssue::Kind::Repairable) {
            bool* flag = is.id == "degenerate" ? &ui.import.options.removeDegenerate
                         : is.id == "duplicate_faces" ? &ui.import.options.removeDuplicateFaces
                         : is.id == "weld" ? &ui.import.options.weld
                         : is.id == "unreferenced" ? &ui.import.options.removeUnreferenced
                         : is.id == "orientation" ? &ui.import.options.fixOrientation
                                                  : nullptr;
            if (flag) ImGui::Checkbox(is.action.c_str(), flag);
            ImGui::Indent();
            mutedWrapped(is.text);
            ImGui::Unindent();
        } else if (is.kind == knit::ImportIssue::Kind::Warning && is.id != "missing_file") {
            glyph(StageStatus::NeedsUpdate);
            ImGui::SameLine();
            ImGui::TextWrapped("%s %s", is.text.c_str(), is.action.c_str());
        }
        ImGui::PopID();
    }
}

// Knit or Fur for the object: the Import step's first choice, for an open file and for an example. Either can be
// chosen again later; the other workflow's settings and results are kept.
void workflowChoice(SceneUi& ui) {
    if (!knit::furEnabled()) return;  // the Fur switch is off: Knit is the only workflow
    ImGui::SeparatorText("What to make");
    const bool fur = knit::furWorkflow(ui.s);
    ImGui::BeginDisabled(ui.job.active);  // a running job's result belongs to the workflow it started in
    if (choiceButton("Knit", ui.workflowChosen && !fur)) chooseWorkflow(ui, knit::Workflow::Knit);
    tooltip("Import -> Knit -> Fibres -> Export: the object is knitted, and wool fibres grow on the yarn.");
    ImGui::SameLine();
    if (choiceButton("Fur", ui.workflowChosen && fur)) chooseWorkflow(ui, knit::Workflow::Fur);
    tooltip("Import -> Fur -> Export: a pile of curled fibres grows on the object's own mesh; nothing is knitted.");
    ImGui::EndDisabled();
    if (ui.job.active) mutedWrapped(waitReason(ui));
    if (!ui.workflowChosen) wrapped(kWarn, std::string(kChooseWorkflow) + " for this object.");
}

void purposeSection(SceneUi& ui, const knit::ImportAnalysis& an) {
    workflowChoice(ui);
    if (!ui.workflowChosen || !knit::furWorkflow(ui.s)) return;
    mutedWrapped("Fur grows on the file's own mesh: no knitting charts are prepared. Unticked components stay plain meshes (eyes, "
                 "a nose): no fur grows on them, and they block the fur. Knit later rebuilds the surface then.");
    bool anyRepair = false;
    for (const knit::ImportIssue& is : an.issues)
        anyRepair |= is.kind == knit::ImportIssue::Kind::Repairable || (is.kind == knit::ImportIssue::Kind::Warning && is.id != "missing_file");
    if (!anyRepair) return;
    importRepairs(ui, an);
    mutedWrapped("Repairs change the prepared copy only; the original file is never written.");
}

// --- components

// One component: its name (a click focuses it in the preview) and its triangles; on the right, whether it is
// knitted or grows fur (the checkbox), or keeps its mesh. Ids are per name, which the
// file makes unique (mesh_read.h). True when the row was scrolled into view in the list for a pick (only
// while the list is on screen: ImGui lays out nothing in a hidden child, so its rows have no position).
bool componentRow(SceneUi& ui, const knit::ImportNode& node, bool listOnScreen) {
    ImGui::PushID(node.name.c_str());
    const bool focused = node.name == ui.import.focusedNode;
    const bool knitted = !nodeKeptAsMesh(ui.import, node.name);
    const ImGuiStyle& st = ImGui::GetStyle();
    const float right = ImGui::CalcTextSize("Keep mesh").x + ImGui::GetFrameHeight() + st.ItemInnerSpacing.x + st.ItemSpacing.x;
    const float x0 = ImGui::GetCursorPosX();
    const float nameWidth = std::max(ImGui::GetContentRegionAvail().x - right, ImGui::GetFontSize() * 4);
    const ImGuiID focusId = ImGui::GetID("##focus");
    if (ImGui::Selectable("##focus", focused, ImGuiSelectableFlags_None, ImVec2(nameWidth, ImGui::GetFrameHeight()))) focusNode(ui, node.name);
    testLabel(focusId, node.name.c_str());
    if (focused) {  // the accent bar at the row's left edge
        const ImVec2 lo = ImGui::GetItemRectMin(), hi = ImGui::GetItemRectMax();
        ImGui::GetWindowDrawList()->AddRectFilled(lo, ImVec2(lo.x + uiPx(2), hi.y), ImGui::GetColorU32(kAccent));
    }
    ImGui::SameLine(x0 + st.FramePadding.x);
    ImGui::AlignTextToFramePadding();
    ImGui::TextUnformatted(node.name.c_str());
    ImGui::SameLine();
    ImGui::TextDisabled("%s triangles", fmtCount(node.triangles).c_str());
    ImGui::SameLine(rightAligned(right - st.ItemSpacing.x));
    bool ticked = knitted;
    const ImGuiID tickId = ImGui::GetID("###knit");
    const char* label = !ticked ? "Keep mesh###knit" : knit::furWorkflow(ui.s) && ui.workflowChosen ? "Fur###knit" : "Knit###knit";
    if (ImGui::Checkbox(label, &ticked)) toggleNode(ui, node.name);
    testLabel(tickId, ("Knit " + node.name).c_str());
    const bool reveal = focused && ui.import.revealFocus && listOnScreen;
    if (reveal) ImGui::SetScrollHereY(0.5f);  // a pick in the preview: the row comes into view in the list
    ImGui::PopID();
    return reveal;
}

void componentsSection(SceneUi& ui, const knit::ImportedMesh& m) {
    ImGui::SeparatorText("Components");
    const int total = static_cast<int>(m.nodes.size()), knitted = knittedNodeCount(ui);
    char b[96];
    const bool fur = ui.workflowChosen && knit::furWorkflow(ui.s);
    std::snprintf(b, sizeof b, fur ? "%d of %d with fur" : "%d of %d to knit", knitted, total);
    ImGui::AlignTextToFramePadding();
    ImGui::TextUnformatted(b);
    {
        const ImGuiStyle& st = ImGui::GetStyle();
        const float buttons = ImGui::CalcTextSize("Knit all").x + ImGui::CalcTextSize("Keep all as mesh").x + 4 * st.FramePadding.x + st.ItemSpacing.x;
        ImGui::SameLine(rightAligned(buttons));
        if (ImGui::Button("Knit all")) selectAllNodes(ui);
        ImGui::SameLine();
        if (ImGui::Button("Keep all as mesh")) selectNoNodes(ui);
        mutedWrapped("Unticked components keep their mesh and textures.");
    }
    const float rows = static_cast<float>(std::min(m.nodes.size(), kNodeRowsShown));
    const float height = rows * ImGui::GetFrameHeightWithSpacing() + ImGui::GetStyle().WindowPadding.y * 2;
    // A pick reveals its row through both scroll containers: the form brings the list on screen (its centre),
    // then the list its row. The request stays until a frame has drawn the row.
    const float listScreenY = ImGui::GetCursorScreenPos().y;
    const bool listOnScreen = ImGui::IsRectVisible(ImVec2(ImGui::GetContentRegionAvail().x, height));
    if (ui.import.revealFocus) ImGui::SetScrollFromPosY(listScreenY + height * 0.5f - ImGui::GetWindowPos().y, 0.5f);
    ImGui::BeginChild("##components", ImVec2(0, height), ImGuiChildFlags_Borders);
    bool revealed = false;
    for (const knit::ImportNode& node : m.nodes) revealed |= componentRow(ui, node, listOnScreen);
    ImGui::EndChild();
    if (revealed) ui.import.revealFocus = false;
    if (!ui.import.focusedNode.empty()) mutedWrapped("Shown in the preview" + std::string(kDot) + ui.import.focusedNode);
    else mutedWrapped("Click a name, or the object, to show a component in the preview.");
}

// A preset's components: the same rows as an import's, the ticks fixed by the example (its mesh selection,
// SceneEntry::glb.plainNodes). Listed once per preset.
void presetComponents(SceneUi& ui, const knit::SceneEntry& e) {
    const std::vector<knit::GlbNode>& nodes = presetNodes(ui, e);
    if (nodes.empty()) return;
    ImGui::SeparatorText("Components");
    int knitted = 0;
    for (const knit::GlbNode& n : nodes) knitted += n.knitted ? 1 : 0;
    char b[96];
    std::snprintf(b, sizeof b, "%d of %d with yarn or fur, fixed by the example", knitted, static_cast<int>(nodes.size()));
    ImGui::TextUnformatted(b);
    const ImGuiStyle& st = ImGui::GetStyle();
    const float right = ImGui::CalcTextSize("Keep mesh").x + ImGui::GetFrameHeight() + st.ItemInnerSpacing.x + st.ItemSpacing.x;
    ImGui::BeginDisabled();
    for (const knit::GlbNode& n : nodes) {
        ImGui::PushID(n.name.c_str());
        ImGui::AlignTextToFramePadding();
        ImGui::TextUnformatted(n.name.c_str());
        ImGui::SameLine();
        ImGui::TextDisabled("%s triangles", fmtCount(n.triangles).c_str());
        ImGui::SameLine(rightAligned(right - st.ItemSpacing.x));
        bool ticked = n.knitted;
        const ImGuiID tickId = ImGui::GetID("###knit");
        ImGui::Checkbox(ticked ? "Knit###knit" : "Keep mesh###knit", &ticked);
        testLabel(tickId, ("Knit " + n.name).c_str());
        ImGui::PopID();
    }
    ImGui::EndDisabled();
    mutedWrapped("Unticked components keep their mesh and textures.");
}

// --- real-world size

// Width, depth and height after the unit and up axis, the height marked: the import's own numbers.
void sizeReadouts(const ImportSize& size) {
    if (!ImGui::BeginTable("##importSize", 3, ImGuiTableFlags_SizingStretchSame)) return;
    const char* names[3] = {"Width", "Depth", "Height"};
    const double values[3] = {size.width, size.depth, size.height};
    for (int i = 0; i < 3; ++i) {
        ImGui::TableNextColumn();
        if (i == 2) ImGui::TableSetBgColor(ImGuiTableBgTarget_CellBg, ImGui::GetColorU32(ImGuiCol_Header));
        ImGui::TextDisabled("%s", names[i]);
        if (size.valid) ImGui::Text("%.1f mm", values[i]);
        else ImGui::TextUnformatted("-");
    }
    ImGui::EndTable();
}

void sizeSection(SceneUi& ui, const knit::ImportAnalysis& an) {
    ImportFlow& im = ui.import;
    ImGui::SeparatorText("Real-world size");
    if (choiceButton("File units", !im.sizeByHeight)) chooseSizeByHeight(ui, false);
    ImGui::SameLine();
    if (choiceButton("Set height", im.sizeByHeight)) chooseSizeByHeight(ui, true);
    if (im.sizeByHeight) {
        double heightMm = im.targetHeightMm;
        if (ImGui::InputDouble("Target height [mm]", &heightMm, 0, 0, "%.1f")) setImportHeight(ui, heightMm);
    } else {
        const char* units[] = {"1 mm", "1 cm", "1 m", "1 in", "Custom"};
        const UnitChoice choices[] = {UnitChoice::Mm, UnitChoice::Cm, UnitChoice::M, UnitChoice::Inch, UnitChoice::Other};
        int shown = 4;
        for (int i = 0; i < 5; ++i)
            if (im.unit == choices[i]) shown = i;
        ImGui::AlignTextToFramePadding();
        ImGui::TextUnformatted("1 file unit =");
        ImGui::SameLine();
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 7);
        if (combo("##fileUnit", &shown, units, 5)) chooseImportUnit(ui, choices[shown]);
        if (im.unit == UnitChoice::Other) {
            double mmPerUnit = im.customUnitMm;
            if (ImGui::InputDouble("mm per file unit", &mmPerUnit, 0, 0, "%.6g")) setCustomImportUnit(ui, mmPerUnit);
        }
    }
    int up = im.options.up == knit::UpAxis::Y ? 0 : 1;
    const char* axes[] = {"Y (OBJ)", "Z"};
    ImGui::AlignTextToFramePadding();
    ImGui::TextUnformatted("Up axis");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 7);
    if (combo("##upAxis", &up, axes, 2)) chooseImportUp(ui, up == 0 ? knit::UpAxis::Y : knit::UpAxis::Z);
    const ImportSize size = importSize(ui);
    sizeReadouts(size);
    if (!size.valid) {
        wrapped(kWarn, sizeProblem(ui));
        return;
    }
    if (!im.unitTouched)
        wrapped(kWarn, suggestedUnitName(an.suggestedUnitMm) + " suggested from the model bounds (" + an.suggestionReason + "): check the height before importing.");
    if (knit::furWorkflow(ui.s) || !(tickedNodeArea(ui) > 0)) return;
    const double stitches = importStitchCount(ui);
    char b[200];
    std::snprintf(b, sizeof b, "Up to about %s stitches at the current stitch size (%.1f x %.1f mm).", fmtCount(static_cast<long long>(stitches)).c_str(),
                  ui.s.wale, ui.s.course);
    mutedWrapped(b);
    if (stitches > 300000) wrapped(kWarn, "That is a lot of knitting: Final quality may take hours. A smaller size or larger stitches knit faster.");
    if (size.height < 5 * ui.s.course) wrapped(kWarn, "The object is only a few stitches tall at this size: most parts would not be knitted.");
}

// --- surface preparation

// Import > Surface preparation: the components marked Knit become one closed surface, built in the app.
void surfaceSection(SceneUi& ui) {
    if ((ui.workflowChosen && knit::furWorkflow(ui.s)) || knittedNodeCount(ui) == 0) return;
    ImGui::SeparatorText("Surface preparation");
    mutedWrapped("The components marked Knit are rebuilt as one closed surface at stitch scale, with new knitting charts and the base "
                 "colour baked onto them.");
    mutedWrapped("Runs in the app; no other software needed. Usually under a minute.");
    if (ImGui::TreeNode("Preparation details")) {
        char b[320];
        std::snprintf(b, sizeof b, "Interiors are removed and gaps under %.0f mm are sealed. The knitting charts are new, and the file's base colour "
                                   "is baked onto them. It takes longer if the surface has to be closed again.",
                      2.9 * ui.import.options.wrapStitchMm);
        mutedWrapped(b);
        ImGui::TreePop();
    }
}

// --- colour and compatibility

// Import > Colour and compatibility: the base colour in one line, what needs attention beside it, and the
// file's details (counts, colour files, what is not used) in a disclosure.
void colourSection(SceneUi& ui, const knit::ImportedMesh& m, const knit::ImportAnalysis& an) {
    ImGui::SeparatorText("Colour and compatibility");
    ImGui::TextWrapped("%s", importColourSummary(ui).c_str());
    for (const std::string& wn : ui.import.session.snapshot().textureWarnings) wrapped(kWarn, wn);
    for (const knit::ImportIssue& is : an.issues)
        if (is.id == "missing_file") wrapped(kWarn, is.text + " " + is.action);
    if (!ImGui::TreeNode("File details & compatibility")) return;
    char b[320];
    std::snprintf(b, sizeof b, "%s vertices%s%s triangles%s%s%d open boundar%s", fmtCount(static_cast<long long>(m.pos.size())).c_str(), kDot,
                  fmtCount(static_cast<long long>(m.tri.size())).c_str(), m.polygons ? (" (from " + fmtCount(m.faces) + " faces)").c_str() : "", kDot,
                  an.boundaryLoops, an.boundaryLoops == 1 ? "y" : "ies");
    mutedWrapped(b);
    if (an.colourSource.rfind("texture ", 0) == 0) mutedWrapped("Base colour files" + std::string(kDot) + an.colourSource.substr(8));
    mutedWrapped("Only the base colour is imported; roughness, metalness, normal and other maps are not used.");
    for (const knit::ImportIssue& is : an.issues) {
        const bool knitWarning = !knit::furWorkflow(ui.s) && is.kind == knit::ImportIssue::Kind::Warning && is.id != "missing_file";
        if (is.kind == knit::ImportIssue::Kind::Feature || knitWarning) ImGui::BulletText("%s %s", is.text.c_str(), is.action.c_str());
    }
    ImGui::TreePop();
}

// --- stage 2

// The file is being read: its stage, a bar, Cancel.
void importReadingView(SceneUi& ui, const knit::ImportSession::Snapshot& file) {
    glyph(StageStatus::Running);
    ImGui::SameLine();
    ImGui::Text("%s ...", file.stage.c_str());
    ImGui::ProgressBar(static_cast<float>(file.fraction), ImVec2(-FLT_MIN, 0));
    if (ImGui::Button("Cancel")) closeImport(ui);  // the session's reset also stops the reading
}

// Why the file is refused: its message, and the way back to step 1.
void importFailedView(SceneUi& ui, const knit::ImportSession::Snapshot& file) {
    glyph(StageStatus::Failed);
    ImGui::SameLine();
    ImGui::TextUnformatted("Cannot import");
    wrapped(kBad, file.error);
    if (ImGui::Button("Close")) closeImport(ui);
}

// What stops the import whatever is chosen, right under the file's name.
void importBlockers(const knit::ImportAnalysis& an) {
    for (const knit::ImportIssue& is : an.issues)
        if (is.kind == knit::ImportIssue::Kind::Fatal) {
            glyph(StageStatus::Failed);
            ImGui::SameLine();
            wrapped(kBad, "Cannot import: " + is.text + " " + is.action);
        }
}

void importContents(SceneUi& ui, const knit::ImportedMesh& m, const knit::ImportAnalysis& an) {
    importBlockers(an);
    purposeSection(ui, an);
    componentsSection(ui, m);
    sizeSection(ui, an);
    surfaceSection(ui);
    colourSection(ui, m, an);
}

// Stage 2 of the Import step: the open file, read and refused or waiting for its choices.
void importFileView(SceneUi& ui) {
    const knit::ImportSession::Snapshot file = ui.import.session.snapshot();
    stepHeading(fs::path(file.path).filename().string().c_str(), importSourceLine(ui));
    if (file.phase == knit::ImportSession::Phase::Reading) {
        importReadingView(ui, file);
    } else if (file.mesh && file.analysis) {
        ImGui::BeginDisabled(file.phase == knit::ImportSession::Phase::Preparing);
        importContents(ui, *file.mesh, *file.analysis);
        ImGui::EndDisabled();
    } else {
        importFailedView(ui, file);
    }
}

// Stage 1: the object so far, and the ways to open another.
void importStartView(SceneUi& ui) {
    ImGui::SeparatorText("Object");
    if (const knit::SceneEntry* e = entryOf(ui)) {
        ImGui::TextUnformatted(e->label.c_str());
        if (ui.s.scene == "import") {
            mutedWrapped(ui.s.importSource);
            char b[200];
            const bool rebuilt = ui.s.importSurface == "wrapped";
            const std::string surface = rebuilt ? "; rebuilt as a closed surface" : "; the file's own mesh";
            std::snprintf(b, sizeof b, "1 file unit = %.4g mm, %s up%s", ui.s.importUnitMm, ui.s.importUp == "y" ? "Y" : "Z", surface.c_str());
            mutedWrapped(b);
            for (const std::string& c : ui.import.changes) mutedWrapped((rebuilt ? "Rebuilt: " : "Repair: ") + c);
            if (!ui.s.importPlainNodes.empty()) mutedWrapped("Plain meshes: " + ui.s.importPlainNodes);
        } else {
            mutedWrapped("Example scene: imported with a fixed size and mesh selection.");
            if (e->isPreset()) mutedWrapped("Made from " + knit::creditLine(e->credit) + ".");
            presetComponents(ui, *e);
        }
        workflowChoice(ui);
        if (!ui.import.note.empty()) wrapped(kWarn, ui.import.note);
    } else {
        mutedWrapped("No object yet. Open an OBJ or PLY file, drop one onto the window, or start from an example.");
    }
    ImGui::Spacing();
    ImGui::BeginDisabled(ui.job.active);
    if (primaryButton("Open mesh (OBJ, PLY) ...", ui.platform.dialogs && !ui.job.active,
                      ui.platform.dialogs ? waitReason(ui).c_str() : "The file dialog is not available here: drop a file onto the window."))
        ui.platform.request = PlatformRequest::OpenMesh;  // main.cpp opens the dialog, then openMeshFile
    if (ui.s.scene == "import" && ImGui::Button("Change the import ...", ImVec2(-FLT_MIN, 0))) {
        changeImport(ui);
    }
    ImGui::EndDisabled();
    mutedWrapped("Or drop an OBJ or PLY file anywhere in the window. Limits: 200 MB of mesh data (textures not counted), 200 components, 2 million triangles.");
    ImGui::SeparatorText("Examples");
    ImGui::BeginDisabled(ui.job.active);
    for (const knit::SceneEntry& e : knit::sceneRegistry()) {
        if (!knit::exampleModelPresent(e)) continue;  // its model is not in this checkout (scene_registry.h)
        if (ImGui::Button(e.label.c_str(), ImVec2(-FLT_MIN, 0))) openExample(ui, e.id);
    }
    ImGui::EndDisabled();
    ImGui::SeparatorText("Method");
    mutedWrapped("Knitting follows each flattened surface piece (UV chart) of the mesh with a regular lattice of stitches. Each chart's "
                 "rows continue a few stitches past its border onto its neighbours, so the seams are covered by overlapping yarn (Seams: "
                 "Overlapped); the seams are not sewn and the stitch count does not shape the fabric. Every import and example is knitted "
                 "this way except the teddy: it is knitted as spiral pieces (one closed spiral per body part, shaped by "
                 "increases and decreases). This is the project's own method, not Stitch Mesh / Stitch Meshing.");
    mutedWrapped("An opened mesh is first rebuilt as one closed surface at stitch scale, with few wide charts and its base colours baked "
                 "onto them, so that any model can be knitted. Components left unticked stay as the file has them.");
}

// --- the footer

// The running preparation: a rebuild reports how far it is, the other preparations show their stage.
void importPreparingFooter(SceneUi& ui, const knit::ImportSession::Snapshot& file) {
    const bool rebuild = ui.import.options.surface == knit::ImportOptions::Surface::Wrapped;
    JobProgressView v;
    v.title = rebuild ? "Rebuilding the surface ..." : "Preparing ...";
    v.fraction = rebuild ? static_cast<float>(file.fraction) : -1.0f;
    v.qualifier = file.stage;
    jobProgressBlock(v, "");
    if (cancelButton("Cancel preparing", false)) ui.import.session.cancel();  // back at Analysed at once
}

} // namespace

void importPanel(SceneUi& ui) {
    if (ui.import.fileOpen) importFileView(ui);
    else importStartView(ui);
}

// Below the choices: the running preparation; or what the import will make, a failure's message, and the
// import and Cancel buttons. A failure keeps every choice, so the import can be corrected and tried again.
void importFooter(SceneUi& ui) {
    const knit::ImportSession::Snapshot file = ui.import.session.snapshot();
    ImGui::PushID("importFooter");
    ImGui::Separator();
    if (file.phase == knit::ImportSession::Phase::Preparing) {
        importPreparingFooter(ui, file);
        ImGui::PopID();
        return;
    }
    if (file.phase == knit::ImportSession::Phase::Failed) {
        // In a box of its own: a long message must not push the buttons out of the window.
        const float height = std::min(ImGui::CalcTextSize(file.error.c_str(), nullptr, false, ImGui::GetContentRegionAvail().x).y + ImGui::GetStyle().ItemSpacing.y, uiPx(160));
        ImGui::BeginChild("##importFailure", ImVec2(0, height), ImGuiChildFlags_None);
        wrapped(kBad, file.error);
        if (ImGui::IsWindowAppearing()) ImGui::SetScrollHereY(1.0f);  // a long message's last line says what went wrong
        ImGui::EndChild();
        if (canKeepFileAfterFailedRebuild(ui)) {  // the user's choice, never automatic
            mutedWrapped("Or keep the file's own surface: the Mark step and the guided stitch mesh work on it; surface tiling needs a "
                         "rebuilt surface.");
            if (ImGui::Button("Keep the file's own surface")) keepFileAfterFailedRebuild(ui);
        }
    }
    const ImportFooterView view = importFooterView(ui);
    ImGui::TextWrapped("%s", view.summary.c_str());
    wrapped(view.whyDisabled.empty() ? muted() : kWarn, view.status);
    // One row, as the space below the form is small: the action takes what Cancel leaves.
    const float cancelWidth = ImGui::CalcTextSize("Cancel").x + ImGui::GetStyle().FramePadding.x * 6;
    const float actionWidth = std::max(ImGui::GetContentRegionAvail().x - cancelWidth - ImGui::GetStyle().ItemSpacing.x, cancelWidth);
    if (primaryButton(view.action.c_str(), view.whyDisabled.empty(), view.whyDisabled.c_str(), actionWidth)) startImport(ui);
    ImGui::SameLine();
    if (ImGui::Button("Cancel", ImVec2(-FLT_MIN, ImGui::GetFrameHeight() * 1.35f))) closeImport(ui);
    ImGui::PopID();
}

} // namespace viewer
