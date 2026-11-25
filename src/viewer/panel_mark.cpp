// The Mark step's panel and footer (interface.h): the imported file's pieces, the Stitch Mesh view's section (the built
// mesh's headline, switches, legend, measured values and Export Stitch Mesh), suggestions for the isolated piece (setups
// to preview, accept or reject), its markings with their meanings, undo and the markings file, the build settings, the
// quick check while drawing and the last layout's findings and measures in four separate groups (validity, layout,
// coverage, cost), each saying whether it was computed. The actions are marking_flow.h's; this file only draws and calls
// them.
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <string>
#include <vector>

#include "imgui.h"
#include "viewer/interface.h"

namespace viewer {

namespace {

constexpr knit::MarkingMeaning kMeanings[6] = {knit::MarkingMeaning::Start,    knit::MarkingMeaning::Finish, knit::MarkingMeaning::RowStartPath,
                                               knit::MarkingMeaning::Boundary, knit::MarkingMeaning::Cut,    knit::MarkingMeaning::Attachment};
const char* const kMeaningLabels[6] = {"Start", "Finish", "Row-start path", "Boundary", "Cut", "Attachment"};
const char* const kMeaningNotes[6] = {
    "Where the knitting begins: a closed curve around the piece (cast-on), or a point (pole start, not available yet).",
    "Where the knitting ends: a closed curve around the piece (bind-off); the cap beyond it waits for a closure.",
    "Where every round begins and the helix climbs: an open path from the start curve to the finish curve. It never cuts.",
    "Where covering stops, e.g. inside the shoulder: a closed curve; the surface beyond it is left out.",
    "A deliberate change of construction (stored and drawn; splitting into panels is not available yet).",
    "Where two components connect later (stored and drawn; joins are not available yet).",
};

int meaningIndex(knit::MarkingMeaning m) {
    for (int i = 0; i < 6; ++i)
        if (kMeanings[i] == m) return i;
    return 0;
}

double secondsNow() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

const char* causeWord(knit::FindingCause c) {
    switch (c) {
    case knit::FindingCause::Marking: return "marking";
    case knit::FindingCause::AlgorithmLimitation: return "not available yet";
    case knit::FindingCause::Input: return "input";
    case knit::FindingCause::SolverLimitation: return "solver limit";
    }
    return "";
}

StageStatus findingStatus(knit::FindingLevel l) {
    switch (l) {
    case knit::FindingLevel::Ok: return StageStatus::Ready;
    case knit::FindingLevel::Warning: return StageStatus::NeedsUpdate;
    case knit::FindingLevel::Problem: return StageStatus::Failed;
    case knit::FindingLevel::Incomplete: return StageStatus::Diagnostic;
    }
    return StageStatus::Diagnostic;
}

std::string markingTitle(const knit::SurfaceMarking& m) {
    std::string t = std::string(knit::meaningLabel(m.meaning)) + " #" + std::to_string(m.id) + "  (" + std::to_string(m.points.size()) + " points)";
    if (m.setup > 0) t += "  suggestion set " + std::to_string(m.setup);
    if (m.origin == knit::MarkingOrigin::Proposal) t += m.proposal == knit::ProposalState::Offered ? "  proposed" : m.proposal == knit::ProposalState::Accepted ? "  accepted" : "  rejected";
    return t;
}

void piecesSection(SceneUi& ui) {
    const MarkingFlow& mk = ui.mark;
    sectionHeading("Pieces", std::to_string(mk.pieceLabels.size()));
    mutedWrapped("Click a piece in the view or here to isolate it.");
    const float rows = static_cast<float>(std::min<size_t>(mk.pieceLabels.size(), 8));
    if (ImGui::BeginChild("##pieces", ImVec2(0, ImGui::GetTextLineHeightWithSpacing() * rows + ImGui::GetStyle().WindowPadding.y), ImGuiChildFlags_Borders)) {
        for (int i = 0; i < static_cast<int>(mk.pieceLabels.size()); ++i) {
            ImGui::PushID(i);
            const int n = i < static_cast<int>(mk.markingsOnPiece.size()) ? mk.markingsOnPiece[i] : 0;
            const std::string label = mk.pieceLabels[i] + (n > 0 ? "  (" + std::to_string(n) + ")" : "");
            if (ImGui::Selectable(label.c_str(), mk.piece == i)) isolatePiece(ui, mk.piece == i ? -1 : i);
            ImGui::PopID();
        }
    }
    ImGui::EndChild();
    ImGui::TextUnformatted("Other pieces");
    const float cell = (ImGui::GetContentRegionAvail().x - 2 * ImGui::GetStyle().ItemSpacing.x) / 3;
    if (choiceButton("See-through", mk.around == NeighbourLook::Transparent, cell)) setNeighbourLook(ui, NeighbourLook::Transparent);
    ImGui::SameLine();
    if (choiceButton("Hidden", mk.around == NeighbourLook::Hidden, cell)) setNeighbourLook(ui, NeighbourLook::Hidden);
    ImGui::SameLine();
    if (choiceButton("Shown", mk.around == NeighbourLook::Shown, cell)) setNeighbourLook(ui, NeighbourLook::Shown);
    if (mk.piece < 0) return;
    // knitting components: one piece knitted as overlapping parts (a collar over a body)
    sectionHeading("Knitting component", std::to_string(mk.knitComponent + 1));
    mutedWrapped("Knit this piece as overlapping parts, such as a collar over a body: each component has its own start, finish and row-start "
                 "path, and lies one yarn layer outside the components before it.");
    const int components = knitComponentsOnPiece(ui);
    for (int c = 0; c < components; ++c) {
        if (c > 0) ImGui::SameLine();
        const std::string label = std::to_string(c + 1) + "##component";
        if (choiceButton(label.c_str(), mk.knitComponent == c)) setKnitComponent(ui, c);
    }
    ImGui::SameLine();
    if (ImGui::SmallButton("Add component")) setKnitComponent(ui, components);
}

// The active marking's details: its meaning and shape, a proposal's purpose and reason, delete.
void markingDetails(SceneUi& ui, const knit::SurfaceMarking& m) {
    ImGui::Indent();
    int meaning = meaningIndex(m.meaning);
    if (combo("##meaning", &meaning, kMeaningLabels, 6)) setMarkingMeaning(ui, m.id, kMeanings[meaning]);
    mutedWrapped(kMeaningNotes[meaningIndex(m.meaning)]);
    const knit::MarkingShape shapes[3] = {knit::MarkingShape::Point, knit::MarkingShape::Path, knit::MarkingShape::ClosedCurve};
    const char* shapeLabels[3] = {"Point", "Open path", "Closed curve"};
    bool first = true;
    for (int s = 0; s < 3; ++s) {
        if (!knit::shapeAllowed(m.meaning, shapes[s])) continue;
        if (!first) ImGui::SameLine();
        first = false;
        if (choiceButton(shapeLabels[s], m.shape == shapes[s])) setMarkingShape(ui, m.id, shapes[s]);
    }
    if (!m.note.empty()) mutedWrapped("Purpose: " + m.note);
    if (!m.reason.empty()) mutedWrapped("Why: " + m.reason);
    if (!m.unsuitable.empty()) wrapped(kWarn, "Unsuitable: " + m.unsuitable);
    if (ImGui::SmallButton("Delete marking")) deleteMarking(ui, m.id);
    if (ui.mark.selectedPoint >= 0) {
        ImGui::SameLine();
        if (ImGui::SmallButton("Delete point")) deleteSelectedPoint(ui);
    }
    ImGui::Unindent();
}

// One card: what the setup is for, why, its limits and its measured numbers, and what can be done with it.
void suggestionCard(SceneUi& ui, int card, const knit::ProposalSetup& setup) {
    MarkingFlow& mk = ui.mark;
    const std::string state = suggestionState(ui, card);
    ImGui::PushID(card);
    ImGui::Separator();
    const char* stateWord = state == "offered" ? "suggested" : state.c_str();
    ImGui::TextWrapped("%d. %s", card, setup.title.c_str());
    ImGui::SameLine();
    ImGui::TextDisabled("(%s)", stateWord);
    if (state == "accepted" || state == "rejected" || state == "edited") {   // decided: the title, and the limits of what is in use
        if (state == "accepted") mutedWrapped("Limits: " + setup.limitation);
        if (setup.refinable && state == "accepted" && ImGui::SmallButton("Refine the start")) refineSuggestion(ui, card);
        ImGui::PopID();
        return;
    }
    mutedWrapped(setup.purpose);
    mutedWrapped("Why: " + setup.reason);
    mutedWrapped("Limits: " + setup.limitation);
    char b[300];
    std::snprintf(b, sizeof b, "About %d rows; about %d stitches round at the start and %d at the finish (estimates from the shape).", setup.rowsAbout, setup.firstRingStitches,
                  setup.lastRingStitches);
    mutedWrapped(b);
    std::snprintf(b, sizeof b, "Left out: %.0f mm2 visible, %.0f mm2 hidden (measured on the piece).", setup.omittedVisibleMm2, setup.omittedHiddenMm2);
    if (setup.omittedVisibleMm2 > 1) wrapped(kWarn, b);
    else mutedWrapped(b);
    if (state == "outdated") wrapped(kWarn, "Made for another gauge, unit or reference: suggest again.");
    const bool offered = state == "offered";
    const bool previewing = offered && mk.previewSetup > 0 && card <= static_cast<int>(mk.cardSetup.size()) && mk.cardSetup[card - 1] == mk.previewSetup;
    if (offered) {
        if (ImGui::SmallButton(previewing ? "Stop preview" : "Preview")) previewSuggestion(ui, previewing ? 0 : card);
        ImGui::SameLine();
        const std::string replaces = suggestionReplaces(ui, card);
        const std::string accept = replaces.empty() ? "Accept" : "Accept and replace your " + replaces;
        if (ImGui::SmallButton(accept.c_str())) acceptSuggestion(ui, card, !replaces.empty());
        if (!replaces.empty() && ImGui::IsItemHovered()) tooltip("Your own markings of the same roles are replaced (one undo step).");
        ImGui::SameLine();
        if (ImGui::SmallButton("Reject")) rejectSuggestion(ui, card);
    }
    if (setup.refinable && (offered || state == "accepted")) {
        if (offered) ImGui::SameLine();
        if (ImGui::SmallButton("Refine the start")) refineSuggestion(ui, card);
        if (ImGui::IsItemHovered())
            tooltip("A deeper calculation in the background: one trial layout with relaxation, nothing holding the start point; the start is moved to where "
                    "it settles and offered as a new suggestion. It is not a verdict on this one, and it never settles yarn.");
    }
    ImGui::PopID();
}

void suggestionsSection(SceneUi& ui) {
    MarkingFlow& mk = ui.mark;
    if (mk.piece < 0) return;
    const knit::MarkingSession::Snapshot snap = mk.session.snapshot();
    const bool forPiece = mk.proposals && mk.proposals->piece == mk.piece;
    sectionHeading("Suggestions", forPiece ? std::to_string(mk.proposals->setups.size()) : std::string());
    const bool working = snap.phase == knit::MarkingSession::Phase::Proposing || snap.phase == knit::MarkingSession::Phase::Refining;
    if (working) {
        ImGui::ProgressBar(snap.progress, ImVec2(-1, 0), snap.progressText.c_str());
        if (ImGui::Button("Cancel##suggest")) cancelSuggestionWork(ui);
    } else if (plainButton(forPiece ? "Suggest again" : "Suggest markings for this piece", !mk.stale && snap.phase == knit::MarkingSession::Phase::Ready,
                           mk.stale ? "the markings are stale" : "the file is being read or a build is running", -FLT_MIN)) {
        suggestMarkings(ui);
    }
    if (!forPiece) {
        mutedWrapped("Suggestions come from this piece's own shape: where knitting starts and finishes and where each round begins. They are offered "
                     "only; nothing is used until you accept a suggestion.");
        return;
    }
    const knit::ProposalSet& ps = *mk.proposals;
    mutedWrapped(ps.shapeText + ".");
    for (const knit::MarkingFinding& f : ps.findings) {
        glyph(findingStatus(f.level));
        ImGui::SameLine();
        ImGui::TextWrapped("%s", f.message.c_str());
    }
    if (mk.previewSetup > 0)
        wrapped(kInfo, "Previewing a suggestion: it is not accepted. Your markings of the same roles are shown grey and are still in use; Knit uses only "
                       "accepted and your own markings.");
    for (int i = 0; i < static_cast<int>(ps.setups.size()); ++i) suggestionCard(ui, i + 1, ps.setups[i]);
}

void markingsSection(SceneUi& ui) {
    const MarkingFlow& mk = ui.mark;
    sectionHeading("Markings", mk.piece >= 0 ? mk.pieceLabels[mk.piece] : std::string());
    if (mk.stale) wrapped(kWarn, "These markings were made on another mesh: they are listed but not applied. Open other markings or start again.");
    if (mk.piece < 0) {
        mutedWrapped("Isolate a piece first.");
        return;
    }
    mutedWrapped("Add a marking, then click points on the piece. Ctrl+click selects a point: the next click moves it, Delete removes it.");
    const float cell = (ImGui::GetContentRegionAvail().x - 2 * ImGui::GetStyle().ItemSpacing.x) / 3;
    for (int i = 0; i < 6; ++i) {
        if (i % 3 != 0) ImGui::SameLine();
        ImGui::PushID(i);
        if (plainButton(kMeaningLabels[i], !mk.stale, "the markings are stale", cell)) addMarking(ui, kMeanings[i]);
        if (ImGui::IsItemHovered()) tooltip(kMeaningNotes[i]);
        ImGui::PopID();
    }
    for (size_t i = 0; i < mk.set.markings.size(); ++i) {
        const knit::SurfaceMarking& m = mk.set.markings[i];
        const bool here = i < mk.tracePiece.size() && (mk.tracePiece[i] == mk.piece || (m.points.empty() && m.id == mk.activeMarking));
        if (!here) continue;
        ImGui::PushID(m.id);
        if (ImGui::Selectable(markingTitle(m).c_str(), m.id == mk.activeMarking)) selectMarking(ui, m.id);
        // a suggestion's markings are accepted together (its card); single proposals (a lab file's) one by one
        if (m.origin == knit::MarkingOrigin::Proposal && m.proposal == knit::ProposalState::Offered && m.setup == 0) {
            if (ImGui::SmallButton("Accept")) acceptProposal(ui, m.id);
            ImGui::SameLine();
            if (ImGui::SmallButton("Reject")) rejectProposal(ui, m.id);
        }
        if (m.id == mk.activeMarking) markingDetails(ui, m);
        ImGui::PopID();
    }
    ImGui::Spacing();
    const float half = (ImGui::GetContentRegionAvail().x - ImGui::GetStyle().ItemSpacing.x) / 2;
    if (plainButton("Undo", !mk.undo.empty(), "nothing to undo", half)) undoMarking(ui);
    ImGui::SameLine();
    if (plainButton("Redo", !mk.redo.empty(), "nothing to redo", half)) redoMarking(ui);
    if (plainButton("Save as ...", ui.platform.dialogs, "no file dialogs here: the markings are saved automatically", half)) ui.platform.request = PlatformRequest::SaveMarkings;
    ImGui::SameLine();
    if (plainButton("Open ...", ui.platform.dialogs, "no file dialogs here", half)) ui.platform.request = PlatformRequest::OpenMarkings;
    if (!mk.savedTo.empty()) mutedWrapped("Saved automatically to " + mk.savedTo);
}

void measureRow(const char* label, const std::string& value) {
    row(label);
    ImGui::TextUnformatted(value.c_str());
}

std::string number(const char* format, double v);

// A measure as shown: its value when computed for today's inputs, marked when older, and never 0 or PASS when not computed.
std::string known(knit::Known k, const std::string& value) {
    switch (k) {
    case knit::Known::No: return "not computed";
    case knit::Known::Outdated: return value + "  (outdated)";
    case knit::Known::Current: return value;
    }
    return value;
}

// The findings grouped by cause, plain words first; a finding with a place has a Show button that highlights it.
void findingList(SceneUi& ui, const knit::MarkingReview& r, bool quick) {
    const knit::FindingCause order[4] = {knit::FindingCause::Marking, knit::FindingCause::AlgorithmLimitation, knit::FindingCause::SolverLimitation, knit::FindingCause::Input};
    const char* const heads[4] = {"Your markings", "Not supported by the mesher (no marking cures it)", "Solver limit (not the markings')", "Input"};
    for (int c = 0; c < 4; ++c) {
        bool any = false;
        for (size_t i = 0; i < r.findings.size(); ++i) {
            const knit::MarkingFinding& f = r.findings[i];
            if (f.cause != order[c]) continue;
            if (!any) ImGui::TextDisabled("%s", heads[c]);
            any = true;
            ImGui::PushID(static_cast<int>(i) + (quick ? 10000 : 0));
            glyph(findingStatus(f.level));
            ImGui::SameLine();
            std::string who = "Piece";
            if (const knit::SurfaceMarking* m = knit::findMarking(ui.mark.set, f.markingId)) who = std::string(knit::meaningLabel(m->meaning)) + " #" + std::to_string(m->id);
            ImGui::TextWrapped("%s: %s", who.c_str(), f.message.c_str());
            if (f.firstTracePoint >= 0 || f.row >= 0 || !f.at.empty()) {
                const bool shown = ui.mark.highlightFinding == static_cast<int>(i) && ui.mark.highlightQuick == quick;
                if (ImGui::SmallButton(shown ? "Shown" : "Show")) highlightFinding(ui, shown ? -1 : static_cast<int>(i), quick);
            }
            ImGui::PopID();
        }
    }
}

void quickSection(SceneUi& ui) {
    const MarkingFlow& mk = ui.mark;
    if (mk.piece < 0) return;
    const knit::MarkingSession::Snapshot snap = mk.session.snapshot();
    const std::shared_ptr<const knit::MarkingSession::QuickView> q = currentQuick(ui);
    sectionHeading("While drawing", q ? "quick check: now" : "checking ...");
    if (mk.stale) {
        mutedWrapped("No checks: these markings were made on another mesh.");
        return;
    }
    if (snap.quickPreparing) mutedWrapped("Preparing the piece for quick checks (once per piece, gauge and reference) ...");
    if (!q) {
        mutedWrapped("Checking the markings ...");
        return;
    }
    findingList(ui, q->review, true);
    const knit::LayoutMeasures& m = q->review.measures;
    if (beginRows("##quickMeasures")) {
        measureRow("Rows", known(m.rowsKnown, std::to_string(m.rows)));
        measureRow("Rows not crossed once", known(m.contourEventsKnown, std::to_string(m.contourEvents)));
        ImGui::EndTable();
    }
    mutedWrapped("Cells, coverage and the layout's checks need a layout build.");
}

// Colour swatches with their labels: what the Stitch Mesh view's colours mean.
void legend(const std::vector<knit::StitchLegendEntry>& entries, const std::vector<int>& counts) {
    const float size = ImGui::GetTextLineHeight();
    for (size_t i = 0; i < entries.size(); ++i) {
        if (i < counts.size() && counts[i] == 0) continue;   // a kind the mesh does not have
        const knit::StitchLegendEntry& e = entries[i];
        ImGui::PushID(static_cast<int>(i));
        ImGui::ColorButton("##swatch", ImVec4(e.rgb[0], e.rgb[1], e.rgb[2], 1.0f),
                           ImGuiColorEditFlags_NoTooltip | ImGuiColorEditFlags_NoPicker | ImGuiColorEditFlags_NoDragDrop, ImVec2(size, size));
        ImGui::SameLine();
        ImGui::TextUnformatted((e.label + (i < counts.size() ? "  " + std::to_string(counts[i]) : std::string())).c_str());
        ImGui::PopID();
    }
}

// The Stitch Mesh view's section: whether a mesh was generated and how well it fits (two facts), the view's switches,
// the legend, the measured values next to their limits, and Export Stitch Mesh. No yarn is made for any of it.
void stitchMeshSection(SceneUi& ui) {
    MarkingFlow& mk = ui.mark;
    if (mk.view != MarkView::StitchMesh || mk.piece < 0) return;
    const knit::MarkedPieceView* b = shownBuild(ui);
    const bool outdated = shownBuildOutdated(ui);
    sectionHeading("Stitch Mesh", !b ? "not built" : outdated ? "outdated" : "current");
    if (!b) {
        mutedWrapped("No stitch mesh for this piece yet: press Build stitch mesh. The view needs no yarn and no Final knit.");
        return;
    }
    const knit::StitchMeshQuality& q = b->quality;
    wrapped(!b->use.drawable ? kBad : q.clean ? kOk : kWarn, knit::stitchMeshHeadline(b->use, q));
    if (outdated) wrapped(kWarn, "Outdated: the markings, the gauge or a preview changed since this build. It is shown and saved as it was built.");
    if (!b->use.drawable) {
        mutedWrapped("The piece is shown as imported. Correct what the findings name and build again.");
        return;
    }
    if (!q.relaxed) mutedWrapped("An unrelaxed layout (Relax was off, or the rows were updated after an edit).");
    ImGui::Checkbox("Source surface", &mk.stitchSource);
    ImGui::SetItemTooltip("The piece's own surface under the stitch faces (coloured by what the build kept or left out).");
    ImGui::SameLine();
    ImGui::Checkbox("Highlight findings", &mk.stitchHighlight);
    ImGui::SetItemTooltip("Colour the faces by what the checks found on them instead of by stitch kind.");
    if (mk.stitchHighlight) {
        int overFit = 0, farOverFit = 0, structural = 0, unsupported = 0;
        for (int f = 0; f < b->stitch.faces(); ++f) {
            const uint8_t finding = b->stitch.finding[f];
            if (finding & knit::kFaceStructural) ++structural;
            else if (finding & knit::kFaceUnsupported) ++unsupported;
            else if ((finding & knit::kFaceOverFit) && b->stitch.deviationMm[f] > 2 * q.fitTargetMm) ++farOverFit;
            else if (finding & knit::kFaceOverFit) ++overFit;
        }
        const int none = q.faces - overFit - farOverFit - structural - unsupported;
        legend(knit::stitchFindingLegend(), {none, overFit, farOverFit, structural, unsupported});
    } else {
        legend(knit::stitchKindLegend(), {q.regular, q.increase, q.decrease, q.special, q.unsupported, q.poleStart});
    }
    for (const std::string& f : knit::stitchMeshFindings(q)) wrapped(kWarn, f + ".");
    if (!knit::stitchMeshRelaxationNote(q).empty()) mutedWrapped(knit::stitchMeshRelaxationNote(q) + ": the solver's limit; the checks above judge the mesh as it is.");
    if (beginRows("##stitchMeasures")) {
        int fewest = 0, most = 0;
        for (size_t r = 0; r < q.facesPerRow.size(); ++r) {
            fewest = r == 0 ? q.facesPerRow[r] : std::min(fewest, q.facesPerRow[r]);
            most = r == 0 ? q.facesPerRow[r] : std::max(most, q.facesPerRow[r]);
        }
        measureRow("Rows", std::to_string(q.rows) + kDot + "stitches per row " + std::to_string(fewest) + " to " + std::to_string(most));
        measureRow("Fit", number("max %.2f mm", q.maxDeviationMm) + kDot + std::to_string(q.facesOverFit) + number(" over %.2f mm", q.fitTargetMm));
        measureRow("Coverage", (q.uncoveredMeasured ? number("%.1f %% uncovered", 100 * q.uncoveredFraction) : std::string("not measurable")) + number(" (limit %.1f %%)", 100 * q.uncoveredLimit) + kDot +
                                   std::to_string(q.coveredTwice) + " covered twice");
        measureRow("Boundary", number("%.2f mm", q.boundaryMm) + number(" (limit %.2f mm)", q.boundaryLimitMm));
        measureRow("Structure", std::to_string(q.inverted) + " inverted" + kDot + std::to_string(q.degenerate) + " degenerate" + kDot + std::to_string(q.intersectingPairs) +
                                    " crossing pairs" + kDot + "connectivity " + (q.connectivityPass ? "PASS" : "FAIL"));
        measureRow("Relaxation", !q.relaxed ? std::string("not requested") : q.relaxRan ? q.relaxStop : "did not run");
        ImGui::EndTable();
    }
    if (detailsNode("Stitches per row", std::to_string(q.facesPerRow.size()) + " rows")) {
        std::string list;
        for (size_t r = 0; r < q.facesPerRow.size(); ++r) list += (r == 0 ? "" : ", ") + std::to_string(q.facesPerRow[r]);
        mutedWrapped(list);
        ImGui::TreePop();
    }
    if (plainButton("Export Stitch Mesh ...", ui.platform.dialogs, "no file dialogs here (script step mark-export-mesh:<path>)", -FLT_MIN))
        ui.platform.request = PlatformRequest::SaveStitchMesh;   // then exportStitchMesh
    ImGui::SetItemTooltip("Saves the stitch polygons shown (not yarn) as an OBJ file, in the imported file's units and axes. Fit findings do not stop it.");
}

void reviewSection(SceneUi& ui) {
    const knit::MarkedPieceView* b = shownBuild(ui);
    if (ui.mark.piece < 0) return;
    if (!b) {
        sectionHeading("Layout", "not built");
        mutedWrapped(ui.mark.layoutRefresh ? "No layout yet: it is built when the markings are complete, or press Build stitch mesh."
                                           : "No layout yet: press Build stitch mesh.");
        return;
    }
    const knit::MarkingReview& r = b->review;
    const bool outdated = shownBuildOutdated(ui);
    char age[120];
    std::snprintf(age, sizeof age, "%s%s, built %.0f s ago", b->previewIds.empty() ? "" : "preview, ", r.depth == knit::ReviewDepth::Relaxed ? "relaxed" : "unrelaxed",
                  secondsNow() - shownBuildArrivedAt(ui));
    sectionHeading("Layout", outdated ? "outdated" : age);
    if (outdated) wrapped(kWarn, "Outdated: the markings, the gauge or a preview changed since this build. Its measures are shown as they were.");
    else if (r.measures.checksKnown == knit::Known::Current && r.measures.checksPass) wrapped(kOk, "The layout passes its checks at this stage; final yarn settling has not been run.");
    else if (b->use.drawable) mutedWrapped(knit::stitchMeshHeadline(b->use, b->quality) + ". See the findings below. Final yarn settling has not been run.");
    else wrapped(kBad, knit::stitchMeshHeadline(b->use, b->quality) + ".");
    if (!b->previewIds.empty()) wrapped(kInfo, "This layout is a preview of a suggestion: it is not saved and not the Knit input.");
    findingList(ui, r, false);
    if (!r.recommendation.empty()) mutedWrapped(r.recommendation);
    const std::string comparison = baselineComparison(ui);
    if (!comparison.empty()) mutedWrapped(comparison);
    const knit::LayoutMeasures& m = r.measures;
    if (!beginRows("##markMeasures")) return;
    char v[160];
    std::snprintf(v, sizeof v, "%d contour %s", m.contourEvents, m.contourEvents == 1 ? "event" : "events");
    measureRow("Validity", known(m.contourEventsKnown, v));
    measureRow("Cells without a construction", known(m.unsupportedKnown, std::to_string(m.unsupportedCells)));
    measureRow("Layout checks", known(m.checksKnown, m.checksPass ? "PASS" : "FAIL"));
    measureRow("Rows", known(m.rowsKnown, std::to_string(m.rows) + kDot + "largest change " + std::to_string(m.maxStitchCountJump) + " stitches"));
    measureRow("Exceptional cells", known(m.exceptionalKnown, std::to_string(m.exceptionalCells)));
    measureRow("Distortion (p95)", known(m.distortionKnown, number("%.3f", m.distortionP95)));
    measureRow("Coverage", known(m.coverageKnown, number("%.2f %%", 100 * m.coveredFraction) + kDot + (m.coveragePass ? "PASS" : "FAIL")));
    measureRow("Left out", known(m.omittedKnown, number("%.0f mm2 visible", m.omittedVisibleMm2) + kDot + number("%.0f mm2 hidden", m.omittedHiddenMm2)));
    measureRow("Time", known(m.secondsKnown, number("%.2f s", m.seconds)));
    ImGui::EndTable();
    if (detailsNode("Details", std::to_string(b->log.size()) + " lines")) {
        for (const std::string& l : b->log) mutedWrapped(l);
        ImGui::TreePop();
    }
}

std::string number(const char* format, double v) {
    char b[64];
    std::snprintf(b, sizeof b, format, v);
    return b;
}

// The settings a build uses (the view itself is chosen in the toolbar: Markings, Rows, Stitch Mesh).
void buildSettingsSection(SceneUi& ui) {
    MarkingFlow& mk = ui.mark;
    sectionHeading("Build settings", "");
    mutedWrapped("The view is chosen above the viewport: Markings, Rows or Stitch Mesh.");
    if (beginRows("##markGauge")) {
        dragRow("Stitch width", &ui.s.wale, kStitchWidthRange, "The stitch width w: the Knit step's gauge (the same setting).");
        dragRow("Row height", &ui.s.course, kRowHeightRange, "The row height h: the Knit step's gauge (the same setting).");
        row("Relax", "Relax the stitch mesh on the surface after building it (Y12 6.1, the experiment's frozen settings).");
        ImGui::Checkbox("##relax", &mk.relax);
        row("Update rows after edits", "After an edit that leaves the rows unknown, an unrelaxed layout is built 0.7 s later, so the row checks become current. "
                                       "Relaxing and knitting stay explicit.");
        bool refresh = mk.layoutRefresh;
        if (ImGui::Checkbox("##refresh", &refresh)) setLayoutRefresh(ui, refresh);
        row("Reference", "The surface the stitch mesh is laid out and relaxed on. Imported facets: the file's surface. Smooth (Loop limit): an experimental "
                         "smooth surface through it; the review states how far it lies from the imported one. Saved with the markings.");
        const bool smooth = mk.set.preparation == "loop_limit";
        ImGui::SetNextItemWidth(-1);
        if (ImGui::BeginCombo("##reference", smooth ? "Smooth (Loop limit)" : "Imported facets")) {
            if (ImGui::Selectable("Imported facets", !smooth)) setMarkingPreparation(ui, "facets");
            if (ImGui::Selectable("Smooth (Loop limit)", smooth)) setMarkingPreparation(ui, "loop_limit");
            ImGui::EndCombo();
        }
        ImGui::EndTable();
    }
}

} // namespace

void markPanel(SceneUi& ui) {
    ImGui::PushID("mark");
    stepHeading("Mark", "Markings for the stitch mesher (optional): where knitting starts and finishes, and where each round begins.");
    if (ui.s.scene != "import" || ui.s.importSource.empty()) {
        mutedWrapped(entryOf(ui) ? "Markings are drawn on an imported file's pieces: import a file (Import)." : "Import an object first (Import).");
        ImGui::PopID();
        return;
    }
    const knit::MarkingSession::Snapshot snap = ui.mark.session.snapshot();
    if (!snap.pieces) {
        if (snap.phase == knit::MarkingSession::Phase::Failed) wrapped(kBad, "Cannot read the file: " + snap.error);
        else mutedWrapped("Reading the file ...");
        ImGui::PopID();
        return;
    }
    if (!ui.mark.note.empty()) mutedWrapped(ui.mark.note);
    piecesSection(ui);
    stitchMeshSection(ui);
    suggestionsSection(ui);
    quickSection(ui);
    reviewSection(ui);
    markingsSection(ui);
    buildSettingsSection(ui);
    ImGui::PopID();
}

bool markFooterShown(const SceneUi& ui) { return ui.s.scene == "import" && ui.mark.session.snapshot().pieces != nullptr; }

void markFooter(SceneUi& ui) {
    ImGui::PushID("markFooter");
    const knit::MarkingSession::Snapshot snap = ui.mark.session.snapshot();
    const bool building = snap.phase == knit::MarkingSession::Phase::Building;
    const bool busy = building || snap.phase == knit::MarkingSession::Phase::Proposing || snap.phase == knit::MarkingSession::Phase::Refining;
    if (building) mutedWrapped("Building the stitch mesh ...");
    else if (!snap.error.empty()) wrapped(kBad, snap.error);
    else if (ui.mark.layoutDueAt > 0) mutedWrapped("The rows update in a moment (an unrelaxed layout) ...");
    else if (shownBuildOutdated(ui)) wrapped(kWarn, "The stitch mesh shown is outdated: the markings, the gauge or a preview changed. Build again.");
    if (ui.mark.previewSetup > 0) wrapped(kInfo, "Building previews the suggestion; it is not saved and not the Knit input.");
    const char* why = ui.mark.piece < 0 ? "isolate a piece first" : ui.mark.stale ? "the markings are stale" : "the session is busy";
    if (primaryButton("Build stitch mesh", ui.mark.piece >= 0 && !ui.mark.stale && !busy, why)) buildMarkedPiece(ui);
    // The guided Knit route (experimental): this piece knitted as a stitch mesh from these markings, the other pieces plain.
    const bool guidedHere = ui.s.knitRoute == "guided" && ui.s.guidedPiece == ui.mark.piece && ui.mark.piece >= 0;
    if (ImGui::Button(guidedHere ? "Knit this piece (guided) - chosen" : "Knit this piece (guided)") && ui.mark.piece >= 0 && !ui.mark.stale) knitThisPiece(ui);
    ImGui::SetItemTooltip("Experimental: the Knit step knits only this piece, as a stockinette stitch mesh from these markings, and settles its yarn; "
                          "the other pieces stay plain meshes. The result says whether the yarn settled.");
    ImGui::PopID();
}

} // namespace viewer
