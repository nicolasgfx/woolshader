// The Knit step (interface.h), in the order a knitter decides: the header (object, scope, Knit this
// object), Construction (four tiles: single tube with its ply texture, ply tubes, fibre strands, fibre maps;
// a yarn of fibres is set up in the Fibres step), Yarn & stitch size with the prepared-surface fit check,
// Stitch pattern (four tiles, each with its flat-swatch verdict at this gauge), Stitch character (preset
// tiles, fine-tuning), Yarn colour, then collapsed details with their current settings in one line (stitch
// layout: seams and row direction; mesh detail). The footer holds the knit's state and its two actions,
// Knit draft and Knit final; the Report window's parts table is here too.
#include "model/choices.h"
#include "model/params.h"
#include "viewer/interface.h"
#include "export/export.h"
#include "jobs/built_scene.h"
#include "jobs/previews.h"
#include "jobs/scene_builder.h"
#include "settings/scene_registry.h"
#include "settings/settings_profiles.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

#include "settings/file_presets.h"
#include "settings/thickness_presets.h"
#include "settings/settings_mapping.h"
#include "settings/yarn_look.h"

namespace viewer {

namespace {

// A twisted-ply style as its cross-section: the maps' radius profile rho(phi) around a unit yarn
// (ply_profile.h), drawn in a small box. Building a profile also tabulates its AO (0.3 to 0.6 s for both
// styles), so both are made once per plies and fill and kept in sketch.
void plyStyleSketch(PlySketch& sketch, int plies, knit::PlyProfileModel model, double fill) {
    if (sketch.plies != plies || sketch.fill != fill) {
        sketch.plies = plies;
        sketch.fill = fill;
        sketch.swept = knit::plyCrossSection(plies, knit::PlyProfileModel::Swept, fill, 160);
        sketch.circles = knit::plyCrossSection(plies, knit::PlyProfileModel::Circles, fill, 160);
    }
    const std::vector<knit::Vec2>& pts = model == knit::PlyProfileModel::Swept ? sketch.swept : sketch.circles;
    const float side = ImGui::GetFrameHeight() * 2.2f, half = side * 0.5f;
    const ImVec2 o = ImGui::GetCursorScreenPos();
    std::vector<ImVec2> screen(pts.size());
    for (size_t i = 0; i < pts.size(); ++i)
        screen[i] = {o.x + half + float(pts[i].x) * half * 0.92f, o.y + half - float(pts[i].y) * half * 0.92f};
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const ImVec2 c{o.x + half, o.y + half};  // rho(phi) is star-shaped about the centre: a fan fills it
    const ImU32 fillCol = ImGui::GetColorU32(ImGuiCol_FrameBg);
    for (size_t i = 0; i < screen.size(); ++i) dl->AddTriangleFilled(c, screen[i], screen[(i + 1) % screen.size()], fillCol);
    dl->AddPolyline(screen.data(), int(screen.size()), ImGui::GetColorU32(ImGuiCol_Text), ImDrawFlags_Closed, 1.0f);
    ImGui::Dummy({side, side});
}

// The step's title, then the object and the settings' scope (Knit or Fur is the Import step's choice).
void knitHeader(SceneUi& ui, const knit::SceneEntry& e) {
    ImGui::AlignTextToFramePadding();
    ImGui::TextUnformatted("Knit");
    const char* kind = ui.s.scene == "import" ? "Imported model" : "Example";
    mutedWrapped(e.label + kDot + kind + kDot + "applies to every part");
}

// --- Construction

// The four constructions (yarn_look.h) as tiles, in the order of their cost.
struct ConstructionTile {
    knit::YarnConstruction construction;
    const char* id;
    const char* label;
    const char* detail;  // the tile's second line
    const char* note;    // the tooltip, and the line under the tiles for the chosen one
};
const ConstructionTile kConstructions[] = {
    {knit::YarnConstruction::SingleTube, "single", "Single tube", "One smooth tube \xC2\xB7 fastest",
     "One validated tube per yarn (the main workflow). Its plies can be shown by the ply texture (below: maps on the tube, no "
     "geometry)."},
    {knit::YarnConstruction::PlyTubes, "ply", "Ply tubes", "Real ply geometry \xC2\xB7 close-ups",
     "Each ply as its own round tube inside the drawn yarn: real geometry for close-ups (more to render). The fuzz grows on the plies."},
    {knit::YarnConstruction::FibreStrands, "strands", "Fibre strands", "Individual fibres \xC2\xB7 detailed close-ups",
     "Each yarn drawn as 1 to 3 plies of long fibres in place of its tube (far more to render), on the same knit: switching knits "
     "nothing again. The fibres are set up in the Fibres step."},
    {knit::YarnConstruction::FibreMaps, "maps", "Fibre maps", "Fibre detail on tubes \xC2\xB7 lower render cost",
     "One tube per yarn under maps baked from those same fibres: about half the memory and render time of the strands; the maps are "
     "baked once per knit and look, then cached. The fibres are set up in the Fibres step."},
};

const float kTwoPi = 6.2831853f;

// A point of a twisted strand across a tile's picture: t in [0, 1] along it, phase around the twist (1.5 turns
// across the picture).
ImVec2 twistPoint(float x0, float x1, float cy, float amplitude, float t, float phase) {
    return {x0 + (x1 - x0) * t, cy + amplitude * std::sin(t * 1.5f * kTwoPi + phase)};
}

void drawTwistedStrand(ImDrawList* dl, float x0, float x1, float cy, float amplitude, float phase, ImU32 col, float width) {
    ImVec2 pts[25];
    for (int i = 0; i <= 24; ++i) pts[i] = twistPoint(x0, x1, cy, amplitude, i / 24.0f, phase);
    dl->AddPolyline(pts, 25, col, ImDrawFlags_None, width);
}

// A construction tile's picture, one yarn across it in the theme's colours: a smooth tube, three twisted
// ply tubes, the same twisted bundle made of thin fibre strands, or the tube with the twist as shaded stripes
// (the maps). An illustration of the construction, not a render.
void constructionSwatch(ImDrawList* dl, ImVec2 lo, ImVec2 hi, knit::YarnConstruction c) {
    dl->AddRectFilled(lo, hi, ImGui::GetColorU32(ImGuiCol_FrameBg), ImGui::GetStyle().FrameRounding);
    const float cy = (lo.y + hi.y) * 0.5f, R = (hi.y - lo.y) * 0.32f, x0 = lo.x + R * 0.6f, x1 = hi.x - R * 0.6f;
    const ImU32 yarn = ImGui::GetColorU32(ImGuiCol_Text, 0.75f), fill = ImGui::GetColorU32(ImGuiCol_Text, 0.5f);
    const ImU32 shade = ImGui::GetColorU32(ImGuiCol_FrameBg, 1.0f);
    dl->PushClipRect(lo, hi, true);
    switch (c) {
    case knit::YarnConstruction::SingleTube:
        dl->AddRectFilled({x0, cy - R}, {x1, cy + R}, fill, R);
        dl->AddLine({x0 + R, cy - R * 0.45f}, {x1 - R, cy - R * 0.45f}, yarn, std::max(1.0f, R * 0.25f));
        break;
    case knit::YarnConstruction::PlyTubes:
        for (int k = 0; k < 3; ++k) drawTwistedStrand(dl, x0, x1, cy, R * 0.6f, k * kTwoPi / 3, k == 1 ? fill : yarn, R * 0.75f);
        break;
    case knit::YarnConstruction::FibreStrands:
        for (int k = 0; k < 14; ++k) drawTwistedStrand(dl, x0, x1, cy, R * (0.35f + 0.35f * float((k * 5) % 7) / 6.0f), k * kTwoPi / 14, yarn, 1.0f);
        break;
    case knit::YarnConstruction::FibreMaps:
        dl->AddRectFilled({x0, cy - R}, {x1, cy + R}, fill, R);
        dl->PushClipRect({x0, cy - R}, {x1, cy + R}, true);
        for (float x = x0 - R; x < x1 + R; x += R * 0.9f) dl->AddLine({x, cy + R}, {x + R * 1.2f, cy - R}, shade, std::max(1.0f, R * 0.22f));
        dl->PopClipRect();
        break;
    }
    dl->PopClipRect();
}

// The tiles: the construction the settings describe is selected; a click applies it (applyYarnLook keeps
// the construction left in ui.panels.yarnMemory).
void constructionTiles(SceneUi& ui) {
    knit::YarnLook look = knit::yarnLookOf(ui.s, ui.panels.yarnMemory);
    const TileRow r = beginTiles(4, true);
    ImDrawList* dl = ImGui::GetWindowDrawList();
    for (int i = 0; i < 4; ++i) {
        const ConstructionTile& t = kConstructions[i];
        TilePicture picture;
        if (tile(r, i, {t.id, t.label, t.detail, t.note, look.construction == t.construction}, &picture) && look.construction != t.construction) {
            look.construction = t.construction;
            knit::applyYarnLook(ui.s, look, ui.panels.yarnMemory);
        }
        constructionSwatch(dl, picture.lo, picture.hi, t.construction);
    }
}

// A yarn of fibres: its look and cost in one line (the full cost line is in the Fibres step), the link to
// its settings there, and a refused body said with the way to its control.
void yarnFibresLines(SceneUi& ui) {
    mutedWrapped(yarnFibresSummary(ui));
    if (const std::string advice = yarnFibresMemoryAdvice(yarnFibresCost(ui)); !advice.empty()) wrapped(kWarn, advice);
    if (textLink("Set up in the Fibres step")) revealYarnFibres(ui, Reveal::Section);
    tooltip("The yarn's fibres are set up in the Fibres step (Yarn fibres): the look, the plies, the fibres per yarn, their structure and "
            "shading. Editing them there knits nothing again.");
    const knit::YarnBodyCheck check = knit::checkYarnBody(ui.s);
    const double longest = ui.results.knit && ui.results.knit->knit ? knit::longestYarnMm(*ui.results.knit->knit) : 0;
    const knit::YarnBodyCheck mapCheck = knit::checkYarnMaps(ui.s, longest);
    if (check.accepted && mapCheck.accepted) return;
    const knit::YarnBodyCheck& refused = check.accepted ? mapCheck : check;
    wrapped(kBad, "The export refuses the yarn's fibres: " + refused.reason + ". " + refused.remedy);
    if (textLink("Correct it in the Fibres step")) revealYarnFibres(ui, check.accepted ? Reveal::Section : revealOfBodyField(check.field));
}

// The single tube's twisted-ply maps (strand.texture_profile): their ply count and the cross-section style
// they depict (sketched). They change no tube, but the profile is a knitting-stage key, so a change marks
// the knit for update.
void plyTextureRows(SceneUi& ui) {
    knit::SceneSettings& s = ui.s;
    knit::YarnLook look = knit::yarnLookOf(s, ui.panels.yarnMemory);
    const knit::YarnLook before = look;
    if (beginRows("##plyTextureRows")) {
        row("Ply texture", "Texture-based plies: height and normal maps of twisted plies on the single yarn tube (photon_path_tracer package). "
                           "The tube itself stays round; no ply geometry is made. Parts knitted at a fallback diameter get a plain tube. A "
                           "change marks the knit for update.");
        const char* details[] = {"Off", "2 plies", "3 plies", "4 plies"};
        int k = look.appearancePlies <= 1 ? 0 : std::min(look.appearancePlies, 4) - 1;
        if (combo("##Ply texture", &k, details, 4)) look.appearancePlies = k == 0 ? 0 : k + 1;
        ImGui::EndTable();
    }
    if (look.appearancePlies > 1) {
        ImGui::TextUnformatted("Ply style");
        const knit::PlyProfileModel styles[] = {knit::PlyProfileModel::Swept, knit::PlyProfileModel::Circles};
        const char* names[] = {"Swept plies", "Circles"};
        const char* tips[] = {"The exact cut of round twisted plies (the geometry of separate ply tubes).", "A closed-form layout of round plies."};
        for (int k = 0; k < 2; ++k) {
            ImGui::PushID(200 + k);
            ImGui::BeginGroup();
            plyStyleSketch(ui.panels.plySketch, look.appearancePlies, styles[k], s.plyFill);
            if (ImGui::RadioButton(names[k], look.style == styles[k])) look.style = styles[k];
            ImGui::EndGroup();
            tooltip(tips[k]);
            ImGui::PopID();
            if (k == 0) ImGui::SameLine(0, ImGui::GetFontSize() * 1.5f);
        }
        mutedWrapped("The maps' cross-section of the plies (not the tube geometry).");
    }
    if (!(look == before)) knit::applyYarnLook(s, look, ui.panels.yarnMemory);
}

// An import's route: the surface tiling (the default), or the guided stitch mesh of one piece (experimental): what will be
// knitted, from which markings, and how the last guided knit ended.
void routeSection(SceneUi& ui) {
    const bool guided = ui.s.knitRoute == "guided";
    sectionHeading("Route", guided ? "experimental" : "");
    const float cell = (ImGui::GetContentRegionAvail().x - ImGui::GetStyle().ItemSpacing.x) / 2;
    if (choiceButton("Surface tiling", !guided, cell)) ui.s.knitRoute = "surface";
    ImGui::SetItemTooltip("Every knitted component on the rebuilt surface's charts (the default).");
    ImGui::SameLine();
    if (choiceButton("Guided stitch mesh", guided, cell)) ui.s.knitRoute = "guided";
    ImGui::SetItemTooltip("Experimental: one piece knitted as a stockinette stitch mesh from its Mark step markings; the other pieces stay plain meshes.");
    if (!guided) return;
    if (ui.s.guidedPiece < 0) {
        mutedWrapped("Choose the piece and its markings in the Mark step: isolate the piece, draw or load markings, then Knit this piece (guided).");
        if (ImGui::Button("Open the Mark step")) ui.stage = Stage::Mark;
        return;
    }
    const MarkingFlow& mk = ui.mark;
    const knit::MarkingSession::Snapshot snap = mk.session.snapshot();
    const bool listed = snap.pieces && snap.path == ui.s.importSource && ui.s.guidedPiece < int(snap.pieces->pieces.size());
    if (beginRows("##guidedRows")) {
        row("Piece", "The one piece the stitch mesh covers. The other pieces of the file stay plain meshes: they are not knitted or joined.");
        std::string label = listed && ui.s.guidedPiece < int(mk.pieceLabels.size()) ? mk.pieceLabels[ui.s.guidedPiece] : "piece " + std::to_string(ui.s.guidedPiece);
        if (ui.s.guidedKnitComponent > 0) label += ", knitting component " + std::to_string(ui.s.guidedKnitComponent + 1);
        ImGui::TextUnformatted(label.c_str());
        row("Markings", "The markings file the Knit step reads (the Mark step saves every edit to it).");
        ImGui::TextUnformatted(std::filesystem::path(ui.s.guidedMarkings).filename().string().c_str());
        row("In use", "The markings the build uses: user markings and accepted proposals on this piece. Offered and rejected proposals never change the knit.");
        if (listed && !mk.stale) {
            int used[6] = {0, 0, 0, 0, 0, 0}, offered = 0, rejected = 0;
            for (const knit::SurfaceMarking& m : mk.set.markings) {
                if (!knit::markingInKnitComponent(*snap.pieces, m, ui.s.guidedPiece, ui.s.guidedKnitComponent)) continue;
                if (m.proposal == knit::ProposalState::Offered) ++offered;
                else if (m.proposal == knit::ProposalState::Rejected) ++rejected;
                else ++used[static_cast<int>(m.meaning)];
            }
            std::string text;
            for (int k = 0; k < 6; ++k)
                if (used[k]) text += (text.empty() ? "" : ", ") + std::to_string(used[k]) + " " + knit::meaningLabel(static_cast<knit::MarkingMeaning>(k));
            if (offered + rejected) text += " (" + std::to_string(offered) + " offered, " + std::to_string(rejected) + " rejected: not used)";
            ImGui::TextWrapped("%s", text.empty() ? "none" : text.c_str());
        } else {
            ImGui::TextDisabled("listed once the Mark step has read the file");
        }
        row("Reference", "The surface the stitch mesh is laid out on (the Mark step's Reference, saved with the markings).");
        if (!listed || mk.stale) ImGui::TextDisabled("listed once the Mark step has read the file");
        else ImGui::TextUnformatted(mk.set.preparation == "loop_limit" ? "Smooth (Loop limit, research: up to 0.7 mm inside the imported shape)" : "Imported facets");
        ImGui::EndTable();
    }
    mutedWrapped("The yarn is generated from the stitch mesh and settled with the experiment's frozen settings. Whether it settled is part of the result: "
                 "yarn that did not settle is shown as not settled and exports only as diagnostic. Stitch pattern and coverage settings do not apply.");
    if (ui.results.knit && ui.results.knit->knit && knit::knitsGuidedPiece(ui.results.knit->knit->settings))
        for (const knit::PartRow& r : ui.results.knit->parts)
            for (const std::string& l : r.lines)
                if (l.rfind("guided stitch mesh", 0) == 0) wrapped(l.find(": settled:") != std::string::npos ? kOk : kWarn, l);
}

// The first choice: how each knitted yarn is built for the export and the render preview.
void constructionSection(SceneUi& ui) {
    sectionHeading("Construction", "");
    constructionTiles(ui);
    knit::YarnLook look = knit::yarnLookOf(ui.s, ui.panels.yarnMemory);
    for (const ConstructionTile& t : kConstructions)
        if (t.construction == look.construction) mutedWrapped(t.note);
    if (look.construction == knit::YarnConstruction::PlyTubes && beginRows("##constructionRows")) {
        int k = std::clamp(look.tubePlies, 2, 4) - 2;
        const char* counts[] = {"2", "3", "4"};
        if (comboRow("Ply tubes", &k, counts, 3)) {
            look.tubePlies = k + 2;
            knit::applyYarnLook(ui.s, look, ui.panels.yarnMemory);
        }
        ImGui::EndTable();
    }
    if (look.construction == knit::YarnConstruction::SingleTube) plyTextureRows(ui);
    if (knit::yarnOfFibres(look.construction)) yarnFibresLines(ui);
}

// --- Yarn & stitch size

// A stitch dimension the Gauge list derives from the yarn, in place of its field: the same row
// height (the rows never shift under the closing Gauge popup), nothing to drag.
void readoutRow(const char* label, double mm, const char* why) {
    row(label);
    ImGui::AlignTextToFramePadding();
    ImGui::Text("%.2f mm", mm);
    tooltip(why);
}

// The yarn diameter, the one Gauge list (gaugeChoices: a preset, a stitch sized for the yarn, or
// Custom), then the stitch: readouts while the yarn sizes it, fields otherwise.
void gaugeRows(knit::SceneSettings& s, const knit::SceneEntry& e) {
    const std::vector<knit::GaugeChoice> choices = knit::gaugeChoices(e);
    const std::string current = knit::currentGaugeChoice(s, e);
    const knit::GaugeChoice* shown = &choices.back();  // Custom
    for (const knit::GaugeChoice& c : choices)
        if (c.id == current) shown = &c;
    if (!beginRows("##gaugeRows")) return;
    const bool fromYarn = s.stitchSize == "yarn";  // Snug, Regular, Loose: w and h are derived
    row("Yarn diameter");
    const bool diameterChanged = dragField("##Yarn diameter", &s.diameter, kYarnDiameterRange);
    tooltip(fromYarn ? "The yarn's thickness. The stitch width and row height follow it."
                     : "The yarn's thickness. Thicker yarn fills the stitches more and needs more room between neighbouring yarns.");
    row("Gauge", shown->note.c_str());
    if (beginList("##Gauge", shown->label.c_str(), ImGuiComboFlags_HeightLargest)) {  // every entry and caption in view
        for (const knit::GaugeChoice& c : choices) {
            if (!c.group.empty()) ImGui::SeparatorText(c.group.c_str());
            if (ImGui::Selectable(c.label.c_str(), c.id == current)) knit::applyGaugeChoice(s, e, c.id);
            tooltip(c.note.c_str());
        }
        ImGui::EndCombo();
    }
    bool stitchChanged = false;
    if (fromYarn) {
        const char* follows = "Follows the yarn diameter and the clearance (Gauge: Snug, Regular or Loose). Choose Custom to set it yourself.";
        readoutRow("Stitch width", s.wale, follows);
        readoutRow("Row height", s.course, follows);
    } else {
        stitchChanged = dragRow("Stitch width", &s.wale, kStitchWidthRange);
        stitchChanged |= dragRow("Row height", &s.course, kRowHeightRange);
    }
    if (diameterChanged || stitchChanged) {
        s.thicknessPreset = "custom";  // an edited preset is Custom; a yarn entry stays one
        knit::applyStitchSize(s);      // yarn mode: w and h follow D, c and the tightness
    }
    ImGui::EndTable();
}

// What the gauge means: the pattern's own gauge (rib and garter knit their relaxed size of the
// stockinette stitch above, patternStitch); one line either way, so the layout stays put.
void gaugeLine(const knit::SceneSettings& s) {
    const knit::StitchPattern pattern = knit::parsePattern(s.pattern).value_or(knit::StitchPattern::Stockinette);
    const knit::StitchSize st = knit::patternStitch(s, pattern);
    const knit::PatternScale k = knit::patternRestScale(pattern);
    if (!(st.w > 0 && st.h > 0)) ImGui::TextDisabled("Stitch width and row height must be positive.");
    else if (k.w == 1 && k.h == 1) ImGui::TextDisabled("%.1f stitches %s %.1f rows / 10 cm", 100 / st.w, kTimes, 100 / st.h);
    else
        ImGui::TextDisabled("%s: %.1f stitches %s %.1f rows / 10 cm (%.2f %s %.2f of this stitch)", patternLabel(pattern), 100 / st.w,
                            kTimes, 100 / st.h, k.w, kTimes, k.h);
}

// Fit check 1, imports with a rebuilt surface only: the stitch width the surface was prepared for
// against the current one (preparedSurfaceCheck). Too wide: the warning and the way back to Import.
void preparedSurfaceLine(SceneUi& ui) {
    const PreparedSurface p = preparedSurfaceCheck(ui.s);
    if (!p.applies) return;
    char b[200];
    if (!p.tooWide) {
        glyph(StageStatus::Ready);
        ImGui::SameLine();
        std::snprintf(b, sizeof b, "Surface prepared for %.1f mm stitches.", p.preparedMm);
        mutedWrapped(b);
        return;
    }
    glyph(StageStatus::NeedsUpdate);
    ImGui::SameLine();
    std::snprintf(b, sizeof b, "Surface prepared for %.1f mm stitches; the current %.2f mm stitches are wider and may tangle.", p.preparedMm, ui.s.wale);
    wrapped(kWarn, b);
    const bool enabled = !ui.job.active && !ui.s.importSource.empty();
    ImGui::BeginDisabled(!enabled);
    if (ImGui::Button("Review Import")) reviewImport(ui);
    ImGui::EndDisabled();
    std::snprintf(b, sizeof b, "Opens the Import step for this object to prepare its surface for %.2f mm stitches. The knitting settings stay as "
                               "they are; the knit then needs an update.", ui.s.wale);
    tooltip(ui.job.active ? waitReason(ui).c_str() : ui.s.importSource.empty() ? "The imported file is not known (no import.source)." : b);
    ImGui::SameLine();
    ImGui::TextDisabled("or make the stitches smaller here.");
}

// Fit check 2: the chosen pattern on a flat swatch (FlatReferences), under the pattern tiles. One short
// line; the explanation is in the help marker. A swatch that does not pass is information, not a warning:
// every gauge preset fails it (thickness_presets_flat_matrix) and knits. A yarn that really does not fit
// warns through thickYarnWarning.
void flatSwatchLine(SceneUi& ui, const knit::SceneEntry& e, const std::vector<knit::PatternReference>* refs) {
    const char* swatch = "The flat swatch is an 8 x 8 stitch patch checked with the same validation as the knit (yarn clearance, curvature, "
                         "crossings), but as the bare stitch template, without settle and repair. It says whether the pattern fits at these "
                         "values on a flat surface; the knitted object's own parts are validated after knitting. Neither result implies the "
                         "other.";
    const knit::StitchPattern pattern = knit::parsePattern(ui.s.pattern).value_or(knit::StitchPattern::Stockinette);
    const FlatSwatchView v = flatSwatchOf(ui, refs, pattern);
    glyph(v.status);
    ImGui::SameLine();
    if (v.status == StageStatus::Optional && knit::knitsSewnPieces(e, ui.s)) {
        const std::string t = "A flat swatch at these values " + v.outcome +
                              ". The teddy's spiral pieces are validated by their own knit: at its calibrated Dense gauge, Final "
                              "passes 14 of 18 pieces (the ears fail).\n\n" + swatch;
        noteWithHelp(muted(), v.text, t.c_str());
    } else {
        noteWithHelp(muted(), v.text, swatch);
    }
}

// How the yarn's tubes touch: its depth and the clearance.
void contactRows(knit::SceneSettings& s) {
    if (beginRows("##shapeRows")) {
        dragRow("Depth", &s.depthAmplitude, {0.3, 1.2, "%.2f D"});
        tooltip("How far the yarn dips in and out of the surface, in yarn diameters.");
        if (dragRow("Clearance", &s.clearance, kClearanceRange)) knit::applyStitchSize(s);  // yarn mode: w and h follow it
        tooltip("The minimum gap the validation requires between neighbouring yarns.");
        ImGui::EndTable();
    }
}

// Gauge fit (yarnFit): a stitch too tight for a round yarn is knitted with the yarn compressed. Drawn
// (disabled) when the yarn fits, so the layout stays put.
void flattenToFitRow(knit::SceneSettings& s) {
    const bool needed = knit::yarnFit(s).need < 1;
    ImGui::BeginDisabled(!needed);
    ImGui::Checkbox("Flatten yarn that is too thick for the stitch", &s.compressToFit);
    ImGui::EndDisabled();
    tooltip(!needed ? "Not needed: this stitch fits a round yarn of this diameter, so nothing is flattened."
                    : "Real wool flattens where yarns cross. The stitch is too small for a round yarn of this diameter, so the yarn is "
                      "knitted at a smaller contact diameter that fits (at least 60 % of it) and drawn at its full diameter, flattened "
                      "where yarns touch. Off: knitted as a round yarn; the parts then fail their checks after long knitting.");
    if (const std::string n = knit::yarnFitNote(s); !n.empty()) mutedWrapped(n);
}

// The yarn's fit on its surface: the compressed fit, its depth, clearance and contact flattening.
void yarnGeometryNode(SceneUi& ui) {
    knit::SceneSettings& s = ui.s;
    if (!detailsNode("Gauge & yarn geometry", yarnGeometrySummary(s))) return;
    flattenToFitRow(s);
    contactRows(s);
    ImGui::BeginDisabled(s.plies > 1);
    ImGui::Checkbox("Contact flattening", &s.flatten);
    tooltip("The rendered tube is inflated and flattened where neighbours touch (a geometric approximation of squashed yarn). The "
            "centreline is still validated at the yarn diameter. Not with separate ply tubes. A compressed yarn (Flatten yarn that is "
            "too thick for the stitch) is always flattened; its ply tubes stay round and cross where yarns touch. Fibres follow the "
            "flattened surface.");
    const bool flattened = knit::renderFlattened(s);
    if ((s.flatten || flattened) && beginRows("##flattenRows")) {
        if (s.flatten) dragRow("Inflation", &s.flattenInflate, {0.0, 0.4, "%.2f"});
        if (flattened)  // the range of FinalSurfaceSpec::valid (kMinMaxCompression .. kMaxMaxCompression)
            percentRow("Max flattening", &s.maxFlattening, 0.05, 0.6,
                       "The most a yarn is squashed where it touches its neighbours, as a share of its radius. Lower keeps thick yarn "
                       "rounder and fuller, so the fabric closes; a change knits again.");
        ImGui::EndTable();
    }
    ImGui::EndDisabled();
    ImGui::TreePop();
}

void yarnStitchSection(SceneUi& ui, const knit::SceneEntry& e) {
    sectionHeading("Yarn & stitch size", gaugeChoiceLabel(ui.s, e));
    gaugeRows(ui.s, e);
    gaugeLine(ui.s);
    preparedSurfaceLine(ui);
    if (const std::string w = knit::thickYarnWarning(ui.s, e); !w.empty()) wrapped(kWarn, w);
    yarnGeometryNode(ui);
}

// --- Stitch character

// A character tile's picture: three rows of knit stitches (V shapes) whose wobble grows with the
// preset's unevenness, column shift and thick-and-thin spots. Drawn from the preset's own values in
// the theme's colours: an illustration of how uneven it is, not a render of the knit.
void characterSwatch(ImDrawList* dl, ImVec2 lo, ImVec2 hi, const knit::HandKnit& h) {
    dl->AddRectFilled(lo, hi, ImGui::GetColorU32(ImGuiCol_FrameBg), ImGui::GetStyle().FrameRounding);
    // The presets' largest values (Rustic: 32 %, 1.0 mm, 40 %) give a wobble of about 0.6.
    const float wobble = std::clamp(float(h.unevenness / 50 * 0.5 + h.columnShift / 2 * 0.3 + h.spotChange / 50 * 0.2), 0.0f, 1.0f);
    const float rowH = (hi.y - lo.y) / 3, cellW = rowH * 1.1f;
    const float lineW = std::max(1.0f, rowH * 0.22f);
    const ImU32 yarn = ImGui::GetColorU32(ImGuiCol_Text, 0.75f);
    dl->PushClipRect(lo, hi, true);
    for (int r = 0; r < 4; ++r) {
        for (int c = 0; c * cellW < hi.x - lo.x + cellW; ++c) {
            const float noise = std::sin(c * 4.5f + r * 3.7f), noise2 = std::cos(c * 2.1f - r * 1.3f);  // a fixed pattern per tile
            const float x = lo.x + c * cellW + noise * wobble * cellW * 0.35f;
            const float y = lo.y + (r - 0.3f) * rowH + noise2 * wobble * rowH * 0.3f;
            const float w = lineW * (1 + noise2 * wobble * 0.6f);
            dl->AddLine({x, y}, {x + cellW * 0.5f, y + rowH}, yarn, w);
            dl->AddLine({x + cellW * 0.5f, y + rowH}, {x + cellW, y}, yarn, w);
        }
    }
    dl->PopClipRect();
}

// The hand-knit presets (presets/hand_knit, in their order) as tiles (beginTiles: 4 across, 2 x 2 or one
// per row when the panel is narrow). The tile of the preset the settings match is selected; none is after
// a manual edit (Custom).
void characterTiles(knit::SceneSettings& s) {
    const std::vector<knit::FilePreset>& list = knit::presets(knit::PresetGroup::HandKnit);
    if (list.empty()) return;
    const knit::FilePreset* match = knit::matchingFilePreset(knit::PresetGroup::HandKnit, s);
    const TileRow r = beginTiles(int(list.size()), false);
    ImDrawList* dl = ImGui::GetWindowDrawList();
    for (size_t i = 0; i < list.size(); ++i) {
        const knit::FilePreset& p = list[i];
        TilePicture picture;
        if (tile(r, int(i), {p.id.c_str(), p.label.c_str(), "", p.note.c_str(), match == &p}, &picture)) knit::applyPreset(knit::PresetGroup::HandKnit, s, p);
        knit::SceneSettings look = s;  // the preset's own values, whatever the settings hold now
        knit::applyPreset(knit::PresetGroup::HandKnit, look, p);
        characterSwatch(dl, picture.lo, picture.hi, look.hand);
    }
}

// How unevenly the object is knitted, value by value: the stitches, the columns, the loops, the
// fabric and the yarn. Editing any value leaves the presets (Custom).
void fineTuneNode(knit::SceneSettings& s) {
    if (!detailsNode("Fine-tune character", characterSummary(s))) return;
    if (beginRows("##handRows")) {
        dragRow("Uneven stitches", &s.hand.unevenness, {0, 50, "%.0f%%"}, "How far each stitch strays from the others: some tighter and thinner, some looser and thicker.");
        dragRow("Column shift", &s.hand.columnShift, {0, 2, "%.2f mm"}, "How far a column of stitches slides up or down against its neighbours, so rows stop lining up.");
        dragRow("Column run", &s.hand.columnRows, {1, 30, "%.0f rows"}, "How many rows a column keeps about the same shift.");
        dragRow("Raised loops", &s.hand.loopRate, {0, 5, "%.2f / cm2"}, "Loops that stand out of the surface, per square centimetre.");
        dragRow("Loop height", &s.hand.loopHeight, {0, 1.2, "%.2f mm"}, "How far the tallest raised loop stands out.");
        dragRow("Thick and thin spots", &s.hand.spotRate, {0, 5, "%.2f / cm2"}, "Places where the yarn is thicker or thinner, per square centimetre.");
        dragRow("Spot strength", &s.hand.spotChange, {0, 50, "%.0f%%"}, "The largest change of the yarn's diameter at a spot.");
        dragRow("Spot length", &s.hand.spotLength, {2, 60, "%.0f mm"}, "How far a spot runs along the yarn.");
        dragRow("Fabric wave", &s.hand.waveHeight, {0, 1, "%.2f mm"}, "A slow rise and fall of the whole surface.");
        dragRow("Variation", &s.hand.seed, 1, 999, "Another number gives another pattern of the same strength.");
        ImGui::EndTable();
    }
    mutedWrapped("Knit draft moves the stitches. The yarn settles into them at Knit final.");
    ImGui::TreePop();
}

void characterSection(knit::SceneSettings& s) {
    const knit::FilePreset* match = knit::matchingFilePreset(knit::PresetGroup::HandKnit, s);
    sectionHeading("Stitch character", match ? "" : "Custom");
    characterTiles(s);
    mutedWrapped(match ? match->note : "Your own values (Fine-tune character).");
    for (const std::string& err : knit::presetErrors(knit::PresetGroup::HandKnit)) wrapped(kWarn, "Hand-made preset not loaded: " + err);
    fineTuneNode(s);
}

// --- Yarn colour

// One colour, per object part, or per stitch (colourwork); the fibres follow the yarn's colour.
void yarnColourSection(knit::SceneSettings& s, const knit::SceneEntry& e) {
    sectionHeading("Yarn colour", "");
    if (beginRows("##colourRows")) {
        row("Source");
        int src = s.yarnColourSource == "texture" ? 2 : s.yarnColourSource == "object" ? 1 : 0;
        const char* sources[] = {"One colour", "From the object", "Colourwork (per stitch)"};
        const char* ids[] = {"single", "object", "texture"};
        // Spiral pieces (the teddy): no colourwork (one colour per bear, or one colour). A yarn of fibres
        // offers it disabled: the export gives a part's fibres one colour (addPartBody).
        const bool pieces = knit::knitsSewnPieces(e, s);
        const bool fibreYarn = s.yarnBody.enabled;
        if (beginList("##Yarn colour", sources[src])) {
            for (int k = 0; k < (pieces ? 2 : 3); ++k) {
                ImGui::BeginDisabled(k == 2 && fibreYarn);
                if (ImGui::Selectable(sources[k], src == k)) s.yarnColourSource = ids[k];
                ImGui::EndDisabled();
                if (k == 2 && fibreYarn) tooltip("Not with a yarn made of fibres: each part's fibres take one colour, the part's yarn colour.");
            }
            ImGui::EndCombo();
        }
        tooltip(pieces ? "From the object: each bear takes its fabric's colour. Colourwork is not available on the teddy's spiral pieces."
                       : "From the object: each knitted part takes the average base colour of the surface it covers (its texture, vertex "
                         "colours or material). Colourwork: every stitch takes the average texture colour over its cell, in at most "
                         "64 yarn colours, so the object's pattern is knitted down to about one stitch (single tube and ply tubes only). "
                         "The fibres take their stitch's colour while they match the yarn.");
        if (src == 0) colourRow("Colour", s.yarnWhite.color);
        ImGui::EndTable();
    }
    if (s.yarnColourSource == "object") mutedWrapped("One colour per part, from the object. Patterns inside a part are not knitted.");
    if (s.yarnColourSource == "texture" && s.yarnBody.enabled)
        wrapped(kWarn, "Colourwork is not exported with a yarn made of fibres: each part's fibres take one colour, so the pattern the "
                       "viewport shows is lost. Choose One colour or From the object, or a single tube for colourwork.");
    else if (s.yarnColourSource == "texture")
        mutedWrapped("Each stitch takes the average texture colour beneath it, in at most 64 yarn colours. Fibres follow it while they "
                     "match the yarn.");
    const std::string conflict = knit::yarnLookConflict(s);
    if (!conflict.empty()) wrapped(kWarn, "Settings file: " + conflict + ".");
}

// --- Stitch pattern

// The patterns the tiles offer, with what each looks like. Rib 2 x 2 is not offered: its overlapping wale
// pairs do not knit cleanly yet; a settings file may still name it.
struct OfferedPattern {
    knit::StitchPattern pattern;
    const char* detail;  // the tile's second line
    const char* note;    // the tooltip, and the line under the tiles for the chosen one
};
const OfferedPattern kOffered[] = {
    {knit::StitchPattern::Stockinette, "Knit V shapes on the front \xC2\xB7 the plain fabric",
     "Knit stitches (V shapes) on the front, purl bumps on the back."},
    {knit::StitchPattern::ReverseStockinette, "Purl bumps on the front \xC2\xB7 inside out",
     "Stockinette turned inside out: purl bumps on the front."},
    {knit::StitchPattern::Garter, "Ridges on both faces \xC2\xB7 shorter rows",
     "Knit and purl rows alternate: horizontal ridges on both faces. A garter row is 0.70 x as tall and 1.10 x as wide as the "
     "stockinette stitch set below, and the fabric folds into ridges (relaxed garter, measured by Hwang et al. 2025)."},
    {knit::StitchPattern::Rib1x1, "Knit and purl columns \xC2\xB7 narrow, stretchy",
     "Knit and purl columns alternate; the knit columns stand in front of the purl columns. A rib stitch is 0.74 x as wide and "
     "1.23 x as tall as the stockinette stitch set below (relaxed rib, measured by Hwang et al. 2025)."},
};
constexpr int kOfferedCount = int(sizeof kOffered / sizeof kOffered[0]);

// A knit stitch (a V) or a purl bump (an arch) in one cell of a pattern tile's picture.
void drawKnitCell(ImDrawList* dl, float x, float y, float w, float h, ImU32 col, float lineW) {
    dl->AddLine({x, y}, {x + w * 0.5f, y + h}, col, lineW);
    dl->AddLine({x + w * 0.5f, y + h}, {x + w, y}, col, lineW);
}

void drawPurlCell(ImDrawList* dl, float x, float y, float w, float h, ImU32 col, float lineW) {
    dl->AddBezierQuadratic({x, y + h * 0.8f}, {x + w * 0.5f, y + h * 0.05f}, {x + w, y + h * 0.8f}, col, lineW * 1.4f);
}

// A pattern tile's picture: rows of stitches, each cell a knit V or a purl bump as the pattern's face rule
// says (stitchFaceIsPurl), regular, in the theme's colours: an illustration of the pattern, not a render.
void patternPicture(ImDrawList* dl, ImVec2 lo, ImVec2 hi, knit::StitchPattern p) {
    dl->AddRectFilled(lo, hi, ImGui::GetColorU32(ImGuiCol_FrameBg), ImGui::GetStyle().FrameRounding);
    const float rowH = (hi.y - lo.y) / 3, cellW = rowH * 1.1f;
    const float lineW = std::max(1.0f, rowH * 0.22f);
    const ImU32 yarn = ImGui::GetColorU32(ImGuiCol_Text, 0.75f);
    dl->PushClipRect(lo, hi, true);
    for (int r = 0; r < 4; ++r)
        for (int c = 0; c * cellW < hi.x - lo.x + cellW; ++c) {
            const float x = lo.x + c * cellW, y = lo.y + (r - 0.3f) * rowH;
            if (knit::stitchFaceIsPurl(p, c, r)) drawPurlCell(dl, x, y, cellW, rowH, yarn, lineW);
            else drawKnitCell(dl, x, y, cellW, rowH, yarn, lineW);
        }
    dl->PopClipRect();
}

// The offered patterns as tiles, each with its flat swatch verdict at the current gauge in the picture's
// corner (the chosen one's verdict is the line below the tiles). On the teddy's spiral pieces every tile but
// Stockinette is disabled (sewn_pieces.cpp: purl increases and decreases are not derived).
void patternTiles(SceneUi& ui, const knit::SceneEntry& e, const std::vector<knit::PatternReference>* refs) {
    knit::SceneSettings& s = ui.s;
    const knit::StitchPattern current = knit::parsePattern(s.pattern).value_or(knit::StitchPattern::Stockinette);
    const bool pieces = knit::knitsSewnPieces(e, s);
    const TileRow r = beginTiles(kOfferedCount, true);
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const float glyphSize = ImGui::GetFontSize() * 0.9f;
    const ImU32 backing = ImGui::GetColorU32(ImGuiCol_WindowBg, 1.0f);  // the glyph on a dark disc: legible over the stitches
    for (int i = 0; i < kOfferedCount; ++i) {
        const OfferedPattern& o = kOffered[i];
        const bool locked = pieces && o.pattern != knit::StitchPattern::Stockinette;
        const FlatSwatchView swatch = flatSwatchOf(ui, refs, o.pattern);
        const std::string note = locked ? std::string("The teddy's spiral pieces knit stockinette only (increases and decreases for purl stitches are not derived).")
                                        : std::string(o.note) + " " + swatch.text;
        TilePicture picture;
        ImGui::BeginDisabled(locked);
        const bool hit = tile(r, i, {knit::toString(o.pattern), patternLabel(o.pattern), o.detail, note.c_str(), o.pattern == current}, &picture);
        ImGui::EndDisabled();
        if (hit && !locked) s.pattern = knit::toString(o.pattern);
        patternPicture(dl, picture.lo, picture.hi, o.pattern);
        const ImVec2 glyphAt{picture.hi.x - glyphSize - 3, picture.lo.y + 2};
        dl->AddCircleFilled({glyphAt.x + glyphSize * 0.5f, glyphAt.y + glyphSize * 0.5f}, glyphSize * 0.62f, backing);
        drawGlyph(dl, glyphAt, glyphSize, swatch.status);
    }
}

// The stitch pattern: the tiles, what the chosen one looks like, and its flat swatch verdict.
void patternSection(SceneUi& ui, const knit::SceneEntry& e) {
    sectionHeading("Stitch pattern", "");
    const std::shared_ptr<const std::vector<knit::PatternReference>> refs = ui.previews.flat.get(ui.s);
    patternTiles(ui, e, refs.get());
    const knit::StitchPattern current = knit::parsePattern(ui.s.pattern).value_or(knit::StitchPattern::Stockinette);
    for (const OfferedPattern& o : kOffered)
        if (o.pattern == current) mutedWrapped(o.note);
    flatSwatchLine(ui, e, refs.get());
}

// --- Details: stitch layout

// The stitch layout: its seams and the direction of its rows.
void layoutRows(knit::SceneSettings& s, bool pieces) {
    // knit.coverage: the uv charts' seams (scene_settings.h).
    row("Seams", "Overlapped: each chart's rows continue a few stitches past its border onto its neighbours, so the seams "
                 "are covered (more stitches, a longer knit). Open: each chart is knitted inside its own border; the seams "
                 "stay bare. UV charts only.");
    // Spiral pieces have no chart seams: their own word, disabled.
    static const char* const kSeams[] = {"Overlapped", "Open", "Sewn pieces"};
    int seams = pieces ? 2 : s.coverage == "overlap" ? 0 : 1;
    ImGui::BeginDisabled(pieces);
    if (combo("##Seams", &seams, kSeams, pieces ? 3 : 2)) s.coverage = seams == 0 ? "overlap" : "component";
    ImGui::EndDisabled();
    if (pieces) tooltip("The teddy's spiral pieces are sewn: no chart seams.");
    row("Row direction");
    ImGui::BeginDisabled(pieces);
    dragField("##Row direction", &s.angleDeg, {-180, 180, "%.0f deg"});
    ImGui::EndDisabled();
    tooltip(pieces ? "Fixed on the teddy's spiral pieces: each piece's rows run around its own axis."
                   : "Rows run across each chart's projected world-up direction, rotated by this angle.");
}

void stitchLayoutNode(SceneUi& ui, const knit::SceneEntry& e) {
    knit::SceneSettings& s = ui.s;
    if (!detailsNode("Stitch layout", stitchLayoutSummary(s, e))) return;
    const bool pieces = knit::knitsSewnPieces(e, s);  // the teddy's spiral pieces
    if (beginRows("##stitchRows")) {
        layoutRows(s, pieces);
        ImGui::EndTable();
    }
    ImGui::TreePop();
}

// --- Details: mesh detail

const char* const kTessellationHelp =
    "Polygons of the yarn tube (sides around, ring spacing): what the viewer draws and the export writes. Changing it rebuilds the tubes from "
    "the knitted yarn (no knitting); fibres must then be generated again (their roots are placed in this tube).";

// Yarn tessellation: the levels are files (presets/strand_mesh/*.preset, strand_mesh_presets.h).
// True when the mesh settings changed.
bool tessellationRows(SceneUi& ui) {
    knit::SceneSettings& s = ui.s;
    bool meshChanged = false;
    const knit::FilePreset* current = knit::matchingFilePreset(knit::PresetGroup::StrandMesh, s);
    const bool custom = ui.panels.customTessellation || current == nullptr;
    row("Yarn tessellation", kTessellationHelp);
    if (beginList("##Yarn tessellation", custom ? "Custom" : current->label.c_str())) {
        for (const knit::FilePreset& l : knit::presets(knit::PresetGroup::StrandMesh)) {
            if (ImGui::Selectable(l.label.c_str(), !custom && current == &l)) {
                if (current != &l) {
                    knit::applyPreset(knit::PresetGroup::StrandMesh, s, l);
                    meshChanged = true;
                }
                ui.panels.customTessellation = false;
            }
            tooltip((l.note + " (" + l.path + ")").c_str());
        }
        if (ImGui::Selectable("Custom", custom)) ui.panels.customTessellation = true;
        tooltip("Your own tessellation (the fields below).");
        ImGui::EndCombo();
    }
    if (custom) {
        meshChanged |= dragRow("Sides around", &s.tubeSides, 3, 32);
        tooltip("Polygon sides around the yarn tube.");
        meshChanged |= dragRow("Ring spacing", &s.tubeRing, {0.05, 2.0, "%.2f mm"});
        tooltip("Distance between the tube's cross-sections along the yarn.");
        if (s.plies > 1) {
            meshChanged |= dragRow("Ply sides", &s.plyTubeSides, 3, 32);
            meshChanged |= dragRow("Ply ring spacing", &s.plyTubeRing, {0.05, 2.0, "%.2f mm"});
        }
    }
    return meshChanged;
}

// The tessellation level, what changing it does, and fibres that it made stale. True when the mesh
// settings changed this frame.
bool meshDetailNode(SceneUi& ui) {
    if (!detailsNode("Mesh detail", meshDetailSummary(ui.s, ui.panels.customTessellation))) return false;
    bool meshChanged = false;
    if (beginRows("##buildRows")) {
        meshChanged = tessellationRows(ui);
        ImGui::EndTable();
    }
    mutedWrapped("Changing it rebuilds the yarn tubes from the knit (no knitting). Generate fibres again so their roots match.");
    if (ui.s.yarnBody.enabled)
        mutedWrapped("The yarn is made of fibres: this tube is what the viewport draws and where the fuzz grows; the export writes the fibres.");
    if (fibresNeedNewTubes(ui)) wrapped(kWarn, "The generated fibres no longer match these yarn tubes: generate them again (Fibres).");
    ImGui::TreePop();
    return meshChanged;
}

// The mesh is rebuilt from the shown knit at once when the level changed; tubes that differ for another
// reason (a settings file) are offered a rebuild. Drawn below the details, so it shows while they are
// collapsed.
void tubeRebuildAction(SceneUi& ui, bool meshChanged) {
    const knit::SceneSettings& s = ui.s;
    for (const std::string& err : knit::presetErrors(knit::PresetGroup::StrandMesh)) wrapped(kWarn, "Yarn tessellation level not loaded: " + err);
    if (ui.results.knit && !ui.results.knit->displayNote.empty()) wrapped(kWarn, ui.results.knit->displayNote);
    std::string why;
    if (!knit::tubeSettingsValid(s, &why)) wrapped(kBad, "Yarn tessellation not applied: " + why + ".");
    if (shownTubesDiffer(ui)) {
        if (meshChanged && !ui.job.active) startStrandMesh(ui);
        else if (primaryButton("Rebuild the yarn tubes (no knitting)", !ui.job.active, waitReason(ui).c_str())) startStrandMesh(ui);
    }
}

void detailsSection(SceneUi& ui, const knit::SceneEntry& e) {
    sectionHeading("Details", "");
    stitchLayoutNode(ui, e);
    const bool meshChanged = meshDetailNode(ui);
    tubeRebuildAction(ui, meshChanged);
}

} // namespace

void knitPanel(SceneUi& ui, const knit::SceneEntry& e) {
    ImGui::PushID("knit");
    knitHeader(ui, e);
    if (!ui.workflowChosen || knit::furWorkflow(ui.s)) {
        ImGui::SeparatorText("Knitting");
        stateCard(knitState(ui), "");
        mutedWrapped(!ui.workflowChosen ? "Choose Knit or Fur for the object in the Import step."
                                        : "Nothing is knitted: the Fur workflow grows fur on the mesh itself (Import -> Fur -> Export).");
        ImGui::Spacing();
        if (primaryButton(ui.workflowChosen ? "Next: Fur" : "Back to Import", true)) ui.stage = ui.workflowChosen ? Stage::Fibres : Stage::Import;
        ImGui::PopID();
        return;
    }
    if (ui.s.scene == "import") routeSection(ui);
    if (ui.import.lacksCharts && !knit::knitsGuidedPiece(ui.s)) {
        ImGui::SeparatorText("Knitting");
        wrapped(kWarn, "This object keeps the file's own mesh (knitting was skipped when it was imported): knitting needs a rebuilt surface.");
        if (primaryButton("Rebuild the surface ...", !ui.job.active && !ui.s.importSource.empty(), waitReason(ui).c_str())) {
            prepareKnittingCharts(ui);
        }
        mutedWrapped("Opens the Import step to rebuild the surface for knitting; its units and its axis stay as they are.");
        ImGui::PopID();
        return;
    }
    constructionSection(ui);
    yarnStitchSection(ui, e);
    patternSection(ui, e);
    characterSection(ui.s);
    yarnColourSection(ui.s, e);
    detailsSection(ui, e);
    ImGui::PopID();
}

namespace {

// The knit's state, once it is not running: its result in one line (a result the settings outdated is
// the previous one), the last attempt, what the export needs, and what this machine cannot run.
void knitStatus(SceneUi& ui, const knit::SceneEntry& e, const StageState& st, bool draftResult, const KnitActionView& draft,
                const KnitActionView& final) {
    const knit::KnitResult* k = knitResult(ui);
    std::string text = st.detail;
    if (k) {
        char b[256];
        std::snprintf(b, sizeof b, "%s#%d: %d part%s, %d PASS%s%s%s%s", st.status == StageStatus::NeedsUpdate ? "previous result " : "", k->id,
                      k->parts, k->parts == 1 ? "" : "s", k->pass, k->fail ? (", " + std::to_string(k->fail) + " FAIL").c_str() : "",
                      k->unresolved ? (", " + std::to_string(k->unresolved) + " UNRESOLVED").c_str() : "",
                      k->plain ? (", " + std::to_string(k->plain) + " plain").c_str() : "", ui.results.knit->cached ? " (from the cache)" : "");
        text = b;
        if (!st.detail.empty()) text += " - " + st.detail;
    } else if (st.status == StageStatus::NotStarted) {
        text = "export needs Final";  // on the state's line: the footer stays short on small windows
    }
    stateCard(st, text);
    // The parts and their checks (with Focus) are in the Report window.
    if (ui.results.knit && !ui.results.partial && !ui.results.knit->parts.empty() && ImGui::TextLink("Parts and checks (Report)")) ui.panels.reportOpen = true;
    // (Once: without a result the attempt's message is already the state's detail.)
    attemptNote(ui.results.knitAttempt, st.detail);
    if (draftResult) wrapped(kWarn, "Export needs a Final knit: this result skipped settle and the surface check.");
    if (!knit::isPreviewProfile(ui.s) && !knit::isFinalProfile(ui.s, e))
        mutedWrapped("Custom profile (from a settings file): not exportable; Knit draft or Knit final replaces it.");
    if (knit::isFinalProfile(ui.s, e) && !knit::thickYarnWarning(ui.s, e).empty())
        wrapped(kWarn, "Yarn too thick for the stitch (Yarn & stitch size): expect FAIL and long knitting.");
    // What cannot run here, for the job each action submits: a note both share is said once.
    for (const std::string& n : final.notes) {
        const bool shared = std::find(draft.notes.begin(), draft.notes.end(), n) != draft.notes.end();
        wrapped(kWarn, shared ? n : "Knit final: " + n);
    }
    for (const std::string& n : draft.notes)
        if (std::find(final.notes.begin(), final.notes.end(), n) == final.notes.end()) wrapped(kWarn, "Knit draft: " + n);
}

// The two actions, each knitting at its quality, or Cancel while the knit runs.
void knitActions(SceneUi& ui, const KnitActionView& draft, const KnitActionView& final) {
    if (runningStage(ui, knit::BuiltScene::Kind::Knit)) {
        if (cancelButton(ui.job.cancelRequested ? "Cancel requested ..." : "Cancel", ui.job.cancelRequested)) cancelJob(ui);
        return;
    }
    const bool enabled = !ui.job.active;
    const float w = (ImGui::GetContentRegionAvail().x - ImGui::GetStyle().ItemSpacing.x) * 0.5f;
    if (plainButton("Knit draft", enabled, waitReason(ui).c_str(), w)) knitAtDraftQuality(ui);
    if (enabled) tooltip(draft.tooltip.c_str());
    ImGui::SameLine();
    if (primaryButton("Knit final", enabled, waitReason(ui).c_str(), w)) knitAtFinalQuality(ui);
    if (enabled) tooltip(final.tooltip.c_str());
}

} // namespace

// The Knit footer: the knit's state and its two actions. The status block grows with its content up
// to a share of the panel and then scrolls; the action row below it always keeps its room.
void knitFooter(SceneUi& ui, const knit::SceneEntry& e) {
    ImGui::PushID("knitFooter");
    const StageState st = knitState(ui);
    const knit::KnitResult* k = knitResult(ui);
    const bool running = runningStage(ui, knit::BuiltScene::Kind::Knit);
    const bool draftResult = !running && k && k->computational && !k->final && st.status != StageStatus::NeedsUpdate;
    const knit::ComputeCaps& caps = knit::computeCaps();
    const KnitActionView draft = knitActionView(ui, KnitQuality::Draft, caps);
    const KnitActionView final = knitActionView(ui, KnitQuality::Final, caps);
    beginFooterStatus("##knitStatus", ImGui::GetFrameHeight() * 1.35f + ImGui::GetStyle().ItemSpacing.y);
    if (running) {
        const JobProgressView v = progressView(ui.job, ui.results.partial != nullptr);
        jobProgressBlock(v, v.count);
    } else {
        knitStatus(ui, e, st, draftResult, draft, final);
    }
    ImGui::EndChild();
    knitActions(ui, draft, final);
    ImGui::PopID();
}

// The Report window's parts table (popOuts): each part's verdict and stitch count, with Focus.
void partsTable(SceneUi& ui, const knit::BuiltScene& b) {
    if (b.parts.empty()) return;
    if (!ImGui::TreeNode("##parts", "Parts (%zu)", b.parts.size())) return;
    if (ImGui::BeginTable("##partsTable", 4, ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_BordersInnerH)) {
        ImGui::TableSetupColumn("Part");
        ImGui::TableSetupColumn("Result");
        ImGui::TableSetupColumn("Stitches", ImGuiTableColumnFlags_WidthFixed);
        ImGui::TableSetupColumn("", ImGuiTableColumnFlags_WidthFixed);
        ImGui::TableHeadersRow();
        for (const knit::PartRow& r : b.parts) {
            ImGui::PushID(&r);
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::TextUnformatted(r.name.c_str());
            if (!r.note.empty() || !r.lines.empty()) {
                std::string t = r.note;
                for (const std::string& l : r.lines) t += (t.empty() ? "" : "\n") + l;
                tooltip(t.c_str());
            }
            ImGui::TableNextColumn();
            const StageStatus st = !r.knitted ? StageStatus::Optional : r.verdict == knit::Status::Pass ? StageStatus::Ready : r.verdict == knit::Status::Fail ? StageStatus::Failed : StageStatus::NeedsUpdate;
            glyph(st, ImGui::GetFontSize() * 0.8f);
            ImGui::SameLine(0, 2);
            ImGui::TextUnformatted(!r.knitted ? "plain" : knit::toString(r.verdict));
            ImGui::TableNextColumn();
            ImGui::Text("%s", fmtCount(r.stitches).c_str());
            ImGui::TableNextColumn();
            if (r.lo.x <= r.hi.x && ImGui::SmallButton("Focus")) {
                ui.camera.focus = Bounds{r.lo, r.hi};
            }
            ImGui::PopID();
        }
        ImGui::EndTable();
    }
    mutedWrapped("PASS: every check holds. FAIL: a check found a defect. UNRESOLVED: a check could not decide. Plain: too small for the "
                 "stitch size, kept as a surface. Details: hover a part, or open the Report (status bar).");
    ImGui::TreePop();
}

} // namespace viewer
