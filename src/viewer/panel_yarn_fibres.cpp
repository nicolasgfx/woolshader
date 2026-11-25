// The Fibres step's Yarn fibres section (interface.h), drawn first when the Knit step's Construction is Fibre
// strands or Fibre maps (yarn_look.h): the settings of the yarn body (yarn_body.*, fibres/yarn_body.h). They
// are Export-stage keys, so an edit here knits nothing again and leaves the fuzz current: the render preview
// and every export build the fibres from the current values, no action needed. Top to bottom: the look (a
// preset of the body's structure and shading, presets/yarn_fibres), the three choices that make the yarn,
// the maps' texel and cavity, its validity and cost (always shown), the maps' state and their bake, then its
// structure and its shading, collapsed. A link from the Knit step or the Export step's fix scrolls here and
// opens the node it names (ui.panels.revealYarnFibres). Arithmetic only per frame.
#include "viewer/interface.h"
#include "export/export.h"
#include "jobs/built_scene.h"

#include <cstdio>
#include <string>

#include "settings/file_presets.h"
#include "settings/yarn_look.h"

namespace viewer {

namespace {

// A refused yarn body's remedy under the control it concerns (value column, bad colour).
void bodyFieldProblem(const knit::YarnBodyCheck& c, knit::YarnBodyField field) {
    if (c.accepted || c.field != field) return;
    ImGui::TableNextRow();
    ImGui::TableNextColumn();
    ImGui::TableNextColumn();
    wrapped(kBad, c.remedy);
}

// The last item drawn comes into view when a link asked for it (ui.panels.revealYarnFibres).
void revealHere(const SceneUi& ui, Reveal what) {
    if (ui.panels.revealYarnFibres == what) ImGui::SetScrollHereY(0.0f);
}

// What the yarn is (strands or maps), where its fibres show, and that nothing here needs an action.
void yarnFibresNotes(bool maps) {
    mutedWrapped(maps ? "The yarn is made of fibres rendered as maps (Construction, Knit step): every export writes one tube per yarn under "
                        "maps baked from those fibres. The viewport shows the tubes under their cavity map once the maps are baked; the "
                        "full maps appear in the render preview and every export. The fuzz below grows on top of them."
                      : "The yarn is made of fibre strands (Construction, Knit step). The viewport shows a simplified yarn (the tube); the "
                        "strands appear in the render preview (Preview, status bar) and in every export, in place of the tube. The fuzz "
                        "below grows on top of them.");
    mutedWrapped("Edits here need no action: the render preview and the export use the current values. The Reference patch (Preview "
                 "target) shows the fuzz on a plain tube, not the yarn's fibres.");
}

// The look, the three choices that make the yarn and, under maps, the texel and the cavity shading.
void yarnFibresRows(SceneUi& ui, const knit::YarnBodyCheck& check, const knit::YarnBodyCheck& mapCheck, bool maps) {
    knit::SceneSettings& s = ui.s;
    knit::SceneSettings::YarnBody& b = s.yarnBody;
    if (!beginRows("##yarnBodyRows")) return;
    presetRow(s, knit::PresetGroup::YarnFibres, "Look", "Your own fibres (editing any value below makes the look Custom).");
    dragRow("Plies", &b.plies, 1, 3, "The yarn's plies, each a bundle of fibres twisted around the yarn's axis (3 is the standard yarn).");
    bodyFieldProblem(check, knit::YarnBodyField::Plies);
    dragRow("Fibres per yarn", &b.fibres, 1, 2000, "Fibres in the plies' shells for the whole yarn, split evenly across its plies "
                                                   "(400 with 3 plies: about 133 per ply).");
    bodyFieldProblem(check, knit::YarnBodyField::Fibres);
    dragRow("Fibre diameter \xC2\xB7 um", &b.fibreDiameterUm, {10, 100, "%.0f"}, "The fibres' diameter; wool is about 18 to 40 um.");
    bodyFieldProblem(check, knit::YarnBodyField::FibreDiameter);
    if (maps) {
        dragRow("Map texel \xC2\xB7 um", &b.mapTexelUm, {10, 200, "%.0f"},
                "The maps' texel along and around the yarn: 40 um resolves the fibres at the article's size; a coarser texel costs less "
                "texture memory and bake time. Each yarn's atlas must fit the renderer's largest texture (4096 texels a side).");
        bodyFieldProblem(mapCheck, knit::YarnBodyField::MapTexel);
        percentRow("Cavity shading", &b.mapCavity, 0, 1, "How much the baked ambient occlusion darkens the valleys between the fibres.");
    }
    ImGui::EndTable();
}

// The viewport's maps: baked into the cache and shown on the tubes by a strand mesh job (startYarnMapBake);
// stale once the body or the texel changed (the shown scene's settings against the current ones). While
// the bake runs, its progress and Cancel are here (the Knit footer shows the same job).
void yarnMapsState(SceneUi& ui) {
    if (!ui.results.knit) return;
    if (ui.job.active && ui.job.bakeMaps) {
        JobProgressView v = progressView(ui.job, false);
        jobProgressBlock(v, "");
        mutedWrapped("Elapsed " + fmtSeconds(ui.job.progress.seconds) + kDot + "the tubes show the maps once every yarn is baked.");
        if (cancelButton(ui.job.cancelRequested ? "Cancel requested ..." : "Cancel the bake", ui.job.cancelRequested)) cancelJob(ui);
        return;
    }
    const bool shown = knit::carriesYarnMaps(*ui.results.knit);
    const bool current = shown && knit::yarnMapInputsText(ui.results.knit->settings) == knit::yarnMapInputsText(ui.s);
    if (current) {
        mutedWrapped("The yarn maps are baked: the viewport shows their cavity map on the tubes.");
        return;
    }
    ImGui::BeginDisabled(ui.job.active);
    if (ImGui::SmallButton(shown ? "Bake the yarn maps again" : "Bake the yarn maps")) startYarnMapBake(ui);
    ImGui::EndDisabled();
    tooltip(ui.job.active ? waitReason(ui).c_str()
                          : "Bakes every yarn's maps from its fibres (minutes once per knit and look, then cached) and shows their cavity "
                            "map on the tubes in the viewport. The render preview and the export bake them too when they are missing.");
}

void structureNode(SceneUi& ui, const knit::YarnBodyCheck& check, bool maps) {
    knit::SceneSettings::YarnBody& b = ui.s.yarnBody;
    char summary[96];
    if (maps) std::snprintf(summary, sizeof summary, "shell %.0f %%", b.shellDepth * 100);
    else std::snprintf(summary, sizeof summary, "shell %.0f %%, %s core", b.shellDepth * 100, b.core.c_str());
    if (ui.panels.revealYarnFibres == Reveal::Structure) ImGui::SetNextItemOpen(true);
    if (!detailsNode("Fibre structure", summary)) return;
    revealHere(ui, Reveal::Structure);
    if (beginRows("##yarnBodyStructure")) {
        percentRow("Shell depth", &b.shellDepth, 0.05, 1, "How deep the fibres reach into a ply, as a share of its radius; 100 % fills it.");
        bodyFieldProblem(check, knit::YarnBodyField::ShellDepth);
        dragRow("Fibre angle", &b.fibreAngleDeg, {0, 60, "%.0f deg"}, "The fibres' twist against their ply's axis.");
        bodyFieldProblem(check, knit::YarnBodyField::FibreAngle);
        ImGui::BeginDisabled(b.plies <= 1);
        dragRow("Ply angle", &b.plyAngleDeg, {0, 60, "%.0f deg"}, "The plies' twist against the yarn's axis.");
        ImGui::EndDisabled();
        if (b.plies <= 1) tooltip("One ply: there is nothing to twist around the yarn's axis.");
        dragRow("Interior fibres", &b.interiorFibres, 0, 2000, "Fibres inside the shells, per yarn: they fill the plies (more to render).");
        bodyFieldProblem(check, knit::YarnBodyField::InteriorFibres);
        if (!maps) {  // the core proxy stands in for the fibres a real yarn has inside; the maps never see it
            int core = b.core == "opaque" ? 1 : b.core == "translucent" ? 2 : 0;
            const char* cores[] = {"None", "Opaque", "Translucent"};
            if (comboRow("Core", &core, cores, 3, "A surface inside each ply that stands for the fibres a real yarn has there: it closes the "
                                                  "gaps between the fibres. Opaque: plain diffuse; translucent: lets light through."))
                b.core = core == 1 ? "opaque" : core == 2 ? "translucent" : "none";
            bodyFieldProblem(check, knit::YarnBodyField::Core);
        }
        dragRow("Irregularity", &b.irregularity, {0, 2, "%.2f"}, "Fibre-scale variation: 0 regular fibres, 2 the most uneven radii, paths and ends.");
        bodyFieldProblem(check, knit::YarnBodyField::Irregularity);
        dragRow("Fibre seed", &b.seed, 1, 1000000, "Another seed gives other fibres of the same kind.");
        ImGui::EndTable();
    }
    ImGui::TreePop();
}

void appearanceNode(SceneUi& ui, bool maps) {
    knit::SceneSettings& s = ui.s;
    knit::SceneSettings::YarnBody& b = s.yarnBody;
    char summary[128];
    const bool shiny = b.mapMaterial == "shiny";
    if (maps)
        std::snprintf(summary, sizeof summary, "maps: %s, roughness %.2f, specular %.2f, anisotropy %.2f", shiny ? "shiny" : "fabric",
                      b.mapRoughness, b.mapSpecular, b.mapAnisotropy);
    else
        std::snprintf(summary, sizeof summary, "roughness %.2f / %.2f, IOR %.2f", s.yarnWhite.fibreLongRoughness, s.yarnWhite.fibreAzimRoughness,
                      s.yarnWhite.fibreIor);
    if (ui.panels.revealYarnFibres == Reveal::Appearance) ImGui::SetNextItemOpen(true);
    if (!detailsNode("Fibre appearance", summary)) return;
    revealHere(ui, Reveal::Appearance);
    if (beginRows("##yarnBodyAppearance")) {
        dragRow("Fibre roughness along", &s.yarnWhite.fibreLongRoughness, {0.01, 1, "%.2f"}, "The fibres' hair shading: the highlight's spread along a fibre.");
        dragRow("Fibre roughness around", &s.yarnWhite.fibreAzimRoughness, {0.01, 1, "%.2f"}, "The highlight's spread around a fibre.");
        dragRow("Fibre IOR", &s.yarnWhite.fibreIor, {1.3, 1.8, "%.2f"}, "The fibres' refractive index: lower gives a softer, less glassy shine.");
        if (maps) {
            int material = shiny ? 1 : 0;
            const char* materials[] = {"Fabric", "Shiny"};
            if (comboRow("Map material", &material, materials, 2,
                         "Fabric: a sheen in front of a soft highlight from the baked maps; lets the wool profile's translucency "
                         "through. Shiny: the highlight alone, brighter; opaque."))
                b.mapMaterial = material == 1 ? "shiny" : "fabric";
            dragRow("Map roughness", &b.mapRoughness, {0.01, 1, "%.2f"}, "The smooth fibre's roughness under the maps; the baked normal spread adds to it.");
            dragRow("Map specular", &b.mapSpecular, {0, 1, "%.2f"}, "The maps' highlight strength (0.35 is a fibre of refractive index 1.4).");
            dragRow("Highlight anisotropy", &b.mapAnisotropy, {0.05, 1, "%.2f"}, "The highlight's roughness along the yarn over across it: lower stretches it along the fibres.");
        }
        ImGui::EndTable();
    }
    ImGui::TreePop();
}

} // namespace

void yarnFibresSection(SceneUi& ui) {
    knit::SceneSettings& s = ui.s;
    const bool maps = knit::yarnBodyRenderOf(s) == knit::YarnBodyRender::Maps;
    const knit::YarnBodyCheck check = knit::checkYarnBody(s);
    const double longest = ui.results.knit && ui.results.knit->knit ? knit::longestYarnMm(*ui.results.knit->knit) : 0;
    const knit::YarnBodyCheck mapCheck = knit::checkYarnMaps(s, longest);
    const knit::FilePreset* match = knit::matchingFilePreset(knit::PresetGroup::YarnFibres, s);
    ImGui::PushID("yarnFibres");
    sectionHeading("Yarn fibres", match ? "" : "Custom");
    revealHere(ui, Reveal::Section);
    yarnFibresNotes(maps);
    yarnFibresRows(ui, check, mapCheck, maps);
    for (const std::string& err : knit::presetErrors(knit::PresetGroup::YarnFibres)) wrapped(kWarn, "Fibre look not loaded: " + err);
    if (!check.accepted) wrapped(kBad, "The export refuses these fibres: " + check.reason + ". " + check.remedy);
    if (!mapCheck.accepted) wrapped(kBad, "The export refuses these maps: " + mapCheck.reason + ". " + mapCheck.remedy);
    const YarnFibresCost cost = yarnFibresCost(ui);
    const std::string advice = yarnFibresMemoryAdvice(cost);
    wrapped(advice.empty() ? muted() : kWarn, yarnFibresCostLine(cost));
    if (!advice.empty()) wrapped(kWarn, advice);
    if (maps) yarnMapsState(ui);
    structureNode(ui, check, maps);
    appearanceNode(ui, maps);
    ui.panels.revealYarnFibres = Reveal::None;  // done this frame
    ImGui::PopID();
}

} // namespace viewer
