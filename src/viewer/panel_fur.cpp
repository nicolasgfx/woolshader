// The Fur step (workflow = fur): where the pile grows, its loops and their variation, the fine fibres over them,
// and how they look. The names follow the Redshift "RS Procedural Fur" panel the owner works from (Density,
// Length, Thickness, Taper, Bend, Curl, Gravity and their variations); their meanings are ours
// (fibres/fur_look.h, settings pile.*). The job and its result are in the shared footer (fibreFooter).
#include <algorithm>
#include <string>
#include <vector>

#include "imgui.h"
#include "settings/settings_mapping.h"
#include "viewer/interface.h"
#include "viewer/workflow.h"

namespace viewer {

namespace {

// Coverage: the whole surface or a mask image, and the parts that grow fur.
void furCoverageSection(SceneUi& ui) {
    knit::SceneSettings::Pile& p = ui.s.pile;
    ImGui::SeparatorText("Coverage");
    int cov = p.coverage == "mask" ? 1 : 0;
    ImGui::RadioButton("Entire surface", &cov, 0);
    ImGui::SameLine();
    ImGui::RadioButton("Mask image", &cov, 1);
    tooltip("A greyscale image in the object's colour texture layout: white grows full density, black none.");
    p.coverage = cov == 1 ? "mask" : "entire";
    if (cov == 1) {
        mutedWrapped(p.maskPath.empty() ? std::string("No mask image loaded.") : p.maskPath);
        if (ImGui::Button("Load ...") && ui.platform.dialogs) ui.platform.request = PlatformRequest::OpenMaskImage;  // then loadMaskImage
        if (!ui.platform.dialogs) tooltip("The file dialog is not available here.");
        ImGui::SameLine();
        ImGui::BeginDisabled(p.maskPath.empty());
        if (ImGui::Button("Reload")) reloadMaskImage(ui);
        tooltip("Read the image again after editing it (the fur is then out of date).");
        ImGui::SameLine();
        if (ImGui::Button("Clear")) {
            p.maskPath.clear();
            p.maskHash.clear();
        }
        ImGui::EndDisabled();
        if (checkMaskImage(ui)) wrapped(kWarn, "The image changed since it was loaded: Reload it.");
    }
    // The parts are known once the fur has grown on the model.
    const knit::FibreResult* fr = fibreResult(ui);
    if (fr && !fr->baseParts.empty()) {
        if (ImGui::TreeNode("Parts")) {
            std::vector<std::string> excluded = knit::furExcludedParts(ui.s);
            bool changed = false;
            for (const std::string& part : fr->baseParts) {
                bool on = std::find(excluded.begin(), excluded.end(), part) == excluded.end();
                if (ImGui::Checkbox(part.c_str(), &on)) {
                    changed = true;
                    if (on) excluded.erase(std::remove(excluded.begin(), excluded.end(), part), excluded.end());
                    else excluded.push_back(part);
                }
            }
            if (changed) {
                std::string t;
                for (const std::string& n : excluded) t += (t.empty() ? "" : ",") + n;
                p.excludedParts = t;
            }
            ImGui::TreePop();
        }
    } else {
        mutedWrapped("Parts can be switched off once the fur has grown on the model.");
    }
    ImGui::Checkbox("Show coverage on the model", &ui.view.coverage);
    tooltip("Tints the areas without fur red (mask below one half, switched-off parts). Model view, after the fur has grown.");
}

// Loops: the thick fibres, each curling round once in its own plane, and their variation.
void furLoopsSection(SceneUi& ui) {
    knit::SceneSettings::Pile& p = ui.s.pile;
    ImGui::SeparatorText("Loops");
    if (beginRows("##furLoops")) {
        dragRow("Density", &p.density, {1, 2000, "%.0f / cm\xC2\xB2"}, "Loop fibres per square centimetre of surface, before the coverage.");
        if (p.curl > 0)
            dragRow("Length", &p.length, {0.05, 10, "%.2f turns"},
                    "How far each fibre curls round, in multiples of one loop of its own size (1: once round, ending beside its root).");
        else
            dragRow("Length##straight", &p.straightLength, {0.1, 50, "%.1f mm"}, "Each fibre's length (Curl 0: straight fibres).");
        dragRow("Thickness", &p.thickness, {0.01, 3, "%.2f mm"}, "Each fibre's diameter at its root.");
        percentRow("Taper", &p.taper, 0.05, 1, "The tip's diameter as a share of the root's.");
        dragRow("Bend", &p.bend, {0.1, 4, "%.2f"}, "How soon a fibre curls fully: 1 after a quarter turn, 2 after half a turn.");
        dragRow("Curl", &p.curl, {0, 10, "%.2f mm"}, "The loop's radius (0: straight fibres).");
        dragRow("Gravity", &p.gravity, {0, 50, "%.2f m/s\xC2\xB2"},
                "Pulls the settled fibres down; fibres shorter than 30 mm take less of it (they stand), a styling rule.");
        ImGui::EndTable();
    }
    if (ImGui::TreeNodeEx("Variation", ImGuiTreeNodeFlags_DefaultOpen)) {
        if (beginRows("##furVariation")) {
            percentRow("Length Var.", &p.varLength, 0, 1, "Each fibre loses up to this share of its length (uniform).");
            percentRow("Thickness Var.", &p.varThickness, 0, 1, "Each fibre's diameter varies by up to this share either way (uniform).");
            percentRow("Direction Var.", &p.varDirection, 0, 1,
                       "100 %: each fibre leaves the surface 0 to 80 deg up, in any direction around the normal; 0 %: straight up.");
            percentRow("Bend Var.", &p.varBend, 0, 1, "How soon each fibre curls varies by up to this share either way.");
            percentRow("Curl Var.", &p.varCurl, 0, 1, "Each loop's radius varies by up to half this share either way.");
            ImGui::EndTable();
        }
        ImGui::TreePop();
    }
}

// Fibres: the fine straight fibres standing over the loops (the silhouette's halo).
void furFibresSection(SceneUi& ui) {
    knit::SceneSettings::Pile& p = ui.s.pile;
    ImGui::SeparatorText("Fibres");
    ImGui::Checkbox("Fine fibres over the loops", &p.fibresEnabled);
    if (!p.fibresEnabled) return;
    if (beginRows("##furFibres")) {
        dragRow("Density##fibres", &p.fibresDensity, {0, 2000, "%.0f / cm\xC2\xB2"}, "Fine fibres per square centimetre of surface.");
        dragRow("Length##fibres", &p.fibresLength, {0.1, 50, "%.1f mm"}, "The typical length; each fibre lies between half and twice it.");
        dragRow("Length spread", &p.fibresLengthSpread, {0, 2, "%.2f"}, "The spread of the log-normal lengths (0: all the same).");
        dragRow("Thickness##fibres", &p.fibresThickness, {0.005, 1, "%.3f mm"}, "Each fine fibre's diameter at its root.");
        percentRow("Taper##fibres", &p.fibresTaper, 0.05, 1, "The tip's diameter as a share of the root's.");
        dragRow("Lift", &p.fibresLiftDeg, {0, 90, "%.0f deg"}, "How steeply the fine fibres leave the surface (90: straight up).");
        percentRow("Direction Var.##fibres", &p.fibresDirection, 0, 1, "100 %: lift +-40 deg, in any direction around the normal.");
        ImGui::EndTable();
    }
}

// Islands: patches of matted, fluffy or thin pile, placed as the wool islands are (fibre_islands.h).
void furIslandsSection(SceneUi& ui) {
    knit::SceneSettings::Pile& p = ui.s.pile;
    ImGui::SeparatorText("Islands");
    if (beginRows("##furIslands")) {
        dragRow("Rate", &p.islandsRate, {0, 10, "%.2f / cm\xC2\xB2"}, "Islands per square centimetre of surface (0: none).");
        dragRow("Size##islands", &p.islandsSize, {1, 100, "%.0f mm"}, "An island's diameter (each varies 0.7 to 1.3 times it).");
        percentRow("Strength##islands", &p.islandsStrength, 0, 1, "How strongly an island changes the pile at its centre (it fades to its edge).");
        dragRow("Matted", &p.islandsMatted, {0, 10, "%.1f"}, "Share of matted islands: loops flatter, shorter and clumped tighter; fine fibres lying down.");
        dragRow("Fluffy", &p.islandsFluffy, {0, 10, "%.1f"}, "Share of fluffy islands: loops more open, more fine fibres.");
        dragRow("Thin", &p.islandsThin, {0, 10, "%.1f"}, "Share of thin islands: both layers sparser.");
        ImGui::EndTable();
    }
}

// Irregularities: variation neighbouring fibres share.
void furIrregularitiesSection(SceneUi& ui) {
    knit::SceneSettings::Pile& p = ui.s.pile;
    ImGui::SeparatorText("Irregularities");
    if (beginRows("##furIrregular")) {
        percentRow("Density Var.", &p.varDensity, 0, 1, "Patches of denser and sparser pile: +- this share of the density.");
        dragRow("Density Var. Scale", &p.varDensityScale, {0.5, 500, "%.1f mm"}, "The patches' size.");
        percentRow("Patch length", &p.patchLength, 0, 1, "The same patches on the loops' length: +- this share.");
        percentRow("Patch lift", &p.patchLift, 0, 1, "The same patches on how steeply the loops stand: +- this share.");
        dragRow("Clump size", &p.clumpSize, {0.5, 20, "%.1f mm"}, "The reach of a nub: loops within it lean together.");
        percentRow("Clump strength", &p.clumpStrength, 0, 1, "How far the loops of a nub lean toward its direction (0: no nubs).");
        percentRow("Wear", &p.wear, 0, 1, "Exposed, convex places flatter, shorter and sparser; creases fuller (a curvature estimate).");
        ImGui::EndTable();
    }
}

// Appearance: the hair material in one colour (the model's own by default), or each fibre the texture's colour at
// its root.
void furAppearanceSection(SceneUi& ui) {
    knit::SceneSettings& s = ui.s;
    ImGui::SeparatorText("Appearance");
    int src = s.pileColourSource == "texture" ? 0 : s.pileColourSource == "object" ? 1 : 2;
    ImGui::RadioButton("Colour from the object", &src, 1);
    ImGui::SameLine();
    ImGui::RadioButton("One colour", &src, 2);
    ImGui::RadioButton("Colour from the texture", &src, 0);
    s.pileColourSource = src == 0 ? "texture" : src == 1 ? "object" : "single";
    if (src == 0) mutedWrapped("Every fibre takes the texture colour at its root, in at most 64 colours. A part without a texture has one colour.");
    if (src == 1) mutedWrapped("One colour per part: its source material's stored colour (with a texture, its average).");
    if (src == 2) colorEdit("Fur colour", s.pileMaterial.color);
    if (beginRows("##furShading")) {
        dragRow("Roughness along", &s.pileMaterial.fibreLongRoughness, {0.01, 1, "%.2f"}, "The hair material's longitudinal roughness.");
        dragRow("Roughness around", &s.pileMaterial.fibreAzimRoughness, {0.01, 1, "%.2f"}, "The hair material's azimuthal roughness.");
        dragRow("Index of refraction", &s.pileMaterial.fibreIor, {1.0, 2.5, "%.2f"});
        ImGui::EndTable();
    }
}

// View: the disc preview's colouring and what it shows (handover section 12: an honest answer per control).
void furViewSection(SceneUi& ui) {
    ImGui::SeparatorText("View");
    if (beginRows("##furView")) {
        const char* items[] = {"Plain", "Turns", "Thickness", "Tilt", "Layer"};
        int c = static_cast<int>(ui.view.discColour);
        if (comboRow("Colour by", &c, items, 5,
                     "Plain: the fur's colour. Turns: each loop's length share (blue short, red full). Thickness: its diameter. Tilt: its "
                     "start direction as RGB. Layer: loops and fine fibres."))
            ui.view.discColour = static_cast<knit::FurDiscColour>(c);
        ImGui::EndTable();
    }
    const std::shared_ptr<const knit::FurPreview>& p = ui.previews.discShown;
    if (ui.view.fibres != FibreView::Discs) {
        mutedWrapped("Choose Discs (Preview target, above the view) for the live disc preview.");
        return;
    }
    if (p) {
        char b[300];
        std::snprintf(b, sizeof b, "%s loop discs, %s fine fibres%s; built in %.0f ms (fields %.0f, buffers %.0f%s).", fmtCount(p->loops).c_str(),
                      fmtCount(p->fibres).c_str(), p->keyBelow < 1 ? " (a quarter, while editing)" : "", p->msField + p->msBuffers + p->msRoots + p->msScene,
                      p->msField, p->msBuffers, p->msRoots > 0 ? (", roots " + std::to_string(int(p->msRoots)) + " ms").c_str() : "");
        mutedWrapped(b);
        mutedWrapped(p->note);
    }
    mutedWrapped("Live: density, curl, length, thickness, direction and their variations. After Grow fur only: gravity, bend's exact shape, "
                 "the settle against the body and the fibres rejected there.");
}

void furAdvancedSection(knit::SceneSettings& s) {
    if (!ImGui::TreeNode("Advanced")) return;
    ImGui::InputInt("Seed", &s.pile.seed);
    long long budget = s.segmentBudget;
    if (ImGui::InputScalar("Segment budget", ImGuiDataType_S64, &budget)) s.segmentBudget = std::max(0LL, budget);
    tooltip("The most fibre segments the export may hold (whole fibres beyond it are dropped before the settle).");
    ImGui::TreePop();
}

} // namespace

void furPanel(SceneUi& ui) {
    ImGui::PushID("fur");
    stepHeading("Fur", std::string("Loops and fine fibres on the mesh") + kDot + "applies to the whole object");
    std::string why;
    if (!knit::pileSettingsValid(ui.s, &why)) wrapped(kWarn, "These settings cannot grow fur: " + why + ".");
    furCoverageSection(ui);
    furLoopsSection(ui);
    furFibresSection(ui);
    furIslandsSection(ui);
    furIrregularitiesSection(ui);
    furAppearanceSection(ui);
    furViewSection(ui);
    furAdvancedSection(ui.s);
    ImGui::PopID();
}

} // namespace viewer
