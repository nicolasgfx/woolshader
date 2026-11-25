// The Draft and Final profiles, the viewer's opening settings, and what this machine can compute.
#include "settings/settings_profiles.h"

#include "settings/settings_mapping.h"
#include "settings/file_presets.h"
#include "settings/thickness_presets.h"

namespace knit {

void applyPreviewProfile(SceneSettings& s) {
    s.settle = false;
    s.repair = false;
    s.diameterFallback.clear();
    s.diameterAttempts = 0;
    s.surfaceGate = s.surfaceHalf = false;
}

void applyFullProfile(SceneSettings& s, const SceneEntry& e) {
    SceneSettings d;
    applySceneDefaults(d, e);
    s.settle = d.settle;
    s.repair = d.repair;
    s.diameterFallback = d.diameterFallback;
    // A thickness preset's own fallback (scaled to its drawn D); Custom: the scene default's.
    if (const ThicknessPreset* p = matchingPreset(s, e)) s.diameterFallback = p->fallback;
    s.surfaceGate = d.surfaceGate;
    s.surfaceHalf = d.surfaceHalf;
    s.diameterAttempts = 1;
}

bool isPreviewProfile(const SceneSettings& s) {
    return !s.settle && !s.repair && s.diameterFallback.empty() && !s.surfaceGate && !s.surfaceHalf;
}

bool isFinalProfile(const SceneSettings& s, const SceneEntry& e) {
    SceneSettings full = s;
    applyFullProfile(full, e);
    return s.settle == full.settle && s.repair == full.repair && s.diameterFallback == full.diameterFallback &&
           s.surfaceGate == full.surfaceGate && s.surfaceHalf == full.surfaceHalf;
}

std::vector<std::string> computeNotes(const SceneSettings& s, const ComputeCaps& c) {
    std::vector<std::string> notes;
    if (!c.device) notes.push_back("No usable GPU (" + c.whyNoDevice + "): knitting and fibres cannot run on this machine.");
    else if (!c.ncp && (s.settle || s.repair)) notes.push_back("The knitting solver is unavailable (" + c.whyNoNcp + "): a Final knit cannot run on this machine.");
    return notes;
}

void applySceneFibreDefaults(SceneSettings& s, const SceneEntry& e) {
    if (e.defaultFibreDensity > 0) s.fibreDensity = e.defaultFibreDensity;
    if (const FilePreset* islands = presetWithId(PresetGroup::FibreIslands, e.defaultFibreIslands))
        applyPreset(PresetGroup::FibreIslands, s, *islands);
}

constexpr const char* kNewObjectHand = "handmade";
constexpr const char* kNewObjectIslands = "worn";

SceneSettings viewerDefaults(const SceneEntry& e) {
    SceneSettings s;
    s.scene = e.id;
    // The CLI's order (scene_cli.cpp): the scene's gauge and materials. No scene sets a workflow.
    applySceneDefaults(s, e);
    applyPreviewProfile(s);
    // New objects: the scene's wool preset, else the first one (presets/fibres, lowest order;
    // fibre_presets.h), if any loaded.
    const std::vector<FilePreset>& wool = presets(PresetGroup::Fibres);
    const FilePreset* use = presetWithId(PresetGroup::Fibres, e.defaultFibrePreset);
    if (!use && !wool.empty()) use = &wool.front();
    if (use) applyPreset(PresetGroup::Fibres, s, *use);
    // New objects look knitted by hand, with worn fuzz (presets/hand_knit, presets/fibre_islands).
    // SceneSettings{} stays the template, so knit_check without --viewer-defaults keeps it.
    if (const FilePreset* hand = presetWithId(PresetGroup::HandKnit, kNewObjectHand)) applyPreset(PresetGroup::HandKnit, s, *hand);
    if (const FilePreset* islands = presetWithId(PresetGroup::FibreIslands, kNewObjectIslands)) applyPreset(PresetGroup::FibreIslands, s, *islands);
    applySceneFibreDefaults(s, e);
    applySceneSettings(s, e);  // the scene's own settings win over the new objects' presets
    return s;
}

} // namespace knit
