#pragma once
// Thickness presets of the Knit step: yarn diameter D and stitch size w x h together. The
// flat-swatch reference (flatPatternReference) is asserted as an expected verdict matrix (tests:
// thickness_presets_flat_matrix); a curved object's own verdicts come from its knit, and one never
// stands in for the other.
//   every scene (imports, examples and the teddy's spiral pieces): the teddy's calibrated stitch
//     5.5 x 3.4 with three drawn yarns, Dense 1.86 and Medium 1.50 (knitted compressed at the
//     contact diameter, yarnFit: the same knit, a chunkier drawn yarn; Dense matches the 1.6x
//     inflated close-up study of 2026-10-28) and Thin 1.16 (the calibrated round yarn: 14 of 18
//     teddy pieces PASS after settle and repair). The flat references are
//     knitted at the contact diameter; at the curved fit limit they FAIL the clearance.
// Editing D, w or h makes the settings Custom.
#include <string>
#include <vector>

#include "settings/scene_registry.h"
#include "settings/scene_settings.h"

namespace knit {

struct ThicknessPreset {
    std::string id, label;    // id: strand.thickness_preset ("thick-dense", "medium", "thin")
    double D = 0, w = 0, h = 0;
    std::string note;         // what it is for, and its flat-swatch evidence
    // strand.diameter_fallback in the Final profile: the calibrated 1.12 / 1.08 / 1.04 / 1.0 scaled
    // by D / 1.16 (knitted at about those diameters whatever the drawn D).
    std::vector<double> fallback;
};

// The presets of a scene (its family); empty ids are never listed.
std::vector<ThicknessPreset> thicknessPresets(const SceneEntry& e);
// The preset whose values the settings hold (nullptr: Custom).
const ThicknessPreset* matchingPreset(const SceneSettings& s, const SceneEntry& e);
// Sets D, w, h and, unless the settings are Draft (no fallback), the preset's fallback.
void applyThicknessPreset(SceneSettings& s, const ThicknessPreset& p);

// The Knit step's one Gauge list: the thickness presets (a fixed stitch and yarn), then the stitch
// sized for the yarn (strand.stitch_size = yarn at a tightness: w and h follow D), then Custom. It
// only maps onto the existing keys strand.thickness_preset, strand.stitch_size and strand.tightness.
struct GaugeChoice {
    std::string id;     // a preset's id, "yarn-tight", "yarn-standard", "yarn-loose" or "custom"
    std::string label;  // shown in the list
    std::string note;   // its tooltip
    std::string group;  // non-empty: a caption drawn before it (the first entry of a group)
};
std::vector<GaugeChoice> gaugeChoices(const SceneEntry& e);
// The id of the entry the settings hold: yarn mode by its tightness, else the matching preset, else Custom.
std::string currentGaugeChoice(const SceneSettings& s, const SceneEntry& e);
// Selects the entry: a preset sets its values, a yarn entry derives w and h from D, Custom keeps the values.
void applyGaugeChoice(SceneSettings& s, const SceneEntry& e, const std::string& id);
// Gauge fit (yarnFit, scene_settings.h), every scene. The note when the yarn is knitted compressed
// (empty: none), and the warning when even the knit diameter does not fit (compression off, ply
// tubes, or below kMinCompression): its parts are expected to fail.
std::string yarnFitNote(const SceneSettings& s);
std::string thickYarnWarning(const SceneSettings& s, const SceneEntry& e);

} // namespace knit
