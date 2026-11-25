#pragma once
// From settings to the numbers the stages use, in the settings' own terms: scene defaults, the yarn
// diameters to try, the gauge fit, the scene's metric mapping and the fibre settings' rules. The
// mapping to knitting's and fibres' types is scene/stage_params.h.
//
// Pure functions of a SceneSettings. Called by the stages, the export, the viewer and knit_check.
// fibreSettingsValid and tubeSettingsValid report a refusal; the rest cannot fail.
#include <string>
#include <vector>

#include "model/choices.h"
#include "model/coordinates.h"
#include "model/params.h"
#include "settings/scene_registry.h"
#include "settings/scene_settings.h"

namespace knit {

// Library mm of a scene about centreMm -> scene metres (Z-up; the exporter maps to Y-up).
SceneMapping sceneMetres(const Vec3& centreMm);

// The fibre settings' own rules (max angle below the fibre model's minimum lift, lengths out of order,
// a curl without a length, segments out of range): false with the reason. A preset is checked at load,
// the fibre stage before it builds its parameters (scene/stage_params.h).
bool fibreSettingsValid(const SceneSettings& s, std::string* err);
// pile.excluded_parts: the part names, in order.
std::vector<std::string> furExcludedParts(const SceneSettings& s);
// The Fur workflow's settings in their bounds (pile.*: every number finite and in its range, a loop wider than
// its fibre): false with the reason naming the key. Checked where the Fur stage builds its look.
bool pileSettingsValid(const SceneSettings& s, std::string* why);
// A texture profile (strand.texture_profile): "none" -> count 1; "<n>-swept" / "<n>-circles" ->
// n plies (2..8) of that model. False with a reason for anything else.
bool parseTextureProfile(const std::string& id, int& plies, PlyProfileModel& model, std::string* err);

// Library defaults a GLB scene starts from (measured flat references): D 1.10 at
// 5.5 x 3.4 (the largest flat PASS; 1.16 passes only on the teddy's curved pieces).
void applySceneDefaults(SceneSettings& s, const SceneEntry& e);
// The scene's own settings (SceneEntry::defaultSettings), applied after everything else that opens it: by
// applySceneDefaults, viewerDefaults and the viewer's scene switch. strand.stitch_size goes first, because every
// setting re-derives a yarn-sized stitch (applyStitchSize) and would replace a gauge set before the size mode. A
// key the table refuses is reported on stderr and skipped (a registry error; the test names it).
void applySceneSettings(SceneSettings& s, const SceneEntry& e);

// The yarn diameters a part tries, in order, on every layout: the nominal one (strand.diameter), then
// the fallbacks below it, until gate (a) passes (a compressed yarn knits each at its contact fraction
// k). knit.diameter_attempts cuts the list; the cut ones are in notAttempted (named in the part's
// note if no attempted diameter passes).
struct DiameterPlan {
    std::vector<double> attempt, notAttempted;
};
DiameterPlan partDiameters(const SceneSettings& s);

// knit.coverage = overlap, in the Knit workflow: under Fur there are no seams, and the parts are whole
// components (an overlapped chart part carries its whole component: fur on it would grow once per chart).
bool overlapCoverage(const SceneSettings& s);
// The strand tube's settings in their bounds (3..64 sides, 0.05..2 mm ring spacing); false with
// the reason otherwise (checked where the tube is built: the computation's boundary).
bool tubeSettingsValid(const SceneSettings& s, std::string* why);

// Gauge fit. A round yarn of diameter D with clearance c fits a stitch
// w x h when w >= kFitW (D + c) and h >= kFitH (D + c): the calibrated teddy (1.16 mm at
// 5.5 x 3.4, c 0.1) passes after settle at 4.37 and 2.70; the static template needs 4.59 and
// 2.76. A tighter gauge has no round-tube solution, and knitting it only runs settle
// rounds that cannot succeed. With strand.compress_to_fit the yarn is knitted (solver, gate a)
// at the contact diameter k D that fits, k >= kMinCompression, and drawn at D: contact
// flattening with inflation (1 + e) / k - 1. Real wool flattens where yarns cross (the
// literature relaxes with soft contacts).
constexpr double kFitW = 4.35, kFitH = 2.68, kMinCompression = 0.6;
// The fit limits of the stockinette stitch. They hold for every pattern: the
// settings' w x h is the stockinette stitch of the yarn, and a pattern knits its own relaxed size of it
// (patternStitch below). The yarn's fit (yarnFit) is therefore the same for every pattern.
struct FitLimits {
    double w = kFitW, h = kFitH;  // [D + c]
};
FitLimits fitLimits();
// The relaxed stitch of a pattern relative to stockinette's, measured on machine-knitted 40 x 48 stitch
// patches (Hwang, Lin, Hsu, Mastripolito, McCann, Yuksel 2025, "Neighbor-Aware Data-Driven Relaxation of
// Stitch Mesh Models", supplemental Table 2): rib draws in across the wales and
// garter shortens along them. Reverse stockinette is stockinette seen from the back (1 x 1).
struct PatternScale {
    double w = 1, h = 1;  // x the stockinette stitch
};
PatternScale patternRestScale(StitchPattern p);
// The stitch a chart knit of pattern p lays out [mm]: the settings' stockinette stitch s.wale x
// s.course scaled by patternRestScale(p). The teddy's spiral pieces knit stockinette (s.wale x s.course).
struct StitchSize {
    double w = 0, h = 0;  // [mm]
};
StitchSize patternStitch(const SceneSettings& s, StitchPattern p);
// The relaxed shape of pattern p on pp: pattern, the stitch w x h (patternStitch)
// and, for a rib, its two layers: loops at stockinette's width (head and leg half widths / the width
// scale; the wales overlap in projection), knit wales kRibLayerOffset in front of the purl wales, and
// for rib 2x2 the pairs spread to stockinette's wale spacing. pp's other fields are the caller's.
void applyPatternShape(PatternParams& pp, const SceneSettings& s, StitchPattern p);
// Rib layer offsets [D] (knit_check flat 8 x 8 at 5.5 x 3.4 / D 1.164 and the rib sizes): rib 1x1
// passes from 1.25 (1.3: the fabric is then 2 x stockinette's thickness, Spencer); rib 2x2's
// pair transitions pass the fold check from 2.5.
constexpr double kRibLayerOffset1x1 = 1.3, kRibLayerOffset2x2 = 2.5;
// Rib 1x1's loop width at most 1.25 x the default head and leg (head side 0.489 < the sinker's 0.5).
constexpr double kRibMaxWiden = 1.25;
// Garter ridges [D] and transition (PatternParams::garterWave, transNeckY / Z, transLegY): the best
// flat margins of a scan over the wave, the depth and the transition at the garter sizes.
constexpr double kGarterWave = 0.3, kGarterNeckY = 0.55, kGarterNeckZ = 0.0, kGarterLegY = 0.30;
// Stitch size from the yarn (strand.stitch_size = yarn): w = f_W (D + c), h = f_H (D + c), rounded
// up to 0.05 mm, with f = fitLimits() x the tightness factor: tight 1 (the fit limit), standard
// the former Medium preset's ratio to the stockinette limit (6.5 x 4.6 at 1.0 mm), loose the former
// Thin preset's (7.5 x 5.5 at 0.8 mm; the presets before 2026-10-28). No tightness ever compresses the yarn (yarnFit need >= 1).
// strand.stitch_size = yarn (above): writes the derived w and h into s.wale / s.course (custom: unchanged); setSetting
// and readSettings call it, the viewer after its edits.
void applyStitchSize(SceneSettings& s);
struct YarnFit {
    double need = 1;       // the contact fraction the gauge needs (>= 1: a round yarn fits)
    double k = 1;          // knit (contact) diameter = k D
    bool tooTight = false; // the knit diameter still does not fit (compression off, or need < kMin)
    bool compressed() const { return k < 1; }
};
YarnFit yarnFit(const SceneSettings& s);
double knitDiameter(const SceneSettings& s, double D);  // k D
bool renderFlattened(const SceneSettings& s);           // contact flattening: strand.flatten or compressed
double renderInflate(const SceneSettings& s);           // e of the rendered D_r = D_knit (1 + e)

} // namespace knit
