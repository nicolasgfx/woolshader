// SceneSettings: the field table that drives write, read, diff and the stage keys (scene_settings.h).
#include "settings/scene_settings.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <optional>
#include <sstream>
#include <utility>

#include "model/env.h"
#include "settings/file_presets.h"
#include "settings/settings_mapping.h"

namespace knit {

const char* toString(Stage s) {
    switch (s) {
    case Stage::Load: return "load";
    case Stage::Knit: return "knit";
    case Stage::Surface: return "surface";
    case Stage::Fibres: return "fibres";
    case Stage::Export: return "export";
    case Stage::Display: return "display";
    }
    return "?";
}

const char* toString(Workflow w) { return w == Workflow::Fur ? "fur" : "knit"; }

Workflow workflowOf(const SceneSettings& s) { return s.workflow == "fur" ? Workflow::Fur : Workflow::Knit; }

namespace {

// The Fur switch: -1 until the environment is read, then 0 (off) or 1 (on). Atomic: jobs read settings on their
// worker threads.
std::atomic<int> furSwitch{-1};

// Why workflow = fur is refused while the switch is off; the same words for a file and for --set.
const char* const kFurSwitchedOff = "workflow = fur: the Fur workflow is switched off (KNIT_FUR=1 switches it on)";

}  // namespace

bool furEnabled() {
    int on = furSwitch.load();
    if (on < 0) {
        const std::optional<std::string> env = envVar("KNIT_FUR");
        on = env ? (*env == "0" ? 0 : 1) : (kFurByDefault ? 1 : 0);  // any value but 0 switches it on
        furSwitch.store(on);
    }
    return on == 1;
}

void setFurEnabled(bool on) { furSwitch.store(on ? 1 : 0); }

SceneSettings::SceneSettings() {
    // Directed fibres: pieceFibres() (fibre.cpp).
    directed.curvature = 1.0;      // bowTheta 0.10 rad, bowRadial 0.03 mm at scale 1
    directed.maxAngleDeg = 24;     // liftMaxDeg
    random = directed;
    random.curvature = 0.3;        // randomBend
    // The strays' deviation from the yarn: the Loose stray preset (presets/strays/loose.preset), with
    // randomSpreadDeg 45 in the header. Older settings files keep theirs (readSettings).
    random.maxAngleDeg = 40;       // randomElevMaxDeg: the strays rise 20 to 40 deg
    random.lengthSpread = 0.5;
    yarnWhite.color = kDefaultYarnColour;
    yarnBrown.color = {0.42, 0.20, 0.08};
    fibreDirected = yarnWhite;
    yarnWhite.sheenTint = yarnBrown.sheenTint = kYarnSheenTint;  // after the copy: the fuzz keeps its own
    fibreRandom.model = MaterialSpec::Model::Fibre;
    fibreRandom.color = {0.55, 0.55, 0.55};  // with fibreRandomInheritsYarn off
    stuffing.color = {0.90, 0.90, 0.88};
    // The pile's hair material (the approved look: the hair model's default roughness 0.3 / 0.3, IOR 1.55,
    // cuticle 2 deg); its colour is the model's own (pileColourSource object) unless one is picked.
    pileMaterial.model = MaterialSpec::Model::Fibre;
    pileMaterial.color = {0.72, 0.45, 0.22};
}

double SceneSettings::Population::thicknessMin() const { return thickness * std::exp(-2 * thicknessSigma); }
double SceneSettings::Population::thicknessMax() const { return thickness * std::exp(2 * thicknessSigma); }
void SceneSettings::Population::setThicknessRange(double mn, double mx) {
    if (!(mn > 0) || !(mx >= mn)) return;
    thickness = std::sqrt(mn * mx);
    thicknessSigma = std::log(mx / mn) / 4;
}

namespace {

constexpr const char* kVersionLine = "knit-settings v1";

// Keys the program no longer uses, removed with what read them. Saved files may still hold them: they
// are accepted and ignored, never written, listed or a value of a stage (a few keep their line in a
// stage's text: kRetiredInputs, result_keys.cpp).
struct RetiredKey {
    const char* key;
    const char* lastDefault;  // the value every saved file holds, not reported; "": every value is reported
};
constexpr RetiredKey kRetiredKeys[] = {
    {"strand.contact", "clearance"},  // 2026-11-01: the compress contact is gone
    {"strand.compression", "0.1"},
    {"strand.core_clearance", "0"},
    {"strand.template", "fixed"},  // 2026-11-01: the relaxed stitch cell is gone
    {"strand.cell_variant", "L"},
    {"strand.cell_tightness", "1"},
    {"knit.layout", "uv-tiling"},  // 2026-11-01: the stitch-mesh layout is gone
    {"strand.teddy_stitches", ""},  // 2026-10-28: the teddy's scale is SceneEntry::gaugeStitches
    {"knit.backend", ""},  // 2026-10-28: scenes knit on the GPU only
    {"knit.solver", ""},  // 2026-10-30: the explicit problem is the one solver
    {"fibres.backend", ""},  // 2026-10-28: scene fibres grow on the GPU
    {"fibres.placement", ""},  // 2026-11-01: every fibre placed all around the yarn
    {"fibres.exposed.widen", ""},
    {"fibres.exposed.cell", ""},
    {"display.ui_scale", ""},  // 2026-10-28: the viewer's scale is ViewerPrefs::uiScale
    {"knit.lattice_anchor", "bbox"},  // 2026-11-01: block gauges are gone
    {"import.charts", "keep"},  // 2026-11-11: knitting charts are never generated on the file's mesh
    // 2026-11-12: the import has one engine, and every wrapped import is rebuilt by its stages. "builtin" asked for
    // those already, so only another engine's name or a recipe is reported.
    {"import.wrap_engine", "builtin"},
    {"import.wrap_recipe", ""},
    {"fibres.wool", "true"},  // 2026-11-01: scalp hair is gone, and with it fibres without wool
    {"hair.enabled", "false"}, {"hair.density", "1.5"}, {"hair.guide_density", "0.05"}, {"hair.root_diameter", "0.07"},
    {"hair.tip_ratio", "0.8"}, {"hair.segment_length", "1.5"}, {"hair.crown", "25"}, {"hair.volume", "30"},
    {"hair.tip_volume", "14"}, {"hair.cut_below_top", "88"}, {"hair.cut_jitter", "2"}, {"hair.lift_near_deg", "35"},
    {"hair.lift_far_deg", "10"}, {"hair.part_wiggle", "1"}, {"hair.clump", "0.55"}, {"hair.wave_amplitude", "1.5"},
    {"hair.frizz", "0.25"}, {"hair.stray_fraction", "0.01"}, {"hair.tones", "3"}, {"hair.tone_spread", "0.04"},
    {"hair.settle_steps", "0"}, {"hair.gravity", "1"}, {"hair.segment_budget", "8000000"}, {"hair.seed", "0"},
    {"material.hair.model", "fibre"}, {"material.hair.color", "0.20000000000000001, 0.5, 0.20000000000000001"},
    {"material.hair.sheen", "1"}, {"material.hair.sheen_tint", "0.5"}, {"material.hair.sheen_roughness", "0.5"},
    {"material.hair.clearcoat", "false"}, {"material.hair.clearcoat_weight", "0.3"}, {"material.hair.clearcoat_ior", "1.5"},
    {"material.hair.clearcoat_roughness", "0.3"}, {"material.hair.fibre_longitudinal_roughness", "0.3"},
    {"material.hair.fibre_azimuthal_roughness", "0.3"}, {"material.hair.fibre_ior", "1.55"},
    {"material.hair.fibre_cuticle_tilt_deg", "2"},
    // 2026-11-19: the Hair mode (fibres.profile = fur) is replaced by the Fur workflow (pile.*). A knitted file
    // holds these at their defaults; knit.enabled and fibres.profile are read by legacyWorkflow.
    {"fur.coverage", "entire"}, {"fur.mask_path", ""}, {"fur.mask_hash", ""}, {"fur.excluded_parts", ""},
    {"fur.density", "1.5"}, {"fur.length_median", "8"}, {"fur.length_spread", "0.25"}, {"fur.length_min", "3"},
    {"fur.length_max", "15"}, {"fur.root_diameter", "0.08"}, {"fur.taper", "0.3"}, {"fur.comb", "0, 0, -1"},
    {"fur.lift_deg", "25"}, {"fur.bend_stiffness", "0.75"}, {"fur.gravity", "0.15"}, {"fur.settle_steps", "60"},
    {"fur.tip_stiffness_ratio", "0.5"}, {"fur.collision_offset", "0.02"}, {"fur.lift_jitter_deg", "8"},
    {"fur.comb_jitter_deg", "15"}, {"fur.segment_length", "0.5"}, {"fur.seed", "0"},
    {"material.fur.model", "fibre"}, {"material.fur.color", "0.71999999999999997, 0.45000000000000001, 0.22"},
    {"material.fur.sheen", "1"}, {"material.fur.sheen_tint", "0.5"}, {"material.fur.sheen_roughness", "0.5"},
    {"material.fur.clearcoat", "false"}, {"material.fur.clearcoat_weight", "0.3"}, {"material.fur.clearcoat_ior", "1.5"},
    {"material.fur.clearcoat_roughness", "0.3"}, {"material.fur.fibre_longitudinal_roughness", "0.3"},
    {"material.fur.fibre_azimuthal_roughness", "0.3"}, {"material.fur.fibre_ior", "1.55"},
    {"material.fur.fibre_cuticle_tilt_deg", "2"}, {"material.fur.translucency", "0"}, {"material.fur_colour_source", "object"},
};
const RetiredKey* findRetiredKey(const std::string& key) {
    const auto it = std::find_if(std::begin(kRetiredKeys), std::end(kRetiredKeys), [&](const RetiredKey& k) { return key == k.key; });
    return it == std::end(kRetiredKeys) ? nullptr : it;
}

bool isRetiredSetting(const std::string& key) { return findRetiredKey(key) != nullptr; }

// Builds the table from an object's members: each overload fixes the kind, so a row cannot disagree
// with the member it names, and the member's place is what every other object's value is found by.
class TableBuilder {
public:
    explicit TableBuilder(SceneSettings& layout) : layout_(layout) {}
    void add(const std::string& key, Stage stage, double& member) { row(key, stage, SettingKind::Real, &member); }
    void add(const std::string& key, Stage stage, int& member) { row(key, stage, SettingKind::Integer, &member); }
    void add(const std::string& key, Stage stage, long long& member) { row(key, stage, SettingKind::LongInteger, &member); }
    void add(const std::string& key, Stage stage, bool& member) { row(key, stage, SettingKind::Flag, &member); }
    void add(const std::string& key, Stage stage, std::string& member) { row(key, stage, SettingKind::Text, &member); }
    void add(const std::string& key, Stage stage, Vec3& member) { row(key, stage, SettingKind::Vector3, &member); }
    void add(const std::string& key, Stage stage, std::vector<double>& member) { row(key, stage, SettingKind::RealList, &member); }
    void add(const std::string& key, Stage stage, MaterialSpec::Model& member) { row(key, stage, SettingKind::MaterialModel, &member); }
    // One of a fixed list of words (any other value is refused).
    void addChoice(const std::string& key, Stage stage, std::string& member, std::vector<std::string> words) {
        row(key, stage, SettingKind::Choice, &member);
        rows_.back().choices = std::move(words);
    }
    void addMaterial(const std::string& prefix, MaterialSpec& m) {
        const Stage E = Stage::Export;
        add(prefix + ".model", E, m.model);
        add(prefix + ".color", E, m.color);
        add(prefix + ".sheen", E, m.sheen);
        add(prefix + ".sheen_tint", E, m.sheenTint);
        add(prefix + ".sheen_roughness", E, m.sheenRoughness);
        add(prefix + ".clearcoat", E, m.clearcoat);
        add(prefix + ".clearcoat_weight", E, m.clearcoatWeight);
        add(prefix + ".clearcoat_ior", E, m.clearcoatIor);
        add(prefix + ".clearcoat_roughness", E, m.clearcoatRoughness);
        add(prefix + ".fibre_longitudinal_roughness", E, m.fibreLongRoughness);
        add(prefix + ".fibre_azimuthal_roughness", E, m.fibreAzimRoughness);
        add(prefix + ".fibre_ior", E, m.fibreIor);
        add(prefix + ".fibre_cuticle_tilt_deg", E, m.fibreCuticleDeg);
        add(prefix + ".translucency", E, m.translucency);
    }
    void addPopulation(const std::string& prefix, SceneSettings::Population& p) {
        const Stage F = Stage::Fibres;
        add(prefix + ".length_min", F, p.lengthMin);
        add(prefix + ".length_max", F, p.lengthMax);
        add(prefix + ".length_median", F, p.lengthMedian);
        add(prefix + ".length_spread", F, p.lengthSpread);
        add(prefix + ".thickness", F, p.thickness);
        add(prefix + ".thickness_sigma", F, p.thicknessSigma);
        add(prefix + ".curvature", F, p.curvature);
        add(prefix + ".max_angle_deg", F, p.maxAngleDeg);
    }
    std::vector<SettingField> take() { return std::move(rows_); }

private:
    void row(const std::string& key, Stage stage, SettingKind kind, const void* member) {
        const size_t offset = static_cast<size_t>(static_cast<const char*>(member) - reinterpret_cast<const char*>(&layout_));
        rows_.push_back({key, stage, kind, offset, {}});
    }
    SceneSettings& layout_;
    std::vector<SettingField> rows_;
};

// The table's order is part of the stage keys: a new key goes where its neighbours' keys stay valid.
std::vector<SettingField> makeTable() {
    SceneSettings s;  // only its members' places are used
    TableBuilder t(s);
    const Stage L = Stage::Load, K = Stage::Knit, S = Stage::Surface, F = Stage::Fibres, E = Stage::Export, D = Stage::Display;
    t.add("scene", L, s.scene);
    t.add("import.source", L, s.importSource);
    t.add("import.identity", L, s.importIdentity);
    t.add("import.unit_mm", L, s.importUnitMm);
    t.addChoice("import.up", L, s.importUp, {"y", "z"});
    t.add("import.repairs", L, s.importRepairs);
    t.addChoice("import.surface", L, s.importSurface, {"file", "wrapped"});
    t.add("import.wrap_stitch_mm", L, s.importWrapStitchMm);
    t.add("import.plain_nodes", L, s.importPlainNodes);
    t.addChoice("workflow", L, s.workflow, {"knit", "fur"});  // the former knit.enabled row (result_keys.cpp)
    t.addChoice("knit.route", K, s.knitRoute, {"surface", "guided"});
    t.add("knit.guided_piece", K, s.guidedPiece);
    t.add("knit.guided_component", K, s.guidedKnitComponent);
    t.add("knit.guided_markings", K, s.guidedMarkings);
    t.add("knit.guided_markings_digest", K, s.guidedMarkingsDigest);
    t.add("strand.thickness_preset", K, s.thicknessPreset);
    t.add("strand.diameter", K, s.diameter);
    t.add("strand.clearance", K, s.clearance);
    t.addChoice("knit.coverage", K, s.coverage, {"component", "overlap"});
    t.add("knit.overlap_stitches", K, s.overlapStitches);
    t.add("knit.gap_patches", K, s.gapPatches);
    t.add("strand.depth_amplitude", K, s.depthAmplitude);
    t.add("strand.wale", K, s.wale);
    t.add("strand.course", K, s.course);
    t.addChoice("strand.stitch_size", K, s.stitchSize, {"custom", "yarn"});
    t.addChoice("strand.tightness", K, s.tightness, {"tight", "standard", "loose"});
    t.add("strand.angle_deg", K, s.angleDeg);
    t.add("strand.pattern", K, s.pattern);
    t.add("strand.diameter_fallback", K, s.diameterFallback);
    t.add("strand.settle", K, s.settle);
    t.add("strand.repair", K, s.repair);
    t.add("knit.hand.seed", K, s.hand.seed);
    t.add("knit.hand.unevenness", K, s.hand.unevenness);
    t.add("knit.hand.column_shift", K, s.hand.columnShift);
    t.add("knit.hand.column_rows", K, s.hand.columnRows);
    t.add("knit.hand.loop_rate", K, s.hand.loopRate);
    t.add("knit.hand.loop_height", K, s.hand.loopHeight);
    t.add("knit.hand.wave_height", K, s.hand.waveHeight);
    t.add("knit.hand.spot_rate", K, s.hand.spotRate);
    t.add("knit.hand.spot_change", K, s.hand.spotChange);
    t.add("knit.hand.spot_length", K, s.hand.spotLength);
    t.add("strand.plies", S, s.plies);
    t.add("strand.ply_pitch", S, s.plyPitch);
    t.add("strand.ply_fill", S, s.plyFill);
    t.add("strand.ply_handedness", S, s.plyHandedness);
    t.add("strand.texture_profile", S, s.textureProfile);
    t.add("strand.flatten", K, s.flatten);
    t.add("strand.flatten_inflate", K, s.flattenInflate);
    t.add("strand.max_flattening", K, s.maxFlattening);
    t.add("strand.compress_to_fit", K, s.compressToFit);
    t.add("strand.surface_gate", K, s.surfaceGate);
    t.add("strand.surface_half", K, s.surfaceHalf);
    t.add("strand.tube_sides", F, s.tubeSides);
    t.add("strand.tube_ring", F, s.tubeRing);
    t.add("strand.ply_tube_sides", F, s.plyTubeSides);
    t.add("strand.ply_tube_ring", F, s.plyTubeRing);
    t.add("knit.diameter_attempts", K, s.diameterAttempts);
    t.add("fibres.enabled", F, s.fibres);
    t.add("fibres.density", F, s.fibreDensity);
    t.add("fibres.random_fraction", F, s.randomFraction);
    t.add("fibres.curvature_bias", F, s.curvatureBias);
    t.add("fibres.preset", F, s.fibrePreset);
    t.add("fibres.islands.rate", F, s.islands.rate);
    t.add("fibres.islands.size", F, s.islands.size);
    t.add("fibres.islands.strength", F, s.islands.strength);
    t.add("fibres.islands.fuzzy", F, s.islands.fuzzy);
    t.add("fibres.islands.matted", F, s.islands.matted);
    t.add("fibres.islands.worn", F, s.islands.worn);
    t.add("fibres.max_segments", F, s.fibreMaxSegments);
    t.add("fibres.segment_length", F, s.fibreSegmentLength);
    t.add("fibres.segment_budget", E, s.segmentBudget);
    t.add("fibres.tail_taper", E, s.tailTaper);
    t.add("fibres.tail_taper_min_diameters", E, s.tailTaperMinDiameters);
    t.add("fibres.tail_taper_tip", E, s.tailTaperTip);
    t.add("strand.ply_phase_jitter", E, s.plyPhaseJitter);
    t.add("strand.ply_pitch_jitter", E, s.plyPitchJitter);
    t.add("strand.ply_crown_tilt_deg", E, s.plyCrownTiltDeg);
    t.add("strand.normal_strength", E, s.normalStrength);
    t.addChoice("fibres.generation", F, s.fibreGeneration, {"budget", "full"});
    t.addChoice("fibres.checks", F, s.fibreChecks, {"light", "full"});
    t.add("fibres.side_band_cos", F, s.sideBandCos);
    t.add("fibres.crimp_amplitude", F, s.crimpAmplitude);
    t.add("fibres.crimp_wavelength", F, s.crimpWavelength);
    t.addChoice("pile.coverage", F, s.pile.coverage, {"entire", "mask"});
    t.add("pile.mask_path", F, s.pile.maskPath);
    t.add("pile.mask_hash", F, s.pile.maskHash);
    t.add("pile.excluded_parts", F, s.pile.excludedParts);
    t.add("pile.density", F, s.pile.density);
    t.add("pile.length", F, s.pile.length);
    t.add("pile.straight_length", F, s.pile.straightLength);
    t.add("pile.thickness", F, s.pile.thickness);
    t.add("pile.taper", F, s.pile.taper);
    t.add("pile.bend", F, s.pile.bend);
    t.add("pile.curl", F, s.pile.curl);
    t.add("pile.gravity", F, s.pile.gravity);
    t.add("pile.gravity_full_length", F, s.pile.gravityFullLength);
    t.add("pile.var.density", F, s.pile.varDensity);
    t.add("pile.var.density_scale", F, s.pile.varDensityScale);
    t.add("pile.var.length", F, s.pile.varLength);
    t.add("pile.var.thickness", F, s.pile.varThickness);
    t.add("pile.var.direction", F, s.pile.varDirection);
    t.add("pile.var.bend", F, s.pile.varBend);
    t.add("pile.var.curl", F, s.pile.varCurl);
    t.add("pile.fibres.enabled", F, s.pile.fibresEnabled);
    t.add("pile.fibres.density", F, s.pile.fibresDensity);
    t.add("pile.fibres.length", F, s.pile.fibresLength);
    t.add("pile.fibres.length_spread", F, s.pile.fibresLengthSpread);
    t.add("pile.fibres.thickness", F, s.pile.fibresThickness);
    t.add("pile.fibres.taper", F, s.pile.fibresTaper);
    t.add("pile.fibres.lift_deg", F, s.pile.fibresLiftDeg);
    t.add("pile.fibres.direction", F, s.pile.fibresDirection);
    t.add("pile.islands.rate", F, s.pile.islandsRate);
    t.add("pile.islands.size", F, s.pile.islandsSize);
    t.add("pile.islands.strength", F, s.pile.islandsStrength);
    t.add("pile.islands.matted", F, s.pile.islandsMatted);
    t.add("pile.islands.fluffy", F, s.pile.islandsFluffy);
    t.add("pile.islands.thin", F, s.pile.islandsThin);
    t.add("pile.irregular.clump_size", F, s.pile.clumpSize);
    t.add("pile.irregular.clump_strength", F, s.pile.clumpStrength);
    t.add("pile.irregular.patch_length", F, s.pile.patchLength);
    t.add("pile.irregular.patch_lift", F, s.pile.patchLift);
    t.add("pile.irregular.wear", F, s.pile.wear);
    t.add("pile.seed", F, s.pile.seed);
    t.addPopulation("fibres.directed", s.directed);
    t.addPopulation("fibres.random", s.random);
    t.add("fibres.random.curl_radius", F, s.random.curlRadius);
    t.add("fibres.random.curl_length", F, s.random.curlLength);
    t.add("fibres.directed.hug_height", F, s.directedHugHeight);
    t.add("fibres.directed.follow_plies", F, s.directedFollowPlies);
    t.add("fibres.directed.drape", F, s.directedDrape);
    t.add("fibres.random.spread_deg", F, s.randomSpreadDeg);
    t.addMaterial("material.yarn_white", s.yarnWhite);
    t.addMaterial("material.yarn_brown", s.yarnBrown);
    t.addMaterial("material.fibre_directed", s.fibreDirected);
    t.addMaterial("material.fibre_random", s.fibreRandom);
    t.addMaterial("material.stuffing", s.stuffing);
    t.addMaterial("material.pile", s.pileMaterial);
    t.addChoice("material.pile_colour_source", E, s.pileColourSource, {"object", "single", "texture"});
    t.add("material.fibre_inherits_yarn", E, s.fibreInheritsYarn);
    t.add("material.fibre_random_inherits_yarn", E, s.fibreRandomInheritsYarn);
    t.addChoice("material.fibre_hair_colour", E, s.fibreHairColour, {"layer", "yarn"});
    t.addChoice("material.yarn_colour_source", E, s.yarnColourSource, {"single", "object", "texture"});
    t.addChoice("material.plain_colour_source", E, s.plainColourSource, {"texture", "average", "material"});
    t.add("filling.enabled", E, s.filling);
    t.add("filling.margin", E, s.fillingMargin);
    t.add("yarn_body.enabled", E, s.yarnBody.enabled);
    t.add("yarn_body.plies", E, s.yarnBody.plies);
    t.add("yarn_body.fibres", E, s.yarnBody.fibres);
    t.add("yarn_body.interior_fibres", E, s.yarnBody.interiorFibres);
    t.add("yarn_body.fibre_diameter_um", E, s.yarnBody.fibreDiameterUm);
    t.add("yarn_body.shell_depth", E, s.yarnBody.shellDepth);
    t.add("yarn_body.fibre_angle_deg", E, s.yarnBody.fibreAngleDeg);
    t.add("yarn_body.ply_angle_deg", E, s.yarnBody.plyAngleDeg);
    t.add("yarn_body.irregularity", E, s.yarnBody.irregularity);
    t.addChoice("yarn_body.core", E, s.yarnBody.core, {"none", "opaque", "translucent"});
    t.add("yarn_body.seed", E, s.yarnBody.seed);
    t.addChoice("yarn_body.render", E, s.yarnBody.render, {"fibres", "maps"});
    t.add("yarn_body.map_texel_um", E, s.yarnBody.mapTexelUm);
    t.add("yarn_body.map_cavity", E, s.yarnBody.mapCavity);
    t.add("yarn_body.map_roughness", E, s.yarnBody.mapRoughness);
    t.add("yarn_body.map_specular", E, s.yarnBody.mapSpecular);
    t.add("yarn_body.map_anisotropy", E, s.yarnBody.mapAnisotropy);
    t.addChoice("yarn_body.map_material", E, s.yarnBody.mapMaterial, {"fabric", "shiny"});
    t.add("knit.time_cap_s", D, s.timeCapS);
    return t.take();
}

bool parseReal(const std::string& t, double& v) {
    char* end = nullptr;
    v = std::strtod(t.c_str(), &end);
    return end && end != t.c_str() && *end == 0 && std::isfinite(v);
}

bool parseInteger(const std::string& t, long long& v) {
    char* end = nullptr;
    v = std::strtoll(t.c_str(), &end, 10);
    return end && end != t.c_str() && *end == 0;
}

// Two values as text: equal, or numbers within 1e-9.
bool sameValueText(const std::string& x, const std::string& y) {
    if (x == y) return true;
    double u = 0, v = 0;
    return parseReal(x, u) && parseReal(y, v) && std::fabs(u - v) < 1e-9;
}

// A retired key's value worth a note: anything but the last default every saved file holds. An empty value
// chose nothing (import.wrap_recipe = , as older settings files hold it), so it is never reported.
bool isReportedRetiredValue(const std::string& key, const std::string& value) {
    const RetiredKey* k = findRetiredKey(key);
    if (value.empty()) return false;
    return !k || !*k->lastDefault || !sameValueText(value, k->lastDefault);
}

std::string withoutBlanks(std::string t) {
    t.erase(0, t.find_first_not_of(' '));
    t.erase(t.find_last_not_of(' ') + 1);
    return t;
}

bool parseVector3(const std::string& t, Vec3& v) {
    double c[3];
    size_t p = 0;
    for (int i = 0; i < 3; ++i) {
        const size_t e = i < 2 ? t.find(',', p) : t.size();
        if (e == std::string::npos) return false;
        if (!parseReal(withoutBlanks(t.substr(p, e - p)), c[i])) return false;
        p = e + 1;
    }
    v = {c[0], c[1], c[2]};
    return true;
}

bool parseRealList(const std::string& t, std::vector<double>& out) {
    std::vector<double> values;
    std::stringstream in(t);
    for (std::string part; std::getline(in, part, ',');) {
        part = withoutBlanks(part);
        if (part.empty()) continue;
        double v;
        if (!parseReal(part, v)) return false;
        values.push_back(v);
    }
    out = std::move(values);
    return true;
}

// The bytes of a row's member in one settings object. Each kind reads them as the type TableBuilder
// took the offset from.
char* memberBytes(SceneSettings& s, const SettingField& f) { return reinterpret_cast<char*>(&s) + f.offset; }
const char* memberBytes(const SceneSettings& s, const SettingField& f) { return reinterpret_cast<const char*>(&s) + f.offset; }

bool setFieldText(SceneSettings& s, const SettingField& f, const std::string& text) {
    char* m = memberBytes(s, f);
    switch (f.kind) {
    case SettingKind::Real: return parseReal(text, *reinterpret_cast<double*>(m));
    case SettingKind::Integer: {
        long long v;
        if (!parseInteger(text, v)) return false;
        *reinterpret_cast<int*>(m) = static_cast<int>(v);
        return true;
    }
    case SettingKind::LongInteger: return parseInteger(text, *reinterpret_cast<long long*>(m));
    case SettingKind::Flag:
        if (text != "true" && text != "false") return false;
        *reinterpret_cast<bool*>(m) = text == "true";
        return true;
    case SettingKind::Text:
        if (text.find_first_of("\r\n") != std::string::npos) return false;
        *reinterpret_cast<std::string*>(m) = text;
        return true;
    case SettingKind::Choice:
        if (std::find(f.choices.begin(), f.choices.end(), text) == f.choices.end()) return false;
        *reinterpret_cast<std::string*>(m) = text;
        return true;
    case SettingKind::Vector3: return parseVector3(text, *reinterpret_cast<Vec3*>(m));
    case SettingKind::RealList: return parseRealList(text, *reinterpret_cast<std::vector<double>*>(m));
    case SettingKind::MaterialModel:
        if (text != "fibre" && text != "diffuse-sheen") return false;
        *reinterpret_cast<MaterialSpec::Model*>(m) = text == "fibre" ? MaterialSpec::Model::Fibre : MaterialSpec::Model::DiffuseSheen;
        return true;
    }
    return false;
}

const SettingField* findField(const std::string& key) {
    for (const SettingField& f : settingFields())
        if (key == f.key) return &f;
    return nullptr;
}

} // namespace

const std::vector<SettingField>& settingFields() {
    static const std::vector<SettingField> table = makeTable();
    return table;
}

std::string settingNumber(double v) {
    char b[40];
    std::snprintf(b, sizeof b, "%.17g", v);
    return b;
}

std::string settingValue(const SceneSettings& s, const SettingField& f) {
    const char* m = memberBytes(s, f);
    switch (f.kind) {
    case SettingKind::Real: return settingNumber(*reinterpret_cast<const double*>(m));
    case SettingKind::Integer: return std::to_string(*reinterpret_cast<const int*>(m));
    case SettingKind::LongInteger: return std::to_string(*reinterpret_cast<const long long*>(m));
    case SettingKind::Flag: return *reinterpret_cast<const bool*>(m) ? "true" : "false";
    case SettingKind::Text:
    case SettingKind::Choice: return *reinterpret_cast<const std::string*>(m);
    case SettingKind::Vector3: {
        const Vec3& v = *reinterpret_cast<const Vec3*>(m);
        return settingNumber(v.x) + ", " + settingNumber(v.y) + ", " + settingNumber(v.z);
    }
    case SettingKind::RealList: {
        std::string o;
        for (double v : *reinterpret_cast<const std::vector<double>*>(m)) o += (o.empty() ? "" : ", ") + settingNumber(v);
        return o;
    }
    case SettingKind::MaterialModel:
        return *reinterpret_cast<const MaterialSpec::Model*>(m) == MaterialSpec::Model::Fibre ? "fibre" : "diffuse-sheen";
    }
    return "";
}

std::string writeSettings(const SceneSettings& s) {
    std::string o = std::string(kVersionLine) + "\n";
    for (const SettingField& f : settingFields()) o += f.key + " = " + settingValue(s, f) + "\n";
    return o;
}

// The fibre values a settings file written before the wool presets (no fibres.preset key) was made
// with, written out: a later change of the built-in values (SceneSettings()) must not change how an
// older file looks. The built-in values of that time, before the strays followed the yarn (2026-11-10).
static void setPrePresetFibreValues(SceneSettings& s) {
    s.fibrePreset = "";
    s.randomFraction = 0.5;
    s.fibreMaxSegments = 3;
    s.directed = SceneSettings::Population{};  // lengths 0.4 / 0.8 / 1.6 mm, spread 0.35; thickness 0.07 mm
    s.directed.curvature = 1.0;
    s.directed.maxAngleDeg = 24;
    s.random = s.directed;
    s.random.curvature = 0.3;
    s.random.maxAngleDeg = 60;
    s.randomSpreadDeg = 180;
}

// The former workflow selectors: knit.enabled (true: knitted) and fibres.profile (wool | fur). They are no
// rows of the table; readSettings reads them before any key is applied.
bool isLegacySelector(const std::string& key) { return key == "knit.enabled" || key == "fibres.profile"; }

// The workflow a file's selector lines ask for (scene_settings.h readSettings has the rules); false with
// the reason: a contradiction, a repeated selector, an unknown value, or the retired Hair mode.
bool legacyWorkflow(const std::vector<KeyValueLine>& lines, std::string& workflow, std::string* err) {
    const KeyValueLine *knitLine = nullptr, *profileLine = nullptr, *workflowLine = nullptr;
    for (const KeyValueLine& l : lines) {
        const KeyValueLine** slot = l.key == "knit.enabled" ? &knitLine : l.key == "fibres.profile" ? &profileLine : l.key == "workflow" ? &workflowLine : nullptr;
        if (!slot) continue;
        if (*slot) {
            if (err) *err = "line " + std::to_string(l.line) + ": " + l.key + " is given twice";
            return false;
        }
        *slot = &l;
    }
    if (knitLine && knitLine->value != "true" && knitLine->value != "false") {
        if (err) *err = "line " + std::to_string(knitLine->line) + ": bad value for knit.enabled: " + knitLine->value;
        return false;
    }
    if (profileLine && profileLine->value != "wool" && profileLine->value != "fur") {
        if (err) *err = "line " + std::to_string(profileLine->line) + ": bad value for fibres.profile: " + profileLine->value;
        return false;
    }
    if (workflowLine && workflowLine->value != "knit" && workflowLine->value != "fur") {
        if (err) *err = "line " + std::to_string(workflowLine->line) + ": bad value for workflow: " + workflowLine->value + " (knit or fur)";
        return false;
    }
    // Each selector as the workflow it meant: knit.enabled = true and fibres.profile = wool knit, the rest is fur.
    const std::string byKnit = !knitLine ? "" : knitLine->value == "true" ? "knit" : "fur";
    const std::string byProfile = !profileLine ? "" : profileLine->value == "wool" ? "knit" : "fur";
    if (workflowLine) {
        for (const KeyValueLine* l : {knitLine, profileLine}) {
            const std::string meant = l == knitLine ? byKnit : byProfile;
            if (l && meant != workflowLine->value) {
                if (err) *err = "workflow = " + workflowLine->value + " contradicts " + l->key + " = " + l->value + " (line " + std::to_string(l->line) + ")";
                return false;
            }
        }
        workflow = workflowLine->value;
        return true;
    }
    if (!byKnit.empty() && !byProfile.empty() && byKnit != byProfile) {
        if (err) *err = "knit.enabled = " + knitLine->value + " contradicts fibres.profile = " + profileLine->value + ": a knitted object grows wool, the Hair mode skipped the knit";
        return false;
    }
    if (byKnit == "fur" && !profileLine) {
        if (err) *err = "knit.enabled = false contradicts fibres.profile = wool (the value it held when absent): wool grows on a knit";
        return false;
    }
    if (byProfile == "fur") {
        if (err)
            *err = "this file was saved by the Hair mode (fibres.profile = fur), which the Fur workflow replaced; its hair settings cannot be carried over: "
                   "open the model and choose Fur";
        return false;
    }
    workflow = "knit";
    return true;
}

bool readSettings(const std::string& textIn, SceneSettings& s, std::string* err, std::vector<std::string>* retiredKeys) {
    std::vector<KeyValueLine> lines;
    if (!parseKeyValueText(textIn, kVersionLine, lines, err)) return false;
    std::string workflow;
    if (!legacyWorkflow(lines, workflow, err)) return false;
    if (workflow == "fur" && !furEnabled()) {
        if (err) *err = kFurSwitchedOff;
        return false;
    }
    SceneSettings out = s;
    out.workflow = workflow;
    // Written before the wool presets (no fibres.preset key): the values a preset owns that the file
    // leaves out take the values of that time (setPrePresetFibreValues), not the caller's (the viewer
    // starts from a preset) nor today's built-in ones, so loading an older file never changes its fibres.
    if (textIn.find("\nfibres.preset") == std::string::npos) setPrePresetFibreValues(out);
    // Written before fibres.random.spread_deg existed (settings files are written whole): its strays
    // left the yarn in any direction.
    if (textIn.find("\nfibres.random.spread_deg") == std::string::npos) out.randomSpreadDeg = 180;
    for (const KeyValueLine& l : lines) {
        if (isLegacySelector(l.key)) continue;  // read by legacyWorkflow
        if (isRetiredSetting(l.key)) {
            if (retiredKeys && isReportedRetiredValue(l.key, l.value)) retiredKeys->push_back(l.key);
            continue;
        }
        if (!setSetting(out, l.key, l.value, err)) {
            if (err) *err = "line " + std::to_string(l.line) + ": " + *err;
            return false;
        }
    }
    // Written before the wool look: every file states the yarns' tint, and 0.5 was its default. No
    // control sets it, so 0.5 here is that default, not a choice; --set applies after the file.
    constexpr double kFormerYarnSheenTint = 0.5;
    if (out.yarnWhite.sheenTint == kFormerYarnSheenTint) out.yarnWhite.sheenTint = kYarnSheenTint;
    if (out.yarnBrown.sheenTint == kFormerYarnSheenTint) out.yarnBrown.sheenTint = kYarnSheenTint;
    s = out;
    return true;
}

bool knitsGuidedPiece(const SceneSettings& s) { return s.scene == "import" && knitWorkflow(s) && s.knitRoute == "guided"; }

bool setSetting(SceneSettings& s, const std::string& key, const std::string& value, std::string* err) {
    if (isLegacySelector(key)) {
        if (err) *err = key + " is replaced by workflow = knit | fur";
        return false;
    }
    if (key == "workflow" && value == "fur" && !furEnabled()) {
        if (err) *err = kFurSwitchedOff;
        return false;
    }
    const SettingField* f = findField(key);
    if (!f && !isRetiredSetting(key)) {
        if (err) *err = "unknown key " + key;
        return false;
    }
    if (f && !setFieldText(s, *f, value)) {
        if (err) *err = "bad value for " + key + ": " + value;
        return false;
    }
    applyStitchSize(s);  // yarn mode: w and h follow D, c and the tightness
    return true;
}

void copySetting(SceneSettings& to, const SceneSettings& from, const std::string& key) {
    const SettingField* f = findField(key);
    if (f) setFieldText(to, *f, settingValue(from, *f));
}

bool sameSetting(const SceneSettings& a, const SceneSettings& b, const std::string& key) {
    const SettingField* f = findField(key);
    if (!f) return isRetiredSetting(key);
    return sameValueText(settingValue(a, *f), settingValue(b, *f));
}

std::vector<std::string> diffSettings(const SceneSettings& a, const SceneSettings& b) {
    std::vector<std::string> d;
    for (const SettingField& f : settingFields())
        if (settingValue(a, f) != settingValue(b, f)) d.push_back(f.key);
    return d;
}

} // namespace knit
