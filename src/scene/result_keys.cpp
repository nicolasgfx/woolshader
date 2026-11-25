#include "scene/result_keys.h"

#include "gpu/knitting.h"
#include "knitting/checks/checks.h"
#include "knitting/stitch/hand_tension.h"
#include "knitting/tube/yarn_surface.h"
#include "model/hash.h"
#include "scene/scene.h"
#include "stitchmesh/yarn_settle.h"
#include "settings/settings_mapping.h"

namespace knit {

namespace {

// The Fur workflow's own keys.
bool isPileKey(const std::string& k) { return k.rfind("pile.", 0) == 0 || k.rfind("material.pile", 0) == 0; }

// The rebuilt surface's own keys. Never inputs of a keep-file import: its saved results must keep their keys.
bool isWrapKey(const std::string& k) {
    return k == "import.surface" || k == "import.wrap_stitch_mm";
}

// Retired keys whose lines stay in the text of `stage` and every later stage, after the key named, with their
// last defaults: saved results are keyed by that text.
struct RetiredInput {
    const char* after;
    Stage stage;
    const char* lines;
};
constexpr RetiredInput kRetiredInputs[] = {
    {"strand.clearance", Stage::Knit, "strand.contact=clearance\nstrand.template=fixed\n"},
    {"knit.coverage", Stage::Knit, "knit.layout=uv-tiling\n"},
    {"strand.angle_deg", Stage::Knit, "knit.lattice_anchor=bbox\n"},
    {"fibres.enabled", Stage::Fibres, "fibres.wool=true\n"},
    // 2026-11-19: the workflow row replaced knit.enabled at the same place; a knitted object always held true.
    {"workflow", Stage::Knit, "knit.enabled=true\n"},
};

// Whether key k can change a stage's result under these settings. A label names values that are
// inputs themselves; a feature's parameters count only while the feature is in use.
bool isEffectiveInput(const SceneSettings& s, const std::string& k) {
    if (k == "strand.thickness_preset" || k == "fibres.preset") return false;      // names the values
    if (k == "strand.stitch_size" || k == "strand.tightness") return false;        // name the derived w and h
    if (k == "strand.compress_to_fit") return false;                               // an input as the contact fraction it gives
    if (k == "strand.max_flattening") return renderFlattened(s) && s.maxFlattening != kDefaultMaxCompression;
    if (k == "import.plain_nodes") return false;                                   // import.identity depends on it
    if (k == "workflow") return false;                                             // selects the keys below; knit.enabled=true above
    if (k == "knit.guided_markings") return false;                                 // where the file is: its content is the digest's
    if (k == "knit.route" || k.rfind("knit.guided_", 0) == 0) return knitsGuidedPiece(s);  // the guided route's rows: no line otherwise
    if (isWrapKey(k)) return s.scene == "import" && s.importSurface == "wrapped";
    if (k.rfind("import.", 0) == 0) return s.scene == "import";
    if (isPileKey(k)) return furWorkflow(s);
    if (k.rfind("knit.hand.", 0) == 0) return s.hand.inUse();
    if (k.rfind("fibres.islands.", 0) == 0) return s.islands.inUse();
    if (k.rfind("yarn_body.", 0) == 0) return s.yarnBody.enabled;
    if (k == "fibres.directed.hug_height") return s.directedHugHeight > 0;
    if (k == "fibres.directed.follow_plies") return s.directedFollowPlies;
    if (k == "fibres.directed.drape") return s.directedDrape;
    if (k == "fibres.random.spread_deg") return s.randomSpreadDeg < 180;
    if (k == "knit.overlap_stitches" || k == "knit.gap_patches") return overlapCoverage(s);
    return true;
}

} // namespace

std::string stageInputs(const SceneSettings& s, Stage stage) {
    std::string inputs;
    const auto upTo = [stage](Stage st) { return st != Stage::Display && static_cast<int>(st) <= static_cast<int>(stage); };
    for (const SettingField& f : settingFields()) {
        if (upTo(f.stage) && isEffectiveInput(s, f.key)) inputs += f.key + "=" + settingValue(s, f) + "\n";
        for (const RetiredInput& r : kRetiredInputs)
            if (f.key == r.after && upTo(r.stage)) inputs += r.lines;
    }
    if (static_cast<int>(stage) < static_cast<int>(Stage::Knit)) return inputs;
    // Derived inputs of the knitting: the crossing attribution of the local signature (checks.h)
    // decides which controls settle frees, and the contact fraction which yarn is knitted.
    inputs += "derived.attribution_window=" + settingNumber(kAttributionWindow) + "\n";
    inputs += "derived.contact_fraction=" + settingNumber(yarnFit(s).k) + "\n";
    // a guided knit runs on its settle's solver options, every other route on the environment's
    const SettleSettings guided;
    inputs += "derived.ncp_solver=" + (knitsGuidedPiece(s) ? gpu::ncpSolverIdentity(guided.contactStorage, guided.curvatureRows) : gpu::ncpSolverIdentity()) + "\n";
    if (s.hand.inUse()) inputs += std::string("derived.hand_rest=") + kHandRestRule + "\n";  // settle's rest shape
    if (knitsGuidedPiece(s)) inputs += "derived.guided_version=" + std::to_string(kGuidedVersion) + "\n";  // its templates, yarn, build and settle settings
    return inputs;
}

uint64_t stageKey(const SceneSettings& s, Stage stage, uint64_t assetHash) {
    // Bump when a generator change alters the output of identical settings (every derived result
    // then misses once; nothing is deleted).
    constexpr uint64_t kGeneratorVersion = 2;
    uint64_t h = fnv1a(&kGeneratorVersion, sizeof kGeneratorVersion);
    h = fnv1a(&assetHash, sizeof assetHash, h);
    const std::string inputs = stageInputs(s, stage);
    return fnv1a(inputs.data(), inputs.size(), h);
}

uint64_t furKey(const SceneSettings& s, uint64_t assetHash) {
    // 2: the Fur workflow's loop pile (fur_field.h, fur_loops.h); 1 was the Hair mode's groom (never reused).
    constexpr uint64_t kFurVersion = 2;  // bump when the fur generator's output changes
    uint64_t h = fnv1a(&kFurVersion, sizeof kFurVersion);
    h = fnv1a(&assetHash, sizeof assetHash, h);
    for (const SettingField& f : settingFields()) {
        const std::string& k = f.key;
        const bool load = f.stage == Stage::Load && (k.rfind("import.", 0) != 0 || s.scene == "import");
        if (!load && !isPileKey(f.key) && k != "fibres.segment_budget") continue;
        if (k == "import.plain_nodes" || (isWrapKey(k) && !isEffectiveInput(s, k))) continue;
        if (k.rfind("material.", 0) == 0) continue;  // appearance: the export only
        const std::string kv = k + "=" + settingValue(s, f) + "\n";
        h = fnv1a(kv.data(), kv.size(), h);
    }
    return h;
}

uint64_t geometryKey(const SceneSettings& s) { return stageKey(s, Stage::Surface, 0); }

uint64_t fibreResultKey(const SceneSettings& s) {
    if (furWorkflow(s)) return furKey(s, 0);
    SceneSettings t = s;
    t.fibres = true;  // fibres.enabled says whether the export includes them, not what they are
    uint64_t h = stageKey(t, Stage::Fibres, 0);
    return fnv1a(&t.segmentBudget, sizeof t.segmentBudget, h);
}

} // namespace knit
