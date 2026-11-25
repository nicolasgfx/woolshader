// S6 tests (scene viewer, settings, GLB scenes, UV-chart tiling, patterns). Run one
// case: knit_tests_s6 <name>; all: knit_tests_s6; list: knit_tests_s6 --list. A skipped
// case (missing asset, changed asset hash) prints [SKIP] and does not count as passed.
#include <algorithm>
#include <random>
#include <cctype>
#include <atomic>
#include <chrono>
#include <cstring>
#include <mutex>
#include <thread>
#include <functional>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <string>

#include "harness.h"
#include "ply_reference.h"
#include "scene/appearance/filling.h"
#include "scene/fibres/reference_samples.h"
#include "scene/stage_params.h"
#include "scene_fixtures.h"
#include "settings/settings_profiles.h"
#include "strand_fixtures.h"
#include "test_compute.h"
#include "cpu_numerics.h"
#include "knitting/knit/mesh_knit.h"
#include "knitting/knit/coverage.h"
#include "scene/knit/knit_cache.h"
#include "knitting/relax/knit_problem.h"
#include "fibres/fibre_checks.h"
#include "model/choices.h"
#include "model/hash.h"
#include "knitting/knit/fixtures.h"
#include "fibres/fibre.h"
#include "fibres/fibre_budget.h"
#include "fibres/fibre_generator.h"
#include "fibres/fibre_shape.h"
#include "knitting/tube/yarn_surface.h"
#include "knitting/checks/checks.h"
#include "knitting/stitch/expected_signature.h"
#include "knitting/tube/param_text.h"
#include "knitting/stitch/stockinette.h"
#include "model/mesh_check.h"
#include "gpu/device.h"
#include "gpu/knitting.h"
#include "io/glb_load.h"
#include "io/image_io.h"
#include "settings/file_presets.h"
#include "io/publish.h"
#include "jobs/scene_builder.h"
#include "jobs/scene_cache.h"
#include "scene/source/sewn_pieces.h"
#include "scene/knit/spiral_cache.h"
#include "settings/scene_registry.h"
#include "settings/scene_settings.h"
#include "scene/knit/gap_patches.h"
#include "settings/thickness_presets.h"
#include "scene/scene.h"
#include "settings/settings_mapping.h"
#include "scene/source/scene_parts.h"
#include "scene/fibres/part_fibre_job.h"
#include "scene/source/scene_solids.h"
#include "scene/knit/knit_stage.h"
#include "scene/knit/knit_report.h"
#include "scene/appearance/rendered_part.h"
#include "scene/fibres/fibre_stage.h"
#include "export/export.h"
#include "scene/appearance/ply_maps.h"
#include "scene/appearance/stitch_colours.h"
#include "model/cancel.h"
#include "model/mesh_query.h"
#include "model/parallel.h"
#include "scene/knit/knit_cache.h"
#include "io/files.h"
#include "scene/result_keys.h"
#include "knitting/knit/pieces.h"

using namespace knit;

namespace {

// ------------------------------------------------------------- settings ----

// Every key round-trips, both at the defaults and with every key changed.
// Another word of a choice key: a new key is covered without naming it here.
std::string otherChoice(const SettingField& k, const std::string& current) {
    for (const std::string& word : k.choices)
        if (word != current) return word;
    return current;
}

// The Fur switch (scene_settings.h furEnabled): while it is off, neither a file nor --set can ask for workflow = fur,
// and a knit's file reads as before; switched on again, the same fur file reads.
TEST(fur_switch_refuses_the_fur_workflow) {
    SceneSettings fur;
    fur.workflow = "fur";
    const std::string furText = writeSettings(fur);
    const std::string knitText = writeSettings(SceneSettings{});
    SceneSettings back;
    std::string err;
    setFurEnabled(false);
    const bool furRead = readSettings(furText, back, &err);
    const std::string readErr = err;
    const bool knitRead = readSettings(knitText, back, &err) && knitWorkflow(back);
    SceneSettings set;
    const bool furSet = setSetting(set, "workflow", "fur", &err);
    const std::string setErr = err;
    const bool knitSet = setSetting(set, "workflow", "knit", &err);
    setFurEnabled(true);
    CHECK(!furRead && readErr.find("switched off") != std::string::npos);
    CHECK(knitRead && knitSet);
    CHECK(!furSet && setErr == readErr);
    CHECK(readSettings(furText, back, &err) && furWorkflow(back));
}

TEST(scene_settings_roundtrip) {
    const SceneSettings d;
    SceneSettings back;
    std::string err;
    CHECK(readSettings(writeSettings(d), back, &err));
    CHECK(back == d);
    CHECK(writeSettings(back) == writeSettings(d));
    // Change every key: numbers by a non-representable-in-short-decimal amount, flags
    // flipped, texts replaced; the re-read must reproduce the exact values.
    SceneSettings m = d;
    const std::vector<SettingField>& keys = settingFields();
    std::set<std::string> seen;
    for (const SettingField& k : keys) {
        CHECK(seen.insert(k.key).second);  // unique keys
        const std::string line = [&] {
            const std::string all = writeSettings(m);
            const size_t p = all.find("\n" + k.key + " = ");
            CHECK(p != std::string::npos);
            const size_t v = p + 1 + k.key.size() + 3;
            return all.substr(v, all.find('\n', v) - v);
        }();
        std::string nv;
        if (!k.choices.empty()) nv = otherChoice(k, line);
        else if (line == "true" || line == "false") nv = line == "true" ? "false" : "true";
        else if (k.key.size() > 6 && k.key.substr(k.key.size() - 6) == ".model") nv = line == "fibre" ? "diffuse-sheen" : "fibre";
        else if (k.key.size() > 8 && k.key.substr(k.key.size() - 8) == ".backend") nv = line == "cpu" ? "gpu" : "cpu";
        else if (k.key == "pile.seed" || k.key == "knit.hand.seed" || k.key == "knit.guided_piece" || k.key == "knit.guided_component")
            nv = "7";
        else if (k.key == "strand.ply_handedness") nv = "-1";
        else if (k.key == "fibres.segment_budget" || k.key == "strand.plies" || k.key == "knit.diameter_attempts" ||
                 k.key == "strand.tube_sides" || k.key == "strand.ply_tube_sides" || k.key == "fibres.max_segments" ||
                 k.key == "knit.overlap_stitches" || k.key == "yarn_body.plies" || k.key == "yarn_body.fibres" ||
                 k.key == "yarn_body.interior_fibres" || k.key == "yarn_body.seed")
            nv = "7";
        else if (line.find(',') != std::string::npos || k.key == "strand.diameter_fallback") nv = "0.1, 0.30000000000000004, 1e-7";
        else if (char* end = nullptr; std::strtod(line.c_str(), &end), end && *end == 0 && !line.empty()) nv = "0.30000000000000004";
        else nv = "changed-" + k.key;
        const bool set = setSetting(m, k.key, nv, &err);
        if (!set) std::printf("    refused: %s = %s (%s)\n", k.key.c_str(), nv.c_str(), err.c_str());
        CHECK(set);
    }
    CHECK(diffSettings(d, m).size() == keys.size());  // every key changed
    SceneSettings back2;
    CHECK(readSettings(writeSettings(m), back2, &err));
    CHECK(back2 == m);
    CHECK(diffSettings(back2, m).empty());
    // The wool look: the yarns' sheen takes their colour, the shared default and the fuzz's stay. A
    // file written before it states the former default for the yarns and reads as the look.
    CHECK(d.yarnWhite.sheenTint == 1.0 && d.yarnBrown.sheenTint == 1.0);
    CHECK(d.fibreDirected.sheenTint == 0.5 && d.fibreRandom.sheenTint == 0.5 && MaterialSpec{}.sheenTint == 0.5);
    SceneSettings saved;
    CHECK(readSettings("knit-settings v1\nmaterial.yarn_white.sheen_tint = 0.5\nmaterial.yarn_brown.sheen_tint = 0.5\n"
                       "material.fibre_directed.sheen_tint = 0.25\n", saved, &err));
    CHECK(saved.yarnWhite.sheenTint == 1.0 && saved.yarnBrown.sheenTint == 1.0 && saved.fibreDirected.sheenTint == 0.25);
    CHECK(readSettings("knit-settings v1\nmaterial.yarn_white.sheen_tint = 0.3\n", saved, &err) && saved.yarnWhite.sheenTint == 0.3);
    CHECK(setSetting(saved, "material.yarn_white.sheen_tint", "0.5", &err) && saved.yarnWhite.sheenTint == 0.5);  // --set applies as stated
}

// The stage inputs say why a key changed: each stage lists its upstream inputs, a feature's
// parameters appear only while it is in use, labels and retired keys never, and the keys hash
// exactly this text (equal inputs and asset: equal keys).
TEST(scene_settings_stage_inputs) {
    auto lists = [](const std::string& inputs, const std::string& key) { return ("\n" + inputs).find("\n" + key + "=") != std::string::npos; };
    SceneSettings s;
    const std::string knit = stageInputs(s, Stage::Knit);
    CHECK(lists(knit, "scene") && lists(knit, "strand.diameter"));
    CHECK(lists(knit, "knit.layout") && lists(knit, "strand.contact") && lists(knit, "strand.template"));  // retired: saved knits keep their keys
    CHECK(lists(knit, "derived.contact_fraction") && !lists(stageInputs(s, Stage::Load), "derived.contact_fraction"));
    CHECK(!lists(knit, "strand.tube_sides") && lists(stageInputs(s, Stage::Fibres), "strand.tube_sides"));
    CHECK(lists(stageInputs(s, Stage::Export), "strand.diameter"));
    for (const char* never : {"strand.thickness_preset", "fibres.preset", "strand.stitch_size", "strand.teddy_stitches", "knit.backend", "knit.time_cap_s"})
        CHECK(!lists(stageInputs(s, Stage::Export), never));
    CHECK(!lists(knit, "strand.compression") && !lists(knit, "knit.overlap_stitches"));
    CHECK(lists(knit, "derived.ncp_solver") && !lists(knit, "knit.solver"));
    SceneSettings u = s;
    u.coverage = "overlap";
    const std::string used = stageInputs(u, Stage::Knit);
    CHECK(lists(used, "knit.overlap_stitches"));
    // Unused parameters leave the text, and with it the key, alone.
    SceneSettings same = s;
    same.fibrePreset = "custom";
    CHECK(stageInputs(same, Stage::Export) == stageInputs(s, Stage::Export));
    CHECK(stageKey(same, Stage::Export, 7) == stageKey(s, Stage::Export, 7));
}

// Retired keys (scene_settings.cpp retiredField): a saved file that still holds them loads, they
// are reported and change nothing, they are never written or listed; a key that never existed is
// still refused (a typo is not a retired key).
TEST(scene_settings_retired_keys) {
    const SceneSettings d;
    std::string text = writeSettings(d);
    CHECK(text.find("strand.teddy_stitches") == std::string::npos && text.find("display.ui_scale") == std::string::npos);
    text += "strand.teddy_stitches = 80\ndisplay.ui_scale = 1.5\nknit.solver = legacy\n";
    SceneSettings back;
    std::string err;
    std::vector<std::string> retired;
    CHECK(readSettings(text, back, &err, &retired));
    CHECK(back == d);
    CHECK(retired == std::vector<std::string>({"strand.teddy_stitches", "display.ui_scale", "knit.solver"}));
    CHECK(setSetting(back, "strand.teddy_stitches", "90", &err) && back == d);
    CHECK(!readSettings(text + "strand.teddy_stiches = 80\n", back, &err));
    std::printf("  a misspelt key: %s\n", err.c_str());
    for (const SettingField& k : settingFields()) CHECK(k.key != "strand.teddy_stitches" && k.key != "display.ui_scale");
}

// The keys of removed knit features: every saved file holds them at their last defaults, which load
// without a note. Any other value is reported and ignored. Either way the Knit key does not move.
TEST(scene_settings_retired_knit_keys) {
    struct RetiredKey { std::string key, lastDefault, otherValue; };
    const std::vector<RetiredKey> keys = {{"strand.contact", "clearance", "compress"},
                                          {"strand.compression", "0.10000000000000001", "0.2"},
                                          {"strand.core_clearance", "0", "0.05"},
                                          {"strand.template", "fixed", "relaxed"},
                                          {"strand.cell_variant", "L", "E"},
                                          {"strand.cell_tightness", "1", "0.9"},
                                          {"knit.layout", "uv-tiling", "stitch-mesh"}};
    const SceneSettings d;
    const std::string saved = writeSettings(d);
    std::string atDefault = saved, changed = saved;
    std::vector<std::string> expected;
    for (const RetiredKey& k : keys) {
        CHECK(saved.find(k.key) == std::string::npos);
        atDefault += k.key + " = " + k.lastDefault + "\n";
        changed += k.key + " = " + k.otherValue + "\n";
        expected.push_back(k.key);
    }
    SceneSettings a, b;
    std::string err;
    std::vector<std::string> notesAtDefault, notesChanged;
    CHECK(readSettings(atDefault, a, &err, &notesAtDefault) && a == d && notesAtDefault.empty());
    CHECK(readSettings(changed, b, &err, &notesChanged) && b == d && notesChanged == expected);
    CHECK(stageKey(a, Stage::Knit, 1) == stageKey(d, Stage::Knit, 1) && stageKey(b, Stage::Knit, 1) == stageKey(d, Stage::Knit, 1));
    for (const SettingField& s : settingFields())
        for (const RetiredKey& k : keys) CHECK(s.key != k.key);
}

// The keys scalp hair left, as every saved file holds them (the text of knit_check 72bd789): they load
// without a note. Another value is reported and ignored, and fibres.wool keeps its line in the Fibres
// text, so saved fibres keep their keys.
TEST(scene_settings_retired_hair_keys) {
    const std::string saved =
        "fibres.wool = true\n"
        "hair.enabled = false\n"
        "hair.density = 1.5\n"
        "hair.guide_density = 0.050000000000000003\n"
        "hair.root_diameter = 0.070000000000000007\n"
        "hair.tip_ratio = 0.80000000000000004\n"
        "hair.segment_length = 1.5\n"
        "hair.crown = 25\n"
        "hair.volume = 30\n"
        "hair.tip_volume = 14\n"
        "hair.cut_below_top = 88\n"
        "hair.cut_jitter = 2\n"
        "hair.lift_near_deg = 35\n"
        "hair.lift_far_deg = 10\n"
        "hair.part_wiggle = 1\n"
        "hair.clump = 0.55000000000000004\n"
        "hair.wave_amplitude = 1.5\n"
        "hair.frizz = 0.25\n"
        "hair.stray_fraction = 0.01\n"
        "hair.tones = 3\n"
        "hair.tone_spread = 0.040000000000000001\n"
        "hair.settle_steps = 0\n"
        "hair.gravity = 1\n"
        "hair.segment_budget = 8000000\n"
        "hair.seed = 0\n"
        "material.hair.model = fibre\n"
        "material.hair.color = 0.20000000000000001, 0.5, 0.20000000000000001\n"
        "material.hair.sheen = 1\n"
        "material.hair.sheen_tint = 0.5\n"
        "material.hair.sheen_roughness = 0.5\n"
        "material.hair.clearcoat = false\n"
        "material.hair.clearcoat_weight = 0.29999999999999999\n"
        "material.hair.clearcoat_ior = 1.5\n"
        "material.hair.clearcoat_roughness = 0.29999999999999999\n"
        "material.hair.fibre_longitudinal_roughness = 0.29999999999999999\n"
        "material.hair.fibre_azimuthal_roughness = 0.29999999999999999\n"
        "material.hair.fibre_ior = 1.55\n"
        "material.hair.fibre_cuticle_tilt_deg = 2\n";
    const SceneSettings d;
    SceneSettings a = d, b = d;
    std::string err;
    std::vector<std::string> notes;
    CHECK(readSettings(writeSettings(d) + saved, a, &err, &notes) && a == d);
    for (const std::string& k : notes) std::printf("    noted: %s\n", k.c_str());
    CHECK(notes.empty());
    CHECK(readSettings(writeSettings(d) + "hair.enabled = true\nhair.density = 3\nfibres.wool = false\nmaterial.hair.color = 0, 0, 0\n", b, &err, &notes) && b == d);
    CHECK(notes == std::vector<std::string>({"hair.enabled", "hair.density", "fibres.wool", "material.hair.color"}));
    const std::string fibres = stageInputs(d, Stage::Fibres);
    CHECK(fibres.find("fibres.wool=true\n") != std::string::npos && stageInputs(d, Stage::Knit).find("fibres.wool") == std::string::npos);
}

// Knitted settings saved before the workflow key (tests/data/legacy_settings, written by knit_check of main
// 4597176, the last build with knit.enabled and fibres.profile) keep every stage key, so their saved knits
// and fibres are reused. The pinned numbers are what main computes for the same files: 4597176's own values until
// main 2a4928e changed the NCP solver identity to v5 (compact contact rows), which moves every Knit and later key;
// with the v4 identity text swapped back, this build gives 4597176's values again (checked 2026-11-19 at the merge).
static std::string legacySettingsText(const std::string& name) {
    std::ifstream in(std::string(KNIT_TEST_DATA) + "/legacy_settings/" + name);
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}
static std::string withoutLine(std::string text, const std::string& key) {
    const size_t at = text.find("\n" + key + " = ");
    if (at != std::string::npos) text.erase(at + 1, text.find('\n', at + 1) - at);
    return text;
}
static std::string withLine(std::string text, const std::string& key, const std::string& value) {
    return withoutLine(text, key) + key + " = " + value + "\n";
}
TEST(legacy_knit_settings_keep_their_keys) {
    const std::string teddy = legacySettingsText("teddy_knit.settings");
    CHECK(!teddy.empty());
    std::string import = withLine(withLine(withLine(teddy, "scene", "import"), "import.identity", "0123456789abcdef"), "import.surface", "file");
    struct Case {
        const char* name;
        std::string text;
        std::vector<uint64_t> keys;  // Load, Knit, Surface, Fibres, Export stage keys (asset 0x1234), fibreResultKey
    };
    const std::vector<Case> cases = {
        {"teddy", teddy,
         {0x1a5bd1293b3c9eb6ull, 0x15dd9ed41515aa88ull, 0xde94d0e7c859e2e6ull, 0xbf2b57dc4b294eb4ull, 0xa68e40dd87ac6a7eull, 0x8e0ab69098acdd58ull}},
        {"teddy before knit.enabled", withoutLine(teddy, "knit.enabled"),
         {0x1a5bd1293b3c9eb6ull, 0x15dd9ed41515aa88ull, 0xde94d0e7c859e2e6ull, 0xbf2b57dc4b294eb4ull, 0xa68e40dd87ac6a7eull, 0x8e0ab69098acdd58ull}},
        {"swatch with fibres", legacySettingsText("swatch_fibres.settings"),
         {0x9c2d09d713548e54ull, 0xa109b71c461bcb2bull, 0x7be1f26ba3d441d4ull, 0x07879a02dd02c5bcull, 0x762923d6f7a651c3ull, 0x8d191da81dd67511ull}},
        {"keep-file import", import,
         {0xda0a8b40247eb08eull, 0x9f27d57951575db0ull, 0x7006cf74c997529eull, 0xf7768dbd19c192acull, 0x1d75137da6537ff6ull, 0xc0f9775c956fa284ull}},
    };
    bool all = true;
    for (const Case& c : cases) {
        SceneSettings s;
        std::string err;
        const bool read = readSettings(c.text, s, &err);
        all &= read;
        if (!read) {
            std::printf("    %s: %s\n", c.name, err.c_str());
            continue;
        }
        std::vector<uint64_t> keys;
        for (Stage st : {Stage::Load, Stage::Knit, Stage::Surface, Stage::Fibres, Stage::Export}) keys.push_back(stageKey(s, st, 0x1234));
        keys.push_back(fibreResultKey(s));
        std::printf("    %s:", c.name);
        for (uint64_t k : keys) std::printf(" 0x%016llxull,", static_cast<unsigned long long>(k));
        std::printf("\n");
        all &= keys == c.keys;
    }
    CHECK(all);
}

TEST(scene_settings_refusals) {
    SceneSettings s;
    std::string err;
    CHECK(!readSettings("scene = man\n", s, &err));                       // no version line
    CHECK(!readSettings("knit-settings v2\nscene = man\n", s, &err));     // other version
    CHECK(!readSettings("knit-settings v1\nno.such.key = 1\n", s, &err));  // unknown key
    CHECK(err.find("unknown key") != std::string::npos);
    CHECK(!readSettings("knit-settings v1\nstrand.diameter = 1.2x\n", s, &err));  // bad number
    CHECK(!readSettings("knit-settings v1\nstrand.settle = yes\n", s, &err));     // bad bool
    CHECK(readSettings("knit-settings v1\nfibres.backend = cuda\n", s, &err));    // retired (GPU only): any value is ignored
    CHECK(!readSettings("knit-settings v1\nstrand.diameter 1.2\n", s, &err));     // no '='
    // A refused file leaves the settings untouched.
    SceneSettings before = s;
    CHECK(!readSettings("knit-settings v1\nstrand.diameter = 2\nbad = 1\n", s, &err));
    CHECK(s == before);
    CHECK(readSettings("# comment\nknit-settings v1\n\nstrand.diameter = 1.25\n", s, &err));
    CHECK(s.diameter == 1.25);
}

// The knit stage key hashes gpu::ncpSolverIdentity(), the revision (v4: the device-controlled
// iteration; v5: the compact contact rows, whose sums round differently), the early-stop defaults, the step budget and settle's one round with its own stopping
// tolerances, so a knit of an earlier solver revision is not reused.
TEST(scene_settings_ncp_key_names_the_solver) {
    const std::string id = gpu::ncpSolverIdentity();
    std::printf("    %s\n", id.c_str());
    CHECK(id.rfind("ncp-gpu v5", 0) == 0);
    CHECK(id.find("settle: one round, progress 0.05, contact tolViol 0.005") != std::string::npos);
    CHECK(stageInputs(SceneSettings{}, Stage::Knit).find(id) != std::string::npos);
}

// The text behind the knit cache keys. A changed line misses every stored knit, so it changes only on
// purpose, and retiring a setting must keep this text.
TEST(knit_params_text_pinned) {
    const uint64_t pattern = fnv1a64(patternParamsText(PatternParams{}));
    const uint64_t surface = fnv1a64(surfaceParamsText(SurfaceParams{}));
    std::printf("  pattern %016llx surface %016llx\n", static_cast<unsigned long long>(pattern), static_cast<unsigned long long>(surface));
    CHECK(pattern == 0x2e164b2901483422ull);
    CHECK(surface == 0xcd7c45132e4f4026ull);
}

// The ply-map set key of the default spec (pinned at 05960b1): the five cached sets under out/cache/plymaps and
// every package's maps are found by it.
TEST(ply_map_key_pinned) {
    const uint64_t key = fnv1a64(plyMapSetKey(PlyMapSpec{}));
    std::printf("  %s\n  %016llx\n", plyMapSetKey(PlyMapSpec{}).c_str(), static_cast<unsigned long long>(key));
    CHECK(key == 0xd55b4a1404eb8b97ull);
}

TEST(knit_stage_key_pinned) {
    SceneSettings plain;
    SceneSettings overlap;
    overlap.coverage = "overlap";
    const uint64_t knit = fnv1a64(stageInputs(plain, Stage::Knit));
    const uint64_t all = fnv1a64(stageInputs(plain, Stage::Export));
    const uint64_t seams = fnv1a64(stageInputs(overlap, Stage::Knit));
    std::printf("  knit %016llx all %016llx seams %016llx\n", static_cast<unsigned long long>(knit), static_cast<unsigned long long>(all),
                static_cast<unsigned long long>(seams));
    // Re-pinned 2026-11-19: the solver identity v5 (compact contact rows); with the v4 text the former pins hold
    // (checked). Before: 2026-11-10, the default scene is "teddy" (the man was removed); 2026-11-08, identity v4.
    CHECK(knit == 0x9a32b06327efe554ull);  // re-pinned 2026-11-25: default scene wool-swatch ("teddy" gives the former 0x347dcae6161bdb72)
    // Re-pinned 2026-11-19: identity v5 (as above). Before: 2026-11-11, the export text holds the fibre inputs,
    // whose built-in strays are now the Loose stray preset (spread 45, rise to 40, length spread 0.5).
    CHECK(all == 0xe39080d9c017d7bdull);  // re-pinned 2026-11-25, as above (former 0x1fcd2a90ca5f1ea7)
    CHECK(seams == 0x03b49e23db2f2f66ull);  // re-pinned 2026-11-25, as above (former 0xe6fcc2547c2ab890)
}

// A hand's values count only while it is in use: the seed of a hand that is off changes no result.
TEST(hand_off_leaves_the_keys_alone) {
    SceneSettings off, otherSeed, on;
    otherSeed.hand.seed = 5;
    on.hand.unevenness = 12;
    CHECK(stageInputs(otherSeed, Stage::Knit) == stageInputs(off, Stage::Knit));
    CHECK(stageInputs(on, Stage::Knit) != stageInputs(off, Stage::Knit));
}

// fibres.directed.hug_height and follow_plies count only while set: off, the Fibres text is the
// default one (fibre_stage_text_pinned); on, each changes it.
TEST(fibre_hug_keys_only_while_set) {
    SceneSettings off, hug, follow;
    hug.directedHugHeight = 0.15;
    follow.directedFollowPlies = true;
    CHECK(stageInputs(off, Stage::Fibres).find("hug_height") == std::string::npos);
    CHECK(stageInputs(off, Stage::Fibres).find("follow_plies") == std::string::npos);
    CHECK(stageInputs(hug, Stage::Fibres) != stageInputs(off, Stage::Fibres));
    CHECK(stageInputs(follow, Stage::Fibres) != stageInputs(off, Stage::Fibres));
}

// fibres.directed.drape is preset-owned and true only in Mohair: choosing another profile after Mohair clears it, and
// the fibre key then is that profile's own (no "drape"). With Mohair the drape is part of the key (saved fibres of the
// cut rule are not reused) and of the Fibres stage's inputs; off, it is in neither.
TEST(fibre_drape_preset_transitions) {
    const FilePreset* mohair = presetWithId(PresetGroup::Fibres, "mohair");
    CHECK(mohair != nullptr);
    if (!mohair) return;
    std::string err;
    SceneSettings m;
    applyPreset(PresetGroup::Fibres, m, *mohair);
    FibreParams mp;
    CHECK(m.directedDrape && fibreParamsFromSettings(m, mp, &err) && mp.drape);
    CHECK(fibreParamsText(mp).find(" drape") != std::string::npos);
    SceneSettings mOff = m;
    mOff.directedDrape = false;
    FibreParams mpOff;
    CHECK(fibreParamsFromSettings(mOff, mpOff, &err) && fibreParamsText(mpOff) != fibreParamsText(mp));
    CHECK(stageInputs(m, Stage::Fibres) != stageInputs(mOff, Stage::Fibres));
    CHECK(stageInputs(SceneSettings{}, Stage::Fibres).find("drape") == std::string::npos);
    int others = 0;
    for (const FilePreset& p : presets(PresetGroup::Fibres)) {
        if (p.id == "mohair") continue;
        ++others;
        SceneSettings after = m, fresh;
        applyPreset(PresetGroup::Fibres, after, p);
        applyPreset(PresetGroup::Fibres, fresh, p);
        FibreParams a, f;
        CHECK(!after.directedDrape && fibreParamsFromSettings(after, a, &err) && fibreParamsFromSettings(fresh, f, &err));
        CHECK(fibreParamsText(a) == fibreParamsText(f) && fibreParamsText(a).find("drape") == std::string::npos);
        CHECK(stageInputs(after, Stage::Fibres) == stageInputs(fresh, Stage::Fibres));
    }
    CHECK(others == 5);
}

// The default Fibres text, which saved fibres are keyed by: a key may leave the settings only if its
// line stays here.
TEST(fibre_stage_text_pinned) {
    const uint64_t fibres = fnv1a64(stageInputs(SceneSettings{}, Stage::Fibres));
    std::printf("  fibres %016llx\n", static_cast<unsigned long long>(fibres));
    // Re-pinned 2026-11-19: the solver identity v5 (the Fibres text includes the knit's inputs); with the v4 text
    // the former pin holds (checked). Before: 2026-11-11, the built-in strays are the Loose stray preset (spread 45,
    // rise to 40, length spread 0.5); 2026-11-10, the default scene is "teddy" (the man was removed).
    CHECK(fibres == 0x37ae0d3e8a1668d6ull);  // re-pinned 2026-11-25: default scene wool-swatch ("teddy" gives the former 0x9f1ec979498bb084)
}

// Stage keys: a key changes exactly the stages at and below its own. The import.*
// keys count for scene = "import" only.
TEST(scene_settings_stage_keys) {
    SceneSettings a;
    const uint64_t asset = 42;
    auto keysOf = [&](const SceneSettings& s, uint64_t h) {
        std::vector<uint64_t> k;
        for (Stage st : {Stage::Load, Stage::Knit, Stage::Surface, Stage::Fibres, Stage::Export}) k.push_back(stageKey(s, st, h));
        return k;
    };
    const std::vector<uint64_t> base = keysOf(a, asset);
    std::string err;
    for (const SettingField& k : settingFields()) {
        SceneSettings b = a;
        // A value valid for every key type that differs from the default ...
        const std::string all = writeSettings(b);
        const size_t p = all.find("\n" + k.key + " = ") + 1 + k.key.size() + 3;
        const std::string cur = all.substr(p, all.find('\n', p) - p);
        const bool numeric = cur.find_first_not_of("-0123456789.e+, ") == std::string::npos;
        std::string nv = cur == "true" ? "false" : cur == "false" ? "true" : cur == "fibre" ? "diffuse-sheen"
                         : cur == "diffuse-sheen" ? "fibre" : cur == "cpu" ? "gpu" : cur == "gpu" ? "cpu" : cur == "legacy" ? "ncp" : cur == "ncp" ? "legacy"
                         : !numeric ? "x" : cur.find(',') != std::string::npos ? "3, 3, 3" : "3";
        if (!k.choices.empty()) nv = otherChoice(k, cur);
        if (k.key == "material.yarn_colour_source") nv = "object";
        if (k.key == "material.plain_colour_source") nv = "material";
        if (k.key == "pile.coverage") nv = "mask";
        if (k.key == "material.pile_colour_source") nv = "single";
        if (k.key == "strand.tightness") nv = "loose";   // custom: a label only
        if (nv == cur) nv = "4";
        const bool set = setSetting(b, k.key, nv, &err);
        if (!set) std::printf("    refused: %s = %s (%s)\n", k.key.c_str(), nv.c_str(), err.c_str());
        CHECK(set);
        const bool importKey = k.key.rfind("import.", 0) == 0;
        // The workflows: the pile.* / material.pile* keys are no inputs under Knit; under Fur they change their stage
        // and the ones below. The workflow itself is no line of a stage's text: under Fur the pile keys enter the
        // Fibres stage and below.
        const bool furKey = k.key.rfind("pile.", 0) == 0 || k.key.rfind("material.pile", 0) == 0;
        const bool workflowKey = k.key == "workflow";
        // fibres.preset and strand.thickness_preset name the values (fibre_presets.h,
        // thickness_presets.h), strand.tightness the derived stitch size (custom here): never inputs.
        const bool labelKey = k.key == "fibres.preset" || k.key == "strand.thickness_preset" || k.key == "strand.tightness";
        // strand.compress_to_fit is hashed as the compression it gives: none at the default gauge (thick_yarn_warning).
        const bool fitKey = k.key == "strand.compress_to_fit";
        // The parameters of overlapped seams (per-component coverage is the default) are inputs only
        // once they are chosen.
        const bool unusedKey = k.key == "knit.overlap_stitches" || k.key == "knit.gap_patches" || k.key == "knit.hand.seed" ||
                               k.key == "knit.hand.loop_height" || k.key == "knit.hand.column_rows" || k.key.rfind("fibres.islands.", 0) == 0 || k.key == "knit.hand.spot_length" || k.key == "knit.hand.spot_change" ||
                               (k.key.rfind("yarn_body.", 0) == 0 && k.key != "yarn_body.enabled") ||  // the yarn body's, while it is off
                               k.key == "strand.max_flattening" ||  // a round yarn's (max_flattening_setting)
                               k.key == "knit.route" || k.key.rfind("knit.guided_", 0) == 0;  // the guided route's: inputs of an import on that route only
        const std::vector<uint64_t> kk = keysOf(b, asset);
        for (int st = 0; st < 5; ++st) {
            const bool affected = workflowKey ? st >= static_cast<int>(Stage::Fibres)
                                              : !importKey && !furKey && !labelKey && !fitKey && !unusedKey && k.stage != Stage::Display &&
                                                    st >= static_cast<int>(k.stage);
            CHECK((kk[st] != base[st]) == affected);
        }
        if (furKey) {
            SceneSettings af = a, bf = b;
            af.workflow = bf.workflow = "fur";
            const std::vector<uint64_t> ka = keysOf(af, asset), kb = keysOf(bf, asset);
            for (int st = 0; st < 5; ++st) CHECK((ka[st] != kb[st]) == (st >= static_cast<int>(k.stage)));
        }
        // import.plain_nodes enters through import.identity, which depends on it (isEffectiveInput).
        if (importKey && k.key != "import.plain_nodes") {  // with scene = import, they enter every stage
            SceneSettings ai = a, bi = b;
            ai.scene = bi.scene = "import";
            if (k.key == "import.wrap_stitch_mm")
                ai.importSurface = bi.importSurface = "wrapped";  // no input of a keep-file import
            const std::vector<uint64_t> ka = keysOf(ai, asset), kb = keysOf(bi, asset);
            for (int st = 0; st < 5; ++st) CHECK(ka[st] != kb[st]);
        }
    }
    // The asset hash enters every stage.
    const std::vector<uint64_t> other = keysOf(a, asset + 1);
    for (int st = 0; st < 5; ++st) CHECK(other[st] != base[st]);
}

TEST(fnv1a_reference) {
    // Published FNV-1a 64 test vectors.
    CHECK(fnv1a("", 0) == 0xcbf29ce484222325ull);
    CHECK(fnv1a("a", 1) == 0xaf63dc4c8601ec8cull);
    CHECK(fnv1a("foobar", 6) == 0x85944171f73967e8ull);
}

TEST(texture_profile_parse) {
    int n = 0;
    PlyProfileModel m = PlyProfileModel::Circles;
    std::string err;
    CHECK(parseTextureProfile("none", n, m, &err) && n == 1);
    CHECK(parseTextureProfile("3-swept", n, m, &err) && n == 3 && m == PlyProfileModel::Swept);
    CHECK(parseTextureProfile("4-circles", n, m, &err) && n == 4 && m == PlyProfileModel::Circles);
    for (const char* bad : {"1-swept", "9-swept", "3", "3-", "x-swept", "3-round", "3-swept-x", ""}) CHECK(!parseTextureProfile(bad, n, m, &err));
}

// Wool presets (fibre_presets.h, viewer plan P7): each preset assigns every value it owns (the same
// settings whatever was chosen before), round-trips through matchingFibrePreset, converts to
// fibre parameters (segments, stray share, the stray population's own lengths), and its label is
// never part of a key; fibres.max_segments takes 2..48; a settings
// file written before the presets keeps its fibre values (the historical defaults), not the
// viewer's preset.
// The two hand-made groups load, an all-zero Off is what default settings match, and only Off is off.
TEST(hand_presets_load_with_an_off) {
    for (const PresetGroup group : {PresetGroup::HandKnit, PresetGroup::FibreIslands}) {
        for (const std::string& e : presetErrors(group)) std::printf("    not loaded: %s\n", e.c_str());
        CHECK(presetErrors(group).empty() && presets(group).size() == 4);
        const FilePreset* off = matchingFilePreset(group, SceneSettings{});
        CHECK(off && off->id == "off" && off == &presets(group).front());
        for (const FilePreset& p : presets(group)) {
            SceneSettings s;
            applyPreset(group, s, p);
            CHECK(matchingFilePreset(group, s) == &p);
            CHECK((s.hand.inUse() || s.islands.inUse()) == (p.id != "off"));
        }
    }
    PresetFile pf;
    FilePreset out;
    std::string err;
    const std::string head = "knit-preset v1\nid = t\nlabel = T\norder = 5\n";
    CHECK(parsePresetFile(head + "fibres.islands.rate = 1\n", pf, &err) && !presetFromFile(PresetGroup::FibreIslands, pf, out, &err) &&
          err.find("missing") != std::string::npos);
    const std::string island = "fibres.islands.rate = 1\nfibres.islands.size = 10\nfibres.islands.strength = 40\nfibres.islands.fuzzy = 1\n"
                               "fibres.islands.matted = 1\nfibres.islands.worn = 1\n";
    CHECK(parsePresetFile(head + island, pf, &err) && presetFromFile(PresetGroup::FibreIslands, pf, out, &err));
    CHECK(parsePresetFile(head + "strays = loose\n" + island, pf, &err) && !presetFromFile(PresetGroup::FibreIslands, pf, out, &err) &&
          err.find("only a fibre preset") != std::string::npos);
}

// The stray presets (presets/strays, owner 2026-11-10): six kinds from Aligned to Wild, ordered by how
// far the strays deviate; Loose is what the built-in settings hold. Each round-trips through
// matchingFilePreset; every wool preset's link resolves and is applied with it; choosing another kind
// leaves the wool preset's name, as editing a wool value leaves the stray preset's.
TEST(stray_presets_load_and_link) {
    const std::vector<FilePreset>& strays = presets(PresetGroup::Strays);
    for (const std::string& e : presetErrors(PresetGroup::Strays)) std::printf("    not loaded: %s\n", e.c_str());
    CHECK(presetErrors(PresetGroup::Strays).empty() && strays.size() == 6);
    const char* ids[] = {"aligned", "along", "loose", "long_tail", "scattered", "wild"};
    for (size_t i = 0; i < strays.size() && i < 6; ++i) CHECK(strays[i].id == ids[i]);
    const FilePreset* loose = matchingFilePreset(PresetGroup::Strays, SceneSettings{});
    CHECK(loose && loose->id == "loose");
    for (size_t i = 0; i < strays.size(); ++i) {
        SceneSettings s;
        applyPreset(PresetGroup::Strays, s, strays[i]);
        CHECK(matchingFilePreset(PresetGroup::Strays, s) == &strays[i]);
        if (i > 0) CHECK(strays[i].values.randomSpreadDeg >= strays[i - 1].values.randomSpreadDeg);
    }
    const std::vector<FilePreset>& wool = presets(PresetGroup::Fibres);
    CHECK(!wool.empty() && presetErrors(PresetGroup::Fibres).empty());
    for (const FilePreset& w : wool) {
        const FilePreset* linked = presetWithId(PresetGroup::Strays, w.strays);
        CHECK(linked != nullptr);
        if (!linked) continue;
        SceneSettings s;
        applyPreset(PresetGroup::Fibres, s, w);
        CHECK(matchingFilePreset(PresetGroup::Fibres, s) == &w && matchingFilePreset(PresetGroup::Strays, s) == linked);
        for (const FilePreset& other : strays) {
            SceneSettings t = s;
            applyPreset(PresetGroup::Strays, t, other);
            CHECK(matchingFilePreset(PresetGroup::Fibres, t) == &w && matchingFilePreset(PresetGroup::Strays, t) == &other);
        }
        SceneSettings u = s;
        u.random.lengthMax += 1;  // a wool value
        CHECK(matchingFilePreset(PresetGroup::Fibres, u) == nullptr && matchingFilePreset(PresetGroup::Strays, u) == linked);
    }
    // The pairs (owner: every wool preset comes with its stray preset).
    for (const FilePreset& w : wool) {
        if (w.id == "standard" || w.id == "soft" || w.id == "look_h") CHECK(w.strays == "loose");
        if (w.id == "hard") CHECK(w.strays == "aligned");
        if (w.id == "fluffy") CHECK(w.strays == "long_tail");
        if (w.id == "mohair") CHECK(w.strays == "loose");
    }
    CHECK(presetWithId(PresetGroup::Fibres, "mohair") != nullptr);
}

// fibres.segment_budget holds per object (a bear): each part's threshold is its object's, so a
// second bear leaves the first one's fibres as they were.
TEST(object_budgets_are_independent) {
    Scene t;
    t.source.prims.resize(3);
    t.source.prims[0].nodeIndex = t.source.prims[1].nodeIndex = 0;  // two pieces of the first bear
    t.source.prims[2].nodeIndex = 1;                                 // the second bear
    t.parts.resize(3);
    for (int i = 0; i < 3; ++i) t.parts[static_cast<size_t>(i)].prim = i;
    CHECK(sceneObjectOf(t.source, 1) == 0 && sceneObjectOf(t.source, 2) == 1);
    // n fibres of 4 segments each, with keys spread over (0, 1).
    const auto fibresOf = [](int n, int seed) {
        FibreSet fs;
        for (int i = 0; i < n; ++i) {
            Fibre f;
            f.key = std::fmod(0.6180339887 * (i + 1) + 0.1 * seed, 1.0);
            fs.addFibre(f, std::vector<Vec3>(5, Vec3{}));
        }
        return fs;
    };
    std::vector<PartFibres> parts(3);
    for (int i = 0; i < 3; ++i) {
        parts[static_cast<size_t>(i)].part = i;
        parts[static_cast<size_t>(i)].fs = fibresOf(10, i);
    }
    const long long budget = 40;  // 10 fibres: the first bear's 20 exceed it
    std::vector<std::pair<float, int>> bear1, bear2;
    for (int i : {0, 1}) addKeySegments(parts[static_cast<size_t>(i)].fs, bear1);
    addKeySegments(parts[2].fs, bear2);
    applyObjectBudgets(t, budget, parts);
    const float q1 = parts[0].keyBelow;
    CHECK(q1 < 1.0f && parts[1].keyBelow == q1 && q1 == budgetThreshold(bear1, budget));
    CHECK(parts[2].keyBelow == budgetThreshold(bear2, budget));
    // A much denser second bear: its threshold drops, the first bear's stays.
    parts[2].fs = fibresOf(40, 2);
    applyObjectBudgets(t, budget, parts);
    CHECK(parts[0].keyBelow == q1 && parts[1].keyBelow == q1 && parts[2].keyBelow < 1.0f);
}

// The yarn fibre looks (presets/yarn_fibres): every shipped file loads; a look sets every value it owns, whatever
// the settings held before; the settings then match it by value; the blue swatch preset is the Smooth look.
TEST(yarn_fibres_presets_assign_all_values) {
    const std::vector<FilePreset>& looks = presets(PresetGroup::YarnFibres);
    for (const std::string& e : presetErrors(PresetGroup::YarnFibres)) std::printf("    not loaded: %s\n", e.c_str());
    CHECK(looks.size() >= 2 && presetErrors(PresetGroup::YarnFibres).empty());
    SceneSettings scrambled;
    scrambled.yarnBody.plies = 1;
    scrambled.yarnBody.fibres = 77;
    scrambled.yarnBody.irregularity = 0.3;
    scrambled.yarnBody.mapRoughness = 0.9;
    scrambled.yarnWhite.fibreIor = 1.7;
    for (const FilePreset& p : looks) {
        SceneSettings a, b = scrambled;
        applyPreset(PresetGroup::YarnFibres, a, p);
        applyPreset(PresetGroup::YarnFibres, b, p);
        CHECK(a.yarnBody == b.yarnBody && a.yarnWhite == b.yarnWhite);
        CHECK(matchingFilePreset(PresetGroup::YarnFibres, a) == &p);
        CHECK(a.yarnBody.plies == 3);  // the standard yarn
        std::printf("    %s: %d fibres of %.0f um, map roughness %.2f\n", p.label.c_str(), a.yarnBody.fibres, a.yarnBody.fibreDiameterUm, a.yarnBody.mapRoughness);
    }
    CHECK(presetWithId(PresetGroup::YarnFibres, "smooth") != nullptr && presetWithId(PresetGroup::YarnFibres, "soft_wool") != nullptr);
    // The map material (yarn_body.map_material) is a look's value: every shipped look is fabric, the default (owner,
    // 2026-11-13); Mohair is the fine, lustrous, loosely twisted look in a deep shell.
    const FilePreset* mohair = presetWithId(PresetGroup::YarnFibres, "mohair");
    const FilePreset* smooth = presetWithId(PresetGroup::YarnFibres, "smooth");
    CHECK(mohair != nullptr && smooth != nullptr && SceneSettings{}.yarnBody.mapMaterial == "fabric");
    for (const FilePreset& p : looks) CHECK(p.values.yarnBody.mapMaterial == "fabric");
    if (mohair && smooth) {
        const SceneSettings::YarnBody &m = mohair->values.yarnBody, &sm = smooth->values.yarnBody;
        CHECK(m.fibreDiameterUm < sm.fibreDiameterUm && m.shellDepth > sm.shellDepth && m.fibreAngleDeg < sm.fibreAngleDeg && m.plyAngleDeg < sm.plyAngleDeg);
        SceneSettings shiny;
        applyPreset(PresetGroup::YarnFibres, shiny, *mohair);
        shiny.yarnBody.mapMaterial = "shiny";
        CHECK(matchingFilePreset(PresetGroup::YarnFibres, shiny) == nullptr);  // the material is part of the look
    }
    {
        // A look file with another map material is refused, naming the key.
        const SceneSettings d;
        std::string body;
        for (const std::string& k : presetKeys(PresetGroup::YarnFibres)) {
            const std::string all = writeSettings(d);
            const size_t at = all.find(k + " = ");
            CHECK(at != std::string::npos);
            if (at != std::string::npos) body += all.substr(at, all.find('\n', at) - at + 1);
        }
        const std::string head = "knit-preset v1\nid = t\nlabel = Test\norder = 5\nnote = n\n";
        PresetFile pf;
        FilePreset fp;
        std::string err;
        CHECK(parsePresetFile(head + body, pf, &err) && presetFromFile(PresetGroup::YarnFibres, pf, fp, &err) && fp.values.yarnBody.mapMaterial == "fabric");
        std::string glossy = body;
        glossy.replace(glossy.find("yarn_body.map_material = fabric"), 31, "yarn_body.map_material = glossy");
        CHECK(parsePresetFile(head + glossy, pf, &err) && !presetFromFile(PresetGroup::YarnFibres, pf, fp, &err) &&
              err.find("yarn_body.map_material") != std::string::npos);
        std::printf("    refused: %s\n", err.c_str());
    }
    const SceneEntry* blue = findScene("wool-swatch-blue");
    CHECK(blue != nullptr);
    if (blue) {
        SceneSettings s;
        applySceneSettings(s, *blue);
        const FilePreset* match = matchingFilePreset(PresetGroup::YarnFibres, s);
        std::printf("    blue swatch: %s\n", match ? match->label.c_str() : "Custom");
        CHECK(match && match->id == "smooth");
    }
}

TEST(fibre_presets_assign_all_values) {
    // The presets are the files of presets/fibres (the user edits them): every shipped file loads.
    const std::vector<FilePreset>& wool = presets(PresetGroup::Fibres);
    for (const std::string& e : presetErrors(PresetGroup::Fibres)) std::printf("    not loaded: %s\n", e.c_str());
    CHECK(!wool.empty() && presetErrors(PresetGroup::Fibres).empty());
    SceneSettings scrambled;
    scrambled.randomFraction = 0.77;
    scrambled.directed.lengthMax = 9;
    scrambled.directed.curvature = 2.5;
    scrambled.random.thicknessSigma = 0.2;
    scrambled.random.maxAngleDeg = 33;
    scrambled.randomSpreadDeg = 77;
    scrambled.fibreMaxSegments = 7;
    scrambled.fibreDirected.fibreLongRoughness = 0.7;
    scrambled.fibreHairColour = scrambled.fibreHairColour == "yarn" ? "layer" : "yarn";
    for (const FilePreset& p : wool) {
        SceneSettings a, b = scrambled;
        applyPreset(PresetGroup::Fibres, a, p);
        applyPreset(PresetGroup::Fibres, b, p);
        CHECK(writeSettings(a) == writeSettings(b));
        CHECK(matchingFilePreset(PresetGroup::Fibres, a) == &p);
        FibreParams fp;
        std::string err;
        CHECK(fibreParamsFromSettings(a, fp, &err));
        CHECK(fp.maxSegments == p.values.fibreMaxSegments && fp.randomFraction == p.values.randomFraction && fp.lengthMax == p.values.directed.lengthMax);
        // The stray population's own lengths (an override only where they differ from the directed ones).
        if (p.values.random.lengthMax != p.values.directed.lengthMax) CHECK(fp.randomLengthMax == p.values.random.lengthMax);
        if (p.values.random.lengthMedian != p.values.directed.lengthMedian) CHECK(fp.randomLengthMean == p.values.random.lengthMedian);
        // The strays deviate as the profile's own stray preset says (strays = <id>, applied with it).
        const FilePreset* strays = presetWithId(PresetGroup::Strays, p.strays);
        CHECK(strays != nullptr && matchingFilePreset(PresetGroup::Strays, a) == strays);
        if (strays)
            CHECK(fp.randomSpreadDeg == strays->values.randomSpreadDeg && fp.randomElevMinDeg == 20 &&
                  fp.randomElevMaxDeg == strays->values.random.maxAngleDeg && a.random.lengthSpread == strays->values.random.lengthSpread);
        // A stray value is the stray preset's: editing it keeps the wool preset and makes the strays Custom.
        SceneSettings c = a;
        c.randomSpreadDeg += 1;
        CHECK(matchingFilePreset(PresetGroup::Fibres, c) == &p && matchingFilePreset(PresetGroup::Strays, c) == nullptr);
        c = a;
        c.random.lengthMedian += 0.01;
        CHECK(matchingFilePreset(PresetGroup::Fibres, c) == nullptr);
        c = a;
        c.fibreDirected.fibreLongRoughness += 0.01;
        CHECK(matchingFilePreset(PresetGroup::Fibres, c) == nullptr);
        c = a;
        c.fibreHairColour = a.fibreHairColour == "yarn" ? "layer" : "yarn";
        CHECK(matchingFilePreset(PresetGroup::Fibres, c) == nullptr);
        c = a;
        c.fibrePreset = "custom";  // the label: no key
        CHECK(fibreResultKey(c) == fibreResultKey(a) && stageKey(c, Stage::Export, 0) == stageKey(a, Stage::Export, 0));
    }
    // The built-in settings (no preset) have the Loose strays (owner, 2026-11-10): within 45 deg of the
    // yarn, rising 20 to 40 deg, length spread 0.5.
    {
        FibreParams fromDefaults;
        std::string err;
        CHECK(fibreParamsFromSettings(SceneSettings{}, fromDefaults, &err));
        CHECK(fromDefaults.randomSpreadDeg == 45 && fromDefaults.randomElevMinDeg == 20 && fromDefaults.randomElevMaxDeg == 40 &&
              SceneSettings{}.random.lengthSpread == 0.5);
    }
    // A preset file (preset_files.h): every owned value, nothing else, values the generator accepts.
    const SceneSettings d;
    std::string body;
    for (const std::string& k : presetKeys(PresetGroup::Fibres)) {
        const std::string all = writeSettings(d);
        const size_t at = all.find(k + " = ");
        body += all.substr(at, all.find('\n', at) - at + 1);
    }
    // A wool preset names its stray preset; the other kinds share the plain header.
    const std::string head = "knit-preset v1\n# test\nid = t\nlabel = Test\norder = 5\nnote = n\n", woolHead = head + "strays = loose\n";
    PresetFile pf;
    FilePreset fpre;
    std::string err;
    CHECK(parsePresetFile(woolHead + body, pf, &err) && pf.id == "t" && pf.order == 5 && pf.strays == "loose" &&
          pf.values.size() == presetKeys(PresetGroup::Fibres).size());
    CHECK(presetFromFile(PresetGroup::Fibres, pf, fpre, &err) && fpre.values.directed == d.directed && fpre.values.random == d.random &&
          fpre.values.fibreMaxSegments == 3 && fpre.strays == "loose");
    CHECK(parsePresetFile(head + body, pf, &err) && !presetFromFile(PresetGroup::Fibres, pf, fpre, &err) && err.find("strays: required") != std::string::npos);
    CHECK(parsePresetFile(head + "strays = nowhere\n" + body, pf, &err) && !presetFromFile(PresetGroup::Fibres, pf, fpre, &err) &&
          err.find("no such stray preset") != std::string::npos);
    CHECK(parsePresetFile(woolHead + body + "fibres.random.spread_deg = 30\n", pf, &err) && !presetFromFile(PresetGroup::Fibres, pf, fpre, &err) &&
          err.find("fibres.random.spread_deg is not a value of a fibre preset") != std::string::npos);
    const std::string missing = body.substr(body.find('\n') + 1);  // without fibres.random_fraction
    CHECK(parsePresetFile(woolHead + missing, pf, &err) && !presetFromFile(PresetGroup::Fibres, pf, fpre, &err) && err.find("missing fibres.random_fraction") != std::string::npos);
    CHECK(parsePresetFile(woolHead + body + "fibres.density = 3\n", pf, &err) && !presetFromFile(PresetGroup::Fibres, pf, fpre, &err) &&
          err.find("fibres.density") != std::string::npos);
    std::string tooMany = body;
    tooMany.replace(tooMany.find("fibres.max_segments = 3"), 23, "fibres.max_segments = 49");
    CHECK(parsePresetFile(woolHead + tooMany, pf, &err) && !presetFromFile(PresetGroup::Fibres, pf, fpre, &err) && err.find("max_segments") != std::string::npos);
    // Every profile but Hard wool has hair-shaded fuzz and curled strays; Hard wool keeps plain strays and
    // matte fuzz. Standard wool and Look H stay the opaque, whole-length-tapered reference
    // coat, Look H with the knit car's stray share. Soft wool is that coat with the study's finish (MAT-22,
    // 2026-11-04): a late taper at the coat's own coverage (thinner fibres), softer relief and a
    // translucent yarn. Fluffy wool's strays reach far and lead its look; Standard's stay rarer and shorter.
    const FilePreset* standard = nullptr;
    const FilePreset* soft = nullptr;
    const FilePreset* fluffy = nullptr;
    const FilePreset* lookH = nullptr;
    const FilePreset* mohair = nullptr;
    double lowestOtherRoughness = 1e9;
    for (const FilePreset& p : wool) {
        const bool hard = p.id == "hard";
        CHECK((p.values.random.curlRadius > 0) != hard && (p.values.fibreDirected.model == MaterialSpec::Model::Fibre) != hard);
        if (p.id == "standard") standard = &p;
        if (p.id == "soft") soft = &p;
        if (p.id == "fluffy") fluffy = &p;
        if (p.id == "look_h") lookH = &p;
        if (p.id == "mohair") mohair = &p;
        else lowestOtherRoughness = std::min(lowestOtherRoughness, p.values.fibreDirected.fibreLongRoughness);
    }
    // Mohair (2026-11-13, from the owner's photos): the fuzz is the halo, longer and finer than Fluffy wool's, with
    // more strays that are only gently wavy (a wider, shallower curl than Fluffy's); the most lustrous hair, in the
    // yarn's own colour (a halo lighter than the layer rule's dark net); fine tips, no crimp; a translucent yarn.
    CHECK(mohair != nullptr && fluffy != nullptr);
    if (mohair && fluffy) {
        const SceneSettings &m = mohair->values, &f = fluffy->values;
        CHECK(m.directed.lengthMedian > f.directed.lengthMedian && m.directed.lengthMax > f.directed.lengthMax &&
              m.directed.thickness < f.directed.thickness && m.randomFraction > f.randomFraction);
        CHECK(m.random.curlRadius > 0 && m.random.curlRadius < f.random.curlRadius && m.random.curlLength > f.random.curlLength);
        CHECK(m.fibreDirected.fibreLongRoughness < lowestOtherRoughness && m.fibreHairColour == "yarn");
        CHECK(m.crimpAmplitude == 0 && m.tailTaper > 0 && m.yarnWhite.translucency >= 0.4 && m.normalStrength < 1);
        // Its longest fibres keep pieces of segment_length (max_segments does not coarsen them).
        CHECK(m.fibreMaxSegments * m.fibreSegmentLength >= std::max(m.directed.lengthMax, m.random.lengthMax));
    }
    CHECK(standard != nullptr && soft != nullptr && fluffy != nullptr && lookH != nullptr);
    if (standard && soft && fluffy && lookH) {
        for (const FilePreset* ref : {standard, lookH})
            CHECK(ref->values.tailTaper == 0 && ref->values.yarnWhite.translucency == 0 && ref->values.normalStrength == 1 &&
                  ref->values.crimpAmplitude == 0);
        CHECK(lookH->values.randomFraction == 0.1 && lookH->values.directed == standard->values.directed &&
              lookH->values.random == standard->values.random);
        CHECK(soft->values.tailTaper > 0 && soft->values.yarnWhite.translucency > 0 && soft->values.normalStrength < 1 &&
              soft->values.directed.thickness < standard->values.directed.thickness &&
              soft->values.directed.lengthMedian == standard->values.directed.lengthMedian);
        // Standard wool's coat keeps the yarn's colour, matte; Fluffy wool's layer keeps the mean rule.
        CHECK(standard->values.fibreHairColour == "yarn" && fluffy->values.fibreHairColour == "layer" &&
              standard->values.fibreDirected.fibreLongRoughness > fluffy->values.fibreDirected.fibreLongRoughness);
        CHECK(fluffy->values.random.lengthMax >= 40);
        CHECK(standard->values.randomFraction < fluffy->values.randomFraction && standard->values.random.lengthMedian < fluffy->values.random.lengthMedian &&
              standard->values.random.lengthMax < fluffy->values.random.lengthMax);
    }
    CHECK(!parsePresetFile("id = x\n", pf, &err) && !parsePresetFile(head + body + "fibres.max_segments = 4\n", pf, &err));
    CHECK(!parsePresetFile("knit-preset v1\nlabel = L\norder = 1\n", pf, &err) && err.find("id") != std::string::npos);
    // The strand mesh levels' files load as well (presets/strand_mesh).
    for (const std::string& e : presetErrors(PresetGroup::StrandMesh)) std::printf("    not loaded: %s\n", e.c_str());
    CHECK(!presets(PresetGroup::StrandMesh).empty() && presetErrors(PresetGroup::StrandMesh).empty());
    FilePreset sm;
    CHECK(parsePresetFile(head + "strand.tube_sides = 2\nstrand.tube_ring = 0.4\nstrand.ply_tube_sides = 6\nstrand.ply_tube_ring = 0.3\n", pf, &err) &&
          !presetFromFile(PresetGroup::StrandMesh, pf, sm, &err) && err.find("tube_sides") != std::string::npos);
    // The density presets' files load (presets/density); Dense is the default density (the teddy
    // renders' 125 per mm); a file sets fibres.density only, within 1..1000.
    for (const std::string& e : presetErrors(PresetGroup::Density)) std::printf("    not loaded: %s\n", e.c_str());
    CHECK(!presets(PresetGroup::Density).empty() && presetErrors(PresetGroup::Density).empty());
    const FilePreset* dense = nullptr;
    for (const FilePreset& p : presets(PresetGroup::Density)) {
        SceneSettings a;
        applyPreset(PresetGroup::Density, a, p);
        CHECK(matchingFilePreset(PresetGroup::Density, a) == &p && a.fibreDensity == p.values.fibreDensity);
        if (p.id == "dense") dense = &p;
    }
    CHECK(dense && dense->values.fibreDensity == SceneSettings{}.fibreDensity && matchingFilePreset(PresetGroup::Density, SceneSettings{}) == dense);
    SceneSettings off;
    off.fibreDensity = 124;
    CHECK(matchingFilePreset(PresetGroup::Density, off) == nullptr);
    FilePreset dp;
    CHECK(parsePresetFile(head + "fibres.density = 60\n", pf, &err) && presetFromFile(PresetGroup::Density, pf, dp, &err) && dp.values.fibreDensity == 60 && dp.id == "t");
    for (const char* v : {"0", "1001"})
        CHECK(parsePresetFile(head + "fibres.density = " + v + "\n", pf, &err) && !presetFromFile(PresetGroup::Density, pf, dp, &err) &&
              err.find("fibres.density") != std::string::npos);
    CHECK(parsePresetFile(head + "fibres.density = 60\nfibres.random_fraction = 0.2\n", pf, &err) && !presetFromFile(PresetGroup::Density, pf, dp, &err) &&
          err.find("fibres.random_fraction") != std::string::npos);
    CHECK(parsePresetFile(head, pf, &err) && !presetFromFile(PresetGroup::Density, pf, dp, &err) && err.find("missing fibres.density") != std::string::npos);
    // Segments per fibre.
    SceneSettings m;
    FibreParams fp;
    for (const int n : {1, 49}) {
        m.fibreMaxSegments = n;
        CHECK(!fibreParamsFromSettings(m, fp, &err) && err.find("max_segments") != std::string::npos);
    }
    m.fibreMaxSegments = 48;
    CHECK(fibreParamsFromSettings(m, fp, &err) && fp.maxSegments == 48);
    m.fibreMaxSegments = 3;
    CHECK(stageKey(m, Stage::Fibres, 0) == stageKey(SceneSettings{}, Stage::Fibres, 0));
    m.fibreMaxSegments = 4;
    CHECK(stageKey(m, Stage::Fibres, 0) != stageKey(SceneSettings{}, Stage::Fibres, 0));
    // New objects take Soft wool (MAT-22, owner 2026-11-04): the lowest order, with its tip 0.3.
    const SceneEntry* teddy = findScene("teddy");
    CHECK(teddy != nullptr);
    SceneSettings v = viewerDefaults(*teddy);
    CHECK(wool[0].id == "soft" && v.fibrePreset == "soft" && v.tailTaperTip == 0.3 && v.tailTaper == 0.2 &&
          v.normalStrength == 0.6 && v.yarnWhite.translucency == 0.3);
    CHECK(matchingFilePreset(PresetGroup::Fibres, v) == &wool[0]);
    // A file saved with Standard keeps Standard after the default moved.
    if (standard) {
        SceneSettings saved = v, reread = viewerDefaults(*teddy);
        applyPreset(PresetGroup::Fibres, saved, *standard);
        CHECK(readSettings(writeSettings(saved), reread, &err) && reread.fibrePreset == "standard" &&
              matchingFilePreset(PresetGroup::Fibres, reread) == standard && reread.tailTaper == 0 &&
              reread.yarnWhite.translucency == 0);
    }
    // A file written before the presets (no fibres.preset key, no population keys) loaded over the
    // viewer's defaults: the values of that time, written out (setPrePresetFibreValues), never the
    // built-in values of today (the Loose strays since 2026-11-10).
    std::string old;
    {
        std::stringstream in(writeSettings(SceneSettings{}));
        std::string line;
        while (std::getline(in, line))
            if (line.rfind("fibres.preset", 0) != 0 && line.rfind("fibres.directed.", 0) != 0 && line.rfind("fibres.random", 0) != 0 &&
                line.rfind("fibres.max_segments", 0) != 0)
                old += line + "\n";
    }
    CHECK(readSettings(old, v, &err));
    CHECK(v.directed == d.directed && v.randomFraction == 0.5 && v.fibreMaxSegments == 3 && matchingFilePreset(PresetGroup::Fibres, v) == nullptr);
    CHECK(v.random.maxAngleDeg == 60 && v.random.lengthSpread == 0.35 && v.random.curvature == 0.3 && v.randomSpreadDeg == 180);
    CHECK(v.random != SceneSettings{}.random && v.randomSpreadDeg != SceneSettings{}.randomSpreadDeg);
    // A current file keeps its preset values, the stray ones too.
    SceneSettings w = viewerDefaults(*teddy), back = SceneSettings{};
    CHECK(readSettings(writeSettings(w), back, &err) && matchingFilePreset(PresetGroup::Fibres, back) == &wool[0] && back.random == w.random &&
          back.randomSpreadDeg == w.randomSpreadDeg);
    // A file written before fibres.random.spread_deg existed: its strays left the yarn in any direction.
    std::string noSpread;
    {
        std::stringstream in(writeSettings(w));
        std::string line;
        while (std::getline(in, line))
            if (line.rfind("fibres.random.spread_deg", 0) != 0) noSpread += line + "\n";
    }
    SceneSettings x = viewerDefaults(*teddy);
    CHECK(readSettings(noSpread, x, &err) && x.randomSpreadDeg == 180 && x.random == w.random);
}

// ------------------------------------------------------- teddy pieces ----

// The teddy (user, 2026-10-28): the same workflow and settings as every scene, knitted as its S1
// spiral pieces (its UV charts knit with holes). The spiral knitting reads the
// scene settings directly: at the defaults the pre-S6 command line's knitting (README S5: --knit
// --settle --repair --surface-half), the common diameter plan, the shared knitting settings in its
// knit cache key (clearance, depth amplitude, stitch size, with the bears at their fixed size at
// every stitch width), and it refuses settings the pieces cannot honour instead of ignoring them.
// The pieces' scene colours each bear by its GLB fabric.
TEST(teddy_options_from_settings) {
    const SceneEntry* e = findScene("teddy");
    CHECK(e != nullptr && e->knitsAsPieces());
    knit_test::requireScene(findProjectFile(e->piecesManifest));
    SceneSettings s;
    s.scene = "teddy";
    applySceneDefaults(s, *e);
    s.fibres = true;
    std::string err;
    CHECK(checkSpiralSettings(*e, s, &err));
    CHECK(s.settle && s.repair && s.surfaceGate && s.surfaceHalf);
    CHECK(s.wale == e->piecesW && s.course == e->piecesH && s.diameter == e->piecesD);
    // The diameters are the common plan (the scene's Dense preset and its fallback chain), each at its
    // knitting diameter (a compressed yarn knits at its contact fraction).
    const DiameterPlan plan = partDiameters(s);
    CHECK(plan.attempt.size() == s.diameterFallback.size() + 1 && plan.attempt[0] == knitDiameter(s, s.diameter));
    for (size_t k = 0; k + 1 < plan.attempt.size(); ++k) CHECK(plan.attempt[k + 1] == knitDiameter(s, s.diameterFallback[k]));
    // The bears keep their size at every stitch width: the gauge stitch count is not rounded, so
    // sigma (scene mm per unit, the scale of every piece) is the same at 5.5, 6.0 and 7.5 mm.
    SpiralSource src0;
    CHECK(loadSpiralSource(*e, s, src0, &err));
    for (double w : {5.5, 6.0, 7.5}) {
        SceneSettings t = s;
        t.wale = w;
        SpiralSource src;
        CHECK(loadSpiralSource(*e, t, src, &err));
        std::printf("    w %.1f mm: sigma %.6f mm per unit\n", w, src.sigma);
        CHECK_NEAR(src.sigma, src0.sigma, 1e-12 * src0.sigma);
    }
    // The shared knitting settings reach the knitting: they are inputs of its knit cache key.
    {
        const std::optional<SpiralKnitCache> c0 = spiralKnitCache(s, src0, "", (std::filesystem::temp_directory_path() / "knit_s6_key").string(), &err);
        SceneSettings t = s;
        t.clearance = 0.15;
        t.depthAmplitude = 0.7;
        const std::optional<SpiralKnitCache> c1 = spiralKnitCache(t, src0, "", (std::filesystem::temp_directory_path() / "knit_s6_key").string(), &err);
        CHECK(c0 && c1 && !c0->key.empty() && c0->key != c1->key && c0->loadDir.find("gpu-f32") != std::string::npos);
    }
    // Refusals (not wired / not meaningful on spiral pieces).
    auto refused = [&](std::function<void(SceneSettings&)> edit) {
        SceneSettings t = s;
        edit(t);
        const bool r = !checkSpiralSettings(*e, t, &err) && !err.empty();
        if (r) std::printf("    refused: %s\n", err.c_str());
        return r;
    };
    CHECK(refused([](SceneSettings& t) { t.pattern = "rib1x1"; }));
    CHECK(refused([](SceneSettings& t) { t.angleDeg = 10; }));
    CHECK(refused([](SceneSettings& t) { t.yarnColourSource = "texture"; }));
    CHECK(refused([](SceneSettings& t) {
        t.flatten = true;
        t.plies = 2;
    }));
    CHECK(refused([](SceneSettings& t) { t.yarnWhite.model = MaterialSpec::Model::Fibre; }));
    CHECK(refused([](SceneSettings& t) { t.diameterAttempts = -1; }));
    CHECK(refused([](SceneSettings& t) { t.directed.maxAngleDeg = 5; }));
    // The pieces' scene: with yarn colour = object each bear's pieces carry its GLB fabric's base colour.
    std::map<std::string, Vec3> fabric;
    CHECK(piecesObjectColours(*e, fabric, &err));
    CHECK(s.yarnColourSource == "object");
    const Vec3 white = fabric.at("fabric_white"), brown = fabric.at("fabric_brown");
    std::printf("    fabric_white %.3f %.3f %.3f, fabric_brown %.3f %.3f %.3f\n", white.x, white.y, white.z, brown.x, brown.y, brown.z);
    CHECK(!(white == brown));
    Scene t;
    SpiralSource src;
    CHECK(loadPiecesScene(*e, s, src, t, &err));
    int whiteParts = 0, brownParts = 0;
    for (const ScenePart& r : t.parts) {
        const SourcePrimitive& g = t.source.prims[static_cast<size_t>(r.prim)];
        CHECK(fabric.count(g.material) && g.baseColor == fabric.at(g.material));
        whiteParts += g.material == "fabric_white";
        brownParts += g.material == "fabric_brown";
    }
    std::printf("    %zu pieces: %d white, %d brown\n", t.parts.size(), whiteParts, brownParts);
    CHECK(whiteParts > 0 && brownParts > 0 && whiteParts + brownParts == static_cast<int>(t.parts.size()));
}

// The teddy's Dense comes from one initialization path: the command line's order (the scene
// defaults) and the viewer's give the same gauge; Final's fallback is the calibration's scaled to the
// drawn 1.86 (from 1.80), Draft has none.
TEST(teddy_defaults_one_path) {
    const SceneEntry* teddy = findScene("teddy");
    CHECK(teddy);
    SceneSettings cli;
    cli.scene = "teddy";
    applySceneDefaults(cli, *teddy);
    const SceneSettings view = viewerDefaults(*teddy);
    std::printf("    cli D %.2f at %.1f x %.1f, fallback %zu; viewer D %.2f at %.1f x %.1f, preset %s\n", cli.diameter, cli.wale, cli.course,
                cli.diameterFallback.size(), view.diameter, view.wale, view.course, view.thicknessPreset.c_str());
    for (const SceneSettings* x : std::initializer_list<const SceneSettings*>{&cli, &view}) {
        CHECK(x->diameter == 1.86 && x->wale == 5.5 && x->course == 3.4 && x->thicknessPreset == "thick-dense");
        CHECK(knitWorkflow(*x) && x->pattern == "stockinette" && x->yarnColourSource == "object");
    }
    const std::vector<double> fallback{1.80, 1.73, 1.67, 1.60};
    CHECK(cli.diameterFallback == fallback);
    CHECK(isPreviewProfile(view) && view.diameterFallback.empty());
    SceneSettings fin = view;
    applyFullProfile(fin, *teddy);
    CHECK(fin.diameterFallback == fallback && isFinalProfile(fin, *teddy));
}

// The teddy's scene cache identity covers what its knit reads: the pieces, the GLB (colours, the kept
// parts, the skipped-knit mesh) and the routing; a changed piece is a new identity.
TEST(teddy_cache_identity) {
    const SceneEntry* e = findScene("teddy");
    CHECK(e != nullptr);
    knit_test::requireScene(findProjectFile(e->piecesManifest));
    std::string err;
    const uint64_t h = sceneAssetHash(*e, &err);
    std::printf("  teddy asset hash %016llx\n", static_cast<unsigned long long>(h));
    // The tracked pieces' bytes (af949d0), the routing v5 (the kept parts from the GLB's selection,
    // 2026-11-11) and bears.glb in teddy cream and brown (2026-10-31): a changed piece, a line-ending
    // conversion or a lost file changes the teddy's identity, and with it every cached teddy knit.
    // Re-pinned 2026-11-11: the pieces re-exported with pole A at the lower vertex index (6 pieces and
    // the manifest changed; teddy_twin_pieces_chart_alike).
    CHECK(h == 0xca283980b1796e54ull);
    // The same GLB without the pieces (the teddy of 4195804): another identity, so its UV-chart
    // entries are never read for the pieces.
    SceneEntry glbOnly = *e;
    glbOnly.piecesManifest.clear();
    CHECK(sceneAssetHash(glbOnly, &err) != h);
    // A copy of the pieces with one piece changed: another identity.
    namespace fs = std::filesystem;
    const fs::path src = fs::path(findProjectFile(e->piecesManifest)).parent_path();
    const fs::path tmp = fs::temp_directory_path() / "knit_teddy_identity";
    std::error_code ec;
    fs::remove_all(tmp, ec);
    fs::copy(src, tmp, fs::copy_options::recursive, ec);
    CHECK(!ec);
    SceneEntry copy = *e;
    copy.piecesManifest = (tmp / "manifest.txt").string();
    CHECK(sceneAssetHash(copy, &err) == h);  // the same bytes: the same identity
    {
        std::ofstream f(tmp / e->gaugePiece, std::ios::app);
        f << "# changed\n";
    }
    CHECK(sceneAssetHash(copy, &err) != h);
    fs::remove_all(tmp, ec);
}

// The two bears are instanced copies of one mesh: each bear.009 piece has the chart of the bear.002 piece
// of its island (same triangles, same (a, b) per corner), so both bears knit alike. Pole A is the lower
// vertex index (export_pieces.py); before that rule bear.009's _04 and _12 had a flipped chart and
// another seam, and its _12 knitted 446 stitches against 443.
TEST(teddy_twin_pieces_chart_alike) {
    const SceneEntry* e = findScene("teddy");
    CHECK(e != nullptr);
    const std::string path = findProjectFile(e->piecesManifest);
    knit_test::requireScene(path);
    Manifest man;
    std::string err;
    CHECK(readManifest(path, man, &err));
    int pairs = 0;
    for (const PieceEntry& b : man.pieces) {
        if (b.object != "bear.009") continue;
        const auto twin = std::find_if(man.pieces.begin(), man.pieces.end(),
                                       [&](const PieceEntry& a) { return a.object == "bear.002" && a.island == b.island; });
        CHECK(twin != man.pieces.end());
        TriMesh ma, mb;
        CHECK(readObj(piecePath(man, twin->file), ma, &err) && readObj(piecePath(man, b.file), mb, &err));
        CHECK(ma.tri == mb.tri);
        bool sameChart = ma.uv.size() == mb.uv.size();
        for (size_t t = 0; sameChart && t < ma.uv.size(); ++t)
            for (int k = 0; k < 3; ++k) sameChart = sameChart && ma.uv[t][k].x == mb.uv[t][k].x && ma.uv[t][k].y == mb.uv[t][k].y;
        if (!sameChart) std::printf("  %s and %s chart apart\n", twin->file.c_str(), b.file.c_str());
        CHECK(sameChart);
        ++pairs;
    }
    CHECK(pairs == 9);
}

// ------------------------------------------------------------- patterns ----

// The pre-S6 template table, frozen here (stockinette.cpp at f42a519), so a change of
// the stockinette generator is caught even if both sides of a comparison change.
std::vector<TemplateControl> frozenTemplate(const PatternParams& pp) {
    const double Xh = pp.headHalfWidth, Xl = pp.legHalfWidth;
    const TemplateControl right[9] = {{0.0, 1.40, -1.0}, {0.64 * Xh, 1.38, -1.0}, {Xh + 0.02, 1.24, -1.0},
                                      {Xh + 0.02, 0.95, +1.0}, {0.80 * Xh, 0.60, +1.0}, {Xl, 0.30, +1.0},
                                      {Xl, -0.05, -1.0}, {0.5 * (Xl + 0.5), -0.30, -1.0}, {0.5, -0.36, -1.0}};
    std::vector<TemplateControl> t(kControlsPerStitch);
    for (int j = 8; j < 16; ++j) t[j] = right[j - 8];
    for (int j = 0; j < 8; ++j) t[j] = {-right[8 - j].X, right[8 - j].Y, right[8 - j].Z};
    return t;
}

// Stockinette output is bit-identical to the pre-S6 generator: every control equals the
// frozen template placed by the (unchanged) templateToLocal, with ==, not a tolerance.
TEST(pattern_stockinette_bit_identical) {
    PatternParams dense;
    dense.waleSpacing = 5.5;
    dense.courseSpacing = 3.4;
    dense.yarnDiameter = 1.16;
    dense.wales = 7;
    dense.courses = 5;
    for (const PatternParams& pp : {PatternParams{}, dense}) {
        const std::vector<TemplateControl> ft = frozenTemplate(pp);
        const std::vector<TemplateControl> t0 = stitchTemplate(pp), t1 = stitchTemplate(pp, StitchFaces{});
        for (int j = 0; j < kControlsPerStitch; ++j) {
            CHECK(t0[j].X == ft[j].X && t0[j].Y == ft[j].Y && t0[j].Z == ft[j].Z);
            CHECK(t1[j].X == ft[j].X && t1[j].Y == ft[j].Y && t1[j].Z == ft[j].Z);
        }
        const CurveSet cs = generateStockinette(pp);
        const int K = kControlsPerStitch;
        for (const YarnPath& y : cs.yarns) {
            std::vector<Vec3> c = y.curve.ctrl;
            if (y.reversed) std::reverse(c.begin(), c.end());
            CHECK(static_cast<int>(c.size()) == pp.wales * K + 3);
            for (int g = -1; g <= pp.wales * K + 1; ++g) {
                const int i = g >= 0 ? g / K : -1, j = g - i * K;
                const Vec3 e = templateToLocal(ft[j], i, y.course, pp);
                CHECK(c[g + 1].x == e.x && c[g + 1].y == e.y && c[g + 1].z == e.z);
            }
        }
        for (const StitchRecord& r : cs.topo.stitches) CHECK(r.face == StitchFace::Knit);
    }
}

// Reverse stockinette is the exact depth mirror (z negated, x and y unchanged), every
// stitch labelled purl; the validator's verdict is recorded (PASS at the expanded preset).
TEST(pattern_reverse_is_mirror) {
    PatternParams pp;
    const CurveSet k = generateStockinette(pp);
    pp.pattern = StitchPattern::ReverseStockinette;
    const CurveSet p = generateStockinette(pp);
    CHECK(k.yarns.size() == p.yarns.size());
    for (size_t y = 0; y < k.yarns.size(); ++y)
        for (size_t c = 0; c < k.yarns[y].curve.ctrl.size(); ++c) {
            const Vec3 &a = k.yarns[y].curve.ctrl[c], &b = p.yarns[y].curve.ctrl[c];
            CHECK(a.x == b.x && a.y == b.y && a.z == -b.z);
        }
    for (const StitchRecord& r : p.topo.stitches) CHECK(r.face == StitchFace::Purl);
    const ValidationReport r = validate(p, pp);
    CHECK(r.verdict == Status::Pass);
}

// The purl expectation is the knit one mirrored, written from the construction (not
// from the template), and is chosen by the CHILD's face label: geometry and labels are
// checked against each other, including face-swapped negative controls.
TEST(purl_signature_expectation) {
    using namespace signature;
    for (int q = 0; q < 4; ++q) {
        CHECK(kPurlNewerAlongNewer[q] != kNewerAlongNewer[q]);
        CHECK(kPurlNewerAlongOlder[q] != kNewerAlongOlder[q]);
    }
    auto signatureWith = [](CurveSet cs, const PatternParams& pp, StitchFace (*relabel)(StitchFace)) {
        for (StitchRecord& r : cs.topo.stitches) r.face = relabel(r.face);
        return validate(cs, pp, {.generatorChecks = false}).signature.status;
    };
    auto keep = [](StitchFace f) { return f; };
    auto allPurl = [](StitchFace) { return StitchFace::Purl; };
    auto allKnit = [](StitchFace) { return StitchFace::Knit; };
    PatternParams pp;
    pp.wales = pp.courses = 4;
    const CurveSet knit = generateStockinette(pp);
    pp.pattern = StitchPattern::ReverseStockinette;
    const CurveSet purl = generateStockinette(pp);
    CHECK(signatureWith(knit, pp, keep) == Status::Pass);
    CHECK(signatureWith(purl, pp, keep) == Status::Pass);
    CHECK(signatureWith(knit, pp, allPurl) == Status::Fail);  // knit geometry, purl labels
    CHECK(signatureWith(purl, pp, allKnit) == Status::Fail);  // purl geometry, knit labels
    // Transition cases: rib (wale transitions inside each row) and garter (course
    // transitions), both faces; one child's label swapped must fail exactly its pair.
    for (StitchPattern pat : {StitchPattern::Rib1x1, StitchPattern::Rib2x2, StitchPattern::Garter}) {
        PatternParams q;
        q.pattern = pat;
        q.wales = q.courses = 6;
        const CurveSet cs = generateStockinette(q);
        CHECK(signatureWith(cs, q, keep) == Status::Pass);
        for (int target : {2 * 6 + 2, 3 * 6 + 3}) {  // one knit and one purl child for each pattern
            CurveSet sw = cs;
            StitchRecord& r = sw.topo.stitches[target];
            r.face = r.face == StitchFace::Knit ? StitchFace::Purl : StitchFace::Knit;
            const ValidationReport v = validate(sw, q, {.generatorChecks = false});
            int failing = 0, failingTarget = 0;
            for (const PairSignature& ps : v.signature.pairs)
                if (ps.status == Status::Fail) {
                    ++failing;
                    failingTarget += ps.child == target;
                }
            CHECK(failing == 1 && failingTarget == 1);
        }
    }
}

// Recorded flat reference verdicts (S6 measurement): verdicts and failure categories,
// not wording. A change in either direction fails and forces a re-review.
TEST(pattern_flat_references) {
    struct Row { StitchPattern p; double w, h, D; Status v; unsigned fail; };
    const unsigned cl = kFailClearance, cr = kFailCrossing;
    const Row rows[] = {
        {StitchPattern::Stockinette, 8, 6.5, 1.0, Status::Pass, 0},
        {StitchPattern::ReverseStockinette, 8, 6.5, 1.0, Status::Pass, 0},
        {StitchPattern::Rib1x1, 8, 6.5, 1.0, Status::Pass, 0},
        {StitchPattern::Rib2x2, 8, 6.5, 1.0, Status::Pass, 0},
        {StitchPattern::Garter, 8, 6.5, 1.0, Status::Pass, 0},          // with the searched transition
        {StitchPattern::Stockinette, 5.5, 3.4, 1.10, Status::Pass, 0},  // thick (dense) preset
        {StitchPattern::Rib1x1, 5.5, 3.4, 1.10, Status::Pass, 0},
        {StitchPattern::Rib2x2, 5.5, 3.4, 1.10, Status::Pass, 0},
        {StitchPattern::ReverseStockinette, 5.5, 3.4, 1.10, Status::Pass, 0},
        {StitchPattern::Stockinette, 5.5, 3.4, 1.16, Status::Fail, cl},  // passes on curved teddy pieces only after settle
        {StitchPattern::Garter, 5.5, 3.4, 1.0, Status::Fail, cl},
        {StitchPattern::Garter, 5.5, 5.0, 1.10, Status::Pass, 0},        // garter needs taller rows
    };
    for (const Row& row : rows) {
        PatternParams pp;
        pp.waleSpacing = row.w;
        pp.courseSpacing = row.h;
        pp.yarnDiameter = row.D;
        const PatternReference ref = flatPatternReference(pp, row.p);
        if (ref.verdict != row.v || ref.fail != row.fail) {
            char b[256];
            std::snprintf(b, sizeof b, "%s at %s: %s (fail %s), recorded %s (fail %s)", toString(row.p), ref.params.c_str(),
                          toString(ref.verdict), failCategoryText(ref.fail).c_str(), toString(row.v), failCategoryText(row.fail).c_str());
            throw std::runtime_error(b);
        }
    }
    // The naive garter transition (the plain template's values) fails at the loose gauge:
    // the depth flip at the leg / parent-crown crossing.
    PatternParams naive;
    naive.transNeckY = 0.60;
    naive.transNeckZ = 1.0;
    naive.transLegY = 0.30;
    const PatternReference g = flatPatternReference(naive, StitchPattern::Garter);
    CHECK(g.verdict == Status::Fail && g.fail == (cl | cr));
}

// The relaxed size of each pattern: the measured ratios (Hwang et al. 2025), the
// chart stitch they make of the settings' stockinette stitch, and a yarn fit that no longer depends
// on the pattern (garter had its own, wrong-way fit limits: rows x 1.53 taller).
TEST(pattern_rest_scale) {
    const PatternScale st = patternRestScale(StitchPattern::Stockinette), rev = patternRestScale(StitchPattern::ReverseStockinette);
    const PatternScale g = patternRestScale(StitchPattern::Garter), r1 = patternRestScale(StitchPattern::Rib1x1),
                       r2 = patternRestScale(StitchPattern::Rib2x2);
    CHECK(st.w == 1 && st.h == 1 && rev.w == 1 && rev.h == 1);
    CHECK(g.w == 1.10 && g.h == 0.70);  // garter: wider and SHORTER than stockinette
    CHECK(r1.w == 0.74 && r1.h == 1.23 && r2.w == 0.65 && r2.h == 1.12);
    SceneSettings s;
    s.wale = 5.5;
    s.course = 3.4;
    s.diameter = 1.86;
    const StitchSize sg = patternStitch(s, StitchPattern::Garter);
    CHECK(std::fabs(sg.w - 6.05) < 1e-12 && std::fabs(sg.h - 2.38) < 1e-12);
    s.pattern = "stockinette";
    const double kStockinette = yarnFit(s).k;
    s.pattern = "garter";
    CHECK(yarnFit(s).k == kStockinette);
    // yarn mode: the stitch size is the stockinette stitch for every pattern.
    s.stitchSize = "yarn";
    s.pattern = "stockinette";
    applyStitchSize(s);
    const double w0 = s.wale, h0 = s.course;
    s.pattern = "garter";
    applyStitchSize(s);
    CHECK(s.wale == w0 && s.course == h0);
    // The shape applied where the knit's PatternParams are made.
    PatternParams pp;
    applyPatternShape(pp, s, StitchPattern::Rib1x1);
    CHECK(pp.pattern == StitchPattern::Rib1x1 && pp.waleSpacing == s.wale * 0.74 && pp.courseSpacing == s.course * 1.23);
    CHECK(pp.ribLayerOffset == kRibLayerOffset1x1 && pp.headHalfWidth == PatternParams{}.headHalfWidth * kRibMaxWiden);
    PatternParams gp;
    applyPatternShape(gp, s, StitchPattern::Garter);
    CHECK(gp.garterWave == kGarterWave && gp.ribLayerOffset == 0 && gp.headHalfWidth == PatternParams{}.headHalfWidth);
    PatternParams sp;
    applyPatternShape(sp, s, StitchPattern::Stockinette);
    PatternParams plain;
    plain.waleSpacing = s.wale;
    plain.courseSpacing = s.course;
    CHECK(sp == plain);  // stockinette: nothing but its stitch
}

// Rib layers (PatternParams::ribLayerOffset): every control of a knit wale moves ribLayerOffset x D to
// the front, of a purl wale to the back, except the sinker between a knit and a purl wale; the crossing
// depths at every interlock keep their sign. Other patterns ignore the field.
TEST(pattern_rib_layers) {
    PatternParams pp;
    pp.pattern = StitchPattern::Rib1x1;
    const double layer = 1.3;
    for (const StitchFaces f : {StitchFaces{StitchFace::Knit, StitchFace::Knit, StitchFace::Purl},
                                StitchFaces{StitchFace::Purl, StitchFace::Purl, StitchFace::Knit},
                                StitchFaces{StitchFace::Knit, StitchFace::Knit, StitchFace::Knit}}) {
        pp.ribLayerOffset = 0;
        const std::vector<TemplateControl> flat = stitchTemplate(pp, f);
        pp.ribLayerOffset = layer;
        const std::vector<TemplateControl> layered = stitchTemplate(pp, f);
        const double sigma = f.own == StitchFace::Knit ? 1 : -1;
        for (int j = 0; j < kControlsPerStitch; ++j) {
            CHECK(layered[j].X == flat[j].X && layered[j].Y == flat[j].Y);
            const double expect = j == 0 && f.left != f.own ? flat[j].Z : flat[j].Z + sigma * layer / pp.depthAmplitude;
            CHECK(std::fabs(layered[j].Z - expect) < 1e-15);
        }
    }
    CHECK(std::fabs(pp.centrelineDepth() - (pp.depthAmplitude + layer) * pp.yarnDiameter) < 1e-15);
    PatternParams g = pp;
    g.pattern = StitchPattern::Garter;  // not a rib: the field changes nothing
    PatternParams g0 = g;
    g0.ribLayerOffset = 0;
    CHECK(generateStockinette(g).yarns[2].curve.ctrl == generateStockinette(g0).yarns[2].curve.ctrl);
    CHECK(g.centrelineDepth() == g.depthAmplitude * g.yarnDiameter);
}

// The garter wave (PatternParams::garterWave) is a depth that depends on the position in the fabric
// only: x and y of every control stay, and every crossing keeps its depth order, so the crossing
// signature is the flat garter's. At the relaxed garter size (6.05 x 2.38 mm, D 1.164) it raises the
// fold margin of the bare template (the scene knit's settle does the rest).
TEST(pattern_garter_wave) {
    PatternParams pp;
    pp.pattern = StitchPattern::Garter;
    pp.wales = pp.courses = 6;
    const CurveSet flat = generateStockinette(pp);
    pp.garterWave = 0.6;
    const CurveSet waved = generateStockinette(pp);
    bool zMoved = false;
    for (size_t y = 0; y < flat.yarns.size(); ++y)
        for (size_t c = 0; c < flat.yarns[y].curve.ctrl.size(); ++c) {
            const Vec3 &a = flat.yarns[y].curve.ctrl[c], &b = waved.yarns[y].curve.ctrl[c];
            CHECK(a.x == b.x && a.y == b.y);
            zMoved = zMoved || a.z != b.z;
        }
    CHECK(zMoved);
    CHECK(validate(waved, pp, {.generatorChecks = false}).signature.status == Status::Pass);
    PatternParams garter;
    garter.pattern = StitchPattern::Garter;
    garter.wales = garter.courses = 8;
    garter.waleSpacing = 6.05;
    garter.courseSpacing = 2.38;
    garter.yarnDiameter = 1.164;
    garter.transNeckY = kGarterNeckY;
    garter.transNeckZ = kGarterNeckZ;
    garter.transLegY = kGarterLegY;
    auto foldMargin = [](const PatternParams& q) {
        const ValidationReport r = validate(generateStockinette(q), q);
        double rho = 1e300;
        for (const CurvatureResult& c : r.curvature) rho = std::min(rho, c.rhoLo);
        return rho;
    };
    const double without = foldMargin(garter);
    garter.garterWave = kGarterWave;
    CHECK(foldMargin(garter) > without);
}

// The rib layer overlap rule (signature.cpp ribLayerNeighbours): with loops wider than the wale
// spacing, a rib's neighbouring knit and purl stitches overlap in the front projection. The crossings
// are expected with the knit loop in front, and fail with the purl loop in front (layers swapped).
TEST(signature_rib_layer_overlap) {
    PatternParams pp;
    pp.pattern = StitchPattern::Rib1x1;
    pp.wales = pp.courses = 8;
    pp.waleSpacing = 4.07;
    pp.courseSpacing = 4.18;
    pp.yarnDiameter = 1.164;
    pp.headHalfWidth = 0.375 / 0.74;  // stockinette's loop width: wider than the rib's wale
    pp.legHalfWidth = 0.125 / 0.74;
    pp.ribLayerOffset = 1.3;
    const ValidationReport front = validate(generateStockinette(pp), pp);
    CHECK(front.verdict == Status::Pass);
    CHECK(!front.signature.layerOverlap.empty() && front.signature.layerWrongOrder == 0 && front.signature.unexpected.empty());
    pp.ribLayerOffset = -1.3;  // the purl wales in front of the knit wales
    const ValidationReport swapped = validate(generateStockinette(pp), pp, {.generatorChecks = false});
    CHECK(swapped.signature.layerWrongOrder > 0 && swapped.signature.status == Status::Fail);
}

// ------------------------------------------------------------ GLB / mesh ----

TEST(glb_load_transforms) {
    const std::string dir = std::filesystem::temp_directory_path().string() + "/knit_s6_tiny";
    const std::string path = writeTinyGltf(dir);
    SceneEntry e;
    e.id = "tiny";
    e.glb.sizeRef = GlbSelection::SizeRef::Height;
    e.glb.sizeMm = 100;
    e.glb.knitAllMaterials = true;
    SourceScene g;
    std::string err;
    CHECK(loadGlb(path, e.glb, g, &err));
    CHECK(g.prims.size() == 2);
    // Independent expectation of the glTF world positions (Y-up).
    auto parentChild = [](Vec3 p) {
        p = p * 2.0 + Vec3{0, 1, 0};                 // child TRS
        return Vec3{p.z, p.y, -p.x} + Vec3{1, 2, 3};  // rotation +90 deg about Y, then translation
    };
    auto mirror = [](Vec3 p) { return Vec3{-p.x, p.y, p.z} + Vec3{5, 0, 0}; };
    const Vec3 src[3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
    std::vector<Vec3> world;
    for (const Vec3& p : src) world.push_back(parentChild(p));
    for (const Vec3& p : src) world.push_back(mirror(p));
    // Library: zUp (x, -z, y), height -> 100 mm, centred, floor at z = 0.
    double zlo = 1e300, zhi = -1e300, xlo = 1e300, xhi = -1e300, ylo = 1e300, yhi = -1e300;
    for (const Vec3& w : world) {
        const Vec3 z{w.x, -w.z, w.y};
        zlo = std::min(zlo, z.z), zhi = std::max(zhi, z.z), xlo = std::min(xlo, z.x), xhi = std::max(xhi, z.x);
        ylo = std::min(ylo, z.y), yhi = std::max(yhi, z.y);
    }
    const double s = 100 / (zhi - zlo);
    CHECK_NEAR(g.scale, s, 1e-12);
    for (int k = 0; k < 2; ++k) {
        const SourcePrimitive& p = g.prims[k];
        CHECK(p.mirrored == (k == 1));
        CHECK(p.welded.pos.size() == 3 && p.welded.tri.size() == 1);
        CHECK(p.baseColor.x == 0.5 && p.baseColor.y == 0.25 && p.baseColor.z == 0.125);
        // Each source vertex is found at its expected library position.
        for (int v = 0; v < 3; ++v) {
            const Vec3& w = world[3 * k + v];
            const Vec3 lib = Vec3{w.x, -w.z, w.y} * s + Vec3{-0.5 * (xlo + xhi), -0.5 * (ylo + yhi), -zlo} * s;
            double best = 1e300;
            for (const Vec3& q : p.welded.pos) best = std::min(best, norm(q - lib));
            CHECK(best < 1e-9);
        }
        // Winding: the face normal is the transformed source normal (+Z in glTF: the
        // normal matrix of the mirror keeps +Z), mapped to Z-up; mirroring must not flip it.
        const Vec3 n = p.welded.faceNormal(0);
        const Vec3 nWorld = k == 0 ? Vec3{1, 0, 0} : Vec3{0, 0, 1};  // R_y(90) maps +Z to +X
        const Vec3 nLib{nWorld.x, -nWorld.z, nWorld.y};
        CHECK(dot(n, nLib) > 1 - 1e-9);
        // Corner uvs travel with the corners (swapped together on the mirrored node).
        for (int c = 0; c < 3; ++c) {
            const Vec2 uv = p.welded.uv[0][c];
            const Vec3 q = p.welded.pos[p.welded.tri[0][c]];
            int sv = -1;
            for (int v = 0; v < 3; ++v) {
                const Vec3& w = world[3 * k + v];
                const Vec3 lib = Vec3{w.x, -w.z, w.y} * s + Vec3{-0.5 * (xlo + xhi), -0.5 * (ylo + yhi), -zlo} * s;
                if (norm(lib - q) < 1e-9) sv = v;
            }
            CHECK(sv >= 0 && uv.x == (sv == 1 ? 1.0 : 0.0) && uv.y == (sv == 2 ? 1.0 : 0.0));
        }
    }
    CHECK_NEAR(g.hi.z - g.lo.z, 100, 1e-9);
    CHECK_NEAR(g.lo.z, 0, 1e-12);
    // An unlisted material is refused (no silent default).
    e.glb.knitAllMaterials = false;
    e.glb.knitMaterials = {"other"};
    CHECK(!loadGlb(path, e.glb, g, &err) && err.find("wool") != std::string::npos);
}

// Welding: only within one source, never across sources; corner uvs kept, so a uv seam
// splits charts while the welded positions join.
// A file that requires an extension is refused by its name. The message used to be built after the
// file was freed, so it read freed memory.
TEST(glb_required_extension_is_named) {
    const std::string dir = std::filesystem::temp_directory_path().string() + "/knit_s6_required_extension";
    std::filesystem::create_directories(dir);
    writeTinyGltf(dir);
    std::string text;
    {
        std::ifstream in(dir + "/tiny.gltf");
        text.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    }
    const std::string required = R"("extensionsUsed":["KHR_draco_mesh_compression"],"extensionsRequired":["KHR_draco_mesh_compression"],)";
    text.insert(text.find('{') + 1, required);
    {
        std::ofstream out(dir + "/tiny.gltf");
        out << text;
    }
    SceneEntry e;
    e.glb.sizeMm = 10;
    e.glb.knitAllMaterials = true;
    SourceScene g;
    std::string err;
    CHECK(!loadGlb(dir + "/tiny.gltf", e.glb, g, &err));
    std::printf("  %s\n", err.c_str());
    CHECK(err.find("requires extension KHR_draco_mesh_compression") != std::string::npos);
}

TEST(glb_weld_keeps_corner_uv) {
    // A 2 x 1 quad strip whose middle column is duplicated with two different uvs (a seam).
    TriMesh m;
    m.pos = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {1, 0, 0}, {2, 0, 0}, {2, 1, 0}, {1, 1, 0}};
    m.tri = {{0, 1, 2}, {0, 2, 3}, {4, 5, 6}, {4, 6, 7}};
    m.uv = {{Vec2{0, 0}, Vec2{1, 0}, Vec2{1, 1}}, {Vec2{0, 0}, Vec2{1, 1}, Vec2{0, 1}},
            {Vec2{5, 0}, Vec2{6, 0}, Vec2{6, 1}}, {Vec2{5, 0}, Vec2{6, 1}, Vec2{5, 1}}};
    const WeldedMesh welded = weldVertices(m, 1e-9);
    const TriMesh& w = welded.mesh;
    const std::vector<int>& remap = welded.remap;
    CHECK(w.pos.size() == 6 && w.tri.size() == 4 && welded.droppedTriangles == 0);
    CHECK(remap[4] == remap[1] && remap[7] == remap[2]);
    for (size_t t = 0; t < 4; ++t)
        for (int c = 0; c < 3; ++c) CHECK(w.uv[t][c].x == m.uv[t][c].x && w.uv[t][c].y == m.uv[t][c].y);
    const MeshCheck mc = checkMesh(w);
    CHECK(mc.components.size() == 1 && mc.components[0].uvCharts == 2);
    // Two coincident sources stay separate: the loader welds per (node, primitive).
    const std::string dir = std::filesystem::temp_directory_path().string() + "/knit_s6_tiny2";
    std::filesystem::create_directories(dir);
    {
        writeTinyGltf(dir);
        std::ifstream in(dir + "/tiny.gltf");
        std::string j((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        // The mirror node becomes a second, identical instance at the child's place.
        const std::string from = R"({"name":"mirror","translation":[5,0,0],"scale":[-1,1,1],"mesh":0})";
        const size_t at = j.find(from);
        CHECK(at != std::string::npos);
        j.replace(at, from.size(),
                  R"({"name":"twin","matrix":[0,0,-1,0, 0,1,0,0, 1,0,0,0, 1,2,3,1],"children":[3]},)"
                  R"({"name":"twinchild","translation":[0,1,0],"scale":[2,2,2],"mesh":0})");
        std::ofstream out(dir + "/tiny.gltf");
        out << j;
    }
    SceneEntry e;
    e.glb.sizeMm = 10;
    e.glb.knitAllMaterials = true;
    SourceScene g;
    std::string err;
    CHECK(loadGlb(dir + "/tiny.gltf", e.glb, g, &err));
    CHECK(g.prims.size() == 2);
    for (const SourcePrimitive& p : g.prims) CHECK(p.welded.pos.size() == 3);  // not merged into one 3-vertex set
}

// KHR_texture_transform reaches the colour uv (the man's textures repeat 1.5 and 2 times); the
// knitting chart keeps the mesh's own uv.
TEST(glb_texture_transform_reaches_the_colour_uv) {
    TextureTransform t;
    t.scale = {2, 3};
    t.offset = {0.25, 0.5};
    Vec2 q = transformTextureUv(t, {1, 1});
    CHECK(std::fabs(q.x - 2.25) < 1e-12 && std::fabs(q.y - 3.5) < 1e-12);
    TextureTransform turn;
    turn.rotation = std::acos(0.0);  // the extension's matrix: u' = cos u + sin v, v' = -sin u + cos v
    q = transformTextureUv(turn, {1, 0});
    CHECK(std::fabs(q.x) < 1e-12 && std::fabs(q.y + 1) < 1e-12);
    const std::string dir = std::filesystem::temp_directory_path().string() + "/knit_s6_texture_transform";
    std::string path = writeTinyGltf(dir);
    CHECK(writePngRgb8(dir + "/tiny.png", 2, 2, std::vector<uint8_t>(12, 128)));  // decoded, its colour does not matter
    std::string j;
    {
        std::ifstream in(path);
        j.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    }
    const std::string from = R"("pbrMetallicRoughness":{"baseColorFactor":[0.5,0.25,0.125,1]}}],)";
    const size_t at = j.find(from);
    CHECK(at != std::string::npos);
    j.replace(at, from.size(),
              R"("pbrMetallicRoughness":{"baseColorTexture":{"index":0,"extensions":{"KHR_texture_transform":{"scale":[2,3],"offset":[0.25,0.5]}}}}}],)"
              R"("textures":[{"source":0}],"images":[{"uri":"tiny.png"}],"extensionsUsed":["KHR_texture_transform"],)");
    {
        std::ofstream out(path);
        out << j;
    }
    SceneEntry e;
    e.glb.sizeMm = 10;
    e.glb.knitAllMaterials = true;
    SourceScene g;
    std::string err;
    CHECK(loadGlb(path, e.glb, g, &err));
    CHECK(!g.prims.empty() && g.prims[0].colourUv.size() == 1);
    const SourcePrimitive& p = g.prims[0];
    // Corner (1, 0): (2.25, 0.5) in glTF's convention, v up in the colour uv; corner (0, 1): (0.25, 3.5).
    CHECK(std::fabs(p.colourUv[0][1].x - 2.25) < 1e-6 && std::fabs(p.colourUv[0][1].y - 0.5) < 1e-6);
    CHECK(std::fabs(p.colourUv[0][2].x - 0.25) < 1e-6 && std::fabs(p.colourUv[0][2].y + 2.5) < 1e-6);
    CHECK(p.welded.uv[0][1].x == 1 && p.welded.uv[0][1].y == 0 && p.welded.uv[0][2].x == 0 && p.welded.uv[0][2].y == 1);
}

TEST(triangle_intersection_cases) {
    const Vec3 a[3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
    const Vec3 apart[3] = {{0, 0, 1}, {1, 0, 1}, {0, 1, 1}};
    const Vec3 pierce[3] = {{0.2, 0.2, -1}, {0.3, 0.2, 1}, {0.2, 0.3, 1}};
    const Vec3 touch[3] = {{0.2, 0.2, 0}, {0.3, 0.2, 1}, {0.2, 0.3, 1}};       // one vertex on a's plane, inside a
    const Vec3 coplanar[3] = {{0.2, 0.2, 0}, {1.2, 0.2, 0}, {0.2, 1.2, 0}};    // overlapping in the plane
    const Vec3 coplanarApart[3] = {{2, 2, 0}, {3, 2, 0}, {2, 3, 0}};
    CHECK(triangleIntersection(a, apart, 1e-9) == Status::Pass);
    CHECK(triangleIntersection(a, pierce, 1e-9) == Status::Fail);
    CHECK(triangleIntersection(a, touch, 1e-9) == Status::Unresolved);
    CHECK(triangleIntersection(a, coplanar, 1e-9) == Status::Fail);
    CHECK(triangleIntersection(a, coplanarApart, 1e-9) == Status::Pass);
}

// A closed, consistently oriented mesh with positive volume is NOT a valid solid if it
// intersects itself: not eligible for inside tests (negative control), while the
// undeformed cube is (positive control).
TEST(solid_eligibility_self_intersecting) {
    TriMesh cube;
    cube.pos = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    cube.tri = {{0, 2, 1}, {0, 3, 2}, {4, 5, 6}, {4, 6, 7}, {0, 1, 5}, {0, 5, 4},
                {1, 2, 6}, {1, 6, 5}, {2, 3, 7}, {2, 7, 6}, {3, 0, 4}, {3, 4, 7}};
    const MeshCheck ok = checkMesh(cube);
    CHECK(ok.components.size() == 1);
    const ComponentCheck& c0 = ok.components[0];
    CHECK(c0.closed && c0.consistent && c0.signedVolume > 0 && c0.selfIntersection == Status::Pass && c0.insideEligible);
    CHECK_NEAR(c0.signedVolume, 1.0, 1e-12);
    TriMesh bad = cube;
    bad.pos[6] = {0.5, 0.5, -0.5};  // the top corner pushed through the bottom face
    const MeshCheck mc = checkMesh(bad);
    const ComponentCheck& c = mc.components[0];
    CHECK(c.closed && c.consistent);
    CHECK(c.signedVolume > 0);
    CHECK(c.selfIntersection == Status::Fail && c.intersectingPairs > 0);
    CHECK(!c.insideEligible && c.reason.find("self-intersection") != std::string::npos);
    // Open and flipped variants are not eligible either.
    TriMesh open = cube;
    open.tri.pop_back();
    CHECK(!checkMesh(open).components[0].insideEligible);
    TriMesh flipped = cube;
    for (auto& t : flipped.tri) std::swap(t[1], t[2]);
    const ComponentCheck& cf = checkMesh(flipped).components[0];
    CHECK(cf.closed && cf.consistent && cf.signedVolume < 0 && !cf.insideEligible);
}

// --------------------------------------------------------- UV-chart tiling ----

std::string failing(const PartKnit& k) {
    std::string o = k.detail;
    for (const CheckLine& l : k.after.lines)
        if (l.status != Status::Pass) o += "\n    [" + std::string(toString(l.status)) + "] " + l.name + ": " + l.detail;
    int shown = 0;
    for (const ClearanceRegion& g : k.after.clearance.regions) {
        if (g.status == Status::Pass || ++shown > 6) continue;
        const KnitTopology& t = k.t.curves.topo;
        const int sa = t.stitchAt(g.yarnA, 0.5 * (g.uA0 + g.uA1)), sb = t.stitchAt(g.yarnB, 0.5 * (g.uB0 + g.uB1));
        char b[300];
        std::snprintf(b, sizeof b, "\n    clearance %s: yarn %d u[%.3f,%.3f] (stitch %d, i %d r %d) vs yarn %d u[%.3f,%.3f] (stitch %d, i %d r %d) d %.4f",
                      toString(g.status), g.yarnA, g.uA0, g.uA1, sa, sa >= 0 ? t.stitches[sa].lattice.i : -1,
                      sa >= 0 ? t.stitches[sa].lattice.r : -1, g.yarnB, g.uB0, g.uB1, sb, sb >= 0 ? t.stitches[sb].lattice.i : -1,
                      sb >= 0 ? t.stitches[sb].lattice.r : -1, g.worstDistance);
        o += b;
    }
    return o;
}

// A bilinear lookup interpolates in linear light (the midpoint of black and white is linear 0.5, sRGB
// 188, not 128), with texel centres at (i + 1/2) / size, repeat wrap and v up.
TEST(texel_bilinear_linear_light) {
    ColourImage across;
    across.width = 2, across.height = 1;
    across.rgba = {0, 0, 0, 255, 255, 255, 255, 255};
    CHECK(bilinearTexel(across, {0.25, 0.5}).x == 0.0);
    CHECK(std::fabs(bilinearTexel(across, {0.75, 0.5}).x - 1.0) < 1e-12);
    const Vec3 middle = bilinearTexel(across, {0.5, 0.5});
    CHECK(std::fabs(middle.x - 0.5) < 1e-12 && stitchColourKey(middle) == 0xbcbcbc);
    CHECK(std::fabs(bilinearTexel(across, {0.0, 0.5}).x - 0.5) < 1e-12);  // between the last texel and the first
    CHECK(bilinearTexel(across, {1.25, 0.5}).x == 0.0);
    CHECK(std::fabs(bilinearTexel(across, {-0.25, 0.5}).x - 1.0) < 1e-12);
    ColourImage down;
    down.width = 1, down.height = 2;
    down.rgba = {255, 255, 255, 255, 0, 0, 0, 255};  // row 0 (the image's top) white
    CHECK(std::fabs(bilinearTexel(down, {0.5, 0.75}).x - 1.0) < 1e-12);
    CHECK(bilinearTexel(down, {0.5, 0.25}).x == 0.0);
    for (uint32_t b = 0; b < 256; ++b) CHECK(stitchColourKey(stitchColourLinear(b << 16 | b << 8 | b)) == (b << 16 | b << 8 | b));
}

// The palette: nothing changes at or below the limit; above it at most the limit, the same on every
// call, and a rare colour far from the rest keeps an entry of its own.
TEST(palette_reduction) {
    const std::map<uint32_t, long long> few{{0xff0000, 10}, {0x00ff00, 3}, {0x123456, 1}};
    const std::map<uint32_t, uint32_t> same = reducePalette(few, 3);
    CHECK(same.size() == 3);
    for (const std::pair<const uint32_t, uint32_t>& entry : same) CHECK(entry.first == entry.second);
    const uint32_t rare = 0x6a008a;
    std::map<uint32_t, long long> many;
    for (uint32_t r = 120; r < 250; r += 3)
        for (uint32_t g = 20; g < 80; g += 3)
            for (uint32_t b = 20; b < 80; b += 3) many[r << 16 | g << 8 | b] = 500;
    many[rare] = 3;
    const std::map<uint32_t, uint32_t> reduced = reducePalette(many, 64);
    CHECK(reduced.size() == many.size() && reduced == reducePalette(many, 64));
    std::set<uint32_t> entries;
    for (const std::pair<const uint32_t, uint32_t>& entry : reduced) entries.insert(entry.second);
    const uint32_t kept = reduced.at(rare);
    auto channelError = [&](int shift) { return std::abs(int(kept >> shift & 0xff) - int(rare >> shift & 0xff)); };
    std::printf("    %zu colours -> %zu, the rare one -> %06x\n", many.size(), entries.size(), kept);
    CHECK(entries.size() <= 64 && entries.size() > 32);
    CHECK(channelError(16) <= 8 && channelError(8) <= 8 && channelError(0) <= 8);
}

// A planar chart with isometric uvs reproduces the flat generator: every row's controls
// equal the flat patch's (unreversed) controls up to ONE translation, and gate (a) passes.
TEST(uv_tiling_plane_equals_flat) {
    const PatternParams pp = expandedD1();
    const TriMesh m = gridMesh(70, 50, 1.0, [](double x, double y) { return Vec3{x, y, 0}; },
                               [](double x, double y) { return Vec2{x / 100, y / 100}; });
    MeshKnitParams p = tiledParams(pp);
    p.settle = false;
    PartKnit k;
    std::string err;
    CHECK(knitMesh(m, p, knit_test::testNumerics(), k, &err));
    const CurveSet& t = k.t.curves;
    CHECK(t.kind == CurveKind::Surface && !t.yarns.empty());
    // All rows have the same length (rectangular chart).
    const size_t n = t.yarns[0].stitches.size();
    for (const YarnPath& y : t.yarns) CHECK(y.stitches.size() == n);
    PatternParams fp = pp;
    fp.wales = static_cast<int>(n);
    fp.courses = static_cast<int>(t.yarns.size());
    const CurveSet f = generateStockinette(fp);
    std::vector<std::vector<Vec3>> fc;
    for (const YarnPath& y : f.yarns) {
        std::vector<Vec3> c = y.curve.ctrl;
        if (y.reversed) std::reverse(c.begin(), c.end());
        fc.push_back(c);
    }
    const Vec3 shift = t.yarns[0].curve.ctrl[0] - fc[0][0];
    double dev = 0, extent = 0;
    for (size_t r = 0; r < t.yarns.size(); ++r) {
        CHECK(t.yarns[r].curve.ctrl.size() == fc[r].size());
        for (size_t c = 0; c < fc[r].size(); ++c) {
            dev = std::max(dev, norm(t.yarns[r].curve.ctrl[c] - fc[r][c] - shift));
            const Vec3& q = t.yarns[r].curve.ctrl[c];
            extent = std::max({extent, std::fabs(q.x), std::fabs(q.y), std::fabs(q.z)});
        }
    }
    // Exact on the CPU. The GPU's finalize stores the controls in float32: then the bound is a few
    // float32 rounding units of the coordinates (2^-24 relative per component, sqrt(3) for 3D).
    const double tol = knit_test::useGpu() ? 4 * std::sqrt(3.0) * std::ldexp(extent, -24) : 1e-9;
    if (!(dev < tol)) throw std::runtime_error("plane vs flat: max control deviation " + std::to_string(dev) + " > " + std::to_string(tol));
    CHECK(k.t.bind.unresolved() == 0);
    CHECK(k.after.verdict == Status::Pass);
    CHECK(k.t.lattice.status == Status::Pass && k.t.layout.status == Status::Pass);
    CHECK(k.boundary.b1 == Status::Pass && k.boundary.b2 == Status::Pass);
    CHECK(k.verdict == Status::Pass);
    CHECK(!k.atlas->charts()[0].upFromWorld);  // the plane's normal is world up: uv' v fallback
}

// Anisotropic affine uvs on a plane: M recovers an isometry (M A orthogonal), every
// singular value is 1, the cells are exactly w x h.
TEST(uv_metric_affine_exact) {
    const double A[4] = {0.02, 0.007, -0.004, 0.031};  // uv = A xy
    const TriMesh m = gridMesh(70, 50, 1.0, [](double x, double y) { return Vec3{x, y, 0}; },
                               [&](double x, double y) { return Vec2{A[0] * x + A[1] * y, A[2] * x + A[3] * y}; });
    const UvAtlas at(m);
    CHECK(at.charts().size() == 1);
    const UvChart& ch = at.charts()[0];
    const Mat2 MA = ch.M * Mat2{A[0], A[1], A[2], A[3]};
    CHECK_NEAR(MA.a * MA.a + MA.c * MA.c, 1, 1e-12);
    CHECK_NEAR(MA.b * MA.b + MA.d * MA.d, 1, 1e-12);
    CHECK_NEAR(MA.a * MA.b + MA.c * MA.d, 0, 1e-12);
    for (int q = 0; q < 5; ++q) {
        CHECK_NEAR(ch.distortion.s1[q], 1, 1e-12);
        CHECK_NEAR(ch.distortion.s2[q], 1, 1e-12);
    }
    CHECK(ch.distortion.outsideFraction == 0);
    MeshKnitParams p = tiledParams(expandedD1());
    p.settle = false;
    PartKnit k;
    std::string err;
    // On the CPU: the 1e-9 below is a double-precision bound, and the device classifies nodes in float32.
    CHECK(knitMesh(m, p, knit_test::cpuNumerics(), k, &err));
    CHECK_NEAR(k.t.layout.minW, 1, 1e-9);
    CHECK_NEAR(k.t.layout.maxW, 1, 1e-9);
    CHECK_NEAR(k.t.layout.minH, 1, 1e-9);
    CHECK_NEAR(k.t.layout.maxH, 1, 1e-9);
}

// Global scale is not distortion: the same uvs on a mesh twice as large give identical
// distortion and the same 3D cell size in mm (twice as many uv' cells).
TEST(uv_scale_whole_chart_invariant) {
    auto make = [](double s) {
        return gridMesh(70, 50, 1.0, [s](double x, double y) { return Vec3{s * x, s * y, 0}; },
                        [](double x, double y) { return Vec2{x / 100 + 0.002 * y / 100, y / 100}; });
    };
    const TriMesh m1 = make(1), m2 = make(2);
    const UvAtlas a1(m1), a2(m2);
    for (int q = 0; q < 5; ++q) {
        CHECK_NEAR(a1.charts()[0].distortion.s1[q], a2.charts()[0].distortion.s1[q], 1e-12);
        CHECK_NEAR(a1.charts()[0].distortion.s2[q], a2.charts()[0].distortion.s2[q], 1e-12);
    }
    MeshKnitParams p = tiledParams(expandedD1());
    p.settle = false;
    PartKnit k1, k2;
    std::string err;
    // On the CPU: the 1e-6 below is a double-precision bound, and the device classifies nodes in float32.
    CHECK(knitMesh(m1, p, knit_test::cpuNumerics(), k1, &err) && knitMesh(m2, p, knit_test::cpuNumerics(), k2, &err));
    for (const PartKnit* k : {&k1, &k2}) {
        CHECK_NEAR(k->t.layout.minW, 1, 1e-6);
        CHECK_NEAR(k->t.layout.maxW, 1, 1e-6);
        CHECK_NEAR(k->t.layout.minH, 1, 1e-6);
        CHECK_NEAR(k->t.layout.maxH, 1, 1e-6);
    }
    CHECK(k2.t.prog.stitches.size() > 3 * k1.t.prog.stitches.size());
}

// Residual local stretch relative to the chart normalization is flagged. Distortion: a
// 24 mm band of a 120 mm plane stretched 1.6x in x has s1 = 1.6 / sqrt(157.44 / 120) =
// 1.397 > 1.3, the rest s2 = 0.873; outside fraction exactly 24 / 120. The layout check
// sees cells inside the band (a band narrower than a stitch is averaged by the cell
// widths: the distortion statistics are the check for that, not the layout).
TEST(uv_local_stretch_flagged) {
    auto u = [](double x) { return x < 48 ? x : x < 72 ? 48 + (x - 48) / 1.6 : 48 + 24 / 1.6 + (x - 72); };
    const TriMesh m = gridMesh(120, 40, 1.0, [](double x, double y) { return Vec3{x, y, 0}; },
                               [&](double x, double y) { return Vec2{u(x) / 100, y / 100}; });
    const UvAtlas at(m);
    const UvChart& ch = at.charts()[0];
    const double g = std::sqrt((96 + 24 * 2.56) / 120);
    CHECK_NEAR(ch.distortion.outsideFraction, 24.0 / 120.0, 1e-12);
    CHECK_NEAR(ch.distortion.s1[4], 1.6 / g, 1e-9);
    CHECK_NEAR(ch.distortion.s1[0], 1.0, 1e-9);      // unstretched y
    CHECK_NEAR(ch.distortion.s2[0], 1.0 / g, 1e-9);  // unstretched x after normalization
    MeshKnitParams p = tiledParams(expandedD1());
    p.settle = false;
    PartKnit k;
    std::string err;
    CHECK(knitMesh(m, p, knit_test::testNumerics(), k, &err));
    CHECK(k.t.layout.status == Status::Fail && k.t.layout.outOfBounds > 0);
    CHECK(k.t.layout.maxW > 1.3);
    CHECK(k.verdict == Status::Fail);
}

// A cylinder (R 20 mm, 48 mm high, axis = world up) cut by a uv seam: one chart, wales
// along the axis (world up), rows around it ending at the seam; gate (a), lattice, layout,
// projection and B1 / B2 pass after settle.
TEST(uv_tiling_cylinder_valid) {
    const double R = 20, H = 48;
    const int na = 160;
    const TriMesh raw = gridMesh(na, 48, 1.0,
                                 [&](double x, double y) {
                                     const double a = 2 * kPi * x / na;
                                     return Vec3{R * std::cos(a), R * std::sin(a), y * H / 48};
                                 },
                                 [&](double x, double y) { return Vec2{2 * kPi * R * x / na / 100, y * H / 48 / 100}; });
    // Weld the seam column (x = 0 and x = na coincide): positions shared, corner uvs differ.
    const TriMesh m = weldVertices(raw, 1e-9).mesh;
    const MeshCheck mc = checkMesh(m);
    CHECK(mc.components.size() == 1 && mc.components[0].uvCharts == 1 && mc.components[0].euler == 0);
    MeshKnitParams p = tiledParams(expandedD1());
    PartKnit k;
    std::string err;
    CHECK(knitMesh(m, p, knit_test::testNumerics(), k, &err));
    const UvChart& ch = k.atlas->charts()[0];
    CHECK(ch.upFromWorld && std::fabs(ch.upDir.x) < 1e-9);  // wales along the axis
    if (k.verdict != Status::Pass) throw std::runtime_error("cylinder: " + failing(k));
    // The seam gap: the rows end on both sides of the cut, nothing crosses it.
    CHECK(k.t.charts[0].runs > 0 && k.t.charts[0].boundary + k.t.charts[0].runEnd > 0);
}

// Two charts on one plane, cut by a uv seam at x = 40: with the footprint pre-filter the
// yarns of the two charts clear each other (gate a PASS, clearance across the seam
// checked); without it (negative control) stitches sit on the seam from both sides and
// collide: clearance FAIL. (A zero inset alone keeps a gap: each chart's lattice starts at
// its own boundary, see TileParams::footprintCheck.)
// An 80 x 50 mm plane split into two uv charts along x = 40 (the seam fixture).
TriMesh seamPlane() {
    const TriMesh raw = gridMesh(80, 50, 1.0, [](double x, double y) { return Vec3{x, y, 0}; },
                                 [](double x, double y) { return Vec2{(x <= 40 ? x : x + 1000) / 100, y / 100}; });
    // gridMesh evaluates the uv per corner from the corner's grid point: the seam column
    // x = 40 gets the left side's uv in both cells. Rebuild the right cells' corners on x = 40.
    TriMesh m = raw;
    for (size_t t = 0; t < m.tri.size(); ++t) {
        double cx = 0;
        for (int c = 0; c < 3; ++c) cx += m.pos[m.tri[t][c]].x / 3;
        if (cx > 40)
            for (int c = 0; c < 3; ++c) {
                const Vec3& q = m.pos[m.tri[t][c]];
                m.uv[t][c] = {(q.x + 1000) / 100, q.y / 100};
            }
    }
    return m;
}

// A device failure while classifying a chart ends that part with a message, never the run.
class FailingClassifier : public knit_test::ForwardingNumerics {
public:
    using ForwardingNumerics::ForwardingNumerics;
    void classifyLattice(const UvAtlas&, const LatticeChartQuery&, const LatticeSolids&, std::vector<LatticeNodeClass>&) override {
        throw PartFailure("device lost");
    }
};

TEST(tiled_classify_failure_ends_the_part) {
    FailingClassifier numerics(knit_test::cpuNumerics());
    PartKnit k;
    std::string err;
    CHECK(!knitMesh(seamPlane(), tiledParams(expandedD1()), numerics, k, &err));
    CHECK(err.find("device lost") != std::string::npos);
}

TEST(uv_tiling_seam_clearance) {
    const TriMesh m = seamPlane();
    const MeshCheck mc = checkMesh(m);
    CHECK(mc.components.size() == 1 && mc.components[0].uvCharts == 2);
    MeshKnitParams p = tiledParams(expandedD1());
    p.settle = false;
    PartKnit k;
    std::string err;
    CHECK(knitMesh(m, p, knit_test::testNumerics(), k, &err));
    CHECK(k.t.charts.size() == 2 && k.t.charts[0].runs > 0 && k.t.charts[1].runs > 0);
    if (k.after.clearance.status != Status::Pass) throw std::runtime_error("seam with inset: " + failing(k));
    CHECK(k.after.verdict == Status::Pass);
    p.tile.footprintCheck = false;
    PartKnit k0;
    CHECK(knitMesh(m, p, knit_test::testNumerics(), k0, &err));
    CHECK(k0.after.clearance.status == Status::Fail);
    for (const PartKnit* x : {&k, &k0})
        std::printf("    %s: tile %.0f, validate %.0f, settle %.0f, boundary %.0f, gate d %.0f ms (%s)\n",
                    x == &k ? "with the inset" : "without", x->msTile, x->msValidate, x->msSettle, x->msBoundary, x->msGateD,
                    knit_test::useGpu() ? "gpu" : "cpu");
}

// ---- coverage (seamless plan M0) ----------------------------------------------------------------

TubeMesh coverageTube(const CurveSet& cs, const PatternParams& pp) {
    SurfaceParams sp;
    sp.compression = false;
    const SurfaceModel sm = buildSurface(cs, pp, sp);
    TubeParams tp;
    tp.ringSpacing = 0.1;
    tp.sides = 12;
    return buildTubeMesh(cs, sm, tp);
}

// The knit without the yarns of lattice rows `rows`: curves, binding states and program.
struct DroppedRows {
    CurveSet cs;
    std::vector<std::vector<ProjState>> states;
    StitchProgram prog;
    int stitches = 0;  // dropped
};
DroppedRows dropRows(const PartKnit& k, const std::set<int>& rows) {
    DroppedRows d;
    std::set<int> turns;
    for (const StitchRecord& s : k.t.prog.stitches)
        if (rows.count(s.lattice.r)) turns.insert(s.turn), ++d.stitches;
    d.cs = k.t.curves;
    d.cs.yarns.clear();
    for (size_t y = 0; y < k.t.curves.yarns.size(); ++y)
        if (!turns.count(static_cast<int>(y))) {
            d.cs.yarns.push_back(k.t.curves.yarns[y]);
            d.states.push_back(k.t.bind.state[y]);
        }
    d.prog = k.t.prog;
    d.prog.stitches.clear();
    for (const StitchRecord& s : k.t.prog.stitches)
        if (!rows.count(s.lattice.r)) d.prog.stitches.push_back(s);
    return d;
}

PartKnit coveragePlaneKnit(const TriMesh& m) {
    MeshKnitParams p = tiledParams(expandedD1());
    p.settle = false;
    PartKnit k;
    std::string err;
    if (!knitMesh(m, p, knit_test::cpuNumerics(), k, &err)) throw std::runtime_error(err);
    CHECK(k.t.curves.yarns.size() == k.t.prog.turns.size());
    return k;
}

TriMesh coveragePlane() {
    return gridMesh(70, 50, 1.0, [](double x, double y) { return Vec3{x, y, 0}; }, [](double x, double y) { return Vec2{x / 100, y / 100}; });
}

// On a plane, the normal-segment coverage equals the front projection's open area (study.h).
TEST(coverage_flat_matches_mesh_open_area) {
    const TriMesh m = coveragePlane();
    const PartKnit k = coveragePlaneKnit(m);
    const TubeMesh tube = coverageTube(k.t.curves, k.t.pp);
    Vec3 lo{1e300, 1e300, 0}, hi{-1e300, -1e300, 0};
    for (const LayoutNode& n : k.t.lay.nodes) {
        lo = {std::min(lo.x, n.p.x), std::min(lo.y, n.p.y), 0};
        hi = {std::max(hi.x, n.p.x), std::max(hi.y, n.p.y), 0};
    }
    const double w = k.t.pp.waleSpacing, h = k.t.pp.courseSpacing;
    const InteriorWindow win{lo.x + 1.5 * w, lo.y + 1.5 * h, hi.x - 1.5 * w, hi.y - 1.5 * h};
    TriMesh wm;
    wm.pos = {{win.x0, win.y0, 0}, {win.x1, win.y0, 0}, {win.x1, win.y1, 0}, {win.x0, win.y1, 0}};
    wm.tri = {{0, 1, 2}, {0, 2, 3}};
    wm.computeNormals();
    CoverageInput in;
    in.mesh = &wm;
    in.own.tube = &tube;  // no feet: every hit counts
    CoverageParams cp;
    cp.spacing = 0.05;
    const CoverageReport r = measureCoverage(in, cp, ExecPolicy{true});
    const double open = meshOpenArea(tube, win, 800);
    std::printf("    window %.1f x %.1f mm: coverage uncovered %.4f (%lld samples), meshOpenArea %.4f\n", win.x1 - win.x0, win.y1 - win.y0,
                r.uncoveredFrac, r.samples, open);
    CHECK(open > 0.05 && open < 0.95);
    CHECK_NEAR(r.uncoveredFrac, open, 0.01);
}

// Dropping stitches never changes the sample set (fixed domain).
TEST(coverage_domain_fixed) {
    const TriMesh m = coveragePlane();
    const PartKnit k = coveragePlaneKnit(m);
    const DroppedRows d = dropRows(k, {3, 4});
    const TubeMesh full = coverageTube(k.t.curves, k.t.pp), part = coverageTube(d.cs, k.t.pp);
    CoverageParams cp;
    CoverageInput a, b;
    a.mesh = b.mesh = &m;
    a.own.tube = &full;
    a.own.vertexFoot = tubeVertexFeet(full, k.t.curves, k.t.bind.state);
    b.own.tube = &part;
    b.own.vertexFoot = tubeVertexFeet(part, d.cs, d.states);
    const CoverageReport ra = measureCoverage(a, cp, ExecPolicy{true}), rb = measureCoverage(b, cp, ExecPolicy{true});
    CHECK(ra.samples == rb.samples && ra.samples > 0);
    CHECK_NEAR(ra.spacing, rb.spacing, 0);
    CHECK(rb.uncovered > ra.uncovered);
}

// What the overlap mode's gap patches read from a coverage report, bit for bit (pinned at 05960b1): the counts
// and every hole region's area, width, box, centre and its triangle, on the plane with two dropped rows.
TEST(coverage_regions_pinned) {
    const TriMesh m = coveragePlane();
    const PartKnit k = coveragePlaneKnit(m);
    const DroppedRows d = dropRows(k, {3, 4});
    const TubeMesh part = coverageTube(d.cs, k.t.pp);
    CoverageInput in;
    in.mesh = &m;
    in.own.tube = &part;
    in.own.vertexFoot = tubeVertexFeet(part, d.cs, d.states);
    const CoverageReport r = measureCoverage(in, CoverageParams{}, ExecPolicy{true});
    std::string text = std::to_string(r.samples) + " " + std::to_string(r.uncovered) + " " + std::to_string(r.unownedHits) + " " +
                       std::to_string(r.holes.count) + " ";
    for (double x : {r.spacing, r.area, r.sampleArea, r.uncoveredFrac, r.holes.maxArea, r.holes.maxWidth, r.holes.totalArea})
        appendExact(text, x);
    for (const CoverageRegion& g : r.holes.items) {
        text += "\n" + std::to_string(g.centreTri) + " ";
        for (double x : {g.area, g.width, g.lo.x, g.lo.y, g.lo.z, g.hi.x, g.hi.y, g.hi.z, g.centre.x, g.centre.y, g.centre.z, g.centreBary[0],
                         g.centreBary[1], g.centreBary[2]})
            appendExact(text, x);
    }
    const uint64_t digest = fnv1a64(text);
    std::printf("    %lld samples, %d hole regions, digest %016llx\n", r.samples, r.holes.count, static_cast<unsigned long long>(digest));
    CHECK(r.holes.count > 0);
    CHECK(digest == 0xac1a1b108d1e9677ull);
}

// One missing course: the uncovered area grows by less than the strip (n w h): this fabric is open
// lace (about half of its interior is uncovered anyway) and the neighbour rows' loops reach into the
// strip; at least 30 % of the strip turns uncovered.
TEST(coverage_missing_course) {
    const TriMesh m = coveragePlane();
    const PartKnit k = coveragePlaneKnit(m);
    const int row = 3;
    const DroppedRows d = dropRows(k, {row});
    CHECK(d.stitches > 4);
    const TubeMesh full = coverageTube(k.t.curves, k.t.pp), part = coverageTube(d.cs, k.t.pp);
    CoverageParams cp;
    cp.spacing = 0.15;
    CoverageInput a, b;
    a.mesh = b.mesh = &m;
    a.own.tube = &full;
    a.own.vertexFoot = tubeVertexFeet(full, k.t.curves, k.t.bind.state);
    b.own.tube = &part;
    b.own.vertexFoot = tubeVertexFeet(part, d.cs, d.states);
    const CoverageReport ra = measureCoverage(a, cp, ExecPolicy{true}), rb = measureCoverage(b, cp, ExecPolicy{true});
    const double w = k.t.pp.waleSpacing, h = k.t.pp.courseSpacing, strip = d.stitches * w * h;
    const double holeGrowth = (rb.uncovered - ra.uncovered) * ra.sampleArea;
    std::printf("    full: %s\n    dropped row %d (%d stitches, %.0f mm2): %s\n    uncovered +%.1f mm2\n", ra.summary().c_str(), row,
                d.stitches, strip, rb.summary().c_str(), holeGrowth);
    CHECK(holeGrowth >= 0.3 * strip);
}

// Thin double sheet: B lies 1.5 mm under A, facing away (a thin ear). Courses missing on B are
// uncovered, although A's yarn crosses B's normal segments: that yarn is not eligible (opposite
// normal, no surface path). Measured on B's dropped strip only (a fixed region: the two sheets'
// mirrored patterns are not independent).
TEST(coverage_thin_double_sheet) {
    const TriMesh A = coveragePlane();
    const PartKnit k = coveragePlaneKnit(A);
    const std::set<int> rows{3, 4};
    const DroppedRows d = dropRows(k, rows);
    auto mirror = [](const Vec3& q) { return Vec3{q.x, 50 - q.y, -1.5 - q.z}; };  // mirrored: B faces -z
    TriMesh M = A;
    const int nv = static_cast<int>(A.pos.size()), ntA = static_cast<int>(A.tri.size());
    for (const Vec3& q : A.pos) M.pos.push_back(mirror(q));
    for (const auto& t : A.tri) M.tri.push_back({t[0] + nv, t[1] + nv, t[2] + nv});
    M.uv.clear();
    M.computeNormals();
    CHECK(M.faceNormal(0).z > 0.99 && M.faceNormal(ntA).z < -0.99);
    const TubeMesh tA = coverageTube(k.t.curves, k.t.pp);
    TubeMesh tB = coverageTube(d.cs, k.t.pp);
    for (Vec3& q : tB.pos) q = mirror(q);
    TubeMesh both = tA;
    for (const Vec3& q : tB.pos) both.pos.push_back(q);
    for (uint32_t i : tB.tri) both.tri.push_back(i + static_cast<uint32_t>(tA.pos.size()));
    std::vector<int> feet = tubeVertexFeet(tA, k.t.curves, k.t.bind.state);
    for (int f : tubeVertexFeet(tB, d.cs, d.states)) feet.push_back(f >= 0 ? f + ntA : -1);
    // The region: sheet B (z = -1.5; A lies at z = 0) inside the dropped rows' cells (A's uv' = A's y,
    // B's y = 50 - A's; 1 mm margins).
    const TemplateExtent extent = templateExtent(k.t.pp);
    const double h = k.t.pp.courseSpacing, ymid = 0.5 * (extent.Ymin + extent.Ymax);
    double yr0 = 1e300, yr1 = -1e300, x0 = 1e300, x1 = -1e300;
    for (const StitchRecord& s : k.t.prog.stitches)
        if (rows.count(s.lattice.r)) {
            const Vec3& p = k.t.lay.nodes[s.id].p;
            yr0 = std::min(yr0, p.y), yr1 = std::max(yr1, p.y), x0 = std::min(x0, p.x), x1 = std::max(x1, p.x);
        }
    const double ya = yr0 + (ymid - 0.5) * h + 1, yb = yr1 + (ymid + 0.5) * h - 1;  // in A's coordinates
    auto run = [&](bool ownership) {
        CoverageInput in;
        in.mesh = &M;
        in.own.tube = &both;
        if (ownership) in.own.vertexFoot = feet;
        in.regionLo = {x0, 50 - yb, -3};
        in.regionHi = {x1, 50 - ya, -0.75};
        CoverageParams cp;
        cp.spacing = 0.1;
        return measureCoverage(in, cp, ExecPolicy{true});
    };
    const CoverageReport r = run(true), naive = run(false);
    const double n = double(r.samples);
    std::printf("    B's dropped strip (%lld samples): %s\n    hit only by A's (ineligible) yarn %.3f; without ownership: uncovered %.3f\n",
                r.samples, r.summary().c_str(), r.unownedHits / n, naive.uncoveredFrac);
    CHECK(r.samples > 1000 && r.samples == naive.samples);
    CHECK(r.uncoveredFrac > 0.6);                          // only the neighbour rows' loops reach in
    CHECK(r.unownedHits / n > 0.25);                       // A's yarn crosses B's normal segments, rejected
    CHECK(naive.uncoveredFrac < r.uncoveredFrac - 0.2);   // without ownership A's yarn would count
}

// A uv seam inside a plane: the per-component tiler leaves a strip open along it. The footprint
// pre-filter keeps every stitch node at least half a wale from the seam (x = 40), and the coverage
// finds an open region running along the seam.
TEST(coverage_chart_seam) {
    const TriMesh m = seamPlane();
    MeshKnitParams p = tiledParams(expandedD1());
    p.settle = false;
    PartKnit k;
    std::string err;
    CHECK(knitMesh(m, p, knit_test::cpuNumerics(), k, &err));
    const double w = k.t.pp.waleSpacing;
    double nearest = 1e300;
    for (const LayoutNode& n : k.t.lay.nodes)
        if (n.p.y > 8 && n.p.y < 42) nearest = std::min(nearest, std::fabs(n.p.x - 40));
    const TubeMesh tube = coverageTube(k.t.curves, k.t.pp);
    CoverageInput in;
    in.mesh = &m;
    in.own.tube = &tube;
    in.own.vertexFoot = tubeVertexFeet(tube, k.t.curves, k.t.bind.state);
    CoverageParams cp;
    cp.spacing = 0.15;
    const CoverageReport r = measureCoverage(in, cp, ExecPolicy{true});
    double alongSeam = 0;  // the longest open region across x = 40, along y
    for (const CoverageRegion& g : r.holes.items)
        if (g.lo.x < 40 && g.hi.x > 40) alongSeam = std::max(alongSeam, g.hi.y - g.lo.y);
    std::printf("    %s\n    nearest node %.2f mm from the seam (w %.2f mm); longest open region across it %.1f mm along y\n",
                r.summary().c_str(), nearest, w, alongSeam);
    CHECK(nearest > 0.5 * w);
    CHECK(alongSeam > 25);
}

// knit.coverage = overlap: each chart of the seam plane knitted alone (TileParams::onlyChart) with
// two stitches / rows of overhang walked beyond it; their tubes together cover the seam strip that
// the per-component tiling leaves open (stitches of chart 0 reach past the seam onto chart 1).
TEST(overlap_seam_covered) {
    const TriMesh m = seamPlane();
    CoverageParams cp;
    cp.spacing = 0.15;
    auto stripCoverage = [&](const std::vector<const PartKnit*>& ks) {
        TubeMesh all;
        std::vector<int> feet;
        for (const PartKnit* k : ks) {
            const TubeMesh tube = coverageTube(k->t.curves, k->t.pp);
            const std::vector<int> f = tubeVertexFeet(tube, k->t.curves, k->t.bind.state);
            const uint32_t base = static_cast<uint32_t>(all.pos.size());
            all.pos.insert(all.pos.end(), tube.pos.begin(), tube.pos.end());
            for (uint32_t v : tube.tri) all.tri.push_back(base + v);
            feet.insert(feet.end(), f.begin(), f.end());
        }
        CoverageInput in;
        in.mesh = &m;
        in.own.tube = &all;
        in.own.vertexFoot = feet;
        in.regionLo = {37, 8, -1e300};  // the seam strip |x - 40| < 3, 8 < y < 42
        in.regionHi = {43, 42, 1e300};
        return measureCoverage(in, cp, ExecPolicy{true});
    };
    MeshKnitParams p = tiledParams(expandedD1());
    p.settle = false;
    PartKnit whole;
    std::string err;
    CHECK(knitMesh(m, p, knit_test::cpuNumerics(), whole, &err));
    const CoverageReport before = stripCoverage({&whole});
    p.tile.footprintCheck = false;
    p.tile.overhang = 2;
    PartKnit k0, k1;
    p.tile.onlyChart = 0;
    CHECK(knitMesh(m, p, knit_test::cpuNumerics(), k0, &err));
    p.tile.onlyChart = 1;
    CHECK(knitMesh(m, p, knit_test::cpuNumerics(), k1, &err));
    const CoverageReport after = stripCoverage({&k0, &k1});
    std::printf("    seam strip: per component %.1f %% uncovered (max region width %.2f mm); overlap %.1f %% (max width %.2f mm)\n",
                100 * before.uncoveredFrac, before.holes.maxWidth, 100 * after.uncoveredFrac, after.holes.maxWidth);
    int walked = 0, skipped = 0;
    for (const ChartTiling& ct : k0.t.charts) walked += ct.overhang, skipped += ct.skipped;
    // Chart 0 is either side of the seam: its nodes cross the seam line by more than a stitch.
    double xMin0 = 1e300, xMax0 = -1e300, xSum = 0;
    for (const LayoutNode& n : k0.t.lay.nodes) xMin0 = std::min(xMin0, n.p.x), xMax0 = std::max(xMax0, n.p.x), xSum += n.p.x;
    const bool left = xSum / k0.t.lay.nodes.size() < 40;
    const double past = left ? xMax0 - 40 : 40 - xMin0;
    std::printf("    chart 0 part (%s of the seam): %zu stitches (%d walked), %d chart(s) skipped, nodes reach %.2f mm past the seam\n",
                left ? "left" : "right", k0.t.prog.stitches.size(), walked, skipped, past);
    CHECK(walked > 0 && skipped == 1);
    CHECK(past > k0.t.pp.waleSpacing);
    CHECK(k0.t.lattice.status == Status::Pass && k1.t.lattice.status == Status::Pass);
    // The strip is covered as the open fabric is (no layout gap left), far better than before.
    CHECK(after.uncoveredFrac < 0.5 * before.uncoveredFrac);
    CHECK(after.holes.maxWidth < k0.t.pp.courseSpacing);
}

// ---- two parallel yarns ---------------------------------------------------------------------------

YarnPath contractLine(double y, int id) {
    YarnPath p;
    p.id = id;
    p.course = id;
    for (int i = -1; i <= 11; ++i) p.curve.ctrl.push_back({double(i), y, 0});
    p.stitches.push_back({0, id, 0.0, p.curve.uMax()});
    p.finalize();
    return p;
}

CurveSet twoLines(double d) {
    CurveSet cs;
    cs.yarns = {contractLine(0, 0), contractLine(d, 1)};
    cs.wales = 1;
    cs.courses = 2;
    fillFlatTopology(cs);
    return cs;
}

// Two parallel yarns drawn at D 1, flattened by at most 0.1 with no widening: at d = 0.92 D the tubes
// flatten by 0.08; at d = 0.88 D they need more than the bound (clamped).
TEST(compress_two_yarns) {
    PatternParams pp;
    pp.yarnDiameter = 1.0;
    SurfaceParams sp;
    sp.compression = true;
    sp.maxExpansion = 0;
    sp.maxCompression = 0.1;
    {
        const SurfaceModel sm = buildSurface(twoLines(0.92), pp, sp);
        std::printf("    d 0.92: clamped %d, min b/r %.4f\n", sm.compression.clamped, sm.compression.minAxisRatio);
        CHECK(sm.compression.clamped == 0);
        CHECK_NEAR(1 - sm.compression.minAxisRatio, 0.08, 0.005);
    }
    {
        const SurfaceModel sm = buildSurface(twoLines(0.88), pp, sp);
        std::printf("    d 0.88: clamped %d, min b/r %.4f\n", sm.compression.clamped, sm.compression.minAxisRatio);
        CHECK(sm.compression.clamped > 0);
    }
}

// Gate (d) on flattened yarn counts the rings that needed more than SurfaceParams::maxCompression,
// and those rings alone are no verdict.
TEST(compress_gate_d_clamped_reported) {
    PartKnit k;
    k.t.curves = twoLines(0.55);
    k.t.pp.yarnDiameter = 1.0;
    k.t.pp.clearance = 0.1;
    k.after.verdict = Status::Pass;
    MeshKnitParams p;
    p.surface.compression = true;
    p.surfaceGate = true;
    runSurfaceGate(k, p, knit_test::cpuNumerics());
    std::printf("    d 0.55: gate d %s, clamped %d\n", toString(k.gateDStatus), k.compressionClamped);
    CHECK(k.compressionClamped > 0);
    CHECK(k.gateDStatus != Status::Fail || k.gateD.back().fold == Status::Fail || k.gateD.back().collisions == Status::Fail);
}

// strand.max_flattening is the flattening limit of the final surface: a knit input only while the
// yarn is flattened and the limit differs from 0.4 (stored keys stay valid), carried by finalSurface
// into gate (d), and refused outside 0.05..0.6. Two yarns at d = 0.8 D need 0.2: within 0.4, past 0.15.
TEST(max_flattening_setting) {
    SceneSettings round;
    SceneSettings limited = round;
    limited.maxFlattening = 0.15;
    CHECK(stageKey(round, Stage::Knit, 1) == stageKey(limited, Stage::Knit, 1));  // a round yarn: no input
    SceneSettings flat = round;
    flat.flatten = true;
    SceneSettings flatLimited = flat;
    flatLimited.maxFlattening = 0.15;
    CHECK(stageKey(flat, Stage::Knit, 1) != stageKey(flatLimited, Stage::Knit, 1));
    SceneSettings flatDefault = flat;
    flatDefault.maxFlattening = 0.4;
    CHECK(stageKey(flat, Stage::Knit, 1) == stageKey(flatDefault, Stage::Knit, 1));
    const FinalSurfaceSpec fs = finalSurface(flatLimited, 1.0);
    CHECK(fs.maxCompression == 0.15);
    std::string why;
    CHECK(fs.valid(&why));
    FinalSurfaceSpec tooMuch = fs;
    tooMuch.maxCompression = 0.7;
    CHECK(!tooMuch.valid(&why));
    int clamped[2] = {0, 0};
    for (int k = 0; k < 2; ++k) {
        PartKnit part;
        part.t.curves = twoLines(0.8);
        part.t.pp.yarnDiameter = 1.0;
        part.t.pp.clearance = 0.1;
        part.after.verdict = Status::Pass;
        MeshKnitParams p;
        p.surface.compression = true;
        p.surface.maxCompression = k == 0 ? 0.4 : 0.15;
        p.surfaceGate = true;
        runSurfaceGate(part, p, knit_test::cpuNumerics());
        clamped[k] = part.compressionClamped;
    }
    std::printf("    d 0.8: clamped %d at 0.4, %d at 0.15\n", clamped[0], clamped[1]);
    CHECK(clamped[0] == 0);
    CHECK(clamped[1] > 0);
}

// A bend of radius 0.475 mm lies between the yarn radii 0.45 and 0.5: the fold check FAILs at 0.5 on
// the CPU and the GPU alike (validate's curvature pass and its fold line), and passes at 0.45.
TEST(fold_radius_cpu_gpu_agree) {
    const double R = 0.475, n = 24, Rc = R * 3 / (2 + std::cos(2 * kPi / n));  // control circle of a B-spline of radius ~R
    YarnPath y;
    for (int i = 0; i <= 14; ++i) y.curve.ctrl.push_back({Rc * std::cos(2 * kPi * i / n), Rc * std::sin(2 * kPi * i / n), 0});
    y.stitches.push_back({0, 0, 0.0, y.curve.uMax()});
    y.finalize();
    CurveSet cs;
    cs.yarns = {y};
    cs.wales = cs.courses = 1;
    cs.kind = CurveKind::Surface;  // the GPU validates surface curve sets (no global flat signature)
    PatternParams vp;
    vp.yarnDiameter = 1.0;
    vp.clearance = 0;
    PatternParams core = vp;
    core.yarnDiameter = 0.9;  // the radius 0.45
    auto fold = [](const ValidationReport& r) {
        for (const CheckLine& l : r.lines)
            if (l.name.rfind("no local fold", 0) == 0) return l.status;
        return Status::Unresolved;
    };
    const ValidationReport cpu = validate(cs, vp, ValidateOptions{.generatorChecks = false});
    const ValidationReport cpuCore = validate(cs, core, ValidateOptions{.generatorChecks = false});
    std::printf("    rho_lo %.4f mm: CPU fold at r 0.5 %s, at r 0.45 %s\n", cpu.curvature[0].rhoLo, toString(fold(cpu)), toString(fold(cpuCore)));
    CHECK(cpu.curvature[0].rhoLo > 0.45 && cpu.curvature[0].rhoLo < 0.5);
    CHECK(fold(cpu) == Status::Fail && cpu.curvature[0].noFold == Status::Fail);
    CHECK(fold(cpuCore) == Status::Pass);
    if (!knit_test::useGpu()) SKIP("no GPU: the CPU half ran");
    std::string err;
    std::unique_ptr<KnitNumerics> g = gpu::makeKnitNumerics(&err);
    CHECK(g != nullptr);
    const ValidationReport gpuR = g->validate(cs, vp, ValidateOptions{.generatorChecks = false});
    std::printf("    GPU fold at r 0.5 %s (noFold %s)\n", toString(fold(gpuR)), toString(gpuR.curvature[0].noFold));
    CHECK(fold(gpuR) == fold(cpu) && gpuR.curvature[0].noFold == cpu.curvature[0].noFold);
}

// The relax call's hard thresholds follow the validation parameters.
TEST(relax_thresholds_follow_validation) {
    PatternParams pp;
    pp.yarnDiameter = 1.2;
    pp.clearance = 0.1;
    const HardThresholds hard = hardThresholds(pp, 0.03, 0.03);
    CHECK_NEAR(hard.contact, 1.2 + 0.1 + 0.03, 1e-12);
    CHECK_NEAR(hard.kappa, 1.0 / (0.6 + 0.05) / 1.03, 1e-12);
}

// capSide's head-side X is the template's own; a purl stitch mirrors its depth.
TEST(capside_fixed_bit_identical) {
    PatternParams pp;
    CHECK(templateHeadSideX(pp) == pp.headHalfWidth + 0.02);
    CHECK(templateHeadSideX(pp) == std::fabs(stitchTemplate(pp)[10].X));
    const std::vector<TemplateControl> t = stitchTemplate(pp);
    StitchFaces purl;
    purl.own = purl.child = purl.left = StitchFace::Purl;
    const std::vector<TemplateControl> tp = stitchTemplate(pp, purl);
    for (size_t j = 0; j < t.size(); ++j) CHECK(tp[j].Z == -t[j].Z && tp[j].X == t[j].X);
}

// Forwards to another KnitNumerics and counts the validate() calls (gate (a) has started).
struct ValidateProbe : knit_test::ForwardingNumerics {
    explicit ValidateProbe(KnitNumerics& n) : ForwardingNumerics(n), inner(n) {}
    ValidationReport validate(const CurveSet& cs, const PatternParams& pp, const ValidateOptions& opt) override {
        ++calls;
        return inner.validate(cs, pp, opt);
    }
    KnitNumerics& inner;
    std::atomic<int> calls{0};
};

// The GPU numerics' accumulated curvature time [ms] (knitNumericsStats: "(upload U, curvature C,").
// It grows only when a curvature launch completes, so it shows which phase a cancel landed in.
double gpuCurvatureMs(KnitNumerics& n) {
    const std::string s = gpu::knitNumericsStats(n);
    const std::string key = ", curvature ";
    const size_t at = s.find("(upload ");
    const size_t k = at == std::string::npos ? at : s.find(key, at);
    if (k == std::string::npos) return -1;
    char* end = nullptr;
    const double c = std::strtod(s.c_str() + k + key.size(), &end);
    return end == s.c_str() + k + key.size() ? -1 : c;
}

// Cancelling during the slowest validation we know: the seam fixture without its inset, where
// gate (a) takes minutes (on the GPU: about 3 s of curvature, then the local signature). It runs
// through the real builder, and each cancel is requested only after gate (a) has started (a
// probe on the numerics): once 0.5 s in (inside the curvature launch) and once 6 s in (after it,
// in the signature). Each must stop within the 1 s target, and the completed result before it
// stays the latest complete one. Afterwards the same
// worker (its numerics and abort flag) runs a normal job: it completes with the validation of the
// first job, and nothing of the cancelled jobs was published.
TEST(tiled_validation_cancel) {
    const TriMesh m = seamPlane();
    std::unique_ptr<ValidateProbe> probe;
    std::atomic<bool> probeReady{false};
    std::mutex digestLock;
    std::map<int, uint64_t> digest;  // revision -> gate (a) report digest of a completed knit
    SceneBuilder b([&](const JobRequest&, const JobControl& j, std::string* err) -> std::shared_ptr<BuiltScene> {
        MeshKnitParams p = tiledParams(expandedD1());
        p.settle = false;
        const bool slow = j.revision == 2 || j.revision == 3;  // the cancelled ones; 1 and 4: the inset (fast, PASS)
        p.tile.footprintCheck = !slow;
        KnitNumerics* numerics = &knit_test::testNumerics();
        if (slow) {
            probe = std::make_unique<ValidateProbe>(*numerics);
            numerics = probe.get();
            probeReady = true;
        }
        PartKnit k;
        if (!knitMesh(m, p, *numerics, k, err)) return nullptr;
        {
            std::lock_guard<std::mutex> lk(digestLock);
            digest[j.revision] = fnv1a64(validationSnapshot(k.after));
        }
        auto r = std::make_shared<BuiltScene>();
        r->revision = j.revision;
        return r;
    });
    const KnitRequest job{std::make_shared<const SceneSettings>()};
    CHECK(b.startKnit(job) == 1);
    while (b.snapshot().busy) std::this_thread::sleep_for(std::chrono::milliseconds(2));
    CHECK(b.snapshot().latestComplete && b.snapshot().latestComplete->revision == 1);
    const int published = b.snapshot().published;
    const bool gpu = knit_test::useGpu();
    for (const bool afterTheLaunch : {false, true}) {
        probeReady = false;
        b.startKnit(job);
        const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(120);
        while (!(probeReady && probe->calls > 0) && std::chrono::steady_clock::now() < deadline)
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        CHECK(probeReady && probe->calls > 0);  // gate (a) is running
        const double curv0 = gpu ? gpuCurvatureMs(probe->inner) : 0;
        // Two phases of gate (a): soon after it starts, inside the curvature launch; or just after
        // that launch ended, which a loaded GPU reaches late, so the test waits for the event itself.
        if (afterTheLaunch && gpu) {
            while (gpuCurvatureMs(probe->inner) == curv0 && std::chrono::steady_clock::now() < deadline)
                std::this_thread::sleep_for(std::chrono::milliseconds(5));
        } else {
            std::this_thread::sleep_for(std::chrono::milliseconds(afterTheLaunch ? 6000 : 500));
        }
        const auto t0 = std::chrono::steady_clock::now();
        b.cancel();
        while (b.snapshot().busy) std::this_thread::sleep_for(std::chrono::milliseconds(1));
        const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        const double curv1 = gpu ? gpuCurvatureMs(probe->inner) : 0;
        std::printf("    cancel %s the curvature launch (%s%s): stopped after %.0f ms\n", afterTheLaunch ? "after" : "during",
                    gpu ? "gpu" : "cpu", !gpu ? "" : curv1 > curv0 ? ", after the curvature launch" : ", inside the curvature launch", ms);
        const JobSnapshot end = b.snapshot();
        CHECK(end.lastOutcome == JobOutcome::Cancelled);
        CHECK(end.latestComplete && end.latestComplete->revision == 1);  // the previous result stays
        CHECK(ms < 1000);
        if (gpu) CHECK(afterTheLaunch ? curv1 > curv0 : curv1 == curv0);  // the phase the cancel landed in
    }
    // Nothing of the cancelled jobs was published, complete or partial.
    const JobSnapshot after = b.snapshot();
    CHECK(after.published == published && after.latest && after.latest->revision == 1);
    CHECK(digest.size() == 1);  // only job 1 completed a knit
    // The same worker again, a normal job: it completes, with job 1's validation (a raised abort
    // flag or a stale device state would end kernels early and change the report).
    CHECK(b.startKnit(job) == 4);
    while (b.snapshot().busy) std::this_thread::sleep_for(std::chrono::milliseconds(2));
    const JobSnapshot last = b.snapshot();
    CHECK(last.lastOutcome == JobOutcome::Done && last.latestComplete && last.latestComplete->revision == 4);
    CHECK(digest.count(4) && digest[4] == digest[1]);
}

// A chart whose uvs are mirrored (u = -x): flagged, frames right-handed, no inverted cell.
TEST(uv_tiling_mirrored_chart) {
    const TriMesh m = gridMesh(70, 50, 1.0, [](double x, double y) { return Vec3{x, y, 0}; },
                               [](double x, double y) { return Vec2{-x / 100, y / 100}; });
    MeshKnitParams p = tiledParams(expandedD1());
    p.settle = false;
    PartKnit k;
    std::string err;
    CHECK(knitMesh(m, p, knit_test::testNumerics(), k, &err));
    CHECK(k.atlas->charts()[0].mirrored);
    CHECK(k.t.layout.inverted == 0 && k.t.layout.status == Status::Pass);
    CHECK(k.after.verdict == Status::Pass);
}

// Degenerate, folded and overlapping uv triangles are unusable and counted; a chart with
// more than 1 % folded + overlap area is UNRESOLVED and not tiled.
TEST(uv_chart_degenerate_folded_overlap) {
    // Fold: the uv of the columns x in [30, 34] runs backwards (u decreases), so those
    // triangles are mirrored against the majority AND overlap their neighbours in uv.
    auto u = [](double x) { return x <= 30 ? x : x <= 34 ? 30 - (x - 30) : 26 + (x - 34) * 3; };
    const TriMesh m = gridMesh(70, 40, 1.0, [](double x, double y) { return Vec3{x, y, 0}; },
                               [&](double x, double y) { return Vec2{u(x) / 100, y / 100}; });
    const UvAtlas at(m);
    CHECK(at.charts().size() == 1);
    const UvChart& ch = at.charts()[0];
    CHECK_NEAR(ch.foldedArea, 4 * 40, 1e-9);
    CHECK(ch.overlapArea > 0);
    CHECK(ch.status == Status::Unresolved);
    // A degenerate (zero uv area) strip on an otherwise clean chart: unusable, still tiled.
    auto u2 = [](double x) { return x <= 30 ? x : x <= 31 ? 30.0 : 30 + (x - 31); };
    const TriMesh m2 = gridMesh(70, 40, 1.0, [](double x, double y) { return Vec3{x, y, 0}; },
                                [&](double x, double y) { return Vec2{u2(x) / 100, y / 100}; });
    const UvAtlas at2(m2);
    CHECK_NEAR(at2.charts()[0].degenerateArea, 40, 1e-9);
    CHECK(at2.charts()[0].status == Status::Pass);
    int unusable = 0;
    for (size_t t = 0; t < m2.tri.size(); ++t) unusable += !at2.usable(static_cast<int>(t));
    CHECK(unusable == 80);
}

// A hole splits rows into two runs; the row above the hole is one run whose stitches over
// the gap have no parent. Parents, frames and the rib parity follow (i, r), never the run.
TEST(lattice_runs_split_merge) {
    PatternParams pp = expandedD1();
    pp.pattern = StitchPattern::Rib1x1;
    const TriMesh m = gridMesh(140, 60, 1.0, [](double x, double y) { return Vec3{x, y, 0}; },
                               [](double x, double y) { return Vec2{x / 100, y / 100}; },
                               [](double x, double y) { return !(x > 62 && x < 76 && y > 14 && y < 22); });
    MeshKnitParams p = tiledParams(pp);
    PartKnit k;
    std::string err;
    CHECK(knitMesh(m, p, knit_test::testNumerics(), k, &err));
    const StitchProgram& prog = k.t.prog;
    // Rows with two runs, and a row with one run spanning columns that are missing below.
    std::map<int, int> runsPerRow;
    for (size_t t = 0; t < prog.turns.size(); ++t) ++runsPerRow[prog.stitches[prog.turnStart[t]].lattice.r];
    int split = 0;
    for (const auto& [r, n] : runsPerRow) split += n == 2;
    CHECK(split >= 1);
    int merged = 0;
    for (const StitchRecord& s : prog.stitches) {
        // Parent and pattern parity from the lattice.
        if (s.nParents) {
            const StitchRecord& q = prog.stitches[s.parents[0]];
            CHECK(q.lattice.i == s.lattice.i && q.lattice.r == s.lattice.r - 1);
        }
        CHECK((s.face == StitchFace::Purl) == stitchFaceIsPurl(pp.pattern, s.lattice.i, s.lattice.r));
        // A parentless stitch inside a run (not in the lattice's first row): over the hole.
        if (!s.nParents && s.lattice.r > 0 && k.t.lay.rowPrev[s.id] >= 0 && k.t.lay.rowNext[s.id] >= 0) {
            const StitchRecord& left = prog.stitches[k.t.lay.rowPrev[s.id]];
            merged += left.turn == s.turn;
        }
    }
    CHECK(merged > 0);
    // Frames from the explicit row neighbours: c along +x everywhere, v along +y.
    for (const LayoutNode& nd : k.t.lay.nodes) {
        CHECK(dot(nd.c, Vec3{1, 0, 0}) > 1 - 1e-9);
        CHECK(dot(nd.v, Vec3{0, 1, 0}) > 1 - 1e-6);
    }
    CHECK(k.t.lattice.status == Status::Pass && k.t.lattice.parentless > 0);
    if (k.verdict != Status::Pass) throw std::runtime_error("split/merge: " + failing(k));
}

// The lattice graph check fails a parent in the wrong column, a loop with two children,
// a non-contiguous run; counts parentless stitches and open loops.
TEST(lattice_program_check) {
    StitchProgram prog;
    auto add = [&](int run, int i, int r, int parent) {
        StitchRecord s;
        s.id = static_cast<int>(prog.stitches.size());
        s.turn = run;
        s.lattice = {0, i, r, run};
        if (parent >= 0) {
            s.parents[0] = parent;
            s.nParents = 1;
        }
        prog.stitches.push_back(s);
    };
    for (int i = 0; i < 3; ++i) add(0, i, 0, -1);
    for (int i = 0; i < 3; ++i) add(1, i, 1, i);
    linkChildren(prog.stitches);
    LatticeCheck c = checkLatticeProgram(prog);
    CHECK(c.status == Status::Pass && c.parentless == 3 && c.uncaught == 3);
    StitchProgram bad = prog;
    bad.stitches[4].parents[0] = 0;  // column 1 drawn through column 0
    linkChildren(bad.stitches);
    CHECK(checkLatticeProgram(bad).status == Status::Fail);
    bad = prog;
    bad.stitches[5].lattice.i = 7;  // run not contiguous
    CHECK(checkLatticeProgram(bad).status == Status::Fail);
}

// Welded normals are continuous across a uv seam; an unwelded copy (the seam's vertices
// duplicated) shows the jump (negative control).
TEST(mesh_surface_normal_continuous_across_seam) {
    const double R = 20;
    const int na = 64;
    const TriMesh raw = gridMesh(na, 10, 1.0,
                                 [&](double x, double y) {
                                     const double a = 2 * kPi * x / na;
                                     return Vec3{R * std::cos(a), R * std::sin(a), y};
                                 },
                                 [&](double x, double y) { return Vec2{x / na, y / 100}; });
    const TriMesh welded = weldVertices(raw, 1e-9).mesh;
    const MeshSurface sw(welded), sr(raw);
    // Points just on either side of the seam (angle 0) at mid height.
    const Vec3 a{R, -1e-6, 5}, b{R, 1e-6, 5};
    const double jumpW = norm(sw.project(a).n - sw.project(b).n), jumpR = norm(sr.project(a).n - sr.project(b).n);
    CHECK(jumpW < 1e-6);
    CHECK(jumpR > 1e-2);
}

// B1 is certified and B2 sampled: a yarn deliberately 0.1 mm too close to an open edge
// FAILS B1; a coarse sampling step is refined (never a PASS by skipping).
TEST(boundary_check_certified_vs_sampled) {
    const TriMesh m = gridMesh(40, 30, 1.0, [](double x, double y) { return Vec3{x, y, 0}; },
                               [](double x, double y) { return Vec2{x / 100, y / 100}; });
    const MeshSurface surf(m);
    const double Rs = 0.5;
    auto line = [](double y, double z) {
        CurveSet cs;
        cs.kind = CurveKind::Surface;
        YarnPath yp;
        for (int k = 0; k < 8; ++k) yp.curve.ctrl.push_back({5.0 + 4 * k, y, z});
        yp.finalize();
        cs.yarns.push_back(yp);
        return cs;
    };
    SurfaceBinding none;
    // Parallel to the edge y = 0 at height 0.3: 3D distance sqrt(y^2 + 0.09).
    const BoundaryCheck far = knit_test::checkBoundary(line(1.0, 0.3), none, surf, Rs);
    CHECK(far.b1 == Status::Pass);
    const BoundaryCheck close = knit_test::checkBoundary(line(std::sqrt(0.4 * 0.4 - 0.09), 0.3), none, surf, Rs);  // 0.1 mm too close
    CHECK(close.b1 == Status::Fail);
    CHECK_NEAR(close.minDistance, 0.4, 1e-9);
    // Coarse step: the same verdicts (the intervals are refined with the enclosed arc length).
    CHECK(knit_test::checkBoundary(line(1.0, 0.3), none, surf, Rs, {5.0}).b1 == Status::Pass);
    CHECK(knit_test::checkBoundary(line(std::sqrt(0.4 * 0.4 - 0.09), 0.3), none, surf, Rs, {5.0}).b1 == Status::Fail);
    // A yarn hanging over the edge: B2 sees the overhang (OpenBoundary) from bound states.
    CurveSet over = line(-0.6, 0.3);
    SurfaceBinding bind;
    bind.surface = &surf;
    const double b[3] = {1.0 / 3, 1.0 / 3, 1.0 / 3};
    bind.state.assign(1, std::vector<ProjState>(over.yarns[0].curve.ctrl.size(), surf.stateAt(0, b, surf.stateAt(0, b, {}).foot)));
    bind.result.assign(1, std::vector<ProjResult>(over.yarns[0].curve.ctrl.size(), ProjResult::Owned));
    bind.cp.lambdaW = 100;  // the states start far from the yarn: generous budget
    const BoundaryCheck oh = knit_test::checkBoundary(over, bind, surf, Rs);
    CHECK(oh.b2 == Status::Fail && oh.b2Unresolved[static_cast<int>(ProjResult::OpenBoundary)] > 0);
}

// A U-folded strip: sheet 1 (z = 0) and sheet 2 (z = 1.5) joined by a half-cylinder bend
// of radius 0.75 at x = L; both normals point into the 1.5 mm gap, so the yarn controls
// (depth +A D = 0.8 mm) are closer to the OTHER sheet. Arc coordinate s along the strip.
struct Fold {
    double L = 40, R = 0.75, W = 30;
    double total() const { return 2 * L + kPi * R; }
    Vec3 at(double s, double v) const {
        if (s <= L) return {s, v, 0};
        if (s <= L + kPi * R) {
            const double a = (s - L) / R;
            return {L + R * std::sin(a), v, R - R * std::cos(a)};
        }
        return {L - (s - L - kPi * R), v, 2 * R};
    }
    // From the arc coordinate s (the mesh carries it in u = s / 100): 1: sheet 1, 2: sheet 2,
    // 0: the bend. (A coarse mesh moves the bend's geometry; the arc coordinate does not.)
    int sheetOfS(double s) const { return s < L ? 1 : s > L + kPi * R ? 2 : 0; }
};

double arcOf(const TriMesh& m, int tri, const double b[3]) {
    return 100 * (m.uv[tri][0].x * b[0] + m.uv[tri][1].x * b[1] + m.uv[tri][2].x * b[2]);
}

TriMesh foldMesh(const Fold& f, double cell) {
    const int nx = static_cast<int>(std::ceil(f.total() / cell)), ny = static_cast<int>(std::ceil(f.W / cell));
    const double cs = f.total() / nx, cv = f.W / ny;
    return gridMesh(nx, ny, 1.0, [&](double x, double y) { return f.at(x * cs, y * cv); },
                    [&](double x, double y) { return Vec2{x * cs / 100, y * cv / 100}; });
}

TEST(two_sheets_fold) {
    const Fold f;
    for (double cell : {0.25, 1.0, 4.0}) {
        const TriMesh m = foldMesh(f, cell);
        MeshKnitParams p = tiledParams(expandedD1());
        PartKnit k;
        std::string err;
        CHECK(knitMesh(m, p, knit_test::testNumerics(), k, &err));
        const SurfaceBinding& b = k.t.bind;
        MeshQuery mq(m);
        int owned = 0, wrongOwned = 0, unresolved = 0, unrestrictedWrong = 0;
        for (size_t y = 0; y < b.state.size(); ++y) {
            const int first = k.t.prog.turnStart[y], n = k.t.prog.turns[y].output;
            for (size_t c = 0; c < b.state[y].size(); ++c) {
                const int st = first + std::clamp(static_cast<int>(c) / kControlsPerStitch, 0, n - 1);
                const ProjState& ns = k.t.lay.nodeState[st];
                const int own = f.sheetOfS(arcOf(m, ns.tri, ns.b));
                if (own == 0) continue;  // a node on the bend itself
                if (b.result[y][c] != ProjResult::Owned) {
                    ++unresolved;
                    continue;
                }
                ++owned;
                const int got = f.sheetOfS(arcOf(m, b.state[y][c].tri, b.state[y][c].b));
                wrongOwned += got != 0 && got != own;
                // Negative control: the unrestricted closest point of the placed control.
                const Vec3 q = k.t.curves.yarns[y].curve.ctrl[c];
                const MeshQuery::Hit fh = mq.closest(q);
                const auto& ft = m.tri[fh.tri];
                const TrianglePoint onTri = closestOnTriangle(q, m.pos[ft[0]], m.pos[ft[1]], m.pos[ft[2]]);
                const int free = f.sheetOfS(arcOf(m, fh.tri, onTri.bary));
                unrestrictedWrong += free != 0 && free != own;
            }
        }
        if (wrongOwned != 0 || owned == 0 || unrestrictedWrong == 0) {
            char m2[200];
            std::snprintf(m2, sizeof m2, "cell %.2f: owned %d, owned on the other sheet %d, unresolved %d, unrestricted wrong %d", cell,
                          owned, wrongOwned, unresolved, unrestrictedWrong);
            throw std::runtime_error(m2);
        }
        // Random queries near the bend from states on sheet 1: never Owned on sheet 2.
        const MeshSurface& s = *k.surface;
        uint64_t h = 1469598103934665603ull;
        auto rnd = [&]() {
            h = h * 6364136223846793005ull + 1442695040888963407ull;
            return double(h >> 11) / double(1ull << 53);
        };
        int bad = 0, turned = 0;
        for (int q = 0; q < 2000; ++q) {
            const double x0 = f.L - 3 * rnd(), y0 = 5 + 20 * rnd();
            const Vec3 start{x0, y0, 0};
            const MeshQuery::Hit hit = mq.closest(start);
            const auto& tv = m.tri[hit.tri];
            const TrianglePoint onTri = closestOnTriangle(start, m.pos[tv[0]], m.pos[tv[1]], m.pos[tv[2]]);
            const ProjState st = s.stateAt(hit.tri, onTri.bary, start);
            if (f.sheetOfS(arcOf(m, st.tri, st.b)) != 1) continue;
            const Vec3 qp{x0 + 2 * rnd() - 0.5, y0 + 2 * rnd() - 1, 1.4 * rnd()};
            ProjState next;
            const ProjResult r = s.continueFrom(st, qp, next, b.cp);
            if (r == ProjResult::Owned && f.sheetOfS(arcOf(m, next.tri, next.b)) == 2) ++bad;
            turned += r == ProjResult::NormalTurn;
        }
        CHECK(bad == 0);
        (void)turned;
    }
}

// Settle keeps the relax result whole, controls and projection states alike: a numerics wrapper flips a
// leg control through the surface in the relaxation (breaking that stitch's signature) and corrupts
// its state; settle stores exactly that result (nothing is restored from the state before the round,
// which put interpenetrating construction controls next to relaxed ones) and its validation reports
// the broken stitch.
TEST(settle_keeps_the_relaxed_state_whole) {
    PatternParams pp;
    pp.waleSpacing = 5.5;
    pp.courseSpacing = 3.4;
    pp.yarnDiameter = 1.16;  // the flat dense gauge FAILS clearance: settle has defects to work on
    const TriMesh m = gridMesh(60, 40, 1.0, [](double x, double y) { return Vec3{x, y, 0}; },
                               [](double x, double y) { return Vec2{x / 100, y / 100}; });
    MeshKnitParams p = tiledParams(pp);
    p.settle = false;
    PartKnit k;
    std::string err;
    CHECK(knitMesh(m, p, knit_test::testNumerics(), k, &err));
    CHECK(k.before.verdict != Status::Pass);
    // The target: every control of a stitch in the middle of a middle run, mirrored
    // through the surface (a purl-like stitch where a knit one was: its signature breaks).
    const int y = static_cast<int>(k.t.prog.turns.size() / 2), local = k.t.prog.turns[y].output / 2;
    const int c0 = local * kControlsPerStitch + 1, c1 = c0 + kControlsPerStitch;
    knit_test::SabotagedRound sabotaged(knit_test::cpuNumerics(), y, c0, c1);
    ValidationReport last;
    const SettleReport rep = settle({k.t.curves, k.t.pp, k.surface.get(), &k.t.bind}, SettleParams{}, sabotaged, &last);
    std::printf("    %d round(s); local signature %d FAIL after settle\n", rep.rounds, last.localSignature.fail);
    CHECK(rep.rounds == 1);
    CHECK(last.localSignature.fail > 0);
    CHECK(sabotaged.relaxed.size() == k.t.curves.yarns.size());
    long long differ = 0;
    for (size_t yy = 0; yy < k.t.curves.yarns.size(); ++yy)
        for (size_t c = 0; c < k.t.curves.yarns[yy].curve.ctrl.size(); ++c) {
            const Vec3 &q = k.t.curves.yarns[yy].curve.ctrl[c], &r = sabotaged.relaxed[yy][c];
            differ += !(q.x == r.x && q.y == r.y && q.z == r.z);
            differ += !(k.t.bind.state[yy][c] == sabotaged.relaxedBinding.state[yy][c] && k.t.bind.result[yy][c] == sabotaged.relaxedBinding.result[yy][c]);
        }
    CHECK(differ == 0);
}

// A rejected certified-repair attempt leaves controls and projection states as they were
// (every attempt is rejected here).
TEST(repair_rejected_attempt_keeps_projection_state) {
    PatternParams pp;
    pp.waleSpacing = 5.5;
    pp.courseSpacing = 3.4;
    pp.yarnDiameter = 1.16;
    const TriMesh m = gridMesh(60, 40, 1.0, [](double x, double y) { return Vec3{x, y, 0}; },
                               [](double x, double y) { return Vec2{x / 100, y / 100}; });
    MeshKnitParams p = tiledParams(pp);
    p.settle = false;
    PartKnit k;
    std::string err;
    CHECK(knitMesh(m, p, knit_test::testNumerics(), k, &err));
    const SurfaceBinding b0 = k.t.bind;
    std::vector<std::vector<Vec3>> c0;
    for (const YarnPath& yp : k.t.curves.yarns) c0.push_back(yp.curve.ctrl);
    knit_test::RejectEveryAttempt rejecting(knit_test::cpuNumerics());
    const RepairReport rr = certifiedRepair({k.t.curves, k.t.pp, k.surface.get(), &k.t.bind}, RepairParams{}, rejecting);
    std::printf("    %s\n", rr.detail.c_str());
    CHECK(rejecting.attempts > 0);
    for (size_t y = 0; y < c0.size(); ++y)
        for (size_t c = 0; c < c0[y].size(); ++c) {
            const Vec3 &a = c0[y][c], &q = k.t.curves.yarns[y].curve.ctrl[c];
            CHECK(a.x == q.x && a.y == q.y && a.z == q.z);
            CHECK(k.t.bind.state[y][c] == b0.state[y][c] && k.t.bind.result[y][c] == b0.result[y][c]);
        }
}

// --------------------------------------------------------------- publish ----

namespace pubfs = std::filesystem;

std::string freshDir(const char* name) {
    const std::string d = (pubfs::temp_directory_path() / name).string();
    std::error_code ec;
    pubfs::remove_all(d, ec);
    for (const auto& e : pubfs::directory_iterator(pubfs::temp_directory_path(), ec))
        if (e.path().filename().string().rfind(std::string(name) + ".", 0) == 0) pubfs::remove_all(e.path(), ec);
    return d;
}

void writeText(const std::string& path, const std::string& text) {
    pubfs::create_directories(pubfs::path(path).parent_path());
    std::ofstream f(path, std::ios::binary);
    f << text;
}

std::string readText(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    return std::string((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
}

int siblingCount(const std::string& target) {
    int n = 0;
    std::error_code ec;
    const std::string prefix = pubfs::path(target).filename().string() + ".";
    for (const auto& e : pubfs::directory_iterator(pubfs::path(target).parent_path(), ec))
        n += e.path().filename().string().rfind(prefix, 0) == 0;
    return n;
}

// Cancel mid-write (abort) and failure (no commit): the previous publication stays byte-
// identical and remains the cache hit; no staging directory is left behind.
TEST(export_cancel_keeps_previous_output) {
    const std::string t = freshDir("knit_s6_pub_a");
    std::string err;
    {
        Publisher p;
        CHECK(p.begin(t, &err));
        writeText(p.staging() + "/scene.json", "one");
        writeText(p.staging() + "/meshes/a.ply", std::string(1000, 'x'));
        CHECK(p.commit("k1", &err));
    }
    CHECK(isPublished(t, "k1"));
    {
        Publisher p;
        CHECK(p.begin(t, &err));
        writeText(p.staging() + "/scene.json", "two");  // cancelled half-way
        p.abort();
    }
    CHECK(readText(t + "/scene.json") == "one" && readText(t + "/meshes/a.ply") == std::string(1000, 'x'));
    CHECK(isPublished(t, "k1") && siblingCount(t) == 0);
}

TEST(export_failure_keeps_previous_output) {
    const std::string t = freshDir("knit_s6_pub_b");
    std::string err;
    {
        Publisher p;
        CHECK(p.begin(t, &err));
        writeText(p.staging() + "/scene.json", "one");
        CHECK(p.commit("k1", &err));
    }
    {
        Publisher p;  // a writer that fails: leaves without commit (destructor aborts)
        CHECK(p.begin(t, &err));
        writeText(p.staging() + "/scene.json", "broken");
    }
    CHECK(readText(t + "/scene.json") == "one" && isPublished(t, "k1") && siblingCount(t) == 0);
    // A successful second publication replaces the first one completely.
    {
        Publisher p;
        CHECK(p.begin(t, &err));
        writeText(p.staging() + "/scene.json", "two");
        CHECK(p.commit("k2", &err));
    }
    CHECK(readText(t + "/scene.json") == "two" && isPublished(t, "k2") && !isPublished(t, "k1") && siblingCount(t) == 0);
}

// An unfinished (cancelled) map set is never a cache hit; neither is a marked set whose
// content changed or grew.
TEST(plymap_cancel_never_cache_hit) {
    const std::string t = freshDir("knit_s6_pub_c");
    std::string err;
    {
        Publisher p;
        CHECK(p.begin(t, &err));
        writeText(p.staging() + "/ply_ao.png", "png");
        // cancelled before commit
    }
    CHECK(!isPublished(t, "maps-3-swept"));
    pubfs::create_directories(t);
    writeText(t + "/ply_ao.png", "png");  // the files without a marker
    CHECK(!isPublished(t, "maps-3-swept"));
}

TEST(plymap_corrupt_marker_not_hit) {
    const std::string t = freshDir("knit_s6_pub_d");
    std::string err;
    {
        Publisher p;
        CHECK(p.begin(t, &err));
        writeText(p.staging() + "/ply_ao.png", "png-bytes");
        writeText(p.staging() + "/ply_normal.png", "normal-bytes");
        CHECK(p.commit("maps-3-swept", &err));
    }
    CHECK(isPublished(t, "maps-3-swept"));
    CHECK(!isPublished(t, "maps-4-swept"));
    writeText(t + "/ply_ao.png", "png-byteS");  // same size, different content
    std::string why;
    CHECK(!isPublished(t, "maps-3-swept", &why) && why.find("checksum") != std::string::npos);
    writeText(t + "/ply_ao.png", "png-bytes");
    CHECK(isPublished(t, "maps-3-swept"));
    writeText(t + "/extra.txt", "x");
    CHECK(!isPublished(t, "maps-3-swept", &why) && why.find("unlisted") != std::string::npos);
}

// A crash between moving the old output aside and renaming the staging directory: the
// startup recovery restores the old output and removes the staging leftovers.
TEST(publish_recovery_from_old_dir) {
    const std::string t = freshDir("knit_s6_pub_e");
    std::string err;
    {
        Publisher p;
        CHECK(p.begin(t, &err));
        writeText(p.staging() + "/scene.json", "one");
        CHECK(p.commit("k1", &err));
    }
    // The publisher marks the previous output as its own before moving it aside.
    CHECK(markPublishOwned(t, t, &err));
    pubfs::rename(t, t + ".old-00000000c0a5ed01");
    pubfs::create_directories(t + ".staging-00000000c0a5ed02");
    writeText(t + ".staging-00000000c0a5ed02/scene.json", "half");
    CHECK(markPublishOwned(t + ".staging-00000000c0a5ed02", t, &err));
    const std::string log = recoverPublished(t);
    CHECK(!log.empty());
    CHECK(readText(t + "/scene.json") == "one" && isPublished(t, "k1") && siblingCount(t) == 0);
}

TEST(staging_cleanup_on_startup) {
    const std::string t = freshDir("knit_s6_pub_f");
    std::string err;
    {
        Publisher p;
        CHECK(p.begin(t, &err));
        writeText(p.staging() + "/a.txt", "a");
        CHECK(p.commit("k", &err));
    }
    pubfs::create_directories(t + ".staging-0000000000000001");
    pubfs::create_directories(t + ".old-0000000000000002");
    CHECK(markPublishOwned(t + ".staging-0000000000000001", t, &err));
    CHECK(markPublishOwned(t + ".old-0000000000000002", t, &err));
    recoverPublished(t);
    CHECK(siblingCount(t) == 0 && isPublished(t, "k"));
}

// ----------------------------------------------------------------- fibres ----

// Fibres on the straight strand demo (as knit_check --strand-demo, without the export's
// serialized checks); bendDeg bends it into an arc.
FibreSet strandFibres(const FibreParams& fp, double bendDeg, CurveSet* curves = nullptr) {
    StrandDemo d = makeStrandDemo(12.0, 1.16);
    if (bendDeg != 0) {
        const double Rb = 12.0 / (bendDeg * kPi / 180);
        YarnPath& y = d.curves.yarns[0];
        for (Vec3& c : y.curve.ctrl) c = {Rb * std::sin(c.x / Rb), Rb * (1 - std::cos(c.x / Rb)), 0};
        y.finalize();
    }
    SurfaceParams sp;
    sp.compression = false;
    const SurfaceModel sm = buildSurface(d.curves, d.pp, sp);
    double kappaMax = 0;
    for (const YarnPath& y : d.curves.yarns)
        kappaMax = std::max(kappaMax, std::min(checkCurvature(y, d.pp.radius(), Allowances::forExtent(d.curves.extent())).kappaHi, 1e3));
    FibreParams f = fp;
    FibreFrame fr;
    FibreSetup fset;
    fset.tubeRadius = d.pp.radius();
    fset.sides = kExportSides;
    fset.ring = kExportRing;
    fset.kappaMax = kappaMax;
    setupPieceFibres(f, fr, d.pp, fset);
    fr.fabric = [](const Vec3& p) { return FabricPoint{{p.x, p.y, 0}, {0, 0, 1}}; };
    if (curves) *curves = d.curves;
    return growYarnFibres(d.curves, sm, d.pp, f, &fr);
}

bool sameFibres(const FibreSet& a, const FibreSet& b) {
    if (a.fibres.size() != b.fibres.size()) return false;
    for (size_t k = 0; k < a.fibres.size(); ++k) {
        const Fibre &x = a.fibres[k], &y = b.fibres[k];
        if (x.radius0 != y.radius0 || x.radius1 != y.radius1 || x.population != y.population || fibrePoints(a, x) != fibrePoints(b, y)) return false;
    }
    return true;
}

// Overrides equal to the shared values draw bit-identically (every channel), and the
// generated set is identical; the legacy fixtures are checked by knit_tests.
TEST(fibres_stay_in_front) {
    // The Fibres step's reference patch faces +z; its back is the lowest point of the yarn tubes. With
    // every wool preset no fibre point from the curve start on (the buried root is exempt), and no
    // segment midpoint, lies behind the back: neither population grows out of the back of the fabric.
    // Fibres inside the fabric's thickness (rooted on yarn seen through the stitch gaps) stay.
    const std::vector<FilePreset>& wool = presets(PresetGroup::Fibres);
    CHECK(!wool.empty());
    int behindTotal = 0;
    for (const FilePreset& p : wool) {
        SceneSettings s;
        applyPreset(PresetGroup::Fibres, s, p);
        std::string err;
        const std::optional<WoolSample> grown = growWoolSample(s, &err);
        CHECK(grown.has_value());
        if (!grown) continue;
        const WoolSample& sample = *grown;
        double back = 1e300;
        for (const YarnSurface& ys : sample.surface.yarns)
            for (const Ring& R : ys.rings) back = std::min(back, R.c.z - std::max(R.a, R.b));
        for (const FibreSet* set : {&sample.along, &sample.stray}) {
            int behind = 0;
            double lowest = 1e300;
            for (const Fibre& f : set->fibres) {
                const std::span<const Vec3f> pts = set->pointsOf(f);
                double low = pts[1].z;
                for (size_t k = 2; k < pts.size(); ++k) low = std::min({low, double(pts[k].z), 0.5 * (double(pts[k - 1].z) + pts[k].z)});
                behind += low < back - 0.01;  // the generator's capsules reach a little below the rings
                lowest = std::min(lowest, low);
            }
            std::printf("  %-14s %-8s %6zu fibres, %4d with a point behind the back at z = %.3f mm (lowest %.3f)\n", p.label.c_str(),
                        set == &sample.along ? "directed" : "stray", set->fibres.size(), behind, back, set->fibres.empty() ? 0.0 : lowest);
            behindTotal += behind;
            CHECK(!set->fibres.empty());
        }
    }
    CHECK(behindTotal == 0);
}

TEST(fibre_random_overrides_inherit_identical) {
    const FibreParams fp = pieceFibres();
    FibreParams eq = fp;
    eq.randomLengthMean = fp.lengthMean;
    eq.randomLengthMin = fp.lengthMin;
    eq.randomLengthMax = fp.lengthMax;
    eq.randomLengthSpread = fp.lengthSpread;
    eq.randomThickness = fp.thickness;
    eq.randomThicknessSigma = fp.thicknessSigma;
    CHECK(hasRandomOverrides(eq) && !hasRandomOverrides(fp));
    for (int i = 0; i < 20000; ++i) {
        const FibreDraw a = drawFibre(fp, 7, i % 3, i), b = drawFibre(eq, 7, i % 3, i);
        CHECK(a.length == b.length && a.radiusFactor == b.radiusFactor && a.population == b.population && a.key == b.key &&
              a.psi == b.psi && a.lift == b.lift && a.elev == b.elev && a.bendAmp == b.bendAmp);
    }
    CHECK(sameFibres(strandFibres(fp, 0), strandFibres(eq, 0)));
}

// Real overrides change the random population only: its lengths stay in the new bounds,
// its radii follow the new thickness; directed draws stay bit-identical.
TEST(fibre_random_overrides_applied) {
    const FibreParams fp = pieceFibres();
    FibreParams ov = fp;
    ov.randomLengthMin = 1.2;
    ov.randomLengthMax = 2.2;
    ov.randomLengthMean = 1.6;
    ov.randomThickness = 0.1;
    ov.randomThicknessSigma = 0;
    int random = 0, directed = 0;
    for (int i = 0; i < 20000; ++i) {
        const FibreDraw a = drawFibre(fp, 3, 0, i), b = drawFibre(ov, 3, 0, i);
        CHECK(a.population == b.population);
        if (b.population == FibrePopulation::Directed) {
            ++directed;
            CHECK(a.length == b.length && a.radiusFactor == b.radiusFactor);
        } else {
            ++random;
            CHECK(b.length >= 1.2 && b.length <= 2.2 && b.radiusFactor == 1.0);
        }
    }
    CHECK(random > 1000 && directed > 1000);
    const FibreSet fs = strandFibres(ov, 0);
    int nr = 0;
    for (const Fibre& f : fs.fibres)
        if (f.population == FibrePopulation::Random) {
            ++nr;
            CHECK_NEAR(f.radius0, 0.05, 1e-7);  // float32
            CHECK(f.targetLength >= 1.2 && f.targetLength <= 2.2);
        } else {
            CHECK(f.radius0 >= 0.5 * fp.thickness * std::exp(-2 * fp.thicknessSigma) - 1e-7);
            CHECK(f.radius0 <= 0.5 * fp.thickness * std::exp(2 * fp.thicknessSigma) + 1e-7);
        }
    CHECK(nr > 0);
}

// The settings' thickness range <-> (t, sigma): exact range, both ways.
TEST(fibre_thickness_minmax_roundtrip) {
    SceneSettings::Population p;
    for (const auto& [mn, mx] : {std::pair{0.05, 0.09}, {0.07, 0.07}, {0.02, 0.2}}) {
        p.setThicknessRange(mn, mx);
        CHECK_NEAR(p.thicknessMin(), mn, 1e-14);
        CHECK_NEAR(p.thicknessMax(), mx, 1e-14);
    }
    // The drawn radii stay inside the range (the factor is clamped at +-2 sigma).
    FibreParams fp = pieceFibres();
    p.setThicknessRange(0.05, 0.09);
    fp.thickness = p.thickness;
    fp.thicknessSigma = p.thicknessSigma;
    for (int i = 0; i < 20000; ++i) {
        const double d = fp.thickness * drawFibre(fp, 11, 0, i).radiusFactor;
        CHECK(d >= 0.05 - 1e-12 && d <= 0.09 + 1e-12);
    }
}

// On a straight yarn without bow and helix, the measured angle of a directed fibre IS its
// nominal lift (the curve is a straight line in yarn coordinates): within [liftMin, liftMax].
TEST(fibre_angle_nominal_equals_physical_straight) {
    FibreParams fp = pieceFibres();
    fp.randomFraction = 0;
    fp.bowTheta = fp.bowRadial = 0;
    fp.helixDeg = fp.helixSigmaDeg = 0;
    fp.minSegments = fp.maxSegments = 1;  // knit_check --fiber-variant straight: one visible chord
    fp.minTurnDeg = 0;
    CurveSet cs;
    const FibreSet fs = strandFibres(fp, 0, &cs);
    const FibreAngles a = measureFibreAngles(cs, fs);
    CHECK(a.pop[0].fibres > 50);
    for (const Fibre& f : fs.fibres) {
        FibreSet one;
        one.addFibre(f, fs.pointsOf(f));
        const FibreAngles x = measureFibreAngles(cs, one);
        CHECK(x.pop[0].emergence[2] >= fp.liftMinDeg - 0.5 && x.pop[0].emergence[2] <= fp.liftMaxDeg + 0.5);
        CHECK_NEAR(x.pop[0].emergence[2], f.lift, 0.5);
    }
}

// On a bent yarn with a strong radial bow the physical angle along the fibre exceeds the
// nominal maximum lift: the setting is a NOMINAL bound, the measurement is reported.
TEST(fibre_angle_measured_on_curved_bowed) {
    FibreParams fp = pieceFibres();
    fp.randomFraction = 0;
    fp.bowRadial = 0.3;
    CurveSet cs;
    const FibreSet fs = strandFibres(fp, 180, &cs);
    const FibreAngles a = measureFibreAngles(cs, fs);
    CHECK(a.pop[0].fibres > 20);
    CHECK(a.pop[0].along[2] > fp.liftMaxDeg);
    CHECK(a.pop[0].along[2] <= 90);
}

namespace {
// Directed fibres only, long enough for their rise to matter (Standard wool's lengths).
FibreParams longDirectedFibres() {
    FibreParams fp = pieceFibres();
    fp.randomFraction = 0;
    fp.lengthMin = 1.2;
    fp.lengthMean = 2.6;
    fp.lengthMax = 5.0;
    fp.maxSegments = 24;
    return fp;
}
// Polar angle about +x: makeStrandDemo's strand is the x axis.
double angleAboutX(const Vec3& p) { return std::atan2(p.z, p.y); }
double wrapAngle(double a) { return std::remainder(a, 2 * kPi); }
} // namespace

// fibres.directed.hug_height: on a straight strand every visible vertex of a hugging fibre stays
// within (1 + H) r of the axis and the last chord lies along the yarn, where the linear rise of the
// same draws ends far out.
TEST(fibre_hug_levels_off) {
    const double r = makeStrandDemo().pp.radius();
    FibreParams fp = longDirectedFibres();
    const FibreSet linear = strandFibres(fp, 0);
    fp.hugHeight = 0.15;
    const FibreSet hug = strandFibres(fp, 0);
    auto maxLambda = [r](const FibreSet& fs) {
        double m = 0;
        for (const Fibre& f : fs.fibres) {
            const auto pts = fs.pointsOf(f);
            for (size_t k = 1; k < pts.size(); ++k) m = std::max(m, std::hypot(pts[k].y, pts[k].z) / r);
        }
        return m;
    };
    std::vector<double> tipLift;  // [deg] the last chord against the yarn axis's radial plane
    for (const Fibre& f : hug.fibres) {
        const auto pts = hug.pointsOf(f);
        const Vec3 a = toVec3(pts[pts.size() - 2]), b = toVec3(pts.back());
        const double radial = std::hypot(b.y, b.z) - std::hypot(a.y, a.z);
        tipLift.push_back(std::asin(std::clamp(radial / norm(b - a), -1.0, 1.0)) * 180 / kPi);
    }
    std::sort(tipLift.begin(), tipLift.end());
    const double medianTip = tipLift.empty() ? 90 : tipLift[tipLift.size() / 2];
    std::printf("  linear: %zu fibres, max lambda %.3f; hugging: %zu fibres, max lambda %.3f, median tip chord %.2f deg\n",
                linear.fibres.size(), maxLambda(linear), hug.fibres.size(), maxLambda(hug), medianTip);
    CHECK(hug.fibres.size() > 50);
    CHECK(maxLambda(linear) > 1.5);
    CHECK(maxLambda(hug) <= 1 + fp.hugHeight + 1e-4);
    CHECK(medianTip < 3);
    // The analytic derivative of the levelled (and ramped) curve against a central difference.
    const StrandDemo d = makeStrandDemo();
    SurfaceParams sp;
    sp.compression = false;
    const SurfaceModel sm = buildSurface(d.curves, d.pp, sp);
    DirectedCurve c;
    c.cs = &d.curves;
    c.ys = &sm.yarns[0];
    c.sE = 6.0;
    c.thetaE = 0.9;
    c.S = 1.5;
    c.p = std::tan(55 * kPi / 180);
    c.q = std::tan(25 * kPi / 180);
    c.bowTheta = 0.09;
    c.r = r;
    c.hug = 0.15;
    for (bool ramp : {false, true})
        for (int sgn : {1, -1}) {
            c.rampWinding = ramp;
            c.sgn = sgn;
            for (double u = 0.05; u < 1; u += 0.1) {
                const double h = 1e-6;
                const Vec3 fd = (c.X(u + h) - c.X(u - h)) / (2 * h), an = c.dX(u);
                CHECK(norm(fd - an) <= 1e-6 * std::max(1.0, norm(an)));
            }
        }
}

// fibres.directed.follow_plies: hugging fibres wind about the strand in the ply's sense and, once
// the winding has ramped in (4 kWindRampLength from the anchor: above 96 %), at the ply's rate
// d theta / ds = h 2 pi / pitch (plies.cpp samplePlyCentres), for both handednesses.
TEST(fibre_follows_ply_helix) {
    const StrandDemo d = makeStrandDemo();
    SurfaceParams sp;
    sp.compression = false;
    const SurfaceModel sm = buildSurface(d.curves, d.pp, sp);
    const double pitch = 4.0;  // 42 deg on this thin strand: 2.5 mm would turn too sharply at 0.45 mm chords
    for (int h : {1, -1}) {
        auto plyAngle = [&](double s) { return angleAboutX(yarnCoordinates(d.curves, sm.yarns[0], s, h * 2 * kPi * s / pitch, 0.5).X); };
        const double plyRate = wrapAngle(plyAngle(6.05) - plyAngle(6.0)) / 0.05;  // [rad / mm] along +x
        FibreParams fp = longDirectedFibres();
        fp.lengthMin = 6;
        fp.lengthMean = 8;
        fp.lengthMax = 10;
        fp.maxSegments = 30;
        fp.hugHeight = 0.15;
        fp.plyPitch = pitch;
        fp.handedness = h;
        fp.helixDeg = 0;
        fp.helixSigmaDeg = 3;
        const FibreSet fs = strandFibres(fp, 0);
        const FibreReport& R = fs.rep;
        std::printf("  %d candidates, %d exposed, %d accepted; refused: %d yarn end, %d no emergence, %d too sharp, %d approx\n",
                    R.candidates, R.exposed, R.accepted, R.yarnEnd, R.noEmergence, R.turnTooSharp, R.approxError);
        std::vector<double> rates;
        int wrongSense = 0;
        for (const Fibre& f : fs.fibres) {
            const auto pts = fs.pointsOf(f);
            size_t k0 = 1;
            while (k0 + 1 < pts.size() && std::fabs(pts[k0].x - pts[1].x) < 4 * kWindRampLength) ++k0;
            double turn = 0;
            for (size_t k = k0 + 1; k < pts.size(); ++k) turn += wrapAngle(angleAboutX(toVec3(pts[k])) - angleAboutX(toVec3(pts[k - 1])));
            const double dx = pts.back().x - pts[k0].x;
            if (std::fabs(dx) < 0.5) continue;  // too little wound length to measure a rate
            rates.push_back(turn / dx);
            wrongSense += rates.back() * plyRate <= 0;
        }
        std::sort(rates.begin(), rates.end());
        const double median = rates.empty() ? 0 : rates[rates.size() / 2];
        std::printf("  h %+d: ply %.3f rad/mm; %zu fibres measured, median %.3f rad/mm, %d winding the other way\n", h, plyRate,
                    rates.size(), median, wrongSense);
        CHECK(rates.size() > 50);
        CHECK(wrongSense == 0);
        CHECK(std::fabs(median / plyRate - 1) < 0.05);
    }
}

// ---------------------------------------------------------------- builder ----

// A synthetic slow stage: `steps` checkpoints 2 ms apart, then a result for the job.
std::shared_ptr<BuiltScene> slowJob(const JobRequest& j, const JobControl& c, int steps, std::vector<int>* order = nullptr,
                                    std::mutex* om = nullptr) {
    for (int k = 0; k < steps; ++k) {
        cancelCheckpoint();
        std::this_thread::sleep_for(std::chrono::milliseconds(2));
    }
    auto b = std::make_shared<BuiltScene>();
    b->revision = c.revision;
    b->kind = jobKind(j);
    b->settings = jobSettings(j);
    if (order) {
        std::lock_guard<std::mutex> lk(*om);
        order->push_back(c.revision);
    }
    return b;
}

void waitIdle(const SceneBuilder& b, double seconds = 20) {
    const auto t0 = std::chrono::steady_clock::now();
    while (b.snapshot().busy && std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() < seconds)
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
    CHECK(!b.snapshot().busy);
}

// Queues a knit (or an export) of `scene`'s default settings; returns its revision.
int startKnit(SceneBuilder& b, const std::string& scene = "teddy") {
    auto s = std::make_shared<SceneSettings>();
    s->scene = scene;
    return b.startKnit(KnitRequest{s});
}
int startExport(SceneBuilder& b) {
    ExportRequest x;
    x.settings = std::make_shared<SceneSettings>();
    return b.startExport(x);
}
int runningRevision(const SceneBuilder& b) { return b.snapshot().progress.revision; }

TEST(builder_rapid_edits_publish_latest_only) {
    SceneBuilder b([](const JobRequest& j, const JobControl& c, std::string*) { return slowJob(j, c, 100); });
    for (int r = 1; r <= 6; ++r) {
        CHECK(startKnit(b) == r);
        std::this_thread::sleep_for(std::chrono::milliseconds(3));
    }
    waitIdle(b);
    CHECK(b.snapshot().latest && b.snapshot().latest->revision == 6);
    CHECK(b.snapshot().published == 1);  // the earlier revisions were cancelled or replaced
    CHECK(b.snapshot().lastOutcome == JobOutcome::Done);
}

TEST(builder_scene_switch_cancels) {
    std::atomic<int> started{0};
    SceneBuilder b([&](const JobRequest& j, const JobControl& c, std::string*) {
        ++started;
        return slowJob(j, c, jobSettings(j).scene == "teddy" ? 2000 : 10);
    });
    CHECK(startKnit(b, "teddy") == 1);
    while (runningRevision(b) != 1) std::this_thread::sleep_for(std::chrono::milliseconds(1));
    const auto t0 = std::chrono::steady_clock::now();
    CHECK(startKnit(b, "wool-swatch") == 2);
    waitIdle(b);
    const double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    CHECK(b.snapshot().latest->revision == 2 && b.snapshot().latest->settings.scene == "wool-swatch");
    CHECK(b.snapshot().published == 1 && started == 2);
    CHECK(s < 2.0);  // the teddy job (4 s) was cancelled at a checkpoint
}

// A job that cannot be cancelled (no checkpoints) completes after a newer request: it is
// published, then replaced by the newer one; publications never go back in revision.
TEST(builder_stale_completion_ignored) {
    std::vector<int> order;
    std::mutex om;
    SceneBuilder b([&](const JobRequest& j, const JobControl& c, std::string*) {
        if (c.revision == 1) {
            std::this_thread::sleep_for(std::chrono::milliseconds(150));  // no checkpoint
            return slowJob(j, c, 0, &order, &om);
        }
        return slowJob(j, c, 5, &order, &om);
    });
    CHECK(startKnit(b) == 1);
    while (runningRevision(b) != 1) std::this_thread::sleep_for(std::chrono::milliseconds(1));
    CHECK(startKnit(b) == 2);
    waitIdle(b);
    CHECK(b.snapshot().latest->revision == 2);
    CHECK(std::is_sorted(order.begin(), order.end()));
}

TEST(builder_shutdown_while_running_joins) {
    auto t0 = std::chrono::steady_clock::now();
    {
        SceneBuilder b([](const JobRequest& j, const JobControl& c, std::string*) { return slowJob(j, c, 100000); });
        CHECK(startKnit(b) == 1);
        while (runningRevision(b) != 1) std::this_thread::sleep_for(std::chrono::milliseconds(1));
        t0 = std::chrono::steady_clock::now();
    }  // destructor: cancel + join
    CHECK(std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() < 1.0);
}

TEST(builder_stage_exception_reported) {
    SceneBuilder b([](const JobRequest& j, const JobControl& c, std::string* err) -> std::shared_ptr<BuiltScene> {
        if (c.revision == 2) throw std::runtime_error("stage exploded");
        if (c.revision == 3) {
            *err = "refused";
            return nullptr;
        }
        return slowJob(j, c, 1);
    });
    CHECK(startKnit(b) == 1);
    waitIdle(b);
    CHECK(b.snapshot().latest->revision == 1);
    CHECK(startKnit(b) == 2);
    waitIdle(b);
    CHECK(b.snapshot().lastOutcome == JobOutcome::Failed && b.snapshot().lastError == "stage exploded");
    CHECK(b.snapshot().latest->revision == 1);  // the previous result stays
    CHECK(startKnit(b) == 3);
    waitIdle(b);
    CHECK(b.snapshot().lastOutcome == JobOutcome::Failed && b.snapshot().lastError == "refused" && b.snapshot().latest->revision == 1);
}

// An export is a user action: an edit during it does not cancel it (it runs from its own
// snapshot); the edit's build follows.
TEST(builder_export_not_cancelled_by_edit) {
    std::vector<int> order;
    std::mutex om;
    SceneBuilder b([&](const JobRequest& j, const JobControl& c, std::string*) { return slowJob(j, c, jobKind(j) == BuiltScene::Kind::Export ? 60 : 5, &order, &om); });
    CHECK(startExport(b) == 1);
    while (runningRevision(b) != 1) std::this_thread::sleep_for(std::chrono::milliseconds(1));
    CHECK(startKnit(b) == 2);
    waitIdle(b);
    CHECK(order.size() == 2 && order[0] == 1 && order[1] == 2);
    CHECK(b.snapshot().published == 2 && b.snapshot().latest->revision == 2);
}

// ------------------------------------------------------------ final surface ----

bool sameControls(const CurveSet& a, const CurveSet& b) {
    if (a.yarns.size() != b.yarns.size()) return false;
    for (size_t y = 0; y < a.yarns.size(); ++y) {
        const auto &p = a.yarns[y].curve.ctrl, &q = b.yarns[y].curve.ctrl;
        if (p.size() != q.size() || std::memcmp(p.data(), q.data(), p.size() * sizeof(Vec3)) != 0) return false;
    }
    return true;
}

// A knitted plane as a tiled scene part (centre at the origin).
std::unique_ptr<ScenePart> planePart(const PatternParams& pp, const FinalSurfaceSpec* surface = nullptr) {
    auto r = std::make_unique<ScenePart>();
    r->name = "plane";
    r->mesh = std::make_shared<const TriMesh>(gridMesh(48, 34, 1.0, [](double x, double y) { return Vec3{x, y, 0}; },
                                                       [](double x, double y) { return Vec2{x / 100, y / 100}; }));
    MeshKnitParams p = tiledParams(pp);
    p.settle = false;
    if (surface) p.surface = *surface;
    r->knit = std::make_unique<PartKnit>();
    std::string err;
    if (!knitMesh(*r->mesh, p, knit_test::testNumerics(), *r->knit, &err)) r->knit.reset();
    r->yarnD = pp.yarnDiameter;
    return r;
}

// e = 0 (no flattening, one strand): the rendered surface IS the centreline model at D, so
// every consumer sees exactly the pre-S6 parameters (the teddy parity run covers the pieces).
TEST(final_surface_e0_bit_identical) {
    const PatternParams pp = expandedD1();
    FinalSurfaceSpec f;
    f.D = pp.yarnDiameter;
    f.c = pp.clearance;
    std::string why;
    CHECK(f.valid(&why));
    CHECK(f.rendered(pp) == pp);
    CHECK(f.supportRadius() == pp.radius());
    CHECK(f.referenceCheckApplicable());
    const SceneSettings s;  // defaults: no flattening, one strand
    const FinalSurfaceSpec g = finalSurface(s, pp.yarnDiameter);
    CHECK(g.inflate == 0 && !g.compression && g.plies == 1 && g.rendered(pp) == pp);
    // Tiling: R_s from the spec equals the explicit pre-S6 default r.
    MeshKnitParams a = tiledParams(pp), b = tiledParams(pp);
    a.settle = b.settle = false;
    b.tile.supportRadius = pp.radius();
    const TriMesh m = gridMesh(48, 34, 1.0, [](double x, double y) { return Vec3{x, y, 0}; },
                               [](double x, double y) { return Vec2{x / 100, y / 100}; });
    PartKnit ka, kb;
    CHECK(knitMesh(m, a, knit_test::testNumerics(), ka, &why) && knitMesh(m, b, knit_test::testNumerics(), kb, &why));
    CHECK(sameControls(ka.t.curves, kb.t.curves) && ka.verdict == kb.verdict);
    // The rendered part is the knitted yarn itself.
    const auto r = planePart(pp);
    CHECK(r->knit != nullptr);
    const std::optional<RenderedPart> rendered = renderedPart(*r, s, &why);
    CHECK(rendered.has_value());
    if (!rendered) return;
    const RenderedPart& rp = *rendered;
    CHECK(rp.strands == &r->knit->t.curves && rp.strandPp == pp && rp.plies == nullptr);
    CHECK(rp.ring == s.tubeRing && rp.sides == s.tubeSides);  // the settings' tube (default: Reduced)
}

// Inflation e with compression on two validated straight yarns at d = D + c: the rendered
// circles (D_r > d) overlap, the M3 flattening brings every ring's minor semi-axis below
// r_R = D_r / 2 and gate (d) runs on that surface. Known gap: the compression is geometric, not
// mechanics.
TEST(inflate_flattening_known_gap) {
    StrandDemo d = makeStrandDemo(12.0, 1.0);
    d.pp.clearance = 0.1;
    YarnPath y2 = d.curves.yarns[0];
    for (Vec3& c : y2.curve.ctrl) c.y += d.pp.yarnDiameter + d.pp.clearance;
    y2.finalize();
    d.curves.yarns.push_back(y2);
    FinalSurfaceSpec f;
    f.D = d.pp.yarnDiameter;
    f.c = d.pp.clearance;
    f.inflate = 0.2;
    std::string why;
    CHECK(!f.valid(&why));  // inflation without compression: overlapping tubes
    f.compression = true;
    CHECK(f.valid(&why));
    {
        FinalSurfaceSpec g = f;
        g.plies = 3;
        CHECK(!g.valid(&why));  // the plies are round
    }
    const PatternParams ppR = f.rendered(d.pp);
    const double rR = ppR.radius();
    CHECK(2 * rR > d.pp.yarnDiameter + d.pp.clearance);
    SurfaceParams sp;
    sp.compression = true;
    const SurfaceModel sm = buildSurface(d.curves, ppR, sp);
    double minB = 1e300, maxA = 0;
    for (const YarnSurface& ys : sm.yarns)
        for (const Ring& r : ys.rings) {
            minB = std::min(minB, std::min(r.a, r.b));
            maxA = std::max(maxA, std::max(r.a, r.b));
        }
    CHECK(minB < rR - 1e-3);
    CHECK(maxA <= f.supportRadius() + 1e-12);  // the B1 support radius bounds the long axis
    const SurfaceCheck c = checkSurface(d.curves, sm);
    CHECK(c.pointsTested > 0);
    std::printf("    flattened: min semi-axis %.4f (r_R %.4f), max %.4f; gate d fold %s collisions %s (%ld penetrating, max %.2e)\n", minB,
                rR, maxA, toString(c.fold), toString(c.collisions), c.penetratingPoints, c.maxDepth);
    // Negative control: the same rendered diameter without compression collides.
    sp.compression = false;
    const SurfaceModel round = buildSurface(d.curves, ppR, sp);
    CHECK(checkSurface(d.curves, round).collisions == Status::Fail);
}

// The round-tube fibre reference check is not applicable to compressed tubes: its verdict
// is UNRESOLVED (never PASS), in the spec and in the tiled pipeline's report.
TEST(reference_check_unresolved_when_compressed) {
    FinalSurfaceSpec f;
    CHECK(f.referenceCheckApplicable());
    f.compression = true;
    CHECK(!f.referenceCheckApplicable());
    const PatternParams pp = expandedD1();
    SceneSettings s;
    s.flatten = true;
    s.flattenInflate = 0.1;
    const FinalSurfaceSpec fs = finalSurface(s, pp.yarnDiameter);
    const auto r = planePart(pp, &fs);
    CHECK(r->knit != nullptr);
    if (!r->knit) return;
    FibreParams fp = pieceFibres();
    fp.density *= 0.25;
    const SceneSolids none;
    FibreRunOptions checked;
    checked.fibreCheck = true;
    std::string err;
    FibreSet fibres;
    std::unique_ptr<PartFibreJob> job = PartFibreJob::prepare({*r, 0, s, fp, none, 0.025, checked}, &err);
    CHECK(job && job->generate(0.0f, 2.0f, fibres, &err));
    const FinishedFibres compressed = job->finish(std::move(fibres));
    CHECK(!compressed.fibres.fibres.empty());
    CHECK(compressed.checks.find("reference check: UNRESOLVED (not applicable") != std::string::npos);
    // Round tubes: the check runs.
    s.flatten = false;
    const auto r0 = planePart(pp);
    FibreSet roundFibres;
    job = PartFibreJob::prepare({*r0, 0, s, fp, none, 0.2, checked}, &err);
    CHECK(job && job->generate(0.0f, 2.0f, roundFibres, &err));
    const FinishedFibres round = job->finish(std::move(roundFibres));
    CHECK(round.checks.find("reference check: UNRESOLVED (not applicable") == std::string::npos);
    CHECK(round.checks.find("reference check: ") != std::string::npos);
}

// Explicit plies on a tiled part: the rendered strands are the plies (count x yarns), at the
// ply export resolution, inside the validated envelope; the fibres stay authored on the yarn.
TEST(tiled_rendered_part_plies) {
    const PatternParams pp = expandedD1();
    SceneSettings s;
    s.plies = 3;
    const auto r = planePart(pp);
    CHECK(r->knit != nullptr);
    if (!r->knit) return;
    std::string err;
    const std::optional<RenderedPart> rendered = renderedPart(*r, s, &err);
    CHECK(rendered.has_value());
    if (!rendered) return;
    const RenderedPart& rp = *rendered;
    CHECK(rp.plies != nullptr && rp.yarn == &r->knit->t.curves);
    CHECK(rp.strands->yarns.size() == 3 * r->knit->t.curves.yarns.size());
    CHECK(rp.strandPp.yarnDiameter < pp.yarnDiameter && rp.ring == s.plyTubeRing && rp.sides == s.plyTubeSides);
    CHECK(rp.fs.supportRadius() == 0.5 * pp.yarnDiameter);
    std::printf("    plies: r_p %.3f mm, %d folds FAIL, %d UNRESOLVED\n", rp.plies->layout.rp, rp.plyFolds, rp.plyFoldsUnresolved);
    s.flatten = true;
    s.flattenInflate = 0.1;
    CHECK(!renderedPart(*r, s, &err));  // flattening with explicit plies: refused
}

// A yarn compressed to fit: the plies knit the same yarn (k D) and wrap its drawn diameter.
TEST(plies_on_compressed_yarn) {
    const PatternParams pp = expandedD1();
    SceneSettings s;
    s.diameter = 2.4;  // too thick for 8 x 6.5 mm (round up to 1.74 mm)
    s.wale = pp.waleSpacing;
    s.course = pp.courseSpacing;
    s.clearance = pp.clearance;
    SceneSettings tube = s;
    s.plies = 3;
    CHECK(yarnFit(s).compressed() && yarnFit(s).k == yarnFit(tube).k && !yarnFit(s).tooTight);
    CHECK(yarnFitNote(s).find("plies") != std::string::npos && yarnFitNote(tube).find("plies") == std::string::npos);
    const auto r = planePart(pp);
    CHECK(r->knit != nullptr);
    if (!r->knit) return;
    std::string err;
    const std::optional<RenderedPart> rendered = renderedPart(*r, s, &err);
    CHECK(rendered.has_value());
    if (!rendered) return;
    const RenderedPart& rp = *rendered;
    const double Dr = pp.yarnDiameter / yarnFit(s).k;
    PlyGeometry g;
    g.count = s.plies;
    g.fill = s.plyFill;
    CHECK(std::fabs(rp.ppR.yarnDiameter - Dr) < 1e-9 && std::fabs(rp.fs.supportRadius() - 0.5 * Dr) < 1e-9);
    CHECK(rp.plies != nullptr && std::fabs(rp.plies->layout.rp - plyLayout(g, 0.5 * Dr).rp) < 1e-9);
    s.flatten = true;
    s.flattenInflate = 0.1;
    CHECK(!renderedPart(*r, s, &err));  // flattening set by hand: still refused
}

// Overlap coverage measures the yarn around separate plies, so they gain no patch of their own
// (2026-10-31: the plies were skipped, read as 100 % uncovered and the car got a whole-car patch).
TEST(overlap_coverage_measures_plied_yarn) {
    if (!gpu::available(nullptr)) SKIP("no usable CUDA GPU");
    auto part = planePart(expandedD1());
    CHECK(part->knit != nullptr);
    if (!part->knit) return;
    part->chart = 0;
    Scene t;
    t.parts.push_back(std::move(*part));
    const ComponentGroup g = componentGroups(t).front();
    SceneSettings s;
    const GroupTube single = groupTube(t, g, s, ExecPolicy{});
    s.plies = 3;
    const GroupTube plied = groupTube(t, g, s, ExecPolicy{});
    CHECK(plied.knitted == 1 && plied.note.empty() && plied.tube.pos.size() == single.tube.pos.size());
    const CoverageReport c = groupCoverage(t, g, plied, ExecPolicy{});
    const CoverageReport c1 = groupCoverage(t, g, single, ExecPolicy{});
    std::printf("    plane: %.2f %% uncovered with plies, %.2f %% with the single tube\n", 100 * c.uncoveredFrac, 100 * c1.uncoveredFrac);
    CHECK(c.uncoveredFrac == c1.uncoveredFrac && c.uncoveredFrac < 1);  // this open gauge: not covered, but measured
    addGapPatches(t, s, KnitRunOptions{});
    const size_t pliedPatches = t.parts.size() - 1;
    t.parts.resize(1);
    s.plies = 1;
    addGapPatches(t, s, KnitRunOptions{});
    CHECK(pliedPatches == t.parts.size() - 1);
}

// Every load makes the scene's solids before any part, so every reader finds them (Scene::solids):
// a GLB scene knitted on its UV charts with and without its knitting setup (bears.glb with the teddy's
// pieces manifest cleared: no registered preset knits on charts since the man went), and the sewn pieces.
TEST(scene_solids_made_at_load) {
    std::string err;
    const SceneEntry* teddy = findScene("teddy");
    CHECK(teddy);
    if (!teddy) return;
    if (findProjectFile(teddy->path).empty()) SKIP("teddy asset not found");
    constexpr size_t kTeddyChartSolids = 18;  // measured 2026-11-10: the closed components of the two bears and their details
    SceneEntry charts = *teddy;
    charts.id = "teddy-charts";
    charts.piecesManifest.clear();
    Scene base, toKnit, pieces;
    CHECK(loadBaseParts(charts, viewerDefaults(charts), SourceOptions{}, base, &err));
    std::printf("    charts: %zu solids at the base load\n", base.solids ? base.solids->queries().size() : size_t(0));
    CHECK(base.solids && base.solids->queries().size() == kTeddyChartSolids);
    CHECK(loadPartsToKnit(charts, viewerDefaults(charts), SourceOptions{}, toKnit, &err).has_value());
    std::printf("    charts: %zu solids to knit\n", toKnit.solids ? toKnit.solids->queries().size() : size_t(0));
    CHECK(toKnit.solids && toKnit.solids->queries().size() == kTeddyChartSolids);
    SpiralSource src;
    CHECK(loadPiecesScene(*teddy, viewerDefaults(*teddy), src, pieces, &err));
    // The 18 pieces. The kept parts are bears.glb's nodes (open eyes, noses and mouths, a non-manifold bow):
    // none is a closed solid, so none hides fibres (the extras' closed bow was the 19th until 2026-11-11).
    CHECK(pieces.solids && pieces.solids->queries().size() == 18);
}

// The wool swatch is an ordinary GLB preset knitted on its own uv chart (no rebuild): one knitted part, the r2
// plane's chart named swatch_yarn_c0_ch0 (the name the re-knit's hand-made randomness was seeded with,
// stage_params.cpp), the needle one plain primitive, 1 mm per mm, and with yarn colour "object" the yarn is the GLB
// material's colour, the white scene's albedo (0.45, 0.44, 0.42), not material.yarn_white.color.
TEST(swatch_preset_loads_on_its_uv) {
    const SceneEntry* e = findScene("wool-swatch");
    CHECK(e != nullptr);
    if (!e) return;
    if (findProjectFile(e->path).empty()) SKIP("scenes/wool-swatch/wool_swatch.glb not found (the scene preparation tools make it)");
    const SceneSettings s = viewerDefaults(*e);
    Scene t;
    std::string err;
    CHECK(loadPartsToKnit(*e, s, SourceOptions{}, t, &err).has_value());
    int plain = 0;
    for (const SourcePrimitive& g : t.source.prims) plain += g.cls == SourcePrimitive::Class::Plain;
    std::printf("    %zu parts (%s), %d plain, scale %.6f mm per unit, height %.3f mm\n", t.parts.size(), t.parts.empty() ? "" : t.parts[0].name.c_str(), plain,
                t.source.scale, t.source.hi.z - t.source.lo.z);
    CHECK(t.parts.size() == 1 && t.parts[0].name == "swatch_yarn_c0_ch0" && t.parts[0].chart == 0 && plain == 1);
    CHECK(std::abs(t.source.scale - 1000) < 1e-3 && std::abs((t.source.hi.z - t.source.lo.z) - e->glb.sizeMm) < 1e-3);
    if (t.parts.empty()) return;
    const Vec3 c = partYarnColour(t, t.parts[0], s);
    std::printf("    yarn colour (linear) %.6f %.6f %.6f\n", c.x, c.y, c.z);
    CHECK(s.yarnColourSource == "object" && std::abs(c.x - 0.45) < 1e-6 && std::abs(c.y - 0.44) < 1e-6 && std::abs(c.z - 0.42) < 1e-6);
}

// Gap coverage reads Knit inputs only: the tube it measures is fixed, whatever the yarn's tessellation
// (strand.tube_*, Fibres rows outside the knit's key), so the tessellation cannot change which gap
// patches exist; an out-of-range tessellation no longer drops every patch (2026-11-06).
TEST(gap_coverage_ignores_tube_settings) {
    if (!gpu::available(nullptr)) SKIP("no usable CUDA GPU");
    auto part = planePart(expandedD1());
    CHECK(part->knit != nullptr);
    if (!part->knit) return;
    part->chart = 0;
    Scene t;
    t.parts.push_back(std::move(*part));
    const ComponentGroup g = componentGroups(t).front();
    const SceneSettings reduced;  // the default tessellation: 0.6 mm rings, 6 sides
    SceneSettings coarser = reduced;
    coarser.tubeRing = 1.2;
    coarser.tubeSides = 12;
    const GroupTube a = groupTube(t, g, reduced, ExecPolicy{}), b = groupTube(t, g, coarser, ExecPolicy{});
    CHECK(a.knitted == 1 && b.knitted == 1 && a.tube.pos == b.tube.pos && a.tube.tri == b.tube.tri);
    SceneSettings outOfRange = reduced;
    outOfRange.tubeSides = 2;  // tubeSettingsValid refuses it; the export refuses such settings, coverage does not read them
    const GroupTube c = groupTube(t, g, outOfRange, ExecPolicy{});
    CHECK(c.knitted == 1 && c.note.empty() && c.tube.pos == a.tube.pos);
}

std::string fileBytes(const std::string& p) {
    std::ifstream f(p, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(f), {});
}

// ------------------------------------------------------ cross-component ----

std::unique_ptr<ScenePart> linePart(const std::string& name, const Vec3& centre, double y, double D) {
    auto r = std::make_unique<ScenePart>();
    r->name = name;
    r->centre = centre;
    r->yarnD = D;
    r->knit = std::make_unique<PartKnit>();
    Polyline pl;
    for (int k = 0; k <= 20; ++k) pl.p.push_back(Vec3{0.5 * k, y, 0} - centre);  // world y, local coordinates
    r->knit->after.polylines.push_back(pl);
    return r;
}

// Two parts' centrelines (in their own re-centred frames) 1.05 mm apart in the world: closer
// than D + c = 1.1, reported FAIL with the pair; 1.2 mm apart: PASS. The minimum distance is
// the world distance (the part centres are applied).
TEST(cross_component_clearance_reported) {
    if (!gpu::available(nullptr)) SKIP("no usable CUDA GPU");
    SceneSettings s;
    std::string err;
    s.clearance = 0.1;
    for (double gap : {1.05, 1.2}) {
        std::vector<ScenePart> parts;
        parts.push_back(std::move(*linePart("a", {3, -7, 2}, 0, 1.0)));
        parts.push_back(std::move(*linePart("b", {-4, 9, 1}, gap, 1.0)));
        CrossComponentReport c;
        CHECK(crossComponentClearance(parts, s, c, &err));
        if (gap < 1.1) {
            CHECK_NEAR(c.minDistance, gap, 1e-12);
            CHECK(c.status == Status::Fail && c.pairs == 1 && c.closeSegments > 0 && c.lines.size() == 1);
            CHECK(c.lines[0].find("a vs b") != std::string::npos);
        } else {
            CHECK(c.status == Status::Pass && c.pairs == 0 && c.lines.empty());
            CHECK(c.minDistance >= 1.1);  // pairs beyond the largest threshold are not measured
        }
    }
    // A thicker part (fallback D) raises its pair's threshold: R_a + R_b + c = 0.5 + 0.6 + 0.1.
    std::vector<ScenePart> parts;
    parts.push_back(std::move(*linePart("a", {}, 0, 1.0)));
    parts.push_back(std::move(*linePart("b", {}, 1.15, 1.2)));
    CrossComponentReport thick;
    CHECK(crossComponentClearance(parts, s, thick, &err) && thick.status == Status::Fail);
}

// ------------------------------------------------------------------ units ----

// The caller's cancellation token holds in every chunk of a parallel range, on the
// executor's workers too, and a range without a token sees none: a cancelled job stops its
// parallel loops; Cancelled reaches the caller.
TEST(parallel_for_propagates_cancel_token) {
    ExecPolicy par;
    par.parallel = true;
    CancelToken tok;
    std::atomic<int> seen{0}, wrong{0};
    {
        ThreadCancelScope scope(&tok);
        parallelFor(256, 1, [&](int64_t, int64_t, int64_t) {
            std::this_thread::sleep_for(std::chrono::microseconds(200));
            (threadCancelToken() == &tok ? seen : wrong)++;
        }, par);
    }
    CHECK(seen == 256 && wrong == 0);
    std::atomic<int> none{0};
    parallelFor(256, 1, [&](int64_t, int64_t, int64_t) { none += threadCancelToken() == nullptr; }, par);
    CHECK(none == 256);
    // Cancelled mid-range: the checkpoints throw, the range stops early, the caller gets it.
    std::atomic<int> ran{0};
    bool threw = false;
    {
        ThreadCancelScope scope(&tok);
        try {
            parallelFor(4096, 1, [&](int64_t b, int64_t, int64_t) {
                if (b == 64) tok.cancel();
                cancelCheckpoint();
                ++ran;
                std::this_thread::sleep_for(std::chrono::microseconds(100));
            }, par);
        } catch (const Cancelled&) {
            threw = true;
        }
    }
    CHECK(threw && ran < 4096);
    std::printf("    cancelled after %d of 4096 chunks\n", ran.load());
}

// ------------------------------------------------------------ scene cache ----

BuiltScene sampleScene() {
    BuiltScene b;
    for (int k = 0; k < 3; ++k) {
        DisplayMesh d;
        d.name = "part " + std::to_string(k);
        for (int i = 0; i < 30; ++i) {
            d.pos.push_back(0.1f * float(i + k));
            d.nrm.push_back(0.5f - 0.01f * float(i));
        }
        d.idx = {0, 1, 2, 2, 3, 4};
        if (k == 1) d.uv = {0.f, 1.f, 0.25f, 0.5f};
        d.color[0] = 0.1f * float(k), d.color[1] = 0.2f, d.color[2] = 0.3f;
        d.kind = k == 2 ? DisplayMesh::Kind::Plain : DisplayMesh::Kind::Yarn;
        b.meshes.push_back(std::move(d));
    }
    PartRow r;
    r.name = "Object_1_02_Body_c2";
    r.knitted = true;
    r.verdict = Status::Fail;
    r.stitches = 1527;
    r.yarnD = 1.1;
    r.note = "D 1.10: gate a FAIL; ";
    r.lines = {"FAIL layout: 2 cells", "UNRESOLVED B2"};
    r.lo = {-1, -2, 0};
    r.hi = {3, 4, 5};
    b.parts.push_back(r);
    b.report = "parts: 22, knitted 5\n";
    b.lo = {-150, -80, 0};
    b.hi = {150, 80, 300};
    b.yarnLength = 98765.4;
    b.ms = 4321;
    return b;
}

// The binary form round-trips exactly; any truncation or a wrong header is rejected; a
// stored entry is read back (cached) and a corrupted file is a miss.
TEST(scene_cache_roundtrip) {
    namespace fs = std::filesystem;
    const BuiltScene b = sampleScene();
    const std::string bytes = serializeBuiltScene(b);
    BuiltScene c;
    CHECK(deserializeBuiltScene(bytes, c));
    CHECK(serializeBuiltScene(c) == bytes);
    CHECK(c.meshes.size() == 3 && c.meshes[1].uv == b.meshes[1].uv && c.meshes[2].kind == DisplayMesh::Kind::Plain);
    CHECK(c.parts.size() == 1 && c.parts[0].lines == b.parts[0].lines && c.parts[0].hi.z == 5 && c.parts[0].verdict == Status::Fail);
    for (size_t n : {size_t(0), size_t(7), bytes.size() / 2, bytes.size() - 1}) {
        BuiltScene t;
        CHECK(!deserializeBuiltScene(bytes.substr(0, n), t));
    }
    {
        BuiltScene t;
        std::string bad = bytes;
        bad[0] = 'X';
        CHECK(!deserializeBuiltScene(bad, t));
        CHECK(!deserializeBuiltScene(bytes + "x", t));
    }
    const std::string root = fs::temp_directory_path().string() + "/knit_s6_scenecache";
    fs::remove_all(root);
    SceneSettings s;
    std::string err;
    CHECK(!loadSceneCache(root, s, 42));
    CHECK(storeSceneCache(root, s, 42, b, &err));
    auto hit = loadSceneCache(root, s, 42);
    CHECK(hit && hit->cached && serializeBuiltScene(*hit) == bytes);
    CHECK(!loadSceneCache(root, s, 43));  // other assets
    for (const auto& d : fs::directory_iterator(root))
        if (fs::is_directory(d)) std::ofstream(d.path() / "scene.bin", std::ios::binary | std::ios::app) << "tail";
    CHECK(!loadSceneCache(root, s, 42));  // checksum mismatch: a miss
    fs::remove_all(root);
}

// The key follows what the display depends on: knit and surface settings (including the
// profile) and the assets; not materials or the UI.
TEST(scene_cache_key) {
    const SceneEntry* e = findScene("teddy");
    CHECK(e != nullptr);
    const SceneSettings a = viewerDefaults(*e);
    CHECK(isPreviewProfile(a));
    const std::string k = sceneCacheKey(a, 7);
    CHECK(k.find(" gen 11 ") != std::string::npos);  // 11: entries made before gap coverage had a fixed tube miss
    auto differs = [&](auto edit, uint64_t asset = 7) {
        SceneSettings t = a;
        edit(t);
        return sceneCacheKey(t, asset) != k;
    };
    CHECK(differs([](SceneSettings& t) { t.diameter = 1.08; }));
    CHECK(differs([](SceneSettings& t) { t.pattern = "rib1x1"; }));
    CHECK(differs([&](SceneSettings& t) { applyFullProfile(t, *e); }));
    CHECK(differs([](SceneSettings& t) { t.textureProfile = "none"; }));
    CHECK(differs([](SceneSettings&) {}, 8));
    CHECK(differs([](SceneSettings& t) { t.scene = "wool-swatch"; }));
    CHECK(!differs([](SceneSettings& t) { t.yarnWhite.color = {0.2, 0.3, 0.4}; }));
    CHECK(!differs([](SceneSettings& t) { t.timeCapS = 60; }));  // run control, never in a key
    CHECK(!differs([](SceneSettings& t) { t.fibreDensity = 7; }));
    // The display's tube follows the strand tessellation.
    CHECK(differs([](SceneSettings& t) { t.tubeSides = 16; }));
    SceneSettings f = a;
    applyFullProfile(f, *e);
    CHECK(!isPreviewProfile(f) && f.settle && f.repair && !f.diameterFallback.empty());
}

// Scenes compute on the GPU; the viewer picks the explicit solver (one diameter attempt in the full
// profile) when the knitting solver is usable, the CLI defaults (SceneSettings) stay legacy / full search.
// A machine without a device gets a note that knitting and fibres cannot run; only fibres on the
// visible sides keep the CPU generator. The preview has no diameter search and no budget.
TEST(viewer_defaults_compute) {
    const SceneEntry* e = findScene("teddy");
    CHECK(e != nullptr);
    const SceneSettings cli;
    CHECK(cli.diameterAttempts == 0);
    const ComputeCaps& caps = computeCaps();
    const SceneSettings v = viewerDefaults(*e);
    CHECK(isPreviewProfile(v) && v.diameterAttempts == 0);
    SceneSettings f = v;
    applyFullProfile(f, *e);
    CHECK(!isPreviewProfile(f) && !f.diameterFallback.empty() && f.diameterAttempts == 1);
    applyPreviewProfile(f);
    CHECK(isPreviewProfile(f) && f.diameterAttempts == 0 && sceneCacheKey(f, 7) == sceneCacheKey(v, 7));
    // Other machines (injected capabilities): no CUDA device; a device whose solver cannot run. What cannot run
    // is said before a job, with the reason; a Draft needs no relaxation, a Final does.
    ComputeCaps none;
    none.whyNoDevice = none.whyNoNcp = "test: no device";
    ComputeCaps noSolver;
    noSolver.device = true;
    noSolver.whyNoNcp = "test: solver refused";
    SceneSettings full = v;
    applyFullProfile(full, *e);
    const std::vector<std::string> n0 = computeNotes(v, none);
    CHECK(n0.size() == 1 && n0[0].find("test: no device") != std::string::npos && n0[0].find("cannot run") != std::string::npos);
    CHECK(computeNotes(v, noSolver).empty());
    const std::vector<std::string> n1 = computeNotes(full, noSolver);
    CHECK(n1.size() == 1 && n1[0].find("test: solver refused") != std::string::npos && n1[0].find("Final") != std::string::npos);
    ComputeCaps ready = noSolver;
    ready.ncp = true;
    CHECK(computeNotes(full, ready).empty());
    std::printf("    this machine: %s; knitting solver %s\n", caps.device ? caps.deviceName.c_str() : caps.whyNoDevice.c_str(),
                caps.ncp ? "available" : caps.whyNoNcp.c_str());
    // The scenes' own fibre defaults (owner, 2026-11-05): the teddy opens at density 135 (neither bear
    // reaches the Reduced export's 7 M segments) with Natural islands; a scene without its own keeps the new
    // objects' (an import's entry).
    auto islandsOf = [](const SceneSettings& x) {
        const FilePreset* p = matchingFilePreset(PresetGroup::FibreIslands, x);
        return p ? p->id : std::string("custom");
    };
    const SceneEntry* teddy = findScene("teddy");
    CHECK(teddy != nullptr);
    const SceneSettings t = viewerDefaults(*teddy);
    SceneEntry import;
    import.id = "import";
    import.kind = SceneEntry::Kind::Mesh;
    import.glb.knitAllMaterials = true;
    const SceneSettings i = viewerDefaults(import);
    std::printf("    teddy: density %g, islands %s; an import: density %g, islands %s\n", t.fibreDensity, islandsOf(t).c_str(),
                i.fibreDensity, islandsOf(i).c_str());
    CHECK(t.fibreDensity == 135 && islandsOf(t) == "natural");
    CHECK(i.fibreDensity == SceneSettings{}.fibreDensity && islandsOf(i) == "worn");
}

// A Knit job whose settings are in the cache is read, not computed: a synthetic entry
// stored under the man scene's key comes back from knitObject at once.
TEST(builder_cache_hit_skips_job) {
    namespace fs = std::filesystem;
    const SceneEntry* e = findScene("teddy");
    CHECK(e != nullptr);
    std::string err;
    const uint64_t asset = sceneAssetHash(*e, &err);
    if (!asset) SKIP("teddy asset not found");
    const std::string root = fs::temp_directory_path().string() + "/knit_s6_scenecache_hit";
    fs::remove_all(root);
    const SceneSettings s = viewerDefaults(*e);
    const BuiltScene b = sampleScene();
    CHECK(storeSceneCache(root, s, asset, b, &err));
    KnitRequest j;
    j.settings = std::make_shared<const SceneSettings>(s);
    j.cacheRoot = root;
    const auto t0 = std::chrono::steady_clock::now();
    auto r = knitObject(j, JobControl{9, nullptr}, &err);
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    CHECK(r && r->cached && r->revision == 9 && r->settings == s);
    CHECK(r->parts.size() == 1 && r->parts[0].stitches == 1527 && r->meshes.size() == 3);
    CHECK(sec < 5);  // the man scene's computation takes minutes
    std::printf("    cache hit in %.3f s\n", sec);
    fs::remove_all(root);
}

// Partial scenes: a job's partial publications are shown while it runs (flagged partial),
// with its progress; the final result replaces them; when a job is cancelled (and its
// successor fails before showing anything) its last partial scene stays shown.
TEST(builder_progressive_publish) {
    std::atomic<int> step{0}, hold{1};
    SceneBuilder b([&](const JobRequest& j, const JobControl& c, std::string* err) -> std::shared_ptr<BuiltScene> {
        if (jobSettings(j).scene == "fail") {
            *err = "synthetic failure";
            return nullptr;
        }
        for (int k = 1; k <= 3; ++k) {
            auto p = std::make_shared<BuiltScene>();
            p->revision = c.revision;
            p->settings = jobSettings(j);
            p->parts.resize(k);
            c.reporter->progress("knitting", k, 3);
            c.reporter->partial(p);
            step = 10 * c.revision + k;
            while (hold == k) {
                cancelCheckpoint();
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
        }
        auto f = std::make_shared<BuiltScene>();
        f->revision = c.revision;
        f->settings = jobSettings(j);
        f->parts.resize(5);
        return f;
    });
    auto waitStep = [&](int v) {
        while (step != v) std::this_thread::sleep_for(std::chrono::milliseconds(1));
    };
    CHECK(startKnit(b) == 1);
    waitStep(11);
    {
        const JobSnapshot j = b.snapshot();
        CHECK(j.latest && j.latest->partial && j.latest->parts.size() == 1);
        CHECK(j.progress.revision == 1 && j.progress.stage == "knitting" && j.progress.done == 1 && j.progress.total == 3);
    }
    hold = 2;
    waitStep(12);
    CHECK(b.snapshot().latest->partial && b.snapshot().latest->parts.size() == 2);
    hold = 4;
    waitIdle(b);
    {
        const JobSnapshot j = b.snapshot();
        CHECK(j.latest && !j.latest->partial && j.latest->parts.size() == 5 && j.published == 1);
        CHECK(j.progress.revision == -1);
    }
    hold = 1;
    CHECK(startKnit(b) == 2);
    waitStep(21);
    CHECK(b.snapshot().latest->revision == 2 && b.snapshot().latest->partial);
    CHECK(startKnit(b, "fail") == 3);  // cancels revision 2; revision 3 fails without publishing
    waitIdle(b);
    const JobSnapshot j = b.snapshot();
    CHECK(j.lastOutcome == JobOutcome::Failed);
    CHECK(j.latest->revision == 2 && j.latest->partial && j.latest->parts.size() == 1);
}

// The real pipeline on the man scene (preview profile, cold): the base meshes are published
// before any part is knitted, progress counts the parts, the final scene is complete; the
// computed scene is stored and the next request reads it.
// Time caps (seamless plan M0): the man knitted on the CPU with a 1.5 s part cap and a 4 s run
// cap. Parts past their deadline are UNRESOLVED "time cap" (not knitted, never cached), the parts
// after the run's deadline are "not started", and the run returns within a bounded overrun (the
// cap is cooperative: a part stops at its next checkpoint).
TEST(deadline_marks_part_unresolved_and_rest_not_started) {
    const SceneEntry* e = findScene("teddy");
    CHECK(e != nullptr);
    std::string err;
    if (!sceneAssetHash(*e, &err)) SKIP("teddy asset not found");
    SceneSettings s;
    s.scene = "teddy";
    applySceneDefaults(s, *e);
    KnitRunOptions ro;
    ro.source.selfIntersection = false;
    ro.capS = 4;
    ro.partCapS = 1.5;
    Scene t;
    const auto t0 = std::chrono::steady_clock::now();
    CHECK(knitScene(*e, s, ro, t, &err));
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    int timedOut = 0, notStarted = 0, knitted = 0;
    for (const ScenePart& r : t.parts) {
        timedOut += r.stop == KnitStop::TimeCap;
        notStarted += r.stop == KnitStop::NotStarted;
        knitted += r.knit != nullptr;
        if (r.stop == KnitStop::TimeCap) CHECK(!r.knit && r.note.rfind("time cap", 0) == 0);
        if (r.stop == KnitStop::NotStarted) CHECK(!r.knit);
    }
    std::printf("    %zu parts: %d knitted, %d time cap, %d not started; run %.1f s (loading included)\n", t.parts.size(), knitted, timedOut,
                notStarted, sec);
    CHECK(timedOut >= 1 && notStarted >= 1);
    CHECK(formatKnitReport(t, s, false).find("UNRESOLVED (time cap") != std::string::npos);
    CHECK(formatKnitReport(t, s, false).find("NOT STARTED (") != std::string::npos);
}

// The parts knit in scene order, one at a time: each part's start is followed by its end, before
// the next part starts, and `done` counts them.
TEST(knit_parts_in_scene_order) {
    const SceneEntry* e = findScene("teddy");
    CHECK(e != nullptr);
    std::string err;
    if (!sceneAssetHash(*e, &err)) SKIP("teddy asset not found");
    SceneSettings s;
    s.scene = "teddy";
    applySceneDefaults(s, *e);
    KnitRunOptions ro;
    ro.source.selfIntersection = false;
    ro.capS = 4;
    ro.partCapS = 1;
    struct Event { size_t part; PartEvent event; int done; };
    std::vector<Event> events;
    ro.onPart = [&](size_t i, PartEvent event, int done, int) { events.push_back({i, event, done}); };
    Scene t;
    CHECK(knitScene(*e, s, ro, t, &err));
    std::printf("    %zu events\n", events.size());
    CHECK(events.size() >= 4 && events.size() % 2 == 0);
    for (size_t k = 0; k + 1 < events.size(); k += 2) {
        CHECK(events[k].event == PartEvent::Started && events[k + 1].event == PartEvent::Finished);
        CHECK(events[k].part == events[k + 1].part && events[k + 1].done == events[k].done + 1);
        if (k > 0) CHECK(events[k].part == events[k - 1].part + 1);
    }
}

} // namespace

std::string fileBytesS6(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    return std::string((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
}

// The knit cache of the spiral layout (spiral_cache.h; knit_check --save-knit / --load-knit). A
// one-piece manifest (the chart fixture) as a scene of sewn pieces, knitted by knitScene: a run
// writes the cache, then runs read it. Same settings: identical curves, verdict, diameter and stitch
// count, and nothing knitted again (no knitting time). Changed surface settings (gate (d) on the
// flattened surface): the cache is used. A changed knitting setting, the piece file edited in place
// (manifest unchanged) or a missing cache: refused, the run fails.
TEST(knit_cache_pipeline) {
    if (!knit_test::useGpu()) SKIP("no GPU: scenes knit on the GPU");
    namespace fs = std::filesystem;
    const fs::path dir = fs::temp_directory_path() / "knit_s6_cache";
    fs::remove_all(dir);
    fs::create_directories(dir);
    const std::string obj = (dir / "piece_fixture.obj").string();
    fs::copy_file(std::string(KNIT_TEST_DATA) + "/chart_fixture.obj", obj);
    {
        std::ofstream m(dir / "manifest.txt");
        m << "piece piece_fixture.obj object=fixture island=0 material=fabric_white rings=5 segments=32 area=1 "
             "centroid=0,0,0 poleA=0,0,-100 poleB=0,0,100\n";
    }
    SceneEntry e;
    e.id = "fixture";
    e.projectName = "fixture";
    e.piecesManifest = (dir / "manifest.txt").string();
    e.gaugePiece = "piece_fixture.obj";
    e.gaugeStitches = 24;
    e.piecesW = 5.5;
    e.piecesH = 3.4;
    e.piecesD = 1.16;
    e.piecesDFallback = {1.12, 1.08, 1.04, 1.0};
    auto settings = [&]() {
        SceneSettings s;
        applySceneDefaults(s, e);
        s.settle = true;
        s.repair = false;
        s.surfaceGate = s.surfaceHalf = false;
        return s;
    };
    // One run: true with the knitted piece, false if the run failed (a refused cache).
    auto knitRun = [&](const SceneSettings& s, const std::string& save, const std::string& load, Scene& t) {
        KnitRunOptions o;
        o.saveKnitCache = save;
        o.loadKnitCache = load;
        std::string err;
        const bool ok = knitScene(e, s, o, t, &err);
        if (!ok) std::printf("    run failed: %s\n", err.c_str());
        return ok && t.parts.size() == 1 && t.parts[0].knit != nullptr;
    };
    auto same = [](const Scene& x, const Scene& y) {
        const PartKnit &p = *x.parts[0].knit, &q = *y.parts[0].knit;
        bool ok = p.verdict == q.verdict && p.t.pp == q.t.pp && p.t.prog.stitches.size() == q.t.prog.stitches.size() &&
                  p.t.curves.yarns.size() == q.t.curves.yarns.size() && x.parts[0].note == y.parts[0].note;
        for (size_t k = 0; ok && k < p.t.curves.yarns.size(); ++k) ok = p.t.curves.yarns[k].curve.ctrl == q.t.curves.yarns[k].curve.ctrl;
        return ok;
    };
    const std::string cache = (dir / "cache").string();
    Scene a, b;
    CHECK(knitRun(settings(), cache, "", a));
    // The cache sits in the producer's subfolder (gpu-f32/ncp).
    CHECK(a.parts.size() == 1 && a.parts[0].name == "piece_fixture" && a.parts[0].knit && a.parts[0].knit->msTile > 0);
    CHECK(fs::exists(dir / "cache" / "gpu-f32" / "ncp" / "piece_fixture.knitcache"));
    CHECK(knitRun(settings(), "", cache, b));
    CHECK(same(a, b) && b.parts[0].knit->msTile == 0);  // from the cache: not knitted again
    SceneSettings s = settings();  // surface settings changed (gate (d) on the flattened surface): the cache is used
    s.surfaceGate = true;
    s.flatten = true;
    CHECK(knitRun(s, "", cache, b));
    CHECK(same(a, b) && b.parts[0].knit->msTile == 0);
    s = settings();  // a knitting setting changed: refused
    s.settle = false;
    CHECK(!knitRun(s, "", cache, b));
    s = settings();
    s.wale = 6.0;
    CHECK(!knitRun(s, "", cache, b));
    CHECK(!knitRun(settings(), "", (dir / "missing").string(), b));  // no cache in the directory: refused, never knitted instead
    {
        // The piece file edited in place (one coordinate's last digit), manifest unchanged: refused.
        std::string text = fileBytesS6(obj);
        const size_t v = text.find("\nv ");
        CHECK(v != std::string::npos);
        const size_t eol = text.find('\n', v + 1);
        text.insert(eol, "1");
        std::ofstream f(obj, std::ios::binary);
        f << text;
    }
    CHECK(!knitRun(settings(), "", cache, b));
}

// knit_check --closeup N on the common export (ExportOptions::closeupPart; the article's teddy
// shots read the cameras named closeup_* and medium_*): scene.json holds "closeup_<part>" 40 mm and
// "medium_<part>" 180 mm (true metres) from one look-at point on the part, with a 20 degree field
// of view. A part index outside the scene is refused.
TEST(export_closeup_cameras) {
    namespace fs = std::filesystem;
    const std::string dir = fs::temp_directory_path().string() + "/knit_s6_closeup";
    const std::string path = writeTinyGltf(dir);
    SceneEntry e;
    e.id = "tiny";
    e.path = path;
    e.projectName = "tiny";
    e.glb.sizeRef = GlbSelection::SizeRef::Height;
    e.glb.sizeMm = 3;
    e.glb.knitAllMaterials = true;
    SceneSettings s;
    applySceneDefaults(s, e);
    s.workflow = "fur";  // the mesh itself is written: the cameras are under test
    Scene t;
    std::string err, log;
    CHECK(loadBaseParts(e, s, SourceOptions{}, t, &err));
    CHECK(!t.parts.empty());
    ExportOptions o;
    o.folder = dir + "/export";
    o.closeupPart = 0;
    fs::remove_all(o.folder);
    CHECK(exportScene({e, t, s}, o, &log, &err));
    const std::string json = fileBytesS6(o.folder + "/scene.json");
    // A camera's "position" and "look_at" (the first two vectors after its name).
    auto camera = [&json](const std::string& name, Vec3& eye, Vec3& at, double& fov) {
        const size_t n = json.find("\"name\": \"" + name + "\"");
        if (n == std::string::npos) return false;
        const size_t p = json.find("\"position\": [", n), l = json.find("\"look_at\": [", n), f = json.find("\"fov_deg\": ", n);
        if (p == std::string::npos || l == std::string::npos || f == std::string::npos) return false;
        // "[x, y, z]": strtod skips the separators' spaces; the commas are stepped over.
        auto vec3 = [&json](size_t at, Vec3& v) {
            char* end = nullptr;
            v.x = std::strtod(json.c_str() + at, &end);
            v.y = std::strtod(end + 1, &end);
            v.z = std::strtod(end + 1, &end);
        };
        vec3(p + 13, eye);
        vec3(l + 12, at);
        fov = std::strtod(json.c_str() + f + 11, nullptr);
        return true;
    };
    const std::string part = t.parts[0].name;
    Vec3 closeEye, closeAt, mediumEye, mediumAt;
    double closeFov = 0, mediumFov = 0;
    CHECK(camera("closeup_" + part, closeEye, closeAt, closeFov));
    CHECK(camera("medium_" + part, mediumEye, mediumAt, mediumFov));
    std::printf("    closeup_%s: eye %.5f %.5f %.5f, look at %.5f %.5f %.5f m\n", part.c_str(), closeEye.x, closeEye.y, closeEye.z,
                closeAt.x, closeAt.y, closeAt.z);
    CHECK_NEAR(norm(closeEye - closeAt), 0.040, 1e-6);
    CHECK_NEAR(norm(mediumEye - mediumAt), 0.180, 1e-6);
    CHECK_NEAR(norm(closeAt - mediumAt), 0.0, 1e-9);
    CHECK(closeFov == 20 && mediumFov == 20);
    // mid (115 mm, 20 degrees) and far (400 mm, 30 degrees) on the same axis: the renderer's bookmarks.
    Vec3 midEye, midAt, farEye, farAt;
    double midFov = 0, farFov = 0;
    CHECK(camera("mid", midEye, midAt, midFov) && camera("far", farEye, farAt, farFov));
    CHECK_NEAR(norm(midEye - midAt), 0.115, 1e-6);
    CHECK_NEAR(norm(farEye - farAt), 0.400, 1e-6);
    CHECK_NEAR(norm(midAt - closeAt), 0.0, 1e-9);
    CHECK_NEAR(norm(farAt - closeAt), 0.0, 1e-9);
    CHECK(norm(normalized(farEye - farAt) - normalized(closeEye - closeAt)) < 1e-4);  // positions are written to 1 um
    CHECK(midFov == 20 && farFov == 30);
    o.closeupPart = static_cast<int>(t.parts.size());
    err.clear();
    CHECK(!exportScene({e, t, s}, o, &log, &err) && err.find("close-up") != std::string::npos);
}

// recompute 2c: the boundary checks (B1 certified, B2 sampled) and a whole tiled part (gate (a),
// settle, boundary, gate (d)) with a parallel policy equal the serial run in every result field.
std::string boundarySnapshot(const BoundaryCheck& b) {
    char s[512];
    std::snprintf(s, sizeof s, "%s %s %a %d %d %ld %d %d %d %d %d %s", toString(b.b1), toString(b.b2), b.minDistance, b.b1Fail,
                  b.b1Unresolved, b.b2Samples, b.b2Unresolved[0], b.b2Unresolved[1], b.b2Unresolved[2], b.b2Unresolved[3],
                  b.b2Unresolved[4], b.detail.c_str());
    return s;
}

std::string knitSnapshot(const PartKnit& k) {
    std::string o = std::string(toString(k.verdict)) + " " + k.detail + "\n" + validationSnapshot(k.before) + validationSnapshot(k.after) +
                    boundarySnapshot(k.boundary) + "\n";
    char b[512];
    std::snprintf(b, sizeof b, "settle %d %d %a %s %d\n", k.settle.rounds, k.settle.freeStitches, k.settle.maxMove, toString(k.gateDStatus),
                  k.projectionUnresolved);
    o += b;
    for (const SurfaceCheck& g : k.gateD) {
        std::snprintf(b, sizeof b, "gate d %a %d %d %a %s %ld %ld %a %d %d %a %a %a %a %a %s\n", g.minFoldRatio, g.foldRings,
                      g.worstFoldYarn, g.worstFoldS, toString(g.fold), g.pointsTested, g.penetratingPoints,
                      g.maxDepth, g.worstYarnA, g.worstYarnB, g.worstSA, g.worstSB, g.worstPoint.x, g.worstPoint.y,
                      g.worstPoint.z, toString(g.collisions));
        o += b;
    }
    for (const YarnPath& y : k.t.curves.yarns)
        for (const Vec3& q : y.curve.ctrl) {
            std::snprintf(b, sizeof b, "%a %a %a\n", q.x, q.y, q.z);
            o += b;
        }
    return o;
}

TEST(tiled_parallel_identical) {
    const TriMesh m = gridMesh(40, 30, 1.0, [](double x, double y) { return Vec3{x, y, 0}; },
                               [](double x, double y) { return Vec2{x / 100, y / 100}; });
    const MeshSurface surf(m);
    auto line = [](double y, double z) {
        CurveSet cs;
        cs.kind = CurveKind::Surface;
        YarnPath yp;
        for (int k = 0; k < 8; ++k) yp.curve.ctrl.push_back({5.0 + 4 * k, y, z});
        yp.finalize();
        cs.yarns.push_back(yp);
        return cs;
    };
    SurfaceBinding none;
    CurveSet over = line(-0.6, 0.3);
    SurfaceBinding bind;
    bind.surface = &surf;
    const double bc[3] = {1.0 / 3, 1.0 / 3, 1.0 / 3};
    bind.state.assign(1, std::vector<ProjState>(over.yarns[0].curve.ctrl.size(), surf.stateAt(0, bc, surf.stateAt(0, bc, {}).foot)));
    bind.result.assign(1, std::vector<ProjResult>(over.yarns[0].curve.ctrl.size(), ProjResult::Owned));
    struct Case { CurveSet cs; const SurfaceBinding* b; double step; };
    const std::vector<Case> cases{{line(1.0, 0.3), &none, 0.1}, {line(std::sqrt(0.4 * 0.4 - 0.09), 0.3), &none, 0.1},
                                  {line(std::sqrt(0.4 * 0.4 - 0.09), 0.3), &none, 5.0}, {over, &bind, 0.1}};
    for (const Case& c : cases) {
        const std::string want = boundarySnapshot(knit_test::checkBoundary(c.cs, *c.b, surf, 0.5, {c.step}));
        for (int conc : {1, 2, 0}) {
            ExecPolicy par;
            par.parallel = true;
            par.maxConcurrency = conc;
            par.chunk = 3;
            CHECK(boundarySnapshot(knit_test::checkBoundary(c.cs, *c.b, surf, 0.5, {c.step, 1e-3, par})) == want);
        }
    }
    // A whole tiled part (the S6 cylinder, open rims: boundary checks; gate (d) both passes).
    const double R = 20, H = 48;
    const int na = 160;
    const TriMesh raw = gridMesh(na, 48, 1.0,
                                 [&](double x, double y) {
                                     const double a = 2 * kPi * x / na;
                                     return Vec3{R * std::cos(a), R * std::sin(a), y * H / 48};
                                 },
                                 [&](double x, double y) { return Vec2{2 * kPi * R * x / na / 100, y * H / 48 / 100}; });
    const TriMesh cyl = weldVertices(raw, 1e-9).mesh;
    MeshKnitParams p = tiledParams(expandedD1());  // CPU reference: the executor's determinism
    p.tile.pp.yarnDiameter = 1.3;  // thicker than the default: settle has defects to work on
    p.surfaceGate = p.surfaceHalf = true;
    std::string err;
    PartKnit serial;
    CHECK(knitMesh(cyl, p, knit_test::cpuNumerics(), serial, &err));
    const std::string want = knitSnapshot(serial);
    for (int conc : {2, 0}) {
        MeshKnitParams q = p;
        q.exec.parallel = true;
        q.exec.maxConcurrency = conc;
        knit_test::CpuNumerics parallel(q.exec);
        PartKnit k;
        CHECK(knitMesh(cyl, q, parallel, k, &err));
        CHECK(knitSnapshot(k) == want);
    }
    // Deferred gate (d) on the kept attempt (the tiled scene's loop) equals knitMesh with the gate.
    MeshKnitParams nogate = p;
    nogate.surfaceGate = false;
    PartKnit deferred;
    CHECK(knitMesh(cyl, nogate, knit_test::cpuNumerics(), deferred, &err));
    runSurfaceGate(deferred, p, knit_test::cpuNumerics());
    CHECK(knitSnapshot(deferred) == want);
    std::printf("  cylinder at D 1.3: %s, %zu gate d passes; identical serial, at 2 and H, and with gate (d) deferred\n",
                toString(serial.verdict), serial.gateD.size());
}

int main(int argc, char** argv) {
    using knit_test::registry;
    if (argc > 1 && std::strcmp(argv[1], "--list") == 0) {
        for (const knit_test::Case& c : registry()) std::printf("%s\n", c.name);
        return 0;
    }
    int failed = 0, ran = 0, skipped = 0;
    for (const knit_test::Case& c : registry()) {
        if (argc > 1 && std::strcmp(argv[1], c.name) != 0) continue;
        ++ran;
        try {
            c.fn();
            std::printf("[ OK ] %s\n", c.name);
        } catch (const knit_test::Skipped& s) {
            ++skipped;
            std::printf("[SKIP] %s: %s\n", c.name, s.why);
        } catch (const std::exception& e) {
            ++failed;
            std::printf("[FAIL] %s: %s\n", c.name, e.what());
        }
    }
    if (ran == 0) { std::printf("no test named %s\n", argc > 1 ? argv[1] : "?"); return 2; }
    if (skipped) std::printf("%d skipped\n", skipped);
    if (failed) return 1;
    return skipped == ran ? 77 : 0;  // 77: every case run was skipped (CTest SKIP_RETURN_CODE)
}

TEST(knit_cache_round_trip) {
    // A knitted and settled sphere piece: written, read and written again, the bytes are equal
    // (every stored field round-trips); the curves, their arc-length tables (queries bit-equal),
    // the rebuilt topology index and gate (a) on the read curves equal the original's. A
    // truncated file and a foreign file are refused.
    FixturePiece f = fixturePiece(25.0);
    PatternParams pp;
    PieceKnitParams kp;
    pp.waleSpacing = kp.g.w;
    pp.courseSpacing = kp.g.h;
    KnitCacheEntry e;
    CHECK(knitPiece({*f.surf, true}, {pp, kp}, knit_test::cpuNumerics(), e.knit));
    ValidationReport r;
    settle({e.knit.curves, pp, f.surf.get()}, SettleParams{}, knit_test::cpuNumerics(), &r);
    e.key = "generator = knit-cache 1\nsigma = 0x1p+0\n";
    e.pp = pp;
    e.poleAInside = true;
    e.attempts = {{1.16, Status::Fail}, {1.12, Status::Pass}};
    e.residue = false;
    e.meshDigest = meshDigest(f.mesh);
    e.validationDigest = fnv1a64(validationSnapshot(r));
    e.log = "     settle: 2 rounds\n";
    e.msKnit = 1234.5;
    const std::string d = scratchDir("knit_cache");
    std::string err;
    CHECK(writeKnitCache(d + "/a.knitcache", e, &err));
    KnitCacheEntry back;
    CHECK(readKnitCache(d + "/a.knitcache", back, &err));
    CHECK(writeKnitCache(d + "/b.knitcache", back, &err));
    const std::string bytes = fileBytes(d + "/a.knitcache");
    CHECK(bytes == fileBytes(d + "/b.knitcache"));
    CHECK(back.key == e.key && back.log == e.log && back.msKnit == e.msKnit && back.pp == e.pp && back.poleAInside && !back.poleBInside);
    CHECK(back.attempts.size() == 2 && back.attempts[1].D == 1.12 && back.attempts[1].verdict == Status::Pass);
    const CurveSet &a = e.knit.curves, &b = back.knit.curves;
    CHECK(a.yarns.size() == b.yarns.size() && a.topo.stitches.size() == b.topo.stitches.size());
    for (size_t y = 0; y < a.yarns.size(); ++y) {
        CHECK(a.yarns[y].curve.ctrl == b.yarns[y].curve.ctrl);
        CHECK(a.yarns[y].length() == b.yarns[y].length());
        for (int q = 0; q <= 200; ++q) {
            const double s = a.yarns[y].length() * q / 200, u = a.yarns[y].uEnd() * q / 200;
            CHECK(a.yarns[y].paramAtArcLength(s) == b.yarns[y].paramAtArcLength(s));
            CHECK(a.yarns[y].arcLengthAt(u) == b.yarns[y].arcLengthAt(u));
            CHECK(a.topo.stitchAt(int(y), u) == b.topo.stitchAt(int(y), u));
            const ArcLengthBounds ba = a.yarns[y].lengthBounds(0, u), bb = b.yarns[y].lengthBounds(0, u);
            CHECK(ba.lo == bb.lo && ba.hi == bb.hi);
        }
    }
    CHECK(fnv1a64(validationSnapshot(validate(b, back.pp, {.generatorChecks = false}))) == e.validationDigest);
    {
        std::ofstream t(d + "/cut.knitcache", std::ios::binary);
        t.write(bytes.data(), static_cast<std::streamsize>(bytes.size() / 2));
    }
    CHECK(!readKnitCache(d + "/cut.knitcache", back, &err) && err.find("truncated") != std::string::npos);
    {
        std::ofstream t(d + "/foreign.knitcache", std::ios::binary);
        t << "something else\n";
    }
    CHECK(!readKnitCache(d + "/foreign.knitcache", back, &err) && err.find("not a knit cache") != std::string::npos);
    CHECK(staleKeyField(e.key, e.key).empty());
    CHECK(staleKeyField(e.key, "generator = knit-cache 1\nsigma = 0x1.1p+0\n") == "sigma");
    CHECK(staleKeyField(e.key, "generator = knit-cache 1\n") == "sigma");
    CHECK(staleKeyField(e.key, e.key + "repair = 1\n") == "repair");
    std::printf("  %zu bytes, %zu stitches, %zu yarns: round trip exact\n", bytes.size(), a.topo.stitches.size(), a.yarns.size());
}

TEST(knit_cache_attempt_selection) {
    // The cached attempts decide another diameter list exactly as the pipeline's loop would:
    // the first passing D is kept; every D failing keeps the requested D (residue).
    KnitCacheEntry e;
    e.attempts = {{1.16, Status::Fail}, {1.12, Status::Fail}, {1.08, Status::Pass}};
    double miss = 0;
    CHECK(selectCachedAttempt(e, {1.16, 1.12, 1.08, 1.04}, 1.16, &miss));
    CHECK(selectCachedAttempt(e, {1.16, 1.08}, 1.16, &miss));            // 1.12 dropped: 1.08 still the first pass
    CHECK(selectCachedAttempt(e, {1.12, 1.08, 1.04}, 1.12, &miss));      // another requested D, same walk
    CHECK(!selectCachedAttempt(e, {1.16, 1.10, 1.08}, 1.16, &miss) && miss == 1.10);  // 1.10 never attempted
    CHECK(!selectCachedAttempt(e, {1.16}, 1.16, &miss) && miss == 1.16);  // kept would be the 1.16 FAIL (not stored)
    KnitCacheEntry res;  // residue: every D failed, the requested D appended and kept
    res.attempts = {{1.16, Status::Fail}, {1.12, Status::Fail}, {1.16, Status::Fail}};
    res.residue = true;
    CHECK(selectCachedAttempt(res, {1.16, 1.12}, 1.16, &miss));
    CHECK(!selectCachedAttempt(res, {1.16, 1.12, 1.08}, 1.16, &miss) && miss == 1.08);
    CHECK(!selectCachedAttempt(res, {1.12}, 1.12, &miss) && miss == 1.12);  // the kept result is at 1.16, not 1.12
    KnitCacheEntry single;
    single.attempts = {{1.16, Status::Fail}};
    CHECK(selectCachedAttempt(single, {1.16}, 1.16, &miss));             // no fallback: the only attempt is kept
    CHECK(!selectCachedAttempt(KnitCacheEntry{}, {1.16}, 1.16, &miss));
    // Knit diameters are the drawn D times k rounded to 0.001: the white teddy's frozen knit holds
    // 1.5 x 0.776 = 1.164, and Dense (1.86 x 0.626 = 1.16436) knits the same stitches, 0.36 um apart.
    KnitCacheEntry dense;
    dense.attempts = {{1.5 * 0.776, Status::Fail}};
    CHECK(selectCachedAttempt(dense, {1.86 * 0.626}, 1.86 * 0.626, &miss));
    CHECK(!selectCachedAttempt(dense, {1.166}, 1.166, &miss) && miss == 1.166);  // 2 um apart: another knit
}

TEST(ply_maps_periodic) {
    PlyGeometry g;
    PlyProfileParams pp;
    pp.aoDirections = 64;
    const double r = 0.56, rep = g.pitch;
    const PlyShadingProfile p = buildPlyShadingProfile(g, r, pp);
    // Analytic map functions: one period apart in u and v.
    for (int k = 0; k < 500; ++k) {
        const double u = std::fmod(0.6180339887 * k, 1.0), v = std::fmod(0.7548776662 * k + 0.1, 1.0);
        const PlyMapSample a = plyMapSample(p, rep, u, v), bu = plyMapSample(p, rep, u + 1, v), bv = plyMapSample(p, rep, u, v + 1);
        for (const PlyMapSample* b : {&bu, &bv}) {
            CHECK_NEAR(a.H, b->H, 1e-9);
            CHECK_NEAR(a.ao, b->ao, 1e-9);
            CHECK(norm(a.nTangent - b->nTangent) < 1e-9);
        }
        CHECK_NEAR(norm(a.nTangent), 1.0, 1e-12);
        CHECK(a.H >= -1e-12 && a.H <= 1 + 1e-12 && a.ao >= 0 && a.ao <= 1);
    }
    // The integral of dh/du and dh/dv over a period is 0 (the height is a scalar function).
    for (int k = 0; k < 8; ++k) {
        const double c = 0.123 * k;
        double iu = 0, iv = 0;
        const int N = 20000;
        for (int i = 0; i < N; ++i) {
            const double t = (i + 0.5) / N;
            const PlyMapSample su = plyMapSample(p, rep, t, c), sv = plyMapSample(p, rep, c, t);
            iu += -su.nTangent.x / su.nTangent.z / N;  // d rho / ds (from the normal), integrated over u
            iv += -sv.nTangent.y / sv.nTangent.z * sv.rho / N;
        }
        CHECK(std::fabs(iu) < 1e-6 && std::fabs(iv) < 1e-6);
    }
    // Texel grids: wrap-filtered (bilinear, wrap) sampling is continuous across the seams, and
    // 8-bit decoded normals are unit within the quantization bound sqrt(3) / 255.
    const int W = 256, H = 256;
    const PlyMaps maps = makePlyMaps(p, rep, W, H);
    double lip = 0;  // largest texel-to-texel step of H (a Lipschitz bound per texel)
    for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x) {
            const float h = maps.H[size_t(y) * W + x];
            lip = std::max({lip, double(std::fabs(h - maps.H[size_t(y) * W + (x + 1) % W])), double(std::fabs(h - maps.H[size_t((y + 1) % H) * W + x]))});
            float c[3];
            for (int i = 0; i < 3; ++i) c[i] = std::round((maps.normal[3 * (size_t(y) * W + x) + i] * 0.5f + 0.5f) * 255) / 255 * 2 - 1;
            CHECK(std::fabs(std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]) - 1) <= std::sqrt(3.0) / 255);
        }
    for (int k = 0; k < 64; ++k) {
        const double t = (k + 0.5) / 64, e = 1e-7;
        CHECK(std::fabs(sampleWrap(maps.H, W, H, 1 - e, t) - sampleWrap(maps.H, W, H, e, t)) <= lip * 2 * e * W + 1e-6);
        CHECK(std::fabs(sampleWrap(maps.H, W, H, t, 1 - e) - sampleWrap(maps.H, W, H, t, e)) <= lip * 2 * e * H + 1e-6);
        CHECK(std::fabs(sampleWrap(maps.normal, W, H, 1 - e, t, 3, 1) - sampleWrap(maps.normal, W, H, e, t, 3, 1)) < 1e-3);
    }
    // Texel layout: row 0 is the top (v = 1 - 0.5 / H).
    CHECK_NEAR(maps.H[3], plyMapSample(p, rep, 3.5 / W, 1 - 0.5 / H).H, 1e-6);
}

// An open sheet gets a backing on the side opposite its knitted +N side, at the offset.
TEST(filling_open_sheet_side) {
    const TriMesh m = gridMesh(10, 8, 1.0, [](double x, double y) { return Vec3{x, y, 0}; },
                               [](double x, double y) { return Vec2{x / 100, y / 100}; });
    const MeshCheck mc = checkMesh(m);
    CHECK(mc.components.size() == 1 && !mc.components[0].closed && mc.components[0].consistent);
    TriMesh mn = m;
    mn.computeNormals();
    CHECK(mn.vnormal[0].z > 0);  // the knitted side
    const Filling f = makeFilling(m, mc.components[0], 0.7);
    CHECK(f.kind == Filling::Kind::Backing);
    for (const Vec3& p : f.mesh.pos) CHECK_NEAR(p.z, -0.7, 1e-12);
    CHECK(f.inverted == 0 && f.selfIntersection == Status::Pass);
    CHECK(f.summary().find("backing") != std::string::npos);
    // Inconsistent orientation: no filling, with the reason.
    TriMesh bad = m;
    std::swap(bad.tri[0][1], bad.tri[0][2]);
    const MeshCheck bc = checkMesh(bad);
    const Filling n = makeFilling(bad, bc.components[0], 0.7);
    CHECK(n.kind == Filling::Kind::None && !n.reason.empty());
}

// A closed thin slab moved inward by more than its thickness: the moved surface turns inside out;
// the report says so (inverted triangles, volume sign). A small offset is clean.
TEST(filling_offset_inversion_reported) {
    TriMesh slab;
    slab.pos = {{0, 0, 0}, {10, 0, 0}, {10, 10, 0}, {0, 10, 0}, {0, 0, 1}, {10, 0, 1}, {10, 10, 1}, {0, 10, 1}};
    slab.tri = {{0, 2, 1}, {0, 3, 2}, {4, 5, 6}, {4, 6, 7}, {0, 1, 5}, {0, 5, 4},
                {1, 2, 6}, {1, 6, 5}, {2, 3, 7}, {2, 7, 6}, {3, 0, 4}, {3, 4, 7}};
    const MeshCheck mc = checkMesh(slab);
    CHECK(mc.components[0].closed && mc.components[0].signedVolume > 0);
    const Filling ok = makeFilling(slab, mc.components[0], 0.1);
    CHECK(ok.kind == Filling::Kind::Stuffing && ok.inverted == 0 && ok.selfIntersection == Status::Pass && ok.signedVolume > 0);
    const Filling inv = makeFilling(slab, mc.components[0], 2.0);
    CHECK(inv.kind == Filling::Kind::Stuffing);
    CHECK(inv.inverted > 0);
    CHECK(inv.signedVolume <= 0);
    CHECK(inv.summary().find("inside out") != std::string::npos);
    std::printf("    %s\n", inv.summary().c_str());
}

// The guided route's rows (knit.route, knit.guided_*): written and read back; a settings file written before them reads
// as the surface route with no piece, so every older file keeps meaning what it meant. (Defect: a missing row refused.)
TEST(scene_settings_guided_rows) {
    SceneSettings g;
    g.scene = "import";
    g.knitRoute = "guided";
    g.guidedPiece = 3;
    g.guidedMarkings = "D:/somewhere/arm.markings";
    g.guidedMarkingsDigest = "5fe3b4d1adf4441c";
    SceneSettings back;
    std::string err;
    CHECK(readSettings(writeSettings(g), back, &err));
    CHECK(back.knitRoute == "guided" && back.guidedPiece == 3 && back.guidedMarkings == g.guidedMarkings && back.guidedMarkingsDigest == g.guidedMarkingsDigest);
    CHECK(knitsGuidedPiece(back));
    std::string older;
    std::istringstream lines(writeSettings(g));
    for (std::string line; std::getline(lines, line);)
        if (line.rfind("knit.route", 0) != 0 && line.rfind("knit.guided_", 0) != 0) older += line + "\n";
    SceneSettings old;
    CHECK(readSettings(older, old, &err));
    CHECK(old.knitRoute == "surface" && old.guidedPiece == -1 && old.guidedMarkings.empty() && old.guidedMarkingsDigest.empty() && !knitsGuidedPiece(old));
    SceneSettings bad = g;
    CHECK(!setSetting(bad, "knit.route", "spiral", &err));  // a route that does not exist is refused
}
