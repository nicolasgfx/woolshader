// The .preset file format and its groups, read from presets/<folder>.
#include "settings/file_presets.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>

#include "io/files.h"
#include "settings/settings_mapping.h"

namespace knit {

namespace {

const char kVersion[] = "knit-preset v1";

// What a group's messages say.
struct PresetKind {
    const char* folder;
    const char* missingFolder;
    const char* what;         // "<key> is not a value of <what>"
    const char* missingNote;  // after "missing <keys>"
};

PresetKind kindOf(PresetGroup group) {
    switch (group) {
        case PresetGroup::Fibres: return {"fibres", "no presets/fibres folder (the presets are files there)", "a fibre preset", " (a preset sets every value it owns)"};
        case PresetGroup::Density: return {"density", "no presets/density folder (the presets are files there)", "a density preset", ""};
        case PresetGroup::StrandMesh: return {"strand_mesh", "no presets/strand_mesh folder (the levels are files there)", "a strand mesh level", ""};
        case PresetGroup::HandKnit: return {"hand_knit", "no presets/hand_knit folder (the presets are files there)", "a hand-knit preset", ""};
        case PresetGroup::FibreIslands: return {"fibre_islands", "no presets/fibre_islands folder (the presets are files there)", "a fibre island preset", ""};
        case PresetGroup::Strays: return {"strays", "no presets/strays folder (the presets are files there)", "a stray preset", " (a preset sets every value it owns)"};
        case PresetGroup::YarnFibres: return {"yarn_fibres", "no presets/yarn_fibres folder (the looks are files there)", "a yarn fibre look", " (a look sets every value it owns)"};
    }
    return {};
}

// The Density control's range (src/viewer/panel_fibres.cpp, Coverage).
constexpr double kMinDensity = 1, kMaxDensity = 1000;

std::vector<std::string> fibreKeys() {
    std::vector<std::string> k = {"fibres.random_fraction", "fibres.max_segments", "fibres.segment_length"};
    for (const char* v : {"length_min", "length_median", "length_max", "length_spread", "thickness", "thickness_sigma", "curvature",
                          "max_angle_deg"})
        k.push_back(std::string("fibres.directed.") + v);
    k.push_back("fibres.directed.drape");
    // The strays' length spread, rise and direction spread belong to the Strays group (presetKeys).
    for (const char* v : {"length_min", "length_median", "length_max", "thickness", "thickness_sigma", "curvature", "curl_radius",
                          "curl_length"})
        k.push_back(std::string("fibres.random.") + v);
    k.push_back("material.fibre_directed.model");
    k.push_back("material.fibre_directed.fibre_longitudinal_roughness");
    k.push_back("material.fibre_hair_colour");
    // The wool look a profile carries into the export (MAT-22, 2026-11-04): the fibres' late taper and
    // crimp, the ply relief's strength and the yarn's translucency. Selecting a profile sets every one,
    // so a profile without crimp or translucency clears what another left. The yarn colour stays the user's.
    for (const char* key : {"fibres.tail_taper", "fibres.tail_taper_min_diameters", "fibres.tail_taper_tip", "fibres.crimp_amplitude",
                            "fibres.crimp_wavelength", "strand.normal_strength", "material.yarn_white.translucency"})
        k.push_back(key);
    return k;
}

std::vector<PresetFile> readGroupFolder(PresetGroup group, std::vector<std::string>& errors) {
    const PresetKind kind = kindOf(group);
    const std::string folder = presetFolder(kind.folder);
    if (folder.empty()) {
        errors.push_back(kind.missingFolder);
        return {};
    }
    return readPresetFiles(folder, &errors);
}

// Applies a file's values: only the group's keys, and all of them.
bool setPresetValues(const PresetFile& f, PresetGroup group, SceneSettings& values, std::string* err) {
    const auto bad = [&](const std::string& why) {
        if (err) *err = why;
        return false;
    };
    const PresetKind kind = kindOf(group);
    const std::vector<std::string>& keys = presetKeys(group);
    std::set<std::string> set;
    for (const auto& [key, value] : f.values) {
        if (std::find(keys.begin(), keys.end(), key) == keys.end()) return bad(key + " is not a value of " + kind.what);
        std::string e;
        if (!setSetting(values, key, value, &e)) return bad(e);
        set.insert(key);
    }
    std::string missing;
    for (const std::string& k : keys)
        if (!set.count(k)) missing += (missing.empty() ? "" : ", ") + k;
    if (!missing.empty()) return bad("missing " + missing + kind.missingNote);
    return true;
}

// A yarn fibre look's ranges: the body's controls (the generator's own refusals need the knit's yarn radius:
// export.h checkYarnBody) and the maps' material.
bool yarnFibresValuesValid(const SceneSettings& v, std::string* err) {
    const SceneSettings::YarnBody& b = v.yarnBody;
    const bool ok = b.plies >= 1 && b.plies <= 3 && b.fibres >= 1 && b.fibreDiameterUm >= 10 && b.shellDepth > 0 && b.shellDepth <= 1 &&
                    b.fibreAngleDeg >= 0 && b.fibreAngleDeg <= 60 && b.plyAngleDeg >= 0 && b.plyAngleDeg <= 60 && b.irregularity >= 0 &&
                    b.irregularity <= 2 && b.mapCavity >= 0 && b.mapCavity <= 1 && b.mapRoughness >= 0.01 && b.mapRoughness <= 1 &&
                    b.mapSpecular >= 0 && b.mapSpecular <= 1 && b.mapAnisotropy >= 0.05 && b.mapAnisotropy <= 1;
    if (ok) return true;
    if (err) *err = "yarn_body: plies 1..3, fibres >= 1, fibre_diameter_um >= 10, shell_depth in (0, 1], angles 0..60, irregularity 0..2, "
                    "map_cavity 0..1, map_roughness 0.01..1, map_specular 0..1, map_anisotropy 0.05..1";
    return false;
}

// The values' own rules, beyond each key's type.
bool presetValuesValid(PresetGroup group, const SceneSettings& values, std::string* err) {
    switch (group) {
        case PresetGroup::Fibres: return fibreSettingsValid(values, err);
        case PresetGroup::Density:
            if (values.fibreDensity >= kMinDensity && values.fibreDensity <= kMaxDensity) return true;
            if (err) *err = "fibres.density = " + std::to_string(values.fibreDensity) + ": 1 to 1000 candidates per mm of yarn";
            return false;
        case PresetGroup::StrandMesh: return tubeSettingsValid(values, err);
        case PresetGroup::HandKnit:
            if (values.hand.columnRows > 0 && values.hand.spotLength > 0) return true;
            if (err) *err = "knit.hand.column_rows and knit.hand.spot_length must be above 0";
            return false;
        case PresetGroup::FibreIslands:
            if (values.islands.size > 0 && values.islands.fuzzy + values.islands.matted + values.islands.worn > 0) return true;
            if (err) *err = "fibres.islands.size and the sum of the three kinds' shares must be above 0";
            return false;
        case PresetGroup::Strays:
            if (!fibreSettingsValid(values, err)) return false;  // the direction spread and the rise
            if (values.random.lengthSpread >= 0) return true;
            if (err) *err = "fibres.random.length_spread must be 0 or more";
            return false;
        case PresetGroup::YarnFibres: return yarnFibresValuesValid(values, err);
    }
    return true;
}

struct LoadedPresets {
    std::vector<FilePreset> presets;
    std::vector<std::string> errors;
};

LoadedPresets loadPresets(PresetGroup group) {
    LoadedPresets r;
    for (const PresetFile& f : readGroupFolder(group, r.errors)) {
        FilePreset p;
        std::string err;
        if (presetFromFile(group, f, p, &err)) r.presets.push_back(std::move(p));
        else r.errors.push_back(f.path + ": " + err);
    }
    return r;
}

// One static per group, each in its own case: loading the Fibres group reads the Strays group (a
// wool preset's link is checked), which must not re-enter an initialisation in progress.
const LoadedPresets& loaded(PresetGroup group) {
    switch (group) {
        case PresetGroup::Fibres: {
            static const LoadedPresets fibres = loadPresets(PresetGroup::Fibres);
            return fibres;
        }
        case PresetGroup::Density: {
            static const LoadedPresets density = loadPresets(PresetGroup::Density);
            return density;
        }
        case PresetGroup::StrandMesh: {
            static const LoadedPresets strandMesh = loadPresets(PresetGroup::StrandMesh);
            return strandMesh;
        }
        case PresetGroup::HandKnit: {
            static const LoadedPresets handKnit = loadPresets(PresetGroup::HandKnit);
            return handKnit;
        }
        case PresetGroup::Strays: {
            static const LoadedPresets strays = loadPresets(PresetGroup::Strays);
            return strays;
        }
        case PresetGroup::YarnFibres: {
            static const LoadedPresets yarnFibres = loadPresets(PresetGroup::YarnFibres);
            return yarnFibres;
        }
        case PresetGroup::FibreIslands: break;
    }
    static const LoadedPresets islands = loadPresets(PresetGroup::FibreIslands);
    return islands;
}

} // namespace

std::string trimmed(const std::string& s) {
    const size_t a = s.find_first_not_of(" \t\r");
    if (a == std::string::npos) return "";
    return s.substr(a, s.find_last_not_of(" \t\r") - a + 1);
}

bool parseKeyValueText(const std::string& text, const std::string& versionLine, std::vector<KeyValueLine>& out, std::string* err) {
    std::stringstream in(text);
    std::string raw;
    bool version = false;
    for (int line = 1; std::getline(in, raw); ++line) {
        const std::string t = trimmed(raw);
        if (t.empty() || t[0] == '#') continue;
        if (!version) {
            if (t != versionLine) return fail(err, "line " + std::to_string(line) + ": expected \"" + versionLine + "\"");
            version = true;
            continue;
        }
        const size_t eq = t.find('=');
        if (eq == std::string::npos) return fail(err, "line " + std::to_string(line) + ": expected key = value");
        out.push_back({line, trimmed(t.substr(0, eq)), trimmed(t.substr(eq + 1))});
    }
    if (!version) return fail(err, "missing \"" + versionLine + "\"");
    return true;
}

bool parsePresetFile(const std::string& text, PresetFile& out, std::string* err) {
    out = PresetFile{};
    std::vector<KeyValueLine> lines;
    if (!parseKeyValueText(text, kVersion, lines, err)) return false;
    bool order = false;
    std::set<std::string> seen;
    for (const KeyValueLine& l : lines) {
        const std::string at = "line " + std::to_string(l.line) + ": ";
        if (l.key.empty()) return fail(err, at + "empty key");
        if (!seen.insert(l.key).second) return fail(err, at + l.key + " appears twice");
        if (l.key == "id") out.id = l.value;
        else if (l.key == "label") out.label = l.value;
        else if (l.key == "strays") out.strays = l.value;
        else if (l.key == "note") out.note = l.value;
        else if (l.key == "order") {
            char* end = nullptr;
            const long v = std::strtol(l.value.c_str(), &end, 10);
            if (l.value.empty() || *end) return fail(err, at + "order = " + l.value + ": a whole number");
            out.order = static_cast<int>(v);
            order = true;
        } else {
            out.values.push_back({l.key, l.value});
        }
    }
    if (out.id.empty() || out.id.find_first_of(" \t") != std::string::npos) return fail(err, "id: one word, required");
    if (out.label.empty()) return fail(err, "label: required");
    if (!order) return fail(err, "order: required");
    return true;
}

std::string presetFolder(const std::string& kind) { return findProjectFile("presets/" + kind); }

std::vector<PresetFile> readPresetFiles(const std::string& folder, std::vector<std::string>* errors) {
    namespace fs = std::filesystem;
    std::vector<PresetFile> out;
    std::error_code ec;
    if (folder.empty() || !fs::is_directory(folder, ec)) return out;
    std::vector<fs::path> files;
    for (const auto& e : fs::directory_iterator(folder, ec))
        if (e.is_regular_file() && e.path().extension() == ".preset") files.push_back(e.path());
    std::sort(files.begin(), files.end());
    std::set<std::string> ids;
    for (const fs::path& p : files) {
        std::ifstream f(p, std::ios::binary);
        std::stringstream ss;
        ss << f.rdbuf();
        PresetFile pf;
        std::string err;
        if (!parsePresetFile(ss.str(), pf, &err)) {
            if (errors) errors->push_back(p.filename().string() + ": " + err);
            continue;
        }
        if (!ids.insert(pf.id).second) {
            if (errors) errors->push_back(p.filename().string() + ": id " + pf.id + " is used by another file");
            continue;
        }
        pf.path = p.string();
        out.push_back(std::move(pf));
    }
    std::stable_sort(out.begin(), out.end(), [](const PresetFile& a, const PresetFile& b) {
        return a.order != b.order ? a.order < b.order : a.id < b.id;
    });
    return out;
}

// --- the groups

const std::vector<std::string>& presetKeys(PresetGroup group) {
    static const std::vector<std::string> fibres = fibreKeys(), density = {"fibres.density"},
                                          strandMesh = {"strand.tube_sides", "strand.tube_ring", "strand.ply_tube_sides", "strand.ply_tube_ring"},
                                          handKnit = {"knit.hand.unevenness", "knit.hand.column_shift", "knit.hand.column_rows", "knit.hand.loop_rate",
                                                      "knit.hand.loop_height", "knit.hand.wave_height", "knit.hand.spot_rate", "knit.hand.spot_change",
                                                      "knit.hand.spot_length"},
                                          islands = {"fibres.islands.rate", "fibres.islands.size", "fibres.islands.strength", "fibres.islands.fuzzy",
                                                     "fibres.islands.matted", "fibres.islands.worn"},
                                          strays = {"fibres.random.spread_deg", "fibres.random.max_angle_deg", "fibres.random.length_spread"},
                                          yarnFibres = {"yarn_body.plies", "yarn_body.fibres", "yarn_body.fibre_diameter_um", "yarn_body.shell_depth",
                                                        "yarn_body.fibre_angle_deg", "yarn_body.ply_angle_deg", "yarn_body.irregularity",
                                                        "yarn_body.map_cavity", "yarn_body.map_roughness", "yarn_body.map_specular",
                                                        "yarn_body.map_anisotropy", "yarn_body.map_material",
                                                        "material.yarn_white.fibre_longitudinal_roughness",
                                                        "material.yarn_white.fibre_azimuthal_roughness", "material.yarn_white.fibre_ior"};
    switch (group) {
        case PresetGroup::Fibres: return fibres;
        case PresetGroup::Density: return density;
        case PresetGroup::StrandMesh: return strandMesh;
        case PresetGroup::HandKnit: return handKnit;
        case PresetGroup::Strays: return strays;
        case PresetGroup::YarnFibres: return yarnFibres;
        case PresetGroup::FibreIslands: break;
    }
    return islands;
}

bool presetFromFile(PresetGroup group, const PresetFile& f, FilePreset& out, std::string* err) {
    SceneSettings values;
    if (!setPresetValues(f, group, values, err) || !presetValuesValid(group, values, err)) return false;
    // A wool preset comes with a stray preset (presets/strays); no other kind names one.
    if (group != PresetGroup::Fibres && !f.strays.empty()) return fail(err, "strays: only a fibre preset names a stray preset");
    if (group == PresetGroup::Fibres) {
        if (f.strays.empty()) return fail(err, "strays: required (the id of a stray preset, presets/strays)");
        if (!presetWithId(PresetGroup::Strays, f.strays)) return fail(err, "strays = " + f.strays + ": no such stray preset (presets/strays)");
    }
    out = FilePreset{f.id, f.label, f.note, f.path, values, f.strays};
    return true;
}

const std::vector<FilePreset>& presets(PresetGroup group) { return loaded(group).presets; }
const std::vector<std::string>& presetErrors(PresetGroup group) { return loaded(group).errors; }

const FilePreset* presetWithId(PresetGroup group, const std::string& id) {
    for (const FilePreset& p : presets(group))
        if (p.id == id) return &p;
    return nullptr;
}

const FilePreset* matchingFilePreset(PresetGroup group, const SceneSettings& s) {
    const std::vector<std::string>& keys = presetKeys(group);
    for (const FilePreset& p : presets(group)) {
        const bool same = std::all_of(keys.begin(), keys.end(), [&](const std::string& k) { return sameSetting(s, p.values, k); });
        if (same) return &p;
    }
    return nullptr;
}

void applyPreset(PresetGroup group, SceneSettings& s, const FilePreset& p) {
    if (group == PresetGroup::Fibres) s.fibrePreset = p.id;
    for (const std::string& k : presetKeys(group)) copySetting(s, p.values, k);
    if (group != PresetGroup::Fibres) return;
    if (const FilePreset* strays = presetWithId(PresetGroup::Strays, p.strays)) applyPreset(PresetGroup::Strays, s, *strays);
}

} // namespace knit
