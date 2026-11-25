// Frozen bytes of everything the file formats promise (part of knit_tests_s6). Each case builds its
// input by hand, needs no GPU, and compares an FNV-1a 64 of exact text or file bytes with a pinned
// value: settings text, stage keys, import identities, what the mesh readers report, every stored
// format, a package's scene.json, an OBJ and report.json. A changed value means stored files, cache keys or exports
// changed: never re-pin to make a case pass, find what moved. All pins of a case are compared at the
// end, so a deliberate re-pin sees every new value at once.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "harness.h"
#include "io/cache_paths.h"
#include "io/glb_load.h"
#include "io/image_io.h"
#include "import/import_store.h"
#include "scene/knit/knit_cache.h"
#include "import/import_analysis.h"
#include "import/prepare_import.h"
#include "io/publish.h"
#include "export/export.h"
#include "export/export_scene.h"
#include "settings/scene_settings.h"
#include "model/coordinates.h"
#include "model/hash.h"
#include "scene/knit/saved_knit.h"
#include "scene/fibres/fibre_files.h"
#include "scene/result_keys.h"
#include "settings/scene_registry.h"
#include "scene/scene.h"

using namespace knit;
namespace fs = std::filesystem;

namespace {

std::string freshDir(const std::string& name) {
    const fs::path p = fs::temp_directory_path() / "knit_io_pins" / name;
    std::error_code ec;
    fs::remove_all(p, ec);
    fs::create_directories(p, ec);
    return p.generic_string();
}

std::string fileBytes(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    return std::string((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
}

void writeBytes(const std::string& path, const std::string& bytes) {
    std::ofstream f(path, std::ios::binary);
    f.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
}

// The pinned hashes, each measured on the code before the IO simplification. The wrapped options and
// identity follow kWrapVersion (now 6, the one engine's): a re-pin needs that kind of reason.
const std::map<std::string, uint64_t>& pinnedHashes() {
    static const std::map<std::string, uint64_t> pinned = {
        {"settings default", 0x1e84e0147d07294bull},  // re-pinned 2026-11-25: the default scene is wool-swatch (with "teddy" the text hashes to the former pin); before: re-pinned 2026-11-22 (hardening/stitchmesh-shirt): the settings text holds knit.guided_component (the guided route's knitting component); built without that row, the text hashes to the former pin; before: re-pinned 2026-11-19 (main 2a4928e merged into the Fur workflow): the Fur workflow's text (workflow, pile.* with the islands and irregularities) plus main's four guided-route rows; without those rows the text hashes to the branch's pin (evidence fur-loops/m5/pins/merge_text_check.py); before: re-pinned 2026-11-19: the settings text holds the guided route's rows knit.route, knit.guided_piece, knit.guided_markings, knit.guided_markings_digest (feature/guided-stockinette-workflow); before: re-pinned 2026-11-14: the settings text holds fibres.directed.drape (the text without that line hashes to the former pin); before: re-pinned 2026-11-14: the settings text holds yarn_body.map_material (fabric | shiny, the baked maps' material); before: re-pinned 2026-11-12: the import has a single engine, so the settings text has no import.wrap_engine and no import.wrap_recipe line; before: re-pinned 2026-11-11: yarn_body.render + the five yarn_body.map_* rows, plies default 3; before: re-pinned 2026-11-11: the settings text holds strand.max_flattening; before: re-pinned 2026-11-11: built-in strays = the Loose stray preset (spread 45, rise to 40, length spread 0.5); before: re-pinned 2026-11-11: import.charts is retired (no line in the settings text); before: re-pinned 2026-11-10: the default scene is "teddy" (the man was removed); before that today: default spread_deg 30 and random max_angle_deg 30 (aligned strays)
        {"settings changed", 0x12def9c2b5934d9dull},  // re-pinned 2026-11-22 (hardening/stitchmesh-shirt): the settings text holds knit.guided_component (the guided route's knitting component); built without that row, the text hashes to the former pin; before: re-pinned 2026-11-19 (main 2a4928e merged into the Fur workflow): the Fur workflow (workflow, pile.* inputs, fur key version 2) on main's guided rows and solver identity v5; before: re-pinned 2026-11-19: the settings text holds the guided route's rows knit.route, knit.guided_piece, knit.guided_markings, knit.guided_markings_digest (feature/guided-stockinette-workflow); before: re-pinned 2026-11-14: the settings text holds fibres.directed.drape (the text without that line hashes to the former pin); before: re-pinned 2026-11-14: the settings text holds yarn_body.map_material (fabric | shiny, the baked maps' material); before: re-pinned 2026-11-12: the import has a single engine, so the settings text has no import.wrap_engine and no import.wrap_recipe line; before: re-pinned 2026-11-11: yarn_body.render + the five yarn_body.map_* rows, plies default 3; before: re-pinned 2026-11-11: the settings text holds strand.max_flattening; before: re-pinned 2026-11-11: built-in strays = the Loose stray preset (spread 45, rise to 40, length spread 0.5); before: re-pinned 2026-11-11: import.charts is retired (no line in the settings text); before: re-pinned 2026-11-10: the default scene is "teddy" (the man was removed); before that today: its changed values derive from the new defaults
        {"keys default", 0x96951b4798e122d5ull},  // re-pinned 2026-11-25: the default scene is wool-swatch (with "teddy" the text hashes to the former pin); before: re-pinned 2026-11-19 (main 2a4928e merged into the Fur workflow): the fur key is the Fur workflow's (version 2, pile.*); the stage texts are main's: with main's fur key they hash to main's pin (evidence fur-loops/m5/pins/merge_pin_check.py); before: re-pinned 2026-11-19: derived.ncp_solver is identity v5 (compact contact rows, fix/contact-memory-real-arm); with the v4 identity text the keys hash to the former pin (checked by a temporary swap); before: re-pinned 2026-11-11: built-in strays = the Loose stray preset (spread 45, rise to 40, length spread 0.5); before: re-pinned 2026-11-10: the default scene is "teddy" (the man was removed); before that today: the fibre key holds ' spread' (default 30) and the new rise
        {"keys changed", 0xc64cbc505f407f0eull},  // re-pinned 2026-11-19 (main 2a4928e merged into the Fur workflow): the Fur workflow (workflow, pile.* inputs, fur key version 2) on main's guided rows and solver identity v5; before: re-pinned 2026-11-19: derived.ncp_solver is identity v5 (compact contact rows, fix/contact-memory-real-arm); with the v4 identity text the keys hash to the former pin (checked by a temporary swap); before: re-pinned 2026-11-15: derived.hand_rest (the hand is on); with that line removed the text hashes to the former pin; before: re-pinned 2026-11-14: fibres.directed.drape changed (on) is a fibre input (with it off the text hashes to the former pin); before: re-pinned 2026-11-14: yarn_body.map_material changed is an export input (every field on); before: re-pinned 2026-11-11: yarn_body.render + the five yarn_body.map_* rows, plies default 3 (every field on: the body's rows count); before: re-pinned 2026-11-11: strand.max_flattening changed with strand.flatten on is a knit input; before: re-pinned 2026-11-11: built-in strays = the Loose stray preset (spread 45, rise to 40, length spread 0.5); before: re-pinned 2026-11-10: the default scene is "teddy" (the man was removed); before that today: the fibre key holds ' spread' (default 30) and the new rise
        {"keys wrapped import", 0xe415020eadee1babull},  // re-pinned 2026-11-19 (main 2a4928e merged into the Fur workflow): the fur key is the Fur workflow's (version 2, pile.*); the stage texts are main's: with main's fur key they hash to main's pin (evidence fur-loops/m5/pins/merge_pin_check.py); before: re-pinned 2026-11-19: derived.ncp_solver is identity v5 (compact contact rows, fix/contact-memory-real-arm); with the v4 identity text the keys hash to the former pin (checked by a temporary swap); before: re-pinned 2026-11-12: the import has a single engine, so import.wrap_engine and import.wrap_recipe are no inputs of a wrapped import; before: re-pinned 2026-11-11: built-in strays = the Loose stray preset (spread 45, rise to 40, length spread 0.5); before: re-pinned 2026-11-11: import.charts is retired, no longer an input of an import; before: re-pinned 2026-11-10: the fibre key holds ' spread' (default 30) and the new rise
        {"keys fur", 0x89090920deb8be37ull},  // re-pinned 2026-11-25: the default scene is wool-swatch (with "teddy" the text hashes to the former pin); before: re-pinned 2026-11-19 (main 2a4928e merged into the Fur workflow): the Fur workflow (workflow, pile.* inputs, fur key version 2) on main's guided rows and solver identity v5; before: re-pinned 2026-11-19: derived.ncp_solver is identity v5 (compact contact rows, fix/contact-memory-real-arm); with the v4 identity text the keys hash to the former pin (checked by a temporary swap); before: re-pinned 2026-11-11: built-in strays = the Loose stray preset (spread 45, rise to 40, length spread 0.5); before: re-pinned 2026-11-10: the default scene is "teddy" (the man was removed); before that today: the fibre key holds ' spread' (default 30) and the new rise
        {"keys features on", 0x1acc2da14895549eull},  // re-pinned 2026-11-25: the default scene is wool-swatch (with "teddy" the text hashes to the former pin); before: re-pinned 2026-11-19 (main 2a4928e merged into the Fur workflow): the fur key is the Fur workflow's (version 2, pile.*); the stage texts are main's: with main's fur key they hash to main's pin (evidence fur-loops/m5/pins/merge_pin_check.py); before: re-pinned 2026-11-19: derived.ncp_solver is identity v5 (compact contact rows, fix/contact-memory-real-arm); with the v4 identity text the keys hash to the former pin (checked by a temporary swap); before: re-pinned 2026-11-11: built-in strays = the Loose stray preset (spread 45, rise to 40, length spread 0.5); before: re-pinned 2026-11-10: the default scene is "teddy" (the man was removed); before that today: the fibre key holds ' spread' (default 30) and the new rise
        {"options default", 0x84e0c379b6d13e9full},
        {"options wrapped", 0x010e8a9bf49dcdaaull},  // re-pinned 2026-11-12: the import has a single engine, so a wrapped import writes no engine, engine build or recipe line and has wrap_version = 6 (the one engine's); before: the pin measured before
        {"options builtin", 0x3d05ed58b0cf823cull},  // re-pinned 2026-11-12: the import has a single engine, so the line wrap_engine = builtin is gone; before: re-pinned 2026-11-11: wrap_version = 6, a chart mostly outside the stretch bounds is split instead of slit; before: re-pinned 2026-11-11: wrap_version = 5, the flatten's solver converges and the surface file holds the flattened triangles, so stored version-4 imports are never reused; before: re-pinned 2026-11-11: the built-in engine has its own stages, wrap_version = 4 (the other engine stays 3); before: re-pinned 2026-11-11: charts are no longer an option (the case set charts = generate)
        {"identity default", 0x8bcd64301d8b6a46ull},  // re-pinned 2026-11-11: PNGs are deflate-compressed by stb_image_write (the fixture's pin.png / colour_pin.png bytes)
        {"identity wrapped", 0x11e5aad84a88840cull},  // re-pinned 2026-11-12: the import has a single engine, so the options text changes as options wrapped; before: re-pinned 2026-11-11: PNGs are deflate-compressed by stb_image_write (the fixture's pin.png / colour_pin.png bytes)
        {"read pin.obj", 0x5dd9e9c5a059cf8dull},
        {"analysis pin.obj", 0x2ff297a1db389040ull},  // re-pinned 2026-11-11: no uv status or chart issues; the features speak of the rebuild
        {"read ascii.ply", 0x86a6b8e8b04a466eull},
        {"analysis ascii.ply", 0x3071f3498794f60aull},  // re-pinned 2026-11-11: as pin.obj
        {"read little.ply", 0x2ca7d0e127355412ull},
        {"analysis little.ply", 0x3071f3498794f60aull},  // re-pinned 2026-11-11: as pin.obj
        {"read big.ply", 0x355ecdb1d99efba4ull},
        {"analysis big.ply", 0x3071f3498794f60aull},  // re-pinned 2026-11-11: as pin.obj
        {"refusals", 0x752ec306b9d41361ull},
        {"kmesh", 0x8fdd377f344f99eaull},
        {"triangle colours", 0x6671e450c04a4acfull},
        {"imported scene", 0x31239c9e99d2e876ull},
        {"knitcache", 0x769cc0e1b503adeeull},  // re-pinned 2026-11-19: knit-cache v2 (85f8a7f, purl patterns: PatternParams::ribLayerOffset, ribPairSpread, garterWave in the file); the pin had not followed that version; before: the pin measured before
        {"tiled knit", 0x5f7beef64a36e436ull},  // re-pinned 2026-11-19: the serialized parts carry PatternParams with the purl fields of knit-cache v2 (85f8a7f); before: the pin measured before
        {"fibrecache", 0x7bd2cc87eed1bc31ull},
        {"mesh digest", 0xc403a6b17aa4ef81ull},
        {"snapshot", 0x4600d8ac489d1cbcull},
        {"scene all", 0x90d02bd4833264b2ull},
        {"textured mesh", 0xe8336b1eef651232ull},
        {"ground mesh", 0x57eeb2c5284996d4ull},
        {"colour texture", 0x1652a091177d310cull},  // re-pinned 2026-11-11: PNGs are deflate-compressed by stb_image_write (the fixture's pin.png / colour_pin.png bytes)
        {"scene other", 0x455f44cdf84355e0ull},
        {"scene plain", 0x772ca9c9e016b3eeull},
        {"obj metres", 0xcc24afcf7705ffe5ull},
        {"mtl", 0x8398efeae61baf49ull},
        {"obj import y up", 0xf1f46e35465ddba1ull},
        {"obj import z up", 0x0d761692de50dbe9ull},
        {"report ready", 0x71281208061d6b43ull},  // re-pinned 2026-11-25: the default scene is wool-swatch (with "teddy" the text hashes to the former pin); before: re-pinned 2026-11-22 (hardening/stitchmesh-shirt): the settings text holds knit.guided_component (the guided route's knitting component); built without that row, the text hashes to the former pin; before: re-pinned 2026-11-19 (main 2a4928e merged into the Fur workflow): the Fur workflow's text (workflow, pile.* with the islands and irregularities) plus main's four guided-route rows; without those rows the text hashes to the branch's pin (evidence fur-loops/m5/pins/merge_text_check.py); before: re-pinned 2026-11-19: the settings text holds the guided route's rows knit.route, knit.guided_piece, knit.guided_markings, knit.guided_markings_digest (feature/guided-stockinette-workflow); before: re-pinned 2026-11-14: the settings text holds fibres.directed.drape (the report without that line hashes to the former pin); before: re-pinned 2026-11-14: the settings text holds yarn_body.map_material (fabric | shiny, the baked maps' material); before: re-pinned 2026-11-12: the import has a single engine, so the settings text has no import.wrap_engine and no import.wrap_recipe line; before: re-pinned 2026-11-11: yarn_body.render + the five yarn_body.map_* rows, plies default 3; before: re-pinned 2026-11-11: the settings text holds strand.max_flattening; before: re-pinned 2026-11-11: built-in strays = the Loose stray preset (spread 45, rise to 40, length spread 0.5); before: re-pinned 2026-11-11: import.charts is retired (no line in the settings text); before: re-pinned 2026-11-10: the default scene is "teddy" (the man was removed); before that today: its settings text has the new defaults
        {"report diagnostic", 0xeffe4d97780e6b58ull},  // re-pinned 2026-11-25: the default scene is wool-swatch (with "teddy" the text hashes to the former pin); before: re-pinned 2026-11-22 (hardening/stitchmesh-shirt): the settings text holds knit.guided_component (the guided route's knitting component); built without that row, the text hashes to the former pin; before: re-pinned 2026-11-19 (main 2a4928e merged into the Fur workflow): the Fur workflow's text (workflow, pile.* with the islands and irregularities) plus main's four guided-route rows; without those rows the text hashes to the branch's pin (evidence fur-loops/m5/pins/merge_text_check.py); before: re-pinned 2026-11-19: the settings text holds the guided route's rows knit.route, knit.guided_piece, knit.guided_markings, knit.guided_markings_digest (feature/guided-stockinette-workflow); before: re-pinned 2026-11-14: the settings text holds fibres.directed.drape (the report without that line hashes to the former pin); before: re-pinned 2026-11-14: the settings text holds yarn_body.map_material (fabric | shiny, the baked maps' material); before: re-pinned 2026-11-12: the import has a single engine, so the settings text has no import.wrap_engine and no import.wrap_recipe line; before: re-pinned 2026-11-11: yarn_body.render + the five yarn_body.map_* rows, plies default 3; before: re-pinned 2026-11-11: the settings text holds strand.max_flattening; before: re-pinned 2026-11-11: built-in strays = the Loose stray preset (spread 45, rise to 40, length spread 0.5); before: re-pinned 2026-11-11: import.charts is retired (no line in the settings text); before: re-pinned 2026-11-10: the default scene is "teddy" (the man was removed); before that today: the same
        {"plain scene json", 0x5b4f4cc94ddcc2f7ull},
        {"plain scene obj", 0xac40eaaf80a1f954ull},
        {"plain scene mtl", 0x86a70e092e65d1f4ull},
    };
    return pinned;
}

// Named values against their pins; every difference is printed before the case fails.
class Pins {
public:
    void check(const std::string& name, const std::string& bytes) {
        const uint64_t got = fnv1a64(bytes);
        std::printf("  pin %s 0x%016llxull\n", name.c_str(), static_cast<unsigned long long>(got));
        if (got != pinnedHashes().at(name)) wrong_.push_back(name);
    }
    void finish() const {
        for (const std::string& name : wrong_) std::printf("  CHANGED: %s\n", name.c_str());
        CHECK(wrong_.empty());
    }

private:
    std::vector<std::string> wrong_;
};

// A value of `key` that is not `value` and that the key accepts.
std::string changedValue(const std::string& key, const std::string& value, const std::vector<std::string>& choices) {
    if (!choices.empty()) return choices.back() != value ? choices.back() : choices.front();
    if (value == "true") return "false";
    if (value == "false") return "true";
    if (key.size() > 6 && key.compare(key.size() - 6, 6, ".model") == 0) return value == "fibre" ? "diffuse-sheen" : "fibre";
    if (value.find(',') != std::string::npos) return "0.25, 0.5, 0.75";
    char* end = nullptr;
    const double v = std::strtod(value.c_str(), &end);
    if (!value.empty() && *end == 0) return std::to_string(std::llround(v) + 1);
    return "x";
}

// Default settings with every key given another value, in table order.
SceneSettings everyFieldChanged() {
    SceneSettings s;
    std::map<std::string, std::string> defaults;
    std::stringstream in(writeSettings(s));
    for (std::string line; std::getline(in, line);) {
        const size_t eq = line.find(" = ");
        if (eq != std::string::npos) defaults[line.substr(0, eq)] = line.substr(eq + 3);
    }
    for (const SettingField& k : settingFields()) {
        std::string err;
        CHECK(setSetting(s, k.key, changedValue(k.key, defaults[k.key], k.choices), &err));
    }
    return s;
}

// Every stage's inputs, and the keys hashed from them.
std::string keysText(const SceneSettings& s) {
    std::string o;
    for (Stage stage : {Stage::Load, Stage::Knit, Stage::Surface, Stage::Fibres, Stage::Export, Stage::Display}) {
        o += std::string(toString(stage)) + ":\n" + stageInputs(s, stage);
        o += std::to_string(stageKey(s, stage, 1234)) + "\n";
    }
    return o + std::to_string(furKey(s, 1234)) + "\n";
}

const char kObjText[] =
    "# a square of two faces, one a quad\n"
    "mtllib pin.mtl\n"
    "o square\n"
    "v 0 0 0\nv 40 0 0\nv 40 40 0 0.5 0.25 0.125\nv 0 40 0\nv 80 0 0\nv 80 40 0\n"
    "vt 0 0\nvt 1 0\nvt 1 1\nvt 0 1\nvt 0.5 0.5\n"
    "vn 0 0 1\n"
    "usemtl dyed\n"
    "f 1/1 2/2 3/3 4/4\n"
    "usemtl plain\n"
    "f -5/5 -2/2 -1/1\n"
    "f 2/2 6/3 3/4\n"
    "unknown statement\n"
    "s off\n";

const char kMtlText[] =
    "newmtl dyed\n"
    "Kd 0.2 0.4 0.6\n"
    "Ns 10\n"
    "map_Kd -s 1 1 1 -bm 0.5 pin.png\n"
    "newmtl plain\n"
    "Kd 2 -1 0.5\n";

const char kPlyHeader[] =
    "format %s 1.0\n"
    "comment TextureFile pin.png\n"
    "obj_info ignored\n"
    "element vertex 4\n"
    "property float x\nproperty float y\nproperty float z\n"
    "property float nx\nproperty float ny\nproperty float nz\n"
    "property uchar red\nproperty uchar green\nproperty uchar blue\n"
    "property float s\nproperty float t\n"
    "property float quality\n"
    "element face 2\n"
    "property list uchar int vertex_indices\n"
    "property list uchar float texcoord\n"
    "property int flags\n"
    "element edge 0\n"
    "property int vertex1\n"
    "end_header\n";

std::string plyHeader(const char* format) {
    char b[1024];
    std::snprintf(b, sizeof b, kPlyHeader, format);
    return "ply\n" + std::string(b);
}

// Little-endian or big-endian bytes of the PLY the ascii text below holds.
class PlyBytes {
public:
    explicit PlyBytes(bool bigEndian) : bigEndian_(bigEndian) {}
    void f32(float v) { raw(&v, 4); }
    void i32(int v) { raw(&v, 4); }
    void u8(unsigned char v) { out_.push_back(static_cast<char>(v)); }
    const std::string& bytes() const { return out_; }

private:
    void raw(const void* p, int n) {
        const char* c = static_cast<const char*>(p);
        for (int k = 0; k < n; ++k) out_.push_back(c[bigEndian_ ? n - 1 - k : k]);
    }
    bool bigEndian_;
    std::string out_;
};

std::string plyBinary(bool bigEndian) {
    PlyBytes b(bigEndian);
    const float pos[4][3] = {{0, 0, 0}, {30, 0, 0}, {30, 30, 0}, {0, 30, 0}};
    const float uv[4][2] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    for (int v = 0; v < 4; ++v) {
        for (float c : pos[v]) b.f32(c);
        for (int k = 0; k < 3; ++k) b.f32(k == 2 ? 1.0f : 0.0f);
        for (int k = 0; k < 3; ++k) b.u8(static_cast<unsigned char>(10 * (v + 1) + k));
        for (float c : uv[v]) b.f32(c);
        b.f32(0.5f);
    }
    b.u8(3); b.i32(0); b.i32(1); b.i32(2);
    b.u8(6);
    for (float c : {0.0f, 0.0f, 1.0f, 0.0f, 1.0f, 1.0f}) b.f32(c);
    b.i32(7);
    b.u8(3); b.i32(0); b.i32(2); b.i32(3);
    b.u8(0);
    b.i32(8);
    return plyHeader(bigEndian ? "binary_big_endian" : "binary_little_endian") + b.bytes();
}

const char kPlyAsciiBody[] =
    "0 0 0 0 0 1 11 12 13 0 0 0.5\n"
    "30 0 0 0 0 1 21 22 23 1 0 0.5\n"
    "30 30 0 0 0 1 31 32 33 1 1 0.5\n"
    "0 30 0 0 0 1 41 42 43 0 1 0.5\n"
    "3 0 1 2 6 0 0 1 0 1 1 7\n"
    "3 0 2 3 0 8\n";

void appendVec(std::string& o, const Vec3& v) {
    char b[120];
    std::snprintf(b, sizeof b, "%a %a %a\n", v.x, v.y, v.z);
    o += b;
}

std::string uvText(const std::vector<std::array<Vec2, 3>>& uv) {
    std::string o;
    for (const auto& corners : uv)
        for (const Vec2& c : corners) appendVec(o, {c.x, c.y, 0});
    return o;
}

// All that a reader reports about a file, with the folder left out of the paths.
std::string meshText(const ImportedMesh& m) {
    std::string o = m.format + "\n" + std::to_string(m.fileBytes) + "\n";
    for (const Vec3& p : m.pos) appendVec(o, p);
    for (const Vec3& c : m.vcolor) appendVec(o, c);
    for (const auto& t : m.tri) o += std::to_string(t[0]) + " " + std::to_string(t[1]) + " " + std::to_string(t[2]) + "\n";
    o += uvText(m.uv);
    for (int t : m.triMaterial) o += std::to_string(t) + "\n";
    for (const ImportMaterial& mat : m.materials) {
        o += mat.name + " " + mat.textureRef + " " + std::to_string(mat.hasKd) + std::to_string(mat.textureMissing) + std::to_string(!mat.texture.empty()) + "\n";
        appendVec(o, mat.kd);
    }
    o += std::to_string(m.faces) + " " + std::to_string(m.polygons) + " " + std::to_string(m.maxPolygon) + " " + std::to_string(m.facesWithoutUv) + "\n";
    for (const std::string& s : m.ignored) o += s + "\n";
    for (const std::string& s : m.dependencies) o += s.substr(s.find_last_of("/\\") + 1) + "\n";
    for (const std::string& s : m.missingDependencies) o += s.substr(s.find_last_of("/\\") + 1) + "\n";
    return o;
}

std::string analysisText(const ImportAnalysis& a) {
    std::string o;
    for (const ImportIssue& i : a.issues)
        o += std::to_string(int(i.kind)) + " " + i.id + " " + std::to_string(i.count) + " " + i.text + " | " + i.action + "\n";
    appendVec(o, a.lo);
    appendVec(o, a.hi);
    char b[320];
    std::snprintf(b, sizeof b, "%a %d %d %d %d %d %lld %lld %lld %lld %a\n", a.area, a.components, a.boundaryLoops, a.boundaryEdges,
                  a.nonManifoldEdges, a.inconsistentEdges, a.degenerate, a.duplicateFaces, a.unreferenced, a.duplicateVertices, a.suggestedUnitMm);
    return o + b + a.colourSource + " | " + a.suggestionReason + "\n";
}

void writeImportFiles(const std::string& dir) {
    writeBytes(dir + "/pin.obj", kObjText);
    writeBytes(dir + "/pin.mtl", kMtlText);
    writePngRgb8(dir + "/pin.png", 2, 2, {255, 0, 0, 0, 255, 0, 0, 0, 255, 255, 255, 255});
}

CurveSet handCurves() {
    CurveSet cs;
    for (int id = 0; id < 2; ++id) {
        YarnPath y;
        y.id = id;
        y.course = id;
        y.name = "yarn " + std::to_string(id);
        y.reversed = id == 1;
        for (int k = 0; k < 8; ++k) y.curve.ctrl.push_back({double(k), 0.25 * id, 0.125 * k * k});
        y.stitches.push_back({0, id, 0.0, y.curve.uMax()});
        y.finalize();
        cs.yarns.push_back(std::move(y));
    }
    cs.wales = 2;
    cs.courses = 2;
    cs.waleSpacing = 5.5;
    cs.courseSpacing = 3.4;
    StitchRecord r;
    r.id = 0;
    r.yarn = 0;
    r.uBegin = 0.5;
    r.uEnd = 1.5;
    r.turn = 1;
    r.index = 2;
    r.origin = {1, 2, 3};
    r.U = {1, 0, 0};
    r.V = {0, 1, 0};
    r.N = {0, 0, 1};
    cs.topo.stitches.push_back(r);
    TerminalRecord t;
    t.yarn = 1;
    t.uBegin = 0.25;
    t.uEnd = 0.75;
    t.stitches = {0};
    t.segments = {{0.25, 0.5}};
    t.centre = {1, 1, 1};
    t.axis = {0, 0, 1};
    t.radius = 2;
    cs.topo.terminals.push_back(t);
    return cs;
}

LayoutNode handNode(double at) {
    LayoutNode n;
    n.p = {at, 1, 2};
    n.n = {0, 0, 1};
    n.c = {1, 0, 0};
    n.v = {0, 1, 0};
    n.w = 5.5;
    n.h = 3.4;
    n.wL = 5.25;
    n.wR = 5.75;
    n.a = at / 2;
    n.b = at / 4;
    n.f = 0.125;
    return n;
}

KnitCacheEntry handKnitEntry() {
    KnitCacheEntry e;
    e.key = "generator = knit-cache 1\nsigma = 0x1p+0\n";
    PieceKnit& k = e.knit;
    k.prog.turns.push_back({6, 6, 1, 0, 7});
    k.prog.stitches = handCurves().topo.stitches;
    k.prog.turnStart = {0, 6};
    k.prog.castRing = true;
    k.lay.nodes = {handNode(1), handNode(2)};
    k.lay.before = handNode(0);
    k.lay.after = handNode(3);
    k.lay.helical = true;
    k.lay.aSign = -1;
    k.lay.virtualNodes = {handNode(4)};
    k.lay.virtualChildOf = {0, -1};
    k.lay.rowPrev = {-1, 0};
    k.lay.rowNext = {1, -1};
    k.lay.waleHint = {{0, 1, 0}};
    ProjState s;
    s.tri = 3;
    s.b[0] = 0.25;
    s.b[1] = 0.25;
    s.b[2] = 0.5;
    s.foot = {1, 2, 3};
    s.n = {0, 0, 1};
    s.query = {1, 2, 4};
    k.lay.nodeState = {s};
    k.curves = handCurves();
    k.startAtA = false;
    k.dir = -1;
    k.counts = {6, 7};
    k.startArc = 1.5;
    k.endArc = 2.5;
    k.detail = "a detail";
    e.pp.yarnDiameter = 1.16;
    e.poleAInside = true;
    e.attempts = {{1.16, Status::Fail}, {1.12, Status::Pass}};
    e.residue = true;
    e.generatorChecks = true;
    e.meshDigest = 0x1234;
    e.validationDigest = 0x5678;
    e.log = "     settle: 2 rounds\n";
    e.msKnit = 1234.5;
    return e;
}

FibreSet handFibres() {
    FibreSet fs;
    for (int k = 0; k < 3; ++k) {
        Fibre f;
        f.candidate = 10 + k;
        f.yarnId = k;
        f.parentYarn = k;
        f.population = k == 2 ? FibrePopulation::Random : FibrePopulation::Directed;
        f.truncated = k == 1;
        f.s = 0.5 * k;
        f.key = 0.25 * k;
        f.theta = 0.125f;
        f.anchor = {1, 2, 3};
        f.normal = {0, 0, 1};
        f.emergeL = 0.25f;
        f.radius0 = 0.05f;
        f.radius1 = 0.02f;
        f.targetLength = 1.5f;
        f.attach = 0.25f;
        f.attachArc = 0.5f;
        f.psi = 1;
        f.lift = 2;
        f.smoothLength = 3;
        f.emergeAngle = 4;
        f.helixActual = 5;
        f.maxDev = 6;
        f.minTurn = 7;
        f.maxTurn = 8;
        fs.addFibre(f, std::vector<Vec3>{{0, 0, 0}, {0.5, 0, 0.5}, {1, 0, 1.5}, {1.5, double(k), 3}});
    }
    fs.rep.candidates = 9;
    fs.rep.accepted = 3;
    fs.rep.totalLength = 4.5;
    fs.rep.ms = 12.5;
    fs.rep.pop.resize(2);
    fs.rep.pop[0].accepted = 2;
    fs.rep.pop[1].accepted = 1;
    return fs;
}

// A scene.json with every kind of thing the file can hold.
Placement placedUnder(const std::string& object, const std::string& group) { return {object, object + " name \"q\"", group}; }

void addHandParts(ExportScene& ex, std::string& err) {
    std::error_code ec;
    fs::create_directories(ex.folder() + "/meshes", ec);
    fs::create_directories(ex.folder() + "/strands", ec);
    TriMesh m;
    m.pos = {{-0.02, -0.02, -0.003}, {0.03, -0.02, -0.003}, {0.03, 0.03, -0.003}, {-0.02, 0.03, -0.003}};
    m.tri = {{0, 1, 2}, {0, 2, 3}};
    const std::vector<std::array<Vec2, 3>> corners = {{Vec2{0, 0}, Vec2{1, 0}, Vec2{1, 1}}, {Vec2{0, 0}, Vec2{1, 1}, Vec2{0.5, 1}}};
    CHECK(ex.addMesh("plate", "wool_white", m, &err));
    CHECK(ex.addTexturedMesh("textured", "tex_plate", m, corners, &err));
    const std::optional<std::string> texture = ex.addColourTexture("pin", ColourImage{"", 1, 1, {10, 20, 30, 255}}, &err);
    CHECK(texture);
    ex.baseMaterials.push_back({"tex_plate", {}});
    ex.albedoTextures["tex_plate"] = texture.value_or("");
    MaterialSpec coated;
    coated.sheen = 0;
    coated.clearcoat = true;
    ex.baseMaterials.push_back({"coated", coated});
    CHECK(ex.addMesh("coated_part", "coated", m, &err));
    CHECK(ex.addMesh("textured_yarn", ExportScene::yarnMaterial("wool_white"), m, &err));
    StrandSetF s;
    s.add({0, -0.01f, 0, 0.004f, -0.01f, 0}, {2e-4f, 2e-4f});
    CHECK(ex.addStrands("fibres_a", ex.orientedFibreMaterial("wool_white"), s, &err));
    CHECK(ex.addStrands("fibres_b", ex.randomFibreMaterial("wool_brown"), s, &err));
    CHECK(ex.addStrands("empty", "wool_white", StrandSetF{}, &err));
    ex.place("plate", placedUnder("bear", ""));
    ex.place("textured", placedUnder("bear", ""));
    ex.place("fibres_a", placedUnder("bear", ExportScene::kOrientedGroup));
    ex.place("fibres_b", placedUnder("hat", ExportScene::kRandomGroup));
    ex.extraCameras.push_back({"zoom", {0, 0, 1}, {0, 0, 0}, 25});
    ex.extraCameras.push_back({"another", {1, 0, 1}, {0, 0, 0}, 12.5});
}

// The wool look with everything on.
void setAllShading(ExportScene& ex) {
    ex.projectName = "pin \"project\"";
    ex.fibreLook.hairMaterials = true;
    ex.fibreLook.hair.model = MaterialSpec::Model::Fibre;
    ex.fibreLook.randomInheritsYarn = true;
    ex.plyShading.on = true;
    ex.plyShading.normal = "textures/n.png";
    ex.plyShading.cavity = "textures/c.png";
    ex.skyIntensity = 0.5;
}

// What the OBJ's colours and names depend on: hair and surface colours at their two precisions, a
// ply-mapped yarn whose spec is a fibre model, a colour texture, a duplicate id, an id no entry has,
// spaces in names, and a mesh named like the ground, which the OBJ leaves out.
void addObjParts(ExportScene& ex, std::string& err) {
    std::error_code ec;
    fs::create_directories(ex.folder() + "/meshes", ec);
    fs::create_directories(ex.folder() + "/strands", ec);
    ex.fibreLook.hairMaterials = true;
    ex.fibreLook.hair.model = MaterialSpec::Model::Fibre;
    ex.plyShading.on = true;
    ex.plyShading.cavity = "textures/c.png";
    MaterialSpec fibreYarn;
    fibreYarn.model = MaterialSpec::Model::Fibre;
    fibreYarn.color = {0.12344996, 0.4567891, 0.7654321};
    ex.baseMaterials[1] = {"wool_brown", fibreYarn};
    MaterialSpec first, second, spaced;
    first.color = {0.111, 0.222, 0.333};
    second.color = {0.444, 0.555, 0.666};
    spaced.color = {0.1239, 0.45649, 0.7};
    ex.baseMaterials.push_back({"dup", first});
    ex.baseMaterials.push_back({"dup", second});
    ex.baseMaterials.push_back({"two words", spaced});
    ex.baseMaterials.push_back({"tex_plate", {}});
    const std::optional<std::string> texture = ex.addColourTexture("pin", ColourImage{"", 1, 1, {10, 20, 30, 255}}, &err);
    CHECK(texture);
    ex.albedoTextures["tex_plate"] = texture.value_or("");
    TriMesh m;
    m.pos = {{-0.02, -0.02, -0.003}, {0.03, -0.02, -0.003}, {0.03, 0.03, -0.003}, {-0.02, 0.03, -0.003}};
    m.tri = {{0, 1, 2}, {0, 2, 3}};
    const std::vector<std::array<Vec2, 3>> corners = {{Vec2{0, 0}, Vec2{1, 0}, Vec2{1, 1}}, {Vec2{0, 0}, Vec2{1, 1}, Vec2{0.5, 1}}};
    CHECK(ex.addMesh("plate", "wool_white", m, &err));
    CHECK(ex.addMesh("fibre model", "wool_brown", m, &err));
    CHECK(ex.addMesh("mapped", ExportScene::yarnMaterial("wool_brown"), m, &err));
    CHECK(ex.addTexturedMesh("textured", "tex_plate", m, corners, &err));
    CHECK(ex.addMesh("twice", "dup", m, &err));
    CHECK(ex.addMesh("spaced", "two words", m, &err));
    CHECK(ex.addMesh("lost", "nowhere", m, &err));
    CHECK(ex.addMesh(ExportScene::kGround, "black", m, &err));
    StrandSetF a, b;
    a.add({0, -0.01f, 0, 0.004f, -0.01f, 0}, {2e-4f, 2e-4f});
    b.add({0, 0, 0, 0.001f, 0, 0, 0.002f, 0.001f, 0}, {1e-4f, 1.2e-4f, 1.4e-4f});
    b.add({0, 0, 0.001f, 0.001f, 0, 0.001f}, {5e-5f, 6e-5f});
    CHECK(ex.addStrands("fibres_a", ex.orientedFibreMaterial("wool_white"), a, &err));
    CHECK(ex.addStrands("fibres_b", ex.randomFibreMaterial("wool_brown"), b, &err));
    ex.sortParts();
}

// The package's OBJ, as an export writes it. The texture goes with it, under the same path.
void writeObjPin(const ExportScene& ex, const std::string& out, std::string& err) {
    std::string log;
    CHECK(writeObjFile(ex, {out, "pin"}, log, &err));
    CHECK(fileBytes(out + "/textures/colour_pin.png") == fileBytes(ex.folder() + "/textures/colour_pin.png"));
}

// A scene that exports without a GPU: one plain primitive and one part that could not be knitted.
void fillPlainScene(SceneEntry& e, Scene& t) {
    e.id = "pin";
    e.projectName = "pin_project";
    TriMesh m;
    m.pos = {{0, 0, 0}, {10, 0, 0}, {10, 10, 0}, {0, 10, 0}};
    m.tri = {{0, 1, 2}, {0, 2, 3}};
    SourcePrimitive plain;
    plain.cls = SourcePrimitive::Class::Plain;
    plain.mesh = "plate";
    plain.material = "metal";
    plain.nodeIndex = 0;
    plain.baseColor = {0.2, 0.3, 0.4};
    plain.welded = m;
    SourcePrimitive body = plain;
    body.cls = SourcePrimitive::Class::Knit;
    body.mesh = "body";
    body.material = "wool";
    body.nodeIndex = 1;
    body.check.components.emplace_back();
    body.check.components[0].triangles = {0, 1};
    t.source.prims = {plain, body};
    ScenePart part;
    part.prim = 1;
    part.component = 0;
    part.name = "body_wool_c0";
    part.mesh = std::make_shared<TriMesh>(m);
    part.note = "too small for the gauge";
    t.parts.push_back(std::move(part));
}

} // namespace

// A corrupt file must never allocate beyond its own size, and a read past the end gives zero.
TEST(io_byte_reader_refuses) {
    ByteWriter w;
    w.u8(200);
    w.i32(-5);
    w.u32(7);
    w.i64(-6);
    w.u64(8);
    w.f64(0.1);
    w.vec2({1, 2});
    w.vec3({3, 4, 5});
    w.string32("short");
    w.string64("long");
    w.doubles({0.5, 0.25});
    w.ints({-1, 2});
    ByteReader r(w.written(), 0);
    CHECK(r.u8() == 200 && r.i32() == -5 && r.u32() == 7 && r.i64() == -6 && r.u64() == 8 && r.f64() == 0.1);
    const Vec2 v2 = r.vec2();
    const Vec3 v3 = r.vec3();
    CHECK(v2.y == 2 && v3.z == 5 && r.string32() == "short" && r.string64() == "long");
    CHECK(r.doubles() == (std::vector<double>{0.5, 0.25}) && r.ints() == (std::vector<int>{-1, 2}));
    CHECK(r.ok() && r.at() == w.written().size());
    CHECK(r.i64() == 0 && !r.ok());  // past the end: zero, and it stays refused
    CHECK(r.string64().empty() && !r.ok());

    ByteWriter huge;
    huge.i64(int64_t(1) << 40);
    huge.f64(0);
    ByteReader tooMany(huge.written(), 0);
    CHECK(tooMany.count(8) == 0 && !tooMany.ok());
    CHECK(tooMany.doubles().empty());  // a refused reader allocates nothing either
    ByteWriter negative;
    negative.i64(-1);
    ByteReader below(negative.written(), 0);
    CHECK(below.count(1) == 0 && !below.ok());
    ByteWriter length;
    length.u32(0xffffffffu);
    ByteReader longString(length.written(), 0);
    CHECK(longString.string32().empty() && !longString.ok());
    ByteWriter exact;
    exact.i64(2);
    exact.f64(1);
    exact.f64(2);
    ByteReader fits(exact.written(), 0);
    CHECK(fits.count(8) == 2 && fits.ok());
}

TEST(io_pin_settings_text) {
    Pins p;
    const SceneSettings changed = everyFieldChanged();
    p.check("settings default", writeSettings(SceneSettings{}));
    p.check("settings changed", writeSettings(changed));
    SceneSettings back;
    std::string err;
    p.finish();
    CHECK(readSettings(writeSettings(changed), back, &err));
    CHECK(writeSettings(back) == writeSettings(changed));
}

TEST(io_pin_stage_keys) {
    Pins p;
    p.check("keys default", keysText(SceneSettings{}));
    p.check("keys changed", keysText(everyFieldChanged()));
    SceneSettings wrapped;
    wrapped.scene = "import";
    wrapped.importSurface = "wrapped";
    p.check("keys wrapped import", keysText(wrapped));
    SceneSettings fur;
    fur.workflow = "fur";
    p.check("keys fur", keysText(fur));
    SceneSettings overlap;
    overlap.coverage = "overlap";
    p.check("keys features on", keysText(overlap));
    p.finish();
}

TEST(io_pin_import_identity) {
    const std::string d = freshDir("identity");
    writeImportFiles(d);
    ImportedMesh m;
    std::string err;
    CHECK(readMeshFile(d + "/pin.obj", ImportLimits{}, m, &err));
    ImportOptions o;
    ImportOptions wrapped = o;
    wrapped.surface = ImportOptions::Surface::Wrapped;
    ImportOptions builtin = wrapped;
    builtin.up = UpAxis::Z;
    builtin.unitMm = 1.0 / 3;     // every digit shows
    builtin.wrapStitchMm = 2.0 / 3;
    Pins p;
    p.check("options default", importOptionsText(o));
    p.check("options wrapped", importOptionsText(wrapped));
    p.check("options builtin", importOptionsText(builtin));
    p.check("identity default", std::to_string(importIdentity(m, o)));
    p.check("identity wrapped", std::to_string(importIdentity(m, wrapped)));
    p.finish();
}

TEST(io_pin_mesh_read) {
    const std::string d = freshDir("mesh_read");
    writeImportFiles(d);
    writeBytes(d + "/ascii.ply", plyHeader("ascii") + kPlyAsciiBody);
    writeBytes(d + "/little.ply", plyBinary(false));
    writeBytes(d + "/big.ply", plyBinary(true));
    Pins p;
    std::string err;
    for (const char* file : {"pin.obj", "ascii.ply", "little.ply", "big.ply"}) {
        ImportedMesh m;
        CHECK(readMeshFile(d + "/" + file, ImportLimits{}, m, &err));
        p.check(std::string("read ") + file, meshText(m));
        p.check(std::string("analysis ") + file, analysisText(analyzeImport(m)));
    }
    // What a user is told, with the folder taken out.
    const std::vector<std::pair<std::string, std::string>> bad = {
        {"no_face.obj", "v 0 0 0\n"},
        {"short_vertex.obj", "v 1 2\n"},
        {"bad_number.obj", "v 1 2 x\n"},
        {"bad_index.obj", "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 9\n"},
        {"partial_uv.obj", "v 0 0 0\nv 1 0 0\nv 0 1 0\nvt 0 0\nf 1/1 2 3\n"},
        {"nan.obj", "v nan 0 0\n"},
        {"magic.ply", "not a ply\n"},
        {"cloud.ply", "ply\nformat ascii 1.0\nelement vertex 1\nproperty float x\nproperty float y\nproperty float z\nend_header\n0 0 0\n"},
        {"truncated.ply", plyBinary(false).substr(0, plyHeader("binary_little_endian").size() + 8)},
        {"too_many.ply", "ply\nformat ascii 1.0\nelement vertex 9000000\nproperty float x\nelement face 1\nproperty list uchar int vertex_indices\nend_header\n"},
    };
    std::string messages;
    for (const auto& [name, text] : bad) {
        writeBytes(d + "/" + name, text);
        ImportedMesh m;
        CHECK(!readMeshFile(d + "/" + name, ImportLimits{}, m, &err));
        std::string message = err;
        for (size_t at = message.find(d + "/"); at != std::string::npos; at = message.find(d + "/")) message.erase(at, d.size() + 1);
        messages += message + "\n";
    }
    p.check("refusals", messages);
    p.finish();
}

TEST(io_pin_prepared_import) {
    const std::string d = freshDir("prepared");
    writeImportFiles(d);
    ImportedMesh m;
    std::string err;
    CHECK(readMeshFile(d + "/pin.obj", ImportLimits{}, m, &err));
    PreparedMesh prepared;
    CHECK(prepareImport(m, ImportOptions{}, prepared, &err));
    CHECK(storeImport(d + "/root", {std::make_shared<const ImportedMesh>(m), ImportOptions{}, {}, 0x1234}, {prepared, {}}, &err));
    const std::string entry = d + "/root/" + identityHex(0x1234);
    Pins p;
    p.check("kmesh", fileBytes(entry + "/mesh.kmesh"));
    std::string colours;
    const std::vector<ColourImage> textures = [&] {
        std::vector<std::string> warnings;
        return loadImportTextures(m, ImportLimits{}, &warnings);
    }();
    for (const Vec3& c : triangleColours(prepared, textures, Vec3{0.5, 0.5, 0.5})) appendVec(colours, c);
    p.check("triangle colours", colours);
    SourceScene scene;
    CHECK(loadImportedScene(entry, scene, &err));
    std::string prims;
    for (const SourcePrimitive& g : scene.prims) {
        prims += g.material + " " + std::to_string(g.welded.pos.size()) + " " + std::to_string(g.welded.tri.size()) + "\n";
        appendVec(prims, g.baseColor);
        for (const Vec3& c : g.triColor) appendVec(prims, c);
    }
    p.check("imported scene", prims);
    p.finish();
}

TEST(io_pin_knit_formats) {
    const std::string d = freshDir("knit_formats");
    std::string err;
    Pins p;
    CHECK(writeKnitCache(d + "/a.knitcache", handKnitEntry(), &err));
    p.check("knitcache", fileBytes(d + "/a.knitcache"));

    Scene t;
    ScenePart knitted;
    knitted.name = "knitted";
    knitted.prim = 0;
    knitted.component = 1;
    knitted.centre = {1, 2, 3};
    auto mesh = std::make_shared<TriMesh>();
    mesh->pos = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
    mesh->tri = {{0, 1, 2}};
    knitted.mesh = mesh;
    knitted.yarnD = 1.16;
    knitted.note = "a note";
    knitted.knit = std::make_unique<PartKnit>();
    knitted.knit->t.curves = handCurves();
    knitted.knit->t.prog.stitches = handCurves().topo.stitches;
    knitted.knit->t.charts.resize(1);
    knitted.knit->t.charts[0].coverage = 0.75;
    knitted.knit->detail = "detail";
    t.parts.push_back(std::move(knitted));
    ScenePart patch;
    patch.name = "patch";
    patch.gapPatch = true;
    patch.mesh = mesh;
    patch.patch.bary[0] = 0.25;
    patch.patch.radius = 12.5;
    t.parts.push_back(std::move(patch));
    t.cross.lines = {"a line"};
    p.check("tiled knit", serializeKnit(t));

    const FibreSet fs = handFibres();
    CHECK(writeFibreCache(d + "/a.fibrecache", {"fibre key = 1\nfibre params = 2\n", 0.5f, 0x77}, fs, &err));
    p.check("fibrecache", fileBytes(d + "/a.fibrecache"));

    TriMesh digestMesh = *mesh;
    p.check("mesh digest", std::to_string(meshDigest(digestMesh)));
    p.finish();

    // Pinned bytes read back into equal bytes, so a file written before any of this still reads.
    KnitCacheEntry back;
    CHECK(readKnitCache(d + "/a.knitcache", back, &err));
    CHECK(writeKnitCache(d + "/b.knitcache", back, &err));
    CHECK(fileBytes(d + "/a.knitcache") == fileBytes(d + "/b.knitcache"));
    const std::optional<FibreCacheEntry> saved = readFibreCache(d + "/a.fibrecache", &err);
    CHECK(saved.has_value());
    if (!saved) return;
    const FibreFileIdentity& id = saved->identity;
    CHECK(id.key == "fibre key = 1\nfibre params = 2\n" && id.keyBound == 0.5f && id.fibreDigest == 0x77);
    CHECK(writeFibreCache(d + "/b.fibrecache", id, saved->fs, &err));
    CHECK(fileBytes(d + "/a.fibrecache") == fileBytes(d + "/b.fibrecache"));
}

// Only fibre-cache v6 is read; an older file says so instead of failing somewhere inside.
TEST(io_old_fibre_files_refused) {
    const std::string d = freshDir("old_fibres");
    for (const char* version : {"fibre-cache v3\n", "fibre-cache v4\n", "fibre-cache v5\n"}) {
        const std::string path = d + "/old.fibrecache";
        writeBytes(path, version + std::string(200, '\0') + version);
        std::string err;
        CHECK(!readFibreCache(path, &err));
        CHECK(err.find("no longer read") != std::string::npos);
    }
}

TEST(io_pin_validation_snapshot) {
    ValidationReport r;
    r.verdict = Status::Fail;
    r.generatorChecksApplied = true;
    r.lines.push_back({"gate", Status::Pass, "a short line"});
    r.lines.push_back({"long", Status::Fail, std::string(700, 'x')});  // the snapshot's buffer cuts it: digests hold that
    r.clearance.minInterYarn = 0.1;
    r.localSignature.detail = "signature";
    Pins p;
    p.check("snapshot", validationSnapshot(r));
    p.finish();
}

TEST(io_pin_scene_json) {
    Pins p;
    std::string err;
    {  // everything at once, in the wool look
        ExportScene ex(freshDir("scene_all"));
        setAllShading(ex);
        addHandParts(ex, err);
        ex.sortParts();
        CHECK(writeSceneJson(ex, &err));
        p.check("scene all", fileBytes(ex.folder() + "/scene.json"));
        p.check("textured mesh", fileBytes(ex.folder() + "/meshes/textured.ply"));
        p.check("ground mesh", fileBytes(ex.folder() + "/meshes/_ground.ply"));
        p.check("colour texture", fileBytes(ex.folder() + "/textures/colour_pin.png"));
        {  // the bytes are stb_image_write's deflate stream; the pixels stay the one colour given
            int w = 0, h = 0;
            std::vector<uint8_t> px;
            const std::vector<uint8_t> want = {10, 20, 30};
            CHECK(readPngRgb8(ex.folder() + "/textures/colour_pin.png", w, h, px) && w == 1 && h == 1 && px == want);
        }
    }
    {  // the other branches: bump map, a dense coat, an oriented override, no ground
        ExportScene ex(freshDir("scene_other"));
        ex.fibreLook.hairMaterials = true;
        ex.fibreLook.hair.model = MaterialSpec::Model::Fibre;
        ex.fibreLook.hairKeepsYarnColour = true;
        ex.plyShading.on = true;
        ex.plyShading.frame = PlyShading::Frame::Bump;
        ex.plyShading.height = "textures/h.exr";
        ex.plyShading.bumpScale = 1.5e-5;
        ex.groundPlane = false;
        addHandParts(ex, err);
        ex.sortParts();
        CHECK(writeSceneJson(ex, &err));
        p.check("scene other", fileBytes(ex.folder() + "/scene.json"));
    }
    {  // plain fibre colours, an override, no placement
        ExportScene ex(freshDir("scene_plain"));
        MaterialSpec own;
        own.color = {0.1, 0.2, 0.3};
        ex.fibreLook.orientedOverride = own;
        std::error_code ec;
        fs::create_directories(ex.folder() + "/meshes", ec);
        fs::create_directories(ex.folder() + "/strands", ec);
        TriMesh m;
        m.pos = {{0, 0, 0}, {0.1, 0, 0}, {0, 0.1, 0}};
        m.tri = {{0, 1, 2}};
        CHECK(ex.addMesh("plate", "black", m, &err));
        StrandSetF s;
        s.add({0, -0.01f, 0, 0.004f, -0.01f, 0}, {2e-4f, 2e-4f});
        CHECK(ex.addStrands("a", ex.orientedFibreMaterial("wool_white"), s, &err));
        CHECK(ex.addStrands("b", ex.randomFibreMaterial("bow"), s, &err));
        CHECK(writeSceneJson(ex, &err));
        p.check("scene plain", fileBytes(ex.folder() + "/scene.json"));
    }
    p.finish();
}

// The OBJ: its bytes in metres, where its numbers pass untouched, and in an import's coordinates.
TEST(io_pin_obj_files) {
    Pins p;
    std::string err;
    const std::string out = freshDir("obj_metres");
    {
        ExportScene ex(freshDir("obj_package"));
        addObjParts(ex, err);
        writeObjPin(ex, out, err);
        p.check("obj metres", fileBytes(out + "/pin.obj"));
        p.check("mtl", fileBytes(out + "/pin.mtl"));
    }
    for (const bool yUp : {true, false}) {
        ExportScene ex(freshDir("obj_package_import"));
        ex.importFrame = ImportFrame{2.5, {1, 2, 3}, yUp ? UpAxis::Y : UpAxis::Z};
        addObjParts(ex, err);
        const std::string there = freshDir("obj_import");
        writeObjPin(ex, there, err);
        p.check(yUp ? "obj import y up" : "obj import z up", fileBytes(there + "/pin.obj"));
    }
    p.finish();
}

// A whole export with no GPU: report.json, plain and diagnostic, and the files beside it.
TEST(io_pin_report_json) {
    Pins p;
    std::string err, log;
    SceneEntry e;
    Scene t;
    fillPlainScene(e, t);
    const SceneSettings s;
    {
        ExportOptions o;
        o.folder = freshDir("report_ready") + "/out";
        CHECK(exportScene({e, t, s}, o, &log, &err));
        p.check("report ready", fileBytes(o.folder + "/report.json"));
        p.check("plain scene json", fileBytes(o.folder + "/scene.json"));
        ExportOptions obj;
        obj.format = "obj";
        obj.folder = freshDir("report_obj") + "/out";
        obj.name = "pin";
        CHECK(exportScene({e, t, s}, obj, &log, &err));
        p.check("plain scene obj", fileBytes(obj.folder + "/pin.obj"));
        p.check("plain scene mtl", fileBytes(obj.folder + "/pin.mtl"));
        CHECK(fileBytes(obj.folder + "/report.json") == fileBytes(o.folder + "/report.json"));
    }
    {
        ExportOptions o;
        o.folder = freshDir("report_diagnostic") + "/out";
        o.status = ExportEligibility::State::Diagnostic;
        CHECK(exportScene({e, t, s}, o, &log, &err));
        p.check("report diagnostic", fileBytes(o.folder + "/report.json"));
    }
    p.finish();
}

// A Windows path in the settings: report.json stays valid JSON.
TEST(export_report_escapes_the_settings) {
    std::string err, log;
    SceneEntry e;
    Scene t;
    fillPlainScene(e, t);
    SceneSettings s;
    CHECK(setSetting(s, "import.source", "C:\\models\\pin.obj", &err));
    ExportOptions o;
    o.folder = freshDir("report_backslash") + "/out";
    CHECK(exportScene({e, t, s}, o, &log, &err));
    CHECK(fileBytes(o.folder + "/report.json").find("C:\\\\models\\\\pin.obj") != std::string::npos);
}

// Every row of the format table is complete and reaches a writer: a format with no case in
// writeFormat fails here, by name, since MSVC does not warn about a missing case.
TEST(export_every_format_writes) {
    std::string err, log;
    SceneEntry e;
    Scene t;
    fillPlainScene(e, t);
    const SceneSettings s;
    std::map<std::string, int> ids;
    for (const ExportFormat& f : exportFormats()) {
        std::printf("  %s\n", f.id.c_str());
        CHECK(++ids[f.id] == 1);
        for (const std::string& text : {f.label, f.name, f.shortName, f.contentsHeading, f.contents, f.summary, f.meshSummary}) CHECK(!text.empty());
        if (!f.available) {
            CHECK(!f.unavailable.empty());
            continue;
        }
        ExportOptions o;
        o.format = f.id;
        o.folder = freshDir("every_format_" + f.id) + "/out";
        CHECK(exportScene({e, t, s}, o, &log, &err));
        CHECK(fs::exists(o.folder + "/report.json") && !fs::exists(o.folder + "/.package"));
        for (const std::string& extension : f.extensions) CHECK(fs::exists(o.folder + "/out" + extension));
    }
    ExportOptions unknown;
    unknown.format = "no_such_format";
    unknown.folder = freshDir("every_format_unknown") + "/out";
    CHECK(!exportScene({e, t, s}, unknown, &log, &err) && !fs::exists(unknown.folder));
}
