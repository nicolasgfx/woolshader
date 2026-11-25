// The knittable import as a whole (import/wrap_surface.h, import/import_pipeline.h), part of knit_tests_import: small generated
// meshes, on the CPU. The properties of each stage are in tests_wrap_properties.cpp and tests_wrap_flatten.cpp.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <sstream>
#include <thread>
#include <utility>
#include <vector>

#include "harness.h"
#include "import/import_pipeline.h"
#include "import/wrap_stages.h"
#include "import/wrap_surface.h"
#include "io/cache_paths.h"
#include "jobs/import_session.h"
#include "scene/knit/knit_stage.h"
#include "scene/scene.h"
#include "scene/appearance/stitch_colours.h"
#include "io/glb_load.h"
#include "io/image_io.h"
#include "import/import_analysis.h"
#include "import/import_store.h"
#include "settings/scene_registry.h"
#include "settings/scene_settings.h"
#include "knitting/layout/uv_atlas.h"
#include "model/cancel.h"
#include "scene/result_keys.h"
#include "model/coordinates.h"
#include "model/parallel.h"
#include "settings/settings_profiles.h"

using namespace knit;
namespace fs = std::filesystem;

namespace {

// --- Fixtures of the colours and the import.

std::string freshDir(const std::string& name) {
    const fs::path p = fs::temp_directory_path() / "knit_wrap_tests" / name;
    std::error_code ec;
    fs::remove_all(p, ec);
    fs::create_directories(p, ec);
    return p.string();
}

void writeText(const std::string& path, const std::string& text) {
    std::ofstream f(path, std::ios::binary);
    f << text;
}

struct ImportCacheAt {
    explicit ImportCacheAt(const std::string& root) { _putenv_s("KNIT_IMPORT_CACHE", root.c_str()); }
    ~ImportCacheAt() { _putenv_s("KNIT_IMPORT_CACHE", ""); }
};

ImportOptions wrappedZUp() {
    ImportOptions o;
    o.unitMm = 1;
    o.up = UpAxis::Z;
    o.surface = ImportOptions::Surface::Wrapped;
    return o;
}

std::string readText(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    std::stringstream text;
    text << file.rdbuf();
    return text.str();
}

const char* kBoxObj =
    "mtllib box.mtl\nv 0 0 0\nv 90 0 0\nv 90 90 0\nv 0 90 0\nv 0 0 90\nv 90 0 90\nv 90 90 90\nv 0 90 90\nusemtl paint\n"
    "f 1 4 3\nf 1 3 2\nf 1 2 6\nf 1 6 5\nf 3 4 8\nf 3 8 7\nf 4 1 5\nf 4 5 8\nf 2 3 7\nf 2 7 6\nf 5 6 7\nf 5 7 8\n";

// The distinct colours of an image's texels that are not black.
std::set<uint32_t> coloursOf(const ColourImage& t) {
    std::set<uint32_t> colours;
    for (size_t i = 0; i + 3 < t.rgba.size(); i += 4) {
        const uint32_t key = uint32_t(t.rgba[i]) << 16 | uint32_t(t.rgba[i + 1]) << 8 | t.rgba[i + 2];
        if (key) colours.insert(key);
    }
    return colours;
}

// A box with one flat colour, through the built-in engine into the store at `root`. Returns its identity.
uint64_t storedBox(const std::string& folder, const std::string& root) {
    writeText(folder + "/box.obj", kBoxObj);
    writeText(folder + "/box.mtl", "newmtl paint\nKd 0.6 0.2 0.1\n");
    ImportedMesh m;
    std::string err;
    CHECK(readMeshFile(folder + "/box.obj", ImportLimits{}, m, &err));
    const ImportOptions o = wrappedZUp();
    const uint64_t identity = importIdentity(m, o);
    std::vector<std::string> changes;
    CHECK(prepareAndStore({std::make_shared<ImportedMesh>(m), o, ImportLimits{}, identity}, changes, &err, {}));
    CHECK(!changes.empty());
    return identity;
}

bool waitForSession(const ImportSession& s, ImportSession::Phase wanted) {
    for (int i = 0; i < 12000; ++i) {
        if (s.snapshot().phase == wanted && !s.busy()) return true;
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    return false;
}

} // namespace

// --- the import

TEST(import_identity_unchanged) {
    // Computed by knit_viewer --import at e4dc1b2, before the wrapped import existed. A saved import
    // refers to its identity: keep-file identities must never move.
    const std::string d = freshDir("identity");
    writeText(d + "/pin.obj", "v 0 0 0\nv 10 0 0\nv 10 10 0\nv 0 10 0\nf 1 2 3\nf 1 3 4\n");
    ImportedMesh m;
    std::string err;
    CHECK(readMeshFile(d + "/pin.obj", ImportLimits{}, m, &err));
    CHECK(identityHex(importIdentity(m, ImportOptions{})) == "9c007be01953941d");
}

TEST(import_wrap_identity_differs) {
    const std::string d = freshDir("identity_wrap");
    writeText(d + "/box.obj", kBoxObj);
    writeText(d + "/box.mtl", "newmtl paint\nKd 0.6 0.2 0.1\n");
    ImportedMesh m;
    std::string err;
    CHECK(readMeshFile(d + "/box.obj", ImportLimits{}, m, &err));
    ImportOptions keep, wrapped = keep, finer = keep;
    wrapped.surface = finer.surface = ImportOptions::Surface::Wrapped;
    finer.wrapStitchMm = 4.5;
    CHECK(importIdentity(m, keep) != importIdentity(m, wrapped));
    CHECK(importIdentity(m, wrapped) != importIdentity(m, finer));
    // The options survive the settings file exactly: otherwise a reopened import is prepared again.
    SceneSettings s;
    s.scene = "import";
    s.importSurface = "wrapped";
    s.importWrapStitchMm = 5.37123456789012;
    wrapped.wrapStitchMm = s.importWrapStitchMm;
    SceneSettings back;
    CHECK(readSettings(writeSettings(s), back, &err));
    CHECK(importIdentity(m, importOptionsFrom(back)) == importIdentity(m, wrapped));
}

TEST(import_settings_old_file) {
    SceneSettings s;
    s.scene = "import";
    std::string text = writeSettings(s), kept;
    std::istringstream lines(text);
    for (std::string line; std::getline(lines, line);)
        if (line.rfind("import.surface", 0) != 0 && line.rfind("import.wrap_stitch_mm", 0) != 0) kept += line + "\n";
    SceneSettings back;
    std::string err;
    CHECK(readSettings(kept, back, &err));
    CHECK(back.importSurface == "file");
    CHECK(importOptionsFrom(back).surface == ImportOptions::Surface::KeepFile);
    // import.charts is retired (2026-11-11): files that hold it still load, the key is never written again, and only a
    // value other than "keep" is reported.
    CHECK(writeSettings(back).find("import.charts") == std::string::npos);
    for (const char* value : {"keep", "generate"}) {
        SceneSettings old;
        std::vector<std::string> retired;
        CHECK(readSettings(text + "import.charts = " + value + "\n", old, &err, &retired));
        CHECK(writeSettings(old) == writeSettings(s) && retired.size() == (std::string(value) == "keep" ? 0u : 1u));
    }
}

TEST(import_stage_keys_unchanged) {
    // The new keys are inputs of a wrapped import only: a keep-file import's saved results keep their keys.
    SceneSettings keep;
    keep.scene = "import";
    const std::string inputs = stageInputs(keep, Stage::Export);
    CHECK(inputs.find("import.surface") == std::string::npos && inputs.find("import.wrap_stitch_mm") == std::string::npos);
    CHECK(inputs.find("import.charts") == std::string::npos);  // retired: charts are never generated on the file's mesh
    SceneSettings wrapped = keep;
    wrapped.importSurface = "wrapped";
    const std::string wrappedInputs = stageInputs(wrapped, Stage::Export);
    CHECK(wrappedInputs.find("import.surface=wrapped") != std::string::npos && wrappedInputs.find("import.wrap_stitch_mm") != std::string::npos);
}

TEST(import_wrap_roundtrip) {
    const std::string d = freshDir("roundtrip"), root = freshDir("roundtrip_cache");
    ImportCacheAt env(root);
    const uint64_t id = storedBox(d, root);
    CHECK(importStored(root, id));
    CHECK(!storedImportChanges(root, id).empty());
    PreparedMesh back;
    std::string source, err;
    CHECK(loadImport(importCacheDir(root, id), back, &source, &err));
    CHECK(!back.mesh.tri.empty() && back.chartsGenerated);
    CHECK(back.materials.size() == 1 && !back.materials[0].texture.empty());
    ImportedMesh holder;
    holder.materials = back.materials;
    std::vector<std::string> warnings;
    const std::vector<ColourImage> textures = loadImportTextures(holder, ImportLimits{}, &warnings);
    CHECK(textures.size() == 1 && textures[0].width == kBakedTextureSide);
    // The working folder is gone: a stored import keeps its own copies.
    for (const fs::directory_entry& entry : fs::directory_iterator(root)) CHECK(entry.path().filename().string().rfind(".wrap_", 0) != 0);
}

TEST(import_wrap_one_primitive) {
    const std::string d = freshDir("one_primitive"), root = freshDir("one_primitive_cache");
    ImportCacheAt env(root);
    const uint64_t id = storedBox(d, root);
    SourceScene g;
    std::string err;
    CHECK(loadImportedScene(importCacheDir(root, id), g, &err));
    CHECK(g.imported && g.prims.size() == 1 && g.prims[0].albedo);
}

TEST(import_wrap_session_cancel) {
    const std::string d = freshDir("session_cancel"), root = freshDir("session_cancel_cache");
    ImportCacheAt env(root);
    // Large enough for the wrap to be running when the cancel arrives.
    std::string obj;
    for (const char* line : {"v 0 0 0\nv 390 0 0\nv 390 390 0\nv 0 390 0\nv 0 0 390\nv 390 0 390\nv 390 390 390\nv 0 390 390\n",
                             "f 1 4 3\nf 1 3 2\nf 1 2 6\nf 1 6 5\nf 3 4 8\nf 3 8 7\nf 4 1 5\nf 4 5 8\nf 2 3 7\nf 2 7 6\nf 5 6 7\nf 5 7 8\n"})
        obj += line;
    writeText(d + "/big.obj", obj);
    using Ph = ImportSession::Phase;
    ImportSession s;
    s.open(d + "/big.obj");
    CHECK(waitForSession(s, Ph::Analysed));
    s.prepare(wrappedZUp());
    for (int i = 0; i < 3000 && s.snapshot().stage.find("wrapping") == std::string::npos; ++i) std::this_thread::sleep_for(std::chrono::milliseconds(1));
    s.cancel();
    CHECK(waitForSession(s, Ph::Analysed));
    CHECK(s.snapshot().identity == 0 && s.snapshot().error.empty());
}

TEST(import_wrap_session_reuses_a_stored_wrap) {
    const std::string d = freshDir("session_reuse"), root = freshDir("session_reuse_cache");
    ImportCacheAt env(root);
    writeText(d + "/box.obj", kBoxObj);
    writeText(d + "/box.mtl", "newmtl paint\nKd 0.6 0.2 0.1\n");
    using Ph = ImportSession::Phase;
    ImportSession s;
    s.open(d + "/box.obj");
    CHECK(waitForSession(s, Ph::Analysed));
    s.prepare(wrappedZUp());
    CHECK(waitForSession(s, Ph::Prepared));
    const std::vector<std::string> first = s.snapshot().changes;
    const uint64_t identity = s.snapshot().identity;
    CHECK(!first.empty());
    // The stored wrap answers: it reports what it did, and a rewritten entry would have a new time.
    std::error_code ec;
    const auto written = fs::last_write_time(fs::path(importCacheDir(root, identity)) / "import.txt", ec);
    s.open(d + "/box.obj");
    CHECK(waitForSession(s, Ph::Analysed));
    s.prepare(wrappedZUp());
    CHECK(waitForSession(s, Ph::Prepared));
    CHECK(s.snapshot().identity == identity && s.snapshot().changes == first);
    CHECK(fs::last_write_time(fs::path(importCacheDir(root, identity)) / "import.txt", ec) == written);
}

// --- Settings from the time of an engine choice

// The import has one engine since 2026-11-12. Settings that still name an engine (any name, builtin included) or a recipe
// load without error and import exactly as the same file without those keys: the keys are accepted and ignored, never
// written again, and the import options and identity hold no engine.
TEST(import_retired_engine_keys_load) {
    ImportedMesh m;
    m.path = "x.obj";
    SceneSettings wrapped;
    wrapped.scene = "import";
    wrapped.importSurface = "wrapped";
    wrapped.importWrapStitchMm = 4.5;
    const std::string text = writeSettings(wrapped);
    CHECK(text.find("import.wrap_engine") == std::string::npos && text.find("import.wrap_recipe") == std::string::npos);
    std::string err;
    SceneSettings plain;
    CHECK(readSettings(text, plain, &err));
    const ImportOptions expected = importOptionsFrom(plain);
    CHECK(importOptionsText(expected).find("engine") == std::string::npos);
    CHECK(importOptionsText(expected).find("wrap_version = 6\n") != std::string::npos);
    const char* retiredLines = "import.wrap_engine = external\nimport.wrap_recipe = \n";
    for (const char* lines : {retiredLines, "import.wrap_engine = builtin\n",
                              "import.wrap_engine = external\nimport.wrap_recipe = finish91\n"}) {
        SceneSettings old;
        std::vector<std::string> retired;
        CHECK(readSettings(text + lines, old, &err, &retired));
        CHECK(writeSettings(old) == text);
        CHECK(importOptionsFrom(old) == expected && importIdentity(m, importOptionsFrom(old)) == importIdentity(m, expected));
        // Only a value that asked for something else is reported: another engine, or a recipe.
        const std::string asked = lines;
        CHECK((std::count(retired.begin(), retired.end(), "import.wrap_engine") == 1) == (asked.find("external") != std::string::npos));
        CHECK((std::count(retired.begin(), retired.end(), "import.wrap_recipe") == 1) == (asked.find("finish91") != std::string::npos));
    }
    // Real saved settings files, as they are and with the retired lines added, give the same settings and import.
    for (const char* file : {"tests/data/swatch_settings/swatch_r1.settings", "tests/data/swatch_settings/swatch_r2.settings",
                             "tests/data/swatch_settings/swatch_blue.settings"}) {
        const std::string onDisk = readText(std::string(KNIT_SOURCE_DIR) + "/" + file);
        CHECK(!onDisk.empty() && onDisk.find("import.wrap_engine") == std::string::npos);
        const std::string withRetired = onDisk + (onDisk.empty() || onDisk.back() == '\n' ? "" : "\n") + retiredLines;
        SceneSettings asSaved, old;
        CHECK(readSettings(onDisk, asSaved, &err));
        CHECK(readSettings(withRetired, old, &err));
        CHECK(writeSettings(old) == writeSettings(asSaved));
        CHECK(importOptionsFrom(old) == importOptionsFrom(asSaved));
        CHECK(importIdentity(m, importOptionsFrom(old)) == importIdentity(m, importOptionsFrom(asSaved)));
    }
}

TEST(import_keep_file_fur_key_unchanged) {
    SceneSettings keep;
    keep.scene = "import";
    const uint64_t key = furKey(keep, 7);
    SceneSettings other = keep;
    other.importWrapStitchMm = 3;
    CHECK(furKey(other, 7) == key);
    CHECK(stageInputs(keep, Stage::Export) == stageInputs(other, Stage::Export));
    // A wrapped import's surface does depend on them.
    SceneSettings wrapped = keep;
    wrapped.importSurface = "wrapped";
    const uint64_t wrappedKey = furKey(wrapped, 7);
    CHECK(wrappedKey != key);
    wrapped.importWrapStitchMm = 3;
    CHECK(furKey(wrapped, 7) != wrappedKey);
}

// Overlapped seams run rows past a chart's border. Those stitches lie in no cell of their chart, and
// must still take the object's colour: without it they were knitted mid-grey along every seam.
TEST(colourwork_past_the_chart_border) {
    const std::string d = freshDir("past_border"), root = freshDir("past_border_cache");
    ImportCacheAt env(root);
    const uint64_t id = storedBox(d, root);
    std::string err;
    SceneSettings s;
    s.scene = "import";
    s.importSource = d + "/box.obj";
    s.importIdentity = identityHex(id);
    s.importUnitMm = 1;
    s.importUp = "z";
    s.importSurface = "wrapped";
    s.yarnColourSource = "texture";
    const SceneEntry* e = sceneEntryFor(s);
    CHECK(e != nullptr);
    applyPreviewProfile(s);
    s.coverage = "overlap";
    Scene t;
    CHECK(knitScene(*e, s, KnitRunOptions{}, t, &err));
    // One flat colour: every stitch, on its chart or past its border, is that colour.
    const SceneStitchColours colours = sceneStitchColours(t);
    CHECK(colours.stitchesOffSurface == 0);
    PreparedMesh back;
    std::string source;
    CHECK(loadImport(importCacheDir(root, id), back, &source, &err));
    ImportedMesh holder;
    holder.materials = back.materials;
    std::vector<std::string> warnings;
    const std::set<uint32_t> baked = coloursOf(loadImportTextures(holder, ImportLimits{}, &warnings).at(0));
    CHECK(baked.size() == 1);
    const uint32_t expected = *baked.begin();
    size_t stitches = 0;
    for (const StitchColours& c : colours.parts) {
        stitches += c.key.size();
        for (uint32_t key : c.key)
            for (int shift : {0, 8, 16}) CHECK(std::abs(int(key >> shift & 255) - int(expected >> shift & 255)) <= 2);
    }
    CHECK(stitches > 100);
}

// --- A selection: the ticked nodes knitted, the others plain meshes in their place.

namespace {

// Two 30 mm boxes 70 mm apart along x: `left` flat red in triangles, `right` textured in quads (quadFirstTriangle, colour uv).
void writeTwoBoxes(const std::string& dir) {
    fs::copy_file(std::string(KNIT_TEST_DATA) + "/wrap_fixtures/torus/checker.png", dir + "/checker.png", fs::copy_options::overwrite_existing);
    writeText(dir + "/boxes.mtl", "newmtl red\nKd 0.8 0.1 0.1\nnewmtl chequered\nmap_Kd checker.png\n");
    std::string obj = "mtllib boxes.mtl\nvt 0 0\nvt 1 0\nvt 1 1\nvt 0 1\n";
    for (const double x : {0.0, 100.0})
        for (const double z : {0.0, 30.0})
            for (const double y : {0.0, 30.0})
                for (const double dx : {0.0, 30.0}) obj += "v " + std::to_string(x + dx) + " " + std::to_string(y) + " " + std::to_string(z) + "\n";
    obj += "o left\nusemtl red\n";
    for (const char* f : {"1 3 4", "1 4 2", "5 6 8", "5 8 7", "1 2 6", "1 6 5", "3 7 8", "3 8 4", "1 5 7", "1 7 3", "2 4 8", "2 8 6"}) obj += std::string("f ") + f + "\n";
    obj += "o right\nusemtl chequered\n";
    for (const char* f : {"9/1 11/2 12/3 10/4", "13/1 14/2 16/3 15/4", "9/1 10/2 14/3 13/4", "11/1 15/2 16/3 12/4", "9/1 13/2 15/3 11/4", "10/1 12/2 16/3 14/4"}) obj += std::string("f ") + f + "\n";
    writeText(dir + "/boxes.obj", obj);
}

ImportedMesh readTwoBoxes(const std::string& dir) {
    writeTwoBoxes(dir);
    ImportedMesh m;
    std::string err;
    CHECK(readMeshFile(dir + "/boxes.obj", ImportLimits{}, m, &err));
    finishNodes(m);
    return m;
}

struct Extent {
    Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
};

Extent extentOf(const std::vector<Vec3>& points) {
    Extent e;
    for (const Vec3& q : points) {
        e.lo = {std::min(e.lo.x, q.x), std::min(e.lo.y, q.y), std::min(e.lo.z, q.z)};
        e.hi = {std::max(e.hi.x, q.x), std::max(e.hi.y, q.y), std::max(e.hi.z, q.z)};
    }
    return e;
}

const SourcePrimitive* primitiveOf(const SourceScene& scene, SourcePrimitive::Class cls) {
    for (const SourcePrimitive& g : scene.prims)
        if (g.cls == cls) return &g;
    return nullptr;
}

// The stored scene of a two-box file with `right` plain: the knitted box left of the origin, the plain one right of it.
SourceScene preparedSelection(ImportOptions o, const std::string& d, const std::string& root) {
    writeTwoBoxes(d);
    using Ph = ImportSession::Phase;
    ImportSession s;
    s.open(d + "/boxes.obj");
    CHECK(waitForSession(s, Ph::Analysed));
    o.plainNodes = {"right"};
    s.prepare(o);
    CHECK(waitForSession(s, Ph::Prepared));
    SourceScene scene;
    std::string err;
    CHECK(loadImportedScene(importCacheDir(root, s.snapshot().identity), scene, &err));
    return scene;
}

} // namespace

TEST(import_split_by_nodes) {
    const ImportedMesh m = readTwoBoxes(freshDir("split"));
    CHECK(m.nodes.size() == 2 && m.nodes[0].name == "left" && m.nodes[1].name == "right");
    CHECK(m.nodes[0].triangles == 12 && m.nodes[1].triangles == 12);
    const NodeSplit split = splitByNodes(m, {"right"});
    CHECK(split.knitted.tri.size() == 12 && split.plain.tri.size() == 12);
    CHECK(split.knitted.pos.size() == 8 && split.plain.pos.size() == 8);  // only the vertices each side uses
    CHECK(extentOf(split.knitted.pos).hi.x == 30 && extentOf(split.plain.pos).lo.x == 100);
    CHECK(split.knitted.triMaterial[0] == 0 && split.plain.triMaterial[0] == 1 && split.plain.materials.size() == 2);
    CHECK(split.knitted.uv.size() == 12 && split.plain.uv.size() == 12 && split.plain.uv[0][1].x == 1);
    CHECK(split.knitted.quadFirstTriangle.empty() && split.plain.quadFirstTriangle.size() == 6);
    for (int first : split.plain.quadFirstTriangle) CHECK(first < 12);  // counted in the side's own triangles
    const NodeSplit unknown = splitByNodes(m, {"gone from the file"});
    CHECK(unknown.knitted.tri.size() == 24 && unknown.plain.tri.empty());
}

TEST(import_selection_keep_file) {
    const std::string d = freshDir("keep_file"), root = freshDir("keep_file_cache");
    ImportCacheAt env(root);
    ImportOptions o;
    o.up = UpAxis::Z;
    const SourceScene scene = preparedSelection(o, d, root);
    const SourcePrimitive* knitted = primitiveOf(scene, SourcePrimitive::Class::Knit);
    const SourcePrimitive* plain = primitiveOf(scene, SourcePrimitive::Class::Plain);
    CHECK(knitted && plain && scene.prims.size() == 2);
    // Both stand where the whole object stands: 130 mm wide, centred, on the floor.
    const Extent k = extentOf(knitted->welded.pos), p = extentOf(plain->welded.pos);
    CHECK(k.lo.x == -65 && k.hi.x == -35 && p.lo.x == 35 && p.hi.x == 65);
    CHECK(k.lo.z == 0 && p.lo.z == 0 && k.lo.y == -15 && p.hi.y == 15);
    CHECK(scene.lo.x == -65 && scene.hi.x == 65);
    // The plain box keeps its texture and colour uv; the export writes both.
    CHECK(plain->albedo && plain->colourUv.size() == 12 && plain->material == "chequered");
    CHECK(!knitted->albedo);
}

// The rebuilt surface is made from the selected box alone, and stays beside the plain one.
TEST(import_selection_built_in) {
    const std::string d = freshDir("built_in"), root = freshDir("built_in_cache");
    ImportCacheAt env(root);
    const SourceScene scene = preparedSelection(wrappedZUp(), d, root);
    const SourcePrimitive* knitted = primitiveOf(scene, SourcePrimitive::Class::Knit);
    const SourcePrimitive* plain = primitiveOf(scene, SourcePrimitive::Class::Plain);
    CHECK(knitted && plain);
    const Extent k = extentOf(knitted->welded.pos), p = extentOf(plain->welded.pos);
    CHECK(std::fabs(k.lo.x + 65) < 3 && std::fabs(k.hi.x + 35) < 3);  // the rebuild grows the surface a little
    CHECK(std::fabs(k.lo.z) < 3);
    CHECK(p.lo.x == 35 && p.hi.x == 65 && p.lo.z == 0);
}

TEST(import_selection_of_nothing_is_refused) {
    const std::string d = freshDir("nothing"), root = freshDir("nothing_cache");
    ImportCacheAt env(root);
    writeTwoBoxes(d);
    using Ph = ImportSession::Phase;
    ImportSession s;
    s.open(d + "/boxes.obj");
    CHECK(waitForSession(s, Ph::Analysed));
    ImportOptions o;
    o.plainNodes = {"left", "right"};
    s.prepare(o);
    CHECK(waitForSession(s, Ph::Failed));
    CHECK(s.snapshot().error.find("every component is left out") != std::string::npos);
}

// A selection changes the identity, and through it every key; the row itself is no input, so no existing key moves.
TEST(import_selection_identity_and_keys) {
    const ImportedMesh m = readTwoBoxes(freshDir("selection_keys"));
    ImportOptions whole, part;
    part.plainNodes = {"right", "left"};
    CHECK(importOptionsText(whole).find("plain_nodes") == std::string::npos);
    CHECK(importOptionsText(part).find("plain_nodes = left|right\n") != std::string::npos);
    CHECK(importIdentity(m, whole) != importIdentity(m, part));
    part.surface = ImportOptions::Surface::Wrapped;
    CHECK(importOptionsText(part).find("plain_nodes = left|right\n") != std::string::npos);
    CHECK(parsePlainNodes(plainNodesText(part.plainNodes)) == part.plainNodes && parsePlainNodes("").empty());

    SceneSettings a;
    a.scene = "import";
    a.importIdentity = "0123456789abcdef";
    SceneSettings b = a;
    b.importPlainNodes = "left|right";
    CHECK(stageKey(a, Stage::Surface, 0) == stageKey(b, Stage::Surface, 0) && furKey(a, 0) == furKey(b, 0));
    SceneSettings read;
    std::string err;
    CHECK(readSettings(writeSettings(b), read, &err, nullptr) && read.importPlainNodes == "left|right");
    CHECK(importOptionsFrom(read).plainNodes == part.plainNodes);
}
