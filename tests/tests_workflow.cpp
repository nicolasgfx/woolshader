// The viewer workflow's backend (knit_tests_workflow): import (mesh_read.h, prepare_import.h, import_store.h), the
// result model and export eligibility (scene_builder.h), the export formats (export.h)
// and the builder's cancellation. Small synthetic meshes: everything runs in seconds on the CPU.
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <map>
#include <mutex>
#include <iterator>
#include <regex>
#include <set>
#include <sstream>
#include <thread>

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#include <psapi.h>
#pragma comment(lib, "psapi.lib")
#endif
#include "harness.h"
#include "io/cache_paths.h"
#include "jobs/results.h"
#include "knitting/tube/param_text.h"
#include "scene/stage_params.h"
#include "settings/settings_profiles.h"
#include "test_compute.h"
#include "io/image_io.h"
#include "fibres/fibre_checks.h"
#include "model/choices.h"
#include "model/colour.h"
#include "model/coordinates.h"
#include "model/hash.h"
#include "model/log.h"
#include "knitting/checks/checks.h"
#include "knitting/layout/uv_atlas.h"
#include "knitting/stitch/hand_tension.h"
#include "export/export.h"
#include "io/blender_process.h"
#include "jobs/import_session.h"
#include "jobs/job_steps.h"
#include "import/import_pipeline.h"
#include "import/import_store.h"
#include "import/import_analysis.h"
#include "import/prepare_import.h"
#include "io/publish.h"
#include "jobs/scene_builder.h"
#include "jobs/scene_cache.h"
#include "scene/appearance/stitch_colours.h"
#include "settings/scene_registry.h"
#include "settings/thickness_presets.h"
#include "scene/knit/saved_knit.h"
#include "scene/scene.h"
#include "settings/settings_mapping.h"
#include "scene/source/scene_parts.h"
#include "scene/fibres/part_fibre_job.h"
#include "scene/source/scene_solids.h"
#include "scene/knit/knit_stage.h"
#include "scene/appearance/rendered_part.h"
#include "scene/fibres/fibre_stage.h"
#include "settings/yarn_look.h"
#include "fibres/fur_field.h"
#include "fibres/fur_loops.h"
#include "jobs/fur_preview.h"
#include "fibres/fibre_strands.h"
#include "jobs/previews.h"
#include "scene/fibres/fur_stage.h"
#include "fibres/fibre.h"
#include "fibres/fur.h"
#include "model/cancel.h"
#include "model/mesh_query.h"
#include "fibres/fibre_strands.h"
#include "io/files.h"
#include "export/export_scene.h"
#include "fibres/fibre_budget.h"
#include "scene/result_keys.h"
#include "scene/fibres/fibre_files.h"

using namespace knit;

// A stage's optional result into `out`; false when it is empty.
static bool grown(SceneFibres& out, std::optional<SceneFibres> result) {
    if (!result) return false;
    out = std::move(*result);
    return true;
}
namespace fs = std::filesystem;

namespace {

// A scene knitted on its UV charts (an import's or a GLB's), for the settings rules that differ from the
// teddy's spiral pieces: no registered preset has been one since the man went (owner, 2026-11-10).
SceneEntry uvChartScene() {
    SceneEntry e;
    e.id = "charts";
    e.label = "a scene knitted on its UV charts";
    e.kind = SceneEntry::Kind::Mesh;
    e.glb.knitAllMaterials = true;
    return e;
}

std::string freshDir(const std::string& name) {
    const fs::path p = fs::temp_directory_path() / "knit_workflow_tests" / name;
    std::error_code ec;
    fs::remove_all(p, ec);
    fs::create_directories(p, ec);
    return p.string();
}

void writeText(const std::string& path, const std::string& text) {
    fs::create_directories(fs::path(path).parent_path());
    std::ofstream f(path, std::ios::binary);
    f << text;
}

std::string readText(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    return std::string((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
}

// An nx x ny grid over w x h (file units) in the xy plane, uv = (x / w, y / h), optional
// per-vertex colour; as OBJ text (with a material) or vertices / faces.
struct Grid {
    std::vector<Vec3> pos;
    std::vector<Vec2> uv;
    std::vector<std::array<int, 3>> tri;
};
Grid grid(int nx, int ny, double w, double h) {
    Grid g;
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i) {
            g.pos.push_back({w * i / nx, h * j / ny, 0});
            g.uv.push_back({double(i) / nx, double(j) / ny});
        }
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            const int a = j * (nx + 1) + i, b = a + 1, c = a + nx + 1, d = c + 1;
            g.tri.push_back({a, b, d});
            g.tri.push_back({a, d, c});
        }
    return g;
}

std::string objText(const Grid& g, const std::string& mtllib, const std::string& mat, bool uv = true, bool colours = false) {
    std::ostringstream o;
    o.precision(17);
    if (!mtllib.empty()) o << "mtllib " << mtllib << "\n";
    for (size_t v = 0; v < g.pos.size(); ++v) {
        o << "v " << g.pos[v].x << " " << g.pos[v].y << " " << g.pos[v].z;
        if (colours) o << " " << g.uv[v].x << " 0.5 " << g.uv[v].y;
        o << "\n";
    }
    if (uv)
        for (const Vec2& t : g.uv) o << "vt " << t.x << " " << t.y << "\n";
    if (!mat.empty()) o << "usemtl " << mat << "\n";
    for (const auto& t : g.tri) {
        o << "f";
        for (int k : t) o << " " << (k + 1) << (uv ? "/" + std::to_string(k + 1) : std::string());
        o << "\n";
    }
    return o.str();
}

void writePlyFile(const std::string& path, const Grid& g, const std::string& format) {
    std::ofstream f(path, std::ios::binary);
    f.precision(9);  // float round trip in ascii
    f << "ply\nformat " << format << " 1.0\ncomment TextureFile tex.png\nelement vertex " << g.pos.size()
      << "\nproperty float x\nproperty float y\nproperty float z\nproperty float s\nproperty float t\nproperty uchar red\nproperty uchar green\n"
         "property uchar blue\nelement face "
      << g.tri.size() << "\nproperty list uchar int vertex_indices\nend_header\n";
    auto put = [&](const void* p, size_t n) {
        std::string b(static_cast<const char*>(p), n);
        if (format == "binary_big_endian") std::reverse(b.begin(), b.end());
        f.write(b.data(), static_cast<std::streamsize>(n));
    };
    for (size_t v = 0; v < g.pos.size(); ++v) {
        const float q[5] = {float(g.pos[v].x), float(g.pos[v].y), float(g.pos[v].z), float(g.uv[v].x), float(g.uv[v].y)};
        const unsigned char c[3] = {static_cast<unsigned char>(v % 256), 128, 7};
        if (format == "ascii") {
            f << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << " " << q[4] << " " << int(c[0]) << " " << int(c[1]) << " " << int(c[2]) << "\n";
        } else {
            for (float x : q) put(&x, 4);
            for (unsigned char x : c) put(&x, 1);
        }
    }
    for (const auto& t : g.tri) {
        if (format == "ascii") {
            f << "3 " << t[0] << " " << t[1] << " " << t[2] << "\n";
        } else {
            const unsigned char n = 3;
            put(&n, 1);
            for (int k : t) {
                const int32_t x = k;
                put(&x, 4);
            }
        }
    }
}

void writeTexture(const std::string& path, uint8_t r, uint8_t g, uint8_t b) {
    std::vector<uint8_t> rgb(16 * 16 * 3);
    for (size_t i = 0; i < rgb.size(); i += 3) rgb[i] = r, rgb[i + 1] = g, rgb[i + 2] = b;
    CHECK(writePngRgb8(path, 16, 16, rgb));
}

// The fixture files written to `dir` (plane.obj and what it names) as an import: read, prepared and
// stored in the import cache rooted at `root`; returns settings for scene = import (Draft, on the GPU
// with the explicit solver when available: test_compute.h).
SceneSettings importWrittenFixture(const std::string& dir, const std::string& root, uint64_t* identity = nullptr) {
    ImportedMesh m;
    std::string err;
    CHECK(readMeshFile(dir + "/plane.obj", ImportLimits{}, m, &err));
    ImportOptions o;
    o.unitMm = 1;
    o.up = UpAxis::Z;
    PreparedMesh p;
    CHECK(prepareImport(m, o, p, &err));
    const uint64_t id = importIdentity(m, o);
    CHECK(storeImport(root, {std::make_shared<const ImportedMesh>(m), o, {}, id}, {p, {}}, &err));
    if (identity) *identity = id;
    SceneSettings s;
    s.scene = "import";
    s.importSource = dir + "/plane.obj";
    s.importIdentity = identityHex(id);
    s.importUnitMm = 1;
    s.importUp = "z";
    const SceneEntry* e = sceneEntryFor(s);
    CHECK(e != nullptr);
    applySceneDefaults(s, *e);
    applyPreviewProfile(s);
    s.textureProfile = "none";
    return s;
}

// The plane fixture: OBJ + MTL + texture.
SceneSettings importPlane(const std::string& dir, const std::string& root, uint64_t* identity = nullptr) {
    writeText(dir + "/plane.mtl", "newmtl wool\nKd 0.2 0.4 0.6\nmap_Kd tex.png\nNs 10\n");
    writeTexture(dir + "/tex.png", 200, 100, 50);
    writeText(dir + "/plane.obj", objText(grid(12, 8, 60, 40), "plane.mtl", "wool"));
    return importWrittenFixture(dir, root, identity);
}

// The plane fixture with a texture of 64 x 64 texels, left half red, right half green.
constexpr double kRedGreenEdgeX = 0;  // library mm: the import centres the 60 x 40 plane
constexpr double kRedGreenTexelMm = 60.0 / 64;
void writeRedGreenTexture(const std::string& path) {
    std::vector<uint8_t> rgb(64 * 64 * 3, 0);
    for (int y = 0; y < 64; ++y)
        for (int x = 0; x < 64; ++x) rgb[3 * size_t(y * 64 + x) + (x < 32 ? 0 : 1)] = 255;
    CHECK(writePngRgb8(path, 64, 64, rgb));
}
SceneSettings importRedGreenPlane(const std::string& dir, const std::string& root) {
    writeText(dir + "/plane.mtl", "newmtl wool\nKd 0.2 0.4 0.6\nmap_Kd tex.png\n");
    writeRedGreenTexture(dir + "/tex.png");
    writeText(dir + "/plane.obj", objText(grid(12, 8, 60, 40), "plane.mtl", "wool"));
    return importWrittenFixture(dir, root);
}

struct ImportCacheEnv {
    explicit ImportCacheEnv(const std::string& root) { _putenv_s("KNIT_IMPORT_CACHE", root.c_str()); }
    ~ImportCacheEnv() { _putenv_s("KNIT_IMPORT_CACHE", ""); }
};

struct StageLog : JobReporter {
    std::mutex m;
    std::vector<std::string> stages;
    std::vector<std::pair<int, int>> values;  // done, total of each report
    void progress(const std::string& stage, int done, int total) override {
        std::lock_guard<std::mutex> lk(m);
        stages.push_back(stage);
        values.push_back({done, total});
    }
    void partial(std::shared_ptr<BuiltScene>) override {}
};

std::shared_ptr<BuiltScene> run(BuiltScene::Kind kind, const SceneSettings& s, std::string* err, JobReporter* sink = nullptr,
                                std::shared_ptr<const KnitResult> knit = nullptr, std::shared_ptr<const FibreResult> fibres = nullptr,
                                const std::string& exportDir = "", const std::string& format = "ppt", bool exportFibres = false,
                                const std::string& cacheRoot = "") {
    const auto settings = std::make_shared<const SceneSettings>(s);
    const JobControl control{0, sink};
    if (kind == BuiltScene::Kind::Knit) return knitObject(KnitRequest{settings, cacheRoot}, control, err);
    if (kind == BuiltScene::Kind::Fibres) return growFibres(FibreRequest{settings, std::move(knit)}, control, err);
    ExportRequest x;
    x.settings = settings;
    x.knit = std::move(knit);
    x.fibres = std::move(fibres);
    x.dir = exportDir;
    x.format = format;
    x.includeFibres = exportFibres;
    return exportResults(x, control, err);
}

std::vector<std::string> filesUnder(const std::string& dir) {
    std::vector<std::string> out;
    for (const auto& e : fs::recursive_directory_iterator(dir))
        if (e.is_regular_file()) out.push_back(fs::relative(e.path(), dir).generic_string());
    std::sort(out.begin(), out.end());
    return out;
}

} // namespace

// ------------------------------------------------------------------ import ----

TEST(import_obj_material_texture_colours) {
    const std::string d = freshDir("obj_basic");
    writeText(d + "/m.mtl", "newmtl red\nKd 0.9 0.1 0.1\nKs 1 1 1\nmap_Kd -s 1 1 1 tex.png\nmap_Bump bump.png\n");
    writeTexture(d + "/tex.png", 255, 0, 0);
    // A quad (polygon) and a triangle, vertex colours, negative indices, a vn (ignored).
    writeText(d + "/a.obj",
              "mtllib m.mtl\nv 0 0 0 1 0 0\nv 1 0 0 0 1 0\nv 1 1 0 0 0 1\nv 0 1 0 1 1 1\nv 2 0 0 0.5 0.5 0.5\nvt 0 0\nvt 1 0\nvt 1 1\nvt 0 1\nvn 0 0 1\n"
              "usemtl red\nf 1/1 2/2 3/3 4/4\nf -4/-4 5/2 -3/-2\n");
    ImportedMesh m;
    std::string err;
    CHECK(readMeshFile(d + "/a.obj", ImportLimits{}, m, &err));
    CHECK(m.format == "OBJ" && m.pos.size() == 5 && m.tri.size() == 3 && m.faces == 2 && m.polygons == 1);
    CHECK(m.uv.size() == 3 && m.vcolor.size() == 5 && m.vcolor[1].y == 1.0);
    CHECK(m.materials.size() == 1 && m.materials[0].hasKd && m.materials[0].kd.x == 0.9 && !m.materials[0].texture.empty());
    CHECK(m.dependencies.size() == 2);  // m.mtl, tex.png
    bool sawKs = false;
    for (const std::string& i : m.ignored) sawKs |= i.find("Ks") != std::string::npos;
    CHECK(sawKs);
    const ImportAnalysis a = analyzeImport(m);
    CHECK(!a.fatal() && a.colourSource.find("tex.png") != std::string::npos);
    std::vector<std::string> warnings;
    const std::vector<ColourImage> tex = loadImportTextures(m, ImportLimits{}, &warnings);
    CHECK(tex.size() == 1 && tex[0].width == 16 && warnings.empty());
    PreparedMesh p;
    CHECK(prepareImport(m, ImportOptions{}, p, &err));
    const std::vector<Vec3> c = triangleColours(p, tex, {0, 0, 0});
    CHECK(c.size() == 3 && std::fabs(c[0].x - 1.0) < 1e-9 && c[0].y == 0);  // the texture's red (linear)
}

TEST(import_ply_three_encodings) {
    const std::string d = freshDir("ply");
    const Grid g = grid(4, 3, 30, 20);
    writeTexture(d + "/tex.png", 0, 255, 0);
    std::vector<ImportedMesh> ms;
    for (const char* fmt : {"ascii", "binary_little_endian", "binary_big_endian"}) {
        const std::string p = d + "/" + fmt + ".ply";
        writePlyFile(p, g, fmt);
        ImportedMesh m;
        std::string err;
        CHECK(readMeshFile(p, ImportLimits{}, m, &err));
        CHECK(m.pos.size() == g.pos.size() && m.tri.size() == g.tri.size() && m.uv.size() == g.tri.size() && m.vcolor.size() == g.pos.size());
        CHECK(m.materials.size() == 1 && !m.materials[0].texture.empty());
        ms.push_back(std::move(m));
    }
    for (size_t k = 1; k < ms.size(); ++k) {
        CHECK(ms[k].tri == ms[0].tri);
        for (size_t v = 0; v < ms[0].pos.size(); ++v) {
            CHECK(norm(ms[k].pos[v] - ms[0].pos[v]) < 1e-5 && norm(ms[k].vcolor[v] - ms[0].vcolor[v]) < 1e-9);
        }
    }
}

TEST(import_units_and_up_axis) {
    const std::string d = freshDir("units");
    writeText(d + "/a.obj", "v 0 0 0\nv 3 0 0\nv 0 5 0\nv 0 0 2\nf 1 2 3\nf 1 3 4\n");
    ImportedMesh m;
    std::string err;
    CHECK(readMeshFile(d + "/a.obj", ImportLimits{}, m, &err));
    ImportOptions o;
    o.unitMm = 10;  // cm
    o.up = UpAxis::Y;
    PreparedMesh p;
    CHECK(prepareImport(m, o, p, &err));
    Vec3 lo, hi;
    p.mesh.bounds(lo, hi);
    // Y-up file: height = file y extent (5 cm = 50 mm); depth = file z (2 cm); x 3 cm.
    CHECK(std::fabs(hi.z - lo.z - 50) < 1e-9 && std::fabs(hi.x - lo.x - 30) < 1e-9 && std::fabs(hi.y - lo.y - 20) < 1e-9);
    CHECK(std::fabs(lo.z) < 1e-12);
    for (size_t v = 0; v < p.mesh.pos.size(); ++v) {
        // The inverse transform returns the file coordinates (OBJ export in the source's units).
        bool found = false;
        const Vec3 back = toFileCoordinates(p, p.mesh.pos[v]);
        for (const Vec3& q : m.pos) found |= norm(q - back) < 1e-9;
        CHECK(found);
    }
    const ImportAnalysis a = analyzeImport(m);
    CHECK(a.suggestedUnitMm == 25.4 || a.suggestedUnitMm == 10 || a.suggestedUnitMm == 1000);  // a suggestion, by the extent (5)
}

TEST(import_oversize_rejected_before_parse) {
    const std::string d = freshDir("oversize");
    const std::string big = d + "/big.obj";
    {
        std::ofstream f(big, std::ios::binary);
        f << "v 0 0 0\n";
    }
    fs::resize_file(big, (200ull << 20) + 1);  // a real file just over the default limit
    ImportedMesh m;
    std::string err;
    const auto t0 = std::chrono::steady_clock::now();
    CHECK(!readMeshFile(big, ImportLimits{}, m, &err));
    const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    CHECK(err.find("import limit is 200.0 MB of mesh data") != std::string::npos && err.find("decimate tool") != std::string::npos);
    CHECK(ms < 200 && m.pos.empty());  // the size is checked before the file is read
    fs::remove(big);
    // An injected limit applies the same rule.
    writeText(d + "/small.obj", objText(grid(2, 2, 1, 1), "", ""));
    ImportLimits lim;
    lim.maxFileBytes = 10;
    CHECK(!readMeshFile(d + "/small.obj", lim, m, &err) && err.find("import limit") != std::string::npos);
}

// The nodes of a file with names: o and g are one rule, a name without faces is no node, a name used twice is one node, and faces
// before any name have a node of their own.
TEST(import_nodes_from_names) {
    const std::string d = freshDir("nodes_named");
    writeText(d + "/named.obj", "v 0 0 0\nv 1 0 0\nv 0 1 0\nv 0 0 1\n"
                                "f 1 2 3\n"
                                "g default\n"
                                "o eyes\nf 1 2 4\nf 1 3 4\n"
                                "g nose|tip\nf 2 3 4\n"
                                "o eyes\nf 1 2 3\n");
    ImportedMesh m;
    std::string err;
    CHECK(readMeshFile(d + "/named.obj", ImportLimits{}, m, &err));
    finishNodes(m);
    CHECK(m.nodes.size() == 3);
    CHECK(m.nodes[0].name == "(unnamed)" && m.nodes[1].name == "eyes" && m.nodes[2].name == "nose_tip");  // | separates names in the settings
    CHECK(m.nodes[0].triangles == 1 && m.nodes[1].triangles == 3 && m.nodes[2].triangles == 1);
    CHECK((m.nodeOfTriangle == std::vector<int>{0, 1, 1, 2, 1}));
    CHECK(std::fabs(m.nodes[0].area - 0.5) < 1e-12 && m.nodes[2].area > 0);
}

// Without two names a file's nodes are its connected pieces on the WELDED mesh: a PLY repeats a corner per face.
TEST(import_nodes_from_loose_parts) {
    const std::string d = freshDir("nodes_loose");
    writeText(d + "/one_object.obj", "o shell\nv 0 0 0\nv 1 0 0\nv 0 1 0\nv 5 0 0\nv 6 0 0\nv 5 1 0\nf 1 2 3\nf 4 5 6\n");
    ImportedMesh m;
    std::string err;
    CHECK(readMeshFile(d + "/one_object.obj", ImportLimits{}, m, &err));
    finishNodes(m);
    CHECK(m.nodes.size() == 2 && m.nodes[0].name == "Part 1" && m.nodes[1].name == "Part 2");
    CHECK((m.nodeOfTriangle == std::vector<int>{0, 1}));

    writeText(d + "/solid.obj", "o shell\nv 0 0 0\nv 1 0 0\nv 0 1 0\nv 1 1 0\nf 1 2 3\nf 2 4 3\n");
    CHECK(readMeshFile(d + "/solid.obj", ImportLimits{}, m, &err));
    finishNodes(m);
    CHECK(m.nodes.size() == 1 && m.nodes[0].name == "shell" && m.nodes[0].triangles == 2);  // one piece keeps its name

    writeText(d + "/corners.ply", "ply\nformat ascii 1.0\nelement vertex 9\nproperty float x\nproperty float y\nproperty float z\n"
                                  "element face 3\nproperty list uchar int vertex_indices\nend_header\n"
                                  "0 0 0\n1 0 0\n0 1 0\n1 0 0\n1 1 0\n0 1 0\n5 0 0\n6 0 0\n5 1 0\n3 0 1 2\n3 3 4 5\n3 6 7 8\n");
    CHECK(readMeshFile(d + "/corners.ply", ImportLimits{}, m, &err));
    finishNodes(m);
    CHECK(m.nodes.size() == 2 && m.nodes[0].triangles == 2 && m.nodes[1].triangles == 1);
}

TEST(import_header_counts_checked_before_allocation) {
    const std::string d = freshDir("counts");
    writeText(d + "/huge.ply", "ply\nformat binary_little_endian 1.0\nelement vertex 3\nproperty float x\nproperty float y\nproperty float z\n"
                               "element face 1000000000\nproperty list uchar int vertex_indices\nend_header\n");
    ImportedMesh m;
    std::string err;
    CHECK(!readMeshFile(d + "/huge.ply", ImportLimits{}, m, &err));
    CHECK(err.find("1,000,000,000 faces") != std::string::npos);
    // Within the count limits but larger than the file: truncated, before reading.
    writeText(d + "/trunc.ply", "ply\nformat binary_little_endian 1.0\nelement vertex 1000\nproperty float x\nproperty float y\nproperty float z\n"
                                "element face 10\nproperty list uchar int vertex_indices\nend_header\nxyz");
    CHECK(!readMeshFile(d + "/trunc.ply", ImportLimits{}, m, &err) && err.find("truncated") != std::string::npos);
    ImportLimits lim;
    lim.maxVertices = 10;
    writeText(d + "/many.obj", objText(grid(4, 4, 1, 1), "", ""));
    CHECK(!readMeshFile(d + "/many.obj", lim, m, &err) && err.find("more than 10 vertices") != std::string::npos);
}

TEST(import_malformed_fatal_messages) {
    const std::string d = freshDir("malformed");
    std::string err;
    ImportedMesh m;
    writeText(d + "/nan.obj", "v 0 0 0\nv nan 0 0\nv 0 1 0\nf 1 2 3\n");
    CHECK(!readMeshFile(d + "/nan.obj", ImportLimits{}, m, &err) && err.find("line 2") != std::string::npos && err.find("finite") != std::string::npos);
    writeText(d + "/idx.obj", "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 7\n");
    CHECK(!readMeshFile(d + "/idx.obj", ImportLimits{}, m, &err) && err.find("line 4") != std::string::npos && err.find("refers to no vertex") != std::string::npos);
    writeText(d + "/garbage.ply", "this is not a mesh");
    CHECK(!readMeshFile(d + "/garbage.ply", ImportLimits{}, m, &err) && err.find("not a PLY") != std::string::npos);
    writeText(d + "/lines.obj", "v 0 0 0\nv 1 0 0\nl 1 2\n");
    CHECK(!readMeshFile(d + "/lines.obj", ImportLimits{}, m, &err) && err.find("no faces") != std::string::npos);
    writeText(d + "/a.stl", "solid x\n");
    CHECK(!readMeshFile(d + "/a.stl", ImportLimits{}, m, &err) && err.find("only OBJ and PLY") != std::string::npos);
    writeText(d + "/flat.obj", "v 1 1 1\nv 1 1 1\nv 1 1 1\nf 1 2 3\n");
    CHECK(readMeshFile(d + "/flat.obj", ImportLimits{}, m, &err));
    CHECK(analyzeImport(m).fatal());  // zero extent
}

// Repairs: the analysis counts and the prepared mesh's changes agree exactly.
TEST(import_repairs_exact_counts) {
    const std::string d = freshDir("repairs");
    // Two triangles of a square with the diagonal's vertices duplicated (split seam), a zero-area
    // face, a duplicate face, a flipped face sharing an edge, and an unreferenced vertex.
    writeText(d + "/r.obj",
              "v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 0 0\nv 1 1 0\nv 0 1 0\nv 5 5 5\nv 2 0 0\nv 2 1 0\n"
              "f 1 2 3\nf 4 5 6\nf 1 2 2\nf 1 2 3\nf 2 3 9\nf 2 9 8\n");
    ImportedMesh m;
    std::string err;
    CHECK(readMeshFile(d + "/r.obj", ImportLimits{}, m, &err));
    const ImportAnalysis a = analyzeImport(m);
    CHECK(a.degenerate == 1 && a.duplicateFaces == 1 && a.duplicateVertices == 2 && a.unreferenced == 1);
    CHECK(a.inconsistentEdges >= 1);
    PreparedMesh p;
    CHECK(prepareImport(m, ImportOptions{}, p, &err));
    CHECK(p.mesh.tri.size() == 4 && p.mesh.pos.size() == 6);
    std::string all;
    for (const std::string& c : p.changes) all += c + "\n";
    CHECK(all.find("merged 2 coincident") != std::string::npos && all.find("removed 1 zero-area") != std::string::npos &&
          all.find("removed 1 duplicate") != std::string::npos && all.find("1 vertices that belong to no face") != std::string::npos &&
          all.find("flipped") != std::string::npos);
    // Consistent afterwards: every interior edge used once in each direction.
    std::map<std::pair<int, int>, int> dir;
    for (const auto& t : p.mesh.tri)
        for (int k = 0; k < 3; ++k) ++dir[{t[k], t[(k + 1) % 3]}];
    for (const auto& [e, n] : dir) CHECK(n == 1);
}

TEST(import_missing_texture_falls_back) {
    const std::string d = freshDir("missing_tex");
    writeText(d + "/m.mtl", "newmtl a\nKd 0.1 0.2 0.3\nmap_Kd gone.png\n");
    writeText(d + "/a.obj", objText(grid(2, 2, 10, 10), "m.mtl", "a"));
    ImportedMesh m;
    std::string err;
    CHECK(readMeshFile(d + "/a.obj", ImportLimits{}, m, &err));
    CHECK(m.materials[0].textureMissing && m.missingDependencies.size() == 1);
    const ImportAnalysis a = analyzeImport(m);
    bool warned = false;
    for (const ImportIssue& i : a.issues) warned |= i.id == "missing_file";
    CHECK(warned && a.colourSource == "material colours (Kd)");
    PreparedMesh p;
    CHECK(prepareImport(m, ImportOptions{}, p, &err));
    const std::vector<Vec3> c = triangleColours(p, std::vector<ColourImage>{}, {1, 1, 1});
    CHECK(c[0].x == 0.1 && c[0].z == 0.3);
}

// The identity covers the MTL and the images: changing only one of them is another import.
TEST(import_identity_dependencies) {
    const std::string d = freshDir("identity");
    writeText(d + "/m.mtl", "newmtl a\nKd 0.1 0.2 0.3\nmap_Kd t.png\n");
    writeTexture(d + "/t.png", 1, 2, 3);
    writeText(d + "/a.obj", objText(grid(2, 2, 10, 10), "m.mtl", "a"));
    auto ident = [&] {
        ImportedMesh m;
        std::string err;
        CHECK(readMeshFile(d + "/a.obj", ImportLimits{}, m, &err));
        return importIdentity(m, ImportOptions{});
    };
    const uint64_t i0 = ident();
    writeText(d + "/unrelated.txt", "x");
    CHECK(ident() == i0);
    writeText(d + "/m.mtl", "newmtl a\nKd 0.1 0.2 0.4\nmap_Kd t.png\n");  // the MTL only
    const uint64_t i1 = ident();
    CHECK(i1 != i0);
    writeTexture(d + "/t.png", 9, 2, 3);  // the image only
    const uint64_t i2 = ident();
    CHECK(i2 != i1);
    fs::remove(d + "/t.png");  // the image gone
    const uint64_t i3 = ident();
    CHECK(i3 != i2);
    ImportOptions o;
    o.unitMm = 10;
    ImportedMesh m;
    std::string err;
    CHECK(readMeshFile(d + "/a.obj", ImportLimits{}, m, &err));
    CHECK(importIdentity(m, o) != i3);  // the options too
}

// The prepared entry: self-contained (texture copied), reads back exactly; the source is untouched.
TEST(import_store_roundtrip) {
    const std::string d = freshDir("store"), root = freshDir("store_cache");
    ImportCacheEnv env(root);
    bool ok = false;
    uint64_t id = 0;
    const SceneSettings s = importPlane(d, root, &id);
    const uint64_t srcHash = fnv1aFile(d + "/plane.obj", &ok);
    CHECK(importStored(root, id));
    PreparedMesh p;
    std::string src, err;
    CHECK(loadImport(importCacheDir(root, id), p, &src, &err));
    CHECK(src == d + "/plane.obj" && p.materials.size() == 1 && p.materials[0].texture.find(importCacheDir(root, id)) == 0);
    CHECK(fs::exists(p.materials[0].texture));
    CHECK(fnv1aFile(d + "/plane.obj", &ok) == srcHash);
    SourceScene g;
    CHECK(loadImportedScene(importCacheDir(root, id), g, &err));
    CHECK(g.imported && g.prims.size() == 1 && g.prims[0].albedo && g.prims[0].triColor.size() == g.prims[0].welded.tri.size());
    // The object's base colour: the texture's (200, 100, 50) sRGB, as linear.
    CHECK(std::fabs(g.prims[0].baseColor.x - std::pow((200 / 255.0 + 0.055) / 1.055, 2.4)) < 1e-9);
    (void)s;
}

// ----------------------------------------------------------- results ----

// Knit (Draft) -> blocked; Final -> eligible; fibres consume the knit without knitting; the
// export writes the reviewed fibres; an edit invalidates downstream.
// Fibre job progress (FibreProgress, job_steps.h): the phases over all parts publish a never-decreasing
// value, below 100 % until finish(); the display lines' stage keeps the cap.
TEST(fibre_progress_aggregator) {
    struct Published : JobReporter {
        std::vector<int> values;
        void progress(const std::string&, int done, int total) override {
            CHECK(total == FibreProgress::kScale);
            values.push_back(done);  // called under the aggregator's lock: in publication order
        }
        void partial(std::shared_ptr<BuiltScene>) override {}
    } reporter;
    std::vector<int>& pub = reporter.values;
    const JobControl control{0, &reporter};
    FibreProgress p(control);
    const int n = 8;
    for (int k = 0; k <= n; ++k) p.phase("fibres: preparing parts", k, n);
    for (int k = 0; k <= n; ++k) p.phase("fibres: generating", k, n);
    p.phase("fibres: preparing parts", 0, n);  // an earlier phase again: never backwards
    for (int k = 0; k <= n; ++k) p.phase("fibres: finishing parts", k, n);
    bool monotonic = true;
    for (size_t i = 1; i < pub.size(); ++i) monotonic = monotonic && pub[i] >= pub[i - 1];
    CHECK(monotonic && pub.back() == FibreProgress::kBeforeResult && p.published() < FibreProgress::kScale);
    p.finishing("building the display lines");
    CHECK(p.published() == FibreProgress::kBeforeResult);
    p.finish();
    CHECK(p.published() == FibreProgress::kScale && pub.back() == FibreProgress::kScale);
}

// The hand in the real pipeline, on the device, with the chart's projection binding: placement alone
// moves stitches within the fabric, and the solver turns tension into height.
TEST(workflow_hand_solver_adds_height) {
    const std::string d = freshDir("hand"), root = freshDir("hand_cache");
    ImportCacheEnv env(root);
    SceneSettings s = importPlane(d, root);
    s.hand.unevenness = 30;
    s.hand.loopRate = 1;
    s.hand.loopHeight = 1;
    std::string err;
    auto draft = run(BuiltScene::Kind::Knit, s, &err);
    CHECK(draft && draft->knit);
    applyFullProfile(s, *sceneEntryFor(s));
    auto fin = run(BuiltScene::Kind::Knit, s, &err);
    CHECK(fin && fin->knit);
    const CurveSet& placed = draft->knit->scene->parts[0].knit->t.curves;
    const CurveSet& settled = fin->knit->scene->parts[0].knit->t.curves;
    const StitchSpread before = measureStitchSpread(placed), after = measureStitchSpread(settled);
    std::printf("    height %.3f -> %.3f mm, shape %.3f -> %.3f mm\n", before.height, after.height, before.shape, after.shape);
    CHECK(after.height > 2 * before.height);
    const ScenePart& part = fin->knit->scene->parts[0];
    const RaisedLoops loops = measureRaisedLoops(settled, drawStitchTensions(settled.topo, handRequest(s, part.name)));
    std::printf("    %d of %d raised loops standing\n", loops.standing, loops.asked);
    CHECK(loops.asked > 0 && loops.standing >= 0.9 * loops.asked);
    CHECK(fin->knit->scene->parts[0].knit->after.localSignature.pass == fin->knit->scene->parts[0].knit->before.localSignature.pass);
}

// Thickness is derived from the stitches, so a saved knit keeps its bytes and gets it back on load.
TEST(workflow_hand_thickness_survives_a_reload) {
    const std::string d = freshDir("handsave"), root = freshDir("handsave_cache");
    ImportCacheEnv env(root);
    SceneSettings s = importPlane(d, root);
    s.hand.spotRate = 3;
    s.hand.spotChange = 30;
    std::string err;
    auto b = run(BuiltScene::Kind::Knit, s, &err);
    CHECK(b && b->knit);
    const CurveSet& knitted = b->knit->scene->parts[0].knit->t.curves;
    CHECK(!knitted.yarns[0].thickness.empty());
    Scene t;
    CHECK(restoreKnit(serializeKnit(*b->knit->scene), *sceneEntryFor(s), s, t, &err));
    for (size_t i = 0; i < t.parts.size(); ++i) {
        if (!t.parts[i].knit) continue;
        CHECK(curvesDigest(t.parts[i].knit->t.curves) == curvesDigest(b->knit->scene->parts[i].knit->t.curves));
    }
}

TEST(workflow_plane_knit_fibres_export) {
    const std::string d = freshDir("flow"), root = freshDir("flow_cache"), out = freshDir("flow_out");
    ImportCacheEnv env(root);
    SceneSettings s = importPlane(d, root);
    std::string err;
    auto draft = run(BuiltScene::Kind::Knit, s, &err);
    CHECK(draft && draft->knit && draft->knit->computational && !draft->knit->final);
    // The chart part and the gap patch that covers what it left open: every part of the import is knitted.
    CHECK(draft->knit->parts == 2 && draft->knit->plain == 0);
    ExportEligibility el = exportEligibility(s, draft->knit.get(), nullptr, FibreChoice::LeftOut);
    CHECK(el.state == ExportEligibility::State::Blocked && el.fix == "Knit at Final quality");
    // Final.
    applyFullProfile(s, *sceneEntryFor(s));
    auto fin = run(BuiltScene::Kind::Knit, s, &err);
    CHECK(fin && fin->knit && fin->knit->final);
    std::printf("    plane: %d stitches, verdict %s\n", fin->parts[0].stitches, toString(fin->parts[0].verdict));
    // Fibres: a job of the GIVEN knit; its stages never knit. With Soft wool's late taper, whose thin
    // ends add vertices the Export step's forecast counts.
    s.fibres = true;
    s.fibreDensity = 20;
    s.tailTaper = 0.2;
    s.tailTaperTip = 0.3;
    StageLog log;
    auto fib = run(BuiltScene::Kind::Fibres, s, &err, &log, fin->knit);
    CHECK(fib && fib->fibreResult && fib->fibreResult->knitId == fin->knit->id && fib->fibreResult->fibres > 0);
    for (const std::string& st : log.stages) CHECK(st.find("knit") == std::string::npos || st.find("restored") != std::string::npos);
    {
        // Work progress (FibreProgress, job_steps.h): never backwards, 100 % only at the end (the result).
        int last = 0, full = 0;
        bool monotonic = true;
        for (const auto& [done, total] : log.values) {
            if (total != FibreProgress::kScale) continue;
            monotonic = monotonic && done >= last;
            last = done;
            full += done == total;
        }
        CHECK(monotonic && last == FibreProgress::kScale && full == 1 && log.values.back().first == FibreProgress::kScale);
    }
    el = exportEligibility(s, fin->knit.get(), fib->fibreResult.get(), FibreChoice::Included);
    CHECK(el.state != ExportEligibility::State::Blocked);
    const bool diagnostic = el.state == ExportEligibility::State::Diagnostic;
    // A material edit keeps everything current; a fibre edit blocks the fibres; a geometry edit the knit.
    SceneSettings m = s;
    m.yarnWhite.color = {0.1, 0.1, 0.1};
    CHECK(exportEligibility(m, fin->knit.get(), fib->fibreResult.get(), FibreChoice::Included).state == el.state);
    m = s;
    m.random.lengthMax = 1.9;
    CHECK(exportEligibility(m, fin->knit.get(), fib->fibreResult.get(), FibreChoice::Included).state == ExportEligibility::State::Blocked);
    CHECK(exportEligibility(m, fin->knit.get(), fib->fibreResult.get(), FibreChoice::LeftOut).state == el.state);  // fibres not included
    m = s;
    m.plies = 3;  // the surface (explicit plies): the knit result is out of date
    CHECK(exportEligibility(m, fin->knit.get(), fib->fibreResult.get(), FibreChoice::Included).state == ExportEligibility::State::Blocked);
    // Export (PPT) with the reviewed fibres, a wool profile's translucency and one yarn per object, as
    // the teddy's (material edits).
    s.yarnWhite.translucency = 0.3;
    s.yarnColourSource = "object";
    const std::string target = out + "/plane_knit";
    auto ex = run(BuiltScene::Kind::Export, s, &err, nullptr, fin->knit, fib->fibreResult, target, "ppt", true);
    CHECK(ex && !ex->exportedTo.empty());
    CHECK(fs::exists(target + "/scene.json") && fs::exists(target + "/report.json"));
    {
        // Every yarn entry of the scene (per object or colour, and its ply-mapped "_yarn" variant) is a
        // copy of material.yarn_white: each carries the translucency, as the renderer's fabric block.
        std::ifstream sj(target + "/scene.json");
        std::string line;
        int yarns = 0, perObject = 0, translucent = 0;
        while (std::getline(sj, line)) {
            const bool object = line.find("{\"id\": \"yarn_") != std::string::npos;
            if (!(object || line.find("{\"id\": \"wool_") != std::string::npos) || line.find("\"type\": \"fabric\"") == std::string::npos) continue;
            ++yarns;
            perObject += object;
            translucent += line.find("\"translucency\": {\"transmit_share\": 0.3}") != std::string::npos;
        }
        std::printf("    yarn materials: %d (%d per object), translucent %d\n", yarns, perObject, translucent);
        CHECK(perObject >= 1 && translucent == yarns);
    }
    bool strands = false;
    for (const std::string& f : filesUnder(target)) strands |= f.find(".pptstrands") != std::string::npos;
    CHECK(strands);
    {
        // The knitted part's source mesh is not written: the meshes are its yarn (and the ground),
        // none in an unknitted part's material (plain_ / mean_ / tex).
        std::ifstream sj(target + "/scene.json");
        std::string line, section, materials;
        int meshes = 0, source = 0;
        while (std::getline(sj, line)) {
            if (line.find("\"meshes\": [") != std::string::npos) section = "meshes";
            else if (line.find("\": [") != std::string::npos) section.clear();
            else if (section == "meshes" && line.find("\"material\": \"") != std::string::npos) {
                const size_t a = line.find("\"material\": \"") + 13;
                const std::string mat = line.substr(a, line.find('"', a) - a);
                materials += " " + mat;
                ++meshes;
                source += mat.rfind("plain_", 0) == 0 || mat.rfind("mean_", 0) == 0 || mat.rfind("tex", 0) == 0;
            }
        }
        std::printf("    exported meshes:%s\n", materials.c_str());
        CHECK(meshes >= 1 && source == 0);
    }
    CHECK(ex->exportLog.find("not generated again") != std::string::npos);
    {
        // The fibre result's forecast of the exported segments is the export's count, taper included.
        long long forecast = 0, grown = 0;
        for (const ObjectSegments& o : fib->fibreResult->objectSegments) forecast += o.exported, grown += o.grown;
        const size_t at = ex->exportLog.find(" exported\n");
        const size_t from = ex->exportLog.rfind("-> ", at);
        const long long exported = at != std::string::npos && from != std::string::npos ? std::stoll(ex->exportLog.substr(from + 3, at - from - 3)) : -1;
        std::printf("    fibre segments: %lld grown, forecast %lld exported, export %lld\n", grown, forecast, exported);
        CHECK(fib->fibreResult->objectSegments.size() == 1 && forecast == exported && forecast > grown);
        CHECK(ex->exportLog.find("renderer GPU memory, render phase: about ") != std::string::npos);
    }
    // The same scene through the CLI path (the fibres grown again with knit_check's diagnostics, then
    // exported): byte-identical files.
    SceneFibres fresh;
    CHECK(grown(fresh, growSceneFibres(*sceneEntryFor(s), *fin->knit->scene, s, FibreRunOptions{}, &err)));
    ExportOptions eo;
    eo.folder = out + "/cli";
    eo.status = diagnostic ? ExportEligibility::State::Diagnostic : ExportEligibility::State::Ready;
    std::string clog;
    CHECK(exportScene({*sceneEntryFor(s), *fin->knit->scene, s, &fresh}, eo, &clog, &err));
    const std::vector<std::string> a = filesUnder(target), b = filesUnder(eo.folder);
    CHECK(a == b);
    for (const std::string& f : a)
        if (f != kPublishMarker) CHECK(readText(target + "/" + f) == readText(eo.folder + "/" + f));
}

// The Export step forecasts the tube triangles as 2 sides per ring spacing along BuiltScene::yarnLength.
// That length is the tubed curves': with plies the plies' (their helix included), so the export's yarn
// tube triangles match the forecast with and without plies, and multiplying by the plies again would not.
TEST(export_tube_triangles_follow_the_tubed_length) {
    const std::string d = freshDir("tubelen"), root = freshDir("tubelen_cache"), out = freshDir("tubelen_out");
    ImportCacheEnv env(root);
    SceneSettings s = importPlane(d, root);
    applyFullProfile(s, *sceneEntryFor(s));
    s.filling = false;
    for (int plies : {1, 3}) {
        s.plies = plies;
        if (plies > 1) s.textureProfile = "none";  // as the viewer exports plies: round tubes
        std::string err;
        auto fin = run(BuiltScene::Kind::Knit, s, &err);
        CHECK(fin && fin->knit && fin->yarnLength > 0);
        if (!fin || !fin->knit) return;
        const std::string target = out + "/plies" + std::to_string(plies);
        auto ex = run(BuiltScene::Kind::Export, s, &err, nullptr, fin->knit, nullptr, target, "ppt", false);
        CHECK(ex && !ex->exportedTo.empty());
        if (!ex || ex->exportedTo.empty()) return;
        // The yarn tubes' triangles, from the PLY headers of the meshes in a yarn material.
        long long triangles = 0;
        std::ifstream sj(target + "/scene.json");
        for (std::string line; std::getline(sj, line);) {
            const size_t path = line.find("\"path\": \"meshes/"), material = line.find("\"material\": \"");
            if (path == std::string::npos || material == std::string::npos) continue;
            const std::string mat = line.substr(material + 13, line.find('"', material + 13) - material - 13);
            if (mat.rfind("yarn", 0) != 0 && mat.rfind("wool", 0) != 0) continue;
            std::ifstream ply(target + "/" + line.substr(path + 9, line.find('"', path + 9) - path - 9), std::ios::binary);
            for (std::string h; std::getline(ply, h) && h != "end_header";)
                if (h.rfind("element face ", 0) == 0) triangles += std::stoll(h.substr(13));
        }
        const double perMm = plies > 1 ? s.plyTubeSides / s.plyTubeRing : s.tubeSides / s.tubeRing;
        const double forecast = 2.0 * perMm * fin->yarnLength;
        std::printf("    plies %d: %.0f mm tubed, forecast %.0f triangles, exported %lld\n", plies, fin->yarnLength, forecast, triangles);
        CHECK(triangles > 0 && std::fabs(double(triangles) - forecast) < 0.05 * forecast);
    }
}

// A fibre job reports where its time goes: per-part progress (done / total), the resolved
// backend, the phase wall times, per-part stage times and (GPU) the generator's run report.
// The stand-ins of a running knit (the parts not knitted yet) lie behind the fabric: each vertex
// is the source vertex moved by fillingOffset along -N (the export filling's depth), and no vertex
// of the finished yarn reaches that depth (the knit is centred on the surface, controls at +-A D).
TEST(partial_base_inside_fabric) {
    const std::string d = freshDir("standin"), root = freshDir("standin_cache");
    ImportCacheEnv env(root);
    SceneSettings s = importPlane(d, root);
    struct First : JobReporter {
        std::mutex m;
        std::shared_ptr<BuiltScene> first;
        void progress(const std::string&, int, int) override {}
        void partial(std::shared_ptr<BuiltScene> p) override {
            std::lock_guard<std::mutex> lk(m);
            if (!first) first = std::move(p);
        }
    } sink;
    std::string err;
    auto b = run(BuiltScene::Kind::Knit, s, &err, &sink);
    CHECK(b && b->knit && sink.first);
    if (!b || !sink.first) return;
    const double inset = fillingOffset(s, s.diameter);
    const double expect = s.depthAmplitude * s.diameter + 0.5 * finalSurface(s, s.diameter).renderedDiameter() + s.clearance + s.fillingMargin;
    CHECK(std::fabs(inset - expect) < 1e-12);
    const DisplayMesh *standIn = nullptr, *source = nullptr;
    for (const DisplayMesh& m : sink.first->meshes)
        if (m.kind == DisplayMesh::Kind::Base) standIn = &m;
    for (const DisplayMesh& m : b->meshes)
        if (m.kind == DisplayMesh::Kind::Base) source = &m;  // the knitted part's "(source)" surface, at offset 0
    CHECK(standIn && source && standIn->pos.size() == source->pos.size() && !source->pos.empty());
    if (!standIn || !source || standIn->pos.size() != source->pos.size() || source->pos.empty()) return;
    double moved = 0;
    for (size_t i = 0; i < source->pos.size(); ++i)
        moved = std::max(moved, double(std::fabs(standIn->pos[i] - (source->pos[i] - source->nrm[i] * float(inset)))));
    // The plane's yarn: its deepest point below the source surface, along the surface normal.
    const Vec3 n{source->nrm[0], source->nrm[1], source->nrm[2]};
    const double h0 = dot(n, Vec3{source->pos[0], source->pos[1], source->pos[2]});
    double deepest = 1e30;
    for (const DisplayMesh& m : b->meshes)
        if (m.kind == DisplayMesh::Kind::Yarn)
            for (size_t i = 0; i + 2 < m.pos.size(); i += 3) deepest = std::min(deepest, dot(n, Vec3{m.pos[i], m.pos[i + 1], m.pos[i + 2]}) - h0);
    std::printf("    stand-in inset %.3f mm (A D %.3f + r_R %.3f + c %.2f + margin %.2f); max vertex error %.2e mm; deepest yarn %.3f mm\n",
                inset, s.depthAmplitude * s.diameter, 0.5 * finalSurface(s, s.diameter).renderedDiameter(), s.clearance, s.fillingMargin,
                moved, deepest);
    CHECK(moved < 1e-4);  // float display positions
    CHECK(deepest > -inset && deepest < 0);
}

TEST(fibre_progress_and_profile_lines) {
    const std::string d = freshDir("fprog"), root = freshDir("fprog_cache");
    ImportCacheEnv env(root);
    SceneSettings s = importPlane(d, root);
    std::string err;
    auto knit = run(BuiltScene::Kind::Knit, s, &err);
    CHECK(knit && knit->knit);
    struct Totals : JobReporter {
        std::mutex m;
        std::vector<std::pair<std::string, int>> seen;  // stage, total
        void progress(const std::string& stage, int, int total) override {
            std::lock_guard<std::mutex> lk(m);
            seen.push_back({stage, total});
        }
        void partial(std::shared_ptr<BuiltScene>) override {}
    } sink;
    s.fibres = true;
    s.fibreDensity = 20;
    auto fib = run(BuiltScene::Kind::Fibres, s, &err, &sink, knit->knit);
    CHECK(fib && fib->fibreResult);
    bool perPart = false;
    // Per part, in the work scale of fibre_progress.h (the part in the label).
    for (const auto& [stage, total] : sink.seen) perPart |= stage == "fibres: finishing parts (part 1 of 2)" && total == FibreProgress::kScale;
    CHECK(perPart);
    const std::string& log = fib->fibreResult->log;
    std::printf("%s", log.c_str());
    CHECK(log.find("generation on the GPU") != std::string::npos);
    CHECK(log.find("fibres: wall time: prepare") != std::string::npos);
    CHECK(log.find("profile: fibres:") != std::string::npos);
    CHECK(log.find("fibres: verify the knit") != std::string::npos);
    CHECK(log.find("GPU fibres:") != std::string::npos);
    CHECK(log.find("budget-ordered generation") != std::string::npos);
}

// The viewer's fibre options leave out diagnostics only (measured angles; on the GPU the audit
// of the omitted root check, the hidden-region states and an eager tube index): every part's
// fibres - candidate ids, populations, polylines, radii, keys, exported strands - are identical
// to the knit_check options.
TEST(fibre_options_do_not_change_strands) {
    const std::string d = freshDir("fopt"), root = freshDir("fopt_cache");
    ImportCacheEnv env(root);
    SceneSettings s = importPlane(d, root);
    std::string err;
    auto knit = run(BuiltScene::Kind::Knit, s, &err);
    CHECK(knit && knit->knit && knit->knit->scene);
    if (!knit_test::useGpu() || !computeCaps().device) SKIP("no usable CUDA GPU");
    s.fibres = true;
    s.fibreDensity = 40;
    SceneSettings b = s;
    FibreRunOptions full;  // knit_check --fiber-check: the checks, the angles and the CPU audit
    full.fibreCheck = true;
    FibreRunOptions viewer;  // the viewer's and knit_check's default: the light checks
    SceneFibres a, v;
    CHECK(grown(a, growWoolFibres(*knit->knit->scene, b, full, &err)));
    CHECK(grown(v, growWoolFibres(*knit->knit->scene, b, viewer, &err)));
    CHECK(a.parts.size() == v.parts.size() && !a.parts.empty());
    size_t fibres = 0;
    for (size_t p = 0; p < a.parts.size(); ++p) {
        const PartFibres &x = a.parts[p], &y = v.parts[p];
        CHECK(x.fs.fibres.size() == y.fs.fibres.size());
        for (size_t i = 0; i < x.fs.fibres.size(); ++i) {
            const Fibre &f = x.fs.fibres[i], &g = y.fs.fibres[i];
            CHECK(f.candidate == g.candidate && f.yarnId == g.yarnId && f.population == g.population && f.key == g.key);
            CHECK(f.radius0 == g.radius0 && f.radius1 == g.radius1 && f.pointCount == g.pointCount);
            CHECK(fibrePoints(x.fs, f) == fibrePoints(y.fs, g));
        }
        fibres += x.fs.fibres.size();
    }
    CHECK(fibres > 0);
    CHECK(a.log.find("angle to the yarn at emergence") != std::string::npos);
    CHECK(v.log.find("angles not measured") != std::string::npos);
    std::printf("    %zu fibres identical with and without the diagnostics\n", fibres);
}

// Budget order (fibre_budget.h): with a segment budget that really cuts (a third of the full run's
// segments), the budget-ordered fibres with float(key) < q are exactly the full run's - candidate
// ids, polylines and radii - so the export keeps the same fibres. The key windows [0, a) + [a, 2)
// of one part reproduce its full run.
TEST(fibre_budget_matches_full_per_backend) {
    const std::string d = freshDir("fbud"), root = freshDir("fbud_cache");
    ImportCacheEnv env(root);
    SceneSettings s = importPlane(d, root);
    std::string err;
    auto knit = run(BuiltScene::Kind::Knit, s, &err);
    CHECK(knit && knit->knit && knit->knit->scene);
    const Scene& t = *knit->knit->scene;
    if (!knit_test::useGpu() || !computeCaps().device) SKIP("no usable CUDA GPU");
    s.fibres = true;
    s.fibreDensity = 60;
    auto same = [](const FibreSet& a, const Fibre& f, const FibreSet& b, const Fibre& g) {
        return f.candidate == g.candidate && f.yarnId == g.yarnId && f.key == g.key && f.radius0 == g.radius0 && f.radius1 == g.radius1 &&
               fibrePoints(a, f) == fibrePoints(b, g);
    };
    FibreRunOptions eo;
    SceneSettings b = s;
    b.fibreGeneration = "full";
    SceneFibres full;
    CHECK(grown(full, growWoolFibres(t, b, eo, &err)));
    long long segs = 0;
    for (const PartFibres& p : full.parts) segs += segmentCount(p.fs);
    CHECK(segs > 300);
    b.fibreGeneration = "budget";
    b.segmentBudget = segs / 3;
    SceneFibres bud;
    CHECK(grown(bud, growWoolFibres(t, b, eo, &err)));
    // The budget order's threshold is the one full generation gives at the same budget.
    std::vector<std::pair<float, int>> fullKeySegs;
    for (const PartFibres& p : full.parts) addKeySegments(p.fs, fullKeySegs);
    // An import is one object: every part shares its threshold.
    CHECK(!bud.parts.empty());
    if (bud.parts.empty()) return;
    const float q = bud.parts[0].keyBelow;
    for (const PartFibres& p : bud.parts) CHECK(p.keyBelow == q);
    CHECK(bud.budgetMode && q < 1.0f && budgetThreshold(fullKeySegs, b.segmentBudget) == q);
    size_t kept = 0;
    for (size_t p = 0; p < full.parts.size(); ++p) {
        std::vector<const Fibre*> a, c;
        for (const Fibre& f : full.parts[p].fs.fibres)
            if (float(f.key) < q) a.push_back(&f);
        for (const Fibre& f : bud.parts[p].fs.fibres)
            if (float(f.key) < q) c.push_back(&f);
        CHECK(a.size() == c.size());
        for (size_t i = 0; i < a.size(); ++i) CHECK(same(full.parts[p].fs, *a[i], bud.parts[p].fs, *c[i]));
        kept += a.size();
    }
    CHECK(kept > 0);
    std::printf("    budget %lld of %lld segments, q %.6f: %zu fibres below q identical to the full run\n", b.segmentBudget, segs, double(q),
                kept);
    // Key windows.
    FibreParams fp;
    CHECK(fibreParamsFromSettings(s, fp, &err));
    const ScenePart* part = nullptr;
    for (const ScenePart& r : t.parts)
        if (r.hasStitches()) part = &r;
    CHECK(part != nullptr);
    auto jg = PartFibreJob::prepare({*part, 0, s, fp, *t.solids, 0.2, FibreRunOptions{}}, &err);
    CHECK(jg != nullptr);
    FibreSet whole, lo, hi;
    CHECK(jg->generate(0.0f, 2.0f, whole, &err));
    const float a = 0.37f;
    CHECK(jg->generate(0.0f, a, lo, &err));
    CHECK(jg->generate(a, 2.0f, hi, &err));
    CHECK(lo.rep.candidates + hi.rep.candidates == whole.rep.candidates);
    CHECK(lo.fibres.size() + hi.fibres.size() == whole.fibres.size() && !lo.fibres.empty() && !hi.fibres.empty());
    for (const Fibre& f : lo.fibres) CHECK(float(f.key) < a);
    for (const Fibre& f : hi.fibres) CHECK(float(f.key) >= a);
    std::map<std::pair<int, int>, const Fibre*> byCand;
    for (const Fibre& f : whole.fibres) byCand[{f.yarnId, f.candidate}] = &f;
    for (const FibreSet* w : {&lo, &hi})
        for (const Fibre& f : w->fibres) {
            auto it = byCand.find({f.yarnId, f.candidate});
            CHECK(it != byCand.end() && same(*w, f, whole, *it->second));
        }
    std::printf("    gpu windows [0, %.2f) + [%.2f, 2): %zu + %zu fibres = the full run's %zu\n", double(a), double(a), lo.fibres.size(),
                hi.fibres.size(), whole.fibres.size());
}

// The Fur workflow's keys (scene_settings.h pile.*): no pile key is an input of a Knit stage, whatever it holds;
// the strand mesh level changes the Fibres and Export keys, never the geometry's; the fur key never depends on a
// knitting key or a material, and changes with a pile key or the asset; pile.* round-trips; the former selectors
// are refused as settings; pileSettingsValid names the key out of bounds.
TEST(pile_keys_and_fur_key) {
    const char* scenes[] = {"teddy", "wool-swatch", "import"};
    const Stage stages[] = {Stage::Load, Stage::Knit, Stage::Surface, Stage::Fibres, Stage::Export};
    for (const char* scene : scenes)
        for (const Stage stage : stages) {
            SceneSettings s;
            s.scene = scene;
            const uint64_t key = stageKey(s, stage, 0x1234);
            // Another strand mesh level: a mesh-only change, the same load, knit and surface keys.
            SceneSettings medium = s;
            medium.tubeSides = kExportSides;
            medium.tubeRing = kExportRing;
            medium.plyTubeSides = kExportPlySides;
            medium.plyTubeRing = kExportPlyRing;
            medium.segmentBudget = 14000000;
            CHECK((stageKey(medium, stage, 0x1234) == key) == (static_cast<int>(stage) < static_cast<int>(Stage::Fibres)));
            // The pile keys do not touch Knit keys.
            s.pile.density = 300;
            s.pile.curl = 1.2;
            s.pileMaterial.color = {0.1, 0.2, 0.3};
            CHECK(stageKey(s, stage, 0x1234) == key);
        }
    // Fur: the workflow changes the fibre key; knitting and material edits do not change the fur key, pile edits do.
    SceneSettings knitted;
    knitted.scene = "teddy";
    SceneSettings f = knitted;
    f.workflow = "fur";
    CHECK(fibreResultKey(f) != fibreResultKey(knitted));
    const uint64_t k0 = furKey(f, 0x1234);
    SceneSettings g = f;
    g.diameter = 0.9;
    g.plies = 3;
    g.pattern = "garter";
    g.pileMaterial.color = {0.5, 0.5, 0.5};
    g.pileColourSource = "single";
    CHECK(furKey(g, 0x1234) == k0);
    g.pile.length = 0.8;
    CHECK(furKey(g, 0x1234) != k0);
    CHECK(furKey(f, 0x1235) != k0);
    // Round trip of the pile keys.
    g.pile.coverage = "mask";
    g.pile.maskPath = "C:/masks/a b.png";
    g.pile.excludedParts = "eyes_c0,nose_c1";
    g.pile.fibresEnabled = false;
    g.pile.varDirection = 0.25;
    g.pile.seed = 17;
    SceneSettings r;
    std::string err;
    CHECK(readSettings(writeSettings(g), r, &err));
    CHECK(r == g && workflowOf(r) == Workflow::Fur);
    CHECK(!setSetting(r, "fibres.profile", "fur", &err) && err.find("workflow") != std::string::npos);
    CHECK(!setSetting(r, "knit.enabled", "false", &err) && err.find("workflow") != std::string::npos);
    CHECK(!setSetting(r, "workflow", "hair", &err));
    // Bounds: the defaults pass; a value out of range, and a loop no wider than its fibre, are refused by key.
    CHECK(pileSettingsValid(f, &err));
    SceneSettings bad = f;
    bad.pile.density = 0;
    CHECK(!pileSettingsValid(bad, &err) && err.find("pile.density") != std::string::npos);
    bad = f;
    bad.pile.curl = 0.4;  // smallest loop 0.2 mm, thickest fibre radius 0.4 mm
    CHECK(!pileSettingsValid(bad, &err) && err.find("pile.curl") != std::string::npos);
    std::printf("    %s\n", err.c_str());
    bad = f;
    bad.pile.bend = 0;
    CHECK(!pileSettingsValid(bad, &err) && err.find("pile.bend") != std::string::npos);
}

// ------------------------------------------------------------------- fur ----

namespace {
TriMesh furSphere(double R, int lat, int lon) {
    TriMesh m;
    for (int i = 0; i <= lat; ++i)
        for (int j = 0; j < lon; ++j) {
            const double th = kPi * i / lat, ph = 2 * kPi * j / lon;
            m.pos.push_back({R * std::sin(th) * std::cos(ph), R * std::sin(th) * std::sin(ph), R * std::cos(th)});
        }
    for (int i = 0; i < lat; ++i)
        for (int j = 0; j < lon; ++j) {
            const int a = i * lon + j, b = i * lon + (j + 1) % lon, c = a + lon, e = b + lon;
            if (i > 0) m.tri.push_back({a, c, b});
            if (i + 1 < lat) m.tri.push_back({b, c, e});
        }
    m.computeNormals();
    return m;
}
// A vertical plane (normal +x) of w x h mm, centred at the origin.
TriMesh furWall(double w, double h, int n) {
    TriMesh m;
    for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= n; ++j) m.pos.push_back({0, -w / 2 + w * j / n, -h / 2 + h * i / n});
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            const int a = i * (n + 1) + j, b = a + 1, c = a + n + 1, d = c + 1;
            m.tri.push_back({a, b, d});
            m.tri.push_back({a, d, c});
        }
    m.computeNormals();
    return m;
}
} // namespace

namespace {
bool sameDisc(const FurDisc& a, const FurDisc& b) {
    return a.root == b.root && a.normal == b.normal && a.d == b.d && a.c == b.c && a.radius == b.radius && a.turns == b.turns &&
           a.rise == b.rise && a.rootRadius == b.rootRadius && a.tipRadius == b.tipRadius && a.straightLength == b.straightLength &&
           a.tri == b.tri && a.index == b.index && a.budgetKey == b.budgetKey && a.layer == b.layer;
}
bool sameField(const FurField& a, const FurField& b) {
    if (a.discs.size() != b.discs.size()) return false;
    for (size_t i = 0; i < a.discs.size(); ++i)
        if (!sameDisc(a.discs[i], b.discs[i])) return false;
    return true;
}
double elevationDeg(const FurDisc& d) { return std::asin(std::clamp(dot(d.d, d.normal), -1.0, 1.0)) * 180 / kPi; }
} // namespace

// The fur field (fur_field.h): area-weighted roots at E[n] = D A; discs independent of the thread count and of the
// cache's capacity (a cache drawn at twice the density gives the same discs: no history); stable identities and
// roots under shape edits; every drawn shape within its variation's bounds; a mask of 0 masks everything; the fine
// layer on and off.
TEST(fur_field_properties) {
    const TriMesh sphere = furSphere(20, 40, 80);
    FurLook look;  // the approved look's defaults
    look.loops.density = 0.5;
    look.fibres.density = 0.2;
    const uint64_t loopsSeed = furPartSeed(look.seed, FurLayer::Loops, 7), fibresSeed = furPartSeed(look.seed, FurLayer::Fibres, 7);
    CHECK(loopsSeed != fibresSeed);
    const FurRootCache loops = sampleFurRoots(sphere, FurLayer::Loops, look.loops.density, loopsSeed, nullptr, ExecPolicy{true});
    const FurRootCache fibres = sampleFurRoots(sphere, FurLayer::Fibres, look.fibres.density, fibresSeed, nullptr, ExecPolicy{true});
    const FurField f = buildFurField({loops, fibres, look, nullptr}, ExecPolicy{true});
    // E[n] = density * area; floor(x + u) has variance <= 1/4 per triangle.
    const double expected = look.loops.density * sphere.area(), sigma = 0.5 * std::sqrt(double(sphere.tri.size()));
    std::printf("    %lld loops (expected %.1f, 4 sigma %.1f), %lld fine fibres; %zu discs\n", static_cast<long long>(f.loops.kept), expected,
                4 * sigma, static_cast<long long>(f.fibres.kept), f.discs.size());
    CHECK(std::fabs(double(f.loops.kept) - expected) <= 4 * sigma + 1);
    CHECK(f.loops.sampled == f.loops.kept && f.loops.hidden == 0 && f.loops.masked == 0 && f.fibres.kept > 0);
    CHECK(static_cast<int64_t>(f.discs.size()) == f.loops.kept + f.fibres.kept);
    // Serial on a cache at twice the capacity: the same discs.
    const FurRootCache bigger = sampleFurRoots(sphere, FurLayer::Loops, 2 * look.loops.density, loopsSeed, nullptr, ExecPolicy{false});
    CHECK(bigger.roots.size() > loops.roots.size());
    CHECK(sameField(f, buildFurField({bigger, fibres, look, nullptr}, ExecPolicy{false})));
    // A shape edit keeps every root and identity (and the start directions: no direction setting changed).
    FurLook shaped = look;
    shaped.loops.curl = 1.2;
    shaped.loops.rootDiameter = 0.4;
    const FurField g = buildFurField({loops, fibres, shaped, nullptr}, ExecPolicy{true});
    CHECK(g.discs.size() == f.discs.size());
    std::set<std::pair<int, int>> ids;
    for (size_t i = 0; i < f.discs.size(); ++i) {
        const FurDisc &a = f.discs[i], &b = g.discs[i];
        CHECK(a.tri == b.tri && a.index == b.index && a.root == b.root && a.d == b.d && a.layer == b.layer);
        if (a.layer == FurLayer::Loops) CHECK(ids.insert({a.tri, a.index}).second);
    }
    // Every drawn value within its variation's bounds (the defaults: direction 1, curl 1, thickness 0.33, bend 0.2,
    // length 0.5; fine fibres lift 85 +- 40 deg, length in [2, 8] mm).
    const double tol = 1e-9;
    for (const FurDisc& d : f.discs) {
        CHECK(std::fabs(norm(d.d) - 1) < 1e-12 && std::fabs(norm(d.normal) - 1) < 1e-12);
        if (d.layer == FurLayer::Loops) {
            CHECK(elevationDeg(d) >= -tol && elevationDeg(d) <= 80 + tol);
            CHECK(d.radius >= 0.45 - tol && d.radius <= 1.35 + tol);
            CHECK(d.rootRadius >= 0.3 * 0.67 - tol && d.rootRadius <= 0.3 * 1.33 + tol && std::fabs(d.tipRadius - 0.8 * d.rootRadius) < tol);
            CHECK(d.turns >= 0.5 - tol && d.turns <= 1 + tol && d.rise >= 0.8 - tol && d.rise <= 1.2 + tol);
            CHECK(std::fabs(norm(d.c) - 1) < 1e-12 && std::fabs(dot(d.c, d.d)) < 1e-12 && dot(d.c, d.normal) >= 0);
        } else {
            CHECK(d.radius == 0 && d.straightLength >= 2 - tol && d.straightLength <= 8 + tol);
            CHECK(elevationDeg(d) >= 45 - tol);
            CHECK(std::fabs(d.rootRadius - 0.03) < tol && std::fabs(d.tipRadius - 0.009) < tol);
        }
    }
    // A coverage of 0: every candidate masked, counted; the fine layer off: loops only.
    const FurField none = buildFurField({loops, fibres, look, [](int, const Vec3&) { return 0.0; }}, ExecPolicy{true});
    CHECK(none.discs.empty() && none.loops.masked == f.loops.kept && none.fibres.masked == f.fibres.kept);
    FurLook loopsOnly = look;
    loopsOnly.fibres.enabled = false;
    const FurField lo = buildFurField({loops, fibres, loopsOnly, nullptr}, ExecPolicy{true});
    CHECK(static_cast<int64_t>(lo.discs.size()) == f.loops.kept && lo.fibres.kept == 0);
}

// The irregularities (fur_field.h) on a ball: all off, the fields equal the
// unmodified ones; density patches follow the noise and keep the expected mean; the same settings give the same discs
// whatever the cache held; thin islands thin both layers and never move a root; matted islands lower the lift;
// fluffy islands raise only the fine fibres' bound; clump frames stay orthonormal and lean further with strength;
// wear on a convex ball shortens the loops.
TEST(fur_irregularities) {
    const TriMesh sphere = furSphere(20, 40, 80);
    FurLook look;
    look.loops.density = 0.5;
    look.fibres.density = 0.3;
    const uint64_t nameHash = 5;
    auto field = [&](const FurLook& l, double capacityScale) {
        const FurRootCache loops = sampleFurRoots(sphere, FurLayer::Loops, furDensityBound(l, FurLayer::Loops) * capacityScale,
                                                  furPartSeed(l.seed, FurLayer::Loops, nameHash), nullptr, ExecPolicy{true});
        const FurRootCache fibres = sampleFurRoots(sphere, FurLayer::Fibres, furDensityBound(l, FurLayer::Fibres) * capacityScale,
                                                   furPartSeed(l.seed, FurLayer::Fibres, nameHash), nullptr, ExecPolicy{true});
        const FurPartModifiers mods = buildFurModifiers(sphere, l, nameHash);
        return buildFurField({loops, fibres, l, nullptr, &mods}, ExecPolicy{true});
    };
    const FurField plain = field(look, 1);
    CHECK(!buildFurModifiers(sphere, look, nameHash).active);
    // Patches: +- 80 % of the density, 8 mm features.
    FurLook patchy = look;
    patchy.variation.density = 0.8;
    patchy.variation.densityScale = 8;
    const FurField p = field(patchy, 1);
    CHECK(sameField(p, field(patchy, 1.6)));  // a larger cache: the same discs
    const FurPartModifiers pm = buildFurModifiers(sphere, patchy, nameHash);
    double expected = 0, high = 0, low = 0, highArea = 0, lowArea = 0;
    for (size_t t = 0; t < sphere.tri.size(); ++t) {
        const std::array<int, 3>& v = sphere.tri[t];
        const double area = 0.5 * norm(cross(sphere.pos[v[1]] - sphere.pos[v[0]], sphere.pos[v[2]] - sphere.pos[v[0]]));
        const double f = (pm.patch[v[0]] + pm.patch[v[1]] + pm.patch[v[2]]) / 3;
        expected += patchy.loops.density * (1 + 0.8 * f) * area;
        if (f > 0.3) highArea += area;
        if (f < -0.3) lowArea += area;
    }
    for (const FurDisc& d : p.discs) {
        if (d.layer != FurLayer::Loops) continue;
        const std::array<int, 3>& v = sphere.tri[static_cast<size_t>(d.tri)];
        const double f = (pm.patch[v[0]] + pm.patch[v[1]] + pm.patch[v[2]]) / 3;
        high += f > 0.3;
        low += f < -0.3;
    }
    std::printf("    patches: %lld loops kept (expected %.0f; plain %lld), %lld thinned; density where f > 0.3: %.3f, f < -0.3: %.3f per mm2\n",
                static_cast<long long>(p.loops.kept), expected, static_cast<long long>(plain.loops.kept), static_cast<long long>(p.loops.thinned),
                high / std::max(1.0, highArea), low / std::max(1.0, lowArea));
    CHECK(std::fabs(double(p.loops.kept) - expected) < 0.05 * expected);
    CHECK(high / std::max(1.0, highArea) > 1.5 * (low / std::max(1.0, lowArea)));
    // Thin islands everywhere: both layers sparser, every kept root one of the unmodified ones, at its place.
    FurLook thin = look;
    thin.islands = {0.5, 20, 0.8, 0, 0, 1};
    const FurField th = field(thin, 1);
    std::map<std::pair<int, int>, Vec3> rootOf;
    for (const FurDisc& d : field(look, 1).discs) rootOf[{d.tri * 2 + int(d.layer), d.index}] = d.root;
    size_t moved = 0;
    for (const FurDisc& d : th.discs) {
        const auto it = rootOf.find({d.tri * 2 + int(d.layer), d.index});
        moved += it == rootOf.end() || !(it->second == d.root);
    }
    std::printf("    thin islands: loops %lld -> %lld, fine fibres %lld -> %lld, %zu moved\n", static_cast<long long>(plain.loops.kept),
                static_cast<long long>(th.loops.kept), static_cast<long long>(plain.fibres.kept), static_cast<long long>(th.fibres.kept), moved);
    CHECK(th.loops.kept < plain.loops.kept && th.fibres.kept < plain.fibres.kept && moved == 0);
    // Matted islands everywhere: the loops lie lower; fluffy: only the fine fibres' bound rises.
    FurLook matted = look;
    matted.islands = {0.5, 20, 0.8, 1, 0, 0};
    double liftPlain = 0, liftMatted = 0;
    for (const FurDisc& d : plain.discs) liftPlain += d.layer == FurLayer::Loops ? dot(d.d, d.normal) : 0;
    const FurField mt = field(matted, 1);
    for (const FurDisc& d : mt.discs) liftMatted += d.layer == FurLayer::Loops ? dot(d.d, d.normal) : 0;
    CHECK(liftMatted / double(mt.loops.kept) < liftPlain / double(plain.loops.kept));
    FurLook fluffy = look;
    fluffy.islands = {0.5, 20, 0.8, 0, 1, 0};
    CHECK(furDensityBound(fluffy, FurLayer::Loops) == look.loops.density && furDensityBound(fluffy, FurLayer::Fibres) > look.fibres.density);
    // Clumps: frames orthonormal; a stronger clump leans further toward its direction.
    double lean[2] = {0, 0};
    for (int k = 0; k < 2; ++k) {
        FurLook clumped = look;
        clumped.irregular.clumpSize = 3;
        clumped.irregular.clumpStrength = k == 0 ? 0.4 : 1.0;
        const FurField cf = field(clumped, 1);
        const FurPartModifiers cm = buildFurModifiers(sphere, clumped, nameHash);
        CHECK(!cm.clumps.empty());
        for (const FurDisc& d : cf.discs) {
            if (d.layer != FurLayer::Loops) continue;
            CHECK(std::fabs(norm(d.d) - 1) < 1e-9 && std::fabs(norm(d.c) - 1) < 1e-9 && std::fabs(dot(d.d, d.c)) < 1e-9);
        }
        // The mean alignment of neighbouring loops' directions (pairs closer than 1 mm): higher with a stronger clump.
        double align = 0;
        int pairs = 0;
        for (size_t i = 0; i < cf.discs.size(); ++i)
            for (size_t j = i + 1; j < cf.discs.size(); ++j)
                if (cf.discs[i].layer == FurLayer::Loops && cf.discs[j].layer == FurLayer::Loops && norm(cf.discs[i].root - cf.discs[j].root) < 1.0) {
                    align += dot(cf.discs[i].d, cf.discs[j].d);
                    ++pairs;
                }
        lean[k] = align / std::max(1, pairs);
    }
    std::printf("    clumps: neighbouring loops' direction alignment %.3f at strength 0.4, %.3f at 1\n", lean[0], lean[1]);
    CHECK(lean[1] > lean[0]);
    // Wear on a convex ball: shorter loops than without.
    FurLook worn = look;
    worn.irregular.wear = 1;
    double turnsPlain = 0, turnsWorn = 0;
    for (const FurDisc& d : plain.discs) turnsPlain += d.layer == FurLayer::Loops ? d.turns : 0;
    const FurField wf = field(worn, 1);
    for (const FurDisc& d : wf.discs) turnsWorn += d.layer == FurLayer::Loops ? d.turns : 0;
    std::printf("    wear: mean turns %.3f -> %.3f\n", turnsPlain / double(plain.loops.kept), turnsWorn / double(wf.loops.kept));
    CHECK(turnsWorn / double(wf.loops.kept) < turnsPlain / double(plain.loops.kept));
}

// The curves (fur_loops.h): a reference loop (variation off, unclamped) has equal chords of L / n, each at least
// kDiscMinSegmentRadii root radii; one nominal loop turns more than 320 deg and ends within the ramp's offset, its
// climb and one diameter of its root. The settle keeps every segment's length, leaves no vertex in the body and does not depend
// on the thread count; gravity droops a straight fibre; the strands keep the Fur export's joint rule and break the
// wool rule (the negative control); a crossing segment is caught, a fibre leaving the surface is not.
TEST(fur_loops_properties) {
    FurDisc ref;
    ref.root = {0, 0, 0};
    ref.normal = {0, 0, 1};
    ref.d = {1, 0, 0};
    ref.c = {0, 0, 1};
    ref.radius = 0.9;
    ref.turns = 1;
    ref.rise = 1;
    ref.rootRadius = 0.1;
    ref.tipRadius = 0.08;
    CHECK(!furLengthClamped(ref));
    const Hair h = furRestShape(ref);
    const int n = furSegmentCount(ref);
    const double L = furFibreLength(ref);
    CHECK(static_cast<int>(h.points.size()) == n + 2);
    double turning = std::acos(std::clamp(dot(ref.d, normalized(h.points[2] - h.points[1])), -1.0, 1.0));
    for (int k = 2; k <= n + 1; ++k) {
        const double seg = norm(h.points[k] - h.points[k - 1]);
        CHECK_NEAR(seg, L / n, 1e-9);
        CHECK(seg >= kDiscMinSegmentRadii * ref.rootRadius);
        if (k > 2)
            turning += std::acos(std::clamp(dot(normalized(h.points[k - 1] - h.points[k - 2]), normalized(h.points[k] - h.points[k - 1])), -1.0, 1.0));
    }
    // The ends cannot meet: the curvature ramp over a = (pi r / 2) rise shifts the loop's circle like a clothoid
    // transition curve, by a / 2 along d and p = a^2 / (24 r) away from it (the transition's shift), so the point where
    // the tangent is d again lies that far from the root; the climb adds kDiscClimbDiameters root diameters across
    // the plane, the equal chords up to one diameter.
    const double a = 0.5 * kPi * ref.radius * ref.rise, p = a * a / (24 * ref.radius), diameter = 2 * ref.rootRadius;
    const double gap = norm(h.points.back() - ref.root), bound = 0.5 * a + p + kDiscClimbDiameters * diameter + diameter;
    std::printf("    reference loop: %d segments of %.3f mm (L %.3f mm), turning %.1f deg, tip %.3f mm from the root (ramp alone %.3f mm, bound "
                "%.3f mm)\n",
                n, L / n, L, turning * 180 / kPi, gap, std::hypot(0.5 * a, p), bound);
    CHECK(turning > 320 * kPi / 180 && gap <= bound);
    // The settle on a ball, at the approved look (both layers).
    const TriMesh sphere = furSphere(20, 40, 80);
    const MeshQuery body(sphere);
    FurLook look;
    look.loops.density = 0.3;
    look.fibres.density = 0.1;
    const FurRootCache loops = sampleFurRoots(sphere, FurLayer::Loops, look.loops.density, furPartSeed(look.seed, FurLayer::Loops, 3), nullptr, ExecPolicy{true});
    const FurRootCache fibres = sampleFurRoots(sphere, FurLayer::Fibres, look.fibres.density, furPartSeed(look.seed, FurLayer::Fibres, 3), nullptr, ExecPolicy{true});
    const FurField field = buildFurField({loops, fibres, look, nullptr}, ExecPolicy{true});
    FurReport rep;
    const std::vector<Hair> hp = growLoops(field, body, look.settle, ExecPolicy{true}, &rep), hs = growLoops(field, body, look.settle, ExecPolicy{false}, nullptr);
    std::printf("    ball: %lld + %lld grown, %lld + %lld stuck, %lld loops clamped\n", static_cast<long long>(rep.loops.grown),
                static_cast<long long>(rep.fibres.grown), static_cast<long long>(rep.loops.stuck), static_cast<long long>(rep.fibres.stuck),
                static_cast<long long>(rep.loops.clamped));
    CHECK(hp.size() == hs.size() && !hp.empty() && rep.loops.grown > 0 && rep.fibres.grown > 0);
    std::map<std::pair<int, int>, const FurDisc*> discOf;
    for (const FurDisc& d : field.discs) discOf[{d.tri * 2 + int(d.layer), d.index}] = &d;
    for (size_t i = 0; i < hp.size(); ++i) {
        CHECK(hp[i].points == hs[i].points);
        if (i % 5) continue;
        const FurDisc* d = discOf[{hp[i].tri * 2 + (hp[i].rootRadius < 0.1 ? 1 : 0), hp[i].index}];
        CHECK(d != nullptr);
        if (!d) continue;
        const Hair rest = furRestShape(*d);
        CHECK(rest.points.size() == hp[i].points.size());
        for (size_t k = 1; k < rest.points.size(); ++k)
            CHECK_NEAR(norm(hp[i].points[k] - hp[i].points[k - 1]), norm(rest.points[k] - rest.points[k - 1]), 1e-9);
        for (size_t k = 2; k < hp[i].points.size(); ++k) CHECK(!body.inside(hp[i].points[k]));
    }
    // The Fur export's joint rule holds; the wool rule (4 r) fails on the thick loops.
    FibreSet set;
    for (const Hair& x : hp) addHairFibre(set, x, {});
    const StrandSetF strands = fibresToStrands(set, SceneMapping{1000, {}}, -1);
    StrandRules furRules;
    furRules.minSegOverR = kDiscMinSegOverR;
    const StrandCheck ok = checkStrands(strands, furRules), wool = checkStrands(strands);
    if (!ok.ok) std::printf("    strand check: %s\n", ok.error.c_str());
    std::printf("    shortest segment %.2f r; the wool rule: %s\n", ok.minSegOverR, wool.ok ? "passes" : wool.error.c_str());
    CHECK(ok.ok && !wool.ok && wool.fault == StrandFault::SegmentShort);
    // Gravity droops a straight fibre on a vertical wall (normal +x); none keeps it level.
    const TriMesh wall = furWall(40, 40, 20);
    const MeshQuery wq(wall);
    FurLook straight;
    straight.loops.curl = 0;
    straight.loops.straightLength = 8;
    straight.loops.rootDiameter = 0.08;
    straight.loops.density = 0.2;
    straight.variation.direction = 0;  // lift 90: along the normal
    straight.variation.length = 0;
    straight.fibres.enabled = false;
    straight.settle.gravityFullLength = 0;
    const FurRootCache wr = sampleFurRoots(wall, FurLayer::Loops, straight.loops.density, 11, nullptr, ExecPolicy{true});
    const FurField wf = buildFurField({wr, FurRootCache{}, straight, nullptr}, ExecPolicy{true});
    double drop[2] = {0, 0};
    for (int gi = 0; gi < 2; ++gi) {
        FurSettle settle = straight.settle;
        settle.strand.gravity = gi == 0 ? 0.0 : 1.0;
        const std::vector<Hair> hw = growLoops(wf, wq, settle, ExecPolicy{true}, nullptr);
        for (const Hair& x : hw) drop[gi] += (x.points[1].z - x.points.back().z) / double(hw.size());
    }
    std::printf("    straight 8 mm fibres on a wall: mean tip drop %.4f mm without gravity, %.3f mm with it\n", drop[0], drop[1]);
    CHECK(std::fabs(drop[0]) < 1e-6 && drop[1] > 0.1);
    // The stuck test's stated geometry: a crossing is caught; a segment leaving the surface, or lying in it, is not.
    CHECK(wq.segmentCrosses({1, 0, 0}, {-1, 0, 0.5}));
    CHECK(!wq.segmentCrosses({0, 0, 0}, {1, 0, 0}, kRootCrossingSkip));
    CHECK(!wq.segmentCrosses({0, -1, 0}, {0, 1, 0}));
    Hair out;
    out.points = {{-0.1, 0, 0}, {0, 0, 0}, {1, 0, 0}, {2, 0, 0.5}};
    out.rootRadius = out.tipRadius = 0.01;
    CHECK(!hairInsideBody(out, wq));
    out.points.push_back({-1, 0, 1});
    CHECK(hairInsideBody(out, wq));
}

// Ready without any knit, and exports its mesh plus the pile; knitting edits keep it current, pile edits block it;
// OBJ holds the mesh and the fibres as lines.
TEST(workflow_plane_fur_export) {
    const std::string d = freshDir("fur"), root = freshDir("fur_cache"), out = freshDir("fur_out");
    ImportCacheEnv env(root);
    // No uv in the file: the file's mesh as it is (fur needs no charts).
    writeText(d + "/plane.obj", objText(grid(12, 8, 60, 40), "", "", false));
    ImportedMesh m;
    std::string err;
    CHECK(readMeshFile(d + "/plane.obj", ImportLimits{}, m, &err));
    ImportOptions o;
    o.unitMm = 1;
    o.up = UpAxis::Z;
    PreparedMesh p;
    CHECK(prepareImport(m, o, p, &err));
    CHECK(p.mesh.uv.empty() && !p.chartsGenerated);
    const uint64_t id = importIdentity(m, o);
    CHECK(storeImport(root, {std::make_shared<const ImportedMesh>(m), o, {}, id}, {p, {}}, &err));
    SceneSettings s;
    s.scene = "import";
    s.importSource = d + "/plane.obj";
    s.importIdentity = identityHex(id);
    s.importUp = "z";
    const SceneEntry* e = sceneEntryFor(s);
    CHECK(e != nullptr);
    s.workflow = "fur";
    s.pile.density = 30;
    s.pile.fibresDensity = 30;
    // A knit refuses this mesh (the Fur workflow has no knit); the pile grows on it.
    std::string kerr;
    CHECK(!run(BuiltScene::Kind::Knit, s, &kerr) && kerr.find("no knit") != std::string::npos);
    auto fur = run(BuiltScene::Kind::Fibres, s, &err);
    if (!fur) std::printf("    %s\n", err.c_str());
    CHECK(fur && fur->fibreResult && fur->fibreResult->workflow == Workflow::Fur && fur->fibreResult->knitId == 0);
    CHECK(fur->fibreResult->fibres > 0 && fur->fibreResult->base);
    std::printf("%s", fur->fibreResult->log.c_str());
    CHECK(exportEligibility(s, nullptr, fur->fibreResult.get(), FibreChoice::Included).state == ExportEligibility::State::Ready);
    SceneSettings k = s;
    k.diameter = 0.9;
    k.plies = 3;
    k.pileMaterial.color = {0.2, 0.3, 0.4};
    CHECK(exportEligibility(k, nullptr, fur->fibreResult.get(), FibreChoice::Included).state == ExportEligibility::State::Ready);
    k = s;
    k.pile.length = 0.8;
    CHECK(exportEligibility(k, nullptr, fur->fibreResult.get(), FibreChoice::Included).state == ExportEligibility::State::Blocked);
    CHECK(exportEligibility(s, nullptr, nullptr, FibreChoice::Included).state == ExportEligibility::State::Blocked);
    // The package: the mesh and the fibres, which read back.
    const std::string target = out + "/plane_fur";
    auto ex = run(BuiltScene::Kind::Export, s, &err, nullptr, nullptr, fur->fibreResult, target, "ppt", true);
    if (!ex) std::printf("    %s\n", err.c_str());
    CHECK(ex && fs::exists(target + "/scene.json"));
    size_t strandFiles = 0, meshFiles = 0, strands = 0, strandVertices = 0;
    for (const std::string& f : filesUnder(target)) {
        if (f.find(".pptstrands") != std::string::npos) {
            ++strandFiles;
            StrandSetF back;
            CHECK(readStrandsFile(target + "/" + f, back, &err));
            CHECK(back.strands() > 0);
            strands += back.strands();
            strandVertices += back.vertices();
        }
        meshFiles += f.rfind("meshes/", 0) == 0;
    }
    CHECK(strandFiles >= 1 && meshFiles >= 1);
    // OBJ: the meshes and the fibres as line elements, one per strand of the package, through vertices of the file.
    auto obj = run(BuiltScene::Kind::Export, s, &err, nullptr, nullptr, fur->fibreResult, out + "/plane_fur_obj", "obj", true);
    CHECK(obj);
    std::string objFile;
    for (const std::string& f : filesUnder(out + "/plane_fur_obj"))
        if (f.size() > 4 && f.compare(f.size() - 4, 4, ".obj") == 0) objFile = out + "/plane_fur_obj/" + f;
    std::ifstream in(objFile);
    size_t v = 0, l = 0, lineVertices = 0, maxIndex = 0;
    std::string line;
    while (std::getline(in, line)) {
        if (line.rfind("v ", 0) == 0) ++v;
        if (line.rfind("l ", 0) != 0) continue;
        ++l;
        std::istringstream is(line.substr(2));
        for (size_t i; is >> i;) {
            ++lineVertices;
            maxIndex = std::max(maxIndex, i);
        }
    }
    std::printf("    OBJ %s: %zu v, %zu lines through %zu vertices; package: %zu strands, %zu vertices\n", objFile.c_str(), v, l, lineVertices,
                strands, strandVertices);
    CHECK(l == strands && lineVertices == strandVertices && maxIndex <= v && maxIndex > 0);
}

// The reference preview (previews.h): built without a knit or a mesh, keyed by the settings
// it depends on; thickness / taper change the drawn tube radii, density the count; the pile on a sphere.
TEST(fibre_preview_reference) {
    SceneSettings s;
    std::string err;
    FibrePreview a;
    CHECK(buildFibrePreview(s, a, &err));
    std::printf("    wool: %s (%.0f ms)\n", a.note.c_str(), a.ms);
    CHECK(a.fibres > 0 && a.meshes.size() == 3 && !a.meshes[0].idx.empty() && !a.meshes[1].idx.empty());
    CHECK(a.ms < 5000);
    // A square patch: its sides differ by less than half a course (5 wales, round(5 w / h - 1) courses).
    std::printf("    patch %.1f x %.1f mm\n", a.hi.x - a.lo.x, a.hi.y - a.lo.y);
    CHECK(std::abs((a.hi.x - a.lo.x) - (a.hi.y - a.lo.y)) < 0.5 * s.course);
    // The two populations in their display colours (yellow along the yarn, blue stray).
    CHECK(a.meshes[1].color[0] == float(kOrientedFibreColour.x) && a.meshes[2].color[2] == float(kStrayFibreColour.z));
    // The key: the fibre and yarn settings, not the scene.
    SceneSettings t = s;
    t.scene = s.scene == "teddy" ? "wool-swatch" : "teddy";
    CHECK(fibrePreviewKey(t) == a.key);
    t.fibreDensity *= 0.5;
    CHECK(fibrePreviewKey(t) != a.key);
    FibrePreview b;
    CHECK(buildFibrePreview(t, b, &err));
    CHECK(b.fibres < a.fibres);
    // Thickness: the tube radii follow (the mesh's extent across a fibre grows with the thickness).
    auto meanRadius = [](const DisplayMesh& d) {
        // Rings of 5 vertices: the mean distance of a ring's vertices from their centroid.
        double sum = 0;
        size_t n = 0;
        for (size_t v = 0; v + 15 <= d.pos.size(); v += 15) {
            Vec3 c{0, 0, 0};
            for (int j = 0; j < 5; ++j) c += Vec3{d.pos[v + 3 * j], d.pos[v + 3 * j + 1], d.pos[v + 3 * j + 2]} / 5.0;
            for (int j = 0; j < 5; ++j) sum += norm(Vec3{d.pos[v + 3 * j], d.pos[v + 3 * j + 1], d.pos[v + 3 * j + 2]} - c);
            n += 5;
        }
        return n ? sum / double(n) : 0.0;
    };
    t = s;
    t.directed.setThicknessRange(2 * s.directed.thicknessMin(), 2 * s.directed.thicknessMax());
    FibrePreview c;
    CHECK(buildFibrePreview(t, c, &err));
    std::printf("    directed fibre tube radius %.4f -> %.4f mm with twice the thickness\n", meanRadius(a.meshes[1]), meanRadius(c.meshes[1]));
    CHECK(meanRadius(c.meshes[1]) > 1.5 * meanRadius(a.meshes[1]));
    // Fur: the pile on the reference sphere; thinner loops draw thinner tubes.
    SceneSettings f;
    f.workflow = "fur";
    FibrePreview fu;
    CHECK(buildFibrePreview(f, fu, &err));
    std::printf("    fur: %s (%.0f ms)\n", fu.note.c_str(), fu.ms);
    CHECK(fu.fibres > 500 && fu.workflow == Workflow::Fur);
    f.pile.thickness = 0.3;
    FibrePreview fw;
    CHECK(buildFibrePreview(f, fw, &err));
    std::printf("    loop tube radius %.4f -> %.4f mm at half the thickness\n", meanRadius(fu.meshes[1]), meanRadius(fw.meshes[1]));
    CHECK(meanRadius(fw.meshes[1]) < 0.75 * meanRadius(fu.meshes[1]));
}

// Fur coverage and budget (fur_stage.h): a mask is used only when chosen; Reload re-hashes it (the fur key changes)
// and a mask file changed since loading is refused; an excluded part grows nothing; a smaller segment budget keeps
// a subset of the same fibres, unmoved.
TEST(fur_coverage_rules) {
    const std::string d = freshDir("furcov"), root = freshDir("furcov_cache");
    ImportCacheEnv env(root);
    SceneSettings s = importPlane(d, root);  // uv present: the mask has colour coordinates
    s.workflow = "fur";
    s.pile.density = 30;
    s.pile.fibresDensity = 30;
    const SceneEntry* e = sceneEntryFor(s);
    CHECK(e != nullptr);
    std::string err;
    Scene t;
    CHECK(loadBaseParts(*e, s, SourceOptions{}, t, &err));
    CHECK(!t.parts.empty());
    // A black mask image that exists but is not chosen: no root is masked.
    const std::string mask = d + "/mask.png";
    {
        std::vector<uint8_t> rgb(16 * 16 * 3, 0);
        CHECK(writePngRgb8(mask, 16, 16, rgb));
    }
    s.pile.maskPath = mask;
    s.pile.maskHash = fileHashHex(mask);
    SceneFibres all;
    CHECK(grown(all, growSceneFur(*e, t, s, &err)));
    CHECK(all.log.find(" 0 outside the coverage") != std::string::npos);
    size_t fibres = 0;
    for (const PartFibres& p : all.parts) fibres += p.fs.fibres.size();
    CHECK(fibres > 0);
    // Chosen: the black mask removes every root.
    s.pile.coverage = "mask";
    SceneFibres none;
    CHECK(grown(none, growSceneFur(*e, t, s, &err)));
    size_t left = 0;
    for (const PartFibres& p : none.parts) left += p.fs.fibres.size();
    CHECK(left == 0);
    // The file changes: refused until reloaded; reloading changes the fur key.
    const uint64_t k0 = furKey(s, 0);
    {
        std::vector<uint8_t> rgb(16 * 16 * 3, 255);
        CHECK(writePngRgb8(mask, 16, 16, rgb));
    }
    SceneFibres stale;
    CHECK(!grown(stale, growSceneFur(*e, t, s, &err)) && err.find("Reload") != std::string::npos);
    s.pile.maskHash = fileHashHex(mask);
    CHECK(furKey(s, 0) != k0);
    SceneFibres white;
    CHECK(grown(white, growSceneFur(*e, t, s, &err)));
    size_t full = 0;
    for (const PartFibres& p : white.parts) full += p.fs.fibres.size();
    CHECK(full == fibres);  // white = no restriction: the same roots as the whole surface
    // An excluded part grows nothing.
    s.pile.coverage = "entire";
    s.pile.excludedParts = t.parts[0].name;
    SceneFibres ex;
    CHECK(grown(ex, growSceneFur(*e, t, s, &err)));
    for (const PartFibres& p : ex.parts) CHECK(p.part != 0);
    CHECK(!setSetting(s, "pile.coverage", "preset", &err));
    // A third of the segments: fewer fibres, each one of the full run's, point for point.
    s.pile.excludedParts.clear();
    long long segments = 0;
    std::map<std::pair<int, int>, std::vector<Vec3f>> fullFibres;
    for (const PartFibres& p : all.parts)
        for (const Fibre& f : p.fs.fibres) {
            const std::span<const Vec3f> pts = p.fs.pointsOf(f);
            segments += static_cast<long long>(pts.size());
            fullFibres[{f.yarnId * 4 + int(f.radius0 < 0.1f), f.candidate}] = std::vector<Vec3f>(pts.begin(), pts.end());
        }
    s.segmentBudget = segments / 3;
    SceneFibres cut;
    CHECK(grown(cut, growSceneFur(*e, t, s, &err)));
    size_t kept = 0, moved = 0;
    for (const PartFibres& p : cut.parts)
        for (const Fibre& f : p.fs.fibres) {
            ++kept;
            const std::span<const Vec3f> pts = p.fs.pointsOf(f);
            const auto it = fullFibres.find({f.yarnId * 4 + int(f.radius0 < 0.1f), f.candidate});
            moved += it == fullFibres.end() || !std::equal(pts.begin(), pts.end(), it->second.begin(), it->second.end(),
                                                           [](const Vec3f& a, const Vec3f& b) { return a.x == b.x && a.y == b.y && a.z == b.z; });
        }
    std::printf("    budget %lld of %lld segments: %zu of %zu fibres, %zu moved\n", s.segmentBudget, segments, kept, fibres, moved);
    CHECK(kept > 0 && kept < fibres && moved == 0);
}

// A fresh knit of a mesh with overlapped seams, then Fibres on that knit in the same session. The
// knit's surface and atlas refer to their part's mesh, and the gap-patch pass appends to the parts
// after the first parts are knitted: moving a part must leave those references valid. The parts are
// moved here on purpose, so a part that owned its mesh inline fails every time, not only when the
// freed memory happens to be overwritten.
TEST(fresh_mesh_knit_then_fibres_same_session) {
    const std::string d = freshDir("session"), root = freshDir("session_cache");
    ImportCacheEnv env(root);
    SceneSettings s = importPlane(d, root);
    applyFullProfile(s, *sceneEntryFor(s));
    CHECK(overlapCoverage(s) && s.gapPatches);
    std::string err;
    auto fin = run(BuiltScene::Kind::Knit, s, &err);
    CHECK(fin && fin->knit && fin->knit->scene);
    if (!fin || !fin->knit || !fin->knit->scene) return;
    Scene& scene = *fin->knit->scene;
    std::vector<ScenePart> moved;
    moved.reserve(scene.parts.size());
    for (ScenePart& r : scene.parts) moved.push_back(std::move(r));
    scene.parts = std::move(moved);
    int knitted = 0, gapPatches = 0;
    for (const ScenePart& r : scene.parts) {
        gapPatches += r.gapPatch ? 1 : 0;
        if (!r.knit) continue;
        ++knitted;
        CHECK(r.mesh && &r.knit->surface->mesh() == r.mesh.get());
        if (r.knit->atlas) CHECK(&r.knit->atlas->mesh() == r.mesh.get());
    }
    std::printf("    %zu parts: %d knitted, %d gap patches\n", scene.parts.size(), knitted, gapPatches);
    CHECK(knitted > 0 && gapPatches > 0);
    s.fibres = true;
    s.fibreDensity = 20;
    auto fib = run(BuiltScene::Kind::Fibres, s, &err, nullptr, fin->knit);
    CHECK(fib && fib->fibreResult && !fib->fibreResult->scene->parts.empty());
}

// A yarn of fibres rendered as maps in the viewport (StrandMeshRequest::bakeMaps): without the entry the rebuild
// keeps the plain tubes and says so; asked to bake, it bakes the maps into the cache (a yarn at a time, reported)
// and shows one mesh per yarn under its cavity atlas; the next rebuild takes them from the cache. The knit result
// stays the same. A fresh seed each run keeps the first rebuild honest (the cache keeps its newest entries).
TEST(strand_mesh_yarn_maps) {
    const std::string d = freshDir("yarn_maps"), root = freshDir("yarn_maps_cache"), scenes = d + "/cache/scenes";
    ImportCacheEnv env(root);
    SceneSettings s = importPlane(d, root);
    applyFullProfile(s, *sceneEntryFor(s));
    s.fibres = false;
    s.filling = false;
    s.yarnBody.enabled = true;
    s.yarnBody.plies = 2;
    s.yarnBody.fibres = 30;
    s.yarnBody.render = "maps";
    s.yarnBody.mapTexelUm = 100;
    s.yarnBody.seed = static_cast<int>(std::chrono::steady_clock::now().time_since_epoch().count() % 1000000) + 1;
    std::string err;
    auto fin = run(BuiltScene::Kind::Knit, s, &err, nullptr, nullptr, nullptr, "", "ppt", false, scenes);  // stored: plain tubes
    CHECK(fin && fin->knit && !fin->cached);
    auto yarnMeshes = [](const BuiltScene& b) {
        int plain = 0, mapped = 0;
        for (const DisplayMesh& m : b.meshes)
            if (m.kind == DisplayMesh::Kind::Yarn) (m.albedo >= 0 ? mapped : plain)++;
        return std::make_pair(plain, mapped);
    };
    StrandMeshRequest j;
    j.settings = std::make_shared<const SceneSettings>(s);
    j.shown = fin;
    StageLog log;
    auto plain = rebuildStrandMesh(j, JobControl{0, &log}, &err);
    if (!plain) std::printf("    %s\n", err.c_str());
    CHECK(plain && plain->knit == fin->knit);
    auto [plainCount, mappedCount] = yarnMeshes(*plain);
    std::printf("    without the maps: %d plain, %d mapped yarn meshes; note: %s\n", plainCount, mappedCount, plain->displayNote.c_str());
    CHECK(plainCount > 0 && mappedCount == 0 && plain->displayNote.find("Bake the yarn maps") != std::string::npos);
    j.bakeMaps = true;
    StageLog bakeLog;
    auto baked = rebuildStrandMesh(j, JobControl{0, &bakeLog}, &err);
    if (!baked) std::printf("    %s\n", err.c_str());
    CHECK(baked && baked->knit == fin->knit && baked->displayNote.empty());
    std::tie(plainCount, mappedCount) = yarnMeshes(*baked);
    bool reported = false;
    for (const std::string& st : bakeLog.stages) reported = reported || st.rfind("baking ", 0) == 0;
    std::printf("    baked: %d plain, %d mapped yarn meshes, %zu textures, bake stages %s\n", plainCount, mappedCount, baked->albedo.size(),
                reported ? "reported" : "missing");
    CHECK(plainCount == 0 && mappedCount > 0 && baked->albedo.size() >= size_t(mappedCount) && reported);
    for (const DisplayMesh& m : baked->meshes)
        if (m.kind == DisplayMesh::Kind::Yarn) CHECK(m.name.find("#y") != std::string::npos && m.albedoUv.size() * 3 == m.pos.size() * 2);
    j.bakeMaps = false;
    auto again = rebuildStrandMesh(j, JobControl{}, &err);
    CHECK(again && yarnMeshes(*again).second == mappedCount && again->displayNote.empty());
    // The Knit job shows the cached maps as well: knitted again, and restored from the scene cache (whose entry
    // holds the plain tubes, never the maps' pixels).
    auto knitAgain = run(BuiltScene::Kind::Knit, s, &err);
    if (!knitAgain) std::printf("    %s\n", err.c_str());
    CHECK(knitAgain && !knitAgain->cached && yarnMeshes(*knitAgain).second == mappedCount);
    auto restored = run(BuiltScene::Kind::Knit, s, &err, nullptr, nullptr, nullptr, "", "ppt", false, scenes);
    if (!restored) std::printf("    %s\n", err.c_str());
    if (restored) std::printf("    restored from the scene cache: %s, %d mapped yarn meshes\n", restored->cached ? "yes" : "no", yarnMeshes(*restored).second);
    CHECK(restored && restored->cached && yarnMeshes(*restored).second == mappedCount && restored->displayNote.empty());
}

// The strand tessellation (strand.tube_*) is never a knitting input (viewer plan P5). Medium ->
// Fine: the knitted curves and their validation are reused (the rebuild job runs no knitting stage
// and returns the same knit result), the shown tubes are rebuilt at the new tessellation, the
// existing fibres become stale (their roots are placed in the exported tube) and regenerate, and
// the export writes the new tube.
TEST(strand_mesh_rebuild_contract) {
    const std::string d = freshDir("tess"), root = freshDir("tess_cache"), out = freshDir("tess_out");
    ImportCacheEnv env(root);
    SceneSettings s = importPlane(d, root);
    applyFullProfile(s, *sceneEntryFor(s));
    std::string err;
    auto fin = run(BuiltScene::Kind::Knit, s, &err);
    CHECK(fin && fin->knit);
    s.fibres = true;
    s.fibreDensity = 20;
    auto fib = run(BuiltScene::Kind::Fibres, s, &err, nullptr, fin->knit);
    CHECK(fib && fib->fibreResult);
    auto yarnTriangles = [](const BuiltScene& b) {
        long long n = 0;
        for (const DisplayMesh& m : b.meshes) n += m.kind == DisplayMesh::Kind::Yarn ? static_cast<long long>(m.idx.size() / 3) : 0;
        return n;
    };
    SceneSettings f = s;  // Fine
    f.tubeSides = 16;
    f.tubeRing = 0.2;
    f.plyTubeSides = 10;
    f.plyTubeRing = 0.15;
    CHECK(geometryKey(f) == geometryKey(s) && fibreResultKey(f) != fibreResultKey(s));
    StageLog log;
    StrandMeshRequest j;
    j.settings = std::make_shared<const SceneSettings>(f);
    j.shown = fin;
    auto fine = rebuildStrandMesh(j, JobControl{0, &log}, &err);
    CHECK(fine && fine->knit == fin->knit && !fine->partial);
    // Its stages: the restored knit's verification (usableKnit, as every downstream job) and the
    // rebuild; nothing knits.
    for (const std::string& st : log.stages) CHECK(st == "rebuilding the strand mesh" || st == "verifying the restored knit (gate a)");
    const long long medium = yarnTriangles(*fin), fineTris = yarnTriangles(*fine);
    std::printf("    yarn triangles: medium %lld, fine %lld (x%.2f); %zu meshes before and after\n", medium, fineTris,
                medium ? double(fineTris) / medium : 0.0, fin->meshes.size());
    CHECK(medium > 0 && fineTris > 3 * medium && fine->meshes.size() == fin->meshes.size() && fine->parts.size() == fin->parts.size());
    // The fibres of the Medium tube are stale for Fine, and regenerate on the same knit.
    const ExportEligibility stale = exportEligibility(f, fin->knit.get(), fib->fibreResult.get(), FibreChoice::Included);
    CHECK(stale.state == ExportEligibility::State::Blocked && stale.fix == "Update fibres");
    auto fib2 = run(BuiltScene::Kind::Fibres, f, &err, nullptr, fin->knit);
    CHECK(fib2 && fib2->fibreResult && fib2->fibreResult->key == fibreResultKey(f) && fib2->fibreResult->knitId == fin->knit->id);
    CHECK(exportEligibility(f, fin->knit.get(), fib2->fibreResult.get(), FibreChoice::Included).state != ExportEligibility::State::Blocked);
    // The export's yarn tube: the package's meshes at both tessellations (vertices per PLY header).
    auto tubeVertices = [&](const SceneSettings& t, const std::string& name) {
        auto ex = run(BuiltScene::Kind::Export, t, &err, nullptr, fin->knit, nullptr, out + "/" + name, "ppt", false);
        CHECK(ex);
        long long v = 0;
        for (const auto& e : fs::directory_iterator(out + "/" + name + "/meshes")) {
            const std::string h = readText(e.path().string()).substr(0, 400);
            const size_t at = h.find("element vertex ");
            if (at != std::string::npos) v += std::atoll(h.c_str() + at + 15);  // fibres off: the tubes only
        }
        return v;
    };
    const long long vm = tubeVertices(s, "medium"), vf = tubeVertices(f, "fine");
    std::printf("    exported tube vertices: medium %lld, fine %lld (x%.2f)\n", vm, vf, vm ? double(vf) / vm : 0.0);
    CHECK(vm > 0 && vf > 3 * vm);
}

// The job writes straight into a folder with files in it: the export's files lie directly in it,
// in no subfolder of its own, under the request's file name whatever the folder is called, and
// the user's files stay. A second export replaces the first's files and adds no other folder.
TEST(export_writes_into_the_chosen_folder) {
    const std::string d = freshDir("objinto"), root = freshDir("objinto_cache"), out = freshDir("objinto_out");
    ImportCacheEnv env(root);
    SceneSettings s = importPlane(d, root);
    applyFullProfile(s, *sceneEntryFor(s));
    std::string err;
    auto fin = run(BuiltScene::Kind::Knit, s, &err);
    CHECK(fin && fin->knit);
    const std::string chosen = out + "/my folder";
    writeText(chosen + "/notes.txt", "mine");
    writeText(chosen + "/old/keep.txt", "mine too");
    ExportRequest x;
    x.settings = std::make_shared<const SceneSettings>(s);
    x.knit = fin->knit;
    x.dir = chosen;
    x.name = "plane";
    x.format = "obj";
    for (int round = 0; round < 2; ++round) {
        auto ex = exportResults(x, JobControl{}, &err);
        std::printf("    export %d: %s\n", round, err.c_str());
        CHECK(ex && ex->exportedTo == fs::absolute(chosen).string());
        CHECK(fs::exists(chosen + "/plane.obj") && fs::exists(chosen + "/plane.mtl") && fs::exists(chosen + "/report.json"));
        CHECK(readText(chosen + "/notes.txt") == "mine" && readText(chosen + "/old/keep.txt") == "mine too");
        for (const auto& e : fs::directory_iterator(chosen)) {
            std::printf("      %s\n", e.path().filename().string().c_str());
            CHECK(e.path().filename().string().find(".staging-") == std::string::npos);
        }
        CHECK(fs::directory_iterator(out)->path().filename() == "my folder");
    }
}

TEST(export_obj_contents_and_coordinates) {
    const std::string d = freshDir("objx"), root = freshDir("objx_cache"), out = freshDir("objx_out");
    ImportCacheEnv env(root);
    SceneSettings s = importPlane(d, root);
    applyFullProfile(s, *sceneEntryFor(s));
    std::string err;
    auto fin = run(BuiltScene::Kind::Knit, s, &err);
    CHECK(fin && fin->knit);
    const std::string target = out + "/plane_obj";
    auto ex = run(BuiltScene::Kind::Export, s, &err, nullptr, fin->knit, nullptr, target, "obj", false);
    CHECK(ex);
    CHECK(filesUnder(target) == (std::vector<std::string>{std::string(kPublishMarker), "plane_obj.mtl", "plane_obj.obj", "report.json"}));
    // Coordinates in the source's units: every vertex within the plane's extent (60 x 40, z near 0
    // up to the yarn's depth) - the knit sits on the surface it was made from.
    std::ifstream f(target + "/plane_obj.obj");
    std::string line;
    Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    int objects = 0;
    while (std::getline(f, line)) {
        if (line.rfind("o ", 0) == 0) ++objects;
        if (line.rfind("v ", 0) != 0) continue;
        Vec3 v;
        std::istringstream ls(line.substr(2));
        ls >> v.x >> v.y >> v.z;
        lo = {std::min(lo.x, v.x), std::min(lo.y, v.y), std::min(lo.z, v.z)};
        hi = {std::max(hi.x, v.x), std::max(hi.y, v.y), std::max(hi.z, v.z)};
    }
    std::printf("    OBJ: %d objects, bounds [%.2f %.2f %.2f] .. [%.2f %.2f %.2f]\n", objects, lo.x, lo.y, lo.z, hi.x, hi.y, hi.z);
    CHECK(objects >= 1 && lo.x > -3 && hi.x < 63 && lo.y > -3 && hi.y < 43 && std::fabs(lo.z) < 5 && std::fabs(hi.z) < 5);
    CHECK(readText(target + "/plane_obj.mtl").find("newmtl wool_white") != std::string::npos);
}

// A new object opens with its texture knitted; other colours per part; without colour, one colour.
TEST(import_default_yarn_colour_source) {
    ImportAnalysis a;
    a.colourSource = "texture bunny.png";
    CHECK(defaultYarnColourSource(a) == "texture");
    a.colourSource = "vertex colours";
    CHECK(defaultYarnColourSource(a) == "object");
    a.colourSource = "material colours (Kd)";
    CHECK(defaultYarnColourSource(a) == "object");
    a.colourSource = "none";
    CHECK(defaultYarnColourSource(a) == "single");
    const std::string d = freshDir("colourdefault");
    writeText(d + "/plane.mtl", "newmtl wool\nKd 0.2 0.4 0.6\nmap_Kd tex.png\n");
    writeTexture(d + "/tex.png", 200, 100, 50);
    writeText(d + "/plane.obj", objText(grid(12, 8, 60, 40), "plane.mtl", "wool"));
    ImportedMesh m;
    std::string err;
    CHECK(readMeshFile(d + "/plane.obj", ImportLimits{}, m, &err));
    CHECK(defaultYarnColourSource(analyzeImport(m)) == "texture");
}

// An example opens with its own colour sources: the teddy (no colour texture) and the finished wool swatch take
// their objects' colours.
TEST(presets_open_with_their_texture_colours) {
    auto sources = [](const char* id) {
        const SceneEntry* e = findScene(id);
        CHECK(e != nullptr);
        SceneSettings s;
        s.scene = id;
        applySceneDefaults(s, *e);
        return s.yarnColourSource + " " + s.pileColourSource;
    };
    CHECK(sources("teddy") == "object object");
    CHECK(sources("wool-swatch").rfind("object", 0) == 0);
}

// Fur on the red-green plane: every fibre takes the texture colour at its root, on screen and in the export, split
// into one strand part and one fur material per colour; the object source keeps one colour per part and the same
// fibres (owner, 2026-11-19: the texture colour stays a choice of the Fur appearance).
TEST(fur_takes_the_texture_colour_at_its_root) {
    const std::string d = freshDir("furcol"), root = freshDir("furcol_cache"), out = freshDir("furcol_out");
    ImportCacheEnv env(root);
    SceneSettings s = importRedGreenPlane(d, root);
    s.workflow = "fur";
    s.pile.density = 30;
    s.pile.fibresDensity = 30;
    s.pileColourSource = "texture";
    std::string err;
    auto fur = run(BuiltScene::Kind::Fibres, s, &err);
    if (!fur) std::printf("    %s\n", err.c_str());
    CHECK(fur && fur->fibreResult);
    if (!fur || !fur->fibreResult) return;
    const DisplayLines& lines = fur->fibres;
    CHECK(lines.ownColor.size() == lines.partColor.size() && lines.population.size() * 3 == lines.pos.size());
    // A fibre reaches at most its length from its root: the fine fibres' 2 x median beats every loop.
    const double margin = 2 * s.pile.fibresLength + kRedGreenTexelMm;
    long red = 0, green = 0, wrong = 0;
    for (size_t v = 0; v < lines.population.size(); ++v) {
        const double r = lines.ownColor[3 * v], g = lines.ownColor[3 * v + 1], x = lines.pos[3 * v];
        if (r == 1 && g == 0) red += 1, wrong += x > kRedGreenEdgeX + margin;
        if (r == 0 && g == 1) green += 1, wrong += x < kRedGreenEdgeX - margin;
        // Away from the colour edge and the plane's sides (the texture repeats: its sides blend): one colour.
        wrong += std::fabs(x) > margin && std::fabs(x) < 30 - margin && r != 1 && g != 1;
    }
    std::printf("    fibre vertices: %ld red, %ld green, %ld wrong\n", red, green, wrong);
    CHECK(red > 0 && green > 0 && wrong == 0);
    auto package = [&](const SceneSettings& settings, const std::string& name) {
        const std::string dir = out + "/" + name;
        auto ex = run(BuiltScene::Kind::Export, settings, &err, nullptr, nullptr, fur->fibreResult, dir, "ppt", true);
        if (!ex) std::printf("    %s: %s\n", name.c_str(), err.c_str());
        CHECK(ex);
        size_t strands = 0;
        for (const std::string& f : filesUnder(dir)) {
            if (f.find(".pptstrands") == std::string::npos) continue;
            StrandSetF back;
            CHECK(readStrandsFile(dir + "/" + f, back, &err));
            strands += back.strands();
        }
        return std::make_pair(readText(dir + "/scene.json"), strands);
    };
    const auto [texture, textureStrands] = package(s, "texture");
    CHECK(texture.find("\"fur_cff0000\"") != std::string::npos && texture.find("\"fur_c00ff00\"") != std::string::npos);
    CHECK(texture.find("_fibres__cff0000") != std::string::npos && texture.find("_fibres__c00ff00") != std::string::npos);
    CHECK(texture.find("_fibres\"") == std::string::npos);  // no fibre left unsplit
    SceneSettings object = s;
    object.pileColourSource = "object";
    const auto [perPart, perPartStrands] = package(object, "object");
    CHECK(perPart.find("fur_c") == std::string::npos && perPart.find("_fibres\"") != std::string::npos);
    CHECK(textureStrands == static_cast<size_t>(fur->fibreResult->fibres) && perPartStrands == textureStrands);
}

// Every stitch takes the mean of the texture over its cell: wholly on one side of the plane's colour
// edge exactly that colour, across it a blend within one sample of the cell's true share.
TEST(stitch_colours_average_each_cell) {
    const std::string d = freshDir("cellcol"), root = freshDir("cellcol_cache");
    ImportCacheEnv env(root);
    SceneSettings s = importRedGreenPlane(d, root);
    s.yarnColourSource = "texture";
    std::string err;
    auto fin = run(BuiltScene::Kind::Knit, s, &err);
    CHECK(fin && fin->knit && fin->knit->scene);
    if (!fin || !fin->knit || !fin->knit->scene) return;
    const Scene& scene = *fin->knit->scene;
    const SceneStitchColours colours = sceneStitchColours(scene);
    CHECK(&sceneStitchColours(scene) == &sceneStitchColours(scene));  // computed once per knit
    CHECK(colours.stitchesOffSurface == 0);
    int red = 0, green = 0, blended = 0, wrong = 0;
    for (size_t i = 0; i < scene.parts.size(); ++i) {
        const ScenePart& r = scene.parts[i];
        if (!r.knit) continue;
        const double w = r.knit->t.pp.waleSpacing, h = r.knit->t.pp.courseSpacing;
        for (const StitchRecord& st : r.knit->t.curves.topo.stitches) {
            const uint32_t key = colours.parts[i].key[st.id];
            const Vec3 c = stitchColourLinear(key);
            // The cell's extent along x: its axes on the plane in either order.
            const double xU = std::fabs(st.U.x) * w / 2, x0 = st.origin.x + r.centre.x, xA = x0 + st.V.x * h;
            const double lo = std::min(x0, xA) - xU, hi = std::max(x0, xA) + xU;
            const double share = std::clamp((kRedGreenEdgeX - lo) / (hi - lo), 0.0, 1.0);
            if (c.z != 0 || std::fabs(c.x + c.y - 1) > 0.02) ++wrong;
            if (hi < kRedGreenEdgeX - kRedGreenTexelMm) {
                ++red;
                if (key != 0xff0000) ++wrong;
            } else if (lo > kRedGreenEdgeX + kRedGreenTexelMm) {
                ++green;
                if (key != 0x00ff00) ++wrong;
            } else {
                ++blended;
                if (std::fabs(c.x - share) > 1.0 / kCellSamples) ++wrong;
            }
        }
    }
    std::printf("    %d red, %d green, %d blended, %d wrong, %zu palette colours\n", red, green, blended, wrong, colours.palette.size());
    CHECK(red > 0 && green > 0 && blended > 0 && wrong == 0);
}

// Overlapped seams: two uv islands meet at the plane's colour edge, the left half red, the right half
// green. A stitch knitted from one chart past the seam lies on the other chart and takes the colour
// there, not its own chart's.
TEST(stitch_colours_across_a_chart_seam) {
    const std::string d = freshDir("seamcol"), root = freshDir("seamcol_cache");
    ImportCacheEnv env(root);
    writeText(d + "/plane.mtl", "newmtl wool\nmap_Kd tex.png\n");
    writeRedGreenTexture(d + "/tex.png");
    const Grid g = grid(12, 8, 60, 40);
    std::ostringstream obj;
    obj << "mtllib plane.mtl\n";
    for (const Vec3& p : g.pos) obj << "v " << p.x << " " << p.y << " 0\n";
    const int vertices = static_cast<int>(g.uv.size());
    for (const Vec2& t : g.uv) obj << "vt " << t.x * 0.9 << " " << t.y << "\n";          // island 1: u in [0, 0.45] at x < 30
    for (const Vec2& t : g.uv) obj << "vt " << 0.1 + t.x * 0.9 << " " << t.y << "\n";    // island 2: u in [0.55, 1] at x > 30
    obj << "usemtl wool\n";
    for (const std::array<int, 3>& t : g.tri) {
        const bool rightHalf = g.pos[t[0]].x + g.pos[t[1]].x + g.pos[t[2]].x > 90;
        obj << "f";
        for (int k : t) obj << " " << k + 1 << "/" << k + 1 + (rightHalf ? vertices : 0);
        obj << "\n";
    }
    writeText(d + "/plane.obj", obj.str());
    SceneSettings s = importWrittenFixture(d, root);
    applyFullProfile(s, *sceneEntryFor(s));
    CHECK(overlapCoverage(s));
    s.yarnColourSource = "texture";
    std::string err;
    auto fin = run(BuiltScene::Kind::Knit, s, &err);
    CHECK(fin && fin->knit && fin->knit->scene);
    if (!fin || !fin->knit || !fin->knit->scene) return;
    const Scene& scene = *fin->knit->scene;
    const SceneStitchColours colours = sceneStitchColours(scene);
    CHECK(colours.stitchesOffSurface == 0);
    int knittedParts = 0, pastTheSeam = 0, wrong = 0;
    for (size_t i = 0; i < scene.parts.size(); ++i) {
        const ScenePart& r = scene.parts[i];
        if (!r.knit) continue;
        ++knittedParts;
        const double w = r.knit->t.pp.waleSpacing, h = r.knit->t.pp.courseSpacing;
        const std::vector<StitchRecord>& stitches = r.knit->t.curves.topo.stitches;
        double meanX = 0;
        for (const StitchRecord& st : stitches) meanX += st.origin.x + r.centre.x;
        const bool leftChart = meanX < 0;
        for (const StitchRecord& st : stitches) {
            const double x0 = st.origin.x + r.centre.x, xA = x0 + st.V.x * h, xU = std::fabs(st.U.x) * w / 2;
            const double lo = std::min(x0, xA) - xU, hi = std::max(x0, xA) + xU;
            const uint32_t key = colours.parts[i].key[st.id];
            if (hi < kRedGreenEdgeX - kRedGreenTexelMm) {
                pastTheSeam += !leftChart;
                wrong += key != 0xff0000;
            } else if (lo > kRedGreenEdgeX + kRedGreenTexelMm) {
                pastTheSeam += leftChart;
                wrong += key != 0x00ff00;
            }
        }
    }
    std::printf("    %d knitted parts, %d stitches wholly past the seam, %d wrong, %zu palette colours\n", knittedParts, pastTheSeam, wrong,
                colours.palette.size());
    CHECK(knittedParts >= 2 && pastTheSeam > 0 && wrong == 0);
}

// A latitude-longitude sphere of radius R about c as a Grid: (lat + 1) x (lon + 1) vertices, uv = (j / lon,
// 1 - i / lat). The import's weld closes the seam column and the poles, so it reads back closed.
Grid sphereGrid(double R, int lat, int lon, const Vec3& c) {
    Grid g;
    for (int i = 0; i <= lat; ++i)
        for (int j = 0; j <= lon; ++j) {
            const double th = kPi * i / lat, ph = 2 * kPi * j / lon;
            g.pos.push_back(c + Vec3{R * std::sin(th) * std::cos(ph), R * std::sin(th) * std::sin(ph), R * std::cos(th)});
            g.uv.push_back({double(j) / lon, 1.0 - double(i) / lat});
        }
    for (int i = 0; i < lat; ++i)
        for (int j = 0; j < lon; ++j) {
            const int a = i * (lon + 1) + j, b = a + 1, d = a + lon + 1, e = d + 1;
            g.tri.push_back({a, d, e});
            g.tri.push_back({a, e, b});
        }
    return g;
}

Grid joinedGrids(const Grid& a, const Grid& b) {
    Grid g = a;
    const int base = static_cast<int>(a.pos.size());
    g.pos.insert(g.pos.end(), b.pos.begin(), b.pos.end());
    g.uv.insert(g.uv.end(), b.uv.begin(), b.uv.end());
    for (const std::array<int, 3>& t : b.tri) g.tri.push_back({base + t[0], base + t[1], base + t[2]});
    return g;
}

int partsWithLayout(const Scene& scene, KnitLayout layout) {
    int n = 0;
    for (const ScenePart& r : scene.parts) n += r.layout == layout;
    return n;
}

// Imports knit on charts only (2026-11-11): two closed spheres sunk into each other, as the
// pieces of a toy are, beside an open plane, all on the file's own mesh (import.surface = file, as scripted runs keep
// it). Before, the spheres knitted as two spirals. Only the parts are loaded; nothing is knitted.
TEST(import_keep_file_knits_on_charts) {
    const std::string d = freshDir("round"), root = freshDir("round_cache");
    ImportCacheEnv env(root);
    const Grid spheres = joinedGrids(sphereGrid(30, 16, 32, {120, 20, 35}), sphereGrid(30, 16, 32, {165, 20, 35}));
    writeText(d + "/plane.obj", objText(joinedGrids(grid(12, 8, 60, 40), spheres), "", ""));
    SceneSettings s = importWrittenFixture(d, root);
    CHECK(s.importSurface == "file");
    const SceneEntry* e = sceneEntryFor(s);
    Scene t;
    std::string err;
    CHECK(e && loadPartsToKnit(*e, s, SourceOptions{}, t, &err));
    if (!err.empty()) std::printf("    %s\n", err.c_str());
    std::printf("    %d spiral, %d chart parts\n", partsWithLayout(t, KnitLayout::Spiral), partsWithLayout(t, KnitLayout::UvTiling));
    CHECK(partsWithLayout(t, KnitLayout::Spiral) == 0 && partsWithLayout(t, KnitLayout::UvTiling) >= 3);
}

// The stitch colours are in the display whatever the source (the viewer shows them when Colourwork is
// chosen, with no job, and the source is in no cache key); a new tessellation keeps them.
TEST(stitch_colours_do_not_depend_on_the_source) {
    const std::string d = freshDir("sourcecol"), root = freshDir("sourcecol_cache");
    ImportCacheEnv env(root);
    SceneSettings s = importRedGreenPlane(d, root);
    applyFullProfile(s, *sceneEntryFor(s));
    s.yarnColourSource = "object";
    std::string err;
    auto fin = run(BuiltScene::Kind::Knit, s, &err);
    CHECK(fin && fin->knit);
    if (!fin || !fin->knit) return;
    auto colouredYarnMeshes = [](const BuiltScene& b) {
        int coloured = 0, yarn = 0;
        for (const DisplayMesh& m : b.meshes) {
            if (m.kind != DisplayMesh::Kind::Yarn) continue;
            ++yarn;
            coloured += !m.vcol.empty() && m.vcol.size() == m.pos.size();
        }
        return yarn > 0 && coloured == yarn;
    };
    CHECK(colouredYarnMeshes(*fin));
    SceneSettings texture = s;
    texture.yarnColourSource = "texture";
    CHECK(sceneCacheKey(texture, 1) == sceneCacheKey(s, 1));
    SceneSettings fine = texture;
    fine.tubeSides = 16;
    fine.tubeRing = 0.2;
    StrandMeshRequest j;
    j.settings = std::make_shared<const SceneSettings>(fine);
    j.shown = fin;
    auto rebuilt = rebuildStrandMesh(j, JobControl{}, &err);
    CHECK(rebuilt && colouredYarnMeshes(*rebuilt));
}

// A display vertex in a stitch colour lies on that colour's side: the stitch lookup along the strand
// is right on the single tube and on separate ply tubes, for the tubes and the oriented fibres.
TEST(stitch_colours_along_tubes_and_plies) {
    const std::string d = freshDir("plycol"), root = freshDir("plycol_cache");
    ImportCacheEnv env(root);
    SceneSettings s = importRedGreenPlane(d, root);
    applyFullProfile(s, *sceneEntryFor(s));
    s.yarnColourSource = "texture";
    s.fibres = true;
    s.fibreDensity = 20;
    for (int plies : {1, 3}) {
        s.plies = plies;
        std::string err;
        auto fin = run(BuiltScene::Kind::Knit, s, &err);
        CHECK(fin && fin->knit);
        if (!fin || !fin->knit) return;
        auto fib = run(BuiltScene::Kind::Fibres, s, &err, nullptr, fin->knit);
        CHECK(fib);
        if (!fib) return;
        const double margin = 3 * s.wale;
        long tubeRed = 0, tubeGreen = 0, tubeWrong = 0;
        for (const DisplayMesh& m : fin->meshes) {
            if (m.kind != DisplayMesh::Kind::Yarn) continue;
            CHECK(m.vcol.size() == m.pos.size());
            for (size_t v = 0; v + 2 < m.pos.size(); v += 3) {
                const double r = m.vcol[v], g = m.vcol[v + 1], x = m.pos[v];
                if (r == 1 && g == 0) tubeRed += 1, tubeWrong += x > kRedGreenEdgeX + margin;
                if (r == 0 && g == 1) tubeGreen += 1, tubeWrong += x < kRedGreenEdgeX - margin;
            }
        }
        const DisplayLines& lines = fib->fibres;
        long fibreRed = 0, fibreGreen = 0, fibreWrong = 0;
        for (size_t v = 0; v < lines.population.size(); ++v) {
            if (lines.population[v] != 0) continue;
            const double r = lines.ownColor[3 * v], g = lines.ownColor[3 * v + 1], x = lines.pos[3 * v];
            if (r == 1 && g == 0) fibreRed += 1, fibreWrong += x > kRedGreenEdgeX + margin;
            if (r == 0 && g == 1) fibreGreen += 1, fibreWrong += x < kRedGreenEdgeX - margin;
        }
        std::printf("    plies %d: tube vertices %ld red, %ld green, %ld wrong; fibre vertices %ld red, %ld green, %ld wrong\n", plies, tubeRed,
                    tubeGreen, tubeWrong, fibreRed, fibreGreen, fibreWrong);
        CHECK(tubeRed > 0 && tubeGreen > 0 && tubeWrong == 0);
        CHECK(fibreRed > 0 && fibreGreen > 0);
        CHECK(fibreWrong * 200 < fibreRed + fibreGreen);  // loose fibres reach past a stitch; under 0.5 %
    }
}

// The export: one yarn material per palette colour, tubes and oriented fibres split by colour, on the
// single tube and on ply tubes; stray fibres one part per object part in the same material as with
// one colour per part.
TEST(export_texture_colourwork) {
    const std::string d = freshDir("xcol"), root = freshDir("xcol_cache"), out = freshDir("xcol_out");
    ImportCacheEnv env(root);
    SceneSettings s = importRedGreenPlane(d, root);
    applyFullProfile(s, *sceneEntryFor(s));
    s.fibres = true;
    s.fibreDensity = 20;
    auto packageOf = [&](const SceneSettings& settings, const std::string& name) {
        std::string err;
        auto fin = run(BuiltScene::Kind::Knit, settings, &err);
        CHECK(fin && fin->knit);
        auto fib = run(BuiltScene::Kind::Fibres, settings, &err, nullptr, fin->knit);
        CHECK(fib && fib->fibreResult);
        const std::string dir = out + "/" + name;
        auto ex = run(BuiltScene::Kind::Export, settings, &err, nullptr, fin->knit, fib->fibreResult, dir, "ppt", true);
        CHECK(ex);
        if (!ex) std::printf("    %s: %s\n", name.c_str(), err.c_str());
        return readText(dir + "/scene.json");
    };
    // The material entry of the stray fibres with this id prefix, from its type on.
    auto strayMaterial = [](const std::string& json, const std::string& prefix) {
        const size_t at = json.find("{\"id\": \"" + prefix);
        if (at == std::string::npos) return std::string();
        const size_t type = json.find("\"type\"", at), endLine = json.find('\n', at);
        return json.substr(type, endLine - type);
    };
    auto count = [](const std::string& json, const std::string& what) {
        int n = 0;
        for (size_t at = json.find(what); at != std::string::npos; at = json.find(what, at + 1)) ++n;
        return n;
    };
    s.yarnColourSource = "object";
    const std::string oneColour = packageOf(s, "object");
    s.yarnColourSource = "texture";
    for (int plies : {1, 3}) {
        s.plies = plies;
        const std::string json = packageOf(s, "texture_" + std::to_string(plies));
        CHECK(json.find("yarn_cff0000") != std::string::npos && json.find("yarn_c00ff00") != std::string::npos);
        CHECK(json.find("fibre_oriented_yarn_cff0000") != std::string::npos && json.find("fibre_oriented_yarn_c00ff00") != std::string::npos);
        CHECK(json.find("_fibres_oriented__cff0000") != std::string::npos && json.find("_fibres_oriented__c00ff00") != std::string::npos);
        CHECK(json.find("_fibres_random__cff0000") != std::string::npos && json.find("_fibres_random__c00ff00") != std::string::npos);
        CHECK(json.find("_fibres_random\"") == std::string::npos);  // no stray part left unsplit
        const std::string red = strayMaterial(json, "fibre_random_yarn_cff0000"), green = strayMaterial(json, "fibre_random_yarn_c00ff00");
        CHECK(red.find("\"type\": \"hair\"") == 0 && red.find("\"color\": [1.000000, 0.000100, 0.000100]") != std::string::npos);
        CHECK(green.find("\"type\": \"hair\"") == 0 && green.find("\"color\": [0.000100, 1.000000, 0.000100]") != std::string::npos);
        std::printf("    plies %d: %d yarn colour mentions, %d stray parts\n", plies, count(json, "\"fibre_oriented_yarn_c"),
                    count(json, "_fibres_random__"));
    }
    // Match yarn colour off: the strays keep their own colour, one part per object part.
    s.fibreRandomInheritsYarn = false;
    const std::string own = packageOf(s, "texture_own");
    CHECK(own.find("_fibres_random__") == std::string::npos && own.find("_fibres_random\"") != std::string::npos);
    const std::string ownStray = strayMaterial(own, "fibre_random_");
    CHECK(ownStray.find("\"type\": \"hair\"") == 0 && ownStray.find("\"color\": [0.550000, 0.550000, 0.550000]") != std::string::npos);
    CHECK(!strayMaterial(oneColour, "fibre_random_").empty());
}

// A part without a texture takes its object colour beside textured parts; with no texture on any
// part the source behaves as one colour per part. Neither stops the export.
TEST(stitch_colours_without_a_texture) {
    const std::string d = freshDir("mixcol"), root = freshDir("mixcol_cache"), out = freshDir("mixcol_out");
    ImportCacheEnv env(root);
    writeText(d + "/plane.mtl", "newmtl wool\nKd 0.2 0.4 0.6\nmap_Kd tex.png\nnewmtl plain\nKd 0.1 0.2 0.3\n");
    writeRedGreenTexture(d + "/tex.png");
    const Grid g = grid(12, 8, 60, 40);
    std::ostringstream obj;
    obj << objText(g, "plane.mtl", "wool");
    const int n = static_cast<int>(g.pos.size());
    for (const Vec3& p : g.pos) obj << "v " << p.x + 80 << " " << p.y << " 0\n";
    obj << "usemtl plain\n";
    for (const std::array<int, 3>& t : g.tri) {
        obj << "f";
        for (int k : t) obj << " " << k + n + 1 << "/" << k + 1;
        obj << "\n";
    }
    writeText(d + "/plane.obj", obj.str());
    SceneSettings s = importWrittenFixture(d, root);
    applyFullProfile(s, *sceneEntryFor(s));
    s.yarnColourSource = "texture";
    std::string err;
    auto fin = run(BuiltScene::Kind::Knit, s, &err);
    CHECK(fin && fin->knit && fin->knit->scene);
    if (!fin || !fin->knit || !fin->knit->scene) return;
    const Scene& scene = *fin->knit->scene;
    const SceneStitchColours colours = sceneStitchColours(scene);
    const uint32_t plainKey = stitchColourKey({0.1, 0.2, 0.3});
    int plainStitches = 0, otherColours = 0;
    for (size_t i = 0; i < scene.parts.size(); ++i) {
        if (scene.source.prims[scene.parts[i].prim].albedo) continue;
        for (uint32_t k : colours.parts[i].key) ++plainStitches, otherColours += k != plainKey;
    }
    std::printf("    %d stitches on the plain part, %d of another colour; %zu palette colours\n", plainStitches, otherColours, colours.palette.size());
    CHECK(plainStitches > 0 && otherColours == 0 && colours.palette.size() > 1);
    CHECK(run(BuiltScene::Kind::Export, s, &err, nullptr, fin->knit, nullptr, out + "/mixed", "ppt", false));
    // No texture at all.
    const std::string d2 = freshDir("mixcol2");
    writeText(d2 + "/plane.mtl", "newmtl plain\nKd 0.1 0.2 0.3\n");
    writeText(d2 + "/plane.obj", objText(grid(12, 8, 60, 40), "plane.mtl", "plain"));
    SceneSettings plain = importWrittenFixture(d2, root);
    applyFullProfile(plain, *sceneEntryFor(plain));
    plain.yarnColourSource = "texture";
    auto fin2 = run(BuiltScene::Kind::Knit, plain, &err);
    CHECK(fin2 && fin2->knit && fin2->knit->scene);
    if (!fin2 || !fin2->knit || !fin2->knit->scene) return;
    CHECK(sceneStitchColours(*fin2->knit->scene).palette.empty());
    CHECK(run(BuiltScene::Kind::Export, plain, &err, nullptr, fin2->knit, nullptr, out + "/plain", "ppt", false));
}

// Unknitted surfaces (material.plain_colour_source; the fur's base, unknitted parts): "texture"
// writes the source texture and the colour uvs (every vertex's uv is its position on the plane
// fixture: uv = (x / w, y / h)), "average" one mean colour, "material" the material colour; the OBJ
// conversion carries the texture as map_Kd.
// An unticked component reaches the export as a mesh with its texture, in the OBJ as in the package.
TEST(export_selection_keeps_the_plain_node) {
    const std::string d = freshDir("selection_export"), root = freshDir("selection_export_cache"), out = freshDir("selection_export_out");
    ImportCacheEnv env(root);
    writeTexture(d + "/picture.png", 200, 100, 50);
    writeText(d + "/two.mtl", "newmtl flat\nKd 0.2 0.3 0.4\nnewmtl picture\nmap_Kd picture.png\n");
    writeText(d + "/two.obj", "mtllib two.mtl\nvt 0 0\nvt 1 0\nvt 1 1\nvt 0 1\n"
                              "v 0 0 0\nv 60 0 0\nv 60 40 0\nv 0 40 0\nv 100 0 0\nv 160 0 0\nv 160 40 0\nv 100 40 0\n"
                              "o near\nusemtl flat\nf 1/1 2/2 3/3 4/4\no far\nusemtl picture\nf 5/1 6/2 7/3 8/4\n");
    auto file = std::make_shared<ImportedMesh>();
    std::string err, log;
    CHECK(readMeshFile(d + "/two.obj", ImportLimits{}, *file, &err));
    finishNodes(*file);
    ImportOptions o;
    o.up = UpAxis::Z;
    o.plainNodes = {"far"};
    const uint64_t id = importIdentity(*file, o);
    std::vector<std::string> changes;
    CHECK(prepareAndStore({file, o, ImportLimits{}, id}, changes, &err, {}));
    SceneSettings s;
    s.scene = "import";
    s.importSource = d + "/two.obj";
    s.importIdentity = identityHex(id);
    s.importUp = "z";
    s.importPlainNodes = "far";
    s.workflow = "fur";  // the quickest way to the writer: the plain mesh takes the same path under a knit
    const SceneEntry* e = sceneEntryFor(s);
    CHECK(e != nullptr);
    Scene t;
    CHECK(loadBaseParts(*e, s, SourceOptions{}, t, &err));
    ExportOptions objExport;
    objExport.format = "obj";
    objExport.folder = out + "/obj";
    objExport.name = "two";
    CHECK(exportScene({*e, t, s}, objExport, &log, &err) || (std::printf("    %s\n", err.c_str()), false));
    const std::string obj = readText(objExport.folder + "/two.obj"), mtl = readText(objExport.folder + "/two.mtl");
    CHECK(obj.find("o two_plain_picture") != std::string::npos);
    CHECK(mtl.find("map_Kd textures/colour_0.png") != std::string::npos && fs::exists(objExport.folder + "/textures/colour_0.png"));
}

TEST(export_unknitted_colour_sources) {
    const std::string d = freshDir("plaincol"), root = freshDir("plaincol_cache"), out = freshDir("plaincol_out");
    ImportCacheEnv env(root);
    SceneSettings s = importPlane(d, root);  // texture (200, 100, 50)
    s.workflow = "fur";
    const SceneEntry* e = sceneEntryFor(s);
    CHECK(e != nullptr);
    std::string err, log;
    Scene t;
    CHECK(loadBaseParts(*e, s, SourceOptions{}, t, &err));  // the Fur path: nothing knitted, every part a body
    CHECK(t.parts.size() == 1);
    const std::string part = "meshes/" + t.parts[0].name + ".ply";
    auto plyFloats = [](const std::string& path, int& props) {
        std::ifstream f(path, std::ios::binary);
        std::string line;
        size_t n = 0;
        props = 0;
        while (std::getline(f, line) && line != "end_header") {
            if (line.rfind("element vertex ", 0) == 0) n = std::stoul(line.substr(15));
            if (line.rfind("property float ", 0) == 0) ++props;
        }
        std::vector<float> v(n * props);
        f.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(v.size() * sizeof(float)));
        return v;
    };
    const double lin[3] = {srgbDecode(200 / 255.0), srgbDecode(100 / 255.0), srgbDecode(50 / 255.0)};
    char mean[96];
    std::snprintf(mean, sizeof mean, "\"albedo\": [%.3f, %.3f, %.3f]", lin[0], lin[1], lin[2]);
    for (const char* src : {"texture", "average", "material"}) {
        s.plainColourSource = src;
        ExportOptions o;
        o.folder = out + "/" + src;
        CHECK(exportScene({*e, t, s}, o, &log, &err));
        const std::string json = readText(o.folder + "/scene.json");
        int props = 0;
        const std::vector<float> v = plyFloats(o.folder + "/" + part, props);
        CHECK(!v.empty());
        if (std::string(src) == "texture") {
            CHECK(json.find("\"albedo\": [1.000, 1.000, 1.000], \"albedo_tex\": \"textures/colour_0.png\"") != std::string::npos);
            int w = 0, h = 0;
            std::vector<uint8_t> rgb;
            CHECK(readPngRgb8(o.folder + "/textures/colour_0.png", w, h, rgb) && w == 16 && h == 16);
            CHECK(rgb.size() == 16 * 16 * 3 && rgb[0] == 200 && rgb[1] == 100 && rgb[2] == 50);
            CHECK(props == 8);
            if (props != 8) continue;
            // Scene axes (x, z, -y) of the library: u runs with x, v against z.
            float x0 = 1e30f, x1 = -1e30f, z0 = 1e30f, z1 = -1e30f;
            for (size_t i = 0; i < v.size(); i += 8) x0 = std::min(x0, v[i]), x1 = std::max(x1, v[i]), z0 = std::min(z0, v[i + 2]), z1 = std::max(z1, v[i + 2]);
            double worst = 0;
            for (size_t i = 0; i < v.size(); i += 8)
                worst = std::max({worst, double(std::fabs(v[i + 6] - (v[i] - x0) / (x1 - x0))),
                                  double(std::fabs(v[i + 7] - (z1 - v[i + 2]) / (z1 - z0)))});
            std::printf("    texture: %zu vertices, worst uv vs position %.2e\n", v.size() / 8, worst);
            CHECK(worst < 1e-5);
            // OBJ: the texture as map_Kd beside the .obj, the uvs as vt.
            ExportOptions objExport;
            objExport.format = "obj";
            objExport.folder = out + "/obj";
            objExport.name = "plane";
            CHECK(exportScene({*e, t, s}, objExport, &log, &err));
            CHECK(readText(objExport.folder + "/plane.mtl").find("map_Kd textures/colour_0.png") != std::string::npos);
            CHECK(fs::exists(objExport.folder + "/textures/colour_0.png"));
            CHECK(readText(objExport.folder + "/plane.obj").find("\nvt ") != std::string::npos);
        } else {
            // average: the mean colour (the uniform texture's); material: an import's material colour
            // is its mean colour too. No texture, no uv.
            CHECK(json.find("albedo_tex\": \"textures/colour_") == std::string::npos);
            CHECK(!fs::exists(o.folder + "/textures/colour_0.png"));
            CHECK(props == 6);
            CHECK(json.find(mean) != std::string::npos);
            CHECK(json.find(std::string("material\": \"") + (std::string(src) == "average" ? "mean_" : "plain_")) != std::string::npos);
        }
    }
}

// A new Export step starts on the Blender file, or on the package where Blender is missing, and then
// the Blender entry says why it is off.
TEST(export_default_format) {
    const ExportFormat* blend = findExportFormat("blend");
    CHECK(blend != nullptr);
    CHECK(exportFormats().front().id == "blend" && !blend->unavailable.empty() && blend->fibres);
    std::printf("    Blender %s, default %s\n", blend->available ? "found" : "not found", defaultExportFormat().c_str());
    CHECK(defaultExportFormat() == (blend->available ? "blend" : "ppt"));
    const std::optional<std::string> chosen = envVar("KNIT_BLENDER");
    if (chosen && !fs::exists(*chosen)) CHECK(!blend->available);
}

// The smallest package the Blender conversion reads: one quad, a camera, the sky; `extra` goes into
// scene.json.
void writeSmallPackage(const std::string& dir, const std::string& extra) {
    writeText(dir + "/scene.json",
              "{\"project_name\": \"small\", \"materials\": [{\"id\": \"m\", \"type\": \"diffuse\", \"diffuse\": {\"albedo\": [0.5, 0.5, 0.5]}}],"
              "\"meshes\": [{\"name\": \"quad\", \"path\": \"meshes/quad.ply\", \"material\": \"m\"}],"
              "\"instances\": [{\"id\": \"i\", \"name\": \"quad\", \"asset_mesh_index\": 0, \"material\": {\"id\": \"m\"},"
              "\"transform\": {\"translate\": [0, 0, 0], \"rotate_euler_deg\": [0, 0, 0], \"scale\": [1, 1, 1]}}],"
              "\"cameras\": [{\"id\": \"c\", \"name\": \"front\", \"position\": [0, 0, 3], \"look_at\": [0, 0, 0], \"up\": [0, 1, 0], \"fov_deg\": 30,"
              "\"width\": 64, \"height\": 64, \"spp\": 1, \"max_bounces\": 2}], \"default_camera_id\": \"c\","
              "\"env\": {\"mode\": \"analytical\", \"sky_color\": [1, 1, 1], \"sky_intensity\": 1}" + extra + "}\n");
    std::string ply = "ply\nformat binary_little_endian 1.0\nelement vertex 4\n";
    for (const char* name : {"x", "y", "z", "nx", "ny", "nz"}) ply += std::string("property float ") + name + "\n";
    ply += "element face 2\nproperty list uchar uint vertex_indices\nend_header\n";
    for (int i = 0; i < 4; ++i) {
        const std::vector<float> v = {float(i & 1), float(i >> 1), 0, 0, 0, 1};
        ply.append(reinterpret_cast<const char*>(v.data()), v.size() * sizeof(float));
    }
    for (const std::vector<uint32_t> t : {std::vector<uint32_t>{0, 1, 3}, std::vector<uint32_t>{0, 3, 2}}) {
        ply += char(3);
        ply.append(reinterpret_cast<const char*>(t.data()), 3 * sizeof(uint32_t));
    }
    writeText(dir + "/meshes/quad.ply", ply);
}

// Real Blender: the .blend is written (zstd compressed), Blender's own output is not left in the
// folder that gets published.
TEST(blender_export_writes_a_blend) {
    if (!findExportFormat("blend")->available) SKIP("no Blender 5.1 on this machine");
    const std::string package = freshDir("blend_pkg"), out = freshDir("blend_out");
    writeSmallPackage(package, "");
    std::string log, err;
    CHECK(writeBlendFile(package, {out, "small"}, log, &err));
    std::printf("    %s%s\n", log.c_str(), err.c_str());
    std::ifstream blend(out + "/small.blend", std::ios::binary);
    unsigned char magic[4] = {};
    blend.read(reinterpret_cast<char*>(magic), 4);
    CHECK(magic[0] == 0x28 && magic[1] == 0xB5 && magic[2] == 0x2F && magic[3] == 0xFD);
    size_t files = 0;
    for (const fs::directory_entry& entry : fs::directory_iterator(out)) files += entry.is_regular_file() ? 1 : 0;
    CHECK(files == 1);
}

// The whole export job with Blender: the package is staged, converted, removed, and only the .blend
// and the verdicts are published.
TEST(blender_export_job_publishes_the_blend) {
    if (!findExportFormat("blend")->available) SKIP("no Blender 5.1 on this machine");
    const std::string d = freshDir("blendjob"), root = freshDir("blendjob_cache"), out = freshDir("blendjob_out");
    ImportCacheEnv env(root);
    SceneSettings s = importPlane(d, root);
    applyFullProfile(s, *sceneEntryFor(s));
    std::string err;
    auto knitted = run(BuiltScene::Kind::Knit, s, &err);
    CHECK(knitted && knitted->knit);
    const std::string target = out + "/plane_blend";
    auto exported = run(BuiltScene::Kind::Export, s, &err, nullptr, knitted->knit, nullptr, target, "blend", false);
    std::printf("    %s\n", err.c_str());
    CHECK(exported != nullptr);
    CHECK(filesUnder(target) == (std::vector<std::string>{std::string(kPublishMarker), "plane_blend.blend", "report.json"}));
}

// A cancelled job stops Blender at once: the Cancelled reaches the job, and nothing holds the folder
// when the publisher removes the staging directory.
TEST(blender_export_cancel_stops_blender) {
    if (!findExportFormat("blend")->available) SKIP("no Blender 5.1 on this machine");
    const std::string package = freshDir("blendcancel_pkg"), out = freshDir("blendcancel_out");
    writeSmallPackage(package, "");
    std::string log;
    CancelToken token;
    token.cancel();
    ThreadCancelScope scope(&token);
    bool cancelled = false;
    std::string err;
    try {
        (void)writeBlendFile(package, {out, "small"}, log, &err);
    } catch (const Cancelled&) {
        cancelled = true;
    }
    CHECK(cancelled && !fs::exists(out + "/small.blend"));
    std::error_code ec;
    fs::remove_all(out, ec);
    CHECK(!ec);
}

// A package the script refuses: the export fails with the script's reason, not a bare exit code.
TEST(blender_export_reports_the_script_message) {
    if (!findExportFormat("blend")->available) SKIP("no Blender 5.1 on this machine");
    const std::string package = freshDir("blendbad_pkg"), out = freshDir("blendbad_out");
    writeSmallPackage(package, ", \"lights\": [{\"id\": \"l\"}]");
    std::string err, log;
    CHECK(!writeBlendFile(package, {out, "small"}, log, &err));
    std::printf("    %s\n", err.c_str());
    CHECK(err.find("area lights") != std::string::npos);
    CHECK(!fs::exists(out + "/small.blend") && !fs::exists(out + "/.blender_output.txt"));
    // The script reads a wool fibre part's shape from a package the exporter wrote: it refuses the
    // one shape it cannot put into the radii. (Radii on a real .blend: checked by reading one back.)
    const std::string fibres = freshDir("blendshape_pkg");
    ExportScene ex(fibres);
    fs::create_directories(fibres + "/strands");
    fs::create_directories(fibres + "/meshes");
    StrandSetF strand;
    strand.add({0, 0, 0, 0.001f, 0, 0}, {1e-5f, 1e-5f});
    CHECK(ex.addStrands("piece_fibres_random", ex.randomFibreMaterial("wool_white"), strand, &err));
    CHECK(writeSceneJson(ex, &err));
    std::string written = readText(fibres + "/scene.json");
    const std::string imported = "\"thickness\": \"imported\"";
    const size_t shape = written.find(imported);
    CHECK(shape != std::string::npos);
    writeText(fibres + "/scene.json", written.replace(shape, imported.size(), "\"thickness\": \"uniform\""));
    err.clear();
    CHECK(!writeBlendFile(fibres, {out, "small"}, log, &err));
    CHECK(err.find("part_fibre_shapes") != std::string::npos);
}

// The OBJ conversion with meshes with and without uv in one package: v / vt / vn keep separate
// running indices. Every face index is parsed, range-checked, and resolved back to the source
// vertex (each vertex's x, u and nx encode its global id, so a wrong offset is caught).
TEST(obj_export_mixed_attributes) {
    for (const std::vector<bool> uv : {std::vector<bool>{false, true}, std::vector<bool>{true, false, true}}) {
        const std::string pkg = freshDir("obj_mixed_pkg"), out = freshDir("obj_mixed_out");
        fs::create_directories(pkg + "/meshes");
        ExportScene ex(pkg);
        MaterialSpec grey;
        grey.color = {0.5, 0.5, 0.5};
        ex.baseMaterials = {{"m", grey}};
        std::string err, log;
        int id = 0;
        for (size_t m = 0; m < uv.size(); ++m) {
            SceneMesh mesh;
            for (int i = 0; i < 4; ++i, ++id) {
                mesh.pos.push_back({double(id), double(i & 1), 0});
                mesh.nrm.push_back({float(id) / 100.f, 0, 1});
                if (uv[m]) mesh.uv.insert(mesh.uv.end(), {float(id), 0.25f});
            }
            mesh.tri = {0, 1, 2, 2, 1, 3};
            CHECK(ex.addMeshPart("p" + std::to_string(m), "m", mesh, &err));
        }
        CHECK(writeObjFile(ex, {out, "mixed"}, log, &err));
        std::ifstream f(out + "/mixed.obj");
        std::string line;
        std::vector<double> vx, vtu, vnx;
        int faces = 0, textured = 0;
        while (std::getline(f, line)) {
            std::istringstream ls(line);
            std::string tag;
            ls >> tag;
            double a = 0;
            if (tag == "v") ls >> a, vx.push_back(a);
            else if (tag == "vt") ls >> a, vtu.push_back(a);
            else if (tag == "vn") ls >> a, vnx.push_back(a);
            else if (tag == "f") {
                ++faces;
                std::string corner;
                int withUv = 0;
                for (int k = 0; k < 3; ++k) {
                    CHECK(static_cast<bool>(ls >> corner));
                    const size_t s1 = corner.find('/'), s2 = corner.find('/', s1 + 1);
                    CHECK(s1 != std::string::npos && s2 != std::string::npos);
                    const long iv = std::stol(corner.substr(0, s1)), in = std::stol(corner.substr(s2 + 1));
                    CHECK(iv >= 1 && iv <= long(vx.size()) && in >= 1 && in <= long(vnx.size()));
                    const double gid = vx[iv - 1];
                    CHECK(std::fabs(vnx[in - 1] * 100 - gid) < 1e-3);
                    if (s2 > s1 + 1) {
                        const long it = std::stol(corner.substr(s1 + 1, s2 - s1 - 1));
                        CHECK(it >= 1 && it <= long(vtu.size()));
                        CHECK(vtu[it - 1] == gid);  // the vt of THIS vertex
                        ++withUv;
                    }
                }
                CHECK(withUv == 0 || withUv == 3);
                textured += withUv == 3;
            }
        }
        size_t nuv = 0;
        for (bool u : uv) nuv += u;
        CHECK(faces == int(2 * uv.size()) && textured == int(2 * nuv));
        CHECK(vx.size() == 4 * uv.size() && vnx.size() == vx.size() && vtu.size() == 4 * nuv);
    }
}

// A settings file that enables only some Final steps (here: no surface check) is not Final: the
// knit result says so and the export (UI rule and backend) refuses it.
TEST(custom_profile_is_not_final) {
    const std::string d = freshDir("custom_profile"), root = freshDir("custom_profile_cache"), out = freshDir("custom_profile_out");
    ImportCacheEnv env(root);
    SceneSettings s = importPlane(d, root);
    const SceneEntry& e = *sceneEntryFor(s);
    applyFullProfile(s, e);
    CHECK(isFinalProfile(s, e) && !isPreviewProfile(s));
    s.surfaceGate = s.surfaceHalf = false;  // custom: settle and fallback, but no gate (d)
    CHECK(!isFinalProfile(s, e) && !isPreviewProfile(s));
    std::string err;
    auto k = run(BuiltScene::Kind::Knit, s, &err);
    CHECK(k && k->knit && !k->knit->final);
    const ExportEligibility el = exportEligibility(s, k->knit.get(), nullptr, FibreChoice::LeftOut);
    CHECK(el.state == ExportEligibility::State::Blocked && el.reasons[0].find("custom profile") != std::string::npos);
    CHECK(!run(BuiltScene::Kind::Export, s, &err, nullptr, k->knit, nullptr, out + "/x", "ppt", false));
    CHECK(err.find("custom profile") != std::string::npos && !fs::exists(out + "/x"));
    // Through the scene cache as well (a restored result keeps the predicate).
    const std::string cache = freshDir("custom_profile_scenes");
    auto a = run(BuiltScene::Kind::Knit, s, &err, nullptr, nullptr, nullptr, "", "ppt", false, cache);
    auto b = run(BuiltScene::Kind::Knit, s, &err, nullptr, nullptr, nullptr, "", "ppt", false, cache);
    CHECK(a && b && b->cached && !a->knit->final && !b->knit->final);
}

TEST(export_eligibility_rows) {
    SceneSettings s;
    s.scene = "teddy";
    KnitResult k;
    k.id = 7;
    k.settings = s;
    k.final = true;
    k.parts = 3;
    FibreResult f;
    f.knitId = 7;
    f.key = fibreResultKey(s);
    using St = ExportEligibility::State;
    CHECK(exportEligibility(s, nullptr, nullptr, FibreChoice::LeftOut).state == St::Blocked);
    CHECK(exportEligibility(s, &k, nullptr, FibreChoice::LeftOut).state == St::Ready);
    k.final = false;
    CHECK(exportEligibility(s, &k, nullptr, FibreChoice::LeftOut).state == St::Blocked);
    k.final = true;
    k.computational = false;
    k.viewOnlyReason = "test";
    CHECK(exportEligibility(s, &k, nullptr, FibreChoice::LeftOut).state == St::Blocked);
    k.computational = true;
    CHECK(exportEligibility(s, &k, nullptr, FibreChoice::Included).state == St::Blocked);   // fibres included, none
    CHECK(exportEligibility(s, &k, &f, FibreChoice::Included).state == St::Ready);
    f.knitId = 6;
    CHECK(exportEligibility(s, &k, &f, FibreChoice::Included).state == St::Blocked);        // fibres of another knit
    f.knitId = 7;
    SceneSettings t = s;
    t.wale = 6;
    CHECK(exportEligibility(t, &k, &f, FibreChoice::Included).state == St::Blocked);        // knit settings changed
    t = s;
    t.flatten = true;
    CHECK(exportEligibility(t, &k, &f, FibreChoice::Included).state == St::Blocked);        // the final surface changed
    t = s;
    t.thicknessPreset = "custom";                                          // a label only
    CHECK(exportEligibility(t, &k, &f, FibreChoice::Included).state == St::Ready);
    t = s;
    t.yarnBrown.color = {0, 0, 0};
    t.fibres = !s.fibres;                                                  // materials, include flag
    CHECK(exportEligibility(t, &k, &f, FibreChoice::Included).state == St::Ready);
    k.failing = {"a (FAIL)"};
    CHECK(exportEligibility(s, &k, &f, FibreChoice::Included).state == St::Diagnostic);
    k.failing.clear();
    k.cross = Status::Unresolved;
    CHECK(exportEligibility(s, &k, &f, FibreChoice::Included).state == St::Diagnostic);
}

// The backend refuses what the UI would block; a diagnostic export is written, and marked so.
TEST(export_job_enforces_eligibility) {
    const std::string d = freshDir("enforce"), root = freshDir("enforce_cache"), out = freshDir("enforce_out");
    ImportCacheEnv env(root);
    SceneSettings s = importPlane(d, root);
    std::string err;
    auto draft = run(BuiltScene::Kind::Knit, s, &err);
    CHECK(draft && draft->knit);
    CHECK(!run(BuiltScene::Kind::Export, s, &err, nullptr, draft->knit, nullptr, out + "/x", "ppt"));
    CHECK(err.find("Draft") != std::string::npos && !fs::exists(out + "/x"));
    auto k = std::make_shared<KnitResult>();
    k->settings = s;
    k->final = true;
    k->parts = 1;
    k->failing = {"p (FAIL)"};
    k->scene = draft->knit->scene;
    CHECK(run(BuiltScene::Kind::Export, s, &err, nullptr, k, nullptr, out + "/y", "ppt", false));
    CHECK(readText(out + "/y/report.json").find("\"export_status\": \"diagnostic\"") != std::string::npos);
    CHECK(!run(BuiltScene::Kind::Fibres, s, &err, nullptr, nullptr));
    CHECK(err.find("knit first") != std::string::npos);
}

// The tiled knit cache: exact round trip and verification; damaged data is refused.
TEST(tiled_knit_cache_roundtrip) {
    const std::string d = freshDir("tkc"), root = freshDir("tkc_cache"), cache = freshDir("tkc_scenes");
    ImportCacheEnv env(root);
    SceneSettings s = importPlane(d, root);
    std::string err;
    auto b = run(BuiltScene::Kind::Knit, s, &err, nullptr, nullptr, nullptr, "", "ppt", false, cache);
    CHECK(b && !b->cached && b->knit->computational);
    auto c = run(BuiltScene::Kind::Knit, s, &err, nullptr, nullptr, nullptr, "", "ppt", false, cache);
    CHECK(c && c->cached && c->knit->computational && c->knit->digest == b->knit->digest && c->knit->id != b->knit->id);
    CHECK(c->knit->scene->restored && !c->knit->scene->verified);
    CHECK(verifyKnit(*c->knit->scene, s, &err));
    // Fibres from the restored knit: the job verifies it first and does not knit.
    StageLog log;
    s.fibreDensity = 10;
    auto f = run(BuiltScene::Kind::Fibres, s, &err, &log, c->knit);
    CHECK(f && f->fibreResult);
    // A damaged knit.bin is refused by the restore itself (the record's checksum), wherever the
    // damage is: in the curves, in report text no verification reads, in the checksum.
    const std::string bytes = serializeKnit(*b->knit->scene);
    Scene t;
    CHECK(restoreKnit(bytes, *sceneEntryFor(s), s, t, &err));
    for (const size_t at : {size_t(20), bytes.size() / 4, bytes.size() / 2, bytes.size() * 3 / 4, bytes.size() - 1}) {
        std::string bad = bytes;
        bad[at] ^= 0x10;
        Scene u;
        CHECK(!restoreKnit(bad, *sceneEntryFor(s), s, u, &err));
        CHECK(err.find("damaged") != std::string::npos || err.find("version") != std::string::npos);
    }
    CHECK(!restoreKnit(bytes.substr(0, bytes.size() - 3), *sceneEntryFor(s), s, t, &err));  // truncated
}

// Caches written before the spiral pieces joined the common path (knit.bin v2) stay usable: a v2
// record (the v3 record without the per-part validator flag) restores and verifies. A cache entry
// whose knit cannot be restored is shown for viewing only when the scene opens, and "Knit again"
// (KnitRequest::knitAgain) knits and replaces it instead of showing the same entry again.
TEST(tiled_knit_cache_previous_version_and_knit_again) {
    const std::string d = freshDir("tkv"), root = freshDir("tkv_cache"), cache = freshDir("tkv_scenes");
    ImportCacheEnv env(root);
    SceneSettings s = importPlane(d, root);
    s.coverage = "component";  // one part, as the knits written in v2 (before knit.coverage = overlap became the default)
    std::string err;
    auto b = run(BuiltScene::Kind::Knit, s, &err, nullptr, nullptr, nullptr, "", "ppt", false, cache);
    CHECK(b && b->knit->computational && b->knit->scene->parts.size() == 1 && b->knit->scene->parts[0].knit);
    const PartKnit& k = *b->knit->scene->parts[0].knit;
    CHECK(!k.gateAOnCpu);
    // The v2 layout: version 2, no gap patch count after it (v4), no validator flag after the gate (a)
    // digest, the record's checksum.
    std::string v2 = serializeKnit(*b->knit->scene);
    const size_t versionAt = std::string("tiled-knit v1\n").size();
    const int64_t two = 2;
    std::memcpy(&v2[versionAt], &two, sizeof two);
    v2.erase(versionAt + sizeof two, sizeof(int64_t));
    const uint64_t after = fnv1a64(validationSnapshot(k.after));
    const size_t digestAt = v2.find(std::string(reinterpret_cast<const char*>(&after), sizeof after));
    CHECK(digestAt != std::string::npos);
    v2.erase(digestAt + sizeof after, sizeof(int64_t));
    v2.resize(v2.size() - sizeof(uint64_t));
    const uint64_t sum = fnv1a64(v2);
    v2.append(reinterpret_cast<const char*>(&sum), sizeof sum);
    Scene t;
    CHECK(restoreKnit(v2, *sceneEntryFor(s), s, t, &err));
    CHECK(t.parts.size() == 1 && t.parts[0].knit && !t.parts[0].knit->gateAOnCpu);
    CHECK(verifyKnit(t, s, &err));
    // An entry whose knit.bin cannot be restored (as one written in a format no longer read): view
    // only on open, knitted again on "Knit again".
    const uint64_t asset = sceneAssetHash(*sceneEntryFor(s));
    const std::string unreadable = std::string(serializeKnit(*b->knit->scene)).substr(0, 20);
    CHECK(storeSceneCache(cache, s, asset, *b, &err, &unreadable));
    auto c = run(BuiltScene::Kind::Knit, s, &err, nullptr, nullptr, nullptr, "", "ppt", false, cache);
    CHECK(c && c->cached && !c->knit->computational);
    KnitRequest j;
    j.settings = std::make_shared<const SceneSettings>(s);
    j.cacheRoot = cache;
    j.knitAgain = true;
    auto again = knitObject(j, JobControl{}, &err);
    CHECK(again && !again->cached && again->knit->computational);
    auto reopened = run(BuiltScene::Kind::Knit, s, &err, nullptr, nullptr, nullptr, "", "ppt", false, cache);
    CHECK(reopened && reopened->cached && reopened->knit->computational && reopened->knit->digest == again->knit->digest);
}

// The viewer's --load-knit (KnitRequest::loadKnitCache) is used or the job fails: a directory that
// holds no knit of these pieces refuses the run (nothing is knitted in its place), and a scene knitted
// on its mesh has no per-piece knit files to load.
TEST(viewer_load_knit_is_used_or_fails) {
    if (!gpu::available(nullptr)) SKIP("no usable CUDA GPU");
    const SceneEntry* teddy = findScene("teddy");
    CHECK(teddy);
    if (!teddy) return;
    knit_test::requireScene(findProjectFile(teddy->piecesManifest));
    KnitRequest j;
    j.settings = std::make_shared<const SceneSettings>(viewerDefaults(*teddy));
    j.loadKnitCache = freshDir("load_knit_empty");
    std::string err;
    CHECK(knitObject(j, JobControl{}, &err) == nullptr && err.find("knit cache refused") != std::string::npos);
    const std::string d = freshDir("load_knit_tiled"), root = freshDir("load_knit_tiled_cache");
    ImportCacheEnv env(root);
    j.settings = std::make_shared<const SceneSettings>(importPlane(d, root));
    err.clear();
    CHECK(knitObject(j, JobControl{}, &err) == nullptr && err.find("per-piece knit cache files") != std::string::npos);
}

// knit.coverage = overlap: gap patches are parts the knitting adds after its first pass, which
// loading the scene does not make. The cache stores their seeds and makes them again, so a knit
// with gap patches restores (user, 2026-10-29: the man's overlapped knit opened "view only").
TEST(tiled_knit_cache_gap_patches) {
    const std::string d = freshDir("tkg"), root = freshDir("tkg_cache"), cache = freshDir("tkg_scenes");
    ImportCacheEnv env(root);
    SceneSettings s = importPlane(d, root);
    std::string err;
    auto b = run(BuiltScene::Kind::Knit, s, &err, nullptr, nullptr, nullptr, "", "ppt", false, cache);
    CHECK(b && b->knit->computational && !b->knit->scene->parts.empty());
    Scene& scene = *b->knit->scene;
    ScenePart first;  // a copy of the identity: the push_back below may move the parts
    first.componentName = scene.parts.front().componentName;
    first.prim = scene.parts.front().prim;
    first.component = scene.parts.front().component;
    first.centre = scene.parts.front().centre;
    first.mesh = scene.parts.front().mesh;
    ScenePart patch;
    patch.name = first.componentName + "_gapx";  // after the knit's own gap patches (if any)
    patch.prim = first.prim;
    patch.component = first.component;
    patch.componentName = first.componentName;
    patch.centre = first.centre;
    patch.mesh = first.mesh;
    patch.gapPatch = true;
    patch.patch.tri = 1;
    patch.patch.bary[0] = 0.2;
    patch.patch.bary[1] = 0.3;
    patch.patch.bary[2] = 0.5;
    patch.patch.wale = Vec3{0, 0, 1};
    patch.patch.radius = 12.5;
    scene.parts.push_back(std::move(patch));
    const std::string bytes = serializeKnit(scene);
    Scene t;
    CHECK(restoreKnit(bytes, *sceneEntryFor(s), s, t, &err));
    CHECK(t.parts.size() == scene.parts.size());
    const ScenePart& got = t.parts.back();
    CHECK(got.gapPatch && got.name == first.componentName + "_gapx" && got.prim == first.prim && got.component == first.component);
    CHECK(got.patch.tri == 1 && got.patch.bary[0] == 0.2 && got.patch.bary[1] == 0.3 && got.patch.bary[2] == 0.5);
    CHECK(got.patch.wale.z == 1 && got.patch.radius == 12.5 && meshDigest(*got.mesh) == meshDigest(*first.mesh));
    CHECK(verifyKnit(t, s, &err));
}

// ------------------------------------------------------------ lifecycle ----

namespace {

// Stops the ImportSession worker once at `point` ("start" / "finish") until release().
struct HookGate {
    explicit HookGate(std::string p) : point(std::move(p)) {}
    std::function<void(const char*)> hook() {
        return [this](const char* p) {
            std::unique_lock<std::mutex> lk(m);
            if (point != p || open) return;
            at = true;
            cv.notify_all();
            cv.wait(lk, [&] { return open; });
        };
    }
    void waitReached() {
        std::unique_lock<std::mutex> lk(m);
        CHECK(cv.wait_for(lk, std::chrono::seconds(30), [&] { return at; }));
    }
    void release() {
        std::lock_guard<std::mutex> lk(m);
        open = true;
        cv.notify_all();
    }
    std::string point;
    std::mutex m;
    std::condition_variable cv;
    bool at = false, open = false;
};

struct ReleaseOnExit {
    std::vector<HookGate*> gates;
    ~ReleaseOnExit() {
        for (HookGate* g : gates) g->release();
    }
};

bool waitUntil(const std::function<bool()>& f, double seconds = 30) {
    const auto end = std::chrono::steady_clock::now() + std::chrono::duration<double>(seconds);
    while (!f()) {
        if (std::chrono::steady_clock::now() > end) return false;
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }
    return true;
}

} // namespace

// reset() while a read runs: neither its progress nor its late result may overwrite Idle.
// (Gates are declared before the session and released by a guard declared after it: a failing
// CHECK unwinds through the release before the session joins its worker - no hang.)
// One triangle per named object: a file of exactly `nodes` components.
std::string objWithNodes(const std::string& dir, int nodes) {
    std::string obj = "v 0 0 0\nv 1 0 0\nv 0 1 0\n";
    for (int n = 0; n < nodes; ++n) obj += "o node" + std::to_string(n) + "\nf 1 2 3\n";
    const std::string path = dir + "/nodes" + std::to_string(nodes) + ".obj";
    writeText(path, obj);
    return path;
}

// The safety limit on components, on load: 200 are fine, 201 are refused with a message and no stage 2.
TEST(import_session_refuses_too_many_nodes) {
    const std::string d = freshDir("node_limit");
    ImportSession s;
    s.open(objWithNodes(d, 200));
    CHECK(waitUntil([&] { return s.snapshot().phase != ImportSession::Phase::Reading; }));
    CHECK(s.snapshot().phase == ImportSession::Phase::Analysed && s.snapshot().mesh->nodes.size() == 200);
    s.open(objWithNodes(d, 201));
    CHECK(waitUntil([&] { return s.snapshot().phase != ImportSession::Phase::Reading; }));
    const ImportSession::Snapshot refused = s.snapshot();
    CHECK(refused.phase == ImportSession::Phase::Failed && !refused.mesh && !refused.preview);
    CHECK(refused.error.find("nodes201.obj holds 201 components; the import limit is 200.") != std::string::npos);
}

TEST(import_session_reset_race) {
    const std::string d = freshDir("session_reset");
    writeText(d + "/plane.obj", objText(grid(12, 8, 60, 40), "", "", false));
    using Ph = ImportSession::Phase;
    HookGate atFinish("finish"), atStart("start");
    ImportSession s;
    ReleaseOnExit guard{{&atFinish, &atStart}};
    // (a) Reset after the work, before its result is published.
    s.setTestHook(atFinish.hook());
    s.open(d + "/plane.obj");
    atFinish.waitReached();
    s.reset();
    const int gen = s.snapshot().generation;
    const Ph afterReset = s.snapshot().phase;
    atFinish.release();
    CHECK(afterReset == Ph::Idle);
    CHECK(waitUntil([&] { return !s.busy(); }));
    const ImportSession::Snapshot stayed = s.snapshot();
    CHECK(stayed.phase == Ph::Idle && !stayed.mesh && !stayed.analysis && !stayed.preview && stayed.generation == gen);
    // (b) Reset before the work starts: its progress reports are dropped too.
    s.setTestHook(atStart.hook());
    s.open(d + "/plane.obj");
    atStart.waitReached();
    s.reset();
    atStart.release();
    CHECK(waitUntil([&] { return !s.busy(); }));
    const ImportSession::Snapshot idle = s.snapshot();
    CHECK(idle.phase == Ph::Idle && idle.stage.empty() && idle.fraction == 0 && !idle.mesh);
    // The session still works afterwards.
    s.setTestHook(nullptr);
    s.open(d + "/plane.obj");
    CHECK(waitUntil([&] { return s.snapshot().phase == Ph::Analysed; }));
    CHECK(s.snapshot().mesh && s.snapshot().preview);
}

// cancel() of a preparation that is still queued leaves Preparing at once (back to Analysed),
// and the retired running one publishes nothing.
TEST(import_session_cancel_queued_prepare) {
    const std::string d = freshDir("session_cancel"), root = freshDir("session_cancel_cache");
    ImportCacheEnv env(root);
    writeText(d + "/plane.obj", objText(grid(12, 8, 60, 40), "", "", false));
    using Ph = ImportSession::Phase;
    HookGate atStart("start");
    ImportSession s;
    ReleaseOnExit guard{{&atStart}};
    s.open(d + "/plane.obj");
    CHECK(waitUntil([&] { return s.snapshot().phase == Ph::Analysed; }));
    ImportOptions o;
    o.unitMm = 1;
    o.up = UpAxis::Z;
    s.setTestHook(atStart.hook());
    s.prepare(o);  // taken by the worker, held at "start"
    atStart.waitReached();
    s.prepare(o);  // queued behind it
    const bool queued = s.snapshot().phase == Ph::Preparing && s.busy();
    const int gen = s.snapshot().generation;
    s.cancel();
    // Immediately, the worker still held.
    const Ph afterCancel = s.snapshot().phase;
    const bool keptMesh = s.snapshot().mesh != nullptr;
    const int genAfter = s.snapshot().generation;
    atStart.release();
    CHECK(queued);
    CHECK(afterCancel == Ph::Analysed && keptMesh && genAfter > gen);
    CHECK(waitUntil([&] { return !s.busy(); }));
    const ImportSession::Snapshot back = s.snapshot();
    CHECK(back.phase == Ph::Analysed && back.identity == 0 && back.error.empty());
    // A later preparation completes normally.
    s.setTestHook(nullptr);
    s.prepare(o);
    CHECK(waitUntil([&] { return s.snapshot().phase == Ph::Prepared; }));
    CHECK(s.snapshot().identity != 0);
}

// Whether an analysed file is the stored import its settings name (the viewer's revalidation): the same
// options and a stored entry match; other options, another identity or a deleted entry do not.
TEST(import_session_stored_match) {
    const std::string d = freshDir("session_stored"), root = freshDir("session_stored_cache");
    ImportCacheEnv env(root);
    writeText(d + "/plane.obj", objText(grid(12, 8, 60, 40), "", "", false));
    using Ph = ImportSession::Phase;
    ImportSession s;
    s.open(d + "/plane.obj");
    CHECK(waitUntil([&] { return s.snapshot().phase == Ph::Analysed; }));
    ImportOptions o;
    o.unitMm = 1;
    o.up = UpAxis::Z;
    CHECK(!s.storedImportMatches(o, ""));  // nothing stored yet
    s.prepare(o);
    CHECK(waitUntil([&] { return s.snapshot().phase == Ph::Prepared; }));
    const uint64_t id = s.snapshot().identity;
    s.open(d + "/plane.obj");  // read again, as the viewer does for an import its settings name
    CHECK(waitUntil([&] { return s.snapshot().phase == Ph::Analysed; }));
    CHECK(s.storedImportMatches(o, identityHex(id)));
    ImportOptions other = o;
    other.unitMm = 2;
    CHECK(!s.storedImportMatches(other, identityHex(id)));
    CHECK(!s.storedImportMatches(o, identityHex(id + 1)));
    std::error_code ec;
    std::filesystem::remove_all(importCacheDir(root, id), ec);
    CHECK(!s.storedImportMatches(o, identityHex(id)));
}

// cancel(): any running job stops, the pending one is dropped, nothing restarts, and the last
// completed result stays the complete one.
TEST(builder_cancel_any_job_drops_pending) {
    std::atomic<int> started{0};
    SceneBuilder b([&](const JobRequest&, const JobControl& c, std::string*) -> std::shared_ptr<BuiltScene> {
        ++started;
        if (c.revision == 1) {
            auto r = std::make_shared<BuiltScene>();
            r->revision = 1;
            return r;
        }
        for (int i = 0; i < 2000; ++i) {
            cancelCheckpoint();
            std::this_thread::sleep_for(std::chrono::milliseconds(5));
        }
        auto r = std::make_shared<BuiltScene>();
        r->revision = c.revision;
        return r;
    });
    const auto settings = std::make_shared<const SceneSettings>();
    auto running = [&b] { return b.snapshot().progress.revision; };
    CHECK(b.startKnit(KnitRequest{settings}) == 1);
    while (b.snapshot().busy) std::this_thread::sleep_for(std::chrono::milliseconds(2));
    CHECK(b.snapshot().latestComplete && b.snapshot().latestComplete->revision == 1);
    CHECK(b.startFibres(FibreRequest{settings}) == 2);
    while (running() != 2) std::this_thread::sleep_for(std::chrono::milliseconds(1));
    // An export queued behind it, then cancel: both go.
    ExportRequest x;
    x.settings = settings;
    CHECK(b.startExport(x) == 3);
    while (running() != 3 && running() != -1) std::this_thread::sleep_for(std::chrono::milliseconds(1));
    const auto t0 = std::chrono::steady_clock::now();
    b.cancel();
    const JobSnapshot cancelling = b.snapshot();
    CHECK(cancelling.cancelRequested || !cancelling.busy);
    while (b.snapshot().busy) std::this_thread::sleep_for(std::chrono::milliseconds(1));
    const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    const JobSnapshot end = b.snapshot();
    CHECK(end.lastOutcome == JobOutcome::Cancelled && !end.cancelRequested && ms < 1000);
    CHECK(end.latestComplete->revision == 1);
    const int n = started;
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    CHECK(started == n && !b.snapshot().busy);  // nothing restarts by itself
}

// The thickness presets' flat-swatch evidence as an expected matrix, the expected FAILs included
// (thickness_presets.h): Dense passes stockinette, reverse and ribs and fails garter (clearance);
// Medium and Thin pass every pattern. The teddy (an ordinary mesh) offers the same presets.
// The Knit step's yarn look (yarn_look.h): the default is one tube with 3-ply appearance maps;
// ply tubes clear the maps and flattening; returning restores both, never a silent loss.
TEST(yarn_look_roundtrip) {
    SceneSettings s;
    YarnLookMemory mem;
    YarnLook l = yarnLookOf(s, mem);
    CHECK(l.construction == YarnConstruction::SingleTube && l.appearancePlies == 3 && l.style == PlyProfileModel::Swept);
    CHECK(s.plies == 1 && s.textureProfile == "3-swept");
    // A single-tube look with flattening and a circles style.
    l.appearancePlies = 4;
    l.style = PlyProfileModel::Circles;
    applyYarnLook(s, l, mem);
    CHECK(s.plies == 1 && s.textureProfile == "4-circles");
    s.flatten = true;
    s.flattenInflate = 0.22;
    // -> 3 ply tubes: no maps, no flattening.
    l = yarnLookOf(s, mem);
    l.construction = YarnConstruction::PlyTubes;
    l.tubePlies = 3;
    applyYarnLook(s, l, mem);
    CHECK(s.plies == 3 && s.textureProfile == "none" && !s.flatten);
    CHECK(yarnLookConflict(s).empty());
    l = yarnLookOf(s, mem);
    CHECK(l.construction == YarnConstruction::PlyTubes && l.tubePlies == 3 && l.appearancePlies == 4 && l.style == PlyProfileModel::Circles);  // remembered
    // Changing the tube count keeps the memory.
    l.tubePlies = 2;
    applyYarnLook(s, l, mem);
    CHECK(s.plies == 2);
    // -> back to the single tube: the maps and the flattening return.
    l = yarnLookOf(s, mem);
    l.construction = YarnConstruction::SingleTube;
    applyYarnLook(s, l, mem);
    CHECK(s.plies == 1 && s.textureProfile == "4-circles" && s.flatten && s.flattenInflate == 0.22);
    CHECK(yarnLookOf(s, mem).tubePlies == 2);
    // Appearance off is the plain tube.
    l = yarnLookOf(s, mem);
    l.appearancePlies = 0;
    applyYarnLook(s, l, mem);
    CHECK(s.textureProfile == "none" && yarnLookOf(s, mem).appearancePlies == 0);
    // A settings file with ply tubes AND a profile is reported, not dropped silently.
    SceneSettings f;
    f.plies = 3;
    f.textureProfile = "2-swept";
    CHECK(yarnLookConflict(f).find("ignored") != std::string::npos);
}

// Fibre strands and fibre maps (yarn_look.h): the yarn body on the single drawn yarn. It keeps the single tube's maps in
// the settings (unused) and the flattening, so switching to it and back changes no knitting key: the knit
// stays current (review 2026-11-11). Ply tubes clear the maps and the flattening into the memory; fibres
// after ply tubes restore both. Conflicting settings files are said.
TEST(yarn_look_fibres_construction) {
    SceneSettings s;
    YarnLookMemory mem;
    s.flatten = true;
    s.flattenInflate = 0.3;
    s.textureProfile = "3-swept";
    const uint64_t singleTubeKey = geometryKey(s);
    YarnLook l = yarnLookOf(s, mem);
    l.construction = YarnConstruction::FibreStrands;
    applyYarnLook(s, l, mem);
    CHECK(s.yarnBody.enabled && s.yarnBody.render == "fibres" && s.plies == 1 && s.textureProfile == "3-swept" && s.flatten && s.flattenInflate == 0.3);
    CHECK(geometryKey(s) == singleTubeKey);  // not Outdated: the knit, its surface and its fuzz stay
    CHECK(yarnConstructionOf(s) == YarnConstruction::FibreStrands && yarnLookConflict(s).empty());
    CHECK(yarnLookOf(s, mem).appearancePlies == 3);
    // -> fibre maps: the same body rendered as maps (yarn_body.render); the knit stays current here too.
    l = yarnLookOf(s, mem);
    l.construction = YarnConstruction::FibreMaps;
    applyYarnLook(s, l, mem);
    CHECK(s.yarnBody.enabled && s.yarnBody.render == "maps" && s.plies == 1 && geometryKey(s) == singleTubeKey);
    CHECK(yarnConstructionOf(s) == YarnConstruction::FibreMaps && yarnOfFibres(YarnConstruction::FibreMaps) && !yarnOfFibres(YarnConstruction::PlyTubes));
    l = yarnLookOf(s, mem);
    l.construction = YarnConstruction::FibreStrands;
    applyYarnLook(s, l, mem);
    CHECK(s.yarnBody.render == "fibres" && yarnConstructionOf(s) == YarnConstruction::FibreStrands);
    // -> ply tubes: the body is off, the maps and the flattening go into the memory.
    l = yarnLookOf(s, mem);
    l.construction = YarnConstruction::PlyTubes;
    applyYarnLook(s, l, mem);
    CHECK(!s.yarnBody.enabled && s.plies == 3 && !s.flatten && s.textureProfile == "none");
    // -> fibres again: one drawn yarn, the flattening and the kept maps back.
    l = yarnLookOf(s, mem);
    l.construction = YarnConstruction::FibreStrands;
    applyYarnLook(s, l, mem);
    CHECK(s.yarnBody.enabled && s.plies == 1 && s.flatten && s.flattenInflate == 0.3 && s.textureProfile == "3-swept");
    // -> single tube: the same settings as at the start.
    l = yarnLookOf(s, mem);
    l.construction = YarnConstruction::SingleTube;
    applyYarnLook(s, l, mem);
    CHECK(!s.yarnBody.enabled && s.plies == 1 && s.textureProfile == "3-swept" && s.flatten);
    CHECK(geometryKey(s) == singleTubeKey);
    // Settings files: a body with ply tubes is refused by the export; a body with maps keeps them, unused.
    SceneSettings f;
    f.yarnBody.enabled = true;
    f.plies = 2;
    CHECK(yarnLookConflict(f).find("refused") != std::string::npos);
    f.plies = 1;
    f.textureProfile = "3-swept";
    CHECK(yarnLookConflict(f).empty());
}

// The thickness presets (thickness_presets.h): one list for every non-block scene, the man and the
// teddy alike. Each preset's flat reference is knitted as the scene knits it (flatPatternReference of
// the settings: the contact diameter of yarnFit, the same for every pattern, and the pattern's relaxed
// shape, applyPatternShape); Dense and Medium knit the same stockinette as Thin's calibrated 1.16
// (within the 0.001 rounding of k). The expected verdicts are the measured ones (2026-11-11): the
// stockinette stitch is the curved fit limit, which a bare flat template does not clear (the curved
// pieces pass after settle), so stockinette and reverse FAIL the clearance; garter at its relaxed 0.70
// row height also folds (it knits on a surface after settle); rib 1x1's two layers
// PASS; rib 2x2's spread pairs FAIL the crossing signature (not offered in the viewer).
TEST(thickness_presets_flat_matrix) {
    const SceneEntry charts = uvChartScene();
    const SceneEntry* man = &charts;
    const SceneEntry* teddy = findScene("teddy");
    CHECK(man && teddy);
    using S = Status;
    struct Expected {
        S verdict;
        const char* fail;  // failCategoryText of a FAIL
    };
    const Expected row[5] = {{S::Fail, "clearance"}, {S::Fail, "clearance"}, {S::Fail, "curvature+clearance"}, {S::Pass, "none"},
                             {S::Fail, "signature+crossing"}};
    const std::vector<ThicknessPreset> mp = thicknessPresets(*man), tp = thicknessPresets(*teddy);
    CHECK(mp.size() == 3 && tp.size() == 3);
    for (size_t k = 0; k < mp.size() && k < tp.size(); ++k)
        CHECK(mp[k].id == tp[k].id && mp[k].label == tp[k].label && mp[k].D == tp[k].D && mp[k].w == tp[k].w && mp[k].h == tp[k].h);
    CHECK(mp[0].D == 1.86 && mp[1].D == 1.50 && mp[2].D == 1.16);
    int checked = 0, unexpected = 0;
    double stockinetteKnit[3] = {0, 0, 0};
    for (size_t i = 0; i < mp.size() && i < 3; ++i) {
        const ThicknessPreset& p = mp[i];
        for (int k = 0; k < 5; ++k) {
            const StitchPattern pattern = static_cast<StitchPattern>(k);
            SceneSettings s;
            applyThicknessPreset(s, p);
            s.pattern = toString(pattern);
            const YarnFit f = yarnFit(s);
            if (pattern == StitchPattern::Stockinette) stockinetteKnit[i] = knitDiameter(s, p.D);
            const PatternReference r = flatPatternReference(s, pattern);
            std::printf("    %s %s: knit %.3f (k %.3f%s): %s (%s)\n", p.id.c_str(), toString(pattern), knitDiameter(s, p.D), f.k,
                        f.tooTight ? ", too tight" : "", toString(r.verdict), failCategoryText(r.fail).c_str());
            // The whole matrix is printed before any verdict is checked.
            if (r.verdict != row[k].verdict || failCategoryText(r.fail) != row[k].fail) ++unexpected;
            ++checked;
        }
    }
    CHECK(checked == 15 && unexpected == 0);
    CHECK(std::fabs(stockinetteKnit[0] - 1.16) < 0.006 && std::fabs(stockinetteKnit[1] - 1.16) < 0.006 && stockinetteKnit[2] == 1.16);
    // Presets apply exactly; editing a value makes the settings Custom.
    SceneSettings s;
    applyThicknessPreset(s, thicknessPresets(*man)[1]);
    CHECK(matchingPreset(s, *man) && matchingPreset(s, *man)->id == "medium");
    s.wale += 0.1;
    CHECK(!matchingPreset(s, *man));
    // Each preset's fallback is the calibrated one scaled to its drawn D, so its knitted fallbacks
    // stay below its knit diameter (a scene-wide list in drawn units left Medium and Thin without
    // any: teddy 12 and 11 of 18 instead of 14). Final restores the matching preset's list; Draft
    // keeps none.
    for (const ThicknessPreset& p : mp) {
        CHECK(p.fallback.size() == 4);
        SceneSettings f = viewerDefaults(*teddy);
        applyFullProfile(f, *teddy);
        applyThicknessPreset(f, p);
        CHECK(f.diameterFallback == p.fallback && isFinalProfile(f, *teddy));
        const double k = knitDiameter(f, p.D) / p.D;
        CHECK(std::fabs(knitDiameter(f, p.fallback.front()) - 1.12) < 0.01 && p.fallback.front() * k < knitDiameter(f, p.D));
        applyPreviewProfile(f);
        applyThicknessPreset(f, p);
        CHECK(f.diameterFallback.empty() && isPreviewProfile(f));
        applyFullProfile(f, *teddy);
        CHECK(f.diameterFallback == p.fallback);
    }
}

// thickYarnWarning: warned only when even the compressed yarn does not fit the stitch (any scene).
TEST(thick_yarn_warning) {
    // Gauge fit (yarnFit): a stitch too tight for a round yarn is knitted at
    // a compressed contact diameter and drawn at D; keys change only when the yarn is compressed.
    const SceneEntry charts = uvChartScene();
    const SceneEntry* man = &charts;
    const SceneEntry* teddy = findScene("teddy");
    CHECK(man && teddy);
    auto at = [](double D, double w, double h) {
        SceneSettings s;
        s.diameter = D;
        s.wale = w;
        s.course = h;
        return s;
    };
    auto key = [](const SceneSettings& s) { return stageKey(s, Stage::Knit, 0); };
    // Gauges that fit: nothing changes (the S1 teddy pieces' 1.16 gauge, the Dense preset).
    for (const SceneSettings& s : {at(1.16, 5.5, 3.4), at(1.10, 5.5, 3.4), at(1.0, 6.5, 4.6)}) {
        CHECK(yarnFit(s).k == 1 && !yarnFit(s).tooTight);
        CHECK(yarnFitNote(s).empty() && thickYarnWarning(s, *teddy).empty());
        CHECK(knitDiameter(s, s.diameter) == s.diameter && !renderFlattened(s) && renderInflate(s) == 0);
        SceneSettings off = s;
        off.compressToFit = false;
        CHECK(key(off) == key(s));
    }
    // The man at the user gauge: compressed to about 66 %, drawn at 1.10 mm, no FAIL warning.
    SceneSettings m = at(1.10, 3.61, 2.52);
    const YarnFit fm = yarnFit(m);
    std::printf("    man 1.10 at 3.61 x 2.52: need %.4f k %.3f; %s\n", fm.need, fm.k, yarnFitNote(m).c_str());
    CHECK(fm.compressed() && fm.k > 0.66 && fm.k < 0.67 && !fm.tooTight);
    CHECK(thickYarnWarning(m, *man).empty() && yarnFitNote(m).find("0.73") != std::string::npos);
    const FinalSurfaceSpec fs = finalSurface(m, knitDiameter(m, m.diameter));
    CHECK(fs.compression && std::fabs(fs.renderedDiameter() - 1.10) < 1e-9 && fs.valid(nullptr));
    m.flatten = true;  // flattening's own inflation stays on top of the compression
    CHECK(std::fabs(finalSurface(m, knitDiameter(m, m.diameter)).renderedDiameter() - 1.10 * 1.15) < 1e-9);
    m.flatten = false;
    // Compression off: knitted round, warned, and a different key.
    SceneSettings mOff = m;
    mOff.compressToFit = false;
    CHECK(yarnFit(mOff).k == 1 && yarnFit(mOff).tooTight && key(mOff) != key(m));
    CHECK(thickYarnWarning(mOff, *man).find("compression is off") != std::string::npos);
    // Separate ply tubes knit the same compressed yarn, drawn round at its full diameter.
    SceneSettings mPly = m;
    mPly.plies = 2;
    CHECK(yarnFit(mPly).k == fm.k && thickYarnWarning(mPly, *man).empty());  // plies knit the same yarn
    // 1.46 mm at 5.5 x 3.4 (the user's run on the teddy): compressed to about 80 %, no FAIL warning.
    const YarnFit ft = yarnFit(at(1.46, 5.5, 3.4));
    CHECK(ft.k > 0.79 && ft.k < 0.81 && !ft.tooTight);
    // Beyond the least compression: k stops at kMinCompression and the warning stays.
    const SceneSettings x = at(1.6, 3.0, 2.0);
    std::printf("    %s\n", thickYarnWarning(x, *man).c_str());
    CHECK(yarnFit(x).k == kMinCompression && yarnFit(x).tooTight && !thickYarnWarning(x, *man).empty());
    // The fallback diameters are compressed by the same factor (pieces options).
    SceneSettings mf = at(1.10, 3.61, 2.52);
    mf.diameterFallback = {1.08, 1.04};
    CHECK(std::fabs(knitDiameter(mf, 1.04) - 1.04 * yarnFit(mf).k) < 1e-12);
}

// The teddy's pieces are drawn (and exported) at the drawn diameter of each preset, not at the knit
// diameter of a compressed yarn: Dense 1.86 knitted at 1.16 is drawn 1.86 thick. A spiral piece is a
// part of the common knit result, rendered by renderedPart as every part.
TEST(teddy_pieces_drawn_diameter) {
    const SceneEntry* teddy = findScene("teddy");
    CHECK(teddy);
    knit_test::requireScene(findProjectFile(teddy->piecesManifest));
    for (const ThicknessPreset& preset : thicknessPresets(*teddy)) {
        SceneSettings s;
        applySceneDefaults(s, *teddy);
        applyThicknessPreset(s, preset);
        s.textureProfile = "none";
        // One straight yarn along x at the knit parameters, as the spiral layout delivers a piece.
        ScenePart piece;
        piece.name = "piece";
        piece.knit = std::make_unique<PartKnit>();
        YarnPath y;
        for (int k = 0; k < 12; ++k) y.curve.ctrl.push_back({4.0 * k, 0, 0});
        y.finalize();
        CurveSet& cs = piece.knit->t.curves;
        cs.yarns.push_back(y);
        cs.kind = CurveKind::Surface;
        PatternParams& pp = piece.knit->t.pp;
        pp.yarnDiameter = knitDiameter(s, s.diameter);
        pp.waleSpacing = s.wale;
        pp.courseSpacing = s.course;
        std::string err;
        const std::optional<RenderedPart> rendered = renderedPart(piece, s, &err);
        CHECK(rendered.has_value());
        if (!rendered) continue;
        const double drawn = rendered->strandPp.yarnDiameter;
        std::printf("    %-6s knit %.3f mm, drawn %.3f mm (preset D %.2f)\n", preset.label.c_str(), pp.yarnDiameter, drawn, preset.D);
        CHECK(std::fabs(drawn - preset.D) < 0.02 * preset.D);
    }
}

// Separate ply tubes (strand.plies > 1): the viewer draws the teddy's pieces as their ply tubes at
// the ply diameter (renderedPart: the tubes the export writes), not the yarn's single tube. The
// smallest piece of a Draft knit with 3 plies, round tubes: every ring of its displayed tube mesh
// (sides + 1 vertices, the seam duplicated; yarn_surface.cpp buildTubeMesh) lies at the ply radius from
// the ring's centre.
TEST(teddy_pieces_ply_tubes) {
    if (!gpu::available(nullptr)) SKIP("no usable CUDA GPU");
    const SceneEntry* teddy = findScene("teddy");
    CHECK(teddy);
    knit_test::requireScene(findProjectFile(teddy->piecesManifest));
    SceneSettings s = viewerDefaults(*teddy);
    s.plies = 3;
    s.textureProfile = "none";  // round tubes: no flattening, no shading-map displacement
    s.flatten = false;
    std::string err;
    auto b = run(BuiltScene::Kind::Knit, s, &err);
    if (!b || !b->knit || !b->knit->scene) {
        std::printf("    knit: %s\n", err.c_str());
        CHECK(false);
        return;
    }
    const Scene& t = *b->knit->scene;
    const ScenePart* smallest = nullptr;
    for (const ScenePart& r : t.parts)
        if (r.knit && (!smallest || r.mesh->area() < smallest->mesh->area())) smallest = &r;
    CHECK(smallest);
    const std::optional<RenderedPart> rendered = renderedPart(*smallest, s, &err);
    CHECK(rendered.has_value());
    if (!rendered) return;
    const RenderedPart& rp = *rendered;
    CHECK(rp.plies && rp.strands->yarns.size() == 3 * smallest->knit->t.curves.yarns.size());
    const double plyR = rp.strandPp.yarnDiameter / 2;
    const DisplayMesh* mesh = nullptr;
    for (const DisplayMesh& d : b->meshes)
        if (d.kind == DisplayMesh::Kind::Yarn && d.name == smallest->name) mesh = &d;
    CHECK(mesh && !mesh->pos.empty());
    const size_t sides = static_cast<size_t>(rp.sides), stride = sides + 1, vertices = mesh->pos.size() / 3;
    CHECK(vertices % stride == 0);
    auto vertex = [mesh](size_t v) { return Vec3{mesh->pos[3 * v], mesh->pos[3 * v + 1], mesh->pos[3 * v + 2]}; };
    double lo = 1e300, hi = 0;
    for (size_t first = 0; first + stride <= vertices; first += stride) {
        Vec3 centre{0, 0, 0};
        for (size_t i = 0; i < sides; ++i) centre = centre + vertex(first + i);
        centre = centre * (1.0 / double(sides));
        for (size_t i = 0; i < sides; ++i) {
            const double r = norm(vertex(first + i) - centre);
            lo = std::min(lo, r);
            hi = std::max(hi, r);
        }
    }
    std::printf("    %s: %zu ply tubes, %zu rings; ring radius %.4f..%.4f mm, ply radius %.4f mm (yarn radius %.3f mm)\n",
                smallest->name.c_str(), rp.strands->yarns.size(), vertices / stride, lo, hi, plyR, smallest->yarnD / 2);
    CHECK(lo > 0.98 * plyR && hi < 1.02 * plyR);
}

// Stitch size from the yarn (strand.stitch_size = yarn, applyStitchSize): w and h follow D, c, the
// pattern and the tightness, never compress the yarn, and leave the stage keys as the same w and h
// set by hand.
TEST(stitch_size_from_yarn) {
    const SceneEntry* teddy = findScene("teddy");
    const SceneEntry charts = uvChartScene();
    const SceneEntry* man = &charts;
    CHECK(teddy && man);
    auto key = [](const SceneSettings& s) { return stageKey(s, Stage::Knit, 0); };
    // Defaults: custom, nothing derived (the labels are no inputs).
    SceneSettings d;
    CHECK(d.stitchSize == "custom" && d.tightness == "standard");
    SceneSettings d2 = d;
    d2.tightness = "loose";
    CHECK(key(d2) == key(d));
    // The numbers: Standard at 1.0 mm stockinette = the Medium preset (6.5 x 4.6), Loose at 0.8 mm
    // = the Thin preset (7.5 x 5.5), Tight = the fit limit 4.35 x 2.68 (D + c) rounded up.
    std::string err;
    SceneSettings s;
    CHECK(setSetting(s, "strand.stitch_size", "yarn", &err) && setSetting(s, "strand.diameter", "1.0", &err));
    CHECK(std::fabs(s.wale - 6.5) < 1e-9 && std::fabs(s.course - 4.6) < 1e-9);
    CHECK(setSetting(s, "strand.tightness", "loose", &err) && setSetting(s, "strand.diameter", "0.8", &err));
    CHECK(std::fabs(s.wale - 7.5) < 1e-9 && std::fabs(s.course - 5.5) < 1e-9);
    CHECK(setSetting(s, "strand.tightness", "tight", &err) && setSetting(s, "strand.diameter", "1.16", &err));
    CHECK(std::fabs(s.wale - 5.5) < 1e-9 && std::fabs(s.course - 3.4) < 1e-9);  // 4.35 x 1.26 = 5.481, 2.68 x 1.26 = 3.377
    // The pattern leaves the stockinette stitch as it is (garter knits its relaxed size of it, patternStitch).
    CHECK(setSetting(s, "strand.pattern", "garter", &err));
    CHECK(std::fabs(s.wale - 5.5) < 1e-9 && std::fabs(s.course - 3.4) < 1e-9);
    // An explicit w in yarn mode is replaced by the derived one (setSetting re-derives).
    const double w0 = s.wale;
    CHECK(setSetting(s, "strand.wale", "9.0", &err) && s.wale == w0);
    // The keys: yarn mode = custom with the same w and h.
    SceneSettings c = s;
    c.stitchSize = "custom";
    CHECK(key(c) == key(s));
    // Never compressed, for every pattern, tightness and diameter of the UI's range (0.5 to 1.6 mm).
    int n = 0;
    for (const char* pat : {"stockinette", "reverse-stockinette", "garter", "rib1x1", "rib2x2"})
        for (const char* t : {"tight", "standard", "loose"})
            for (double D = 0.5; D <= 1.6 + 1e-9; D += 0.01) {
                SceneSettings y;
                y.stitchSize = "yarn";
                y.pattern = pat;
                y.tightness = t;
                y.diameter = D;
                applyStitchSize(y);
                const YarnFit f = yarnFit(y);
                CHECK(!f.compressed() && !f.tooTight && f.need >= 1 - 1e-9);
                ++n;
            }
    CHECK(n == 5 * 3 * 111);
    // A preset is a stitch size of its own: back to custom (the teddy's Dense).
    SceneSettings t;
    t.stitchSize = "yarn";
    applyStitchSize(t);
    applyThicknessPreset(t, thicknessPresets(*teddy)[0]);
    CHECK(t.stitchSize == "custom" && t.wale == 5.5 && t.course == 3.4 && t.diameter == 1.86);
    // Settings files: the mode and the derived values survive a write / read.
    SceneSettings back;
    CHECK(readSettings(writeSettings(s), back, &err) && back == s);
}

// The Knit step's one Gauge list (gaugeChoices): every entry, once applied, is the entry the
// settings show (currentGaugeChoice), from any starting entry; Custom keeps the values; the yarn
// entries never compress the yarn, the compressed presets do.
TEST(gauge_choices_round_trip) {
    const SceneEntry* teddy = findScene("teddy");
    CHECK(teddy);
    const std::vector<GaugeChoice> list = gaugeChoices(*teddy);
    std::vector<std::string> ids;
    for (const GaugeChoice& c : list) ids.push_back(c.id);
    CHECK((ids == std::vector<std::string>{"thick-dense", "medium", "thin", "yarn-tight", "yarn-standard", "yarn-loose", "custom"}));
    CHECK(!list[0].group.empty() && !list[3].group.empty() && !list[6].group.empty() && list[1].group.empty());
    int checked = 0;
    for (const std::string& from : ids)
        for (const std::string& to : ids) {
            SceneSettings s;
            applySceneDefaults(s, *teddy);
            applyGaugeChoice(s, *teddy, from);
            const SceneSettings before = s;
            applyGaugeChoice(s, *teddy, to);
            CHECK(currentGaugeChoice(s, *teddy) == to);
            if (to == "custom") CHECK(s.diameter == before.diameter && s.wale == before.wale && s.course == before.course);
            if (to.rfind("yarn-", 0) == 0) CHECK(!yarnFit(s).compressed() && yarnFit(s).need >= 1 - 1e-9);
            ++checked;
        }
    CHECK(checked == 49);
    // Dense and Medium flatten their yarn (k < 1); Thin is round.
    for (const char* id : {"thick-dense", "medium", "thin"}) {
        SceneSettings s;
        applyGaugeChoice(s, *teddy, id);
        CHECK(yarnFit(s).compressed() == (std::string(id) != "thin"));
    }
}

// Diagnostic (not in ctest): flat-swatch verdicts over candidate preset parameters; prints the
// grid and the time per reference. Used to choose the Medium / Thin presets (README).
// The Fur baseline's measured shape (fur-loops M2; not registered, run by hand from the main checkout, which holds
// scenes/): the brown teddy at the approved look with the export's 14 M segment budget. Per layer: the total turning
// of each fibre's visible polyline (the sum of its joint angles), the loops' segments per full turn, and the share of fibres whose centreline comes closer to the body than its own radius
// (sampled at the vertices and segment midpoints from the second visible segment on: the buried root and the
// emergence segment are left out; between samples up to half a segment goes unmeasured). A measurement, not a
// rejection: the stuck test catches only a centreline inside or crossing the body.
TEST(diag_fur_baseline_measures) {
    const SceneEntry* e = findScene("teddy-brown");
    CHECK(e != nullptr);
    SceneSettings s;
    s.scene = "teddy-brown";
    applySceneDefaults(s, *e);
    s.workflow = "fur";
    s.segmentBudget = 14000000;
    Scene t;
    std::string err;
    if (!loadBaseParts(*e, s, SourceOptions{}, t, &err)) SKIP(err.c_str());
    SceneFibres grownFur;
    CHECK(grown(grownFur, growSceneFur(*e, t, s, &err)));
    std::printf("%s", grownFur.log.c_str());
    TriMesh body;
    for (const SourcePrimitive& g : t.source.prims) {
        if (g.cls == SourcePrimitive::Class::Skip) continue;
        const int base = static_cast<int>(body.pos.size());
        body.pos.insert(body.pos.end(), g.welded.pos.begin(), g.welded.pos.end());
        for (const auto& tr : g.welded.tri) body.tri.push_back({tr[0] + base, tr[1] + base, tr[2] + base});
    }
    body.computeNormals();
    const MeshQuery query(body);
    struct Layer {
        long long fibres = 0, close = 0;
        std::vector<double> turning, perTurn;
    } layers[2];
    for (const PartFibres& pf : grownFur.parts) {
        const Vec3 centre = t.parts[static_cast<size_t>(pf.part)].centre;
        const std::vector<Fibre>& fibres = pf.fs.fibres;
        std::vector<char> close(fibres.size(), 0);
        std::vector<double> turning(fibres.size(), 0);
        parallelFor(static_cast<int64_t>(fibres.size()), 512, [&](int64_t b, int64_t end, int64_t) {
            for (int64_t q = b; q < end; ++q) {
                const Fibre& f = fibres[static_cast<size_t>(q)];
                const std::span<const Vec3f> pts = pf.fs.pointsOf(f);
                std::vector<Vec3> p;
                for (const Vec3f& v : pts) p.push_back(Vec3{v.x, v.y, v.z} + centre);
                const int nv = static_cast<int>(p.size());
                double turn = 0;  // the joints of the visible polyline
                for (int i = 3; i < nv; ++i) turn += std::acos(std::clamp(dot(normalized(p[i - 1] - p[i - 2]), normalized(p[i] - p[i - 1])), -1.0, 1.0));
                turning[static_cast<size_t>(q)] = turn;
                for (int i = 2; i < nv && !close[static_cast<size_t>(q)]; ++i) {
                    const double r = taperedRadius(f.radius0, f.radius1, i, nv - 1);
                    // The vertex, and from the second visible segment on its midpoint (the emergence segment leaves the skin).
                    const Vec3 samples[2] = {p[i], (p[i] + p[i - 1]) * 0.5};
                    for (int k = 0; k < (i >= 3 ? 2 : 1); ++k) {
                        const MeshQuery::Hit h = query.closest(samples[k], r + 1e-9);
                        if (h.tri >= 0 && (h.inside || h.dist < r)) close[static_cast<size_t>(q)] = 1;
                    }
                }
            }
        }, ExecPolicy{true});
        for (size_t q = 0; q < fibres.size(); ++q) {
            Layer& l = layers[fibres[q].radius0 < 0.1f ? 1 : 0];
            ++l.fibres;
            l.close += close[q];
            l.turning.push_back(turning[q] * 180 / kPi);
            const int segments = static_cast<int>(pf.fs.pointsOf(fibres[q]).size()) - 2;
            if (fibres[q].radius0 >= 0.1f && turning[q] > 0) l.perTurn.push_back(segments / (turning[q] / (2 * kPi)));
        }
    }
    auto percentile = [](std::vector<double> v, double share) {
        if (v.empty()) return 0.0;
        std::sort(v.begin(), v.end());
        return v[static_cast<size_t>(share * double(v.size() - 1))];
    };
    for (int li = 0; li < 2; ++li) {
        const Layer& l = layers[li];
        std::printf("%s: %lld fibres; turning [deg] p5 %.0f, p50 %.0f, p95 %.0f; centreline closer to the body than its radius: %lld (%.3f %%)\n",
                    li == 0 ? "loops" : "fine fibres", l.fibres, percentile(l.turning, 0.05), percentile(l.turning, 0.5), percentile(l.turning, 0.95),
                    l.close, 100.0 * double(l.close) / double(std::max(1LL, l.fibres)));
    }
    // The loops' segments per full turn, in bins of one.
    std::map<int, long long> bins;
    for (double v : layers[0].perTurn) ++bins[std::min(12, static_cast<int>(v))];
    std::printf("loops, segments per full turn (polyline):");
    for (const auto& [b, n] : bins) std::printf(" [%d%s): %lld", b, b == 12 ? "+" : "", n);
    std::printf("\nloops turning more than 320 deg: %.1f %%\n",
                100.0 * double(std::count_if(layers[0].turning.begin(), layers[0].turning.end(), [](double d) { return d > 320; })) /
                    double(std::max<size_t>(1, layers[0].turning.size())));
    CHECK(layers[0].fibres > 0 && layers[1].fibres > 0);
}

// The disc preview's latency on the brown teddy at the approved look (fur-loops M4; not registered, run by hand
// from the main checkout): request to a published preview, end to end through the worker thread, and the parts of
// it (scene, roots, fields, buffers). Single edits are built in full; a drag (a request every 16 ms) is drawn as
// drafts (a quarter of the fibres) and in full once it stops. Targets (handover section 8): disc rebuild <= 100 ms,
// root cache rebuild <= 5 s.
TEST(diag_fur_preview_latency) {
    using Clock = std::chrono::steady_clock;
    const SceneEntry* e = findScene("teddy-brown");
    CHECK(e != nullptr);
    SceneSettings s;
    s.scene = "teddy-brown";
    applySceneDefaults(s, *e);
    s.workflow = "fur";
    FurPreviewWorker worker;
    const auto ms = [](Clock::time_point a, Clock::time_point b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
    // A request, then the preview of it in full (keyBelow >= 1).
    auto full = [&](const SceneSettings& t, const char* what) {
        const Clock::time_point t0 = Clock::now();
        const uint64_t g = worker.request(t, FurDiscColour::Plain);
        std::shared_ptr<const FurPreview> p;
        while (!(p = worker.latest()) || p->generation < g || p->keyBelow < 1) {
            if (!worker.error().empty() && !worker.busy()) break;
            std::this_thread::yield();  // not sleep_for(1 ms): Windows rounds it up to its 15.6 ms timer
        }
        const double total = ms(t0, Clock::now());
        if (p)
            std::printf("    %-34s %7.0f ms end to end: scene %.0f, roots %.0f, fields %.0f, buffers %.0f ms; %lld discs, %lld needles\n", what, total,
                        p->msScene, p->msRoots, p->msField, p->msBuffers, p->loops, p->fibres);
        else std::printf("    %s: %s\n", what, worker.error().c_str());
        std::this_thread::sleep_for(std::chrono::milliseconds(400));  // the next request is a single edit, not a drag
        return total;
    };
    const auto percentile = [](std::vector<double> v, double share) {
        std::sort(v.begin(), v.end());
        return v.empty() ? 0.0 : v[static_cast<size_t>(share * double(v.size() - 1))];
    };
    full(s, "first (scene and roots)");
    std::vector<double> single;
    for (int k = 0; k < 12; ++k) {
        SceneSettings t = s;
        t.pile.curl = 0.9 + 0.02 * (k % 6 + 1);
        single.push_back(full(t, "curl edit"));
    }
    std::printf("    single curl edits in full: p50 %.0f ms, p95 %.0f ms (n %zu)\n", percentile(single, 0.5), percentile(single, 0.95), single.size());
    // A drag: 60 requests 16 ms apart; every preview published on the way is timed from its own request.
    std::map<uint64_t, Clock::time_point> asked;
    std::vector<double> drafts;
    uint64_t seen = 0;
    const auto watch = [&] {
        const std::shared_ptr<const FurPreview> p = worker.latest();
        if (!p || p->generation == seen || !asked.count(p->generation)) return p;
        seen = p->generation;
        if (p->keyBelow < 1) drafts.push_back(ms(asked[p->generation], Clock::now()));
        return p;
    };
    Clock::time_point last;
    for (int k = 0; k < 60; ++k) {
        SceneSettings t = s;
        t.pile.curl = 0.8 + 0.005 * k;
        last = Clock::now();
        asked[worker.request(t, FurDiscColour::Plain)] = last;
        for (const Clock::time_point until = last + std::chrono::milliseconds(16); Clock::now() < until;) {
            watch();
            std::this_thread::yield();
        }
    }
    std::shared_ptr<const FurPreview> p;
    while (!(p = worker.latest()) || p->generation < asked.rbegin()->first || p->keyBelow < 1) {
        watch();
        std::this_thread::yield();
    }
    std::printf("    drag of 60 requests (16 ms apart): %zu drafts shown, request to draft p50 %.0f ms, p95 %.0f ms; last request to "
                "the full preview %.0f ms\n",
                drafts.size(), percentile(drafts, 0.5), percentile(drafts, 0.95), ms(last, Clock::now()));
    SceneSettings lower = s;
    lower.pile.density = 300;
    full(lower, "density 300 (within the capacity)");
    SceneSettings within = s;
    within.pile.density = 480;
    full(within, "density 480 (within 1.25 x 400)");
    SceneSettings beyond = s;
    beyond.pile.density = 600;
    full(beyond, "density 600 (resampled)");
    p = worker.latest();
    CHECK(p != nullptr);
    if (!p) return;
    const double mb = double(p->discPos.size() + p->discNormal.size() + p->discLocal.size() + p->discColour.size() + p->markPos.size() +
                             p->markKey.size() + p->markColour.size() + p->needlePos.size() + p->needleKey.size() + p->needleColour.size()) *
                      4 / 1048576.0;
    std::printf("    buffers at density 600: %.1f MB; sizeof FurRoot %zu B, FurDisc %zu B\n", mb, sizeof(FurRoot), sizeof(FurDisc));
    // M5: every irregularity on (patches, clumps, islands, wear): their fields are built once, then edits reuse them.
    SceneSettings irregular = s;
    irregular.pile.varDensity = 0.6;
    irregular.pile.varDensityScale = 12;
    irregular.pile.patchLength = irregular.pile.patchLift = 0.3;
    irregular.pile.clumpSize = 2.5;
    irregular.pile.clumpStrength = 0.7;
    irregular.pile.islandsRate = 0.15;
    irregular.pile.islandsSize = 20;
    irregular.pile.islandsStrength = 0.7;
    irregular.pile.wear = 0.8;
    full(irregular, "irregularities on (fields built)");
    std::vector<double> irregularEdits;
    for (int k = 0; k < 8; ++k) {
        SceneSettings t = irregular;
        t.pile.curl = 0.9 + 0.02 * (k % 4 + 1);
        irregularEdits.push_back(full(t, "curl edit, irregularities on"));
    }
    std::printf("    curl edits with the irregularities on, in full: p50 %.0f ms, p95 %.0f ms (n %zu)\n", percentile(irregularEdits, 0.5),
                percentile(irregularEdits, 0.95), irregularEdits.size());
}

TEST(diag_preset_sweep) {
    const char* names[] = {"st", "rev", "garter", "rib1", "rib2"};
    struct C { double D, w, h; };
    std::vector<C> cs;
    for (double w : {6.0, 6.5, 7.0})
        for (double h : {4.2, 4.6, 5.0}) cs.push_back({1.0, w, h});
    for (double D : {0.8, 0.9})
        for (double w : {7.0, 7.5, 8.0})
            for (double h : {5.0, 5.5, 6.0}) cs.push_back({D, w, h});
    for (const C& c : cs) {
        PatternParams pp;
        pp.yarnDiameter = c.D;
        pp.waleSpacing = c.w;
        pp.courseSpacing = c.h;
        const auto t0 = std::chrono::steady_clock::now();
        std::string line;
        bool all = true;
        for (int k = 0; k < 5; ++k) {
            const PatternReference r = flatPatternReference(pp, static_cast<StitchPattern>(k));
            line += std::string(" ") + names[k] + "=" + toString(r.verdict);
            all &= r.verdict == Status::Pass;
        }
        const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        std::printf("    D %.2f w %.2f h %.2f:%s%s (%.0f ms)\n", c.D, c.w, c.h, line.c_str(), all ? "  ALL PASS" : "", ms);
    }
}

// ------------------------------------------------------- unified workflow ----

// The workflow key and the former selectors (scene_settings.h readSettings): knitted files read as workflow = knit;
// the Hair mode's files are refused (owner, 2026-11-19) and the settings they were read into stay as they were;
// contradictions, repeated selectors and unknown workflows are refused whatever the line order; workflow = fur
// round-trips; no preset carries a workflow; the builder refuses a knit job under Fur.
TEST(workflow_key_and_migration) {
    std::string err;
    SceneSettings a;
    a.scene = "teddy";
    SceneSettings b = a;
    b.workflow = "fur";
    SceneSettings r;
    CHECK(readSettings(writeSettings(b), r, &err) && r == b && furWorkflow(r));
    const std::string text = writeSettings(a);
    CHECK(text.find("\nworkflow = knit\n") != std::string::npos && text.find("knit.enabled") == std::string::npos &&
          text.find("fibres.profile") == std::string::npos);
    // Knit stages do not see the workflow; the fur result key does.
    for (Stage st : {Stage::Load, Stage::Knit, Stage::Surface}) CHECK(stageKey(b, st, 7) == stageKey(a, st, 7));
    CHECK(stageInputs(a, Stage::Knit).find("knit.enabled=true\n") != std::string::npos);
    CHECK(fibreResultKey(b) != fibreResultKey(a));
    // The migration table, each in both line orders where two selectors meet.
    struct Row {
        const char* lines;
        const char* workflow;   // "" refused
        const char* mentions;   // the refusal names it
    };
    const Row rows[] = {
        {"fibres.profile = wool\n", "knit", ""},
        {"knit.enabled = true\n", "knit", ""},
        {"knit.enabled = true\nfibres.profile = wool\n", "knit", ""},
        {"", "knit", ""},
        {"fibres.profile = fur\n", "", "Hair mode"},
        {"fibres.profile = fur\nknit.enabled = false\n", "", "Hair mode"},
        {"knit.enabled = false\n", "", "fibres.profile"},
        {"knit.enabled = false\nfibres.profile = wool\n", "", "fibres.profile = wool"},
        {"fibres.profile = wool\nknit.enabled = false\n", "", "fibres.profile = wool"},
        {"knit.enabled = true\nfibres.profile = fur\n", "", "knit.enabled = true"},
        {"fibres.profile = fur\nknit.enabled = true\n", "", "knit.enabled = true"},
        {"workflow = fur\nknit.enabled = false\n", "fur", ""},
        {"knit.enabled = false\nworkflow = fur\n", "fur", ""},
        {"workflow = fur\nfibres.profile = fur\n", "fur", ""},
        {"workflow = knit\nknit.enabled = false\n", "", "contradicts"},
        {"knit.enabled = false\nworkflow = knit\n", "", "contradicts"},
        {"fibres.profile = wool\nworkflow = fur\n", "", "contradicts"},
        {"workflow = knit\nworkflow = knit\n", "", "twice"},
        {"knit.enabled = true\nknit.enabled = true\n", "", "twice"},
        {"workflow = hair\n", "", "workflow"},
        {"knit.enabled = maybe\n", "", "knit.enabled"},
    };
    for (const Row& row : rows) {
        SceneSettings s;
        s.workflow = "fur";  // a refused file leaves the settings it was read into as they were
        s.pile.density = 123;
        const SceneSettings before = s;
        const bool read = readSettings(std::string("knit-settings v1\nscene = teddy\n") + row.lines, s, &err);
        std::printf("    %-46s -> %s\n", std::string(row.lines).substr(0, 44).c_str(), read ? s.workflow.c_str() : err.c_str());
        CHECK(read == (*row.workflow != 0));
        if (read) CHECK(s.workflow == row.workflow);
        else CHECK(s == before && err.find(row.mentions) != std::string::npos);
    }
    // The Hair mode's saved file as main's last build wrote it.
    {
        std::ifstream in(std::string(KNIT_TEST_DATA) + "/legacy_settings/teddy_hair.settings");
        const std::string hair((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        SceneSettings s;
        CHECK(!hair.empty() && !readSettings(hair, s, &err) && err.find("Hair mode") != std::string::npos);
    }
    // No preset carries a workflow: every one opens with the settings' default, knit (the viewer's draft asks first).
    for (const SceneEntry& e : sceneRegistry()) {
        const SceneSettings v = viewerDefaults(e);
        std::printf("    %-17s workflow = %s\n", e.id.c_str(), v.workflow.c_str());
        CHECK(knitWorkflow(v));
    }
    // A knit job under Fur is refused.
    SceneSettings fur;
    fur.scene = "teddy";
    fur.workflow = "fur";
    CHECK(!run(BuiltScene::Kind::Knit, fur, &err) && err.find("no knit") != std::string::npos);
}

namespace {
// knit_check (the command-line tool): KNIT_TEST_KNIT_CHECK overrides the built one (to run these
// checks against another build).
std::string knitCheckExe() {
    char* env = nullptr;
    size_t envLen = 0;
    if (_dupenv_s(&env, &envLen, "KNIT_TEST_KNIT_CHECK") == 0 && env) {
        const std::string p = env;
        std::free(env);
        return p;
    }
#ifdef KNIT_CHECK_EXE
    return KNIT_CHECK_EXE;
#else
    return "";
#endif
}
int runLogged(const std::string& args, const std::string& log) {
    const std::string cmd = "\"\"" + knitCheckExe() + "\" " + args + " > \"" + log + "\" 2>&1\"";
    return std::system(cmd.c_str());
}
} // namespace

// knit_check under Fur: the package holds the base meshes (no more vertices than the source), the pile's strands
// and nothing knitted: no ply maps, no yarn or wool fibre material, a finite camera (the leak of 2026-10-25, when a
// hair-only scene was knitted under its hair, stays checked). The former selectors and a pile value out of bounds
// are refused before anything runs (exit 2).
TEST(cli_fur_export_has_no_knit) {
    if (knitCheckExe().empty() || !fs::exists(knitCheckExe())) SKIP("knit_check not configured");
    const std::string d = freshDir("cli_fur"), root = freshDir("cli_fur_cache");
    ImportCacheEnv env(root);
    SceneSettings s = importPlane(d, root);
    s.workflow = "fur";
    s.pile.density = 100;
    s.pile.fibresDensity = 50;
    const SceneEntry* e = sceneEntryFor(s);
    CHECK(e != nullptr);
    if (!e) return;
    {
        std::ofstream f(d + "/fur.settings");
        f << writeSettings(s);
    }
    const std::string out = d + "/plane_fur", log = d + "/export.log";
    const int rc = runLogged("--scene import --settings \"" + d + "/fur.settings\" --export-ppt \"" + out + "\"", log);
    const std::string text = readText(log);
    if (rc != 0) std::printf("%s\n", text.c_str());
    CHECK(rc == 0);
    CHECK(fs::exists(out + "/scene.json"));
    // Nothing was knitted.
    for (const char* knitLine : {"gate a", "stitches,", "texture profile"})
        if (text.find(knitLine) != std::string::npos) {
            std::printf("    the log knits: found \"%s\"\n", knitLine);
            CHECK(!"knit_check knitted a Fur scene");
        }
    // The base meshes: each file is one part with exactly its source triangles (the texture's uv seams may split
    // vertices; yarn tubes would have far more faces).
    Scene t;
    std::string err;
    CHECK(loadBaseParts(*e, s, SourceOptions{}, t, &err));
    std::map<std::string, size_t> partTris;
    for (const ScenePart& r : t.parts) partTris["meshes/" + r.name + ".ply"] = r.mesh->tri.size();
    size_t meshes = 0, strands = 0;
    for (const std::string& f : filesUnder(out)) {
        if (f.rfind("textures/ply_", 0) == 0) {
            std::printf("    ply map %s\n", f.c_str());
            CHECK(!"ply maps in a Fur export");
        }
        if (f.find("_stuffing") != std::string::npos || f.find("_backing") != std::string::npos) {
            std::printf("    filling %s\n", f.c_str());
            CHECK(!"filling in a Fur export");
        }
        if (f.rfind("strands/", 0) == 0) ++strands;
        if (f.rfind("meshes/", 0) != 0 || f.find("_ground") != std::string::npos) continue;
        ++meshes;
        std::ifstream pf(out + "/" + f, std::ios::binary);
        size_t verts = 0, faces = 0;
        for (std::string line; std::getline(pf, line) && line != "end_header";) {
            if (line.rfind("element vertex ", 0) == 0) verts = std::stoul(line.substr(15));
            if (line.rfind("element face ", 0) == 0) faces = std::stoul(line.substr(13));
        }
        const auto it = partTris.find(f);
        std::printf("    %s: %zu vertices, %zu faces (source part: %zu triangles)\n", f.c_str(), verts, faces,
                    it == partTris.end() ? size_t(0) : it->second);
        CHECK(it != partTris.end() && verts > 0 && faces == it->second);
    }
    CHECK(meshes == t.parts.size() && strands >= 1);
    const std::string json = readText(out + "/scene.json");
    for (const char* mat : {"\"wool_white\"", "\"wool_brown\"", "\"stuffing\"", "\"fibre_directed\"", "\"yarn_"})
        if (json.find(mat) != std::string::npos) {
            std::printf("    scene.json lists %s\n", mat);
            CHECK(!"a knitting material in a Fur export");
        }
    // A finite camera: no "inf" / "nan" number token (as a word: "provenance" contains "nan").
    auto numberToken = [&](const std::string& w) {
        for (size_t p = json.find(w); p != std::string::npos; p = json.find(w, p + 1)) {
            const bool before = p == 0 || !std::isalpha(static_cast<unsigned char>(json[p - 1]));
            const bool after = p + w.size() >= json.size() || !std::isalpha(static_cast<unsigned char>(json[p + w.size()]));
            if (before && after) return true;
        }
        return false;
    };
    CHECK(!numberToken("inf") && !numberToken("nan") && json.find("\"cameras\"") != std::string::npos);
    // Refused before anything runs.
    CHECK(runLogged("--scene teddy --set fibres.profile=fur --export-ppt \"" + d + "/refused1\"", d + "/r1.log") != 0);
    CHECK(readText(d + "/r1.log").find("replaced by workflow") != std::string::npos);
    CHECK(runLogged("--scene teddy --set workflow=fur --set pile.curl=0.2 --export-ppt \"" + d + "/refused2\"", d + "/r2.log") != 0);
    CHECK(readText(d + "/r2.log").find("pile.curl") != std::string::npos);
    CHECK(!fs::exists(d + "/refused1") && !fs::exists(d + "/refused2"));
}

// knit_check --save-knit / --load-knit on a scene knitted on its mesh (an import): the knit is kept
// as one file with its settings key, and a load exports that knit without knitting again (the
// solver is not deterministic: a frozen knit is what later comparisons need). Other knitting
// settings are refused.
TEST(cli_tiled_knit_save_load) {
    if (knitCheckExe().empty() || !fs::exists(knitCheckExe())) SKIP("knit_check not configured");
    const std::string d = freshDir("cli_tiled_knit"), root = freshDir("cli_tiled_knit_cache");
    ImportCacheEnv env(root);
    SceneSettings s = importPlane(d, root);
    s.fibres = false;
    {
        std::ofstream f(d + "/plane.settings");
        f << writeSettings(s);
    }
    const std::string settings = "--scene import --settings \"" + d + "/plane.settings\"";
    const std::string knitDir = d + "/frozen";
    const int saved = runLogged(settings + " --save-knit \"" + knitDir + "\" --export-ppt \"" + d + "/first\"", d + "/save.log");
    if (saved != 0) std::printf("%s\n", readText(d + "/save.log").c_str());
    CHECK(saved == 0 && fs::exists(knitDir + "/knit.bin") && fs::exists(knitDir + "/knit.key"));
    const int loaded = runLogged(settings + " --load-knit \"" + knitDir + "\" --export-ppt \"" + d + "/second\"", d + "/load.log");
    const std::string loadLog = readText(d + "/load.log");
    if (loaded != 0) std::printf("%s\n", loadLog.c_str());
    CHECK(loaded == 0 && loadLog.find("knit: loaded from") != std::string::npos);
    // The same knit: every yarn tube file is byte-identical.
    size_t tubes = 0;
    for (const std::string& f : filesUnder(d + "/first")) {
        if (f.rfind("meshes/", 0) != 0 || f.find("_ground") != std::string::npos) continue;
        ++tubes;
        CHECK(readText(d + "/first/" + f) == readText(d + "/second/" + f));
    }
    CHECK(tubes >= 1);
    // Other knitting settings: refused, nothing exported.
    const int other = runLogged(settings + " --set strand.course=" + std::to_string(s.course * 1.1) + " --load-knit \"" + knitDir +
                                "\" --export-ppt \"" + d + "/third\"", d + "/other.log");
    CHECK(other == 2 && readText(d + "/other.log").find("other knitting settings") != std::string::npos && !fs::exists(d + "/third"));
}

// The yarn body (fibres/yarn_body.h, a prototype) in a ppt export: each knitted part's yarn as structural
// fibres in place of its tube, in a hair material of the yarn's colour with no wool tip taper, without ply
// maps; the core proxy as a mesh when asked; explicit ply tubes refused. The renderer loads and renders it.
TEST(cli_yarn_body_export) {
    if (knitCheckExe().empty() || !fs::exists(knitCheckExe())) SKIP("knit_check not configured");
    const std::string d = freshDir("cli_yarn_body"), root = freshDir("cli_yarn_body_cache");
    ImportCacheEnv env(root);
    SceneSettings s = importPlane(d, root);
    s.fibres = false;
    s.textureProfile = "2-swept";
    s.yarnBody.enabled = true;
    s.yarnBody.plies = 2;
    s.yarnBody.fibres = 40;
    s.yarnBody.core = "opaque";
    s.filling = false;  // no backing mesh: the knitted part's only mesh is its core
    {
        std::ofstream f(d + "/body.settings");
        f << writeSettings(s);
    }
    const std::string out = d + "/body", log = d + "/body.log";
    const int rc = runLogged("--scene import --settings \"" + d + "/body.settings\" --export-ppt \"" + out + "\"", log);
    const std::string text = readText(log);
    if (rc != 0) std::printf("%s\n", text.c_str());
    CHECK(rc == 0 && text.find("yarn body (2 plies, opaque core)") != std::string::npos);
    size_t bodies = 0, cores = 0;
    long long bodySegments = 0;
    for (const std::string& f : filesUnder(out)) {
        CHECK(f.rfind("textures/ply_", 0) != 0);  // no ply maps
        if (f.rfind("strands/", 0) == 0 && f.find("_body") != std::string::npos) {
            ++bodies;
            StrandSetF back;
            std::string err;
            CHECK(readStrandsFile(out + "/" + f, back, &err) && back.strands() > 0);
            bodySegments += static_cast<long long>(back.vertices() - back.strands());
            StrandRules rules;
            rules.minSegOverRMax = 4;
            const StrandCheck check = checkStrands(back, rules);  // in metres, as written
            if (!check.ok) std::printf("    %s\n", check.error.c_str());
            CHECK(check.ok);
        }
        if (f.rfind("meshes/", 0) == 0 && f.find("_core") != std::string::npos) ++cores;
        // No yarn tube: a knitted part's mesh is its core or the ground only.
        if (f.rfind("meshes/", 0) == 0 && f.find("_core") == std::string::npos && f.find("_ground") == std::string::npos)
            CHECK(!"a yarn tube in a yarn-body export");
    }
    CHECK(bodies >= 1 && cores == bodies);
    const std::string json = readText(out + "/scene.json");
    CHECK(json.find("\"fibre_body_") != std::string::npos && json.find("\"type\": \"hair\"") != std::string::npos);
    CHECK(json.find("part_fibre_shapes") == std::string::npos);  // the body is no fuzz: no wool tip taper
    CHECK(json.find("\"core_") != std::string::npos);
    // Split into strand files of at most 20000 segments (knit_check --body-file-segments; the default is
    // 2.5 M): whole strands, the same segments in all, every file in the object's yarn_body group; the
    // OBJ export writes one object per file.
    const long long kFileSegments = 20000;
    const std::string split = d + "/split", obj = d + "/obj";
    const std::string splitArgs = "--scene import --settings \"" + d + "/body.settings\" --body-file-segments " + std::to_string(kFileSegments);
    CHECK(runLogged(splitArgs + " --export-ppt \"" + split + "\"", d + "/split.log") == 0);
    size_t files = 0;
    long long splitSegments = 0;
    for (const std::string& f : filesUnder(split)) {
        if (f.rfind("strands/", 0) != 0 || f.find("_body") == std::string::npos) continue;
        ++files;
        CHECK(std::regex_search(f, std::regex("_body_[0-9]+\\.pptstrands$")));
        StrandSetF back;
        std::string err;
        CHECK(readStrandsFile(split + "/" + f, back, &err));
        const long long own = static_cast<long long>(back.vertices() - back.strands());
        CHECK(own <= kFileSegments || back.strands() == 1);
        splitSegments += own;
    }
    std::printf("    body: %lld segments in %zu file(s); split at %lld: %zu files\n", bodySegments, bodies, kFileSegments, files);
    CHECK(files > bodies && splitSegments == bodySegments);
    const std::string splitJson = readText(split + "/scene.json");
    std::set<std::string> bodyGroups;
    const std::regex groupRow("\"id\": \"(g[0-9]+)\", \"name\": \"yarn_body\"");
    for (std::sregex_iterator it(splitJson.begin(), splitJson.end(), groupRow), end; it != end; ++it) bodyGroups.insert((*it)[1]);
    size_t placed = 0;
    const std::regex bodyInstance("\"name\": \"[^\"]*_body_[0-9]+\", \"group_id\": \"(g[0-9]+)\"");
    for (std::sregex_iterator it(splitJson.begin(), splitJson.end(), bodyInstance), end; it != end; ++it) {
        ++placed;
        CHECK(bodyGroups.count((*it)[1]) == 1);
    }
    CHECK(bodyGroups.size() == 1 && placed == files);
    CHECK(runLogged(splitArgs + " --export-format obj --export-ppt \"" + obj + "\"", d + "/obj.log") == 0);
    size_t objects = 0;
    const std::regex bodyObject("\no [^\n]*_body_[0-9]+\n");
    for (const std::string& f : filesUnder(obj)) {
        if (f.size() < 4 || f.substr(f.size() - 4) != ".obj") continue;
        const std::string objText = readText(obj + "/" + f);
        for (std::sregex_iterator it(objText.begin(), objText.end(), bodyObject), end; it != end; ++it) ++objects;
    }
    std::printf("    OBJ: %zu body objects\n", objects);
    CHECK(objects == files);
#ifdef KNIT_PPT_EXE
    if (fs::exists(KNIT_PPT_EXE)) {
        const std::string render = d + "/render.log";
        const std::string cmd = std::string("\"\"") + KNIT_PPT_EXE + "\" render \"" + out + "/scene.json\" --output \"" + d +
                                "/out.png\" --width 64 --height 64 --spp 1 > \"" + render + "\" 2>&1\"";
        const int rr = std::system(cmd.c_str());
        if (rr != 0) std::printf("%s\n", readText(render).c_str());
        CHECK(rr == 0 && fs::exists(d + "/out.png"));
    }
#endif
    // Explicit ply tubes with the body: refused before anything is written.
    s.plies = 3;
    {
        std::ofstream f(d + "/refused.settings");
        f << writeSettings(s);
    }
    CHECK(runLogged("--scene import --settings \"" + d + "/refused.settings\" --export-ppt \"" + d + "/refused\"", d + "/refused.log") == 2);
    CHECK(readText(d + "/refused.log").find("needs strand.plies = 1") != std::string::npos);
}

// The yarn of fibres rendered as maps (yarn_body.render = maps; fibres/yarn_bake.h, scene/appearance/baked_maps.h):
// no body strands and no cores, one textured tube mesh per yarn with uv in the yarn body's group, five maps per
// yarn under textures/; the renderer loads the package when it is configured; a second export of the same
// configuration takes the maps from the cache and writes the same bytes. The map material (yarn_body.map_material):
// the default fabric is the yarn's fabric under the cavity and normal maps with the yarn's translucency and, under its
// sheen, the highlight of the roughness and specular maps (the renderer's R1); shiny is a plastic with the same
// highlight block and no translucency, from the same cached maps (a material value never bakes).
TEST(cli_yarn_body_maps_export) {
    if (knitCheckExe().empty() || !fs::exists(knitCheckExe())) SKIP("knit_check not configured");
    const std::string d = freshDir("cli_yarn_body_maps"), root = freshDir("cli_yarn_body_maps_cache");
    ImportCacheEnv env(root);
    SceneSettings s = importPlane(d, root);
    s.fibres = false;
    s.yarnBody.enabled = true;
    s.yarnBody.plies = 2;
    s.yarnBody.fibres = 40;
    s.yarnBody.core = "opaque";  // not an input of the maps: no core is written
    s.yarnBody.render = "maps";
    s.yarnBody.mapTexelUm = 80;  // coarse: a fast test
    s.filling = false;
    s.yarnWhite.translucency = s.yarnBrown.translucency = 0.4;  // a translucent yarn (a wool profile's value)
    CHECK(s.yarnBody.mapMaterial == "fabric");                  // the default
    {
        std::ofstream f(d + "/maps.settings");
        f << writeSettings(s);
        SceneSettings shiny = s;
        shiny.yarnBody.mapMaterial = "shiny";
        std::ofstream g(d + "/shiny.settings");
        g << writeSettings(shiny);
    }
    std::vector<std::string> logs;
    for (const char* run : {"first", "second", "shiny"}) {
        const std::string out = d + "/" + run, log = d + "/" + run + ".log";
        const std::string settings = d + (std::string(run) == "shiny" ? "/shiny.settings" : "/maps.settings");
        const int rc = runLogged("--scene import --settings \"" + settings + "\" --export-ppt \"" + out + "\"", log);
        const std::string text = readText(log);
        if (rc != 0) std::printf("%s\n", text.c_str());
        CHECK(rc == 0 && text.find("yarn maps (2 plies)") != std::string::npos);
        logs.push_back(text);
    }
    // The first export bakes unless an earlier run left the entry in the tool's cache; the second always reads it.
    CHECK(logs[1].find("from the cache") != std::string::npos && logs[2].find("from the cache") != std::string::npos);
    const std::string out = d + "/first";
    const auto bytesOf = [](const std::string& path) {
        std::ifstream f(path, std::ios::binary);
        return std::string(std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>());
    };
    size_t tubes = 0, maps = 0;
    std::string err;
    for (const std::string& f : filesUnder(out)) {
        CHECK(f.find("_body") == std::string::npos && f.find("_core") == std::string::npos);  // the maps replace both
        if (f.rfind("textures/baked_", 0) == 0) {
            ++maps;
            int w = 0, h = 0;
            std::vector<uint8_t> rgb;
            CHECK(readPngRgb8(out + "/" + f, w, h, rgb) && w == 4096 && h > 0 && h <= 4096 && rgb.size() == size_t(3) * w * h);
            CHECK(bytesOf(out + "/" + f) == bytesOf(d + "/second/" + f));  // the cache, byte for byte
        }
        if (f.rfind("meshes/", 0) == 0 && f.find("_baked_y") != std::string::npos) {
            ++tubes;
            const std::optional<StoredMesh> mesh = readPartMesh(out + "/" + f, &err);
            CHECK(mesh && mesh->stride == 8 && !mesh->tri.empty());
            if (mesh) {
                bool inRange = true;
                for (size_t v = 0; v + 7 < mesh->vertices.size(); v += 8)
                    inRange = inRange && mesh->vertices[v + 6] >= 0 && mesh->vertices[v + 6] <= 1 && mesh->vertices[v + 7] >= 0 && mesh->vertices[v + 7] <= 1;
                CHECK(inRange);
            }
        }
    }
    std::printf("    %zu mapped tubes, %zu maps\n", tubes, maps);
    CHECK(tubes >= 1 && maps == 5 * tubes);
    const std::string json = readText(out + "/scene.json");
    // The baked yarns' material entries (one line each in scene.json).
    const auto bakedMaterials = [](const std::string& text) {
        std::vector<std::string> lines;
        std::stringstream in(text);
        std::string line;
        while (std::getline(in, line))
            if (line.find("{\"id\": \"baked_") != std::string::npos) lines.push_back(line);
        return lines;
    };
    const std::vector<std::string> fabric = bakedMaterials(json), shiny = bakedMaterials(readText(d + "/shiny/scene.json"));
    std::printf("    %zu fabric and %zu shiny baked materials\n", fabric.size(), shiny.size());
    CHECK(!fabric.empty() && fabric.size() == shiny.size());
    // The highlight block of a material line: from "specular": { to its closing brace.
    const auto highlightOf = [](const std::string& m) {
        const size_t at = m.find("\"specular\": {");
        return at == std::string::npos ? std::string() : m.substr(at, m.find('}', at) + 1 - at);
    };
    for (size_t i = 0; i < fabric.size(); ++i) {
        const std::string& m = fabric[i];
        const std::string block = highlightOf(m);
        CHECK(m.find("\"type\": \"fabric\"") != std::string::npos && m.find("\"albedo_tex\": \"textures/baked_") != std::string::npos &&
              m.find("\"normal_tex\": \"textures/baked_") != std::string::npos && m.find("\"sheen\": ") != std::string::npos &&
              m.find("\"translucency\": {\"transmit_share\": 0.4}") != std::string::npos);
        CHECK(block.find("\"reflectance_tex\": \"textures/baked_") != std::string::npos && block.find("\"roughness_tex\": \"textures/baked_") != std::string::npos &&
              block.find("\"aniso_rotation\": ") != std::string::npos && block.find("\"roughness_y\": 1.0") != std::string::npos);
        CHECK(i < shiny.size() && block == highlightOf(shiny[i]));  // one highlight, written the same in both materials
    }
    for (const std::string& m : shiny)
        CHECK(m.find("\"type\": \"plastic\"") != std::string::npos && m.find("\"reflectance_tex\": \"textures/baked_") != std::string::npos &&
              m.find("\"roughness_tex\": \"textures/baked_") != std::string::npos && m.find("\"normal_tex\": \"textures/baked_") != std::string::npos &&
              m.find("\"albedo_tex\": \"textures/baked_") != std::string::npos && m.find("\"aniso_rotation\": ") != std::string::npos &&
              m.find("translucency") == std::string::npos);
    CHECK(json.find("\"fibre_body_") == std::string::npos && json.find("\"core_") == std::string::npos);
    // One group per construction: every mapped tube in the yarn body's group.
    const std::regex group("\"name\": \"yarn_body\"");
    CHECK(std::regex_search(json, group));
#ifdef KNIT_PPT_EXE
    if (fs::exists(KNIT_PPT_EXE)) {
        for (const char* run : {"first", "shiny"}) {
            const std::string check = d + "/dump_" + run + ".log";
            const std::string cmd = std::string("\"\"") + KNIT_PPT_EXE + "\" dump-material \"" + d + "/" + run + "/scene.json\" > \"" + check + "\" 2>&1\"";
            const int rr = std::system(cmd.c_str());
            if (rr != 0) std::printf("%s\n", readText(check).c_str());
            CHECK(rr == 0);
        }
        // The renderer reads the fabric's highlight with its maps (photon_path_tracer dd8870fda; an older build
        // ignores the block silently, which this catches).
        const size_t idAt = fabric.empty() ? std::string::npos : fabric[0].find("{\"id\": \"") + 8;
        if (idAt != std::string::npos) {
            const std::string id = fabric[0].substr(idAt, fabric[0].find('"', idAt) - idAt), check = d + "/dump_fabric.log";
            const std::string cmd = std::string("\"\"") + KNIT_PPT_EXE + "\" dump-material \"" + out + "/scene.json\" " + id + " > \"" + check + "\" 2>&1\"";
            CHECK(std::system(cmd.c_str()) == 0);
            const std::string dump = readText(check);
            const size_t specularAt = dump.find("\nspecular:");
            CHECK(dump.find("type: fabric") != std::string::npos && specularAt != std::string::npos &&
                  dump.find("specular (tex): textures/baked_", specularAt) != std::string::npos &&
                  dump.find("roughness (tex): textures/baked_", specularAt) != std::string::npos);
        }
    }
#endif
}

// A preset under Fur exports like an import, its mesh under the pile: no rule of its own blocks it, only the fur it
// has not grown yet.
TEST(preset_fur_export_like_an_import) {
    for (const char* id : {"teddy", "teddy-brown"}) {
        const SceneEntry* e = findScene(id);
        CHECK(e != nullptr && e->isPreset());
        if (!e) continue;
        SceneSettings fur = viewerDefaults(*e);
        fur.workflow = "fur";
        const ExportEligibility el = exportEligibility(fur, nullptr, nullptr, FibreChoice::Included);
        std::printf("    %s: fix \"%s\"\n", id, el.fix.c_str());
        CHECK(el.state == ExportEligibility::State::Blocked && el.fix == "Grow fur" && el.reasons.size() == 1);
    }
}

// Every preset names its model in full, what CC BY and CC BY-NC ask of a work made from it: title, author, source,
// licence and licence URL, all in the credit text its exports carry (CREDITS.txt), with what each workflow changed.
TEST(every_preset_names_its_model) {
    for (const SceneEntry& e : sceneRegistry()) {
        const ModelCredit& c = e.credit;
        CHECK(e.isPreset() && !c.author.empty() && !c.source.empty() && !c.licence.empty() && c.licenceUrl.rfind("https://", 0) == 0);
        std::printf("    %s: %s\n", e.id.c_str(), creditLine(c).c_str());
        for (const Workflow w : {Workflow::Knit, Workflow::Fur}) {
            const std::string text = creditText(c, w);
            for (const std::string& part : {c.title, c.author, c.source, c.licence, c.licenceUrl, std::string("Changed: ")})
                CHECK(text.find(part) != std::string::npos);
            CHECK(text.find(w == Workflow::Knit ? "knitted yarn" : "fur") != std::string::npos);
        }
    }
}

// The blue swatch (scene_registry.cpp): the wool swatch's knit with its yarn made of fibres by default (the
// re-knit's blue look: 400 fibres of 36 um per yarn across 3 plies, an opaque core, flattened by at most 0.15, the
// fitted blue dye), the same knit keys as the wool swatch, and fibres the export accepts at its drawn yarn.
TEST(blue_swatch_preset) {
    const SceneEntry* blue = findScene("wool-swatch-blue");
    const SceneEntry* cream = findScene("wool-swatch");
    CHECK(blue && cream);
    if (!blue || !cream) return;
    const SceneSettings s = viewerDefaults(*blue), c = viewerDefaults(*cream);
    CHECK(yarnConstructionOf(s) == YarnConstruction::FibreStrands && yarnLookConflict(s).empty());
    CHECK(s.yarnBody.plies == 3 && s.yarnBody.fibres == 400 && s.yarnBody.fibreDiameterUm == 36 && s.yarnBody.core == "opaque");
    CHECK(s.maxFlattening == 0.15 && !s.filling && s.yarnColourSource == "single");
    CHECK(s.yarnWhite.color.x == 0.0005 && s.yarnWhite.color.y == 0.1523 && s.yarnWhite.color.z == 0.4479);
    const YarnBodyCheck check = checkYarnBody(s);
    if (!check.accepted) std::printf("    %s\n", check.reason.c_str());
    CHECK(check.accepted);
    // The same stitches as the wool swatch: besides the scene id, only the flattening limit is a knitting difference.
    SceneSettings creamLimited = c;
    creamLimited.scene = s.scene;
    creamLimited.maxFlattening = s.maxFlattening;
    CHECK(stageKey(s, Stage::Knit, 1) != stageKey(c, Stage::Knit, 1));
    CHECK(stageKey(s, Stage::Knit, 1) == stageKey(creamLimited, Stage::Knit, 1));
    CHECK(!c.yarnBody.enabled);
}

// A yarn of fibres the generator refuses (review 2026-11-11: an opaque core at shell depth 1) names the
// control and a remedy, and blocks the export before it starts, with that remedy among the reasons. The
// yarn body's keys are no knitting keys: the knit stays current throughout.
TEST(yarn_body_refusal_blocks_export) {
    const SceneEntry* blue = findScene("wool-swatch-blue");
    CHECK(blue != nullptr);
    if (!blue) return;
    SceneSettings s = viewerDefaults(*blue);
    KnitResult k;
    k.id = 3;
    k.settings = s;
    k.final = true;
    k.parts = 1;
    using St = ExportEligibility::State;
    CHECK(exportEligibility(s, &k, nullptr, FibreChoice::LeftOut).state == St::Ready);
    s.yarnBody.shellDepth = 1;
    const YarnBodyCheck c = checkYarnBody(s);
    std::printf("    refused: %s / %s\n", c.reason.c_str(), c.remedy.c_str());
    CHECK(!c.accepted && c.field == YarnBodyField::Core && c.remedy.find("Core: None") != std::string::npos);
    const ExportEligibility el = exportEligibility(s, &k, nullptr, FibreChoice::LeftOut);
    CHECK(el.state == St::Blocked && el.fix == "Correct the yarn's fibres");
    bool remedyListed = false;
    for (const std::string& r : el.reasons) remedyListed = remedyListed || r.find(c.remedy) != std::string::npos;
    CHECK(remedyListed);
    // Baked maps write no core and never bake one (H4, 2026-11-14): the same opaque core at shell depth 1 is no refusal
    // under yarn_body.render = maps, and the parameters carry no core there.
    SceneSettings maps = s;
    maps.yarnBody.render = "maps";
    YarnBodyParams mapped;
    CHECK(checkYarnBody(maps).accepted && yarnBodyParams(maps, mapped, nullptr) && mapped.core == BodyCore::None);
    CHECK(exportEligibility(maps, &k, nullptr, FibreChoice::LeftOut).state == St::Ready);
    s.yarnBody.core = "none";
    CHECK(checkYarnBody(s).accepted && exportEligibility(s, &k, nullptr, FibreChoice::LeftOut).state == St::Ready);
    s.yarnBody.fibres = 2000;  // more than the shells hold: the count is the control named
    CHECK(checkYarnBody(s).field == YarnBodyField::Fibres);
}

namespace {
// Vertex positions of a package's binary little-endian PLY (x y z nx ny nz [u v] float32).
std::vector<Vec3> packagePlyPositions(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    std::string line;
    size_t n = 0;
    int props = 0;
    while (std::getline(f, line) && line != "end_header") {
        if (line.rfind("element vertex ", 0) == 0) n = std::stoul(line.substr(15));
        if (line.rfind("property float ", 0) == 0) ++props;
    }
    std::vector<Vec3> out;
    std::vector<float> v(props);
    for (size_t i = 0; i < n && f.read(reinterpret_cast<char*>(v.data()), props * sizeof(float)); ++i) out.push_back({v[0], v[1], v[2]});
    return out;
}

// Every source vertex of the scene as the package writes it (metres, Y up), sorted by x.
std::vector<Vec3> sourcePointsMetres(const Scene& t) {
    std::vector<Vec3> out;
    for (const SourcePrimitive& g : t.source.prims)
        for (const Vec3& p : g.welded.pos) out.push_back(sceneMetres({}).toExport(p));
    std::sort(out.begin(), out.end(), [](const Vec3& a, const Vec3& b) { return a.x < b.x; });
    return out;
}

// How many of `points` lie within 1 um of a source vertex (float32 positions of a 10 cm model are
// exact to about 10 nm; a tube's vertices lie a yarn radius away from the surface).
size_t pointsOnSource(const std::vector<Vec3>& points, const std::vector<Vec3>& source) {
    const double tol = 1e-6;
    size_t hits = 0;
    for (const Vec3& p : points) {
        auto it = std::lower_bound(source.begin(), source.end(), p.x - tol, [](const Vec3& a, double x) { return a.x < x; });
        for (; it != source.end() && it->x <= p.x + tol; ++it)
            if (norm(*it - p) <= tol) {
                ++hits;
                break;
            }
    }
    return hits;
}

// The positions of every mesh of a package except the ground, then those of its OBJ.
std::vector<Vec3> packageMeshPoints(const std::string& dir) {
    std::vector<Vec3> out;
    for (const std::string& f : filesUnder(dir)) {
        if (f.rfind("meshes/", 0) != 0 || f.find("_ground") != std::string::npos) continue;
        const std::vector<Vec3> p = packagePlyPositions(dir + "/" + f);
        out.insert(out.end(), p.begin(), p.end());
    }
    return out;
}
std::vector<Vec3> objPoints(const std::string& path) {
    std::vector<Vec3> out;
    std::ifstream f(path);
    for (std::string line; std::getline(f, line);) {
        if (line.rfind("v ", 0) != 0) continue;
        std::istringstream ls(line.substr(2));
        Vec3 p;
        ls >> p.x >> p.y >> p.z;
        out.push_back(p);
    }
    return out;
}

// One spiral piece that needs no scene files: the chart fixture and its manifest in dir.
void writeFixturePiece(const std::string& dir) {
    fs::copy_file(std::string(KNIT_TEST_DATA) + "/chart_fixture.obj", dir + "/piece_fixture.obj");
    std::ofstream m(dir + "/manifest.txt");
    m << "piece piece_fixture.obj object=fixture island=0 material=fabric_white rings=5 segments=32 area=1 centroid=0,0,0 poleA=0,0,-100 poleB=0,0,100\n";
}
} // namespace

// One spiral piece that needs no scene files (the committed chart fixture, 24 stitches around): a Final
// knit as the viewer runs it, Standard wool fibres on it, the package with fibres (the piece's tube
// and fibre parts, the ground), then the Wavefront OBJ without fibres, read back by the importer
// (import_analysis.h): no fatal problem, every face of the file, every material's colour from the .mtl,
// no fibre lines, the size of the source, and in neither file a vertex of the source or a credit
// (this scene is no preset).
TEST(pieces_knit_fibres_export) {
    if (!gpu::available(nullptr)) SKIP("no usable CUDA GPU");
    const std::string pieces = freshDir("pieces_scene"), cache = freshDir("pieces_scenes") + "/scenes", out = freshDir("pieces_out");
    writeFixturePiece(pieces);
    std::string err;
    const SceneEntry* e = manifestScene(pieces + "/manifest.txt", "piece_fixture.obj", 24, &err);
    CHECK(e != nullptr);
    if (!e) return;
    SceneSettings s = viewerDefaults(*e);
    applyFullProfile(s, *e);
    s.timeCapS = 900;
    auto b = run(BuiltScene::Kind::Knit, s, &err, nullptr, nullptr, nullptr, "", "ppt", false, cache);
    if (!b || !b->knit || !b->knit->scene) {
        std::printf("    knit: %s\n", err.c_str());
        CHECK(false);
        return;
    }
    const Scene& t = *b->knit->scene;
    CHECK(b->knit->final && t.parts.size() == 1 && t.parts[0].layout == KnitLayout::Spiral && t.parts[0].knit != nullptr);
    std::printf("    knit: %d stitches, verdict %s, %.1f s\n", b->parts.empty() ? 0 : b->parts[0].stitches,
                b->parts.empty() ? "-" : toString(b->parts[0].verdict), b->ms / 1000);
    // Fibres: the Standard wool preset (a manifest scene's default).
    s.fibres = true;
    CHECK(s.fibrePreset == "standard");
    auto f = run(BuiltScene::Kind::Fibres, s, &err, nullptr, b->knit);
    if (!f || !f->fibreResult) {
        std::printf("    fibres: %s\n", err.c_str());
        CHECK(false);
        return;
    }
    CHECK(f->fibreResult->fibres > 0);
    std::printf("    fibres: %lld\n", static_cast<long long>(f->fibreResult->fibres));
    const ExportEligibility el = exportEligibility(s, b->knit.get(), f->fibreResult.get(), FibreChoice::Included);
    CHECK(el.state != ExportEligibility::State::Blocked);
    for (const std::string& r : el.reasons) std::printf("    %s\n", r.c_str());
    // The package, with the fibres.
    auto x = run(BuiltScene::Kind::Export, s, &err, nullptr, b->knit, f->fibreResult, out + "/package", "ppt", true);
    if (!x) std::printf("    package: %s\n", err.c_str());
    CHECK(x != nullptr);
    const std::string json = readText(out + "/package/scene.json");
    auto has = [&](const std::string& name) { return json.find("\"name\": \"" + name + "\"") != std::string::npos; };
    CHECK(has("piece_fixture") && has("piece_fixture_fibres_oriented") && has("Ground"));
    CHECK(has("piece_fixture_stuffing") || has("piece_fixture_backing"));
    const std::vector<Vec3> source = sourcePointsMetres(t);
    const std::vector<Vec3> mesh = packageMeshPoints(out + "/package");
    CHECK(!mesh.empty() && pointsOnSource(mesh, source) == 0);
    CHECK(!fs::exists(out + "/package/CREDITS.txt"));
    // The OBJ, without fibres.
    s.fibres = false;
    const std::string target = out + "/knitted_fixture";
    auto o = run(BuiltScene::Kind::Export, s, &err, nullptr, b->knit, nullptr, target, "obj", false);
    if (!o) std::printf("    obj: %s\n", err.c_str());
    CHECK(o != nullptr);
    const std::string base = fs::path(target).filename().string();
    const std::string objPath = target + "/" + base + ".obj";
    CHECK(fs::exists(objPath) && fs::exists(target + "/" + base + ".mtl") && !fs::exists(target + "/CREDITS.txt"));
    CHECK(pointsOnSource(objPoints(objPath), source) == 0);
    // Read back as an import (the fibres' "l" lines would be listed as ignored: there are none).
    int faceLines = 0, lineElements = 0;
    {
        std::ifstream objFile(objPath);
        for (std::string line; std::getline(objFile, line);) {
            faceLines += line.rfind("f ", 0) == 0;
            lineElements += line.rfind("l ", 0) == 0;
        }
    }
    ImportLimits lim;
    lim.maxFileBytes = 4ull << 30;
    lim.maxTriangles = lim.maxVertices = 50000000;
    ImportedMesh m;
    CHECK(readMeshFile(objPath, lim, m, &err));
    const ImportAnalysis an = analyzeImport(m);
    const Vec3 ext = an.hi - an.lo;
    std::printf("    OBJ %s: %.1f MB, %zu vertices, %d faces (%zu triangles), %zu materials, %d components; extent %.4f x %.4f x %.4f m\n",
                objPath.c_str(), double(m.fileBytes) / 1048576.0, m.pos.size(), m.faces, m.tri.size(), m.materials.size(), an.components,
                ext.x, ext.y, ext.z);
    for (const ImportIssue& i : an.issues)
        if (i.kind == ImportIssue::Kind::Fatal) std::printf("    fatal: %s\n", i.text.c_str());
    CHECK(!an.fatal() && m.faces == faceLines && faceLines > 0 && lineElements == 0);
    CHECK(m.missingDependencies.empty() && !m.materials.empty());
    for (const ImportMaterial& mat : m.materials) CHECK(mat.hasKd);
    // The size of the source: the yarn lies on the surface, a few mm either way.
    Vec3 lo{1e30, 1e30, 1e30}, hi{-1e30, -1e30, -1e30};
    for (const Vec3& q : source) {
        lo = {std::min(lo.x, q.x), std::min(lo.y, q.y), std::min(lo.z, q.z)};
        hi = {std::max(hi.x, q.x), std::max(hi.y, q.y), std::max(hi.z, q.z)};
    }
    const Vec3 size = hi - lo;
    std::printf("    source extent %.4f x %.4f x %.4f m\n", size.x, size.y, size.z);
    CHECK(std::fabs(ext.x - size.x) < 0.012 && std::fabs(ext.y - size.y) < 0.012 && std::fabs(ext.z - size.z) < 0.012);
}

// Two pieces on the chart fixture, the first piece knitted at Draft (a manifest without a GLB: no kept parts).
const SceneEntry* knitTwoPieceScene(const std::string& d, SceneSettings& s, Scene& t, std::string& err) {
    fs::copy_file(std::string(KNIT_TEST_DATA) + "/chart_fixture.obj", d + "/piece_a.obj");
    fs::copy_file(std::string(KNIT_TEST_DATA) + "/chart_fixture.obj", d + "/piece_b.obj");
    {
        std::ofstream m(d + "/manifest.txt");
        for (const char* piece : {"piece_a.obj", "piece_b.obj"})
            m << "piece " << piece << " object=fixture island=0 material=fabric_white rings=5 segments=32 area=1 centroid=0,0,0 poleA=0,0,-100 poleB=0,0,100\n";
    }
    const SceneEntry* e = manifestScene(d + "/manifest.txt", "piece_a.obj", 24, &err);
    if (!e) return nullptr;
    applySceneDefaults(s, *e);
    applyPreviewProfile(s);
    KnitRunOptions ro;
    ro.only = "piece_a";
    if (!knitScene(*e, s, ro, t, &err)) std::printf("    knit: %s\n", err.c_str());
    return e;
}

// Two bears of one piece each: the fixture piece, and beside it the same piece moved by three of its
// widths (so neither hides the other's fibres). Both knitted (Draft).
const SceneEntry* knitTwoObjectScene(const std::string& d, SceneSettings& s, Scene& t, std::string& err) {
    fs::copy_file(std::string(KNIT_TEST_DATA) + "/chart_fixture.obj", d + "/piece_a.obj");
    {
        std::ifstream in(std::string(KNIT_TEST_DATA) + "/chart_fixture.obj");
        std::vector<std::string> lines;
        double x0 = 1e300, x1 = -1e300;
        for (std::string line; std::getline(in, line);) {
            lines.push_back(line);
            if (line.rfind("v ", 0) == 0) {
                const double x = std::stod(line.substr(2));
                x0 = std::min(x0, x), x1 = std::max(x1, x);
            }
        }
        std::ofstream out(d + "/piece_b.obj");
        for (const std::string& line : lines) {
            if (line.rfind("v ", 0) != 0) {
                out << line << "\n";
                continue;
            }
            std::istringstream v(line.substr(2));
            double x = 0, y = 0, z = 0;
            v >> x >> y >> z;
            out << "v " << x + 3 * (x1 - x0) << " " << y << " " << z << "\n";
        }
        std::ofstream m(d + "/manifest.txt");
        m << "piece piece_a.obj object=bear_a island=0 material=fabric_white rings=5 segments=32 area=1 centroid=0,0,0 poleA=0,0,-100 poleB=0,0,100\n";
        m << "piece piece_b.obj object=bear_b island=0 material=fabric_white rings=5 segments=32 area=1 centroid=0,0,0 poleA=0,0,-100 poleB=0,0,100\n";
    }
    const SceneEntry* e = manifestScene(d + "/manifest.txt", "piece_a.obj", 24, &err);
    if (!e) return nullptr;
    applySceneDefaults(s, *e);
    applyPreviewProfile(s);
    if (!knitScene(*e, s, KnitRunOptions{}, t, &err)) std::printf("    knit: %s\n", err.c_str());
    return e;
}

// fibres.segment_budget holds per object: with a budget below either bear's fibres, each bear keeps
// about the whole budget, and the first bear's fibres are the same whether the second bear has stitches
// or not (its stitches taken away on the same knit: the solver is not deterministic).
TEST(fibre_budget_per_object) {
    if (!gpu::available(nullptr)) SKIP("no usable CUDA GPU");
    const std::string d = freshDir("budget_objects");
    std::string err;
    SceneSettings s;
    Scene t;
    const SceneEntry* e = knitTwoObjectScene(d, s, t, err);
    CHECK(e != nullptr && t.parts.size() == 2 && t.parts[0].hasStitches() && t.parts[1].hasStitches());
    if (!e || t.parts.size() != 2 || !t.parts[0].hasStitches() || !t.parts[1].hasStitches()) return;
    CHECK(sceneObjectOf(t.source, t.parts[0].prim) != sceneObjectOf(t.source, t.parts[1].prim));
    s.fibres = true;
    s.fibreDensity = 40;
    s.fibreGeneration = "budget";
    s.segmentBudget = 0;  // unlimited: every fibre, to choose the budget from
    SceneFibres all;
    CHECK(grown(all, growWoolFibres(t, s, FibreRunOptions{}, &err)));
    CHECK(all.parts.size() == 2);
    if (all.parts.size() != 2) return;
    const long long segmentsA = segmentCount(all.parts[0].fs), segmentsB = segmentCount(all.parts[1].fs);
    s.segmentBudget = std::min(segmentsA, segmentsB) / 2;
    SceneFibres both, alone;
    CHECK(grown(both, growWoolFibres(t, s, FibreRunOptions{}, &err)));
    std::unique_ptr<PartKnit> held;
    std::swap(held, t.parts[1].knit);  // the second bear without stitches
    CHECK(grown(alone, growWoolFibres(t, s, FibreRunOptions{}, &err)));
    std::swap(held, t.parts[1].knit);
    CHECK(both.parts.size() == 2 && alone.parts.size() == 1);
    if (both.parts.size() != 2 || alone.parts.size() != 1) return;
    const auto keptBelow = [](const PartFibres& pf) {
        std::vector<const Fibre*> kept;
        for (const Fibre& f : pf.fs.fibres)
            if (float(f.key) < pf.keyBelow) kept.push_back(&f);
        return kept;
    };
    const auto segmentsOf = [](const std::vector<const Fibre*>& fibres) {
        long long n = 0;
        for (const Fibre* f : fibres) n += f->segments();
        return n;
    };
    const std::vector<const Fibre*> a2 = keptBelow(both.parts[0]), a1 = keptBelow(alone.parts[0]), b2 = keptBelow(both.parts[1]);
    std::printf("    budget %lld per object (bears %lld and %lld segments): kept %lld and %lld; the first bear alone %lld\n", s.segmentBudget,
                segmentsA, segmentsB, segmentsOf(a2), segmentsOf(b2), segmentsOf(a1));
    CHECK(both.parts[0].keyBelow == alone.parts[0].keyBelow && a1.size() == a2.size());
    for (size_t i = 0; i < std::min(a1.size(), a2.size()); ++i)
        CHECK(a1[i]->candidate == a2[i]->candidate && fibrePoints(alone.parts[0].fs, *a1[i]) == fibrePoints(both.parts[0].fs, *a2[i]));
    CHECK(segmentsOf(a2) > 0.9 * s.segmentBudget && segmentsOf(b2) > 0.9 * s.segmentBudget);
}

// The same scene as a preset (a credit), whose detail is the unknitted piece.
SceneEntry presetOfTwoPieceScene(const SceneEntry& e) {
    SceneEntry preset = e;
    preset.credit = {"Fixture", "Someone", "https://example.org/fixture"};
    return preset;
}

// A preset's .blend holds its credit as a text block, because the file may travel without
// CREDITS.txt: the credit must be in the package folder when Blender runs. Blender reads it back.
TEST(blender_export_keeps_the_credit) {
    if (!findExportFormat("blend")->available) SKIP("no Blender 5.1 on this machine");
    if (!gpu::available(nullptr)) SKIP("no usable CUDA GPU");
    const std::string d = freshDir("blend_credit");
    std::string err, log;
    SceneSettings s;
    Scene t;
    const SceneEntry* e = knitTwoPieceScene(d, s, t, err);
    CHECK(e != nullptr);
    if (!e) return;
    const SceneEntry preset = presetOfTwoPieceScene(*e);
    ExportOptions o;
    o.format = "blend";
    o.folder = d + "/blend";
    o.name = "credit";
    CHECK(exportScene({preset, t, s}, o, &log, &err));
    std::printf("    %s\n", err.c_str());
    CHECK(fs::exists(o.folder + "/CREDITS.txt"));
    writeText(d + "/read_credit.py", "import bpy, sys\n"
                                     "bpy.ops.wm.open_mainfile(filepath=sys.argv[sys.argv.index('--') + 1])\n"
                                     "sys.exit(0 if 'CREDITS.txt' in bpy.data.texts else 3)\n");
    BlenderScriptRun run;
    run.blender = findBlender();
    run.script = d + "/read_credit.py";
    run.arguments = {fs::absolute(o.folder + "/credit.blend").string()};
    run.logPath = d + "/read_credit.log";
    run.timeCapSeconds = 120;
    BlenderScriptResult result;
    CHECK(runBlenderScript(run, result, &err));
    std::printf("    Blender exit %d\n", result.exitCode);
    CHECK(result.exitCode == 0);
}

// A wool profile's translucent yarn (Soft: 0.3) reaches the Blender file as a Translucent BSDF in the
// fabric group, carrying the transmit share, and the close-up cameras (mid, far) are camera objects.
// Before, the conversion stopped on the package's unknown "translucency" key.
TEST(blender_export_translucent_yarn_and_cameras) {
    if (!findExportFormat("blend")->available) SKIP("no Blender 5.1 on this machine");
    if (!gpu::available(nullptr)) SKIP("no usable CUDA GPU");
    const std::string d = freshDir("blend_translucent");
    std::string err, log;
    SceneSettings s;
    Scene t;
    const SceneEntry* e = knitTwoPieceScene(d, s, t, err);
    CHECK(e != nullptr);
    if (!e) return;
    s.yarnWhite.translucency = 0.3;
    ExportOptions o;
    o.format = "blend";
    o.folder = d + "/blend";
    o.name = "translucent";
    const bool written = exportScene({*e, t, s}, o, &log, &err);
    if (!written) std::printf("    %s\n", err.c_str());
    CHECK(written);
    if (!written) return;
    writeText(d + "/read_blend.py", "import bpy, sys\n"
                                    "bpy.ops.wm.open_mainfile(filepath=sys.argv[sys.argv.index('--') + 1])\n"
                                    "cams = {o.name for o in bpy.data.objects if o.type == 'CAMERA'}\n"
                                    "groups = [g for g in bpy.data.node_groups if g.name.endswith(' t0.3')]\n"
                                    "lucent = any(n.bl_idname == 'ShaderNodeBsdfTranslucent' for g in groups for n in g.nodes)\n"
                                    "sys.exit(0 if lucent and {'mid', 'far'} <= cams else 3)\n");
    BlenderScriptRun run;
    run.blender = findBlender();
    run.script = d + "/read_blend.py";
    run.arguments = {fs::absolute(o.folder + "/translucent.blend").string()};
    run.logPath = d + "/read_blend.log";
    run.timeCapSeconds = 120;
    BlenderScriptResult result;
    CHECK(runBlenderScript(run, result, &err));
    std::printf("    Blender exit %d\n", result.exitCode);
    CHECK(result.exitCode == 0);
}

// A yarn of fibres rendered as maps reaches Blender in both map materials (yarn_body.map_material); the yarn tubes
// stay meshes. Blender opens each .blend and counts:
//   fabric (the default, here translucent): the package's fabric node group with its translucency and highlight
//     ("PPT Fabric ... t0.4 h", a Principled BSDF inside) under a normal map, four image textures (the cavity map on
//     the albedo, the normal map, the roughness and reflectance maps on the highlight) and the uv tangent;
//   shiny: a Principled BSDF under a normal map with the cavity, roughness and reflectance maps on their inputs
//     (four image textures) and the uv tangent.
TEST(blender_export_baked_yarn) {
    if (!findExportFormat("blend")->available) SKIP("no Blender 5.1 on this machine");
    if (!gpu::available(nullptr)) SKIP("no usable CUDA GPU");
    const std::string d = freshDir("blend_baked"), root = freshDir("blend_baked_cache");
    ImportCacheEnv env(root);
    SceneSettings s = importPlane(d, root);
    applyFullProfile(s, *sceneEntryFor(s));  // an export needs a Final knit
    s.fibres = false;
    s.filling = false;
    s.yarnBody.enabled = true;
    s.yarnBody.plies = 2;
    s.yarnBody.fibres = 30;
    s.yarnBody.render = "maps";
    s.yarnBody.mapTexelUm = 100;  // coarse: a fast bake
    s.yarnWhite.translucency = s.yarnBrown.translucency = 0.4;
    std::string err;
    auto fin = run(BuiltScene::Kind::Knit, s, &err);
    CHECK(fin && fin->knit);
    if (!fin) return;
    // The script counts the baked yarns' materials (named baked_*) of the asked kind; exit 3 when there are none or
    // one of another kind.
    writeText(d + "/read_blend.py",
              "import bpy, sys\n"
              "args = sys.argv[sys.argv.index('--') + 1:]\n"
              "bpy.ops.wm.open_mainfile(filepath=args[0])\n"
              "good = bad = 0\n"
              "for m in bpy.data.materials:\n"
              "    if not m.use_nodes or not m.name.startswith('baked_'): continue\n"
              "    kinds = [n.bl_idname for n in m.node_tree.nodes]\n"
              "    trees = [n.node_tree for n in m.node_tree.nodes if n.bl_idname == 'ShaderNodeGroup']\n"
              "    groups = [t.name for t in trees]\n"
              "    tangent = any(n.bl_idname == 'ShaderNodeTangent' and n.direction_type == 'UV_MAP' for n in m.node_tree.nodes)\n"
              "    if args[1] == 'shiny':\n"
              "        ok = 'ShaderNodeBsdfPrincipled' in kinds and 'ShaderNodeNormalMap' in kinds and kinds.count('ShaderNodeTexImage') == 4\n"
              "    else:\n"
              "        ok = len(groups) == 1 and groups[0].startswith('PPT Fabric') and groups[0].endswith(' t0.4 h') and "
              "any(n.bl_idname == 'ShaderNodeBsdfPrincipled' for n in trees[0].nodes) and "
              "'ShaderNodeNormalMap' in kinds and kinds.count('ShaderNodeTexImage') == 4\n"
              "    ok = ok and tangent\n"
              "    print(m.name, 'ok' if ok else 'unexpected', kinds, groups)\n"
              "    good, bad = good + ok, bad + (not ok)\n"
              "print('baked materials', good, 'unexpected', bad)\n"
              "sys.exit(0 if good > 0 and bad == 0 else 3)\n");
    for (const char* material : {"fabric", "shiny"}) {
        SceneSettings e = s;
        e.yarnBody.mapMaterial = material;
        const std::string dir = d + "/blend_" + material;
        auto exported = run(BuiltScene::Kind::Export, e, &err, nullptr, fin->knit, nullptr, dir, "blend");
        if (!exported) std::printf("    %s\n", err.c_str());
        CHECK(exported != nullptr);
        if (!exported) return;
        std::string blend;
        for (const std::string& f : filesUnder(dir))
            if (f.size() > 6 && f.compare(f.size() - 6, 6, ".blend") == 0) blend = dir + "/" + f;
        std::printf("    %s: %s\n", material, blend.c_str());
        CHECK(!blend.empty());
        if (blend.empty()) return;
        BlenderScriptRun br;
        br.blender = findBlender();
        br.script = d + "/read_blend.py";
        br.arguments = {fs::absolute(blend).string(), material};
        br.logPath = d + "/read_blend_" + material + ".log";
        br.timeCapSeconds = 120;
        BlenderScriptResult result;
        CHECK(runBlenderScript(br, result, &err));
        std::printf("    Blender exit %d (%s)\n", result.exitCode, material);
        if (result.exitCode != 0) std::printf("%s\n", readText(br.logPath).c_str());
        CHECK(result.exitCode == 0);
    }
}

namespace {
// The distance [m] from each point to its nearest source vertex (both in export metres).
std::vector<double> distancesToSource(const std::vector<Vec3>& points, const std::vector<Vec3>& source) {
    std::vector<double> out;
    for (const Vec3& p : points) {
        double best = 1e30;
        for (const Vec3& q : source) best = std::min(best, norm(q - p));
        out.push_back(best);
    }
    return out;
}
} // namespace

// A preset exports like an import, with its own details. Two pieces, one knitted and one outside
// the part filter (the fallback every failed, too small or time-capped part takes), and a plain extra:
// the unknitted piece and the extra are their source meshes vertex for vertex; the knitted piece has
// its filling behind the yarn at fillingOffset; report.json lists nothing left out; CREDITS.txt names
// the model. Without the filling, no filling mesh or material is written.
TEST(preset_export_writes_its_details) {
    if (!gpu::available(nullptr)) SKIP("no usable CUDA GPU");
    const std::string d = freshDir("preset_details");
    std::string err, log;
    SceneSettings s;
    Scene t;
    const SceneEntry* e = knitTwoPieceScene(d, s, t, err);
    CHECK(e != nullptr);
    if (!e) return;
    CHECK(t.parts.size() == 2 && t.parts[0].hasStitches() && !t.parts[1].hasStitches());
    if (t.parts.size() != 2 || !t.parts[0].hasStitches()) return;
    const SceneEntry preset = presetOfTwoPieceScene(*e);
    CHECK(s.filling);
    ExportOptions o;
    o.folder = d + "/package";
    CHECK(exportScene({preset, t, s}, o, &log, &err));
    std::printf("    %s\n", err.c_str());
    CHECK(log.find("left out") == std::string::npos);
    const std::vector<Vec3> source = sourcePointsMetres(t);
    for (const char* plain : {"piece_b"}) {
        const std::vector<Vec3> p = packagePlyPositions(o.folder + "/meshes/" + plain + ".ply");
        std::printf("    %s: %zu of %zu vertices on the source\n", plain, pointsOnSource(p, source), p.size());
        CHECK(!p.empty() && pointsOnSource(p, source) == p.size());
    }
    const std::string fillingName = fs::exists(o.folder + "/meshes/piece_a_stuffing.ply") ? "piece_a_stuffing" : "piece_a_backing";
    const double offset = fillingOffset(s, t.parts[0].knit->t.pp.yarnDiameter) * 1e-3;
    const std::vector<double> filling = distancesToSource(packagePlyPositions(o.folder + "/meshes/" + fillingName + ".ply"), source);
    CHECK(!filling.empty());
    if (!filling.empty()) {
        const auto [nearest, farthest] = std::minmax_element(filling.begin(), filling.end());
        std::printf("    %s: %zu vertices, %.4f to %.4f mm from the source (offset %.4f mm)\n", fillingName.c_str(), filling.size(),
                    *nearest * 1e3, *farthest * 1e3, offset * 1e3);
        CHECK(*nearest > 1e-6 && *farthest < offset + 1e-6);
    }
    const std::string json = readText(o.folder + "/scene.json");
    for (const char* name : {"\"piece_b\"", "\"stuffing\""}) CHECK(json.find(name) != std::string::npos);
    CHECK(readText(o.folder + "/report.json").find("\"omitted\"") == std::string::npos);
    CHECK(readText(o.folder + "/CREDITS.txt").find("Fixture") != std::string::npos);
    // The OBJ writes the same surfaces.
    ExportOptions objExport;
    objExport.format = "obj";
    objExport.folder = d + "/obj";
    objExport.name = "knit";
    CHECK(exportScene({preset, t, s}, objExport, &log, &err));
    CHECK(pointsOnSource(objPoints(objExport.folder + "/knit.obj"), source) > 0);
    // No filling: neither its mesh nor its material.
    s.filling = false;
    ExportOptions bare;
    bare.folder = d + "/bare";
    CHECK(exportScene({preset, t, s}, bare, &log, &err));
    const std::string bareJson = readText(bare.folder + "/scene.json");
    CHECK(bareJson.find("\"stuffing\"") == std::string::npos && bareJson.find(fillingName) == std::string::npos);
    CHECK(bareJson.find("\"piece_b\"") != std::string::npos);
}

// With knitting on and no part knitted, a preset exports like an import: its parts as plain meshes.
TEST(everything_unknitted_exports_plain_meshes) {
    const std::string d = freshDir("nothing_knitted");
    const SceneEntry* man = findScene("teddy");
    CHECK(man != nullptr);
    if (!man) return;
    knit_test::requireScene(findProjectFile(man->path));
    const SceneSettings s = viewerDefaults(*man);
    KnitResult k;
    k.settings = s;
    k.final = true;
    k.parts = k.plain = 2;
    const ExportEligibility el = exportEligibility(s, &k, nullptr, FibreChoice::LeftOut);
    for (const std::string& r : el.reasons) std::printf("    %s\n", r.c_str());
    CHECK(el.state != ExportEligibility::State::Blocked);
    Scene t;
    std::string err, log;
    CHECK(loadBaseParts(*man, s, SourceOptions{}, t, &err));
    ExportOptions o;
    o.folder = d + "/out";
    CHECK(exportScene({*man, t, s}, o, &log, &err));
    const std::vector<Vec3> mesh = packageMeshPoints(o.folder);
    std::printf("    %zu of %zu exported vertices on the source\n", pointsOnSource(mesh, sourcePointsMetres(t)), mesh.size());
    CHECK(!mesh.empty() && pointsOnSource(mesh, sourcePointsMetres(t)) == mesh.size());
}

// The saved fibres' key check: names the first stale line, values to rounding only (a last-bit
// difference is not a changed input; a changed setting is), a line cut by the former writer's buffer
// still matches, a longer stored line matches, a file without a placement line is legacy.
TEST(saved_fibre_key_check) {
    const std::string have = "strands = 00ff 1\nexport tube = 0x1p-2 0x1p-1 6 0x1p-9 metres 0 0 0\nchecks = 0 0x1.8afibre params = 0x1.4p+3 7\n";
    CHECK(staleFibreKeyField(have, "strands = 00ff 1\nexport tube = 0x1p-2 0x1p-1 6 \nfibre params = 0x1.4p+3 7\nfibre placement = legacy\n").empty());
    CHECK(staleFibreKeyField(have, "strands = 00fe 1\n") == "strands");
    CHECK(staleFibreKeyField(have, "export tube = 0x1p-2 0x1p-1 8 \n") == "export tube");
    CHECK(staleFibreKeyField(have, "fibre params = 0x1.4000000000001p+3 7\n").empty());
    CHECK(staleFibreKeyField(have, "fibre params = 0x1.41p+3 7\n") == "fibre params");
    CHECK(staleFibreKeyField(have, "fibre params = 0x1.4p+3 7 8\n") == "fibre params");
    CHECK(staleFibreKeyField(have, "fibre placement = exposed\n") == "fibre placement");
    CHECK(staleFibreKeyField(have, "fibre backend = gpu\n") == "fibre backend");
    // Fibres saved by the v2 generator (every far sample searched) are refused by today's key (v3: far samples
    // bounded by a slope); the line is the last of every saved key (partFibreKey).
    const std::string v2 = have + "fibre placement = legacy\nfibre backend = gpu float32 (gpu fibres v2: none behind the fabric, no serialized root check)\n";
    CHECK(staleFibreKeyField(v2, have + "fibre placement = legacy\n" + fibreBackendLine()) == "fibre backend");
    CHECK(staleFibreKeyField(have + "fibre placement = legacy\n" + fibreBackendLine(), have + "fibre placement = legacy\n" + fibreBackendLine()).empty());
}

// knit_check --pieces with --save-fibres / --load-fibres: a manifest of no registered scene is a
// scene of spiral pieces on the common path (exported, as any scene of the user's), a second run
// loads its fibres instead of growing them and exports the same fibre files; changed settings
// refuse the saved fibres.
TEST(cli_pieces_manifest_and_saved_fibres) {
    if (knitCheckExe().empty() || !fs::exists(knitCheckExe())) SKIP("knit_check not configured");
    if (!gpu::available(nullptr)) SKIP("no usable CUDA GPU");
    const std::string d = freshDir("cli_pieces");
    fs::create_directories(d + "/own_pieces");
    writeFixturePiece(d + "/own_pieces");
    const std::string pieces = "--pieces \"" + d + "/own_pieces/manifest.txt\" --gauge piece_fixture.obj 24 --set fibres.enabled=true ";
    CHECK(runLogged(pieces + "--stage-profile --save-knit \"" + d + "/knit\" --save-fibres \"" + d + "/fibres\" --export-ppt \"" + d + "/a\"",
                    d + "/a.log") == 0);
    CHECK(fs::exists(d + "/a/scene.json") && fs::exists(d + "/fibres/piece_fixture.fibrecache"));
    // The light checks by default (as the viewer): no CPU audit of the GPU fibres, no measured angles.
    const std::string light = readText(d + "/a.log");
    CHECK(light.find("audit of the omitted root check") == std::string::npos && light.find("angles not measured") != std::string::npos);
    CHECK(runLogged(pieces + "--fiber-check --load-knit \"" + d + "/knit\" --export-ppt \"" + d + "/checked\"", d + "/checked.log") == 0);
    const std::string checked = readText(d + "/checked.log");
    CHECK(checked.find("GPU fibres audit (CPU, not applied)") != std::string::npos && checked.find("angle to the yarn at emergence") != std::string::npos);
    CHECK(runLogged(pieces + "--load-knit \"" + d + "/knit\" --load-fibres \"" + d + "/fibres\" --export-ppt \"" + d + "/b\"", d + "/b.log") == 0);
    CHECK(readText(d + "/b.log").find("from 1 saved parts") != std::string::npos);
    int compared = 0;
    for (const auto& e : fs::recursive_directory_iterator(d + "/a")) {
        if (!e.is_regular_file() || e.path().string().find("fibres") == std::string::npos) continue;
        const std::string other = d + "/b/" + fs::relative(e.path(), d + "/a").string();
        CHECK(fs::exists(other) && readText(other) == readText(e.path().string()));
        ++compared;
    }
    CHECK(compared > 0);
    const int stale = runLogged(pieces + "--set fibres.density=7 --load-knit \"" + d + "/knit\" --load-fibres \"" + d + "/fibres\"", d + "/stale.log");
    CHECK(stale == 2 && readText(d + "/stale.log").find("stale: fibre params") != std::string::npos);
    const int noGauge = runLogged("--pieces \"" + d + "/own_pieces/manifest.txt\"", d + "/nogauge.log");
    CHECK(noGauge == 2 && readText(d + "/nogauge.log").find("--gauge") != std::string::npos);
}

// The teddy's own research route was retired in S2 (its knitting is the common spiral layout, its
// fibres and export the common ones): its other fibre and export flags are refused with what
// replaces them, never ignored; --closeup needs the export its cameras belong to. A refused run
// writes nothing.
TEST(cli_retired_research_flags) {
    if (knitCheckExe().empty() || !fs::exists(knitCheckExe())) SKIP("knit_check not configured");
    const std::string d = freshDir("cli_retired_flags");
    const std::vector<std::string> gone = {"--coverage-report x", "--mesh-route", "--xatlas-charts", "--jobs 4"};
    for (size_t i = 0; i < gone.size(); ++i) {
        const std::string log = d + "/gone" + std::to_string(i) + ".log";
        const int rc = runLogged("--scene teddy " + gone[i], log);
        CHECK(rc == 2 && readText(log).find("retired") != std::string::npos);
    }
    const std::string lightLog = d + "/light.log";
    CHECK(runLogged("--strand-demo \"" + d + "/strand\" --light 0 0 1 5", lightLog) == 2 && readText(lightLog).find("retired") != std::string::npos);
    CHECK(!fs::exists(d + "/strand/scene.json"));
    const std::vector<std::string> retired = {"--piece 7",        "--fibre-trace x",      "--fibre-domain-dump x", "--fibre-domain-box 5",
                                              "--no-root-check",  "--fiber-samples 0",    "--no-attachment-check"};
    for (size_t i = 0; i < retired.size(); ++i) {
        const std::string log = d + "/r" + std::to_string(i) + ".log";
        const int rc = runLogged("--scene teddy " + retired[i] + " --export-ppt \"" + d + "/out\"", log);
        const std::string text = readText(log);
        if (rc != 2 || text.find("retired (S2)") == std::string::npos) std::printf("    %s: exit %d: %s\n", retired[i].c_str(), rc, text.c_str());
        CHECK(rc == 2 && text.find("retired (S2)") != std::string::npos);
    }
    for (const std::string& flag : {std::string("--fiber-profile"), std::string("--fibre-regions on"), std::string("--surface-diag"),
                                    std::string("--viewer-checks"), std::string("--no-gpu-audit"), std::string("--no-fibre-angles")}) {
        CHECK(runLogged("--scene teddy " + flag, d + "/fibre_flag.log") == 2);
        CHECK(readText(d + "/fibre_flag.log").find("retired") != std::string::npos);
    }
    for (const std::string& flag : {"--contract compress", "--k 0.1", "--core-c 0", "--relax-cell", "--relaxed", "--variant E", "--tau 0.9", "--cell-time 5"}) {
        CHECK(runLogged(flag, d + "/knit_flag.log") == 2);
        CHECK(readText(d + "/knit_flag.log").find("retired") != std::string::npos);
    }
    for (const std::string& flag : {std::string("--dump-relax x"), std::string("--dump-filter phase=settle"), std::string("--dump-exit"), std::string("--dump-trace")}) {
        CHECK(runLogged("--scene teddy " + flag, d + "/relax_flag.log") == 2);
        CHECK(readText(d + "/relax_flag.log").find("retired") != std::string::npos);
    }
    CHECK(runLogged("--relax-frames in out", d + "/relax_frames.log") == 2);
    CHECK(readText(d + "/relax_frames.log").find("retired") != std::string::npos);
    CHECK(runLogged("--scene teddy --closeup 7", d + "/closeup.log") == 2);
    CHECK(readText(d + "/closeup.log").find("--export-ppt") != std::string::npos);
    CHECK(!fs::exists(d + "/out"));
}

// knit.time_cap_s stops the spiral pieces as every scene's parts (the common part loop): the piece
// running at the deadline is timed out, the later ones are not started, the result is returned
// (shown) and the knit job does not cache it.
TEST(teddy_spiral_pieces_time_cap) {
    if (!gpu::available(nullptr)) SKIP("no usable CUDA GPU");
    const SceneEntry* e = findScene("teddy");
    CHECK(e);
    knit_test::requireScene(findProjectFile(e->piecesManifest));
    SceneSettings s = viewerDefaults(*e);
    s.timeCapS = 1;
    std::string err;
    const std::string cache = freshDir("teddy_cap_scenes") + "/scenes";
    auto b = run(BuiltScene::Kind::Knit, s, &err, nullptr, nullptr, nullptr, "", "ppt", false, cache);
    if (!b || !b->knit || !b->knit->scene) {
        std::printf("    knit: %s\n", err.c_str());
        CHECK(false);
        return;
    }
    int timedOut = 0, notStarted = 0;
    for (const ScenePart& r : b->knit->scene->parts) {
        CHECK(r.knit || r.stop == KnitStop::TimeCap || r.stop == KnitStop::NotStarted);
        timedOut += r.stop == KnitStop::TimeCap;
        notStarted += r.stop == KnitStop::NotStarted;
        CHECK(partIncomplete(r, "") == !r.knit);  // a piece without its knit is incomplete
    }
    std::printf("    %d pieces timed out, %d not started\n", timedOut, notStarted);
    CHECK(timedOut + notStarted > 0);
    CHECK(b->report.find("scene cache not written") != std::string::npos);
    CHECK(!fs::exists(sceneCacheEntryDir(cache, s, sceneAssetHash(*e))));
}

// The teddy (spiral pieces) takes the path of every scene. The Knit job delivers the common knit
// result (layout Spiral: one knitted part per piece, the kept GLB nodes as plain primitives, no cache
// directories); the scene cache restores it from knit.bin and verifies gate (a) with the CPU
// validator it was made with; the Fibres job grows the fibres on the restored knit; the Export job
// writes yarn, filling and fibres, the nodes the teddy's selection keeps (eyes, noses, mouths, the bow) as
// their source meshes, and the credit travels with it.
// SYNTHETIC Final: the knit is Draft (minutes less) and is marked Final by hand before the export.
TEST(teddy_spiral_pieces_result_path_synthetic_final) {
    if (!gpu::available(nullptr)) SKIP("no usable CUDA GPU");
    const SceneEntry* e = findScene("teddy");
    CHECK(e);
    Manifest man;
    std::string err;
    knit_test::requireScene(findProjectFile(e->piecesManifest));
    CHECK(readManifest(findProjectFile(e->piecesManifest), man, &err));
    SceneSettings s = viewerDefaults(*e);  // Draft: placement and gate (a) only
    s.fibres = true;
    s.segmentBudget = 200000;  // the budget order evaluates only the fibres the export keeps
    const std::string cache = freshDir("teddy_scenes") + "/scenes", out = freshDir("teddy_out");
    auto b = run(BuiltScene::Kind::Knit, s, &err, nullptr, nullptr, nullptr, "", "ppt", false, cache);
    if (!b || !b->knit || !b->knit->scene) {
        std::printf("    knit: %s\n", err.c_str());
        CHECK(false);
        return;
    }
    const Scene& t = *b->knit->scene;
    int knitted = 0, plain = 0;
    for (const ScenePart& r : t.parts) knitted += r.knit != nullptr;
    for (const SourcePrimitive& g : t.source.prims) plain += g.cls == SourcePrimitive::Class::Plain;
    CHECK(!hasChartParts(t) && t.parts.size() == man.pieces.size() && knitted == static_cast<int>(man.pieces.size()));
    CHECK(plain == static_cast<int>(e->glb.plainNodes.size()) && b->knit->parts == static_cast<int>(man.pieces.size()));
    for (const ScenePart& r : t.parts) CHECK(!r.knit || r.knit->gateAOnCpu);
    // The second knit job reads the cache: the same knit, restored and verified before the fibres.
    auto c = run(BuiltScene::Kind::Knit, s, &err, nullptr, nullptr, nullptr, "", "ppt", false, cache);
    if (!c || !c->knit || !c->knit->scene) {
        CHECK(false);
        return;
    }
    CHECK(c->cached && c->knit->computational && c->knit->digest == b->knit->digest);
    CHECK(!hasChartParts(*c->knit->scene) && c->knit->scene->restored && !c->knit->scene->verified);
    auto f = run(BuiltScene::Kind::Fibres, s, &err, nullptr, c->knit);
    if (!f || !f->fibreResult) {
        std::printf("    fibres: %s\n", err.c_str());
        CHECK(false);
        return;
    }
    CHECK(c->knit->scene->verified && f->fibreResult->knitId == c->knit->id && f->fibreResult->fibres > 0);
    // The export, of the restored Draft knit marked Final (synthetic).
    auto k = std::make_shared<KnitResult>();
    k->id = c->knit->id;
    k->settings = c->knit->settings;
    k->final = true;
    k->parts = c->knit->parts;
    k->failing = c->knit->failing;
    k->scene = c->knit->scene;
    auto x = run(BuiltScene::Kind::Export, s, &err, nullptr, k, f->fibreResult, out + "/teddy", "ppt", true);
    if (!x) std::printf("    export: %s\n", err.c_str());
    CHECK(x != nullptr);
    if (!x) return;
    // The eyes, noses and mouths of both bears and the bow (bears.glb's kept nodes), as their source meshes;
    // a filling per piece.
    const std::vector<Vec3> source = sourcePointsMetres(*c->knit->scene);
    const std::vector<std::string> plainNames = plainPrimitiveNames(c->knit->scene->source);
    int keptWritten = 0;
    for (const std::string& kept : plainNames) {
        if (kept.empty()) continue;
        const std::vector<Vec3> p = packagePlyPositions(out + "/teddy/meshes/" + kept + ".ply");
        std::printf("    %s: %zu of %zu vertices on the source\n", kept.c_str(), pointsOnSource(p, source), p.size());
        CHECK(!p.empty() && pointsOnSource(p, source) == p.size());
        ++keptWritten;
    }
    CHECK(keptWritten == static_cast<int>(e->glb.plainNodes.size()));
    int fillings = 0;
    for (const std::string& name : filesUnder(out + "/teddy/meshes"))
        fillings += name.find("_stuffing.ply") != std::string::npos || name.find("_backing.ply") != std::string::npos;
    std::printf("    %d fillings for %zu pieces\n", fillings, man.pieces.size());
    CHECK(fillings == static_cast<int>(man.pieces.size()));
    const std::string report = readText(out + "/teddy/report.json");
    CHECK(report.find("\"omitted\"") == std::string::npos && report.find("hectopod") != std::string::npos);
    CHECK(readText(out + "/teddy/CREDITS.txt").find("hectopod") != std::string::npos);
}

// Evidence (not a CTest case: knit_tests_workflow viewer_chain_evidence): a scene's viewer chain - the
// Knit (Draft), Fibres and Export jobs as the viewer issues them; hair scenes without the knit - with
// wall times, counts and the peak memory after each job, written under KNIT_EVIDENCE_DIR (scene:
// KNIT_EVIDENCE_SCENE, default man). Run in two checkouts for a before/after comparison.
TEST(viewer_chain_evidence) {
    char* env = nullptr;
    size_t envLen = 0;
    const bool set = _dupenv_s(&env, &envLen, "KNIT_EVIDENCE_DIR") == 0 && env != nullptr;
    const std::string dir = set ? env : "";
    std::free(env);
    if (!set) SKIP("KNIT_EVIDENCE_DIR not set");
    char* sceneEnv = nullptr;
    const bool named = _dupenv_s(&sceneEnv, &envLen, "KNIT_EVIDENCE_SCENE") == 0 && sceneEnv != nullptr;
    const std::string scene = named ? sceneEnv : "teddy";
    std::free(sceneEnv);
    const SceneEntry* e = findScene(scene);
    CHECK(e);
    if (!e) return;
    SceneSettings s = viewerDefaults(*e);
    s.fibres = true;
    // KNIT_EVIDENCE_PROFILE=final: the Final profile (exported as it is); else the viewer's Draft.
    char* profileEnv = nullptr;
    const bool profiled = _dupenv_s(&profileEnv, &envLen, "KNIT_EVIDENCE_PROFILE") == 0 && profileEnv != nullptr;
    const bool finalProfile = profiled && std::string(profileEnv) == "final";
    std::free(profileEnv);
    if (finalProfile) applyFullProfile(s, *e);
    // KNIT_EVIDENCE_CACHE: the scene cache root to open (another checkout's entries); else a fresh one.
    char* cacheEnv = nullptr;
    const bool cached = _dupenv_s(&cacheEnv, &envLen, "KNIT_EVIDENCE_CACHE") == 0 && cacheEnv != nullptr;
    const std::string cacheRoot = cached ? cacheEnv : dir + "/cache/scenes";
    std::free(cacheEnv);
    std::string err;
    using Clock = std::chrono::steady_clock;
    auto secondsSince = [](Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); };
    std::string memory;  // the peak working set after each job
    auto notePeak = [&memory](const char* after) {
#ifdef _WIN32
        PROCESS_MEMORY_COUNTERS pmc;
        if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof pmc))
            memory += std::string(after) + " peak " + std::to_string(static_cast<int>(pmc.PeakWorkingSetSize / 1048576.0)) + " MB, now " +
                      std::to_string(static_cast<int>(pmc.WorkingSetSize / 1048576.0)) + " MB\n";
#endif
    };
    Clock::time_point t0 = Clock::now();
    // Fur: no knit job.
    auto b = knitWorkflow(s) ? run(BuiltScene::Kind::Knit, s, &err, nullptr, nullptr, nullptr, "", "ppt", false, cacheRoot)
                             : std::make_shared<BuiltScene>();
    if (furWorkflow(s)) b->knit = std::make_shared<KnitResult>();
    const double knitS = secondsSince(t0);
    notePeak("knit");
    if (!b || !b->knit) {
        std::printf("knit failed: %s\n", err.c_str());
        CHECK(false);
        return;
    }
    t0 = Clock::now();
    auto f = run(BuiltScene::Kind::Fibres, s, &err, nullptr, knitWorkflow(s) ? b->knit : nullptr);
    const double fibreS = secondsSince(t0);
    notePeak("fibres");
    if (!f || !f->fibreResult) {
        std::printf("fibres failed: %s\n", err.c_str());
        CHECK(false);
        return;
    }
    // The export job needs a Final knit: a Final knit as it is; a Draft knit marked Final (synthetic),
    // exported as a diagnostic.
    const bool synthetic = knitWorkflow(s) && !b->knit->final;
    auto k = std::make_shared<KnitResult>();
    k->id = b->knit->id;
    k->digest = b->knit->digest;
    k->settings = b->knit->settings;
    k->final = true;
    k->parts = b->knit->parts;
    k->pass = b->knit->pass;
    k->fail = b->knit->fail;
    k->unresolved = b->knit->unresolved;
    k->plain = b->knit->plain;
    k->cross = b->knit->cross;
    k->failing = b->knit->failing;
    k->scene = b->knit->scene;
    t0 = Clock::now();
    const ExportEligibility el = exportEligibility(s, knitWorkflow(s) ? (synthetic ? k.get() : b->knit.get()) : nullptr, f->fibreResult.get(), FibreChoice::Included);
    auto x = run(BuiltScene::Kind::Export, s, &err, nullptr, knitWorkflow(s) ? (synthetic ? k : b->knit) : nullptr, f->fibreResult, dir + "/export", "ppt", true);
    const double exportS = secondsSince(t0);
    notePeak("export");
    if (!x) std::printf("export failed: %s\n", err.c_str());
    CHECK(x != nullptr);
    double peakMB = 0;
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof pmc)) peakMB = pmc.PeakWorkingSetSize / 1048576.0;
#endif
    char state[512];
    std::snprintf(state, sizeof state, "knit state: %s, %s, %s%s; export eligibility %s\n",
                  b->cached ? "from the cache" : "knitted", b->knit->computational ? "computational" : ("view only: " + b->knit->viewOnlyReason).c_str(),
                  b->knit->final ? "Final" : "Draft", synthetic ? " (exported marked Final: synthetic)" : "",
                  el.state == ExportEligibility::State::Ready ? "ready" : el.state == ExportEligibility::State::Diagnostic ? "diagnostic" : "blocked");
    char line[1024];
    std::snprintf(line, sizeof line,
                  "knit %.1f s: %d parts, %d PASS, %d FAIL, %d UNRESOLVED, %d plain, digest %016llx, %zu display meshes\n"
                  "fibres %.1f s: %lld fibres, %lld segments, %lld displayed, %zu display line vertices\n"
                  "export %.1f s: %s\npeak working set %.0f MB\n",
                  knitS, b->knit->parts, b->knit->pass, b->knit->fail, b->knit->unresolved, b->knit->plain,
                  static_cast<unsigned long long>(b->knit->digest), b->meshes.size(), fibreS, f->fibreResult->fibres,
                  f->fibreResult->segments, f->fibreResult->displayedFibres, f->fibres.pos.size() / 3, exportS,
                  x ? x->exportedTo.c_str() : err.c_str(), peakMB);
    std::printf("%s%s%s", state, line, memory.c_str());
    writeText(dir + "/summary.txt", state + std::string(line) + memory + std::string("\n--- fibre log\n") + f->fibreResult->log + "--- export log\n" +
                                        (x ? x->exportLog : err) + "--- knit report\n" + b->report);
}

int main(int argc, char** argv) {
    int failed = 0, run = 0, skipped = 0;
    const std::string only = argc > 1 ? argv[1] : "";
    if (only == "--list") {
        for (const auto& c : knit_test::registry()) std::printf("%s\n", c.name);
        return 0;
    }
    for (const auto& c : knit_test::registry()) {
        if (!only.empty() && only != c.name) continue;
        ++run;
        try {
            c.fn();
            std::printf("[PASS] %s\n", c.name);
        } catch (const knit_test::Skipped& s) {
            ++skipped;
            std::printf("[SKIP] %s: %s\n", c.name, s.why);
        } catch (const std::exception& e) {
            ++failed;
            std::printf("[FAIL] %s: %s\n", c.name, e.what());
        }
    }
    if (run == 0) {
        std::printf("no test named %s\n", only.c_str());
        return 1;
    }
    std::printf("%d run, %d failed, %d skipped\n", run, failed, skipped);
    if (failed) return 1;
    return skipped == run ? 77 : 0;  // 77: every case run was skipped (CTest SKIP_RETURN_CODE)
}

namespace {

// A textured square of two triangles at x0, its own image of w x h texels.
SourcePrimitive texturedSquare(const std::string& name, SourcePrimitive::Class cls, double x0, int w, int h) {
    SourcePrimitive g;
    g.cls = cls;
    g.mesh = name;
    g.material = name;
    g.welded.pos = {{x0, 0, 0}, {x0 + 10, 0, 0}, {x0 + 10, 10, 0}, {x0, 10, 0}};
    g.welded.tri = {{0, 1, 2}, {0, 2, 3}};
    g.welded.uv = {{Vec2{0, 0}, Vec2{1, 0}, Vec2{1, 1}}, {Vec2{0, 0}, Vec2{1, 1}, Vec2{0, 1}}};
    g.welded.computeNormals();
    g.check = checkMesh(g.welded);
    auto image = std::make_shared<ColourImage>();
    image->width = w;
    image->height = h;
    image->rgba.assign(static_cast<size_t>(w) * h * 4, 200);
    g.albedo = image;
    return g;
}

} // namespace

// The Knit job's partial scene after the second part list (gap patches): every textured mesh indexes a
// texture of that publication, its own primitive's, and the parts already knitted stay (2026-11-06: the
// texture index outlived the first onLoaded, so the second publication's meshes pointed past its textures).
TEST(knit_stand_ins_survive_the_second_part_list) {
    struct Keep : JobReporter {
        std::vector<std::shared_ptr<BuiltScene>> scenes;
        void progress(const std::string&, int, int) override {}
        void partial(std::shared_ptr<BuiltScene> scene) override { scenes.push_back(std::move(scene)); }
    } sink;
    const JobControl control{0, &sink};
    Scene t;
    t.source.prims.push_back(texturedSquare("eye", SourcePrimitive::Class::Plain, 0, 4, 4));
    t.source.prims.push_back(texturedSquare("body", SourcePrimitive::Class::Knit, 20, 8, 2));
    ScenePart body;
    body.prim = 1;
    body.component = 0;
    body.name = "body_c0";
    t.parts.push_back(std::move(body));
    PartialKnitScenes partials(control, std::make_shared<const BuiltScene>());
    const SceneSettings s;
    partials.showStandIns(s, t);
    DisplayMesh yarn;
    yarn.name = "body_c0 yarn";
    partials.showPart(0, yarn);
    ScenePart patch;
    patch.prim = 1;
    patch.component = 0;
    patch.name = "body_c0_gap0";
    patch.gapPatch = true;
    t.parts.push_back(std::move(patch));
    partials.showStandIns(s, t);
    CHECK(sink.scenes.size() == 2);
    if (sink.scenes.size() != 2) return;
    const BuiltScene& last = *sink.scenes.back();
    std::set<std::string> names;
    int textured = 0;
    for (const DisplayMesh& m : last.meshes) {
        names.insert(m.name);
        if (m.albedo < 0) continue;
        ++textured;
        CHECK(m.albedo < static_cast<int>(last.albedo.size()));
        if (m.albedo >= static_cast<int>(last.albedo.size())) continue;
        const int expectedWidth = m.name == "eye" ? 4 : 8;  // each mesh shows its own primitive's image
        CHECK(last.albedo[static_cast<size_t>(m.albedo)]->width == expectedWidth);
    }
    CHECK(textured == 2);  // the eye and the gap patch's stand-in; the knitted part shows its yarn
    CHECK(names.count("eye") && names.count("body_c0 yarn") && names.count("body_c0_gap0") && !names.count("body_c0"));
}
