// The export's tests (knit_tests_export, with tests_io_pins.cpp): the neutral scene and the package, the
// formats' materials, the export's tubes. Run one case: knit_tests_export <name>; list: --list.
#include <algorithm>
#include <chrono>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <regex>
#include <stdexcept>
#include <set>
#include <string>
#include <vector>
#include <unordered_map>
#include <thread>

#include "cpu_numerics.h"
#include "scene/appearance/baked_maps.h"
#include "knitting/tube/param_text.h"
#include "harness.h"
#include "export/collect_scene.h"
#include "fibres/fibre.h"
#include "fibres/fibre_budget.h"
#include "fibres/fibre_checks.h"
#include "fibres/fibre_generator.h"
#include "fibres/fibre_shape.h"
#include "fibres/fibre_strands.h"
#include "io/image_io.h"
#include "io/publish.h"
#include "knitting/checks/checks.h"
#include "knitting/knit/fixtures.h"
#include "knitting/knit/piece_knit.h"
#include "knitting/layout/stitch_layout.h"
#include "knitting/relax/settle.h"
#include "knitting/stitch/hand_tension.h"
#include "knitting/stitch/stockinette.h"
#include "knitting/tube/plies.h"
#include "knitting/tube/ply_profile.h"
#include "knitting/tube/yarn_surface.h"
#include "model/cancel.h"
#include "model/choices.h"
#include "model/colour.h"
#include "model/coordinates.h"
#include "model/hash.h"
#include "model/parallel.h"
#include "model/spatial_hash.h"
#include "scene/appearance/ply_maps.h"
#include "scene/source/scene_parts.h"
#include "scene_fixtures.h"
#include "settings/settings_mapping.h"
#include "strand_fixtures.h"

using namespace knit;

namespace {

// Vertex positions of a binary little-endian PLY (x y z nx ny nz [u v] float32).
std::vector<Vec3> readPlyPositions(const std::string& path) {
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

// A yarn tube as one part of the package.
bool addExportTube(ExportScene& package, const std::string& name, const std::string& material, const CurveSet& cs, const PatternParams& pp,
                   const SceneMapping& map, std::string* err, const TubeStyle& style = TubeStyle{}) {
    const std::optional<ExportTube> tube = buildExportTube(cs, pp, map, style, err);
    return tube && package.addMeshPart(name, material, tube->mesh, err);
}

std::string sceneJsonText(const std::string& dir) {
    std::ifstream f(dir + "/scene.json");
    return std::string((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
}

// No inf / nan value (a bare literal after [ , : or a space; the renderer's JSON reader rejects them).
bool sceneJsonFinite(const std::string& dir) {
    const std::string t = sceneJsonText(dir);
    static const std::regex bad(R"([\[,: ]-?(inf|nan|INF|NAN))");
    return !t.empty() && !std::regex_search(t, bad);
}

// The vertices of an exporter PLY without UVs: x y z nx ny nz per vertex.
std::vector<float> plyVertices(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    std::string line;
    size_t n = 0;
    while (std::getline(f, line) && line != "end_header")
        if (line.rfind("element vertex ", 0) == 0) n = std::stoul(line.substr(15));
    std::vector<float> v(6 * n);
    f.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(v.size() * sizeof(float)));
    return f ? v : std::vector<float>{};
}

// The share of a PLY's triangles whose winding normal (v1 - v0) x (v2 - v0) agrees with the sum of
// their vertex normals (-1 if unreadable).
double plyWindingAgreement(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    std::string line;
    size_t nv = 0, nf = 0, stride = 6;
    while (std::getline(f, line) && line != "end_header") {
        if (line.rfind("element vertex ", 0) == 0) nv = std::stoul(line.substr(15));
        if (line.rfind("element face ", 0) == 0) nf = std::stoul(line.substr(13));
        if (line == "property float u") stride = 8;
    }
    std::vector<float> v(stride * nv);
    f.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(v.size() * sizeof(float)));
    size_t agree = 0;
    for (size_t t = 0; t < nf; ++t) {
        unsigned char n = 0;
        uint32_t i[3];
        f.read(reinterpret_cast<char*>(&n), 1);
        f.read(reinterpret_cast<char*>(i), sizeof i);
        if (!f || n != 3) return -1;
        auto P = [&](uint32_t k) { return Vec3{v[stride * k], v[stride * k + 1], v[stride * k + 2]}; };
        auto N = [&](uint32_t k) { return Vec3{v[stride * k + 3], v[stride * k + 4], v[stride * k + 5]}; };
        agree += dot(cross(P(i[1]) - P(i[0]), P(i[2]) - P(i[0])), N(i[0]) + N(i[1]) + N(i[2])) > 0;
    }
    return nf ? double(agree) / double(nf) : -1;
}

} // namespace

TEST(pptstrands_refusals) {
    const std::string d = scratchDir("refusals");
    std::string err;
    {  // a one-vertex strand
        StrandSetF s = tinyStrands();
        s.add({0.f, 0.f, 0.f}, {1e-4f});
        CHECK(!writeStrandsFile(d + "/b.pptstrands", s, StrandRules{}, &err));
        CHECK(!std::filesystem::exists(d + "/b.pptstrands"));
    }
    {  // a radius below 1 um after float conversion
        StrandSetF s = tinyStrands(1e-3);
        CHECK(!writeStrandsFile(d + "/b.pptstrands", s, StrandRules{}, &err));
    }
    {  // a radius step above 50 %
        StrandSetF s;
        s.add({0, 0, 0, 0.01f, 0, 0, 0.02f, 0, 0}, {2e-4f, 2e-4f, 1e-4f});
        CHECK(!writeStrandsFile(d + "/b.pptstrands", s, StrandRules{}, &err));
    }
    {  // a segment shorter than 4 r, and one longer than 200 r
        StrandSetF s;
        s.add({0, 0, 0, 5e-4f, 0, 0}, {2e-4f, 2e-4f});
        CHECK(!writeStrandsFile(d + "/b.pptstrands", s, StrandRules{}, &err));
        StrandSetF l;
        l.add({0, 0, 0, 0.05f, 0, 0}, {2e-4f, 2e-4f});
        CHECK(!writeStrandsFile(d + "/b.pptstrands", l, StrandRules{}, &err));
    }
    {  // the renderer's curvature rule (circumradius >= r) is implied by the segment
       // rule: a triple's circumradius is at least half its longest side >= 2 r, so
       // even a hairpin of admissible segments passes, with the bound reported
        StrandSetF s;
        s.add({0, 0, 0, 0.002f, 0, 0, 0.0f, 1e-5f, 0}, {2e-4f, 2e-4f, 2e-4f});
        StrandCheck c = checkStrands(s);
        CHECK(c.ok);
        CHECK(c.minCurvOverR >= 2.0);
    }
    {  // empty: refused by the writer, omitted by the exporter
        StrandSetF s;
        CHECK(!writeStrandsFile(d + "/b.pptstrands", s, StrandRules{}, &err));
        ExportScene ex(d);
        CHECK(ex.addStrands("none", ex.orientedFibreMaterial("wool_white"), s, &err));
        CHECK(ex.strandParts().empty());
    }
}

// The late taper (StrandTaper): the body keeps the root radius, the tail narrows by a cosine to tip x root,
// the path is the generated polyline (every generated vertex kept, every new one on a chord), and the strand
// passes checkStrands. A fibre too thick for its length keeps fibresToStrands' strand and is counted.
TEST(late_taper_keeps_the_body_and_the_path) {
    const std::vector<Vec3> path = {{0, 0, -0.2}, {0, 0, 0}, {0.8, 0, 0.1}, {1.6, 0, 0.1}, {2.4, 0, 0.3}};
    FibreSet fs;
    Fibre f{};
    f.radius0 = 0.02f;
    f.radius1 = 0.01f;
    f.emergeL = 0.1f;
    fs.addFibre(f, path);
    Fibre stub{};
    stub.radius0 = 0.1f;
    stub.radius1 = 0.05f;
    fs.addFibre(stub, std::vector<Vec3>{{0, 0, -0.1}, {0, 0, 0.3}});
    const SceneMapping map;  // sigma 1: scene units are library mm, the export's axes
    const StrandTaper taper{0.1, 8, 0.1};
    const TaperedStrands t = fibresToTaperedStrands(fs, map, -1, taper);
    CHECK(t.tapered == 1 && t.untapered == 1);
    CHECK(t.strands.strands() == 2);
    CHECK(checkStrands(t.strands).ok);
    // The tapered fibre: arc length along the exported points.
    const uint32_t a = t.strands.first[0], b = t.strands.first[1];
    std::vector<Vec3> out;
    for (uint32_t v = a; v < b; ++v) out.push_back({t.strands.xyz[3 * v], t.strands.xyz[3 * v + 1], t.strands.xyz[3 * v + 2]});
    double L = 0;
    for (size_t k = 1; k < path.size(); ++k) L += norm(path[k] - path[k - 1]);
    const double visible = L - 0.1, tail = std::max(0.1 * visible, 8 * 2 * 0.02), x0 = L - tail;
    CHECK(out.size() > path.size());
    double x = 0;
    for (uint32_t v = a; v < b; ++v) {
        if (v > a) x += norm(out[v - a] - out[v - a - 1]);
        if (x < x0 - 1e-6) CHECK(std::fabs(t.strands.radius[v] - 0.02f) < 1e-7f);
    }
    CHECK(std::fabs(t.strands.radius[b - 1] - 0.002f) < 1e-6f);
    CHECK(std::fabs(x - L) < 1e-5);  // new points lie on the chords: the length is the polyline's
    for (const Vec3& p : path) {
        const Vec3 q = map.toExport(p);
        bool kept = false;
        for (const Vec3& o : out) kept = kept || norm(o - q) < 1e-6;
        CHECK(kept);
    }
    // The stub, as fibresToStrands writes it.
    FibreSet only;
    only.addFibre(stub, std::vector<Vec3>{{0, 0, -0.1}, {0, 0, 0.3}});
    const StrandSetF plain = fibresToStrands(only, map);
    CHECK(std::equal(plain.radius.begin(), plain.radius.end(), t.strands.radius.begin() + b));
    // tail 0: no late taper, fibresToStrands' strands.
    const StrandSetF none = fibresToStrands(fs, map);
    CHECK(none.radius.size() == fs.points.size());
    // knit_check --retaper: the same taper on the written strands. With the emergence at the curve start
    // (the file holds no other), the strand is the export's own taper of that fibre.
    FibreSet atStart;
    Fibre g = f;
    g.emergeL = 0.2f;  // the root segment's length: the curve start
    atStart.addFibre(g, path);
    const TaperedStrands exported = fibresToTaperedStrands(atStart, map, -1, taper);
    const TaperedStrands retapered = taperStrandSet(fibresToStrands(atStart, map), taper);
    CHECK(exported.tapered == 1 && retapered.tapered == 1 && retapered.untapered == 0);
    CHECK(exported.strands.radius.size() == retapered.strands.radius.size());
    for (size_t v = 0; v < exported.strands.radius.size() && v < retapered.strands.radius.size(); ++v)
        CHECK(std::fabs(exported.strands.radius[v] - retapered.strands.radius[v]) < 1e-6f);
    for (size_t k = 0; k < exported.strands.xyz.size() && k < retapered.strands.xyz.size(); ++k)
        CHECK(std::fabs(exported.strands.xyz[k] - retapered.strands.xyz[k]) < 1e-5f);
    // The stub keeps its strand, counted.
    CHECK(taperStrandSet(fibresToStrands(fs, map), taper).untapered == 1);
}

// The Export step's forecast counts the segments the export writes, the late taper's vertices included,
// without building the strands: taperedSegmentCount equals the strands of fibresToTaperedStrands for
// a tapered fibre, a stub that keeps its strand and a stray, below any key.
TEST(late_taper_segment_count_exact) {
    const std::vector<Vec3> path = {{0, 0, -0.2}, {0, 0, 0}, {0.8, 0, 0.1}, {1.6, 0, 0.1}, {2.4, 0, 0.3}};
    FibreSet fs;
    Fibre f{};
    f.radius0 = 0.02f;
    f.radius1 = 0.01f;
    f.emergeL = 0.1f;
    f.key = 0.2;
    fs.addFibre(f, path);
    Fibre stub{};
    stub.radius0 = 0.1f;
    stub.radius1 = 0.05f;
    stub.key = 0.5;
    fs.addFibre(stub, std::vector<Vec3>{{0, 0, -0.1}, {0, 0, 0.3}});
    Fibre stray = f;
    stray.population = FibrePopulation::Random;
    stray.key = 0.8;
    std::vector<Vec3> longer = path;
    longer.push_back({3.2, 0, 0.5});
    fs.addFibre(stray, longer);
    const SceneMapping map;
    const auto segments = [](const StrandSetF& s) { return static_cast<long long>(s.vertices() - s.strands()); };
    long long grown = 0;
    for (const Fibre& g : fs.fibres) grown += g.segments();
    for (double tip : {0.1, 0.3}) {
        const StrandTaper taper{0.2, 8, tip};
        const long long all = taperedSegmentCount(fs, map, taper, 2.0f);
        CHECK(all == segments(fibresToTaperedStrands(fs, map, -1, taper).strands));
        CHECK(all > grown);  // the tail's vertices: counting the grown segments falls short
        FibreSet below;  // the fibres with key below 0.6: the tapered one and the stub
        for (size_t i = 0; i < 2; ++i) below.addFibre(fs.fibres[i], fs.pointsOf(fs.fibres[i]));
        CHECK(taperedSegmentCount(fs, map, taper, 0.6f) == segments(fibresToTaperedStrands(below, map, -1, taper).strands));
        std::printf("    tip %.1f: %lld segments exported for %lld grown\n", tip, all, grown);
        for (TaperPoints points : {TaperPoints::Existing, TaperPoints::Onset, TaperPoints::Refined}) {
            const StrandTaper sparse{0.2, 8, tip, points, 0.05};  // the forecast of every mode is the export's
            CHECK(taperedSegmentCount(fs, map, sparse, 2.0f) == segments(fibresToTaperedStrands(fs, map, -1, sparse).strands));
            CHECK(taperedSegmentCount(fs, map, sparse, 0.6f) == segments(fibresToTaperedStrands(below, map, -1, sparse).strands));
            CHECK(measureTaper(fs, map, sparse, 2.0f).finalSegments == taperedSegmentCount(fs, map, sparse, 2.0f));
        }
    }
    CHECK(taperedSegmentCount(fs, map, StrandTaper{0, 8, 0.1}, 2.0f) == grown);  // no late taper
}

// --- The sparse late taper (TaperPoints Existing, Onset, Refined)

// The intended profile, written out again here so the tests do not share the code under test.
double intendedTaperRadius(double s, double x0, double L, double r0, double tip) {
    if (s <= x0) return r0;
    const double w = std::min(1.0, (s - x0) / (L - x0));
    return r0 * (tip + (1 - tip) * 0.5 * (1 + std::cos(kPi * w)));
}

// One fibre on a straight line along x: `xs` its vertices (xs[0] the buried root), radius r0.
FibreSet straightFibre(const std::vector<double>& xs, double r0, double emergeL) {
    std::vector<Vec3> path;
    for (double x : xs) path.push_back({x, 0, 0});
    FibreSet fs;
    Fibre f{};
    f.radius0 = float(r0);
    f.radius1 = float(r0 / 2);
    f.emergeL = float(emergeL);
    fs.addFibre(f, path);
    return fs;
}

// The written strand k: arc lengths from its stored positions, and its stored radii.
struct WrittenStrand {
    std::vector<double> s, r;
    std::vector<Vec3> p;
};
WrittenStrand writtenStrand(const StrandSetF& st, size_t k) {
    WrittenStrand w;
    for (uint32_t v = st.first[k]; v < st.first[k + 1]; ++v) {
        const Vec3 p{st.xyz[3 * v], st.xyz[3 * v + 1], st.xyz[3 * v + 2]};
        w.s.push_back(w.p.empty() ? 0.0 : w.s.back() + norm(p - w.p.back()));
        w.p.push_back(p);
        w.r.push_back(st.radius[v]);
    }
    return w;
}

// The largest |r(s) - r_linear(s)| of a written strand, sampled 10000 times per interval.
double sampledTaperError(const WrittenStrand& w, double x0, double L, double r0, double tip) {
    double worst = 0;
    for (size_t v = 0; v + 1 < w.s.size(); ++v)
        for (int i = 0; i <= 10000; ++i) {
            const double t = i / 10000.0, s = w.s[v] + (w.s[v + 1] - w.s[v]) * t;
            worst = std::max(worst, std::fabs(intendedTaperRadius(s, x0, L, r0, tip) - (w.r[v] + (w.r[v + 1] - w.r[v]) * t)));
        }
    return worst;
}

// r_linear at arc length x of a written strand.
double writtenRadiusAt(const WrittenStrand& w, double x) {
    for (size_t v = 0; v + 1 < w.s.size(); ++v)
        if (x <= w.s[v + 1]) return w.r[v] + (w.r[v + 1] - w.r[v]) * (x - w.s[v]) / (w.s[v + 1] - w.s[v]);
    return w.r.back();
}

// The tail of a straightFibre from its own float32 values, as the export computes it (16 r0 the shortest).
struct FibreTail {
    double L = 0, x0 = 0, r0 = 0;
};
FibreTail tailOfFibre(const FibreSet& one, double tailFraction) {
    const std::span<const Vec3f> points = one.pointsOf(one.fibres[0]);
    FibreTail t;
    for (size_t k = 1; k < points.size(); ++k) t.L += norm(toVec3(points[k]) - toVec3(points[k - 1]));
    t.r0 = one.fibres[0].radius0;
    t.x0 = t.L - std::max(tailFraction * (t.L - double(one.fibres[0].emergeL)), 16 * t.r0);
    return t;
}

bool holdsVertex(const WrittenStrand& w, const Vec3& q) {
    for (const Vec3& p : w.p)
        if (norm(p - q) < 1e-6) return true;
    return false;
}

constexpr TaperPoints kAllTaperPoints[] = {TaperPoints::Dense, TaperPoints::Existing, TaperPoints::Onset, TaperPoints::Refined};

// Every mode keeps every generated vertex, puts new ones on the chords, ends at tip x r0 and gives the
// same strand twice.
TEST(late_taper_modes_keep_the_path) {
    const std::vector<Vec3> path = {{0, 0, -0.2}, {0, 0, 0}, {0.8, 0, 0.1}, {1.6, 0, 0.1}, {2.4, 0, 0.3}};
    const SceneMapping map;
    for (double tip : {0.3, 0.1}) {
        const double r0 = tip < 0.2 ? 0.05 : 0.02;  // a tip of 0.1 x 0.02 makes a 0.8 chord longer than 200 r
        FibreSet fs;
        Fibre f{};
        f.radius0 = float(r0);
        f.radius1 = float(r0 / 2);
        f.emergeL = 0.1f;
        fs.addFibre(f, path);
        for (TaperPoints points : kAllTaperPoints) {
            const StrandTaper taper{0.2, 8, tip, points, 0.05};
            const TaperedStrands t = fibresToTaperedStrands(fs, map, -1, taper);
            CHECK(t.tapered == 1 && t.untapered == 0);
            const WrittenStrand w = writtenStrand(t.strands, 0);
            for (const Vec3& p : path) CHECK(holdsVertex(w, map.toExport(p)));
            double L = 0;
            for (size_t k = 1; k < path.size(); ++k) L += norm(path[k] - path[k - 1]);
            CHECK(std::fabs(w.s.back() - L) < 1e-5);  // new vertices on the chords: the length is the polyline's
            CHECK(w.r.back() == float(tip * r0));
            const TaperedStrands again = fibresToTaperedStrands(fs, map, -1, taper);
            CHECK(again.strands.xyz == t.strands.xyz && again.strands.radius == t.strands.radius);
            if (points == TaperPoints::Existing) CHECK(w.p.size() == path.size());
        }
    }
}

// Existing (B0) writes the generated vertices with the profile's radii there, so the cone of the
// segment holding x0 already thins before x0. The report measures that thinning as it is.
TEST(late_taper_existing_is_b0) {
    const double r0 = 0.02, tip = 0.3;
    const FibreSet fs = straightFibre({-0.2, 0, 0.8, 1.6, 2.4}, r0, 0.2);
    const StrandTaper taper{0.2, 8, tip, TaperPoints::Existing, 0.05};
    const double L = 2.6, visible = 2.4, x0 = L - std::max(0.2 * visible, 16 * r0);  // x0 = 2.12, inside [1.8, 2.6]
    const TaperedStrands t = fibresToTaperedStrands(fs, SceneMapping{}, -1, taper);
    const WrittenStrand w = writtenStrand(t.strands, 0);
    CHECK(w.s.size() == 5);
    for (size_t v = 0; v < w.s.size(); ++v) CHECK(std::fabs(w.r[v] - intendedTaperRadius(w.s[v], x0, L, r0, tip)) < 1e-7 * r0);
    const double thinning = (r0 - writtenRadiusAt(w, x0)) / r0;  // the cone from (1.8, r0) to (2.6, tip r0) at x0
    CHECK(std::fabs(thinning - (1 - tip) * (x0 - 1.8) / 0.8) < 1e-5);
    const TaperReport report = measureTaper(fs, SceneMapping{}, taper, 2.0f);
    CHECK(std::fabs(report.maxEarlyThinningOverR0 - thinning) < 1e-5);
    CHECK(report.finalSegments == report.originalSegments);
}

// Onset (B1): the body keeps r0 up to x0. A vertex too close to x0 moves the new vertex away by the
// shortest legal segment; a segment too short for any legal place gets none. No vertex moves.
TEST(late_taper_onset_keeps_the_body) {
    const double r0 = 0.02, tip = 0.3, L = 2.6, visible = 2.4, slack = 1.01;
    struct Case {
        const char* name;
        std::vector<double> xs;
        double x0;
    };
    const Case cases[] = {
        {"exact", {-0.2, 0, 0.8, 1.6, 2.4}, 2.12},
        {"vertex just before x0", {-0.2, 0, 0.8, 1.6, 2.4}, 1.83},
        {"vertex just after x0", {-0.2, 0, 0.8, 1.6, 2.4}, 1.77},
        {"segment too short", {-0.2, 0, 0.8, 1.6, 1.7, 2.4}, 1.85},
    };
    for (const Case& c : cases) {
        const FibreSet fs = straightFibre(c.xs, r0, 0.2);
        const StrandTaper taper{(L - c.x0) / visible, 8, tip, TaperPoints::Onset, 0.05};
        const TaperedStrands t = fibresToTaperedStrands(fs, SceneMapping{}, -1, taper);
        CHECK(t.tapered == 1);
        const WrittenStrand w = writtenStrand(t.strands, 0);
        for (double x : c.xs) CHECK(holdsVertex(w, SceneMapping{}.toExport({x, 0, 0})));
        const TaperReport report = measureTaper(fs, SceneMapping{}, taper, 2.0f);
        const double r = intendedTaperRadius(1.8, c.x0, L, r0, tip);  // at the polyline vertex 1.8
        if (std::string(c.name) == "exact") {
            CHECK(report.onsetExact == 1);
            bool atX0 = false;
            for (double s : w.s) atX0 = atX0 || std::fabs(s - c.x0) < 1e-5;
            CHECK(atX0);
            for (size_t v = 0; v < w.s.size(); ++v)
                if (w.s[v] <= c.x0 + 1e-6) CHECK(w.r[v] == float(r0));
            CHECK(report.maxEarlyThinningOverR0 < 1e-6);
        } else if (std::string(c.name) == "vertex just before x0") {
            CHECK(report.onsetShifted == 1);
            const double s = 1.8 + slack * 4 * r0, rs = intendedTaperRadius(s, c.x0, L, r0, tip);  // after x0
            const double expected = (r0 - (r0 + (rs - r0) * (c.x0 - 1.8) / (s - 1.8))) / r0;
            CHECK(std::fabs(report.maxEarlyThinningOverR0 - expected) < 1e-5);
            CHECK(s - c.x0 <= slack * 4 * r0);
        } else if (std::string(c.name) == "vertex just after x0") {
            CHECK(report.onsetShifted == 1);
            const double s = 1.8 - slack * 4 * std::max(r, r0);  // before x0: radius r0
            const double expected = (r0 - (r0 + (r - r0) * (c.x0 - s) / (1.8 - s))) / r0;
            CHECK(std::fabs(report.maxEarlyThinningOverR0 - expected) < 1e-5);
            for (size_t v = 0; v < w.s.size(); ++v)
                if (w.s[v] <= s + 1e-6) CHECK(w.r[v] == float(r0));
        } else {
            CHECK(report.onsetImpossible == 1);
            CHECK(w.p.size() == c.xs.size());  // no new vertex: the generated ones with the profile's radii
        }
    }
}

// x0 = L - tail with visible = L - emergeL; no visible length or no tail: the generated strand.
TEST(late_taper_onset_relative_to_emergence) {
    const double r0 = 0.02, tip = 0.3, L = 2.6;
    for (double emergeL : {0.2, 0.6}) {  // x0 = 2.0 and 2.1: both legal, so both exact
        const FibreSet fs = straightFibre({-0.2, 0, 0.8, 1.6, 2.4}, r0, emergeL);
        const StrandTaper taper{0.25, 8, tip, TaperPoints::Onset, 0.05};
        const double visible = L - emergeL, x0 = L - std::max(0.25 * visible, 16 * r0);
        const WrittenStrand w = writtenStrand(fibresToTaperedStrands(fs, SceneMapping{}, -1, taper).strands, 0);
        bool atX0 = false;
        for (double s : w.s) atX0 = atX0 || std::fabs(s - x0) < 1e-5;
        CHECK(atX0);
    }
    const FibreSet hidden = straightFibre({-0.2, 0, 0.8, 1.6, 2.4}, r0, 3.0);  // emerges past its tip: nothing visible
    CHECK(measureTaper(hidden, SceneMapping{}, StrandTaper{0.2, 8, tip, TaperPoints::Onset, 0.05}, 2.0f).fallbackNoTail == 1);
    const FibreSet fs = straightFibre({-0.2, 0, 0.8, 1.6, 2.4}, r0, 0.2);
    CHECK(measureTaper(fs, SceneMapping{}, StrandTaper{0, 0, tip, TaperPoints::Onset, 0.05}, 2.0f).fallbackNoTail == 1);
}

// The error is the exact maximum over each interval. The single cone of Onset is the planted case:
// its error is antisymmetric, zero at the middle and largest at w = asin(2/pi)/pi, so a mid-point
// check would report nothing. Refined meets its tolerance or names the strand nonconforming.
TEST(late_taper_refined_error_bound) {
    const double tip = 0.3;
    const FibreSet fs = straightFibre({-0.2, 0, 0.8, 1.6, 2.4}, 0.02, 0.2);
    // The tail from the fibre's own float32 values, as the export computes it: a tail start from the
    // decimal literals sits about 1e-7 away and moves the error by more than the 1e-9 compared here.
    const FibreTail t = tailOfFibre(fs, 0.2);
    const double L = t.L, x0 = t.x0, r0 = t.r0;
    const StrandTaper onset{0.2, 8, tip, TaperPoints::Onset, 0.05};
    const TaperReport cone = measureTaper(fs, SceneMapping{}, onset, 2.0f);
    const double w = std::asin(2 / kPi) / kPi, peak = (1 - tip) * (0.5 * (1 + std::cos(kPi * w)) - (1 - w));
    CHECK(std::fabs(cone.maxErrorOverR0 - peak) < 1e-4);  // 0.0737 at tip 0.3
    const WrittenStrand wc = writtenStrand(fibresToTaperedStrands(fs, SceneMapping{}, -1, onset).strands, 0);
    CHECK(sampledTaperError(wc, x0, L, r0, tip) <= cone.maxErrorOverR0 * r0 + 1e-9);
    for (double tolerance : {0.10, 0.05, 0.02}) {
        const StrandTaper refined{0.2, 8, tip, TaperPoints::Refined, tolerance};
        const TaperReport report = measureTaper(fs, SceneMapping{}, refined, 2.0f);
        const WrittenStrand wr = writtenStrand(fibresToTaperedStrands(fs, SceneMapping{}, -1, refined).strands, 0);
        const double sampled = sampledTaperError(wr, x0, L, r0, tip);
        CHECK(sampled <= report.maxErrorOverR0 * r0 + 1e-9);
        CHECK(sampled >= 0.99 * report.maxErrorOverR0 * r0);
        CHECK(report.nonconforming == 0 && report.maxErrorOverR0 <= tolerance);
        std::printf("    tolerance %.2f: %zu vertices, max error %.4f r0\n", tolerance, wr.s.size(), report.maxErrorOverR0);
    }
    // The shortest tail (16 r0) leaves no legal room for 0.5 %: nonconforming, and named.
    const FibreSet thick = straightFibre({-0.5, 0, 1.0, 2.0, 3.0}, 0.1, 0.5);  // every chord >= 4 r0
    const TaperReport tight = measureTaper(thick, SceneMapping{}, StrandTaper{0.01, 8, tip, TaperPoints::Refined, 0.005}, 2.0f);
    CHECK(tight.nonconforming == 1 && tight.nonconformingFirst.size() == 1 && tight.nonconformingFirst[0] == 0);
    CHECK(tight.maxErrorOverR0 > 0.005);
    // The report reads the stored float32 radii, not the profile in double. At tip 0.9 the tip radius's
    // float spacing is about 1e-7 r0; for a root radius whose tip radius rounds by more than 2e-8 r0, the
    // cone's largest error moves by at least 0.22 of that (4.4e-9 r0), above the 2e-9 compared here.
    double r0f = 0;
    for (int i = 0; i < 1000 && r0f == 0; ++i) {
        const float r = float(0.02 + i * 1e-7);
        if (std::fabs(double(float(0.9 * double(r))) - 0.9 * double(r)) / double(r) > 2e-8) r0f = r;
    }
    CHECK(r0f > 0);
    const FibreSet planted = straightFibre({-0.2, 0, 0.8, 1.6, 2.4}, r0f, 0.2);
    const StrandTaper nearlyFlat{0.2, 8, 0.9, TaperPoints::Onset, 0.05};
    const FibreTail tp = tailOfFibre(planted, 0.2);
    const double Lp = tp.L, x0p = tp.x0, r0p = tp.r0;
    const TaperReport stored = measureTaper(planted, SceneMapping{}, nearlyFlat, 2.0f);
    const WrittenStrand wp = writtenStrand(fibresToTaperedStrands(planted, SceneMapping{}, -1, nearlyFlat).strands, 0);
    CHECK(std::fabs(sampledTaperError(wp, x0p, Lp, r0p, 0.9) / r0p - stored.maxErrorOverR0) < 2e-9);
}

// Two radii rounded to float32 change their ratio by about 2u; taperRules allows 4u above 1 / tip.
TEST(late_taper_float_ratio_margin) {
    for (double tip : {0.3, 0.1}) {
        StrandTaper taper{0.2, 8, tip, TaperPoints::Onset, 0.05};
        StrandRules noMargin = taperRules(taper);
        noMargin.maxRadiusStep = 1 / tip - 1;
        long long above = 0, below = 0, refusedWithout = 0, refusedWith = 0;
        for (int i = 0; i <= 5500; ++i) {
            const double r0 = 2e-5 + i * 1e-8;  // tip x r0 stays above the 1e-6 floor
            const float fr0 = float(r0), ftip = float(r0 * tip);
            const double ratio = double(fr0) / double(ftip);
            above += ratio > 1 / tip;
            below += ratio < 1 / tip;
            StrandSetF s;
            s.add({0, 0, 0, float(15 * r0), 0, 0}, {fr0, ftip});
            refusedWith += !checkStrands(s, taperRules(taper)).ok;
            refusedWithout += !checkStrands(s, noMargin).ok;
        }
        CHECK(above > 0 && below > 0);
        CHECK(refusedWith == 0);
        CHECK(refusedWithout > 0);  // the defect the margin removes
        std::printf("    tip %.1f: %lld ratios round above 1/tip, %lld below; %lld refused without the margin\n", tip, above,
                    below, refusedWithout);
        StrandSetF far;
        far.add({0, 0, 0, float(15 * 2e-5), 0, 0}, {float(2e-5), float(2e-5 * tip / (1 + 1e-5))});
        CHECK(!checkStrands(far, taperRules(taper)).ok);
    }
    constexpr double r0 = 11e-6;  // the owner's reproduction: 3.333333471 stored
    constexpr double stored = double(float(r0)) / double(float(r0 * 0.3));
    static_assert(stored > 1 / 0.3);
    StrandSetF s;
    s.add({0, 0, 0, float(15 * r0), 0, 0}, {float(r0), float(r0 * 0.3)});
    CHECK(checkStrands(s, taperRules(StrandTaper{0.2, 8, 0.3, TaperPoints::Onset, 0.05})).ok);
}

// The sparse rules are scoped to the sparse modes; the default writer still refuses a 3.33x step.
TEST(late_taper_sparse_rules_scoped) {
    const StrandRules dense = taperRules(StrandTaper{0.2, 8, 0.3, TaperPoints::Dense, 0.05});
    CHECK(dense.maxRadiusStep == StrandRules{}.maxRadiusStep && dense.minSegOverRMax == 0);
    const StrandRules sparse = taperRules(StrandTaper{0.2, 8, 0.3, TaperPoints::Existing, 0.05});
    StrandSetF step;
    step.add({0, 0, 0, 0.004f, 0, 0}, {2e-4f, 6e-5f});
    CHECK(!checkStrands(step).ok && checkStrands(step).fault == StrandFault::RadiusStep);
    CHECK(checkStrands(step, sparse).ok);
    StrandSetF shortForLarger;  // 5 x the smaller radius, 1.5 x the larger: only minSegOverRMax (4) refuses it
    shortForLarger.add({0, 0, 0, 3e-4f, 0, 0}, {2e-4f, 6e-5f});
    const StrandCheck c = checkStrands(shortForLarger, sparse);
    CHECK(!c.ok && c.fault == StrandFault::SegmentShort);
    CHECK(taperRules(StrandTaper{0.2, 8, 0.9, TaperPoints::Onset, 0.05}).maxRadiusStep == 0.5);  // never below the default
    std::string err;
    const std::string d = scratchDir("sparse_rules");
    CHECK(!writeStrandsFile(d + "/a.pptstrands", step, StrandRules{}, &err));
    CHECK(writeStrandsFile(d + "/b.pptstrands", step, sparse, &err));
}

// A fibre whose sparse taper breaks a rule keeps the generated strand, bit for bit, counted by reason.
TEST(late_taper_fallback_profile) {
    const FibreSet fs = straightFibre({-0.2, 0, 0.8, 1.6, 2.4}, 0.005, 0.2);  // tip radius 0.0015: the last 0.8 chord is 533 r
    const StrandTaper taper{0.2, 8, 0.3, TaperPoints::Existing, 0.05};
    const TaperedStrands t = fibresToTaperedStrands(fs, SceneMapping{}, -1, taper);
    CHECK(t.untapered == 1);
    const StrandSetF plain = fibresToStrands(fs, SceneMapping{});
    CHECK(t.strands.xyz == plain.xyz && t.strands.radius == plain.radius);
    const TaperReport report = measureTaper(fs, SceneMapping{}, taper, 2.0f);
    CHECK(report.fallbackRule[static_cast<int>(StrandFault::SegmentLong)] == 1);
    const FibreSet thin = straightFibre({-2e-4, 0, 8e-4, 1.6e-3, 2.4e-3}, 1.5e-6, 2e-4);  // tip 4.5e-7: under the floor
    CHECK(measureTaper(thin, SceneMapping{}, taper, 2.0f).fallbackRule[static_cast<int>(StrandFault::RadiusFloor)] == 1);
}

// The close-up cameras (closeup, mid, medium, far) are written by default: they look at the knitted
// part of the largest mesh area unless a part is asked for; with nothing knitted there are none.
TEST(closeup_part_default_is_the_largest_knitted) {
    const auto triangle = [](double size) {
        auto m = std::make_shared<TriMesh>();
        m->pos = {{0, 0, 0}, {size, 0, 0}, {0, size, 0}};
        m->tri = {{0, 1, 2}};
        return std::shared_ptr<const TriMesh>(m);
    };
    std::vector<ScenePart> parts(3);
    parts[0].mesh = triangle(1);
    parts[1].mesh = triangle(3);
    parts[2].mesh = triangle(5);  // the largest, not knitted
    for (int i : {0, 1}) {
        parts[static_cast<size_t>(i)].knit = std::make_unique<PartKnit>();
        parts[static_cast<size_t>(i)].knit->t.prog.stitches.push_back(StitchRecord{});
    }
    CHECK(closeupPartOf(parts, -1) == 1);
    CHECK(closeupPartOf(parts, 2) == 2);  // asked for: any part
    for (ScenePart& p : parts) p.knit.reset();
    CHECK(closeupPartOf(parts, -1) == -1);
}

// The renderer's GPU memory forecast reproduces photon_path_tracer's measured process peaks of the
// white teddy's packages (2026-11-05, 1600 x 1000) within 1 %, and the cap is the owner's 9 GiB.
TEST(renderer_memory_matches_the_teddy) {
    struct Measured {
        long long segments, triangles;
        double mib;
    };
    const Measured runs[] = {{0, 4023450, 1512}, {0, 7886506, 1758}, {33970516, 4023450, 7832}, {49010279, 4023450, 10636}};
    for (const Measured& m : runs) {
        const double forecast = rendererMemoryGiB(m.segments, m.triangles) * 1024;
        std::printf("    %lld segments, %lld triangles: forecast %.0f MiB, measured %.0f MiB\n", m.segments, m.triangles, forecast, m.mib);
        CHECK(std::fabs(forecast - m.mib) < 0.01 * m.mib);
    }
    static_assert(kRendererMemoryCapGiB == 9.0);
}

// Exported tubes wind counter-clockwise about their vertex normals (buildTubeMesh's radial, outward
// normals; its strips wind the other way, so the renderer's import warned "winding opposes normals"
// on every yarn tube), on the CPU tube path; the ground and a plain mesh as they were.
TEST(ppt_tube_winding_outward) {
    ExportScene ex(scratchDir("ppt_winding"));
    const StrandDemo sd = makeStrandDemo(6.0, 1.12);
    SceneMapping map;
    map.sigma = 1000;
    std::string err;
    CHECK(addExportTube(ex, "strand", "wool_white", sd.curves, sd.pp, map, &err));
    TriMesh m;
    m.pos = {{-0.02, -0.02, -0.003}, {0.03, -0.02, -0.003}, {0.03, 0.03, -0.003}};
    m.tri = {{0, 1, 2}};
    CHECK(ex.addMesh("plate", "black", m, &err));
    CHECK(writeSceneJson(ex, &err));
    const double tube = plyWindingAgreement(ex.folder() + "/meshes/strand.ply"), plate = plyWindingAgreement(ex.folder() + "/meshes/plate.ply"),
                 ground = plyWindingAgreement(ex.folder() + "/meshes/_ground.ply");
    std::printf("  winding agrees with the normals: tube %.4f, plate %.4f, ground %.4f\n", tube, plate, ground);
    CHECK(tube == 1.0 && plate == 1.0 && ground == 1.0);
}

// The exported scene frames every part except lights (meshes and strands too, not only
// tubes: an export without tubes once wrote an inf camera), puts the ground quad just under
// the lowest part point, and refuses what it cannot frame.
TEST(ppt_export_bounds_ground) {
    TriMesh m;  // library z = -0.003 -> scene y = -0.003
    m.pos = {{-0.02, -0.02, -0.003}, {0.03, -0.02, -0.003}, {0.03, 0.03, -0.003}, {-0.02, 0.03, -0.003}};
    m.tri = {{0, 1, 2}, {0, 2, 3}};
    std::string err;
    auto checkGround = [](const ExportScene& ex) {
        const std::vector<float> v = plyVertices(ex.folder() + "/meshes/_ground.ply");
        CHECK(v.size() == 24);
        if (v.size() != 24) return;
        const double extent = norm(ex.bounds().hi - ex.bounds().lo), cx = 0.5 * (ex.bounds().lo.x + ex.bounds().hi.x);
        for (int i = 0; i < 4; ++i) {
            CHECK(v[6 * i + 1] == float(ex.bounds().lo.y - 1e-4 * extent));
            CHECK(v[6 * i + 3] == 0.f && v[6 * i + 4] == 1.f && v[6 * i + 5] == 0.f);
            const double h = ExportScene::kGroundHalfSize * extent;
            CHECK(std::fabs(std::fabs(v[6 * i] - cx) - h) < 1e-6 * h);
        }
        const std::string t = sceneJsonText(ex.folder());
        CHECK(t.find("\"path\": \"meshes/_ground.ply\", \"material\": \"_ground\"") != std::string::npos);
        CHECK(t.find("{\"id\": \"_ground\"") != std::string::npos);
        CHECK(t.find("\"albedo\": [0.500, 0.500, 0.500]") != std::string::npos);
    };
    {  // a mesh only (no tubes)
        ExportScene ex(scratchDir("ppt_ground_mesh"));
        CHECK(ex.addMesh("plate", "wool_white", m, &err));
        CHECK(ex.bounds().lo.x == -0.02 && ex.bounds().lo.y == -0.003 && ex.bounds().lo.z == -0.03 && ex.bounds().hi.x == 0.03 && ex.bounds().hi.y == -0.003 && ex.bounds().hi.z == 0.02);
        CHECK(writeSceneJson(ex, &err));
        CHECK(sceneJsonFinite(ex.folder()));
        checkGround(ex);
    }
    {  // a strand below the mesh lowers the ground
        ExportScene ex(scratchDir("ppt_ground_strand"));
        CHECK(ex.addMesh("plate", "wool_white", m, &err));
        StrandSetF s;
        s.add({0, -0.01f, 0, 0.004f, -0.01f, 0}, {2e-4f, 2e-4f});
        CHECK(ex.addStrands("fibres", "wool_white", s, &err));
        CHECK(ex.bounds().lo.y == double(-0.01f));
        CHECK(writeSceneJson(ex, &err));
        CHECK(sceneJsonFinite(ex.folder()));
        checkGround(ex);
    }
    {  // ground off: neither the mesh nor the material
        ExportScene ex(scratchDir("ppt_ground_off"));
        ex.groundPlane = false;
        CHECK(ex.addMesh("plate", "wool_white", m, &err));
        CHECK(writeSceneJson(ex, &err));
        CHECK(sceneJsonText(ex.folder()).find("_ground") == std::string::npos);
        CHECK(!std::filesystem::exists(ex.folder() + "/meshes/_ground.ply"));
    }
    {  // nothing to frame: refused, nothing written
        ExportScene ex(scratchDir("ppt_ground_empty"));
        err.clear();
        CHECK(!writeSceneJson(ex, &err) && !err.empty());
        CHECK(!std::filesystem::exists(ex.folder() + "/scene.json"));
    }
    {  // a non-finite extra camera: refused
        ExportScene ex(scratchDir("ppt_ground_badcam"));
        CHECK(ex.addMesh("plate", "wool_white", m, &err));
        ex.extraCameras.push_back({"bad", {0, HUGE_VAL, 0}, {0, 0, 0}, 20});
        err.clear();
        CHECK(!writeSceneJson(ex, &err) && !err.empty());
        CHECK(!std::filesystem::exists(ex.folder() + "/scene.json"));
    }
}

// Tiny scene (one mesh + one strand part, fibre material) loaded and rendered by
// the actual renderer. Skipped when ppt.exe is not configured.
TEST(ppt_load_tiny) {
#ifndef KNIT_PPT_EXE
    SKIP("KNIT_PPT_EXE not configured");
#else
    if (!std::filesystem::exists(KNIT_PPT_EXE)) SKIP("the configured ppt.exe is missing");
    const std::string d = scratchDir("ppt_tiny");
    ExportScene ex(d);
    ex.fibreLook.hairMaterials = true;  // the hair-BSDF option (fibre_oriented_wool_brown as a hair BSDF)
    TriMesh m;
    m.pos = {{-0.02, -0.02, -0.003}, {0.03, -0.02, -0.003}, {0.03, 0.03, -0.003}, {-0.02, 0.03, -0.003}};
    m.tri = {{0, 1, 2}, {0, 2, 3}};
    std::string err;
    CHECK(ex.addMesh("plate", "wool_white", m, &err));
    CHECK(ex.addStrands("fibres", ex.orientedFibreMaterial("wool_brown"), tinyStrands(), &err));
    CHECK(writeSceneJson(ex, &err));  // no tubes: the camera frames the mesh and the strands (the shiba export)
    CHECK(sceneJsonFinite(d));
    const std::string log = d + "/render.log";
    const std::string cmd = std::string("\"\"") + KNIT_PPT_EXE + "\" render \"" + d + "/scene.json\" --output \"" + d +
                            "/out.png\" --width 64 --height 64 --spp 1 > \"" + log + "\" 2>&1\"";
    const int rc = std::system(cmd.c_str());
    std::ifstream f(log);
    std::string text((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    if (rc != 0) std::printf("%s\n", text.c_str());
    CHECK(rc == 0);
    CHECK(std::filesystem::exists(d + "/out.png"));
    for (const char* bad : {"strand.", "instance.", "error", "Error", "refus"})
        if (text.find(bad) != std::string::npos) {
            std::printf("%s\n", text.c_str());
            CHECK(!"renderer log mentions a strand/instance problem or an error");
        }
    // S5: a textured tube (normal map + cavity PNGs; then the float EXR height map as a
    // bump) beside untextured fibre parts loads and renders.
    for (PlyModeOptions::Shading shading : {PlyModeOptions::Shading::Normal, PlyModeOptions::Shading::Bump}) {
        const bool bump = shading == PlyModeOptions::Shading::Bump;
        const std::string dt = scratchDir(bump ? "ppt_tiny_bump" : "ppt_tiny_normal");
        ExportScene et(dt);
        PlyModeOptions o;
        o.shading = shading;
        o.cavity = true;
        o.spec.size = 64;
        const std::optional<PlyShading> shading2 = setupPlyTextures(o, et.folder(), 68.103, &err);
        CHECK(shading2);
        et.plyShading = *shading2;
        TubeStyle style;
        style.uvRepeat = shading2->uvRepeat;
        style.uvWrap = shading2->uvWrap;
        const StrandDemo sd = makeStrandDemo(6.0, 1.12);
        SceneMapping map;
        map.sigma = 68.103;
        CHECK(addExportTube(et, "strand", ExportScene::yarnMaterial("wool_white"), sd.curves, sd.pp, map, &err, style));
        CHECK(et.addStrands("fibres", et.orientedFibreMaterial("wool_white"), tinyStrands(0.1), &err));
        CHECK(writeSceneJson(et, &err));
        CHECK(sceneJsonFinite(dt));
        const std::string lt = dt + "/render.log";
        const std::string ct = std::string("\"\"") + KNIT_PPT_EXE + "\" render \"" + dt + "/scene.json\" --output \"" + dt +
                               "/out.png\" --width 64 --height 64 --spp 1 > \"" + lt + "\" 2>&1\"";
        const int rt = std::system(ct.c_str());
        std::ifstream ft(lt);
        const std::string tt((std::istreambuf_iterator<char>(ft)), std::istreambuf_iterator<char>());
        if (rt != 0) std::printf("%s\n", tt.c_str());
        CHECK(rt == 0);
        CHECK(std::filesystem::exists(dt + "/out.png"));
        for (const char* bad : {"strand.", "instance.", "error", "Error", "refus", "texture load failed", "cannot"})
            if (tt.find(bad) != std::string::npos) {
                std::printf("%s\n", tt.c_str());
                CHECK(!"renderer log mentions a texture/strand problem or an error");
            }
    }
#endif
}

// A tiny grouped scene (two objects, both fibre populations, a diffuse and a hair fibre
// material) loaded and rendered by the actual renderer with its strict scene checks (unknown
// keys fatal). Skipped when ppt.exe is not configured.
TEST(ppt_load_grouped_strict) {
#ifndef KNIT_PPT_EXE
    SKIP("KNIT_PPT_EXE not configured");
#else
    if (!std::filesystem::exists(KNIT_PPT_EXE)) SKIP("the configured ppt.exe is missing");
    const std::string d = scratchDir("ppt_grouped");
    ExportScene ex(d);
    ex.fibreLook.random.model = MaterialSpec::Model::Fibre;  // random fibres a hair BSDF, oriented diffuse
    TriMesh m;
    m.pos = {{-0.02, -0.02, -0.003}, {0.03, -0.02, -0.003}, {0.03, 0.03, -0.003}, {-0.02, 0.03, -0.003}};
    m.tri = {{0, 1, 2}, {0, 2, 3}};
    std::string err;
    CHECK(ex.addMesh("plate", "wool_white", m, &err));
    ex.place("plate", {"a", "object a", ""});
    for (const auto& [object, yarn] : {std::pair<std::string, std::string>{"a", "wool_white"}, {"b", "wool_brown"}}) {
        const double shift = object == "a" ? 0.0 : 0.5;
        CHECK(ex.addStrands(object + "_fibres_oriented", ex.orientedFibreMaterial(yarn), tinyStrands(1.0 + shift), &err));
        CHECK(ex.addStrands(object + "_fibres_random", ex.randomFibreMaterial(yarn), tinyStrands(0.5 + shift), &err));
        ex.place(object + "_fibres_oriented", {object, "object " + object, ExportScene::kOrientedGroup});
        ex.place(object + "_fibres_random", {object, "object " + object, ExportScene::kRandomGroup});
    }
    ex.sortParts();
    CHECK(writeSceneJson(ex, &err));
    CHECK(sceneJsonFinite(d));
    const std::string j = sceneJsonText(d);
    CHECK(j.find("\"groups\": [") != std::string::npos && j.find("\"type\": \"hair\"") != std::string::npos);
    _putenv_s("PPT_STRICT_SCENE", "1");
    const std::string log = d + "/render.log";
    const std::string cmd = std::string("\"\"") + KNIT_PPT_EXE + "\" render \"" + d + "/scene.json\" --output \"" + d +
                            "/out.png\" --width 64 --height 64 --spp 1 > \"" + log + "\" 2>&1\"";
    const int rc = std::system(cmd.c_str());
    _putenv_s("PPT_STRICT_SCENE", "");
    std::ifstream f(log);
    const std::string text((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    if (rc != 0) std::printf("%s\n", text.c_str());
    CHECK(rc == 0);
    CHECK(std::filesystem::exists(d + "/out.png"));
    for (const char* bad : {"strand.", "instance.", "group.", "unknown", "error", "Error", "refus", "warn"})
        if (text.find(bad) != std::string::npos) {
            std::printf("%s\n", text.c_str());
            CHECK(!"renderer log mentions a group/strand/instance problem, an unknown key or an error");
        }
#endif
}

TEST(export_mix_partition) {
    // After the ONE global budget threshold, each population is its own strand part: the
    // parts hold exactly the kept fibres of their population (disjoint, none lost);
    // random fibres use their own material (neutral grey), which is written only when used.
    const StraightSet s = straightSet({0});
    SceneMapping map;
    map.sigma = 68.103;
    for (double frac : {0.5, 0.0}) {
        const FibreSet fs = growYarnFibres(s.cs, s.sm, s.pp, mixParams(frac, 30), nullptr);
        std::vector<std::pair<float, int>> all;
        std::vector<float> keys[2];
        long long segs = 0;
        for (const Fibre& f : fs.fibres) {
            all.push_back({float(f.key), f.segments()});
            keys[static_cast<int>(f.population)].push_back(float(f.key));
            segs += f.segments();
        }
        const float q = budgetThreshold(all, segs / 2);
        const StrandSetF d = fibresToStrands(fs, map, 0), r = fibresToStrands(fs, map, 1);
        CHECK(d.strands() + r.strands() == fs.fibres.size());
        const StrandSetF dk = keepStrands(d, keys[0], q), rk = keepStrands(r, keys[1], q);
        size_t wantD = 0, wantR = 0;
        for (const Fibre& f : fs.fibres)
            if (float(f.key) < q) (f.population == FibrePopulation::Random ? wantR : wantD) += 1;
        std::printf("  random %.1f: %zu + %zu fibres -> kept %zu directed, %zu random (q %.4f)\n", frac, d.strands(), r.strands(),
                    dk.strands(), rk.strands(), double(q));
        CHECK(dk.strands() == wantD && rk.strands() == wantR);
        CHECK(wantD + wantR < fs.fibres.size());  // the budget did thin
        for (bool hair : {false, true}) {
            const std::string dir = scratchDir("mix_partition");
            ExportScene ex(dir);
            ex.fibreLook.hairMaterials = hair;
            std::string err;
            CHECK(ex.addStrands("p_fibres_oriented", ex.orientedFibreMaterial("wool_white"), dk, &err));
            CHECK(ex.addStrands("p_fibres_random", ex.randomFibreMaterial("wool_white"), rk, &err));
            CHECK(ex.strandParts().size() == (rk.strands() ? 2u : 1u));
            CHECK(writeSceneJson(ex, &err));
            std::ifstream f(dir + "/scene.json");
            const std::string j((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
            const std::string grey = "fibre_random_wool_white";
            CHECK((j.find("\"id\": \"" + grey + "\"") != std::string::npos) == (rk.strands() > 0));
            if (rk.strands()) {
                CHECK(j.find(hair ? "[0.550000, 0.550000, 0.550000]" : "[0.550, 0.550, 0.550]") != std::string::npos);
                CHECK(j.find("\"path\": \"strands/p_fibres_random.pptstrands\", \"material\": \"" + grey + "\"") != std::string::npos);
                // Hair mode makes it a hair BSDF in its own colour.
                const size_t p = j.find("\"id\": \"" + grey + "\"");
                CHECK((j.substr(p, j.find('\n', p) - p).find("\"type\": \"hair\"") != std::string::npos) == hair);
            }
        }
    }
}

TEST(ppt_fibre_material_is_yarn) {
    // The oriented fibres' own material (edited apart from the yarn's) takes the yarn entry's
    // appearance by default (diffuse + sheen, the colour authored once); the hair-BSDF option
    // reads the yarn colour.
    const std::string d = scratchDir("fibre_material");
    std::string err;
    for (bool fibre : {false, true}) {
        ExportScene ex(d);
        ex.fibreLook.hairMaterials = fibre;
        TriMesh m;
        m.pos = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
        m.tri = {{0, 1, 2}};
        CHECK(ex.addMesh("piece", "wool_white", m, &err));
        CHECK(ex.addStrands("piece_fibres_oriented", ex.orientedFibreMaterial("wool_white"), tinyStrands(), &err));
        CHECK(writeSceneJson(ex, &err));
        std::ifstream f(d + "/scene.json");
        const std::string j((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
        const size_t p = j.find("{\"id\": \"fibre_oriented_wool_white\"");
        CHECK(p != std::string::npos);
        const std::string line = j.substr(p, j.find('\n', p) - p);
        CHECK((line.find("\"type\": \"hair\"") != std::string::npos) == fibre);
        // Hair: the mean of the yarn's colour and that colour as one fibre shows it.
        CHECK(line.find(fibre ? "\"color\": [0.334537, 0.312873, 0.257308]" : "\"albedo\": [0.800, 0.780, 0.720]") != std::string::npos);
        CHECK(j.find("\"path\": \"strands/piece_fibres_oriented.pptstrands\", \"material\": \"fibre_oriented_wool_white\"") !=
              std::string::npos);
    }
}

// Each fibre material id resolves by one rule (ExportScene::fibreMaterialSpec): every rule's
// case, its model and colour checked.
TEST(ppt_fibre_material_resolution) {
    using Model = MaterialSpec::Model;
    const Vec3 white{0.80, 0.78, 0.72}, brown{0.42, 0.20, 0.08}, grey{0.55, 0.55, 0.55};
    MaterialSpec hairParams;  // material.fibre_directed with model = fibre
    hairParams.model = Model::Fibre;
    hairParams.color = {0.9, 0.1, 0.1};
    hairParams.fibreLongRoughness = 0.21;
    hairParams.fibreAzimRoughness = 0.42;
    hairParams.fibreIor = 1.6;
    hairParams.fibreCuticleDeg = 3;
    auto resolve = [](ExportScene& ex, bool random, const std::string& yarn) {
        std::string err;
        const std::optional<MaterialSpec> m = ex.fibreMaterialSpec(random ? ex.randomFibreMaterial(yarn) : ex.orientedFibreMaterial(yarn), &err);
        CHECK(m);
        return m.value_or(MaterialSpec{});
    };
    auto hairOf = [&](const MaterialSpec& m, const Vec3& colour) {
        return m.model == Model::Fibre && norm(m.color - colour) < 1e-8 && m.fibreLongRoughness == hairParams.fibreLongRoughness &&
               m.fibreAzimRoughness == hairParams.fibreAzimRoughness && m.fibreIor == hairParams.fibreIor &&
               m.fibreCuticleDeg == hairParams.fibreCuticleDeg;
    };
    {  // (a) defaults: oriented = the yarn entry (diffuse + sheen), random = neutral grey
        ExportScene ex("");
        CHECK(resolve(ex, false, "wool_white") == ex.baseMaterials[0].spec);
        CHECK(resolve(ex, false, "wool_brown").color == brown);
        const MaterialSpec r = resolve(ex, true, "wool_white");
        CHECK(r.model == Model::DiffuseSheen && r.color == grey && r == defaultRandomFibre());
    }
    {  // (b) a custom random colour, diffuse
        ExportScene ex("");
        ex.fibreLook.random.color = {0.1, 0.2, 0.3};
        CHECK(resolve(ex, true, "wool_brown") == ex.fibreLook.random);
    }
    {  // (c) the random fibres' own hair model beside diffuse oriented fibres
        ExportScene ex("");
        ex.fibreLook.random = hairParams;
        ex.fibreLook.random.color = {0.3, 0.3, 0.2};
        CHECK(resolve(ex, true, "wool_white") == ex.fibreLook.random);
        CHECK(resolve(ex, false, "wool_white").model == Model::DiffuseSheen);
    }
    {  // (d) hair mode: oriented the mean of the yarn colour and one fibre of it (literals: the fuzz
       // colours must not move), a diffuse random spec gets the hair params in its own colour
        ExportScene ex("");
        ex.fibreLook.hairMaterials = true;
        ex.fibreLook.hair = hairParams;
        CHECK(hairOf(resolve(ex, false, "wool_white"), Vec3{0.344826837, 0.323036184, 0.266941747}) &&
              hairOf(resolve(ex, false, "wool_brown"), Vec3{0.098957305, 0.034579436, 0.011452300}));
        ex.fibreLook.random.color = {0.1, 0.2, 0.3};
        CHECK(hairOf(resolve(ex, true, "wool_white"), Vec3{0.1, 0.2, 0.3}));
    }
    {  // (e) fibre_inherits_yarn = false: the own spec (colour and model), one id per yarn colour
        ExportScene ex("");
        ex.fibreLook.orientedOverride = hairParams;
        CHECK(resolve(ex, false, "wool_white") == hairParams && resolve(ex, false, "wool_brown") == hairParams);
        CHECK(ex.orientedFibreMaterial("wool_white") != ex.orientedFibreMaterial("wool_brown"));
        CHECK(resolve(ex, true, "wool_white").model == Model::DiffuseSheen);  // the random spec is its own
    }
    {  // (e2) fibre_random_inherits_yarn: the random hair spec in the yarn's colour as it stands,
       // under either material.fibre_hair_colour, a zero channel at the renderer's floor
        ExportScene ex("");
        ex.fibreLook.random = hairParams;
        ex.fibreLook.random.fibreAzimRoughness = 0.3;
        ex.fibreLook.randomInheritsYarn = true;
        ex.baseMaterials[0].spec.color = {0.651, 0.054, 0.0};
        for (bool keepsYarnColour : {false, true}) {
            ex.fibreLook.hairKeepsYarnColour = keepsYarnColour;
            const MaterialSpec red = resolve(ex, true, "wool_white");
            CHECK(red.model == Model::Fibre && red.fibreLongRoughness == hairParams.fibreLongRoughness);
            CHECK(red.color == (Vec3{0.651, 0.054, 1e-4}));
            CHECK(resolve(ex, true, "wool_brown").color == brown);
        }
        ex.baseMaterials[0].spec.color = {1.0, 1.0, 1.0};
        CHECK(resolve(ex, true, "wool_white").color == (Vec3{1.0, 1.0, 1.0}));
        ex.fibreLook.random = MaterialSpec{};
        ex.baseMaterials[0].spec.color = {0.651, 0.054, 0.0};
        CHECK(resolve(ex, true, "wool_white").color == (Vec3{0.651, 0.054, 0.0}));  // diffuse: no floor
        ex.fibreLook.randomInheritsYarn = false;
        ex.fibreLook.random = hairParams;
        CHECK(resolve(ex, true, "wool_white").color == hairParams.color);  // off: its own colour
        ex.fibreLook.randomInheritsYarn = true;
        // (e3) hairKeepsYarnColour: a dense coat. The fuzz keeps the colour less the shine, at most the
        // lowest channel, so the three differences and the hue stay; strays do not depend on it.
        ExportScene coat("");
        coat.fibreLook.hairMaterials = true;
        coat.fibreLook.hair = hairParams;
        coat.fibreLook.hair.fibreAzimRoughness = 0.3;
        coat.fibreLook.random = coat.fibreLook.hair;
        coat.fibreLook.randomInheritsYarn = true;
        coat.baseMaterials[0].spec.color = {0.578, 0.376, 0.171};
        const Vec3 layerFuzz = resolve(coat, false, "wool_white").color;
        CHECK(norm(layerFuzz - Vec3{0.162757910, 0.078230239, 0.025998074}) < 1e-8);  // must not move
        CHECK(resolve(coat, true, "wool_white").color == coat.baseMaterials[0].spec.color);
        coat.fibreLook.hairKeepsYarnColour = true;
        const Vec3 tan = resolve(coat, false, "wool_white").color;
        CHECK(std::abs(tan.x - 0.518) < 1e-12 && std::abs(tan.y - 0.316) < 1e-12 && std::abs(tan.z - 0.111) < 1e-12);
        CHECK(resolve(coat, true, "wool_white").color == coat.baseMaterials[0].spec.color);
        coat.baseMaterials[0].spec.color = {0.60, 0.05, 0.02};  // the lowest channel limits the shine
        const Vec3 deep = resolve(coat, false, "wool_white").color;
        CHECK(std::abs(deep.x - 0.58) < 1e-12 && std::abs(deep.y - 0.03) < 1e-12 && deep.z == 1e-4);
        coat.baseMaterials[0].spec.color = {0.03, 0.03, 0.03};  // grey stays grey: no tint out of a dark yarn
        const Vec3 dark = resolve(coat, false, "wool_white").color;
        CHECK(dark.x == dark.y && dark.y == dark.z);
        std::string err;
        CHECK(!ex.fibreMaterialSpec(ex.randomFibreMaterial("no_such_yarn"), &err) && err.find("no_such_yarn") != std::string::npos);
    }
    {  // (f) the pieces CLI's hair mode: fibreParams still says diffuse, the fibres are hair all the same
        ExportScene ex(scratchDir("fibre_resolution_cli"));
        ex.fibreLook.hairMaterials = true;
        CHECK(ex.fibreLook.hair.model == Model::DiffuseSheen);
        const MaterialSpec o = resolve(ex, false, "wool_white");
        CHECK(o.model == Model::Fibre && o.color.x > 0 && o.color.x < white.x);  // the fuzz's mean, darker than the yarn's dense mass
        std::string err;
        CHECK(ex.addStrands("s_fibres_oriented", ex.orientedFibreMaterial("wool_white"), tinyStrands(), &err));
        CHECK(writeSceneJson(ex, &err));
        const std::string j = sceneJsonText(ex.folder());
        const size_t p = j.find("{\"id\": \"fibre_oriented_wool_white\"");
        CHECK(p != std::string::npos && j.substr(p, j.find('\n', p) - p).find("\"type\": \"hair\"") != std::string::npos);
    }
    {  // refused: an id never registered, a yarn without a base entry (writeScene fails, writes nothing)
        ExportScene ex(scratchDir("fibre_resolution_bad"));
        std::string err;
        CHECK(!ex.fibreMaterialSpec("fibre_oriented_wool_white", &err) && !err.empty());
        CHECK(ex.addStrands("s_fibres_oriented", ex.orientedFibreMaterial("yarn_missing"), tinyStrands(), &err));
        err.clear();
        CHECK(!writeSceneJson(ex, &err) && err.find("yarn_missing") != std::string::npos);
        CHECK(!std::filesystem::exists(ex.folder() + "/scene.json"));
    }
}

// The scene graph (renderer groups[]): a group per source object, a child group per fibre
// type holding the per-piece fibre parts, each type x yarn colour its own material; a
// population emptied by the budget leaves neither a group nor a material.
TEST(ppt_export_fibre_groups) {
    ExportScene ex(scratchDir("fibre_groups"));
    TriMesh m;
    m.pos = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
    m.tri = {{0, 1, 2}};
    std::string err;
    for (const auto& [object, yarn] : {std::pair<std::string, std::string>{"bear.002", "wool_white"}, {"bear.009", "wool_brown"}}) {
        const std::string stem = "piece_" + object;
        CHECK(ex.addMesh(stem, yarn, m, &err));
        ex.place(stem, {object, object, ""});
        CHECK(ex.addStrands(stem + "_fibres_oriented", ex.orientedFibreMaterial(yarn), tinyStrands(), &err));
        // The brown piece's random fibres were all dropped by the budget (an empty part is omitted).
        const StrandSetF random = yarn == "wool_white" ? tinyStrands() : StrandSetF{};
        CHECK(ex.addStrands(stem + "_fibres_random", ex.randomFibreMaterial(yarn), random, &err));
        ex.place(stem + "_fibres_oriented", {object, object, ExportScene::kOrientedGroup});
        ex.place(stem + "_fibres_random", {object, object, ExportScene::kRandomGroup});
    }
    CHECK(ex.addMesh("loose", "black", m, &err));  // not placed: a root instance
    // Fur and scene hair: strand parts of their own materials, which are no wool fibres.
    MaterialSpec furHair;
    furHair.model = MaterialSpec::Model::Fibre;
    ex.baseMaterials.push_back({"fur_body", furHair});
    ex.baseMaterials.push_back({"hair_c102030", furHair});
    CHECK(ex.addStrands("body_fur", "fur_body", tinyStrands(), &err));
    CHECK(ex.addStrands("scalp_hair", "hair_c102030", tinyStrands(), &err));
    ex.sortParts();
    CHECK(writeSceneJson(ex, &err));
    CHECK(sceneJsonFinite(ex.folder()));
    const std::string j = sceneJsonText(ex.folder());
    auto lineWith = [&](const std::string& a, const std::string& b) {
        for (size_t p = 0, e; p < j.size(); p = e + 1) {
            e = j.find('\n', p);
            if (e == std::string::npos) e = j.size();
            const std::string line = j.substr(p, e - p);
            if (line.find(a) != std::string::npos && line.find(b) != std::string::npos) return line;
        }
        return std::string();
    };
    // The scene node (projectName), the objects by key under it, then (object, group) pairs that
    // hold a part, then the ground under the scene node.
    const char* groups[][3] = {{"g0001", "knitted_bears", ""},          {"g0002", "bear.002", "g0001"},
                               {"g0003", "bear.009", "g0001"},          {"g0004", "fibres_oriented", "g0002"},
                               {"g0005", "fibres_random", "g0002"},     {"g0006", "fibres_oriented", "g0003"},
                               {"g0007", "Ground", "g0001"}};
    for (const auto& g : groups) {
        const std::string line = lineWith(std::string("{\"id\": \"") + g[0] + "\", \"name\": \"" + g[1] + "\"", "\"origin\": \"User\"");
        CHECK(!line.empty());
        CHECK((line.find("\"parent_id\"") != std::string::npos) == (*g[2] != 0));
        if (*g[2]) CHECK(line.find(std::string("\"parent_id\": \"") + g[2] + "\"") != std::string::npos);
    }
    // Pivot: the object's bottom centre over all its parts; its child groups share it.
    auto pivotOf = [&](const char* id) {
        const std::string line = lineWith(std::string("{\"id\": \"") + id + "\"", "\"pivot\": [");
        const size_t p = line.find("\"pivot\": [");
        return p == std::string::npos ? std::string() : line.substr(p + 10, line.find(']', p) - p - 10);
    };
    {
        Vec3 l{1e300, 1e300, 1e300}, h{-1e300, -1e300, -1e300};
        std::vector<PartFile> everyPart = ex.meshParts();
        everyPart.insert(everyPart.end(), ex.strandParts().begin(), ex.strandParts().end());
        for (const PartFile& part : everyPart)
            if (part.name.rfind("piece_bear.002", 0) == 0) {
                const Box& b = ex.partBounds().at(part.name);
                l = {std::min(l.x, b.lo.x), std::min(l.y, b.lo.y), std::min(l.z, b.lo.z)};
                h = {std::max(h.x, b.hi.x), std::max(h.y, b.hi.y), std::max(h.z, b.hi.z)};
            }
        char want[128];
        std::snprintf(want, sizeof want, "%.6f, %.6f, %.6f", 0.5 * (l.x + h.x), l.y, 0.5 * (l.z + h.z));
        CHECK(pivotOf("g0002") == want);
    }
    CHECK(pivotOf("g0004") == pivotOf("g0002") && pivotOf("g0005") == pivotOf("g0002") && pivotOf("g0006") == pivotOf("g0003"));
    CHECK(j.find("\"id\": \"g0008\"") == std::string::npos);
    CHECK(j.find("\"next_group_seq\": 8,") != std::string::npos);
    auto groupOf = [&](const std::string& name) {
        const std::string line = lineWith("\"name\": \"" + name + "\"", "\"asset_mesh_index\"");
        CHECK(!line.empty());
        const size_t p = line.find("\"group_id\": \"");
        return p == std::string::npos ? std::string() : line.substr(p + 13, 5);
    };
    CHECK(groupOf("piece_bear.002") == "g0002" && groupOf("piece_bear.009") == "g0003");
    CHECK(groupOf("piece_bear.002_fibres_oriented") == "g0004" && groupOf("piece_bear.002_fibres_random") == "g0005");
    CHECK(groupOf("piece_bear.009_fibres_oriented") == "g0006");
    CHECK(groupOf("loose").empty());  // unplaced: root
    CHECK(groupOf(ExportScene::kGround) == "g0007");
    // The wool look's thin ends: every wool fibre instance has one shape under its own part index,
    // nothing else has any.
    int shaped = 0;
    for (size_t p = 0; p < ex.meshParts().size() + 1 + ex.strandParts().size(); ++p) {  // + 1: the ground
        const std::string index = "\"asset_mesh_index\": " + std::to_string(p) + ",";
        const std::string line = lineWith("\"id\": \"inst_", index);
        CHECK(!line.empty());
        const bool woolFibre = line.find("\"id\": \"fibre_oriented_") != std::string::npos || line.find("\"id\": \"fibre_random_") != std::string::npos;
        CHECK((line.find("part_fibre_shapes") != std::string::npos) == woolFibre);
        if (!woolFibre) continue;
        ++shaped;
        const std::string key = "\"part_fibre_shapes\": {\"" + std::to_string(p) + "\": {\"thickness\": \"imported\", \"diameter_scale\": 1.0, \"tip_scale\": 0.2, ";
        CHECK(line.find(key) != std::string::npos);
        CHECK(line.find("\"seed\": " + std::to_string(1000 + p) + "}}") != std::string::npos);
    }
    CHECK(shaped == 3);
    // The yarns' sheen takes the yarn's colour; written here from the default table.
    CHECK(lineWith("{\"id\": \"wool_white\"", "\"tint\": 1.0").size() > 0 && lineWith("{\"id\": \"wool_brown\"", "\"tint\": 1.0").size() > 0);
    // Each written fibre material once; the emptied population's material not at all.
    auto count = [&](const std::string& s) {
        size_t n = 0;
        for (size_t p = j.find(s); p != std::string::npos; p = j.find(s, p + 1)) ++n;
        return n;
    };
    for (const char* id : {"fibre_oriented_wool_white", "fibre_random_wool_white", "fibre_oriented_wool_brown"})
        CHECK(count(std::string("{\"id\": \"") + id + "\"") == 1);
    CHECK(j.find("fibre_random_wool_brown") == std::string::npos);
    CHECK(lineWith("{\"id\": \"fibre_random_wool_white\"", "\"albedo\": [0.550, 0.550, 0.550]").size() > 0);
    CHECK(lineWith("{\"id\": \"fibre_oriented_wool_brown\"", "\"albedo\": [0.420, 0.200, 0.080]").size() > 0);
}

// The wool look a profile carries (file_presets.cpp fibreKeys): a translucent yarn writes the renderer
// fabric's translucency block (only when set), and an export whose fibres carry a late taper writes
// tip_scale 1, so the renderer does not narrow them a second time.
TEST(ppt_export_wool_look_keys) {
    ExportScene ex(scratchDir("wool_look_keys"));
    TriMesh m;
    m.pos = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
    m.tri = {{0, 1, 2}};
    std::string err;
    MaterialSpec translucent = ex.baseMaterials.front().spec;
    translucent.translucency = 0.3;
    ex.baseMaterials.front().spec = translucent;     // wool_white
    ex.fibreTipScale = 1.0;                         // collectScene, under fibres.tail_taper > 0
    CHECK(ex.addMesh("piece_a", "wool_white", m, &err));
    CHECK(ex.addMesh("piece_b", "wool_brown", m, &err));
    ex.place("piece_a", {"a", "a", ""});
    CHECK(ex.addStrands("piece_a_fibres_oriented", ex.orientedFibreMaterial("wool_white"), tinyStrands(), &err));
    ex.place("piece_a_fibres_oriented", {"a", "a", ExportScene::kOrientedGroup});
    ex.sortParts();
    CHECK(writeSceneJson(ex, &err));
    const std::string j = sceneJsonText(ex.folder());
    auto lineOf = [&](const std::string& id) {
        const size_t at = j.find("{\"id\": \"" + id + "\"");
        return at == std::string::npos ? std::string() : j.substr(at, j.find('\n', at) - at);
    };
    CHECK(lineOf("wool_white").find("\"type\": \"fabric\"") != std::string::npos);
    CHECK(lineOf("wool_white").find("\"translucency\": {\"transmit_share\": 0.3}") != std::string::npos);
    CHECK(!lineOf("wool_brown").empty() && lineOf("wool_brown").find("translucency") == std::string::npos);
    CHECK(j.find("\"tip_scale\": 1, ") != std::string::npos && j.find("\"tip_scale\": 0.2") == std::string::npos);
}

// Without placements and fibre materials, scene.json is byte-identical to the exporter before
// the scene graph, re-pinned for material schema 12, for the table listing only referenced
// materials and for the yarn's sheen tint 1.0 (FNV-1a 64 of the text, as read back).
// All three are compared at the end, so a re-pin sees every new value at once.
TEST(ppt_unplaced_bytes_unchanged) {
    uint64_t got[3] = {};
    const uint64_t want[3] = {0x7d0c9837b8ab0dfdull, 0x66893d626646f0c5ull, 0xe970827ba4805336ull};
    TriMesh m;
    m.pos = {{-0.02, -0.02, -0.003}, {0.03, -0.02, -0.003}, {0.03, 0.03, -0.003}, {-0.02, 0.03, -0.003}};
    m.tri = {{0, 1, 2}, {0, 2, 3}};
    std::string err;
    {  // a mesh only
        ExportScene ex(scratchDir("ppt_bytes_mesh"));
        CHECK(ex.addMesh("plate", "wool_white", m, &err));
        CHECK(writeSceneJson(ex, &err));
        got[0] = fnv1a64(sceneJsonText(ex.folder()));
    }
    {  // ground off
        ExportScene ex(scratchDir("ppt_bytes_off"));
        ex.groundPlane = false;
        CHECK(ex.addMesh("plate", "wool_white", m, &err));
        CHECK(writeSceneJson(ex, &err));
        got[1] = fnv1a64(sceneJsonText(ex.folder()));
    }
    {  // a strand part on a base material
        ExportScene ex(scratchDir("ppt_bytes_strand"));
        CHECK(ex.addMesh("plate", "wool_white", m, &err));
        StrandSetF s;
        s.add({0, -0.01f, 0, 0.004f, -0.01f, 0}, {2e-4f, 2e-4f});
        CHECK(ex.addStrands("fibres", "wool_white", s, &err));
        CHECK(writeSceneJson(ex, &err));
        got[2] = fnv1a64(sceneJsonText(ex.folder()));
    }
    for (int k = 0; k < 3; ++k)
        if (got[k] != want[k]) std::printf("  scene %d: FNV-1a 0x%016llxull, pinned 0x%016llxull\n", k, (unsigned long long)got[k], (unsigned long long)want[k]);
    for (int k = 0; k < 3; ++k) CHECK(got[k] == want[k]);
}

// Renderer material schema 12 (one concrete type per material, photon_path_tracer RDR-320):
// each MaterialSpec maps to the type migrate_materials.py gives its old layer stack; nothing
// writes a layer, an AO map or path guiding.
TEST(ppt_material_schema12) {
    TriMesh m;
    m.pos = {{-0.02, -0.02, -0.003}, {0.03, -0.02, -0.003}, {0.03, 0.03, -0.003}, {-0.02, 0.03, -0.003}};
    m.tri = {{0, 1, 2}, {0, 2, 3}};
    auto spec = [](double sheen, bool coat, bool fibre) {
        MaterialSpec s;
        s.sheen = sheen;
        s.clearcoat = coat;
        if (fibre) s.model = MaterialSpec::Model::Fibre;
        return s;
    };
    ExportScene ex(scratchDir("ppt_schema12"));
    ex.baseMaterials = {{"m_hair", spec(1, false, true)},   {"m_fabric", spec(1, false, false)}, {"m_fabric_half", spec(0.5, false, false)},
                        {"m_coat", spec(0, true, false)},   {"m_coat_sheen", spec(1, true, false)}, {"m_plain", spec(0, false, false)}};
    std::string err;
    for (const BaseMaterial& b : ex.baseMaterials) CHECK(ex.addMesh("plate_" + b.id, b.id, m, &err));
    CHECK(writeSceneJson(ex, &err));
    const std::string j = sceneJsonText(ex.folder());
    auto line = [&](const std::string& id) {
        const size_t p = j.find("{\"id\": \"" + id + "\"");
        CHECK(p != std::string::npos);
        return p == std::string::npos ? std::string() : j.substr(p, j.find('\n', p) - p);
    };
    auto has = [](const std::string& l, const std::string& s) { return l.find(s) != std::string::npos; };
    struct Want { const char* id; const char* type; bool sheen, coat; const char* coverage; };
    const Want wants[] = {{"m_hair", "hair", false, false, nullptr},
                          {"m_fabric", "fabric", true, false, nullptr},
                          {"m_fabric_half", "fabric", true, false, "\"coverage\": 0.5}"},
                          {"m_coat", "coated_diffuse", false, true, "\"coverage\": 0.3}"},
                          {"m_coat_sheen", "fabric", true, false, nullptr},  // no type has a coat over a sheen
                          {"m_plain", "diffuse", false, false, nullptr},
                          {ExportScene::kGround, "diffuse", false, false, nullptr}};
    for (const Want& w : wants) {
        const std::string l = line(w.id);
        CHECK(has(l, "\"schema_version\": 12, ") && has(l, std::string("\"type\": \"") + w.type + "\""));
        CHECK(has(l, "\"fibre\": {\"color\": ") == (std::string(w.type) == "hair"));
        CHECK(has(l, "\"diffuse\": {\"albedo\": ") == (std::string(w.type) != "hair"));
        CHECK(has(l, "\"sheen\": {") == w.sheen && has(l, "\"coat\": {") == w.coat);
        CHECK(has(l, "coverage") == (w.coverage != nullptr));
        if (w.coverage) CHECK(has(l, w.coverage));
    }
    size_t n = 0;
    for (size_t p = j.find("\"schema_version\""); p != std::string::npos; p = j.find("\"schema_version\"", p + 1)) ++n;
    CHECK(n == std::size(wants));
    for (const char* bad : {"\"layers\"", "\"weight\"", "\"params\"", "ao_tex", "guiding"}) CHECK(j.find(bad) == std::string::npos);
}

// scene.json lists what a part references: a base entry stays the colour source of its textured
// and fibre materials without being written itself.
TEST(ppt_material_table_lists_only_referenced) {
    TriMesh m;
    m.pos = {{-0.02, -0.02, -0.003}, {0.03, -0.02, -0.003}, {0.03, 0.03, -0.003}, {-0.02, 0.03, -0.003}};
    m.tri = {{0, 1, 2}, {0, 2, 3}};
    MaterialSpec used;
    used.color = {0.1, 0.2, 0.3};
    ExportScene ex(scratchDir("ppt_referenced"));
    ex.baseMaterials = {{"used", used}, {"unused", MaterialSpec{}}};
    ex.yarnEntries = 2;
    ex.plyShading.on = true;
    std::string err;
    CHECK(ex.addMesh("tube", ExportScene::yarnMaterial("used"), m, &err));
    CHECK(ex.addStrands("fibres", ex.orientedFibreMaterial("used"), tinyStrands(), &err));
    ex.orientedFibreMaterial("unused");  // registered, but no part uses it
    CHECK(writeSceneJson(ex, &err));
    const std::string j = sceneJsonText(ex.folder());
    auto line = [&](const std::string& id) {
        const size_t p = j.find("{\"id\": \"" + id + "\"");
        return p == std::string::npos ? std::string() : j.substr(p, j.find('\n', p) - p);
    };
    for (const char* id : {"used_yarn", "fibre_oriented_used"})
        CHECK(line(id).find("\"albedo\": [0.100, 0.200, 0.300]") != std::string::npos);
    for (const char* id : {"used", "unused", "unused_yarn", "fibre_oriented_unused"}) CHECK(line(id).empty());
}

// Export with plies as maps: strand parts reference untextured materials only; the tube
// material references the map files, which decode as specified (EXR float height, linear
// normal, sRGB-encoded cavity multiplier); no AO map (renderer schema 12 has none).
TEST(ppt_fibre_parts_untextured) {
    for (PlyModeOptions::Shading shading : {PlyModeOptions::Shading::Normal, PlyModeOptions::Shading::Bump}) {
        const bool bump = shading == PlyModeOptions::Shading::Bump;
        const std::string d = scratchDir("ppt_textured");
        ExportScene ex(d);
        PlyModeOptions o;
        o.shading = shading;
        o.cavity = true;
        o.spec.size = 64;
        const PlyGeometry g;
        std::string err;
        o.spec.g = g;
        o.spec.sigma = 68.103;
        const PlyMapSet set = plyMapSetOf(o.spec);
        const std::optional<PlyShading> plyMaps = setupPlyTextures(o, ex.folder(), 68.103, &err);
        CHECK(plyMaps);
        ex.plyShading = *plyMaps;
        TubeStyle style;
        style.uvRepeat = plyMaps->uvRepeat;
        style.uvWrap = plyMaps->uvWrap;
        CHECK(style.uvRepeat == g.pitch && style.uvWrap > 0 && ex.plyShading.on);
        const StrandDemo sd = makeStrandDemo(6.0, 1.12);
        SceneMapping map;
        map.sigma = 68.103;
        CHECK(addExportTube(ex, "strand", ExportScene::yarnMaterial("wool_white"), sd.curves, sd.pp, map, &err, style));
        CHECK(ex.addStrands("strand_fibres_oriented", ex.orientedFibreMaterial("wool_white"), tinyStrands(), &err));
        CHECK(ex.addStrands("strand_fibres_random", ex.randomFibreMaterial("wool_white"), tinyStrands(), &err));
        ex.sortParts();
        CHECK(writeSceneJson(ex, &err));
        std::ifstream f(d + "/scene.json");
        const std::string j((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
        auto matLine = [&](const std::string& id) {
            const size_t p = j.find("{\"id\": \"" + id + "\"");
            CHECK(p != std::string::npos);
            return j.substr(p, j.find('\n', p) - p);
        };
        const char* maps[] = {"normal_tex", "bump_tex", "ao_tex", "albedo_tex", "displacement_tex"};
        for (const PartFile& e : ex.strandParts())
            for (const char* k : maps) CHECK(matLine(e.material).find(k) == std::string::npos);
        auto texPath = [&](const std::string& line, const char* key) {
            const size_t p = line.find(std::string("\"") + key + "\": \"");
            CHECK(p != std::string::npos);
            const size_t s = p + std::strlen(key) + 5;
            const std::string rel = line.substr(s, line.find('"', s) - s);
            CHECK(std::filesystem::exists(d + "/" + rel));
            return d + "/" + rel;
        };
        CHECK(ex.meshParts().size() == 1 && ex.meshParts()[0].material == "wool_white_yarn");
        const std::string tube = matLine("wool_white_yarn");
        CHECK(tube.find("\"albedo\": [0.800, 0.780, 0.720]") != std::string::npos);  // the entry's colour
        CHECK((tube.find("normal_tex") != std::string::npos) == !bump);
        CHECK(tube.find("\"type\": \"fabric\"") != std::string::npos && tube.find("ao_tex") == std::string::npos);
        const PlyMaps ref = makePlyMaps(set.prof, set.uvRepeat, 64, 64);
        int w, h;
        if (bump) {
            std::vector<float> hf;
            CHECK(readExrRgbF32(texPath(tube, "bump_tex"), w, h, hf));
            CHECK(w == 64 && h == 64);
            for (size_t i = 0; i < ref.H.size(); ++i) CHECK(hf[3 * i] == ref.H[i] && hf[3 * i + 1] == ref.H[i]);
            char b[64];
            std::snprintf(b, sizeof b, "\"bump_scale\": %.8g", set.bumpScale);
            CHECK(tube.find(b) != std::string::npos);
            CHECK_NEAR(set.bumpScale, (set.prof.rhoMax - set.prof.rhoMin) / 68.103, 1e-12);
        } else {
            std::vector<uint8_t> n;
            CHECK(readPngRgb8(texPath(tube, "normal_tex"), w, h, n));
            CHECK(w == 64 && h == 64);
            for (size_t i = 0; i < ref.H.size(); ++i) {
                double c[3], l = 0;
                for (int k = 0; k < 3; ++k) {
                    c[k] = n[3 * i + k] / 255.0 * 2 - 1;
                    l += c[k] * c[k];
                    CHECK(std::fabs(c[k] - ref.normal[3 * i + k]) <= 1.0 / 255 + 1e-9);
                }
                CHECK(std::fabs(std::sqrt(l) - 1) <= std::sqrt(3.0) / 255);
            }
        }
        std::vector<uint8_t> cav;
        CHECK(readPngRgb8(texPath(tube, "albedo_tex"), w, h, cav));
        for (size_t i = 0; i < ref.H.size(); ++i) {
            // sRGB 8-bit: the decoded multiplier is the intended linear one within half a code
            // step times the decode slope (<= 1 / 12.92 ... 2.4 / 1.055: < 2.3).
            CHECK(std::fabs(srgbDecode(cav[3 * i] / 255.0) - (1 - o.spec.cavity * (1 - ref.ao[i]))) <= 2.3 * 0.5 / 255);
        }
    }
}

// The grouped tube export partitions the tubes: the group meshes hold exactly the triangles of the
// ungrouped tube, and each group's part is written with its own material.
TEST(tube_split_partition) {
    const PatternParams pp = expandedD1();
    const TriMesh m = gridMesh(40, 30, 1.0, [](double x, double y) { return Vec3{x, y, 0}; },
                               [](double x, double y) { return Vec2{x / 100, y / 100}; });
    MeshKnitParams p = tiledParams(pp);
    p.settle = false;
    PartKnit k;
    std::string err;
    CHECK(knitMesh(m, p, knit_test::cpuNumerics(), k, &err));
    const CurveSet& cs = k.t.curves;
    CHECK(!cs.topo.empty());
    namespace fs = std::filesystem;
    const std::string dir = fs::temp_directory_path().string() + "/knit_s6_tube_split";
    fs::remove_all(dir);
    fs::create_directories(dir + "/meshes");
    const SceneMapping map{};
    ExportScene one(dir), many(dir);
    const std::optional<ExportTube> built = buildExportTube(cs, k.t.pp, map, TubeStyle{}, &err);
    CHECK(built);
    if (!built) return;
    const ExportTube& tube = *built;
    CHECK(one.addMeshPart("whole", "wool_white", tube.mesh, &err));
    const StitchAlongYarn along(cs);
    std::map<int, long> calls;
    TubeGroups groups;
    for (int g = 0; g < 3; ++g) groups.groups.push_back({"_g" + std::to_string(g), "m" + std::to_string(g)});
    for (size_t t = 0; t + 2 < tube.mesh.tri.size(); t += 3) {
        const uint32_t a = tube.mesh.tri[t], b = tube.mesh.tri[t + 1], c = tube.mesh.tri[t + 2];
        const int st = along.at(tube.yarn[a], (tube.s[a] + tube.s[b] + tube.s[c]) / 3);
        const int g = st >= 0 ? st % 3 : -1;
        ++calls[g];
        groups.groupOfTriangle.push_back(g);
    }
    const std::optional<std::vector<std::string>> split = many.addTubeGroups("split", tube, groups, &err);
    CHECK(split);
    const std::vector<std::string> parts = split.value_or(std::vector<std::string>{});
    std::printf("    %zu triangles whole, %zu split over %zu groups\n", one.triangleCount(), many.triangleCount(), parts.size());
    CHECK(one.triangleCount() == many.triangleCount() && one.triangleCount() > 0);
    CHECK(parts.size() == 3 && many.meshParts().size() == 3 && calls.count(-1) == 0);
    long sum = 0;
    for (const auto& [g, n] : calls) sum += n;
    CHECK(sum == static_cast<long>(one.triangleCount()));
    for (size_t i = 0; i < many.meshParts().size(); ++i) CHECK(many.meshParts()[i].material == "m" + std::to_string(i));
    fs::remove_all(dir);
}

// Units contract: a GLB scene exported through the tiled pipeline lands in metres,
// Y-up, with the registry size: the exported vertices span exactly the library bbox / 1000
// (axes (x, z, -y)), every part's re-centring undone.
TEST(glb_units_roundtrip) {
    namespace fs = std::filesystem;
    const std::string dir = fs::temp_directory_path().string() + "/knit_s6_units";
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
    s.workflow = "fur";  // the mesh itself is written: the mapping is under test
    Scene t;
    std::string err, log;
    CHECK(loadBaseParts(e, s, SourceOptions{}, t, &err));
    Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    for (const SourcePrimitive& g : t.source.prims)
        for (const Vec3& p : g.welded.pos) {
            const Vec3 q{p.x / 1000, p.z / 1000, -p.y / 1000};
            lo = {std::min(lo.x, q.x), std::min(lo.y, q.y), std::min(lo.z, q.z)};
            hi = {std::max(hi.x, q.x), std::max(hi.y, q.y), std::max(hi.z, q.z)};
        }
    CHECK_NEAR(hi.y - lo.y, 0.003, 1e-12);
    CHECK_NEAR(lo.y, 0.0, 1e-12);
    ExportOptions o;
    o.folder = dir + "/export";
    fs::remove_all(o.folder);
    CHECK(exportScene({e, t, s}, o, &log, &err));
    Vec3 elo{1e300, 1e300, 1e300}, ehi{-1e300, -1e300, -1e300};
    int files = 0;
    for (const auto& f : fs::directory_iterator(o.folder + "/meshes")) {
        if (f.path().extension() != ".ply") continue;
        if (f.path().stem() == ExportScene::kGround) {  // the staging ground: just under the parts, not a part
            for (const Vec3& q : readPlyPositions(f.path().string())) CHECK(q.y < lo.y && q.y > lo.y - 1e-5);
            continue;
        }
        ++files;
        for (const Vec3& q : readPlyPositions(f.path().string())) {
            elo = {std::min(elo.x, q.x), std::min(elo.y, q.y), std::min(elo.z, q.z)};
            ehi = {std::max(ehi.x, q.x), std::max(ehi.y, q.y), std::max(ehi.z, q.z)};
        }
    }
    if (files < 2) {
        std::printf("%s\nparts %zu prims %zu\n", log.c_str(), t.parts.size(), t.source.prims.size());
        for (const auto& f : fs::recursive_directory_iterator(o.folder)) std::printf("  %s\n", f.path().string().c_str());
    }
    CHECK(files >= 2);
    const double tol = 1e-8;  // float32 at 3 mm in metres
    CHECK_NEAR(elo.x, lo.x, tol);
    CHECK_NEAR(elo.y, lo.y, tol);
    CHECK_NEAR(elo.z, lo.z, tol);
    CHECK_NEAR(ehi.x, hi.x, tol);
    CHECK_NEAR(ehi.y, hi.y, tol);
    CHECK_NEAR(ehi.z, hi.z, tol);
    CHECK_NEAR(ehi.y - elo.y, 0.003, tol);
    fs::remove_all(dir);
}

// The export's object of a primitive is its source node, not a name: two nodes sharing one
// mesh are two objects (two groups), and an unnamed node with an unnamed mesh still gets one.
TEST(glb_object_identity) {
    const std::string dir = std::filesystem::temp_directory_path().string() + "/knit_s6_identity";
    const std::string path = writeTinyGltf(dir);
    {
        std::ifstream in(path);
        std::string j((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        auto replace = [&](const std::string& from, const std::string& to) {
            const size_t at = j.find(from);
            CHECK(at != std::string::npos);
            j.replace(at, from.size(), to);
        };
        const std::string prim = R"({"attributes":{"POSITION":0,"TEXCOORD_0":1},"indices":2,"material":0})";
        replace(R"("nodes":[0,2]})", R"("nodes":[0,2,3]})");
        replace(R"("mesh":0}],)", R"("mesh":0},{"translation":[9,0,0],"mesh":1}],)");
        replace(R"("primitives":[)" + prim + "]}]", R"("primitives":[)" + prim + R"(]},{"primitives":[)" + prim + "]}]");
        std::ofstream out(path);
        out << j;
    }
    SceneEntry e;
    e.glb.sizeMm = 10;
    e.glb.knitAllMaterials = true;
    SourceScene g;
    std::string err;
    CHECK(loadGlb(path, e.glb, g, &err));
    CHECK(g.prims.size() == 3);
    if (g.prims.size() != 3) return;
    // Scene order: the parent's child (node 1), the mirror (node 2), the unnamed node (3).
    CHECK(g.prims[0].nodeIndex == 1 && g.prims[1].nodeIndex == 2 && g.prims[2].nodeIndex == 3);
    CHECK(g.prims[0].mesh == "tri" && g.prims[1].mesh == "tri" && g.prims[2].mesh.empty() && g.prims[2].node.empty());
    const Placement a = exportObjectOf(g, 0), b = exportObjectOf(g, 1), c = exportObjectOf(g, 2, ExportScene::kRandomGroup);
    CHECK(a.object != b.object && b.object != c.object && a.object != c.object);
    CHECK(a.objectName == "child" && b.objectName == "mirror" && c.objectName == "node_3" && c.group == "fibres_random");
    // Exported: three object groups (and the fibre group under the unnamed one).
    ExportScene ex(dir + "/export");
    std::filesystem::remove_all(ex.folder());
    std::filesystem::create_directories(ex.folder() + "/meshes");
    TriMesh tri;
    tri.pos = {{0, 0, 0}, {0.01, 0, 0}, {0, 0.01, 0}};
    tri.tri = {{0, 1, 2}};
    for (int k = 0; k < 3; ++k) {
        const std::string name = "p" + std::to_string(k);
        CHECK(ex.addMesh(name, "wool_white", tri, &err));
        ex.place(name, exportObjectOf(g, k));
    }
    StrandSetF s;
    s.add({0, 0, 0, 0.001f, 0, 0}, {2e-5f, 2e-5f});
    CHECK(ex.addStrands("p2_fibres_random", ex.randomFibreMaterial("wool_white"), s, &err));
    ex.place("p2_fibres_random", c);
    CHECK(writeSceneJson(ex, &err));
    std::ifstream in(ex.folder() + "/scene.json");
    const std::string j((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    for (const char* group : {R"({"id": "g0001", "name": "knitted_bears", "origin")",
                              R"({"id": "g0002", "name": "child", "parent_id": "g0001")",
                              R"({"id": "g0003", "name": "mirror", "parent_id": "g0001")",
                              R"({"id": "g0004", "name": "node_3", "parent_id": "g0001")",
                              R"({"id": "g0005", "name": "fibres_random", "parent_id": "g0004")",
                              R"({"id": "g0006", "name": "Ground", "parent_id": "g0001")"})
        CHECK(j.find(group) != std::string::npos);
    CHECK(j.find("\"next_group_seq\": 7,") != std::string::npos);
}

// --------------------------------------------------------- map set cache ----
// (from knit_tests_s6, its CTest name kept: attaching a set to a package is the export's)

namespace {

std::string mapFileBytes(const std::string& p) {
    std::ifstream f(p, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(f), {});
}

PlyMapSpec smallMapSpec() {
    PlyMapSpec s;
    s.g.count = 3;
    s.size = 32;
    s.pp.aoDirections = 16;
    s.yarnD = 1.1;
    return s;
}

} // namespace
// The baked maps' cache (scene/appearance/baked_maps.h): a miss builds the body on the given surface and publishes
// one atlas per yarn, reporting a yarn at a time; a hit returns it; the key follows the geometry (surface, body,
// seed, texel), not the core proxy; the entry's plans fit the tubes the export meshes; a package gets its five
// material maps from an entry; a cancelled bake leaves no entry, a finish-once bake publishes regardless.
TEST(baked_map_set_cache) {
    namespace fs = std::filesystem;
    const std::string root = fs::temp_directory_path().string() + "/knit_s6_bakes";
    fs::remove_all(root);
    const CurveSet cs = makeSet({makeYarn(lineCtrl({0, 0, 0}, {12, 0, 0}, 24), 0), makeYarn(lineCtrl({0, 3, 0}, {8, 3, 0}, 16), 1)});
    PatternParams pp;
    pp.yarnDiameter = 1.5;
    pp.clearance = 0.1;
    SurfaceParams sp;
    sp.compression = false;
    sp.ringSpacing = 0.05;
    const SurfaceModel sm = buildSurface(cs, pp, sp);
    BakedMapSpec spec;
    spec.surfaceDigest = surfaceBuildDigest(cs, pp, sp);
    spec.body.plies = 2;
    spec.body.fibres = 40;
    spec.body.irregularity = 0;
    spec.partSeed = 7;
    spec.texel = 0.05;
    spec.atlasWidth = 256;
    {
        BakedMapSpec t = spec;
        t.texel = 0.08;
        CHECK(bakedMapId(t) != bakedMapId(spec));
        t = spec;
        t.body.fibres = 41;
        CHECK(bakedMapId(t) != bakedMapId(spec));
        t = spec;
        t.partSeed = 8;
        CHECK(bakedMapId(t) != bakedMapId(spec));
        t = spec;
        t.body.core = BodyCore::Opaque;
        CHECK(bakedMapId(t) == bakedMapId(spec));  // the core proxy is no input of the maps
    }
    BakedMapEntry a, b;
    std::string err;
    int reports = 0;
    BakeProgress progress;
    progress.report = [&](int, int) { ++reports; };
    const bool baked = ensureBakedMapSet(root, spec, cs, sm, BakeInterrupt::Cancellable, progress, a, &err);
    if (!baked) std::printf("    %s\n", err.c_str());
    CHECK(baked && !a.hit && a.yarns.size() == 2 && reports == 2 && isPublished(a.dir, bakedMapKey(spec)));
    std::printf("    baked 2 yarns in %.0f ms: %d x %d and %d x %d\n", a.ms, a.yarns[0].width, a.yarns[0].height, a.yarns[1].width, a.yarns[1].height);
    for (int y = 0; y < 2; ++y) CHECK(bakedPlanMatches(a.yarns[y], sm.yarns[y].length, sm.yarns[y].radius, spec));
    CHECK(!bakedPlanMatches(a.yarns[0], 2 * sm.yarns[0].length, sm.yarns[0].radius, spec));
    CHECK(ensureBakedMapSet(root, spec, cs, sm, BakeInterrupt::Cancellable, BakeProgress{}, b, &err) && b.hit && b.dir == a.dir && b.yarns.size() == 2);
    // A package's material maps from the entry, five files of the atlas's size.
    const std::string pkg = root + "/package";
    BakedYarnFiles files;
    CHECK(attachBakedYarn(a, 1, pkg, "textures/baked_y001", 0.6, 0.3, files, &err));
    for (const std::string& f : {files.cavity, files.normal, files.specular, files.roughness, files.height}) {
        int w = 0, h = 0;
        std::vector<uint8_t> px;
        CHECK(readPngRgb8(pkg + "/" + f, w, h, px) && w == a.yarns[1].width && h == a.yarns[1].height);
    }
    CHECK(!attachBakedYarn(a, 2, pkg, "textures/x", 0.6, 0.3, files, &err));
    // Cancelled: no entry. Finish once: the entry, the cancelled token notwithstanding.
    fs::remove_all(root);
    CancelToken token;
    token.cancel();
    const ThreadCancelScope scope(&token);
    bool cancelled = false;
    try {
        (void)ensureBakedMapSet(root, spec, cs, sm, BakeInterrupt::Cancellable, BakeProgress{}, a, &err);
    } catch (const Cancelled&) {
        cancelled = true;
    }
    CHECK(cancelled && !isPublished(root + "/" + bakedMapId(spec), bakedMapKey(spec)));
    CHECK(ensureBakedMapSet(root, spec, cs, sm, BakeInterrupt::FinishOnce, BakeProgress{}, b, &err) && !b.hit && isPublished(b.dir, bakedMapKey(spec)));
}


// The cache: generated once, then a verified hit; the key covers the D / pitch ratio but not
// sigma; a corrupted file is not a hit (regenerated); a cancelled generation publishes
// nothing; an attached set equals setupPlyTextures' output for the same spec.
TEST(ply_map_set_cache) {
    namespace fs = std::filesystem;
    const std::string root = fs::temp_directory_path().string() + "/knit_s6_plymaps";
    fs::remove_all(root);
    const PlyMapSpec spec = smallMapSpec();
    {
        PlyMapSpec t = spec;
        t.sigma = 3;
        CHECK(plyMapSetId(t) == plyMapSetId(spec));  // bump_scale is per scene
        t = spec;
        t.yarnD = 1.0;
        CHECK(plyMapSetId(t) != plyMapSetId(spec));
        t = spec;
        t.g.pitch = 3.0;
        CHECK(plyMapSetId(t) != plyMapSetId(spec));
        t = spec;
        t.pp.model = PlyProfileModel::Circles;
        CHECK(plyMapSetId(t) != plyMapSetId(spec));
    }
    PlyMapSetResult a, b;
    std::string err;
    CHECK(ensurePlyMapSet(root, spec, a, &err));
    CHECK(!a.hit && isPublished(a.dir, plyMapSetKey([&] {
              PlyMapSpec c = spec;
              c.sigma = 1;
              return c;
          }())));
    CHECK(ensurePlyMapSet(root, spec, b, &err));
    CHECK(b.hit && b.dir == a.dir);
    // Corrupt one map: no hit, regenerated byte-identical.
    const std::string png = a.dir + "/" + a.set.normal;
    const std::string good = mapFileBytes(png);
    std::ofstream(png, std::ios::binary) << "garbage";
    PlyMapSetResult c;
    CHECK(ensurePlyMapSet(root, spec, c, &err));
    CHECK(!c.hit && mapFileBytes(png) == good);
    // Cancelled generation of another key: nothing published, no staging left.
    {
        PlyMapSpec t = spec;
        t.g.count = 2;
        CancelToken tok;
        tok.cancel();
        bool threw = false;
        {
            ThreadCancelScope scope(&tok);
            PlyMapSetResult r;
            try {
                ensurePlyMapSet(root, t, r, &err);
            } catch (const Cancelled&) {
                threw = true;
            }
        }
        CHECK(threw);
        PlyMapSpec tc = t;
        tc.sigma = 1;
        CHECK(!isPublished(root + "/" + plyMapSetId(tc), plyMapSetKey(tc)));
        int staging = 0;
        for (const auto& d : fs::directory_iterator(root)) staging += d.path().filename().string().find(".staging-") != std::string::npos;
        CHECK(staging == 0);
    }
    // A package holds the maps its material can use, byte for byte what the generator writes (writePlyMaps
    // called directly: the oracle): height, depth, both normal maps and ply_maps.json; never the AO map, and the
    // cavity map only when the material uses it. Both paths into a package (attachPlySet, setupPlyTextures).
    const std::string e1 = root + "_exp1", e2 = root + "_exp2", e3 = root + "_oracle", e4 = root + "_cavity";
    for (const std::string& d : {e1, e2, e3, e4}) fs::remove_all(d);
    PlyModeOptions o;
    o.spec = spec;
    PlyMapSpec atScene = spec;
    atScene.sigma = 1000.0;
    PlyMapSet oracle;
    CHECK(writePlyMaps(e3, atScene, oracle, &err));
    const std::optional<PlyShading> x1 = attachPlySet(b, o, e1, 1000.0, &err), x2 = setupPlyTextures(o, e2, 1000.0, &err);
    CHECK(x1 && x2);
    for (const std::string& pkg : {e1, e2}) {
        for (const std::string& f : {b.set.height, b.set.depth, b.set.normal, b.set.normalFlipGreen, std::string("ply_maps.json")})
            CHECK(!mapFileBytes(e3 + "/" + f).empty() && mapFileBytes(pkg + "/textures/" + f) == mapFileBytes(e3 + "/" + f));
        CHECK(!fs::exists(pkg + "/textures/" + b.set.ao) && !fs::exists(pkg + "/textures/" + b.set.cavity));
    }
    CHECK(x1->on && x1->bumpScale == x2->bumpScale && x1->uvRepeat == x2->uvRepeat && x1->normal == x2->normal);
    CHECK(x1->cavity.empty() && x2->cavity.empty());
    o.cavity = true;
    const std::optional<PlyShading> x4 = attachPlySet(b, o, e4, 1000.0, &err);
    CHECK(x4 && !x4->cavity.empty() && mapFileBytes(e4 + "/textures/" + b.set.cavity) == mapFileBytes(e3 + "/" + b.set.cavity));
    CHECK(!fs::exists(e4 + "/textures/" + b.set.ao));
    for (const std::string& d : {root, e1, e2, e3, e4}) fs::remove_all(d);
}
