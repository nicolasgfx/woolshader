#pragma once
// Yarns and strands the test programs build by hand (knit_tests, knit_tests_export, knit_tests_s6, knit_tests_gpu):
// straight and fixture yarns, the chart fixture piece and the exact text of its placement, a scratch folder, two
// tiny strands.
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "fibres/fibre.h"
#include "fibres/strands_file.h"
#include "knitting/knit/piece_knit.h"
#include "knitting/layout/chart.h"
#include "knitting/stitch/hand_tension.h"
#include "knitting/tube/param_text.h"
#include "knitting/tube/yarn_surface.h"
#include "model/mesh.h"
#include "model/params.h"
#include "model/yarn_path.h"

namespace knit {

inline YarnPath makeYarn(std::vector<Vec3> ctrl, int id) {
    YarnPath y;
    y.id = id;
    y.course = id;
    y.name = "test " + std::to_string(id);
    y.curve.ctrl = std::move(ctrl);
    y.stitches.push_back({0, id, 0.0, y.curve.uMax()});
    y.finalize();
    return y;
}

inline CurveSet makeSet(std::vector<YarnPath> ys) {
    CurveSet cs;
    cs.yarns = std::move(ys);
    cs.wales = 1;
    cs.courses = static_cast<int>(cs.yarns.size());
    fillFlatTopology(cs);
    return cs;
}

// Straight line from a to b with n spans (collinear, evenly spaced controls).
inline std::vector<Vec3> lineCtrl(Vec3 a, Vec3 b, int n) {
    std::vector<Vec3> c;
    Vec3 step = (b - a) / n;
    for (int i = -1; i <= n + 1; ++i) c.push_back(a + step * i);
    return c;
}

// The chart fixture (a distorted UV sphere, exported like the scene pieces) scaled to about
// 30 mm radius: a small closed piece with two poles.
struct FixturePiece {
    TriMesh mesh;
    std::unique_ptr<Chart> chart;
    std::unique_ptr<PieceSurface> surf;
};

inline FixturePiece fixturePiece(double scale) {
    FixturePiece f;
    std::string err;
    if (!readObj(std::string(KNIT_TEST_DATA) + "/chart_fixture.obj", f.mesh, &err)) throw std::runtime_error(err);
    f.mesh.transform(scale, f.mesh.centroid());
    f.mesh.computeNormals();
    f.chart = std::make_unique<Chart>(f.mesh, 5);
    f.surf = std::make_unique<PieceSurface>(*f.chart);
    return f;
}

// A number as exact text: a hex float, so equal text means equal bits.
inline void appendExact(std::string& text, double v) {
    char b[40];
    std::snprintf(b, sizeof b, "%a ", v);
    text += b;
}

inline void appendNode(std::string& text, const LayoutNode& n) {
    for (const Vec3& v : {n.p, n.n, n.c, n.v})
        for (double x : {v.x, v.y, v.z}) appendExact(text, x);
    for (double x : {n.w, n.h, n.wL, n.wR, n.a, n.b, n.f}) appendExact(text, x);
}

// Everything a spiral piece's placement and hand produce before settle, as exact text: every layout node and
// virtual node, the controls and stitch spans (curvesDigest), the hand's records and every yarn's thickness.
// The placement pins hash it.
inline std::string placementText(const PieceKnit& pk, const std::vector<StitchTension>& tensions) {
    std::string text = "nodes\n";
    for (const LayoutNode& n : pk.lay.nodes) appendNode(text, n);
    appendNode(text, pk.lay.before);
    appendNode(text, pk.lay.after);
    text += "\nvirtual\n";
    for (const LayoutNode& n : pk.lay.virtualNodes) appendNode(text, n);
    char b[64];
    std::snprintf(b, sizeof b, "\ncurves %016llx\nhand\n", static_cast<unsigned long long>(curvesDigest(pk.curves)));
    text += b;
    for (const StitchTension& t : tensions)
        for (float x : {t.tension, t.lean, t.raise, t.spot}) appendExact(text, x);
    text += "\nthickness\n";
    for (const YarnPath& y : pk.curves.yarns)
        for (float x : y.thickness) appendExact(text, x);
    return text;
}

// A closed octahedron of radius r around c, faces outward: a cover whose winding number is 1 at c.
inline TriMesh octahedronAround(const Vec3& c, double r) {
    TriMesh m;
    m.pos = {c + Vec3{r, 0, 0}, c + Vec3{-r, 0, 0}, c + Vec3{0, r, 0}, c + Vec3{0, -r, 0}, c + Vec3{0, 0, r}, c + Vec3{0, 0, -r}};
    m.tri = {{0, 2, 4}, {2, 1, 4}, {1, 3, 4}, {3, 0, 4}, {2, 0, 5}, {1, 2, 5}, {3, 1, 5}, {0, 3, 5}};
    m.computeNormals();
    return m;
}

// A hand with every channel on (spiral pieces ignore the column shift: it slides lattice columns only).
inline HandRequest everyChannelHand() {
    HandRequest r;
    r.settings.unevenness = 12;
    r.settings.columnShift = 0.4;
    r.settings.loopRate = 2;
    r.settings.loopHeight = 0.6;
    r.settings.waveHeight = 0.3;
    r.settings.spotRate = 1.5;
    r.settings.spotChange = 25;
    r.seed = 0x5eedf00dull;
    return r;
}

inline std::string scratchDir(const char* name) {
    std::filesystem::path d = std::filesystem::temp_directory_path() / "knit_tests" / name;
    std::filesystem::remove_all(d);
    std::filesystem::create_directories(d / "meshes");
    return d.string();
}

// Two tapered strands: a gentle arc (6 vertices) and a short straight one (3).
inline StrandSetF tinyStrands(double scale = 1.0) {
    StrandSetF s;
    std::vector<float> p, r;
    for (int i = 0; i < 6; ++i) {
        const double t = i / 5.0;
        p.insert(p.end(), {float(scale * 0.01 * t), float(scale * (0.002 * std::sin(t))), 0.f});
        r.push_back(float(scale * 2e-4 * (1 - 0.5 * t)));
    }
    s.add(p, r);
    p.clear();
    r.clear();
    for (int i = 0; i < 3; ++i) {
        p.insert(p.end(), {0.f, float(scale * 0.004 * i), float(scale * 0.001)});
        r.push_back(float(scale * 1.5e-4 * (1 - 0.2 * i)));
    }
    s.add(p, r);
    return s;
}

// Straight round yarns along +x at the given y (fibre front +z, flat frame), D = 1.
struct StraightSet {
    CurveSet cs;
    PatternParams pp;
    SurfaceModel sm;
};

inline StraightSet straightSet(std::vector<double> ys, std::vector<double> lengths = {}) {
    StraightSet s;
    std::vector<YarnPath> v;
    for (size_t k = 0; k < ys.size(); ++k) {
        const double L = k < lengths.size() ? lengths[k] : 20.0;
        v.push_back(makeYarn(lineCtrl({0, ys[k], 0}, {L, ys[k], 0}, static_cast<int>(std::ceil(L))), static_cast<int>(k)));
    }
    s.cs = makeSet(std::move(v));
    SurfaceParams sp;
    sp.compression = false;
    s.sm = buildSurface(s.cs, s.pp, sp);
    return s;
}

inline FibreParams mixParams(double frac, double density) {
    FibreParams fp = pieceFibres();
    fp.randomFraction = frac;
    fp.density = density;
    fp.rootDepth = 0.126;
    return fp;
}

} // namespace knit
