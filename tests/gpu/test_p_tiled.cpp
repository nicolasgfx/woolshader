// GPU port tests: the tiled pipeline's device parts (lattice node classification, mesh surface
// continuations, boundary checks B1 / B2, gate (d), the knitting numerics) against the CPU on
// synthetic meshes; deviations are printed, the checks are structural.
#include <chrono>
#include <cmath>
#include <cstdio>
#include <functional>

#include "gpu_common.h"
#include "../cpu_numerics.h"
#include "model/mesh_check.h"
#include "model/mesh_query.h"
#include "knitting/knit/mesh_knit.h"
#include "gpu/knitting.h"

using namespace knit;

namespace {

TriMesh gridMeshP(int nx, int ny, const std::function<Vec3(double, double)>& pos, const std::function<Vec2(double, double)>& uv) {
    TriMesh m;
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i) m.pos.push_back(pos(i, j));
    auto id = [&](int i, int j) { return j * (nx + 1) + i; };
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            const int a = id(i, j), b = id(i + 1, j), d = id(i + 1, j + 1), e = id(i, j + 1);
            auto corner = [&](int v) { return uv(v % (nx + 1), v / (nx + 1)); };
            m.tri.push_back({a, b, d});
            m.uv.push_back({corner(a), corner(b), corner(d)});
            m.tri.push_back({a, d, e});
            m.uv.push_back({corner(a), corner(d), corner(e)});
        }
    m.computeNormals();
    return m;
}

TriMesh cylinder() {
    const double R = 20, H = 48;
    const int na = 160;
    const TriMesh raw = gridMeshP(na, 48,
                                  [&](double x, double y) {
                                      const double a = 2 * kPi * x / na;
                                      return Vec3{R * std::cos(a), R * std::sin(a), y * H / 48};
                                  },
                                  [&](double x, double y) { return Vec2{2 * kPi * R * x / na / 100, y * H / 48 / 100}; });
    return weldVertices(raw, 1e-9).mesh;
}

// A plane in two charts (a uv seam at x = 40, the right half's uv shifted far away).
TriMesh twoCharts() {
    TriMesh m = gridMeshP(80, 50, [](double x, double y) { return Vec3{x, y, 0}; },
                          [](double x, double y) { return Vec2{(x <= 40 ? x : x + 1000) / 100, y / 100}; });
    for (size_t t = 0; t < m.tri.size(); ++t) {
        double cx = 0;
        for (int c = 0; c < 3; ++c) cx += m.pos[m.tri[t][c]].x / 3;
        if (cx > 40)
            for (int c = 0; c < 3; ++c) {
                const Vec3 p = m.pos[m.tri[t][c]];
                m.uv[t][c] = {(p.x + 1000) / 100, p.y / 100};
            }
    }
    return m;
}

// A closed axis-aligned box with outward normals: a solid for the lattice's inside test.
TriMesh closedBox(const Vec3& lo, const Vec3& hi) {
    TriMesh m;
    for (int k = 0; k < 8; ++k) m.pos.push_back({k & 1 ? hi.x : lo.x, k & 2 ? hi.y : lo.y, k & 4 ? hi.z : lo.z});
    m.tri = {{0, 2, 1}, {1, 2, 3}, {4, 5, 6}, {5, 7, 6}, {0, 1, 4}, {1, 5, 4}, {2, 6, 3}, {3, 6, 7}, {0, 4, 2}, {2, 4, 6}, {1, 3, 5}, {3, 7, 5}};
    m.computeNormals();
    return m;
}

double msSince(std::chrono::steady_clock::time_point t) { return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t).count(); }

} // namespace

TEST(tiled_gpu_vs_cpu) {
    REQUIRE_GPU();
    struct Case {
        const char* name;
        TriMesh m;
    };
    for (const Case& c : {Case{"cylinder", cylinder()}, Case{"two charts", twoCharts()}}) {
        MeshKnitParams p;
        p.tile.pp = PatternParams{};
        p.tile.g.w = p.tile.pp.waleSpacing;
        p.tile.g.h = p.tile.pp.courseSpacing;
        p.surfaceGate = true;
        PartKnit a, b;
        std::string err;
        auto t0 = std::chrono::steady_clock::now();
        CHECK(knitMesh(c.m, p, knit_test::cpuNumerics(), a, &err));
        const double msCpu = msSince(t0);
        std::unique_ptr<KnitNumerics> num = gpu::makeKnitNumerics(&err);
        CHECK(num != nullptr);
        t0 = std::chrono::steady_clock::now();
        CHECK(knitMesh(c.m, p, *num, b, &err));
        const double msGpu = msSince(t0);
        int sameCharts = 0;
        for (size_t k = 0; k < a.t.charts.size() && k < b.t.charts.size(); ++k) {
            const ChartTiling &x = a.t.charts[k], &y = b.t.charts[k];
            sameCharts += x.candidates == y.candidates && x.kept == y.kept && x.outside == y.outside && x.ambiguous == y.ambiguous &&
                          x.boundary == y.boundary && x.insideOther == y.insideOther && x.runEnd == y.runEnd && x.runs == y.runs;
        }
        std::printf("  %s: charts with equal lattice counts %d / %zu; stitches %zu / %zu; gate a %s / %s (settle %d / %d rounds); "
                    "B1 %s / %s (min %.4f / %.4f mm), B2 %s / %s; gate d %s / %s; verdict %s / %s; CPU %.0f ms, GPU %.0f ms; %s\n",
                    c.name, sameCharts, a.t.charts.size(), a.t.prog.stitches.size(), b.t.prog.stitches.size(), toString(a.after.verdict),
                    toString(b.after.verdict), a.settle.rounds, b.settle.rounds, toString(a.boundary.b1), toString(b.boundary.b1),
                    a.boundary.minDistance, b.boundary.minDistance, toString(a.boundary.b2), toString(b.boundary.b2), toString(a.gateDStatus),
                    toString(b.gateDStatus), toString(a.verdict), toString(b.verdict), msCpu, msGpu, gpu::knitPlacementStats(*num).c_str());
        CHECK(sameCharts == static_cast<int>(a.t.charts.size()));
        CHECK(a.t.prog.stitches.size() == b.t.prog.stitches.size());
    }
}

// The port: a small plane knitted through the CPU double and through the device numerics gives the
// same stitches, the same yarn lengths and the same verdicts, without and with a solid to keep
// stitches out of (in another frame: the offset must be applied). Settle is off: nothing here
// depends on the solver, so a difference is a difference in what the port's methods compute.
TEST(knit_mesh_port_cpu_gpu) {
    REQUIRE_GPU();
    const TriMesh plane = gridMeshP(60, 40, [](double x, double y) { return Vec3{x, y, 0}; }, [](double x, double y) { return Vec2{x / 100, y / 100}; });
    const Vec3 offset{100, 0, 0};
    // x 20..35 in the plane's frame; no lattice node on a face, where float and double may differ.
    const TriMesh box = closedBox({120, 12, -4}, {135, 25, 4});
    const MeshQuery boxQuery(box);
    LatticeSolids none, solid;
    solid.meshes = {&boxQuery};
    solid.offset = offset;
    CHECK(solid.contains({27, 18, 0}) && !solid.contains({127, 18, 0}) && !none.contains({27, 18, 0}));
    std::string err;
    std::unique_ptr<KnitNumerics> device = gpu::makeKnitNumerics(&err);
    CHECK(device != nullptr);
    for (const LatticeSolids* solids : {&none, &solid}) {
        MeshKnitParams p;
        p.tile.pp = PatternParams{};
        p.tile.g.w = p.tile.pp.waleSpacing;
        p.tile.g.h = p.tile.pp.courseSpacing;
        p.tile.solids = *solids;
        p.settle = false;
        p.surfaceGate = true;
        PartKnit a, b;
        CHECK(knitMesh(plane, p, knit_test::cpuNumerics(), a, &err));
        CHECK(knitMesh(plane, p, *device, b, &err));
        std::printf("  %zu solids: stitches %zu / %zu, inside %d / %d, verdict %s / %s\n", solids->meshes.size(), a.t.prog.stitches.size(),
                    b.t.prog.stitches.size(), a.t.charts[0].insideOther, b.t.charts[0].insideOther, toString(a.verdict), toString(b.verdict));
        CHECK(a.t.prog.stitches.size() == b.t.prog.stitches.size());
        CHECK(a.t.charts.size() == 1 && b.t.charts.size() == 1);
        CHECK(a.t.charts[0].kept == b.t.charts[0].kept && a.t.charts[0].insideOther == b.t.charts[0].insideOther);
        CHECK((a.t.charts[0].insideOther > 0) == !solids->meshes.empty());
        CHECK(a.verdict == b.verdict && a.after.verdict == b.after.verdict);
        CHECK(a.boundary.b1 == b.boundary.b1 && a.boundary.b2 == b.boundary.b2 && a.gateDStatus == b.gateDStatus);
        CHECK(a.verdict == Status::Pass);
        CHECK(a.t.curves.yarns.size() == b.t.curves.yarns.size());
        double worst = 0;
        for (size_t y = 0; y < a.t.curves.yarns.size() && y < b.t.curves.yarns.size(); ++y)
            worst = std::max(worst, std::fabs(a.t.curves.yarns[y].length() - b.t.curves.yarns[y].length()));
        std::printf("  largest difference in a yarn's length %.2e mm\n", worst);
        CHECK(worst < 1e-3);
    }
}
