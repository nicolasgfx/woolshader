// GPU tests: every public entry point survives a device failure at any of its checks (plan D10). A test hook
// (runtime.h failCheckForTests) makes the n-th next check of this thread throw DeviceError("<what>: injected
// failure"). After each injected failure: the public failure contract (false + the reason, or the port's
// exception with it), the tracked device memory back at its baseline, no lease held or queued, the abort flag
// clear, no allocation cap hit; then the uninjected call succeeds with the reference's exact results.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <functional>
#include <random>
#include <string>
#include <vector>

#include "gpu_common.h"
#include "gpu/device/runtime.h"
#include "gpu/fibre_generator.h"
#include "gpu/knitting.h"
#include "gpu/surfaces.h"
#include "fibres/fibre_generator.h"
#include "knitting/knit/fixtures.h"
#include "knitting/knit/mesh_knit.h"
#include "knitting/stitch/stockinette.h"
#include "knitting/tube/plies.h"
#include "knitting/tube/yarn_surface.h"
#include "scene/knit/knit_cache.h"

using namespace knit;

namespace {

// One entry point: runs it, returns success and the reason, and appends its results as numbers.
using EntryPoint = std::function<bool(std::string& err, std::vector<double>& results)>;

void addVec(std::vector<double>& r, const Vec3& v) {
    r.push_back(v.x);
    r.push_back(v.y);
    r.push_back(v.z);
}

// The checks of `what` at up to 12 positions (the first, the last and evenly spread): each injected failure is
// reported and leaves nothing behind; the next call gives the reference results (`exact`: compared exactly).
void injectEverywhere(const char* what, const EntryPoint& run, bool exact) {
    gpu::releaseResident();
    const size_t baseline = gpu::memStats().current;
    std::string err;
    std::vector<double> reference;
    const long long before = gpu::checksForTests();
    CHECK(run(err, reference));
    const long long checks = gpu::checksForTests() - before;
    gpu::releaseResident();
    CHECK(gpu::memStats().current == baseline);
    std::vector<long long> at;
    for (int k = 0; k < 12 && checks > 0; ++k) at.push_back(1 + (checks - 1) * k / 11);
    at.erase(std::unique(at.begin(), at.end()), at.end());
    std::printf("    %-28s %lld checks; injected at %zu of them\n", what, checks, at.size());
    CHECK(checks > 0);
    for (long long n : at) {
        std::vector<double> results;
        err.clear();
        const long long first = gpu::checksForTests();
        gpu::failCheckForTests(n);
        const bool ok = run(err, results);
        gpu::failCheckForTests(0);
        const bool reached = gpu::checksForTests() - first >= n;
        const bool reported = !ok && err.find("injected failure") != std::string::npos;
        gpu::releaseResident();
        const bool released = gpu::memStats().current == baseline;
        const bool idle = gpu::leaseDepth() == 0 && gpu::queuedRequests() == 0 && gpu::abortFlagValue() == 0 && !gpu::AllocationCap::capHit();
        results.clear();
        const bool again = run(err, results);
        gpu::releaseResident();
        const bool same = !exact || results == reference;
        if (!reached || !reported || !released || !idle || !again || !same)
            std::printf("      at check %lld: reached %d, reported %d, memory released %d (%zu vs %zu bytes), idle %d, rerun %d, same results %d (%s)\n", n,
                        reached, reported, released, gpu::memStats().current, baseline, idle, again, same, err.c_str());
        CHECK(reached);
        CHECK(reported);
        CHECK(released);
        CHECK(idle);
        CHECK(again);
        CHECK(same);
    }
}

// The port's calls report a device failure by exception.
bool portCall(const std::function<void(std::vector<double>&)>& call, std::string& err, std::vector<double>& results) {
    try {
        call(results);
        return true;
    } catch (const std::exception& e) {
        err = e.what();
        return false;
    }
}

TriMesh plane(int nx, int ny) {
    TriMesh m;
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i) m.pos.push_back({double(i), double(j), 0});
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            const int a = j * (nx + 1) + i, b = a + 1, c = a + nx + 2, d = a + nx + 1;
            m.tri.push_back({a, b, c});
            m.uv.push_back({Vec2{i / 100.0, j / 100.0}, Vec2{(i + 1) / 100.0, j / 100.0}, Vec2{(i + 1) / 100.0, (j + 1) / 100.0}});
            m.tri.push_back({a, c, d});
            m.uv.push_back({Vec2{i / 100.0, j / 100.0}, Vec2{(i + 1) / 100.0, (j + 1) / 100.0}, Vec2{i / 100.0, (j + 1) / 100.0}});
        }
    m.computeNormals();
    return m;
}

YarnPath strand(double y, int id) {
    YarnPath p;
    p.id = id;
    p.course = id;
    for (int i = -1; i <= 11; ++i) p.curve.ctrl.push_back({1.0 * i, y, 0});
    p.stitches.push_back({0, id, 0.0, p.curve.uMax()});
    p.finalize();
    return p;
}

} // namespace

TEST(gpu_injected_failures_release_everything) {
    REQUIRE_GPU();
    // The same checks in every run: launch sizes from the pacer's growth cap, not from measured times (a loaded GPU
    // would otherwise split a kernel into more or fewer launches than the reference run).
    _putenv_s("KNIT_GPU_LAUNCH_MS", "1e9");
    const PatternParams pp = findPreset("tightest")->params;
    const CurveSet patch = generateStockinette(pp);
    const SurfaceModel cpuModel = buildSurface(patch, pp, SurfaceParams{});

    injectEverywhere("buildSurface", [&](std::string& err, std::vector<double>& r) {
        SurfaceModel sm;
        if (!gpu::buildSurface({patch, pp, SurfaceParams{}}, sm, nullptr, &err)) return false;
        for (const YarnSurface& ys : sm.yarns)
            for (const Ring& R : ys.rings) {
                addVec(r, R.c);
                addVec(r, R.n1);
                r.push_back(R.a);
                r.push_back(R.b);
            }
        return true;
    }, true);
    injectEverywhere("checkSurface", [&](std::string& err, std::vector<double>& r) {
        SurfaceCheck check;
        if (!gpu::checkSurface({patch, pp, SurfaceParams{}}, {}, check, nullptr, &err)) return false;
        r = {check.minFoldRatio, double(check.foldRings), double(check.pointsTested), double(check.penetratingPoints), check.maxDepth};
        return true;
    }, true);
    injectEverywhere("buildTubeMesh", [&](std::string& err, std::vector<double>& r) {
        TubeMesh m;
        if (!gpu::buildTubeMesh({patch, pp, SurfaceParams{}}, TubeParams{}, m, nullptr, &err)) return false;
        for (const Vec3& p : m.pos) addVec(r, p);
        for (uint32_t t : m.tri) r.push_back(t);
        return true;
    }, true);
    injectEverywhere("buildPlies", [&](std::string& err, std::vector<double>& r) {
        PlyGeometry g;
        g.count = 3;
        PlySet plies;
        if (!gpu::buildPlies({patch, cpuModel, g, pp.radius()}, plies, nullptr, &err)) return false;
        for (const YarnPath& y : plies.curves.yarns)
            for (const Vec3& c : y.curve.ctrl) addVec(r, c);
        r.push_back(plies.containment);
        return true;
    }, true);
    injectEverywhere("checkSceneClearance", [&](std::string& err, std::vector<double>& r) {
        std::mt19937 random(5);
        std::uniform_real_distribution<double> position(0, 20), step(-1, 1);
        std::vector<gpu::ClearanceSegment> segments;
        for (int i = 0; i < 2000; ++i) {
            const Vec3 a{position(random), position(random), position(random)};
            segments.push_back({a, a + Vec3{step(random), step(random), step(random)}, 0.5, i % 5, i % 5});
        }
        gpu::SceneClearance out;
        if (!gpu::checkSceneClearance({segments, 0.1}, out, &err)) return false;
        for (const gpu::ClosePartPair& p : out.pairs) r.insert(r.end(), {double(p.partA), double(p.partB), double(p.segments), p.minDistance});
        r.push_back(out.minDistance);
        return true;
    }, true);
    injectEverywhere("generateFibres", [&](std::string& err, std::vector<double>& r) {
        FibreParams fp = pieceFibres();
        fp.density = 5;
        fp.lightChecks = true;
        FibreSet set;
        gpu::FibreRunReport rep;
        if (!gpu::generateFibres({patch, cpuModel, pp, fp, gpu::frameFrom(FibreFrame{})}, set, &rep, nullptr, &err)) return false;
        for (const Vec3f& p : set.points) r.insert(r.end(), {p.x, p.y, p.z});
        r.push_back(double(set.fibres.size()));
        return true;
    }, true);

    std::string err;
    std::unique_ptr<KnitNumerics> port = gpu::makeKnitNumerics(&err);
    CHECK(port != nullptr);
    if (!port) return;
    // Gate (a) on a knitted teddy piece (the port refuses flat patches), when the fixture is there.
    KnitCacheEntry piece;
    const std::string dir = knit_test::knitCacheFixtureDir();
    std::string readErr;
    if (!dir.empty() && std::filesystem::exists(dir + "/piece_bear_009_12.knitcache") && readKnitCache(dir + "/piece_bear_009_12.knitcache", piece, &readErr)) {
        injectEverywhere("port validate", [&](std::string& e, std::vector<double>& r) {
            return portCall([&](std::vector<double>& out) {
                const ValidationReport v = port->validate(piece.knit.curves, piece.pp, ValidateOptions{});
                out.push_back(static_cast<int>(v.verdict));
                for (const CurvatureResult& c : v.curvature) out.push_back(c.rhoLo);
                out.push_back(v.clearance.minInterYarn);
                out.push_back(double(v.clearance.issues.size()));
                out.push_back(v.localSignature.pass);
            }, e, r);
        }, true);
    } else {
        std::printf("    port validate: skipped (KNIT_TEST_KNITCACHE without piece_bear_009_12)\n");
    }
    // The relaxation: two parallel strands closer than the hard contact (results not compared: the solver is
    // not bit-reproducible).
    injectEverywhere("port relax", [&](std::string& e, std::vector<double>& r) {
        return portCall([&](std::vector<double>& out) {
            CurveSet cs;
            cs.yarns = {strand(0, 0), strand(0.7, 1)};
            cs.wales = 1;
            cs.courses = 2;
            Controls P{cs.yarns[0].curve.ctrl, cs.yarns[1].curve.ctrl};
            const Controls rest = P;
            ControlMask freeCtrl(2, std::vector<char>(P[0].size(), 0));
            for (size_t c = 2; c + 2 < P[0].size(); ++c) freeCtrl[0][c] = freeCtrl[1][c] = 1;
            RelaxMoveParams mp;
            mp.thr = 1.0;
            mp.kappaMax = 2.0;
            mp.window = 3.0;
            mp.maxDisplacement = 0.5;
            port->relax({P, rest, freeCtrl, cs, nullptr, mp, nullptr});
            for (const Vec3& c : P[0]) addVec(out, c);
        }, e, r);
    }, false);
    // The tiled route through the port, settle off: the mesh surface, the lattice, B1 / B2, gate (d), gate (a).
    const TriMesh mesh = plane(40, 30);
    injectEverywhere("knitMesh through the port", [&](std::string& e, std::vector<double>& r) {
        return portCall([&](std::vector<double>& out) {
            MeshKnitParams p;
            p.tile.pp = PatternParams{};
            p.tile.g.w = p.tile.pp.waleSpacing;
            p.tile.g.h = p.tile.pp.courseSpacing;
            p.settle = false;
            p.surfaceGate = true;
            PartKnit k;
            std::string knitErr;
            if (!knitMesh(mesh, p, *port, k, &knitErr)) throw std::runtime_error(knitErr);
            out.push_back(static_cast<int>(k.verdict));
            out.push_back(k.boundary.minDistance);
            for (const YarnPath& y : k.t.curves.yarns)
                for (const Vec3& c : y.curve.ctrl) addVec(out, c);
        }, e, r);
    }, true);
}
