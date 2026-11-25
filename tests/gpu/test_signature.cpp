// GPU port tests (WP-S): the local signature on the device (signature.cu) against the CPU
// checkLocalSignature, on synthetic spiral / shaping / tiled pattern sets, broken sets, and the
// teddy knit caches (read only; SKIP if absent).
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <functional>
#include <map>
#include <string>

#include "gpu_common.h"
#include "../cpu_numerics.h"
#include "knitting/checks/checks.h"
#include "knitting/knit/fixtures.h"
#include "model/env.h"
#include "scene/knit/knit_cache.h"
#include "knitting/layout/stitch_layout.h"
#include "knitting/stitch/stitch_program.h"
#include "knitting/knit/mesh_knit.h"
#include "gpu/validation/checks.h"

using namespace knit;

namespace {

struct SigCompare {
    int pairs = 0, statusAgree = 0, crossAgree = 0;
    int conf[3][3] = {};  // [cpu][gpu]: Pass, Fail, Unresolved
    double msCpu = 0, msGpu = 0;
    std::string firstDiff;
};

int idx(Status s) { return s == Status::Pass ? 0 : s == Status::Fail ? 1 : 2; }

SigCompare compareSig(const CurveSet& cs, const PatternParams& pp, const char* label, bool print = true) {
    SigCompare c;
    ExecPolicy ex;
    ex.parallel = true;
    auto t0 = std::chrono::steady_clock::now();
    // the device port computes the ownership-span attribution (window 0); the attribution
    // correction is shared CPU code (applyCrossingAttribution), applied after it in validate
    const LocalSignatureResult cpu = checkLocalSignature(cs, pp, ex, 0.0);
    c.msCpu = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    t0 = std::chrono::steady_clock::now();
    const LocalSignatureResult g = gpu::checkSignatureOnDevice(gpu::uploadCurves(cs));
    c.msGpu = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    CHECK(cpu.pairs.size() == g.pairs.size());
    c.pairs = static_cast<int>(cpu.pairs.size());
    for (size_t k = 0; k < cpu.pairs.size(); ++k) {
        const LocalPairResult &a = cpu.pairs[k], &b = g.pairs[k];
        CHECK(a.child == b.child && a.parent == b.parent && a.terminal == b.terminal);
        ++c.conf[idx(a.status)][idx(b.status)];
        c.statusAgree += a.status == b.status;
        const bool cross = a.crossings == b.crossings && a.unexpected == b.unexpected && a.overlap == b.overlap;
        c.crossAgree += cross;
        if ((a.status != b.status || !cross) && c.firstDiff.empty()) {
            char buf[512];
            std::snprintf(buf, sizeof buf, "pair %zu (child %d%s): CPU %s %d/%d/%d '%s' vs GPU %s %d/%d/%d '%s'", k, a.child,
                          a.terminal ? " CLOSE" : "", toString(a.status), a.crossings, a.unexpected, a.overlap, a.detail.c_str(),
                          toString(b.status), b.crossings, b.unexpected, b.overlap, b.detail.c_str());
            c.firstDiff = buf;
        }
    }
    if (print) {
        std::printf("  %s: %d pairs, status agree %d, crossings agree %d; CPU P/F/U -> GPU [P %d %d %d | F %d %d %d | U %d %d %d]; "
                    "CPU %.1f ms (all threads), GPU %.1f ms; CPU %s / GPU %s\n",
                    label, c.pairs, c.statusAgree, c.crossAgree, c.conf[0][0], c.conf[0][1], c.conf[0][2], c.conf[1][0], c.conf[1][1],
                    c.conf[1][2], c.conf[2][0], c.conf[2][1], c.conf[2][2], c.msCpu, c.msGpu, cpu.detail.c_str(), g.detail.c_str());
        if (!c.firstDiff.empty()) std::printf("    first difference: %s\n", c.firstDiff.c_str());
    }
    return c;
}

struct Cyl {
    PatternParams pp;
    CurveSet cs;
};
Cyl cylinder(int n, int turns) {
    Cyl c;
    GaugeParams g;
    c.pp = findPreset("tightest")->params;
    c.pp.waleSpacing = g.w;
    c.pp.courseSpacing = g.h;
    const double R = n * g.w / (2 * kPi), H = (turns + 2) * g.h;
    ConeSurface cyl(R, R, H);
    const StitchProgram prog = buildProgram(std::vector<int>(turns, n));
    const PieceLayout lay = helicalLayout(cyl, prog, g, evenHelix(prog, g.h, H));
    c.cs = placeSpiral(cyl, lay, prog, c.pp);
    finalizeCurves(c.cs);
    return c;
}

TriMesh planeMesh(int nx, int ny) {
    TriMesh m;
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i) m.pos.push_back({double(i), double(j), 0});
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            const int a = j * (nx + 1) + i, b = a + 1, d = a + nx + 2, e = a + nx + 1;
            auto uv = [&](int v) { return Vec2{(v % (nx + 1)) / 100.0, (v / (nx + 1)) / 100.0}; };
            m.tri.push_back({a, b, d});
            m.uv.push_back({uv(a), uv(b), uv(d)});
            m.tri.push_back({a, d, e});
            m.uv.push_back({uv(a), uv(d), uv(e)});
        }
    m.computeNormals();
    return m;
}

} // namespace

TEST(sig_gpu_cylinder) {
    REQUIRE_GPU();
    const Cyl c = cylinder(40, 8);
    const SigCompare r = compareSig(c.cs, c.pp, "cylinder 40 x 8");
    CHECK(r.statusAgree == r.pairs && r.crossAgree == r.pairs);
    CHECK(r.conf[0][0] == 40 * 7);
}

TEST(sig_gpu_shaping) {
    REQUIRE_GPU();
    for (const char* kind : {"cone"}) {
        ShapingDemo d;
        CHECK(makeShapingDemo(kind, {}, d));
        const SigCompare r = compareSig(d.curves, d.pp, kind);
        CHECK(r.statusAgree == r.pairs && r.crossAgree == r.pairs);
    }
}

TEST(sig_gpu_broken) {
    // A stitch pulled through to the front of its parent: the CPU's FAIL pairs are FAIL on the GPU.
    REQUIRE_GPU();
    Cyl c = cylinder(30, 5);
    const StitchRecord& st = c.cs.topo.stitches[30 * 2 + 7];
    const int k0 = static_cast<int>(std::floor(st.uBegin)) + 1, k1 = static_cast<int>(std::floor(st.uEnd));
    for (int k = k0; k <= k1; ++k) c.cs.yarns[st.yarn].curve.ctrl[k] += st.N * 1.5;
    c.cs.yarns[st.yarn].finalize();
    const SigCompare r = compareSig(c.cs, c.pp, "cylinder, one stitch moved to the front");
    CHECK(r.conf[1][1] >= 1);
    CHECK(r.conf[1][0] == 0 && r.conf[2][0] == 0);  // nothing the CPU rejects passes on the GPU
    CHECK(r.statusAgree == r.pairs);
}

TEST(sig_gpu_tiled_patterns) {
    REQUIRE_GPU();
    const TriMesh m = planeMesh(70, 50);
    for (StitchPattern pat : {StitchPattern::Stockinette, StitchPattern::ReverseStockinette, StitchPattern::Rib1x1, StitchPattern::Garter}) {
        MeshKnitParams p;
        PatternParams pp;  // expanded preset 8 x 6.5, D 1
        pp.pattern = pat;
        p.tile.pp = pp;
        p.tile.g.w = pp.waleSpacing;
        p.tile.g.h = pp.courseSpacing;
        p.settle = false;
        PartKnit k;
        std::string err;
        CHECK(knitMesh(m, p, knit_test::cpuNumerics(), k, &err));
        char label[64];
        std::snprintf(label, sizeof label, "tiled plane, pattern %d", static_cast<int>(pat));
        const SigCompare r = compareSig(k.t.curves, k.t.pp, label);
        CHECK(r.statusAgree == r.pairs && r.crossAgree == r.pairs);
    }
}

TEST(sig_gpu_teddy_caches) {
    REQUIRE_GPU();
    const std::string dir = knit_test::knitCacheFixtureDir();
    if (dir.empty()) SKIP("KNIT_TEST_KNITCACHE (folder of teddy knit caches) not set");
    if (!std::filesystem::exists(dir)) SKIP("teddy knit caches not found at KNIT_TEST_KNITCACHE");
    int pairs = 0, agree = 0, cross = 0, conf[3][3] = {};
    double cpu = 0, gpuMs = 0;
    for (const auto& e : std::filesystem::directory_iterator(dir)) {
        if (e.path().extension() != ".knitcache") continue;
        if (const auto only = envVar("KNIT_SIG_PIECE"); only && e.path().stem().string() != *only) continue;
        KnitCacheEntry ce;
        std::string err;
        CHECK(readKnitCache(e.path().string(), ce, &err));
        compareSig(ce.knit.curves, ce.pp, e.path().stem().string().c_str(), false);  // warm-up (context, first-use)
        const SigCompare r = compareSig(ce.knit.curves, ce.pp, e.path().stem().string().c_str());
        pairs += r.pairs;
        agree += r.statusAgree;
        cross += r.crossAgree;
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) conf[a][b] += r.conf[a][b];
        cpu += r.msCpu;
        gpuMs += r.msGpu;
    }
    std::printf("  teddy total: %d pairs, status agree %d (%.4f), crossings agree %d; CPU P/F/U -> GPU [P %d %d %d | F %d %d %d | U %d %d %d]; "
                "CPU %.0f ms, GPU %.0f ms\n",
                pairs, agree, pairs ? double(agree) / pairs : 0.0, cross, conf[0][0], conf[0][1], conf[0][2], conf[1][0], conf[1][1], conf[1][2],
                conf[2][0], conf[2][1], conf[2][2], cpu, gpuMs);
    CHECK(pairs > 0 && agree >= pairs * 0.99);
}
