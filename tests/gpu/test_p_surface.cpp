// GPU port tests: the surface model (buildSurface) and gate (d) (checkSurface) against the CPU
// on cached teddy curves; deviations and times are printed, the checks are structural.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <vector>

#include "gpu_common.h"
#include "scene/knit/knit_cache.h"
#include "knitting/tube/yarn_surface.h"
#include "gpu/surfaces.h"

using namespace knit;

namespace {

bool loadCached(const char* stem, KnitCacheEntry& e) {
    std::string err;
    const std::string dir = knit_test::knitCacheFixtureDir();
    if (dir.empty()) return false;
    const std::string p = dir + "/" + stem + ".knitcache";
    return std::filesystem::exists(p) && readKnitCache(p, e, &err);
}
double msSince(std::chrono::steady_clock::time_point t) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t).count();
}

} // namespace

TEST(surface_model_gpu_vs_cpu) {
    REQUIRE_GPU();
    for (const char* stem : {"piece_bear_009_12", "piece_bear_002_00"}) {
        KnitCacheEntry e;
        if (!loadCached(stem, e)) SKIP("teddy knit cache not found");
        const CurveSet& cs = e.knit.curves;
        SurfaceParams sp;
        sp.compression = false;
        ExecPolicy ex;
        ex.parallel = true;
        auto t0 = std::chrono::steady_clock::now();
        const SurfaceModel a = buildSurface(cs, e.pp, sp, ex);
        const double msCpu = msSince(t0);
        SurfaceModel b;
        gpu::SurfaceGpuStats st;
        std::string err;
        t0 = std::chrono::steady_clock::now();
        CHECK(gpu::buildSurface({cs, e.pp, sp}, b, &st, &err));
        const double msGpu = msSince(t0);
        CHECK(a.yarns.size() == b.yarns.size());
        double dc = 0, dT = 0, dn = 0, du = 0;
        size_t rings = 0;
        int nonFinite = 0;
        for (size_t y = 0; y < a.yarns.size(); ++y) {
            CHECK(a.yarns[y].rings.size() == b.yarns[y].rings.size());
            for (size_t k = 0; k < a.yarns[y].rings.size(); ++k) {
                const Ring &A = a.yarns[y].rings[k], &B = b.yarns[y].rings[k];
                nonFinite += !std::isfinite(B.c.x + B.n1.x + B.n2.x + B.T.x);
                dc = std::max(dc, norm(A.c - B.c));
                dT = std::max(dT, norm(A.T - B.T));
                dn = std::max(dn, norm(A.n1 - B.n1));
                du = std::max(du, std::fabs(A.u - B.u));
                ++rings;
            }
        }
        std::printf("  %s: %zu rings; max |dc| %.2e mm, |dT| %.2e, |dn1| %.2e, |du| %.2e; CPU %.0f ms, GPU %.0f ms (%s)\n", stem, rings, dc, dT,
                    dn, du, msCpu, msGpu, st.line().c_str());
        CHECK(nonFinite == 0);
        CHECK(dc < 1e-2);
    }
}

// Rings start at r times the yarn's thickness on both sides, and the fit never grows one: the model's
// extent is known before any ring is computed, and no ring has to come back to the host for it.
TEST(surface_thickness_gpu_vs_cpu) {
    REQUIRE_GPU();
    KnitCacheEntry e;
    if (!loadCached("piece_bear_009_12", e)) SKIP("teddy knit cache not found");
    CurveSet cs = e.knit.curves;
    for (YarnPath& y : cs.yarns) {
        y.thickness.resize(y.curve.ctrl.size());
        for (size_t c = 0; c < y.thickness.size(); ++c) y.thickness[c] = 1.0f + 0.3f * static_cast<float>(std::sin(0.7 * double(c)));
    }
    PatternParams ppR = e.pp;
    ppR.yarnDiameter *= 1.25;
    for (const bool compression : {false, true}) {
        SurfaceParams sp;
        sp.compression = compression;
        ExecPolicy ex;
        ex.parallel = true;
        const SurfaceModel a = buildSurface(cs, ppR, sp, ex);
        SurfaceModel b;
        std::string err;
        CHECK(gpu::buildSurface({cs, ppR, sp}, b, nullptr, &err));
        double largest = 0, worst = 0, sum = 0;
        long long rings = 0;
        for (size_t y = 0; y < a.yarns.size(); ++y) {
            const double limit = ppR.radius() * cs.yarns[y].largestThickness();
            CHECK(a.yarns[y].rings.size() == b.yarns[y].rings.size());
            for (size_t k = 0; k < a.yarns[y].rings.size(); ++k) {
                const Ring& A = a.yarns[y].rings[k];
                const Ring& B = b.yarns[y].rings[k];
                largest = std::max({largest, B.a / limit, B.b / limit});
                worst = std::max(worst, std::fabs(A.a - B.a));
                sum += std::fabs(A.a - B.a);
                ++rings;
            }
        }
        std::printf("    compression %d: largest axis %.4f of the limit, a differs by mean %.5f, worst %.4f mm\n", int(compression), largest, sum / double(rings), worst);
        CHECK(largest <= 1 + 1e-5);
        CHECK(sum / double(rings) < (compression ? 0.01 : 1e-5));
    }
}

TEST(surface_compression_gpu_vs_cpu) {
    REQUIRE_GPU();
    for (const char* stem : {"piece_bear_009_12", "piece_bear_002_00"}) {
        KnitCacheEntry e;
        if (!loadCached(stem, e)) SKIP("teddy knit cache not found");
        const CurveSet& cs = e.knit.curves;
        // An inflated rendered tube (as FinalSurfaceSpec::inflate) so that neighbours overlap.
        PatternParams ppR = e.pp;
        ppR.yarnDiameter *= 1.25;
        SurfaceParams sp;
        sp.compression = true;
        ExecPolicy ex;
        ex.parallel = true;
        auto t0 = std::chrono::steady_clock::now();
        const SurfaceModel a = buildSurface(cs, ppR, sp, ex);
        const double msCpu = msSince(t0);
        SurfaceModel b;
        gpu::SurfaceGpuStats st;
        std::string err;
        t0 = std::chrono::steady_clock::now();
        CHECK(gpu::buildSurface({cs, ppR, sp}, b, &st, &err));
        const double msGpu = msSince(t0);
        CHECK(a.yarns.size() == b.yarns.size());
        const double r = ppR.radius();
        const double W = kPi * (1 + sp.maxExpansion) * r;
        // Brute force (ring centres, every yarn): the true half distance to the nearest centreline
        // outside the own half-turn window, capped at r.
        auto trueH = [&](size_t y, const Ring& R) {
            double best = 2 * r;
            for (size_t y2 = 0; y2 < a.yarns.size(); ++y2)
                for (const Ring& Q : a.yarns[y2].rings) {
                    if (y2 == y && std::fabs(Q.s - R.s) <= W) continue;
                    best = std::min(best, norm(Q.c - R.c));
                }
            return 0.5 * best;
        };
        double da = 0;
        std::vector<double> dAxis;
        size_t rings = 0, flatA = 0, flatB = 0, flatOne = 0;
        int nonFinite = 0;
        std::vector<double> dbs;
        std::vector<std::pair<size_t, size_t>> differ;  // |db| > 0.01 mm
        for (size_t y = 0; y < a.yarns.size() && y < b.yarns.size(); ++y)
            for (size_t k = 0; k < a.yarns[y].rings.size() && k < b.yarns[y].rings.size(); ++k) {
                const Ring &A = a.yarns[y].rings[k], &B = b.yarns[y].rings[k];
                nonFinite += !std::isfinite(B.a + B.b + B.phi);
                const bool fa = A.b < r - 1e-5, fb = B.b < r - 1e-5;  // float r differs from r by ~1e-8
                flatA += fa;
                flatB += fb;
                flatOne += fa != fb;
                dbs.push_back(std::fabs(A.b - B.b));
                if (dbs.back() > 0.01) differ.push_back({y, k});
                da = std::max(da, std::fabs(A.a - B.a));
                // The a-axis direction (phi mod pi) where both are flattened by more than 1 %.
                if (A.b < 0.99 * r && B.b < 0.99 * r && dbs.back() <= 0.01) dAxis.push_back(std::fabs(std::sin(A.phi - B.phi)));
                ++rings;
            }
        // The rings where the backends disagree, referred to the brute force (at most 48, evenly).
        int gpuCloser = 0, refereed = 0;
        const size_t stride = std::max<size_t>(1, differ.size() / 48);
        for (size_t j = 0; j < differ.size(); j += stride) {
            const auto [y, k] = differ[j];
            const Ring &A = a.yarns[y].rings[k], &B = b.yarns[y].rings[k];
            const double h = trueH(y, A);
            ++refereed;
            gpuCloser += std::fabs(B.b - h) < std::fabs(A.b - h);
            if (refereed <= 3)
                std::printf("    differ y %zu k %zu: brute force h %.4f | CPU b %.4f | GPU b %.4f\n", y, k, h, A.b, B.b);
        }
        std::sort(dbs.begin(), dbs.end());
        std::sort(dAxis.begin(), dAxis.end());
        dAxis.push_back(0);  // never empty
        std::printf("  %s: %zu rings; flattened CPU %zu, GPU %zu, on one side only %zu; |db| p50 %.2e, p99.9 %.2e, max %.2e mm; "
                    "|da| max %.2e mm; |sin dphi| p99.9 %.2e, max %.2e (both > 1 %%, agreeing rings); rings with |db| > 0.01 mm: %zu, refereed %d, "
                    "GPU closer to the brute force %d; clamped CPU %d / GPU %d; CPU %.0f ms, GPU %.0f ms (%s)\n",
                    stem, rings, flatA, flatB, flatOne, dbs[dbs.size() / 2], dbs[dbs.size() * 999 / 1000], dbs.back(), da, dAxis[dAxis.size() * 999 / 1000], dAxis[dAxis.size() - 2 < dAxis.size() ? dAxis.size() - 2 : 0],
                    differ.size(), refereed, gpuCloser, a.compression.clamped, b.compression.clamped, msCpu, msGpu, st.line().c_str());
        CHECK(nonFinite == 0);
        CHECK(flatA > 0);
        CHECK(dbs[dbs.size() * 999 / 1000] < 1e-4);
        // A former CPU defect (fixed with the Multi fit): yarn_surface.cpp accepted a
        // Newton convergence back onto the query point when its residual exceeded 1e-12, and clamped
        // the ring. Both now test the converged parameter against the own window; any remaining
        // disagreement must be in the GPU's favour.
        CHECK(gpuCloser == refereed);
    }
}

TEST(gate_d_gpu_vs_cpu) {
    REQUIRE_GPU();
    for (const char* stem : {"piece_bear_009_12", "piece_bear_002_00"}) {
        KnitCacheEntry e;
        if (!loadCached(stem, e)) SKIP("teddy knit cache not found");
        const CurveSet& cs = e.knit.curves;
        for (double spacing : {0.025, 0.0125}) {
            for (double Dscale : {1.0, 1.3}) {  // 1.3 x D: penetrations
                PatternParams pp = e.pp;
                pp.yarnDiameter *= Dscale;
                SurfaceParams sp;
                sp.compression = false;
                sp.ringSpacing = spacing;
                ExecPolicy ex;
                ex.parallel = true;
                auto t0 = std::chrono::steady_clock::now();
                const SurfaceModel sm = buildSurface(cs, pp, sp, ex);
                const SurfaceCheck a = checkSurface(cs, sm, {}, ex);
                const double msCpu = msSince(t0);
                SurfaceCheck b;
                gpu::SurfaceGpuStats st;
                std::string err;
                t0 = std::chrono::steady_clock::now();
                CHECK(gpu::checkSurface({cs, pp, sp}, {}, b, &st, &err));
                const double msGpu = msSince(t0);
                std::printf("  %s %.4f D x%.1f: folds %s/%s (rings %d/%d, min ratio %.4f/%.4f), collisions %s/%s (%ld/%ld of %ld/%ld points, "
                            "max depth %.4f/%.4f mm); CPU %.0f ms, GPU %.0f ms (%s)\n",
                            stem, spacing, Dscale, toString(a.fold), toString(b.fold), a.foldRings, b.foldRings, a.minFoldRatio, b.minFoldRatio,
                            toString(a.collisions), toString(b.collisions), a.penetratingPoints, b.penetratingPoints, a.pointsTested,
                            b.pointsTested, a.maxDepth, b.maxDepth, msCpu, msGpu, st.line().c_str());
                CHECK(a.pointsTested == b.pointsTested);
                CHECK(std::isfinite(b.minFoldRatio));
            }
        }
    }
}

TEST(tube_mesh_gpu_vs_cpu) {
    REQUIRE_GPU();
    for (const char* stem : {"piece_bear_009_12", "piece_bear_002_00"}) {
        KnitCacheEntry e;
        if (!loadCached(stem, e)) SKIP("teddy knit cache not found");
        SurfaceParams sp;
        sp.compression = false;
        sp.ringSpacing = 0.025;  // the export's model (kExportModelRing with compression off is 0.2; both tested by the pipeline)
        TubeParams tp;
        tp.ringSpacing = 0.4;
        tp.sides = 8;
        tp.uvRepeat = 2.5;
        tp.uvWrap = 64;
        auto t0 = std::chrono::steady_clock::now();
        const SurfaceModel sm = buildSurface(e.knit.curves, e.pp, sp);
        const TubeMesh a = buildTubeMesh(e.knit.curves, sm, tp);
        const double msCpu = msSince(t0);
        TubeMesh b;
        std::string err;
        gpu::SurfaceGpuStats st;
        t0 = std::chrono::steady_clock::now();
        CHECK(gpu::buildTubeMesh({e.knit.curves, e.pp, sp}, tp, b, &st, &err));
        const double msGpu = msSince(t0);
        CHECK(a.pos.size() == b.pos.size() && a.tri == b.tri && a.yarnFirstVertex == b.yarnFirstVertex);
        double dp = 0, dn = 0, duv = 0;
        for (size_t v = 0; v < a.pos.size(); ++v) {
            dp = std::max(dp, norm(a.pos[v] - b.pos[v]));
            dn = std::max(dn, norm(a.nrm[v] - b.nrm[v]));
            duv = std::max(duv, double(std::fabs(a.uv[2 * v] - b.uv[2 * v])));
        }
        std::printf("  %s: %zu vertices, %zu triangles identical; max |dp| %.2e mm, |dn| %.2e, |du| %.2e; CPU (serial, as the export) %.0f ms, "
                    "GPU %.0f ms\n", stem, a.pos.size(), a.tri.size() / 3, dp, dn, duv, msCpu, msGpu);
        CHECK(dp < 1e-2);
    }
}

TEST(plies_gpu_vs_cpu) {
    REQUIRE_GPU();
    for (const char* stem : {"piece_bear_009_12", "piece_bear_002_00"}) {
        KnitCacheEntry e;
        if (!loadCached(stem, e)) SKIP("teddy knit cache not found");
        SurfaceParams sp;
        sp.compression = false;
        ExecPolicy ex;
        ex.parallel = true;
        const SurfaceModel sm = buildSurface(e.knit.curves, e.pp, sp, ex);
        PlyGeometry g;
        g.count = 3;
        auto t0 = std::chrono::steady_clock::now();
        const PlySet a = buildPlies(e.knit.curves, sm, g, e.pp.radius(), ex);
        const double msCpu = msSince(t0);
        PlySet b;
        gpu::PliesGpuStats st;
        std::string err;
        t0 = std::chrono::steady_clock::now();
        CHECK(gpu::buildPlies({e.knit.curves, sm, g, e.pp.radius()}, b, &st, &err));
        const double msGpu = msSince(t0);
        CHECK(a.curves.yarns.size() == b.curves.yarns.size() && a.yarnOf == b.yarnOf && a.plyOf == b.plyOf);
        double dctrl = 0, dlen = 0;
        for (size_t p = 0; p < a.curves.yarns.size(); ++p) {
            const auto &P = a.curves.yarns[p].curve.ctrl, &Q = b.curves.yarns[p].curve.ctrl;
            CHECK(P.size() == Q.size());
            for (size_t k = 0; k < P.size(); ++k) dctrl = std::max(dctrl, norm(P[k] - Q[k]));
            dlen = std::max(dlen, std::fabs(a.curves.yarns[p].length() - b.curves.yarns[p].length()));
        }
        std::printf("  %s: %zu plies; max |dctrl| %.2e mm, |dlength| %.2e mm; containment %.3e / %.3e mm, min centre distance %.4f / %.4f; "
                    "CPU %.0f ms, GPU %.0f ms (%s)\n", stem, a.curves.yarns.size(), dctrl, dlen, a.containment, b.containment,
                    a.minCentreDistance, b.minCentreDistance, msCpu, msGpu, st.line().c_str());
        CHECK(dctrl < 1e-2);
    }
}
