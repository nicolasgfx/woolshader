// GPU port tests (WP-V): curvature, sampling, clearance and finalize on the device against the
// CPU code. Deviations are recorded (printed), not required to vanish; the checks demand
// finite results, no lost work, equal structure where the algorithm fixes it, and agreeing
// verdicts on synthetic cases far from their thresholds.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

#include "gpu_common.h"
#include "knitting/checks/checks.h"
#include "scene/knit/knit_cache.h"
#include "model/parallel.h"
#include "gpu/validation/checks.h"
#include "gpu/device.h"
#include "gpu/surfaces.h"

using namespace knit;

namespace {

ExecPolicy par() {
    ExecPolicy e;
    e.parallel = true;
    return e;
}

double msSince(std::chrono::steady_clock::time_point t) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t).count();
}

// A helix (radius R, pitch p per turn) as one uniform B-spline, `turns` turns, 24 controls per turn.
CurveSet helix(double R, double pitch, int turns) {
    CurveSet cs;
    cs.kind = CurveKind::Surface;
    YarnPath y;
    y.id = 0;
    const int per = 24, n = turns * per + 3;
    for (int i = 0; i < n; ++i) {
        const double a = 2 * kPi * i / per;
        y.curve.ctrl.push_back({R * std::cos(a), R * std::sin(a), pitch * i / per});
    }
    y.finalize();
    cs.yarns.push_back(std::move(y));
    return cs;
}

bool loadCache(const std::string& stem, KnitCacheEntry& e) {
    const std::string dir = knit_test::knitCacheFixtureDir();
    if (dir.empty()) return false;
    const std::string p = dir + "/" + stem + ".knitcache";
    std::string err;
    return std::filesystem::exists(p) && readKnitCache(p, e, &err);
}

struct Q {
    std::vector<double> v;
    void add(double x) { v.push_back(x); }
    double q(double p) {
        if (v.empty()) return 0;
        std::sort(v.begin(), v.end());
        return v[std::min(v.size() - 1, static_cast<size_t>(p * (v.size() - 1) + 0.5))];
    }
};

// The four parts on one curve set: CPU and GPU, compared; returns nothing, prints a report.
void compareParts(const char* label, const CurveSet& cs, const PatternParams& pp, bool timing) {
    using namespace knit::gpu;
    const Allowances al = Allowances::forExtent(cs.extent());
    const ExecPolicy ex = par();
    const UploadedCurves d = uploadCurves(cs);
    // ---- curvature
    std::vector<CurvatureResult> cc;
    auto t0 = std::chrono::steady_clock::now();
    for (const YarnPath& y : cs.yarns) cc.push_back(checkCurvature(y, pp.radius(), al, ex));
    const double cpuCurv = msSince(t0);
    checkCurvatureOnDevice(d, pp.radius(), al);  // warm-up
    t0 = std::chrono::steady_clock::now();
    const std::vector<CurvatureResult> gc = checkCurvatureOnDevice(d, pp.radius(), al);
    const double gpuCurv = msSince(t0);
    CHECK(gc.size() == cc.size());
    Q kdev;
    int statusAgree = 0, statusTotal = 0;
    double rhoLoC = 1e300, rhoLoG = 1e300;
    for (size_t y = 0; y < cc.size(); ++y) {
        CHECK(gc[y].spanKappaHi.size() == cc[y].spanKappaHi.size());
        for (size_t k = 0; k < cc[y].spanKappaHi.size(); ++k) {
            const double a = cc[y].spanKappaHi[k], b = gc[y].spanKappaHi[k];
            CHECK(std::isfinite(b));
            if (a < 1e299 && b < 1e299 && a > 0) kdev.add(std::fabs(b - a) / a);
        }
        statusAgree += (gc[y].regular == cc[y].regular) + (gc[y].noFold == cc[y].noFold);
        statusTotal += 2;
        rhoLoC = std::min(rhoLoC, cc[y].rhoLo);
        rhoLoG = std::min(rhoLoG, gc[y].rhoLo);
    }
    std::printf("  [%s] curvature: span kappa_hi rel. dev p50 %.2e p99 %.2e max %.2e; statuses agree %d/%d; rho_lo CPU %.6f GPU %.6f\n",
                label, kdev.q(0.5), kdev.q(0.99), kdev.q(1.0), statusAgree, statusTotal, rhoLoC, rhoLoG);
    // ---- sampling
    std::vector<Polyline> cp;
    t0 = std::chrono::steady_clock::now();
    for (const YarnPath& y : cs.yarns) cp.push_back(samplePath(y, pp.sampleTolerance, 512, ex));
    const double cpuSamp = msSince(t0);
    std::vector<Polyline> gp;
    t0 = std::chrono::steady_clock::now();
    const gpu::DevicePolylines dpl = samplePolylinesOnDevice(d, pp.sampleTolerance, 512, gp);
    const double gpuSamp = msSince(t0);
    CHECK(gp.size() == cp.size());
    long sameCount = 0, yarnsSameCount = 0;
    double maxVert = 0, maxS = 0, maxCum = 0;
    for (size_t y = 0; y < cp.size(); ++y) {
        if (gp[y].p.size() != cp[y].p.size()) continue;
        ++yarnsSameCount;
        sameCount += static_cast<long>(gp[y].p.size());
        for (size_t v = 0; v < cp[y].p.size(); ++v) {
            if (gp[y].u[v] != cp[y].u[v]) continue;
            maxVert = std::max(maxVert, norm(gp[y].p[v] - cp[y].p[v]));
            maxS = std::max(maxS, std::fabs(gp[y].s[v] - cp[y].s[v]));
            CHECK(std::isfinite(gp[y].s[v]));
        }
        for (size_t v = 0; v < cp[y].cumHi.size(); ++v) maxCum = std::max(maxCum, std::fabs(gp[y].cumHi[v] - cp[y].cumHi[v]));
    }
    std::printf("  [%s] sampling: %ld/%zu yarns with equal vertex counts (%ld vertices); max |dp| %.2e mm, max |ds| %.2e mm, max |dcumHi| %.2e mm\n",
                label, yarnsSameCount, cp.size(), sameCount, maxVert, maxS, maxCum);
    // ---- clearance (on the GPU polylines; CPU on its own)
    t0 = std::chrono::steady_clock::now();
    const ClearanceResult cr = checkClearance({cs, cp, cc}, pp, al, ex);
    const double cpuClear = msSince(t0);
    t0 = std::chrono::steady_clock::now();
    const ClearanceResult gr = checkClearanceOnDevice({d, dpl, gp, gc, pp, al});
    const double gpuClear = msSince(t0);
    auto regs = [](const ClearanceResult& r, Status s) {
        int n = 0;
        for (const ClearanceRegion& g : r.regions) n += g.status == s;
        return n;
    };
    std::printf("  [%s] clearance: status CPU %s GPU %s; FAIL regions %d / %d, UNRESOLVED %d / %d; issues %zu / %zu; min inter %.5g / %.5g, "
                "min self %.5g / %.5g; segment pairs %ld / %ld, Schur %ld / %ld; delta %.2e / %.2e mm\n",
                label, toString(cr.status), toString(gr.status), regs(cr, Status::Fail), regs(gr, Status::Fail), regs(cr, Status::Unresolved),
                regs(gr, Status::Unresolved), cr.issues.size(), gr.issues.size(), cr.minInterYarn, gr.minInterYarn, cr.minSelfNonlocal,
                gr.minSelfNonlocal, cr.segmentPairs, gr.segmentPairs, cr.schurCertified, gr.schurCertified, cr.deltaD, gr.deltaD);
    // ---- finalize
    CurveSet a = cs, b = cs;
    t0 = std::chrono::steady_clock::now();
    for (YarnPath& y : a.yarns) y.finalize(1e-9, ex);
    const double cpuFin = msSince(t0);
    t0 = std::chrono::steady_clock::now();
    computeArcLengthTables(b);
    const double gpuFin = msSince(t0);
    double maxRelL = 0, maxDu = 0, maxDs = 0;
    size_t leavesA = 0, leavesB = 0;
    for (size_t y = 0; y < a.yarns.size(); ++y) {
        const YarnPath &ya = a.yarns[y], &yb = b.yarns[y];
        CHECK(std::isfinite(yb.length()));
        maxRelL = std::max(maxRelL, std::fabs(yb.length() - ya.length()) / std::max(1e-12, ya.length()));
        leavesA += ya.tables().tu.size();
        leavesB += yb.tables().tu.size();
        for (int k = 0; k <= 1000; ++k) {
            const double s = ya.length() * k / 1000;
            maxDu = std::max(maxDu, std::fabs(yb.paramAtArcLength(s) - ya.paramAtArcLength(s)));
            const double u = ya.uEnd() * k / 1000;
            maxDs = std::max(maxDs, std::fabs(yb.arcLengthAt(u) - ya.arcLengthAt(u)));
        }
    }
    std::printf("  [%s] finalize: length rel. dev max %.2e; paramAtArcLength max |du| %.2e; arcLengthAt max |ds| %.2e mm; table nodes %zu / %zu\n",
                label, maxRelL, maxDu, maxDs, leavesA, leavesB);
    if (timing)
        std::printf("  [%s] timing CPU (%d threads) vs GPU [ms]: curvature %.1f / %.1f, sampling %.1f / %.1f, clearance %.1f / %.1f, "
                    "finalize %.1f / %.1f\n",
                    label, threadBudget(), cpuCurv, gpuCurv, cpuSamp, gpuSamp, cpuClear, gpuClear, cpuFin, gpuFin);
}

} // namespace

// Synthetic helices: far from every threshold, the verdicts must agree.
// The failure model (src/gpu/ARCHITECTURE.md): a device failure is an exception below the public
// headers and false plus the reason at them, and it leaves the device usable.
TEST(gpu_device_failure_is_false_and_err) {
    REQUIRE_GPU();
    const CurveSet cs = helix(5, 3, 2);
    PatternParams pp;
    const SurfaceParams sp;
    SurfaceModel model;
    std::string err;
    {
        gpu::AllocationCap cap(gpu::memStats().current + 1);
        CHECK(!gpu::buildSurface({cs, pp, sp}, model, nullptr, &err));
        CHECK(err.find("device memory cap") != std::string::npos);
        bool thrown = false;
        try {
            gpu::uploadCurves(cs);
        } catch (const gpu::DeviceBudgetExceeded&) {
            thrown = true;
        }
        CHECK(thrown);
    }
    CHECK(gpu::buildSurface({cs, pp, sp}, model, nullptr, &err));
    CHECK(!model.yarns.empty());
}

TEST(gpu_validate_helix_parts) {
    REQUIRE_GPU();
    using namespace knit::gpu;
    PatternParams pp;
    pp.yarnDiameter = 1.0;
    for (double pitch : {3.0, 0.8}) {  // turns 3 mm apart (PASS) or 0.8 mm (FAIL: < D + c)
        const CurveSet cs = helix(4.0, pitch, 12);
        char label[64];
        std::snprintf(label, sizeof label, "helix pitch %.1f", pitch);
        compareParts(label, cs, pp, false);
        const Allowances al = Allowances::forExtent(cs.extent());
        const UploadedCurves d = uploadCurves(cs);
        const std::vector<CurvatureResult> gc = checkCurvatureOnDevice(d, pp.radius(), al);
        CHECK(gc[0].regular == Status::Pass && gc[0].noFold == Status::Pass);
        CHECK(std::fabs(gc[0].rhoEst - checkCurvature(cs.yarns[0], pp.radius(), al).rhoEst) < 1e-3 * gc[0].rhoEst);
        std::vector<Polyline> gp;
        const gpu::DevicePolylines dpl = samplePolylinesOnDevice(d, pp.sampleTolerance, 512, gp);
        CHECK(gp[0].p.size() == samplePath(cs.yarns[0], pp.sampleTolerance).p.size());
        const ClearanceResult gr = checkClearanceOnDevice({d, dpl, gp, gc, pp, al});
        const ClearanceResult cr = checkClearance({cs, std::vector<Polyline>{samplePath(cs.yarns[0], pp.sampleTolerance)},
                                                   std::vector<CurvatureResult>{checkCurvature(cs.yarns[0], pp.radius(), al)}}, pp, al);
        CHECK(gr.status == cr.status);
        CHECK(gr.status == (pitch > 2 ? Status::Pass : Status::Fail));
        CHECK(gr.regions.size() == cr.regions.size());
        CHECK(std::fabs(gr.minSelfNonlocal - cr.minSelfNonlocal) < 1e-4);
    }
}

// Two yarns (inter-yarn pairs): parallel helices at a known distance.
TEST(gpu_validate_two_yarns_clearance) {
    REQUIRE_GPU();
    using namespace knit::gpu;
    PatternParams pp;
    pp.yarnDiameter = 1.0;
    for (double gap : {1.3, 0.9}) {
        CurveSet cs = helix(4.0, 6.0, 4);
        CurveSet other = helix(4.0 + gap, 6.0, 4);
        other.yarns[0].id = 1;
        cs.yarns.push_back(other.yarns[0]);
        const Allowances al = Allowances::forExtent(cs.extent());
        const UploadedCurves d = uploadCurves(cs);
        const std::vector<CurvatureResult> gc = checkCurvatureOnDevice(d, pp.radius(), al);
        std::vector<Polyline> gp;
        const gpu::DevicePolylines dpl = samplePolylinesOnDevice(d, pp.sampleTolerance, 512, gp);
        const ClearanceResult gr = checkClearanceOnDevice({d, dpl, gp, gc, pp, al});
        std::vector<CurvatureResult> cc;
        std::vector<Polyline> cp;
        for (const YarnPath& y : cs.yarns) cc.push_back(checkCurvature(y, pp.radius(), al)), cp.push_back(samplePath(y, pp.sampleTolerance));
        const ClearanceResult cr = checkClearance({cs, cp, cc}, pp, al);
        std::printf("  gap %.1f: status CPU %s GPU %s, regions %zu / %zu, min inter %.6f / %.6f\n", gap, toString(cr.status), toString(gr.status),
                    cr.regions.size(), gr.regions.size(), cr.minInterYarn, gr.minInterYarn);
        CHECK(gr.status == cr.status);
        CHECK(gr.regions.size() == cr.regions.size());
        CHECK(std::fabs(gr.minInterYarn - cr.minInterYarn) < 1e-4);
    }
}

// The teddy pieces (CPU knit cache, final curves): every part, deviations and timings.
TEST(gpu_validate_teddy_parts) {
    REQUIRE_GPU();
    KnitCacheEntry e;
    for (const char* stem : {"piece_bear_009_12", "piece_bear_002_00"}) {
        if (!loadCache(stem, e)) SKIP("teddy knit cache (knitcache_tex) not found");
        compareParts(stem, e.knit.curves, e.pp, true);
        // Perturbed: the same centreline at a larger diameter (D + 0.3 mm), so clearance FAILs.
        PatternParams big = e.pp;
        big.yarnDiameter += 0.3;
        compareParts((std::string(stem) + " D+0.3").c_str(), e.knit.curves, big, std::string(stem) == "piece_bear_002_00");
    }
}

// Every cached teddy piece: verdict agreement of curvature and clearance (record only).
TEST(gpu_validate_teddy_all_pieces) {
    REQUIRE_GPU();
    using namespace knit::gpu;
    const std::string cacheDir = knit_test::knitCacheFixtureDir();
    if (cacheDir.empty()) SKIP("KNIT_TEST_KNITCACHE (folder of teddy knit caches) not set");
    if (!std::filesystem::exists(cacheDir)) SKIP("teddy knit caches not found at KNIT_TEST_KNITCACHE");
    int pieces = 0, agree = 0;
    for (const auto& f : std::filesystem::directory_iterator(cacheDir)) {
        if (f.path().extension() != ".knitcache") continue;
        KnitCacheEntry e;
        std::string err;
        if (!readKnitCache(f.path().string(), e, &err)) continue;
        const CurveSet& cs = e.knit.curves;
        const Allowances al = Allowances::forExtent(cs.extent());
        const ExecPolicy ex = par();
        const UploadedCurves d = uploadCurves(cs);
        const std::vector<CurvatureResult> gc = checkCurvatureOnDevice(d, e.pp.radius(), al);
        std::vector<CurvatureResult> cc;
        std::vector<Polyline> gp, cp;
        const gpu::DevicePolylines dpl = samplePolylinesOnDevice(d, e.pp.sampleTolerance, 512, gp);
        const ClearanceResult gr = checkClearanceOnDevice({d, dpl, gp, gc, e.pp, al});
        for (const YarnPath& y : cs.yarns) cc.push_back(checkCurvature(y, e.pp.radius(), al, ex)), cp.push_back(samplePath(y, e.pp.sampleTolerance, 512, ex));
        const ClearanceResult cr = checkClearance({cs, cp, cc}, e.pp, al, ex);
        const bool same = gr.status == cr.status && gc[0].noFold == cc[0].noFold;
        ++pieces;
        agree += same;
        std::printf("  %-22s noFold %s/%s rho_lo %.4f/%.4f clearance %s/%s regions %zu/%zu min self %.5g/%.5g\n",
                    f.path().stem().string().c_str(), toString(cc[0].noFold), toString(gc[0].noFold), cc[0].rhoLo, gc[0].rhoLo,
                    toString(cr.status), toString(gr.status), cr.regions.size(), gr.regions.size(), cr.minSelfNonlocal, gr.minSelfNonlocal);
    }
    std::printf("  %d/%d pieces with equal curvature and clearance verdicts (CPU / GPU)\n", agree, pieces);
    CHECK(pieces > 0);
}
