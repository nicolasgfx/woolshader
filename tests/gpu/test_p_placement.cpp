// GPU port tests: the placement's surface queries on the device (makePieceQuery) against the
// CPU PieceSurface, and knitPiece with GPU queries against the CPU knitPiece on teddy pieces.
// Deviations are measured and printed (float32 geometry); the checks are structural.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <map>
#include <random>

#include "gpu_common.h"
#include "../cpu_numerics.h"
#include "../strand_fixtures.h"
#include "model/env.h"
#include "model/hash.h"
#include "knitting/knit/piece_knit.h"
#include "knitting/knit/pieces.h"
#include "knitting/stitch/stitch_program.h"
#include "gpu/knitting.h"
#include "gpu/surface/placement.h"

using namespace knit;

namespace {

// The teddy pieces live in the main checkout (generated, untracked): KNIT_TEDDY_MANIFEST, or
// the worktree's / main checkout's scenes/teddy-bears/pieces/manifest.txt.
std::string teddyManifest() {
    if (const auto e = envVar("KNIT_TEDDY_MANIFEST")) return *e;
    const std::string p = KNIT_SOURCE_DIR "/scenes/teddy-bears/pieces/manifest.txt";
    return std::filesystem::exists(p) ? p : "";
}

struct Teddy {
    Manifest man;
    double sigma = 0;
    std::vector<TriMesh> scene;
};
bool loadTeddy(Teddy& t) {
    const std::string path = teddyManifest();
    std::string err;
    if (path.empty() || !readManifest(path, t.man, &err)) return false;
    Gauge g;
    if (!computeGauge(t.man, {"piece_bear_002_00.obj", 65, 5.5}, g, &err)) return false;
    t.sigma = g.sigma;
    t.scene.resize(t.man.pieces.size());
    for (size_t k = 0; k < t.man.pieces.size(); ++k)
        if (!readObj(piecePath(t.man, t.man.pieces[k].file), t.scene[k], &err)) return false;
    return true;
}

double ms(std::chrono::steady_clock::time_point t) { return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t).count(); }

} // namespace

// Random chart coordinates (also in the pole caps) and random points near the surface: the
// device at / project against PieceSurface::at / project.
TEST(placement_queries_match_cpu) {
    REQUIRE_GPU();
    Teddy t;
    knit_test::requireScene(teddyManifest());
    CHECK(loadTeddy(t));
    for (int piece : {17, 0}) {
        LoadedPiece lp;
        std::string err;
        CHECK(loadPiece(t.man, piece, t.sigma, lp, &err));
        Chart ch(lp.mesh, lp.entry.rings);
        PieceSurface cpu(ch);
        gpu::PlacementStats st;
        std::unique_ptr<BaseSurface> q = gpu::makePieceQuery(ch, &st, &err);
        CHECK(q != nullptr);
        std::mt19937_64 rng(7);
        std::uniform_real_distribution<double> U(0, 1);
        std::vector<Vec2> ab;
        for (int i = 0; i < 20000; ++i) ab.push_back({U(rng), U(rng)});
        std::vector<SurfaceSample> g, c;
        q->atBatch(ab, g);
        cpu.atBatch(ab, c);
        double dp = 0, dn = 0;
        int nonFinite = 0;
        for (size_t i = 0; i < ab.size(); ++i) {
            nonFinite += !std::isfinite(g[i].p.x + g[i].p.y + g[i].p.z);
            dp = std::max(dp, norm(g[i].p - c[i].p));
            dn = std::max(dn, norm(g[i].n - c[i].n));
        }
        std::vector<Vec3> pts;
        for (size_t i = 0; i < ab.size(); ++i) pts.push_back(c[i].p + c[i].n * (4 * U(rng) - 2));
        std::vector<SurfaceSample> gp, cp;
        q->projectBatch(pts, gp);
        cpu.projectBatch(pts, cp);
        double pp = 0, pn = 0;
        int far = 0;
        for (size_t i = 0; i < pts.size(); ++i) {
            const double d = norm(gp[i].p - cp[i].p);
            pp = std::max(pp, d);
            far += d > 1e-3;  // a different closest triangle (ties / float) shows as a jump
            pn = std::max(pn, norm(gp[i].n - cp[i].n));
        }
        std::printf("  piece %d: at: max |dp| %.2e mm, |dn| %.2e; project: max |dp| %.2e mm (%d of %zu > 1e-3 mm), |dn| %.2e; %s\n",
                    piece, dp, dn, pp, far, pts.size(), pn, st.line().c_str());
        CHECK(nonFinite == 0);
        CHECK(dp < 1e-2);
    }
}

// knitPiece with the device queries against the CPU knitPiece: program, layout and controls.
TEST(placement_knit_piece_gpu_vs_cpu) {
    REQUIRE_GPU();
    Teddy t;
    knit_test::requireScene(teddyManifest());
    CHECK(loadTeddy(t));
    // Every piece in library mm where it stands, for the pole tests (sewn_pieces.cpp places them the same way).
    std::vector<LoadedPiece> loaded(t.man.pieces.size());
    std::vector<PlacedPiece> placed(t.man.pieces.size());
    std::map<std::string, int> objects;
    for (size_t k = 0; k < t.man.pieces.size(); ++k) {
        std::string err;
        CHECK(loadPiece(t.man, static_cast<int>(k), t.sigma, loaded[k], &err));
        const int object = objects.emplace(loaded[k].entry.object, static_cast<int>(objects.size())).first->second;
        placed[k] = {&loaded[k].mesh, loaded[k].map.centre * t.sigma, object};
    }
    for (int piece : {17, 0}) {
        PatternParams pp;
        pp.yarnDiameter = 1.16;
        PieceKnitParams kp;
        std::string err;
        ScenePieceKnit a, b;
        const auto t0 = std::chrono::steady_clock::now();
        const ScenePieceRequest request{placed[piece], loaded[piece].entry.rings, placed};
        CHECK(knitScenePiece(request, {pp, kp}, knit_test::cpuNumerics(), a, &err));
        const double msCpu = ms(t0);
        std::unique_ptr<KnitNumerics> numerics = gpu::makeKnitNumerics(&err);
        CHECK(numerics != nullptr);
        const auto t1 = std::chrono::steady_clock::now();
        CHECK(knitScenePiece(request, {pp, kp}, *numerics, b, &err));
        const double msGpu = ms(t1);
        const PieceKnit &A = a.knit, &B = b.knit;
        const bool sameProgram = A.counts == B.counts && A.prog.stitches.size() == B.prog.stitches.size();
        double maxD = 0, sumD = 0;
        size_t n = 0;
        int nonFinite = 0;
        if (sameProgram && A.curves.yarns.size() == B.curves.yarns.size())
            for (size_t y = 0; y < A.curves.yarns.size(); ++y) {
                const auto &P = A.curves.yarns[y].curve.ctrl, &Q = B.curves.yarns[y].curve.ctrl;
                if (P.size() != Q.size()) continue;
                for (size_t c = 0; c < P.size(); ++c) {
                    nonFinite += !std::isfinite(Q[c].x + Q[c].y + Q[c].z);
                    const double d = norm(P[c] - Q[c]);
                    maxD = std::max(maxD, d);
                    sumD += d;
                    ++n;
                }
            }
        std::printf("  piece %d: program %s (%zu / %zu stitches, %zu / %zu turns), dir %+.0f / %+.0f; controls %zu: mean |d| %.2e, max %.2e mm; "
                    "CPU %.0f ms, GPU-query %.0f ms (%s)\n",
                    piece, sameProgram ? "identical" : "DIFFERENT", A.prog.stitches.size(), B.prog.stitches.size(), A.prog.turns.size(),
                    B.prog.turns.size(), A.dir, B.dir, n, n ? sumD / n : 0.0, maxD, msCpu, msGpu, gpu::knitPlacementStats(*numerics).c_str());
        CHECK(nonFinite == 0);
        CHECK(checkProgram(B.prog).status == Status::Pass);
    }
}

// The device's placement and the hand of a spiral piece before settle, bit for bit (pinned at 05960b1): the
// device's batched chart queries and store, on the fixture piece of spiral_placement_pinned, gathered and
// threaded ends.
TEST(gpu_spiral_placement_pinned) {
    REQUIRE_GPU();
    std::string err;
    std::unique_ptr<KnitNumerics> numerics = gpu::makeKnitNumerics(&err);
    CHECK(numerics != nullptr);
    FixturePiece f = fixturePiece(25.0);
    const TriMesh cover = octahedronAround(f.chart->poleB(), 3.0);
    const PlacedPiece piece{&f.mesh, {0, 0, 0}, 0};
    const std::vector<PlacedPiece> alone = {piece}, covered = {piece, {&cover, {0, 0, 0}, 0}};
    PieceKnitParams kp;
    kp.hand = everyChannelHand();
    const uint64_t pinned[2] = {0x36c96fe7621b039eull, 0x9679e1bd967b2c50ull};
    uint64_t digest[2] = {0, 0};
    int k = 0;
    for (const std::vector<PlacedPiece>* others : {&alone, &covered}) {
        ScenePieceKnit spk;
        CHECK(knitScenePiece({piece, 5, *others}, {PatternParams{}, kp}, *numerics, spk, &err));
        digest[k] = fnv1a64(placementText(spk.knit, spk.hand.tensions));
        std::printf("    %s end: %zu stitches, close ring %.3f mm, digest %016llx\n", k == 0 ? "gathered" : "threaded",
                    spk.knit.prog.stitches.size(), spk.knit.closeRing, static_cast<unsigned long long>(digest[k]));
        ++k;
    }
    CHECK(digest[0] == pinned[0] && digest[1] == pinned[1]);
}
