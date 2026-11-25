// Fibre CPU plan (2026-10-23) tests: candidate sources, budget-ordered generation, fibre cache
// coverage, region queries and exposed-surface placement. A separate executable
// (knit_tests_fp, CTest names fp_<case>) so the knitting work in tests.cpp merges without
// conflicts. Run one case: knit_tests_fp <name>.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "harness.h"
#include "cpu_numerics.h"
#include "fibres/fur.h"
#include "fibres/fur_field.h"
#include "fibres/fur_loops.h"
#include "knitting/tube/param_text.h"
#include "model/coordinates.h"
#include "model/hash.h"
#include "model/cancel.h"
#include "knitting/layout/chart.h"
#include "fibres/fibre_budget.h"
#include "model/mesh_check.h"
#include "model/mesh_query.h"
#include "fibres/fibre.h"
#include "fibres/fibre_checks.h"
#include "fibres/fibre_generator.h"
#include "knitting/knit/piece_knit.h"
#include "knitting/relax/settle.h"
#include "knitting/tube/yarn_surface.h"

using namespace knit;
using knit_test::registry;

namespace {

// ---- Synthetic budget pieces: candidates with chosen keys, a deterministic accept/reject and
// segment count, and an evaluator that returns the accepted ones as fibres (key, candidate,
// yarnId, pts of segs + 1 points).
struct SynthPiece {
    CandidateSource src;
    std::vector<std::vector<float>> key;   // [strand][i]
    std::vector<std::vector<int>> segs;    // 0: rejected
    int evaluations = 0;                   // candidates evaluated
    FibreSet all() const {
        FibreSet fs;
        for (int s = 0; s < src.strands(); ++s)
            for (int i = 0; i < src.count[s]; ++i)
                if (segs[s][i]) add(fs, s, i);
        return fs;
    }
    void add(FibreSet& fs, int s, int i) const {
        Fibre f;
        f.strand = s;
        f.yarnId = src.yarnId[s];
        f.candidate = i;
        f.key = key[s][i];  // exactly representable: float(f.key) == key
        std::vector<Vec3> pts;
        for (int k = 0; k <= segs[s][i]; ++k) pts.push_back(Vec3{double(s), double(i), double(k)});
        fs.addFibre(f, pts);
    }
};

uint64_t mix(uint64_t x) {
    x += 0x9e3779b97f4a7c15ull;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
    return x ^ (x >> 31);
}
double u01(uint64_t a, uint64_t b) { return double(mix(mix(a) ^ b) >> 11) * (1.0 / 9007199254740992.0); }

// `quant` > 0 quantizes keys to multiples of 1 / quant (ties within and across pieces); keys of
// `oneEvery` > 0 candidates round to 1.0f.
SynthPiece synthPiece(uint64_t seed, int strands, int perStrand, double accept, int quant, int oneEvery) {
    SynthPiece p;
    for (int s = 0; s < strands; ++s) {
        const int n = perStrand + int(u01(seed, 1000 + s) * perStrand);
        p.src.yarnId.push_back(s);
        p.src.count.push_back(n);
        p.key.emplace_back();
        p.segs.emplace_back();
        for (int i = 0; i < n; ++i) {
            double k = u01(seed, uint64_t(s) << 32 | uint64_t(i));
            if (quant > 0) k = std::floor(k * quant) / quant;
            float fk = float(k);
            if (oneEvery > 0 && (i % oneEvery) == oneEvery - 1) fk = float(1.0 - 1e-12);  // rounds to 1.0f
            p.key.back().push_back(fk);
            p.segs.back().push_back(u01(seed ^ 77, uint64_t(s) << 32 | uint64_t(i)) < accept ? 2 + int(u01(seed ^ 99, i) * 3) : 0);
        }
    }
    return p;
}

int polySegs(const Fibre& f) { return f.segments(); }

// A piece's candidates per bucket, from its chosen keys (ties, keys that round to 1.0f).
std::vector<int64_t> perBucket(const SynthPiece& p) {
    std::vector<int64_t> count(kKeyBuckets, 0);
    for (int s = 0; s < p.src.strands(); ++s)
        for (int i = 0; i < p.src.count[s]; ++i) ++count[keyBucket(p.key[s][i])];
    return count;
}

// The piece's candidates with a key in the window, evaluated: its accepted ones as fibres.
FibreSet evaluate(SynthPiece& p, const KeyWindow& window) {
    FibreSet fs;
    for (int s = 0; s < p.src.strands(); ++s)
        for (int i = 0; i < p.src.count[s]; ++i) {
            if (!(p.key[s][i] >= window.keyLo && p.key[s][i] < window.keyHi)) continue;
            ++p.evaluations;
            ++fs.rep.candidates;
            if (p.segs[s][i]) {
                p.add(fs, s, i);
                ++fs.rep.accepted;
            }
        }
    return fs;
}

// The schedule driven as the scene drives it: a window, every piece's fibres in it, until finished.
BudgetRun runSchedule(const std::vector<SynthPiece*>& ps, long long budget) {
    std::vector<std::vector<int64_t>> buckets;
    for (const SynthPiece* p : ps) buckets.push_back(perBucket(*p));
    BudgetSchedule schedule(std::move(buckets), budget);
    while (!schedule.finished()) {
        const KeyWindow& window = schedule.nextWindow();
        std::vector<FibreSet> fibres(ps.size());
        for (size_t k = 0; k < ps.size(); ++k)
            if (window.candidates[k] > 0) fibres[k] = evaluate(*ps[k], window);
        schedule.addWindow(std::move(fibres));
    }
    return schedule.finish();
}

// Full-generation reference: q over every accepted fibre, and the kept fibres per piece.
float fullThreshold(const std::vector<SynthPiece*>& ps, long long budget) {
    std::vector<std::pair<float, int>> ks;
    for (const SynthPiece* p : ps)
        for (const Fibre& f : p->all().fibres) ks.push_back({float(f.key), polySegs(f)});
    return budgetThreshold(ks, budget);
}

// The budget run reproduces full generation for `budget`: same q, the fibres below q are
// exactly the full set's (same order), every piece's output is exactly the accepted fibres
// with key < keyBound, and the counters cover exactly those candidates.
void checkAgainstFull(std::vector<SynthPiece*> ps, long long budget, BudgetRun* out = nullptr) {
    const BudgetRun run = runSchedule(ps, budget);
    const float q = fullThreshold(ps, budget);
    if (run.q != q) std::printf("  budget %lld: q %.9g, full %.9g\n", budget, double(run.q), double(q));
    CHECK(run.q == q);
    CHECK(run.keyBound >= run.q || run.complete);
    for (size_t k = 0; k < ps.size(); ++k) {
        const BudgetPartResult& r = run.parts[k];
        std::vector<const Fibre*> want;
        const FibreSet all = ps[k]->all();
        for (const Fibre& f : all.fibres)
            if (float(f.key) < r.keyBound) want.push_back(&f);
        CHECK(r.fs.fibres.size() == want.size());
        for (size_t i = 0; i < want.size(); ++i) {
            CHECK(r.fs.fibres[i].yarnId == want[i]->yarnId && r.fs.fibres[i].candidate == want[i]->candidate);
            CHECK(r.fs.fibres[i].key == want[i]->key);
        }
        int covered = 0;
        for (int s = 0; s < ps[k]->src.strands(); ++s)
            for (int i = 0; i < ps[k]->src.count[s]; ++i) covered += ps[k]->key[s][i] < r.keyBound;
        CHECK(r.covered == covered && r.fs.rep.candidates == covered);
        CHECK(r.fs.rep.accepted == static_cast<int>(want.size()));
    }
    if (q < 2.0f) CHECK(run.segmentsBelowQ <= budget);
    if (out) *out = run;
}

} // namespace

TEST(budget_sched_matches_full_synthetic) {
    // Three pieces, ties within and across pieces (keys quantized to 1/4096), keys rounding to
    // 1.0f; budgets from 1 segment to beyond the total, exact exhaustion, unlimited.
    SynthPiece a = synthPiece(1, 7, 300, 0.4, 4096, 97), b = synthPiece(2, 3, 900, 0.25, 4096, 0), c = synthPiece(3, 11, 50, 0.9, 4096, 13);
    std::vector<SynthPiece*> ps{&a, &b, &c};
    long long total = 0;
    for (SynthPiece* p : ps)
        for (const Fibre& f : p->all().fibres) total += polySegs(f);
    for (long long budget : {1LL, 7LL, 100LL, total / 10, total / 4, total / 2, total - 3, total, total + 5, 0LL, -1LL}) {
        for (SynthPiece* p : ps) p->evaluations = 0;
        BudgetRun run;
        checkAgainstFull(ps, budget, &run);
        int evals = 0;
        for (SynthPiece* p : ps) evals += p->evaluations;
        std::printf("  budget %lld of %lld: q %.6f, %d batches, evaluated %d of %lld%s\n", budget, total, double(run.q), run.batches, evals,
                    static_cast<long long>(run.candidates), run.complete ? ", complete" : "");
        CHECK(evals == run.evaluated);
        if (budget <= 0 || budget >= total) CHECK(run.complete && run.q == 2.0f);
        if (budget > 0 && budget <= total / 4) CHECK(run.evaluated < run.candidates);  // the point of it
    }
}

TEST(budget_sched_exact_exhaustion) {
    // A budget equal to the total accepted segments is not crossed (strict >): q = 2, and every
    // candidate is evaluated before that is known.
    SynthPiece a = synthPiece(11, 5, 200, 0.5, 0, 0);
    long long total = 0;
    for (const Fibre& f : a.all().fibres) total += polySegs(f);
    BudgetRun run;
    checkAgainstFull({&a}, total, &run);
    CHECK(run.q == 2.0f && run.complete && run.evaluated == run.candidates);
    checkAgainstFull({&a}, total - 1, &run);
    CHECK(run.q < 2.0f);
}

TEST(budget_sched_ties_across_pieces) {
    // Every candidate of two pieces has one of 8 keys: whole equal-key groups straddle the budget,
    // across pieces; the group at q is excluded entirely (key < q), as in full generation.
    SynthPiece a = synthPiece(21, 4, 100, 0.6, 8, 0), b = synthPiece(22, 4, 100, 0.6, 8, 0);
    long long total = 0;
    for (SynthPiece* p : {&a, &b})
        for (const Fibre& f : p->all().fibres) total += polySegs(f);
    for (long long budget = 1; budget < total; budget += std::max(1LL, total / 37)) checkAgainstFull({&a, &b}, budget);
}

TEST(budget_sched_key_rounds_to_one) {
    // Keys that round to 1.0f sit in their own top bucket: a budget that reaches them keeps
    // everything below 1.0f and, as full generation, q = 1.0f excludes them.
    SynthPiece a = synthPiece(31, 3, 200, 1.0, 0, 2);  // every second key is 1.0f
    long long below = 0, total = 0;
    for (const Fibre& f : a.all().fibres) {
        total += polySegs(f);
        if (float(f.key) < 1.0f) below += polySegs(f);
    }
    BudgetRun run;
    checkAgainstFull({&a}, below + 1, &run);
    CHECK(run.q == 1.0f);
    CHECK(keyBucket(1.0f) == kKeyBuckets - 1 && bucketLowerKey(kKeyBuckets - 1) == 1.0f);
    checkAgainstFull({&a}, total, &run);
    CHECK(run.q == 2.0f);
}

TEST(budget_sched_empty_and_rejection_heavy) {
    // Keys only in a few distant buckets, and a piece that accepts almost nothing (most batches
    // add no segments): batches skip empty buckets and still stop exactly.
    SynthPiece a = synthPiece(41, 2, 400, 0.02, 4, 0), b = synthPiece(42, 6, 300, 0.01, 0, 0);
    long long total = 0;
    for (SynthPiece* p : {&a, &b})
        for (const Fibre& f : p->all().fibres) total += polySegs(f);
    for (long long budget : {1LL, total / 3, total / 2, total - 1, total}) checkAgainstFull({&a, &b}, budget);
}

// Acceptance falls with the key, so the pilot overestimates the segments per candidate and later
// windows undershoot: the schedule needs several windows, and the doubling rule decides their size.
SynthPiece skewedPiece(uint64_t seed, int strands, int perStrand) {
    SynthPiece p;
    for (int s = 0; s < strands; ++s) {
        p.src.yarnId.push_back(s);
        p.src.count.push_back(perStrand);
        p.key.emplace_back();
        p.segs.emplace_back();
        for (int i = 0; i < perStrand; ++i) {
            const double k = u01(seed, uint64_t(s) << 32 | uint64_t(i));
            p.key.back().push_back(float(k));
            const double accept = std::pow(1 - k, 6);
            p.segs.back().push_back(u01(seed ^ 77, uint64_t(s) << 32 | uint64_t(i)) < accept ? 2 + int(u01(seed ^ 99, i) * 3) : 0);
        }
    }
    return p;
}

TEST(budget_sched_windows_pinned) {
    // The windows, q and coverage of the budget schedule before it became data (bcecbaf's budgetSchedule,
    // its log on these pieces): the covered fibres depend on the windows, and the cache and the digests
    // store them. 20000 never crosses and takes four windows; the third is sized by the doubling rule.
    struct Pinned {
        long long budget;
        std::vector<std::pair<int, int>> windows;
        float q, keyBound;
    };
    const Pinned pinned[] = {
        {3000, {{0, 15870}}, 0x1.1318b0p-5f, 0x1.fcp-4f},
        {9000, {{0, 15870}, {15870, 15994}}, 0x1.099682p-3f, 0x1.f4p-3f},
        {20000, {{0, 15870}, {15870, 16024}, {16024, 16164}, {16164, 16257}}, 2.0f, 2.0f},
    };
    for (const Pinned& want : pinned) {
        SynthPiece a = skewedPiece(61, 8, 3000), b = skewedPiece(62, 5, 2000);
        BudgetSchedule schedule({perBucket(a), perBucket(b)}, want.budget);
        std::vector<std::pair<int, int>> windows;
        while (!schedule.finished()) {
            const KeyWindow& window = schedule.nextWindow();
            windows.push_back({window.firstBucket, window.endBucket});
            schedule.addWindow({evaluate(a, window), evaluate(b, window)});
        }
        const BudgetRun run = schedule.finish();
        std::printf("  budget %lld: %zu windows, q %a, keyBound %a\n", want.budget, windows.size(), double(run.q), double(run.keyBound));
        CHECK(windows == want.windows);
        CHECK(run.q == want.q && run.keyBound == want.keyBound);
    }
}

TEST(budget_sched_protocol_enforced) {
    // The schedule refuses what would change q silently: a window added twice, a missing part's fibres,
    // a part that evaluated other candidates than the window holds, and a second finish.
    SynthPiece a = synthPiece(51, 4, 2000, 0.5, 0, 0), b = synthPiece(52, 3, 500, 0.5, 0, 0);
    BudgetSchedule schedule({perBucket(a), perBucket(b)}, 1000);
    const KeyWindow window = schedule.nextWindow();
    bool refused = false;
    try {
        schedule.nextWindow();  // a window is pending
    } catch (const std::logic_error&) {
        refused = true;
    }
    CHECK(refused);
    refused = false;
    try {
        schedule.addWindow({evaluate(a, window)});  // one set for two parts
    } catch (const std::invalid_argument&) {
        refused = true;
    }
    CHECK(refused);
    refused = false;
    try {
        schedule.addWindow({FibreSet{}, evaluate(b, window)});  // a's candidates missing
    } catch (const std::runtime_error&) {
        refused = true;
    }
    CHECK(refused);
    refused = false;
    try {
        schedule.addWindow({evaluate(a, window), evaluate(b, window)});  // the window was consumed by the refusal
    } catch (const std::logic_error&) {
        refused = true;
    }
    CHECK(refused);
    BudgetSchedule whole({perBucket(a)}, 0);  // no budget: one window holds everything
    const KeyWindow all = whole.nextWindow();
    CHECK(all.firstBucket == 0 && all.endBucket == kKeyBuckets && all.candidates[0] == a.src.total());
    whole.addWindow({evaluate(a, all)});
    CHECK(whole.finished());
    whole.finish();
    refused = false;
    try {
        whole.finish();
    } catch (const std::logic_error&) {
        refused = true;
    }
    CHECK(refused);
}

// ---- The real generator: a knitted, settled sphere (tests.cpp's fixture, rebuilt here).
namespace {
struct Sphere {
    TriMesh mesh;
    std::unique_ptr<Chart> chart;
    std::unique_ptr<PieceSurface> surf;
    PatternParams pp;
    PieceKnit pk;
    SurfaceModel sm;
    SceneMapping map;
    std::unique_ptr<ExportTubeIndex> tubes;
    double eSer = 1e-4;
};
const Sphere& sphere() {
    static std::unique_ptr<Sphere> c;
    if (c) return *c;
    c = std::make_unique<Sphere>();
    std::string err;
    if (!readObj(std::string(KNIT_TEST_DATA) + "/chart_fixture.obj", c->mesh, &err)) throw std::runtime_error(err);
    c->mesh.transform(25.0, c->mesh.centroid());
    c->mesh.computeNormals();
    c->chart = std::make_unique<Chart>(c->mesh, 5);
    c->surf = std::make_unique<PieceSurface>(*c->chart);
    PieceKnitParams kp;
    c->pp.waleSpacing = kp.g.w;
    c->pp.courseSpacing = kp.g.h;
    CHECK(knitPiece({*c->surf, true}, {c->pp, kp}, knit_test::cpuNumerics(), c->pk));
    settle({c->pk.curves, c->pp, c->surf.get()}, SettleParams{}, knit_test::cpuNumerics());
    SurfaceParams sp;
    sp.compression = false;
    c->sm = buildSurface(c->pk.curves, c->pp, sp);
    c->map.sigma = 68.103;
    c->map.centre = {0.3, -1.2, 0.8};
    const SceneMapping map = c->map;
    c->tubes = std::make_unique<ExportTubeIndex>(c->pk.curves, c->pp, ExportTubeSpec{kExportModelRing, kExportRing, kExportSides, false,
                                                 [map](const Vec3& p) { return map.toExportF(p); }, c->eSer});
    return *c;
}
// The pieces' frame and parameters on the sphere, piece seed `piece`.
void sphereSetup(const Sphere& s, size_t piece, double density, FibreParams& fp, FibreFrame& fr) {
    fp = pieceFibres();
    fp.density = density;
    fp.serializedMargin = s.eSer;
    fp.exportChordFactor = kExportRing * kExportRing / 8;
    FibreSetup fset;
    fset.tubeRadius = s.pp.radius();
    fset.sides = kExportSides;
    fset.ring = kExportRing;
    fset.kappaMax = 1.0;
    fset.pieceIndex = piece;
    setupPieceFibres(fp, fr, s.pp, fset);
    const PieceSurface* surf = s.surf.get();
    fr.fabric = [surf](const Vec3& p) {
        const SurfaceSample q = surf->project(p);
        return FabricPoint{q.p, q.n};
    };
    fr.exec.parallel = true;
}
} // namespace

// ---- P2: exact query reuse, position-only evaluation, region states, shared support.
TEST(yarn_position_equals_coordinates) {
    // yarnPosition is bit-identical to yarnCoordinates(...).X (ring boundaries, ends, lambda != 1).
    const Sphere& s = sphere();
    long n = 0;
    for (const YarnSurface& ys : s.sm.yarns)
        for (int k = 0; k < 400; ++k) {
            double sa = k < 4 ? (k < 2 ? 0.0 : ys.length) : k < 40 ? ys.rings[static_cast<size_t>(k) % ys.rings.size()].s : u01(ys.yarnId, k) * ys.length;
            const double th = 2 * kPi * u01(ys.yarnId + 7, k), lam = 0.8 + 0.6 * u01(ys.yarnId + 9, k);
            const Vec3 a = yarnCoordinates(s.pk.curves, ys, sa, th, lam).X, b = yarnPosition(s.pk.curves, ys, sa, th, lam);
            CHECK(std::memcmp(&a, &b, sizeof a) == 0);
            ++n;
        }
    std::printf("  %ld points bit-identical\n", n);
}

namespace {
// A closed, consistently oriented UV sphere (radius R, centre c).
TriMesh uvSphere(double R, const Vec3& c, int nu, int nv) {
    TriMesh m;
    m.pos.push_back(c + Vec3{0, 0, R});
    for (int i = 1; i < nv; ++i)
        for (int j = 0; j < nu; ++j) {
            const double th = kPi * i / nv, ph = 2 * kPi * j / nu;
            m.pos.push_back(c + Vec3{R * std::sin(th) * std::cos(ph), R * std::sin(th) * std::sin(ph), R * std::cos(th)});
        }
    m.pos.push_back(c + Vec3{0, 0, -R});
    const int S = static_cast<int>(m.pos.size()) - 1;
    auto v = [&](int i, int j) { return 1 + (i - 1) * nu + ((j % nu) + nu) % nu; };
    for (int j = 0; j < nu; ++j) m.tri.push_back({0, v(1, j), v(1, j + 1)});
    for (int i = 1; i < nv - 1; ++i)
        for (int j = 0; j < nu; ++j) {
            m.tri.push_back({v(i, j), v(i + 1, j), v(i + 1, j + 1)});
            m.tri.push_back({v(i, j), v(i + 1, j + 1), v(i, j + 1)});
        }
    for (int j = 0; j < nu; ++j) m.tri.push_back({S, v(nv - 1, j + 1), v(nv - 1, j)});
    return m;
}
} // namespace

TEST(support_inputs_equal_reuse) {
    // ExportTubeIndex takes a caller's surface model only when the digests of the build inputs are
    // equal, and then answers exactly as with its own model; other inputs (ring spacing,
    // compression, curves) are not reused.
    const Sphere& s = sphere();
    SurfaceParams sp;
    sp.compression = false;
    sp.ringSpacing = 0.025;
    const uint64_t d = surfaceBuildDigest(s.pk.curves, s.pp, sp);
    const SceneMapping map = s.map;
    auto toE = [map](const Vec3& p) { return map.toExportF(p); };
    const ExportTubeSpec spec{0.025, kExportRing, kExportSides, false, toE, s.eSer};
    ExportTubeSpec coarseSpec = spec, compressedSpec = spec;
    coarseSpec.modelRingSpacing = 0.2;
    compressedSpec.compression = true;
    const ExportTubeIndex own(s.pk.curves, s.pp, spec);
    const ExportTubeIndex shared(s.pk.curves, s.pp, spec, {}, &s.sm, d);
    const ExportTubeIndex coarse(s.pk.curves, s.pp, coarseSpec, {}, &s.sm, d);
    SurfaceParams spc = sp;
    spc.compression = true;
    const ExportTubeIndex comp(s.pk.curves, s.pp, compressedSpec, {}, &s.sm, d);
    CHECK(shared.sharedSurface() && !own.sharedSurface() && !coarse.sharedSurface() && !comp.sharedSurface());
    CHECK(surfaceBuildDigest(s.pk.curves, s.pp, spc) != d);
    PatternParams pp2 = s.pp;
    pp2.yarnDiameter *= 1.01;
    CHECK(surfaceBuildDigest(s.pk.curves, pp2, sp) != d);
    CurveSet cs2 = s.pk.curves;
    cs2.yarns[0].curve.ctrl[3].x += 1e-9;
    CHECK(surfaceBuildDigest(cs2, s.pp, sp) != d);
    // Same answers on real fibres (accepts and depth).
    FibreParams fp;
    FibreFrame fr;
    sphereSetup(s, 0, 4.0, fp, fr);
    const FibreSet fs = growYarnFibres(s.pk.curves, s.sm, s.pp, fp, &fr);
    for (const Fibre& f : fs.fibres) {
        const Vec3 root = toVec3(fs.pointsOf(f)[0]);
        CHECK(own.accepts(fs, f, 1e-4) == shared.accepts(fs, f, 1e-4));
        CHECK(own.depth(f.yarnId, f.s - 0.5, f.s + 0.5, root) == shared.depth(f.yarnId, f.s - 0.5, f.s + 0.5, root));
    }
    std::printf("  %zu fibres: identical answers with the shared model\n", fs.fibres.size());
}

// ---- Pins: what must keep its exact bits while the fibre code is rewritten.
namespace {
uint64_t pointsHash(const std::vector<Vec3>& points, uint64_t h) {
    char b[96];
    for (const Vec3& p : points) {
        std::snprintf(b, sizeof b, "%a %a %a\n", p.x, p.y, p.z);
        h = fnv1a64(b, h);
    }
    return h;
}
} // namespace

TEST(fibre_params_text_pinned) {
    // The text behind every saved fibre result's key: a change here refuses all saved fibres. Light
    // checks and curled strays add their own words, so all three spellings are pinned.
    FibreParams fp = pieceFibres();
    const uint64_t plain = fnv1a64(fibreParamsText(fp));
    fp.lightChecks = true;
    const uint64_t light = fnv1a64(fibreParamsText(fp));
    fp.randomCurlRadius = 1.5;
    const uint64_t curled = fnv1a64(fibreParamsText(fp));
    std::printf("  %016llx %016llx %016llx\n", static_cast<unsigned long long>(plain), static_cast<unsigned long long>(light),
                static_cast<unsigned long long>(curled));
    CHECK(plain == 0xf3260c710d965b0eull);
    CHECK(light == 0x3a433e6e67b5d7aeull);
    CHECK(curled == 0x268b3dc98c696c35ull);
}

// Every coordinate of a small pile (both layers, settled) on a ball, at the approved look's shape: the Fur
// workflow's curves (fur_field.h, fur_loops.h) pinned. The scene-level reproduction of the frozen baseline
// (gate/top15, SHA-256 of every strand file) is the M2 gate; this pins the unit path.
TEST(fur_strands_pinned) {
    TriMesh ball = uvSphere(20.0, {0, 0, 0}, 40, 30);
    ball.computeNormals();
    const MeshQuery ballQuery(ball);
    FurLook look;
    look.loops.density = 0.2;
    look.fibres.density = 0.1;
    const FurRootCache loops = sampleFurRoots(ball, FurLayer::Loops, look.loops.density, furPartSeed(look.seed, FurLayer::Loops, 0), nullptr, ExecPolicy{true});
    const FurRootCache fibres = sampleFurRoots(ball, FurLayer::Fibres, look.fibres.density, furPartSeed(look.seed, FurLayer::Fibres, 0), nullptr, ExecPolicy{true});
    const std::vector<Hair> coat = growLoops(buildFurField({loops, fibres, look, nullptr}, ExecPolicy{true}), ballQuery, look.settle, ExecPolicy{true}, nullptr);
    uint64_t furHash = fnv1a64("fur");
    for (const Hair& h : coat) furHash = pointsHash(h.points, furHash);
    std::printf("  %zu fibres %016llx\n", coat.size(), static_cast<unsigned long long>(furHash));
    CHECK(coat.size() > 500);
    CHECK(furHash == 0xf142dbb9b1804e80ull);
}

int main(int argc, char** argv) {
    if (argc > 1 && std::strcmp(argv[1], "--list") == 0) {
        for (const knit_test::Case& c : registry()) std::printf("%s\n", c.name);
        return 0;
    }
    int failed = 0, ran = 0, skipped = 0;
    for (const knit_test::Case& c : registry()) {
        if (argc > 1 && std::strcmp(argv[1], c.name) != 0) continue;
        ++ran;
        try {
            c.fn();
            std::printf("[ OK ] %s\n", c.name);
        } catch (const knit_test::Skipped& s) {
            ++skipped;
            std::printf("[SKIP] %s: %s\n", c.name, s.why);
        } catch (const std::exception& e) {
            ++failed;
            std::printf("[FAIL] %s: %s\n", c.name, e.what());
        }
    }
    if (ran == 0) {
        std::printf("no test named %s\n", argc > 1 ? argv[1] : "?");
        return 2;
    }
    if (failed) return 1;
    return skipped == ran ? 77 : 0;  // 77: every case run was skipped (CTest SKIP_RETURN_CODE)
}
