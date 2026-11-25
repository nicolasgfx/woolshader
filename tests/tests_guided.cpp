// The guided Knit route (knit.route = guided; scene/source/guided_source.h, scene/knit/guided_part.cpp), part of
// knit_tests_workflow: a two-capsule OBJ imported as the file, markings on its first capsule, the route's refusals, its
// keys and stage boundaries, a refused construction and the cancellation of a running settle. Each case fails with its
// defect back (the comments say which).
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <thread>

#include "harness.h"
#include "import/import_analysis.h"
#include "import/import_store.h"
#include "import/prepare_import.h"
#include "io/cache_paths.h"
#include "io/marking_file.h"
#include "jobs/job_steps.h"
#include "model/cancel.h"
#include "model/hash.h"
#include "model/surface_path.h"
#include "model/surface_pieces.h"
#include "gpu/knitting.h"
#include "scene/knit/knit_report.h"
#include "scene/knit/knit_run.h"
#include "scene/knit/knit_stage.h"
#include "scene/result_keys.h"
#include "scene/scene.h"
#include "scene/source/scene_parts.h"
#include "settings/scene_registry.h"
#include "settings/settings_profiles.h"
#include "test_compute.h"

using namespace knit;
namespace fs = std::filesystem;

namespace {

constexpr double kTestPi = 3.14159265358979323846;
constexpr double kRadius = 10.5;  // 12 stitches of 5.5 mm round

std::string guidedDir(const std::string& name) {
    const fs::path p = fs::temp_directory_path() / "knit_guided_tests" / name;
    std::error_code ec;
    fs::remove_all(p, ec);
    fs::create_directories(p, ec);
    return p.generic_string();
}

// Two capsules round y (Y up, mm), the second 40 mm along x: two objects, so two pieces.
void writeTwoCapsules(const std::string& path) {
    std::ofstream f(path, std::ios::binary);
    const int segments = 48, capRings = 8;
    const double length = 24;
    int base = 0;
    for (int object = 0; object < 2; ++object) {
        const double x0 = 40.0 * object;
        f << "o " << (object == 0 ? "left" : "right") << "\n";
        std::vector<double> ringY, ringR;
        for (int i = 1; i < capRings; ++i) {
            const double a = -kTestPi / 2 + kTestPi / 2 * i / capRings;
            ringY.push_back(kRadius * std::sin(a));
            ringR.push_back(kRadius * std::cos(a));
        }
        for (int i = 0; i <= 16; ++i) {
            ringY.push_back(length * i / 16);
            ringR.push_back(kRadius);
        }
        for (int i = 1; i < capRings; ++i) {
            const double a = kTestPi / 2 * i / capRings;
            ringY.push_back(length + kRadius * std::sin(a));
            ringR.push_back(kRadius * std::cos(a));
        }
        f << "v " << x0 << " " << -kRadius << " 0\n";
        for (size_t j = 0; j < ringY.size(); ++j)
            for (int i = 0; i < segments; ++i) {
                const double a = 2 * kTestPi * i / segments;
                f << "v " << x0 + ringR[j] * std::cos(a) << " " << ringY[j] << " " << -ringR[j] * std::sin(a) << "\n";
            }
        f << "v " << x0 << " " << length + kRadius << " 0\n";
        const int rings = static_cast<int>(ringY.size()), top = 1 + rings * segments;
        auto v = [&](int j, int i) { return base + 2 + j * segments + (i % segments); };  // OBJ indices start at 1
        for (int i = 0; i < segments; ++i) f << "f " << base + 1 << " " << v(0, i + 1) << " " << v(0, i) << "\n";
        for (int j = 0; j + 1 < rings; ++j)
            for (int i = 0; i < segments; ++i) {
                f << "f " << v(j, i) << " " << v(j, i + 1) << " " << v(j + 1, i + 1) << "\n";
                f << "f " << v(j, i) << " " << v(j + 1, i + 1) << " " << v(j + 1, i) << "\n";
            }
        for (int i = 0; i < segments; ++i) f << "f " << v(rings - 1, i) << " " << v(rings - 1, i + 1) << " " << base + 1 + top << "\n";
        base += top + 1;
    }
}

// A marking of `meaning` through file points (piece 0's nearest surface points).
SurfaceMarking markingThrough(const SurfaceFile& file, int id, MarkingMeaning meaning, MarkingShape shape, const std::vector<Vec3>& points) {
    const SurfacePiece& piece = file.pieces->pieces[0];
    SurfaceMarking m;
    m.id = id;
    m.meaning = meaning;
    m.shape = shape;
    for (const Vec3& q : points) {
        SurfacePoint p = closestSurfacePoint(piece.mesh, q);
        p.tri = piece.fileTriangle[p.tri];
        m.points.push_back(p);
    }
    return m;
}

std::vector<Vec3> ring(double y, int n = 6) {
    std::vector<Vec3> p;
    for (int i = 0; i < n; ++i) {
        const double a = 2 * kTestPi * i / n + 0.031;
        p.push_back({kRadius * std::cos(a), y, -kRadius * std::sin(a)});
    }
    return p;
}

// The band route's markings on the left capsule: a start ring, a finish ring and a row-start path.
MarkingSet capsuleMarkings(const SurfaceFile& file, const std::string& source) {
    MarkingSet set;
    set.source = source;
    set.digest = file.digest;
    set.markings.push_back(markingThrough(file, 1, MarkingMeaning::Start, MarkingShape::ClosedCurve, ring(4)));
    set.markings.push_back(markingThrough(file, 2, MarkingMeaning::Finish, MarkingShape::ClosedCurve, ring(18)));
    const double a = 0.013;
    set.markings.push_back(markingThrough(file, 3, MarkingMeaning::RowStartPath, MarkingShape::Path,
                                          {{kRadius * std::cos(a), 1, -kRadius * std::sin(a)}, {kRadius * std::cos(a), 23, -kRadius * std::sin(a)}}));
    set.nextId = 4;
    return set;
}

struct ImportCacheScope {
    explicit ImportCacheScope(const std::string& root) { _putenv_s("KNIT_IMPORT_CACHE", root.c_str()); }
    ~ImportCacheScope() { _putenv_s("KNIT_IMPORT_CACHE", ""); }
};

// The two-capsule file imported as it is (Keep file), with markings written beside it: settings on the guided route.
struct GuidedFixture {
    std::string dir, obj, markings;
    SurfaceFile file;
    MarkingSet set;
    SceneSettings s;
};

GuidedFixture guidedFixture(const std::string& name) {
    GuidedFixture g;
    g.dir = guidedDir(name);
    g.obj = g.dir + "/capsules.obj";
    writeTwoCapsules(g.obj);
    ImportedMesh m;
    std::string err;
    CHECK(readMeshFile(g.obj, ImportLimits{}, m, &err));
    ImportOptions o;
    o.unitMm = 1;
    o.up = UpAxis::Y;
    o.surface = ImportOptions::Surface::KeepFile;
    PreparedMesh p;
    CHECK(prepareImport(m, o, p, &err));
    const uint64_t id = importIdentity(m, o);
    CHECK(storeImport(g.dir + "/imports", {std::make_shared<const ImportedMesh>(m), o, {}, id}, {p, {}}, &err));
    CHECK(readSurfaceFile(g.obj, g.file, &err));
    CHECK(g.file.pieces->pieces.size() == 2);
    g.set = capsuleMarkings(g.file, g.obj);
    g.markings = g.dir + "/capsule.markings";
    CHECK(saveMarkingFile(g.markings, g.set, &err));
    g.s.scene = "import";
    g.s.importSource = g.obj;
    g.s.importIdentity = identityHex(id);
    g.s.importUnitMm = 1;
    g.s.importUp = "y";
    g.s.importSurface = "file";
    g.s.workflow = "knit";
    g.s.knitRoute = "guided";
    g.s.guidedPiece = 0;
    g.s.guidedMarkings = g.markings;
    g.s.guidedMarkingsDigest = identityHex(effectiveMarkingsDigest(g.set, *g.file.pieces, 0));
    return g;
}

bool loadGuided(const SceneSettings& s, Scene& scene, std::string& err) {
    const SceneEntry* e = sceneEntryFor(s);
    if (!e) {
        err = "no scene entry";
        return false;
    }
    return loadPartsToKnit(*e, s, SourceOptions{}, scene, &err).has_value();
}

} // namespace

// The guided scene: one Guided part (the left capsule) placed exactly where the import placed it, the right capsule a
// plain primitive. Refused, with the reason, for markings of another file, a piece that does not exist and markings
// edited after they were chosen; never rebound. (Defects: a placement without the import's offset moves the part; a
// missing digest check knits stale markings.)
TEST(guided_source_refuses_other_file) {
    GuidedFixture g = guidedFixture("source");
    ImportCacheScope cache(g.dir + "/imports");
    Scene scene;
    std::string err;
    CHECK(loadGuided(g.s, scene, err));
    std::printf("  %s\n", err.c_str());
    CHECK(scene.parts.size() == 1 && scene.parts[0].layout == KnitLayout::Guided);
    int plain = 0;
    for (const SourcePrimitive& p : scene.source.prims) plain += p.cls == SourcePrimitive::Class::Plain;
    CHECK(plain == 1 && scene.source.prims.size() == 2);
    // the import's own placement: the same file loaded as base parts puts the left capsule at the same library bounds
    Scene base;
    CHECK(loadBaseParts(*sceneEntryFor(g.s), g.s, SourceOptions{}, base, &err));
    const ScenePart* left = nullptr;  // the import's left capsule (the KeepFile import: one part per component)
    for (const ScenePart& r : base.parts)
        if (r.centre.x < 0) left = &r;
    CHECK(left != nullptr && base.parts.size() == 2);
    if (!left) return;
    Vec3 lo, hi;
    left->mesh->bounds(lo, hi);
    lo = lo + left->centre;
    hi = hi + left->centre;
    Vec3 glo, ghi;
    scene.parts[0].mesh->bounds(glo, ghi);
    glo = glo + scene.parts[0].centre;
    ghi = ghi + scene.parts[0].centre;
    std::printf("  guided part (%.4f %.4f %.4f)-(%.4f %.4f %.4f), import (%.4f %.4f %.4f)-(%.4f %.4f %.4f)\n", glo.x, glo.y, glo.z, ghi.x, ghi.y, ghi.z, lo.x, lo.y, lo.z,
                hi.x, hi.y, hi.z);
    CHECK(norm(glo - lo) < 1e-6 && norm(ghi - hi) < 1e-6);

    // markings of another file (another digest): refused, not rebound
    MarkingSet other = g.set;
    other.digest ^= 1;
    CHECK(saveMarkingFile(g.dir + "/other.markings", other, &err));
    SceneSettings s = g.s;
    s.guidedMarkings = g.dir + "/other.markings";
    Scene refused;
    CHECK(!loadGuided(s, refused, err) && err.find("another version") != std::string::npos);
    // a piece the file does not have
    s = g.s;
    s.guidedPiece = 2;
    CHECK(!loadGuided(s, refused, err) && err.find("does not exist") != std::string::npos);
    // the markings file edited after it was chosen (a user marking added on the knitted piece)
    MarkingSet edited = g.set;
    edited.markings.push_back(markingThrough(g.file, 4, MarkingMeaning::Boundary, MarkingShape::ClosedCurve, ring(12)));
    CHECK(saveMarkingFile(g.markings, edited, &err));
    CHECK(!loadGuided(g.s, refused, err) && err.find("changed after it was chosen") != std::string::npos);
}

// What the build uses decides the digest: offered and rejected proposals and markings on the other piece change nothing;
// accepting a proposal changes it, and with it the Knit key. The markings file's path is not a stage input (its digest
// is), and on the surface route the guided rows are no input at all. (Defect: hashing every marking.)
TEST(guided_ignores_offered_and_rejected) {
    GuidedFixture g = guidedFixture("proposals");
    const uint64_t base = effectiveMarkingsDigest(g.set, *g.file.pieces, 0);
    MarkingSet set = g.set;
    SurfaceMarking offered = markingThrough(g.file, 10, MarkingMeaning::Finish, MarkingShape::ClosedCurve, ring(14));
    offered.origin = MarkingOrigin::Proposal;
    offered.proposal = ProposalState::Offered;
    SurfaceMarking rejected = offered;
    rejected.id = 11;
    rejected.proposal = ProposalState::Rejected;
    set.markings.push_back(offered);
    set.markings.push_back(rejected);
    SurfaceMarking otherPiece;
    otherPiece.id = 12;
    otherPiece.meaning = MarkingMeaning::Start;
    otherPiece.shape = MarkingShape::Point;
    otherPiece.points.push_back(SurfacePoint{g.file.pieces->pieces[1].fileTriangle[0], {1.0 / 3, 1.0 / 3, 1.0 / 3}});
    set.markings.push_back(otherPiece);
    CHECK(effectiveMarkingsDigest(set, *g.file.pieces, 0) == base);
    set.markings[set.markings.size() - 3].proposal = ProposalState::Accepted;
    const uint64_t accepted = effectiveMarkingsDigest(set, *g.file.pieces, 0);
    CHECK(accepted != base);
    SceneSettings a = g.s, b = g.s;
    b.guidedMarkingsDigest = identityHex(accepted);
    CHECK(stageKey(a, Stage::Knit, 7) != stageKey(b, Stage::Knit, 7));
    b = a;
    b.guidedMarkings = g.dir + "/elsewhere.markings";
    CHECK(stageKey(a, Stage::Knit, 7) == stageKey(b, Stage::Knit, 7));
    SceneSettings surface = a, surfaceEdited = a;
    surface.knitRoute = surfaceEdited.knitRoute = "surface";
    surfaceEdited.guidedPiece = 1;
    surfaceEdited.guidedMarkingsDigest = identityHex(accepted);
    CHECK(stageKey(surface, Stage::Knit, 7) == stageKey(surfaceEdited, Stage::Knit, 7));
}

// A proposal's setup number and basis (file v4) are bookkeeping for the Mark step: neither enters the effective digest,
// so setting, changing or clearing them leaves the guided route's digest and the Knit key as they are. (Defect: the setup
// or basis hashed into effectiveMarkingsDigest, so regenerating proposals outdates a knit made from the same intent.)
TEST(guided_setup_fields_not_input) {
    GuidedFixture g = guidedFixture("setup_fields");
    const uint64_t base = effectiveMarkingsDigest(g.set, *g.file.pieces, 0);
    MarkingSet set = g.set;
    for (SurfaceMarking& m : set.markings) {
        m.setup = 3;
        m.basis = 0x1234;
    }
    CHECK(effectiveMarkingsDigest(set, *g.file.pieces, 0) == base);
    for (SurfaceMarking& m : set.markings) {
        m.setup = 9;
        m.basis = 0x5678;
    }
    CHECK(effectiveMarkingsDigest(set, *g.file.pieces, 0) == base);
    SceneSettings a = g.s, b = g.s;
    b.guidedMarkingsDigest = identityHex(effectiveMarkingsDigest(set, *g.file.pieces, 0));
    CHECK(stageKey(a, Stage::Knit, 7) == stageKey(b, Stage::Knit, 7));
    // control: an edit of a point does change it
    set.markings[0].points[0].bary[0] += 1e-9;
    set.markings[0].points[0].bary[1] -= 1e-9;
    CHECK(effectiveMarkingsDigest(set, *g.file.pieces, 0) != base);
}

// The stage boundaries of the guided route: a gauge, diameter, piece or markings edit makes the knit stale; a fibre-only
// edit leaves the knit and changes only the fibres' key. (Defect: a guided row missing from the stage inputs.)
TEST(guided_stage_boundaries) {
    GuidedFixture g = guidedFixture("boundaries");
    const SceneSettings s = g.s;
    const uint64_t knitKey = stageKey(s, Stage::Knit, 7);
    for (int k = 0; k < 5; ++k) {
        SceneSettings t = s;
        if (k == 0) t.wale += 0.5;
        if (k == 1) t.course += 0.4;
        if (k == 2) t.diameter += 0.1;
        if (k == 3) t.guidedPiece = 1;
        if (k == 4) t.guidedMarkingsDigest = "0123456789abcdef";
        CHECK(stageKey(t, Stage::Knit, 7) != knitKey);
    }
    SceneSettings fibre = s;
    fibre.fibreDensity += 10;
    CHECK(stageKey(fibre, Stage::Knit, 7) == knitKey);
    CHECK(fibreResultKey(fibre) != fibreResultKey(s));
    SceneSettings surface = s;
    surface.knitRoute = "surface";
    CHECK(stageKey(surface, Stage::Knit, 7) != knitKey);  // the route itself is an input
}

// Markings the stitch mesher cannot use (a finish ring and nothing to start from): the part is refused with the build's
// findings, not knitted, incomplete, and nothing is settled. (Defect: falling back to another route, or a knit anyway.)
TEST(guided_refused_construction) {
    if (!gpu::available(nullptr)) SKIP("no usable CUDA GPU");
    GuidedFixture g = guidedFixture("refused");
    ImportCacheScope cache(g.dir + "/imports");
    MarkingSet set = g.set;
    set.markings.erase(set.markings.begin());  // no start
    std::string err;
    CHECK(saveMarkingFile(g.markings, set, &err));
    SceneSettings s = g.s;
    s.guidedMarkingsDigest = identityHex(effectiveMarkingsDigest(set, *g.file.pieces, 0));
    Scene scene;
    KnitRunOptions o;
    CHECK(knitScene(*sceneEntryFor(s), s, o, scene, &err));
    CHECK(scene.parts.size() == 1);
    const ScenePart& part = scene.parts[0];
    CHECK(part.guided != nullptr && !part.knit);
    if (!part.guided) return;
    std::printf("  %s: %s\n", toString(part.guided->outcome), part.guided->reason.c_str());
    CHECK(part.guided->outcome == GuidedOutcome::Refused && !part.guided->reason.empty());
    CHECK(part.guided->calls == 0 && partIncomplete(part, ""));
}

// A running guided knit stops at the cancellation: the token fires while the settle runs and knitScene throws Cancelled
// within seconds (the builder then publishes nothing: builder_cancel_any_job_drops_pending). (Defect: a settle loop that
// never reaches a checkpoint.)
TEST(guided_cancel_publishes_nothing) {
    if (!gpu::available(nullptr)) SKIP("no usable CUDA GPU");
    GuidedFixture g = guidedFixture("cancel");
    ImportCacheScope cache(g.dir + "/imports");
    std::atomic<bool> settling{false}, cancelled{false}, finished{false};
    CancelToken token;
    std::thread worker([&] {
        ThreadCancelScope scope(&token);
        Scene scene;
        KnitRunOptions o;
        o.onStep = [&](const std::string& step, int, int) {
            if (step == "settling") settling = true;
        };
        std::string err;
        try {
            (void)knitScene(*sceneEntryFor(g.s), g.s, o, scene, &err);
        } catch (const Cancelled&) {
            cancelled = true;
        }
        finished = true;
    });
    const auto t0 = std::chrono::steady_clock::now();
    while (!settling && !finished && std::chrono::steady_clock::now() - t0 < std::chrono::seconds(120)) std::this_thread::sleep_for(std::chrono::milliseconds(5));
    CHECK(settling.load());
    const auto t1 = std::chrono::steady_clock::now();
    token.cancel();
    worker.join();
    const double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();
    std::printf("  cancelled %.1f s after the request\n", s);
    CHECK(cancelled.load() && s < 60);
}

// The guided knit's key follows the guided pipeline's version (scene.h kGuidedVersion, the stage input
// derived.guided_version): version 2 changed the layout, Draft and Final for the same markings and settings, so a knit kept
// under version 1 is not reused and a knit under the current version is. The tiling route's inputs name no guided version,
// so its keys are unchanged (io_pin_stage_keys pins them). A guided knit never reads the persistent knit cache
// (knit_job.cpp), so the in-memory knit by its key is the only one to reuse. (Defect: the bump missing, a reopened project
// showing the old behaviour; or the version line on every route, every tiled knit stale.)
TEST(guided_knit_key_follows_the_guided_version) {
    GuidedFixture g = guidedFixture("version_key");
    const std::string inputs = stageInputs(g.s, Stage::Knit);
    const std::string line = "derived.guided_version=" + std::to_string(kGuidedVersion) + "\n";
    CHECK(kGuidedVersion == 6 && inputs.find(line) != std::string::npos);
    CHECK(stageInputs(SceneSettings{}, Stage::Knit).find("derived.guided_version") == std::string::npos);
    SceneSettings surface = g.s;
    surface.knitRoute = "surface";
    CHECK(stageInputs(surface, Stage::Knit).find("derived.guided_version") == std::string::npos);
    // stageKey's hash, made here from the inputs (checked against stageKey), then the same inputs at version 1
    const uint64_t asset = 7, generator = 2;
    const auto keyOf = [&](const std::string& text) {
        uint64_t h = fnv1a(&generator, sizeof generator);
        h = fnv1a(&asset, sizeof asset, h);
        return fnv1a(text.data(), text.size(), h);
    };
    const uint64_t current = stageKey(g.s, Stage::Knit, asset);
    CHECK(keyOf(inputs) == current);
    std::string old = inputs;
    old.replace(old.find(line), line.size(), "derived.guided_version=1\n");
    const uint64_t version1 = keyOf(old);
    CHECK(version1 != current);
    auto kept = std::make_shared<Scene>();
    rememberKnit(version1, kept);
    CHECK(rememberedKnit(current) == nullptr);   // the old knit is not reused: the part is knitted again
    auto fresh = std::make_shared<Scene>();
    rememberKnit(current, fresh);
    CHECK(rememberedKnit(current) == fresh);     // a knit under the current version is
}

// Draft on the guided route: the generated yarn as curves with zero solver calls, the outcome Drafted (not "not settled"),
// a knit to show and export as a diagnostic, the target gate UNRESOLVED (it did not run) and the report saying what is
// shown and exported. (Defect at 4745fdc's parent: Draft ran the whole settle.)
TEST(guided_draft_without_a_solve) {
    if (!gpu::available(nullptr)) SKIP("no usable CUDA GPU");
    GuidedFixture g = guidedFixture("draft");
    ImportCacheScope cache(g.dir + "/imports");
    SceneSettings s = g.s;
    s.settle = false;
    Scene scene;
    KnitRunOptions o;
    bool settling = false;
    o.onStep = [&](const std::string& step, int, int) { settling = settling || step == "settling"; };
    std::string err;
    const auto t0 = std::chrono::steady_clock::now();
    CHECK(knitScene(*sceneEntryFor(s), s, o, scene, &err));
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    CHECK(scene.parts.size() == 1);
    if (scene.parts.empty() || !scene.parts[0].guided) return;
    const ScenePart& part = scene.parts[0];
    const GuidedPart& gp = *part.guided;
    std::printf("  %s: %s (%.1f s)\n", toString(gp.outcome), gp.reason.c_str(), seconds);
    CHECK(gp.outcome == GuidedOutcome::Drafted && gp.calls == 0 && !settling);
    CHECK(part.knit != nullptr && part.knit && part.knit->verdict != Status::Pass && part.knit->after.verdict == Status::Unresolved);
    CHECK(gp.yarnShown.find("not settled") != std::string::npos && gp.yarnExported.find("not settled") != std::string::npos);
    bool shownLine = false;
    for (const std::string& l : guidedPartLines(gp)) shownLine = shownLine || l.rfind("shown: ", 0) == 0;
    CHECK(shownLine && seconds < 120);
}

// The shirt trial's cuffs (2026-11-22, at 100 mm): the settle passed the target gate with the construction kept, but the
// stitch mesh's coverage check failed, and the part was headed "not settled". The yarn settled; only the mesh fails.
TEST(guided_settled_yarn_with_failing_mesh_checks_is_not_called_unsettled) {
    CHECK(outcomeAfterSettle(Status::Pass, Status::Pass, Status::Pass) == GuidedOutcome::Settled);
    CHECK(outcomeAfterSettle(Status::Pass, Status::Pass, Status::Fail) == GuidedOutcome::SettledMeshFails);
    CHECK(outcomeAfterSettle(Status::Pass, Status::Pass, Status::Unresolved) == GuidedOutcome::SettledMeshFails);
    CHECK(outcomeAfterSettle(Status::Fail, Status::Pass, Status::Pass) == GuidedOutcome::NotSettled);
    CHECK(outcomeAfterSettle(Status::Pass, Status::Unresolved, Status::Fail) == GuidedOutcome::NotSettled);
    CHECK(outcomeAfterSettle(Status::Unresolved, Status::Pass, Status::Pass) == GuidedOutcome::NotSettled);
    CHECK(std::string(toString(GuidedOutcome::SettledMeshFails)).find("not settled") == std::string::npos);
    // the report keeps the settle's own lines for such a part
    GuidedPart gp;
    gp.outcome = GuidedOutcome::SettledMeshFails;
    gp.reason = "the target gate passes";
    bool settleLine = false, headline = false;
    for (const std::string& l : guidedPartLines(gp)) {
        settleLine = settleLine || l.rfind("settle: D ", 0) == 0;
        headline = headline || l.rfind("guided stitch mesh (experimental): settled (", 0) == 0;
    }
    CHECK(settleLine && headline);
}

// The shirt trial at 2x (2026-11-22): with a long markings path the identity line filled its 600-character buffer and cut
// off the end of the solver identity, so a knit with KNIT_NCP_CONTACT_STORAGE=live did not say so in its report.
TEST(guided_inputs_keep_the_whole_solver_identity) {
    SceneSettings s;
    s.guidedMarkings = "D:/" + std::string(400, 'm') + "/shirt_pieces_v2-7f0c266ca025f74e-guided-516b75f66e475204.markings";
    s.guidedMarkingsDigest = "516b75f66e475204";
    GuidedInputs in;
    in.piece = 0;
    in.pieceName = "torso_lower piece 1";
    in.markings.preparation = "facets";
    const std::string text = guidedInputsText(s, in, 1.164);
    const SettleSettings settle = guidedSettleSettings(s);
    const std::string id = gpu::ncpSolverIdentity(settle.contactStorage, settle.curvatureRows);  // the settle's own options, whatever the environment says
    CHECK(settle.contactStorage == ContactStorage::LiveRows && settle.curvatureRows == CurvatureRows::SamplesAndSpanMaximum);
    CHECK(text.size() > 600);
    CHECK(text.find(s.guidedMarkings) != std::string::npos);
    CHECK(text.size() >= id.size() && text.compare(text.size() - id.size(), id.size(), id) == 0);
}

// The knitting component is part of what a guided knit is made from: its row enters the guided route's Knit stage inputs
// (component 1 and component 2 of a piece are different knits), never the default route's, and the identity line names it.
// (Defect: the component left out of the key, so the collar reuses the body's knit.)
TEST(guided_knit_component_enters_the_knit_key) {
    GuidedFixture g = guidedFixture("component_key");
    SceneSettings collar = g.s;
    collar.guidedKnitComponent = 1;
    CHECK(stageInputs(g.s, Stage::Knit).find("knit.guided_component=0\n") != std::string::npos);
    CHECK(stageInputs(collar, Stage::Knit).find("knit.guided_component=1\n") != std::string::npos);
    CHECK(stageKey(g.s, Stage::Knit, 7) != stageKey(collar, Stage::Knit, 7));
    SceneSettings surface = collar;
    surface.knitRoute = "surface";
    CHECK(stageInputs(surface, Stage::Knit).find("guided_component") == std::string::npos);
    SceneSettings surfaceBody = g.s;
    surfaceBody.knitRoute = "surface";
    CHECK(stageKey(surface, Stage::Knit, 7) == stageKey(surfaceBody, Stage::Knit, 7));   // the default route ignores the row
    GuidedInputs in;
    in.piece = 0;
    in.knitComponent = 1;
    CHECK(guidedInputsText(collar, in, 1.164).find("piece 0 knitting component 2 (") != std::string::npos);
}
