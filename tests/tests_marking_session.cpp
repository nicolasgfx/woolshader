// The Mark step's background work (src/jobs/marking_session.h), part of knit_tests_workflow: the quick lane never waits
// behind a long build, a superseded or cancelled build publishes nothing, and only the newest quick request of the open
// file is published; a generated stitch mesh is kept with its findings, a failed build replaces it with its reason, and
// its OBJ file holds the built polygons. Each case names the defect it catches.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <thread>
#include <vector>

#include "harness.h"
#include "io/stitch_mesh_obj.h"
#include "jobs/marking_session.h"
#include "model/surface_pieces.h"

using namespace knit;
namespace fs = std::filesystem;
using Clock = std::chrono::steady_clock;

namespace {

constexpr double kTwoPi = 6.283185307179586;

std::string sessionDir(const std::string& name) {
    const fs::path p = fs::temp_directory_path() / "knit_marking_session_tests" / name;
    std::error_code ec;
    fs::remove_all(p, ec);
    fs::create_directories(p, ec);
    return p.generic_string();
}

// An open tube round y (mm): `segments` round, rings - 1 a prime so that no row level lies on a ring of vertices.
void writeTube(const std::string& path, double radius, double length, int segments, int rings) {
    std::ofstream f(path, std::ios::binary);
    f << "o tube\n";
    for (int j = 0; j < rings; ++j)
        for (int i = 0; i < segments; ++i) {
            const double a = kTwoPi * i / segments;
            f << "v " << radius * std::cos(a) << " " << length * j / (rings - 1) << " " << -radius * std::sin(a) << "\n";
        }
    for (int j = 0; j + 1 < rings; ++j)
        for (int i = 0; i < segments; ++i) {
            const int a = j * segments + i + 1, b = j * segments + (i + 1) % segments + 1, c = a + segments, d = b + segments;
            f << "f " << a << " " << b << " " << d << "\nf " << a << " " << d << " " << c << "\n";
        }
}

bool waitFor(MarkingSession& s, double seconds, bool (*done)(const MarkingSession::Snapshot&)) {
    const Clock::time_point end = Clock::now() + std::chrono::milliseconds(static_cast<int>(seconds * 1000));
    while (Clock::now() < end) {
        if (done(s.snapshot())) return true;
        std::this_thread::yield();
    }
    return done(s.snapshot());
}

bool ready(const MarkingSession::Snapshot& s) { return s.phase == MarkingSession::Phase::Ready; }
bool proposed(const MarkingSession::Snapshot& s) { return s.phase == MarkingSession::Phase::Ready && s.proposals; }

// The first proposed setup as a set whose markings are accepted (the intent a build and a review take).
MarkingSet acceptedSetup(const MarkingSession::Snapshot& s) {
    MarkingSet set;
    set.digest = s.digest;
    for (const SurfaceMarking& m : s.proposals->setups.front().markings) {
        SurfaceMarking a = m;
        a.proposal = ProposalState::Accepted;
        set.markings.push_back(a);
        set.nextId = std::max(set.nextId, a.id + 1);
    }
    return set;
}

// A file opened, proposals made, setup 1 accepted: the session ready for builds and reviews of piece 0.
bool openWithSetup(MarkingSession& s, const std::string& obj, MarkingSet& set) {
    s.open(obj);
    if (!waitFor(s, 30, ready)) return false;
    s.propose(0, 1, 5.5, 3.4, "facets");
    if (!waitFor(s, 60, proposed)) return false;
    set = acceptedSetup(s.snapshot());
    return !set.markings.empty();
}

} // namespace

// A quick review requested while a long relaxed build runs is published before that build ends: the quick lane has its
// own worker and never waits behind the deep lane. (Defect: quick reviews queued on the build's worker.)
TEST(marking_quick_not_behind_build) {
    const std::string dir = sessionDir("not_behind");
    writeTube(dir + "/tube.obj", 30, 260, 160, 192);
    MarkingSession s;
    MarkingSet set;
    CHECK(openWithSetup(s, dir + "/tube.obj", set));
    const uint64_t first = s.review(set, 0, 1, 5.5, 3.4, {});   // warm the prepared piece
    CHECK(waitFor(s, 30, [](const MarkingSession::Snapshot& x) { return x.quick != nullptr; }));
    CHECK(s.snapshot().quick->serial == first);
    const Clock::time_point t0 = Clock::now();
    s.build(set, 0, 1, 5.5, 3.4, ReviewDepth::Relaxed);
    std::this_thread::sleep_for(std::chrono::milliseconds(50));
    const uint64_t serial = s.review(set, 0, 1, 5.6, 3.4, {});
    bool buildingWhenPublished = false;
    double quickMs = 0;
    while (std::chrono::duration<double>(Clock::now() - t0).count() < 120) {
        const MarkingSession::Snapshot x = s.snapshot();
        if (x.quick && x.quick->serial == serial) {
            buildingWhenPublished = x.phase == MarkingSession::Phase::Building;
            quickMs = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
            break;
        }
        std::this_thread::yield();
    }
    CHECK(waitFor(s, 300, [](const MarkingSession::Snapshot& x) { return x.phase == MarkingSession::Phase::Ready && x.built != nullptr; }));
    const double buildMs = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
    std::printf("  quick review published %.0f ms after the build started; the build took %.0f ms\n", quickMs, buildMs);
    CHECK(buildingWhenPublished && quickMs > 0 && quickMs < buildMs);
}

// A build superseded by a newer one stops (Cancelled) and publishes nothing: the shown build is the newer one's; cancelDeep
// ends a running build without a result. (Defect: a build that runs to its end before its result is dropped, holding the
// newer one back: cancelledTasks stays unchanged.) The times are printed, not checked: under a parallel test run they
// measure the machine's load as much as the cancellation.
TEST(marking_cancel_superseded) {
    const std::string dir = sessionDir("cancel");
    writeTube(dir + "/tube.obj", 30, 260, 160, 192);
    MarkingSession s;
    MarkingSet set;
    CHECK(openWithSetup(s, dir + "/tube.obj", set));
    // one build alone, for its time
    Clock::time_point t0 = Clock::now();
    s.build(set, 0, 1, 5.5, 3.4, ReviewDepth::Relaxed);
    CHECK(waitFor(s, 300, [](const MarkingSession::Snapshot& x) { return x.phase == MarkingSession::Phase::Ready && x.built != nullptr; }));
    const double alone = std::chrono::duration<double>(Clock::now() - t0).count();
    // a build superseded half-way by another gauge
    t0 = Clock::now();
    const int cancelledBefore = s.snapshot().cancelledTasks;
    s.build(set, 0, 1, 5.5, 3.4, ReviewDepth::Relaxed);
    std::this_thread::sleep_for(std::chrono::milliseconds(static_cast<int>(alone * 500)));
    s.build(set, 0, 1, 5.7, 3.4, ReviewDepth::Relaxed);
    const uint64_t newer = markBuildInputs(set, *s.snapshot().pieces, 0, 1, 5.7, 3.4, ReviewDepth::Relaxed);
    CHECK(waitFor(s, 300, [](const MarkingSession::Snapshot& x) { return x.phase == MarkingSession::Phase::Ready; }));
    const double both = std::chrono::duration<double>(Clock::now() - t0).count();
    std::printf("  one build %.2f s; superseded half-way and rebuilt: %.2f s in all\n", alone, both);
    CHECK(s.snapshot().built && s.snapshot().built->inputs == newer);
    CHECK(s.snapshot().cancelledTasks == cancelledBefore + 1);
    // cancelDeep during a build: Ready again at once, the earlier build stays shown, nothing new published
    const uint64_t shown = s.snapshot().built->inputs;
    s.build(set, 0, 1, 5.9, 3.4, ReviewDepth::Relaxed);
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    s.cancelDeep();
    CHECK(s.snapshot().phase == MarkingSession::Phase::Ready);
    std::this_thread::sleep_for(std::chrono::milliseconds(static_cast<int>(alone * 1500)));
    CHECK(s.snapshot().built && s.snapshot().built->inputs == shown);
}

// Only the newest quick request publishes: with the reviews held back, request A then B gives B's result and never A's
// (even though A finishes first); a review still held when another file is opened never publishes over it. (Defect: a
// rule that publishes anything newer than the last published result, so A lands while B is pending.)
TEST(marking_quick_latest_only) {
    const std::string dir = sessionDir("latest");
    writeTube(dir + "/tube.obj", 10.5, 40, 64, 30);
    writeTube(dir + "/other.obj", 12, 40, 64, 30);
    MarkingSession s;
    MarkingSet set;
    CHECK(openWithSetup(s, dir + "/tube.obj", set));
    const uint64_t warm = s.review(set, 0, 1, 5.5, 3.4, {});
    CHECK(waitFor(s, 30, [](const MarkingSession::Snapshot& x) { return x.quick != nullptr; }) && s.snapshot().quick->serial == warm);
    s.delayQuickForTest(300);
    const uint64_t a = s.review(set, 0, 1, 5.5, 3.4, {});
    std::this_thread::sleep_for(std::chrono::milliseconds(20));   // A is running (held), B pending
    MarkingSet edited = set;
    edited.markings.back().points.pop_back();
    const uint64_t b = s.review(edited, 0, 1, 5.5, 3.4, {});
    bool sawA = false;
    const Clock::time_point end = Clock::now() + std::chrono::seconds(5);
    while (Clock::now() < end) {
        const MarkingSession::Snapshot x = s.snapshot();
        if (x.quick && x.quick->serial == a) sawA = true;
        if (x.quick && x.quick->serial == b) break;
        std::this_thread::yield();
    }
    CHECK(!sawA && s.snapshot().quick && s.snapshot().quick->serial == b);
    CHECK(s.snapshot().quick->inputs == quickReviewInputs(edited, *s.snapshot().pieces, 0, 1, 5.5, 3.4, {}));
    // a held review, then another file: nothing of the first file is published
    const int generation = s.snapshot().fileGeneration;
    s.review(set, 0, 1, 5.5, 3.4, {});
    std::this_thread::sleep_for(std::chrono::milliseconds(20));
    s.open(dir + "/other.obj");
    std::this_thread::sleep_for(std::chrono::milliseconds(600));
    const MarkingSession::Snapshot x = s.snapshot();
    CHECK(x.fileGeneration == generation + 1 && x.quick == nullptr);
}

// --- the Stitch Mesh view's data: a generated mesh, its findings, its file

namespace {

// A quad strip of `faces` unit squares in the x-y plane, rows along y.
StitchMeshData quadStrip(int faces) {
    StitchMeshData d;
    for (int i = 0; i <= faces; ++i) {
        d.pos.push_back(Vec3{double(i), 0, 0});
        d.pos.push_back(Vec3{double(i), 1, 0});
    }
    d.faceStart.push_back(0);
    for (int f = 0; f < faces; ++f) {
        for (const int v : {2 * f, 2 * f + 2, 2 * f + 3, 2 * f + 1}) d.corner.push_back(v);
        d.faceStart.push_back(static_cast<int>(d.corner.size()));
        d.kind.push_back(f == 1 ? StitchKind::Increase : StitchKind::Regular);
        d.finding.push_back(0);
        d.deviationMm.push_back(0.1f);
        d.row.push_back(0);
    }
    return d;
}

// An open tube round the y axis through `centre`, in file units of `unitMm` mm (`radius` and `length` in mm).
void writeTubeAt(const std::string& path, const Vec3& centre, double unitMm, double radius, double length, int segments, int rings) {
    std::ofstream f(path, std::ios::binary);
    f.precision(17);
    f << "o tube\n";
    for (int j = 0; j < rings; ++j)
        for (int i = 0; i < segments; ++i) {
            const double a = kTwoPi * i / segments;
            const Vec3 mm{radius * std::cos(a), length * j / (rings - 1), -radius * std::sin(a)};
            const Vec3 p = centre + mm * (1.0 / unitMm);
            f << "v " << p.x << " " << p.y << " " << p.z << "\n";
        }
    for (int j = 0; j + 1 < rings; ++j)
        for (int i = 0; i < segments; ++i) {
            const int a = j * segments + i + 1, b = j * segments + (i + 1) % segments + 1, c = a + segments, d = b + segments;
            f << "f " << a << " " << b << " " << d << "\nf " << a << " " << d << " " << c << "\n";
        }
}

// The build of piece 0 that answers request `request` (null until it is published).
std::shared_ptr<const MarkedPieceView> buildOf(MarkingSession& s, uint64_t request, double seconds) {
    const Clock::time_point end = Clock::now() + std::chrono::milliseconds(static_cast<int>(seconds * 1000));
    while (Clock::now() < end) {
        const MarkingSession::Snapshot x = s.snapshot();
        const std::shared_ptr<const MarkedPieceView> b = findBuild(x, 0, 0);
        if (b && b->request == request && x.phase == MarkingSession::Phase::Ready) return b;
        std::this_thread::sleep_for(std::chrono::milliseconds(2));
    }
    return nullptr;
}

} // namespace

// The draw-and-save gate refuses only a mesh that cannot be drawn, each case with its own reason; fit findings never make a
// mesh undrawable, and the headline states generation and quality apart. (Defect: a mesh with findings hidden, or an
// unsafe one drawn.)
TEST(stitch_mesh_use_gate) {
    CHECK(!stitchMeshUse(StitchMeshData{}).drawable && stitchMeshUse(StitchMeshData{}).reason == "the build made no stitch faces");
    StitchMeshData nan = quadStrip(3);
    nan.pos[4].y = std::nan("");
    CHECK(!stitchMeshUse(nan).drawable && stitchMeshUse(nan).reason == "vertex 4 has a non-finite coordinate");
    StitchMeshData index = quadStrip(3);
    index.corner[5] = 99;
    CHECK(!stitchMeshUse(index).drawable && stitchMeshUse(index).reason == "face 1 names vertex 99 of 8");
    StitchMeshData twoCorners = quadStrip(3);
    twoCorners.corner.erase(twoCorners.corner.begin() + 2, twoCorners.corner.begin() + 4);
    for (size_t f = 1; f < twoCorners.faceStart.size(); ++f) twoCorners.faceStart[f] -= 2;
    CHECK(!stitchMeshUse(twoCorners).drawable && stitchMeshUse(twoCorners).reason == "face 0 has 2 corners");
    // findings everywhere, still drawable
    StitchMeshData found = quadStrip(3);
    found.finding.assign(3, kFaceOverFit | kFaceInverted);
    const StitchMeshUse use = stitchMeshUse(found);
    CHECK(use.drawable && use.reason.empty());
    StitchMeshQuality q;
    q.faces = 563;
    q.rows = 25;
    q.fitTargetMm = 0.55;
    q.maxDeviationMm = 2.758;
    q.facesOverFit = 117;
    q.connectivityPass = true;
    q.uncoveredLimit = 0.005;
    q.areaLow = 0.98;
    q.areaHigh = 1.02;
    q.areaRatio = 1;
    q.boundaryLimitMm = 0.85;
    const std::string headline = stitchMeshHeadline(use, q);
    std::printf("  %s\n", headline.c_str());
    CHECK(headline == "Mesh generated \xC2\xB7 563 faces \xC2\xB7 25 rows \xC2\xB7 117 faces exceed the 0.55 mm fit target");
    q.facesOverFit = 0;
    q.clean = true;
    CHECK(stitchMeshHeadline(use, q) == "Mesh generated \xC2\xB7 563 faces \xC2\xB7 25 rows \xC2\xB7 passes its checks");
    q.relaxed = true;
    q.relaxRan = true;
    q.relaxStop = "stalled";   // every relaxation on this branch stalls: a note, never a reason to read a passing mesh as failing
    CHECK(stitchMeshHeadline(use, q) == "Mesh generated \xC2\xB7 563 faces \xC2\xB7 25 rows \xC2\xB7 passes its checks");
    CHECK(stitchMeshRelaxationNote(q) == "The relaxation stopped without converging (stalled)");
    q.inverted = 1;
    q.facesOverFit = 117;
    const std::vector<std::string> findings = stitchMeshFindings(q);
    CHECK(findings.size() == 2 && findings[0] == "1 face is wound against the piece (inverted)");   // structure first
    CHECK(stitchMeshHeadline(StitchMeshUse{false, "face 1 names vertex 99 of 8"}, q) == "No stitch mesh: face 1 names vertex 99 of 8");
}

// A build that fails replaces its piece's entry with a view that has no mesh and says why: an error (the backend returns
// false), and a build that returns without stitch faces (the row-start path removed) with the backend's own finding. The
// earlier mesh is never kept for the piece. (Defect: the previous build left on screen after a failed rebuild.)
TEST(stitch_mesh_failed_rebuild_replaces) {
    const std::string dir = sessionDir("failed_rebuild");
    writeTube(dir + "/tube.obj", 10.5, 40, 64, 30);
    MarkingSession s;
    MarkingSet set;
    CHECK(openWithSetup(s, dir + "/tube.obj", set));
    const std::shared_ptr<const MarkedPieceView> good = buildOf(s, s.build(set, 0, 1, 5.5, 3.4, ReviewDepth::Built), 120);
    CHECK(good && good->use.drawable && good->stitch.faces() > 0 && good->scene->meshes.size() == 3);
    std::printf("  good build: %s\n", stitchMeshHeadline(good->use, good->quality).c_str());
    // an unusable request (a stitch width of 0): buildFromMarkings returns false
    const std::shared_ptr<const MarkedPieceView> failed = buildOf(s, s.build(set, 0, 1, 0.0, 3.4, ReviewDepth::Built), 60);
    CHECK(failed && !failed->use.drawable && failed->use.reason.rfind("the build failed: ", 0) == 0 && failed->stitch.faces() == 0);
    CHECK(failed->scene && failed->scene->meshes.size() == 1 && failed->stitchEdges.empty());
    CHECK(!s.snapshot().error.empty() && s.snapshot().builds.size() == 1 && findBuild(s.snapshot(), 0, 0) == failed);
    std::printf("  failed build: %s\n", stitchMeshHeadline(failed->use, failed->quality).c_str());
    // a good build again, then one without the row-start path: it returns, without faces, and names the missing marking
    CHECK(buildOf(s, s.build(set, 0, 1, 5.5, 3.4, ReviewDepth::Built), 120) != nullptr);
    MarkingSet noPath = set;
    noPath.markings.erase(std::remove_if(noPath.markings.begin(), noPath.markings.end(), [](const SurfaceMarking& m) { return m.meaning == MarkingMeaning::RowStartPath; }),
                          noPath.markings.end());
    CHECK(noPath.markings.size() + 1 == set.markings.size());
    const std::shared_ptr<const MarkedPieceView> empty = buildOf(s, s.build(noPath, 0, 1, 5.5, 3.4, ReviewDepth::Built), 60);
    CHECK(empty && !empty->review.valid && !empty->use.drawable && empty->scene->meshes.size() == 1 && s.snapshot().error.empty());
    bool named = false;
    for (const MarkingFinding& f : empty->review.findings) named = named || f.message == empty->use.reason;
    CHECK(named);
    std::printf("  build without a row-start path: %s\n", stitchMeshHeadline(empty->use, empty->quality).c_str());
}

// The OBJ file holds exactly the built polygons: the file read back has the same vertices (bit for bit), the same corner
// lists, rows and kinds; the vertices lie on the imported tube in the file's own units and place (off the origin, 68.1 mm
// per unit), and the drawn faces are the same corners lifted 0.15 mm. (Defect: an export scaled, moved, re-ordered,
// triangulated or lifted like the display.)
TEST(stitch_mesh_obj_matches_build) {
    const std::string dir = sessionDir("obj_matches");
    const Vec3 centre{37.5, -12.25, 104};
    const double unitMm = 68.1, radius = 10.5, length = 40;
    writeTubeAt(dir + "/tube.obj", centre, unitMm, radius, length, 64, 30);
    MarkingSession s;
    s.open(dir + "/tube.obj");
    CHECK(waitFor(s, 30, ready));
    s.propose(0, unitMm, 5.5, 3.4, "facets");
    CHECK(waitFor(s, 60, proposed));
    const MarkingSet set = acceptedSetup(s.snapshot());
    const std::shared_ptr<const MarkedPieceView> b = buildOf(s, s.build(set, 0, unitMm, 5.5, 3.4, ReviewDepth::Relaxed), 180);
    CHECK(b && b->use.drawable);
    if (!b || !b->use.drawable) return;
    const StitchMeshData& d = b->stitch;
    std::printf("  %s\n", stitchMeshHeadline(b->use, b->quality).c_str());
    // on the tube, in file units and place
    double worstRadial = 0, lowY = 1e300, highY = -1e300;
    for (const Vec3& p : d.pos) {
        const double r = std::hypot(p.x - centre.x, p.z - centre.z) * unitMm;
        worstRadial = std::max(worstRadial, std::fabs(r - radius));
        lowY = std::min(lowY, (p.y - centre.y) * unitMm);
        highY = std::max(highY, (p.y - centre.y) * unitMm);
    }
    std::printf("  vertices %zu: radius off by at most %.4f mm, height %.3f .. %.3f mm\n", d.pos.size(), worstRadial, lowY, highY);
    CHECK(worstRadial < 0.1 && lowY > -0.1 && highY < length + 0.1);
    // the per-face findings add up to the check's counts, and the largest per-face deviation is the check's maximum
    int overFit = 0, inverted = 0, degenerate = 0;
    double largest = 0;
    for (int f = 0; f < d.faces(); ++f) {
        overFit += (d.finding[f] & kFaceOverFit) != 0;
        inverted += (d.finding[f] & kFaceInverted) != 0;
        degenerate += (d.finding[f] & kFaceDegenerate) != 0;
        largest = std::max(largest, double(d.deviationMm[f]));
    }
    const StitchMeshQuality& q = b->quality;
    CHECK(overFit == q.facesOverFit && inverted == q.inverted && degenerate == q.degenerate && std::fabs(largest - q.maxDeviationMm) < 1e-5);
    int rowFaces = 0;
    for (int n : q.facesPerRow) rowFaces += n;
    CHECK(rowFaces == q.faces && q.rows > 0);
    // the file read back
    const std::string path = dir + "/stitch_mesh.obj";
    std::string err;
    CHECK(writeStitchMeshObj(path, d, {"test"}, "tube", &err));
    StitchMeshData back;
    CHECK(readPolygonObj(path, back, &err));
    CHECK(back.pos == d.pos && back.faceStart == d.faceStart && back.corner == d.corner && back.kind == d.kind && back.row == d.row);
    // the display: the first face's first triangle is its corners 0, 1, 2 lifted along the face normal by 0.15 mm
    const DisplayMesh& byKind = b->scene->meshes[1];
    const Vec3 shown{byKind.pos[0], byKind.pos[1], byKind.pos[2]};
    const double lift = norm(shown - d.pos[d.corner[0]]) * unitMm;
    std::printf("  display lift %.4f mm\n", lift);
    CHECK(std::fabs(lift - 0.15) < 0.01);
}
