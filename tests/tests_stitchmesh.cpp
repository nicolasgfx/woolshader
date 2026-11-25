// knit_tests_stitchmesh: the stockinette stitch-mesh experiment (src/stitchmesh).
// Checkpoint 1: the reference surface, the time field, the rings and the orientation field. Checkpoint 2: the
// construction fixtures and negative controls, the virtual shift, Method A's typed mesh and coverage. Checkpoint 3:
// the relaxation (forces, rest state, perturbation, the virtual helix, the stall rule). Checkpoint 4: the helix boundary
// contract. Run one case:
// knit_tests_stitchmesh <name>; knit_tests_stitchmesh --list.
#include <algorithm>
#include <filesystem>
#include <cmath>
#include <map>
#include <set>

#include "harness.h"
#include "model/hash.h"
#include "model/mesh_check.h"
#include "model/mesh_query.h"
#include "stitchmesh/border_placement.h"
#include "stitchmesh/construction_fixtures.h"
#include "stitchmesh/coverage.h"
#include "stitchmesh/guided_rows.h"
#include "stitchmesh/reference_markings.h"
#include "stitchmesh/orientation_field.h"
#include "stitchmesh/relax_mesh.h"
#include "stitchmesh/revolution_surface.h"
#include "stitchmesh/row_rings.h"
#include "stitchmesh/time_field.h"
#include "stitchmesh/typed_mesh.h"
#include "stitchmesh/virtual_shift.h"

using namespace knit;
using namespace knit::refmark;

namespace {

RevolutionMesh tube(double jitter = 0) {
    RevolutionMeshParams p;
    p.jitter = jitter;
    p.jitterSeed = jitter > 0 ? 1 : 0;
    RevolutionMesh rm;
    std::string err;
    CHECK(makeRevolutionMesh(RevolutionProfile{}, p, rm, &err));
    return rm;
}

std::vector<double> timeFromBottom(const RevolutionMesh& rm) {
    std::vector<double> t;
    TimeFieldReport tr;
    std::string err;
    CHECK(geodesicDistance(rm.mesh, rm.bottomLoop, t, tr, &err));
    return t;
}

} // namespace

// The simple reference surface: one manifold, consistently oriented component with two boundary loops, every vertex
// on the exact surface, the frozen tessellation's counts.
TEST(ssm_revolution_mesh_valid) {
    const RevolutionMesh rm = tube();
    CHECK(rm.rings == 152);
    CHECK(rm.mesh.tri.size() == (size_t) 2 * 192 * 151);
    const MeshCheck mc = checkMesh(rm.mesh, SelfIntersectionTest::Skip);
    CHECK(mc.components.size() == 1);
    CHECK(mc.components[0].nonManifoldEdges == 0 && mc.components[0].consistent);
    CHECK(mc.components[0].boundaryEdges == 2 * 192);
    CHECK(mc.components[0].euler == 0);   // an annulus
    const RevolutionProfile prof;
    for (const Vec3& q : rm.mesh.pos) CHECK_NEAR(std::hypot(q.x, q.y), prof.radius(q.z), 1e-9);
    CHECK_NEAR(rm.meridianLength, 151.0152, 1e-4);
    // outward normals (counter-clockwise triangles seen from outside)
    for (size_t f = 0; f < rm.mesh.tri.size(); f += 97) {
        const Vec3 c = (rm.mesh.pos[rm.mesh.tri[f][0]] + rm.mesh.pos[rm.mesh.tri[f][1]] + rm.mesh.pos[rm.mesh.tri[f][2]]) / 3.0;
        CHECK(dot(rm.mesh.faceNormal((int) f), Vec3{c.x, c.y, 0}) > 0);
    }
}

// Variant J keeps the boundary rings, moves every other vertex by at most the jitter, and stays on the exact surface.
TEST(ssm_revolution_variant_j) {
    const RevolutionMesh base = tube(), j = tube(0.25);
    CHECK(j.mesh.pos.size() == base.mesh.pos.size());
    const RevolutionProfile prof;
    int moved = 0;
    for (size_t v = 0; v < j.mesh.pos.size(); ++v) {
        CHECK_NEAR(std::hypot(j.mesh.pos[v].x, j.mesh.pos[v].y), prof.radius(j.mesh.pos[v].z), 1e-9);
        CHECK(std::fabs(j.arc[v] - base.arc[v]) <= 0.25 + 1e-12);
        if (v < 192 || v >= j.mesh.pos.size() - 192) CHECK(j.mesh.pos[v] == base.mesh.pos[v]);
        else moved += j.mesh.pos[v] == base.mesh.pos[v] ? 0 : 1;
    }
    CHECK(moved == (int) j.mesh.pos.size() - 2 * 192);
}

// The planar-front update is exact for a straight front on a plane of acute triangles: an equilateral lattice and the
// distance from its first column (a segment along the lattice vector (0.5, 0.866), which edge updates alone would
// overestimate). Exact only inside the scheme's domain of dependence of the segment: each update averages two corners,
// one of them a step towards the segment's end, so the circular front beyond the end (a point source, where a
// first-order update is not exact) leaks in and decays by one half per step along i + j = const. Vertices with
// i + j <= ny never see it. On the tube the source is a closed boundary loop with no end, so this does not arise.
TEST(ssm_time_field_planar_front_exact) {
    TriMesh m;
    const int nx = 30, ny = 20;
    const double s3 = std::sqrt(3.0) / 2.0;
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i) m.pos.push_back(Vec3{1.0 * i + 0.5 * j, s3 * j, 0});
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            const int a = j * (nx + 1) + i, b = a + 1, c = a + nx + 1, d = c + 1;
            m.tri.push_back({a, b, c});
            m.tri.push_back({b, d, c});
        }
    std::vector<int> sources;
    for (int j = 0; j <= ny; ++j) sources.push_back(j * (nx + 1));
    std::vector<double> t;
    TimeFieldReport tr;
    std::string err;
    CHECK(geodesicDistance(m, sources, t, tr, &err));
    const Vec3 n{s3, -0.5, 0};   // unit normal of the first column's line
    int checked = 0;
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i + j <= ny && i <= nx; ++i) {
            const size_t v = (size_t) (j * (nx + 1) + i);
            CHECK_NEAR(t[v], dot(m.pos[v], n), 1e-9);
            ++checked;
        }
    CHECK(checked == (ny + 1) * (ny + 2) / 2);
    CHECK(tr.triangleUpdates > 0);
}

// On the tube the field from the bottom boundary is the meridian arc length. Bound: 0.05 mm, so consecutive rings at
// the rectangular gauge (h' = 3.43 mm) keep their spacing within the frozen +-3 % (2 x 0.05 < 0.03 x 3.43).
TEST(ssm_time_field_tube_matches_arc_length) {
    const RevolutionMesh rm = tube();
    const std::vector<double> t = timeFromBottom(rm);
    double worst = 0;
    for (size_t v = 0; v < t.size(); ++v) worst = std::max(worst, std::fabs(t[v] - rm.arc[v]));
    CHECK(worst <= 0.05);
}

// Every interior level on the tube is one closed contour, crossed once by the cut half-plane, oriented in +phi (growth
// on the left seen from outside), and as long as the exact parallel within the chord deficit (0.5 %).
TEST(ssm_rings_one_closed_contour_per_level) {
    const RevolutionMesh rm = tube();
    const std::vector<double> t = timeFromBottom(rm);
    const RevolutionProfile prof;
    MeridianTable table;
    table.build(prof, 600001);
    for (double level = 3.4; level < rm.meridianLength - 1; level += 3.4) {
        const RingLevel lv = extractLevel(rm.mesh, t, level);
        CHECK(lv.contours.size() == 1);
        const Contour& c = lv.contours.front();
        CHECK(c.closed);
        int crossings = 0;
        CHECK(cutCrossing(c, Vec3{0, 0, 0}, Vec3{0, 1, 0}, Vec3{1, 0, 0}, &crossings) >= 0);
        CHECK(crossings == 1);
        const double exact = 2 * kPi * prof.radius(table.heightAt(level));
        CHECK(std::fabs(c.length / exact - 1) < 0.005);
        // +phi: the cross product of consecutive radial directions points up
        const Vec3 p0 = c.points[0].p, p1 = c.points[1].p;
        CHECK(p0.x * p1.y - p0.y * p1.x > 0);
    }
}

// A hole in the tube: the levels through it are not one closed contour, and that is reported (open contours ending on
// the hole's boundary), while levels away from it stay single and closed.
TEST(ssm_rings_report_open_contours_at_a_hole) {
    RevolutionMesh rm = tube();
    TriMesh holed;
    holed.pos = rm.mesh.pos;
    const Vec3 centre{-20 - 2.5, 0, 90};   // on the taper, opposite the cut
    for (const auto& tr : rm.mesh.tri) {
        const Vec3 c = (rm.mesh.pos[tr[0]] + rm.mesh.pos[tr[1]] + rm.mesh.pos[tr[2]]) / 3.0;
        if (norm(c - centre) > 6) holed.tri.push_back(tr);
    }
    std::vector<double> t;
    TimeFieldReport tr;
    std::string err;
    CHECK(geodesicDistance(holed, rm.bottomLoop, t, tr, &err));
    int eventLevels = 0, cleanLevels = 0;
    for (double level = 3.4; level < rm.meridianLength - 1; level += 3.4) {
        const RingLevel lv = extractLevel(holed, t, level);
        const bool clean = lv.contours.size() == 1 && lv.contours.front().closed;
        if (clean) ++cleanLevels;
        else {
            ++eventLevels;
            bool open = false;
            for (const Contour& c : lv.contours) open = open || !c.closed;
            CHECK(open);
        }
    }
    CHECK(eventLevels >= 2);
    CHECK(cleanLevels >= 30);
}

// Method B with the axis guide: RIM's fixed orientation reproduces the guide (no sweeps), so both methods see one
// row axis; RIM subdivides nothing at the coarse edge 2 w.
TEST(ssm_orientation_guided_equals_guide) {
    const RevolutionMesh rm = tube();
    OrientationFieldRequest r;
    r.spacing = 11;
    r.guide = ringAxisGuide(rm.mesh, Vec3{0, 0, 1});
    OrientationField of;
    std::string err;
    CHECK(computeOrientationField(rm.mesh, r, of, &err));
    CHECK(of.finestVertices == (int) rm.mesh.pos.size());
    double worst = 0;
    for (size_t v = 0; v < of.q.size(); ++v) worst = std::max(worst, directionDeviation(of.q[v], r.guide[v], rm.mesh.vnormal[v], 180));
    CHECK(worst < 1e-6);
}

namespace {

TypedCheckOptions fixtureOptions() {
    TypedCheckOptions o;
    o.boundaryTolerance = 1e-6;
    o.bindingTolerance = 0.01;
    o.faceDeviation = 1e9;
    o.degenerateArea = 0.01;
    return o;
}

const FixtureCase& findCase(const std::vector<FixtureCase>& all, const std::string& name) {
    for (const FixtureCase& c : all)
        if (c.name == name) return c;
    CHECK(false);
    return all.front();
}

TypedMeshCheck checkCase(const std::vector<FixtureCase>& all, const std::string& name) {
    const FixtureCase& c = findCase(all, name);
    return checkTypedMesh(c.mesh, c.reference, c.declared, fixtureOptions());
}

VirtualShiftCheck shiftOf(const FixtureCase& c) {
    VirtualShift v;
    std::string err;
    CHECK(buildVirtualShift(c.mesh, c.reading, v, &err));
    return validateVirtualShift(c.mesh, v);
}

} // namespace

// Every fixture and negative control comes out as written in construction_fixtures.cpp: the positive ones pass, each
// negative control is refused (FAIL) or listed (UNRESOLVED) by the independent check or the virtual-shift validator.
TEST(ssm_fixtures_as_expected) {
    const std::vector<FixtureCase> all = constructionFixtureCases();
    CHECK(all.size() == 16);
    for (const FixtureCase& c : all) {
        const TypedMeshCheck k = checkTypedMesh(c.mesh, c.reference, c.declared, fixtureOptions());
        CHECK(k.connectivity == c.expectConnectivity);
        if (c.runShift) CHECK(shiftOf(c).status == c.expectShift);
        if (!c.negative) CHECK(k.geometry == Status::Pass);
    }
}

// The open tube 8 x 3: 24 K, every ring-0 edge CastOn and every ring-3 edge BindOff, 16 interlocks, three closed strands
// stored; the virtual helix traced on the re-read polygons is ONE open strand through all 24 faces, from the start step
// (strip 0, slot 0) to the end step (strip 2, last slot), with the same 16 interlocks and boundary edges.
TEST(ssm_fixture_open_tube_virtual_helix) {
    const std::vector<FixtureCase> all = constructionFixtureCases();
    const FixtureCase& c = findCase(all, "k-open-tube-8x3");
    const TypedMeshCheck k = checkTypedMesh(c.mesh, c.reference, c.declared, fixtureOptions());
    CHECK(k.faces == 24 && k.regular == 24 && k.interlocks == 16 && k.castOnEdges == 8 && k.bindOffEdges == 8);
    CHECK(k.strands.strands.size() == 3 && k.closedStrands == 3 && k.strands.yarnEnds == 0 && k.waleAgree == 24);
    const VirtualShiftCheck v = shiftOf(c);
    CHECK(v.status == Status::Pass && v.strands == 1 && v.longestStrand == 24 && v.startSteps == 1 && v.endSteps == 1);
    CHECK(c.mesh.strip[v.startFace] == 0 && c.mesh.slot[v.startFace] == 0 && c.mesh.strip[v.endFace] == 2 && c.mesh.slot[v.endFace] == 7);
    CHECK(v.interlocksStored == 16 && v.interlocksVirtual == 16 && v.interlocksMissing == 0 && v.interlocksExtra == 0);
    // the helix visits strip 0 slots 0..7, then strip 1 slots 0..7, then strip 2
    for (int i = 0; i < 24; ++i) CHECK(c.mesh.strip[v.strandOrder[i]] == i / 8 && c.mesh.slot[v.strandOrder[i]] == i % 8);
}

// kp 8 -> 9 and d12k 9 -> 8: one shaping face each, the count invariant per strip, the interfaces, and the virtual helix
// through all 25 faces.
TEST(ssm_fixture_kp_and_d12k) {
    const std::vector<FixtureCase> all = constructionFixtureCases();
    for (const char* name : {"kp-ring-8-9", "d12k-ring-9-8"}) {
        const FixtureCase& c = findCase(all, name);
        const TypedMeshCheck k = checkTypedMesh(c.mesh, c.reference, c.declared, fixtureOptions());
        const bool inc = std::string(name) == "kp-ring-8-9";
        CHECK(k.connectivity == Status::Pass && k.faces == 25 && k.increase == (inc ? 1 : 0) && k.decrease == (inc ? 0 : 1));
        CHECK(k.strips.size() == 3 && k.stripsBroken == 0);
        CHECK(k.strips[1].topEdges - k.strips[1].bottomEdges == (inc ? 1 : -1));
        CHECK(k.interlocks == 17 && k.bothTop == 0 && k.bothBottom == 0 && k.waleDisagree == 0);
        const VirtualShiftCheck v = shiftOf(c);
        CHECK(v.status == Status::Pass && v.longestStrand == 25);
    }
}

// The special-quad neighbourhood (dbf5f41 coarse fixture 2 through the ported subdivision): two special quads, and the
// yarn passes THROUGH them: five strands, all closed, no yarn end. A special quad ends a stitch row, not a strand: the
// strand of the short row (14 faces) enters each special quad through one wale edge and leaves through the other.
TEST(ssm_fixture_special_quad_strands_stay_closed) {
    const std::vector<FixtureCase> all = constructionFixtureCases();
    const FixtureCase& c = findCase(all, "special-quad-neighbourhood");
    const TypedMeshCheck k = checkTypedMesh(c.mesh, c.reference, c.declared, fixtureOptions());
    CHECK(k.connectivity == Status::Pass && k.regular == 44 && k.special == 2 && k.unsupported == 0);
    CHECK(k.strands.strands.size() == 5 && k.closedStrands == 5 && k.openStrands == 0 && k.strands.yarnEnds == 0);
    CHECK(k.interlocks == 38 && k.castOnEdges == 8 && k.bindOffEdges == 8);
    int strandOfSpecial = -1;
    for (int f = 0; f < c.mesh.mesh.faces(); ++f)
        if (c.mesh.kind[f] == StitchFaceKind::SpecialQuad) {
            if (strandOfSpecial < 0) strandOfSpecial = k.strands.strandOfFace[f];
            CHECK(k.strands.strandOfFace[f] == strandOfSpecial);
        }
    CHECK(strandOfSpecial >= 0 && k.strands.strands[strandOfSpecial].size() == 14 && k.strands.closed[strandOfSpecial]);
}

// The negative controls are caught by the specific count, not only by the verdict.
TEST(ssm_negative_controls_name_their_defect) {
    const std::vector<FixtureCase> all = constructionFixtureCases();
    CHECK(checkCase(all, "neg-cast-on-on-interior-ring").boundaryRoleOnInteriorEdge == 1);
    CHECK(checkCase(all, "neg-bind-off-with-face-above").boundaryRoleOnInteriorEdge == 1);
    CHECK(checkCase(all, "neg-first-strip-dec").unsupported == 1);
    CHECK(checkCase(all, "neg-shaping-beside-path").unsupported == 1);
    CHECK(checkCase(all, "neg-kp-roles-out-of-order").kindMismatches == 1);
    const TypedMeshCheck flipped = checkCase(all, "neg-face-flipped");
    CHECK(flipped.bothTop == 1 && flipped.bothBottom == 1 && flipped.waleDisagree == 2);
    CHECK(checkCase(all, "neg-missing-face").stripsBroken >= 1);
    CHECK(checkCase(all, "neg-special-quad-wale-edges-apart").kindMismatches == 1);
    const VirtualShiftCheck none = shiftOf(findCase(all, "neg-shift-not-applied"));
    CHECK(none.strands == 3 && none.closedStrands == 3 && !none.oneStrandThroughAll);
    const VirtualShiftCheck up = shiftOf(findCase(all, "neg-shift-read-upward"));
    CHECK(!up.startsInFirstStrip);
    CHECK(shiftOf(findCase(all, "neg-shift-path-skips-ring")).stripsWithoutPathEdge == 1);
    CHECK(shiftOf(findCase(all, "neg-shift-path-two-edges-in-ring")).stripsWithSeveralPathEdges == 1);
}

// One vertex of the tube pushed 6 mm outward (a spike, as a photogrammetry or remeshed import has): its tip lies
// farther from the cast-on than its neighbours above it, so the levels just below the tip's value close into a tiny
// contour round it beside the row. Too small for a round of 3 stitches, it is knitted over: no contour event, the bump recorded
// with its ring, the rows as without the spike. Controls: the same spike with a stitch so narrow that the contour would
// hold a round is a contour event again; the helper alone keeps a closed contour longer than the limit and an open one.
TEST(ssm_guided_rows_knit_over_a_bump_too_small_for_a_round) {
    const RevolutionMesh plain = tube();
    RevolutionMesh rm = tube();
    const std::vector<double> t0 = timeFromBottom(plain);
    int tip = -1;
    for (int v = 0; v < static_cast<int>(rm.mesh.pos.size()); ++v)
        if (tip < 0 || std::fabs(t0[v] - 0.5 * rm.meridianLength) < std::fabs(t0[tip] - 0.5 * rm.meridianLength)) tip = v;
    rm.mesh.computeNormals();
    rm.mesh.pos[tip] += rm.mesh.vnormal[tip] * 6.0;
    rm.mesh.computeNormals();
    const std::vector<double> t = timeFromBottom(rm);
    GuidedRowsRequest r;
    r.w = 5.5;
    r.h = 3.4;
    GuidedRowsProduct p, q;
    std::string err;
    CHECK(buildGuidedRows(plain.mesh, t0, plain.bottomLoop, plain.topLoop, r, q, &err));
    CHECK(buildGuidedRows(rm.mesh, t, rm.bottomLoop, rm.topLoop, r, p, &err));
    CHECK(p.contourEvents == 0 && p.status != Status::Fail && p.rows == q.rows);
    CHECK(!p.knittedOver.empty());
    for (const KnittedOverBump& b : p.knittedOver) CHECK(b.length < kKnitOverStitches * r.w && b.ring > 0 && b.ring < p.rows && norm(b.at - rm.mesh.pos[tip]) < 6.5);
    CHECK(std::abs(p.mesh.mesh.faces() - q.mesh.mesh.faces()) <= 2);
    // control: a stitch narrower than the bump's contour knits nothing over, and the spike is a contour event
    GuidedRowsRequest narrow = r;
    narrow.w = 0.3 * p.knittedOver.front().length;
    GuidedRowsProduct e;
    CHECK(buildGuidedRows(rm.mesh, t, rm.bottomLoop, rm.topLoop, narrow, e, &err));
    CHECK(e.contourEvents >= 1 && e.status == Status::Fail && e.knittedOver.empty());
    // control: the helper keeps a short open contour and a closed one at the limit or longer
    RingLevel level;
    level.contours.resize(4);
    level.contours[0].closed = true, level.contours[0].length = 40;
    level.contours[1].closed = true, level.contours[1].length = 2;
    level.contours[2].closed = false, level.contours[2].length = 1;
    level.contours[3].closed = true, level.contours[3].length = 5.5;
    const std::vector<Contour> dropped = dropContoursShorterThan(level, 5.5);
    CHECK(dropped.size() == 1 && dropped[0].length == 2 && level.contours.size() == 3 && level.contours[0].length == 40);
}

// Method A on the tube (case D, gauge R): every ring one closed contour, counts 34 -> 23 with 11 DEC in 11 strips, all
// clear of the path, connectivity and coverage PASS, the virtual helix through every face; removing one strip's faces
// makes coverage FAIL (the coverage check can fail).
TEST(ssm_guided_rows_tube_case_d) {
    const RevolutionMesh rm = tube();
    const std::vector<double> t = timeFromBottom(rm);
    GuidedRowsRequest r;
    r.w = 5.5;
    r.h = 3.4;
    GuidedRowsProduct p;
    std::string err;
    CHECK(buildGuidedRows(rm.mesh, t, rm.bottomLoop, rm.topLoop, r, p, &err));
    CHECK(p.status == Status::Pass && p.rows == 44 && p.contourEvents == 0);
    CHECK(p.rings.front().n == 34 && p.rings.back().n == 23 && p.shaping.size() == 11);
    for (const ShapingPlacement& s : p.shaping) CHECK(s.placed && s.slot >= 1 && s.slot <= s.slotsInStrip - 2 && s.kind == StitchFaceKind::Decrease);
    std::vector<DeclaredBoundary> declared(2);
    declared[0].role = EdgeBoundary::CastOn;
    declared[1].role = EdgeBoundary::BindOff;
    for (int v : rm.bottomLoop) declared[0].loop.push_back(rm.mesh.pos[v]);
    for (int v : rm.topLoop) declared[1].loop.push_back(rm.mesh.pos[v]);
    TypedCheckOptions o;
    o.w = 5.5;
    o.h = 3.4;
    o.boundaryTolerance = 0.85;
    o.requireStrips = true;
    const TypedMeshCheck k = checkTypedMesh(p.mesh, rm.mesh, declared, o);
    CHECK(k.connectivity == Status::Pass && k.decrease == 11 && k.stripsBroken == 0 && k.castOnEdges == 34 && k.bindOffEdges == 23);
    VirtualShift v;
    CHECK(buildVirtualShift(p.mesh, ShiftReading::PathSideDown, v, &err));
    const VirtualShiftCheck vc = validateVirtualShift(p.mesh, v);
    CHECK(vc.status == Status::Pass && vc.longestStrand == p.mesh.mesh.faces());
    CoverageOptions co;
    co.coverDistance = 0.55;
    co.boundaryBand = 1.7;
    co.hausdorffLimit = 0.85;
    const CoverageReport cov = measureCoverage(rm.mesh, {rm.bottomLoop, rm.topLoop}, p.mesh.mesh, co);
    CHECK(cov.status == Status::Pass && cov.coveredTwice == 0);
    PolygonMesh holed;
    holed.pos = p.mesh.mesh.pos;
    for (int f = 0; f < p.mesh.mesh.faces(); ++f) {
        if (p.mesh.strip[f] == 20) continue;
        std::vector<int> c;
        for (int i = 0; i < p.mesh.mesh.sides(f); ++i) c.push_back(p.mesh.mesh.vertexAt(f, i));
        holed.addFace(c);
    }
    const CoverageReport gap = measureCoverage(rm.mesh, {rm.bottomLoop, rm.topLoop}, holed, co);
    CHECK(gap.status == Status::Fail && gap.uncoveredFraction > 0.005 && gap.outputLoops == 4);
    // a boundary band wider than the piece leaves no sample to measure: not measurable, never a pass and never "100 %"
    // (defect: an empty considered area read as everything uncovered, a false failure on one-row rim strips)
    CoverageOptions wide = co;
    wide.boundaryBand = 1000;
    const CoverageReport none = measureCoverage(rm.mesh, {rm.bottomLoop, rm.topLoop}, p.mesh.mesh, wide);
    CHECK(!none.uncoveredMeasured && none.uncoveredFraction == 0 && none.status == Status::Unresolved && none.detail.find("not measurable") != std::string::npos);
    CHECK(cov.uncoveredMeasured);
}

namespace {

// A tube of 16 around at the radius whose chord is exactly w, `rows` strips of height h, on a reference cylinder that
// ends at the border rings.
struct RelaxTube {
    TypedStitchMesh mesh;
    TriMesh reference;
    std::vector<std::vector<int>> loops;
    double radius = 0;
};

RelaxTube relaxTube(const std::vector<int>& counts, const std::vector<std::vector<int>>& shaping) {
    RelaxTube t;
    t.radius = 5.5 / (2 * std::sin(kPi / 16));
    TriMesh unused;
    std::vector<DeclaredBoundary> declared;
    std::string err;
    CHECK(ringTubeFixture(counts, shaping, t.radius, 5.5, t.mesh, unused, declared, &err));
    const int R = static_cast<int>(counts.size()) - 1;
    t.reference = cylinderReference(t.radius, 0, R * 5.5, 5.5 / 8);
    t.loops = triangleBoundaryLoops(t.reference);
    return t;
}

std::vector<char> borderOf(const RelaxTube& t) {
    std::vector<char> b(t.mesh.mesh.pos.size(), 0);
    for (const std::vector<int>& loop : polygonBorderLoops(t.mesh.mesh))
        for (int v : loop) b[v] = 1;
    return b;
}

void perturb(RelaxTube& t, double amount) {
    const double top = (t.mesh.rings.size() - 1) * 5.5;
    for (size_t v = 0; v < t.mesh.mesh.pos.size(); ++v) {
        const Vec3 p = t.mesh.mesh.pos[v];
        if (p.z < 1e-9 || p.z > top - 1e-9) continue;
        const double u = hash01(7, static_cast<int>(v), 0, 0) - 0.5, s = hash01(7, static_cast<int>(v), 1, 0) - 0.5;
        const double phi = std::atan2(p.y, p.x) + 2 * amount * 5.5 * u / t.radius;
        t.mesh.mesh.pos[v] = Vec3{t.radius * std::cos(phi), t.radius * std::sin(phi), p.z + 2 * amount * 5.5 * s};
    }
}

double worstRatioError(const TypedStitchMesh& s) {
    double worst = 0;
    for (const CellMeasure& c : measureCells(s, 5.5, 5.5))
        if (c.kind == StitchFaceKind::Regular) worst = std::max(worst, std::max(std::fabs(c.widthRatio - 1), std::fabs(c.heightRatio - 1)));
    return worst;
}

} // namespace

// The model's forces conserve momentum (every term acts in equal and opposite parts) and vanish on the rest state.
TEST(ssm_relax_forces_balance) {
    RelaxTube t = relaxTube(std::vector<int>(7, 16), std::vector<std::vector<int>>(6));
    MeshRelaxRequest r;
    const std::vector<RelaxFace> faces = relaxFacesOf(t.mesh.mesh, t.mesh.role);
    std::vector<Vec3> x;
    for (const Vec3& p : t.mesh.mesh.pos) x.push_back(p / 5.5);
    double restTangential = 0;   // on the rest tube every force is radial: no tangential part
    for (size_t v = 0; v < x.size(); ++v) {
        const Vec3 f = relaxForces(faces, x, r)[v];
        const Vec3 radial = normalized(Vec3{x[v].x, x[v].y, 0});
        restTangential = std::max(restTangential, norm(f - radial * dot(f, radial)));
    }
    perturb(t, 0.15);
    x.clear();
    for (const Vec3& p : t.mesh.mesh.pos) x.push_back(p / 5.5);
    Vec3 sum{0, 0, 0};
    double largest = 0;
    for (const Vec3& f : relaxForces(faces, x, r)) {
        sum += f;
        largest = std::max(largest, norm(f));
    }
    CHECK(largest > 0.1 && norm(sum) < 1e-9 * largest * x.size());
    CHECK(restTangential < 1e-9);
}

// The rest state is kept (0 steps), a perturbed tube relaxes back to w x h cells, every vertex stays on the reference and
// every border vertex on its border ring.
TEST(ssm_relax_rest_and_perturbation) {
    RelaxTube rest = relaxTube(std::vector<int>(7, 16), std::vector<std::vector<int>>(6));
    MeshRelaxRequest r;
    RelaxProduct p;
    std::string err;
    CHECK(relaxStitchMesh(rest.reference, rest.loops, relaxFacesOf(rest.mesh.mesh, rest.mesh.role), rest.mesh.mesh.pos, borderOf(rest), r, p, &err));
    CHECK(p.stopReason == "converged" && p.steps == 0);
    RelaxTube t = relaxTube(std::vector<int>(7, 16), std::vector<std::vector<int>>(6));
    perturb(t, 0.15);
    CHECK(worstRatioError(t.mesh) > 0.2);
    CHECK(relaxStitchMesh(t.reference, t.loops, relaxFacesOf(t.mesh.mesh, t.mesh.role), t.mesh.mesh.pos, borderOf(t), r, p, &err));
    CHECK(p.stopReason == "converged" && p.steps < 100 && p.walkClamps == 0);
    TypedStitchMesh after = t.mesh;
    after.mesh.pos = p.pos;
    CHECK(worstRatioError(after) < 0.01);
    const MeshQuery q(t.reference);
    const std::vector<char> border = borderOf(t);
    for (size_t v = 0; v < p.pos.size(); ++v) {
        CHECK(q.closest(p.pos[v]).dist < 1e-9);
        if (border[v]) CHECK(std::fabs(p.pos[v].z) < 1e-9 || std::fabs(p.pos[v].z - 6 * 5.5) < 1e-9);
    }
}

// On Method A's virtual faces the rings turn into a helix (each stitch climbs; the flat border rings limit it to about
// half of one row per round on this short tube), on the stored faces they stay flat; the virtual cells meet the frozen
// statistics; kp and d12k converge on both.
TEST(ssm_relax_virtual_helix) {
    for (const bool virt : {false, true}) {
        RelaxTube t = relaxTube(std::vector<int>(7, 16), std::vector<std::vector<int>>(6));
        std::vector<RelaxFace> faces = relaxFacesOf(t.mesh.mesh, t.mesh.role);
        VirtualShift vs;
        std::string err;
        if (virt) {
            CHECK(buildVirtualShift(t.mesh, ShiftReading::PathSideDown, vs, &err) && validateVirtualShift(t.mesh, vs).status == Status::Pass);
            for (int f = 0; f < t.mesh.mesh.faces(); ++f) faces[f].corners = vs.corners[f];
        }
        RelaxProduct p;
        CHECK(relaxStitchMesh(t.reference, t.loops, faces, t.mesh.mesh.pos, borderOf(t), MeshRelaxRequest{}, p, &err));
        CHECK(p.stopReason == "converged");
        double rise = 0;
        int n = 0;
        for (int k = 1; k < 6; ++k)
            for (int j = 1; j + 1 < 16; ++j) {
                rise += p.pos[t.mesh.rings[k][j + 1]].z - p.pos[t.mesh.rings[k][j]].z;
                ++n;
            }
        rise /= n;
        if (virt) {
            CHECK(rise > 0.1 && rise < 5.5 / 16 + 1e-9);
            TypedStitchMesh after = t.mesh;
            after.mesh.pos = p.pos;
            const TypedStitchMesh va = virtualTypedMesh(after, vs);
            std::vector<double> err2;
            for (const CellMeasure& c : measureCells(va, 5.5, 5.5))
                if (c.kind == StitchFaceKind::Regular) {
                    err2.push_back(std::fabs(c.widthRatio - 1));
                    err2.push_back(std::fabs(c.heightRatio - 1));
                }
            std::sort(err2.begin(), err2.end());
            CHECK(err2[static_cast<size_t>(0.95 * (err2.size() - 1))] <= 0.15);
        } else {
            CHECK(std::fabs(rise) < 1e-6);
        }
    }
    for (const bool inc : {true, false}) {
        RelaxTube t = inc ? relaxTube({16, 16, 16, 17, 17, 17, 17}, {{}, {}, {8}, {}, {}, {}}) : relaxTube({17, 17, 17, 16, 16, 16, 16}, {{}, {}, {8}, {}, {}, {}});
        RelaxProduct p;
        std::string err;
        CHECK(relaxStitchMesh(t.reference, t.loops, relaxFacesOf(t.mesh.mesh, t.mesh.role), t.mesh.mesh.pos, borderOf(t), MeshRelaxRequest{}, p, &err));
        CHECK(p.stopReason == "converged" && p.subQuads == t.mesh.mesh.faces() + 1);
    }
}

// The stall rule stops a solve that makes no progress (a timestep of 1e-6) at its window, not at the step cap.
TEST(ssm_relax_stall_rule_stops) {
    RelaxTube t = relaxTube(std::vector<int>(7, 16), std::vector<std::vector<int>>(6));
    perturb(t, 0.15);
    MeshRelaxRequest r;
    r.timestep = 1e-6;
    RelaxProduct p;
    std::string err;
    CHECK(relaxStitchMesh(t.reference, t.loops, relaxFacesOf(t.mesh.mesh, t.mesh.role), t.mesh.mesh.pos, borderOf(t), r, p, &err));
    CHECK(p.stopReason == "stalled" && p.steps == r.stallWindow && p.status == Status::Unresolved);
}

namespace {

std::set<std::pair<int, int>> borderEdgeSet(const PolygonMesh& m) {
    std::map<std::pair<int, int>, int> count;
    for (int f = 0; f < m.faces(); ++f)
        for (int k = 0; k < m.sides(f); ++k) {
            const int a = m.vertexAt(f, k), b = m.vertexAt(f, (k + 1) % m.sides(f));
            if (a != b) count[{std::min(a, b), std::max(a, b)}] += 1;
        }
    std::set<std::pair<int, int>> out;
    for (const auto& [e, n] : count)
        if (n == 1) out.insert(e);
    return out;
}

} // namespace

// The helix boundary contract on the 8 x 3 fixture:
// the shifted faces' border equals the stored border except for exactly three named edges: S (the stored top edge of
// slot 0 in the last strip) leaves it, E (the last path edge) and D (the slanted top edge of slot 0 in the last strip)
// enter it; every cast-on edge is kept; the start step is the only face whose shifted polygon is a triangle.
TEST(ssm_helix_boundary_contract) {
    const std::vector<FixtureCase> all = constructionFixtureCases();
    const FixtureCase& c = findCase(all, "k-open-tube-8x3");
    VirtualShift v;
    std::string err;
    CHECK(buildVirtualShift(c.mesh, ShiftReading::PathSideDown, v, &err));
    const TypedStitchMesh vt = virtualTypedMesh(c.mesh, v);
    const std::set<std::pair<int, int>> stored = borderEdgeSet(c.mesh.mesh), shifted = borderEdgeSet(vt.mesh);
    const std::vector<std::vector<int>>& rings = c.mesh.rings;
    const int R = static_cast<int>(rings.size()) - 1;
    auto key = [](int a, int b) { return std::make_pair(std::min(a, b), std::max(a, b)); };
    const std::pair<int, int> S = key(rings[R][0], rings[R][1]), E = key(rings[R - 1][0], rings[R][0]), D = key(rings[R - 1][0], rings[R][1]);
    std::set<std::pair<int, int>> removed, added;
    for (const auto& e : stored)
        if (!shifted.count(e)) removed.insert(e);
    for (const auto& e : shifted)
        if (!stored.count(e)) added.insert(e);
    CHECK((removed == std::set<std::pair<int, int>>{S}));
    CHECK((added == std::set<std::pair<int, int>>{E, D}));
    for (size_t j = 0; j < rings[0].size(); ++j) CHECK(shifted.count(key(rings[0][j], rings[0][(j + 1) % rings[0].size()])) == 1);
    int triangles = 0;
    for (int f = 0; f < vt.mesh.faces(); ++f)
        if (vt.mesh.sides(f) == 3) {
            ++triangles;
            CHECK(c.mesh.strip[f] == 0 && c.mesh.slot[f] == 0 && vt.reason[f] == "helix-start-step");
        }
    CHECK(triangles == 1);
}

// ---------------------------------------------------------------------------------------------------------------------
// The imported arm: marking records, the pole start, adjacent shaping, Method A
// on markings, the relaxation's marking constraints.

namespace {

SurfacePoint vertexOn(const TriMesh& m, int v) {
    SurfacePoint p;
    for (int t = 0; t < static_cast<int>(m.tri.size()); ++t)
        for (int k = 0; k < 3; ++k)
            if (m.tri[t][k] == v) {
                p.tri = t;
                p.bary[0] = p.bary[1] = p.bary[2] = 0;
                p.bary[k] = 1;
                return p;
            }
    return p;
}

// A flat disc reference (the pole-disc fixture's) with markings: start at the centre, a straight row-start path along +x
// through the disc's vertices, the finish point on that path at radius `finish`.
struct DiscMarked {
    TriMesh ref;
    MarkingSet set;
    std::vector<int> pathVertices;
};

DiscMarked discMarked(double finish) {
    DiscMarked d;
    TypedStitchMesh unused;
    std::vector<DeclaredBoundary> decl;
    std::string err;
    CHECK(poleDiscFixture({5, 7, 7, 7, 7, 7}, {{1, 3}, {}, {}, {}, {}}, 3.4, unused, d.ref, decl, &err));
    // the disc's vertices on the +x axis, from the centre outwards
    std::vector<std::pair<double, int>> onAxis;
    for (int v = 0; v < static_cast<int>(d.ref.pos.size()); ++v)
        if (std::fabs(d.ref.pos[v].y) < 1e-9 && d.ref.pos[v].x >= -1e-9) onAxis.push_back({d.ref.pos[v].x, v});
    std::sort(onAxis.begin(), onAxis.end());
    d.set.component = ComponentIdentity{"disc", "fnv1a64:0", "0000000000000000", "disc", "fabric", 0};
    d.set.reference.digest = referenceDigest(d.ref);
    d.set.reference.vertices = static_cast<int>(d.ref.pos.size());
    d.set.reference.triangles = static_cast<int>(d.ref.tri.size());
    d.set.reference.frame = "test disc, mm";
    Marking st;
    st.id = "start";
    st.type = MarkingType::StartPoint;
    st.points.push_back(vertexOn(d.ref, onAxis[0].second));
    Marking pa;
    pa.id = "row_start";
    pa.type = MarkingType::RowStartPath;
    for (const auto& [x, v] : onAxis) {
        pa.points.push_back(vertexOn(d.ref, v));
        d.pathVertices.push_back(v);
    }
    Marking fi;
    fi.id = "finish";
    fi.type = MarkingType::FinishBoundary;
    // the finish point: on the axis at x = finish, inside the triangle above the axis edge that contains it
    for (size_t i = 0; i + 1 < onAxis.size(); ++i)
        if (onAxis[i].first <= finish && onAxis[i + 1].first > finish) {
            const SurfacePoint a = vertexOn(d.ref, onAxis[i].second), b = vertexOn(d.ref, onAxis[i + 1].second);
            const int t = commonTriangle(d.ref, a, b);
            const double u = (finish - onAxis[i].first) / (onAxis[i + 1].first - onAxis[i].first);
            fi.points.push_back(SurfacePoint{t, {0, 0, 0}});
            for (int k = 0; k < 3; ++k) {
                if (d.ref.tri[t][k] == onAxis[i].second) fi.points.back().bary[k] = 1 - u;
                if (d.ref.tri[t][k] == onAxis[i + 1].second) fi.points.back().bary[k] = u;
            }
        }
    Marking gd;
    gd.id = "growth";
    gd.type = MarkingType::GrowthDirection;
    gd.settings = {{"from", "start"}, {"to", "finish"}};
    d.set.markings = {st, fi, gd, pa};
    return d;
}

} // namespace

// A marking record survives the file unchanged, and the checks refuse what the contract forbids: another reference, another
// component, a point outside its triangle, a path that leaves the surface, a closed row-start path, an attachment claiming
// a join. A physical cut is listed as not consumed; automatic proposals are listed as not reviewed.
TEST(ssm_markings_round_trip_and_refusals) {
    DiscMarked d = discMarked(14.0);
    const std::string file = (std::filesystem::temp_directory_path() / "ssm_markings_test.markings").string();
    std::string err;
    CHECK(writeMarkings(file, d.set, &err));
    MarkingSet back;
    CHECK(readMarkings(file, back, &err));
    CHECK(back.markings.size() == d.set.markings.size());
    CHECK(back.component == d.set.component && back.reference.digest == d.set.reference.digest);
    for (size_t i = 0; i < back.markings.size(); ++i) {
        CHECK(back.markings[i].points.size() == d.set.markings[i].points.size());
        for (size_t j = 0; j < back.markings[i].points.size(); ++j) {
            CHECK(back.markings[i].points[j].tri == d.set.markings[i].points[j].tri);
            CHECK(back.markings[i].points[j].bary[1] == d.set.markings[i].points[j].bary[1]);
        }
    }
    const MarkingCheck ok = validateMarkings(back, d.ref, d.set.component);
    CHECK(ok.status == Status::Unresolved && ok.refusals.empty() && ok.proposals == 4);
    auto refused = [&](MarkingSet s, const ComponentIdentity& expected) { return validateMarkings(s, d.ref, expected).status == Status::Fail; };
    {
        MarkingSet s = back;
        s.reference.digest = "0123456789abcdef";
        CHECK(refused(s, d.set.component));
    }
    {
        ComponentIdentity other = d.set.component;
        other.piece = 1;
        CHECK(refused(back, other));
    }
    {
        MarkingSet s = back;
        s.markings[0].points[0].bary[0] = 1.5;
        s.markings[0].points[0].bary[1] = -0.5;
        CHECK(refused(s, d.set.component));
    }
    {
        MarkingSet s = back;
        Marking& p = s.markings[3];
        p.points.erase(p.points.begin() + 2, p.points.begin() + 6);   // a gap: the path no longer lies on the surface
        CHECK(refused(s, d.set.component));
    }
    {
        MarkingSet s = back;
        s.markings[3].closed = true;
        CHECK(refused(s, d.set.component));
    }
    {
        MarkingSet s = back;
        Marking at;
        at.id = "attach";
        at.type = MarkingType::AttachmentRegion;
        at.triangles = {0, 1, 2};
        at.settings = {{"other", "component_9"}, {"join", "established"}};
        s.markings.push_back(at);
        CHECK(refused(s, d.set.component));
        s.markings.back().settings[1].second = "unresolved";
        CHECK(!refused(s, d.set.component));
    }
    {
        MarkingSet s = back;
        Marking cut = s.markings[3];
        cut.id = "cut";
        cut.type = MarkingType::Cut;
        s.markings.push_back(cut);
        const MarkingCheck c = validateMarkings(s, d.ref, d.set.component);
        CHECK(c.status == Status::Unresolved);
        bool listed = false;
        for (const std::string& o : c.open) listed = listed || o.find("does not consume cuts") != std::string::npos;
        CHECK(listed);
    }
    std::filesystem::remove(file);
}

// A changed triangulation: the markings are refused on it as they are; a checked rebinding moves every point (none farther
// than the tolerance), densifies the path across the finer triangles, and the result validates on the new reference.
TEST(ssm_markings_rebind_onto_a_refinement) {
    DiscMarked coarse = discMarked(14.0);
    TypedStitchMesh unused;
    std::vector<DeclaredBoundary> decl;
    std::string err;
    TriMesh fine;
    CHECK(poleDiscFixture({5, 7, 7, 7, 7, 7}, {{1, 3}, {}, {}, {}, {}}, 3.4, unused, fine, decl, &err));
    // refine: every triangle into four at its edge midpoints (the same plane)
    {
        TriMesh r;
        r.pos = fine.pos;
        std::map<std::pair<int, int>, int> mid;
        auto midpoint = [&](int a, int b) {
            const std::pair<int, int> k{std::min(a, b), std::max(a, b)};
            if (!mid.count(k)) {
                mid[k] = static_cast<int>(r.pos.size());
                r.pos.push_back(0.5 * (fine.pos[a] + fine.pos[b]));
            }
            return mid[k];
        };
        for (const auto& t : fine.tri) {
            const int a = midpoint(t[0], t[1]), b = midpoint(t[1], t[2]), c = midpoint(t[2], t[0]);
            r.tri.push_back({t[0], a, c});
            r.tri.push_back({a, t[1], b});
            r.tri.push_back({c, b, t[2]});
            r.tri.push_back({a, b, c});
        }
        r.computeNormals();
        fine = r;
    }
    CHECK(validateMarkings(coarse.set, fine, coarse.set.component).status == Status::Fail);
    MarkingSet moved;
    RebindReport rep;
    CHECK(rebindMarkings(coarse.set, coarse.ref, fine, 1e-9, moved, rep));
    CHECK(rep.refused == 0 && rep.maxDistance < 1e-9 && rep.inserted > 0);
    CHECK(validateMarkings(moved, fine, coarse.set.component).refusals.empty());
    // a different surface: the disc lifted by 1 mm, so every point is 1 mm away
    TriMesh lifted = fine;
    for (Vec3& p : lifted.pos) p.z += 1.0;
    MarkingSet none;
    CHECK(!rebindMarkings(coarse.set, coarse.ref, lifted, 1e-3, none, rep));
}

// The arm's construction fixtures (pole start, adjacent shaping) and their negative controls come out as written.
TEST(ssm_arm_fixtures_as_expected) {
    const std::vector<FixtureCase> all = armConstructionFixtureCases();
    CHECK(all.size() == 11);
    for (const FixtureCase& c : all) {
        const TypedMeshCheck k = checkTypedMesh(c.mesh, c.reference, c.declared, fixtureOptions());
        CHECK(k.connectivity == c.expectConnectivity);
        if (c.runShift) CHECK(shiftOf(c).status == c.expectShift);
        if (!c.negative) CHECK(k.geometry == Status::Pass);
    }
    // the pole round: one closed strand of 5 pole-start faces; the helix through the other 12, starting in strip 1
    const VirtualShiftCheck v = shiftOf(findCase(all, "pole-start-5-7-7"));
    CHECK(v.poleStart && v.poleRoundFaces == 5 && v.longestStrand == 12 && v.startsInFirstStrip);
}

// Method A on markings: the row origin follows the marked path (ring vertex 0 on it), the pole is the start point, the
// bind-off ring is the finish level; a path that turns back is a guide problem (no mesh, the crossings listed).
TEST(ssm_marked_rows_follow_the_path) {
    const DiscMarked d = discMarked(14.0);
    std::vector<double> t;
    TimeFieldReport rep;
    std::string err;
    CHECK(geodesicDistance(d.ref, {d.pathVertices[0]}, t, rep, &err));
    auto request = [&](const std::vector<int>& vertices) {
        MarkedRowsRequest r;
        r.w = 5.5;
        r.h = 3.4;
        r.pole = d.ref.pos[d.pathVertices[0]];
        r.poleBinding.tri = d.set.markings[0].points[0].tri;
        for (int k = 0; k < 3; ++k) r.poleBinding.bary[k] = d.set.markings[0].points[0].bary[k];
        r.finishLevel = 14.0;
        r.adjacentShaping = true;
        for (size_t i = 0; i < vertices.size(); ++i) {
            MarkedPathPoint p;
            p.p = d.ref.pos[vertices[i]];
            p.t = t[vertices[i]];
            p.segmentTri = i + 1 < vertices.size() ? commonTriangle(d.ref, vertexOn(d.ref, vertices[i]), vertexOn(d.ref, vertices[i + 1])) : -1;
            r.path.push_back(p);
        }
        return r;
    };
    MarkedRowsProduct ok;
    CHECK(buildMarkedRows(d.ref, t, request(d.pathVertices), ok, &err));
    CHECK(ok.guideProblems == 0 && ok.rows.contourEvents == 0 && ok.rows.mesh.mesh.faces() > 0);
    CHECK(ok.rows.rows == 4);   // round(14 / 3.4)
    for (size_t k = 1; k < ok.rows.mesh.rings.size(); ++k) {
        const Vec3 origin = ok.rows.mesh.mesh.pos[ok.rows.mesh.rings[k][0]];
        CHECK(std::fabs(origin.y) < 1e-6 && origin.x > 0);   // on the +x path
    }
    VirtualShift vs;
    CHECK(buildVirtualShift(ok.rows.mesh, ShiftReading::PathSideDown, vs, &err));
    const VirtualShiftCheck vc = validateVirtualShift(ok.rows.mesh, vs);
    CHECK(vc.poleStart && vc.status == Status::Pass);
    // the path turned back across the level 7 mm (ring 2, h' = 3.5): out to 8.5 mm (vertex 10, 0.85 mm apart), back to
    // 5.1 mm, out again
    std::vector<int> back(d.pathVertices.begin(), d.pathVertices.begin() + 11);
    for (int i = 9; i >= 6; --i) back.push_back(d.pathVertices[i]);
    back.insert(back.end(), d.pathVertices.begin() + 7, d.pathVertices.end());
    MarkedRowsProduct bad;
    CHECK(buildMarkedRows(d.ref, t, request(back), bad, &err));
    CHECK(bad.rows.status == Status::Fail && bad.guideProblems > 0 && bad.rows.mesh.mesh.faces() == 0);
}

// The relaxation keeps what the markings fix: a pinned vertex does not move, a guided vertex stays on its line.
TEST(ssm_relax_guides_hold) {
    TypedStitchMesh s;
    TriMesh ref;
    std::vector<DeclaredBoundary> decl;
    std::string err;
    CHECK(poleDiscFixture({5, 7, 7}, {{1, 3}, {}}, 3.4, s, ref, decl, &err));
    std::vector<Vec3> start = s.mesh.pos;
    for (size_t v = 1; v < start.size(); ++v) start[v] = start[v] + Vec3{0.3 * std::sin(3.0 * v), 0.3 * std::cos(5.0 * v), 0};   // perturbed in the plane
    GuideLine line;   // the +x axis from the centre to the edge of the disc
    for (int i = 0; i <= 40; ++i) line.points.push_back({0.4 * i, 0, 0});
    const MeshQuery q(ref);
    for (int i = 0; i < 40; ++i) line.segmentTri.push_back(q.closest(0.5 * (line.points[i] + line.points[i + 1]) + Vec3{0, 1e-6, 0}).tri);
    RelaxGuides g;
    g.lines = {line};
    g.lineOf.assign(start.size(), -1);
    g.pinned.assign(start.size(), 0);
    g.pinned[0] = 1;
    g.lineOf[s.rings[1][0]] = 0;
    MeshRelaxRequest r;
    r.w = 5.5;
    r.h = 3.4;
    r.maxSteps = 200;
    RelaxProduct p;
    std::vector<RelaxFace> faces = relaxFacesOf(s.mesh, s.role);
    CHECK(relaxStitchMesh(ref, {}, faces, start, std::vector<char>(start.size(), 0), r, p, &err, &g));
    CHECK(norm(p.pos[0] - start[0]) < 1e-12);
    CHECK(std::fabs(p.pos[s.rings[1][0]].y) < 1e-9);
    CHECK(p.pinnedVertices == 1 && p.guidedVertices == 1);
}

// The border ring on a loop with bends (border_placement.h).
// A square loop of 25 stitch widths around, the ring's corners evenly spaced from a first corner 5 mm past a bend: the
// sides across the other three bends cut them off by 1.7 mm. Placed, the same 25 corners keep within 0.85 mm, the first
// one stays, the order around the loop stays, and no side is shorter or longer than 30 % off the mean. With the ring left
// as it is (the defect) the first check fails.
TEST(ssm_border_ring_is_placed_on_the_bends) {
    const double w = 5.5, side = 25 * w / 4;
    const std::vector<Vec3> loop = {{0, 0, 0}, {side, 0, 0}, {side, side, 0}, {0, side, 0}};
    auto onLoop = [&](double arc) {
        arc = std::fmod(arc + 4 * side, 4 * side);
        const int edge = std::min(3, static_cast<int>(arc / side));
        const Vec3 a = loop[edge], b = loop[(edge + 1) % 4];
        return a + (b - a) * ((arc - edge * side) / side);
    };
    std::vector<Vec3> ring;
    for (int k = 0; k < 25; ++k) ring.push_back(onLoop(5.0 + k * w));
    BorderRingRequest request;
    request.limit = 0.85;
    const BorderRingPlacement p = placeBorderRing(loop, ring, request);
    std::printf("  border ring: %.3f -> %.3f mm, sides %.2f to %.2f mm, spread %.0f %%\n", p.distanceBefore, p.distance, p.shortestSide, p.longestSide, 100 * p.spread);
    CHECK(p.placed && p.inOrder && p.distanceBefore > 1.5);
    CHECK(p.distance <= 0.85);
    CHECK(norm(p.corners[0] - ring[0]) <= 0.26);                       // the yarn's end stays (to the loop's sampling step)
    CHECK(p.shortestSide >= 0.7 * w && p.longestSide <= 1.3 * w);
    // the order around the loop: the arc position of the corners grows from the first one
    auto arcOf = [&](const Vec3& q) {
        double best = 1e300, arc = 0;
        for (int e = 0; e < 4; ++e) {
            const Vec3 a = loop[e], d = loop[(e + 1) % 4] - a;
            const double t = std::clamp(dot(q - a, d) / dot(d, d), 0.0, 1.0);
            const double off = norm(q - (a + d * t));
            if (off < best) {
                best = off;
                arc = (e + t) * side;
            }
        }
        return arc;
    };
    double previous = arcOf(p.corners[0]), turned = 0;
    for (int k = 1; k <= 25; ++k) {
        const double arc = arcOf(p.corners[k % 25]);
        const double step = std::fmod(arc - previous + 4 * side, 4 * side);
        CHECK(step > 0 && step < 2 * w);
        turned += step;
        previous = arc;
    }
    CHECK(std::fabs(turned - 4 * side) < 1e-6);
    // a ring whose two neighbours have changed places is still placed, and reported as out of order
    std::swap(ring[7], ring[8]);
    const BorderRingPlacement tangled = placeBorderRing(loop, ring, request);
    CHECK(tangled.placed && !tangled.inOrder && tangled.distance <= 0.85);
}

// The teddy mesh repair: what the relaxation keeps of a valid mesh.
namespace {

// Corners of the tube's faces that are turned over: (b - a) x (c - b) against the outward radial direction at b.
int turnedCorners(const std::vector<RelaxFace>& faces, const std::vector<Vec3>& pos) {
    int turned = 0;
    for (const RelaxFace& f : faces) {
        const int n = static_cast<int>(f.corners.size());
        for (int k = 0; k < n; ++k) {
            const Vec3 a = pos[f.corners[k]], b = pos[f.corners[(k + 1) % n]], c = pos[f.corners[(k + 2) % n]];
            if (dot(cross(b - a, c - b), normalized(Vec3{b.x, b.y, 0})) < 0) ++turned;
        }
    }
    return turned;
}

MeshRelaxRequest relaxAsBefore() {
    MeshRelaxRequest r;
    r.borderOrder = BorderOrder::Free;
    r.faceCorners = FaceCorners::Free;
    r.stepControl = StepControl::FixedTimeStep;
    return r;
}

} // namespace

// The first vertex of the bind-off ring is held and the rows below stand two stitches further round, so the ring is drawn
// back past it. The faces of the solver are the virtual ones: the edge from the first to the second vertex of the ring
// belongs to none of them. As before (no spring there, each border vertex on its own) the second vertex passes the held
// one. With the spring of that edge and the ring moving as one, no vertex passes another and the ring keeps its order.
TEST(ssm_relax_border_ring_keeps_its_order) {
    RelaxTube t = relaxTube(std::vector<int>(7, 16), std::vector<std::vector<int>>(6));
    std::vector<RelaxFace> faces = relaxFacesOf(t.mesh.mesh, t.mesh.role);
    VirtualShift vs;
    std::string err;
    CHECK(buildVirtualShift(t.mesh, ShiftReading::PathSideDown, vs, &err));
    for (int f = 0; f < t.mesh.mesh.faces(); ++f) faces[f].corners = vs.corners[f];
    const std::vector<int> top = t.mesh.rings.back();
    std::vector<char> onTop(t.mesh.mesh.pos.size(), 0);
    for (int v : top) onTop[v] = 1;
    std::vector<Vec3> start = t.mesh.mesh.pos;
    for (size_t v = 0; v < start.size(); ++v) {
        if (onTop[v]) continue;
        const double phi = std::atan2(start[v].y, start[v].x) - 2 * 5.5 / t.radius;
        start[v] = Vec3{t.radius * std::cos(phi), t.radius * std::sin(phi), start[v].z};
    }
    RelaxGuides guides;
    guides.pinned.assign(start.size(), 0);
    guides.lineOf.assign(start.size(), -1);
    guides.pinned[top[0]] = 1;
    RelaxProduct before;
    CHECK(relaxStitchMesh(t.reference, t.loops, faces, start, borderOf(t), relaxAsBefore(), before, &err, &guides));
    std::printf("  as before: %d order violations; %s\n", before.borderOrderViolations, before.firstBorderOrderViolation.c_str());
    CHECK(before.borderOrderViolations > 0);
    guides.courseEdges.push_back({top[0], top[1]});
    RelaxProduct now;
    CHECK(relaxStitchMesh(t.reference, t.loops, faces, start, borderOf(t), MeshRelaxRequest{}, now, &err, &guides));
    CHECK(now.addedCourseSprings == 1 && now.borderOrderViolations == 0);
    double turned = 0;
    for (size_t j = 0; j < top.size(); ++j) {
        const Vec3 a = now.pos[top[j]], b = now.pos[top[(j + 1) % top.size()]];
        double step = std::atan2(b.y, b.x) - std::atan2(a.y, a.x);
        if (step < -kPi) step += 2 * kPi;
        CHECK(step > 0.25 * 5.5 / t.radius);   // every vertex still ahead of the one before it, by a quarter stitch or more
        turned += step;
    }
    CHECK(std::fabs(turned - 2 * kPi) < 1e-6);
}

// The shear force is applied explicitly. With a shear stiffness of 5 (the paper: 0.2) a fixed time step no longer lowers
// the forces and the run ends as stalled at its first window. With steps that raise the forces taken back and the time
// step halved, the same start converges; no corner of a face is turned over on the way or at the end.
TEST(ssm_relax_takes_back_steps_that_raise_the_forces) {
    RelaxTube t = relaxTube(std::vector<int>(7, 16), std::vector<std::vector<int>>(6));
    perturb(t, 0.3);
    const std::vector<RelaxFace> faces = relaxFacesOf(t.mesh.mesh, t.mesh.role);
    CHECK(turnedCorners(faces, t.mesh.mesh.pos) == 0);
    MeshRelaxRequest stiff = relaxAsBefore();
    stiff.kShear = 5;
    RelaxProduct before, now;
    std::string err;
    CHECK(relaxStitchMesh(t.reference, t.loops, faces, t.mesh.mesh.pos, borderOf(t), stiff, before, &err));
    MeshRelaxRequest kept;
    kept.kShear = 5;
    CHECK(relaxStitchMesh(t.reference, t.loops, faces, t.mesh.mesh.pos, borderOf(t), kept, now, &err));
    std::printf("  as before: %s after %d steps, RMS force %.3g; now: %s after %d steps, %.3g, %d steps taken back\n", before.stopReason.c_str(), before.steps, before.rmsForceEnd,
                now.stopReason.c_str(), now.steps, now.rmsForceEnd, now.stepsTakenBack);
    CHECK(before.stopReason == "stalled");
    CHECK(now.stopReason == "converged" && now.stepsTakenBack > 0);
    CHECK(turnedCorners(faces, now.pos) == 0);
}

// The first corner of a border ring may move by up to the freedom the request gives: the fit is then at least as close,
// and the corner stays within that distance of its place.
TEST(ssm_border_ring_first_corner_may_move) {
    const double w = 5.5, side = 25 * w / 4;
    const std::vector<Vec3> loop = {{0, 0, 0}, {side, 0, 0}, {side, side, 0}, {0, side, 0}};
    std::vector<Vec3> ring;
    for (int k = 0; k < 25; ++k) {
        const double arc = std::fmod(2.0 + k * w, 4 * side);
        const int edge = std::min(3, static_cast<int>(arc / side));
        ring.push_back(loop[edge] + (loop[(edge + 1) % 4] - loop[edge]) * ((arc - edge * side) / side));
    }
    BorderRingRequest held, free;
    free.anchorFreedom = 0.5 * w;
    const BorderRingPlacement a = placeBorderRing(loop, ring, held), b = placeBorderRing(loop, ring, free);
    std::printf("  first corner held: %.3f mm; free within %.2f mm: %.3f mm, moved %.2f mm\n", a.distance, free.anchorFreedom, b.distance, b.firstCornerMoved);
    CHECK(a.placed && b.placed && a.firstCornerMoved == 0);
    CHECK(b.distance <= a.distance + 1e-9 && std::fabs(b.firstCornerMoved) <= free.anchorFreedom + 1e-9);
    CHECK(norm(b.corners[0] - ring[0]) <= free.anchorFreedom + 0.26);
}
