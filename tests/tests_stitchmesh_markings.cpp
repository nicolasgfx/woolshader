// knit_tests_stitchmesh, surface markings (src/stitchmesh/surface_trim.h, guided_rows.h's SurfacePathGuide,
// marked_build.h): the trim along closed curves, the drawn row-start path, Method A on a marked piece. Each case has a
// control that must fail.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <set>
#include <utility>

#include "harness.h"
#include "model/surface_path.h"
#include "model/surface_pieces.h"
#include "stitchmesh/construction_fixtures.h"
#include "stitchmesh/coverage.h"
#include "stitchmesh/guided_rows.h"
#include "stitchmesh/marked_build.h"
#include "stitchmesh/prepared_reference.h"
#include "stitchmesh/relax_mesh.h"
#include "stitchmesh/revolution_surface.h"
#include "stitchmesh/row_rings.h"
#include "stitchmesh/surface_trim.h"
#include "stitchmesh/time_field.h"
#include "stitchmesh/typed_mesh.h"
#include "stitchmesh/virtual_shift.h"
#include "stitchmesh/yarn_embed.h"

using namespace knit;

namespace {

RevolutionMesh markingTube() {
    RevolutionMesh rm;
    std::string err;
    CHECK(makeRevolutionMesh(RevolutionProfile{}, RevolutionMeshParams{}, rm, &err));
    return rm;
}

std::vector<double> fieldFromBottom(const RevolutionMesh& rm) {
    std::vector<double> t;
    TimeFieldReport tr;
    std::string err;
    CHECK(geodesicDistance(rm.mesh, rm.bottomLoop, t, tr, &err));
    return t;
}

// One triangle's piece of the meridian half-plane.
struct MeridianPiece {
    Vec3 lo, hi;
    int tri = -1;
};

bool lowerFirst(const MeridianPiece& a, const MeridianPiece& b) { return a.lo.z < b.lo.z; }

// The meridian half-plane at azimuth phi between two heights, drawn as a path: its segment in every triangle it
// crosses, chained by height.
SurfacePathGuide meridianGuide(const TriMesh& m, double phi, double zFrom, double zTo) {
    const Vec3 normal{-std::sin(phi), std::cos(phi), 0}, side{std::cos(phi), std::sin(phi), 0};
    std::vector<MeridianPiece> pieces;
    for (int t = 0; t < static_cast<int>(m.tri.size()); ++t) {
        std::vector<Vec3> x;
        for (int k = 0; k < 3; ++k) {
            const Vec3& p = m.pos[m.tri[t][k]];
            const Vec3& q = m.pos[m.tri[t][(k + 1) % 3]];
            const double gp = dot(p, normal), gq = dot(q, normal);
            if ((gp < 0) == (gq < 0)) continue;
            const Vec3 c = p + (q - p) * (gp / (gp - gq));
            if (dot(c, side) > 0) x.push_back(c);
        }
        if (x.size() != 2) continue;
        if (x[0].z > x[1].z) std::swap(x[0], x[1]);
        if (x[1].z <= zFrom || x[0].z >= zTo) continue;
        pieces.push_back(MeridianPiece{x[0], x[1], t});
    }
    std::sort(pieces.begin(), pieces.end(), lowerFirst);
    SurfacePathGuide g;
    for (const MeridianPiece& p : pieces) {
        if (g.points.empty()) g.points.push_back(p.lo);
        g.points.push_back(p.hi);
        g.segmentTri.push_back(p.tri);
    }
    g.segmentTri.push_back(-1);  // the last point continues nowhere
    return g;
}

int capsuleVertex(int segments, int ring, int i) { return 1 + ring * segments + (i % segments); }

// A closed capsule round z: a cylinder of radius r from z = 0 to z = length with hemispherical ends, aligned rings.
TriMesh capsule(double r, double length, int segments) {
    TriMesh m;
    std::vector<double> ringZ, ringR;
    const int capRings = 16;
    for (int i = 1; i < capRings; ++i) {  // the bottom cap, pole excluded
        const double a = -kPi / 2 + kPi / 2 * i / capRings;
        ringZ.push_back(r * std::sin(a));
        ringR.push_back(r * std::cos(a));
    }
    const int bodyRings = static_cast<int>(length / 1.5);
    for (int i = 0; i <= bodyRings; ++i) {
        ringZ.push_back(length * i / bodyRings);
        ringR.push_back(r);
    }
    for (int i = 1; i < capRings; ++i) {
        const double a = kPi / 2 * i / capRings;
        ringZ.push_back(length + r * std::sin(a));
        ringR.push_back(r * std::cos(a));
    }
    m.pos.push_back(Vec3{0, 0, -r});
    for (size_t j = 0; j < ringZ.size(); ++j)
        for (int i = 0; i < segments; ++i) {
            const double a = 2 * kPi * i / segments;
            m.pos.push_back(Vec3{ringR[j] * std::cos(a), ringR[j] * std::sin(a), ringZ[j]});
        }
    m.pos.push_back(Vec3{0, 0, length + r});
    const int rings = static_cast<int>(ringZ.size()), top = static_cast<int>(m.pos.size()) - 1;
    for (int i = 0; i < segments; ++i) m.tri.push_back({0, capsuleVertex(segments, 0, i + 1), capsuleVertex(segments, 0, i)});
    for (int j = 0; j + 1 < rings; ++j)
        for (int i = 0; i < segments; ++i) {
            const int a = capsuleVertex(segments, j, i), b = capsuleVertex(segments, j, i + 1);
            const int c = capsuleVertex(segments, j + 1, i), d = capsuleVertex(segments, j + 1, i + 1);
            m.tri.push_back({a, b, d});
            m.tri.push_back({a, d, c});
        }
    for (int i = 0; i < segments; ++i) m.tri.push_back({capsuleVertex(segments, rings - 1, i), capsuleVertex(segments, rings - 1, i + 1), top});
    m.computeNormals();
    return m;
}

// A closed marking round z at height z (n control points) on m's triangles.
SurfaceMarking ringMarking(const TriMesh& m, int id, MarkingMeaning meaning, double r, double z, int n) {
    SurfaceMarking k;
    k.id = id;
    k.meaning = meaning;
    k.shape = MarkingShape::ClosedCurve;
    for (int i = 0; i < n; ++i) {
        const double a = 2 * kPi * i / n + 0.031;
        k.points.push_back(closestSurfacePoint(m, Vec3{r * std::cos(a), r * std::sin(a), z}));
    }
    return k;
}

SurfaceMarking pathMarking(const TriMesh& m, int id, double r, double phi, double zFrom, double zTo) {
    SurfaceMarking k;
    k.id = id;
    k.meaning = MarkingMeaning::RowStartPath;
    k.shape = MarkingShape::Path;
    k.points.push_back(closestSurfacePoint(m, Vec3{r * std::cos(phi), r * std::sin(phi), zFrom}));
    k.points.push_back(closestSurfacePoint(m, Vec3{r * std::cos(phi), r * std::sin(phi), zTo}));
    return k;
}

bool hasFinding(const MarkingReview& r, FindingLevel level, FindingCause cause) {
    for (const MarkingFinding& f : r.findings)
        if (f.level == level && f.cause == cause) return true;
    return false;
}

double meshArea(const TriMesh& m) {
    double a = 0;
    for (const std::array<int, 3>& c : m.tri) a += 0.5 * norm(cross(m.pos[c[1]] - m.pos[c[0]], m.pos[c[2]] - m.pos[c[0]]));
    return a;
}

} // namespace

// Two closed curves round the tube's cylinder part split it into three regions; the band touches both curves, is a tube
// (Euler characteristic 0) and has the cylinder's area between the heights. Control: an open path cannot trim.
TEST(ssm_trim_tube_plane_curve) {
    const RevolutionMesh rm = markingTube();
    const TriangleNeighbours nb = triangleNeighbours(rm.mesh);
    const SurfaceMarking low = ringMarking(rm.mesh, 1, MarkingMeaning::Start, 30, 20, 6);
    const SurfaceMarking high = ringMarking(rm.mesh, 2, MarkingMeaning::Finish, 30, 50, 6);
    const MarkingTrace a = traceMarking(rm.mesh, nb, low.points, true), b = traceMarking(rm.mesh, nb, high.points, true);
    CHECK(!a.edgePathUsed && !b.edgePathUsed);
    TrimmedSurface s;
    std::string err;
    CHECK(splitAlongCurves(rm.mesh, {a, b}, s, &err));
    CHECK(s.regions == 3);
    int band = -1;
    for (int g = 0; g < s.regions; ++g)
        if (s.touches[g][0] && s.touches[g][1]) band = g;
    CHECK(band >= 0);
    const RegionMesh kept = regionMesh(s, band);
    const double exact = 2 * kPi * 30 * 30;
    std::printf("  band area %.3f mm2, cylinder %.3f\n", meshArea(kept.mesh), exact);
    CHECK(std::fabs(meshArea(kept.mesh) - exact) < 0.005 * exact);
    std::set<std::pair<int, int>> edges;
    for (const std::array<int, 3>& t : kept.mesh.tri)
        for (int k = 0; k < 3; ++k) edges.insert({std::min(t[k], t[(k + 1) % 3]), std::max(t[k], t[(k + 1) % 3])});
    CHECK(static_cast<int>(kept.mesh.pos.size()) - static_cast<int>(edges.size()) + static_cast<int>(kept.mesh.tri.size()) == 0);
    CHECK(kept.loops[0].size() > 100 && kept.loops[1].size() > 100);
    const MarkingTrace open = traceMarking(rm.mesh, nb, low.points, false);
    CHECK(!splitAlongCurves(rm.mesh, {open}, s, &err));
}

// The drawn row-start path along the meridian half-plane gives Method A exactly the half-plane's product (case D): the
// path guide plays the same role, it is not a new construction.
TEST(ssm_path_guide_equals_half_plane) {
    const RevolutionMesh rm = markingTube();
    const std::vector<double> t = fieldFromBottom(rm);
    const double phi = kPi / 384;
    GuidedRowsRequest half;
    half.w = 5.5;
    half.h = 3.4;
    half.cutNormal = Vec3{-std::sin(phi), std::cos(phi), 0};
    half.cutSide = Vec3{std::cos(phi), std::sin(phi), 0};
    GuidedRowsProduct a, b;
    std::string err;
    CHECK(buildGuidedRows(rm.mesh, t, rm.bottomLoop, rm.topLoop, half, a, &err));
    const SurfacePathGuide g = meridianGuide(rm.mesh, phi, -1, 1000);
    GuidedRowsRequest drawn = half;
    drawn.path = &g;
    CHECK(buildGuidedRows(rm.mesh, t, rm.bottomLoop, rm.topLoop, drawn, b, &err));
    CHECK(a.status == Status::Pass && b.status == Status::Pass && a.rows == b.rows && a.rows == 44);
    CHECK(a.mesh.mesh.faces() == b.mesh.mesh.faces() && a.mesh.mesh.pos.size() == b.mesh.mesh.pos.size());
    double worst = 0;
    for (size_t v = 0; v < a.mesh.mesh.pos.size() && v < b.mesh.mesh.pos.size(); ++v) worst = std::max(worst, norm(a.mesh.mesh.pos[v] - b.mesh.mesh.pos[v]));
    std::printf("  largest vertex difference %.3g mm\n", worst);
    CHECK(worst < 1e-9);
    for (int f = 0; f < a.mesh.mesh.faces() && f < b.mesh.mesh.faces(); ++f) CHECK(a.mesh.kind[f] == b.mesh.kind[f]);
    for (const RingSample& r : b.rings) CHECK(r.guideCrossings == 1);
}

// A path that turns back (up the meridian at 0.1 rad to 80 mm, down at 0.3 rad to 60 mm, up again at 0.5 rad) crosses
// the rings between 60 and 80 mm three times: those rows have no origin, the product is FAIL with no mesh, and the log
// names the row-start path. Control: the straight path of the previous case passes.
TEST(ssm_path_guide_double_crossing) {
    const RevolutionMesh rm = markingTube();
    const std::vector<double> t = fieldFromBottom(rm);
    const SurfacePathGuide up = meridianGuide(rm.mesh, 0.1, -1, 80), down = meridianGuide(rm.mesh, 0.3, 60, 80), again = meridianGuide(rm.mesh, 0.5, 60, 1000);
    SurfacePathGuide g;
    g.points = up.points;
    g.segmentTri = up.segmentTri;  // its last entry (-1): the jump to the next meridian runs in no triangle
    for (size_t i = down.points.size(); i-- > 0;) {
        g.points.push_back(down.points[i]);
        g.segmentTri.push_back(i > 0 ? down.segmentTri[i - 1] : -1);
    }
    g.points.insert(g.points.end(), again.points.begin(), again.points.end());
    g.segmentTri.insert(g.segmentTri.end(), again.segmentTri.begin(), again.segmentTri.end());
    GuidedRowsRequest r;
    r.w = 5.5;
    r.h = 3.4;
    r.path = &g;
    GuidedRowsProduct p;
    std::string err;
    CHECK(buildGuidedRows(rm.mesh, t, rm.bottomLoop, rm.topLoop, r, p, &err));
    CHECK(p.status == Status::Fail && p.mesh.mesh.faces() == 0 && p.contourEvents > 0);
    int three = 0, once = 0;
    for (const RingSample& rs : p.rings) {
        if (rs.guideCrossings == 3) ++three;
        if (rs.guideCrossings == 1) ++once;
    }
    std::printf("  rings crossed three times %d, once %d, of %zu\n", three, once, p.rings.size());
    CHECK(three > 0 && once > 0 && three + once == static_cast<int>(p.rings.size()));
    bool named = false;
    for (const std::string& l : p.log) named = named || l.find("row-start path 3 times") != std::string::npos;
    CHECK(named);
}

// Method A on a marked closed capsule: a start curve, a finish curve and a row-start path give a relaxed stitch mesh of
// the band between the curves, the caps omitted on purpose (their area reported), the path crossing every row once.
// Controls: a start point needs the pole start the mesher lacks (an algorithm limitation, no mesh); without the
// row-start path there is no mesh and the finding names the missing marking.
TEST(ssm_marked_capsule_build) {
    const TriMesh m = capsule(30, 90, 96);
    const std::vector<int> oneNode(m.tri.size(), 0);
    const SurfacePieces pieces = findSurfacePieces(m.pos, m.tri, oneNode, 1e-9);
    CHECK(pieces.pieces.size() == 1 && pieces.pieces[0].mesh.tri.size() == m.tri.size());
    MarkingSet set;
    set.markings = {ringMarking(m, 1, MarkingMeaning::Start, 30, 15, 6), ringMarking(m, 2, MarkingMeaning::Finish, 30, 75, 6),
                    pathMarking(m, 3, 30, 0.013, 5, 85)};
    MarkedBuildRequest rq;
    rq.pieces = &pieces;
    rq.piece = 0;
    rq.markings = &set;
    MarkedBuildProduct p;
    std::string err;
    CHECK(buildFromMarkings(rq, p, &err));
    for (const MarkingFinding& f : p.review.findings)
        std::printf("  finding on %d, level %d, cause %d: %s\n", f.markingId, static_cast<int>(f.level), static_cast<int>(f.cause), f.message.c_str());
    for (const std::string& l : p.log) std::printf("  %s\n", l.c_str());
    const LayoutMeasures& ms = p.review.measures;
    std::printf("  rows %d, cells %d, exceptional %d, covered %.4f, omitted %.1f mm2, p95 %.4f, %.3f s\n", ms.rows, p.stitchMesh.mesh.faces(), ms.exceptionalCells,
                ms.coveredFraction, ms.omittedAreaMm2, ms.distortionP95, ms.seconds);
    CHECK(p.review.valid && p.stitchMesh.mesh.faces() > 0);
    CHECK(ms.rows == 11 && ms.contourEvents == 0 && ms.unsupportedCells == 0);
    CHECK(hasFinding(p.review, FindingLevel::Ok, FindingCause::Marking));
    CHECK(!hasFinding(p.review, FindingLevel::Problem, FindingCause::Marking));
    CHECK(ms.coveredFraction > 0.99);
    const double caps = 2 * (2 * kPi * 30 * 15 + 2 * kPi * 30 * 30);
    CHECK(std::fabs(ms.omittedAreaMm2 - caps) < 0.02 * caps);
    int kept = 0, start = 0, finish = 0;
    for (MarkedRegion g : p.region) {
        kept += g == MarkedRegion::Kept ? 1 : 0;
        start += g == MarkedRegion::OmittedStart ? 1 : 0;
        finish += g == MarkedRegion::OmittedFinish ? 1 : 0;
    }
    CHECK(kept > 0 && start > 0 && finish > 0);
    // a start point takes the (experimental) pole route: a gathered pole start, the finish curve as the bind-off loop
    MarkingSet point = set;
    point.markings[0].shape = MarkingShape::Point;
    point.markings[0].points.resize(1);
    rq.markings = &point;
    CHECK(buildFromMarkings(rq, p, &err));
    CHECK(p.evaluation.route == "pole");
    CHECK(hasFinding(p.review, FindingLevel::Warning, FindingCause::AlgorithmLimitation));   // flagged experimental
    // control: a finish point after a start curve is an unsupported request, not a marking problem
    MarkingSet finishPoint = set;
    finishPoint.markings[1].shape = MarkingShape::Point;
    finishPoint.markings[1].points.resize(1);
    rq.markings = &finishPoint;
    CHECK(buildFromMarkings(rq, p, &err));
    CHECK(!p.review.valid && hasFinding(p.review, FindingLevel::Problem, FindingCause::AlgorithmLimitation));
    // control: no row-start path
    MarkingSet noPath = set;
    noPath.markings.pop_back();
    rq.markings = &noPath;
    CHECK(buildFromMarkings(rq, p, &err));
    CHECK(!p.review.valid && hasFinding(p.review, FindingLevel::Problem, FindingCause::Marking));
}

namespace {

// The lowest and highest centroid height of the piece triangles a build kept.
void keptHeights(const TriMesh& piece, const MarkedBuildProduct& p, double& lo, double& hi) {
    lo = 1e30;
    hi = -1e30;
    for (size_t t = 0; t < p.region.size(); ++t) {
        if (p.region[t] != MarkedRegion::Kept) continue;
        const double z = (piece.pos[piece.tri[t][0]].z + piece.pos[piece.tri[t][1]].z + piece.pos[piece.tri[t][2]].z) / 3;
        lo = std::min(lo, z);
        hi = std::max(hi, z);
    }
}

// The mean distance from the z axis of a reference's vertices on the capsule's cylinder part (10 < z < 80).
double meanCylinderRadius(const TriMesh& reference) {
    double sum = 0;
    int n = 0;
    for (const Vec3& x : reference.pos)
        if (x.z > 10 && x.z < 80) {
            sum += std::hypot(x.x, x.y);
            ++n;
        }
    return n > 0 ? sum / n : 0;
}

} // namespace

// One piece knitted as two overlapping knitting components, as a collar over a body: the body (component 0) from z 15 to
// 55 and the collar (component 1) from z 45 to 75, each with its own start, finish and row-start path. Each build takes only
// its own component's markings and keeps only its own band (45 to 55 lies in both), and the collar is laid out one yarn
// layer (D + c = 1.264 mm) outside the piece. (Defect: markings selected by piece alone, so each build sees two starts, two
// finishes and two paths and neither is valid.)
TEST(ssm_knit_components_overlap_on_one_piece) {
    const TriMesh m = capsule(30, 90, 96);
    const std::vector<int> oneNode(m.tri.size(), 0);
    const SurfacePieces pieces = findSurfacePieces(m.pos, m.tri, oneNode, 1e-9);
    MarkingSet set;
    set.markings = {ringMarking(m, 1, MarkingMeaning::Start, 30, 15, 6), ringMarking(m, 2, MarkingMeaning::Finish, 30, 55, 6),
                    pathMarking(m, 3, 30, 0.013, 5, 65),
                    ringMarking(m, 4, MarkingMeaning::Start, 30, 45, 6), ringMarking(m, 5, MarkingMeaning::Finish, 30, 75, 6),
                    pathMarking(m, 6, 30, 0.013, 35, 85)};
    for (size_t i = 3; i < set.markings.size(); ++i) set.markings[i].knitComponent = 1;
    set.nextId = 7;
    MarkedBuildRequest rq;
    rq.pieces = &pieces;
    rq.piece = 0;
    rq.markings = &set;
    MarkedBuildProduct body, collar;
    std::string err;
    CHECK(buildFromMarkings(rq, body, &err));
    rq.knitComponent = 1;
    rq.layerOffsetMm = 1.264;
    CHECK(buildFromMarkings(rq, collar, &err));
    for (const MarkedBuildProduct* p : {&body, &collar})
        for (const MarkingFinding& f : p->review.findings)
            if (f.level == FindingLevel::Problem) std::printf("  problem on %d: %s\n", f.markingId, f.message.c_str());
    CHECK(body.review.valid && collar.review.valid);
    CHECK(body.markingIds == (std::vector<int>{1, 2, 3}) && collar.markingIds == (std::vector<int>{4, 5, 6}));
    double bodyLo = 0, bodyHi = 0, collarLo = 0, collarHi = 0;
    keptHeights(pieces.pieces[0].mesh, body, bodyLo, bodyHi);
    keptHeights(pieces.pieces[0].mesh, collar, collarLo, collarHi);
    std::printf("  body kept z %.1f .. %.1f, collar kept z %.1f .. %.1f; rows %d and %d\n", bodyLo, bodyHi, collarLo, collarHi, body.review.measures.rows,
                collar.review.measures.rows);
    CHECK(bodyLo > 13 && bodyHi < 57 && collarLo > 43 && collarHi < 77);
    CHECK(bodyHi > collarLo + 5);   // the overlap
    const double rBody = meanCylinderRadius(body.evaluation.reference), rCollar = meanCylinderRadius(collar.evaluation.reference);
    std::printf("  reference radius: body %.4f, collar %.4f mm\n", rBody, rCollar);
    CHECK(std::fabs(rBody - 30) < 0.05 && std::fabs(rCollar - rBody - 1.264) < 0.01);
}

// ---------------------------------------------------------------------------------------------------------------------
// The integration of the Mark step with the imported-arm experiment.

namespace {

// A start point near the capsule's lower pole, a finish point on the cylinder at z, a row-start path over the cap and up
// a meridian.
MarkingSet poleCapsuleSet(const TriMesh& m, double r, double z) {
    MarkingSet set;
    SurfaceMarking s;
    s.id = 1;
    s.meaning = MarkingMeaning::Start;
    s.shape = MarkingShape::Point;
    s.points.push_back(closestSurfacePoint(m, Vec3{0.31, 0.17, -r}));
    SurfaceMarking f;
    f.id = 2;
    f.meaning = MarkingMeaning::Finish;
    f.shape = MarkingShape::Point;
    f.points.push_back(closestSurfacePoint(m, Vec3{r * std::cos(0.013), r * std::sin(0.013), z}));
    SurfaceMarking p;
    p.id = 3;
    p.meaning = MarkingMeaning::RowStartPath;
    p.shape = MarkingShape::Path;
    p.points.push_back(s.points[0]);
    for (int i = 1; i <= 8; ++i) {
        const double a = -kPi / 2 + kPi / 2 * i / 8;
        p.points.push_back(closestSurfacePoint(m, Vec3{r * std::cos(a) * std::cos(0.013), r * std::cos(a) * std::sin(0.013), r * std::sin(a)}));
    }
    p.points.push_back(closestSurfacePoint(m, Vec3{r * std::cos(0.013), r * std::sin(0.013), z + 10}));
    set.markings = {s, f, p};
    set.nextId = 4;
    return set;
}

} // namespace

// A start point takes the pole route; a finish point asks for the bind-off at the field's level through it. The layout is
// typed throughout (gathered pole start, the bind-off strip's decreases admitted by the boundary-shaping variant), the
// virtual shift is the closed pole round plus one helix, the area above the level is omitted (not uncovered), and the
// relaxation holds the pole and the row origins (the markings' constraints).
TEST(ssm_marked_pole_route_finish_level) {
    const TriMesh m = capsule(30, 90, 96);
    const std::vector<int> oneNode(m.tri.size(), 0);
    const SurfacePieces pieces = findSurfacePieces(m.pos, m.tri, oneNode, 1e-9);
    const MarkingSet set = poleCapsuleSet(m, 30, 50);
    MarkedBuildRequest rq;
    rq.pieces = &pieces;
    rq.piece = 0;
    rq.markings = &set;
    rq.w = 5.5;
    rq.h = 3.4;
    MarkedBuildProduct p;
    std::string err;
    CHECK(buildFromMarkings(rq, p, &err));
    for (const MarkingFinding& f : p.review.findings) std::printf("  finding %d: %s\n", f.markingId, f.message.c_str());
    CHECK(p.review.valid && p.evaluation.route == "pole");
    CHECK(p.review.measures.unsupportedCells == 0 && p.review.measures.contourEvents == 0);
    CHECK(p.evaluation.shiftBefore.status == Status::Pass && p.evaluation.shiftBefore.poleStart);
    CHECK(p.evaluation.storedBefore.connectivity == Status::Pass);
    const double expectedLevel = 30 * kPi / 2 + 50;   // the geodesic distance from the pole to z = 50
    CHECK(std::fabs(p.evaluation.finishLevel - expectedLevel) < 1.0);
    CHECK(p.review.measures.omittedAreaMm2 > 0.9 * (2 * kPi * 30 * 40 + 2 * kPi * 30 * 30));   // above z = 50: cylinder rest + top cap
    // the defaults (RelaxAnchors::None, BorderRings::PlacedFirst): the bind-off ring is placed on its loop with its first
    // vertex within half a stitch of where the layout put it and held; the rows and the pole are free (the pole stays
    // within a third of a stitch here), and no cell is turned over
    std::printf("  pole moved %.4f mm\n", p.evaluation.poleMoved);
    CHECK(p.evaluation.relaxed && p.evaluation.poleMoved < 0.3 * rq.w);   // measured 0.99 mm: the start point is drawn 0.35 mm off the axis
    const int yarnEnd = p.evaluation.after.rings.back()[0];
    CHECK(norm(p.evaluation.after.mesh.pos[yarnEnd] - p.evaluation.before.mesh.pos[yarnEnd]) < 0.5 * rq.w + 0.26);
    CHECK(p.evaluation.storedAfter.invertedFaces == 0);
    // the first integration's anchors hold every row's first vertex on the row-start path instead
    rq.options.anchors = RelaxAnchors::RowOriginsOnPath;
    MarkedBuildProduct onPath;
    CHECK(buildFromMarkings(rq, onPath, &err));
    CHECK(onPath.evaluation.relaxed && onPath.evaluation.poleMoved < 1e-9 && onPath.evaluation.rowOriginMean < 0.2);
}

// The yarn generated on a marked piece with a pole start and a closed bind-off ring: every face crossed as its
// construction requires, and the binding thread through every cap exactly once. On a closed ring the first cap follows the
// last one, so the thread's end must stop short of it.
TEST(ssm_marked_pole_route_yarn_is_clean) {
    for (double radius : {30.0, 11.0}) {   // a wide ring, and one of 13 stitches (the teddy arm's), where the caps sit close
    const TriMesh m = capsule(radius, 3 * radius, 96);
    const std::vector<int> oneNode(m.tri.size(), 0);
    const SurfacePieces pieces = findSurfacePieces(m.pos, m.tri, oneNode, 1e-9);
    const MarkingSet set = poleCapsuleSet(m, radius, 1.6 * radius);
    MarkedBuildRequest rq;
    rq.pieces = &pieces;
    rq.piece = 0;
    rq.markings = &set;
    rq.w = 5.5;
    rq.h = 3.4;
    MarkedBuildProduct p;
    std::string err;
    CHECK(buildFromMarkings(rq, p, &err));
    const MarkedEvaluation& ev = p.evaluation;
    YarnEmbedding e;
    const YarnSettings st;
    CHECK(embedYarn(ev.after, ev.shift.corners, ev.reference, rq.w, rq.h, st, e, &err));
    const YarnCheck c = checkYarn(ev.after, ev.shift.corners, ev.reference, e, st);
    std::printf("  %s", c.text().c_str());
    for (size_t i = 0; i < c.defects.size() && i < 6; ++i) std::printf("  defect: %s\n", c.defects[i].c_str());
    std::printf("  radius %.0f mm: %d caps\n", radius, e.caps);
    CHECK(e.caps > 0);
    CHECK(c.status == Status::Pass && c.facesCrossingsOk == c.facesChecked);
    }
}

// Only user markings and accepted proposals are input. Offered proposals are input only when named for a research
// preview (labelled in the review); rejected proposals never, not even when named.
TEST(ssm_marked_offered_preview_rejected) {
    const TriMesh m = capsule(30, 90, 96);
    const std::vector<int> oneNode(m.tri.size(), 0);
    const SurfacePieces pieces = findSurfacePieces(m.pos, m.tri, oneNode, 1e-9);
    MarkingSet set = poleCapsuleSet(m, 30, 50);
    for (SurfaceMarking& k : set.markings) {
        k.origin = MarkingOrigin::Proposal;
        k.proposal = ProposalState::Offered;
        k.proposalId = k.id;
    }
    MarkedBuildRequest rq;
    rq.pieces = &pieces;
    rq.piece = 0;
    rq.markings = &set;
    rq.h = 3.4;
    rq.depth = ReviewDepth::Built;
    MarkedBuildProduct p;
    std::string err;
    CHECK(buildFromMarkings(rq, p, &err));
    CHECK(!p.review.valid && p.markingIds.empty() && hasFinding(p.review, FindingLevel::Ok, FindingCause::Input));
    rq.previewOffered = {1, 2, 3};
    CHECK(buildFromMarkings(rq, p, &err));
    CHECK(p.review.valid && p.previewIds.size() == 3 && hasFinding(p.review, FindingLevel::Warning, FindingCause::Input));
    set.markings[2].proposal = ProposalState::Rejected;   // the path rejected, still named for the preview
    CHECK(buildFromMarkings(rq, p, &err));
    CHECK(!p.review.valid && std::find(p.markingIds.begin(), p.markingIds.end(), 3) == p.markingIds.end());
    for (SurfaceMarking& k : set.markings) k.proposal = ProposalState::Accepted;   // accepted: input without a preview
    rq.previewOffered.clear();
    CHECK(buildFromMarkings(rq, p, &err));
    CHECK(p.review.valid && p.previewIds.empty());
}

// The prepared reference: midpoint refinement keeps the surface, maps every point there and back exactly, and a marking
// traced on the refined mesh is the same curve as on the piece (the plane section does not depend on the refinement).
TEST(ssm_prepared_reference_maps_exactly) {
    const TriMesh m = capsule(30, 90, 48);
    const RefinedMesh r = refineMidpoint(m, 2);
    CHECK(r.mesh.tri.size() == 16 * m.tri.size());
    double worst = 0;
    for (int t = 0; t < static_cast<int>(m.tri.size()); t += 7) {
        const SurfacePoint p{t, {0.2, 0.3, 0.5}};
        const SurfacePoint q = toRefined(r, p);
        const SurfacePoint back = fromRefined(r, q);
        worst = std::max(worst, norm(surfacePosition(m, p) - surfacePosition(r.mesh, q)));
        CHECK(back.tri == t);
        for (int k = 0; k < 3; ++k) worst = std::max(worst, std::fabs(back.bary[k] - p.bary[k]));
    }
    CHECK(worst < 1e-12);
    const SurfaceMarking ring = ringMarking(m, 1, MarkingMeaning::Finish, 30, 45, 6);
    std::vector<SurfacePoint> refinedPoints;
    for (const SurfacePoint& p : ring.points) refinedPoints.push_back(toRefined(r, p));
    const MarkingTrace a = traceMarking(m, triangleNeighbours(m), ring.points, true);
    const MarkingTrace b = traceMarking(r.mesh, triangleNeighbours(r.mesh), refinedPoints, true);
    CHECK(!a.edgePathUsed && !b.edgePathUsed && b.points.size() > a.points.size());
    double far = 0;   // every point of the refined trace lies on the piece trace's polyline
    for (const TracePoint& x : b.points) {
        double d = 1e300;
        for (size_t i = 0; i < a.points.size(); ++i)
            d = std::min(d, closestSegmentPoints(x.pos, x.pos, a.points[i].pos, a.points[(i + 1) % a.points.size()].pos).distance);
        far = std::max(far, d);
    }
    CHECK(far < 1e-9);
}

// The faceted stationarity measure on a minimal crease (two triangles meeting along y, each sloping down by 15 degrees):
// a vertex on the crease pushed into both facets is at rest (0); pushed along the crease, or pulled into one facet, it is
// not; inside a facet it is the force's tangential part.
TEST(ssm_walk_residual_on_a_crease) {
    TriMesh roof;
    const double s = std::sin(kPi / 12), c = std::cos(kPi / 12);
    roof.pos = {Vec3{0, -1, 0}, Vec3{0, 1, 0}, Vec3{-c, 0, -s}, Vec3{c, 0, -s}};
    roof.tri = {{0, 1, 2}, {1, 0, 3}};
    roof.computeNormals();
    const std::vector<Vec3> pos{Vec3{0, 0, 0}};
    const std::vector<char> border{0};
    auto residual = [&](const Vec3& f) { return walkResidualOf(roof, pos, {f}, border, {}, nullptr, 1e-6).max; };
    CHECK(residual(Vec3{0, 0, 1}) < 1e-12);                     // into the ridge from both sides
    CHECK(std::fabs(residual(Vec3{0, 0.5, 0}) - 0.5) < 1e-12);   // along the crease
    CHECK(residual(Vec3{1, 0, 0}) > 0.9);                       // pulled into one facet
    const std::vector<Vec3> inside{Vec3{0.3, 0.1, -0.3 * s / c}};
    CHECK(std::fabs(walkResidualOf(roof, inside, {Vec3{0, 1, 0}}, border, {}, nullptr, 1e-6).max - 1) < 1e-9);
}

// The smooth-reference preparation (prepared_reference.h, smooth-reference report): Loop subdivision on the same split as
// the midpoint refinement (identical parents and corner coordinates, so a file point maps to the same refined triangle and
// weights); each original vertex's limit point the same whatever level the limit is taken from (Loop's interior rule and
// its limit mask, the cubic B-spline border rule and its limit mask, are consistent); the border stays in its plane.
// Control: the midpoint refinement does not move the surface, so its vertices are not the limit points.
TEST(ssm_loop_reference_limit_points) {
    // an open prism with mixed valences (alternating diagonals) and two border rings, and a closed octahedron
    TriMesh prism;
    const int sides = 7, rings = 3;
    for (int k = 0; k <= rings; ++k)
        for (int j = 0; j < sides; ++j) prism.pos.push_back({10 * std::cos(2 * kPi * j / sides), 10 * std::sin(2 * kPi * j / sides), 6.0 * k});
    auto id = [&](int k, int j) { return k * sides + ((j % sides) + sides) % sides; };
    for (int k = 0; k < rings; ++k)
        for (int j = 0; j < sides; ++j) {
            if ((k + j) % 2 == 0) {
                prism.tri.push_back({id(k, j), id(k, j + 1), id(k + 1, j + 1)});
                prism.tri.push_back({id(k, j), id(k + 1, j + 1), id(k + 1, j)});
            } else {
                prism.tri.push_back({id(k, j), id(k, j + 1), id(k + 1, j)});
                prism.tri.push_back({id(k, j + 1), id(k + 1, j + 1), id(k + 1, j)});
            }
        }
    TriMesh octa;
    octa.pos = {Vec3{10, 0, 0}, Vec3{-10, 0, 0}, Vec3{0, 10, 0}, Vec3{0, -10, 0}, Vec3{0, 0, 10}, Vec3{0, 0, -10}};
    octa.tri = {{0, 2, 4}, {2, 1, 4}, {1, 3, 4}, {3, 0, 4}, {2, 0, 5}, {1, 2, 5}, {3, 1, 5}, {0, 3, 5}};
    for (const TriMesh* m : {&prism, &octa}) {
        const RefinedMesh mid = refineMidpoint(*m, 2), l0 = refineLoop(*m, 0), l1 = refineLoop(*m, 1), l2 = refineLoop(*m, 2), l3 = refineLoop(*m, 3);
        CHECK(mid.parent == l2.parent && mid.mesh.tri == l2.mesh.tri);
        bool sameBary = true;
        for (size_t t = 0; t < mid.cornerBary.size(); ++t) sameBary = sameBary && mid.cornerBary[t] == l2.cornerBary[t];
        CHECK(sameBary);
        double drift = 0, midDrift = 0;
        for (size_t v = 0; v < m->pos.size(); ++v) {
            drift = std::max(drift, std::max(norm(l1.mesh.pos[v] - l0.mesh.pos[v]), std::max(norm(l2.mesh.pos[v] - l0.mesh.pos[v]), norm(l3.mesh.pos[v] - l0.mesh.pos[v]))));
            midDrift = std::max(midDrift, norm(mid.mesh.pos[v] - l0.mesh.pos[v]));
        }
        std::printf("  limit points: drift over levels %.3g mm (midpoint control %.3g mm)\n", drift, midDrift);
        CHECK(drift < 1e-9);
        CHECK(midDrift > 0.1);
        // a point of the input keeps its refined triangle and weights between the preparations
        const SurfacePoint p{1, {0.2, 0.3, 0.5}};
        const SurfacePoint q = toRefined(l2, p), r = mapBetween(mid, l2, toRefined(mid, p));
        CHECK(q.tri == r.tri && std::fabs(q.bary[0] - r.bary[0]) < 1e-12 && std::fabs(q.bary[1] - r.bary[1]) < 1e-12);
        // across levels the map goes through the piece: the same input point
        const SurfacePoint s3 = mapBetween(l2, l3, q), back = fromRefined(l3, s3);
        CHECK(back.tri == 1 && std::fabs(back.bary[0] - 0.2) < 1e-9 && std::fabs(back.bary[2] - 0.5) < 1e-9);
        // every level-1 vertex is a limit point already: the same position at levels 1, 2 and 3
        double d12 = 0, d23 = 0;
        for (size_t v = 0; v < l1.mesh.pos.size(); ++v) d12 = std::max(d12, norm(l1.mesh.pos[v] - l2.mesh.pos[v]));
        for (size_t v = 0; v < l2.mesh.pos.size(); ++v) d23 = std::max(d23, norm(l2.mesh.pos[v] - l3.mesh.pos[v]));
        std::printf("  shared vertices across levels: 1-2 %.3g, 2-3 %.3g mm\n", d12, d23);
        CHECK(d12 < 1e-9 && d23 < 1e-9);
    }
    const RefinedMesh lp = refineLoop(prism, 2);
    for (const std::vector<int>& loop : borderLoops(lp.mesh)) {
        double zlo = 1e300, zhi = -1e300;
        for (int v : loop) {
            zlo = std::min(zlo, lp.mesh.pos[v].z);
            zhi = std::max(zhi, lp.mesh.pos[v].z);
        }
        CHECK(zhi - zlo < 1e-9);
    }
}

// The marked build on the smooth-reference variant and the controlled comparison: the review names the preparation and its
// distance from the imported surface; the markings keep their file points (reported shift); the controlled comparison
// keeps the facet layout's connectivity exactly and relaxes it on the Loop reference. Control: the facets build reports
// no change of surface.
TEST(ssm_marked_loop_preparation) {
    const TriMesh m = capsule(30, 90, 24);
    const std::vector<int> oneNode(m.tri.size(), 0);
    const SurfacePieces pieces = findSurfacePieces(m.pos, m.tri, oneNode, 1e-9);
    const MarkingSet set = poleCapsuleSet(m, 30, 50);
    MarkedBuildRequest rq;
    rq.pieces = &pieces;
    rq.piece = 0;
    rq.markings = &set;
    rq.w = 5.5;
    rq.h = 3.4;
    MarkedBuildProduct facets, loop, controlled;
    std::string err;
    CHECK(buildFromMarkings(rq, facets, &err));
    rq.options.preparation = ReferencePreparation::LoopLimit;
    CHECK(buildFromMarkings(rq, loop, &err));
    rq.options.preparation = ReferencePreparation::Facets;
    rq.options.relaxOnOther = true;
    CHECK(buildFromMarkings(rq, controlled, &err));
    CHECK(facets.review.valid && loop.review.valid && controlled.review.valid);
    CHECK(facets.evaluation.preparation == "facets" && facets.evaluation.changeMax == 0 && facets.evaluation.markingShift.empty());
    std::printf("  loop: change %.3f mm (back %.3f), markings moved up to %.3f mm, faces %d (facets %d)\n", loop.evaluation.changeMax, loop.evaluation.changeBackMax,
                loop.evaluation.markingShift.empty() ? -1.0 : loop.evaluation.markingShift[0].controlMax, loop.stitchMesh.mesh.faces(), facets.stitchMesh.mesh.faces());
    CHECK(loop.evaluation.preparation == "loop_limit" && loop.evaluation.changeMax > 0.01 && loop.evaluation.changeMax < 2.0);
    CHECK(!loop.evaluation.markingShift.empty() && loop.evaluation.referenceDigest != facets.evaluation.referenceDigest);
    bool named = false;
    for (const MarkingFinding& f : loop.review.findings) named = named || f.message.find("Smooth reference (Loop limit)") != std::string::npos;
    CHECK(named);
    // controlled: the same faces, corners, roles and kinds as the facets layout; judged against the Loop reference
    const TypedStitchMesh& a = facets.evaluation.before;
    const TypedStitchMesh& b = controlled.evaluation.before;
    CHECK(a.mesh.faces() == b.mesh.faces() && a.mesh.pos.size() == b.mesh.pos.size());
    bool same = a.mesh.faces() == b.mesh.faces();
    for (int f = 0; same && f < a.mesh.faces(); ++f) {
        same = a.kind[f] == b.kind[f] && a.mesh.sides(f) == b.mesh.sides(f);
        for (int k = 0; same && k < a.mesh.sides(f); ++k) same = a.mesh.vertexAt(f, k) == b.mesh.vertexAt(f, k) && a.role[f][k] == b.role[f][k];
    }
    CHECK(same);
    CHECK(controlled.evaluation.relaxPreparation == "loop_limit" && controlled.evaluation.referenceDigest == loop.evaluation.referenceDigest);
    double moved = 0;
    for (size_t v = 0; v < a.mesh.pos.size(); ++v) moved = std::max(moved, norm(a.mesh.pos[v] - b.mesh.pos[v]));
    CHECK(moved > 0 && moved < 2.0);   // carried onto the Loop surface before relaxing
}

// The shaping-cell assessment (typed_mesh.h measureShapingCells): an isolated kp relaxed by the model is within it although
// the frozen pentagon bounds reject it (their trapezoid is not the model's rest shape); the same cell folded (its shared top
// vertex below the bottom edge) is rejected by both. The frozen bounds are kept and reported.
TEST(ssm_shaping_assessment_model_and_defect) {
    const double w = 5.5, h = 3.4;
    TypedStitchMesh m;
    TriMesh unusedRef;
    std::vector<DeclaredBoundary> decl;
    std::string err;
    CHECK(ringTubeFixture({24, 24, 25, 25}, {{}, {12}, {}}, 24.5 * w / (2 * kPi), h, m, unusedRef, decl, &err));
    const TriMesh ref = cylinderReference(24.5 * w / (2 * kPi), 0, 3 * h, h / 4);
    VirtualShift vs;
    CHECK(buildVirtualShift(m, ShiftReading::PathSideDown, vs, &err));
    std::vector<RelaxFace> faces;
    for (int f = 0; f < m.mesh.faces(); ++f) faces.push_back(RelaxFace{vs.corners[f], m.role[f]});
    std::vector<char> border(m.mesh.pos.size(), 0);
    for (const std::vector<int>& l : polygonBorderLoops(m.mesh))
        for (int v : l) border[v] = 1;
    MeshRelaxRequest r;
    r.w = w;
    r.h = h;
    RelaxProduct p;
    CHECK(relaxStitchMesh(ref, triangleBoundaryLoops(ref), faces, m.mesh.pos, border, r, p, &err));
    m.mesh.pos = p.pos;
    const std::vector<ShapingCellMeasure> sc = measureShapingCells(m, w, h);
    const std::vector<CellMeasure> old = measureCells(m, w, h);
    CHECK(sc.size() == 1 && sc[0].kind == StitchFaceKind::Increase && sc[0].subQuads == 2);
    if (sc.empty()) return;
    std::printf("  relaxed kp: course-course %.1f deg (frozen bounds %s), sub-quad corners %.1f-%.1f, wale/h %.3f-%.3f (assessment %s)\n", old[sc[0].face].courseCourseAngle,
                old[sc[0].face].withinBounds ? "within" : old[sc[0].face].outside.c_str(), sc[0].minCorner, sc[0].maxCorner, sc[0].minWale, sc[0].maxWale,
                sc[0].within ? "within" : sc[0].outside.c_str());
    CHECK(sc[0].within && sc[0].foldedSubQuads == 0);
    CHECK(!old[sc[0].face].withinBounds);   // the finding: the trapezoid bounds reject the model's own equilibrium
    const int f = sc[0].face, n = m.mesh.sides(f);
    int mid = -1, b0 = -1, b1 = -1;
    for (int k = 0; k < n; ++k) {
        if (m.role[f][(k + n - 1) % n] == StitchEdgeRole::Top && m.role[f][k] == StitchEdgeRole::Top) mid = m.mesh.vertexAt(f, k);
        if (m.role[f][k] == StitchEdgeRole::Bottom) {
            b0 = m.mesh.vertexAt(f, k);
            b1 = m.mesh.vertexAt(f, (k + 1) % n);
        }
    }
    CHECK(mid >= 0 && b0 >= 0);
    if (mid < 0 || b0 < 0) return;
    const Vec3 bm = 0.5 * (m.mesh.pos[b0] + m.mesh.pos[b1]);
    TypedStitchMesh folded = m;
    folded.mesh.pos[mid] = bm - (m.mesh.pos[mid] - bm) * 0.3;
    const std::vector<ShapingCellMeasure> bad = measureShapingCells(folded, w, h);
    CHECK(!bad[0].within && bad[0].foldedSubQuads > 0);
    CHECK(!measureCells(folded, w, h)[f].withinBounds);
}

// The rows of a marked piece are turned for upright wales (RowPhase::WalesUpright). With the first
// vertex of every ring on the row-start path, the wale edges of a shaping strip all lean one way; the mean lean of a strip
// is what the relaxation then has to turn the ring by. Turned, no strip keeps a mean lean.
namespace {

// Per strip, the mean of (top - bottom) of its wale edges along the direction of the top ring at the top vertex [mm]; the
// largest of their sizes.
double largestMeanLean(const TypedStitchMesh& s) {
    std::vector<int> ringOf(s.mesh.pos.size(), -1), indexOf(s.mesh.pos.size(), -1);
    for (size_t k = 0; k < s.rings.size(); ++k)
        for (size_t j = 0; j < s.rings[k].size(); ++j) {
            ringOf[s.rings[k][j]] = static_cast<int>(k);
            indexOf[s.rings[k][j]] = static_cast<int>(j);
        }
    std::vector<double> sum(s.rings.size(), 0);
    std::vector<int> count(s.rings.size(), 0);
    for (int f = 0; f < s.mesh.faces(); ++f)
        for (int k = 0; k < s.mesh.sides(f); ++k) {
            if (s.role[f][k] != StitchEdgeRole::Wale) continue;
            int a = s.mesh.vertexAt(f, k), b = s.mesh.vertexAt(f, (k + 1) % s.mesh.sides(f));
            if (ringOf[a] > ringOf[b]) std::swap(a, b);
            const std::vector<int>& ring = s.rings[ringOf[b]];
            if (ringOf[a] < 1 || ring.size() < 3) continue;
            const int n = static_cast<int>(ring.size()), j = indexOf[b];
            const Vec3 along = normalized(s.mesh.pos[ring[(j + 1) % n]] - s.mesh.pos[ring[(j + n - 1) % n]]);
            sum[ringOf[b]] += dot(s.mesh.pos[b] - s.mesh.pos[a], along);
            ++count[ringOf[b]];
        }
    double largest = 0;
    for (size_t k = 0; k < sum.size(); ++k)
        if (count[k] > 0) largest = std::max(largest, std::fabs(sum[k] / count[k]));
    return largest;
}

} // namespace

TEST(ssm_marked_rows_turned_for_upright_wales) {
    const TriMesh m = capsule(30, 90, 96);
    const std::vector<int> oneNode(m.tri.size(), 0);
    const SurfacePieces pieces = findSurfacePieces(m.pos, m.tri, oneNode, 1e-9);
    const MarkingSet set = poleCapsuleSet(m, 30, 50);
    MarkedBuildRequest rq;
    rq.pieces = &pieces;
    rq.piece = 0;
    rq.markings = &set;
    rq.w = 5.5;
    rq.h = 3.4;
    MarkedBuildProduct onPath, upright;
    std::string err;
    rq.options.rowPhase = RowPhase::OnPath;
    CHECK(buildFromMarkings(rq, onPath, &err));
    rq.options.rowPhase = RowPhase::WalesUpright;
    CHECK(buildFromMarkings(rq, upright, &err));
    const double leanOnPath = largestMeanLean(onPath.evaluation.before), leanUpright = largestMeanLean(upright.evaluation.before);
    std::printf("  largest mean lean of a strip: on the path %.3f mm, turned %.3f mm\n", leanOnPath, leanUpright);
    CHECK(leanOnPath > 1.0);
    CHECK(leanUpright < 0.25);
    CHECK(upright.evaluation.storedBefore.connectivity == Status::Pass && upright.evaluation.shiftBefore.status == Status::Pass);
    CHECK(upright.evaluation.storedAfter.invertedFaces == 0);
}

// The bind-off ring takes the stitch count of the round below (BindOffCount::OfRoundBelow). A
// capsule finished on its far cap narrows towards the finish. Counted from the loop's length, the bind-off ring has fewer
// stitches than the round below and the last strip holds decreases, each with its top edge on the border. With the count
// of the round below the last strip is plain, and the mesh and its helical reading are still valid.
TEST(ssm_marked_bind_off_strip_is_plain) {
    const TriMesh m = capsule(30, 90, 96);
    const std::vector<int> oneNode(m.tri.size(), 0);
    const SurfacePieces pieces = findSurfacePieces(m.pos, m.tri, oneNode, 1e-9);
    const MarkingSet set = poleCapsuleSet(m, 30, 108);
    MarkedBuildRequest rq;
    rq.pieces = &pieces;
    rq.piece = 0;
    rq.markings = &set;
    rq.w = 5.5;
    rq.h = 3.4;
    auto shapingInLastStrip = [](const TypedStitchMesh& s) {
        const int last = static_cast<int>(s.rings.size()) - 2;
        int shaping = 0;
        for (int f = 0; f < s.mesh.faces(); ++f)
            if (s.strip[f] == last && (s.kind[f] == StitchFaceKind::Decrease || s.kind[f] == StitchFaceKind::Increase)) ++shaping;
        return shaping;
    };
    MarkedBuildProduct fromLength, below;
    std::string err;
    rq.options.bindOffCount = BindOffCount::FromLoopLength;
    CHECK(buildFromMarkings(rq, fromLength, &err));
    rq.options.bindOffCount = BindOffCount::OfRoundBelow;
    CHECK(buildFromMarkings(rq, below, &err));
    const TypedStitchMesh& a = fromLength.evaluation.before;
    const TypedStitchMesh& b = below.evaluation.before;
    std::printf("  bind-off ring: from the loop %zu stitches over %zu, %d shaping faces in the last strip; of the round below %zu over %zu, %d\n", a.rings.back().size(),
                a.rings[a.rings.size() - 2].size(), shapingInLastStrip(a), b.rings.back().size(), b.rings[b.rings.size() - 2].size(), shapingInLastStrip(b));
    CHECK(a.rings.back().size() < a.rings[a.rings.size() - 2].size() && shapingInLastStrip(a) > 0);
    CHECK(b.rings.back().size() == b.rings[b.rings.size() - 2].size() && shapingInLastStrip(b) == 0);
    CHECK(below.evaluation.storedBefore.connectivity == Status::Pass && below.evaluation.shiftBefore.status == Status::Pass);
    CHECK(below.evaluation.storedAfter.invertedFaces == 0);
}
