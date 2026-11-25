// Surface markings below the mesher: paths on a mesh (model/surface_path.h), the pieces of an imported file
// (model/surface_pieces.h) and the markings file (io/marking_file.h). Each case has a control that must fail.
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "harness.h"
#include "io/marking_file.h"
#include "model/surface_marking.h"
#include "model/surface_path.h"
#include "model/surface_pieces.h"

using namespace knit;

namespace {

// An n x n grid of unit squares in the plane z = 0, each split into two triangles, offset by (x0, 0, 0).
TriMesh grid(int n, double x0 = 0) {
    TriMesh m;
    for (int j = 0; j <= n; ++j)
        for (int i = 0; i <= n; ++i) m.pos.push_back(Vec3{x0 + i, double(j), 0});
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
            const int a = j * (n + 1) + i, b = a + 1, c = a + n + 1, d = c + 1;
            m.tri.push_back({a, b, d});
            m.tri.push_back({a, d, c});
        }
    m.computeNormals();
    return m;
}

// An open cylinder of radius r around z, from z = 0 to z = height, segments x rings quads.
TriMesh cylinder(double r, double height, int segments, int rings) {
    TriMesh m;
    for (int j = 0; j <= rings; ++j)
        for (int i = 0; i < segments; ++i) {
            const double a = 2 * kPi * i / segments;
            m.pos.push_back(Vec3{r * std::cos(a), r * std::sin(a), height * j / rings});
        }
    for (int j = 0; j < rings; ++j)
        for (int i = 0; i < segments; ++i) {
            const int a = j * segments + i, b = j * segments + (i + 1) % segments, c = a + segments, d = b + segments;
            m.tri.push_back({a, b, d});
            m.tri.push_back({a, d, c});
        }
    m.computeNormals();
    return m;
}

// The 12 triangles of an axis-aligned cube at `origin`, as a triangle soup (every triangle its own three vertices).
void addCubeSoup(const Vec3& origin, std::vector<Vec3>& pos, std::vector<std::array<int, 3>>& tri) {
    const int faces[6][4] = {{0, 2, 3, 1}, {4, 5, 7, 6}, {0, 1, 5, 4}, {2, 6, 7, 3}, {0, 4, 6, 2}, {1, 3, 7, 5}};
    Vec3 c[8];
    for (int k = 0; k < 8; ++k) c[k] = origin + Vec3{double(k & 1), double((k >> 1) & 1), double((k >> 2) & 1)};
    for (const int* f : faces)
        for (const std::array<int, 3>& t : {std::array<int, 3>{f[0], f[1], f[2]}, std::array<int, 3>{f[0], f[2], f[3]}}) {
            const int base = static_cast<int>(pos.size());
            for (int k = 0; k < 3; ++k) pos.push_back(c[t[k]]);
            tri.push_back({base, base + 1, base + 2});
        }
}

double pointToSegment(const Vec3& p, const Vec3& a, const Vec3& b) {
    const Vec3 d = b - a;
    const double s = std::min(1.0, std::max(0.0, dot(p - a, d) / dot(d, d)));
    return norm(p - (a + d * s));
}

MarkingSet sampleSet() {
    MarkingSet s;
    s.source = "scenes/teddy-bears/obj/bears.obj";
    s.digest = 0x0123456789abcdefull;
    SurfaceMarking a;
    a.id = 1;
    a.meaning = MarkingMeaning::Start;
    a.shape = MarkingShape::ClosedCurve;
    a.note = "around the paw";
    a.points.push_back(SurfacePoint{12, {1.0 / 3, 1.0 / 7, 1 - 1.0 / 3 - 1.0 / 7}});
    a.points.push_back(SurfacePoint{40, {0.5, 0.25, 0.25}});
    SurfaceMarking b;
    b.id = 3;
    b.meaning = MarkingMeaning::RowStartPath;
    b.shape = MarkingShape::Path;
    b.origin = MarkingOrigin::Proposal;
    b.proposal = ProposalState::Offered;
    b.alternativeGroup = 2;
    b.reason = "steepest descent from the end pole";
    b.unsuitable = "";
    b.points.push_back(SurfacePoint{7, {0.1, 0.2, 0.7}});
    s.markings = {a, b};
    s.nextId = 4;
    return s;
}

} // namespace

// On a plane the path is the straight segment: every point on it, in order, and its length the distance. A control:
// two grids apart cannot be joined at all.
TEST(markings_trace_flat_grid) {
    const TriMesh m = grid(10);
    const TriangleNeighbours nb = triangleNeighbours(m);
    const SurfacePoint a = closestSurfacePoint(m, Vec3{0.13, 0.27, 0}), b = closestSurfacePoint(m, Vec3{7.71, 5.39, 0});
    std::vector<TracePoint> path;
    bool edges = true;
    CHECK(traceSurfacePath(m, nb, a, b, path, edges));
    CHECK(!edges && path.size() > 10);
    const Vec3 pa = surfacePosition(m, a), pb = surfacePosition(m, b);
    double length = 0, along = -1;
    for (size_t i = 0; i < path.size(); ++i) {
        CHECK(pointToSegment(path[i].pos, pa, pb) < 1e-12);
        const double t = dot(path[i].pos - pa, pb - pa);
        CHECK(t >= along);  // in order
        along = t;
        if (i > 0) length += norm(path[i].pos - path[i - 1].pos);
        if (i > 0 && i + 1 < path.size()) CHECK(path[i].edgeFrom >= 0);  // every inner point lies on an edge
    }
    CHECK_NEAR(length, norm(pb - pa), 1e-9);
    // control: two separate grids
    TriMesh two = grid(4);
    const TriMesh other = grid(4, 10);
    const int base = static_cast<int>(two.pos.size());
    two.pos.insert(two.pos.end(), other.pos.begin(), other.pos.end());
    for (const std::array<int, 3>& t : other.tri) two.tri.push_back({t[0] + base, t[1] + base, t[2] + base});
    two.computeNormals();
    std::vector<TracePoint> none;
    CHECK(!traceSurfacePath(two, triangleNeighbours(two), closestSurfacePoint(two, Vec3{1, 1, 0}), closestSurfacePoint(two, Vec3{11, 1, 0}), none, edges));
}

// A closed curve of three points round a cylinder at one height: every point at that height (the plane through two
// points and their mean normal is horizontal there), closed, and as long as the section polygon.
TEST(markings_trace_cylinder_wrap) {
    const int segments = 64;
    const TriMesh m = cylinder(10, 20, segments, 10);
    const TriangleNeighbours nb = triangleNeighbours(m);
    std::vector<SurfacePoint> pts;
    for (int k = 0; k < 3; ++k) {
        const double a = 2 * kPi * k / 3 + 0.05;
        pts.push_back(closestSurfacePoint(m, Vec3{10 * std::cos(a), 10 * std::sin(a), 7.3}));
    }
    const MarkingTrace tr = traceMarking(m, nb, pts, true);
    CHECK(!tr.edgePathUsed && tr.failedSegments.empty() && tr.controlAt.size() == 3);
    double length = 0;
    for (size_t i = 0; i < tr.points.size(); ++i) {
        CHECK_NEAR(tr.points[i].pos.z, 7.3, 1e-9);
        if (i > 0) length += norm(tr.points[i].pos - tr.points[i - 1].pos);
    }
    CHECK(norm(tr.points.back().pos - tr.points.front().pos) < 1e-12);
    const double polygon = 2 * segments * 10 * std::sin(kPi / segments);
    CHECK(std::fabs(length - polygon) < 0.01 * polygon);
    // control: the same points as an open path stop at the last point (two thirds of the way round)
    const MarkingTrace open = traceMarking(m, nb, pts, false);
    double openLength = 0;
    for (size_t i = 1; i < open.points.size(); ++i) openLength += norm(open.points[i].pos - open.points[i - 1].pos);
    CHECK(std::fabs(openLength - polygon * 2 / 3) < 0.02 * polygon);
}

// A closed latitude-longitude sphere of radius r: `rings` latitude rings between the poles, `segments` vertices per ring,
// each quad split into two triangles, a fan at each pole. Every vertex lies on the grid, so a plane through two of
// them often passes through further vertices: the degenerate case of the surface tracer.
TriMesh latLongSphere(double r, int segments, int rings) {
    TriMesh m;
    m.pos.push_back(Vec3{0, 0, r});
    for (int j = 1; j <= rings; ++j) {
        const double polar = kPi * j / (rings + 1);
        for (int i = 0; i < segments; ++i) {
            const double a = 2 * kPi * i / segments;
            m.pos.push_back(Vec3{r * std::sin(polar) * std::cos(a), r * std::sin(polar) * std::sin(a), r * std::cos(polar)});
        }
    }
    m.pos.push_back(Vec3{0, 0, -r});
    const int south = static_cast<int>(m.pos.size()) - 1;
    auto at = [&](int j, int i) { return 1 + (j - 1) * segments + (i % segments); };
    for (int i = 0; i < segments; ++i) m.tri.push_back({0, at(1, i), at(1, i + 1)});
    for (int j = 1; j < rings; ++j)
        for (int i = 0; i < segments; ++i) {
            m.tri.push_back({at(j, i), at(j + 1, i), at(j + 1, i + 1)});
            m.tri.push_back({at(j, i), at(j + 1, i + 1), at(j, i + 1)});
        }
    for (int i = 0; i < segments; ++i) m.tri.push_back({south, at(rings, i + 1), at(rings, i)});
    m.computeNormals();
    return m;
}

// The surface point exactly at vertex v, bound to the first triangle that uses it (as a click snapped to a vertex).
SurfacePoint atVertex(const TriMesh& m, int v) {
    for (int t = 0; t < static_cast<int>(m.tri.size()); ++t)
        for (int k = 0; k < 3; ++k)
            if (m.tri[t][k] == v) {
                SurfacePoint p;
                p.tri = t;
                p.bary[0] = p.bary[1] = p.bary[2] = 0;
                p.bary[k] = 1;
                return p;
            }
    return SurfacePoint{};
}

// Control points exactly on vertices (and on edges) of a closed mesh: each traced segment is the short arc between its
// ends, never the long way round the section. The degenerate start (the section leaves a vertex through its own corner)
// used to send the forward and the backward walk the same arbitrary way: on this sphere some neighbouring vertices
// were joined by a path more than ten times their distance. Control: an arc from a vertex to the vertex opposite it
// across the sphere is as long as half the section either way (no short way exists), and still traced.
TEST(markings_trace_vertex_points_take_the_short_way) {
    const TriMesh m = latLongSphere(10, 24, 11);
    const TriangleNeighbours nb = triangleNeighbours(m);
    auto traced = [&](const SurfacePoint& a, const SurfacePoint& b, double& straight) {
        std::vector<TracePoint> path;
        bool edges = false;
        straight = norm(surfacePosition(m, b) - surfacePosition(m, a));
        if (!traceSurfacePath(m, nb, a, b, path, edges)) return -1.0;
        double L = 0;
        for (size_t i = 1; i < path.size(); ++i) L += norm(path[i].pos - path[i - 1].pos);
        return L;
    };
    int pairs = 0, longWay = 0;
    double worst = 0;
    // neighbours along a ring, along a meridian, diagonally, and two apart, from every vertex of three rings
    for (int j : {2, 6, 9})
        for (int i = 0; i < 24; ++i) {
            const int v = 1 + (j - 1) * 24 + i;
            for (int w : {1 + (j - 1) * 24 + (i + 1) % 24, 1 + j * 24 + i, 1 + j * 24 + (i + 1) % 24, 1 + (j - 1) * 24 + (i + 2) % 24}) {
                double straight = 0;
                const double L = traced(atVertex(m, v), atVertex(m, w), straight);
                ++pairs;
                CHECK(L > 0);
                worst = std::max(worst, L / straight);
                if (L > 1.5 * straight) ++longWay;
            }
        }
    CHECK(pairs == 288);
    CHECK(longWay == 0);
    CHECK(worst < 1.2);
    // points on edges: the midpoint of each ring edge to the midpoint of the edge below it
    for (int i = 0; i < 24; ++i) {
        const int a0 = 1 + 4 * 24 + i, a1 = 1 + 4 * 24 + (i + 1) % 24, b0 = a0 + 24, b1 = a1 + 24;
        const SurfacePoint pa = closestSurfacePoint(m, (m.pos[a0] + m.pos[a1]) * 0.5), pb = closestSurfacePoint(m, (m.pos[b0] + m.pos[b1]) * 0.5);
        double straight = 0;
        const double L = traced(pa, pb, straight);
        CHECK(L > 0 && L < 1.2 * straight);
    }
    // control: a vertex to its opposite vertex (same ring, half a turn away): both ways are half the section
    double straight = 0;
    const double half = traced(atVertex(m, 1 + 5 * 24), atVertex(m, 1 + 5 * 24 + 12), straight);
    CHECK(half > 1.4 * straight);
}

// Two cubes written as triangle soups (every triangle its own vertices, as a file with texture seams everywhere) are
// two pieces of eight vertices each, numbered by their lowest triangle, and every file triangle maps to its piece with
// the same corner order. Controls: no welding (tolerance 0) leaves 24 pieces; different components split pieces.
TEST(markings_pieces_weld_uv_seams) {
    std::vector<Vec3> pos;
    std::vector<std::array<int, 3>> tri;
    addCubeSoup(Vec3{0, 0, 0}, pos, tri);
    addCubeSoup(Vec3{3, 0, 0}, pos, tri);
    const std::vector<int> oneNode(tri.size(), 0);
    const SurfacePieces p = findSurfacePieces(pos, tri, oneNode, 1e-9);
    CHECK(p.pieces.size() == 2);
    CHECK(p.pieces[0].mesh.pos.size() == 8 && p.pieces[1].mesh.pos.size() == 8);
    CHECK(p.pieces[0].ordinal == 1 && p.pieces[1].ordinal == 2);
    CHECK(p.pieceOfTriangle[0] == 0 && p.pieceOfTriangle[12] == 1 && p.localTriangle[13] == 1);
    int piece = -1;
    SurfacePoint local;
    CHECK(pieceLocalPoint(p, SurfacePoint{14, {0.2, 0.3, 0.5}}, piece, local));
    CHECK(piece == 1 && local.tri == 2);
    const std::array<int, 3>& ft = tri[14];
    const std::array<int, 3>& lt = p.pieces[1].mesh.tri[2];
    for (int k = 0; k < 3; ++k) CHECK(norm(pos[ft[k]] - p.pieces[1].mesh.pos[lt[k]]) == 0);
    const SurfacePieces soup = findSurfacePieces(pos, tri, oneNode, 0);
    CHECK(soup.pieces.size() == 24);
    std::vector<int> halves(tri.size(), 0);
    for (int t = 6; t < 12; ++t) halves[t] = 1;  // the first cube's second half is another component
    const SurfacePieces split = findSurfacePieces(pos, tri, halves, 1e-9);
    CHECK(split.pieces.size() == 3 && split.pieces[1].node == 1 && split.pieces[1].ordinal == 1);
}

// A set written and read back is the same set, weights bit for bit, and writes the same text again.
TEST(io_markings_roundtrip) {
    const MarkingSet s = sampleSet();
    const std::string text = writeMarkings(s);
    MarkingSet r;
    std::string err;
    CHECK(readMarkings(text, r, &err));
    CHECK(r.source == s.source && r.digest == s.digest && r.nextId == 4 && r.markings.size() == 2);
    for (size_t i = 0; i < s.markings.size(); ++i) {
        const SurfaceMarking& a = s.markings[i];
        const SurfaceMarking& b = r.markings[i];
        CHECK(a.id == b.id && a.meaning == b.meaning && a.shape == b.shape && a.origin == b.origin && a.proposal == b.proposal);
        CHECK(a.alternativeGroup == b.alternativeGroup && a.note == b.note && a.reason == b.reason && a.points.size() == b.points.size());
        for (size_t k = 0; k < a.points.size(); ++k) {
            CHECK(a.points[k].tri == b.points[k].tri);
            for (int c = 0; c < 3; ++c) CHECK(a.points[k].bary[c] == b.points[k].bary[c]);
        }
    }
    CHECK(writeMarkings(r) == text);
}

// The reader refuses what it does not know, with the line: an unknown key, a missing header, a short point, a
// repeated id, a marking without its end. Control: the same text without the defect reads.
TEST(io_markings_refuses_unknown_key) {
    const std::string good = "knit-markings v1\nsource = a.obj\ndigest = 00000000000000ff\nmarking 1\nmeaning = cut\nshape = path\npoint 3 0.5 0.25 0.25\nend\n";
    MarkingSet s;
    std::string err;
    CHECK(readMarkings(good, s, &err));
    CHECK(s.digest == 0xff && s.markings[0].meaning == MarkingMeaning::Cut && s.nextId == 2);
    std::string bad = good;
    bad.replace(bad.find("shape = path"), 12, "colour = red");
    CHECK(!readMarkings(bad, s, &err) && err.find("line 6") != std::string::npos && err.find("unknown key") != std::string::npos);
    CHECK(!readMarkings("knit-settings v1\n", s, &err));
    bad = good;
    bad.replace(bad.find("point 3 0.5 0.25 0.25"), 21, "point 3 0.5 0.25");
    CHECK(!readMarkings(bad, s, &err));
    CHECK(!readMarkings(good + "marking 1\nend\n", s, &err) && err.find("used twice") != std::string::npos);
    CHECK(!readMarkings(good + "marking 2\n", s, &err));
    bad = good;
    bad.replace(bad.find("meaning = cut"), 13, "meaning = seam");
    CHECK(!readMarkings(bad, s, &err));
}

// The format's text is pinned: a change of it is a version change (io/ARCHITECTURE.md). Re-pinned for v5 (knitting
// components): the header, and a marking's knit_component line, written only when > 0; before: re-pinned for v4 (assisted
// stockinette): the header, and a proposal's setup and basis lines, written only when set.
TEST(io_pin_markings_text) {
    MarkingSet s = sampleSet();
    s.markings[0].points.resize(1);
    s.markings[0].points[0] = SurfacePoint{12, {0.5, 0.25, 0.25}};
    s.markings[1].setup = 2;
    s.markings[1].basis = 0x00ff00ff00ff00ffull;
    s.markings[1].knitComponent = 1;
    const std::string expected =
        "knit-markings v5\n"
        "source = scenes/teddy-bears/obj/bears.obj\n"
        "digest = 0123456789abcdef\n"
        "next_id = 4\n"
        "preparation = facets\n"
        "marking 1\n"
        "meaning = start\n"
        "shape = closed\n"
        "origin = user\n"
        "proposal = none\n"
        "proposal_id = 0\n"
        "group = 0\n"
        "note = around the paw\n"
        "point 12 0.5 0.25 0.25\n"
        "end\n"
        "marking 3\n"
        "meaning = row_start\n"
        "shape = path\n"
        "origin = proposal\n"
        "proposal = offered\n"
        "proposal_id = 0\n"
        "group = 2\n"
        "setup = 2\n"
        "knit_component = 1\n"
        "basis = 00ff00ff00ff00ff\n"
        "reason = steepest descent from the end pole\n"
        "point 7 0.10000000000000001 0.20000000000000001 0.69999999999999996\n"
        "end\n";
    CHECK(writeMarkings(s) == expected);
}

// File version 2 adds per-marking settings; a v1 file loads (version 1 recorded) and is rewritten as v2; a setting in a v1
// file and an unknown newer version are refused with a message. A finish may be a point (a bind-off level) in v2.
TEST(io_markings_reads_v1_writes_v2) {
    const std::string v1 = "knit-markings v1\nsource = a.obj\ndigest = 00000000000000ff\nnext_id = 2\nmarking 1\nmeaning = start\nshape = point\npoint 3 0.5 0.25 0.25\nend\n";
    MarkingSet s;
    std::string err;
    CHECK(readMarkings(v1, s, &err) && s.version == 1 && s.markings.size() == 1);
    s.markings[0].settings = {{"ring1", "5"}};
    const std::string v2 = writeMarkings(s);
    CHECK(v2.rfind("knit-markings v5\n", 0) == 0 && v2.find("setting ring1 = 5\n") != std::string::npos);
    MarkingSet back;
    CHECK(readMarkings(v2, back, &err) && back.version == 5 && markingSetting(back.markings[0], "ring1") == "5");
    const std::string v1WithSetting = "knit-markings v1\nsource = a.obj\ndigest = 00000000000000ff\nmarking 1\nmeaning = start\nshape = point\nsetting ring1 = 5\nend\n";
    CHECK(!readMarkings(v1WithSetting, back, &err) && err.find("v2") != std::string::npos);
    CHECK(!readMarkings("knit-markings v6\n", back, &err) && err.find("newer") != std::string::npos);
    CHECK(shapeAllowed(MarkingMeaning::Finish, MarkingShape::Point) && !shapeAllowed(MarkingMeaning::Boundary, MarkingShape::Point));
}

// File version 3 adds the set's preparation (the reference the mesher prepares). A v2 file loads with "facets" (what it
// meant) and is rewritten as v3 with the line; the line in a v2 file, or a value other than facets / loop_limit, refuses
// the file; a v4 file is newer than this program.
TEST(io_markings_reads_v2_writes_v3) {
    const std::string v2 = "knit-markings v2\nsource = a.obj\ndigest = 00000000000000ff\nnext_id = 2\nmarking 1\nmeaning = start\nshape = point\nsetting ring1 = 5\npoint 3 0.5 0.25 0.25\nend\n";
    MarkingSet s;
    std::string err;
    CHECK(readMarkings(v2, s, &err) && s.version == 2 && s.preparation == "facets" && markingSetting(s.markings[0], "ring1") == "5");
    s.preparation = "loop_limit";
    const std::string v3 = writeMarkings(s);
    CHECK(v3.rfind("knit-markings v5\n", 0) == 0 && v3.find("preparation = loop_limit\n") != std::string::npos);
    MarkingSet back;
    CHECK(readMarkings(v3, back, &err) && back.version == 5 && back.preparation == "loop_limit" && back.markings.size() == 1);
    CHECK(back.markings[0].points[0].tri == 3 && back.markings[0].points[0].bary[1] == 0.25);
    std::string inV2 = v2;
    inV2.insert(inV2.find("marking 1"), "preparation = loop_limit\n");
    CHECK(!readMarkings(inV2, back, &err) && err.find("v3") != std::string::npos);
    std::string bad = v3;
    bad.replace(bad.find("loop_limit"), 10, "catmull");
    CHECK(!readMarkings(bad, back, &err) && err.find("not facets or loop_limit") != std::string::npos);
    CHECK(!readMarkings("knit-markings v6\n", back, &err) && err.find("newer") != std::string::npos);
}

// File version 4 adds a proposal's setup and basis. A v3 file loads unchanged (setup 0, basis 0) and is rewritten as v4; a
// setup or basis line in a v3 file refuses it (with the line); both round-trip; v5 is newer than this program.
// (Defect: the reader takes "setup" from a v3 file, or the writer drops the basis.)
TEST(io_markings_reads_v3_writes_v4) {
    const std::string v3 = "knit-markings v3\nsource = a.obj\ndigest = 00000000000000ff\nnext_id = 3\npreparation = facets\n"
                           "marking 1\nmeaning = start\nshape = closed\norigin = proposal\nproposal = offered\nproposal_id = 1\ngroup = 1\n"
                           "point 3 0.5 0.25 0.25\npoint 4 0.5 0.25 0.25\npoint 5 0.5 0.25 0.25\nend\n";
    MarkingSet s;
    std::string err;
    CHECK(readMarkings(v3, s, &err) && s.version == 3 && s.markings[0].setup == 0 && s.markings[0].basis == 0);
    s.markings[0].setup = 7;
    s.markings[0].basis = 0xfedcba9876543210ull;
    const std::string v4 = writeMarkings(s);
    CHECK(v4.rfind("knit-markings v5\n", 0) == 0 && v4.find("setup = 7\n") != std::string::npos && v4.find("basis = fedcba9876543210\n") != std::string::npos);
    MarkingSet back;
    CHECK(readMarkings(v4, back, &err) && back.version == 5 && back.markings[0].setup == 7 && back.markings[0].basis == 0xfedcba9876543210ull);
    CHECK(writeMarkings(back) == v4);
    std::string inV3 = v3;
    inV3.insert(inV3.find("point 3"), "setup = 2\n");
    CHECK(!readMarkings(inV3, back, &err) && err.find("v4") != std::string::npos && err.find("line") != std::string::npos);
    std::string bad = v4;
    bad.replace(bad.find("setup = 7"), 9, "setup = x");
    CHECK(!readMarkings(bad, back, &err));
    CHECK(!readMarkings("knit-markings v6\n", back, &err) && err.find("newer") != std::string::npos);
}

// File version 5 adds a marking's knitting component. A v4 file loads with every marking in component 0 and is rewritten as
// v5 without a component line (0 is not written); a component round-trips; a knit_component line in a v4 file, or a
// negative one, refuses the file; v6 is newer. (Defect: the reader takes "knit_component" from a v4 file, or the writer
// drops it.)
TEST(io_markings_reads_v4_writes_v5) {
    const std::string v4 = "knit-markings v4\nsource = a.obj\ndigest = 00000000000000ff\nnext_id = 2\npreparation = facets\n"
                           "marking 1\nmeaning = start\nshape = point\norigin = user\nproposal = none\nproposal_id = 0\ngroup = 0\n"
                           "point 3 0.5 0.25 0.25\nend\n";
    MarkingSet s;
    std::string err;
    CHECK(readMarkings(v4, s, &err) && s.version == 4 && s.markings[0].knitComponent == 0);
    const std::string asV5 = writeMarkings(s);
    CHECK(asV5.rfind("knit-markings v5\n", 0) == 0 && asV5.find("knit_component") == std::string::npos);
    s.markings[0].knitComponent = 2;
    const std::string v5 = writeMarkings(s);
    CHECK(v5.find("knit_component = 2\n") != std::string::npos);
    MarkingSet back;
    CHECK(readMarkings(v5, back, &err) && back.version == 5 && back.markings[0].knitComponent == 2 && writeMarkings(back) == v5);
    std::string inV4 = v4;
    inV4.insert(inV4.find("point 3"), "knit_component = 1\n");
    CHECK(!readMarkings(inV4, back, &err) && err.find("v5") != std::string::npos);
    std::string negative = v5;
    negative.replace(negative.find("knit_component = 2"), 18, "knit_component = -1");
    CHECK(!readMarkings(negative, back, &err));
    CHECK(!readMarkings("knit-markings v6\n", back, &err) && err.find("newer") != std::string::npos);
}

namespace {

// Two pieces of one file: two grids apart (piece 0 from triangle 0, piece 1 after it). A marking on a piece is placed by
// its first triangle.
struct TwoPieces {
    SurfacePieces pieces;
    int firstTriOf[2] = {0, 0};
};

TwoPieces twoPieces() {
    const TriMesh a = grid(4), b = grid(4, 10);
    std::vector<Vec3> pos = a.pos;
    std::vector<std::array<int, 3>> tri = a.tri;
    const int base = static_cast<int>(pos.size());
    for (const Vec3& p : b.pos) pos.push_back(p);
    for (const std::array<int, 3>& t : b.tri) tri.push_back({t[0] + base, t[1] + base, t[2] + base});
    TwoPieces t;
    t.pieces = findSurfacePieces(pos, tri, std::vector<int>(tri.size(), 0), 1e-9);
    t.firstTriOf[1] = static_cast<int>(a.tri.size());
    return t;
}

SurfaceMarking markingOn(int id, int tri, MarkingMeaning meaning, MarkingOrigin origin, ProposalState state, int setup, uint64_t basis) {
    SurfaceMarking m;
    m.id = id;
    m.meaning = meaning;
    m.shape = defaultShape(meaning);
    m.origin = origin;
    m.proposal = state;
    m.setup = setup;
    m.basis = basis;
    m.points.push_back(SurfacePoint{tri, {1.0 / 3, 1.0 / 3, 1.0 / 3}});
    return m;
}

int countState(const MarkingSet& s, ProposalState p) {
    int n = 0;
    for (const SurfaceMarking& m : s.markings) n += m.origin == MarkingOrigin::Proposal && m.proposal == p ? 1 : 0;
    return n;
}

} // namespace

// Accepting a setup takes exactly its own start, finish and path, never another setup's; an accepted setup of the same
// roles is set aside (rejected); user markings only when the user asked to replace them; an outdated basis or another
// piece refuses with nothing changed; rejecting a setup rejects its offered markings only. (Defect: acceptance by
// alternative group, which mixes one start with another start's finish.)
TEST(markings_setup_accept_reject) {
    const TwoPieces tp = twoPieces();
    const uint64_t basis = 0xabc;
    MarkingSet s;
    s.markings = {markingOn(1, 0, MarkingMeaning::Start, MarkingOrigin::Proposal, ProposalState::Offered, 1, basis),
                  markingOn(2, 1, MarkingMeaning::Finish, MarkingOrigin::Proposal, ProposalState::Offered, 1, basis),
                  markingOn(3, 2, MarkingMeaning::RowStartPath, MarkingOrigin::Proposal, ProposalState::Offered, 1, basis),
                  markingOn(4, 3, MarkingMeaning::Start, MarkingOrigin::Proposal, ProposalState::Offered, 2, basis),
                  markingOn(5, 4, MarkingMeaning::Finish, MarkingOrigin::Proposal, ProposalState::Offered, 2, basis),
                  markingOn(6, 5, MarkingMeaning::RowStartPath, MarkingOrigin::Proposal, ProposalState::Offered, 2, basis)};
    for (SurfaceMarking& m : s.markings) m.alternativeGroup = static_cast<int>(m.meaning) + 1;   // the old grouping: by meaning across setups
    s.nextId = 7;
    CHECK(setupMarkings(s, 1) == (std::vector<int>{1, 2, 3}));
    CHECK(acceptSetup(s, tp.pieces, 0, 1, basis, false));
    CHECK(findMarking(s, 1)->proposal == ProposalState::Accepted && findMarking(s, 2)->proposal == ProposalState::Accepted && findMarking(s, 3)->proposal == ProposalState::Accepted);
    CHECK(findMarking(s, 4)->proposal == ProposalState::Offered && findMarking(s, 5)->proposal == ProposalState::Offered);
    // the second setup replaces the first: the first's markings are rejected, not mixed with the second's
    CHECK(acceptSetup(s, tp.pieces, 0, 2, basis, false));
    CHECK(countState(s, ProposalState::Accepted) == 3 && findMarking(s, 4)->proposal == ProposalState::Accepted);
    CHECK(findMarking(s, 1)->proposal == ProposalState::Rejected && findMarking(s, 2)->proposal == ProposalState::Rejected);
    // a user start on the piece: refused without replace (nothing changes), removed with it
    MarkingSet u = s;
    u.markings.push_back(markingOn(7, 6, MarkingMeaning::Start, MarkingOrigin::User, ProposalState::None, 0, 0));
    u.markings.push_back(markingOn(8, 7, MarkingMeaning::Start, MarkingOrigin::Proposal, ProposalState::Offered, 3, basis));
    u.markings.push_back(markingOn(9, 7, MarkingMeaning::Finish, MarkingOrigin::Proposal, ProposalState::Offered, 3, basis + 1));   // outdated
    const MarkingSet before = u;
    CHECK(!acceptSetup(u, tp.pieces, 0, 3, basis, true));            // one of its markings is outdated
    u.markings.back().basis = basis;
    CHECK(setupConflicts(u, tp.pieces, 0, 3).userMarkings == std::vector<int>{7});
    CHECK(!acceptSetup(u, tp.pieces, 0, 3, basis, false) && u.markings.size() == before.markings.size() && findMarking(u, 8)->proposal == ProposalState::Offered);
    CHECK(acceptSetup(u, tp.pieces, 0, 3, basis, true) && !findMarking(u, 7) && findMarking(u, 8)->proposal == ProposalState::Accepted);
    CHECK(!acceptSetup(u, tp.pieces, 1, 3, basis, true));            // not on piece 1
    // reject: the offered markings of the setup only
    MarkingSet r;
    r.markings = {markingOn(1, 0, MarkingMeaning::Start, MarkingOrigin::Proposal, ProposalState::Offered, 5, basis),
                  markingOn(2, 1, MarkingMeaning::Finish, MarkingOrigin::User, ProposalState::None, 5, basis),   // an edit of the setup's finish
                  markingOn(3, 2, MarkingMeaning::Start, MarkingOrigin::Proposal, ProposalState::Offered, 6, basis)};
    rejectSetup(r, 5);
    CHECK(findMarking(r, 1)->proposal == ProposalState::Rejected && findMarking(r, 2)->origin == MarkingOrigin::User && findMarking(r, 3)->proposal == ProposalState::Offered);
}

// Conflicts are by construction role: a user Boundary is the same role as a setup's Finish (one termination per piece),
// so accepting with replace removes it and without replace nothing changes; a cut and an attachment on the piece and
// markings on the other piece are never touched. (Defect: conflicts by meaning, which keeps the Boundary beside the
// accepted Finish: two terminations.)
TEST(markings_setup_role_conflicts) {
    const TwoPieces tp = twoPieces();
    CHECK(markingRole(MarkingMeaning::Finish) == MarkingRole::Termination && markingRole(MarkingMeaning::Boundary) == MarkingRole::Termination);
    CHECK(markingRole(MarkingMeaning::Cut) == MarkingRole::None && markingRole(MarkingMeaning::Attachment) == MarkingRole::None);
    MarkingSet s;
    s.markings = {markingOn(1, 0, MarkingMeaning::Boundary, MarkingOrigin::User, ProposalState::None, 0, 0),
                  markingOn(2, 1, MarkingMeaning::Cut, MarkingOrigin::User, ProposalState::None, 0, 0),
                  markingOn(3, 2, MarkingMeaning::Attachment, MarkingOrigin::User, ProposalState::None, 0, 0),
                  markingOn(4, tp.firstTriOf[1], MarkingMeaning::Boundary, MarkingOrigin::User, ProposalState::None, 0, 0),   // the other piece
                  markingOn(5, 3, MarkingMeaning::Start, MarkingOrigin::Proposal, ProposalState::Offered, 1, 9),
                  markingOn(6, 4, MarkingMeaning::Finish, MarkingOrigin::Proposal, ProposalState::Offered, 1, 9),
                  markingOn(7, 5, MarkingMeaning::RowStartPath, MarkingOrigin::Proposal, ProposalState::Offered, 1, 9)};
    const SetupConflicts c = setupConflicts(s, tp.pieces, 0, 1);
    CHECK(c.userMarkings == std::vector<int>{1} && c.acceptedProposals.empty());
    MarkingSet keep = s;
    CHECK(!acceptSetup(keep, tp.pieces, 0, 1, 9, false) && findMarking(keep, 1) && findMarking(keep, 6)->proposal == ProposalState::Offered);
    CHECK(acceptSetup(s, tp.pieces, 0, 1, 9, true));
    CHECK(!findMarking(s, 1) && findMarking(s, 2) && findMarking(s, 3) && findMarking(s, 4));
}

// A preview of setup B with setup A accepted is reviewed on a set where B's start, termination and path stand alone for
// their roles: A's are left out of the preview set, everything else (a cut, an attachment, the other piece) stays, B's
// markings stay offered and are named for the preview, and the source set is unchanged byte for byte. (Defect: a
// preview that only adds B through previewOffered, giving two starts, two terminations and two paths.)
TEST(markings_preview_shadows_roles) {
    const TwoPieces tp = twoPieces();
    MarkingSet s;
    s.markings = {markingOn(1, 0, MarkingMeaning::Start, MarkingOrigin::Proposal, ProposalState::Accepted, 1, 9),
                  markingOn(2, 1, MarkingMeaning::Boundary, MarkingOrigin::User, ProposalState::None, 1, 9),   // A's finish, edited into a boundary
                  markingOn(3, 2, MarkingMeaning::RowStartPath, MarkingOrigin::Proposal, ProposalState::Accepted, 1, 9),
                  markingOn(4, 3, MarkingMeaning::Cut, MarkingOrigin::User, ProposalState::None, 0, 0),
                  markingOn(5, tp.firstTriOf[1], MarkingMeaning::Start, MarkingOrigin::User, ProposalState::None, 0, 0),
                  markingOn(6, 4, MarkingMeaning::Start, MarkingOrigin::Proposal, ProposalState::Offered, 2, 9),
                  markingOn(7, 5, MarkingMeaning::Finish, MarkingOrigin::Proposal, ProposalState::Offered, 2, 9),
                  markingOn(8, 6, MarkingMeaning::RowStartPath, MarkingOrigin::Proposal, ProposalState::Offered, 2, 9)};
    const std::string before = writeMarkings(s);
    const uint64_t digest = effectiveMarkingsDigest(s, tp.pieces, 0);
    std::vector<int> ids;
    const MarkingSet p = previewMarkings(s, tp.pieces, 0, 2, ids);
    CHECK(ids == (std::vector<int>{6, 7, 8}));
    CHECK(!findMarking(p, 1) && !findMarking(p, 2) && !findMarking(p, 3));
    CHECK(findMarking(p, 4) && findMarking(p, 5) && findMarking(p, 6)->proposal == ProposalState::Offered);
    int starts = 0, terminations = 0, paths = 0;
    for (const SurfaceMarking& m : p.markings) {
        if (pieceOfMarking(tp.pieces, m) != 0) continue;
        starts += markingRole(m.meaning) == MarkingRole::Start ? 1 : 0;
        terminations += markingRole(m.meaning) == MarkingRole::Termination ? 1 : 0;
        paths += markingRole(m.meaning) == MarkingRole::RowStart ? 1 : 0;
    }
    CHECK(starts == 1 && terminations == 1 && paths == 1);
    CHECK(writeMarkings(s) == before && effectiveMarkingsDigest(s, tp.pieces, 0) == digest);
}

// New proposals for a piece replace its offered and rejected proposals only: user markings, accepted proposals (an
// accepted setup stays accepted) and the other piece's proposals stay; every new marking gets a fresh id above all used
// ones and its setup a fresh number, the setups' grouping kept. (Defect: regeneration that clears every proposal of the
// piece, losing the accepted setup.)
TEST(markings_merge_keeps_user) {
    const TwoPieces tp = twoPieces();
    MarkingSet s;
    s.markings = {markingOn(1, 0, MarkingMeaning::Start, MarkingOrigin::Proposal, ProposalState::Accepted, 1, 9),
                  markingOn(2, 1, MarkingMeaning::Finish, MarkingOrigin::Proposal, ProposalState::Offered, 2, 9),
                  markingOn(3, 2, MarkingMeaning::Finish, MarkingOrigin::Proposal, ProposalState::Rejected, 3, 9),
                  markingOn(4, 3, MarkingMeaning::RowStartPath, MarkingOrigin::User, ProposalState::None, 1, 9),
                  markingOn(5, tp.firstTriOf[1], MarkingMeaning::Start, MarkingOrigin::Proposal, ProposalState::Offered, 4, 9)};
    s.nextId = 6;
    const std::vector<SurfaceMarking> offered = {markingOn(100, 4, MarkingMeaning::Start, MarkingOrigin::Proposal, ProposalState::Offered, 1, 10),
                                                 markingOn(101, 5, MarkingMeaning::Finish, MarkingOrigin::Proposal, ProposalState::Offered, 1, 10),
                                                 markingOn(102, 6, MarkingMeaning::Start, MarkingOrigin::Proposal, ProposalState::Offered, 2, 10)};
    mergeProposals(s, tp.pieces, 0, offered);
    CHECK(findMarking(s, 1) && findMarking(s, 1)->proposal == ProposalState::Accepted);
    CHECK(!findMarking(s, 2) && !findMarking(s, 3));
    CHECK(findMarking(s, 4) && findMarking(s, 5));
    CHECK(findMarking(s, 6) && findMarking(s, 7) && findMarking(s, 8) && s.nextId == 9);
    CHECK(findMarking(s, 6)->setup == 5 && findMarking(s, 7)->setup == 5 && findMarking(s, 8)->setup == 6);
    CHECK(findMarking(s, 6)->proposalId == 6 && findMarking(s, 6)->basis == 10);
}

// Two knitting components on one piece are kept apart: each component's effective digest depends on its own markings only
// (component 0's is the digest the set had before the second component was added); accepting a setup of component 1 sets
// aside only component 1's markings of its roles; new suggestions for component 1 replace only component 1's offered ones.
// (Defect: markings selected by piece alone, so a collar's start counts as a second start of the body.)
TEST(markings_knit_components_apart) {
    const TwoPieces tp = twoPieces();
    MarkingSet s;
    s.markings = {markingOn(1, 0, MarkingMeaning::Start, MarkingOrigin::Proposal, ProposalState::Accepted, 1, 9),
                  markingOn(2, 1, MarkingMeaning::Finish, MarkingOrigin::User, ProposalState::None, 0, 0),
                  markingOn(3, 2, MarkingMeaning::RowStartPath, MarkingOrigin::User, ProposalState::None, 0, 0)};
    s.nextId = 4;
    const uint64_t bodyBefore = effectiveMarkingsDigest(s, tp.pieces, 0);
    SurfaceMarking collarStart = markingOn(4, 3, MarkingMeaning::Start, MarkingOrigin::User, ProposalState::None, 0, 0);
    collarStart.knitComponent = 1;
    s.markings.push_back(collarStart);
    s.nextId = 5;
    CHECK(effectiveMarkingsDigest(s, tp.pieces, 0) == bodyBefore);       // the body's knit does not go stale
    CHECK(effectiveMarkingsDigest(s, tp.pieces, 0, 1) != bodyBefore);
    CHECK(effectiveMarkingsDigest(s, tp.pieces, 0, 1) != effectiveMarkingsDigest(s, tp.pieces, 0, 2));
    CHECK(markingInKnitComponent(tp.pieces, *findMarking(s, 4), 0, 1) && !markingInKnitComponent(tp.pieces, *findMarking(s, 4), 0, 0));
    // a suggested setup for the collar: its start replaces the collar's own start (on request), never the body's
    std::vector<SurfaceMarking> offered = {markingOn(100, 4, MarkingMeaning::Start, MarkingOrigin::Proposal, ProposalState::Offered, 1, 9),
                                           markingOn(101, 5, MarkingMeaning::Finish, MarkingOrigin::Proposal, ProposalState::Offered, 1, 9)};
    for (SurfaceMarking& m : offered) m.knitComponent = 1;
    std::vector<int> numbers;
    mergeProposals(s, tp.pieces, 0, offered, &numbers);
    CHECK(numbers.size() == 1);
    const SetupConflicts c = setupConflicts(s, tp.pieces, 0, numbers[0]);
    CHECK(c.userMarkings == std::vector<int>{4} && c.acceptedProposals.empty());
    CHECK(acceptSetup(s, tp.pieces, 0, numbers[0], 9, true));
    CHECK(findMarking(s, 1)->proposal == ProposalState::Accepted && findMarking(s, 2) && findMarking(s, 3) && !findMarking(s, 4));
    CHECK(effectiveMarkingsDigest(s, tp.pieces, 0) == bodyBefore);
    // new suggestions for the body leave the collar's accepted setup alone, and the collar's offered ones too
    std::vector<SurfaceMarking> collarOffer = {markingOn(200, 6, MarkingMeaning::Start, MarkingOrigin::Proposal, ProposalState::Offered, 1, 9)};
    collarOffer[0].knitComponent = 1;
    mergeProposals(s, tp.pieces, 0, collarOffer);
    const int collarOffered = s.nextId - 1;
    mergeProposals(s, tp.pieces, 0, {markingOn(300, 7, MarkingMeaning::Start, MarkingOrigin::Proposal, ProposalState::Offered, 1, 9)});
    CHECK(findMarking(s, collarOffered) && findMarking(s, collarOffered)->knitComponent == 1);
}

// A closed marking whose points all lie on one border loop is that loop: the points at the bottom rim's edge midpoints of
// an open cylinder give the rim (one vertex per segment), walked the way the points go (reversed points: the reversed
// loop), and its trace runs along the rim through every control point. Refused: a point off the rim, points on both rims,
// a chord between two rim vertices, too few points. (Defect: a near-border point taken as on it, or the loop's direction
// ignored.)
TEST(markings_border_loop_detection) {
    const int segments = 24;
    const TriMesh m = cylinder(2.0, 4.0, segments, 4);
    // the bottom rim's edges: vertices 0 .. segments - 1 at z = 0; triangle {a, b, d} walks a -> b
    auto rimPoint = [&](int i) {
        const int a = i, b = (i + 1) % segments;
        for (int t = 0; t < static_cast<int>(m.tri.size()); ++t)
            for (int k = 0; k < 3; ++k)
                if (m.tri[t][k] == a && m.tri[t][(k + 1) % 3] == b) {
                    SurfacePoint p{t, {0, 0, 0}};
                    p.bary[k] = 0.5;
                    p.bary[(k + 1) % 3] = 0.5;
                    return p;
                }
        return SurfacePoint{};
    };
    std::vector<SurfacePoint> pts;
    for (int i = 0; i < segments; i += 4) pts.push_back(rimPoint(i));
    std::vector<int> loop;
    CHECK(pointsOnOneBorderLoop(m, pts, loop) && static_cast<int>(loop.size()) == segments);
    for (int v : loop) CHECK(v < segments);
    // the direction follows the points: the walk visits 0 then 1 here, the reversed points the other way
    const auto at = [&](const std::vector<int>& l, int v) { return static_cast<int>(std::find(l.begin(), l.end(), v) - l.begin()); };
    CHECK((at(loop, 0) + 1) % segments == at(loop, 1));
    std::vector<SurfacePoint> back(pts.rbegin(), pts.rend());
    std::vector<int> loopBack;
    CHECK(pointsOnOneBorderLoop(m, back, loopBack) && (at(loopBack, 1) + 1) % segments == at(loopBack, 0));
    const MarkingTrace tr = borderLoopTrace(m, loop, pts);
    CHECK(tr.controlAt.size() == pts.size() && static_cast<int>(tr.points.size()) == segments + static_cast<int>(pts.size()) + 1);
    for (size_t k = 1; k < tr.controlAt.size(); ++k) CHECK(tr.controlAt[k] > tr.controlAt[k - 1]);
    for (const TracePoint& p : tr.points) CHECK(std::fabs(p.pos.z) < 1e-12);
    for (size_t k = 0; k < pts.size(); ++k) CHECK(norm(tr.points[tr.controlAt[k]].pos - surfacePosition(m, pts[k])) < 1e-12);
    // refusals
    std::vector<SurfacePoint> off = pts;
    off[1].bary[0] = 0.25;   // the weights no longer have an exact zero: inside the triangle
    off[1].bary[1] = 0.5;
    off[1].bary[2] = 0.25;
    CHECK(!pointsOnOneBorderLoop(m, off, loop));
    std::vector<SurfacePoint> both = pts;
    const SurfacePoint top = closestSurfacePoint(m, Vec3{2.0, 0.0, 4.0});
    both.push_back(top);
    CHECK(top.tri >= 0 && !pointsOnOneBorderLoop(m, both, loop));
    std::vector<SurfacePoint> chord = pts;
    chord[0] = SurfacePoint{0, {0.5, 0, 0.5}};   // triangle {0, 1, 25}: the edge 25 -> 0 is inside the cylinder, not on the rim
    CHECK(!pointsOnOneBorderLoop(m, chord, loop));
    CHECK(!pointsOnOneBorderLoop(m, {pts[0], pts[1]}, loop));
}
