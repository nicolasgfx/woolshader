// Stage 9 of the built-in import, flattenFaces (import/wrap_stages.h), part of knit_tests_import: small generated quad meshes
// on the CPU, each case well under a second.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <thread>
#include <vector>

#include "harness.h"
#include "import/wrap_stages.h"
#include "model/coordinates.h"


using namespace knit;

namespace {

// A quad mesh with the face graph flattenFaces reads (faces of each edge, edge of each corner). Open meshes are allowed: an
// edge with one face has -1 as its second. The graph points into the mesh, so the pair lives behind a unique_ptr.
struct TestSurface {
    QuadMesh mesh;
    FaceGraph graph;
    std::vector<int> label;  // a number per face, for seams between differently numbered faces
};

std::unique_ptr<TestSurface> surfaceOf(QuadMesh mesh, std::vector<int> label) {
    std::unique_ptr<TestSurface> s = std::make_unique<TestSurface>();
    s->mesh = std::move(mesh);
    s->label = std::move(label);
    s->mesh.edge = sortedEdges(s->mesh.quad);
    std::map<std::array<int, 2>, int> edgeOfPair;
    for (int e = 0; e < static_cast<int>(s->mesh.edge.size()); ++e) edgeOfPair[s->mesh.edge[e]] = e;
    const int faces = static_cast<int>(s->mesh.quad.size());
    FaceGraph& g = s->graph;
    g.mesh = &s->mesh;
    g.start.assign(faces + 1, 0);
    g.facesOfEdge.assign(s->mesh.edge.size(), {-1, -1});
    g.edgeOfCorner.resize(faces);
    for (int f = 0; f < faces; ++f)
        for (int k = 0; k < 4; ++k) {
            const int a = s->mesh.quad[f][k];
            const int b = s->mesh.quad[f][(k + 1) % 4];
            const int e = edgeOfPair.at({std::min(a, b), std::max(a, b)});
            g.edgeOfCorner[f][k] = e;
            std::array<int, 2>& pair = g.facesOfEdge[e];
            pair[pair[0] < 0 ? 0 : 1] = f;
        }
    return s;
}

// A seam on every edge between faces of different labels.
std::vector<char> seamsBetweenLabels(const TestSurface& s) {
    std::vector<char> seam(s.graph.facesOfEdge.size(), 0);
    for (size_t e = 0; e < seam.size(); ++e) {
        const std::array<int, 2>& pair = s.graph.facesOfEdge[e];
        seam[e] = pair[1] >= 0 && s.label[pair[0]] != s.label[pair[1]];
    }
    return seam;
}

// An n x n grid of unit quads in the plane z = 0, inner vertices moved by up to `jitter` along each axis: a disc. Counter-
// clockwise seen from +z.
QuadMesh planarGrid(int n, double jitter, unsigned seed) {
    std::mt19937 random(seed);
    std::uniform_real_distribution<double> unit(-1.0, 1.0);
    QuadMesh m;
    for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= n; ++j) {
            const bool inner = i > 0 && i < n && j > 0 && j < n;
            const double dx = jitter * unit(random), dy = jitter * unit(random);
            m.pos.push_back({i + (inner ? dx : 0.0), j + (inner ? dy : 0.0), 0.0});
        }
    const auto at = [n](int i, int j) { return i * (n + 1) + j; };
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) m.quad.push_back({at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)});
    return m;
}

// The upper half of a cube-sphere of radius 10 mm: the top face and the upper halves of the four sides of a cube with 2n quads
// a side, pushed out onto the sphere. Each face's label is its cube side (0 top, 1 to 4 the sides). Outward winding.
std::unique_ptr<TestSurface> hemisphere(int n) {
    QuadMesh m;
    std::vector<int> label;
    std::map<std::array<int, 3>, int> vertexOf;
    const auto vertexAt = [&](const std::array<int, 3>& p) {
        const std::map<std::array<int, 3>, int>::const_iterator it = vertexOf.find(p);
        if (it != vertexOf.end()) return it->second;
        const int v = static_cast<int>(m.pos.size());
        vertexOf[p] = v;
        m.pos.push_back(normalized(Vec3{double(p[0]), double(p[1]), double(p[2])}) * 10.0);
        return v;
    };
    // origin, first and second axis (their cross product points outward), counts along each
    struct Side {
        std::array<int, 3> origin, first, second;
        int count1, count2;
    };
    const std::array<Side, 5> sides = {{
        {{-n, -n, n}, {1, 0, 0}, {0, 1, 0}, 2 * n, 2 * n},
        {{n, -n, 0}, {0, 1, 0}, {0, 0, 1}, 2 * n, n},
        {{-n, -n, 0}, {0, 0, 1}, {0, 1, 0}, n, 2 * n},
        {{-n, n, 0}, {0, 0, 1}, {1, 0, 0}, n, 2 * n},
        {{-n, -n, 0}, {1, 0, 0}, {0, 0, 1}, 2 * n, n},
    }};
    for (int s = 0; s < 5; ++s) {
        const Side& side = sides[s];
        const auto corner = [&](int a, int b) {
            std::array<int, 3> p;
            for (int d = 0; d < 3; ++d) p[d] = side.origin[d] + a * side.first[d] + b * side.second[d];
            return vertexAt(p);
        };
        for (int a = 0; a < side.count1; ++a)
            for (int b = 0; b < side.count2; ++b) {
                m.quad.push_back({corner(a, b), corner(a + 1, b), corner(a + 1, b + 1), corner(a, b + 1)});
                label.push_back(s);
            }
    }
    return surfaceOf(std::move(m), std::move(label));
}

// A planar ring of quads between radii 1 and 2: two boundaries.
QuadMesh planarAnnulus(int rings, int sectors) {
    QuadMesh m;
    for (int r = 0; r <= rings; ++r)
        for (int s = 0; s < sectors; ++s) {
            const double radius = 1.0 + double(r) / rings;
            const double angle = 2 * kPi * s / sectors;
            m.pos.push_back({radius * std::cos(angle), radius * std::sin(angle), 0.0});
        }
    const auto at = [sectors](int r, int s) { return r * sectors + (s % sectors); };
    for (int r = 0; r < rings; ++r)
        for (int s = 0; s < sectors; ++s) m.quad.push_back({at(r, s), at(r + 1, s), at(r + 1, s + 1), at(r, s + 1)});
    return m;
}

std::vector<char> allFaces(const TestSurface& s) { return std::vector<char>(s.mesh.quad.size(), 1); }

// The uv triangles of each flagged face as quadSplit cuts its 3D corners, with their signed areas: none may be flipped.
int flippedTriangles(const TestSurface& s, const std::vector<char>& faces, const QuadUv& uv) {
    int flipped = 0;
    for (size_t f = 0; f < s.mesh.quad.size(); ++f) {
        if (!faces[f]) continue;
        const std::array<int, 4>& q = s.mesh.quad[f];
        const std::array<Vec3, 4> c = {s.mesh.pos[q[0]], s.mesh.pos[q[1]], s.mesh.pos[q[2]], s.mesh.pos[q[3]]};
        for (const std::array<int, 3>& t : quadSplit(c)) {
            const double area = cross(uv[f][t[1]] - uv[f][t[0]], uv[f][t[2]] - uv[f][t[0]]);
            if (!(area > 0)) ++flipped;
        }
    }
    return flipped;
}

bool allFinite(const QuadUv& uv) {
    for (const std::array<Vec2, 4>& q : uv)
        for (const Vec2& p : q)
            if (!std::isfinite(p.x) || !std::isfinite(p.y)) return false;
    return true;
}

// The singular values of the 2x2 map taking the planar triangle (p0, p1, p2) to (u0, u1, u2), and the sign of its determinant.
struct Stretch {
    double larger = 0, smaller = 0;
    bool orientationKept = false;
};

Stretch stretchOf(const std::array<Vec2, 3>& p, const std::array<Vec2, 3>& u) {
    const Vec2 a = p[1] - p[0], b = p[2] - p[0];
    const Vec2 c = u[1] - u[0], d = u[2] - u[0];
    const double det = a.x * b.y - a.y * b.x;
    // J = [c d] [a b]^-1
    const double j00 = (c.x * b.y - d.x * a.y) / det, j01 = (d.x * a.x - c.x * b.x) / det;
    const double j10 = (c.y * b.y - d.y * a.y) / det, j11 = (d.y * a.x - c.y * b.x) / det;
    const double e = 0.5 * (j00 + j11), f = 0.5 * (j00 - j11), g = 0.5 * (j10 + j01), h = 0.5 * (j10 - j01);
    const double q = std::hypot(e, h), r = std::hypot(f, g);
    return {q + r, std::fabs(q - r), j00 * j11 - j01 * j10 > 0};
}

// While alive, every chart's solver fails at its first iteration.
struct BrokenSolver {
    BrokenSolver() { breakFlattenSolverForTests(true); }
    ~BrokenSolver() { breakFlattenSolverForTests(false); }
};

} // namespace

// The largest deviation of a singular value from 1 over the triangles (0, 1, 2) and (0, 2, 3) of every quad of a planar chart
// (z = 0), and whether every triangle keeps its orientation.
double isometryError(const TestSurface& s, const QuadUv& uv, bool& orientationKept) {
    double worst = 0;
    orientationKept = true;
    for (size_t f = 0; f < s.mesh.quad.size(); ++f) {
        const std::array<int, 4>& q = s.mesh.quad[f];
        for (const std::array<int, 3>& t : {std::array<int, 3>{0, 1, 2}, std::array<int, 3>{0, 2, 3}}) {
            std::array<Vec2, 3> p, u;
            for (int c = 0; c < 3; ++c) {
                p[c] = {s.mesh.pos[q[t[c]]].x, s.mesh.pos[q[t[c]]].y};
                u[c] = uv[f][t[c]];
            }
            const Stretch st = stretchOf(p, u);
            orientationKept = orientationKept && st.orientationKept;
            worst = std::max({worst, std::fabs(st.larger - 1), std::fabs(st.smaller - 1)});
        }
    }
    return worst;
}

// A planar disc flattens isometrically: the rest operators, the start map and the common scale add no distortion, and the
// solver converges to the isometry. It stops once an iteration lowers the energy by less than 1e-4 of it; near the isometry
// the energy exceeds 4 by about 4 times the sum of the squared deviations of the singular values from 1, so what is left is a
// deviation of 1e-3 to 1e-2 (measured 3.5e-3 on the regular 2 x 2 grid, 1.9e-2 on the jittered 6 x 6 disc), far inside the
// stretch bounds [0.75, 1.3]. A tolerance of 1e-5 would lower that about threefold at about twice the iterations on whole imports.
TEST(flatten_planar_disc_is_isometric) {
    struct Disc {
        int n;
        double jitter;
        double tolerance;
    };
    for (const Disc& disc : {Disc{2, 0.0, 1e-2}, Disc{6, 0.25, 3e-2}}) {
        const std::unique_ptr<TestSurface> s = surfaceOf(planarGrid(disc.n, disc.jitter, 7u), std::vector<int>(disc.n * disc.n, 0));
        QuadUv uv(s->mesh.quad.size());
        FlattenReport report;
        std::string err;
        CHECK(flattenFaces(s->graph, std::vector<char>(s->mesh.edge.size(), 0), allFaces(*s), uv, report, &err));
        CHECK(report.charts == 1);
        CHECK(report.stoppedEarly == 0);
        CHECK(report.solves.size() == 1 && report.solves[0].stop == SolveStop::Converged);
        bool orientationKept = false;
        CHECK_NEAR(isometryError(*s, uv, orientationKept), 0.0, disc.tolerance);
        CHECK(orientationKept);
    }
}

TEST(flatten_hemisphere_has_no_flip) {
    const std::unique_ptr<TestSurface> s = hemisphere(4);
    const std::vector<char> faces = allFaces(*s);
    QuadUv uv(s->mesh.quad.size());
    FlattenReport report;
    std::string err;
    CHECK(flattenFaces(s->graph, std::vector<char>(s->mesh.edge.size(), 0), faces, uv, report, &err));
    CHECK(report.charts == 1);
    CHECK(allFinite(uv));
    CHECK(flippedTriangles(*s, faces, uv) == 0);
}

// The chart's vertices are numbered by first appearance in face order, so shuffled faces, vertices and corner starts put the
// boundary vertex numbers out of cycle order; the boundary is still walked as one oriented cycle.
TEST(flatten_boundary_out_of_order) {
    const std::unique_ptr<TestSurface> source = hemisphere(3);
    std::mt19937 random(11u);
    std::vector<int> newVertex(source->mesh.pos.size());
    for (size_t v = 0; v < newVertex.size(); ++v) newVertex[v] = static_cast<int>(v);
    std::shuffle(newVertex.begin(), newVertex.end(), random);
    std::vector<int> faceOrder(source->mesh.quad.size());
    for (size_t f = 0; f < faceOrder.size(); ++f) faceOrder[f] = static_cast<int>(f);
    std::shuffle(faceOrder.begin(), faceOrder.end(), random);
    QuadMesh m;
    m.pos.resize(source->mesh.pos.size());
    for (size_t v = 0; v < newVertex.size(); ++v) m.pos[newVertex[v]] = source->mesh.pos[v];
    for (int f : faceOrder) {
        const std::array<int, 4>& q = source->mesh.quad[f];
        const int turn = static_cast<int>(random() % 4);
        std::array<int, 4> moved;
        for (int k = 0; k < 4; ++k) moved[k] = newVertex[q[(k + turn) % 4]];
        m.quad.push_back(moved);
    }
    const std::unique_ptr<TestSurface> s = surfaceOf(std::move(m), std::vector<int>(faceOrder.size(), 0));
    const std::vector<char> faces = allFaces(*s);
    QuadUv uv(s->mesh.quad.size());
    FlattenReport report;
    std::string err;
    CHECK(flattenFaces(s->graph, std::vector<char>(s->mesh.edge.size(), 0), faces, uv, report, &err));
    CHECK(report.charts == 1);
    CHECK(allFinite(uv));
    CHECK(flippedTriangles(*s, faces, uv) == 0);
}

// An inner vertex moved onto its neighbour: one edge of length zero and two quads with a collapsed side. The mollified rest
// shape gives every triangle a positive area, and the chart flattens without flips.
TEST(flatten_zero_length_edge_mollified) {
    QuadMesh m = planarGrid(6, 0.25, 3u);
    const int moved = 3 * 7 + 3, onto = 3 * 7 + 2;
    m.pos[moved] = m.pos[onto];
    const std::unique_ptr<TestSurface> s = surfaceOf(std::move(m), std::vector<int>(36, 0));
    const std::vector<char> faces = allFaces(*s);
    QuadUv uv(s->mesh.quad.size());
    FlattenReport report;
    std::string err;
    CHECK(flattenFaces(s->graph, std::vector<char>(s->mesh.edge.size(), 0), faces, uv, report, &err));
    CHECK(report.charts == 1);
    CHECK(allFinite(uv));
    CHECK(flippedTriangles(*s, faces, uv) == 0);
}

TEST(flatten_annulus_fails) {
    const std::unique_ptr<TestSurface> s = surfaceOf(planarAnnulus(2, 12), std::vector<int>(24, 0));
    const std::vector<char> faces = allFaces(*s);
    QuadUv uv(s->mesh.quad.size());
    const std::array<Vec2, 4> before = {Vec2{7, 7}, Vec2{7, 7}, Vec2{7, 7}, Vec2{7, 7}};
    uv[0] = before;
    FlattenReport report;
    std::string err;
    CHECK(!flattenFaces(s->graph, std::vector<char>(s->mesh.edge.size(), 0), faces, uv, report, &err));
    CHECK(err.find("chart 0") != std::string::npos);
    CHECK(err.find("not a disc") != std::string::npos);
    CHECK(std::memcmp(&uv[0], &before, sizeof before) == 0);  // nothing written on failure
}

// The solver's first global solve fails: the chart keeps the start map, which is free of flips, and counts as stopped early.
TEST(flatten_forced_solver_failure_keeps_start_map) {
    const std::unique_ptr<TestSurface> s = hemisphere(3);
    const std::vector<char> faces = allFaces(*s);
    QuadUv uv(s->mesh.quad.size());
    FlattenReport report;
    std::string err;
    {
        const BrokenSolver broken;
        CHECK(flattenFaces(s->graph, std::vector<char>(s->mesh.edge.size(), 0), faces, uv, report, &err));
    }
    CHECK(report.charts == 1);
    CHECK(report.stoppedEarly == 1);
    CHECK(allFinite(uv));
    CHECK(flippedTriangles(*s, faces, uv) == 0);
    // the start map puts the boundary (the equator, z = 0) on a circle: all its corners at one distance from the centre of
    // the circle through three of them
    std::vector<Vec2> boundary;
    for (size_t f = 0; f < s->mesh.quad.size(); ++f)
        for (int k = 0; k < 4; ++k)
            if (s->mesh.pos[s->mesh.quad[f][k]].z == 0) boundary.push_back(uv[f][k]);
    CHECK(boundary.size() > 3);
    const Vec2 a = boundary[0];
    Vec2 b = a, c = a;
    for (const Vec2& p : boundary)
        if (norm(p - a) > norm(b - a)) b = p;
    for (const Vec2& p : boundary)
        if (std::fabs(cross(b - a, p - a)) > std::fabs(cross(b - a, c - a))) c = p;
    // circumcentre of a, b, c
    const Vec2 ab = b - a, ac = c - a;
    const double d = 2 * cross(ab, ac);
    const Vec2 centre = a + Vec2{ac.y * dot(ab, ab) - ab.y * dot(ac, ac), ab.x * dot(ac, ac) - ac.x * dot(ab, ab)} * (1.0 / d);
    const double radius = norm(a - centre);
    for (const Vec2& p : boundary) CHECK_NEAR(norm(p - centre), radius, 1e-9 * radius);
}

// Five charts (the cube sides), one face left unflagged: the same bytes from the calling thread and from eight threads at
// once, and the unflagged face keeps its uv.
TEST(flatten_same_bytes_at_any_thread_count) {
    const std::unique_ptr<TestSurface> s = hemisphere(4);
    const std::vector<char> seam = seamsBetweenLabels(*s);
    std::vector<char> faces = allFaces(*s);
    faces[5] = 0;
    QuadUv start(s->mesh.quad.size());
    start[5] = {Vec2{1, 2}, Vec2{3, 4}, Vec2{5, 6}, Vec2{7, 8}};
    QuadUv reference = start;
    FlattenReport report;
    std::string err;
    CHECK(flattenFaces(s->graph, seam, faces, reference, report, &err));
    CHECK(report.charts == 5);
    CHECK(flippedTriangles(*s, faces, reference) == 0);
    CHECK(std::memcmp(&reference[5], &start[5], sizeof start[5]) == 0);

    std::vector<QuadUv> results(8, start);
    std::vector<char> ok(8, 0);
    std::vector<std::thread> threads;
    for (int t = 0; t < 8; ++t)
        threads.emplace_back([&, t]() {
            FlattenReport r;
            std::string e;
            ok[t] = flattenFaces(s->graph, seam, faces, results[t], r, &e) && r.charts == 5;
        });
    for (std::thread& t : threads) t.join();
    for (int t = 0; t < 8; ++t) {
        CHECK(ok[t]);
        CHECK(std::memcmp(results[t].data(), reference.data(), reference.size() * sizeof reference[0]) == 0);
    }
}

namespace {

// A tube of the given radius closed at its top by a half cube-sphere and open at the bottom: the cube-sphere cap of `hemisphere`
// (2n quads a side) above z = 0, and below it `rings` rows of quads of height radius / n going straight down. The only
// boundary is the bottom ring, so Tutte's embedding squeezes the tube's far end and the cap exponentially towards the middle,
// as on a gun barrel ending in its muzzle. One label. Outward winding.
std::unique_ptr<TestSurface> cappedTube(int n, int rings, double radius) {
    QuadMesh m;
    std::map<std::array<int, 3>, int> vertexOf;
    const double rowHeight = radius / n;
    const auto vertexAt = [&](const std::array<int, 3>& p) {
        const std::map<std::array<int, 3>, int>::const_iterator it = vertexOf.find(p);
        if (it != vertexOf.end()) return it->second;
        const int v = static_cast<int>(m.pos.size());
        vertexOf[p] = v;
        if (p[2] >= 0) {
            m.pos.push_back(normalized(Vec3{double(p[0]), double(p[1]), double(p[2])}) * radius);
        } else {
            m.pos.push_back(normalized(Vec3{double(p[0]), double(p[1]), 0.0}) * radius + Vec3{0, 0, p[2] * rowHeight});
        }
        return v;
    };
    struct Side {
        std::array<int, 3> origin, first, second;
        int count1, count2;
    };
    const std::array<Side, 5> sides = {{
        {{-n, -n, n}, {1, 0, 0}, {0, 1, 0}, 2 * n, 2 * n},
        {{n, -n, -rings}, {0, 1, 0}, {0, 0, 1}, 2 * n, n + rings},
        {{-n, -n, -rings}, {0, 0, 1}, {0, 1, 0}, n + rings, 2 * n},
        {{-n, n, -rings}, {0, 0, 1}, {1, 0, 0}, n + rings, 2 * n},
        {{-n, -n, -rings}, {1, 0, 0}, {0, 0, 1}, 2 * n, n + rings},
    }};
    for (const Side& side : sides) {
        const auto corner = [&](int a, int b) {
            std::array<int, 3> p;
            for (int d = 0; d < 3; ++d) p[d] = side.origin[d] + a * side.first[d] + b * side.second[d];
            return vertexAt(p);
        };
        for (int a = 0; a < side.count1; ++a)
            for (int b = 0; b < side.count2; ++b) m.quad.push_back({corner(a, b), corner(a + 1, b), corner(a + 1, b + 1), corner(a, b + 1)});
    }
    const size_t faces = m.quad.size();
    return surfaceOf(std::move(m), std::vector<int>(faces, 0));
}

// A cube-sphere of the given radius (2n quads a side) with one quad missing at the centre of its bottom side: a closed surface
// whose only boundary is that small hole. One label. Outward winding.
std::unique_ptr<TestSurface> sphereWithHole(int n, double radius) {
    QuadMesh m;
    std::map<std::array<int, 3>, int> vertexOf;
    const auto vertexAt = [&](const std::array<int, 3>& p) {
        const std::map<std::array<int, 3>, int>::const_iterator it = vertexOf.find(p);
        if (it != vertexOf.end()) return it->second;
        const int v = static_cast<int>(m.pos.size());
        vertexOf[p] = v;
        m.pos.push_back(normalized(Vec3{double(p[0]), double(p[1]), double(p[2])}) * radius);
        return v;
    };
    // origin, first and second axis (their cross product points outward); the first side is the bottom
    struct Side {
        std::array<int, 3> origin, first, second;
    };
    const std::array<Side, 6> sides = {{
        {{-n, -n, -n}, {0, 1, 0}, {1, 0, 0}},
        {{-n, -n, n}, {1, 0, 0}, {0, 1, 0}},
        {{n, -n, -n}, {0, 1, 0}, {0, 0, 1}},
        {{-n, -n, -n}, {0, 0, 1}, {0, 1, 0}},
        {{-n, n, -n}, {0, 0, 1}, {1, 0, 0}},
        {{-n, -n, -n}, {1, 0, 0}, {0, 0, 1}},
    }};
    for (size_t s = 0; s < sides.size(); ++s) {
        const Side& side = sides[s];
        const auto corner = [&](int a, int b) {
            std::array<int, 3> p;
            for (int d = 0; d < 3; ++d) p[d] = side.origin[d] + a * side.first[d] + b * side.second[d];
            return vertexAt(p);
        };
        for (int a = 0; a < 2 * n; ++a)
            for (int b = 0; b < 2 * n; ++b) {
                if (s == 0 && a == n && b == n) continue;  // the hole
                m.quad.push_back({corner(a, b), corner(a + 1, b), corner(a + 1, b + 1), corner(a, b + 1)});
            }
    }
    const size_t faces = m.quad.size();
    return surfaceOf(std::move(m), std::vector<int>(faces, 0));
}

// The symmetric Dirichlet energy per unit area of the flagged faces' map, from the uv as written (after the common scale):
// sum A (s1^2 + s1^-2 + s2^2 + s2^-2) / sum A over the quadSplit triangles, each laid out in its own plane.
double symmetricDirichletOfUv(const TestSurface& s, const QuadUv& uv) {
    double energy = 0, area = 0;
    for (size_t f = 0; f < s.mesh.quad.size(); ++f) {
        const std::array<int, 4>& q = s.mesh.quad[f];
        const std::array<Vec3, 4> c = {s.mesh.pos[q[0]], s.mesh.pos[q[1]], s.mesh.pos[q[2]], s.mesh.pos[q[3]]};
        for (const std::array<int, 3>& t : quadSplit(c)) {
            const Vec3 e1 = c[t[1]] - c[t[0]], e2 = c[t[2]] - c[t[0]];
            const double length = norm(e1);
            const double triangleArea = 0.5 * norm(cross(e1, e2));
            const std::array<Vec2, 3> plane = {Vec2{0, 0}, Vec2{length, 0}, Vec2{dot(e2, e1) / length, 2 * triangleArea / length}};
            const Stretch st = stretchOf(plane, {uv[f][t[0]], uv[f][t[1]], uv[f][t[2]]});
            energy += triangleArea * (st.larger * st.larger + 1 / (st.larger * st.larger) + st.smaller * st.smaller + 1 / (st.smaller * st.smaller));
            area += triangleArea;
        }
    }
    return energy / area;
}

// Per chart (by lowest face, as flattenFaces numbers them), the uv area of its faces' quadSplit triangles over their 3D area.
std::vector<double> uvOverSurfaceAreaPerChart(const TestSurface& s, const std::vector<int>& chartOfFace, int charts, const QuadUv& uv) {
    std::vector<double> areaUv(charts, 0.0), area3d(charts, 0.0);
    for (size_t f = 0; f < s.mesh.quad.size(); ++f) {
        const std::array<int, 4>& q = s.mesh.quad[f];
        const std::array<Vec3, 4> c = {s.mesh.pos[q[0]], s.mesh.pos[q[1]], s.mesh.pos[q[2]], s.mesh.pos[q[3]]};
        for (const std::array<int, 3>& t : quadSplit(c)) {
            areaUv[chartOfFace[f]] += 0.5 * std::fabs(cross(uv[f][t[1]] - uv[f][t[0]], uv[f][t[2]] - uv[f][t[0]]));
            area3d[chartOfFace[f]] += 0.5 * norm(cross(c[t[1]] - c[t[0]], c[t[2]] - c[t[0]]));
        }
    }
    std::vector<double> ratio(charts);
    for (int k = 0; k < charts; ++k) ratio[k] = areaUv[k] / area3d[k];
    return ratio;
}

// While alive, the solver stops after `cap` iterations.
struct IterationCap {
    explicit IterationCap(int cap) { capFlattenIterationsForTests(cap); }
    ~IterationCap() { capFlattenIterationsForTests(0); }
};

} // namespace

// Two planar discs in one call, of about 1 mm^2 and 1e5 mm^2. The start map is scaled to each chart's own rest area, so with
// the solver broken (each chart keeps its start map) every chart's uv area still equals its 3D area after the common scale;
// and the solver, run from there, keeps it so (both discs flatten isometrically).
TEST(flatten_start_map_scaled_to_the_chart) {
    const double smallSide = 0.25, largeSide = 79.0;  // 16 quads: 1 mm^2 and 99856 mm^2
    QuadMesh m = planarGrid(4, 0.25, 5u);
    const size_t firstCount = m.pos.size();
    const QuadMesh other = planarGrid(4, 0.25, 6u);
    for (Vec3& p : m.pos) p = p * smallSide;
    for (const Vec3& p : other.pos) m.pos.push_back(p * largeSide + Vec3{1000, 0, 0});
    for (const std::array<int, 4>& q : other.quad)
        m.quad.push_back({q[0] + int(firstCount), q[1] + int(firstCount), q[2] + int(firstCount), q[3] + int(firstCount)});
    std::vector<int> chartOfFace(32, 0);
    for (int f = 16; f < 32; ++f) chartOfFace[f] = 1;
    const std::unique_ptr<TestSurface> s = surfaceOf(std::move(m), chartOfFace);
    const std::vector<char> faces = allFaces(*s);
    const std::vector<char> noSeams(s->mesh.edge.size(), 0);

    QuadUv start(s->mesh.quad.size());
    FlattenReport report;
    std::string err;
    {
        const BrokenSolver broken;
        CHECK(flattenFaces(s->graph, noSeams, faces, start, report, &err));
    }
    CHECK(report.charts == 2 && report.stoppedEarly == 2);
    CHECK(flippedTriangles(*s, faces, start) == 0);
    for (double ratio : uvOverSurfaceAreaPerChart(*s, chartOfFace, 2, start)) CHECK_NEAR(ratio, 1.0, 1e-9);

    QuadUv solved(s->mesh.quad.size());
    CHECK(flattenFaces(s->graph, noSeams, faces, solved, report, &err));
    CHECK(report.charts == 2 && report.stoppedEarly == 0 && report.reachedCap == 0);
    for (double ratio : uvOverSurfaceAreaPerChart(*s, chartOfFace, 2, solved)) CHECK_NEAR(ratio, 1.0, 1e-3);
}

// A sphere whose only boundary is one small hole: Tutte's embedding puts the hole on the circle and squeezes the far side of
// the sphere exponentially into the middle. The solve converges from there, without flips, to a map within the stretch bounds
// on half the chart (measured 0.50; flattening a whole sphere cannot do much better). A property check: the fixed 20
// iterations from the unscaled start reached the same here (0.49), so this case does not tell the two apart.
TEST(flatten_sphere_with_hole_within_bounds) {
    const std::unique_ptr<TestSurface> s = sphereWithHole(3, 10.0);
    const std::vector<char> faces = allFaces(*s);
    QuadUv uv(s->mesh.quad.size());
    FlattenReport report;
    std::string err;
    CHECK(flattenFaces(s->graph, std::vector<char>(s->mesh.edge.size(), 0), faces, uv, report, &err));
    CHECK(report.charts == 1 && report.stoppedEarly == 0 && report.reachedCap == 0);
    CHECK(report.solves.size() == 1 && report.solves[0].stop == SolveStop::Converged);
    CHECK(flippedTriangles(*s, faces, uv) == 0);
    const ChartQuality quality = measureCharts(s->mesh, s->label, uv);
    CHECK(quality.unusable[0] == 0);
    CHECK(quality.stretched[0] < 0.6);
}

// A long tube closed by a cap, open at one end (a gun barrel and its muzzle): its map stays outside the stretch bounds
// whatever the solver does (the chart rules must cut it), but from Tutte's exponentially squeezed start the solve needs more
// than 20 iterations to converge (measured 31), and it gets there: stopped as converged, a finite energy, no flip, and the
// energy of the written uv near its optimum (measured 23.5).
TEST(flatten_capped_tube_converges) {
    const std::unique_ptr<TestSurface> s = cappedTube(2, 36, 10.0);
    const std::vector<char> faces = allFaces(*s);
    QuadUv uv(s->mesh.quad.size());
    FlattenReport report;
    std::string err;
    CHECK(flattenFaces(s->graph, std::vector<char>(s->mesh.edge.size(), 0), faces, uv, report, &err));
    CHECK(report.charts == 1 && report.stoppedEarly == 0 && report.reachedCap == 0);
    CHECK(report.solves.size() == 1 && report.solves[0].stop == SolveStop::Converged);
    CHECK(report.solves[0].iterations > 20);
    CHECK(std::isfinite(report.solves[0].energy));
    CHECK(flippedTriangles(*s, faces, uv) == 0);
    CHECK(symmetricDirichletOfUv(*s, uv) < 25.0);
}

// With a cap of 3 iterations the solve stops there, says so, and keeps a flip-free map; a cap is not a failure.
TEST(flatten_iteration_cap_reported) {
    const std::unique_ptr<TestSurface> s = cappedTube(3, 9, 10.0);
    const std::vector<char> faces = allFaces(*s);
    QuadUv uv(s->mesh.quad.size());
    FlattenReport report;
    std::string err;
    {
        const IterationCap cap(3);
        CHECK(flattenFaces(s->graph, std::vector<char>(s->mesh.edge.size(), 0), faces, uv, report, &err));
    }
    CHECK(report.charts == 1 && report.reachedCap == 1 && report.stoppedEarly == 0);
    CHECK(report.solves.size() == 1);
    CHECK(report.solves[0].stop == SolveStop::IterationCap);
    CHECK(report.solves[0].iterations == 3);
    CHECK(std::isfinite(report.solves[0].energy));
    CHECK(allFinite(uv));
    CHECK(flippedTriangles(*s, faces, uv) == 0);
}
