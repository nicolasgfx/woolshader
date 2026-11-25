// Stage 9 of the built-in import: the flagged faces, split along the seams into charts, each flattened into the plane with
// least distortion. Per chart, in this order:
//   1. its quads become triangles, its vertices are split along the seams (one chart vertex per fan of corners)
//   2. it must be a disc: every edge used once or twice, connected, V - E + F = 1, one boundary cycle
//   3. its edge lengths are mollified (Sharp and Crane 2020, section 4.5) so every triangle is safely non-degenerate
//   4. the rest operators: per triangle the gradients of the hat functions, laid out from the lengths alone
//   5. the start map: Tutte's embedding with mean value weights, the boundary on a circle (Tutte 1963; Floater 2003), scaled
//      to the chart's rest area
//   6. SLIM minimises the symmetric Dirichlet energy from there without flips until it converges (Rabinovich et al. 2017,
//      third_party/slim), measured at first against a progressive reference of bounded distortion (after Liu et al. 2018)
// Then all charts of the call get one common scale.
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <limits>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include "import/wrap_stages.h"
#include "model/cancel.h"
#include "model/parallel.h"
#include "slim.h"

namespace knit {

namespace {

constexpr double kMollifyMargin = 1e-3;  // eta = kMollifyMargin * the mean edge length (Sharp and Crane 2020, 4.5)
// SLIM runs until one iteration lowers the energy by less than this share of it. Symmetric Dirichlet is at least 4 (per unit
// area, at an isometry), and SLIM converges linearly near its optimum, so the remaining decrease is a small multiple of the
// last step's: at 1e-4 the map's stretch no longer changes visibly in the chart quality (measured on the tank, skull, car,
// heart and bunny), while most charts stop well before the cap.
constexpr double kEnergyTolerance = 1e-4;
constexpr int kMaxSolverIterations = 200;  // a cap, so that a chart converging very slowly cannot stall the import
// The progressive reference allows singular values in [1 / kReferenceBound, kReferenceBound] (progressiveReference).
constexpr double kReferenceBound = 4.0;

// Tests set this to hand the solver a system it cannot solve (breakFlattenSolverForTests).
std::atomic<bool> solverBrokenForTests{false};
// Tests set this to lower the iteration cap (capFlattenIterationsForTests); 0 means kMaxSolverIterations.
std::atomic<int> iterationCapForTests{0};

// --- 1. a chart as triangles

// Union-find over integers; the root of a set is its lowest member, so the result does not depend on the order of joins.
struct DisjointSets {
    std::vector<int> parent;
    explicit DisjointSets(int count) : parent(count) {
        for (int i = 0; i < count; ++i) parent[i] = i;
    }
    int root(int i) {
        while (parent[i] != i) {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        return i;
    }
    void join(int a, int b) {
        const int ra = root(a);
        const int rb = root(b);
        if (ra < rb) parent[rb] = ra;
        else if (rb < ra) parent[ra] = rb;
    }
};

// One chart's faces as triangles on its own vertices. A surface vertex on a seam is one chart vertex per side.
struct ChartTriangles {
    std::vector<int> faces;                    // ascending
    std::vector<int> vertexOfCorner;           // 4 per face of `faces`: the chart vertex of each corner
    std::vector<Vec3> position;                // per chart vertex, in millimetres
    std::vector<std::array<int, 3>> triangles; // chart vertices; two per face, in face order
};

// Each face's chart (-1 for faces not flagged): flagged faces joined across edges that are no seam, numbered by their lowest
// face. Returns the number of charts.
int labelCharts(const FaceGraph& graph, const std::vector<char>& seam, const std::vector<char>& faces, std::vector<int>& chartOfFace) {
    const int faceCount = graph.faceCount();
    chartOfFace.assign(faceCount, -1);
    int charts = 0;
    std::vector<int> stack;
    for (int start = 0; start < faceCount; ++start) {
        if (!faces[start] || chartOfFace[start] >= 0) continue;
        chartOfFace[start] = charts;
        stack.push_back(start);
        while (!stack.empty()) {
            const int f = stack.back();
            stack.pop_back();
            for (int k = 0; k < 4; ++k) {
                const int e = graph.edgeOfCorner[f][k];
                if (seam[e]) continue;
                for (int g : graph.facesOfEdge[e])
                    if (g >= 0 && faces[g] && chartOfFace[g] < 0) {
                        chartOfFace[g] = charts;
                        stack.push_back(g);
                    }
            }
        }
        ++charts;
    }
    return charts;
}

// The corner of `face` at surface vertex `vertex`, or -1.
int cornerAt(const QuadMesh& mesh, int face, int vertex) {
    for (int k = 0; k < 4; ++k)
        if (mesh.quad[face][k] == vertex) return k;
    return -1;
}

// The chart's triangles. Corners are joined across the chart's edges that are no seam (union-find on corners), and each
// joined fan becomes one chart vertex, numbered in order of first appearance: faces ascending, corners 0 to 3.
ChartTriangles trianglesOfChart(const FaceGraph& graph, const std::vector<char>& seam, std::vector<int> chartFaces) {
    const QuadMesh& mesh = *graph.mesh;
    ChartTriangles chart;
    chart.faces = std::move(chartFaces);
    const int faceCount = static_cast<int>(chart.faces.size());
    const auto localIndexOf = [&](int face) {
        const std::vector<int>::const_iterator it = std::lower_bound(chart.faces.begin(), chart.faces.end(), face);
        return (it != chart.faces.end() && *it == face) ? static_cast<int>(it - chart.faces.begin()) : -1;
    };

    DisjointSets corners(4 * faceCount);
    for (int i = 0; i < faceCount; ++i) {
        const int f = chart.faces[i];
        for (int k = 0; k < 4; ++k) {
            const int e = graph.edgeOfCorner[f][k];
            if (seam[e]) continue;
            const std::array<int, 2>& pair = graph.facesOfEdge[e];
            const int g = pair[0] == f ? pair[1] : pair[0];
            const int j = g >= 0 && g != f ? localIndexOf(g) : -1;
            if (j < 0) continue;
            const int from = mesh.quad[f][k];
            const int to = mesh.quad[f][(k + 1) % 4];
            const int cornerFrom = cornerAt(mesh, g, from);
            const int cornerTo = cornerAt(mesh, g, to);
            if (cornerFrom >= 0) corners.join(4 * i + k, 4 * j + cornerFrom);
            if (cornerTo >= 0) corners.join(4 * i + (k + 1) % 4, 4 * j + cornerTo);
        }
    }

    chart.vertexOfCorner.assign(4 * faceCount, -1);
    std::vector<int> vertexOfRoot(4 * faceCount, -1);
    for (int c = 0; c < 4 * faceCount; ++c) {
        const int r = corners.root(c);
        if (vertexOfRoot[r] < 0) {
            vertexOfRoot[r] = static_cast<int>(chart.position.size());
            chart.position.push_back(mesh.pos[mesh.quad[chart.faces[c / 4]][c % 4]]);
        }
        chart.vertexOfCorner[c] = vertexOfRoot[r];
    }

    chart.triangles.reserve(2 * faceCount);
    for (int i = 0; i < faceCount; ++i) {
        const std::array<int, 4>& q = mesh.quad[chart.faces[i]];
        const std::array<Vec3, 4> c = {mesh.pos[q[0]], mesh.pos[q[1]], mesh.pos[q[2]], mesh.pos[q[3]]};
        for (const std::array<int, 3>& t : quadSplit(c))
            chart.triangles.push_back({chart.vertexOfCorner[4 * i + t[0]], chart.vertexOfCorner[4 * i + t[1]],
                                       chart.vertexOfCorner[4 * i + t[2]]});
    }
    return chart;
}

// --- 2. the disc check

// The chart's edges, and its boundary as one cycle.
struct ChartEdges {
    std::vector<std::array<int, 2>> vertices;  // per edge: (lower, higher) chart vertex
    std::vector<std::array<int, 3>> ofSide;    // per triangle: the edge of side s, from corner s to corner s + 1
    std::vector<int> boundaryCycle;            // boundary vertices in cycle order, interior on the left
};

// A triangle's side as it appears in the triangle, for sorting into edges.
struct TriangleSide {
    int lower = 0, higher = 0;
    int from = 0, to = 0;
    int triangle = 0, side = 0;
};

// Numbers the chart's edges and checks that the chart is a disc: every edge used by one or two triangles, in opposite
// directions when two; the triangles connected; V - E + F = 1; the edges used once form exactly one cycle in which every
// boundary vertex has two boundary edges. False with the reason otherwise.
bool checkDisc(const ChartTriangles& chart, ChartEdges& edges, std::string& reason) {
    const int triangleCount = static_cast<int>(chart.triangles.size());
    const int vertexCount = static_cast<int>(chart.position.size());
    std::vector<TriangleSide> sides;
    sides.reserve(3 * triangleCount);
    for (int t = 0; t < triangleCount; ++t)
        for (int s = 0; s < 3; ++s) {
            const int a = chart.triangles[t][s];
            const int b = chart.triangles[t][(s + 1) % 3];
            if (a == b) {
                reason = "a triangle has two corners on one vertex (a collapsed quad)";
                return false;
            }
            sides.push_back({std::min(a, b), std::max(a, b), a, b, t, s});
        }
    std::sort(sides.begin(), sides.end(), [](const TriangleSide& p, const TriangleSide& q) {
        if (p.lower != q.lower) return p.lower < q.lower;
        if (p.higher != q.higher) return p.higher < q.higher;
        if (p.triangle != q.triangle) return p.triangle < q.triangle;
        return p.side < q.side;
    });

    edges.vertices.clear();
    edges.ofSide.assign(triangleCount, {-1, -1, -1});
    DisjointSets connected(triangleCount);
    std::vector<int> boundaryNext(vertexCount, -1);
    std::vector<int> boundaryDegree(vertexCount, 0);
    int boundaryEdges = 0;
    for (size_t first = 0; first < sides.size();) {
        size_t end = first + 1;
        while (end < sides.size() && sides[end].lower == sides[first].lower && sides[end].higher == sides[first].higher) ++end;
        const int e = static_cast<int>(edges.vertices.size());
        edges.vertices.push_back({sides[first].lower, sides[first].higher});
        for (size_t i = first; i < end; ++i) edges.ofSide[sides[i].triangle][sides[i].side] = e;
        if (end - first > 2) {
            reason = "an edge is used by more than two triangles";
            return false;
        }
        if (end - first == 2) {
            if (sides[first].from == sides[first + 1].from) {
                reason = "two triangles run an edge in the same direction (inconsistent winding)";
                return false;
            }
            connected.join(sides[first].triangle, sides[first + 1].triangle);
        } else {
            const TriangleSide& b = sides[first];
            ++boundaryDegree[b.from];
            ++boundaryDegree[b.to];
            if (boundaryNext[b.from] >= 0) {
                reason = "a boundary vertex has more than two boundary edges";
                return false;
            }
            boundaryNext[b.from] = b.to;
            ++boundaryEdges;
        }
        first = end;
    }

    for (int t = 1; t < triangleCount; ++t)
        if (connected.root(t) != connected.root(0)) {
            reason = "its triangles are not connected";
            return false;
        }

    const long long euler = static_cast<long long>(vertexCount) - static_cast<long long>(edges.vertices.size()) + triangleCount;
    if (euler != 1) {
        char text[160];
        std::snprintf(text, sizeof text, "V - E + F = %lld, not 1: it is not a disc%s", euler,
                      euler < 1 ? " (more than one boundary, or a handle)" : " (closed)");
        reason = text;
        return false;
    }

    int start = -1;
    for (int v = 0; v < vertexCount; ++v) {
        if (boundaryDegree[v] == 0) continue;
        if (boundaryDegree[v] != 2) {
            reason = "a boundary vertex does not have exactly two boundary edges";
            return false;
        }
        if (start < 0) start = v;
    }
    if (start < 0) {
        reason = "it has no boundary";
        return false;
    }
    edges.boundaryCycle.clear();
    int v = start;
    do {
        edges.boundaryCycle.push_back(v);
        v = boundaryNext[v];
    } while (v >= 0 && v != start && static_cast<int>(edges.boundaryCycle.size()) <= boundaryEdges);
    if (v != start || static_cast<int>(edges.boundaryCycle.size()) != boundaryEdges) {
        reason = "its boundary edges do not form exactly one cycle";
        return false;
    }
    return true;
}

// --- 3. the mollified rest shape

// The three sides of a triangle sorted to largest, middle, smallest.
std::array<double, 3> sortedDescending(double a, double b, double c) {
    std::array<double, 3> s = {a, b, c};
    std::sort(s.begin(), s.end(), [](double x, double y) { return x > y; });
    return s;
}

// Intrinsic mollification (Sharp and Crane 2020, section 4.5): every edge grows by the least delta that gives each triangle
// a margin of eta = kMollifyMargin * (mean edge length) on its triangle inequality. Returns the mollified lengths; false
// if a triangle still misses half that margin, which would be a bug.
bool mollifiedLengths(const ChartEdges& edges, const std::vector<double>& length3d, std::vector<double>& mollified, std::string& reason) {
    double sum = 0;
    for (double l : length3d) sum += l;
    const double eta = kMollifyMargin * sum / static_cast<double>(length3d.size());
    double delta = 0;
    for (const std::array<int, 3>& side : edges.ofSide) {
        const std::array<double, 3> l = sortedDescending(length3d[side[0]], length3d[side[1]], length3d[side[2]]);
        delta = std::max(delta, eta + l[0] - l[1] - l[2]);
    }
    mollified.resize(length3d.size());
    for (size_t e = 0; e < length3d.size(); ++e) mollified[e] = length3d[e] + delta;
    for (const std::array<int, 3>& side : edges.ofSide) {
        const std::array<double, 3> l = sortedDescending(mollified[side[0]], mollified[side[1]], mollified[side[2]]);
        if (!(l[1] + l[2] - l[0] >= 0.5 * eta)) {
            reason = "bug: a mollified triangle misses the margin on its triangle inequality";
            return false;
        }
    }
    return true;
}

// --- 4. the rest operators

// A triangle's area from its side lengths by Kahan's stable form of Heron's formula (Kahan 2014).
double kahanArea(double a, double b, double c) {
    const std::array<double, 3> s = sortedDescending(a, b, c);
    const double x = s[0], y = s[1], z = s[2];
    const double product = (x + (y + z)) * (z - (x - y)) * (z + (x - y)) * (x + (y - z));
    return 0.25 * std::sqrt(std::max(0.0, product));
}

// Per triangle, the gradients of its corners' hat functions in its own plane, as rows of Dx and Dy, and its area.
struct RestOperators {
    Eigen::SparseMatrix<double> dx, dy;
    Eigen::VectorXd area;
};

// Each triangle (i, j, k) laid out from its mollified lengths alone: p_i = (0, 0), p_j = (l_ij, 0), p_k by the law of
// cosines above the x axis (height 2A / l_ij). The gradient of corner v's hat function is rot90(e_v) / 2A, e_v the
// opposite edge counter-clockwise.
RestOperators restOperatorsOf(const ChartTriangles& chart, const ChartEdges& edges, const std::vector<double>& mollified) {
    const int triangleCount = static_cast<int>(chart.triangles.size());
    const int vertexCount = static_cast<int>(chart.position.size());
    RestOperators rest;
    rest.area.resize(triangleCount);
    std::vector<Eigen::Triplet<double>> xs, ys;
    xs.reserve(3 * triangleCount);
    ys.reserve(3 * triangleCount);
    for (int t = 0; t < triangleCount; ++t) {
        const double lij = mollified[edges.ofSide[t][0]];
        const double ljk = mollified[edges.ofSide[t][1]];
        const double lki = mollified[edges.ofSide[t][2]];
        const double area = kahanArea(lij, ljk, lki);
        const Vec2 pi = {0, 0};
        const Vec2 pj = {lij, 0};
        const Vec2 pk = {(lij * lij + lki * lki - ljk * ljk) / (2 * lij), 2 * area / lij};
        const std::array<Vec2, 3> opposite = {pk - pj, pi - pk, pj - pi};
        for (int c = 0; c < 3; ++c) {
            const Vec2 g = {-opposite[c].y / (2 * area), opposite[c].x / (2 * area)};
            xs.emplace_back(t, chart.triangles[t][c], g.x);
            ys.emplace_back(t, chart.triangles[t][c], g.y);
        }
        rest.area(t) = area;
    }
    rest.dx.resize(triangleCount, vertexCount);
    rest.dy.resize(triangleCount, vertexCount);
    rest.dx.setFromTriplets(xs.begin(), xs.end());
    rest.dy.setFromTriplets(ys.begin(), ys.end());
    return rest;
}

// --- 5. the start map

// The angle opposite side a in a triangle of sides a, b, c, by the law of cosines.
double angleOpposite(double a, double b, double c) {
    const double cosine = (b * b + c * c - a * a) / (2 * b * c);
    return std::acos(std::clamp(cosine, -1.0, 1.0));
}

// Tutte's embedding (Tutte 1963) with mean value weights (Floater 2003): the boundary on the unit circle at angles
// proportional to its 3D arc length, each inner vertex where sum_j w_ij (u_i - u_j) = 0 with
// w_ij = (tan(alpha / 2) + tan(beta / 2)) / |x_i - x_j|, alpha and beta the angles at i beside edge ij. Free of flips.
bool tutteStartMap(const ChartTriangles& chart, const ChartEdges& edges, const std::vector<double>& mollified, Eigen::MatrixXd& start,
                   std::string& reason) {
    const int vertexCount = static_cast<int>(chart.position.size());
    start = Eigen::MatrixXd::Zero(vertexCount, 2);
    std::vector<int> innerIndex(vertexCount, 0);

    const std::vector<int>& cycle = edges.boundaryCycle;
    const int boundaryCount = static_cast<int>(cycle.size());
    std::vector<double> arc(boundaryCount + 1, 0.0);
    for (int k = 0; k < boundaryCount; ++k) arc[k + 1] = arc[k] + norm(chart.position[cycle[(k + 1) % boundaryCount]] - chart.position[cycle[k]]);
    for (int k = 0; k < boundaryCount; ++k) {
        const double theta = 2 * kPi * arc[k] / arc[boundaryCount];
        start(cycle[k], 0) = std::cos(theta);
        start(cycle[k], 1) = std::sin(theta);
        innerIndex[cycle[k]] = -1;
    }
    int innerCount = 0;
    for (int v = 0; v < vertexCount; ++v)
        if (innerIndex[v] >= 0) innerIndex[v] = innerCount++;
    if (innerCount == 0) return true;

    std::vector<Eigen::Triplet<double>> entries;
    entries.reserve(7 * innerCount);
    Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(innerCount, 2);
    for (size_t t = 0; t < chart.triangles.size(); ++t) {
        const std::array<double, 3> sideLength = {mollified[edges.ofSide[t][0]], mollified[edges.ofSide[t][1]],
                                                  mollified[edges.ofSide[t][2]]};
        for (int c = 0; c < 3; ++c) {
            const int i = chart.triangles[t][c];
            if (innerIndex[i] < 0) continue;
            // corner c lies between side c (to corner c + 1) and side c + 2 (from corner c + 2); side c + 1 is opposite
            const double toNext = sideLength[c];
            const double toPrevious = sideLength[(c + 2) % 3];
            const double halfTan = std::tan(0.5 * angleOpposite(sideLength[(c + 1) % 3], toNext, toPrevious));
            const std::array<int, 2> neighbour = {chart.triangles[t][(c + 1) % 3], chart.triangles[t][(c + 2) % 3]};
            const std::array<double, 2> distance = {toNext, toPrevious};
            for (int n = 0; n < 2; ++n) {
                const double w = halfTan / distance[n];
                const int j = neighbour[n];
                entries.emplace_back(innerIndex[i], innerIndex[i], w);
                if (innerIndex[j] >= 0) {
                    entries.emplace_back(innerIndex[i], innerIndex[j], -w);
                } else {
                    rhs(innerIndex[i], 0) += w * start(j, 0);
                    rhs(innerIndex[i], 1) += w * start(j, 1);
                }
            }
        }
    }
    Eigen::SparseMatrix<double> system(innerCount, innerCount);
    system.setFromTriplets(entries.begin(), entries.end());
    system.makeCompressed();
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(system);
    if (solver.info() != Eigen::Success) {
        reason = "the start map's system cannot be factorised";
        return false;
    }
    const Eigen::MatrixXd inner = solver.solve(rhs);
    if (solver.info() != Eigen::Success || !inner.allFinite()) {
        reason = "the start map's system cannot be solved";
        return false;
    }
    for (int v = 0; v < vertexCount; ++v)
        if (innerIndex[v] >= 0) {
            start(v, 0) = inner(innerIndex[v], 0);
            start(v, 1) = inner(innerIndex[v], 1);
        }
    return true;
}

// Symmetric Dirichlet is not scale invariant, and SLIM without pins has no global scale factor, so a start map far from the
// chart's size starts with a huge energy and huge weights in the global solve. Both Smith and Schaefer 2015 (Bijective
// parameterization with free boundaries) and Rabinovich et al. 2017 (SLIM) start from Tutte's embedding scaled to the
// surface's area: the map is scaled about the origin so that its area equals the chart's rest (mollified) area.
bool scaleToRestArea(const ChartTriangles& chart, const RestOperators& rest, Eigen::MatrixXd& start, std::string& reason) {
    double areaUv = 0;
    for (const std::array<int, 3>& t : chart.triangles) {
        const Vec2 a = {start(t[0], 0), start(t[0], 1)};
        const Vec2 b = {start(t[1], 0), start(t[1], 1)};
        const Vec2 c = {start(t[2], 0), start(t[2], 1)};
        areaUv += 0.5 * cross(b - a, c - a);
    }
    const double areaRest = rest.area.sum();
    if (!(areaUv > 0) || !(areaRest > 0) || !std::isfinite(areaUv) || !std::isfinite(areaRest)) {
        reason = "the start map or the rest shape has no positive area";
        return false;
    }
    start *= std::sqrt(areaRest / areaUv);
    return true;
}

// --- 6. the solve

// The symmetric Dirichlet energy of a map per unit rest area, sum A (s1^2 + s1^-2 + s2^2 + s2^-2) / sum A, with s1, s2 the
// singular values of each triangle's Jacobian J = [Dx u, Dy u; Dx v, Dy v], written as A |J|^2 (1 + 1 / det(J)^2).
// det(J) is the triangle's signed uv area over its rest area, from the uv corners: the flip-avoiding line search needs every
// triangle's corner determinant positive at its start, so a triangle whose corners are collinear in floating point counts as
// degenerate even where the gradients would give a tiny nonzero det(J). Infinite if a triangle is flipped or degenerate.
double symmetricDirichletOf(const ChartTriangles& chart, const RestOperators& rest, const Eigen::MatrixXd& map) {
    const Eigen::VectorXd ux = rest.dx * map.col(0), uy = rest.dy * map.col(0);
    const Eigen::VectorXd vx = rest.dx * map.col(1), vy = rest.dy * map.col(1);
    double energy = 0;
    for (Eigen::Index t = 0; t < ux.size(); ++t) {
        const std::array<int, 3>& corner = chart.triangles[t];
        const Vec2 a = {map(corner[0], 0), map(corner[0], 1)};
        const Vec2 b = {map(corner[1], 0), map(corner[1], 1)};
        const Vec2 c = {map(corner[2], 0), map(corner[2], 1)};
        const double det = 0.5 * cross(b - a, c - a) / rest.area(t);
        if (!(det > 0)) return std::numeric_limits<double>::infinity();
        const double frobenius = ux(t) * ux(t) + uy(t) * uy(t) + vx(t) * vx(t) + vy(t) * vy(t);
        energy += rest.area(t) * frobenius * (1 + 1 / (det * det));
    }
    energy /= rest.area.sum();
    return std::isfinite(energy) ? energy : std::numeric_limits<double>::infinity();
}

// The progressive reference, after Liu et al. 2018 (Progressive Parameterizations). SLIM's weights grow with a triangle's
// distortion, so from a start map with exponentially squeezed triangles (Tutte's embedding of a long chart) the global solve
// spends itself on those few and the line search takes tiny steps. Instead the map is measured against a reference in which
// no triangle is distorted beyond kReferenceBound. Per triangle, with J = U diag(s1, s2) V^T the map's Jacobian against the
// rest shape: where a singular value lies outside [1 / bound, bound], the triangle's reference is moved from its rest shape
// towards the map, R = V diag(s1^a, s2^a) V^T, with the least a that brings the map's singular values against the
// reference, s^(1 - a), within the bound. Gradients on the reference are R^-1 times those on the rest shape, so the returned
// operators are the rows of Dx and Dy mixed by R^-1. False when every triangle is within the bound: then the reference is
// the rest shape itself.
bool progressiveReference(const RestOperators& rest, const Eigen::MatrixXd& map, Eigen::SparseMatrix<double>& dx,
                          Eigen::SparseMatrix<double>& dy) {
    const Eigen::VectorXd ux = rest.dx * map.col(0), uy = rest.dy * map.col(0);
    const Eigen::VectorXd vx = rest.dx * map.col(1), vy = rest.dy * map.col(1);
    const Eigen::Index triangleCount = ux.size();
    const double logBound = std::log(kReferenceBound);
    Eigen::VectorXd r00 = Eigen::VectorXd::Ones(triangleCount);
    Eigen::VectorXd r01 = Eigen::VectorXd::Zero(triangleCount);
    Eigen::VectorXd r11 = Eigen::VectorXd::Ones(triangleCount);
    bool moved = false;
    for (Eigen::Index t = 0; t < triangleCount; ++t) {
        // J^T J = V diag(s1^2, s2^2) V^T, a symmetric 2 x 2 matrix [a b; b c]
        const double a = ux(t) * ux(t) + vx(t) * vx(t);
        const double b = ux(t) * uy(t) + vx(t) * vy(t);
        const double c = uy(t) * uy(t) + vy(t) * vy(t);
        const double det = ux(t) * vy(t) - uy(t) * vx(t);
        const double larger = 0.5 * (a + c) + std::sqrt(0.25 * (a - c) * (a - c) + b * b);
        if (!(det > 0) || !(larger > 0)) continue;
        const double smaller = det * det / larger;  // s1^2 s2^2 = det^2, free of the cancellation in the difference form
        const double logLarger = 0.5 * std::log(larger);
        const double logSmaller = 0.5 * std::log(smaller);
        const double worst = std::max(std::fabs(logLarger), std::fabs(logSmaller));
        if (worst <= logBound) continue;
        const double share = 1 - logBound / worst;  // the exponent a
        // the unit eigenvector e of `larger`, from the longer of the two rows of J^T J - larger I
        double ex = b, ey = larger - a;
        if (std::fabs(larger - c) > std::fabs(ey)) {
            ex = larger - c;
            ey = b;
        }
        const double length = std::hypot(ex, ey);
        if (length > 0) {
            ex /= length;
            ey /= length;
        } else {
            ex = 1;
            ey = 0;
        }
        // R^-1 = s2^-a I + (s1^-a - s2^-a) e e^T
        const double inverseSmaller = std::exp(-share * logSmaller);
        const double inverseLarger = std::exp(-share * logLarger);
        r00(t) = inverseSmaller + (inverseLarger - inverseSmaller) * ex * ex;
        r01(t) = (inverseLarger - inverseSmaller) * ex * ey;
        r11(t) = inverseSmaller + (inverseLarger - inverseSmaller) * ey * ey;
        moved = true;
    }
    if (!moved) return false;
    dx = r00.asDiagonal() * rest.dx + r01.asDiagonal() * rest.dy;
    dy = r01.asDiagonal() * rest.dx + r11.asDiagonal() * rest.dy;
    return true;
}

// One SLIM iteration (Rabinovich et al. 2017) on the symmetric Dirichlet energy against the reference given by `dx`, `dy`
// and the face areas: local step, global solve, flip-avoiding line search. False if the global solve failed, reported or
// thrown; then `map` is unchanged.
bool slimIteration(Eigen::MatrixXi& faces, Eigen::SparseMatrix<double>& dx, Eigen::SparseMatrix<double>& dy, Eigen::VectorXd& area,
                   Eigen::MatrixXd& map) {
    slim::SLIMData data;
    try {
        slim::slim_precompute_rest(faces, dx, dy, area, map, data, slim::SLIMData::SYMMETRIC_DIRICHLET);
        slim::slim_solve(data, 1);
    } catch (...) {
        return false;
    }
    if (data.failed_iteration >= 0 || data.V_o.rows() != map.rows() || data.V_o.cols() != 2) return false;
    map = data.V_o;
    return true;
}

// SLIM from the start map, one iteration at a time: against the progressive reference while one is needed, against the rest
// shape from then on. Once an iteration against a moved reference lowers the true energy by less than kEnergyTolerance of
// it, the reference stays at the rest shape for good; against the rest shape the same rule means converged. A failed global
// solve, or a map whose energy is not finite, stops the solve with the last map of finite energy, which is free of flips.
ChartSolve solveSlim(const ChartTriangles& chart, const RestOperators& rest, Eigen::MatrixXd& map) {
    const int triangleCount = static_cast<int>(chart.triangles.size());
    Eigen::MatrixXi faces(triangleCount, 3);
    for (int t = 0; t < triangleCount; ++t)
        for (int c = 0; c < 3; ++c) faces(t, c) = chart.triangles[t][c];
    Eigen::SparseMatrix<double> restDx = rest.dx, restDy = rest.dy;
    Eigen::VectorXd area = rest.area;
    if (solverBrokenForTests.load()) area.setConstant(std::nan(""));
    const int cap = iterationCapForTests.load() > 0 ? iterationCapForTests.load() : kMaxSolverIterations;

    ChartSolve solve;
    solve.energy = symmetricDirichletOf(chart, rest, map);
    if (!std::isfinite(solve.energy)) {
        solve.stop = SolveStop::NonFiniteEnergy;
        return solve;
    }
    bool referenceAtRest = false;
    Eigen::SparseMatrix<double> referenceDx, referenceDy;
    for (;;) {
        if (solve.iterations == cap) {
            solve.stop = SolveStop::IterationCap;
            return solve;
        }
        cancelCheckpoint();
        const bool moved = !referenceAtRest && progressiveReference(rest, map, referenceDx, referenceDy);
        Eigen::MatrixXd next = map;
        const bool solved = moved ? slimIteration(faces, referenceDx, referenceDy, area, next)
                                  : slimIteration(faces, restDx, restDy, area, next);
        ++solve.iterations;
        if (!solved) {
            solve.stop = SolveStop::SolverFailed;
            return solve;
        }
        const double energy = symmetricDirichletOf(chart, rest, next);
        if (!std::isfinite(energy)) {
            solve.stop = SolveStop::NonFiniteEnergy;
            return solve;
        }
        map = next;
        const bool converged = !(solve.energy - energy > kEnergyTolerance * solve.energy);
        solve.energy = energy;
        if (!converged) continue;
        if (!moved) {
            solve.stop = SolveStop::Converged;
            return solve;
        }
        referenceAtRest = true;
    }
}

double triangleArea3d(const ChartTriangles& chart, const std::array<int, 3>& t) {
    return 0.5 * norm(cross(chart.position[t[1]] - chart.position[t[0]], chart.position[t[2]] - chart.position[t[0]]));
}

double signedAreaUv(const std::vector<Vec2>& uv, const std::array<int, 3>& t) {
    return 0.5 * cross(uv[t[1]] - uv[t[0]], uv[t[2]] - uv[t[0]]);
}

// --- 7. one chart, start to end

struct FlattenedChart {
    ChartTriangles triangles;
    std::vector<Vec2> uv;   // per chart vertex
    double area3d = 0;      // the true 3D area of its triangles
    double areaUv = 0;
    ChartSolve solve;
    std::string error;      // empty on success
};

// Flattens one chart. On failure, `error` says why (without the chart's name).
FlattenedChart flattenChart(const FaceGraph& graph, const std::vector<char>& seam, std::vector<int> chartFaces) {
    FlattenedChart result;
    result.triangles = trianglesOfChart(graph, seam, std::move(chartFaces));
    const ChartTriangles& chart = result.triangles;

    ChartEdges edges;
    if (!checkDisc(chart, edges, result.error)) return result;
    std::vector<double> length3d(edges.vertices.size());
    for (size_t e = 0; e < edges.vertices.size(); ++e)
        length3d[e] = norm(chart.position[edges.vertices[e][1]] - chart.position[edges.vertices[e][0]]);
    double boundaryLength = 0;
    for (size_t k = 0; k < edges.boundaryCycle.size(); ++k)
        boundaryLength += norm(chart.position[edges.boundaryCycle[(k + 1) % edges.boundaryCycle.size()]] - chart.position[edges.boundaryCycle[k]]);
    if (!(boundaryLength > 0) || !std::isfinite(boundaryLength)) {
        result.error = "its boundary has no positive length";
        return result;
    }

    std::vector<double> mollified;
    if (!mollifiedLengths(edges, length3d, mollified, result.error)) return result;
    const RestOperators rest = restOperatorsOf(chart, edges, mollified);
    Eigen::MatrixXd map;
    if (!tutteStartMap(chart, edges, mollified, map, result.error)) return result;
    if (!scaleToRestArea(chart, rest, map, result.error)) return result;
    result.solve = solveSlim(chart, rest, map);
    if (!map.allFinite()) {
        result.error = "bug: the flattened map is not finite";
        return result;
    }

    result.uv.resize(chart.position.size());
    for (size_t v = 0; v < chart.position.size(); ++v) result.uv[v] = {map(v, 0), map(v, 1)};
    for (const std::array<int, 3>& t : chart.triangles) {
        result.area3d += triangleArea3d(chart, t);
        result.areaUv += std::fabs(signedAreaUv(result.uv, t));
    }
    return result;
}

bool fail(std::string* err, const std::string& message) {
    if (err) *err = message;
    return false;
}

} // namespace

// Test seam: while set, every chart's solver gets face areas that are not numbers, so its first global solve fails.
void breakFlattenSolverForTests(bool broken) { solverBrokenForTests.store(broken); }

void capFlattenIterationsForTests(int cap) { iterationCapForTests.store(cap); }

bool flattenFaces(const FaceGraph& graph, const std::vector<char>& seam, const std::vector<char>& faces, QuadUv& uv,
                  FlattenReport& report, std::string* err) {
    report = FlattenReport{};
    if (!graph.mesh) return fail(err, "flatten: the face graph has no mesh");
    const int faceCount = graph.faceCount();
    if (faceCount != static_cast<int>(graph.mesh->quad.size()) || static_cast<int>(graph.edgeOfCorner.size()) != faceCount)
        return fail(err, "flatten: the face graph does not match its mesh");
    if (static_cast<int>(faces.size()) != faceCount || static_cast<int>(uv.size()) != faceCount)
        return fail(err, "flatten: one face flag and one uv quad per face are needed");
    if (seam.size() != graph.facesOfEdge.size()) return fail(err, "flatten: one seam flag per edge is needed");

    std::vector<int> chartOfFace;
    const int chartCount = labelCharts(graph, seam, faces, chartOfFace);
    std::vector<std::vector<int>> facesOfChart(chartCount);
    for (int f = 0; f < faceCount; ++f)
        if (chartOfFace[f] >= 0) facesOfChart[chartOfFace[f]].push_back(f);

    // Each chart is flattened on its own and written to its own slot: the result does not depend on the thread count.
    std::vector<FlattenedChart> flattened(chartCount);
    ExecPolicy policy;
    policy.parallel = true;
    parallelFor(chartCount, 1, [&](int64_t begin, int64_t end, int64_t) {
        for (int64_t c = begin; c < end; ++c) flattened[c] = flattenChart(graph, seam, std::move(facesOfChart[c]));
    }, policy);

    double area3d = 0, areaUv = 0;
    for (int c = 0; c < chartCount; ++c) {
        const FlattenedChart& chart = flattened[c];
        if (!chart.error.empty()) {
            char name[96];
            std::snprintf(name, sizeof name, "flatten: chart %d (lowest face %d, %d faces): ", c, chart.triangles.faces.front(),
                          static_cast<int>(chart.triangles.faces.size()));
            return fail(err, name + chart.error);
        }
        area3d += chart.area3d;
        areaUv += chart.areaUv;
    }
    if (chartCount > 0 && !(areaUv > 0 && area3d > 0 && std::isfinite(area3d) && std::isfinite(areaUv)))
        return fail(err, "flatten: the charts have no positive area");

    // One scale for all charts: their uv area equals their 3D area. Each chart is scaled about its uv bounding box's centre.
    const double scale = chartCount > 0 ? std::sqrt(area3d / areaUv) : 1.0;
    for (int c = 0; c < chartCount; ++c) {
        const FlattenedChart& chart = flattened[c];
        Vec2 low = chart.uv[0], high = chart.uv[0];
        for (const Vec2& p : chart.uv) {
            low = {std::min(low.x, p.x), std::min(low.y, p.y)};
            high = {std::max(high.x, p.x), std::max(high.y, p.y)};
        }
        const Vec2 centre = (low + high) * 0.5;
        const std::vector<int>& chartFaces = chart.triangles.faces;
        for (size_t i = 0; i < chartFaces.size(); ++i)
            for (int k = 0; k < 4; ++k)
                uv[chartFaces[i]][k] = centre + (chart.uv[chart.triangles.vertexOfCorner[4 * i + k]] - centre) * scale;
        const SolveStop stop = chart.solve.stop;
        report.stoppedEarly += stop == SolveStop::SolverFailed || stop == SolveStop::NonFiniteEnergy ? 1 : 0;
        report.reachedCap += stop == SolveStop::IterationCap ? 1 : 0;
        report.solves.push_back(chart.solve);
    }
    report.charts = chartCount;
    return true;
}

} // namespace knit
