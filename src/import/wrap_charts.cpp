// The chart stage, in doubles. A heap of (distance, face) pairs breaks equal distances by the lower face index, and every
// "farthest" or "largest" is the first of equals: the same charts on every machine.
#include "import/wrap_stages.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <queue>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include "model/cancel.h"

namespace knit {

namespace {

constexpr int kGrowRounds = 5;
constexpr double kStretchLow = 0.75;     // what knitting accepts (src/knitting/layout/uv_atlas.h)
constexpr double kStretchHigh = 1.3;
constexpr double kTargetOutside = 0.03;  // area share allowed outside those bounds
constexpr double kMaxChartOutside = 0.15;
constexpr double kMaxUnusable = 0.005;   // a chart folding or overlapping more than this is split
constexpr double kSplitStretched = 0.5;  // a chart with more than this outside the bounds is split, not slit
constexpr double kOverlapPixelMm = 1.5;  // points of a flattened chart are written into pixels of this size,
constexpr double kOverlapSpreadMm = 8.0; // and surface points this far apart in one pixel mean two sheets
constexpr double kDegenerateDet = 1e-14;

struct WalkSource {
    int face = 0;
    int label = 0;
};

struct Walk {
    std::vector<double> distance;
    std::vector<int> label;
};

using FaceQueue = std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>>;

// Dijkstra from labelled faces. `inside`, if given, limits the walk. `directions`, if given, holds each label's mean normal, and
// faces turning away from it cost more, which keeps a chart from folding around a rim onto a plate's other side.
Walk walkFaces(const FaceGraph& graph, const std::vector<WalkSource>& sources, const std::vector<char>* inside,
               const std::vector<Vec3>* directions) {
    Walk walk;
    walk.distance.assign(graph.faceCount(), std::numeric_limits<double>::infinity());
    walk.label.assign(graph.faceCount(), -1);
    std::vector<std::pair<double, int>> first;
    for (const WalkSource& source : sources) {
        walk.distance[source.face] = 0.0;
        walk.label[source.face] = source.label;
        first.push_back({0.0, source.face});
    }
    FaceQueue queue(std::greater<std::pair<double, int>>(), std::move(first));
    while (!queue.empty()) {
        const auto [d, f] = queue.top();
        queue.pop();
        if (d > walk.distance[f]) continue;
        for (int k = graph.start[f]; k < graph.start[f + 1]; ++k) {
            const FaceNeighbour& next = graph.neighbour[k];
            if (inside && !(*inside)[next.face]) continue;
            double cost = next.cost;
            if (directions) cost *= 1.0 + kTurnCost * (1.0 - dot(graph.normal[next.face], (*directions)[walk.label[f]]));
            if (d + cost < walk.distance[next.face]) {
                walk.distance[next.face] = d + cost;
                walk.label[next.face] = walk.label[f];
                queue.push({d + cost, next.face});
            }
        }
    }
    return walk;
}

int firstLargest(const std::vector<double>& values) {
    return static_cast<int>(std::max_element(values.begin(), values.end()) - values.begin());
}

std::vector<char> facesOfChart(const std::vector<int>& label, int chart) {
    std::vector<char> inside(label.size());
    for (size_t f = 0; f < label.size(); ++f) inside[f] = label[f] == chart;
    return inside;
}

// The chart's face farthest from its border.
int centreOfChart(const FaceGraph& graph, const std::vector<int>& label, int chart) {
    const std::vector<char> inside = facesOfChart(label, chart);
    std::vector<WalkSource> border;
    int firstFace = -1;
    for (int f = 0; f < graph.faceCount(); ++f) {
        if (!inside[f]) continue;
        if (firstFace < 0) firstFace = f;
        for (int k = graph.start[f]; k < graph.start[f + 1]; ++k)
            if (label[graph.neighbour[k].face] != chart) {
                border.push_back({f, chart});
                break;
            }
    }
    if (border.empty()) return firstFace;
    Walk walk = walkFaces(graph, border, &inside, nullptr);
    for (int f = 0; f < graph.faceCount(); ++f)
        if (!inside[f]) walk.distance[f] = -1.0;
    return firstLargest(walk.distance);
}

// The chart's edges that a breadth-first spanning tree of its faces does not cross, as flags.
std::vector<char> edgesTheTreeLeaves(const FaceGraph& graph, const std::vector<int>& label, int chart) {
    const int root = centreOfChart(graph, label, chart);
    std::vector<char> visited(graph.faceCount(), 0);
    std::vector<char> crossed(graph.facesOfEdge.size(), 0);
    visited[root] = 1;
    std::vector<int> frontier = {root};
    while (!frontier.empty()) {  // breadth first: fronts meet far from the root, near the borders
        std::vector<int> following;
        for (int f : frontier)
            for (int k = graph.start[f]; k < graph.start[f + 1]; ++k) {
                const FaceNeighbour& next = graph.neighbour[k];
                if (label[next.face] != chart || visited[next.face]) continue;
                visited[next.face] = 1;
                crossed[next.edge] = 1;
                following.push_back(next.face);
            }
        frontier.swap(following);
    }
    return crossed;
}

// Adds to `seam` what is left of a chart's inner edges after its dangling branches are pruned: the cuts its tunnels need.
void addTunnelCuts(const FaceGraph& graph, const std::vector<int>& label, int chart, std::vector<char>& seam) {
    const std::vector<char> crossed = edgesTheTreeLeaves(graph, label, chart);
    std::vector<int> kept, border;
    std::vector<char> seen(graph.facesOfEdge.size(), 0);
    for (int e = 0; e < static_cast<int>(graph.facesOfEdge.size()); ++e) {
        const bool first = label[graph.facesOfEdge[e][0]] == chart;
        const bool second = label[graph.facesOfEdge[e][1]] == chart;
        if (first && second && !crossed[e]) kept.push_back(e);
        else if (first != second) border.push_back(e);
    }
    std::unordered_map<int, int> degree;
    std::unordered_map<int, std::vector<int>> ends;
    const auto addEnds = [&](int e) {
        for (int v : graph.mesh->edge[e]) {
            ++degree[v];
            ends[v].push_back(e);
        }
    };
    for (int e : kept) addEnds(e);
    for (int e : border) addEnds(e);
    std::vector<char> stillKept(graph.facesOfEdge.size(), 0);
    for (int e : kept) stillKept[e] = 1;
    std::vector<int> leaves;
    for (const std::pair<const int, int>& entry : degree)
        if (entry.second == 1) leaves.push_back(entry.first);
    while (!leaves.empty()) {  // only inner edges go: the border stays, so the order of pruning does not matter
        const int v = leaves.back();
        leaves.pop_back();
        for (int e : ends[v]) {
            if (!stillKept[e]) continue;
            stillKept[e] = 0;
            for (int w : graph.mesh->edge[e])
                if (--degree[w] == 1) leaves.push_back(w);
        }
    }
    for (int e : kept)
        if (stillKept[e]) seam[e] = 1;
}

std::vector<int> chartsInUse(const std::vector<int>& label) {
    std::vector<int> charts = label;
    std::sort(charts.begin(), charts.end());
    charts.erase(std::unique(charts.begin(), charts.end()), charts.end());
    return charts;
}

// The triangles lying over another part of their own chart: four points of each flattened triangle go into pixels, and two
// distant surface points in one pixel mean two sheets. `flat` is in the chart's isometric plane, three corners per triangle.
std::vector<char> trianglesOverlapping(const std::vector<std::array<Vec3, 3>>& surface, const std::vector<std::array<Vec2, 3>>& flat) {
    static const double kWeight[4][3] = {{1.0 / 3, 1.0 / 3, 1.0 / 3}, {0.8, 0.1, 0.1}, {0.1, 0.8, 0.1}, {0.1, 0.1, 0.8}};
    double lowU = std::numeric_limits<double>::infinity(), lowV = lowU;
    std::vector<Vec2> plane(flat.size() * 4);
    std::vector<Vec3> point(flat.size() * 4);
    for (size_t t = 0; t < flat.size(); ++t)
        for (int s = 0; s < 4; ++s) {
            const double* w = kWeight[s];
            plane[t * 4 + s] = {w[0] * flat[t][0].x + w[1] * flat[t][1].x + w[2] * flat[t][2].x,
                                w[0] * flat[t][0].y + w[1] * flat[t][1].y + w[2] * flat[t][2].y};
            point[t * 4 + s] = {w[0] * surface[t][0].x + w[1] * surface[t][1].x + w[2] * surface[t][2].x,
                                w[0] * surface[t][0].y + w[1] * surface[t][1].y + w[2] * surface[t][2].y,
                                w[0] * surface[t][0].z + w[1] * surface[t][1].z + w[2] * surface[t][2].z};
            lowU = std::min(lowU, plane[t * 4 + s].x);
            lowV = std::min(lowV, plane[t * 4 + s].y);
        }
    struct Pixel {
        Vec3 low, high;
    };
    std::unordered_map<int64_t, Pixel> pixels;
    std::vector<int64_t> keyOfSample(plane.size());
    for (size_t i = 0; i < plane.size(); ++i) {
        const int64_t key = static_cast<int64_t>(std::floor((plane[i].x - lowU) / kOverlapPixelMm)) * 1000003 +
                            static_cast<int64_t>(std::floor((plane[i].y - lowV) / kOverlapPixelMm));
        keyOfSample[i] = key;
        const Vec3& p = point[i];
        const auto found = pixels.find(key);
        if (found == pixels.end()) {
            pixels.insert({key, Pixel{p, p}});
            continue;
        }
        Pixel& pixel = found->second;
        pixel.low = {std::min(pixel.low.x, p.x), std::min(pixel.low.y, p.y), std::min(pixel.low.z, p.z)};
        pixel.high = {std::max(pixel.high.x, p.x), std::max(pixel.high.y, p.y), std::max(pixel.high.z, p.z)};
    }
    std::vector<char> over(flat.size(), 0);
    for (size_t i = 0; i < plane.size(); ++i) {
        const Pixel& pixel = pixels.at(keyOfSample[i]);
        const double spread = std::max({pixel.high.x - pixel.low.x, pixel.high.y - pixel.low.y, pixel.high.z - pixel.low.z});
        if (spread > kOverlapSpreadMm) over[i / 4] = 1;
    }
    return over;
}

// A cut from the seams to the stretched point of `chart` farthest from them. Stretch gathers around cone-like tips, and a short
// cut there lets the chart relax without adding a chart. Returns the new seam edges.
std::vector<int> slitFromSeam(const FaceGraph& graph, const std::vector<int>& label, int chart, const std::vector<int>& stretchedFaces,
                              const std::vector<char>& seam) {
    const QuadMesh& mesh = *graph.mesh;
    const double unreached = std::numeric_limits<double>::infinity();
    std::vector<double> distance(mesh.pos.size(), unreached);
    std::vector<int> parentVertex(mesh.pos.size(), -1);
    std::vector<int> parentEdge(mesh.pos.size(), -1);
    std::vector<std::pair<double, int>> first;
    for (int f = 0; f < graph.faceCount(); ++f) {
        if (label[f] != chart) continue;
        for (int v : mesh.quad[f]) {
            if (distance[v] == 0.0) continue;
            for (int k = graph.edgeStart[v]; k < graph.edgeStart[v + 1]; ++k)
                if (seam[graph.edgeOfVertex[k]]) {
                    distance[v] = 0.0;
                    first.push_back({0.0, v});
                    break;
                }
        }
    }
    FaceQueue queue(std::greater<std::pair<double, int>>(), std::move(first));
    while (!queue.empty()) {
        const auto [d, v] = queue.top();
        queue.pop();
        if (d > distance[v]) continue;
        for (int k = graph.edgeStart[v]; k < graph.edgeStart[v + 1]; ++k) {
            const int e = graph.edgeOfVertex[k];
            if (label[graph.facesOfEdge[e][0]] != chart || label[graph.facesOfEdge[e][1]] != chart) continue;
            const int w = mesh.edge[e][0] == v ? mesh.edge[e][1] : mesh.edge[e][0];
            const double reached = d + graph.edgeLength[e];
            if (reached < distance[w]) {
                distance[w] = reached;
                parentVertex[w] = v;
                parentEdge[w] = e;
                queue.push({reached, w});
            }
        }
    }
    int tip = -1;
    double farthest = -std::numeric_limits<double>::infinity();
    for (int f : stretchedFaces)
        for (int v : mesh.quad[f]) {
            const double d = distance[v] == unreached ? -1.0 : distance[v];
            if (d > farthest) {
                farthest = d;
                tip = v;
            }
        }
    std::vector<int> cut;
    while (parentVertex[tip] >= 0) {
        cut.push_back(parentEdge[tip]);
        tip = parentVertex[tip];
    }
    return cut;
}

// Divides a chart between its two most distant faces. The second half gets a new number, which is returned.
int splitChart(const FaceGraph& graph, std::vector<int>& label, int chart) {
    const std::vector<char> inside = facesOfChart(label, chart);
    int first = static_cast<int>(std::find(inside.begin(), inside.end(), char(1)) - inside.begin());
    int second = first;
    for (int pass = 0; pass < 2; ++pass) {
        Walk walk = walkFaces(graph, {{first, 0}}, &inside, nullptr);
        for (double& d : walk.distance)
            if (!std::isfinite(d)) d = -1.0;
        second = first;
        first = firstLargest(walk.distance);
    }
    const Walk halves = walkFaces(graph, {{first, 0}, {second, 1}}, &inside, nullptr);
    const int other = *std::max_element(label.begin(), label.end()) + 1;
    for (int f = 0; f < graph.faceCount(); ++f)
        if (inside[f] && halves.label[f] == 1) label[f] = other;
    return other;
}

// The chart of the first entry with the largest value.
int chartWithLargest(const std::vector<int>& charts, const std::vector<double>& values) {
    return charts[firstLargest(values)];
}

} // namespace

FaceGraph buildFaceGraph(const QuadMesh& mesh) {
    FaceGraph graph;
    graph.mesh = &mesh;
    const int faces = static_cast<int>(mesh.quad.size());
    graph.area.resize(faces);
    graph.normal.resize(faces);
    std::vector<Vec3> centre(faces);
    for (int f = 0; f < faces; ++f) {
        // The diagonals' cross product: half its length is the area, its direction the normal. The centre is the corners' mean.
        const std::array<int, 4>& q = mesh.quad[f];
        const Vec3 diagonals = cross(mesh.pos[q[0]] - mesh.pos[q[2]], mesh.pos[q[1]] - mesh.pos[q[3]]);
        graph.area[f] = 0.5 * norm(diagonals);
        graph.normal[f] = normalized(diagonals);
        centre[f] = (mesh.pos[q[0]] + mesh.pos[q[1]] + mesh.pos[q[2]] + mesh.pos[q[3]]) * 0.25;
    }

    std::unordered_map<uint64_t, int> edgeOfPair;
    edgeOfPair.reserve(mesh.edge.size() * 2);
    for (int e = 0; e < static_cast<int>(mesh.edge.size()); ++e)
        edgeOfPair[(uint64_t(uint32_t(mesh.edge[e][0])) << 32) | uint32_t(mesh.edge[e][1])] = e;
    graph.facesOfEdge.assign(mesh.edge.size(), {-1, -1});
    graph.edgeOfCorner.resize(faces);
    for (int f = 0; f < faces; ++f)
        for (int k = 0; k < 4; ++k) {
            const uint32_t a = mesh.quad[f][k];
            const uint32_t b = mesh.quad[f][(k + 1) % 4];
            const int e = edgeOfPair.at((uint64_t(std::min(a, b)) << 32) | std::max(a, b));
            graph.edgeOfCorner[f][k] = e;
            std::array<int, 2>& pair = graph.facesOfEdge[e];
            pair[pair[0] < 0 ? 0 : 1] = f;
        }

    graph.start.assign(faces + 1, 0);
    for (const std::array<int, 2>& pair : graph.facesOfEdge) {
        ++graph.start[pair[0] + 1];
        ++graph.start[pair[1] + 1];
    }
    for (int f = 0; f < faces; ++f) graph.start[f + 1] += graph.start[f];
    graph.neighbour.resize(graph.start[faces]);
    std::vector<int> fill(graph.start.begin(), graph.start.end() - 1);
    for (int e = 0; e < static_cast<int>(graph.facesOfEdge.size()); ++e) {
        const int a = graph.facesOfEdge[e][0];
        const int b = graph.facesOfEdge[e][1];
        const double crease = 1.0 - dot(graph.normal[a], graph.normal[b]);
        const double cost = norm(centre[a] - centre[b]) * (1.0 + kCreaseCost * crease);
        graph.neighbour[fill[a]++] = {b, cost, e};
        graph.neighbour[fill[b]++] = {a, cost, e};
    }

    graph.edgeStart.assign(mesh.pos.size() + 1, 0);
    graph.edgeLength.resize(mesh.edge.size());
    for (size_t e = 0; e < mesh.edge.size(); ++e) {
        ++graph.edgeStart[mesh.edge[e][0] + 1];
        ++graph.edgeStart[mesh.edge[e][1] + 1];
        graph.edgeLength[e] = norm(mesh.pos[mesh.edge[e][0]] - mesh.pos[mesh.edge[e][1]]);
    }
    for (size_t v = 0; v < mesh.pos.size(); ++v) graph.edgeStart[v + 1] += graph.edgeStart[v];
    graph.edgeOfVertex.resize(graph.edgeStart.back());
    std::vector<int> fillEdges(graph.edgeStart.begin(), graph.edgeStart.end() - 1);
    for (int e = 0; e < static_cast<int>(mesh.edge.size()); ++e) {
        graph.edgeOfVertex[fillEdges[mesh.edge[e][0]]++] = e;
        graph.edgeOfVertex[fillEdges[mesh.edge[e][1]]++] = e;
    }
    return graph;
}

std::vector<int> spreadSeeds(const FaceGraph& graph, int count) {
    std::vector<int> seeds = {firstLargest(walkFaces(graph, {{0, 0}}, nullptr, nullptr).distance)};
    while (static_cast<int>(seeds.size()) < count) {
        cancelCheckpoint();
        std::vector<WalkSource> sources;
        for (int s : seeds) sources.push_back({s, 0});
        seeds.push_back(firstLargest(walkFaces(graph, sources, nullptr, nullptr).distance));
    }
    return seeds;
}

std::vector<int> growCharts(const FaceGraph& graph, std::vector<int> seeds) {
    std::vector<Vec3> directions;
    std::vector<int> label;
    for (int round = 0; round < kGrowRounds; ++round) {
        cancelCheckpoint();
        std::vector<int> distinct;  // two seeds on one face would leave a patch empty
        for (int s : seeds)
            if (std::find(distinct.begin(), distinct.end(), s) == distinct.end()) distinct.push_back(s);
        seeds.swap(distinct);
        std::vector<WalkSource> sources;
        for (int i = 0; i < static_cast<int>(seeds.size()); ++i) sources.push_back({seeds[i], i});
        label = walkFaces(graph, sources, nullptr, directions.empty() ? nullptr : &directions).label;

        directions.assign(seeds.size(), Vec3());
        for (int f = 0; f < graph.faceCount(); ++f) {
            Vec3& mean = directions[label[f]];
            mean = {mean.x + graph.normal[f].x * graph.area[f], mean.y + graph.normal[f].y * graph.area[f],
                    mean.z + graph.normal[f].z * graph.area[f]};
        }
        for (Vec3& mean : directions) {
            const double length = std::max(norm(mean), 1e-9);
            mean = {mean.x / length, mean.y / length, mean.z / length};
        }
        for (int i = 0; i < static_cast<int>(seeds.size()); ++i) seeds[i] = centreOfChart(graph, label, i);
    }
    return label;
}

std::vector<char> cutToDiscs(const FaceGraph& graph, const std::vector<int>& label) {
    std::vector<char> seam(graph.facesOfEdge.size(), 0);
    for (size_t e = 0; e < graph.facesOfEdge.size(); ++e)
        seam[e] = label[graph.facesOfEdge[e][0]] != label[graph.facesOfEdge[e][1]];
    for (int chart : chartsInUse(label)) {
        cancelCheckpoint();
        addTunnelCuts(graph, label, chart, seam);
    }
    return seam;
}

ChartQuality measureCharts(const QuadMesh& mesh, const std::vector<int>& label, const QuadUv& uv) {
    struct TriangleMeasure {
        std::array<Vec3, 3> p;
        std::array<Vec2, 3> q;
        double area = 0, det = 0, signOfDet = 0;
        std::array<std::array<double, 2>, 3> jacobian;  // 3 x 2
        bool good = false;
    };
    ChartQuality quality;
    std::vector<TriangleMeasure> triangles;
    triangles.reserve(mesh.quad.size() * 2);
    for (int f = 0; f < static_cast<int>(mesh.quad.size()); ++f) {
        const std::array<int, 4>& q = mesh.quad[f];
        const std::array<Vec3, 4> corners = {mesh.pos[q[0]], mesh.pos[q[1]], mesh.pos[q[2]], mesh.pos[q[3]]};
        for (const std::array<int, 3>& split : quadSplit(corners)) {
            TriangleMeasure t;
            for (int k = 0; k < 3; ++k) {
                t.p[k] = corners[split[k]];
                t.q[k] = uv[f][split[k]];
            }
            const Vec3 e1 = t.p[1] - t.p[0];
            const Vec3 e2 = t.p[2] - t.p[0];
            const double a = t.q[1].x - t.q[0].x, b = t.q[2].x - t.q[0].x;
            const double c = t.q[1].y - t.q[0].y, d = t.q[2].y - t.q[0].y;
            const Vec3 normal = {e1.y * e2.z - e1.z * e2.y, e1.z * e2.x - e1.x * e2.z, e1.x * e2.y - e1.y * e2.x};
            t.area = 0.5 * norm(normal);
            t.det = a * d - c * b;
            t.signOfDet = t.det > 0 ? 1.0 : (t.det < 0 ? -1.0 : 0.0);
            t.good = std::fabs(t.det) > kDegenerateDet;
            if (t.good) {
                const double i00 = d / t.det, i01 = -b / t.det, i10 = -c / t.det, i11 = a / t.det;
                t.jacobian = {{{e1.x * i00 + e2.x * i10, e1.x * i01 + e2.x * i11},
                               {e1.y * i00 + e2.y * i10, e1.y * i01 + e2.y * i11},
                               {e1.z * i00 + e2.z * i10, e1.z * i01 + e2.z * i11}}};
            }
            triangles.push_back(t);
            quality.faceOfTriangle.push_back(f);
        }
    }

    const int count = static_cast<int>(triangles.size());
    quality.worst.assign(count, std::numeric_limits<double>::infinity());
    quality.chart = chartsInUse(label);
    std::map<int, std::vector<int>> membersOfChart;
    for (int t = 0; t < count; ++t) membersOfChart[label[quality.faceOfTriangle[t]]].push_back(t);

    for (int chart : quality.chart) {
        const std::vector<int>& members = membersOfChart[chart];
        double chartArea = 0, signedArea = 0;
        for (int t : members) {
            chartArea += triangles[t].area;
            signedArea += triangles[t].area * triangles[t].signOfDet;
        }
        const double majority = signedArea > 0 ? 1.0 : (signedArea < 0 ? -1.0 : 0.0);
        std::vector<int> usable;
        double foldedArea = 0;
        for (int t : members) {
            if (triangles[t].det * majority <= 0) foldedArea += triangles[t].area;
            else if (triangles[t].good) usable.push_back(t);
        }
        if (usable.empty()) {  // flattening failed: nothing of this chart is usable, so it is split
            quality.stretched.push_back(1.0);
            quality.unusable.push_back(1.0);
            continue;
        }
        // The mean metric, its root and inverse root (uv_atlas.h).
        double g00 = 0, g01 = 0, g11 = 0, usableArea = 0;
        for (int t : usable) {
            const std::array<std::array<double, 2>, 3>& j = triangles[t].jacobian;
            const double area = triangles[t].area;
            g00 += area * (j[0][0] * j[0][0] + j[1][0] * j[1][0] + j[2][0] * j[2][0]);
            g01 += area * (j[0][0] * j[0][1] + j[1][0] * j[1][1] + j[2][0] * j[2][1]);
            g11 += area * (j[0][1] * j[0][1] + j[1][1] * j[1][1] + j[2][1] * j[2][1]);
            usableArea += area;
        }
        g00 /= usableArea;
        g01 /= usableArea;
        g11 /= usableArea;
        const double rootDet = std::sqrt(g00 * g11 - g01 * g01);
        const double trace = std::sqrt(g00 + g11 + 2.0 * rootDet);
        const double r00 = (g00 + rootDet) / trace, r01 = g01 / trace, r11 = (g11 + rootDet) / trace;
        const double rDet = r00 * r11 - r01 * r01;
        const double i00 = r11 / rDet, i01 = -r01 / rDet, i11 = r00 / rDet;

        double outsideArea = 0;
        std::vector<std::array<Vec3, 3>> surface;
        std::vector<std::array<Vec2, 3>> flat;
        for (int t : usable) {
            const std::array<std::array<double, 2>, 3>& j = triangles[t].jacobian;
            double k0[3], k1[3];
            for (int row = 0; row < 3; ++row) {
                k0[row] = j[row][0] * i00 + j[row][1] * i01;
                k1[row] = j[row][0] * i01 + j[row][1] * i11;
            }
            const double a = k0[0] * k0[0] + k0[1] * k0[1] + k0[2] * k0[2];
            const double b = k0[0] * k1[0] + k0[1] * k1[1] + k0[2] * k1[2];
            const double c = k1[0] * k1[0] + k1[1] * k1[1] + k1[2] * k1[2];
            const double larger = 0.5 * (a + c + std::sqrt((a - c) * (a - c) + 4.0 * b * b));
            const double smaller = larger > 0 ? std::max(0.0, (a * c - b * b) / larger) : 0.0;
            const double s0 = std::sqrt(larger);
            const double s1 = std::sqrt(smaller);
            const double worst = std::max(s0 / kStretchHigh, kStretchLow / std::max(s1, 1e-9));
            quality.worst[t] = worst;
            if (worst > 1) outsideArea += triangles[t].area;
            surface.push_back(triangles[t].p);
            std::array<Vec2, 3> plane;
            for (int k = 0; k < 3; ++k)
                plane[k] = {triangles[t].q[k].x * r00 + triangles[t].q[k].y * r01, triangles[t].q[k].x * r01 + triangles[t].q[k].y * r11};
            flat.push_back(plane);
        }
        const std::vector<char> over = trianglesOverlapping(surface, flat);
        double overArea = 0;
        for (size_t i = 0; i < usable.size(); ++i)
            if (over[i]) overArea += triangles[usable[i]].area;
        quality.stretched.push_back(outsideArea / chartArea);
        quality.unusable.push_back((foldedArea + overArea) / chartArea);
    }

    double totalArea = 0, outside = 0;
    for (int t = 0; t < count; ++t) {
        totalArea += triangles[t].area;
        if (quality.worst[t] > 1) outside += triangles[t].area;
    }
    quality.totalOutside = outside / totalArea;
    return quality;
}

namespace {

// The chart split in two, the seams cut again (the slits kept), and both halves to be flattened.
ChartStep splitStep(const FaceGraph& graph, int chart, ChartRound& round) {
    ChartStep step;
    const int other = splitChart(graph, round.label, chart);
    round.seam = cutToDiscs(graph, round.label);
    for (size_t e = 0; e < round.seam.size(); ++e) round.seam[e] = round.seam[e] || round.slit[e];
    step.decision = ChartDecision::FlattenChanged;
    step.facesToFlatten.resize(round.label.size());
    for (size_t f = 0; f < round.label.size(); ++f) step.facesToFlatten[f] = round.label[f] == chart || round.label[f] == other;
    return step;
}

} // namespace

ChartStep decideNextRound(const FaceGraph& graph, const ChartQuality& quality, const ChartSettings& settings, ChartRound& round) {
    ChartStep step;
    const bool chartsLeft = static_cast<int>(quality.chart.size()) < settings.maxCharts;
    std::vector<int> broken;
    std::vector<double> brokenShare;
    for (size_t i = 0; i < quality.chart.size(); ++i)
        if (quality.unusable[i] > kMaxUnusable) {
            broken.push_back(quality.chart[i]);
            brokenShare.push_back(quality.unusable[i]);
        }
    if (!broken.empty() && chartsLeft) return splitStep(graph, chartWithLargest(broken, brokenShare), round);
    const double mostStretched = *std::max_element(quality.stretched.begin(), quality.stretched.end());
    if (quality.totalOutside <= kTargetOutside && mostStretched <= kMaxChartOutside) return step;

    // The worst chart, not the largest: one spike skews a chart's whole scale, so every stitch on it is sized wrong.
    const int chart = chartWithLargest(quality.chart, quality.stretched);
    // A slit relieves a spike, the tip of a cone. A chart mostly outside the bounds has no spike to relieve: its shape does not
    // lie flat (a tube closed at one end, say), so it is split, as a folded chart is.
    if (mostStretched > kSplitStretched && chartsLeft) return splitStep(graph, chart, round);
    if (round.slitCount >= settings.maxSlits) return step;
    double largestFinite = -std::numeric_limits<double>::infinity();
    for (size_t t = 0; t < quality.worst.size(); ++t)
        if (round.label[quality.faceOfTriangle[t]] == chart && std::isfinite(quality.worst[t])) largestFinite = std::max(largestFinite, quality.worst[t]);
    if (!std::isfinite(largestFinite)) throw std::runtime_error("a chart with no usable triangle cannot be slit");
    const double threshold = std::max(1.0, 0.5 * largestFinite);
    std::vector<int> tipFaces;
    for (size_t t = 0; t < quality.worst.size(); ++t)
        if (round.label[quality.faceOfTriangle[t]] == chart && quality.worst[t] >= threshold) tipFaces.push_back(quality.faceOfTriangle[t]);
    std::sort(tipFaces.begin(), tipFaces.end());
    tipFaces.erase(std::unique(tipFaces.begin(), tipFaces.end()), tipFaces.end());
    const std::vector<int> cut = slitFromSeam(graph, round.label, chart, tipFaces, round.seam);
    if (cut.empty()) return step;
    for (int e : cut) {
        round.slit[e] = 1;
        round.seam[e] = 1;
    }
    ++round.slitCount;
    step.decision = ChartDecision::FlattenChanged;
    step.facesToFlatten.resize(round.label.size());
    for (size_t f = 0; f < round.label.size(); ++f) step.facesToFlatten[f] = round.label[f] == chart;
    return step;
}

} // namespace knit
