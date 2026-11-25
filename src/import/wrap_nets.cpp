// Surface extraction and smoothing of the built-in import (surface nets, then Taubin 1995). The width of each number is fixed:
// a crossing's position along its edge is a float division, and the sum of the crossings is in doubles.
#include "import/wrap_stages.h"

#include <algorithm>
#include <cstdint>
#include <utility>

#include "model/cancel.h"
#include "model/parallel.h"

namespace knit {

namespace {

constexpr double kTaubinShrink = 0.5;
constexpr double kTaubinInflate = -0.53;

// The eight corners of a cell, x as the slowest index, and its twelve edges in a fixed order.
constexpr int kCorner[8][3] = {{0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1}, {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1}};
constexpr int kCellEdge[12][2] = {{0, 1}, {0, 2}, {0, 4}, {1, 3}, {1, 5}, {2, 3}, {2, 6}, {3, 7}, {4, 5}, {4, 6}, {5, 7}, {6, 7}};

// The mean of the crossings on a cell's edges, in grid coordinates.
Vec3 meanOfEdgeCrossings(const FieldGrid& field, int x, int y, int z) {
    Vec3 sum;
    int crossings = 0;
    for (const auto& edge : kCellEdge) {
        const int* a = kCorner[edge[0]];
        const int* b = kCorner[edge[1]];
        const float f0 = field.values[voxelIndex(field.size, x + a[0], y + a[1], z + a[2])];
        const float f1 = field.values[voxelIndex(field.size, x + b[0], y + b[1], z + b[2])];
        if ((f0 > 0) == (f1 > 0)) continue;
        const double t = f0 / (f0 - f1);
        sum.x += (x + a[0]) + t * (b[0] - a[0]);
        sum.y += (y + a[1]) + t * (b[1] - a[1]);
        sum.z += (z + a[2]) + t * (b[2] - a[2]);
        ++crossings;
    }
    return {sum.x / crossings, sum.y / crossings, sum.z / crossings};
}

// Quads around the grid edges of one axis that the surface crosses, in C order over those edges.
void appendQuadsOfAxis(const FieldGrid& field, const std::vector<int>& vertexOfCell, int axis, QuadMesh& mesh) {
    const std::array<int, 3>& size = field.size;
    const int uAxis = (axis + 1) % 3;
    const int vAxis = (axis + 2) % 3;
    std::array<int, 3> low = {1, 1, 1};
    std::array<int, 3> high = size;
    low[axis] = 0;
    high[axis] = size[axis] - 1;
    std::array<int, 3> p;
    for (p[0] = low[0]; p[0] < high[0]; ++p[0])
        for (p[1] = low[1]; p[1] < high[1]; ++p[1])
            for (p[2] = low[2]; p[2] < high[2]; ++p[2]) {
                std::array<int, 3> next = p;
                ++next[axis];
                const bool here = field.values[voxelIndex(size, p[0], p[1], p[2])] > 0;
                const bool there = field.values[voxelIndex(size, next[0], next[1], next[2])] > 0;
                if (here == there) continue;
                std::array<int, 3> pu = p;  // p - u
                std::array<int, 3> pv = p;  // p - v
                std::array<int, 3> puv = p;  // p - u - v
                --pu[uAxis];
                --pv[vAxis];
                --puv[uAxis];
                --puv[vAxis];
                std::array<int, 4> quad = {vertexOfCell[voxelIndex(size, puv[0], puv[1], puv[2])],
                                           vertexOfCell[voxelIndex(size, pv[0], pv[1], pv[2])],
                                           vertexOfCell[voxelIndex(size, p[0], p[1], p[2])],
                                           vertexOfCell[voxelIndex(size, pu[0], pu[1], pu[2])]};
                if (!here) std::reverse(quad.begin(), quad.end());
                mesh.quad.push_back(quad);
            }
}

// For each vertex its neighbours along quad edges, in ascending order, as compressed lists.
struct VertexNeighbours {
    std::vector<int> start;  // one entry more than there are vertices
    std::vector<int> neighbour;
};

VertexNeighbours neighboursOf(const QuadMesh& mesh) {
    std::vector<std::pair<int, int>> edges;
    edges.reserve(mesh.quad.size() * 8);
    for (const std::array<int, 4>& q : mesh.quad)
        for (int k = 0; k < 4; ++k) {
            edges.push_back({q[k], q[(k + 1) % 4]});
            edges.push_back({q[(k + 1) % 4], q[k]});
        }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    VertexNeighbours lists;
    lists.start.assign(mesh.pos.size() + 1, 0);
    lists.neighbour.reserve(edges.size());
    for (const std::pair<int, int>& e : edges) {
        ++lists.start[e.first + 1];
        lists.neighbour.push_back(e.second);
    }
    for (size_t v = 0; v < mesh.pos.size(); ++v) lists.start[v + 1] += lists.start[v];
    return lists;
}

// One step of the filter: every vertex moves by `factor` towards the mean of its neighbours' old positions.
void moveTowardsNeighbours(const VertexNeighbours& lists, double factor, std::vector<Vec3>& pos, std::vector<Vec3>& next) {
    const int count = static_cast<int>(pos.size());
    parallelFor(count, 4096,
                [&](int64_t begin, int64_t end, int64_t) {
                    for (int v = static_cast<int>(begin); v < end; ++v) {
                        Vec3 sum;
                        for (int k = lists.start[v]; k < lists.start[v + 1]; ++k) sum += pos[lists.neighbour[k]];
                        const double degree = lists.start[v + 1] - lists.start[v];
                        next[v] = {pos[v].x + factor * (sum.x / degree - pos[v].x), pos[v].y + factor * (sum.y / degree - pos[v].y),
                                   pos[v].z + factor * (sum.z / degree - pos[v].z)};
                    }
                },
                ExecPolicy{.parallel = true});
    pos.swap(next);
}

} // namespace

QuadMesh extractSurfaceNets(const FieldGrid& field, const GridPlacement& grid) {
    const std::array<int, 3>& size = field.size;
    std::vector<int> vertexOfCell(field.values.size(), -1);
    QuadMesh mesh;
    for (int x = 0; x < size[0] - 1; ++x) {
        cancelCheckpoint();
        for (int y = 0; y < size[1] - 1; ++y)
            for (int z = 0; z < size[2] - 1; ++z) {
                int inside = 0;
                for (const int* c : kCorner) inside += field.values[voxelIndex(size, x + c[0], y + c[1], z + c[2])] > 0;
                if (inside == 0 || inside == 8) continue;
                vertexOfCell[voxelIndex(size, x, y, z)] = static_cast<int>(mesh.pos.size());
                const Vec3 p = meanOfEdgeCrossings(field, x, y, z);
                mesh.pos.push_back({p.x * grid.voxelMm + grid.origin.x, p.y * grid.voxelMm + grid.origin.y,
                                    p.z * grid.voxelMm + grid.origin.z});
            }
    }
    for (int axis = 0; axis < 3; ++axis) appendQuadsOfAxis(field, vertexOfCell, axis, mesh);
    return mesh;
}

void smoothTaubin(QuadMesh& mesh, int iterations) {
    const VertexNeighbours lists = neighboursOf(mesh);
    std::vector<Vec3> next(mesh.pos.size());
    for (int i = 0; i < iterations; ++i) {
        cancelCheckpoint();
        moveTowardsNeighbours(lists, kTaubinShrink, mesh.pos, next);
        moveTowardsNeighbours(lists, kTaubinInflate, mesh.pos, next);
    }
}

} // namespace knit
