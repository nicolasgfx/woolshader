// Connected parts by vertices, their areas, and the removal of the small ones.
#include "import/wrap_stages.h"

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <unordered_map>
#include <vector>

namespace knit {

namespace {

constexpr double kFragmentShare = 0.01;  // a part above this share of the area is no fragment

// The connected parts, each vertex's part number in order of first appearance.
std::vector<int> partOfVertices(const QuadMesh& mesh, int& partCount) {
    std::vector<int> parent(mesh.pos.size());
    std::iota(parent.begin(), parent.end(), 0);
    const auto rootOf = [&](int v) {
        while (parent[v] != v) {
            parent[v] = parent[parent[v]];
            v = parent[v];
        }
        return v;
    };
    for (const std::array<int, 4>& q : mesh.quad)
        for (int k = 1; k < 4; ++k) {
            const int a = rootOf(q[0]);
            const int b = rootOf(q[k]);
            parent[std::max(a, b)] = std::min(a, b);
        }
    std::vector<int> part(mesh.pos.size(), -1);
    partCount = 0;
    for (size_t v = 0; v < mesh.pos.size(); ++v) {
        const int root = rootOf(static_cast<int>(v));
        if (part[root] < 0) part[root] = partCount++;
        part[v] = part[root];
    }
    return part;
}

// A quad's area: half the length of the cross product of its diagonals, (c0 - c2) x (c1 - c3).
double quadArea(const QuadMesh& mesh, const std::array<int, 4>& q) {
    return 0.5 * norm(cross(mesh.pos[q[0]] - mesh.pos[q[2]], mesh.pos[q[1]] - mesh.pos[q[3]]));
}

std::vector<double> areaOfParts(const QuadMesh& mesh, const std::vector<int>& part, int partCount) {
    std::vector<double> area(partCount, 0.0);
    for (const std::array<int, 4>& q : mesh.quad) area[part[q[0]]] += quadArea(mesh, q);
    return area;
}

} // namespace

double secondPartShare(const QuadMesh& mesh) {
    int partCount = 0;
    const std::vector<int> part = partOfVertices(mesh, partCount);
    if (partCount < 2) return 0.0;
    std::vector<double> area = areaOfParts(mesh, part, partCount);
    std::sort(area.begin(), area.end(), std::greater<double>());
    return area[1] / std::accumulate(area.begin(), area.end(), 0.0);
}

bool dropFragments(QuadMesh& mesh, int& genus, std::string* err) {
    int partCount = 0;
    const std::vector<int> part = partOfVertices(mesh, partCount);
    const std::vector<double> area = areaOfParts(mesh, part, partCount);
    const double total = std::accumulate(area.begin(), area.end(), 0.0);
    const int largest = static_cast<int>(std::max_element(area.begin(), area.end()) - area.begin());
    for (int p = 0; p < partCount; ++p)
        if (p != largest && area[p] > kFragmentShare * total) {
            if (err) *err = "a separate part of " + std::to_string(int(area[p])) + " of " + std::to_string(int(total)) +
                            " mm^2 remains: it needs a wider closing";
            return false;
        }

    if (partCount > 1) {
        std::vector<int> newVertex(mesh.pos.size(), -1);
        std::vector<Vec3> pos;
        for (size_t v = 0; v < mesh.pos.size(); ++v)
            if (part[v] == largest) {
                newVertex[v] = static_cast<int>(pos.size());
                pos.push_back(mesh.pos[v]);
            }
        std::vector<std::array<int, 4>> quads;
        for (const std::array<int, 4>& q : mesh.quad)
            if (part[q[0]] == largest) quads.push_back({newVertex[q[0]], newVertex[q[1]], newVertex[q[2]], newVertex[q[3]]});
        std::vector<std::array<int, 2>> edges;
        for (const std::array<int, 2>& e : mesh.edge)
            if (part[e[0]] == largest) edges.push_back({newVertex[e[0]], newVertex[e[1]]});
        mesh.pos.swap(pos);
        mesh.quad.swap(quads);
        mesh.edge.swap(edges);
    }

    std::unordered_map<uint64_t, int> facesOnEdge;
    facesOnEdge.reserve(mesh.quad.size() * 2);
    for (const std::array<int, 4>& q : mesh.quad)
        for (int k = 0; k < 4; ++k) {
            const uint32_t a = q[k];
            const uint32_t b = q[(k + 1) % 4];
            ++facesOnEdge[(uint64_t(std::min(a, b)) << 32) | std::max(a, b)];
        }
    long long openEdges = 0, branchingEdges = 0;
    for (const std::pair<const uint64_t, int>& entry : facesOnEdge) {
        openEdges += entry.second == 1;
        branchingEdges += entry.second > 2;
    }
    if (openEdges || branchingEdges) {
        if (err) *err = "not a closed manifold: " + std::to_string(openEdges) + " open edges, " + std::to_string(branchingEdges) + " non-manifold edges";
        return false;
    }
    const long long euler = static_cast<long long>(mesh.pos.size()) - static_cast<long long>(mesh.edge.size()) + static_cast<long long>(mesh.quad.size());
    genus = static_cast<int>((2 - euler) / 2);
    return true;
}

} // namespace knit
