// The rules the repairs apply and the analysis counts by (repair_rules.h).
#include "import/repair_rules.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_map>

#include "model/cancel.h"

namespace knit {

double triangleArea(const std::vector<Vec3>& p, const std::array<int, 3>& t) {
    return 0.5 * norm(cross(p[t[1]] - p[t[0]], p[t[2]] - p[t[0]]));
}

uint64_t edgeKey(int a, int b) { return (uint64_t(uint32_t(std::min(a, b))) << 32) | uint32_t(std::max(a, b)); }

std::vector<int> coincidentVertices(const std::vector<Vec3>& pos, double eps) {
    std::vector<int> remap(pos.size());
    if (!(eps > 0)) {
        std::iota(remap.begin(), remap.end(), 0);
        return remap;
    }
    auto cell = [&](double x) { return static_cast<long long>(std::floor(x / eps)); };
    auto key = [](long long x, long long y, long long z) {
        return (static_cast<uint64_t>(x) * 73856093u) ^ (static_cast<uint64_t>(y) * 19349663u) ^ (static_cast<uint64_t>(z) * 83492791u);
    };
    std::unordered_map<uint64_t, std::vector<int>> grid;
    grid.reserve(pos.size());
    for (size_t v = 0; v < pos.size(); ++v) {
        if ((v & 0xffff) == 0) cancelCheckpoint();
        const Vec3& p = pos[v];
        const long long cx = cell(p.x), cy = cell(p.y), cz = cell(p.z);
        int rep = -1;
        for (long long dx = -1; dx <= 1 && rep < 0; ++dx)
            for (long long dy = -1; dy <= 1 && rep < 0; ++dy)
                for (long long dz = -1; dz <= 1 && rep < 0; ++dz) {
                    auto it = grid.find(key(cx + dx, cy + dy, cz + dz));
                    if (it == grid.end()) continue;
                    for (int r : it->second)
                        if (norm(pos[r] - p) <= eps) {
                            rep = r;
                            break;
                        }
                }
        if (rep < 0) {
            rep = static_cast<int>(v);
            grid[key(cx, cy, cz)].push_back(rep);
        }
        remap[v] = rep;
    }
    return remap;
}

// Cancellable: a checkpoint per 64k triangles; the sort itself is about 0.5 s for 2 M triangles.
std::vector<std::pair<uint64_t, int>> sortedEdges(const std::vector<std::array<int, 3>>& tri) {
    std::vector<std::pair<uint64_t, int>> e;
    e.reserve(tri.size() * 3);
    for (size_t t = 0; t < tri.size(); ++t) {
        if ((t & 0xffff) == 0) cancelCheckpoint();
        for (int k = 0; k < 3; ++k) e.push_back({edgeKey(tri[t][k], tri[t][(k + 1) % 3]), static_cast<int>(3 * t + k)});
    }
    cancelCheckpoint();
    std::sort(e.begin(), e.end());
    cancelCheckpoint();
    return e;
}

BoundingBox boundingBoxOf(const std::vector<Vec3>& p) {
    BoundingBox box;
    for (const Vec3& q : p) {
        box.lo = {std::min(box.lo.x, q.x), std::min(box.lo.y, q.y), std::min(box.lo.z, q.z)};
        box.hi = {std::max(box.hi.x, q.x), std::max(box.hi.y, q.y), std::max(box.hi.z, q.z)};
    }
    box.diagonal = p.empty() ? 0 : norm(box.hi - box.lo);
    return box;
}

bool isZeroArea(const std::vector<Vec3>& pos, const std::array<int, 3>& t, double diagonal) {
    return t[0] == t[1] || t[1] == t[2] || t[0] == t[2] || triangleArea(pos, t) <= kZeroAreaRelative * diagonal * diagonal;
}

std::array<int, 3> sortedCorners(const std::array<int, 3>& t) {
    std::array<int, 3> s = t;
    std::sort(s.begin(), s.end());
    return s;
}

long long unreferencedVertices(const std::vector<std::array<int, 3>>& tri, size_t vertexCount) {
    std::vector<char> used(vertexCount, 0);
    for (const std::array<int, 3>& t : tri)
        for (int v : t) used[v] = 1;
    long long unused = 0;
    for (char u : used) unused += !u;
    return unused;
}

} // namespace knit
