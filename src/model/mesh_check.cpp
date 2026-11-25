// Is a mesh fit to knit on: welding, manifoldness, orientation, closedness, self-intersection (mesh_check.h).
#include "model/mesh_check.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <numeric>

#include "model/spatial_hash.h"

namespace knit {

WeldedMesh weldVertices(const TriMesh& in, double eps) {
    WeldedMesh welded;
    TriMesh& out = welded.mesh;
    std::vector<int>& map = welded.remap;
    map.assign(in.pos.size(), -1);
    SpatialHash3 hash(eps > 0 ? 4 * eps : 1.0);
    const Vec3 e{eps, eps, eps};
    for (size_t v = 0; v < in.pos.size(); ++v) {
        const Vec3& p = in.pos[v];
        int hit = -1;
        if (eps > 0)
            hash.query(p - e, p + e, [&](int id) {
                if (hit < 0 && norm(out.pos[id] - p) <= eps) hit = id;
            });
        if (hit < 0) {
            hit = static_cast<int>(out.pos.size());
            out.pos.push_back(p);
            hash.insert(p, p);
        }
        map[v] = hit;
    }
    for (size_t t = 0; t < in.tri.size(); ++t) {
        const std::array<int, 3> w = {map[in.tri[t][0]], map[in.tri[t][1]], map[in.tri[t][2]]};
        if (w[0] == w[1] || w[1] == w[2] || w[0] == w[2]) {
            ++welded.droppedTriangles;
            continue;
        }
        out.tri.push_back(w);
        if (!in.uv.empty()) out.uv.push_back(in.uv[t]);
    }
    out.computeNormals();
    return welded;
}

namespace {

// The gap between the triangles' projections on a unit axis: positive iff they are apart along it.
double separationAlong(const Vec3 a[3], const Vec3 b[3], const Vec3& axis) {
    double a0 = 1e300, a1 = -1e300, b0 = 1e300, b1 = -1e300;
    for (int k = 0; k < 3; ++k) {
        a0 = std::min(a0, dot(a[k], axis));
        a1 = std::max(a1, dot(a[k], axis));
        b0 = std::min(b0, dot(b[k], axis));
        b1 = std::max(b1, dot(b[k], axis));
    }
    return std::max(b0 - a1, a0 - b1);
}

} // namespace

Status triangleIntersection(const Vec3 a[3], const Vec3 b[3], double tol) {
    // Separating axis test over the 2 face normals, the 9 edge-edge cross products and
    // the 6 in-plane edge normals (coplanar case). The largest normalized gap over these
    // axes is positive iff the triangles are disjoint; it is a true separation distance
    // along that axis, so gap > tol decides PASS, every gap <= -tol decides FAIL.
    const Vec3 ea[3] = {a[1] - a[0], a[2] - a[1], a[0] - a[2]}, eb[3] = {b[1] - b[0], b[2] - b[1], b[0] - b[2]};
    const Vec3 na = cross(ea[0], a[2] - a[0]), nb = cross(eb[0], b[2] - b[0]);
    const double la = norm(na), lb = norm(nb);
    if (la > 0 && lb > 0 && norm(cross(na, nb)) <= 1e-12 * la * lb) {
        // Coplanar within tol: the normal axis has zero thickness, so overlap there is
        // never "by a margin"; decide in the plane (the 6 in-plane edge normals).
        const Vec3 u = na / la;
        const double d = std::fabs(dot(b[0] - a[0], u));
        if (d > tol) return Status::Pass;  // parallel planes apart
        double maxGap = -1e300;
        for (const Vec3* e : {ea, eb})
            for (int k = 0; k < 3; ++k) {
                const Vec3 ax = cross(u, e[k]);
                const double l = norm(ax);
                if (l < 1e-300) continue;
                maxGap = std::max(maxGap, separationAlong(a, b, ax / l));
            }
        return maxGap > tol ? Status::Pass : maxGap <= -tol ? Status::Fail : Status::Unresolved;
    }
    std::vector<Vec3> axes = {na, nb};
    for (const Vec3& x : ea)
        for (const Vec3& y : eb) axes.push_back(cross(x, y));
    for (const Vec3& x : ea) axes.push_back(cross(na, x));
    for (const Vec3& y : eb) axes.push_back(cross(nb, y));
    double maxGap = -1e300;
    for (const Vec3& ax : axes) {
        const double l = norm(ax);
        if (l < 1e-300) continue;
        maxGap = std::max(maxGap, separationAlong(a, b, ax / l));
        if (maxGap > tol) return Status::Pass;
    }
    return maxGap <= -tol ? Status::Fail : Status::Unresolved;
}

std::string MeshCheck::summary() const {
    std::string s;
    char b[400];
    for (size_t c = 0; c < components.size(); ++c) {
        const ComponentCheck& k = components[c];
        std::snprintf(b, sizeof b,
                      "component %zu: %.1f x %.1f x %.1f, area %.0f; %zu tris, %d verts, chi %d, %d boundary / %d non-manifold / %d inconsistent edges, "
                      "%s%s, %d uv charts, self-intersection %s (%d pairs, %d unresolved); inside tests %s%s%s\n",
                      c, k.hi.x - k.lo.x, k.hi.y - k.lo.y, k.hi.z - k.lo.z, k.area, k.triangles.size(), k.vertices, k.euler, k.boundaryEdges, k.nonManifoldEdges, k.inconsistentEdges,
                      k.closed ? "closed" : "open", k.closed ? (k.signedVolume > 0 ? ", volume > 0" : ", volume <= 0") : "",
                      k.uvCharts, toString(k.selfIntersection), k.intersectingPairs, k.unresolvedPairs,
                      k.insideEligible ? "eligible" : "NOT eligible", k.reason.empty() ? "" : ": ", k.reason.c_str());
        s += b;
    }
    return s;
}

namespace {

// Two corner uvs are the same within this, for chart continuity.
constexpr double kUvTolerance = 1e-7;

class UnionFind {
public:
    explicit UnionFind(size_t n) : parent_(n) { std::iota(parent_.begin(), parent_.end(), 0); }
    int find(int x) {
        while (parent_[x] != x) x = parent_[x] = parent_[parent_[x]];
        return x;
    }
    void join(int a, int b) { parent_[find(a)] = find(b); }

private:
    std::vector<int> parent_;
};

// Undirected edge (lower vertex, higher vertex) -> its uses: (triangle, the corner of the lower
// vertex k if the triangle runs lower to higher, else -k - 1).
using EdgeUses = std::map<std::pair<int, int>, std::vector<std::pair<int, int>>>;

// Components by shared vertices, largest first.
std::vector<ComponentCheck> splitIntoComponents(const TriMesh& m) {
    UnionFind vertices(m.pos.size());
    for (const auto& t : m.tri) {
        vertices.join(t[0], t[1]);
        vertices.join(t[1], t[2]);
    }
    std::map<int, std::vector<int>> byRoot;
    for (int t = 0; t < static_cast<int>(m.tri.size()); ++t) byRoot[vertices.find(m.tri[t][0])].push_back(t);
    std::vector<ComponentCheck> components;
    for (auto& [root, tris] : byRoot) {
        ComponentCheck c;
        c.triangles = std::move(tris);
        components.push_back(std::move(c));
    }
    std::sort(components.begin(), components.end(),
              [](const ComponentCheck& a, const ComponentCheck& b) { return a.triangles.size() > b.triangles.size(); });
    return components;
}

EdgeUses collectEdgeUses(const TriMesh& m, const ComponentCheck& c) {
    EdgeUses uses;
    for (int t : c.triangles)
        for (int k = 0; k < 3; ++k) {
            const int a = m.tri[t][k], b = m.tri[t][(k + 1) % 3];
            uses[{std::min(a, b), std::max(a, b)}].push_back({t, a < b ? k : -k - 1});
        }
    return uses;
}

// Bounds, area, the Euler characteristic, edge classes and the signed volume.
void measureComponent(const TriMesh& m, const EdgeUses& uses, ComponentCheck& c) {
    std::vector<int> verts;
    c.lo = {1e300, 1e300, 1e300};
    c.hi = {-1e300, -1e300, -1e300};
    for (int t : c.triangles) {
        const std::array<int, 3>& f = m.tri[t];
        for (int v : f) {
            verts.push_back(v);
            c.lo = componentMin(c.lo, m.pos[v]);
            c.hi = componentMax(c.hi, m.pos[v]);
        }
        c.area += 0.5 * norm(cross(m.pos[f[1]] - m.pos[f[0]], m.pos[f[2]] - m.pos[f[0]]));
    }
    std::sort(verts.begin(), verts.end());
    c.vertices = static_cast<int>(std::unique(verts.begin(), verts.end()) - verts.begin());
    c.edges = static_cast<int>(uses.size());
    for (const auto& [edge, u] : uses) {
        if (u.size() == 1) ++c.boundaryEdges;
        else if (u.size() > 2) ++c.nonManifoldEdges;
        else if ((u[0].second >= 0) == (u[1].second >= 0)) ++c.inconsistentEdges;  // same direction twice
    }
    c.euler = c.vertices - c.edges + static_cast<int>(c.triangles.size());
    c.closed = c.boundaryEdges == 0 && c.nonManifoldEdges == 0;
    c.consistent = c.inconsistentEdges == 0 && c.nonManifoldEdges == 0;
    if (c.closed)
        for (int t : c.triangles) {
            const std::array<int, 3>& f = m.tri[t];
            c.signedVolume += dot(m.pos[f[0]], cross(m.pos[f[1]], m.pos[f[2]])) / 6.0;
        }
}

// UV charts: triangles joined across edges whose two corner uvs agree on both sides.
int countUvCharts(const TriMesh& m, const ComponentCheck& c, const EdgeUses& uses) {
    std::map<int, int> local;
    for (size_t q = 0; q < c.triangles.size(); ++q) local[c.triangles[q]] = static_cast<int>(q);
    UnionFind charts(c.triangles.size());
    const auto uvAt = [&](int t, int vtx) {
        for (int k = 0; k < 3; ++k)
            if (m.tri[t][k] == vtx) return m.uv[t][k];
        return Vec2{};
    };
    const auto same = [&](const Vec2& p, const Vec2& q) {
        return std::fabs(p.x - q.x) <= kUvTolerance && std::fabs(p.y - q.y) <= kUvTolerance;
    };
    for (const auto& [edge, u] : uses) {
        if (u.size() != 2) continue;
        const int t0 = u[0].first, t1 = u[1].first;
        if (same(uvAt(t0, edge.first), uvAt(t1, edge.first)) && same(uvAt(t0, edge.second), uvAt(t1, edge.second)))
            charts.join(local[t0], local[t1]);
    }
    int count = 0;
    for (size_t q = 0; q < c.triangles.size(); ++q) count += charts.find(static_cast<int>(q)) == static_cast<int>(q);
    return count;
}

// Every pair of triangles that share no vertex, through a hash of their boxes.
void checkSelfIntersection(const TriMesh& m, double tol, ComponentCheck& c) {
    double meanEdge = 0;
    for (int t : c.triangles) meanEdge += norm(m.pos[m.tri[t][1]] - m.pos[m.tri[t][0]]);
    meanEdge /= std::max<size_t>(1, c.triangles.size());
    SpatialHash3 hash(2 * std::max(meanEdge, 1e-12));
    const Vec3 pad{tol, tol, tol};
    for (int t : c.triangles) {
        const std::array<int, 3>& f = m.tri[t];
        const Vec3 lo = componentMin(m.pos[f[0]], componentMin(m.pos[f[1]], m.pos[f[2]]));
        const Vec3 hi = componentMax(m.pos[f[0]], componentMax(m.pos[f[1]], m.pos[f[2]]));
        hash.insert(lo - pad, hi + pad);
    }
    for (size_t q = 0; q < c.triangles.size(); ++q) {
        const int t = c.triangles[q];
        const Vec3 A[3] = {m.pos[m.tri[t][0]], m.pos[m.tri[t][1]], m.pos[m.tri[t][2]]};
        hash.query(hash.lo(static_cast<int>(q)), hash.hi(static_cast<int>(q)), [&](int id) {
            if (id <= static_cast<int>(q)) return;
            const int s = c.triangles[id];
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    if (m.tri[t][i] == m.tri[s][j]) return;  // adjacent (shares a vertex): not tested
            const Vec3 B[3] = {m.pos[m.tri[s][0]], m.pos[m.tri[s][1]], m.pos[m.tri[s][2]]};
            const Status st = triangleIntersection(A, B, tol);
            c.intersectingPairs += st == Status::Fail;
            c.unresolvedPairs += st == Status::Unresolved;
        });
    }
    c.selfIntersection = c.intersectingPairs ? Status::Fail : c.unresolvedPairs ? Status::Unresolved : Status::Pass;
}

void decideEligibility(ComponentCheck& c) {
    std::string why;
    if (!c.closed) why += "open";
    if (!c.consistent) why += std::string(why.empty() ? "" : ", ") + "inconsistent orientation or non-manifold";
    if (c.closed && !(c.signedVolume > 0)) why += std::string(why.empty() ? "" : ", ") + "volume <= 0";
    if (c.selfIntersection != Status::Pass)
        why += std::string(why.empty() ? "" : ", ") + "self-intersection " + toString(c.selfIntersection);
    c.insideEligible = why.empty();
    c.reason = why;
}

} // namespace

MeshCheck checkMesh(const TriMesh& m, SelfIntersectionTest test) {
    MeshCheck out;
    out.components = splitIntoComponents(m);
    Vec3 lo, hi;
    m.bounds(lo, hi);
    const double tol = 1e-9 * std::max(1e-12, norm(hi - lo));
    for (ComponentCheck& c : out.components) {
        const EdgeUses uses = collectEdgeUses(m, c);
        measureComponent(m, uses, c);
        if (!m.uv.empty()) c.uvCharts = countUvCharts(m, c, uses);
        if (test == SelfIntersectionTest::Run) checkSelfIntersection(m, tol, c);
        decideEligibility(c);
    }
    return out;
}

} // namespace knit
