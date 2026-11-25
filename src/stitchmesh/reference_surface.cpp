#include "stitchmesh/reference_surface.h"

#include <cmath>
#include <set>
#include <utility>

namespace knit {

bool prepareReferenceSurface(const TriMesh& m, ReferenceSurface& out, std::string* why) {
    out = ReferenceSurface{};
    if (m.tri.empty()) {
        if (why) *why = "the reference surface has no triangles";
        return false;
    }
    out.mesh = &m;
    out.check = checkMesh(m);
    out.closed = true;
    for (size_t k = 0; k < out.check.components.size(); ++k) {
        const ComponentCheck& c = out.check.components[k];
        if (c.nonManifoldEdges > 0 || !c.consistent) {
            if (why) *why = "component " + std::to_string(k) + ": " + std::to_string(c.nonManifoldEdges) + " non-manifold edges, " +
                            std::to_string(c.inconsistentEdges) + " inconsistently oriented edges";
            return false;
        }
        if (!c.closed) out.closed = false;
    }
    out.faceNormal.resize(m.tri.size());
    double edgeSum = 0;
    std::set<std::pair<int, int>> edges;
    for (size_t t = 0; t < m.tri.size(); ++t) {
        const Vec3 n = cross(m.pos[m.tri[t][1]] - m.pos[m.tri[t][0]], m.pos[m.tri[t][2]] - m.pos[m.tri[t][0]]);
        const double twiceArea = norm(n);
        if (!(twiceArea > 0)) ++out.degenerateTriangles;
        out.faceNormal[t] = twiceArea > 0 ? n / twiceArea : Vec3{};
        out.area += 0.5 * twiceArea;
        for (int k = 0; k < 3; ++k) {
            const int a = m.tri[t][k], b = m.tri[t][(k + 1) % 3];
            if (edges.insert({std::min(a, b), std::max(a, b)}).second) edgeSum += norm(m.pos[a] - m.pos[b]);
        }
    }
    out.meanEdge = edges.empty() ? 0 : edgeSum / edges.size();
    out.digest = meshDigest(m);
    out.query = std::make_unique<MeshQuery>(m);
    return true;
}

std::vector<VertexBinding> bindToSurface(const ReferenceSurface& s, const std::vector<Vec3>& points) {
    std::vector<VertexBinding> out(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        const MeshQuery::Hit hit = s.query->closest(points[i]);
        if (hit.tri < 0) continue;
        const std::array<int, 3>& t = s.mesh->tri[hit.tri];
        const TrianglePoint foot = closestOnTriangle(points[i], s.mesh->pos[t[0]], s.mesh->pos[t[1]], s.mesh->pos[t[2]]);
        out[i].tri = hit.tri;
        for (int k = 0; k < 3; ++k) out[i].bary[k] = foot.bary[k];
        out[i].distance = norm(points[i] - foot.p);
    }
    return out;
}

} // namespace knit
