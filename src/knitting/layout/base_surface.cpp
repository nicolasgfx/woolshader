// The surfaces stitches are placed on: cone, plane, and MeshSurface with its ownership-preserving projection (base_surface.h).
#include "knitting/layout/base_surface.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include "model/mesh_query.h"

namespace knit {

void BaseSurface::atBatch(const std::vector<Vec2>& ab, std::vector<SurfaceSample>& out) const {
    out.resize(ab.size());
    for (size_t i = 0; i < ab.size(); ++i) out[i] = at(ab[i].x, ab[i].y);
}

void BaseSurface::projectBatch(const std::vector<Vec3>& q, std::vector<SurfaceSample>& out) const {
    out.resize(q.size());
    for (size_t i = 0; i < q.size(); ++i) out[i] = project(q[i]);
}

SurfaceSample ConeSurface::at(double a, double b) const {
    const double t = 2 * kPi * b;
    Vec3 radial{std::cos(t), std::sin(t), 0};
    // Outward normal of the cone: radial tilted against the slope dR/dz.
    const double slope = (R1_ - R0_) / H_;
    Vec3 n = normalized(radial - Vec3{0, 0, slope});
    return {radial * radiusAt(a) + Vec3{0, 0, a * H_}, n};
}

SurfaceSample ConeSurface::project(const Vec3& q) const {
    // Closest point on the generating line in the (rho, z) half-plane of q.
    Vec3 r{q.x, q.y, 0};
    const double rho = norm(r);
    Vec3 radial = rho > 0 ? r / rho : Vec3{1, 0, 0};
    const Vec3 A{R0_, 0, 0}, d{R1_ - R0_, 0, H_};  // (rho, -, z) line, a in R
    const Vec3 Q{rho, 0, q.z};
    double a = dot(Q - A, d) / dot(d, d);
    const double slope = (R1_ - R0_) / H_;
    Vec3 n = normalized(radial - Vec3{0, 0, slope});
    return {radial * radiusAt(a) + Vec3{0, 0, a * H_}, n};
}

bool ProjState::operator==(const ProjState& o) const {
    return tri == o.tri && b[0] == o.b[0] && b[1] == o.b[1] && b[2] == o.b[2] && foot.x == o.foot.x && foot.y == o.foot.y &&
           foot.z == o.foot.z && n.x == o.n.x && n.y == o.n.y && n.z == o.n.z && query.x == o.query.x && query.y == o.query.y &&
           query.z == o.query.z;
}

const char* toString(ProjResult r) {
    switch (r) {
    case ProjResult::Owned: return "owned";
    case ProjResult::Budget: return "budget";
    case ProjResult::NormalTurn: return "normal-turn";
    case ProjResult::OpenBoundary: return "open-boundary";
    case ProjResult::DepthBand: return "depth-band";
    }
    return "?";
}

int SurfaceBinding::unresolved() const {
    int n = 0;
    for (const auto& y : result)
        for (ProjResult r : y) n += r != ProjResult::Owned;
    return n;
}

MeshSurface::MeshSurface(const TriMesh& m) : m_(m), q_(m) {
    const int T = static_cast<int>(m.tri.size());
    nbr_.assign(T, {-1, -1, -1});
    vertexTris_.assign(m.pos.size(), {});
    std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> edges;
    for (int t = 0; t < T; ++t)
        for (int k = 0; k < 3; ++k) {
            const int a = m.tri[t][k], b = m.tri[t][(k + 1) % 3];
            edges[{std::min(a, b), std::max(a, b)}].push_back({t, k});
            vertexTris_[a].push_back(t);
        }
    for (const auto& [e, u] : edges)
        if (u.size() == 2) {
            nbr_[u[0].first][u[0].second] = u[1].first;
            nbr_[u[1].first][u[1].second] = u[0].first;
        }
}

SurfaceSample MeshSurface::at(double, double) const {
    throw std::logic_error("MeshSurface::at: a mesh has no global chart (use LayoutNode wale hints)");
}

Vec3 MeshSurface::normalAt(int tri, const double b[3]) const {
    Vec3 n = m_.vnormal[m_.tri[tri][0]] * b[0] + m_.vnormal[m_.tri[tri][1]] * b[1] + m_.vnormal[m_.tri[tri][2]] * b[2];
    if (norm(n) < 1e-12) n = m_.faceNormal(tri);
    return normalized(n);
}

void MeshSurface::continueBatch(const std::vector<ProjState>& prev, const std::vector<Vec3>& q, const ContinuationParams& cp,
                                std::vector<ProjState>& next, std::vector<ProjResult>& res) const {
    next.resize(prev.size());
    res.resize(prev.size());
    for (size_t i = 0; i < prev.size(); ++i) res[i] = continueFrom(prev[i], q[i], next[i], cp);
}

ProjState MeshSurface::stateAt(int tri, const double b[3], const Vec3& query) const {
    ProjState s;
    s.tri = tri;
    for (int k = 0; k < 3; ++k) s.b[k] = b[k];
    s.foot = m_.pos[m_.tri[tri][0]] * b[0] + m_.pos[m_.tri[tri][1]] * b[1] + m_.pos[m_.tri[tri][2]] * b[2];
    s.n = normalAt(tri, b);
    s.query = query;
    return s;
}

SurfaceSample MeshSurface::project(const Vec3& q) const {
    const MeshQuery::Hit h = q_.closest(q);
    const std::array<int, 3>& t = m_.tri[h.tri];
    const TrianglePoint onTri = closestOnTriangle(q, m_.pos[t[0]], m_.pos[t[1]], m_.pos[t[2]]);
    return {onTri.p, normalAt(h.tri, onTri.bary)};
}

ProjResult MeshSurface::continueFrom(const ProjState& prev, const Vec3& q, ProjState& next, const ContinuationParams& cp) const {
    next = prev;
    if (prev.tri < 0) return ProjResult::Budget;
    const double budget = norm(q - prev.query) + cp.lambdaW;
    auto closest = [&](int t, double b[3]) {
        const std::array<int, 3>& v = m_.tri[t];
        const TrianglePoint onTri = closestOnTriangle(q, m_.pos[v[0]], m_.pos[v[1]], m_.pos[v[2]]);
        std::copy(onTri.bary, onTri.bary + 3, b);
        return onTri.p;
    };
    int t = prev.tri;
    double b[3];
    Vec3 x = closest(t, b);
    double d = norm(q - x), walked = norm(x - prev.foot);
    if (dot(normalAt(t, b), prev.n) < cp.minNormalCos) return ProjResult::NormalTurn;
    if (walked > budget) return ProjResult::Budget;
    const double zero = 1e-12;
    int steps = 0;
    for (; steps < cp.maxSteps; ++steps) {
        // Neighbours that touch the current closest point: across an edge (one zero
        // weight) or around a vertex (two zero weights).
        int nz = 0, corner = -1;
        for (int k = 0; k < 3; ++k)
            if (b[k] <= zero) ++nz;
            else corner = corner < 0 ? k : corner;
        std::vector<int> cand;
        if (nz == 1) {
            for (int k = 0; k < 3; ++k)
                if (b[k] <= zero) {
                    const int e = (k + 1) % 3;  // edge opposite corner k: corners k + 1, k + 2
                    if (nbr_[t][e] >= 0) cand.push_back(nbr_[t][e]);
                }
        } else if (nz == 2) {
            for (int u : vertexTris_[m_.tri[t][corner]])
                if (u != t) cand.push_back(u);
        }
        int best = -1;
        double bestD = d - 1e-12 * (1 + d), bb[3] = {0, 0, 0};
        Vec3 bx;
        for (int c : cand) {
            double bc[3];
            const Vec3 xc = closest(c, bc);
            const double dc = norm(q - xc);
            if (dc < bestD) {
                best = c;
                bestD = dc;
                bx = xc;
                for (int k = 0; k < 3; ++k) bb[k] = bc[k];
            }
        }
        if (best < 0) break;
        walked += norm(bx - x);
        t = best;
        x = bx;
        d = bestD;
        for (int k = 0; k < 3; ++k) b[k] = bb[k];
        if (dot(normalAt(t, b), prev.n) < cp.minNormalCos) return ProjResult::NormalTurn;
        if (walked > budget) return ProjResult::Budget;
    }
    if (steps >= cp.maxSteps) return ProjResult::Budget;
    // A foot on an open-boundary edge (or vertex) of the mesh: q lies beyond the edge.
    for (int k = 0; k < 3; ++k)
        if (b[k] <= zero && nbr_[t][(k + 1) % 3] < 0) return ProjResult::OpenBoundary;
    const Vec3 n = normalAt(t, b);
    if (std::fabs(dot(q - x, n)) > cp.band) return ProjResult::DepthBand;
    next = stateAt(t, b, q);
    return ProjResult::Owned;
}

} // namespace knit
