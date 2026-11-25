// MeshQuery: the closest point on a mesh and the inside test (mesh_query.h).
#include "model/mesh_query.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <unordered_map>

namespace knit {

namespace {
uint64_t edgeKey(int a, int b) {
    if (a > b) std::swap(a, b);
    return (uint64_t(uint32_t(a)) << 32) | uint32_t(b);
}

double boxDistanceSquared(const Vec3& q, const Vec3& lo, const Vec3& hi) {
    const double dx = std::max({lo.x - q.x, 0.0, q.x - hi.x}), dy = std::max({lo.y - q.y, 0.0, q.y - hi.y}),
                 dz = std::max({lo.z - q.z, 0.0, q.z - hi.z});
    return dx * dx + dy * dy + dz * dz;
}
} // namespace

MeshQuery::MeshQuery(const TriMesh& m, double cell) : m_(m) {
    double edges = 0;
    for (const auto& t : m.tri)
        for (int k = 0; k < 3; ++k) edges += norm(m.pos[t[(k + 1) % 3]] - m.pos[t[k]]);
    const double mean = m.tri.empty() ? 1.0 : edges / (3.0 * m.tri.size());
    hash_ = SpatialHash3(cell > 0 ? cell : 2 * mean);
    m.bounds(lo_, hi_);
    faceN_.resize(m.tri.size());
    vertN_.assign(m.pos.size(), Vec3{0, 0, 0});
    std::unordered_map<uint64_t, Vec3> edgeSum;
    for (size_t f = 0; f < m.tri.size(); ++f) {
        const std::array<int, 3>& t = m.tri[f];
        const Vec3 &a = m.pos[t[0]], &b = m.pos[t[1]], &c = m.pos[t[2]];
        faceN_[f] = normalized(cross(b - a, c - a));
        Vec3 lo{std::min({a.x, b.x, c.x}), std::min({a.y, b.y, c.y}), std::min({a.z, b.z, c.z})};
        Vec3 hi{std::max({a.x, b.x, c.x}), std::max({a.y, b.y, c.y}), std::max({a.z, b.z, c.z})};
        hash_.insert(lo, hi);
        for (int k = 0; k < 3; ++k) {
            // Angle-weighted vertex pseudo-normal; edge pseudo-normal = sum of the faces.
            const Vec3 e1 = normalized(m.pos[t[(k + 1) % 3]] - m.pos[t[k]]);
            const Vec3 e2 = normalized(m.pos[t[(k + 2) % 3]] - m.pos[t[k]]);
            vertN_[t[k]] += faceN_[f] * std::acos(std::clamp(dot(e1, e2), -1.0, 1.0));
            edgeSum[edgeKey(t[k], t[(k + 1) % 3])] += faceN_[f];
        }
    }
    hash_.finalize();
    cornerEdgeN_.resize(3 * m.tri.size());
    for (size_t f = 0; f < m.tri.size(); ++f)
        for (int k = 0; k < 3; ++k) cornerEdgeN_[3 * f + k] = edgeSum.at(edgeKey(m.tri[f][(k + 1) % 3], m.tri[f][(k + 2) % 3]));
}

MeshQuery::Hit MeshQuery::closest(const Vec3& q, double maxDist) const {
    Hit best;
    double bb[3] = {1, 0, 0};
    // Expanding search: a triangle within distance r of q has a box overlapping the
    // query box of half size r, so the best hit is final once best.dist <= r.
    const double toBox = std::sqrt(boxDistanceSquared(q, lo_, hi_));
    const double rMax = std::min(maxDist, toBox + norm(hi_ - lo_));
    // The rounds' boxes are nested (same centre, growing half size): a triangle was
    // evaluated in an earlier round iff its box overlaps the previous round's box, and then
    // dist >= best.dist already, so `dist < best.dist` cannot hold again.
    // A triangle whose box is farther from q than the current best (by more than the
    // rounding margin of both distances) cannot pass `dist < best.dist` either.
    const double margin = 1e-9 * (1 + std::max({std::fabs(q.x), std::fabs(q.y), std::fabs(q.z)}) + norm(hi_ - lo_));
    bool havePrev = false;
    Vec3 pLo, pHi;
    for (double r = std::max(hash_.cell(), toBox); ; r *= 2) {
        const double rr = std::min(r, rMax);
        const Vec3 qLo = q - Vec3{rr, rr, rr}, qHi = q + Vec3{rr, rr, rr};
        const auto visit = [&](int f) {
            const Vec3 &a = hash_.lo(f), &b = hash_.hi(f);
            if (havePrev && !(a.x > pHi.x || b.x < pLo.x || a.y > pHi.y || b.y < pLo.y || a.z > pHi.z || b.z < pLo.z)) return;
            if (best.tri >= 0) {
                const double lim = best.dist + margin;
                if (boxDistanceSquared(q, a, b) > lim * lim) return;
            }
            const std::array<int, 3>& t = m_.tri[f];
            const TrianglePoint onTri = closestOnTriangle(q, m_.pos[t[0]], m_.pos[t[1]], m_.pos[t[2]]);
            const double dist = norm(q - onTri.p);
            if (dist < best.dist) {
                best = {f, onTri.p, dist, false};
                std::copy(onTri.bary, onTri.bary + 3, bb);
            }
        };
        hash_.query(qLo, qHi, std::ref(visit));  // by reference: too large for std::function's inline buffer
        if (best.dist <= rr || rr >= rMax) break;
        havePrev = true;
        pLo = qLo;
        pHi = qHi;
    }
    if (best.tri < 0 || best.dist > maxDist) return Hit{};
    const std::array<int, 3>& t = m_.tri[best.tri];
    int zeros = 0;
    for (double w : bb) zeros += w == 0.0;
    Vec3 n = faceN_[best.tri];
    if (zeros == 2) {  // vertex
        for (int k = 0; k < 3; ++k)
            if (bb[k] == 1.0) n = vertN_[t[k]];
    } else if (zeros == 1) {  // edge
        for (int k = 0; k < 3; ++k)
            if (bb[k] == 0.0) n = cornerEdgeN_[3 * best.tri + k];
    }
    best.inside = best.dist > 0 && dot(q - best.p, n) < 0;
    return best;
}

bool MeshQuery::inside(const Vec3& q) const {
    if (q.x < lo_.x || q.y < lo_.y || q.z < lo_.z || q.x > hi_.x || q.y > hi_.y || q.z > hi_.z) return false;
    return closest(q).inside;
}

bool MeshQuery::segmentCrosses(const Vec3& a, const Vec3& b, double skip) const {
    const Vec3 d = b - a;
    const double length = norm(d);
    if (!(length > skip)) return false;
    const Vec3 lo{std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z)};
    const Vec3 hi{std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)};
    const double tMin = skip / length;
    bool hit = false;
    hash_.query(lo, hi, [&](int f) {
        if (hit) return;
        // Moeller-Trumbore: a + t d = v0 + u e1 + v e2, with u, v >= 0, u + v <= 1, t in (tMin, 1].
        const std::array<int, 3>& t = m_.tri[static_cast<size_t>(f)];
        const Vec3& v0 = m_.pos[t[0]];
        const Vec3 e1 = m_.pos[t[1]] - v0, e2 = m_.pos[t[2]] - v0;
        const Vec3 pv = cross(d, e2);
        const double det = dot(e1, pv);
        if (std::fabs(det) < 1e-18) return;  // parallel to the triangle
        const double inv = 1 / det;
        const Vec3 tv = a - v0;
        const double u = dot(tv, pv) * inv;
        if (u < 0 || u > 1) return;
        const Vec3 qv = cross(tv, e1);
        const double v = dot(d, qv) * inv;
        if (v < 0 || u + v > 1) return;
        const double tt = dot(e2, qv) * inv;
        hit = tt > tMin && tt <= 1;
    });
    return hit;
}

} // namespace knit
