// TriMesh and the closest-point primitives the mesh queries and checks are built on (mesh.h).
#include "model/mesh.h"
#include "model/hash.h"

#include <algorithm>
#include <cstdint>
#include <string>

namespace knit {

void TriMesh::computeNormals() {
    vnormal.assign(pos.size(), Vec3{});
    for (const auto& t : tri) {
        Vec3 n = cross(pos[t[1]] - pos[t[0]], pos[t[2]] - pos[t[0]]);  // |n| = 2 area
        for (int k = 0; k < 3; ++k) vnormal[t[k]] += n;
    }
    for (Vec3& n : vnormal) n = normalized(n);
}

void TriMesh::transform(double s, const Vec3& origin) {
    for (Vec3& p : pos) p = (p - origin) * s;
}

double TriMesh::area() const {
    double a = 0;
    for (const auto& t : tri) a += 0.5 * norm(cross(pos[t[1]] - pos[t[0]], pos[t[2]] - pos[t[0]]));
    return a;
}

Vec3 TriMesh::centroid() const {
    Vec3 c{};
    double a = 0;
    for (const auto& t : tri) {
        double w = 0.5 * norm(cross(pos[t[1]] - pos[t[0]], pos[t[2]] - pos[t[0]]));
        c += (pos[t[0]] + pos[t[1]] + pos[t[2]]) * (w / 3.0);
        a += w;
    }
    return a > 0 ? c / a : c;
}

Vec3 TriMesh::faceNormal(int t) const {
    const std::array<int, 3>& f = tri[t];
    return normalized(cross(pos[f[1]] - pos[f[0]], pos[f[2]] - pos[f[0]]));
}

void TriMesh::bounds(Vec3& lo, Vec3& hi) const {
    lo = hi = pos.empty() ? Vec3{} : pos[0];
    for (const Vec3& p : pos) {
        lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
        hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
    }
}

ClosestSegmentPoints closestSegmentPoints(const Vec3& p1, const Vec3& q1, const Vec3& p2, const Vec3& q2) {
    const double EPS = 1e-300;
    double s = 0, t = 0;
    Vec3 d1 = q1 - p1, d2 = q2 - p2, r = p1 - p2;
    double a = dot(d1, d1), e = dot(d2, d2), f = dot(d2, r);
    if (a <= EPS && e <= EPS) { s = t = 0; return {norm(p1 - p2), s, t}; }
    if (a <= EPS) {
        s = 0; t = std::clamp(f / e, 0.0, 1.0);
    } else {
        double c = dot(d1, r);
        if (e <= EPS) {
            t = 0; s = std::clamp(-c / a, 0.0, 1.0);
        } else {
            double b = dot(d1, d2), denom = a * e - b * b;
            s = denom > 0 ? std::clamp((b * f - c * e) / denom, 0.0, 1.0) : 0.0;
            t = (b * s + f) / e;
            if (t < 0) { t = 0; s = std::clamp(-c / a, 0.0, 1.0); }
            else if (t > 1) { t = 1; s = std::clamp((b - c) / a, 0.0, 1.0); }
        }
    }
    Vec3 c1 = p1 + d1 * s, c2 = p2 + d2 * t;
    return {norm(c1 - c2), s, t};
}

TrianglePoint closestOnTriangle(const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c) {
    TrianglePoint r;
    double* bary = r.bary;
    Vec3 ab = b - a, ac = c - a, ap = p - a;
    double d1 = dot(ab, ap), d2 = dot(ac, ap);
    if (d1 <= 0 && d2 <= 0) { bary[0] = 1; bary[1] = bary[2] = 0; r.p = a; return r; }
    Vec3 bp = p - b;
    double d3 = dot(ab, bp), d4 = dot(ac, bp);
    if (d3 >= 0 && d4 <= d3) { bary[1] = 1; bary[0] = bary[2] = 0; r.p = b; return r; }
    double vc = d1 * d4 - d3 * d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) {
        double v = d1 / (d1 - d3);
        bary[0] = 1 - v; bary[1] = v; bary[2] = 0;
        r.p = a + ab * v;
        return r;
    }
    Vec3 cp = p - c;
    double d5 = dot(ab, cp), d6 = dot(ac, cp);
    if (d6 >= 0 && d5 <= d6) { bary[2] = 1; bary[0] = bary[1] = 0; r.p = c; return r; }
    double vb = d5 * d2 - d1 * d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) {
        double w = d2 / (d2 - d6);
        bary[0] = 1 - w; bary[1] = 0; bary[2] = w;
        r.p = a + ac * w;
        return r;
    }
    double va = d3 * d6 - d5 * d4;
    if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
        double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        bary[0] = 0; bary[1] = 1 - w; bary[2] = w;
        r.p = b + (c - b) * w;
        return r;
    }
    double denom = 1.0 / (va + vb + vc);
    double v = vb * denom, w = vc * denom;
    bary[0] = 1 - v - w; bary[1] = v; bary[2] = w;
    r.p = a + ab * v + ac * w;
    return r;
}

uint64_t meshDigest(const TriMesh& m) {
    std::string s;
    s.append(reinterpret_cast<const char*>(m.pos.data()), m.pos.size() * sizeof(Vec3));
    for (const auto& t : m.tri) s.append(reinterpret_cast<const char*>(t.data()), sizeof t);
    return fnv1a64(s);
}

} // namespace knit
