#include "stitchmesh/row_rings.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <utility>

namespace knit {

namespace {

struct Crossing {
    int a = -1, b = -1;   // a < b
    double f = 0;
    Vec3 p;
};

// The level's crossing on the edge (u, v), always computed from the lower index so both triangles of the edge agree
// bit for bit.
Crossing crossingOn(const TriMesh& m, const std::vector<double>& t, int u, int v, double level) {
    Crossing c;
    c.a = std::min(u, v);
    c.b = std::max(u, v);
    const double ta = t[(size_t) c.a], tb = t[(size_t) c.b];
    c.f = (level - ta) / (tb - ta);
    c.p = m.pos[c.a] + (m.pos[c.b] - m.pos[c.a]) * c.f;
    return c;
}

struct Segment {
    Crossing from, to;
    int tri = -1;
};

std::pair<int, int> key(const Crossing& c) { return {c.a, c.b}; }

} // namespace

std::vector<Contour> dropContoursShorterThan(RingLevel& level, double minLength) {
    std::vector<Contour> dropped;
    if (level.contours.size() < 2) return dropped;
    std::vector<Contour> kept;
    for (const Contour& c : level.contours) {
        if (c.closed && c.length < minLength) dropped.push_back(c);
        else kept.push_back(c);
    }
    if (kept.empty() || dropped.empty()) return {};   // all short (nothing to knit over) or none short: left as found
    std::stable_sort(kept.begin(), kept.end(), [](const Contour& a, const Contour& b) { return a.length > b.length; });
    level.contours = kept;
    return dropped;
}

RingLevel extractLevel(const TriMesh& m, const std::vector<double>& t, double level) {
    RingLevel out;
    out.level = level;
    std::vector<Segment> segs;
    for (int f = 0; f < (int) m.tri.size(); ++f) {
        const std::array<int, 3>& tr = m.tri[f];
        bool above[3];
        int nAbove = 0;
        for (int k = 0; k < 3; ++k) {
            above[k] = t[(size_t) tr[k]] >= level;
            nAbove += above[k] ? 1 : 0;
        }
        if (nAbove == 0 || nAbove == 3) continue;
        // Orient: larger values on the left of the segment, seen from outside (the triangle's winding). Walking the
        // corners in winding order, the segment starts on the edge that goes from above to below and ends on the edge that
        // goes from below to above. This is decided from the corners' classes, not from positions: when the level passes
        // through a vertex (a regular grid whose rings fall on the rows), crossings coincide with that vertex or with each
        // other and a geometric left/right test has no sign, which broke the chaining into false split contours.
        Segment s;
        s.tri = f;
        for (int k = 0; k < 3; ++k) {
            const int k1 = (k + 1) % 3;
            if (above[k] && !above[k1]) s.from = crossingOn(m, t, tr[k], tr[k1], level);
            if (!above[k] && above[k1]) s.to = crossingOn(m, t, tr[k], tr[k1], level);
        }
        segs.push_back(s);
    }
    out.crossedTriangles = (int) segs.size();
    std::map<std::pair<int, int>, int> startsAt;
    std::map<std::pair<int, int>, int> endsAt;
    for (int i = 0; i < (int) segs.size(); ++i) {
        startsAt[key(segs[i].from)] = i;
        endsAt[key(segs[i].to)] = i;
    }
    std::vector<char> used(segs.size(), 0);
    // Open contours first: they start where no segment ends (an edge on a boundary or a hole).
    for (int pass = 0; pass < 2; ++pass)
        for (int i = 0; i < (int) segs.size(); ++i) {
            if (used[(size_t) i]) continue;
            if (pass == 0 && endsAt.count(key(segs[i].from))) continue;
            Contour c;
            int s = i;
            while (s >= 0 && !used[(size_t) s]) {
                used[(size_t) s] = 1;
                RingPoint rp;
                rp.p = segs[(size_t) s].from.p;
                rp.a = segs[(size_t) s].from.a;
                rp.b = segs[(size_t) s].from.b;
                rp.f = segs[(size_t) s].from.f;
                rp.tri = segs[(size_t) s].tri;
                c.points.push_back(rp);
                const auto next = startsAt.find(key(segs[(size_t) s].to));
                if (next == startsAt.end()) {   // the contour ends on a boundary edge
                    RingPoint last;
                    last.p = segs[(size_t) s].to.p;
                    last.a = segs[(size_t) s].to.a;
                    last.b = segs[(size_t) s].to.b;
                    last.f = segs[(size_t) s].to.f;
                    c.points.push_back(last);
                    s = -1;
                } else {
                    s = next->second;
                    if (used[(size_t) s]) c.closed = (s == i);
                }
            }
            for (size_t k = 1; k < c.points.size(); ++k) c.length += norm(c.points[k].p - c.points[k - 1].p);
            if (c.closed && c.points.size() > 1) c.length += norm(c.points.front().p - c.points.back().p);
            out.contours.push_back(std::move(c));
        }
    return out;
}

double cutCrossing(const Contour& c, const Vec3& origin, const Vec3& normal, const Vec3& side, int* crossings) {
    int count = 0;
    double where = -1;
    const size_t n = c.points.size();
    const size_t pairs = c.closed ? n : (n > 0 ? n - 1 : 0);
    for (size_t i = 0; i < pairs; ++i) {
        const Vec3& p = c.points[i].p;
        const Vec3& q = c.points[(i + 1) % n].p;
        const double gp = dot(p - origin, normal), gq = dot(q - origin, normal);
        if ((gp < 0) == (gq < 0)) continue;
        const double f = gp / (gp - gq);
        const Vec3 x = p + (q - p) * f;
        if (dot(x - origin, side) <= 0) continue;
        ++count;
        where = (double) i + f;
    }
    if (crossings) *crossings = count;
    return count == 1 ? where : -1;
}

namespace {

// 2D coordinates in the plane of triangle f (corner 0 at the origin).
struct TrianglePlane {
    Vec3 origin, e1, e2;
};

TrianglePlane trianglePlane(const TriMesh& m, int f) {
    TrianglePlane p;
    p.origin = m.pos[m.tri[f][0]];
    const Vec3 n = m.faceNormal(f);
    p.e1 = normalized(m.pos[m.tri[f][1]] - p.origin);
    p.e2 = cross(n, p.e1);
    return p;
}

Vec2 inPlane(const TrianglePlane& p, const Vec3& x) { return Vec2{dot(x - p.origin, p.e1), dot(x - p.origin, p.e2)}; }

double cross2(const Vec2& a, const Vec2& b) { return a.x * b.y - a.y * b.x; }

} // namespace

double pathCrossing(const TriMesh& m, const Contour& c, const SurfacePathGuide& g, int* crossings, std::vector<double>* pathParameters) {
    // the path's segments by the triangle they lie in
    std::map<int, std::vector<int>> segmentsOf;
    for (int i = 0; i + 1 < static_cast<int>(g.points.size()); ++i) segmentsOf[g.segmentTri[i]].push_back(i);
    int count = 0;
    double where = -1;
    const size_t n = c.points.size();
    const size_t pairs = c.closed ? n : (n > 0 ? n - 1 : 0);
    for (size_t i = 0; i < pairs; ++i) {
        const int f = c.points[i].tri;
        if (f < 0) continue;
        const int parent = g.parentTri ? (*g.parentTri)[f] : f;
        const std::map<int, std::vector<int>>::const_iterator it = segmentsOf.find(parent);
        if (it == segmentsOf.end()) continue;
        const TrianglePlane pl = trianglePlane(m, f);
        const Vec2 p = inPlane(pl, c.points[i].p), q = inPlane(pl, c.points[(i + 1) % n].p);
        const Vec2 r = q - p;
        for (int k : it->second) {
            const Vec2 u = inPlane(pl, g.points[k]), v = inPlane(pl, g.points[k + 1]);
            const Vec2 s = v - u;
            const double den = cross2(r, s);
            if (std::fabs(den) < 1e-18) continue;  // parallel
            const double a = cross2(u - p, s) / den, b = cross2(u - p, r) / den;
            if (a < 0 || a >= 1 || b < 0 || b >= 1) continue;  // half-open: a shared end point counts once
            ++count;
            where = static_cast<double>(i) + a;
            if (pathParameters) pathParameters->push_back(static_cast<double>(k) + b);
        }
    }
    if (crossings) *crossings = count;
    return count == 1 ? where : -1;
}

double pathNearest(const Contour& c, const SurfacePathGuide& g) {
    double best = 1e300, where = -1;
    const size_t n = c.points.size();
    const size_t pairs = c.closed ? n : (n > 0 ? n - 1 : 0);
    for (size_t i = 0; i < pairs; ++i)
        for (size_t k = 0; k + 1 < g.points.size(); ++k) {
            const ClosestSegmentPoints cs = closestSegmentPoints(c.points[i].p, c.points[(i + 1) % n].p, g.points[k], g.points[k + 1]);
            if (cs.distance < best) {
                best = cs.distance;
                where = static_cast<double>(i) + std::min(cs.s, 0.999999);
            }
        }
    return where;
}

} // namespace knit
