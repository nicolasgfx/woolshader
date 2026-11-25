// Splitting a piece along closed marking curves into regions (surface_trim.h).
#include "stitchmesh/surface_trim.h"

#include <algorithm>
#include <cstdio>
#include <map>
#include <set>
#include <utility>

namespace knit {

namespace {

using EdgeKey = std::pair<int, int>;  // undirected: smaller vertex first

EdgeKey edgeKey(int a, int b) { return {std::min(a, b), std::max(a, b)}; }

struct Crossing {
    int from = -1, to = -1;
    double s = 0;
};

// A crossing's position along its edge, kept away from the corners so that no split triangle is a sliver
// (a shift below 5 % of one edge).
constexpr double kCornerMargin = 0.05;

// The curve's edge crossings in order, without the back-and-forth pairs a sharp turn inside one triangle leaves
// (two consecutive crossings of the same edge cancel).
std::vector<Crossing> curveCrossings(const MarkingTrace& c) {
    std::vector<Crossing> x;
    for (const TracePoint& p : c.points)
        if (p.edgeFrom >= 0) x.push_back(Crossing{p.edgeFrom, p.edgeTo, p.s});
    bool removed = true;
    while (removed && x.size() >= 2) {
        removed = false;
        for (size_t i = 0; i < x.size(); ++i) {
            const size_t j = (i + 1) % x.size();
            if (edgeKey(x[i].from, x[i].to) == edgeKey(x[j].from, x[j].to)) {
                const size_t hi = std::max(i, j), lo = std::min(i, j);
                x.erase(x.begin() + hi);
                x.erase(x.begin() + lo);
                removed = true;
                break;
            }
        }
    }
    return x;
}

struct CrossingAt {
    int curve = -1, index = -1, vertex = -1;
};

} // namespace

bool splitAlongCurves(const TriMesh& piece, const std::vector<MarkingTrace>& curves, TrimmedSurface& out, std::string* err) {
    out = TrimmedSurface{};
    out.mesh.pos = piece.pos;
    std::map<EdgeKey, CrossingAt> onEdge;
    std::vector<int> crossingsOfCurve;
    char b[200];
    for (int c = 0; c < static_cast<int>(curves.size()); ++c) {
        const MarkingTrace& tr = curves[c];
        if (tr.edgePathUsed || !tr.failedSegments.empty()) {
            std::snprintf(b, sizeof b, "curve %d: a segment could not be drawn across the surface; add a control point in between", c + 1);
            if (err) *err = b;
            return false;
        }
        const std::vector<Crossing> x = curveCrossings(tr);
        if (x.size() < 3) {
            std::snprintf(b, sizeof b, "curve %d is too small to cut the surface (it crosses fewer than 3 mesh edges)", c + 1);
            if (err) *err = b;
            return false;
        }
        std::vector<int> loop;
        for (int i = 0; i < static_cast<int>(x.size()); ++i) {
            const EdgeKey k = edgeKey(x[i].from, x[i].to);
            if (onEdge.count(k)) {
                std::snprintf(b, sizeof b, "curve %d crosses a mesh edge that is already cut (a curve crosses itself or another curve)", c + 1);
                if (err) *err = b;
                return false;
            }
            const double s = std::min(1 - kCornerMargin, std::max(kCornerMargin, x[i].s));
            const int v = static_cast<int>(out.mesh.pos.size());
            out.mesh.pos.push_back(piece.pos[x[i].from] + (piece.pos[x[i].to] - piece.pos[x[i].from]) * s);
            onEdge[k] = CrossingAt{c, i, v};
            loop.push_back(v);
        }
        out.loops.push_back(loop);
        crossingsOfCurve.push_back(static_cast<int>(x.size()));
    }
    // split the crossed triangles; the chord of each split is a barrier between regions
    std::map<EdgeKey, int> barrierCurve;
    for (int t = 0; t < static_cast<int>(piece.tri.size()); ++t) {
        const std::array<int, 3>& c = piece.tri[t];
        int cut[3];
        int count = 0;
        for (int k = 0; k < 3; ++k) {
            const std::map<EdgeKey, CrossingAt>::const_iterator it = onEdge.find(edgeKey(c[k], c[(k + 1) % 3]));
            cut[k] = it == onEdge.end() ? -1 : k;
            count += it == onEdge.end() ? 0 : 1;
        }
        if (count == 0) {
            out.mesh.tri.push_back(c);
            out.parentTri.push_back(t);
            continue;
        }
        if (count != 2) {
            if (err) *err = "a triangle is crossed by a curve on one or three of its edges (curves too close to each other or to themselves)";
            return false;
        }
        // edges i and i + 1 are cut: their shared corner is c[i + 1]
        int i = 0;
        while (!(cut[i] >= 0 && cut[(i + 1) % 3] >= 0)) ++i;
        const int a = c[(i + 1) % 3], next = c[(i + 2) % 3], prev = c[i];
        const CrossingAt P = onEdge.at(edgeKey(prev, a)), Q = onEdge.at(edgeKey(a, next));
        const int n = crossingsOfCurve[P.curve];
        if (P.curve != Q.curve || !((P.index + 1) % n == Q.index || (Q.index + 1) % n == P.index)) {
            if (err) *err = "a curve passes the same triangle twice; move its control points apart";
            return false;
        }
        // (prev, a, next) counter-clockwise: corner triangle (P, a, Q) and the quad (prev, P, Q, next)
        out.mesh.tri.push_back({P.vertex, a, Q.vertex});
        out.mesh.tri.push_back({prev, P.vertex, Q.vertex});
        out.mesh.tri.push_back({prev, Q.vertex, next});
        out.parentTri.insert(out.parentTri.end(), {t, t, t});
        barrierCurve[edgeKey(P.vertex, Q.vertex)] = P.curve;
    }
    // regions: flood fill across edges that are not a curve's chord
    const int nt = static_cast<int>(out.mesh.tri.size());
    std::map<EdgeKey, std::vector<int>> trianglesOfEdge;
    for (int t = 0; t < nt; ++t)
        for (int k = 0; k < 3; ++k) trianglesOfEdge[edgeKey(out.mesh.tri[t][k], out.mesh.tri[t][(k + 1) % 3])].push_back(t);
    out.region.assign(nt, -1);
    for (int seed = 0; seed < nt; ++seed) {
        if (out.region[seed] >= 0) continue;
        const int r = out.regions++;
        std::vector<int> stack{seed};
        out.region[seed] = r;
        while (!stack.empty()) {
            const int t = stack.back();
            stack.pop_back();
            for (int k = 0; k < 3; ++k) {
                const EdgeKey e = edgeKey(out.mesh.tri[t][k], out.mesh.tri[t][(k + 1) % 3]);
                if (barrierCurve.count(e)) continue;
                for (int o : trianglesOfEdge[e])
                    if (out.region[o] < 0) {
                        out.region[o] = r;
                        stack.push_back(o);
                    }
            }
        }
    }
    out.touches.assign(out.regions, std::vector<char>(curves.size(), 0));
    for (const std::pair<const EdgeKey, int>& e : barrierCurve)
        for (int t : trianglesOfEdge[e.first]) out.touches[out.region[t]][e.second] = 1;
    out.mesh.computeNormals();
    return true;
}

RegionMesh regionMesh(const TrimmedSurface& s, int region) {
    RegionMesh r;
    std::vector<int> local(s.mesh.pos.size(), -1);
    for (int t = 0; t < static_cast<int>(s.mesh.tri.size()); ++t) {
        if (s.region[t] != region) continue;
        std::array<int, 3> tri{};
        for (int k = 0; k < 3; ++k) {
            const int v = s.mesh.tri[t][k];
            if (local[v] < 0) {
                local[v] = static_cast<int>(r.mesh.pos.size());
                r.mesh.pos.push_back(s.mesh.pos[v]);
            }
            tri[k] = local[v];
        }
        r.mesh.tri.push_back(tri);
        r.parentTri.push_back(s.parentTri[t]);
    }
    for (const std::vector<int>& loop : s.loops) {
        std::vector<int> l;
        for (int v : loop)
            if (local[v] >= 0) l.push_back(local[v]);
        r.loops.push_back(l);
    }
    r.local = local;
    r.mesh.computeNormals();
    return r;
}

} // namespace knit
