// Plane-section walks and edge paths between surface points (surface_path.h).
#include "model/surface_path.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <queue>
#include <utility>

namespace knit {

TriangleNeighbours triangleNeighbours(const TriMesh& m) {
    TriangleNeighbours nb;
    const int nt = static_cast<int>(m.tri.size());
    nb.across.assign(nt, {-1, -1, -1});
    // undirected edge -> (triangle, corner) of its first use; a second use links both; a third makes it non-manifold
    std::map<std::pair<int, int>, std::pair<int, int>> first;
    std::map<std::pair<int, int>, int> uses;
    for (int t = 0; t < nt; ++t)
        for (int k = 0; k < 3; ++k) {
            const int a = m.tri[t][k], b = m.tri[t][(k + 1) % 3];
            const std::pair<int, int> key{std::min(a, b), std::max(a, b)};
            const int n = ++uses[key];
            if (n == 1) {
                first[key] = {t, k};
            } else if (n == 2) {
                const std::pair<int, int> f = first[key];
                nb.across[t][k] = f.first;
                nb.across[f.first][f.second] = t;
            } else {
                const std::pair<int, int> f = first[key];
                nb.across[f.first][f.second] = -1;  // non-manifold: no walk across it
                nb.across[t][k] = -1;
            }
        }
    // a non-manifold edge may have linked its first two uses before the third appeared: unlink every such edge
    for (int t = 0; t < nt; ++t)
        for (int k = 0; k < 3; ++k) {
            const int a = m.tri[t][k], b = m.tri[t][(k + 1) % 3];
            if (uses[{std::min(a, b), std::max(a, b)}] > 2) nb.across[t][k] = -1;
        }
    return nb;
}

Vec3 surfacePosition(const TriMesh& m, const SurfacePoint& p) {
    const std::array<int, 3>& t = m.tri[p.tri];
    return m.pos[t[0]] * p.bary[0] + m.pos[t[1]] * p.bary[1] + m.pos[t[2]] * p.bary[2];
}

SurfacePoint closestSurfacePoint(const TriMesh& m, const Vec3& p) {
    SurfacePoint best;
    double bestDistance = 1e300;
    for (int t = 0; t < static_cast<int>(m.tri.size()); ++t) {
        const TrianglePoint tp = closestOnTriangle(p, m.pos[m.tri[t][0]], m.pos[m.tri[t][1]], m.pos[m.tri[t][2]]);
        const double d = norm(tp.p - p);
        if (d >= bestDistance) continue;
        bestDistance = d;
        best.tri = t;
        for (int k = 0; k < 3; ++k) best.bary[k] = tp.bary[k];
    }
    return best;
}

namespace {

double pathLength(const std::vector<TracePoint>& p) {
    double L = 0;
    for (size_t i = 1; i < p.size(); ++i) L += norm(p[i].pos - p[i - 1].pos);
    return L;
}

// The corner k of triangle t whose edge (k, k + 1) joins a and b in either direction (-1: none).
int edgeIndex(const TriMesh& m, int t, int a, int b) {
    for (int k = 0; k < 3; ++k) {
        const int u = m.tri[t][k], v = m.tri[t][(k + 1) % 3];
        if ((u == a && v == b) || (u == b && v == a)) return k;
    }
    return -1;
}

// Whether p lies on triangle t (within tol).
bool triangleHolds(const TriMesh& m, int t, const Vec3& p, double tol) {
    const TrianglePoint tp = closestOnTriangle(p, m.pos[m.tri[t][0]], m.pos[m.tri[t][1]], m.pos[m.tri[t][2]]);
    return norm(tp.p - p) <= tol;
}

// Whether the plane's section crosses triangle t (some edge has corners on both sides).
bool sectionCrosses(const TriMesh& m, int t, const Vec3& pa, const Vec3& q) {
    for (int k = 0; k < 3; ++k) {
        const double du = dot(m.pos[m.tri[t][k]] - pa, q), dv = dot(m.pos[m.tri[t][(k + 1) % 3]] - pa, q);
        if ((du >= 0) != (dv >= 0)) return true;
    }
    return false;
}

// The triangles that hold surface point a: a's own triangle for an inner point; also the triangle across the edge for
// a point on an edge; the whole fan round the corner for a point on a vertex (walked through the neighbours both ways
// round, stopping at an open edge).
std::vector<int> trianglesAtPoint(const TriMesh& m, const TriangleNeighbours& nb, const SurfacePoint& a) {
    constexpr double kOnBoundary = 1e-9;   // a barycentric coordinate this close to 0 (or 1) puts a on an edge (or a corner)
    std::vector<int> out = {a.tri};
    int corner = -1, offEdge = -1;
    for (int k = 0; k < 3; ++k) {
        if (a.bary[k] > 1 - kOnBoundary) corner = k;
        if (a.bary[k] < kOnBoundary) offEdge = k;
    }
    if (corner >= 0) {
        const int v = m.tri[a.tri][corner];
        // round the fan: first across the edge leaving v, then (from a's triangle again) across the edge arriving at v
        for (int way = 0; way < 2; ++way) {
            int t = a.tri;
            for (int step = 0; step < 64; ++step) {
                int cv = 0;
                while (m.tri[t][cv] != v) ++cv;
                const int next = nb.across[t][way == 0 ? cv : (cv + 2) % 3];
                if (next < 0 || std::find(out.begin(), out.end(), next) != out.end()) break;
                out.push_back(next);
                t = next;
            }
        }
    } else if (offEdge >= 0) {
        const int across = nb.across[a.tri][(offEdge + 1) % 3];   // edge (offEdge + 1, offEdge + 2): the one a lies on
        if (across >= 0) out.push_back(across);
    }
    return out;
}

// The triangle to start a walk in: among the triangles that hold a, the one where the plane's section leaves a furthest
// towards b (forward) or away from it (backward); crossings at a itself are not counted. A point on a vertex or an edge
// lies in several triangles and the section leaves it through more than one of them: starting in an arbitrary one sent
// the forward and the backward walk the same way round the section, the long way for one of them. -1 when the section
// leaves a in no triangle (the caller keeps a's own triangle).
int leavingTriangle(const TriMesh& m, const TriangleNeighbours& nb, const SurfacePoint& a, const Vec3& pa, const Vec3& q, const Vec3& dir,
                    bool forward, double tol) {
    int best = -1;
    double bestAhead = 0;
    for (int t : trianglesAtPoint(m, nb, a))
        for (int k = 0; k < 3; ++k) {
            const Vec3 pu = m.pos[m.tri[t][k]], pv = m.pos[m.tri[t][(k + 1) % 3]];
            const double du = dot(pu - pa, q), dv = dot(pv - pa, q);
            if ((du >= 0) == (dv >= 0)) continue;
            const Vec3 x = pu + (pv - pu) * (du / (du - dv));
            if (norm(x - pa) <= tol) continue;
            const double ahead = dot(x - pa, dir);
            if (forward ? ahead > bestAhead : ahead < bestAhead) {
                best = t;
                bestAhead = ahead;
            }
        }
    return best;
}

// Walks the section of the surface with the plane {x : dot(x - pa, q) = 0} from a, leaving a's triangle through the
// crossing that lies forward (towards b) or backward. Signed distances >= 0 count as positive, so the section never
// passes exactly through a vertex. Ends in the first triangle that holds b. When a lies on an edge and the plane
// leaves a's triangle on one side (its corners in the plane), the walk starts from the triangle across that edge.
bool walkPlane(const TriMesh& m, const TriangleNeighbours& nb, const SurfacePoint& a, const SurfacePoint& b, const Vec3& q, bool forward,
               std::vector<TracePoint>& out) {
    const Vec3 pa = surfacePosition(m, a), pb = surfacePosition(m, b);
    const Vec3 dir = pb - pa;
    const double tol = 1e-9 * (1 + norm(dir));
    int t = leavingTriangle(m, nb, a, pa, q, dir, forward, 1e-7 * (1 + norm(dir))), entry = -1;
    if (t < 0) t = a.tri;
    if (!sectionCrosses(m, t, pa, q))
        for (int k = 0; k < 3; ++k) {
            const int across = nb.across[a.tri][k];
            if (across >= 0 && triangleHolds(m, across, pa, tol) && sectionCrosses(m, across, pa, q)) t = across;
        }
    out.clear();
    TracePoint start;
    start.pos = pa;
    start.tri = t;
    out.push_back(start);
    const int maxSteps = static_cast<int>(m.tri.size()) + 1;
    for (int step = 0; step < maxSteps; ++step) {
        if (t == b.tri || triangleHolds(m, t, pb, tol)) {
            TracePoint end;
            end.pos = pb;
            end.tri = t;
            out.push_back(end);
            return true;
        }
        int chosen = -1;
        double chosenAhead = 0, chosenS = 0;
        for (int k = 0; k < 3; ++k) {
            if (k == entry) continue;
            const int u = m.tri[t][k], v = m.tri[t][(k + 1) % 3];
            const double du = dot(m.pos[u] - pa, q), dv = dot(m.pos[v] - pa, q);
            if ((du >= 0) == (dv >= 0)) continue;
            const double s = du / (du - dv);
            const Vec3 x = m.pos[u] + (m.pos[v] - m.pos[u]) * s;
            const double ahead = dot(x - pa, dir);
            const bool better = chosen < 0 || (forward ? ahead > chosenAhead : ahead < chosenAhead);
            if (entry >= 0 || better) {
                chosen = k;
                chosenAhead = ahead;
                chosenS = s;
            }
        }
        if (chosen < 0) return false;
        const int next = nb.across[t][chosen];
        if (next < 0) return false;  // an open or non-manifold edge
        TracePoint p;
        p.edgeFrom = m.tri[t][chosen];
        p.edgeTo = m.tri[t][(chosen + 1) % 3];
        p.s = chosenS;
        p.pos = m.pos[p.edgeFrom] + (m.pos[p.edgeTo] - m.pos[p.edgeFrom]) * chosenS;
        p.tri = next;
        out.push_back(p);
        entry = edgeIndex(m, next, p.edgeFrom, p.edgeTo);
        t = next;
    }
    return false;
}

int nearestCorner(const TriMesh& m, const SurfacePoint& p) {
    int best = 0;
    for (int k = 1; k < 3; ++k)
        if (p.bary[k] > p.bary[best]) best = k;
    return m.tri[p.tri][best];
}

// The shortest path along mesh edges from a's nearest corner to b's nearest corner, with a and b at its ends.
bool edgePath(const TriMesh& m, const SurfacePoint& a, const SurfacePoint& b, std::vector<TracePoint>& out) {
    const int nv = static_cast<int>(m.pos.size());
    std::vector<std::vector<int>> adj(nv);
    std::vector<int> triOfEdgeStart(nv, -1);
    for (int t = 0; t < static_cast<int>(m.tri.size()); ++t)
        for (int k = 0; k < 3; ++k) {
            const int u = m.tri[t][k], v = m.tri[t][(k + 1) % 3];
            adj[u].push_back(v);
            triOfEdgeStart[u] = t;
        }
    const int s = nearestCorner(m, a), g = nearestCorner(m, b);
    std::vector<double> dist(nv, 1e300);
    std::vector<int> prev(nv, -1);
    using Item = std::pair<double, int>;
    std::priority_queue<Item, std::vector<Item>, std::greater<Item>> open;
    dist[s] = 0;
    open.push({0, s});
    while (!open.empty()) {
        const Item it = open.top();
        open.pop();
        if (it.first > dist[it.second]) continue;
        if (it.second == g) break;
        for (int v : adj[it.second]) {
            const double d = it.first + norm(m.pos[v] - m.pos[it.second]);
            if (d < dist[v]) {
                dist[v] = d;
                prev[v] = it.second;
                open.push({d, v});
            }
        }
    }
    if (dist[g] >= 1e300) return false;
    std::vector<int> chain;
    for (int v = g; v >= 0; v = prev[v]) chain.push_back(v);
    std::reverse(chain.begin(), chain.end());
    out.clear();
    TracePoint start;
    start.pos = surfacePosition(m, a);
    start.tri = a.tri;
    out.push_back(start);
    for (int v : chain) {
        TracePoint p;
        p.pos = m.pos[v];
        p.tri = triOfEdgeStart[v];
        out.push_back(p);
    }
    TracePoint end;
    end.pos = surfacePosition(m, b);
    end.tri = b.tri;
    out.push_back(end);
    return true;
}

} // namespace

bool traceSurfacePath(const TriMesh& m, const TriangleNeighbours& nb, const SurfacePoint& a, const SurfacePoint& b, std::vector<TracePoint>& out,
                      bool& usedEdgePath) {
    usedEdgePath = false;
    std::vector<TracePoint> path;
    const Vec3 pa = surfacePosition(m, a), pb = surfacePosition(m, b);
    bool found = false;
    if (a.tri == b.tri) {
        TracePoint p0, p1;
        p0.pos = pa;
        p0.tri = a.tri;
        p1.pos = pb;
        p1.tri = b.tri;
        path = {p0, p1};
        found = true;
    } else {
        const Vec3 n = m.faceNormal(a.tri) + m.faceNormal(b.tri);
        Vec3 q = cross(pb - pa, n);
        const double ql = norm(q);
        if (ql > 1e-12 * std::max(1.0, norm(pb - pa))) {
            q = q * (1.0 / ql);
            std::vector<TracePoint> ahead, behind;
            const bool okAhead = walkPlane(m, nb, a, b, q, true, ahead);
            const bool okBehind = walkPlane(m, nb, a, b, q, false, behind);
            if (okAhead && (!okBehind || pathLength(ahead) <= pathLength(behind))) {
                path = ahead;
                found = true;
            } else if (okBehind) {
                path = behind;
                found = true;
            }
        }
        if (!found) {
            if (!edgePath(m, a, b, path)) return false;
            usedEdgePath = true;
        }
    }
    const size_t from = out.empty() ? 0 : 1;
    if (!out.empty()) out.back().tri = path.front().tri;  // the previous segment's end continues into this one
    out.insert(out.end(), path.begin() + from, path.end());
    return true;
}

MarkingTrace traceMarking(const TriMesh& m, const TriangleNeighbours& nb, const std::vector<SurfacePoint>& points, bool closed) {
    MarkingTrace tr;
    if (points.empty()) return tr;
    if (points.size() == 1) {
        TracePoint p;
        p.pos = surfacePosition(m, points[0]);
        p.tri = points[0].tri;
        tr.points.push_back(p);
        tr.controlAt.push_back(0);
        return tr;
    }
    const size_t segments = closed && points.size() > 2 ? points.size() : points.size() - 1;
    for (size_t k = 0; k < segments; ++k) {
        const SurfacePoint& a = points[k];
        const SurfacePoint& b = points[(k + 1) % points.size()];
        tr.controlAt.push_back(tr.points.empty() ? 0 : static_cast<int>(tr.points.size()) - 1);  // control point k
        bool edges = false;
        if (!traceSurfacePath(m, nb, a, b, tr.points, edges)) {
            tr.failedSegments.push_back(static_cast<int>(k));
            // keep the chain going: jump to b
            TracePoint p;
            p.pos = surfacePosition(m, b);
            p.tri = b.tri;
            tr.points.push_back(p);
        }
        tr.edgePathUsed = tr.edgePathUsed || edges;
    }
    if (!closed || points.size() <= 2) tr.controlAt.push_back(static_cast<int>(tr.points.size()) - 1);
    return tr;
}

namespace {

// A control point on the border: the loop edge it lies on (a -> b in the loop's walk) and its fraction from a, or the
// loop vertex it lies on (fraction 0 from that vertex). False when it is not on the border.
struct BorderPlace {
    int loop = -1, index = -1;   // the loop, and the position of the edge's first vertex in it
    double fraction = 0;
};

// The border loops of m as walked by its triangles, with each border vertex's (loop, position).
struct BorderLoops {
    std::vector<std::vector<int>> loops;
    std::map<int, std::pair<int, int>> where;   // vertex -> (loop, position)
};

BorderLoops findBorderLoops(const TriMesh& m) {
    std::map<std::pair<int, int>, int> uses;
    for (const std::array<int, 3>& t : m.tri)
        for (int k = 0; k < 3; ++k) ++uses[{std::min(t[k], t[(k + 1) % 3]), std::max(t[k], t[(k + 1) % 3])}];
    std::map<int, int> next;
    for (const std::array<int, 3>& t : m.tri)
        for (int k = 0; k < 3; ++k) {
            const int a = t[k], b = t[(k + 1) % 3];
            if (uses[{std::min(a, b), std::max(a, b)}] == 1) next[a] = b;
        }
    BorderLoops out;
    for (const std::pair<const int, int>& start : next) {
        if (out.where.count(start.first)) continue;
        std::vector<int> loop;
        for (int v = start.first; !out.where.count(v) && next.count(v); v = next[v]) {
            out.where[v] = {static_cast<int>(out.loops.size()), static_cast<int>(loop.size())};
            loop.push_back(v);
        }
        out.loops.push_back(loop);
    }
    return out;
}

bool borderPlace(const TriMesh& m, const BorderLoops& b, const SurfacePoint& p, BorderPlace& out) {
    if (p.tri < 0 || p.tri >= static_cast<int>(m.tri.size())) return false;
    const std::array<int, 3>& t = m.tri[p.tri];
    for (int k = 0; k < 3; ++k)
        if (p.bary[k] == 1) {   // on a vertex
            const auto it = b.where.find(t[k]);
            if (it == b.where.end()) return false;
            out.loop = it->second.first;
            out.index = it->second.second;
            out.fraction = 0;
            return true;
        }
    for (int z = 0; z < 3; ++z) {
        if (p.bary[z] != 0) continue;
        const int a = t[(z + 1) % 3], c = t[(z + 2) % 3];   // the edge opposite the zero weight
        const auto ia = b.where.find(a), ic = b.where.find(c);
        if (ia == b.where.end() || ic == b.where.end() || ia->second.first != ic->second.first) return false;
        const std::vector<int>& loop = b.loops[ia->second.first];
        const int n = static_cast<int>(loop.size());
        const int pa = ia->second.second, pc = ic->second.second;
        if ((pa + 1) % n == pc) {   // a -> c is a border edge of the walk
            out.loop = ia->second.first;
            out.index = pa;
            out.fraction = p.bary[(z + 2) % 3];
            return true;
        }
        if ((pc + 1) % n == pa) {
            out.loop = ia->second.first;
            out.index = pc;
            out.fraction = p.bary[(z + 1) % 3];
            return true;
        }
        return false;   // both on the loop but not an edge of it (a chord)
    }
    return false;
}

} // namespace

bool pointsOnOneBorderLoop(const TriMesh& m, const std::vector<SurfacePoint>& points, std::vector<int>& loop) {
    loop.clear();
    if (points.size() < 3) return false;
    const BorderLoops b = findBorderLoops(m);
    std::vector<double> along;   // each point's place along the loop (edge index + fraction)
    int which = -1;
    for (const SurfacePoint& p : points) {
        BorderPlace place;
        if (!borderPlace(m, b, p, place) || (which >= 0 && place.loop != which)) return false;
        which = place.loop;
        along.push_back(place.index + place.fraction);
    }
    // the direction the points go round: the walk's when its forward steps are shorter in total than its backward ones
    const double n = static_cast<double>(b.loops[which].size());
    double forward = 0, backward = 0;
    for (size_t i = 0; i < along.size(); ++i) {
        const double d = std::fmod(along[(i + 1) % along.size()] - along[i] + 2 * n, n);
        forward += d;
        backward += d > 0 ? n - d : 0;
    }
    loop = b.loops[which];
    if (backward < forward) std::reverse(loop.begin(), loop.end());
    return true;
}

MarkingTrace borderLoopTrace(const TriMesh& m, const std::vector<int>& loop, const std::vector<SurfacePoint>& points) {
    MarkingTrace tr;
    const int n = static_cast<int>(loop.size());
    if (n < 3 || points.empty()) return tr;
    // every control point's place along this loop (in the loop's own order)
    std::map<int, int> position;
    for (int i = 0; i < n; ++i) position[loop[i]] = i;
    std::vector<double> along;
    for (const SurfacePoint& p : points) {
        const std::array<int, 3>& t = m.tri[p.tri];
        double best = 0;
        bool found = false;
        for (int k = 0; k < 3 && !found; ++k)
            if (p.bary[k] == 1 && position.count(t[k])) {
                best = position[t[k]];
                found = true;
            }
        for (int z = 0; z < 3 && !found; ++z) {
            if (p.bary[z] != 0) continue;
            const int a = t[(z + 1) % 3], c = t[(z + 2) % 3];
            if (!position.count(a) || !position.count(c)) continue;
            const int pa = position[a], pc = position[c];
            if ((pa + 1) % n == pc) {
                best = pa + p.bary[(z + 2) % 3];
                found = true;
            } else if ((pc + 1) % n == pa) {
                best = pc + p.bary[(z + 1) % 3];
                found = true;
            }
        }
        along.push_back(best);
    }
    auto edgeTriangle = [&](int a, int b) {
        for (int t = 0; t < static_cast<int>(m.tri.size()); ++t)
            for (int k = 0; k < 3; ++k)
                if ((m.tri[t][k] == a && m.tri[t][(k + 1) % 3] == b) || (m.tri[t][k] == b && m.tri[t][(k + 1) % 3] == a)) return t;
        return -1;
    };
    // walk once round from the first control point's edge; insert the control points in order as their places pass
    const int startEdge = static_cast<int>(along[0]);
    size_t nextPoint = 0;
    for (int s = 0; s <= n; ++s) {
        const int e = (startEdge + s) % n;
        const int a = loop[e], b = loop[(e + 1) % n];
        const int tri = edgeTriangle(a, b);
        if (s > 0) {
            TracePoint v;
            v.pos = m.pos[a];
            v.tri = tri;
            tr.points.push_back(v);
        }
        if (s == n) {
            // back on the first edge: a last control point before the first one's fraction
            for (; nextPoint < along.size(); ++nextPoint) {
                const double f = along[nextPoint] - static_cast<int>(along[nextPoint]);
                TracePoint c;
                c.pos = (1 - f) * m.pos[a] + f * m.pos[b];
                c.tri = tri;
                tr.controlAt.push_back(static_cast<int>(tr.points.size()));
                tr.points.push_back(c);
            }
            break;
        }
        // the control points on this edge, by fraction (a point on the first edge before its own fraction waits for the end)
        while (nextPoint < along.size()) {
            const double place = along[nextPoint];
            const int pe = static_cast<int>(place);
            if (pe != e) break;
            TracePoint c;
            c.pos = (1 - (place - pe)) * m.pos[a] + (place - pe) * m.pos[b];
            c.tri = tri;
            tr.controlAt.push_back(static_cast<int>(tr.points.size()));
            tr.points.push_back(c);
            ++nextPoint;
        }
    }
    if (!tr.points.empty()) tr.points.push_back(tr.points[0]);   // closed: back to the first control point
    return tr;
}

} // namespace knit
