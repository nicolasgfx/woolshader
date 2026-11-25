#include "stitchmesh/relax_mesh.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <utility>

#include "model/cancel.h"
#include "model/mesh_query.h"
#include "stitchmesh/typed_mesh.h"

namespace knit {

namespace {

using EdgeKey = std::pair<int, int>;
EdgeKey keyOf(int a, int b) { return {std::min(a, b), std::max(a, b)}; }

struct Spring {
    int i = -1, j = -1;
    double rest = 1, k = 2;
};
struct Corner {
    int i = -1, j = -1, k = -1;   // shear at j between i and k
};
struct Strut {
    int i = -1, j = -1, k = -1;   // wale edges i-j and j-k
};

struct Model {
    std::vector<Spring> springs;
    std::vector<Corner> corners;
    std::vector<Strut> struts;
    double rWale = 1;
};

// The elements of Eqs. 2-4 in units of r_course (S3).
Model buildModel(const std::vector<RelaxFace>& faces, const MeshRelaxRequest& r, RelaxProduct* counts, const std::vector<std::array<int, 2>>* courseEdges = nullptr) {
    Model m;
    const double rw = r.h / r.w;
    m.rWale = rw;
    std::map<EdgeKey, StitchEdgeRole> edges;
    int conflicts = 0;
    for (const RelaxFace& f : faces) {
        const int n = static_cast<int>(f.corners.size());
        for (int k = 0; k < n; ++k) {
            const int a = f.corners[k], b = f.corners[(k + 1) % n];
            if (a == b) continue;
            const bool wale = f.roles[k] == StitchEdgeRole::Wale;
            const auto it = edges.find(keyOf(a, b));
            if (it == edges.end()) edges[keyOf(a, b)] = wale ? StitchEdgeRole::Wale : StitchEdgeRole::Bottom;
            else if ((it->second == StitchEdgeRole::Wale) != wale) ++conflicts;
        }
    }
    int added = 0;
    if (courseEdges)
        for (const std::array<int, 2>& e : *courseEdges)
            if (e[0] != e[1] && edges.find(keyOf(e[0], e[1])) == edges.end()) {
                edges[keyOf(e[0], e[1])] = StitchEdgeRole::Bottom;
                ++added;
            }
    int course = 0, waleN = 0;
    for (const auto& [e, role] : edges) {
        const bool wale = role == StitchEdgeRole::Wale;
        m.springs.push_back(Spring{e.first, e.second, wale ? rw : 1.0, r.kStretch});
        (wale ? waleN : course) += 1;
    }
    // sub-quads: bottom course edges paired with top course edges (paper Fig. 15)
    const double diag = std::sqrt(1.0 + rw * rw);
    int subQuads = 0, without = 0, diagonals = 0;
    for (const RelaxFace& f : faces) {
        const int n = static_cast<int>(f.corners.size());
        std::vector<std::pair<int, int>> bottoms, tops;
        for (int k = 0; k < n; ++k) {
            const int a = f.corners[k], b = f.corners[(k + 1) % n];
            if (a == b) continue;
            if (f.roles[k] == StitchEdgeRole::Bottom) bottoms.push_back({a, b});
            if (f.roles[k] == StitchEdgeRole::Top) tops.push_back({a, b});
        }
        int made = 0;
        if (!bottoms.empty() && !tops.empty() && (bottoms.size() == 1 || tops.size() == 1)) {
            for (const auto& bt : bottoms)
                for (const auto& tp : tops) {
                    const int p = bt.first, q = bt.second, rr = tp.first, s = tp.second;
                    if (p == rr || p == s || q == rr || q == s) continue;   // the course edges meet: a triangle, no sub-quad
                    m.springs.push_back(Spring{p, rr, diag, r.kStretch});
                    m.springs.push_back(Spring{q, s, diag, r.kStretch});
                    diagonals += 2;
                    m.corners.push_back(Corner{s, p, q});
                    m.corners.push_back(Corner{p, q, rr});
                    m.corners.push_back(Corner{q, rr, s});
                    m.corners.push_back(Corner{rr, s, p});
                    ++subQuads;
                    ++made;
                }
        }
        if (made == 0) ++without;
    }
    // wale struts: wale edges oriented bottom -> top by the faces that have their own course direction
    std::set<EdgeKey> directed;   // (bottom, top)
    for (const RelaxFace& f : faces) {
        const int n = static_cast<int>(f.corners.size());
        const int ex = exitWaleCorner(f.roles), en = entryWaleCorner(f.roles);
        if (ex >= 0 && f.corners[ex] != f.corners[(ex + 1) % n]) directed.insert({f.corners[ex], f.corners[(ex + 1) % n]});
        if (en >= 0 && f.corners[en] != f.corners[(en + 1) % n]) directed.insert({f.corners[(en + 1) % n], f.corners[en]});
    }
    std::map<int, std::vector<int>> below, above;   // vertex -> the other end
    for (const EdgeKey& e : directed) {
        below[e.second].push_back(e.first);
        above[e.first].push_back(e.second);
    }
    int noStrut = 0;
    for (const auto& [j, lows] : below) {
        const auto it = above.find(j);
        if (it == above.end()) continue;
        if (lows.size() == 1 && it->second.size() == 1) m.struts.push_back(Strut{lows[0], j, it->second[0]});
        else ++noStrut;
    }
    if (counts) {
        counts->courseSprings = course;
        counts->addedCourseSprings = added;
        counts->waleSprings = waleN;
        counts->diagonalSprings = diagonals;
        counts->shearCorners = static_cast<int>(m.corners.size());
        counts->struts = static_cast<int>(m.struts.size());
        counts->subQuads = subQuads;
        counts->facesWithoutSubQuads = without;
        counts->waleJunctionsWithoutStrut = noStrut;
        counts->conflictingEdgeRoles = conflicts;
    }
    return m;
}

void forces(const Model& m, const MeshRelaxRequest& r, const std::vector<Vec3>& x, std::vector<Vec3>& F) {
    F.assign(x.size(), Vec3{0, 0, 0});
    for (const Spring& s : m.springs) {
        const Vec3 d = x[s.i] - x[s.j];
        const double L = norm(d);
        if (L < 1e-12) continue;
        const Vec3 f = d * (s.k * (s.rest / L - 1.0) / L);
        F[s.i] += f;
        F[s.j] -= f;
    }
    for (const Corner& c : m.corners) {
        const Vec3 a = x[c.i] - x[c.j], b = x[c.k] - x[c.j];
        const Vec3 fj = (x[c.j] - (x[c.i] + x[c.k]) * 0.5) * (-r.kShear * dot(a, b));
        F[c.j] += fj;
        F[c.i] -= fj * 0.5;
        F[c.k] -= fj * 0.5;
    }
    for (const Strut& s : m.struts) {
        const Vec3 d = x[s.i] - x[s.k];
        const double L = norm(d);
        if (L < 1e-12) continue;
        const double rr = std::max(m.rWale, norm(x[s.i] - x[s.j]) + norm(x[s.k] - x[s.j]));
        const Vec3 f = d * (-r.kWale * (L / rr - 1.0) / L);
        F[s.i] += f;
        F[s.k] -= f;
    }
}

// A symmetric 3 x 3 block.
struct Block {
    double a[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    Vec3 apply(const Vec3& v) const {
        return Vec3{a[0][0] * v.x + a[0][1] * v.y + a[0][2] * v.z, a[1][0] * v.x + a[1][1] * v.y + a[1][2] * v.z, a[2][0] * v.x + a[2][1] * v.y + a[2][2] * v.z};
    }
};
struct Pair {
    int i = -1, j = -1;
    Block K;   // df_i / dx_i of the pair's force (negative semidefinite, S4)
};

Block outer(const Vec3& u, double c) {
    Block b;
    const double v[3] = {u.x, u.y, u.z};
    for (int p = 0; p < 3; ++p)
        for (int q = 0; q < 3; ++q) b.a[p][q] = c * v[p] * v[q];
    return b;
}

std::vector<Pair> jacobian(const Model& m, const MeshRelaxRequest& r, const std::vector<Vec3>& x) {
    std::vector<Pair> P;
    for (const Spring& s : m.springs) {
        const Vec3 d = x[s.i] - x[s.j];
        const double L = norm(d);
        if (L < 1e-12) continue;
        const Vec3 u = d / L;
        Block K = outer(u, -s.k * s.rest / (L * L));
        const double c2 = std::min(0.0, s.k * (s.rest / L - 1.0) / L);   // transverse term only while stretched (S4)
        for (int p = 0; p < 3; ++p) K.a[p][p] += c2;
        const Block uu = outer(u, -c2);
        for (int p = 0; p < 3; ++p)
            for (int q = 0; q < 3; ++q) K.a[p][q] += uu.a[p][q];
        P.push_back(Pair{s.i, s.j, K});
    }
    for (const Strut& s : m.struts) {
        const Vec3 d = x[s.i] - x[s.k];
        const double L = norm(d);
        if (L < 1e-12) continue;
        const double rr = std::max(m.rWale, norm(x[s.i] - x[s.j]) + norm(x[s.k] - x[s.j]));
        P.push_back(Pair{s.i, s.k, outer(d / L, -r.kWale / rr)});
    }
    return P;
}

// --- the surface: triangle adjacency, the directional walk, the border loops ------------------------------------------

struct Surface {
    const TriMesh* m = nullptr;
    std::map<EdgeKey, std::vector<int>> edgeTris;
    std::vector<Vec3> normal;
};

Surface surfaceOf(const TriMesh& m) {
    Surface s;
    s.m = &m;
    for (int t = 0; t < static_cast<int>(m.tri.size()); ++t) {
        for (int k = 0; k < 3; ++k) s.edgeTris[keyOf(m.tri[t][k], m.tri[t][(k + 1) % 3])].push_back(t);
        s.normal.push_back(m.faceNormal(t));
    }
    return s;
}

Vec3 pointOf(const TriMesh& m, int tri, const double b[3]) {
    return m.pos[m.tri[tri][0]] * b[0] + m.pos[m.tri[tri][1]] * b[1] + m.pos[m.tri[tri][2]] * b[2];
}

// Moves the bound point by d [mm] along the surface: straight in each triangle's plane, the remainder carried across an
// edge into the neighbour's plane. Stops at the reference border (a clamp).
void walk(const Surface& s, VertexBinding& b, Vec3 d, int& clamps) {
    const TriMesh& m = *s.m;
    for (int crossing = 0; crossing < 1000; ++crossing) {
        const Vec3 n = s.normal[b.tri];
        d = d - n * dot(n, d);
        if (norm(d) < 1e-14) return;
        const Vec3 V0 = m.pos[m.tri[b.tri][0]], V1 = m.pos[m.tri[b.tri][1]], V2 = m.pos[m.tri[b.tri][2]];
        const Vec3 e1 = V1 - V0, e2 = V2 - V0;
        const double g11 = dot(e1, e1), g12 = dot(e1, e2), g22 = dot(e2, e2), r1 = dot(e1, d), r2 = dot(e2, d);
        const double det = g11 * g22 - g12 * g12;
        if (std::fabs(det) < 1e-30) return;
        const double al = (r1 * g22 - r2 * g12) / det, be = (r2 * g11 - r1 * g12) / det;
        const double db[3] = {-al - be, al, be};
        double smin = 1.0;
        int exitIdx = -1;
        for (int i = 0; i < 3; ++i)
            if (db[i] < -1e-18) {
                const double sv = -b.bary[i] / db[i];
                if (sv < smin) {
                    smin = std::max(0.0, sv);
                    exitIdx = i;
                }
            }
        if (exitIdx < 0) {
            for (int i = 0; i < 3; ++i) b.bary[i] += db[i];
            return;
        }
        for (int i = 0; i < 3; ++i) b.bary[i] = std::max(0.0, b.bary[i] + smin * db[i]);
        b.bary[exitIdx] = 0;
        const double sum = b.bary[0] + b.bary[1] + b.bary[2];
        for (int i = 0; i < 3; ++i) b.bary[i] /= sum;
        const Vec3 rem = d * (1.0 - smin);
        const int va = m.tri[b.tri][(exitIdx + 1) % 3], vb = m.tri[b.tri][(exitIdx + 2) % 3];
        const auto it = s.edgeTris.find(keyOf(va, vb));
        int nt = -1;
        if (it != s.edgeTris.end())
            for (int t : it->second)
                if (t != b.tri) nt = t;
        if (nt < 0) {
            ++clamps;
            return;
        }
        const Vec3 P = pointOf(m, b.tri, b.bary);
        const Vec3 e = normalized(m.pos[vb] - m.pos[va]);
        int other = m.tri[nt][0];
        for (int k = 0; k < 3; ++k)
            if (m.tri[nt][k] != va && m.tri[nt][k] != vb) other = m.tri[nt][k];
        Vec3 into = m.pos[other] - m.pos[va];
        into = normalized(into - e * dot(into, e));
        const Vec3 perp = rem - e * dot(rem, e);
        d = e * dot(rem, e) + into * norm(perp);
        const TrianglePoint tp = closestOnTriangle(P, m.pos[m.tri[nt][0]], m.pos[m.tri[nt][1]], m.pos[m.tri[nt][2]]);
        b.tri = nt;
        for (int k = 0; k < 3; ++k) b.bary[k] = tp.bary[k];
    }
}

struct BorderSlot {
    int loop = -1;
    double s = 0;   // arc parameter along the loop [mm]
};

struct Loop {
    std::vector<int> v;
    std::vector<double> cum;   // cum[i]: arc length at vertex i; cum[n]: the loop length
};

void placeOnLoop(const Surface& sf, const Loop& L, double s, Vec3& p, Vec3& tangent, VertexBinding& b) {
    const TriMesh& m = *sf.m;
    const int n = static_cast<int>(L.v.size());
    const double len = L.cum[n];
    s = std::fmod(s, len);
    if (s < 0) s += len;
    int i = static_cast<int>(std::upper_bound(L.cum.begin(), L.cum.end(), s) - L.cum.begin()) - 1;
    i = std::clamp(i, 0, n - 1);
    const int a = L.v[i], c = L.v[(i + 1) % n];
    const double seg = L.cum[i + 1] - L.cum[i];
    const double t = seg > 0 ? (s - L.cum[i]) / seg : 0;
    p = m.pos[a] + (m.pos[c] - m.pos[a]) * t;
    tangent = normalized(m.pos[c] - m.pos[a]);
    const auto it = sf.edgeTris.find(keyOf(a, c));
    b.tri = it == sf.edgeTris.end() ? -1 : it->second.front();
    b.distance = 0;
    if (b.tri >= 0)
        for (int k = 0; k < 3; ++k) b.bary[k] = m.tri[b.tri][k] == a ? 1 - t : (m.tri[b.tri][k] == c ? t : 0);
}

// A guide line's arc parameterization and the point at arc parameter s (an open line clamps s to its ends).
struct GuideTrack {
    const GuideLine* line = nullptr;
    std::vector<double> cum;   // cum[i]: arc length at point i; the last entry the line length
};

GuideTrack trackOf(const GuideLine& g) {
    GuideTrack t;
    t.line = &g;
    const int n = static_cast<int>(g.points.size());
    const int segs = g.closed ? n : n - 1;
    t.cum.push_back(0);
    for (int i = 0; i < segs; ++i) t.cum.push_back(t.cum.back() + norm(g.points[(i + 1) % n] - g.points[i]));
    return t;
}

void placeOnGuide(const TriMesh& m, const GuideTrack& t, double& s, Vec3& p, Vec3& tangent, VertexBinding& b) {
    const GuideLine& g = *t.line;
    const int n = static_cast<int>(g.points.size());
    const double len = t.cum.back();
    if (g.closed) {
        s = std::fmod(s, len);
        if (s < 0) s += len;
    } else {
        s = std::clamp(s, 0.0, len);
    }
    const int segs = static_cast<int>(t.cum.size()) - 1;
    int i = static_cast<int>(std::upper_bound(t.cum.begin(), t.cum.end(), s) - t.cum.begin()) - 1;
    i = std::clamp(i, 0, segs - 1);
    const Vec3 a = g.points[i], c = g.points[(i + 1) % n];
    const double seg = t.cum[i + 1] - t.cum[i];
    const double u = seg > 0 ? (s - t.cum[i]) / seg : 0;
    p = a + (c - a) * u;
    tangent = normalized(c - a);
    b.tri = g.segmentTri[i];
    b.distance = 0;
    if (b.tri >= 0) {
        const TrianglePoint tp = closestOnTriangle(p, m.pos[m.tri[b.tri][0]], m.pos[m.tri[b.tri][1]], m.pos[m.tri[b.tri][2]]);
        for (int k = 0; k < 3; ++k) b.bary[k] = tp.bary[k];
    }
}

} // namespace

std::vector<RelaxFace> relaxFacesOf(const PolygonMesh& m, const std::vector<std::vector<StitchEdgeRole>>& role) {
    std::vector<RelaxFace> faces;
    for (int f = 0; f < m.faces(); ++f) {
        RelaxFace rf;
        for (int k = 0; k < m.sides(f); ++k) rf.corners.push_back(m.vertexAt(f, k));
        rf.roles = role[f];
        faces.push_back(rf);
    }
    return faces;
}

std::vector<Vec3> relaxForces(const std::vector<RelaxFace>& faces, const std::vector<Vec3>& x, const MeshRelaxRequest& r) {
    const Model m = buildModel(faces, r, nullptr);
    std::vector<Vec3> F;
    forces(m, r, x, F);
    return F;
}

bool relaxStitchMesh(const TriMesh& reference, const std::vector<std::vector<int>>& refLoops, const std::vector<RelaxFace>& faces,
                     const std::vector<Vec3>& start, const std::vector<char>& border, const MeshRelaxRequest& r, RelaxProduct& out, std::string* err,
                     const RelaxGuides* guides) {
    out = RelaxProduct{};
    const auto t0 = std::chrono::steady_clock::now();
    const int nV = static_cast<int>(start.size());
    if (static_cast<int>(border.size()) != nV || faces.empty()) {
        if (err) *err = "relaxStitchMesh: sizes disagree";
        return false;
    }
    const Model model = buildModel(faces, r, &out, guides ? &guides->courseEdges : nullptr);
    const Surface sf = surfaceOf(reference);
    std::vector<Loop> loops;
    for (const std::vector<int>& l : refLoops) {
        Loop L;
        L.v = l;
        L.cum.push_back(0);
        for (size_t i = 0; i < l.size(); ++i) L.cum.push_back(L.cum.back() + norm(reference.pos[l[(i + 1) % l.size()]] - reference.pos[l[i]]));
        loops.push_back(L);
    }
    // initial bindings: interior vertices by the closest point, border vertices snapped onto the nearest loop
    const MeshQuery q(reference);
    std::vector<Vec3> pos = start;
    std::vector<VertexBinding> bind(nV);
    std::vector<BorderSlot> slot(nV);
    std::vector<Vec3> tangent(nV);
    // the markings' constraints: pinned vertices and vertices on guide lines (their arc parameter in guideS)
    std::vector<GuideTrack> tracks;
    std::vector<int> lineOf(nV, -1);
    std::vector<char> pinned(nV, 0);
    std::vector<double> guideS(nV, 0);
    if (guides) {
        for (const GuideLine& g : guides->lines) tracks.push_back(trackOf(g));
        for (int v = 0; v < nV; ++v) {
            if (v < static_cast<int>(guides->pinned.size()) && guides->pinned[v]) pinned[v] = 1;
            else if (v < static_cast<int>(guides->lineOf.size())) lineOf[v] = guides->lineOf[v];
        }
    }
    // the place of a border vertex: the closest point of the reference loops, as (loop, arc parameter); returns the distance
    auto nearestSlot = [&](int v) {
        double best = 1e300;
        for (int li = 0; li < static_cast<int>(loops.size()); ++li) {
            const Loop& L = loops[li];
            const int n = static_cast<int>(L.v.size());
            for (int i = 0; i < n; ++i) {
                const Vec3 a = reference.pos[L.v[i]], c = reference.pos[L.v[(i + 1) % n]];
                const ClosestSegmentPoints cp = closestSegmentPoints(start[v], start[v], a, c);
                if (cp.distance < best) {
                    best = cp.distance;
                    slot[v].loop = li;
                    slot[v].s = L.cum[i] + cp.t * (L.cum[i + 1] - L.cum[i]);
                }
            }
        }
        return best;
    };
    for (int v = 0; v < nV; ++v) {
        if (lineOf[v] >= 0) {
            // the closest point of the line
            const GuideTrack& t = tracks[lineOf[v]];
            const GuideLine& g = *t.line;
            const int n = static_cast<int>(g.points.size());
            double best = 1e300;
            for (int i = 0; i + 1 < static_cast<int>(t.cum.size()); ++i) {
                const ClosestSegmentPoints cp = closestSegmentPoints(start[v], start[v], g.points[i], g.points[(i + 1) % n]);
                if (cp.distance < best) {
                    best = cp.distance;
                    guideS[v] = t.cum[i] + cp.t * (t.cum[i + 1] - t.cum[i]);
                }
            }
            placeOnGuide(reference, t, guideS[v], pos[v], tangent[v], bind[v]);
            ++out.guidedVertices;
            out.maxGuideSnap = std::max(out.maxGuideSnap, best);
            continue;
        }
        if (pinned[v]) {
            ++out.pinnedVertices;
            if (border[v]) nearestSlot(v);   // a held border vertex keeps its place in the order around its loop
            const MeshQuery::Hit h = q.closest(start[v]);
            bind[v].tri = h.tri;
            const TrianglePoint tp = closestOnTriangle(start[v], reference.pos[reference.tri[h.tri][0]], reference.pos[reference.tri[h.tri][1]],
                                                       reference.pos[reference.tri[h.tri][2]]);
            for (int k = 0; k < 3; ++k) bind[v].bary[k] = tp.bary[k];
            pos[v] = tp.p;
            continue;
        }
        if (border[v]) {
            ++out.borderVertices;
            const double best = nearestSlot(v);
            if (slot[v].loop < 0) {
                if (err) *err = "relaxStitchMesh: a border vertex without a reference loop";
                return false;
            }
            placeOnLoop(sf, loops[slot[v].loop], slot[v].s, pos[v], tangent[v], bind[v]);
            out.maxBorderSnap = std::max(out.maxBorderSnap, best);
            if (best > 1e-9) ++out.borderVerticesSnapped;
        } else {
            const MeshQuery::Hit h = q.closest(start[v]);
            bind[v].tri = h.tri;
            const TrianglePoint tp = closestOnTriangle(start[v], reference.pos[reference.tri[h.tri][0]], reference.pos[reference.tri[h.tri][1]],
                                                       reference.pos[reference.tri[h.tri][2]]);
            for (int k = 0; k < 3; ++k) bind[v].bary[k] = tp.bary[k];
            pos[v] = tp.p;
        }
    }
    // the order of the border vertices around each loop, as the layout gave it: nextOnLoop[v] is the next vertex in the
    // direction of growing arc parameter (-1: not on a loop, or fewer than three vertices on it)
    std::vector<int> nextOnLoop(nV, -1), prevOnLoop(nV, -1);
    for (int li = 0; li < static_cast<int>(loops.size()); ++li) {
        std::vector<int> on;
        for (int v = 0; v < nV; ++v)
            if (border[v] && lineOf[v] < 0 && slot[v].loop == li) on.push_back(v);
        std::sort(on.begin(), on.end(), [&](int a, int b) { return slot[a].s < slot[b].s; });
        for (size_t k = 0; on.size() > 2 && k < on.size(); ++k) {
            nextOnLoop[on[k]] = on[(k + 1) % on.size()];
            prevOnLoop[on[(k + 1) % on.size()]] = on[k];
        }
    }
    // the gap from v to its next vertex along the loop [mm]: positive while the two are in order
    auto gapToNext = [&](int v) {
        const double len = loops[slot[v].loop].cum.back();
        double gap = std::fmod(slot[nextOnLoop[v]].s - slot[v].s, len);
        if (gap < -0.5 * len) gap += len;
        if (gap > 0.5 * len) gap -= len;
        return gap;
    };
    std::vector<double> gapAhead(nV, 0);
    // every corner (a, b, c) of the solver's faces, a repeated corner of a start step left out
    std::vector<std::array<int, 3>> corners;
    for (const RelaxFace& f : faces) {
        std::vector<int> ring;
        for (int v : f.corners)
            if (ring.empty() || (ring.back() != v && ring.front() != v)) ring.push_back(v);
        const int n = static_cast<int>(ring.size());
        for (int k = 0; n >= 3 && k < n; ++k) corners.push_back({ring[k], ring[(k + 1) % n], ring[(k + 2) % n]});
    }
    // the filter of each vertex (Baraff and Witkin): the tangent plane, or the border tangent
    std::vector<Vec3> vnormal = reference.vnormal;
    if (vnormal.size() != reference.pos.size()) {
        TriMesh copy = reference;
        copy.computeNormals();
        vnormal = copy.vnormal;
    }
    auto normalAt = [&](const VertexBinding& b) {
        if (!r.smoothNormals) return sf.normal[b.tri];
        const auto& t = reference.tri[b.tri];
        const Vec3 n = vnormal[t[0]] * b.bary[0] + vnormal[t[1]] * b.bary[1] + vnormal[t[2]] * b.bary[2];
        const double len = norm(n);
        return len > 0 ? n / len : sf.normal[b.tri];
    };
    std::vector<char> kinkMode(nV, 0);
    std::vector<Vec3> kinkDir(nV);
    std::vector<std::array<int, 2>> kinkEdge(nV, std::array<int, 2>{-1, -1});
    auto filter = [&](int v, const Vec3& a) {
        if (pinned[v] || kinkMode[v] == 2) return Vec3{0, 0, 0};
        if (kinkMode[v] == 1) return kinkDir[v] * dot(kinkDir[v], a);
        if (lineOf[v] >= 0 || border[v]) return tangent[v] * dot(tangent[v], a);
        const Vec3 n = normalAt(bind[v]);
        return a - n * dot(n, a);
    };
    int dofs = 0;
    for (int v = 0; v < nV; ++v) dofs += pinned[v] ? 0 : (lineOf[v] >= 0 || border[v]) ? 1 : 2;
    double hstep = r.timestep;
    std::vector<Vec3> keptPos, keptTangent;       // StepControl::ForcesMayNotGrow: the last state that was kept
    std::vector<VertexBinding> keptBind;
    std::vector<BorderSlot> keptSlot;
    std::vector<double> keptGuideS;
    double keptRms = 0;
    std::vector<Vec3> x(nV), F, rhs(nV), dx(nV), res(nV), z(nV), p(nV), qv(nV);
    std::vector<double> rmsHistory;
    out.stopReason = "step cap";
    // R2: per step, the vertices held at a kink (1: on a crease edge, moving along kinkDir; 2: held at a polyline corner)
    auto polyCorner = [&](const std::vector<Vec3>& pts, bool closed, const Vec3& p, const Vec3& f) {
        const int n = static_cast<int>(pts.size());
        for (int i = 0; i < n; ++i) {
            if (norm(pts[i] - p) > r.kinkTolerance * r.w) continue;
            if (!closed && (i == 0 || i == n - 1)) return false;
            const Vec3 tNext = normalized(pts[(i + 1) % n] - pts[i]), tPrev = normalized(pts[i] - pts[(i + n - 1) % n]);
            return dot(f, tNext) <= 0 && dot(f, tPrev) >= 0;
        }
        return false;
    };
    auto chooseKinks = [&](const std::vector<Vec3>& Fn) {
        out.kinkHeldLast = out.kinkOnEdgeLast = 0;
        for (int v = 0; v < nV; ++v) {
            kinkMode[v] = 0;
            if (!r.kinkConstraints || pinned[v]) continue;
            const Vec3& f = Fn[v];
            if (lineOf[v] >= 0) {
                if (polyCorner(tracks[lineOf[v]].line->points, tracks[lineOf[v]].line->closed, pos[v], f)) kinkMode[v] = 2;
            } else if (border[v]) {
                std::vector<Vec3> pts;
                for (int i : loops[slot[v].loop].v) pts.push_back(reference.pos[i]);
                if (polyCorner(pts, true, pos[v], f)) kinkMode[v] = 2;
            } else if (bind[v].tri >= 0) {
                const auto& tr = reference.tri[bind[v].tri];
                double best = r.kinkTolerance * r.w;
                for (int k = 0; k < 3; ++k) {
                    const int a = tr[k], b = tr[(k + 1) % 3];
                    const double d = closestSegmentPoints(pos[v], pos[v], reference.pos[a], reference.pos[b]).distance;
                    if (d > best) continue;
                    const auto it = sf.edgeTris.find(keyOf(a, b));
                    if (it == sf.edgeTris.end() || it->second.size() != 2) continue;
                    const Vec3 l = normalized(reference.pos[b] - reference.pos[a]);
                    bool both = true;
                    for (int t : it->second) {
                        int c = -1;
                        for (int q = 0; q < 3; ++q)
                            if (reference.tri[t][q] != a && reference.tri[t][q] != b) c = reference.tri[t][q];
                        Vec3 into = reference.pos[c] - reference.pos[a];
                        into = normalized(into - l * dot(l, into));
                        both = both && dot(f, into) <= 0;
                    }
                    if (!both) continue;
                    best = d;
                    kinkMode[v] = 1;
                    kinkDir[v] = l;
                    kinkEdge[v] = {a, b};
                }
            }
            if (kinkMode[v] == 1) ++out.kinkOnEdgeLast;
            if (kinkMode[v] == 2) ++out.kinkHeldLast;
        }
    };
    for (int step = 0; step <= r.maxSteps; ++step) {
        cancelCheckpoint();   // a superseded build stops within one step (without a token: nothing)
        double mx = 0;
        int mxV = -1;
        // the forces of the present state, the kinks they choose, and the RMS of the filtered force
        auto evaluate = [&]() {
            for (int v = 0; v < nV; ++v) x[v] = pos[v] / r.w;
            forces(model, r, x, F);
            chooseKinks(F);
            double sum = 0;
            mx = 0;
            mxV = -1;
            for (int v = 0; v < nV; ++v) {
                const Vec3 g = filter(v, F[v]);
                sum += dot(g, g);
                if (norm(g) > mx) {
                    mx = norm(g);
                    mxV = v;
                }
            }
            return std::sqrt(sum / std::max(1, dofs));
        };
        double rms = evaluate();
        if (r.stepControl == StepControl::ForcesMayNotGrow) {
            // a step after which the forces are larger is taken back and the time step halved (down to a 64th); a step
            // that lowers them is kept and the time step grows again by a quarter, up to the request's
            if (step > 0 && rms > keptRms && hstep > r.timestep / 64) {
                pos = keptPos;
                bind = keptBind;
                slot = keptSlot;
                guideS = keptGuideS;
                tangent = keptTangent;
                hstep *= 0.5;
                ++out.stepsTakenBack;
                rms = evaluate();
            } else {
                keptPos = pos;
                keptBind = bind;
                keptSlot = slot;
                keptGuideS = guideS;
                keptTangent = tangent;
                keptRms = rms;
                if (step > 0) hstep = std::min(r.timestep, 1.25 * hstep);
            }
        }
        if (step == 0) out.rmsForceStart = rms;
        out.rmsForceEnd = rms;
        RelaxStep rec;
        rec.step = step;
        rec.rmsForce = rms;
        rec.maxForce = mx;
        rec.maxForceVertex = mxV;
        if (mxV >= 0 && bind[mxV].tri >= 0) {
            rec.maxForceTri = bind[mxV].tri;
            const auto& tr = reference.tri[bind[mxV].tri];
            double e = 1e300;
            for (int k = 0; k < 3; ++k) e = std::min(e, closestSegmentPoints(pos[mxV], pos[mxV], reference.pos[tr[k]], reference.pos[tr[(k + 1) % 3]]).distance);
            rec.maxForceEdgeDistance = e;
        }
        rmsHistory.push_back(rms);
        const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        bool stop = false;
        if (!std::isfinite(rms) || rms > 100 * out.rmsForceStart + 1) {
            out.stopReason = "diverged";
            stop = true;
        } else if (rms < r.forceTolerance) {
            out.stopReason = "converged";
            stop = true;
        } else if (step >= r.stallWindow && rms > r.stallRatio * rmsHistory[step - r.stallWindow]) {
            out.stopReason = "stalled";
            stop = true;
        } else if (step == r.maxSteps) {
            out.stopReason = "step cap";
            stop = true;
        } else if (elapsed > r.capSeconds) {
            out.stopReason = "time cap";
            stop = true;
        }
        if (stop) {
            out.history.push_back(rec);
            out.steps = step;
            break;
        }
        // (I - h J) dx = h f, filtered PCG
        const std::vector<Pair> J = jacobian(model, r, x);
        std::vector<Vec3> diagA(nV, Vec3{1, 1, 1});
        for (const Pair& pr : J)
            for (int c = 0; c < 3; ++c) {
                const double d = -hstep * pr.K.a[c][c];
                (c == 0 ? diagA[pr.i].x : c == 1 ? diagA[pr.i].y : diagA[pr.i].z) += d;
                (c == 0 ? diagA[pr.j].x : c == 1 ? diagA[pr.j].y : diagA[pr.j].z) += d;
            }
        auto applyA = [&](const std::vector<Vec3>& in, std::vector<Vec3>& o) {
            for (int v = 0; v < nV; ++v) o[v] = in[v];
            for (const Pair& pr : J) {
                const Vec3 t = pr.K.apply(in[pr.i] - in[pr.j]) * hstep;
                o[pr.i] -= t;
                o[pr.j] += t;
            }
            for (int v = 0; v < nV; ++v) o[v] = filter(v, o[v]);
        };
        auto precond = [&](const std::vector<Vec3>& in, std::vector<Vec3>& o) {
            for (int v = 0; v < nV; ++v) o[v] = filter(v, Vec3{in[v].x / diagA[v].x, in[v].y / diagA[v].y, in[v].z / diagA[v].z});
        };
        double bnorm = 0;
        for (int v = 0; v < nV; ++v) {
            rhs[v] = filter(v, F[v] * hstep);
            dx[v] = Vec3{0, 0, 0};
            res[v] = rhs[v];
            bnorm += dot(rhs[v], rhs[v]);
        }
        bnorm = std::sqrt(bnorm);
        precond(res, z);
        p = z;
        double rz = 0;
        for (int v = 0; v < nV; ++v) rz += dot(res[v], z[v]);
        int it = 0;
        for (; it < r.cgMaxIterations; ++it) {
            applyA(p, qv);
            double pq = 0;
            for (int v = 0; v < nV; ++v) pq += dot(p[v], qv[v]);
            if (pq <= 0) break;
            const double alpha = rz / pq;
            double rn = 0;
            for (int v = 0; v < nV; ++v) {
                dx[v] += p[v] * alpha;
                res[v] -= qv[v] * alpha;
                rn += dot(res[v], res[v]);
            }
            if (std::sqrt(rn) <= r.cgTolerance * bnorm) {
                ++it;
                break;
            }
            precond(res, z);
            double rz2 = 0;
            for (int v = 0; v < nV; ++v) rz2 += dot(res[v], z[v]);
            const double beta = rz2 / rz;
            rz = rz2;
            for (int v = 0; v < nV; ++v) p[v] = filter(v, z[v] + p[v] * beta);
        }
        rec.cgIterations = it;
        // move on the surface
        for (int v = 0; v < nV; ++v)
            if (nextOnLoop[v] >= 0) gapAhead[v] = gapToNext(v);
        double maxMove = 0;
        for (int v = 0; v < nV; ++v) maxMove = std::max(maxMove, std::min(norm(dx[v]), r.maxStep));
        // BorderOrder::RingSlidesAsOne: the common step of each loop's ring [mm]
        std::vector<double> ringStep(loops.size(), 0);
        if (r.borderOrder == BorderOrder::RingSlidesAsOne) {
            std::vector<int> count(loops.size(), 0);
            std::vector<char> held(loops.size(), 0);
            for (int v = 0; v < nV; ++v) {
                if (nextOnLoop[v] < 0) continue;
                const int li = slot[v].loop;
                if (pinned[v] || kinkMode[v] == 2) held[li] = 1;
                Vec3 d = dx[v];
                const double len = norm(d);
                if (len > r.maxStep) d = d * (r.maxStep / len);
                ringStep[li] += dot(d, tangent[v]) * r.w;
                ++count[li];
            }
            for (size_t li = 0; li < loops.size(); ++li) ringStep[li] = held[li] || count[li] == 0 ? 0 : ringStep[li] / count[li];
        }
        // one vertex by `share` of its solved step, from its state at the start of the step
        auto moveVertex = [&](int v, double share) {
            Vec3 d = dx[v];
            const double len = norm(d);
            if (len > r.maxStep) d = d * (r.maxStep / len);
            d = d * share;
            if (pinned[v] || kinkMode[v] == 2) return;
            if (kinkMode[v] == 1) {
                // onto the crease edge, then along it
                const Vec3 A = reference.pos[kinkEdge[v][0]], B = reference.pos[kinkEdge[v][1]];
                const double u = std::clamp(dot(pos[v] - A, B - A) / dot(B - A, B - A), 0.0, 1.0);
                const int t = bind[v].tri;
                for (int k = 0; k < 3; ++k)
                    bind[v].bary[k] = reference.tri[t][k] == kinkEdge[v][0] ? 1 - u : reference.tri[t][k] == kinkEdge[v][1] ? u : 0;
                walk(sf, bind[v], kinkDir[v] * (dot(d, kinkDir[v]) * r.w), out.walkClamps);
                pos[v] = pointOf(reference, bind[v].tri, bind[v].bary);
                return;
            }
            if (lineOf[v] >= 0) {
                guideS[v] += dot(d, tangent[v]) * r.w;
                placeOnGuide(reference, tracks[lineOf[v]], guideS[v], pos[v], tangent[v], bind[v]);
            } else if (border[v]) {
                double along = dot(d, tangent[v]) * r.w;
                if (r.borderOrder == BorderOrder::RingSlidesAsOne && nextOnLoop[v] >= 0) along = ringStep[slot[v].loop] * share;
                if (r.borderOrder == BorderOrder::Kept && nextOnLoop[v] >= 0) {
                    const double ahead = 0.45 * std::max(0.0, gapAhead[v]), behind = 0.45 * std::max(0.0, gapAhead[prevOnLoop[v]]);
                    if (share == 1 && (along > ahead || along < -behind)) ++out.borderMovesLimited;
                    along = std::clamp(along, -behind, ahead);
                }
                slot[v].s += along;
                placeOnLoop(sf, loops[slot[v].loop], slot[v].s, pos[v], tangent[v], bind[v]);
            } else {
                walk(sf, bind[v], d * r.w, out.walkClamps);
                pos[v] = pointOf(reference, bind[v].tri, bind[v].bary);
            }
        };
        if (r.faceCorners == FaceCorners::Free) {
            for (int v = 0; v < nV; ++v) moveVertex(v, 1.0);
        } else {
            // FaceCorners::Kept: the state at the start of the step, the size of every corner then, and the step taken in
            // shares: the three vertices of a clearly convex corner (a quarter of w h or more) that would turn over (end
            // under a twentieth of w h) take half their share, up to six times, then none
            const std::vector<Vec3> pos0 = pos, tangent0 = tangent;
            const std::vector<VertexBinding> bind0 = bind;
            const std::vector<BorderSlot> slot0 = slot;
            const std::vector<double> guideS0 = guideS;
            std::vector<double> before(corners.size());
            std::vector<Vec3> cornerNormal(corners.size());
            auto sizeOf = [&](size_t c) { return dot(cross(pos[corners[c][1]] - pos[corners[c][0]], pos[corners[c][2]] - pos[corners[c][1]]), cornerNormal[c]); };
            for (size_t c = 0; c < corners.size(); ++c) {
                const VertexBinding& b = bind[corners[c][1]];
                cornerNormal[c] = b.tri >= 0 ? sf.normal[b.tri] : Vec3{0, 0, 0};
                before[c] = sizeOf(c);
            }
            std::vector<double> share(nV, 1.0);
            std::vector<char> redo(nV, 1);
            const double judged = 0.25 * r.w * r.h, least = 0.05 * r.w * r.h;
            for (int pass = 0; pass < 8; ++pass) {
                for (int v = 0; v < nV; ++v) {
                    if (!redo[v]) continue;
                    pos[v] = pos0[v];
                    tangent[v] = tangent0[v];
                    bind[v] = bind0[v];
                    slot[v] = slot0[v];
                    guideS[v] = guideS0[v];
                    if (share[v] > 0) moveVertex(v, share[v]);
                    redo[v] = 0;
                }
                bool any = false;
                for (size_t c = 0; c < corners.size(); ++c) {
                    if (before[c] < judged || sizeOf(c) >= least) continue;
                    if (out.firstCornerLimit.empty()) {
                        char text[400];
                        std::snprintf(text, sizeof text, "step %d: corner %d %d %d, size %.3f -> %.3f mm2; filtered forces %.3g %.3g %.3g; solved steps %.3f %.3f %.3f mm",
                                      step, corners[c][0], corners[c][1], corners[c][2], before[c], sizeOf(c), norm(filter(corners[c][0], F[corners[c][0]])),
                                      norm(filter(corners[c][1], F[corners[c][1]])), norm(filter(corners[c][2], F[corners[c][2]])), norm(dx[corners[c][0]]) * r.w,
                                      norm(dx[corners[c][1]]) * r.w, norm(dx[corners[c][2]]) * r.w);
                        out.firstCornerLimit = text;
                    }
                    for (int v : corners[c]) {
                        if (share[v] == 0) continue;
                        share[v] = pass >= 6 ? 0.0 : 0.5 * share[v];
                        redo[v] = 1;
                        any = true;
                        // a ring that slides as one takes the smaller share as one
                        if (r.borderOrder == BorderOrder::RingSlidesAsOne && nextOnLoop[v] >= 0)
                            for (int u = nextOnLoop[v]; u != v; u = nextOnLoop[u]) {
                                share[u] = share[v];
                                redo[u] = 1;
                            }
                    }
                }
                if (!any) break;
            }
            for (int v = 0; v < nV; ++v)
                if (share[v] < 1) ++out.cornerMovesLimited;
        }
        // a border vertex that has passed its neighbour (the order the layout gave is an invariant of a valid mesh)
        for (int v = 0; v < nV; ++v) {
            if (nextOnLoop[v] < 0 || gapAhead[v] <= 0) continue;
            const double after = gapToNext(v);
            if (after > 0) continue;
            if (out.borderOrderViolations++ == 0) {
                const int u = nextOnLoop[v];
                char text[400];
                std::snprintf(text, sizeof text,
                              "step %d: vertices %d and %d on a loop of %.1f mm: gap %.3f -> %.3f mm; filtered force on them %.3g and %.3g, "
                              "solved steps along the loop %.3f and %.3f mm, held %d and %d",
                              step, v, u, loops[slot[v].loop].cum.back(), gapAhead[v], after, norm(filter(v, F[v])), norm(filter(u, F[u])),
                              dot(dx[v], tangent[v]) * r.w, dot(dx[u], tangent[u]) * r.w, int(pinned[v]), int(pinned[u]));
                out.firstBorderOrderViolation = text;
            }
        }
        rec.maxMove = maxMove;
        out.history.push_back(rec);
    }
    out.pos = pos;
    out.binding = bind;
    for (VertexBinding& b : out.binding) b.distance = 0;
    out.ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    out.status = out.stopReason == "converged" ? Status::Pass : out.stopReason == "diverged" ? Status::Fail : Status::Unresolved;
    char b[300];
    std::snprintf(b, sizeof b, "%s after %d steps, RMS force %.4g -> %.4g, %.0f ms", out.stopReason.c_str(), out.steps, out.rmsForceStart, out.rmsForceEnd, out.ms);
    out.detail = b;
    return true;
}

namespace {

// The force's part along the feasible motions on a polyline at p: |F . t| on a segment; at a corner (within tol of a
// polyline vertex) only what leads out along either segment, max(0, F . t_next) and max(0, -F . t_prev).
double lineResidual(const std::vector<Vec3>& pts, bool closed, const Vec3& p, const Vec3& f, double tol) {
    const int n = static_cast<int>(pts.size());
    for (int i = 0; i < n; ++i) {
        if (norm(pts[i] - p) > tol) continue;
        const bool hasPrev = closed || i > 0, hasNext = closed || i + 1 < n;
        double s = 0;
        if (hasNext) {
            const Vec3 t = normalized(pts[(i + 1) % n] - pts[i]);
            s += std::pow(std::max(0.0, dot(f, t)), 2);
        }
        if (hasPrev) {
            const Vec3 t = normalized(pts[i] - pts[(i + n - 1) % n]);
            s += std::pow(std::max(0.0, -dot(f, t)), 2);
        }
        return std::sqrt(s);
    }
    const int segs = closed ? n : n - 1;
    double best = 1e300;
    Vec3 t{1, 0, 0};
    for (int i = 0; i < segs; ++i) {
        const Vec3 a = pts[i], b = pts[(i + 1) % n];
        const double d = closestSegmentPoints(p, p, a, b).distance;
        if (d < best && norm(b - a) > 0) {
            best = d;
            t = normalized(b - a);
        }
    }
    return std::fabs(dot(f, t));
}

// The unit tangent of the closest segment of a polyline (closed or open) to p.
Vec3 lineTangent(const std::vector<Vec3>& pts, bool closed, const Vec3& p) {
    const int n = static_cast<int>(pts.size());
    const int segs = closed ? n : n - 1;
    double best = 1e300;
    Vec3 t{1, 0, 0};
    for (int i = 0; i < segs; ++i) {
        const Vec3 a = pts[i], b = pts[(i + 1) % n];
        const double d = closestSegmentPoints(p, p, a, b).distance;
        if (d < best && norm(b - a) > 0) {
            best = d;
            t = normalized(b - a);
        }
    }
    return t;
}

} // namespace

WalkResidual walkResidualOf(const TriMesh& reference, const std::vector<Vec3>& pos, const std::vector<Vec3>& F, const std::vector<char>& border,
                            const std::vector<std::vector<int>>& refLoops, const RelaxGuides* guides, double edgeTolerance) {
    WalkResidual w;
    const int nV = static_cast<int>(pos.size());
    w.perVertex.assign(nV, 0);
    const MeshQuery q(reference);
    std::map<EdgeKey, std::vector<int>> edgeTris;
    std::vector<std::vector<int>> trisOf(reference.pos.size());
    for (int t = 0; t < static_cast<int>(reference.tri.size()); ++t)
        for (int k = 0; k < 3; ++k) {
            edgeTris[keyOf(reference.tri[t][k], reference.tri[t][(k + 1) % 3])].push_back(t);
            trisOf[reference.tri[t][k]].push_back(t);
        }
    std::vector<std::vector<Vec3>> loopPts;
    for (const std::vector<int>& l : refLoops) {
        std::vector<Vec3> p;
        for (int v : l) p.push_back(reference.pos[v]);
        loopPts.push_back(p);
    }
    double sum = 0;
    for (int v = 0; v < nV; ++v) {
        const Vec3& f = F[v];
        double r = 0;
        if (guides && v < static_cast<int>(guides->pinned.size()) && guides->pinned[v]) {
            ++w.constrained;
        } else if (guides && v < static_cast<int>(guides->lineOf.size()) && guides->lineOf[v] >= 0) {
            const GuideLine& g = guides->lines[guides->lineOf[v]];
            r = lineResidual(g.points, g.closed, pos[v], f, edgeTolerance);
            ++w.constrained;
        } else if (v < static_cast<int>(border.size()) && border[v] && !loopPts.empty()) {
            double best = 1e300;
            for (const std::vector<Vec3>& l : loopPts) {
                double d = 1e300;
                for (size_t i = 0; i < l.size(); ++i) d = std::min(d, closestSegmentPoints(pos[v], pos[v], l[i], l[(i + 1) % l.size()]).distance);
                if (d < best) {
                    best = d;
                    r = lineResidual(l, true, pos[v], f, edgeTolerance);
                }
            }
            ++w.constrained;
        } else {
            const MeshQuery::Hit h = q.closest(pos[v]);
            const auto& tr = reference.tri[h.tri];
            // edges within the tolerance
            std::vector<int> near;
            for (int k = 0; k < 3; ++k) {
                const Vec3 a = reference.pos[tr[k]], b = reference.pos[tr[(k + 1) % 3]];
                if (closestSegmentPoints(h.p, h.p, a, b).distance <= edgeTolerance) near.push_back(k);
            }
            if (near.empty()) {
                const Vec3 n = reference.faceNormal(h.tri);
                r = norm(f - n * dot(n, f));
                ++w.interior;
            } else if (near.size() == 1) {
                const int a = tr[near[0]], b = tr[(near[0] + 1) % 3];
                const Vec3 l = normalized(reference.pos[b] - reference.pos[a]);
                double s = dot(f, l) * dot(f, l);
                for (int t : edgeTris[keyOf(a, b)]) {
                    int c = -1;
                    for (int k = 0; k < 3; ++k)
                        if (reference.tri[t][k] != a && reference.tri[t][k] != b) c = reference.tri[t][k];
                    Vec3 d = reference.pos[c] - reference.pos[a];
                    d = normalized(d - l * dot(l, d));   // into facet t, across the edge
                    const double pull = std::max(0.0, dot(f, d));
                    s += pull * pull;
                }
                r = std::sqrt(s);
                ++w.onEdge;
            } else {
                // at a mesh vertex: the corner shared by the near edges
                const int k0 = near[0], k1 = near[1];
                const int corner = (k0 + 1) % 3 == k1 ? tr[k1] : tr[k0];
                for (int t : trisOf[corner]) {
                    int i = 0;
                    while (reference.tri[t][i] != corner) ++i;
                    const Vec3 e1 = normalized(reference.pos[reference.tri[t][(i + 1) % 3]] - reference.pos[corner]);
                    const Vec3 e2 = normalized(reference.pos[reference.tri[t][(i + 2) % 3]] - reference.pos[corner]);
                    for (int sIdx = 0; sIdx <= 12; ++sIdx) {
                        const Vec3 d = normalized(e1 * (1.0 - sIdx / 12.0) + e2 * (sIdx / 12.0));
                        r = std::max(r, dot(f, d));
                    }
                }
                ++w.onVertex;
            }
        }
        w.perVertex[v] = r;
        sum += r * r;
        if (r > w.max) {
            w.max = r;
            w.maxVertex = v;
        }
    }
    w.rms = std::sqrt(sum / std::max(1, nV));
    return w;
}

WalkResidual walkStationarity(const TriMesh& reference, const std::vector<RelaxFace>& faces, const std::vector<Vec3>& pos, const std::vector<char>& border,
                              const std::vector<std::vector<int>>& refLoops, const MeshRelaxRequest& r, const RelaxGuides* guides, double edgeTolerance) {
    std::vector<Vec3> x(pos.size());
    for (size_t v = 0; v < pos.size(); ++v) x[v] = pos[v] / r.w;
    const std::vector<Vec3> F = relaxForces(faces, x, r);
    return walkResidualOf(reference, pos, F, border, refLoops, guides, edgeTolerance);
}

} // namespace knit
