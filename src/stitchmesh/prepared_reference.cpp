#include "stitchmesh/prepared_reference.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <utility>

namespace knit {

namespace {

constexpr double kPi = 3.14159265358979323846;

// Loop's weight of each neighbour of an interior vertex of valence n (Loop 1987): beta = (5/8 - (3/8 + cos(2 pi / n) / 4)^2) / n.
double loopBeta(int n) {
    const double c = 0.375 + 0.25 * std::cos(2 * kPi / n);
    return (0.625 - c * c) / n;
}

// The neighbourhood of every vertex of a triangle mesh, as Loop's rules need it.
struct LoopNeighbourhood {
    std::vector<std::set<int>> ring;          // all neighbours
    std::vector<std::vector<int>> borderNb;   // neighbours along border edges (edges with one triangle)
    std::map<std::pair<int, int>, std::vector<int>> opposite;   // per edge: the third corner of each triangle on it
};

LoopNeighbourhood neighbourhood(const TriMesh& m) {
    LoopNeighbourhood n;
    n.ring.resize(m.pos.size());
    n.borderNb.resize(m.pos.size());
    for (const auto& t : m.tri)
        for (int k = 0; k < 3; ++k) {
            const int a = t[k], b = t[(k + 1) % 3], c = t[(k + 2) % 3];
            n.ring[a].insert(b);
            n.ring[b].insert(a);
            n.opposite[{std::min(a, b), std::max(a, b)}].push_back(c);
        }
    for (const auto& [e, opp] : n.opposite)
        if (opp.size() == 1) {
            n.borderNb[e.first].push_back(e.second);
            n.borderNb[e.second].push_back(e.first);
        }
    return n;
}

// Loop's limit position of every vertex of m (the point of the limit surface the vertex converges to): interior
// (1 - n chi) v + chi sum(neighbours), chi = 1 / (n + 3 / (8 beta)); on a smooth border 2/3 v + 1/6 (b0 + b1); a corner (a
// border vertex without exactly two border neighbours) stays.
std::vector<Vec3> loopLimit(const TriMesh& m, const LoopNeighbourhood& nb) {
    std::vector<Vec3> out(m.pos.size());
    for (size_t v = 0; v < m.pos.size(); ++v) {
        const std::vector<int>& bn = nb.borderNb[v];
        if (!bn.empty()) {
            out[v] = bn.size() == 2 ? (2.0 / 3.0) * m.pos[v] + (1.0 / 6.0) * (m.pos[bn[0]] + m.pos[bn[1]]) : m.pos[v];
            continue;
        }
        const int n = static_cast<int>(nb.ring[v].size());
        if (n < 3) {
            out[v] = m.pos[v];
            continue;
        }
        const double chi = 1.0 / (n + 3.0 / (8.0 * loopBeta(n)));
        Vec3 sum{0, 0, 0};
        for (int u : nb.ring[v]) sum = sum + m.pos[u];
        out[v] = (1 - n * chi) * m.pos[v] + chi * sum;
    }
    return out;
}

// One level of the 1-to-4 split (every triangle into four at its edge midpoints, in a fixed order, so the connectivity,
// parent and corner coordinates are the same under either rule). loop false: the new vertices at the midpoints, the old
// ones stay (the surface does not move). loop true: Loop's rules (Loop 1987), with the cubic B-spline rule on borders:
//   old vertex   interior (1 - n beta) v + beta sum(neighbours); smooth border 3/4 v + 1/8 (b0 + b1); corner stays
//   edge vertex  interior 3/8 (a + b) + 1/8 (c + d), c and d the opposite corners; border (or non-manifold) (a + b) / 2
RefinedMesh splitOnce(const RefinedMesh& r, bool loop) {
    RefinedMesh next;
    next.levels = r.levels;
    LoopNeighbourhood nb;
    if (loop) nb = neighbourhood(r.mesh);
    next.mesh.pos = r.mesh.pos;
    if (loop)
        for (size_t v = 0; v < r.mesh.pos.size(); ++v) {
            const std::vector<int>& bn = nb.borderNb[v];
            if (!bn.empty()) {
                if (bn.size() == 2) next.mesh.pos[v] = 0.75 * r.mesh.pos[v] + 0.125 * (r.mesh.pos[bn[0]] + r.mesh.pos[bn[1]]);
                continue;
            }
            const int n = static_cast<int>(nb.ring[v].size());
            if (n < 3) continue;
            const double beta = loopBeta(n);
            Vec3 sum{0, 0, 0};
            for (int u : nb.ring[v]) sum = sum + r.mesh.pos[u];
            next.mesh.pos[v] = (1 - n * beta) * r.mesh.pos[v] + beta * sum;
        }
    std::map<std::pair<int, int>, int> mid;
    auto midpoint = [&](int a, int b) {
        const std::pair<int, int> k{std::min(a, b), std::max(a, b)};
        const auto it = mid.find(k);
        if (it != mid.end()) return it->second;
        const int v = static_cast<int>(next.mesh.pos.size());
        Vec3 p = 0.5 * (r.mesh.pos[a] + r.mesh.pos[b]);
        if (loop) {
            const std::vector<int>& opp = nb.opposite.at(k);
            if (opp.size() == 2) p = 0.375 * (r.mesh.pos[a] + r.mesh.pos[b]) + 0.125 * (r.mesh.pos[opp[0]] + r.mesh.pos[opp[1]]);
        }
        next.mesh.pos.push_back(p);
        mid[k] = v;
        return v;
    };
    for (size_t t = 0; t < r.mesh.tri.size(); ++t) {
        const auto& tr = r.mesh.tri[t];
        const int m01 = midpoint(tr[0], tr[1]), m12 = midpoint(tr[1], tr[2]), m20 = midpoint(tr[2], tr[0]);
        const auto& cb = r.cornerBary[t];
        std::array<double, 3> b01{}, b12{}, b20{};
        for (int k = 0; k < 3; ++k) {
            b01[k] = 0.5 * (cb[0][k] + cb[1][k]);
            b12[k] = 0.5 * (cb[1][k] + cb[2][k]);
            b20[k] = 0.5 * (cb[2][k] + cb[0][k]);
        }
        const std::array<std::array<int, 3>, 4> kids{{{tr[0], m01, m20}, {m01, tr[1], m12}, {m20, m12, tr[2]}, {m01, m12, m20}}};
        const std::array<std::array<std::array<double, 3>, 3>, 4> bary{{{cb[0], b01, b20}, {b01, cb[1], b12}, {b20, b12, cb[2]}, {b01, b12, b20}}};
        for (int k = 0; k < 4; ++k) {
            next.mesh.tri.push_back(kids[k]);
            next.parent.push_back(r.parent[t]);
            next.cornerBary.push_back(bary[k]);
        }
    }
    return next;
}

RefinedMesh unrefined(const TriMesh& in, int levels) {
    RefinedMesh r;
    r.levels = levels;
    r.mesh.pos = in.pos;
    r.mesh.tri = in.tri;
    for (size_t t = 0; t < in.tri.size(); ++t) {
        r.parent.push_back(static_cast<int>(t));
        r.cornerBary.push_back({std::array<double, 3>{1, 0, 0}, std::array<double, 3>{0, 1, 0}, std::array<double, 3>{0, 0, 1}});
    }
    return r;
}

} // namespace

const char* toString(ReferencePreparation p) { return p == ReferencePreparation::LoopLimit ? "loop_limit" : "facets"; }

bool parsePreparation(const std::string& text, ReferencePreparation& out) {
    if (text == "facets") out = ReferencePreparation::Facets;
    else if (text == "loop_limit") out = ReferencePreparation::LoopLimit;
    else return false;
    return true;
}

RefinedMesh refineMidpoint(const TriMesh& in, int levels) {
    RefinedMesh r = unrefined(in, levels);
    for (int level = 0; level < levels; ++level) r = splitOnce(r, false);
    r.mesh.computeNormals();
    return r;
}

RefinedMesh refineLoop(const TriMesh& in, int levels) {
    RefinedMesh r = unrefined(in, levels);
    for (int level = 0; level < levels; ++level) r = splitOnce(r, true);
    r.mesh.pos = loopLimit(r.mesh, neighbourhood(r.mesh));
    r.mesh.computeNormals();
    return r;
}

RefinedMesh prepareReference(const TriMesh& in, int levels, ReferencePreparation p) {
    return p == ReferencePreparation::LoopLimit ? refineLoop(in, levels) : refineMidpoint(in, levels);
}

SurfacePoint mapBetween(const RefinedMesh& from, const RefinedMesh& to, const SurfacePoint& p) {
    if (p.tri < 0) return SurfacePoint{};
    if (from.mesh.tri.size() == to.mesh.tri.size() && from.parent == to.parent) return p;   // the same split: the same triangle
    return toRefined(to, fromRefined(from, p));
}

int refinementLevels(const TriMesh& in, double maxEdge) {
    double sum = 0;
    int n = 0;
    for (const auto& t : in.tri)
        for (int k = 0; k < 3; ++k) {
            sum += norm(in.pos[t[(k + 1) % 3]] - in.pos[t[k]]);
            ++n;
        }
    double mean = n ? sum / n : 0;
    int levels = 0;
    while (mean > maxEdge && levels < 4) {
        mean *= 0.5;
        ++levels;
    }
    return levels;
}

namespace {

// Solves for the weights w of p = sum w_k c_k with c the corners' barycentric coordinates in the parent (2 unknowns).
bool weightsIn(const std::array<std::array<double, 3>, 3>& c, const double p[3], double w[3]) {
    // p - c2 = w0 (c0 - c2) + w1 (c1 - c2), least squares over the 3 coordinates
    double a00 = 0, a01 = 0, a11 = 0, r0 = 0, r1 = 0;
    for (int k = 0; k < 3; ++k) {
        const double u = c[0][k] - c[2][k], v = c[1][k] - c[2][k], d = p[k] - c[2][k];
        a00 += u * u;
        a01 += u * v;
        a11 += v * v;
        r0 += u * d;
        r1 += v * d;
    }
    const double det = a00 * a11 - a01 * a01;
    if (!(std::fabs(det) > 0)) return false;
    w[0] = (r0 * a11 - r1 * a01) / det;
    w[1] = (a00 * r1 - a01 * r0) / det;
    w[2] = 1 - w[0] - w[1];
    return true;
}

} // namespace

SurfacePoint toRefined(const RefinedMesh& r, const SurfacePoint& input) {
    SurfacePoint best;
    double bestNeg = -1e300;
    for (size_t t = 0; t < r.parent.size(); ++t) {
        if (r.parent[t] != input.tri) continue;
        double w[3];
        if (!weightsIn(r.cornerBary[t], input.bary, w)) continue;
        const double worst = std::min(w[0], std::min(w[1], w[2]));
        if (worst > bestNeg) {
            bestNeg = worst;
            best.tri = static_cast<int>(t);
            for (int k = 0; k < 3; ++k) best.bary[k] = w[k];
        }
        if (worst >= -1e-12) break;
    }
    if (best.tri >= 0) {   // clamp the rounding and renormalize
        double sum = 0;
        for (double& x : best.bary) sum += (x = std::max(0.0, x));
        for (double& x : best.bary) x /= sum;
    }
    return best;
}

SurfacePoint fromRefined(const RefinedMesh& r, const SurfacePoint& refined) {
    SurfacePoint p;
    p.tri = r.parent[refined.tri];
    for (int k = 0; k < 3; ++k) {
        p.bary[k] = 0;
        for (int c = 0; c < 3; ++c) p.bary[k] += refined.bary[c] * r.cornerBary[refined.tri][c][k];
    }
    return p;
}

TriMesh clipBelowLevel(const TriMesh& m, const std::vector<double>& t, double level, std::vector<int>& parent) {
    TriMesh out;
    out.pos = m.pos;
    std::map<std::pair<int, int>, int> cut;
    auto edgePoint = [&](int a, int b) {
        const std::pair<int, int> k{std::min(a, b), std::max(a, b)};
        const auto it = cut.find(k);
        if (it != cut.end()) return it->second;
        const double f = (level - t[a]) / (t[b] - t[a]);
        const int v = static_cast<int>(out.pos.size());
        out.pos.push_back((1 - f) * m.pos[a] + f * m.pos[b]);
        cut[k] = v;
        return v;
    };
    parent.clear();
    for (int tt = 0; tt < static_cast<int>(m.tri.size()); ++tt) {
        const auto& tr = m.tri[tt];
        int below = 0;
        for (int k = 0; k < 3; ++k) below += t[tr[k]] < level ? 1 : 0;
        if (below == 0) continue;
        if (below == 3) {
            out.tri.push_back(tr);
            parent.push_back(tt);
            continue;
        }
        int r0 = 0;
        for (int k = 0; k < 3; ++k) {
            const bool b = t[tr[k]] < level;
            if ((below == 1 && b) || (below == 2 && !b)) r0 = k;
        }
        const int a = tr[r0], b = tr[(r0 + 1) % 3], c = tr[(r0 + 2) % 3];
        const int ab = edgePoint(a, b), ac = edgePoint(a, c);
        if (below == 1) {
            out.tri.push_back({a, ab, ac});
            parent.push_back(tt);
        } else {
            out.tri.push_back({ab, b, c});
            out.tri.push_back({ab, c, ac});
            parent.push_back(tt);
            parent.push_back(tt);
        }
    }
    std::vector<int> map(out.pos.size(), -1);
    TriMesh compact;
    for (auto& tr : out.tri)
        for (int& v : tr) {
            if (map[v] < 0) {
                map[v] = static_cast<int>(compact.pos.size());
                compact.pos.push_back(out.pos[v]);
            }
            v = map[v];
        }
    compact.tri = out.tri;
    compact.computeNormals();
    return compact;
}

std::vector<std::vector<int>> borderLoops(const TriMesh& m) {
    std::map<std::pair<int, int>, int> use;
    for (const auto& t : m.tri)
        for (int k = 0; k < 3; ++k) use[{std::min(t[k], t[(k + 1) % 3]), std::max(t[k], t[(k + 1) % 3])}] += 1;
    std::map<int, int> next;
    for (const auto& t : m.tri)
        for (int k = 0; k < 3; ++k) {
            const int x = t[k], y = t[(k + 1) % 3];
            if (use[{std::min(x, y), std::max(x, y)}] == 1) next[x] = y;
        }
    std::vector<std::vector<int>> loops;
    std::set<int> seen;
    for (const auto& [s, unused] : next) {
        (void) unused;
        if (seen.count(s)) continue;
        std::vector<int> loop;
        for (int v = s; !seen.count(v) && next.count(v); v = next[v]) {
            seen.insert(v);
            loop.push_back(v);
        }
        loops.push_back(loop);
    }
    return loops;
}

} // namespace knit
