#include "stitchmesh/labels.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <set>

#include "Highs.h"

namespace knit {

namespace {

using EdgeKey = std::pair<int, int>;
EdgeKey edgeKey(int a, int b) { return a < b ? EdgeKey(a, b) : EdgeKey(b, a); }

bool quadValid(int l0, int l1, int l2, int l3) { return l0 == l2 && l1 == l3 && l0 != l1; }
bool triangleValid(int courses) { return courses >= 1 && courses <= 2; }

// Faces as vertex lists with their directed edges, killed in place and added at the end: the shape the operations of
// sections 5.1 to 5.3 edit. live() gives the surviving faces in order.
struct FaceList {
    std::vector<std::vector<int>> faces;
    std::vector<char> alive;
    std::map<EdgeKey, int> faceOfDirected;  // u -> v: the face on its left

    int add(const std::vector<int>& f) {
        const int id = static_cast<int>(faces.size());
        faces.push_back(f);
        alive.push_back(1);
        for (size_t k = 0; k < f.size(); ++k) faceOfDirected[{f[k], f[(k + 1) % f.size()]}] = id;
        return id;
    }
    void kill(int id) {
        const std::vector<int>& f = faces[id];
        for (size_t k = 0; k < f.size(); ++k) {
            const auto it = faceOfDirected.find({f[k], f[(k + 1) % f.size()]});
            if (it != faceOfDirected.end() && it->second == id) faceOfDirected.erase(it);
        }
        alive[id] = 0;
    }
    int across(int u, int v) const {  // the face on the other side of the directed edge u -> v, or -1
        const auto it = faceOfDirected.find({v, u});
        return it == faceOfDirected.end() ? -1 : it->second;
    }
    bool isEdge(int a, int b) const { return faceOfDirected.count({a, b}) || faceOfDirected.count({b, a}); }
    std::vector<std::vector<int>> live() const {
        std::vector<std::vector<int>> out;
        for (size_t i = 0; i < faces.size(); ++i)
            if (alive[i]) out.push_back(faces[i]);
        return out;
    }
};

// f rotated to start at vertex u.
std::vector<int> startAt(const std::vector<int>& f, int u) {
    const auto it = std::find(f.begin(), f.end(), u);
    std::vector<int> r(it, f.end());
    r.insert(r.end(), f.begin(), it);
    return r;
}

// The quad of triangles a (containing u -> v) and b (containing v -> u): [v, a's third vertex, u, b's third vertex].
std::vector<int> mergedQuad(const std::vector<int>& a, const std::vector<int>& b, int u, int v) {
    const std::vector<int> ra = startAt(a, u), rb = startAt(b, v);
    return {v, ra[2], u, rb[2]};
}

int courseCount(const std::map<EdgeKey, char>& roles, const std::vector<int>& f) {
    int n = 0;
    for (size_t k = 0; k < f.size(); ++k) n += roles.at(edgeKey(f[k], f[(k + 1) % f.size()])) == 1;
    return n;
}

bool faceValid(const std::map<EdgeKey, char>& roles, const std::vector<int>& f) {
    if (f.size() == 4) {
        int l[4];
        for (int k = 0; k < 4; ++k) l[k] = roles.at(edgeKey(f[k], f[(k + 1) % 4]));
        return quadValid(l[0], l[1], l[2], l[3]);
    }
    if (f.size() == 3) return triangleValid(courseCount(roles, f));
    return false;
}

struct Topology {
    std::set<int> vertices;
    std::set<EdgeKey> boundary;
    long euler = 0;
    int overusedEdges = 0, degenerateFaces = 0;
};

Topology topologyOf(const std::vector<std::vector<int>>& faces) {
    Topology t;
    std::map<EdgeKey, int> count;
    for (const std::vector<int>& f : faces) {
        t.degenerateFaces += std::set<int>(f.begin(), f.end()).size() != f.size();
        for (size_t k = 0; k < f.size(); ++k) {
            t.vertices.insert(f[k]);
            ++count[edgeKey(f[k], f[(k + 1) % f.size()])];
        }
    }
    for (const auto& [e, n] : count) {
        if (n == 1) t.boundary.insert(e);
        t.overusedEdges += n > 2;
    }
    t.euler = static_cast<long>(t.vertices.size()) - static_cast<long>(count.size()) + static_cast<long>(faces.size());
    return t;
}

std::vector<std::vector<int>> facesOf(const PolygonMesh& m) {
    std::vector<std::vector<int>> faces(m.faces());
    for (int f = 0; f < m.faces(); ++f)
        for (int k = 0; k < m.sides(f); ++k) faces[f].push_back(m.vertexAt(f, k));
    return faces;
}

// --- 5.1: the quads around interior vertices of odd valence are triangulated through that vertex

std::vector<std::vector<int>> triangulateOddVertices(const std::vector<std::vector<int>>& faces, const Topology& topology, LabelProduct& p) {
    std::map<int, std::set<int>> neighbours;
    for (const std::vector<int>& f : faces)
        for (size_t k = 0; k < f.size(); ++k) {
            neighbours[f[k]].insert(f[(k + 1) % f.size()]);
            neighbours[f[(k + 1) % f.size()]].insert(f[k]);
        }
    std::set<int> boundaryVertices;
    for (const EdgeKey& e : topology.boundary) {
        boundaryVertices.insert(e.first);
        boundaryVertices.insert(e.second);
    }
    std::set<int> odd;
    for (const auto& [v, around] : neighbours)
        if (!boundaryVertices.count(v) && around.size() % 2 == 1) odd.insert(v);
    p.oddVertices = static_cast<int>(odd.size());
    std::set<EdgeKey> edges;  // a diagonal that is an edge elsewhere already would make the mesh non-manifold
    for (const std::vector<int>& f : faces)
        for (size_t k = 0; k < f.size(); ++k) edges.insert(edgeKey(f[k], f[(k + 1) % f.size()]));
    std::vector<std::vector<int>> out;
    for (const std::vector<int>& f : faces) {
        int at = -1;
        if (f.size() == 4)
            for (int v : odd)  // the smallest odd vertex of the quad: deterministic
                if (std::find(f.begin(), f.end(), v) != f.end()) {
                    at = v;
                    break;
                }
        std::vector<int> r = at >= 0 ? startAt(f, at) : f;  // the diagonal r0 - r2
        if (at >= 0 && edges.count(edgeKey(r[0], r[2]))) {
            r = startAt(f, r[1]);  // the other diagonal
            if (edges.count(edgeKey(r[0], r[2]))) {
                at = -1;
                ++p.untriangulatedQuads;
            }
        }
        if (at < 0) {
            out.push_back(f);
            continue;
        }
        edges.insert(edgeKey(r[0], r[2]));
        out.push_back({r[0], r[1], r[2]});
        out.push_back({r[2], r[3], r[0]});
        ++p.triangulatedQuads;
    }
    return out;
}

// --- the half-edge program

struct HalfEdges {
    std::vector<std::vector<int>> ofFace;     // face -> its half-edge ids, in corner order
    std::vector<int> edgeOfHalf;              // half-edge -> undirected edge id
    std::vector<EdgeKey> edgeKeys;            // undirected edge id -> (min, max)
    std::vector<std::pair<int, int>> interior;  // undirected edges with two half-edges: their ids
    int count = 0;
};

HalfEdges halfEdgesOf(const std::vector<std::vector<int>>& faces) {
    HalfEdges h;
    h.ofFace.resize(faces.size());
    std::map<EdgeKey, int> idOf;
    std::vector<std::vector<int>> halvesOfEdge;
    for (size_t f = 0; f < faces.size(); ++f)
        for (size_t k = 0; k < faces[f].size(); ++k) {
            const EdgeKey e = edgeKey(faces[f][k], faces[f][(k + 1) % faces[f].size()]);
            auto it = idOf.find(e);
            if (it == idOf.end()) {
                it = idOf.emplace(e, static_cast<int>(h.edgeKeys.size())).first;
                h.edgeKeys.push_back(e);
                halvesOfEdge.emplace_back();
            }
            h.ofFace[f].push_back(h.count);
            h.edgeOfHalf.push_back(it->second);
            halvesOfEdge[it->second].push_back(h.count);
            ++h.count;
        }
    for (const std::vector<int>& halves : halvesOfEdge)
        if (halves.size() == 2) h.interior.push_back({halves[0], halves[1]});
    return h;
}

// Union-find with parity over the edges: opposite edges of a quad equal, adjacent ones different.
struct ParityClasses {
    std::vector<int> parent, relation;
    int conflicts = 0;
    explicit ParityClasses(int n) : parent(n), relation(n, 0) {
        for (int e = 0; e < n; ++e) parent[e] = e;
    }
    int find(int e, int& parity) {
        parity = 0;
        int root = e;
        while (parent[root] != root) {
            parity ^= relation[root];
            root = parent[root];
        }
        int q = e, pq = parity;  // path compression, every node straight to the root
        while (parent[q] != q) {
            const int next = parent[q], rq = relation[q];
            parent[q] = root;
            relation[q] = pq;
            pq ^= rq;
            q = next;
        }
        return root;
    }
    void unite(int a, int b, int relate) {
        int pa = 0, pb = 0;
        int ra = find(a, pa), rb = find(b, pb);
        if (ra == rb) {
            conflicts += (pa ^ pb) != relate;
            return;
        }
        if (rb < ra) std::swap(ra, rb);
        parent[rb] = ra;
        relation[rb] = pa ^ pb ^ relate;
    }
};

// A feasible start for HiGHS: edge roles from the quads' parity classes, the classes' values by a deterministic
// local search on the triangles (no triangle with three equal roles), then per face the nearest valid half-edge
// roles. The program's LP bound is 0, so a zero-conflict start is all the work where one exists; the search is
// seeded, so the start and the solution are the same on every run.
std::vector<char> heuristicStart(const std::vector<std::vector<int>>& faces, const HalfEdges& h, const std::vector<char>& preferCourse,
                                 int& parityConflicts, int& startConflicts) {
    const int nEdges = static_cast<int>(h.edgeKeys.size());
    ParityClasses classes(nEdges);
    std::vector<std::array<int, 3>> triangles;
    for (size_t f = 0; f < faces.size(); ++f) {
        int e[4];
        for (size_t k = 0; k < faces[f].size(); ++k) e[k] = h.edgeOfHalf[h.ofFace[f][k]];
        if (faces[f].size() == 4) {
            classes.unite(e[0], e[2], 0);
            classes.unite(e[1], e[3], 0);
            classes.unite(e[0], e[1], 1);
        } else {
            triangles.push_back({e[0], e[1], e[2]});
        }
    }
    parityConflicts = classes.conflicts;
    std::vector<char> value(nEdges, 0);
    {
        std::vector<int> agree0(nEdges, 0), agree1(nEdges, 0);  // per root: edges that get their preferred role at value 0 / 1
        for (int e = 0; e < nEdges; ++e) {
            int parity = 0;
            const int root = classes.find(e, parity);
            (preferCourse[e] == parity ? agree0 : agree1)[root]++;
        }
        for (int r = 0; r < nEdges; ++r) value[r] = agree1[r] > agree0[r] ? 1 : 0;
    }
    auto roleOf = [&](int e) {
        int parity = 0;
        const int root = classes.find(e, parity);
        return value[root] ^ parity;
    };
    auto violated = [&](const std::array<int, 3>& t) { return roleOf(t[0]) == roleOf(t[1]) && roleOf(t[1]) == roleOf(t[2]); };
    std::map<int, std::vector<int>> trianglesOfRoot;
    for (size_t t = 0; t < triangles.size(); ++t) {
        std::set<int> roots;
        for (int e : triangles[t]) {
            int parity = 0;
            roots.insert(classes.find(e, parity));
        }
        for (int r : roots) trianglesOfRoot[r].push_back(static_cast<int>(t));
    }
    uint64_t rng = 0x9e3779b97f4a7c15ull;
    auto next = [&]() { return (rng = rng * 6364136223846793005ull + 1442695040888963407ull) >> 33; };
    const size_t maxFlips = 50 * triangles.size() + 100;
    for (size_t it = 0; it < maxFlips; ++it) {
        int bad = -1, nBad = 0;
        for (size_t t = 0; t < triangles.size(); ++t)
            if (violated(triangles[t])) {
                ++nBad;
                if (bad < 0 || next() % static_cast<uint64_t>(nBad) == 0) bad = static_cast<int>(t);
            }
        if (bad < 0) break;
        int roots[3];
        for (int k = 0; k < 3; ++k) {
            int parity = 0;
            roots[k] = classes.find(triangles[bad][k], parity);
        }
        int pick = roots[0];
        if (next() % 5 == 0) {
            pick = roots[next() % 3];  // noise
        } else {
            int best = 1 << 30;
            for (int r : roots) {
                int before = 0, after = 0;
                for (int t : trianglesOfRoot[r]) before += violated(triangles[t]);
                value[r] ^= 1;
                for (int t : trianglesOfRoot[r]) after += violated(triangles[t]);
                value[r] ^= 1;
                if (after - before < best || (after - before == best && r < pick)) {
                    best = after - before;
                    pick = r;
                }
            }
        }
        value[pick] ^= 1;
    }
    std::vector<char> x;
    for (size_t f = 0; f < faces.size(); ++f) {
        const size_t n = faces[f].size();
        int l[4];
        for (size_t k = 0; k < n; ++k) l[k] = roleOf(h.edgeOfHalf[h.ofFace[f][k]]);
        if (n == 4 && !quadValid(l[0], l[1], l[2], l[3])) {  // the alternation matching most edges
            const int m0 = (l[0] == 0) + (l[1] == 1) + (l[2] == 0) + (l[3] == 1);
            const int v = m0 >= 2 ? 0 : 1;
            l[0] = l[2] = v;
            l[1] = l[3] = 1 - v;
        } else if (n == 3 && !triangleValid(l[0] + l[1] + l[2])) {
            l[0] ^= 1;
        }
        for (size_t k = 0; k < n; ++k) x.push_back(static_cast<char>(l[k]));
    }
    startConflicts = 0;
    for (const std::pair<int, int>& e : h.interior) startConflicts += x[e.first] != x[e.second];
    return x;
}

// One binary per half-edge (1 course, 0 wale) and one continuous d_e in [0, 1] per interior edge, d_e >= +-(x_a - x_b),
// sum of d_e minimized: the number of edges whose two half-edges disagree.
bool solveHalfEdgeRoles(const std::vector<std::vector<int>>& faces, const HalfEdges& h, const std::vector<Vec3>& pos, const LabelOptions& o,
                        LabelProduct& p, std::vector<char>& x, std::string* err) {
    const int nH = h.count, nE = static_cast<int>(h.interior.size()), nCol = nH + nE;
    HighsLp lp;
    lp.num_col_ = nCol;
    lp.sense_ = ObjSense::kMinimize;
    int triangles = 0;
    for (const std::vector<int>& f : faces) triangles += f.size() == 3;
    // one disagreeing edge costs more than every triangle's course half-edges together: the preference never buys a disagreement
    const double disagreementCost = o.preferShaping && triangles > 0 ? 2.0 * triangles + 1.0 : 1.0;
    lp.col_cost_.assign(nCol, 0.0);
    for (int e = 0; e < nE; ++e) lp.col_cost_[nH + e] = disagreementCost;
    if (o.preferShaping)
        for (size_t f = 0; f < faces.size(); ++f)
            if (faces[f].size() == 3)
                for (int hf : h.ofFace[f]) lp.col_cost_[hf] = 1.0;
    lp.col_lower_.assign(nCol, 0.0);
    lp.col_upper_.assign(nCol, 1.0);
    lp.integrality_.assign(nCol, HighsVarType::kContinuous);
    for (int i = 0; i < nH; ++i) lp.integrality_[i] = HighsVarType::kInteger;
    std::vector<HighsInt> start{0}, index;
    std::vector<double> value;
    auto addRow = [&](std::initializer_list<std::pair<int, double>> terms, double lo, double hi) {
        for (const std::pair<int, double>& t : terms) {
            index.push_back(t.first);
            value.push_back(t.second);
        }
        start.push_back(static_cast<HighsInt>(index.size()));
        lp.row_lower_.push_back(lo);
        lp.row_upper_.push_back(hi);
    };
    for (size_t f = 0; f < faces.size(); ++f) {
        const std::vector<int>& hf = h.ofFace[f];
        if (faces[f].size() == 4) {
            addRow({{hf[0], 1}, {hf[2], -1}}, 0, 0);
            addRow({{hf[1], 1}, {hf[3], -1}}, 0, 0);
            addRow({{hf[0], 1}, {hf[1], 1}}, 1, 1);
        } else {
            addRow({{hf[0], 1}, {hf[1], 1}, {hf[2], 1}}, 1, 2);
        }
    }
    for (int e = 0; e < nE; ++e) {
        const std::pair<int, int> ab = h.interior[e];
        addRow({{nH + e, 1}, {ab.first, -1}, {ab.second, 1}}, 0, kHighsInf);
        addRow({{nH + e, 1}, {ab.first, 1}, {ab.second, -1}}, 0, kHighsInf);
    }
    // Section 5.3's last rule inside the program: two triangles with two wale edges each may share a wale edge only if
    // removing it leaves a valid quad (the post-processing merges them); otherwise the paper flips the shared edge to
    // course, which turns both into row ends. With the shaping preference the program otherwise chooses exactly those
    // fans at odd vertices (sphere, 2026-11-13: 40 row ends after the flips). Triangle t1 = (u, w, x) with the shared
    // half-edge u -> w then p: w -> x, q: x -> u; t2 = (w, u, y) with r: u -> y, s: y -> w; the merged quad runs p, q, r, s,
    // valid iff p = r and q = s differ. Forbidden: shared wale, one course edge each, p = s: (p, q, r, s) = (1, 0, 0, 1)
    // or (0, 1, 1, 0). A shared course edge (either half) loosens both rows by one.
    p.pairRules = 0;
    if (o.preferShaping) {
        std::vector<int> faceOfHalf(nH, -1), cornerOfHalf(nH, -1);
        for (size_t f = 0; f < faces.size(); ++f)
            for (size_t k = 0; k < h.ofFace[f].size(); ++k) {
                faceOfHalf[h.ofFace[f][k]] = static_cast<int>(f);
                cornerOfHalf[h.ofFace[f][k]] = static_cast<int>(k);
            }
        for (int e = 0; e < nE; ++e) {
            const std::pair<int, int> ab = h.interior[e];
            const int t1 = faceOfHalf[ab.first], t2 = faceOfHalf[ab.second];
            if (faces[t1].size() != 3 || faces[t2].size() != 3) continue;
            const int k1 = cornerOfHalf[ab.first], k2 = cornerOfHalf[ab.second];
            const int hp = h.ofFace[t1][(k1 + 1) % 3], hq = h.ofFace[t1][(k1 + 2) % 3];
            const int hr = h.ofFace[t2][(k2 + 1) % 3], hs = h.ofFace[t2][(k2 + 2) % 3];
            addRow({{hp, 1}, {hs, 1}, {hq, -1}, {hr, -1}, {ab.first, -1}, {ab.second, -1}}, -kHighsInf, 1);
            addRow({{hq, 1}, {hr, 1}, {hp, -1}, {hs, -1}, {ab.first, -1}, {ab.second, -1}}, -kHighsInf, 1);
            ++p.pairRules;
        }
        // A row end whose two course edges border the same row forces a wale-direction mismatch in stage 3: that row
        // has both edges on its side towards the apex, so both are top or both bottom, while the special quad of
        // section 7 has one of each (sphere, 2026-11-13: every one of its 18 mismatches). The two neighbours are in one
        // row iff the row reaches around the apex through wale edges only: for a triangle corner at an interior vertex
        // with edges e, e' and the other edges E at that vertex, sum(E) >= x_e + x_e' - 1 (at least one of E is course
        // when both of the corner's edges are). The others' roles are read on one half-edge each.
        std::vector<std::vector<int>> halvesAtVertex(h.edgeKeys.size() > 0 ? 1 + static_cast<size_t>(std::max_element(h.edgeKeys.begin(), h.edgeKeys.end(),
                                                                                                            [](const EdgeKey& a, const EdgeKey& b) { return a.second < b.second; })->second)
                                                                             : 0);
        std::vector<char> boundaryVertex(halvesAtVertex.size(), 0);
        std::vector<int> halvesOfEdge(h.edgeKeys.size(), 0);
        for (int i = 0; i < nH; ++i) ++halvesOfEdge[h.edgeOfHalf[i]];
        for (int i = 0; i < nH; ++i) {
            const EdgeKey& e = h.edgeKeys[h.edgeOfHalf[i]];
            halvesAtVertex[e.first].push_back(i);
            halvesAtVertex[e.second].push_back(i);
            if (halvesOfEdge[h.edgeOfHalf[i]] == 1) boundaryVertex[e.first] = boundaryVertex[e.second] = 1;
        }
        for (size_t f = 0; f < faces.size(); ++f) {
            if (faces[f].size() != 3) continue;
            for (int c = 0; c < 3; ++c) {
                const int v = faces[f][c];
                if (boundaryVertex[v]) continue;
                const int he = h.ofFace[f][c], heIn = h.ofFace[f][(c + 2) % 3];  // the corner's two edges
                const int edgeOut = h.edgeOfHalf[he], edgeIn = h.edgeOfHalf[heIn];
                std::vector<std::pair<int, double>> terms{{he, -1.0}, {heIn, -1.0}};
                std::vector<char> seen(h.edgeKeys.size(), 0);
                for (int other : halvesAtVertex[v]) {
                    const int edge = h.edgeOfHalf[other];
                    if (edge == edgeOut || edge == edgeIn || seen[edge]) continue;
                    seen[edge] = 1;
                    terms.push_back({other, 1.0});
                }
                for (const std::pair<int, double>& t : terms) {
                    index.push_back(t.first);
                    value.push_back(t.second);
                }
                start.push_back(static_cast<HighsInt>(index.size()));
                lp.row_lower_.push_back(-1.0);
                lp.row_upper_.push_back(kHighsInf);
                ++p.apexRules;
            }
        }
    }
    lp.num_row_ = static_cast<HighsInt>(lp.row_lower_.size());
    lp.a_matrix_.format_ = MatrixFormat::kRowwise;
    lp.a_matrix_.num_col_ = nCol;
    lp.a_matrix_.num_row_ = lp.num_row_;
    lp.a_matrix_.start_ = start;
    lp.a_matrix_.index_ = index;
    lp.a_matrix_.value_ = value;
    Highs highs;
    highs.setOptionValue("output_flag", false);
    highs.setOptionValue("threads", 1);
    highs.setOptionValue("random_seed", 0);
    highs.setOptionValue("mip_rel_gap", 0.0);
    highs.setOptionValue("mip_abs_gap", 0.99);  // the objective is integral
    highs.setOptionValue("time_limit", o.timeLimitS);
    const auto t0 = std::chrono::steady_clock::now();
    if (highs.passModel(lp) == HighsStatus::kError) {
        if (err) *err = "labels: HiGHS refused the half-edge program";
        return false;
    }
    if (o.heuristicStart) {
        std::vector<char> preferCourse;
        const Vec3 up = norm(o.up) > 0 ? normalized(o.up) : Vec3{0, 0, 1};
        for (const EdgeKey& e : h.edgeKeys) {
            const Vec3 d = pos[e.second] - pos[e.first];
            preferCourse.push_back(std::fabs(dot(d, up)) < 0.70710678 * norm(d) ? 1 : 0);
        }
        int parityConflicts = 0, startConflicts = 0;
        const std::vector<char> xs = heuristicStart(faces, h, preferCourse, parityConflicts, startConflicts);
        HighsSolution sol;
        sol.col_value.assign(nCol, 0.0);
        for (int i = 0; i < nH; ++i) sol.col_value[i] = xs[i];
        for (int e = 0; e < nE; ++e) sol.col_value[nH + e] = xs[h.interior[e].first] != xs[h.interior[e].second] ? 1.0 : 0.0;
        sol.value_valid = true;
        highs.setSolution(sol);
        p.report.counts.push_back({"start conflicts", startConflicts});
        p.report.counts.push_back({"start parity conflicts", parityConflicts});
    }
    highs.run();
    p.msSolve = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    const HighsModelStatus status = highs.getModelStatus();
    p.optimal = status == HighsModelStatus::kOptimal;
    p.solverStatus = p.optimal ? "optimal" : status == HighsModelStatus::kTimeLimit ? "time limit" : highs.modelStatusToString(status);
    const HighsInfo& info = highs.getInfo();
    if (info.primal_solution_status != kSolutionStatusFeasible) {
        if (err) *err = "labels: no feasible half-edge labelling (" + p.solverStatus + ")";
        return false;
    }
    p.bound = info.mip_dual_bound / disagreementCost;
    p.gap = info.mip_gap;
    const std::vector<double>& col = highs.getSolution().col_value;
    x.assign(nH, 0);
    for (int i = 0; i < nH; ++i) {
        if (col[i] > 0.01 && col[i] < 0.99) {  // never a rounded fraction
            if (err) *err = "labels: a fractional half-edge value in the solution";
            return false;
        }
        x[i] = col[i] > 0.5 ? 1 : 0;
    }
    // the objective reported is the paper's: the disagreeing edges; the triangles' preference is counted apart
    p.objective = 0;
    for (const std::pair<int, int>& e : h.interior) p.objective += x[e.first] != x[e.second];
    p.rowEndTriangles = 0;
    for (size_t f = 0; f < faces.size(); ++f) {
        if (faces[f].size() != 3) continue;
        int courses = 0;
        for (int hf : h.ofFace[f]) courses += x[hf];
        p.rowEndTriangles += courses == 2;
    }
    return true;
}

// --- 5.2: the disagreeing edges get one role by splitting a quad or rotating a triangle pair

bool convexQuad(const std::vector<Vec3>& pos, const std::vector<int>& q);   // section 5.3 below

void resolveConflicts(FaceList& P, std::map<EdgeKey, char>& directed, const HalfEdges& h, const std::vector<Vec3>& pos, LabelProduct& p) {
    std::vector<char> mustSplit(P.faces.size(), 0);
    std::vector<EdgeKey> rotate;
    for (const std::pair<int, int>& halves : h.interior) {
        // the directed edges of the two half-edges: the face lists are the same as when the half-edges were numbered
        int fa = -1, ka = -1, fb = -1, kb = -1;
        for (size_t f = 0; f < h.ofFace.size() && (fa < 0 || fb < 0); ++f)
            for (size_t k = 0; k < h.ofFace[f].size(); ++k) {
                if (h.ofFace[f][k] == halves.first) fa = static_cast<int>(f), ka = static_cast<int>(k);
                if (h.ofFace[f][k] == halves.second) fb = static_cast<int>(f), kb = static_cast<int>(k);
            }
        const std::vector<int>& A = P.faces[fa];
        const int u = A[ka], v = A[(ka + 1) % A.size()];
        (void) kb;
        const char luv = directed.at({u, v}), lvu = directed.at({v, u});
        if (luv == lvu) continue;
        ++p.conflicts;
        const bool qa = P.faces[fa].size() == 4, qb = P.faces[fb].size() == 4;
        if (qa && qb) {  // course; the quad whose half-edge said wale is split
            ++p.conflictQuadQuad;
            mustSplit[luv == 0 ? fa : fb] = 1;
            directed[{u, v}] = directed[{v, u}] = 1;
        } else if (qa != qb) {  // the triangle's role; the quad is split
            ++p.conflictQuadTriangle;
            const char l = qa ? lvu : luv;
            mustSplit[qa ? fa : fb] = 1;
            directed[{u, v}] = directed[{v, u}] = l;
        } else {
            ++p.conflictTriangleTriangle;
            rotate.push_back(edgeKey(u, v));
        }
    }
    auto roleOf = [&](int a, int b) { return static_cast<int>(directed.at({a, b})); };
    const size_t before = P.faces.size();
    for (size_t fi = 0; fi < before; ++fi) {
        if (!mustSplit[fi]) continue;
        const std::vector<int> q = P.faces[fi];
        const int L[4] = {roleOf(q[0], q[1]), roleOf(q[1], q[2]), roleOf(q[2], q[3]), roleOf(q[3], q[0])};
        if (quadValid(L[0], L[1], L[2], L[3])) {  // every role changed around it: still a valid quad
            ++p.keptQuads;
            continue;
        }
        // the diagonal whose two triangles can both be valid (Fig. 6: one always can); both: the shorter
        auto diagonalRole = [&](int sideA, int sideB, int& role) {
            for (int c : {0, 1})
                if (triangleValid(sideA + c) && triangleValid(sideB + c)) {
                    role = c;
                    return true;
                }
            return false;
        };
        int d02 = 0, d13 = 0;
        const bool ok02 = diagonalRole(L[0] + L[1], L[2] + L[3], d02) && !P.isEdge(q[0], q[2]);
        const bool ok13 = diagonalRole(L[1] + L[2], L[3] + L[0], d13) && !P.isEdge(q[1], q[3]);
        if (!ok02 && !ok13) {
            ++p.unsplittableQuads;
            continue;
        }
        const bool use02 = ok02 && (!ok13 || norm(pos[q[0]] - pos[q[2]]) <= norm(pos[q[1]] - pos[q[3]]));
        P.kill(static_cast<int>(fi));
        if (use02) {
            P.add({q[0], q[1], q[2]});
            P.add({q[2], q[3], q[0]});
            directed[{q[0], q[2]}] = directed[{q[2], q[0]}] = static_cast<char>(d02);
        } else {
            P.add({q[1], q[2], q[3]});
            P.add({q[3], q[0], q[1]});
            directed[{q[1], q[3]}] = directed[{q[3], q[1]}] = static_cast<char>(d13);
        }
        ++p.splitQuads;
    }
    for (const EdgeKey& e : rotate) {  // [u, v, a] | [v, u, b] -> [a, u, b] | [b, v, a], the new edge a - b course
        const int u = e.first, v = e.second;
        const auto ia = P.faceOfDirected.find({u, v}), ib = P.faceOfDirected.find({v, u});
        if (ia == P.faceOfDirected.end() || ib == P.faceOfDirected.end()) {
            ++p.unrotatablePairs;
            continue;
        }
        const int fa = ia->second, fb = ib->second;
        if (P.faces[fa].size() != 3 || P.faces[fb].size() != 3 || roleOf(u, v) == roleOf(v, u)) {
            if (roleOf(u, v) != roleOf(v, u)) ++p.unrotatablePairs;
            continue;
        }
        const int a = startAt(P.faces[fa], u)[2], b = startAt(P.faces[fb], v)[2];
        // the rotated diagonal must lie inside the pair's quad (a, u, b, v): in a dart it would cross outside and invert a
        // triangle (as section 5.3's merges refuse darts); refused, the conflict stays and the validator reports it
        if (a == b || P.isEdge(a, b) || !convexQuad(pos, {a, u, b, v})) {
            ++p.unrotatablePairs;
            continue;
        }
        P.kill(fa);
        P.kill(fb);
        directed.erase({u, v});
        directed.erase({v, u});
        P.add({a, u, b});
        P.add({b, v, a});
        directed[{a, b}] = directed[{b, a}] = 1;
        ++p.rotatedPairs;
    }
}

// --- 5.3: triangle pairs merge into valid quads; a course edge between triangles flips when that lets them merge;
// a wale edge between two two-wale triangles that do not merge becomes course (section 7 needs it)

// A merged quad must be convex about its own normal: two triangles whose third corners lie on the same side of their
// shared edge make a dart (a reflex corner), which stage 4's subdivision turns into an inverted cell (the car's
// coarse face 923, 2026-11-13). Corners turning the other way by more than a hundredth of the quad's area refuse it.
bool convexQuad(const std::vector<Vec3>& pos, const std::vector<int>& q) {
    Vec3 n{0, 0, 0};
    for (size_t k = 0; k < 4; ++k) n += cross(pos[q[k]], pos[q[(k + 1) % 4]]);   // Newell
    const double area2 = norm(n);
    if (area2 <= 0) return false;
    for (size_t k = 0; k < 4; ++k) {
        const Vec3& a = pos[q[(k + 3) % 4]];
        const Vec3& c = pos[q[k]];
        const Vec3& d = pos[q[(k + 1) % 4]];
        if (dot(cross(c - a, d - c), n) <= 0.01 * area2 * area2) return false;
    }
    return true;
}

void mergeTriangles(FaceList& P, std::map<EdgeKey, char>& roles, const std::vector<Vec3>& pos, LabelProduct& p) {
    auto triangleAcross = [&](int fi, int a, int b) {
        const int g = P.across(a, b);
        return g >= 0 && P.alive[g] && P.faces[g].size() == 3 && g != fi ? g : -1;
    };
    auto tryMerge = [&](int fa, int u, int v, bool apply) {
        const int fb = triangleAcross(fa, u, v);
        if (fb < 0 || P.faces[fa].size() != 3) return false;
        const std::vector<int> q = mergedQuad(P.faces[fa], P.faces[fb], u, v);
        if (q[1] == q[3] || !faceValid(roles, q)) return false;  // a folded pair is no quad
        if (!convexQuad(pos, q)) {
            ++p.dartsRefused;
            return false;
        }
        if (apply) {
            P.kill(fa);
            P.kill(fb);
            roles.erase(edgeKey(u, v));
            P.add(q);
        }
        return true;
    };
    auto mergePass = [&]() {
        int n = 0;
        for (bool changed = true; changed;) {
            changed = false;
            for (size_t fi = 0; fi < P.faces.size(); ++fi) {
                if (!P.alive[fi] || P.faces[fi].size() != 3) continue;
                const std::vector<int> t = P.faces[fi];
                for (int k = 0; k < 3; ++k)
                    if (tryMerge(static_cast<int>(fi), t[k], t[(k + 1) % 3], true)) {
                        ++n;
                        changed = true;
                        break;
                    }
            }
        }
        return n;
    };
    p.merged += mergePass();
    std::vector<EdgeKey> courseEdges;
    for (const auto& [e, l] : roles)
        if (l == 1) courseEdges.push_back(e);
    for (const EdgeKey& e : courseEdges) {  // Fig. 10
        const auto ia = P.faceOfDirected.find({e.first, e.second}), ib = P.faceOfDirected.find({e.second, e.first});
        if (ia == P.faceOfDirected.end() || ib == P.faceOfDirected.end()) continue;
        const int fa = ia->second, fb = ib->second;
        if (P.faces[fa].size() != 3 || P.faces[fb].size() != 3 || roles.at(e) != 1) continue;
        auto otherTriangle = [&](int f, int notThis) {
            const std::vector<int>& t = P.faces[f];
            for (int k = 0; k < 3; ++k) {
                const int g = triangleAcross(f, t[k], t[(k + 1) % 3]);
                if (g >= 0 && g != notThis) return true;
            }
            return false;
        };
        if (!otherTriangle(fa, fb) && !otherTriangle(fb, fa)) continue;
        roles[e] = 0;
        auto mergeEdge = [&](int f, int notThis, int& eu, int& ev) {
            const std::vector<int>& t = P.faces[f];
            for (int k = 0; k < 3; ++k) {
                const int a = t[k], b = t[(k + 1) % 3];
                if (edgeKey(a, b) == e) continue;
                const int g = triangleAcross(f, a, b);
                if (g < 0 || g == notThis) continue;
                const std::vector<int> q = mergedQuad(t, P.faces[g], a, b);
                if (q[1] != q[3] && faceValid(roles, q)) {
                    eu = a;
                    ev = b;
                    return true;
                }
            }
            return false;
        };
        int au = -1, av = -1, bu = -1, bv = -1;
        const bool ma = mergeEdge(fa, fb, au, av), mb = mergeEdge(fb, fa, bu, bv);
        const bool aOk = ma || triangleValid(courseCount(roles, P.faces[fa])), bOk = mb || triangleValid(courseCount(roles, P.faces[fb]));
        if (!(ma || mb) || !aOk || !bOk) {
            roles[e] = 1;
            continue;
        }
        int done = 0;
        if (ma && tryMerge(fa, au, av, true)) ++done;
        if (mb && P.alive[fb] && tryMerge(fb, bu, bv, true)) ++done;
        p.flipMerged += done;
    }
    p.merged += mergePass();
    std::vector<EdgeKey> waleEdges;
    for (const auto& [e, l] : roles)
        if (l == 0) waleEdges.push_back(e);
    for (const EdgeKey& e : waleEdges) {
        const auto ia = P.faceOfDirected.find({e.first, e.second}), ib = P.faceOfDirected.find({e.second, e.first});
        if (ia == P.faceOfDirected.end() || ib == P.faceOfDirected.end()) continue;
        const int fa = ia->second, fb = ib->second;
        if (P.faces[fa].size() != 3 || P.faces[fb].size() != 3) continue;
        if (3 - courseCount(roles, P.faces[fa]) != 2 || 3 - courseCount(roles, P.faces[fb]) != 2) continue;
        // every pair the merge pass could take is gone; what is left here was not a valid quad or was refused as a
        // dart, and stage 4 cannot subdivide it either way: the paper's flip
        roles[e] = 1;
        ++p.walePairFlips;
    }
}

} // namespace

// --- EdgeLabels, the validator, the stage

bool EdgeLabels::has(int a, int b) const { return role.count(edgeKey(a, b)) > 0; }
EdgeRole EdgeLabels::at(int a, int b) const { return role.at(edgeKey(a, b)); }
int EdgeLabels::courseEdges() const {
    int n = 0;
    for (const auto& [e, r] : role) n += r == EdgeRole::Course;
    return n;
}

LabelCheck validateLabels(const PolygonMesh& m, const EdgeLabels& labels) {
    LabelCheck c;
    c.faces = m.faces();
    std::map<EdgeKey, std::vector<int>> facesOfEdge;
    for (int f = 0; f < m.faces(); ++f) {
        const int s = m.sides(f);
        int courses = 0, l[4] = {0, 0, 0, 0};
        bool complete = true;
        for (int k = 0; k < s; ++k) {
            const int a = m.vertexAt(f, k), b = m.vertexAt(f, (k + 1) % s);
            facesOfEdge[edgeKey(a, b)].push_back(f);
            if (!labels.has(a, b)) {
                complete = false;
                ++c.unlabeledEdges;
                continue;
            }
            const int r = labels.at(a, b) == EdgeRole::Course ? 1 : 0;
            if (k < 4) l[k] = r;
            courses += r;
        }
        if (!complete) {
            ++c.invalidFaces;
            continue;
        }
        const bool ok = s == 4 ? quadValid(l[0], l[1], l[2], l[3]) : s == 3 ? triangleValid(courses) : false;
        if (!ok) ++c.invalidFaces;
    }
    for (const auto& [e, fs] : facesOfEdge) {
        if (fs.size() != 2 || !labels.has(e.first, e.second) || labels.at(e.first, e.second) != EdgeRole::Wale) continue;
        bool twoWale[2] = {false, false};
        for (int i = 0; i < 2; ++i) {
            const int f = fs[i];
            if (m.sides(f) != 3) continue;
            int wales = 0;
            for (int k = 0; k < 3; ++k) wales += labels.has(m.vertexAt(f, k), m.vertexAt(f, (k + 1) % 3)) && labels.at(m.vertexAt(f, k), m.vertexAt(f, (k + 1) % 3)) == EdgeRole::Wale;
            twoWale[i] = wales == 2;
        }
        if (twoWale[0] && twoWale[1]) ++c.twoWalePairs;
    }
    c.status = c.invalidFaces == 0 && c.unlabeledEdges == 0 ? Status::Pass : Status::Fail;
    if (c.status != Status::Pass) c.detail = std::to_string(c.invalidFaces) + " faces without a valid configuration, " + std::to_string(c.unlabeledEdges) + " unlabeled edge uses";
    return c;
}

bool labelEdges(const PolygonMesh& m, const LabelOptions& o, LabelProduct& out, std::string* err) {
    const auto t0 = std::chrono::steady_clock::now();
    out = LabelProduct{};
    out.report.stage = "labels";
    out.report.algorithmVersion = 1;
    out.report.inputDigest = m.digest();
    char settings[200];
    std::snprintf(settings, sizeof settings, "    time limit = %.0f s\n    triangulate odd vertices = %s\n    post-process = %s\n    heuristic start = %s\n    prefer shaping = %s\n",
                  o.timeLimitS, o.triangulateOdd ? "yes" : "no", o.postProcess ? "yes" : "no", o.heuristicStart ? "yes" : "no", o.preferShaping ? "yes" : "no");
    out.report.settings = settings;
    auto fail = [&](const std::string& why) {
        if (err) *err = "labels: " + why;
        return false;
    };
    for (int f = 0; f < m.faces(); ++f)
        if (m.sides(f) < 3 || m.sides(f) > 4) return fail("face " + std::to_string(f) + " has " + std::to_string(m.sides(f)) + " sides: normalize the mesh first");
    std::vector<std::vector<int>> faces = facesOf(m);
    const Topology before = topologyOf(faces);
    if (before.overusedEdges) return fail(std::to_string(before.overusedEdges) + " edges with more than two faces");
    {
        std::set<EdgeKey> directed;
        for (const std::vector<int>& f : faces)
            for (size_t k = 0; k < f.size(); ++k)
                if (!directed.insert({f[k], f[(k + 1) % f.size()]}).second) return fail("a directed edge appears twice: inconsistent winding");
    }
    if (o.triangulateOdd) faces = triangulateOddVertices(faces, before, out);
    const HalfEdges h = halfEdgesOf(faces);
    out.halfEdges = h.count;
    out.interiorEdges = static_cast<int>(h.interior.size());
    std::vector<char> x;
    if (!solveHalfEdgeRoles(faces, h, m.pos, o, out, x, err)) return false;
    std::map<EdgeKey, char> directed;
    for (size_t f = 0; f < faces.size(); ++f)
        for (size_t k = 0; k < faces[f].size(); ++k) directed[{faces[f][k], faces[f][(k + 1) % faces[f].size()]}] = x[h.ofFace[f][k]];
    FaceList P;
    for (const std::vector<int>& f : faces) P.add(f);
    resolveConflicts(P, directed, h, m.pos, out);
    std::map<EdgeKey, char> roles;  // undirected; a pair still disagreeing (unrotatable) counts course and shows as invalid
    for (const auto& [de, l] : directed) {
        const EdgeKey k = edgeKey(de.first, de.second);
        const auto it = roles.find(k);
        if (it == roles.end()) roles[k] = l;
        else if (it->second != l) it->second = 1;
    }
    if (o.postProcess) mergeTriangles(P, roles, m.pos, out);
    out.mesh.pos = m.pos;
    for (const std::vector<int>& f : P.live()) out.mesh.addFace(f);
    for (int f = 0; f < out.mesh.faces(); ++f)
        for (int k = 0; k < out.mesh.sides(f); ++k) {
            const EdgeKey e = edgeKey(out.mesh.vertexAt(f, k), out.mesh.vertexAt(f, (k + 1) % out.mesh.sides(f)));
            const auto it = roles.find(e);
            out.labels.role[e] = it != roles.end() && it->second == 0 ? EdgeRole::Wale : EdgeRole::Course;
        }
    out.check = validateLabels(out.mesh, out.labels);
    const Topology after = topologyOf(facesOf(out.mesh));
    out.sameVertices = after.vertices == before.vertices;
    out.sameBoundary = after.boundary == before.boundary;
    out.sameEuler = after.euler == before.euler;
    out.eulerBefore = before.euler;
    out.eulerAfter = after.euler;
    out.report.outputDigest = out.mesh.digest();
    out.report.counts.push_back({"faces in", m.faces()});
    out.report.counts.push_back({"odd interior vertices", out.oddVertices});
    out.report.counts.push_back({"quads triangulated", out.triangulatedQuads});
    out.report.counts.push_back({"half-edges", out.halfEdges});
    out.report.counts.push_back({"interior edges", out.interiorEdges});
    out.report.counts.push_back({"objective (disagreeing edges)", out.objective});
    out.report.counts.push_back({"row-end triangles (program)", out.rowEndTriangles});
    out.report.counts.push_back({"triangle pair rules", out.pairRules});
    out.report.counts.push_back({"row-end apex rules", out.apexRules});
    out.report.counts.push_back({"best bound", out.bound});
    out.report.counts.push_back({"ms solve", out.msSolve});
    out.report.counts.push_back({"conflicts resolved", out.conflicts});
    out.report.counts.push_back({"quads split", out.splitQuads});
    out.report.counts.push_back({"pairs rotated", out.rotatedPairs});
    out.report.counts.push_back({"unrotatable pairs", out.unrotatablePairs});
    out.report.counts.push_back({"triangles merged", out.merged});
    out.report.counts.push_back({"flip-merged", out.flipMerged});
    out.report.counts.push_back({"merges refused as darts", out.dartsRefused});
    out.report.counts.push_back({"wale-pair flips", out.walePairFlips});
    out.report.counts.push_back({"faces out", out.mesh.faces()});
    out.report.counts.push_back({"course edges", out.labels.courseEdges()});
    out.report.counts.push_back({"invalid faces", out.check.invalidFaces});
    out.report.counts.push_back({"two-wale pairs", out.check.twoWalePairs});
    const bool topologyKept = out.sameVertices && out.sameBoundary && out.sameEuler && after.overusedEdges == 0 && after.degenerateFaces == 0;
    if (out.check.status != Status::Pass) {
        out.report.status = Status::Fail;
        out.report.detail = "validator: " + out.check.detail;
    } else if (!topologyKept) {
        out.report.status = Status::Fail;
        out.report.detail = "the topology changed: vertices " + std::string(out.sameVertices ? "kept" : "CHANGED") + ", boundary " + (out.sameBoundary ? "kept" : "CHANGED") +
                            ", Euler " + std::to_string(out.eulerBefore) + " -> " + std::to_string(out.eulerAfter);
    } else if (!out.optimal) {
        out.report.status = Status::Unresolved;
        out.report.detail = "feasible, not proven optimal: " + out.solverStatus;
    } else {
        out.report.status = Status::Pass;
        out.report.detail = "every face valid, topology kept, program " + out.solverStatus;
    }
    out.report.elapsedMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    return true;
}

double courseToWaleLength(const PolygonMesh& m, const EdgeLabels& labels) {
    double course = 0, wale = 0;
    int courses = 0, wales = 0;
    for (const auto& [edge, role] : labels.role) {
        const double l = norm(m.pos[edge.second] - m.pos[edge.first]);
        if (role == EdgeRole::Course) { course += l; ++courses; }
        else { wale += l; ++wales; }
    }
    if (!courses || !wales || wale <= 0) return 0;
    return (course / courses) / (wale / wales);
}

} // namespace knit
