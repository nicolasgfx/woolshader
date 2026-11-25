// Stage 4: the typed stitch mesh (stitch_mesh.h). In the order it happens: the subdivision with its roles, the apex
// rewrite into pentagons, the assembly, the strands, the loops; then the independent validator.
#include "stitchmesh/stitch_mesh.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <map>
#include <string>
#include <utility>

#include "model/mesh.h"
#include "model/mesh_check.h"

namespace knit {

const char* toString(StitchFaceKind k) {
    switch (k) {
    case StitchFaceKind::Regular: return "regular";
    case StitchFaceKind::Increase: return "increase";
    case StitchFaceKind::Decrease: return "decrease";
    case StitchFaceKind::SpecialQuad: return "special";
    case StitchFaceKind::PoleStart: return "pole-start";
    default: return "unsupported";
    }
}

const char* toString(StitchEdgeRole r) {
    switch (r) {
    case StitchEdgeRole::Wale: return "wale";
    case StitchEdgeRole::Bottom: return "bottom";
    default: return "top";
    }
}

namespace {

using EdgeKey = std::pair<int, int>;

EdgeKey edgeKey(int a, int b) { return {std::min(a, b), std::max(a, b)}; }

StitchEdgeRole oppositeSide(StitchEdgeRole r) { return r == StitchEdgeRole::Bottom ? StitchEdgeRole::Top : StitchEdgeRole::Bottom; }

StitchEdgeRole roleOfFlag(char top) {
    if (top < 0) return StitchEdgeRole::Wale;
    return top ? StitchEdgeRole::Top : StitchEdgeRole::Bottom;
}

// The kind a face's roles (in corner order) support: the paper's four, else unsupported.
StitchFaceKind kindOfRoles(const std::vector<StitchEdgeRole>& r) {
    const int n = static_cast<int>(r.size());
    int wale = 0, bottom = 0, top = 0;
    bool waleAdjacent = false;
    for (int k = 0; k < n; ++k) {
        if (r[k] == StitchEdgeRole::Wale) ++wale;
        else if (r[k] == StitchEdgeRole::Bottom) ++bottom;
        else ++top;
        if (r[k] == StitchEdgeRole::Wale && r[(k + 1) % n] == StitchEdgeRole::Wale) waleAdjacent = true;
    }
    if (wale != 2) return StitchFaceKind::Unsupported;
    if (n == 4 && bottom == 1 && top == 1) return waleAdjacent ? StitchFaceKind::SpecialQuad : StitchFaceKind::Regular;
    if (n == 5 && bottom == 1 && top == 2 && !waleAdjacent) return StitchFaceKind::Increase;
    if (n == 5 && bottom == 2 && top == 1 && !waleAdjacent) return StitchFaceKind::Decrease;
    return StitchFaceKind::Unsupported;
}

// The wale edge the yarn leaves through when it follows the face's own course direction: the wale edge right after the
// (last) bottom edge, counter-clockwise. -1 for a special quad (its wale edges meet) and an unsupported face.
int exitCorner(const std::vector<StitchEdgeRole>& r, StitchFaceKind kind) {
    if (kind == StitchFaceKind::SpecialQuad || kind == StitchFaceKind::Unsupported) return -1;
    const int n = static_cast<int>(r.size());
    for (int k = 0; k < n; ++k)
        if (r[k] == StitchEdgeRole::Wale && r[(k + n - 1) % n] == StitchEdgeRole::Bottom) return k;
    return -1;
}

// The two wale corners of a face (fewer or more: an unsupported face).
std::vector<int> waleCorners(const std::vector<StitchEdgeRole>& r) {
    std::vector<int> w;
    for (int k = 0; k < static_cast<int>(r.size()); ++k)
        if (r[k] == StitchEdgeRole::Wale) w.push_back(k);
    return w;
}

// A stitch face under construction.
struct DraftFace {
    std::vector<int> corners;
    std::vector<StitchEdgeRole> roles;
    int coarseFace = -1;
    bool removed = false;
};

// --- the subdivision (section 7) ---------------------------------------------------------------------------------------

struct Subdivision {
    std::vector<Vec3> pos;                     // coarse vertices, midpoints, centres
    std::map<EdgeKey, int> midpoint;           // coarse edge -> vertex
    std::vector<int> centre;                   // coarse face -> vertex
    std::vector<DraftFace> faces;              // every sub-face, in (coarse face, corner) order
    std::vector<std::vector<int>> subFaceAt;   // coarse face, corner -> sub-face
};

// The role of the centre edge perpendicular to a half course edge: the other side of the row (the centre line is the top
// of the lower stitch row and the bottom of the upper). For a corner between two wale edges (the apex of a triangle with
// two wale edges) both centre edges are course edges on the side of the coarse face's one course edge.
void subdivide(const PolygonMesh& coarse, const std::vector<std::vector<char>>& top, Subdivision& s) {
    const int nV = static_cast<int>(coarse.pos.size());
    s.pos = coarse.pos;
    for (int f = 0; f < coarse.faces(); ++f) {
        const int n = coarse.sides(f);
        for (int k = 0; k < n; ++k) {
            const EdgeKey e = edgeKey(coarse.vertexAt(f, k), coarse.vertexAt(f, (k + 1) % n));
            if (s.midpoint.count(e)) continue;
            s.midpoint[e] = static_cast<int>(s.pos.size());
            s.pos.push_back((coarse.pos[e.first] + coarse.pos[e.second]) * 0.5);
        }
    }
    s.centre.resize(coarse.faces());
    for (int f = 0; f < coarse.faces(); ++f) {
        s.centre[f] = static_cast<int>(s.pos.size());
        s.pos.push_back(coarse.faceCentroid(f));
    }
    (void) nV;
    s.subFaceAt.resize(coarse.faces());
    for (int f = 0; f < coarse.faces(); ++f) {
        const int n = coarse.sides(f);
        // the coarse face's one course side, for an apex corner
        StitchEdgeRole apexSide = StitchEdgeRole::Bottom;
        int courseEdges = 0;
        for (int k = 0; k < n; ++k)
            if (top[f][k] >= 0) {
                apexSide = roleOfFlag(top[f][k]);
                ++courseEdges;
            }
        s.subFaceAt[f].resize(n);
        for (int k = 0; k < n; ++k) {
            const int prev = (k + n - 1) % n;
            const int v = coarse.vertexAt(f, k);
            const int mk = s.midpoint.at(edgeKey(v, coarse.vertexAt(f, (k + 1) % n)));
            const int mPrev = s.midpoint.at(edgeKey(coarse.vertexAt(f, prev), v));
            const StitchEdgeRole r0 = roleOfFlag(top[f][k]);     // v -> mk, half of coarse edge k
            const StitchEdgeRole r3 = roleOfFlag(top[f][prev]);  // mPrev -> v, half of coarse edge k - 1
            const bool wale0 = r0 == StitchEdgeRole::Wale, wale3 = r3 == StitchEdgeRole::Wale;
            StitchEdgeRole r1 = StitchEdgeRole::Wale, r2 = StitchEdgeRole::Wale;  // mk -> centre, centre -> mPrev
            if (wale0 != wale3) {
                const StitchEdgeRole other = oppositeSide(wale0 ? r3 : r0);
                if (wale0) r1 = other;
                else r2 = other;
            } else if (wale0 && wale3) {
                r1 = r2 = courseEdges == 1 ? apexSide : StitchEdgeRole::Bottom;  // two course edges at a wale-wale corner (a quad) stay unsupported
            }
            DraftFace d;
            d.corners = {v, mk, s.centre[f], mPrev};
            d.roles = {r0, r1, r2, r3};
            d.coarseFace = f;
            s.subFaceAt[f][k] = static_cast<int>(s.faces.size());
            s.faces.push_back(d);
        }
    }
}

// --- the apex rewrite: two triangles merged across the apex's wale edges into two pentagons ---------------------------

struct ApexRewrite {
    int merges = 0, unmerged = 0;
    std::vector<std::string> reasons;
};

// The sub-face of the coarse neighbour across half-edge h that holds the other half of the apex's wale edge, or -1.
int neighbourSubFace(const PolygonMesh& coarse, const HalfEdgeMesh& half, const Subdivision& s, int h, int apexVertex) {
    const int twin = half.twin[h];
    if (twin < 0) return -1;
    const int g = half.face[twin];
    const int n = coarse.sides(g);
    for (int k = 0; k < n; ++k)
        if (coarse.vertexAt(g, k) == apexVertex) return s.subFaceAt[g][k];
    return -1;
}

bool isPlainQuad(const DraftFace& d) { return !d.removed && d.corners.size() == 4 && kindOfRoles(d.roles) == StitchFaceKind::Regular; }

void rewriteApexFaces(const PolygonMesh& coarse, const HalfEdgeMesh& half, const std::vector<std::vector<char>>& top,
                      Subdivision& s, ApexRewrite& out) {
    for (int f = 0; f < coarse.faces(); ++f) {
        if (coarse.sides(f) != 3) continue;
        int waleEdges = 0, apexCorner = -1;
        for (int k = 0; k < 3; ++k)
            if (top[f][k] < 0) ++waleEdges;
        if (waleEdges != 2) continue;
        for (int k = 0; k < 3; ++k)
            if (top[f][k] < 0 && top[f][(k + 2) % 3] < 0) apexCorner = k;  // both edges at corner k are wale
        if (apexCorner < 0) continue;
        const int a = coarse.vertexAt(f, apexCorner);
        DraftFace& apex = s.faces[s.subFaceAt[f][apexCorner]];
        // corners (a, mk, c, mPrev): roles (wale, course, course, wale)
        const int mk = apex.corners[1], c = apex.corners[2], mPrev = apex.corners[3];
        const int hNext = coarse.faceStart[f] + apexCorner;                 // a -> v(k+1), the wale edge after the apex
        const int hPrev = coarse.faceStart[f] + (apexCorner + 2) % 3;       // v(k-1) -> a, the wale edge before the apex
        const int t1 = neighbourSubFace(coarse, half, s, hNext, a);
        const int t2 = neighbourSubFace(coarse, half, s, hPrev, a);
        std::string why;
        if (t1 < 0 || t2 < 0) why = "an apex wale edge lies on the mesh border";
        else if (t1 == t2) why = "both apex wale edges lead to the same neighbouring stitch face";
        else if (!isPlainQuad(s.faces[t1]) || !isPlainQuad(s.faces[t2])) why = "a neighbouring stitch face is already rewritten or not regular";
        else if (s.faces[t1].corners[0] != a || s.faces[t1].corners[3] != mk || s.faces[t2].corners[0] != a || s.faces[t2].corners[1] != mPrev)
            why = "the neighbouring stitch faces do not share the apex's half edges as expected";
        if (!why.empty()) {
            ++out.unmerged;
            out.reasons.push_back("coarse triangle " + std::to_string(f) + " at vertex " + std::to_string(a) + ": " + why);
            continue;
        }
        DraftFace& n1 = s.faces[t1];  // (a, mN, cN, mk): its last edge mk -> a is the shared half
        DraftFace& n2 = s.faces[t2];  // (a, mPrev, cN2, mN2): its first edge a -> mPrev is the shared half
        DraftFace p1;
        p1.corners = {a, n1.corners[1], n1.corners[2], mk, c};
        p1.roles = {n1.roles[0], n1.roles[1], n1.roles[2], apex.roles[1], StitchEdgeRole::Wale};
        p1.coarseFace = n1.coarseFace;
        DraftFace p2;
        p2.corners = {c, mPrev, n2.corners[2], n2.corners[3], a};
        p2.roles = {apex.roles[2], n2.roles[1], n2.roles[2], n2.roles[3], StitchEdgeRole::Wale};
        p2.coarseFace = n2.coarseFace;
        n1 = p1;
        n2 = p2;
        apex.removed = true;
        out.merges += 2;
    }
}

// --- the strands: faces joined across their wale edges ---------------------------------------------------------------

struct WaleStep {
    int face = -1, entryCorner = -1;   // the face entered and the corner edge it was entered through
};

WaleStep stepAcross(const HalfEdgeMesh& half, const PolygonMesh& m, int face, int corner) {
    WaleStep s;
    const int twin = half.twin[m.faceStart[face] + corner];
    if (twin < 0) return s;
    s.face = half.face[twin];
    s.entryCorner = twin - m.faceStart[s.face];
    return s;
}

void traceStrands(const HalfEdgeMesh& half, StitchMesh& s) {
    const PolygonMesh& m = s.mesh;
    const int nF = m.faces();
    std::vector<std::vector<int>> wale(nF);
    std::vector<int> exitOf(nF, -1);
    s.strandOfFace.assign(nF, -1);
    s.reversed.assign(nF, 0);
    s.yarnEnds = 0;
    s.unsupportedEnds = 0;
    for (int f = 0; f < nF; ++f) {
        wale[f] = waleCorners(s.role[f]);
        exitOf[f] = exitCorner(s.role[f], s.kind[f]);
    }
    auto supported = [&](int f) { return f >= 0 && wale[f].size() == 2; };
    auto otherWale = [&](int f, int corner) { return wale[f][0] == corner ? wale[f][1] : wale[f][0]; };
    for (int f = 0; f < nF; ++f)
        for (int k : wale[f]) {
            const WaleStep n = stepAcross(half, m, f, k);
            if (n.face < 0) ++s.yarnEnds;
            else if (supported(f) != supported(n.face) && f < n.face) ++s.unsupportedEnds;
            else if (!supported(f) && !supported(n.face) && f < n.face) ++s.unsupportedEnds;
        }
    std::vector<char> visited(nF, 0);
    for (int f0 = 0; f0 < nF; ++f0) {
        if (visited[f0] || !supported(f0)) continue;
        // walk back through the first wale edge to the strand's start (a border, an unsupported face, or f0 again: a cycle)
        int start = f0, entry = wale[f0][0];
        bool closed = false;
        for (int cur = f0, via = wale[f0][0];;) {
            const WaleStep n = stepAcross(half, m, cur, via);
            if (n.face < 0 || !supported(n.face)) { start = cur; entry = via; break; }
            if (n.face == f0) {  // a cycle: start at f0, entering so that f0 runs in its own course direction when it has one
                start = f0;
                entry = exitOf[f0] >= 0 ? otherWale(f0, exitOf[f0]) : wale[f0][0];
                closed = true;
                break;
            }
            cur = n.face;
            via = otherWale(cur, n.entryCorner);
        }
        Strand strand;
        strand.closed = closed;
        const int id = static_cast<int>(s.strands.size());
        int cur = start, exit = otherWale(start, entry);
        for (;;) {
            strand.faces.push_back(cur);
            visited[cur] = 1;
            s.strandOfFace[cur] = id;
            s.reversed[cur] = exitOf[cur] >= 0 && exit != exitOf[cur];
            const WaleStep n = stepAcross(half, m, cur, exit);
            if (n.face < 0 || !supported(n.face) || n.face == start) break;
            cur = n.face;
            exit = otherWale(cur, n.entryCorner);
        }
        // an open strand has two equivalent ends: run it so that its first face goes in its own course direction
        if (!closed && exitOf[start] >= 0 && s.reversed[start]) {
            std::reverse(strand.faces.begin(), strand.faces.end());
            for (int f : strand.faces) s.reversed[f] = !s.reversed[f];
        }
        s.strands.push_back(strand);
    }
}

// --- the loops: every course edge seen from both sides ----------------------------------------------------------------

std::vector<LoopPort> classifyLoops(const HalfEdgeMesh& half, const PolygonMesh& m, const std::vector<std::vector<StitchEdgeRole>>& role) {
    std::vector<LoopPort> loops;
    for (int f = 0; f < m.faces(); ++f) {
        const int n = m.sides(f);
        for (int k = 0; k < n; ++k) {
            if (role[f][k] == StitchEdgeRole::Wale) continue;
            const int h = m.faceStart[f] + k, twin = half.twin[h];
            const int g = twin < 0 ? -1 : half.face[twin];
            if (g >= 0 && g < f) continue;  // once per edge, from its lower face
            LoopPort p;
            const EdgeKey e = edgeKey(m.vertexAt(f, k), m.vertexAt(f, (k + 1) % n));
            p.a = e.first;
            p.b = e.second;
            p.faceA = f;
            p.faceB = g;
            const bool topA = role[f][k] == StitchEdgeRole::Top;
            if (g < 0) {
                p.kind = topA ? LoopPort::Kind::BoundaryTop : LoopPort::Kind::BoundaryBottom;
                (topA ? p.offering : p.receiving) = f;
            } else {
                const int kg = twin - m.faceStart[g];
                const bool topB = role[g][kg] == StitchEdgeRole::Top;
                if (role[g][kg] == StitchEdgeRole::Wale) {
                    p.kind = LoopPort::Kind::BoundaryTop;  // a course edge against a wale edge: the labels disagree
                    (topA ? p.offering : p.receiving) = f;
                } else if (topA == topB) {
                    p.kind = topA ? LoopPort::Kind::BothTop : LoopPort::Kind::BothBottom;
                } else {
                    p.kind = LoopPort::Kind::Interlock;
                    p.offering = topA ? f : g;
                    p.receiving = topA ? g : f;
                }
            }
            loops.push_back(p);
        }
    }
    return loops;
}

} // namespace

// --- the stage ------------------------------------------------------------------------------------------------------

bool buildStitchMesh(const PolygonMesh& coarse, const std::vector<std::vector<char>>& top, StitchMeshProduct& out, std::string* err) {
    const auto t0 = std::chrono::steady_clock::now();
    out = StitchMeshProduct{};
    if (static_cast<int>(top.size()) != coarse.faces()) {
        if (err) *err = "stitch mesh: top flags for " + std::to_string(top.size()) + " faces, mesh has " + std::to_string(coarse.faces());
        return false;
    }
    for (int f = 0; f < coarse.faces(); ++f) {
        if (coarse.sides(f) < 3 || coarse.sides(f) > 4) {
            if (err) *err = "stitch mesh: face " + std::to_string(f) + " has " + std::to_string(coarse.sides(f)) + " sides";
            return false;
        }
        if (static_cast<int>(top[f].size()) != coarse.sides(f)) {
            if (err) *err = "stitch mesh: face " + std::to_string(f) + " has flags for " + std::to_string(top[f].size()) + " edges";
            return false;
        }
        if (coarse.sides(f) == 4) ++out.coarseQuads;
        else {
            ++out.coarseTriangles;
            int wale = 0;
            for (int k = 0; k < 3; ++k) wale += top[f][k] < 0;
            if (wale == 1) ++out.rowEndTriangles;
            else if (wale == 2) ++out.apexTriangles;
        }
    }
    HalfEdgeMesh half;
    std::string why;
    if (!buildHalfEdges(coarse, half, &why)) {
        if (err) *err = "stitch mesh: " + why;
        return false;
    }
    Subdivision sub;
    subdivide(coarse, top, sub);
    ApexRewrite rewrite;
    rewriteApexFaces(coarse, half, top, sub, rewrite);
    out.merges = rewrite.merges;
    out.unmergedApexFaces = rewrite.unmerged;

    StitchMesh& s = out.stitch;
    s.mesh.pos = sub.pos;
    for (const DraftFace& d : sub.faces) {
        if (d.removed) continue;
        s.mesh.addFace(d.corners);
        s.role.push_back(d.roles);
        s.kind.push_back(kindOfRoles(d.roles));
        s.coarseFace.push_back(d.coarseFace);
    }
    HalfEdgeMesh stitchHalf;
    if (buildHalfEdges(s.mesh, stitchHalf, &why)) {
        traceStrands(stitchHalf, s);
        s.loops = classifyLoops(stitchHalf, s.mesh, s.role);
    }
    out.check = validateStitchMesh(s, coarse);

    StageReport& r = out.report;
    r.stage = "stitch";
    r.inputDigest = coarse.digest();
    r.outputDigest = s.mesh.digest();
    r.elapsedMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    r.counts = {{"coarse quads", out.coarseQuads}, {"coarse triangles", out.coarseTriangles}, {"row-end triangles", out.rowEndTriangles},
                {"apex triangles", out.apexTriangles}, {"pentagons", out.merges}, {"unmerged apex faces", out.unmergedApexFaces},
                {"faces", out.check.faces}, {"regular", out.check.regular}, {"increase", out.check.increase}, {"decrease", out.check.decrease},
                {"special", out.check.special}, {"unsupported", out.check.unsupported}, {"strands", out.check.strands},
                {"closed strands", out.check.closedStrands}, {"yarn ends", out.check.yarnEnds}, {"interlocks", out.check.interlocks},
                {"mismatch both top", out.check.bothTop}, {"mismatch both bottom", out.check.bothBottom}, {"shared edge pairs", out.check.sharedEdgePairs},
                {"boundary bottom", out.check.boundaryBottom}, {"boundary top", out.check.boundaryTop}};
    if (out.check.status == Status::Fail) {
        r.status = Status::Fail;
        r.detail = "the stitch mesh fails its check: " + out.check.detail;
    } else if (out.check.status == Status::Unresolved) {
        r.status = Status::Unresolved;
        r.detail = std::to_string(out.check.unsupported) + " unsupported faces, " + std::to_string(out.check.bothTop + out.check.bothBottom) +
                   " direction mismatches, " + std::to_string(out.check.sharedEdgePairs) +
                   " face pairs sharing two edges (listed; the yarn stage connects them or names them unsupported)";
        for (const std::string& reason : rewrite.reasons) r.detail += "\n    " + reason;
    } else {
        r.status = Status::Pass;
        r.detail = std::to_string(out.check.strands) + " strands (" + std::to_string(out.check.closedStrands) + " closed), " +
                   std::to_string(out.check.interlocks) + " interlocks, every port traced once";
    }
    return true;
}

// --- the surface checks of the validator ---------------------------------------------------------------------------------

// The faces around a vertex walked as a fan through next and twin: one closed cycle (interior) or one open chain
// (boundary) must reach every face incident to the vertex; two fans mean a pinched surface.
int brokenVertexLinks(const PolygonMesh& m, const HalfEdgeMesh& half) {
    const int nV = static_cast<int>(m.pos.size());
    std::vector<int> incident(nV, 0);
    for (int h = 0; h < static_cast<int>(half.origin.size()); ++h) ++incident[half.origin[h]];
    int broken = 0;
    for (int v = 0; v < nV; ++v) {
        if (incident[v] == 0) continue;
        // start at a boundary half-edge out of v when there is one (the fan's open end), else anywhere
        int start = half.halfEdgeOfVertex[v];
        for (int h = 0; h < static_cast<int>(half.origin.size()); ++h)
            if (half.origin[h] == v && half.twin[h] < 0) { start = h; break; }
        int reached = 0;
        for (int h = start, guard = 0; h >= 0 && guard <= incident[v]; ++guard) {
            ++reached;
            // the next face counter-clockwise around v: the twin of the half-edge INTO v, i.e. of prev(h)
            int into = half.next[h];
            while (half.next[into] != h) into = half.next[into];
            const int t = half.twin[into];
            if (t < 0 || t == start) break;
            h = t;
        }
        if (reached != incident[v]) ++broken;
    }
    return broken;
}

// A face's corners in its own plane: positive area with the coarse face's normal, and the corner turns.
void faceGeometry(const StitchMesh& s, const PolygonMesh& coarse, int& inverted, int& nonConvex, std::vector<int>& invertedIds) {
    inverted = nonConvex = 0;
    const PolygonMesh& m = s.mesh;
    for (int f = 0; f < m.faces(); ++f) {
        const Vec3 n = m.faceNormal(f);
        const Vec3 reference = s.coarseFace[f] >= 0 && s.coarseFace[f] < coarse.faces() ? coarse.faceNormal(s.coarseFace[f]) : n;
        if (m.faceArea(f) <= 0 || dot(n, reference) <= 0) {
            ++inverted;
            if (invertedIds.size() < 64) invertedIds.push_back(f);
            continue;
        }
        const int k = m.sides(f);
        bool convex = true;
        for (int i = 0; i < k; ++i) {
            const Vec3& a = m.pos[m.vertexAt(f, i)];
            const Vec3& b = m.pos[m.vertexAt(f, (i + 1) % k)];
            const Vec3& c = m.pos[m.vertexAt(f, (i + 2) % k)];
            if (dot(cross(b - a, c - b), n) < 0) convex = false;
        }
        if (!convex) ++nonConvex;
    }
}

// The faces fanned from their centroids as triangles (right for the star-shaped faces of a stitch mesh), through the
// triangle mesh's self-intersection test.
// A face triangulated in its own plane by ear clipping (the corners projected along the face's normal; an ear is a
// convex corner whose triangle holds no other corner), validated: n - 2 triangles, every one turning the face's way.
// A centroid fan is not used: a face that is not star-shaped about its centroid gets a reversed fan triangle, and an
// intersection found with such a fan says nothing about the surface. False when no valid ear clipping exists (the
// polygon is self-intersecting or degenerate in its plane): the caller leaves the verdict unresolved.
bool triangulateFace(const PolygonMesh& m, int f, std::vector<std::array<int, 3>>& out) {
    const int n = m.sides(f);
    const Vec3 normal = m.faceNormal(f);
    if (n < 3 || norm(normal) == 0) return false;
    // a basis of the plane
    Vec3 u = cross(normal, std::fabs(normal.x) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0});
    u = normalized(u);
    const Vec3 v = cross(normal, u);
    std::vector<double> px(n), py(n);
    std::vector<int> ids(n);
    for (int k = 0; k < n; ++k) {
        const Vec3 d = m.pos[m.vertexAt(f, k)] - m.pos[m.vertexAt(f, 0)];
        px[k] = dot(d, u);
        py[k] = dot(d, v);
        ids[k] = k;
    }
    auto turn = [&](int a, int b, int c) { return (px[b] - px[a]) * (py[c] - py[a]) - (py[b] - py[a]) * (px[c] - px[a]); };
    auto inside = [&](int a, int b, int c, int q) {  // q strictly inside triangle abc (ccw)
        return turn(a, b, q) > 0 && turn(b, c, q) > 0 && turn(c, a, q) > 0;
    };
    std::vector<std::array<int, 3>> tris;
    std::vector<int> ring = ids;
    while (ring.size() > 3) {
        bool clipped = false;
        for (size_t i = 0; i < ring.size() && !clipped; ++i) {
            const int a = ring[(i + ring.size() - 1) % ring.size()], b = ring[i], c = ring[(i + 1) % ring.size()];
            if (turn(a, b, c) <= 0) continue;
            bool empty = true;
            for (int q : ring)
                if (q != a && q != b && q != c && inside(a, b, c, q)) empty = false;
            if (!empty) continue;
            tris.push_back({m.vertexAt(f, a), m.vertexAt(f, b), m.vertexAt(f, c)});
            ring.erase(ring.begin() + static_cast<long>(i));
            clipped = true;
        }
        if (!clipped) return false;
    }
    if (turn(ring[0], ring[1], ring[2]) <= 0) return false;
    tris.push_back({m.vertexAt(f, ring[0]), m.vertexAt(f, ring[1]), m.vertexAt(f, ring[2])});
    out.insert(out.end(), tris.begin(), tris.end());
    return static_cast<int>(tris.size()) == n - 2;
}

// Every face's validated triangles, per face; faces without a valid triangulation are counted and left out.
int triangulateFaces(const PolygonMesh& m, std::vector<std::vector<std::array<int, 3>>>& perFace, std::vector<int>& unvalidatedFaces) {
    int unvalidated = 0;
    perFace.assign(m.faces(), {});
    for (int f = 0; f < m.faces(); ++f)
        if (!triangulateFace(m, f, perFace[f])) {
            perFace[f].clear();
            ++unvalidated;
            if (unvalidatedFaces.size() < 64) unvalidatedFaces.push_back(f);
        }
    return unvalidated;
}

// The faces of a FAILed test that intersect, by face bounding boxes and the triangle test on the validated triangles
// (pairs sharing a vertex are adjacent, not tested), for inspection (at most 64 pairs).
std::vector<std::pair<int, int>> intersectingFacePairs(const PolygonMesh& m, const std::vector<std::vector<std::array<int, 3>>>& tris) {
    std::vector<std::pair<int, int>> out;
    const int nF = m.faces();
    std::vector<Vec3> lo(nF), hi(nF);
    double extent = 0;
    for (int f = 0; f < nF; ++f) {
        lo[f] = hi[f] = m.pos[m.vertexAt(f, 0)];
        for (int k = 1; k < m.sides(f); ++k) {
            const Vec3& p = m.pos[m.vertexAt(f, k)];
            lo[f] = {std::min(lo[f].x, p.x), std::min(lo[f].y, p.y), std::min(lo[f].z, p.z)};
            hi[f] = {std::max(hi[f].x, p.x), std::max(hi[f].y, p.y), std::max(hi[f].z, p.z)};
        }
        extent = std::max(extent, norm(hi[f] - lo[f]));
    }
    const double tol = 1e-7 * extent;
    std::vector<int> order(nF);
    for (int f = 0; f < nF; ++f) order[f] = f;
    std::sort(order.begin(), order.end(), [&](int a, int b) { return lo[a].x < lo[b].x; });
    for (int i = 0; i < nF && out.size() < 64; ++i) {
        const int f = order[i];
        for (int j = i + 1; j < nF && lo[order[j]].x <= hi[f].x && out.size() < 64; ++j) {
            const int g = order[j];
            if (lo[g].y > hi[f].y || hi[g].y < lo[f].y || lo[g].z > hi[f].z || hi[g].z < lo[f].z) continue;
            bool hit = false;
            for (const std::array<int, 3>& ta : tris[f])
                for (const std::array<int, 3>& tb : tris[g]) {
                    if (hit) break;
                    bool share = false;
                    for (int a : ta)
                        for (int b : tb)
                            if (a == b) share = true;
                    if (share) continue;
                    const Vec3 pa[3] = {m.pos[ta[0]], m.pos[ta[1]], m.pos[ta[2]]};
                    const Vec3 pb[3] = {m.pos[tb[0]], m.pos[tb[1]], m.pos[tb[2]]};
                    if (triangleIntersection(pa, pb, tol) == Status::Fail) hit = true;
                }
            if (hit) out.push_back({std::min(f, g), std::max(f, g)});
        }
    }
    return out;
}

// The validated triangles of every face through the triangle mesh's self-intersection test. A face without a valid
// triangulation leaves the verdict UNRESOLVED (not certified), never FAIL.
void intersections(const PolygonMesh& m, const std::vector<std::vector<std::array<int, 3>>>& tris, int unvalidated, Status& status, int& pairs, int& unresolved) {
    TriMesh tri;
    tri.pos = m.pos;
    for (const std::vector<std::array<int, 3>>& face : tris)
        for (const std::array<int, 3>& t : face) tri.tri.push_back(t);
    const MeshCheck check = checkMesh(tri, SelfIntersectionTest::Run);
    status = unvalidated ? Status::Unresolved : Status::Pass;
    pairs = unresolved = 0;
    for (const ComponentCheck& c : check.components) {
        status = combine(status, c.selfIntersection);
        pairs += c.intersectingPairs;
        unresolved += c.unresolvedPairs;
    }
}

// --- the validator ---------------------------------------------------------------------------------------------------

StitchMeshCheck validateStitchMesh(const StitchMesh& s, const PolygonMesh& coarse) {
    StitchMeshCheck c;
    const PolygonMesh& m = s.mesh;
    c.mesh = checkPolygonMesh(m);
    c.faces = m.faces();
    const PolygonMeshCheck before = checkPolygonMesh(coarse);
    c.eulerBefore = before.euler;
    c.eulerAfter = c.mesh.euler;
    c.boundaryLoopsBefore = before.boundaryLoops;
    c.boundaryLoopsAfter = c.mesh.boundaryLoops;
    std::string detail;
    if (c.mesh.status != Status::Pass) detail += "polygon mesh: " + c.mesh.detail + "; ";
    if (s.role.size() != static_cast<size_t>(m.faces()) || s.kind.size() != static_cast<size_t>(m.faces())) {
        c.status = Status::Fail;
        c.detail = detail + "roles or kinds missing for some faces";
        return c;
    }
    HalfEdgeMesh half;
    std::string why;
    const bool haveHalf = buildHalfEdges(m, half, &why);
    // faces: two wale edges, the kind the roles support
    for (int f = 0; f < m.faces(); ++f) {
        if (static_cast<int>(s.role[f].size()) != m.sides(f)) { ++c.kindMismatches; continue; }
        const StitchFaceKind k = kindOfRoles(s.role[f]);
        if (k != s.kind[f]) ++c.kindMismatches;
        if (waleCorners(s.role[f]).size() != 2) ++c.facesWithoutTwoWaleEdges;
        switch (s.kind[f]) {
        case StitchFaceKind::Regular: ++c.regular; break;
        case StitchFaceKind::Increase: ++c.increase; break;
        case StitchFaceKind::Decrease: ++c.decrease; break;
        case StitchFaceKind::SpecialQuad: ++c.special; break;
        default: ++c.unsupported; break;
        }
    }
    if (haveHalf) {
        // the surface: one fan per vertex, faces wound with their coarse face, no two faces intersecting
        c.brokenVertexLinks = brokenVertexLinks(m, half);
        faceGeometry(s, coarse, c.invertedFaces, c.nonConvexFaces, c.invertedFaceIds);
        std::vector<std::vector<std::array<int, 3>>> tris;
        c.unvalidatedTriangulations = triangulateFaces(m, tris, c.unvalidatedFaces);
        intersections(m, tris, c.unvalidatedTriangulations, c.intersections, c.intersectingPairs, c.unresolvedPairs);
        if (c.intersections == Status::Fail) c.intersectingFaces = intersectingFacePairs(m, tris);
        // two faces sharing two edges (a degree-two interior vertex between them): listed for the yarn stage
        std::vector<int> valence(m.pos.size(), 0);
        std::map<EdgeKey, int> edgeFaces;
        for (int f = 0; f < m.faces(); ++f)
            for (int k = 0; k < m.sides(f); ++k) ++edgeFaces[edgeKey(m.vertexAt(f, k), m.vertexAt(f, (k + 1) % m.sides(f)))];
        for (const auto& [e, n] : edgeFaces) { ++valence[e.first]; ++valence[e.second]; }
        for (size_t v = 0; v < m.pos.size(); ++v)
            if (valence[v] == 2 && !half.boundaryVertex[v]) ++c.sharedEdgePairs;
        // every wale edge traced exactly once by the strands; every supported face in exactly one strand
        std::map<EdgeKey, int> waleUse;    // interior wale edge -> times a strand crosses it
        std::vector<int> inStrands(m.faces(), 0);
        c.strands = static_cast<int>(s.strands.size());
        auto crossings = [&](int f, int g) {  // the wale edges shared by two consecutive faces of a strand
            int n = 0;
            for (int k : waleCorners(s.role[f])) {
                const int twin = half.twin[m.faceStart[f] + k];
                if (twin >= 0 && half.face[twin] == g) {
                    ++waleUse[edgeKey(m.vertexAt(f, k), m.vertexAt(f, (k + 1) % m.sides(f)))];
                    ++n;
                }
            }
            return n;
        };
        for (const Strand& st : s.strands) {
            if (st.closed) ++c.closedStrands;
            else ++c.openStrands;
            bool valid = true;
            for (int f : st.faces)
                if (f < 0 || f >= m.faces()) valid = false;
            if (!valid) { ++c.waleEdgesNotTracedOnce; continue; }  // a strand naming a face the mesh does not have
            for (size_t i = 0; i < st.faces.size(); ++i) {
                ++inStrands[st.faces[i]];
                if (i + 1 < st.faces.size()) crossings(st.faces[i], st.faces[i + 1]);
            }
            if (st.closed && st.faces.size() > 1) crossings(st.faces.back(), st.faces.front());
            if (st.closed && st.faces.size() == 1) ++waleUse[EdgeKey{-1, -1}];  // one face closing on itself: both its wale edges are one edge
        }
        for (int f = 0; f < m.faces(); ++f) {
            const bool supported = waleCorners(s.role[f]).size() == 2 && s.kind[f] != StitchFaceKind::Unsupported;
            if (supported && inStrands[f] != 1) ++c.waleEdgesNotTracedOnce;
            if (!supported && inStrands[f] != 0) ++c.waleEdgesNotTracedOnce;
            for (int k : waleCorners(s.role[f])) {
                const int twin = half.twin[m.faceStart[f] + k];
                if (twin < 0) { ++c.yarnEnds; continue; }
                const int g = half.face[twin];
                if (g < f) continue;
                const bool gSupported = waleCorners(s.role[g]).size() == 2 && s.kind[g] != StitchFaceKind::Unsupported;
                const int used = waleUse.count(edgeKey(m.vertexAt(f, k), m.vertexAt(f, (k + 1) % m.sides(f)))) ?
                                 waleUse.at(edgeKey(m.vertexAt(f, k), m.vertexAt(f, (k + 1) % m.sides(f)))) : 0;
                if (supported && gSupported) {
                    if (used == 1) ++c.waleEdgesTracedOnce;
                    else ++c.waleEdgesNotTracedOnce;
                } else if (used != 0) ++c.waleEdgesNotTracedOnce;
            }
        }
        // loops: the classification from the roles, compared with the builder's list
        const std::vector<LoopPort> loops = classifyLoops(half, m, s.role);
        std::map<EdgeKey, LoopPort::Kind> listed;
        for (const LoopPort& p : s.loops) listed[edgeKey(p.a, p.b)] = p.kind;
        for (const LoopPort& p : loops) {
            switch (p.kind) {
            case LoopPort::Kind::Interlock: ++c.interlocks; break;
            case LoopPort::Kind::BothTop: ++c.bothTop; break;
            case LoopPort::Kind::BothBottom: ++c.bothBottom; break;
            case LoopPort::Kind::BoundaryBottom: ++c.boundaryBottom; break;
            default: ++c.boundaryTop; break;
            }
            const auto it = listed.find(edgeKey(p.a, p.b));
            if (it == listed.end() || it->second != p.kind) ++c.loopsNotClassified;
        }
        if (loops.size() != s.loops.size()) c.loopsNotClassified += static_cast<int>(std::abs(static_cast<long>(loops.size()) - static_cast<long>(s.loops.size())));
    } else {
        detail += "no half-edge structure: " + why + "; ";
    }
    if (c.brokenVertexLinks) detail += std::to_string(c.brokenVertexLinks) + " vertices whose faces form more than one fan; ";
    if (c.invertedFaces) detail += std::to_string(c.invertedFaces) + " faces wound against their coarse face or of no area; ";
    if (c.intersections == Status::Fail) detail += std::to_string(c.intersectingPairs) + " pairs of faces intersect; ";
    if (c.unvalidatedTriangulations) detail += std::to_string(c.unvalidatedTriangulations) + " faces without a valid triangulation in their plane (intersections not certified); ";
    if (c.kindMismatches) detail += std::to_string(c.kindMismatches) + " faces whose kind the roles do not support; ";
    if (c.waleEdgesNotTracedOnce) detail += std::to_string(c.waleEdgesNotTracedOnce) + " wale edges or faces not traced exactly once; ";
    if (c.loopsNotClassified) detail += std::to_string(c.loopsNotClassified) + " course edges not classified as the roles say; ";
    if (c.eulerBefore != c.eulerAfter || c.boundaryLoopsBefore != c.boundaryLoopsAfter)
        detail += "topology changed: Euler " + std::to_string(c.eulerBefore) + " -> " + std::to_string(c.eulerAfter) + ", boundary loops " +
                  std::to_string(c.boundaryLoopsBefore) + " -> " + std::to_string(c.boundaryLoopsAfter) + "; ";
    const bool surfaceFault = c.mesh.status != Status::Pass || !haveHalf || c.brokenVertexLinks || c.invertedFaces ||
                              c.intersections == Status::Fail || c.eulerBefore != c.eulerAfter || c.boundaryLoopsBefore != c.boundaryLoopsAfter;
    c.surface = surfaceFault ? Status::Fail : c.intersections == Status::Unresolved ? Status::Unresolved : Status::Pass;
    const bool recordFault = c.kindMismatches || c.waleEdgesNotTracedOnce || c.loopsNotClassified;
    if (surfaceFault || recordFault) c.status = Status::Fail;
    else if (c.unsupported || c.bothTop || c.bothBottom || c.sharedEdgePairs || c.surface == Status::Unresolved) {
        c.status = Status::Unresolved;
        detail += std::to_string(c.unsupported) + " unsupported faces, " + std::to_string(c.bothTop) + " course edges offered by both sides, " +
                  std::to_string(c.bothBottom) + " received by both, " + std::to_string(c.sharedEdgePairs) +
                  " face pairs sharing two edges (yarn support open)";
        if (c.surface == Status::Unresolved) detail += "; " + std::to_string(c.unresolvedPairs) + " face pairs within tolerance of touching";
    } else {
        c.status = Status::Pass;
        detail = "every face typed, every wale edge traced once, every course edge an interlock or a declared end";
    }
    c.detail = detail;
    return c;
}

std::string StitchMeshCheck::text() const {
    std::string t = "stitch mesh check " + std::string(toString(status)) + ": " + detail + "\n";
    t += "    faces " + std::to_string(faces) + " (regular " + std::to_string(regular) + ", increase " + std::to_string(increase) + ", decrease " +
         std::to_string(decrease) + ", special " + std::to_string(special) + ", unsupported " + std::to_string(unsupported) + ")\n";
    t += "    strands " + std::to_string(strands) + " (closed " + std::to_string(closedStrands) + ", open " + std::to_string(openStrands) +
         "), yarn ends " + std::to_string(yarnEnds) + ", wale edges traced once " + std::to_string(waleEdgesTracedOnce) + "\n";
    t += "    loops: interlocks " + std::to_string(interlocks) + ", both top " + std::to_string(bothTop) + ", both bottom " + std::to_string(bothBottom) +
         ", boundary bottom " + std::to_string(boundaryBottom) + ", boundary top " + std::to_string(boundaryTop) + "\n";
    t += "    surface " + std::string(toString(surface)) + ": fans broken " + std::to_string(brokenVertexLinks) + ", inverted faces " +
         std::to_string(invertedFaces) + ", non-convex faces " + std::to_string(nonConvexFaces) + ", intersections " + toString(intersections) +
         " (" + std::to_string(intersectingPairs) + " pairs, " + std::to_string(unresolvedPairs) + " unresolved, " + std::to_string(unvalidatedTriangulations) + " faces not triangulated); Euler " + std::to_string(eulerBefore) +
         " -> " + std::to_string(eulerAfter) + ", boundary loops " + std::to_string(boundaryLoopsBefore) + " -> " + std::to_string(boundaryLoopsAfter) + "\n";
    t += "    face pairs sharing two edges " + std::to_string(sharedEdgePairs) + "\n";
    return t;
}

} // namespace knit
