#include "stitchmesh/typed_mesh.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <utility>

#include "model/hash.h"
#include "model/mesh_check.h"
#include "model/mesh_query.h"

namespace knit {

const char* toString(EdgeBoundary b) {
    switch (b) {
    case EdgeBoundary::Interior: return "interior";
    case EdgeBoundary::CastOn: return "cast-on";
    case EdgeBoundary::BindOff: return "bind-off";
    case EdgeBoundary::Side: return "side";
    }
    return "?";
}

namespace {

using EdgeKey = std::pair<int, int>;
EdgeKey keyOf(int a, int b) { return {std::min(a, b), std::max(a, b)}; }

struct CornerUse {
    int face = -1, corner = -1;
};

// Every undirected edge with the face corners that use it.
std::map<EdgeKey, std::vector<CornerUse>> edgeUses(const PolygonMesh& m) {
    std::map<EdgeKey, std::vector<CornerUse>> uses;
    for (int f = 0; f < m.faces(); ++f) {
        const int n = m.sides(f);
        for (int k = 0; k < n; ++k) uses[keyOf(m.vertexAt(f, k), m.vertexAt(f, (k + 1) % n))].push_back(CornerUse{f, k});
    }
    return uses;
}

// The other use of the edge at (f, k): {-1, -1} on the border or at a non-manifold edge.
CornerUse across(const std::map<EdgeKey, std::vector<CornerUse>>& uses, const PolygonMesh& m, int f, int k) {
    const int n = m.sides(f);
    const auto it = uses.find(keyOf(m.vertexAt(f, k), m.vertexAt(f, (k + 1) % n)));
    if (it == uses.end() || it->second.size() != 2) return CornerUse{};
    const CornerUse& a = it->second[0];
    const CornerUse& b = it->second[1];
    return (a.face == f && a.corner == k) ? b : a;
}

double pointPolylineDistance(const Vec3& p, const std::vector<Vec3>& loop) {
    double best = 1e300;
    for (size_t i = 0; i < loop.size(); ++i) {
        const ClosestSegmentPoints c = closestSegmentPoints(p, p, loop[i], loop[(i + 1) % loop.size()]);
        best = std::min(best, c.distance);
    }
    return best;
}

// Interior angle at corner c of face f in degrees, in [0, 360), measured with the face's Newell normal.
double interiorAngle(const PolygonMesh& m, int f, int c, const Vec3& n) {
    const int s = m.sides(f);
    const Vec3 P = m.pos[m.vertexAt(f, c)];
    const Vec3 a = m.pos[m.vertexAt(f, (c + s - 1) % s)] - P;
    const Vec3 b = m.pos[m.vertexAt(f, (c + 1) % s)] - P;
    double ang = std::atan2(dot(n, cross(b, a)), dot(b, a)) * 180.0 / kPi;
    if (ang < 0) ang += 360.0;
    return ang;
}

double edgeLength(const PolygonMesh& m, int f, int k) {
    const int s = m.sides(f);
    return norm(m.pos[m.vertexAt(f, (k + 1) % s)] - m.pos[m.vertexAt(f, k)]);
}

} // namespace

StitchFaceKind kindFromRoles(const std::vector<StitchEdgeRole>& r) {
    const int n = static_cast<int>(r.size());
    int wale = 0, bottom = 0, top = 0;
    bool waleAdjacent = false, bottomAdjacent = false, topAdjacent = false;
    for (int k = 0; k < n; ++k) {
        const StitchEdgeRole next = r[(k + 1) % n];
        if (r[k] == StitchEdgeRole::Wale) ++wale;
        else if (r[k] == StitchEdgeRole::Bottom) ++bottom;
        else ++top;
        if (r[k] == StitchEdgeRole::Wale && next == StitchEdgeRole::Wale) waleAdjacent = true;
        if (r[k] == StitchEdgeRole::Bottom && next == StitchEdgeRole::Bottom) bottomAdjacent = true;
        if (r[k] == StitchEdgeRole::Top && next == StitchEdgeRole::Top) topAdjacent = true;
    }
    if (wale != 2) return StitchFaceKind::Unsupported;
    if (n == 3 && top == 1 && bottom == 0) return StitchFaceKind::PoleStart;
    if (n == 4 && bottom == 1 && top == 1) return waleAdjacent ? StitchFaceKind::SpecialQuad : StitchFaceKind::Regular;
    if (n == 5 && bottom == 1 && top == 2 && topAdjacent && !waleAdjacent) return StitchFaceKind::Increase;
    if (n == 5 && bottom == 2 && top == 1 && bottomAdjacent && !waleAdjacent) return StitchFaceKind::Decrease;
    return StitchFaceKind::Unsupported;
}

int exitWaleCorner(const std::vector<StitchEdgeRole>& r) {
    const StitchFaceKind kind = kindFromRoles(r);
    if (kind == StitchFaceKind::SpecialQuad || kind == StitchFaceKind::Unsupported) return -1;
    const int n = static_cast<int>(r.size());
    if (kind == StitchFaceKind::PoleStart) {
        for (int k = 0; k < n; ++k)
            if (r[k] == StitchEdgeRole::Wale && r[(k + 1) % n] == StitchEdgeRole::Top) return k;
        return -1;
    }
    for (int k = 0; k < n; ++k)
        if (r[k] == StitchEdgeRole::Wale && r[(k + n - 1) % n] == StitchEdgeRole::Bottom) return k;
    return -1;
}

int entryWaleCorner(const std::vector<StitchEdgeRole>& r) {
    const StitchFaceKind kind = kindFromRoles(r);
    if (kind == StitchFaceKind::SpecialQuad || kind == StitchFaceKind::Unsupported) return -1;
    const int n = static_cast<int>(r.size());
    if (kind == StitchFaceKind::PoleStart) {
        for (int k = 0; k < n; ++k)
            if (r[k] == StitchEdgeRole::Wale && r[(k + n - 1) % n] == StitchEdgeRole::Top) return k;
        return -1;
    }
    for (int k = 0; k < n; ++k)
        if (r[k] == StitchEdgeRole::Wale && r[(k + 1) % n] == StitchEdgeRole::Bottom) return k;
    return -1;
}

StrandTrace traceYarnStrands(const PolygonMesh& m, const std::vector<std::vector<StitchEdgeRole>>& role) {
    StrandTrace t;
    const int nF = m.faces();
    const std::map<EdgeKey, std::vector<CornerUse>> uses = edgeUses(m);
    std::vector<std::vector<int>> wale(nF);
    std::vector<int> exitOf(nF, -1);
    for (int f = 0; f < nF; ++f) {
        for (int k = 0; k < static_cast<int>(role[f].size()); ++k)
            if (role[f][k] == StitchEdgeRole::Wale) wale[f].push_back(k);
        exitOf[f] = exitWaleCorner(role[f]);
    }
    auto supported = [&](int f) { return f >= 0 && wale[f].size() == 2; };
    auto otherWale = [&](int f, int corner) { return wale[f][0] == corner ? wale[f][1] : wale[f][0]; };
    std::map<EdgeKey, int> crossings;
    for (int f = 0; f < nF; ++f)
        for (int k : wale[f])
            if (across(uses, m, f, k).face < 0) ++t.yarnEnds;
    t.strandOfFace.assign(nF, -1);
    std::vector<char> visited(nF, 0);
    for (int f0 = 0; f0 < nF; ++f0) {
        if (visited[f0] || !supported(f0)) continue;
        int start = f0, entry = wale[f0][0];
        bool closed = false;
        int cur = f0, via = wale[f0][0];
        for (int step = 0; step <= nF; ++step) {   // walk back to the strand's start
            const CornerUse n = across(uses, m, cur, via);
            if (n.face < 0 || !supported(n.face)) {
                start = cur;
                entry = via;
                break;
            }
            if (n.face == f0) {
                closed = true;
                start = f0;
                entry = exitOf[f0] >= 0 ? otherWale(f0, exitOf[f0]) : wale[f0][0];
                break;
            }
            cur = n.face;
            via = otherWale(cur, n.corner);
        }
        std::vector<int> strand;
        const int id = static_cast<int>(t.strands.size());
        cur = start;
        int exit = otherWale(start, entry);
        for (int step = 0; step <= nF; ++step) {
            strand.push_back(cur);
            visited[cur] = 1;
            t.strandOfFace[cur] = id;
            const CornerUse n = across(uses, m, cur, exit);
            const int s = m.sides(cur);
            if (n.face >= 0) crossings[keyOf(m.vertexAt(cur, exit), m.vertexAt(cur, (exit + 1) % s))] += 1;
            if (n.face < 0 || !supported(n.face) || n.face == start || visited[n.face]) break;
            cur = n.face;
            exit = otherWale(cur, n.corner);
        }
        t.strands.push_back(strand);
        t.closed.push_back(closed ? 1 : 0);
    }
    for (const auto& [key, list] : uses) {
        if (list.size() != 2) continue;
        const CornerUse& a = list[0];
        const CornerUse& b = list[1];
        if (role[a.face][a.corner] != StitchEdgeRole::Wale || role[b.face][b.corner] != StitchEdgeRole::Wale) continue;
        if (!supported(a.face) || !supported(b.face)) continue;
        const auto c = crossings.find(key);
        if (c != crossings.end() && c->second == 1) ++t.waleEdgesTracedOnce;
        else ++t.waleEdgesNotTracedOnce;
    }
    return t;
}

bool earClipFace(const PolygonMesh& m, int f, std::vector<std::array<int, 3>>& out) {
    out.clear();
    const int n = m.sides(f);
    std::vector<int> idx;
    for (int k = 0; k < n; ++k) idx.push_back(m.vertexAt(f, k));
    if (n == 3) {
        out.push_back({idx[0], idx[1], idx[2]});
        return true;
    }
    const Vec3 nrm = m.faceNormal(f);
    Vec3 u = cross(nrm, Vec3{1, 0, 0});
    if (norm(u) < 0.5) u = cross(nrm, Vec3{0, 1, 0});
    u = normalized(u);
    const Vec3 v = cross(nrm, u);
    std::vector<double> px, py;
    for (int k = 0; k < n; ++k) {
        px.push_back(dot(m.pos[idx[k]], u));
        py.push_back(dot(m.pos[idx[k]], v));
    }
    std::vector<int> ring;
    for (int k = 0; k < n; ++k) ring.push_back(k);
    auto cross2 = [&](int a, int b, int c) { return (px[b] - px[a]) * (py[c] - py[a]) - (py[b] - py[a]) * (px[c] - px[a]); };
    while (ring.size() > 3) {
        bool clipped = false;
        const int r = static_cast<int>(ring.size());
        for (int i = 0; i < r && !clipped; ++i) {
            const int a = ring[(i + r - 1) % r], b = ring[i], c = ring[(i + 1) % r];
            if (cross2(a, b, c) <= 1e-12) continue;   // reflex or flat corner
            bool inside = false;
            for (int j = 0; j < r && !inside; ++j) {
                const int q = ring[j];
                if (q == a || q == b || q == c) continue;
                inside = cross2(a, b, q) >= 0 && cross2(b, c, q) >= 0 && cross2(c, a, q) >= 0;
            }
            if (inside) continue;
            out.push_back({idx[a], idx[b], idx[c]});
            ring.erase(ring.begin() + i);
            clipped = true;
        }
        if (!clipped) {
            out.clear();
            for (int k = 1; k + 1 < n; ++k) out.push_back({idx[0], idx[k], idx[k + 1]});
            return false;
        }
    }
    out.push_back({idx[ring[0]], idx[ring[1]], idx[ring[2]]});
    return true;
}

TypedMeshCheck checkTypedMesh(const TypedStitchMesh& s, const TriMesh& reference, const std::vector<DeclaredBoundary>& declared,
                              const TypedCheckOptions& o) {
    TypedMeshCheck c;
    const PolygonMesh& m = s.mesh;
    const int nF = m.faces();
    c.faces = nF;
    std::string why;
    if (static_cast<int>(s.kind.size()) != nF || static_cast<int>(s.role.size()) != nF || static_cast<int>(s.boundary.size()) != nF) {
        c.connectivity = Status::Fail;
        c.geometry = Status::Unresolved;
        c.detail = "record sizes disagree with the faces";
        return c;
    }
    c.polygon = checkPolygonMesh(m);
    // kinds from the roles alone
    std::map<std::string, int> reasons;
    for (int f = 0; f < nF; ++f) {
        if (static_cast<int>(s.role[f].size()) != m.sides(f) || static_cast<int>(s.boundary[f].size()) != m.sides(f)) {
            c.connectivity = Status::Fail;
            c.detail = "roles or boundary roles do not match a face's corners";
            return c;
        }
        const StitchFaceKind derived = kindFromRoles(s.role[f]);
        int wale = 0;
        for (StitchEdgeRole r : s.role[f]) wale += r == StitchEdgeRole::Wale ? 1 : 0;
        if (wale != 2) ++c.facesWithoutTwoWaleEdges;
        switch (s.kind[f]) {
        case StitchFaceKind::Regular: ++c.regular; break;
        case StitchFaceKind::Increase: ++c.increase; break;
        case StitchFaceKind::Decrease: ++c.decrease; break;
        case StitchFaceKind::SpecialQuad: ++c.special; break;
        case StitchFaceKind::PoleStart: ++c.poleStart; break;
        case StitchFaceKind::Unsupported:
            ++c.unsupported;
            reasons[f < static_cast<int>(s.reason.size()) && !s.reason[f].empty() ? s.reason[f] : std::string("unsupported-roles")] += 1;
            break;
        }
        // a face stored Unsupported (a slot the lookup refused) may carry valid-looking roles; any other kind must agree
        if (s.kind[f] != StitchFaceKind::Unsupported && derived != s.kind[f]) ++c.kindMismatches;
    }
    for (const auto& [r, n] : reasons) c.unsupportedReasons.push_back({r, n});
    // edges
    const std::map<EdgeKey, std::vector<CornerUse>> uses = edgeUses(m);
    for (const auto& [key, list] : uses) {
        if (list.size() == 2) {
            const CornerUse& a = list[0];
            const CornerUse& b = list[1];
            const StitchEdgeRole ra = s.role[a.face][a.corner], rb = s.role[b.face][b.corner];
            if (s.boundary[a.face][a.corner] != EdgeBoundary::Interior) ++c.boundaryRoleOnInteriorEdge;
            if (s.boundary[b.face][b.corner] != EdgeBoundary::Interior) ++c.boundaryRoleOnInteriorEdge;
            const bool courseA = ra != StitchEdgeRole::Wale, courseB = rb != StitchEdgeRole::Wale;
            if (courseA && courseB) {
                if (ra != rb) ++c.interlocks;
                else if (ra == StitchEdgeRole::Top) ++c.bothTop;
                else ++c.bothBottom;
            } else if (courseA != courseB) {
                ++c.courseAgainstWale;
            } else {
                const int exA = exitWaleCorner(s.role[a.face]), enA = entryWaleCorner(s.role[a.face]);
                const int exB = exitWaleCorner(s.role[b.face]), enB = entryWaleCorner(s.role[b.face]);
                if (exA < 0 || exB < 0) continue;   // a special quad or unsupported face: no own course direction
                const bool agree = (a.corner == exA && b.corner == enB) || (a.corner == enA && b.corner == exB);
                if (agree) ++c.waleAgree; else ++c.waleDisagree;
            }
        } else if (list.size() == 1) {
            const CornerUse& a = list[0];
            const StitchEdgeRole r = s.role[a.face][a.corner];
            const EdgeBoundary b = s.boundary[a.face][a.corner];
            if (r == StitchEdgeRole::Wale) {
                if (b == EdgeBoundary::Side) ++c.sideEdges; else ++c.borderWaleEdges;
                continue;
            }
            if (b == EdgeBoundary::Interior) {
                ++c.borderCourseWithoutRole;
                continue;
            }
            if (b == EdgeBoundary::Side || (b == EdgeBoundary::CastOn && r != StitchEdgeRole::Bottom) || (b == EdgeBoundary::BindOff && r != StitchEdgeRole::Top)) {
                ++c.boundaryRoleWrongSide;
                continue;
            }
            if (b == EdgeBoundary::CastOn) ++c.castOnEdges; else ++c.bindOffEdges;
            const Vec3 p = m.pos[key.first], q = m.pos[key.second];
            bool onDeclared = false;
            for (const DeclaredBoundary& d : declared)
                if (d.role == b && pointPolylineDistance(p, d.loop) <= o.boundaryTolerance && pointPolylineDistance(q, d.loop) <= o.boundaryTolerance) onDeclared = true;
            if (!onDeclared) ++c.boundaryRoleOffDeclared;
        }
    }
    c.strands = traceYarnStrands(m, s.role);
    for (char cl : c.strands.closed) (cl ? c.closedStrands : c.openStrands) += 1;
    // the count invariant per strip, from the faces alone, against the rings of the record
    const bool hasStrips = !s.strip.empty();
    if (o.requireStrips && (!hasStrips || static_cast<int>(s.strip.size()) != nF || s.rings.size() < 2)) {
        c.stripsBroken = 1;
    } else if (hasStrips) {
        int strips = 0;
        for (int f = 0; f < nF; ++f) strips = std::max(strips, s.strip[f] + 1);
        c.strips.resize(strips);
        for (int k = 0; k < strips; ++k) {
            c.strips[k].strip = k;
            if (k + 1 < static_cast<int>(s.rings.size())) {
                c.strips[k].ringBelow = static_cast<int>(s.rings[k].size());
                c.strips[k].ringAbove = static_cast<int>(s.rings[k + 1].size());
            }
        }
        for (int f = 0; f < nF; ++f) {
            if (s.strip[f] < 0) {
                ++c.stripsBroken;
                continue;
            }
            StripCount& sc = c.strips[s.strip[f]];
            for (StitchEdgeRole r : s.role[f]) {
                if (r == StitchEdgeRole::Bottom) ++sc.bottomEdges;
                if (r == StitchEdgeRole::Top) ++sc.topEdges;
            }
            const StitchFaceKind derived = kindFromRoles(s.role[f]);
            if (derived == StitchFaceKind::Increase) ++sc.increases;
            if (derived == StitchFaceKind::Decrease) ++sc.decreases;
            if (derived == StitchFaceKind::PoleStart && s.kind[f] == StitchFaceKind::PoleStart) ++sc.poleStarts;
            if (s.kind[f] == StitchFaceKind::Unsupported) ++sc.unsupported;
        }
        std::vector<int> facesInStrip(c.strips.size(), 0);
        for (int f = 0; f < nF; ++f)
            if (s.strip[f] >= 0) facesInStrip[s.strip[f]] += 1;
        for (StripCount& sc : c.strips) {
            if (sc.ringBelow == 1)   // a pole strip: the ring below is the pole vertex
                sc.invariant = sc.bottomEdges == 0 && sc.topEdges == sc.ringAbove && sc.poleStarts == facesInStrip[sc.strip] && sc.poleStarts == sc.ringAbove;
            else
                sc.invariant = sc.bottomEdges == sc.ringBelow && sc.topEdges == sc.ringAbove && sc.topEdges - sc.bottomEdges == sc.increases - sc.decreases;
            if (!sc.invariant) ++c.stripsBroken;
            if (sc.increases > 0 && sc.decreases > 0) ++c.mixedShapingStrips;
        }
    }
    // the connectivity verdict
    std::string fails, opens;
    auto add = [](std::string& list, const char* what, int n) {
        if (n <= 0) return;
        char b[160];
        std::snprintf(b, sizeof b, "%s%s %d", list.empty() ? "" : "; ", what, n);
        list += b;
    };
    if (c.polygon.status != Status::Pass) fails += "polygon mesh: " + c.polygon.detail;
    add(fails, "kind mismatches", c.kindMismatches);
    add(fails, "boundary role on an interior edge", c.boundaryRoleOnInteriorEdge);
    add(fails, "boundary role on the wrong side", c.boundaryRoleWrongSide);
    add(fails, "boundary role away from its declared loop", c.boundaryRoleOffDeclared);
    add(fails, "strips breaking the count invariant", c.stripsBroken);
    add(opens, "unsupported faces", c.unsupported);
    add(opens, "both-top course edges", c.bothTop);
    add(opens, "both-bottom course edges", c.bothBottom);
    add(opens, "course edges against wale edges", c.courseAgainstWale);
    add(opens, "wale edges whose yarn directions disagree", c.waleDisagree);
    add(opens, "border course edges without a boundary role", c.borderCourseWithoutRole);
    add(opens, "yarn ends on undeclared border wale edges", c.borderWaleEdges);
    add(opens, "strips with both INC and DEC", c.mixedShapingStrips);
    c.connectivity = !fails.empty() ? Status::Fail : !opens.empty() ? Status::Unresolved : Status::Pass;
    // geometry
    const MeshQuery query(reference);
    TriMesh fan;
    fan.pos = m.pos;
    c.perFace.assign(nF, FaceCheck{});
    for (int f = 0; f < nF; ++f) {
        FaceCheck& fc = c.perFace[f];
        const Vec3 n = m.faceNormal(f);
        const MeshQuery::Hit hit = query.closest(m.faceCentroid(f));
        fc.inverted = hit.tri < 0 || dot(n, reference.faceNormal(hit.tri)) <= 0;
        fc.degenerate = m.faceArea(f) < o.degenerateArea * o.w * o.h;
        c.invertedFaces += fc.inverted;
        c.degenerateFaces += fc.degenerate;
        double dev = hit.tri >= 0 ? hit.dist : 1e300;
        for (int k = 0; k < m.sides(f); ++k) {
            const MeshQuery::Hit e = query.closest(0.5 * (m.pos[m.vertexAt(f, k)] + m.pos[m.vertexAt(f, (k + 1) % m.sides(f))]));
            dev = std::max(dev, e.tri >= 0 ? e.dist : 1e300);
        }
        c.maxFaceDeviation = std::max(c.maxFaceDeviation, dev);
        fc.deviation = static_cast<float>(std::min(dev, 1e30));
        fc.overDeviation = dev > o.faceDeviation;
        c.facesOverDeviation += fc.overDeviation;
        std::vector<std::array<int, 3>> tris;
        fc.untriangulable = !earClipFace(m, f, tris);
        c.untriangulableFaces += fc.untriangulable;
        for (const std::array<int, 3>& t : tris) fan.tri.push_back(t);
    }
    const MeshCheck mc = checkMesh(fan);
    c.intersections = Status::Pass;
    for (const ComponentCheck& cc : mc.components) {
        c.intersectingPairs += cc.intersectingPairs;
        if (cc.selfIntersection == Status::Fail) c.intersections = Status::Fail;
        else if (cc.selfIntersection == Status::Unresolved && c.intersections == Status::Pass) c.intersections = Status::Unresolved;
    }
    for (const VertexBinding& b : s.binding) {
        if (b.tri < 0) ++c.unboundVertices;
        else c.maxBindingDistance = std::max(c.maxBindingDistance, b.distance);
    }
    if (s.binding.size() != m.pos.size()) c.unboundVertices += static_cast<int>(m.pos.size()) - static_cast<int>(s.binding.size());
    const bool geomFail = c.invertedFaces > 0 || c.degenerateFaces > 0 || c.facesOverDeviation > 0 || c.untriangulableFaces > 0 || c.intersections == Status::Fail || c.unboundVertices > 0 ||
                          c.maxBindingDistance > o.bindingTolerance;
    c.geometry = geomFail ? Status::Fail : c.intersections == Status::Unresolved ? Status::Unresolved : Status::Pass;
    c.detail = "connectivity " + std::string(toString(c.connectivity)) + (fails.empty() ? "" : " [" + fails + "]") + (opens.empty() ? "" : " {" + opens + "}");
    return c;
}

std::string TypedMeshCheck::text() const {
    char b[2048];
    std::snprintf(b, sizeof b,
                  "connectivity %s: faces %d (K %d, INC %d, DEC %d, special %d, unsupported %d); kind mismatches %d; faces without two wale edges %d\n"
                  "  course edges: interlocks %d, both-top %d, both-bottom %d, against wale %d; wale edges: yarn agrees %d, disagrees %d\n"
                  "  boundary: cast-on %d, bind-off %d, side %d; role on interior edge %d, border course without role %d, wrong side %d, off declared %d; "
                  "border wale edges %d\n"
                  "  strands %zu (closed %d, open %d), yarn ends %d, wale edges traced once %d / not once %d; strips %zu, broken %d, mixed shaping %d\n"
                  "geometry %s: inverted %d, degenerate %d, intersections %s (%d pairs), max binding distance %.3g mm, unbound %d, max face deviation %.4g mm (over %d), untriangulable %d\n",
                  toString(connectivity), faces, regular, increase, decrease, special, unsupported, kindMismatches, facesWithoutTwoWaleEdges, interlocks, bothTop,
                  bothBottom, courseAgainstWale, waleAgree, waleDisagree, castOnEdges, bindOffEdges, sideEdges, boundaryRoleOnInteriorEdge, borderCourseWithoutRole,
                  boundaryRoleWrongSide, boundaryRoleOffDeclared, borderWaleEdges, strands.strands.size(), closedStrands, openStrands, strands.yarnEnds,
                  strands.waleEdgesTracedOnce, strands.waleEdgesNotTracedOnce, strips.size(), stripsBroken, mixedShapingStrips, toString(geometry), invertedFaces,
                  degenerateFaces, toString(intersections), intersectingPairs, maxBindingDistance, unboundVertices, maxFaceDeviation, facesOverDeviation, untriangulableFaces);
    std::string t = b;
    if (poleStart > 0) t += "  pole-start faces " + std::to_string(poleStart) + "\n";
    for (const auto& [r, n] : unsupportedReasons) t += "  unsupported reason " + r + ": " + std::to_string(n) + "\n";
    return t;
}

std::vector<CellMeasure> measureCells(const TypedStitchMesh& s, double w, double h) {
    std::vector<CellMeasure> out;
    const PolygonMesh& m = s.mesh;
    const double slant = std::sqrt(0.25 * w * w + h * h);
    for (int f = 0; f < m.faces(); ++f) {
        CellMeasure c;
        c.face = f;
        c.kind = s.kind[f];
        c.area = m.faceArea(f) / (w * h);
        const int n = m.sides(f);
        const Vec3 normal = m.faceNormal(f);
        std::vector<double> course, wale;
        double minA = 360, maxA = 0, ccA = 0;
        for (int k = 0; k < n; ++k) {
            const double len = edgeLength(m, f, k);
            if (s.role[f][k] == StitchEdgeRole::Wale) wale.push_back(len); else course.push_back(len);
            // the corner at vertex k lies between edge k - 1 and edge k
            const StitchEdgeRole before = s.role[f][(k + n - 1) % n], after = s.role[f][k];
            const double ang = interiorAngle(m, f, k, normal);
            const bool beforeCourse = before != StitchEdgeRole::Wale, afterCourse = after != StitchEdgeRole::Wale;
            if (beforeCourse && afterCourse) ccA = ang;
            else if (beforeCourse != afterCourse) {
                minA = std::min(minA, ang);
                maxA = std::max(maxA, ang);
            }
        }
        c.minAngle = minA;
        c.maxAngle = maxA;
        c.courseCourseAngle = ccA;
        std::string why;
        auto check = [&](const char* name, double v, double lo, double hi) {
            if (v < lo || v > hi) {
                char b[96];
                std::snprintf(b, sizeof b, "%s%s %.3f", why.empty() ? "" : ", ", name, v);
                why += b;
            }
        };
        auto range = [](const std::vector<double>& v, double scale, double& lo, double& hi) {
            lo = 1e300;
            hi = 0;
            for (double x : v) {
                lo = std::min(lo, x / scale);
                hi = std::max(hi, x / scale);
            }
            if (v.empty()) lo = hi = 0;
        };
        switch (c.kind) {
        case StitchFaceKind::Regular: {
            double sc = 0, sw = 0;
            for (double x : course) sc += x;
            for (double x : wale) sw += x;
            c.widthRatio = course.empty() ? 0 : sc / course.size() / w;
            c.heightRatio = wale.empty() ? 0 : sw / wale.size() / h;
            range(course, w, c.minCourse, c.maxCourse);
            range(wale, h, c.minWale, c.maxWale);
            check("width", c.widthRatio, 0.75, 1.30);
            check("height", c.heightRatio, 0.75, 1.30);
            check("min angle", c.minAngle, 60, 120);
            check("max angle", c.maxAngle, 60, 120);
            break;
        }
        case StitchFaceKind::Increase:
        case StitchFaceKind::Decrease:
            range(course, w, c.minCourse, c.maxCourse);
            range(wale, slant, c.minWale, c.maxWale);
            c.widthRatio = c.minCourse;
            c.heightRatio = c.minWale;
            check("min course", c.minCourse, 0.75, 1.30);
            check("max course", c.maxCourse, 0.75, 1.30);
            check("min wale/slant", c.minWale, 0.75, 1.30);
            check("max wale/slant", c.maxWale, 0.75, 1.30);
            check("min wale-course angle", c.minAngle, 45, 135);
            check("max wale-course angle", c.maxAngle, 45, 135);
            check("course-course angle", c.courseCourseAngle, 150, 210);
            break;
        case StitchFaceKind::SpecialQuad:
            range(course, w, c.minCourse, c.maxCourse);
            range(wale, h, c.minWale, c.maxWale);
            check("min course", c.minCourse, 0.5, 1.5);
            check("max course", c.maxCourse, 0.5, 1.5);
            check("min wale", c.minWale, 0.5, 1.5);
            check("max wale", c.maxWale, 0.5, 1.5);
            check("min angle", std::min(c.minAngle, c.courseCourseAngle > 0 ? c.courseCourseAngle : 90.0), 30, 150);
            check("max angle", std::max(c.maxAngle, c.courseCourseAngle), 30, 150);
            break;
        case StitchFaceKind::PoleStart: {
            // proposed bounds (untested): the top course edge and both wale edges as for K, the
            // apex at the pole between 30 and 150 degrees (the ideal triangle with sides h and base w has 2 asin(w / 2h))
            range(course, w, c.minCourse, c.maxCourse);
            range(wale, h, c.minWale, c.maxWale);
            c.widthRatio = c.minCourse;
            c.heightRatio = 0.5 * (c.minWale + c.maxWale);
            double apex = 0;
            for (int k = 0; k < n; ++k)
                if (s.role[f][(k + n - 1) % n] == StitchEdgeRole::Wale && s.role[f][k] == StitchEdgeRole::Wale) apex = interiorAngle(m, f, k, normal);
            c.courseCourseAngle = apex;   // for a pole start: the apex (wale-wale) angle
            check("top course", c.minCourse, 0.75, 1.30);
            check("min wale", c.minWale, 0.75, 1.30);
            check("max wale", c.maxWale, 0.75, 1.30);
            check("apex angle", apex, 30, 150);
            break;
        }
        case StitchFaceKind::Unsupported:
            why = "unsupported";
            break;
        }
        c.withinBounds = why.empty();
        c.outside = why;
        out.push_back(c);
    }
    return out;
}

uint64_t typedDigest(const TypedStitchMesh& s) {
    uint64_t h = s.mesh.digest();
    for (StitchFaceKind k : s.kind) h = fnv1a64Bytes(&k, sizeof k, h);
    for (const std::vector<StitchEdgeRole>& r : s.role)
        for (StitchEdgeRole x : r) h = fnv1a64Bytes(&x, sizeof x, h);
    for (const std::vector<EdgeBoundary>& r : s.boundary)
        for (EdgeBoundary x : r) h = fnv1a64Bytes(&x, sizeof x, h);
    return h;
}

TypedStitchMesh typedFromStitchMesh(const StitchMesh& sm, const std::vector<DeclaredBoundary>& declared, double tolerance) {
    TypedStitchMesh t;
    t.mesh = sm.mesh;
    t.kind = sm.kind;
    t.role = sm.role;
    const int nF = t.mesh.faces();
    t.boundary.assign(nF, {});
    t.reason.assign(nF, "");
    const std::map<EdgeKey, std::vector<CornerUse>> uses = edgeUses(t.mesh);
    for (int f = 0; f < nF; ++f) {
        const int n = t.mesh.sides(f);
        t.boundary[f].assign(n, EdgeBoundary::Interior);
        if (t.kind[f] == StitchFaceKind::Unsupported) t.reason[f] = "unsupported-roles";
        for (int k = 0; k < n; ++k) {
            if (t.role[f][k] == StitchEdgeRole::Wale) continue;
            const int a = t.mesh.vertexAt(f, k), b = t.mesh.vertexAt(f, (k + 1) % n);
            if (uses.at(keyOf(a, b)).size() != 1) continue;
            for (const DeclaredBoundary& d : declared) {
                if (pointPolylineDistance(t.mesh.pos[a], d.loop) > tolerance || pointPolylineDistance(t.mesh.pos[b], d.loop) > tolerance) continue;
                const bool fits = (d.role == EdgeBoundary::CastOn && t.role[f][k] == StitchEdgeRole::Bottom) ||
                                  (d.role == EdgeBoundary::BindOff && t.role[f][k] == StitchEdgeRole::Top);
                if (fits) t.boundary[f][k] = d.role;
            }
        }
    }
    return t;
}

std::vector<ShapingCellMeasure> measureShapingCells(const TypedStitchMesh& s, double w, double h) {
    std::vector<ShapingCellMeasure> out;
    const PolygonMesh& m = s.mesh;
    for (int f = 0; f < m.faces(); ++f) {
        if (s.kind[f] != StitchFaceKind::Increase && s.kind[f] != StitchFaceKind::Decrease) continue;
        ShapingCellMeasure c;
        c.face = f;
        c.kind = s.kind[f];
        const int n = m.sides(f);
        const Vec3 normal = m.faceNormal(f);
        std::vector<std::pair<int, int>> bottoms, tops;
        c.minCourse = c.minWale = 1e300;
        for (int k = 0; k < n; ++k) {
            const int a = m.vertexAt(f, k), b = m.vertexAt(f, (k + 1) % n);
            const double len = norm(m.pos[b] - m.pos[a]);
            if (s.role[f][k] == StitchEdgeRole::Wale) {
                c.minWale = std::min(c.minWale, len / h);
                c.maxWale = std::max(c.maxWale, len / h);
            } else {
                c.minCourse = std::min(c.minCourse, len / w);
                c.maxCourse = std::max(c.maxCourse, len / w);
                (s.role[f][k] == StitchEdgeRole::Bottom ? bottoms : tops).push_back({a, b});
            }
        }
        // the sub-quads as the solver pairs them (relax_mesh.cpp buildModel): p q bottom, r s top, cycle p q r s
        c.minCorner = 360;
        c.minSubArea = 1e300;
        for (const auto& bt : bottoms)
            for (const auto& tp : tops) {
                const int q4[4] = {bt.first, bt.second, tp.first, tp.second};
                if (q4[0] == q4[2] || q4[0] == q4[3] || q4[1] == q4[2] || q4[1] == q4[3]) continue;
                ++c.subQuads;
                Vec3 P[4];
                for (int i = 0; i < 4; ++i) P[i] = m.pos[q4[i]];
                for (int i = 0; i < 4; ++i) {
                    const Vec3 a = P[(i + 3) % 4] - P[i], b = P[(i + 1) % 4] - P[i];
                    const double la = norm(a), lb = norm(b);
                    const double ang = la > 0 && lb > 0 ? std::acos(std::max(-1.0, std::min(1.0, dot(a, b) / (la * lb)))) * 180.0 / kPi : 0;
                    c.minCorner = std::min(c.minCorner, ang);
                    c.maxCorner = std::max(c.maxCorner, ang);
                }
                // both diagonal splits must turn like the face (a folded or bow-tie sub-quad has a triangle turning back)
                const Vec3 t0 = cross(P[1] - P[0], P[2] - P[0]), t1 = cross(P[2] - P[0], P[3] - P[0]);
                const Vec3 t2 = cross(P[2] - P[1], P[3] - P[1]), t3 = cross(P[3] - P[1], P[0] - P[1]);
                if (dot(t0, normal) <= 0 || dot(t1, normal) <= 0 || dot(t2, normal) <= 0 || dot(t3, normal) <= 0) ++c.foldedSubQuads;
                c.minSubArea = std::min(c.minSubArea, 0.25 * (norm(t0) + norm(t1) + norm(t2) + norm(t3)) / (w * h));
            }
        if (c.subQuads == 0) c.minCorner = c.minSubArea = 0;
        std::string why;
        auto check = [&](const char* name, double v, double lo, double hi) {
            if (v < lo || v > hi) {
                char b[96];
                std::snprintf(b, sizeof b, "%s%s %.3f", why.empty() ? "" : ", ", name, v);
                why += b;
            }
        };
        check("min course/w", c.minCourse, 0.75, 1.30);
        check("max course/w", c.maxCourse, 0.75, 1.30);
        check("min wale/h", c.minWale, 0.75, 1.30);
        check("max wale/h", c.maxWale, 0.75, 1.30);
        c.interfacesWithin = why.empty();
        const size_t before = why.size();
        check("sub-quads", c.subQuads, 2, 2);
        check("min sub-quad corner", c.minCorner, 45, 135);
        check("max sub-quad corner", c.maxCorner, 45, 135);
        check("min sub-quad area/wh", c.minSubArea, 0.1, 1e9);
        check("folded sub-quads", c.foldedSubQuads, 0, 0);
        c.subQuadsWithin = why.size() == before;
        c.within = c.interfacesWithin && c.subQuadsWithin;
        c.outside = why;
        out.push_back(c);
    }
    return out;
}

ShapingSummary summarizeShaping(const std::vector<ShapingCellMeasure>& cells) {
    ShapingSummary s;
    for (const ShapingCellMeasure& c : cells) {
        ++s.cells;
        s.interfacesWithin += c.interfacesWithin;
        s.subQuadWithin += c.subQuadsWithin;
        s.within += c.within;
        s.invertedSubQuads += c.foldedSubQuads;
    }
    return s;
}


} // namespace knit
