#include "stitchmesh/virtual_shift.h"

#include <algorithm>
#include <cstdio>
#include <map>
#include <set>
#include <utility>

namespace knit {

const char* toString(ShiftReading r) {
    switch (r) {
    case ShiftReading::PathSideDown: return "path-side-down";
    case ShiftReading::PathSideUp: return "path-side-up";
    case ShiftReading::None: return "none";
    }
    return "?";
}

namespace {

using EdgeKey = std::pair<int, int>;
EdgeKey keyOf(int a, int b) { return {std::min(a, b), std::max(a, b)}; }

// The path edges of each strip.
std::vector<std::vector<const ShiftPathEdge*>> pathByStrip(const TypedStitchMesh& s, int strips) {
    std::vector<std::vector<const ShiftPathEdge*>> by(strips);
    for (const ShiftPathEdge& e : s.shiftPath)
        if (e.strip >= 0 && e.strip < strips) by[e.strip].push_back(&e);
    return by;
}

int stripCount(const TypedStitchMesh& s) {
    int n = 0;
    for (int k : s.strip) n = std::max(n, k + 1);
    return n;
}

// Course edges with two faces: (offering face, receiving face) per undirected edge; course edges with one face counted
// by role. Degenerate edges (one vertex twice) are skipped. Re-derived from a corner list and the roles.
struct CourseIndex {
    std::set<std::pair<int, int>> interlocks;   // (face whose Top it is, face whose Bottom it is)
    int boundaryBottom = 0, boundaryTop = 0, nonManifold = 0, mismatched = 0;
};

CourseIndex courseIndex(const std::vector<std::vector<int>>& corners, const std::vector<std::vector<StitchEdgeRole>>& role) {
    struct Use {
        int face, corner;
    };
    std::map<EdgeKey, std::vector<Use>> uses;
    for (int f = 0; f < static_cast<int>(corners.size()); ++f) {
        const int n = static_cast<int>(corners[f].size());
        for (int k = 0; k < n; ++k) {
            const int a = corners[f][k], b = corners[f][(k + 1) % n];
            if (a == b || role[f][k] == StitchEdgeRole::Wale) continue;
            uses[keyOf(a, b)].push_back(Use{f, k});
        }
    }
    CourseIndex ci;
    for (const auto& [key, list] : uses) {
        if (list.size() == 1) {
            (role[list[0].face][list[0].corner] == StitchEdgeRole::Top ? ci.boundaryTop : ci.boundaryBottom) += 1;
        } else if (list.size() == 2) {
            const bool topA = role[list[0].face][list[0].corner] == StitchEdgeRole::Top;
            const bool topB = role[list[1].face][list[1].corner] == StitchEdgeRole::Top;
            if (topA == topB) ++ci.mismatched;
            else ci.interlocks.insert(topA ? std::make_pair(list[0].face, list[1].face) : std::make_pair(list[1].face, list[0].face));
        } else {
            ++ci.nonManifold;
        }
    }
    return ci;
}

} // namespace

bool buildVirtualShift(const TypedStitchMesh& s, ShiftReading reading, VirtualShift& out, std::string* err) {
    const int nF = s.mesh.faces();
    if (static_cast<int>(s.strip.size()) != nF) {
        if (err) *err = "virtual shift: the record has no strips";
        return false;
    }
    out = VirtualShift{};
    out.reading = reading;
    out.corners.resize(nF);
    out.reRead.assign(nF, 0);
    for (int f = 0; f < nF; ++f)
        for (int k = 0; k < s.mesh.sides(f); ++k) out.corners[f].push_back(s.mesh.vertexAt(f, k));
    switch (reading) {
    case ShiftReading::PathSideDown:
        out.rule = "slot 0 of strip k (its W- edge on the path) reads its path-side corners one ring lower: ring k+1's path vertex as ring k's, ring k's as "
                   "ring k-1's; strip 0 has no ring below, so its W- edge collapses to the cast-on path vertex";
        break;
    case ShiftReading::PathSideUp:
        out.rule = "negative control: slot 0's path-side corners read one ring higher";
        break;
    case ShiftReading::None:
        out.rule = "negative control: the stored closed rings";
        return true;
    }
    const int strips = stripCount(s);
    const std::vector<std::vector<const ShiftPathEdge*>> path = pathByStrip(s, strips);
    for (int f = 0; f < nF; ++f) {
        const int k = s.strip[f];
        if (k < 0 || path[k].size() != 1) continue;
        const ShiftPathEdge& e = *path[k][0];
        // the face whose W- edge is the path edge: its entry wale corner runs from the top path vertex down to the bottom one
        const int entry = entryWaleCorner(s.role[f]);
        if (entry < 0) continue;
        const int n = s.mesh.sides(f);
        const int a = s.mesh.vertexAt(f, entry), b = s.mesh.vertexAt(f, (entry + 1) % n);
        if (!(a == e.top && b == e.bottom)) continue;
        int newBottom = e.bottom, newTop = e.top;
        if (reading == ShiftReading::PathSideDown) {
            newTop = e.bottom;
            newBottom = (k > 0 && path[k - 1].size() == 1) ? path[k - 1][0]->bottom : e.bottom;
        } else {
            newBottom = e.top;
            newTop = (k + 1 < strips && path[k + 1].size() == 1) ? path[k + 1][0]->top : e.top;
        }
        for (int c = 0; c < n; ++c) {
            const int v = s.mesh.vertexAt(f, c);
            if (v == e.bottom) out.corners[f][c] = newBottom;
            else if (v == e.top) out.corners[f][c] = newTop;
        }
        out.reRead[f] = 1;
    }
    return true;
}

VirtualShiftCheck validateVirtualShift(const TypedStitchMesh& s, const VirtualShift& v) {
    VirtualShiftCheck c;
    const int nF = s.mesh.faces();
    c.strips = stripCount(s);
    // the path: one wale edge per strip, a connected chain, each a wale edge of the stored mesh
    const std::vector<std::vector<const ShiftPathEdge*>> path = pathByStrip(s, c.strips);
    c.pathEdges = static_cast<int>(s.shiftPath.size());
    std::set<EdgeKey> waleEdges;
    for (int f = 0; f < nF; ++f)
        for (int k = 0; k < s.mesh.sides(f); ++k)
            if (s.role[f][k] == StitchEdgeRole::Wale) waleEdges.insert(keyOf(s.mesh.vertexAt(f, k), s.mesh.vertexAt(f, (k + 1) % s.mesh.sides(f))));
    // a gathered pole start: every face of strip 0 is a pole start (and strip 0 has no path edge)
    int strip0 = 0, strip0Pole = 0;
    for (int f = 0; f < nF; ++f)
        if (s.strip[f] == 0) {
            ++strip0;
            if (s.kind[f] == StitchFaceKind::PoleStart) ++strip0Pole;
        }
    c.poleStart = strip0 > 0 && strip0Pole == strip0;
    const int firstHelixStrip = c.poleStart ? 1 : 0;
    for (int k = 0; k < c.strips; ++k) {
        if (c.poleStart && k == 0) {
            if (!path[0].empty()) ++c.stripsWithSeveralPathEdges;   // a pole round has no path edge
            continue;
        }
        if (path[k].empty()) ++c.stripsWithoutPathEdge;
        if (path[k].size() > 1) ++c.stripsWithSeveralPathEdges;
        for (const ShiftPathEdge* e : path[k])
            if (!waleEdges.count(keyOf(e->bottom, e->top))) ++c.pathEdgesNotWale;
        if (k > 0 && path[k].size() == 1 && path[k - 1].size() == 1 && path[k - 1][0]->top != path[k][0]->bottom) ++c.pathBreaks;
    }
    // the virtual polygons as a polygon mesh (repeated corners kept: a collapsed edge has one vertex twice)
    PolygonMesh vm;
    vm.pos = s.mesh.pos;
    bool sameSizes = static_cast<int>(v.corners.size()) == nF;
    for (int f = 0; f < nF && sameSizes; ++f) {
        if (static_cast<int>(v.corners[f].size()) != s.mesh.sides(f)) sameSizes = false;
        else vm.addFace(v.corners[f]);
    }
    c.countsUnchanged = sameSizes;
    if (!sameSizes) {
        c.status = Status::Fail;
        c.detail = "the virtual polygons do not match the stored faces";
        return c;
    }
    // start steps: a collapsed wale edge; faces without two proper wale edges otherwise
    for (int f = 0; f < nF; ++f) {
        int proper = 0, collapsed = 0;
        const int n = vm.sides(f);
        for (int k = 0; k < n; ++k) {
            if (s.role[f][k] != StitchEdgeRole::Wale) continue;
            if (vm.vertexAt(f, k) == vm.vertexAt(f, (k + 1) % n)) ++collapsed; else ++proper;
        }
        if (collapsed == 1 && proper == 1) {
            ++c.startSteps;
            c.startFace = f;
        } else if (proper != 2) {
            ++c.virtualFacesWithoutTwoWaleEdges;
        }
    }
    // the yarn traced on the virtual polygons (a collapsed edge has no face across: the strand starts there)
    const StrandTrace t = traceYarnStrands(vm, s.role);
    c.strands = static_cast<int>(t.strands.size());
    if (c.poleStart)
        for (size_t i = 0; i < t.strands.size(); ++i) {
            bool allStrip0 = !t.strands[i].empty();
            for (int f : t.strands[i]) allStrip0 = allStrip0 && s.strip[f] == 0;
            if (allStrip0 && t.closed[i]) c.poleRoundFaces = static_cast<int>(t.strands[i].size());
        }
    size_t longest = 0;
    for (size_t i = 0; i < t.strands.size(); ++i) {
        (t.closed[i] ? c.closedStrands : c.openStrands) += 1;
        if (t.strands[i].size() > t.strands[longest].size()) longest = i;
    }
    if (!t.strands.empty()) {
        std::vector<int> order = t.strands[longest];
        // orient along the yarn: the first face's W+ edge leads to the second face
        if (order.size() >= 2) {
            const int f = order[0], g = order[1];
            const int ex = exitWaleCorner(s.role[f]);
            bool forward = false;
            if (ex >= 0) {
                const EdgeKey e = keyOf(vm.vertexAt(f, ex), vm.vertexAt(f, (ex + 1) % vm.sides(f)));
                for (int k = 0; k < vm.sides(g); ++k)
                    if (keyOf(vm.vertexAt(g, k), vm.vertexAt(g, (k + 1) % vm.sides(g))) == e) forward = true;
            }
            if (!forward) std::reverse(order.begin(), order.end());
        }
        c.strandOrder = order;
        c.longestStrand = static_cast<int>(order.size());
        if (c.poleStart)
            c.oneStrandThroughAll = c.strands == 2 && c.openStrands == 1 && c.closedStrands == 1 && c.poleRoundFaces == strip0 && c.longestStrand == nF - strip0;
        else
            c.oneStrandThroughAll = c.strands == 1 && c.openStrands == 1 && c.longestStrand == nF;
        if (!order.empty()) {
            c.startsInFirstStrip = s.strip[order.front()] == firstHelixStrip && order.front() == c.startFace;
            c.endsInLastStrip = s.strip[order.back()] == c.strips - 1;
        }
    }
    // end steps: a face whose W+ edge has no face across in the virtual mesh
    {
        std::map<EdgeKey, int> uses;
        for (int f = 0; f < nF; ++f)
            for (int k = 0; k < vm.sides(f); ++k) {
                const int a = vm.vertexAt(f, k), b = vm.vertexAt(f, (k + 1) % vm.sides(f));
                if (a != b) uses[keyOf(a, b)] += 1;
            }
        for (const auto& [key, n] : uses)
            if (n > 2) ++c.nonManifoldVirtualEdges;
        for (int f = 0; f < nF; ++f) {
            const int ex = exitWaleCorner(s.role[f]);
            if (ex < 0) continue;
            const int a = vm.vertexAt(f, ex), b = vm.vertexAt(f, (ex + 1) % vm.sides(f));
            if (a != b && uses[keyOf(a, b)] == 1) {
                ++c.endSteps;
                c.endFace = f;
            }
        }
    }
    if (c.endsInLastStrip && !c.strandOrder.empty()) c.endsInLastStrip = c.strandOrder.back() == c.endFace;
    // interlocks and boundary course edges, stored against virtual
    std::vector<std::vector<int>> stored(nF);
    for (int f = 0; f < nF; ++f)
        for (int k = 0; k < s.mesh.sides(f); ++k) stored[f].push_back(s.mesh.vertexAt(f, k));
    const CourseIndex cs = courseIndex(stored, s.role), cv = courseIndex(v.corners, s.role);
    c.interlocksStored = static_cast<int>(cs.interlocks.size());
    c.interlocksVirtual = static_cast<int>(cv.interlocks.size());
    for (const auto& p : cs.interlocks)
        if (!cv.interlocks.count(p)) ++c.interlocksMissing;
    for (const auto& p : cv.interlocks)
        if (!cs.interlocks.count(p)) ++c.interlocksExtra;
    c.boundaryBottomStored = cs.boundaryBottom;
    c.boundaryTopStored = cs.boundaryTop;
    c.boundaryBottomVirtual = cv.boundaryBottom;
    c.boundaryTopVirtual = cv.boundaryTop;
    c.nonManifoldVirtualEdges += cv.nonManifold;
    // verdict
    std::string fails;
    auto add = [&](bool bad, const char* what) {
        if (!bad) return;
        if (!fails.empty()) fails += "; ";
        fails += what;
    };
    add(c.stripsWithoutPathEdge > 0, "a strip without a path edge");
    add(c.stripsWithSeveralPathEdges > 0, "a strip with several path edges");
    add(c.pathBreaks > 0, "the path is not a connected chain");
    add(c.pathEdgesNotWale > 0, "a path edge that is no wale edge");
    add(!c.oneStrandThroughAll, c.poleStart ? "not one closed pole round plus one open strand through every other face" : "not one open strand through every face");
    add(c.startSteps != 1 || !c.startsInFirstStrip, "the strand does not start at exactly one start step in the first strip");
    add(c.endSteps != 1 || !c.endsInLastStrip, "the strand does not end at exactly one end step in the last strip");
    add(c.virtualFacesWithoutTwoWaleEdges > 0, "virtual faces without two wale edges");
    add(c.nonManifoldVirtualEdges > 0, "non-manifold virtual edges");
    add(c.interlocksMissing > 0 || c.interlocksExtra > 0, "course interlocks changed");
    add(c.boundaryBottomVirtual != c.boundaryBottomStored || c.boundaryTopVirtual != c.boundaryTopStored, "boundary course edges changed");
    c.status = fails.empty() ? Status::Pass : Status::Fail;
    c.detail = fails;
    return c;
}

TypedStitchMesh virtualTypedMesh(const TypedStitchMesh& s, const VirtualShift& v) {
    TypedStitchMesh t;
    t.mesh.pos = s.mesh.pos;
    t.binding = s.binding;
    t.strip = s.strip;
    t.slot = s.slot;
    for (int f = 0; f < s.mesh.faces(); ++f) {
        const std::vector<int>& c = v.corners[f];
        const int n = static_cast<int>(c.size());
        std::vector<int> corners;
        std::vector<StitchEdgeRole> roles;
        std::vector<EdgeBoundary> bounds;
        for (int k = 0; k < n; ++k) {
            if (c[k] == c[(k + 1) % n]) continue;   // a collapsed edge
            corners.push_back(c[k]);
            roles.push_back(s.role[f][k]);
            bounds.push_back(s.boundary[f][k]);
        }
        t.mesh.addFace(corners);
        t.role.push_back(roles);
        t.boundary.push_back(bounds);
        const bool changed = static_cast<int>(corners.size()) != n;
        t.kind.push_back(changed ? StitchFaceKind::Unsupported : s.kind[f]);
        t.reason.push_back(changed ? "helix-start-step" : (f < static_cast<int>(s.reason.size()) ? s.reason[f] : ""));
    }
    return t;
}

std::string VirtualShiftCheck::text() const {
    char b[1600];
    std::snprintf(b, sizeof b,
                  "virtual shift %s: path edges %d for %d strips (missing %d, several %d, breaks %d, not wale %d); strands %d (open %d, closed %d), longest %d, "
                  "one through all %s; start steps %d (face %d, first strip %s), end steps %d (face %d, last strip %s); faces without two wale edges %d; "
                  "non-manifold %d; interlocks stored %d virtual %d (missing %d, extra %d); boundary bottom %d/%d, top %d/%d%s%s%s\n",
                  toString(status), pathEdges, strips, stripsWithoutPathEdge, stripsWithSeveralPathEdges, pathBreaks, pathEdgesNotWale, strands, openStrands,
                  closedStrands, longestStrand, oneStrandThroughAll ? "yes" : "no", startSteps, startFace, startsInFirstStrip ? "yes" : "no", endSteps, endFace,
                  endsInLastStrip ? "yes" : "no", virtualFacesWithoutTwoWaleEdges, nonManifoldVirtualEdges, interlocksStored, interlocksVirtual, interlocksMissing,
                  interlocksExtra, boundaryBottomStored, boundaryBottomVirtual, boundaryTopStored, boundaryTopVirtual,
                  poleStart ? ("; pole start: closed round of " + std::to_string(poleRoundFaces) + " faces, helix from strip 1").c_str() : "", detail.empty() ? "" : ": ",
                  detail.c_str());
    return b;
}

} // namespace knit
