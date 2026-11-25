// Connected pieces of an imported file's components (surface_pieces.h).
#include "model/surface_pieces.h"

#include <algorithm>
#include <numeric>

#include "model/hash.h"
#include "model/mesh_check.h"

namespace knit {

namespace {

int findRoot(std::vector<int>& parent, int i) {
    while (parent[i] != i) {
        parent[i] = parent[parent[i]];  // path halving
        i = parent[i];
    }
    return i;
}

} // namespace

SurfacePieces findSurfacePieces(const std::vector<Vec3>& pos, const std::vector<std::array<int, 3>>& tri, const std::vector<int>& nodeOfTriangle,
                                double weldTolerance) {
    SurfacePieces out;
    const int nt = static_cast<int>(tri.size());
    out.pieceOfTriangle.assign(nt, -1);
    out.localTriangle.assign(nt, -1);
    TriMesh all;
    all.pos = pos;
    all.tri = tri;
    const WeldedMesh welded = weldVertices(all, weldTolerance);
    const std::vector<int>& weld = welded.remap;
    const std::vector<Vec3>& weldedPos = welded.mesh.pos;
    // union of triangles sharing a welded vertex within the same component
    std::vector<int> parent(nt);
    std::iota(parent.begin(), parent.end(), 0);
    // per welded vertex: the first triangle of each component that uses it (usually one entry)
    std::vector<std::vector<std::pair<int, int>>> firstOfVertex(weldedPos.size());  // (component, triangle)
    std::vector<char> degenerate(nt, 0);
    for (int t = 0; t < nt; ++t) {
        const int a = weld[tri[t][0]], b = weld[tri[t][1]], c = weld[tri[t][2]];
        if (a == b || b == c || a == c) {
            degenerate[t] = 1;
            continue;
        }
        const int node = nodeOfTriangle.empty() ? 0 : nodeOfTriangle[t];
        for (int v : {a, b, c}) {
            int other = -1;
            for (const std::pair<int, int>& e : firstOfVertex[v])
                if (e.first == node) other = e.second;
            if (other < 0) firstOfVertex[v].push_back({node, t});
            else parent[findRoot(parent, t)] = findRoot(parent, other);
        }
    }
    // pieces in order of their lowest triangle
    std::vector<int> pieceOfRoot(nt, -1);
    std::vector<int> ordinalOfNode;
    for (int t = 0; t < nt; ++t) {
        if (degenerate[t]) continue;
        const int r = findRoot(parent, t);
        if (pieceOfRoot[r] < 0) {
            pieceOfRoot[r] = static_cast<int>(out.pieces.size());
            SurfacePiece p;
            p.node = nodeOfTriangle.empty() ? 0 : nodeOfTriangle[t];
            if (p.node >= static_cast<int>(ordinalOfNode.size())) ordinalOfNode.resize(p.node + 1, 0);
            p.ordinal = ++ordinalOfNode[p.node];
            out.pieces.push_back(p);
        }
        const int pi = pieceOfRoot[r];
        out.pieceOfTriangle[t] = pi;
        out.localTriangle[t] = static_cast<int>(out.pieces[pi].fileTriangle.size());
        out.pieces[pi].fileTriangle.push_back(t);
    }
    // vertices numbered per piece (a welded vertex may touch pieces of other components)
    std::vector<int> localOf(weldedPos.size(), -1);
    for (SurfacePiece& p : out.pieces) {
        for (int t : p.fileTriangle) {
            std::array<int, 3> local{};
            for (int k = 0; k < 3; ++k) {
                const int w = weld[tri[t][k]];
                if (localOf[w] < 0) {
                    localOf[w] = static_cast<int>(p.mesh.pos.size());
                    p.mesh.pos.push_back(weldedPos[w]);
                }
                local[k] = localOf[w];
            }
            p.mesh.tri.push_back(local);
        }
        for (int t : p.fileTriangle)
            for (int k = 0; k < 3; ++k) localOf[weld[tri[t][k]]] = -1;
    }
    for (SurfacePiece& p : out.pieces) p.mesh.computeNormals();
    return out;
}

double pieceWeldTolerance(const std::vector<Vec3>& pos) {
    if (pos.empty()) return 0;
    Vec3 lo = pos[0], hi = pos[0];
    for (const Vec3& p : pos) {
        lo = componentMin(lo, p);
        hi = componentMax(hi, p);
    }
    return 1e-6 * norm(hi - lo);
}

uint64_t surfaceDigest(const std::vector<Vec3>& pos, const std::vector<std::array<int, 3>>& tri) {
    TriMesh m;
    m.pos = pos;
    m.tri = tri;
    return meshDigest(m);
}

bool pieceLocalPoint(const SurfacePieces& p, const SurfacePoint& filePoint, int& piece, SurfacePoint& local) {
    if (filePoint.tri < 0 || filePoint.tri >= static_cast<int>(p.pieceOfTriangle.size())) return false;
    piece = p.pieceOfTriangle[filePoint.tri];
    if (piece < 0) return false;
    local = filePoint;
    local.tri = p.localTriangle[filePoint.tri];
    return true;
}

int pieceOfMarking(const SurfacePieces& p, const SurfaceMarking& m) {
    if (m.points.empty()) return -1;
    const int t = m.points.front().tri;
    if (t < 0 || t >= static_cast<int>(p.pieceOfTriangle.size())) return -1;
    return p.pieceOfTriangle[t];
}

bool markingInKnitComponent(const SurfacePieces& p, const SurfaceMarking& m, int piece, int knitComponent) {
    return m.knitComponent == knitComponent && pieceOfMarking(p, m) == piece;
}

uint64_t effectiveMarkingsDigest(const MarkingSet& set, const SurfacePieces& pieces, int piece, int knitComponent) {
    uint64_t h = fnv1a64Bytes(&set.digest, sizeof set.digest);
    h = fnv1a64(set.preparation, h);
    for (const SurfaceMarking& m : set.markings) {
        if (m.proposal == ProposalState::Offered || m.proposal == ProposalState::Rejected) continue;
        if (!markingInKnitComponent(pieces, m, piece, knitComponent)) continue;
        const int head[3] = {m.id, static_cast<int>(m.meaning), static_cast<int>(m.shape)};
        h = fnv1a64Bytes(head, sizeof head, h);
        for (const std::pair<std::string, std::string>& kv : m.settings) h = fnv1a64(kv.first + "=" + kv.second + ";", h);
        for (const SurfacePoint& p : m.points) {
            h = fnv1a64Bytes(&p.tri, sizeof p.tri, h);
            h = fnv1a64Bytes(p.bary, sizeof p.bary, h);
        }
    }
    return h;
}

uint64_t proposalBasis(uint64_t fileDigest, int piece, const std::string& preparation, double mmPerUnit, double w, double h) {
    uint64_t hash = fnv1a64Bytes(&fileDigest, sizeof fileDigest);
    hash = fnv1a64(preparation, hash);
    const double numbers[3] = {mmPerUnit, w, h};
    hash = fnv1a64Bytes(numbers, sizeof numbers, hash);
    return fnv1a64Bytes(&piece, sizeof piece, hash);
}

namespace {

// A marking the mesher takes as input: the user's own, or an accepted proposal.
bool inUse(const SurfaceMarking& m) { return m.origin == MarkingOrigin::User || m.proposal == ProposalState::Accepted; }

bool roleTaken(const std::vector<MarkingRole>& roles, MarkingRole r) {
    for (MarkingRole x : roles)
        if (x == r) return true;
    return false;
}

bool listed(const std::vector<int>& ids, int id) {
    for (int x : ids)
        if (x == id) return true;
    return false;
}

} // namespace

SetupConflicts setupConflicts(const MarkingSet& s, const SurfacePieces& pieces, int piece, int setup) {
    SetupConflicts c;
    std::vector<MarkingRole> roles;
    int knitComponent = 0;   // the setup's knitting component: its markings all carry it
    for (const SurfaceMarking& m : s.markings)
        if (setup > 0 && m.setup == setup && markingRole(m.meaning) != MarkingRole::None) {
            roles.push_back(markingRole(m.meaning));
            knitComponent = m.knitComponent;
        }
    for (const SurfaceMarking& m : s.markings) {
        if (m.setup == setup || !inUse(m) || !markingInKnitComponent(pieces, m, piece, knitComponent) || !roleTaken(roles, markingRole(m.meaning))) continue;
        if (m.origin == MarkingOrigin::Proposal) c.acceptedProposals.push_back(m.id);
        else c.userMarkings.push_back(m.id);
    }
    return c;
}

bool acceptSetup(MarkingSet& s, const SurfacePieces& pieces, int piece, int setup, uint64_t currentBasis, bool replaceUser) {
    const std::vector<int> ids = setupMarkings(s, setup);
    if (ids.empty()) return false;
    for (int id : ids) {
        const SurfaceMarking* m = findMarking(s, id);
        if (m->origin != MarkingOrigin::Proposal || m->proposal != ProposalState::Offered || m->basis != currentBasis || pieceOfMarking(pieces, *m) != piece)
            return false;
    }
    const SetupConflicts c = setupConflicts(s, pieces, piece, setup);
    if (!c.userMarkings.empty() && !replaceUser) return false;
    for (int id : c.acceptedProposals) findMarking(s, id)->proposal = ProposalState::Rejected;
    std::vector<SurfaceMarking> kept;
    for (const SurfaceMarking& m : s.markings)
        if (!listed(c.userMarkings, m.id)) kept.push_back(m);
    s.markings = kept;
    for (int id : ids) findMarking(s, id)->proposal = ProposalState::Accepted;
    return true;
}

MarkingSet previewMarkings(const MarkingSet& s, const SurfacePieces& pieces, int piece, int setup, std::vector<int>& previewIds) {
    const SetupConflicts c = setupConflicts(s, pieces, piece, setup);
    MarkingSet out = s;
    out.markings.clear();
    previewIds.clear();
    for (const SurfaceMarking& m : s.markings) {
        if (listed(c.acceptedProposals, m.id) || listed(c.userMarkings, m.id)) continue;   // shadowed by the setup
        if (setup > 0 && m.setup == setup && m.origin == MarkingOrigin::Proposal && m.proposal == ProposalState::Offered) previewIds.push_back(m.id);
        out.markings.push_back(m);
    }
    return out;
}

void mergeProposals(MarkingSet& s, const SurfacePieces& pieces, int piece, const std::vector<SurfaceMarking>& offered, std::vector<int>* newSetups) {
    const int knitComponent = offered.empty() ? 0 : offered.front().knitComponent;   // the caller stamps every offered proposal alike
    std::vector<SurfaceMarking> kept;
    for (const SurfaceMarking& m : s.markings) {
        const bool replaced = m.origin == MarkingOrigin::Proposal && (m.proposal == ProposalState::Offered || m.proposal == ProposalState::Rejected) &&
                              markingInKnitComponent(pieces, m, piece, knitComponent);
        if (!replaced) kept.push_back(m);
    }
    s.markings = kept;
    appendProposals(s, offered, newSetups);
}

void appendProposals(MarkingSet& s, const std::vector<SurfaceMarking>& offered, std::vector<int>* newSetups) {
    int highestSetup = 0;
    for (const SurfaceMarking& m : s.markings) highestSetup = std::max(highestSetup, m.setup);
    if (newSetups) newSetups->clear();
    std::vector<std::pair<int, int>> setupNumber;   // the proposal's own setup -> its number in this set
    for (const SurfaceMarking& o : offered) {
        SurfaceMarking m = o;
        m.id = s.nextId++;
        m.proposalId = m.id;
        m.origin = MarkingOrigin::Proposal;
        m.proposal = ProposalState::Offered;
        if (o.setup > 0) {
            int number = 0;
            for (const std::pair<int, int>& sn : setupNumber)
                if (sn.first == o.setup) number = sn.second;
            if (number == 0) {
                number = ++highestSetup;
                setupNumber.push_back({o.setup, number});
                if (newSetups) newSetups->push_back(number);
            }
            m.setup = number;
        }
        s.markings.push_back(m);
    }
}

} // namespace knit
