// Stitch records of a stitch-mesh yarn: one per loop, with parents, ownership and frame (guided_records.h).
#include "stitchmesh/guided_records.h"

#include <algorithm>
#include <map>
#include <utility>

#include "model/mesh_query.h"

namespace knit {

namespace {

// A course edge by its two vertex ids, smaller first: both faces of the edge name it the same way.
std::pair<int, int> edgeOf(const YarnPort& p) { return {std::min(p.a, p.b), std::max(p.a, p.b)}; }

StitchOp opOfFace(StitchFaceKind kind, bool& known) {
    known = true;
    switch (kind) {
    case StitchFaceKind::Regular: return StitchOp::K;
    case StitchFaceKind::Increase: return StitchOp::Inc;
    case StitchFaceKind::Decrease: return StitchOp::Dec;
    case StitchFaceKind::PoleStart: return StitchOp::K;  // a first-turn loop of the gathered pole round
    case StitchFaceKind::SpecialQuad:
    case StitchFaceKind::Unsupported: break;
    }
    known = false;
    return StitchOp::K;
}

Vec3 centroidOf(const TypedStitchMesh& s, const std::vector<int>& corners) {
    Vec3 c{0, 0, 0};
    for (int v : corners) c = c + s.mesh.pos[v];
    return corners.empty() ? c : c * (1.0 / static_cast<double>(corners.size()));
}

} // namespace

bool guidedStitchRecords(const TypedStitchMesh& s, const std::vector<std::vector<int>>& faces, const YarnEmbedding& e, const YarnConversion& cv,
                         const TriMesh& reference, int expectedInterlocks, CurveSet& curves, GuidedRecordCounts& counts, std::string* err) {
    auto fail = [&](const std::string& why) {
        if (err) *err = "stitch records: " + why;
        return false;
    };
    const int nFaces = static_cast<int>(faces.size());
    if (static_cast<int>(s.kind.size()) != nFaces || cv.u0.size() != e.pieces.size()) return fail("the mesh, the embedding and the conversion disagree in size");
    // Each face's row pieces and heads, in parameter order (row pieces) and head order (heads).
    std::vector<std::vector<int>> rows(nFaces), heads(nFaces);
    for (size_t i = 0; i < e.pieces.size(); ++i) {
        const YarnPiece& p = e.pieces[i];
        if (p.face < 0 || p.face >= nFaces) continue;  // the binding thread
        if (cv.u0[i] >= cv.u1[i]) return fail("piece " + std::to_string(i) + " runs against the strand's parameter");
        if (p.kind == YarnPieceKind::Row) rows[p.face].push_back(static_cast<int>(i));
        else if (p.kind == YarnPieceKind::Head) heads[p.face].push_back(static_cast<int>(i));
    }
    TriMesh ref = reference;
    if (ref.vnormal.size() != ref.pos.size()) ref.computeNormals();
    const MeshQuery query(ref);

    KnitTopology topo;
    std::map<std::pair<int, int>, int> loopOfTopEdge;  // course edge -> the loop passing it as its top
    std::vector<std::vector<int>> loopsOfFace(nFaces);
    std::map<int, int> nextIndexInTurn;
    counts = GuidedRecordCounts{};
    for (int f = 0; f < nFaces; ++f) {
        std::vector<int>& r = rows[f];
        if (r.size() < 2) return fail("face " + std::to_string(f) + " has " + std::to_string(r.size()) + " row pieces (a loop needs two)");
        std::sort(r.begin(), r.end(), [&cv](int a, int b) { return cv.u0[a] < cv.u0[b]; });
        const int yarn = cv.strandOfPiece[r.front()];
        for (int i : r)
            if (cv.strandOfPiece[i] != yarn) return fail("face " + std::to_string(f) + "'s row pieces lie on two strands");
        bool known = false;
        const StitchOp op = opOfFace(s.kind[f], known);
        if (!known) return fail("face " + std::to_string(f) + " is " + toString(s.kind[f]) + ": no stitch of production corresponds");
        const Vec3 origin = centroidOf(s, faces[f]);
        const MeshQuery::Hit hit = query.closest(origin);
        Vec3 N{0, 0, 1};
        if (hit.tri >= 0) {
            const std::array<int, 3>& t = ref.tri[hit.tri];
            N = normalized(ref.vnormal[t[0]] + ref.vnormal[t[1]] + ref.vnormal[t[2]]);
        }
        const BSpline& curve = curves.yarns[yarn].curve;
        const int turn = s.strip.empty() ? 0 : s.strip[f];
        const int nLoops = static_cast<int>(r.size()) - 1;
        for (int j = 0; j < nLoops; ++j) {
            StitchRecord rec;
            rec.id = static_cast<int>(topo.stitches.size());
            rec.yarn = yarn;
            rec.uBegin = j == 0 ? cv.u0[r[0]] : 0.5 * (cv.u0[r[j]] + cv.u1[r[j]]);
            rec.uEnd = j == nLoops - 1 ? cv.u1[r[j + 1]] : 0.5 * (cv.u0[r[j + 1]] + cv.u1[r[j + 1]]);
            rec.turn = turn;
            rec.index = nextIndexInTurn[turn]++;
            rec.op = op;
            rec.incSlot = op == StitchOp::Inc ? j : 0;
            rec.origin = origin;
            const Vec3 chord = curve.eval(rec.uEnd) - curve.eval(rec.uBegin);
            const Vec3 U = chord - N * dot(chord, N);
            rec.N = N;
            rec.U = norm(U) > 1e-12 ? normalized(U) : Vec3{1, 0, 0};
            rec.V = cross(rec.N, rec.U);
            loopOfTopEdge[edgeOf(e.pieces[r[j]].end)] = rec.id;  // row j ends at the top edge this loop passes
            loopsOfFace[f].push_back(rec.id);
            topo.stitches.push_back(rec);
            ++counts.loops;
            counts.increases += op == StitchOp::Inc;
            counts.decreases += op == StitchOp::Dec;
        }
    }
    // Parents: the heads in each face, through its bottom edges.
    for (int f = 0; f < nFaces; ++f) {
        std::vector<int>& h = heads[f];
        std::sort(h.begin(), h.end(), [&e](int a, int b) { return e.pieces[a].index < e.pieces[b].index; });
        if (h.size() > 2) return fail("face " + std::to_string(f) + " holds " + std::to_string(h.size()) + " heads");
        for (int id : loopsOfFace[f]) {
            StitchRecord& rec = topo.stitches[id];
            for (int piece : h) {
                const auto it = loopOfTopEdge.find(edgeOf(e.pieces[piece].start));
                if (it == loopOfTopEdge.end()) return fail("the head in face " + std::to_string(f) + " passes no loop's top edge");
                rec.parents[rec.nParents++] = it->second;
            }
            counts.firstTurn += rec.nParents == 0;
        }
        counts.parentLinks += static_cast<int>(h.size());
    }
    if (counts.parentLinks != expectedInterlocks)
        return fail(std::to_string(counts.parentLinks) + " parent links against " + std::to_string(expectedInterlocks) + " interlocks of the yarn check");
    linkChildren(topo.stitches);
    // Ownership on each strand: sorted, non-empty, no overlap.
    std::vector<std::vector<int>> byYarn(curves.yarns.size());
    for (const StitchRecord& rec : topo.stitches) byYarn[rec.yarn].push_back(rec.id);
    std::vector<double> yarnEnd;
    for (size_t y = 0; y < curves.yarns.size(); ++y) {
        std::vector<int>& ids = byYarn[y];
        std::sort(ids.begin(), ids.end(), [&topo](int a, int b) { return topo.stitches[a].uBegin < topo.stitches[b].uBegin; });
        YarnPath& path = curves.yarns[y];
        path.stitches.clear();
        for (size_t k = 0; k < ids.size(); ++k) {
            const StitchRecord& rec = topo.stitches[ids[k]];
            if (!(rec.uBegin < rec.uEnd)) return fail("loop " + std::to_string(rec.id) + " is empty on its strand");
            if (k > 0 && topo.stitches[ids[k - 1]].uEnd > rec.uBegin + 1e-9) return fail("loops " + std::to_string(ids[k - 1]) + " and " + std::to_string(rec.id) + " overlap");
            path.stitches.push_back({rec.index, rec.turn, rec.uBegin, rec.uEnd, rec.id});
        }
        path.finalize(1e-6);
        yarnEnd.push_back(path.curve.uMax());
    }
    topo.buildIndex(yarnEnd);
    curves.topo = std::move(topo);
    return true;
}

} // namespace knit
