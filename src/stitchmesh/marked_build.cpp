// Method A on a marked piece (marked_build.h).
#include "stitchmesh/marked_build.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <map>
#include <memory>
#include <set>
#include <utility>

#include "model/cancel.h"
#include "model/hash.h"
#include "model/mesh_query.h"
#include "stitchmesh/border_placement.h"
#include "stitchmesh/coverage.h"
#include "stitchmesh/guided_rows.h"
#include "stitchmesh/prepared_reference.h"
#include "stitchmesh/quantiles.h"
#include "stitchmesh/relax_mesh.h"
#include "stitchmesh/row_rings.h"
#include "stitchmesh/surface_trim.h"
#include "stitchmesh/time_field.h"
#include "stitchmesh/virtual_shift.h"

namespace knit {

const char* toString(MarkedRegion r) {
    switch (r) {
    case MarkedRegion::Kept: return "kept";
    case MarkedRegion::OmittedStart: return "omitted_start";
    case MarkedRegion::OmittedFinish: return "omitted_finish";
    case MarkedRegion::OutsideBoundary: return "outside_boundary";
    }
    return "kept";
}

uint64_t preparedPieceKey(uint64_t fileDigest, int piece, double mmPerUnit, double h, ReferencePreparation p) {
    uint64_t hash = fnv1a64Bytes(&fileDigest, sizeof fileDigest);
    const double numbers[2] = {mmPerUnit, h};
    hash = fnv1a64Bytes(numbers, sizeof numbers, hash);
    const int ints[2] = {piece, static_cast<int>(p)};
    return fnv1a64Bytes(ints, sizeof ints, hash);
}

void preparePiece(const SurfacePieces& pieces, int piece, uint64_t fileDigest, double mmPerUnit, double h, ReferencePreparation p, PreparedPiece& out) {
    const std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
    TriMesh pieceMm = pieces.pieces[piece].mesh;
    pieceMm.transform(mmPerUnit, Vec3{0, 0, 0});
    out.key = preparedPieceKey(fileDigest, piece, mmPerUnit, h, p);
    out.piece = piece;
    out.refined = prepareReference(pieceMm, refinementLevels(pieceMm, 0.5 * h), p);
    out.nb = triangleNeighbours(out.refined.mesh);
    out.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

namespace {

void addFinding(MarkingReview& r, int markingId, FindingLevel level, FindingCause cause, const std::string& message, int first = -1, int last = -1) {
    MarkingFinding f;
    f.markingId = markingId;
    f.level = level;
    f.cause = cause;
    f.message = message;
    f.firstTracePoint = first;
    f.lastTracePoint = last;
    r.findings.push_back(f);
}

// One warning per bump the rows knitted over (a closed contour beside a row too small for a round of its own), with its
// row and place, on the start marking.
void knittedOverFindings(const std::vector<KnittedOverBump>& bumps, int startId, double toFile, MarkingReview& r) {
    for (const KnittedOverBump& bump : bumps) {
        MarkingFinding f;
        f.markingId = startId;
        f.level = FindingLevel::Warning;
        f.cause = FindingCause::Input;
        char b[240];
        std::snprintf(b, sizeof b, "Row %d passes a bump or dip too small for a round of its own (a closed row of %.1f mm beside it); it is knitted over.", bump.ring,
                      bump.length);
        f.message = b;
        f.row = bump.ring;
        f.at.push_back(bump.at * toFile);
        r.findings.push_back(f);
    }
}

// A problem, or a marking still being drawn: either way no layout is built.
bool hasProblem(const MarkingReview& r) {
    for (const MarkingFinding& f : r.findings)
        if (f.level == FindingLevel::Problem || f.level == FindingLevel::Incomplete) return true;
    return false;
}

// "3, 4, 7-9": ring numbers for a message.
std::string ringList(const std::vector<int>& k) {
    std::string s;
    for (size_t i = 0; i < k.size();) {
        size_t j = i;
        while (j + 1 < k.size() && k[j + 1] == k[j] + 1) ++j;
        if (!s.empty()) s += ", ";
        s += std::to_string(k[i]);
        if (j > i) s += "-" + std::to_string(k[j]);
        i = j + 1;
    }
    return s;
}

std::vector<Vec3> scaled(const std::vector<Vec3>& p, double s) {
    std::vector<Vec3> out;
    out.reserve(p.size());
    for (const Vec3& x : p) out.push_back(x * s);
    return out;
}

double triangleArea(const TriMesh& m, int t) {
    const Vec3& a = m.pos[m.tri[t][0]];
    return 0.5 * norm(cross(m.pos[m.tri[t][1]] - a, m.pos[m.tri[t][2]] - a));
}

int eulerCharacteristic(const TriMesh& m) {
    std::set<std::pair<int, int>> edges;
    for (const std::array<int, 3>& t : m.tri)
        for (int k = 0; k < 3; ++k) edges.insert({std::min(t[k], t[(k + 1) % 3]), std::max(t[k], t[(k + 1) % 3])});
    return static_cast<int>(m.pos.size()) - static_cast<int>(edges.size()) + static_cast<int>(m.tri.size());
}

SurfacePoint locate(const MeshQuery& q, const TriMesh& m, const Vec3& p) {
    SurfacePoint s;
    const MeshQuery::Hit h = q.closest(p);
    if (h.tri < 0) return s;
    const TrianglePoint tp = closestOnTriangle(p, m.pos[m.tri[h.tri][0]], m.pos[m.tri[h.tri][1]], m.pos[m.tri[h.tri][2]]);
    s.tri = h.tri;
    for (int k = 0; k < 3; ++k) s.bary[k] = tp.bary[k];
    return s;
}

double valueAt(const TriMesh& m, const std::vector<double>& t, const SurfacePoint& p) {
    const auto& tr = m.tri[p.tri];
    return p.bary[0] * t[tr[0]] + p.bary[1] * t[tr[1]] + p.bary[2] * t[tr[2]];
}

// The row field from a point inside a triangle: its corners seeded with their distances from the point.
bool fieldFromPoint(const TriMesh& m, const SurfacePoint& p, std::vector<double>& t, std::string* err) {
    const Vec3 x = surfacePosition(m, p);
    std::vector<int> seeds;
    std::vector<double> values;
    for (int k = 0; k < 3; ++k) {
        seeds.push_back(m.tri[p.tri][k]);
        values.push_back(norm(m.pos[m.tri[p.tri][k]] - x));
    }
    TimeFieldReport rep;
    return geodesicDistance(m, seeds, t, rep, err, &values);
}

// The other closed pieces of the same component (mm): a point inside one is hidden.
struct HiddenTest {
    std::vector<TriMesh> meshes;
    std::vector<std::unique_ptr<MeshQuery>> queries;
    bool inside(const Vec3& p) const {
        for (const std::unique_ptr<MeshQuery>& q : queries)
            if (q->inside(p)) return true;
        return false;
    }
};

void hiddenTest(const SurfacePieces& pieces, int piece, double mm, HiddenTest& out) {
    for (int i = 0; i < static_cast<int>(pieces.pieces.size()); ++i) {
        if (i == piece || pieces.pieces[i].node != pieces.pieces[piece].node) continue;
        TriMesh m = pieces.pieces[i].mesh;
        std::set<std::pair<int, int>> once, twice;
        for (const auto& t : m.tri)
            for (int k = 0; k < 3; ++k) {
                const std::pair<int, int> e{std::min(t[k], t[(k + 1) % 3]), std::max(t[k], t[(k + 1) % 3])};
                if (!once.insert(e).second) twice.insert(e);
            }
        if (once.size() != twice.size()) continue;   // open: no inside
        m.transform(mm, Vec3{0, 0, 0});
        out.meshes.push_back(m);
    }
    for (const TriMesh& m : out.meshes) out.queries.push_back(std::make_unique<MeshQuery>(m));
}

// The row-start path as the guide of the marked builder: positions, field values, the triangle of each segment.
std::vector<MarkedPathPoint> pathOnField(const MarkingTrace& tr, const TriMesh& onMesh, const MeshQuery& q, const std::vector<double>& t) {
    std::vector<MarkedPathPoint> out;
    for (size_t i = 0; i < tr.points.size(); ++i) {
        MarkedPathPoint p;
        p.p = tr.points[i].pos;
        const SurfacePoint s = locate(q, onMesh, p.p);
        p.t = s.tri >= 0 ? valueAt(onMesh, t, s) : 0;
        if (i + 1 < tr.points.size()) p.segmentTri = locate(q, onMesh, 0.5 * (tr.points[i].pos + tr.points[i + 1].pos)).tri;
        out.push_back(p);
    }
    return out;
}

// The row-start path's trace as the rows builder's guide: each segment in the triangle the trace names; parentTri maps the
// rows' mesh triangles to the traced mesh's (null: the rows run on the traced mesh itself).
SurfacePathGuide guideFromTrace(const MarkingTrace& path, const std::vector<int>* parentTri) {
    SurfacePathGuide guide;
    for (const TracePoint& p : path.points) {
        guide.points.push_back(p.pos);
        guide.segmentTri.push_back(p.tri);
    }
    guide.parentTri = parentTri;
    return guide;
}

GuideLine guideFromPath(const MarkingTrace& tr, const TriMesh& onMesh, const MeshQuery& q) {
    GuideLine g;
    for (size_t i = 0; i < tr.points.size(); ++i) {
        g.points.push_back(tr.points[i].pos);
        if (i + 1 < tr.points.size()) g.segmentTri.push_back(locate(q, onMesh, 0.5 * (tr.points[i].pos + tr.points[i + 1].pos)).tri);
    }
    return g;
}

double distanceToPolyline(const Vec3& x, const std::vector<Vec3>& line) {
    double best = 1e300;
    for (size_t i = 0; i + 1 < line.size(); ++i) best = std::min(best, closestSegmentPoints(x, x, line[i], line[i + 1]).distance);
    return best;
}

// A closed marking on a border loop of the piece (pieceLoop: piece vertex ids, walked as the control points go), on the
// prepared reference: refinement appends midpoints, so the piece's border vertices keep their ids and each piece border
// edge becomes 2^levels reference edges. refLoop: the reference's loop through them, walked the same way; refPoints: each
// control point on its reference edge with exact weights (for borderLoopTrace). False when the loop is not found.
bool borderOnReference(const RefinedMesh& refined, const TriMesh& pieceMesh, const std::vector<int>& pieceLoop, const std::vector<SurfacePoint>& pieceLocal,
                       std::vector<int>& refLoop, std::vector<SurfacePoint>& refPoints) {
    refLoop.clear();
    refPoints.clear();
    const TriMesh& m = refined.mesh;
    for (const std::vector<int>& loop : borderLoops(m))
        if (std::find(loop.begin(), loop.end(), pieceLoop[0]) != loop.end()) refLoop = loop;
    if (refLoop.empty() || pieceLoop.size() < 2) return false;
    const int n = static_cast<int>(refLoop.size()), steps = 1 << refined.levels;
    auto at = [&](int v) { return static_cast<int>(std::find(refLoop.begin(), refLoop.end(), v) - refLoop.begin()); };
    if (at(pieceLoop[1]) != (at(pieceLoop[0]) + steps) % n) std::reverse(refLoop.begin(), refLoop.end());
    if (at(pieceLoop[1]) != (at(pieceLoop[0]) + steps) % n) return false;
    const int pn = static_cast<int>(pieceLoop.size());
    auto onEdge = [&](int a, int b, SurfacePoint& out) {
        for (int t = 0; t < static_cast<int>(m.tri.size()); ++t)
            for (int k = 0; k < 3; ++k)
                if ((m.tri[t][k] == a && m.tri[t][(k + 1) % 3] == b) || (m.tri[t][k] == b && m.tri[t][(k + 1) % 3] == a)) {
                    out.tri = t;
                    out.bary[0] = out.bary[1] = out.bary[2] = 0;
                    return k;
                }
        return -1;
    };
    for (const SurfacePoint& p : pieceLocal) {
        // the piece edge (a -> b in the loop's walk) the point lies on, and its fraction from a
        const std::array<int, 3>& tr = pieceMesh.tri[p.tri];
        int a = -1;
        double f = 0;
        for (int k = 0; k < 3 && a < 0; ++k)
            if (p.bary[k] == 1) a = tr[k];
        for (int z = 0; z < 3 && a < 0; ++z) {
            if (p.bary[z] != 0) continue;
            const int u = tr[(z + 1) % 3], v = tr[(z + 2) % 3];
            const int iu = static_cast<int>(std::find(pieceLoop.begin(), pieceLoop.end(), u) - pieceLoop.begin());
            const int iv = static_cast<int>(std::find(pieceLoop.begin(), pieceLoop.end(), v) - pieceLoop.begin());
            if (iu < pn && iv < pn && (iu + 1) % pn == iv) {
                a = u;
                f = p.bary[(z + 2) % 3];
            } else if (iu < pn && iv < pn && (iv + 1) % pn == iu) {
                a = v;
                f = p.bary[(z + 1) % 3];
            }
        }
        if (a < 0) return false;
        const double sub = f * steps;
        const int s = static_cast<int>(sub);
        const double g = sub - s;
        const int u = refLoop[(at(a) + s) % n], v = refLoop[(at(a) + s + 1) % n];
        SurfacePoint q;
        const int k = onEdge(u, v, q);
        if (k < 0) return false;
        const bool forward = m.tri[q.tri][k] == u;
        if (g == 0) q.bary[forward ? k : (k + 1) % 3] = 1;   // on the vertex u
        else {
            q.bary[forward ? k : (k + 1) % 3] = 1 - g;
            q.bary[forward ? (k + 1) % 3 : k] = g;
        }
        refPoints.push_back(q);
    }
    return true;
}

// The region of a trimmed surface holding the edge (a, b) of the reference (a border edge is never cut: one trim triangle has it).
int regionOfEdge(const TrimmedSurface& trim, int a, int b) {
    for (int t = 0; t < static_cast<int>(trim.mesh.tri.size()); ++t)
        for (int k = 0; k < 3; ++k) {
            const int x = trim.mesh.tri[t][k], y = trim.mesh.tri[t][(k + 1) % 3];
            if ((x == a && y == b) || (x == b && y == a)) return trim.region[t];
        }
    return -1;
}

// Everything after the stored tiling exists: the virtual shift, the checks, the relaxation with the markings'
// constraints, coverage of the knit region, the review's measures. surface: where vertices walk during relaxation.
struct Tail {
    const MarkedBuildRequest* r = nullptr;
    MarkedBuildProduct* out = nullptr;
    TypedStitchMesh s;                          // mm
    const TriMesh* surface = nullptr;
    std::vector<std::vector<int>> surfaceLoops; // border loops of surface the bind-off / cast-on slide on (band, pole-curve)
    RelaxGuides guides;                         // pole route / anchoring
    std::vector<Vec3> pathLine;                 // the row-start path (row-origin measure)
    int pole = -1;                              // stitch vertex of the pole (pole route)
    int finishId = 0;
    std::vector<DeclaredBoundary> declared;
    const TriMesh* relaxTo = nullptr;           // controlled comparison: another preparation of the reference, the same split
    std::vector<GuideLine> relaxToLines;        // ... its constraint lines derived again there (empty: carried)
    std::vector<Vec3> relaxToPath;              // ... the row-start path traced there (row-origin measure)
    TriMesh relaxToRegion;                      // ... the knit region clipped there at the same level (empty: carried)
};

// The controlled comparison (MarkedBuildOptions::relaxOnOther): everything the layout placed on the reference keeps its
// refined triangle and weights and takes that point of the other preparation. The split is the same, so a trimmed
// region's triangle (inside one reference triangle) stays inside the same triangle of the other preparation.
void carryToOtherPreparation(Tail& T, TriMesh& carriedSurface) {
    MarkedEvaluation& ev = T.out->evaluation;
    const TriMesh& from = ev.reference;
    const TriMesh& to = *T.relaxTo;
    const MeshQuery q(from);
    auto carry = [&](Vec3& x) {
        const SurfacePoint p = locate(q, from, x);
        if (p.tri >= 0) x = surfacePosition(to, p);
    };
    for (Vec3& x : T.s.mesh.pos) carry(x);
    for (VertexBinding& vb : T.s.binding) vb.distance = 0;
    if (T.surface != &ev.reference) {
        carriedSurface = *T.surface;
        for (Vec3& x : carriedSurface.pos) carry(x);
        carriedSurface.computeNormals();
        T.surface = &carriedSurface;
    }
    for (GuideLine& g : T.guides.lines)
        for (Vec3& x : g.points) carry(x);
    for (Vec3& x : ev.knitRegion.pos) carry(x);
    ev.knitRegion.computeNormals();
    for (DeclaredBoundary& d : T.declared)
        for (Vec3& x : d.loop) carry(x);
    for (Vec3& x : T.pathLine) carry(x);
    if (!T.relaxToLines.empty() && T.relaxToLines.size() == T.guides.lines.size()) {
        T.guides.lines = T.relaxToLines;
        T.pathLine = T.relaxToPath;
        if (!T.relaxToRegion.tri.empty()) {
            ev.knitRegion = T.relaxToRegion;
            ev.regionLoops = borderLoops(ev.knitRegion);
            if (!T.declared.empty() && !T.guides.lines.empty()) T.declared.back().loop = T.guides.lines.front().points;   // the bind-off loop
        }
    }
    ev.reference = to;   // a walk on the whole reference now walks on the other preparation
}

// BorderRings::PlacedOnLoop: the first and the last ring, where all their vertices lie on the stitch mesh's border, placed
// on the nearest boundary loop of the knit region (placeBorderRing) and marked as held. One log line per ring.
void placeBorderRings(TypedStitchMesh& s, const MarkedEvaluation& ev, double limit, double anchorFreedom, std::vector<char>& held, std::vector<std::string>& log) {
    if (s.rings.empty() || ev.regionLoops.empty()) return;
    std::vector<const std::vector<int>*> rings = {&s.rings.back()};
    if (s.rings.size() > 1) rings.push_back(&s.rings.front());
    for (const std::vector<int>* ring : rings) {
        bool onBorder = ring->size() >= 3;
        for (int v : *ring) onBorder = onBorder && v < static_cast<int>(ev.border.size()) && ev.border[v];
        if (!onBorder) continue;
        std::vector<Vec3> corners;
        for (int v : *ring) corners.push_back(s.mesh.pos[v]);
        // the loop the ring lies on: the one nearest to its vertices
        int nearest = -1;
        double nearestSum = 1e300;
        for (size_t l = 0; l < ev.regionLoops.size(); ++l) {
            double sum = 0;
            for (const Vec3& c : corners) {
                double d = 1e300;
                for (int v : ev.regionLoops[l]) d = std::min(d, norm(ev.knitRegion.pos[v] - c));
                sum += d;
            }
            if (sum < nearestSum) {
                nearestSum = sum;
                nearest = static_cast<int>(l);
            }
        }
        std::vector<Vec3> loop;
        for (int v : ev.regionLoops[nearest]) loop.push_back(ev.knitRegion.pos[v]);
        BorderRingRequest request;
        request.limit = limit;
        request.anchorFreedom = anchorFreedom;
        const BorderRingPlacement placed = placeBorderRing(loop, corners, request);
        char b[300];
        if (!placed.placed) {
            std::snprintf(b, sizeof b, "border ring of %zu vertices: no placement on its loop found, left to the relaxation (distance %.3f mm)", ring->size(), placed.distanceBefore);
            log.push_back(b);
            continue;
        }
        for (size_t k = 0; k < ring->size(); ++k) {
            s.mesh.pos[(*ring)[k]] = placed.corners[k];
            held[(*ring)[k]] = 1;
        }
        std::snprintf(b, sizeof b, "border ring of %zu vertices placed on its loop: distance to the loop %.3f -> %.3f mm%s, sides %.2f to %.2f mm (within %.0f %% of the mean)",
                      ring->size(), placed.distanceBefore, placed.distance, placed.inOrder ? "" : " (the relaxed ring was out of order)", placed.shortestSide,
                      placed.longestSide, 100 * placed.spread);
        log.push_back(b);
    }
}

void finishBuild(Tail& T) {
    const MarkedBuildRequest& r = *T.r;
    MarkedBuildProduct& out = *T.out;
    MarkingReview& review = out.review;
    MarkedEvaluation& ev = out.evaluation;
    TypedStitchMesh& s = T.s;
    char b[400];
    TriMesh carriedSurface;
    if (T.relaxTo) carryToOtherPreparation(T, carriedSurface);
    TypedCheckOptions o;
    o.w = r.w;
    o.h = r.h;
    o.requireStrips = true;
    o.faceDeviation = 0.1 * r.w;
    o.bindingTolerance = 0.001 * r.w;
    o.boundaryTolerance = 0.05;
    TypedCheckOptions ov = o;
    ov.requireStrips = false;
    CoverageOptions co;
    co.coverDistance = 0.1 * r.w;
    co.boundaryBand = 0.5 * r.h;
    co.hausdorffLimit = 0.25 * r.h;
    ev.checkOptions = o;
    ev.coverOptions = co;
    std::string shiftErr;
    const bool shifted = buildVirtualShift(s, ShiftReading::PathSideDown, ev.shift, &shiftErr);
    ev.before = s;
    auto virtualOf = [&](const TypedStitchMesh& x) {
        TypedStitchMesh v = virtualTypedMesh(x, ev.shift);
        v.strip.clear();   // the count invariant is the stored mesh's
        v.slot.clear();
        return v;
    };
    ev.storedBefore = checkTypedMesh(s, ev.reference, T.declared, o);
    cancelCheckpoint();
    ev.coverStoredBefore = measureCoverage(ev.knitRegion, ev.regionLoops, s.mesh, co);
    if (shifted) {
        ev.shiftBefore = validateVirtualShift(s, ev.shift);
        const TypedStitchMesh v = virtualOf(s);
        ev.virtualBefore = checkTypedMesh(v, ev.reference, T.declared, ov);
        cancelCheckpoint();
        ev.coverVirtualBefore = measureCoverage(ev.knitRegion, ev.regionLoops, v.mesh, co);
    }
    cancelCheckpoint();
    if (r.depth == ReviewDepth::Relaxed) {
        for (int f = 0; f < s.mesh.faces(); ++f) {
            std::vector<int> corners;
            for (int k = 0; k < s.mesh.sides(f); ++k) corners.push_back(s.mesh.vertexAt(f, k));
            ev.solverFaces.push_back(RelaxFace{shifted ? ev.shift.corners[f] : corners, s.role[f]});
        }
        MeshRelaxRequest rq;  // the frozen CP3 settings, scaled by the gauge
        rq.w = r.w;
        rq.h = r.h;
        rq.smoothNormals = r.options.smoothFilterNormals;
        rq.kinkConstraints = r.options.kinkConstraints;
        rq.borderOrder = r.options.borderOrder;
        rq.faceCorners = r.options.faceCorners;
        rq.stepControl = r.options.stepControl;
        if (r.options.borderSprings == BorderSprings::StoredBorderEdges) {
            // the stored mesh's border edges (one stored face each); the relaxation adds a spring where no solver face has one
            std::map<std::pair<int, int>, int> uses;
            for (int f = 0; f < s.mesh.faces(); ++f)
                for (int k = 0; k < s.mesh.sides(f); ++k) {
                    const int a = s.mesh.vertexAt(f, k), c = s.mesh.vertexAt(f, (k + 1) % s.mesh.sides(f));
                    ++uses[{std::min(a, c), std::max(a, c)}];
                }
            for (const auto& [e, n] : uses)
                if (n == 1) T.guides.courseEdges.push_back({e.first, e.second});
        }
        std::string relaxErr;
        if (r.options.borderRings == BorderRings::PlacedFirst && r.options.anchors != RelaxAnchors::RowOriginsOnPath && !T.relaxTo) {
            T.guides.pinned.resize(s.mesh.pos.size(), 0);
            placeBorderRings(s, ev, co.hausdorffLimit, 0.5 * r.w, T.guides.pinned, out.log);
        }
        const bool useGuides = !T.guides.lines.empty() || T.pole >= 0 || !T.guides.courseEdges.empty();
        if (relaxStitchMesh(*T.surface, T.surfaceLoops, ev.solverFaces, s.mesh.pos, ev.border, rq, ev.relax, &relaxErr, useGuides ? &T.guides : nullptr)) {
            ev.relaxed = true;
            if (T.pole >= 0) ev.poleMoved = norm(ev.relax.pos[T.pole] - s.mesh.pos[T.pole]);
            const std::vector<Vec3> laidOut = s.mesh.pos;
            s.mesh.pos = ev.relax.pos;
            // BorderRings::PlacedOnLoop: the relaxed rows have found their places around the piece; now the border rings are
            // placed on their loops for the fit, held, and the other vertices relaxed once more around them
            if (r.options.borderRings == BorderRings::PlacedOnLoop && r.options.anchors != RelaxAnchors::RowOriginsOnPath && !T.relaxTo) {
                RelaxGuides held = T.guides;
                held.pinned.resize(s.mesh.pos.size(), 0);
                const std::vector<char> heldBefore = held.pinned;
                placeBorderRings(s, ev, co.hausdorffLimit, r.options.anchors == RelaxAnchors::YarnEnds ? 0.0 : 0.5 * r.w, held.pinned, out.log);
                if (held.pinned != heldBefore) {
                    RelaxProduct second;
                    const bool ran = relaxStitchMesh(*T.surface, T.surfaceLoops, ev.solverFaces, s.mesh.pos, ev.border, rq, second, &relaxErr, &held);
                    const bool kept = ran && second.status != Status::Fail && second.rmsForceEnd < second.rmsForceStart;
                    std::snprintf(b, sizeof b, "relaxation around the held border rings: %s, RMS force %.3g -> %.3g, %d steps: %s", ran ? second.stopReason.c_str() : "did not run",
                                  second.rmsForceStart, second.rmsForceEnd, second.steps, kept ? "kept" : "dropped (the forces did not go down)");
                    out.log.push_back(b);
                    if (kept) {
                        second.steps += ev.relax.steps;
                        second.ms += ev.relax.ms;
                        second.rmsForceStart = ev.relax.rmsForceStart;
                        ev.relax = second;
                        s.mesh.pos = ev.relax.pos;
                    }
                }
            }
            if (T.pole >= 0) ev.poleMoved = norm(s.mesh.pos[T.pole] - laidOut[T.pole]);
            // bindings on the prepared reference (the walk ran on `surface`, which may be a trimmed copy of it)
            const MeshQuery q(ev.reference);
            s.binding.assign(s.mesh.pos.size(), VertexBinding{});
            for (size_t v = 0; v < s.mesh.pos.size(); ++v) {
                const SurfacePoint p = locate(q, ev.reference, s.mesh.pos[v]);
                s.binding[v].tri = p.tri;
                for (int k = 0; k < 3; ++k) s.binding[v].bary[k] = p.bary[k];
                s.binding[v].distance = p.tri >= 0 ? norm(surfacePosition(ev.reference, p) - s.mesh.pos[v]) : 0;
            }
            if (ev.relax.status != Status::Pass)
                addFinding(review, 0, FindingLevel::Warning, FindingCause::SolverLimitation,
                           "The relaxation stopped without converging (" + ev.relax.stopReason + "); the solver's limit, not the markings'.");
            std::snprintf(b, sizeof b, "relaxation %s: %s, %d steps, %.0f ms, filter normals %s, kinks %s (%d on a crease, %d held at a corner), %d pinned, %d guided",
                          toString(ev.relax.status), ev.relax.stopReason.c_str(), ev.relax.steps, ev.relax.ms, rq.smoothNormals ? "smooth (R1)" : "facet",
                          rq.kinkConstraints ? "held (R2)" : "free", ev.relax.kinkOnEdgeLast, ev.relax.kinkHeldLast, ev.relax.pinnedVertices, ev.relax.guidedVertices);
            out.log.push_back(b);
            if (!ev.relax.firstCornerLimit.empty()) out.log.push_back("first corner step limited: " + ev.relax.firstCornerLimit);
            std::snprintf(b, sizeof b, "steps taken back: %d", ev.relax.stepsTakenBack);
            out.log.push_back(b);
            std::snprintf(b, sizeof b, "border springs added: %d; corner steps limited: %d; border order: %d steps along a loop limited, %d order violations%s%s",
                          ev.relax.addedCourseSprings, ev.relax.cornerMovesLimited, ev.relax.borderMovesLimited, ev.relax.borderOrderViolations,
                          ev.relax.borderOrderViolations ? "; first: " : "", ev.relax.firstBorderOrderViolation.c_str());
            out.log.push_back(b);
        } else {
            addFinding(review, 0, FindingLevel::Warning, FindingCause::SolverLimitation, "The relaxation did not run: " + relaxErr);
        }
        cancelCheckpoint();
    }
    ev.after = s;
    ev.storedAfter = checkTypedMesh(s, ev.reference, T.declared, o);
    if (ev.relaxPreparation != "facets") {
        ev.storedAfterVsOriginal = checkTypedMesh(s, ev.original, T.declared, o);
        std::snprintf(b, sizeof b, "Against the imported surface (not the reference the cells are judged on), the stitch faces deviate by up to %.2f mm (%d faces over %.2f mm).",
                      ev.storedAfterVsOriginal.maxFaceDeviation, ev.storedAfterVsOriginal.facesOverDeviation, o.faceDeviation);
        addFinding(review, 0, FindingLevel::Ok, FindingCause::Input, b);
        out.log.push_back(b);
    }
    cancelCheckpoint();
    ev.coverStoredAfter = measureCoverage(ev.knitRegion, ev.regionLoops, s.mesh, co);
    if (shifted) {
        cancelCheckpoint();
        ev.shiftAfter = validateVirtualShift(s, ev.shift);
        const TypedStitchMesh v = virtualOf(s);
        ev.virtualAfter = checkTypedMesh(v, ev.reference, T.declared, ov);
        ev.coverVirtualAfter = measureCoverage(ev.knitRegion, ev.regionLoops, v.mesh, co);
    }
    if (!T.pathLine.empty() && !s.rings.empty()) {
        double sum = 0;
        int n = 0;
        for (size_t k = 1; k < s.rings.size(); ++k) {
            const double d = distanceToPolyline(s.mesh.pos[s.rings[k][0]], T.pathLine);
            sum += d;
            ++n;
            ev.rowOriginMax = std::max(ev.rowOriginMax, d);
        }
        ev.rowOriginMean = n ? sum / n : 0;
    }
    // the review
    int unsupported = 0, exceptional = 0;
    std::set<std::string> reasons;
    for (int f = 0; f < s.mesh.faces(); ++f) {
        if (s.kind[f] != StitchFaceKind::Regular) ++exceptional;
        if (s.kind[f] == StitchFaceKind::Unsupported) {
            ++unsupported;
            reasons.insert(s.reason[f]);
        }
    }
    review.measures.unsupportedCells = unsupported;
    review.measures.exceptionalCells = exceptional;
    if (unsupported > 0) {
        std::string list;
        for (const std::string& x : reasons) list += (list.empty() ? "" : ", ") + x;
        std::snprintf(b, sizeof b, "%d %s no construction yet (%s).", unsupported, unsupported == 1 ? "cell has" : "cells have", list.c_str());
        addFinding(review, 0, FindingLevel::Warning, FindingCause::AlgorithmLimitation, b);
    }
    std::vector<double> distortion;
    for (const CellMeasure& c : measureCells(s, r.w, r.h)) {
        if (c.kind != StitchFaceKind::Regular) continue;
        distortion.push_back(std::fabs(c.widthRatio - 1));
        distortion.push_back(std::fabs(c.heightRatio - 1));
    }
    review.measures.distortionP95 = quantiles(distortion).p95;
    const CoverageReport& cov = ev.coverStoredAfter;
    review.measures.coveredFraction = cov.uncoveredMeasured ? 1 - cov.uncoveredFraction : 0;
    review.measures.coveragePass = cov.status == Status::Pass;
    out.log.push_back("coverage of the stored tiling: " + cov.text());
    if (shifted) out.log.push_back("coverage of the virtual (solver) faces: " + ev.coverVirtualAfter.text());
    if (!review.measures.coveragePass)
        addFinding(review, T.finishId, FindingLevel::Warning, FindingCause::AlgorithmLimitation,
                   std::string(cov.status == Status::Fail ? "The coverage check fails: " : "The coverage check is unresolved: ") + cov.detail);
    if (ev.storedAfter.geometry == Status::Fail) {
        const std::string text = ev.storedAfter.text();
        std::string line = text.substr(text.find("geometry"));
        line = line.substr(0, line.find('\n'));
        addFinding(review, 0, FindingLevel::Warning, FindingCause::AlgorithmLimitation, "Geometry: " + line);
    }
    if (ev.storedAfter.connectivity == Status::Fail)
        addFinding(review, 0, FindingLevel::Warning, FindingCause::AlgorithmLimitation, "Connectivity: " + ev.storedAfter.detail);
    TypedStitchMesh fileUnits = s;
    fileUnits.mesh.pos = scaled(s.mesh.pos, 1.0 / r.mmPerUnit);
    out.stitchMesh = fileUnits;
    review.valid = true;
    // a built mesh is not yet a layout that passes: its own checks say that
    LayoutMeasures& ms = review.measures;
    ms.checksPass = ev.storedAfter.connectivity == Status::Pass && ev.storedAfter.geometry == Status::Pass && ms.coveragePass && ms.unsupportedCells == 0 &&
                    ms.contourEvents == 0 && !hasProblem(review);
    ms.unsupportedKnown = ms.exceptionalKnown = ms.distortionKnown = ms.coverageKnown = ms.checksKnown = Known::Current;
}

// A knitting component that lies over another (a collar over a body) is laid out on its own copy of the reference, moved
// out along the vertex normals by `offsetMm`, so its yarn clears the yarn underneath. Only positions change: the markings'
// refined triangles and weights, the parents and the neighbours stay valid.
std::shared_ptr<const PreparedPiece> offsetReference(const PreparedPiece& prepared, double offsetMm) {
    std::shared_ptr<PreparedPiece> moved = std::make_shared<PreparedPiece>(prepared);
    TriMesh& mesh = moved->refined.mesh;
    mesh.computeNormals();   // area-weighted; triangles are counter-clockwise seen from outside
    for (size_t v = 0; v < mesh.pos.size(); ++v) mesh.pos[v] = mesh.pos[v] + mesh.vnormal[v] * offsetMm;
    mesh.computeNormals();
    moved->key = 0;   // no longer the prepared piece of preparedPieceKey: never shared with a quick review
    return moved;
}

} // namespace

void resolvePieceMarkings(const MarkedBuildRequest& r, const RefinedMesh& refined, const TriangleNeighbours& nb, MarkingReview& review, ResolvedMarkings& out) {
    out = ResolvedMarkings{};
    const TriMesh& mm = refined.mesh;
    const double toFile = 1.0 / r.mmPerUnit;
    char b[400];
    MeaningChoice& choice = out.choice;
    for (const SurfaceMarking& m : r.markings->markings) {
        if (!markingInKnitComponent(*r.pieces, m, r.piece, r.knitComponent)) continue;   // another piece, or another knitting component
        if (m.origin == MarkingOrigin::Proposal && m.proposal == ProposalState::Rejected) continue;  // a rejected proposal is never input
        const bool preview = m.origin == MarkingOrigin::Proposal && m.proposal == ProposalState::Offered &&
                             std::find(r.previewOffered.begin(), r.previewOffered.end(), m.id) != r.previewOffered.end();
        if (m.origin == MarkingOrigin::Proposal && m.proposal == ProposalState::Offered && !preview) {
            ++out.offeredSkipped;   // an offered proposal is input only once accepted (or named for a research preview)
            continue;
        }
        ResolvedMarking u;
        u.marking = &m;
        u.preview = preview;
        u.resolved = true;
        std::vector<SurfacePoint> pieceLocal;   // on the piece's own triangles (exact weights)
        for (const SurfacePoint& p : m.points) {
            int pi = -1;
            SurfacePoint lp;
            if (!pieceLocalPoint(*r.pieces, p, pi, lp) || pi != r.piece) u.resolved = false;
            else {
                pieceLocal.push_back(lp);
                u.local.push_back(toRefined(refined, lp));
            }
        }
        const int index = static_cast<int>(out.used.size());
        out.markingIds.push_back(m.id);
        if (preview) out.previewIds.push_back(m.id);
        // a closed start or termination whose points all lie on one border loop of the piece is that opening itself
        std::vector<int> pieceLoop, refLoop;
        std::vector<SurfacePoint> refPoints;
        const bool opening = u.resolved && m.shape == MarkingShape::ClosedCurve &&
                             (markingRole(m.meaning) == MarkingRole::Start || markingRole(m.meaning) == MarkingRole::Termination) &&
                             pointsOnOneBorderLoop(r.pieces->pieces[r.piece].mesh, pieceLocal, pieceLoop) &&
                             borderOnReference(refined, r.pieces->pieces[r.piece].mesh, pieceLoop, pieceLocal, refLoop, refPoints);
        if (opening) {
            u.borderLoop = refLoop;
            u.trace = borderLoopTrace(mm, refLoop, refPoints);
        } else {
            u.trace = traceMarking(mm, nb, u.local, m.shape == MarkingShape::ClosedCurve);
        }
        MarkingTrace pieceTrace = u.trace;   // file units, on the piece's own triangles (the viewer's contract)
        for (TracePoint& p : pieceTrace.points) {
            p.pos = p.pos * toFile;
            if (p.tri >= 0) p.tri = refined.parent[p.tri];
            p.edgeFrom = p.edgeTo = -1;
        }
        out.traces.push_back(pieceTrace);
        out.used.push_back(u);
        const MarkingTrace& tr = out.used.back().trace;
        if (!u.resolved) {
            addFinding(review, m.id, FindingLevel::Problem, FindingCause::Input, "Some points of this marking lie on another piece or a degenerate triangle; keep a marking on one piece.");
            continue;
        }
        if (!shapeAllowed(m.meaning, m.shape)) {
            std::snprintf(b, sizeof b, "A %s cannot be a %s.", meaningLabel(m.meaning), toString(m.shape));
            addFinding(review, m.id, FindingLevel::Problem, FindingCause::Marking, b);
            continue;
        }
        if (!markingComplete(m)) {
            // still being drawn: incomplete, not wrong; its role counts as taken, so no "add a start" beside it
            addFinding(review, m.id, FindingLevel::Incomplete, FindingCause::Marking, "Add more points: a path needs 2, a closed curve 3.");
            if (markingRole(m.meaning) == MarkingRole::Start) choice.incompleteStart = true;
            if (markingRole(m.meaning) == MarkingRole::Termination) choice.incompleteFinish = true;
            if (markingRole(m.meaning) == MarkingRole::RowStart) choice.incompleteRowStart = true;
            continue;
        }
        for (int k : tr.failedSegments) {
            const int first = tr.controlAt[k];
            const int last = k + 1 < static_cast<int>(tr.controlAt.size()) ? tr.controlAt[k + 1] : static_cast<int>(tr.points.size()) - 1;
            std::snprintf(b, sizeof b, "Points %d and %d cannot be joined across the surface.", k + 1, k + 2);
            addFinding(review, m.id, FindingLevel::Problem, FindingCause::Marking, b, first, last);
        }
        if (tr.edgePathUsed && m.shape == MarkingShape::ClosedCurve)
            addFinding(review, m.id, FindingLevel::Problem, FindingCause::Marking,
                       "Part of this curve follows mesh edges and cannot trim the surface; add a point where it bends.");
        for (const auto& [k, v] : m.settings)
            if (!(m.meaning == MarkingMeaning::Start && m.shape == MarkingShape::Point && k == "ring1"))
                addFinding(review, m.id, FindingLevel::Warning, FindingCause::Input, "Setting \"" + k + "\" is not used by this marking's meaning and shape; it is kept, not applied.");
        switch (m.meaning) {
        case MarkingMeaning::Start:
            if (choice.start >= 0) addFinding(review, m.id, FindingLevel::Problem, FindingCause::Marking, "Only one start per knitting component; delete the other, or give it its own component.");
            else choice.start = index;
            if (m.shape == MarkingShape::Point)
                addFinding(review, m.id, FindingLevel::Warning, FindingCause::AlgorithmLimitation,
                           "A start point uses the gathered pole start, an experimental construction (contract and fixtures exist; not yet checked at yarn level).");
            break;
        case MarkingMeaning::Finish:
        case MarkingMeaning::Boundary:
            if (choice.finish >= 0) {
                addFinding(review, m.id, FindingLevel::Problem, FindingCause::Marking, "Only one finish or boundary per knitting component; delete the other, or give it its own component.");
            } else {
                choice.finish = index;
                choice.finishIsBoundary = m.meaning == MarkingMeaning::Boundary;
            }
            break;
        case MarkingMeaning::RowStartPath:
            if (choice.rowStart >= 0) addFinding(review, m.id, FindingLevel::Problem, FindingCause::Marking, "Only one row-start path per knitting component; delete the other, or give it its own component.");
            else choice.rowStart = index;
            break;
        case MarkingMeaning::Cut:
            addFinding(review, m.id, FindingLevel::Warning, FindingCause::AlgorithmLimitation,
                       "Cuts are kept and drawn, but the mesher does not split pieces into panels yet.");
            break;
        case MarkingMeaning::Attachment:
            addFinding(review, m.id, FindingLevel::Warning, FindingCause::AlgorithmLimitation,
                       "Attachments are kept for joining components later; the mesher does not use them yet, and they establish no join.");
            break;
        }
    }
}

void reviewMeaningChoice(const ResolvedMarkings& m, MarkingReview& review) {
    char b[400];
    const MeaningChoice& choice = m.choice;
    if (!m.previewIds.empty()) {
        std::string ids;
        for (int id : m.previewIds) ids += (ids.empty() ? "" : ", ") + std::to_string(id);
        addFinding(review, 0, FindingLevel::Warning, FindingCause::Input,
                   "Research preview: offered proposals " + ids + " are evaluated as if accepted. They are not accepted; accept or reject them in the Mark step.");
    }
    if (m.offeredSkipped > 0) {
        std::snprintf(b, sizeof b, "%d offered %s on this piece %s not input until accepted.", m.offeredSkipped, m.offeredSkipped == 1 ? "proposal" : "proposals",
                      m.offeredSkipped == 1 ? "is" : "are");
        addFinding(review, 0, FindingLevel::Ok, FindingCause::Input, b);
    }
    if (m.used.empty()) addFinding(review, 0, FindingLevel::Problem, FindingCause::Input, "This piece has no markings in use (accept a proposal or draw one).");
    else {
        if (choice.start < 0 && !choice.incompleteStart) addFinding(review, 0, FindingLevel::Problem, FindingCause::Marking, "Add a start: a curve (cast-on) or a point (pole start).");
        if (choice.finish < 0 && !choice.incompleteFinish)
            addFinding(review, 0, FindingLevel::Problem, FindingCause::Marking, "Add a finish: a curve (bind-off loop) or, after a start point, a point (the bind-off level through it).");
        if (choice.rowStart < 0 && !choice.incompleteRowStart) addFinding(review, 0, FindingLevel::Problem, FindingCause::Marking, "Add a row-start path from the start to the finish.");
    }
    const bool poleRoute = choice.start >= 0 && m.used[choice.start].marking->shape == MarkingShape::Point;
    const bool finishPoint = choice.finish >= 0 && m.used[choice.finish].marking->shape == MarkingShape::Point;
    if (choice.start >= 0 && choice.finish >= 0 && !poleRoute && finishPoint)
        addFinding(review, m.used[choice.finish].marking->id, FindingLevel::Problem, FindingCause::AlgorithmLimitation,
                   "A finish point (a bind-off level) after a start curve is not supported yet; draw a finish curve.");
}

uint64_t markedFieldKey(const ResolvedMarkings& m, uint64_t fileDigest, const std::string& preparation, int piece, double mmPerUnit, double h) {
    uint64_t hash = fnv1a64Bytes(&fileDigest, sizeof fileDigest);
    hash = fnv1a64(preparation, hash);
    const double numbers[2] = {mmPerUnit, h};
    hash = fnv1a64Bytes(numbers, sizeof numbers, hash);
    const int head[2] = {piece, m.choice.finishIsBoundary ? 1 : 0};
    hash = fnv1a64Bytes(head, sizeof head, hash);
    for (const int index : {m.choice.start, m.choice.finish}) {
        if (index < 0) {
            const int none = -1;
            hash = fnv1a64Bytes(&none, sizeof none, hash);
            continue;
        }
        const SurfaceMarking& k = *m.used[index].marking;
        const int ints[3] = {k.id, static_cast<int>(k.meaning), static_cast<int>(k.shape)};
        hash = fnv1a64Bytes(ints, sizeof ints, hash);
        for (const SurfacePoint& p : k.points) {
            hash = fnv1a64Bytes(&p.tri, sizeof p.tri, hash);
            hash = fnv1a64Bytes(p.bary, sizeof p.bary, hash);
        }
    }
    return hash;
}

std::vector<RowCrossing> rowCrossings(const MarkedField& f, const MarkingTrace& path) {
    std::vector<RowCrossing> rows;
    auto split = [&](int k) { return k < static_cast<int>(f.contoursPerRing.size()) && f.contoursPerRing[k] != 1; };
    // Both routes count as their rows builder does: the path's field values against each level (levelCrossings).
    const bool band = f.route == "band";
    const std::vector<MarkedPathPoint> pts = pathOnField(path, f.mesh, *f.query, f.t);
    for (int k = 1; k <= f.rows; ++k) {
        if ((band || f.loopFinish) && k == f.rows) continue;   // a finish loop: its origin is the point nearest the path
        RowCrossing rc;
        rc.row = k;
        rc.split = split(k) || (band && (k >= static_cast<int>(f.contours.size()) || f.contours[k].points.empty()));
        if (!rc.split) rc.count = levelCrossings(pts, k * f.rowHeight, &rc.pathParameters);
        for (double u : rc.pathParameters) rc.at.push_back(pathPointAt(pts, u));
        rows.push_back(rc);
    }
    return rows;
}

void crossingFindings(const std::vector<RowCrossing>& rows, const std::string& route, int startId, int pathId, double toFile, MarkingReview& review) {
    std::vector<int> split, missed, twice;
    for (const RowCrossing& rc : rows) {
        if (rc.split) split.push_back(rc.row);
        else if (rc.count == 0) missed.push_back(rc.row);
        else if (rc.count > 1) twice.push_back(rc.row);
    }
    auto place = [&](MarkingFinding& f, const RowCrossing& rc) {
        for (const Vec3& x : rc.at) f.at.push_back(x * toFile);
        for (double u : rc.pathParameters) {
            const int lo = static_cast<int>(std::floor(u)), hi = static_cast<int>(std::ceil(u));
            f.firstTracePoint = f.firstTracePoint < 0 ? lo : std::min(f.firstTracePoint, lo);
            f.lastTracePoint = std::max(f.lastTracePoint, std::max(hi, lo + 1));
        }
    };
    if (!split.empty()) {
        MarkingFinding f;
        f.markingId = startId;
        f.level = FindingLevel::Problem;
        f.cause = FindingCause::Marking;
        f.row = split.front();
        f.message = route == "band" ? "Rows " + ringList(split) + " split into several pieces between these curves. Move the start or finish curve (splitting the piece is not available yet)."
                                    : "Rows " + ringList(split) + " split into several pieces around this start. Move the start or the finish.";
        review.findings.push_back(f);
    }
    if (!twice.empty() && twice.size() <= 3) {
        for (const RowCrossing& rc : rows) {
            if (rc.split || rc.count <= 1) continue;
            MarkingFinding f;
            f.markingId = pathId;
            f.level = FindingLevel::Problem;
            f.cause = FindingCause::Marking;
            f.row = rc.row;
            char b[200];
            std::snprintf(b, sizeof b, "This path crosses row %d %s; straighten it so that it crosses every row once.", rc.row, rc.count == 2 ? "twice" : (std::to_string(rc.count) + " times").c_str());
            f.message = b;
            place(f, rc);
            review.findings.push_back(f);
        }
    } else if (!twice.empty()) {
        MarkingFinding f;
        f.markingId = pathId;
        f.level = FindingLevel::Problem;
        f.cause = FindingCause::Marking;
        f.row = twice.front();
        f.message = "The row-start path crosses rows " + ringList(twice) + " more than once. Straighten it so that it runs from start to finish without turning back.";
        for (const RowCrossing& rc : rows)
            if (!rc.split && rc.count > 1) place(f, rc);
        review.findings.push_back(f);
    }
    if (!missed.empty()) {
        MarkingFinding f;
        f.markingId = pathId;
        f.level = FindingLevel::Problem;
        f.cause = FindingCause::Marking;
        f.row = missed.front();
        f.message = "The row-start path misses rows " + ringList(missed) + ". Extend it from the start to the finish.";
        review.findings.push_back(f);
    }
    // where it turns back: along the path the rows it crosses must keep rising
    std::vector<std::pair<double, int>> order;
    for (const RowCrossing& rc : rows)
        for (size_t i = 0; i < rc.pathParameters.size(); ++i) order.push_back({rc.pathParameters[i], rc.row});
    std::sort(order.begin(), order.end());
    int turns = 0;
    for (size_t i = 0; i + 1 < order.size() && turns < 3; ++i) {
        if (order[i + 1].second > order[i].second) continue;   // the same row again, or a lower one: it turned back
        if (i > 0 && order[i].second <= order[i - 1].second) continue;   // still the same turn (one finding per turn)
        ++turns;
        MarkingFinding f;
        f.markingId = pathId;
        f.level = FindingLevel::Warning;
        f.cause = FindingCause::Marking;
        f.row = order[i].second;
        char b[160];
        std::snprintf(b, sizeof b, "The path turns back after row %d (it comes to row %d again).", order[i].second, order[i + 1].second);
        f.message = b;
        f.firstTracePoint = static_cast<int>(std::floor(order[i].first));
        f.lastTracePoint = static_cast<int>(std::ceil(order[i + 1].first));
        for (const RowCrossing& rc : rows)
            for (size_t j = 0; j < rc.pathParameters.size(); ++j)
                if (rc.pathParameters[j] == order[i].first || rc.pathParameters[j] == order[i + 1].first) f.at.push_back(rc.at[j] * toFile);
        review.findings.push_back(f);
    }
    if (split.empty() && twice.empty() && missed.empty()) addFinding(review, pathId, FindingLevel::Ok, FindingCause::Marking, "This path crosses every row once.");
}

bool buildFromMarkings(const MarkedBuildRequest& r, MarkedBuildProduct& out, std::string* err) {
    out = MarkedBuildProduct{};
    if (!r.pieces || !r.markings || r.piece < 0 || r.piece >= static_cast<int>(r.pieces->pieces.size()) || r.mmPerUnit <= 0 || r.w <= 0 || r.h <= 0) {
        if (err) *err = "buildFromMarkings: no piece, no markings, or a gauge or unit that is not positive";
        return false;
    }
    const std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
    MarkingReview& review = out.review;
    MarkedEvaluation& ev = out.evaluation;
    review.depth = r.depth;
    const SurfacePiece& piece = r.pieces->pieces[r.piece];
    TriMesh pieceMm = piece.mesh;
    pieceMm.transform(r.mmPerUnit, Vec3{0, 0, 0});
    // the prepared reference: the piece refined so the field is sampled finely; Facets does not move the surface, LoopLimit
    // is the smooth-reference variant (prepared_reference.h)
    review.measures.secondsKnown = Known::Current;
    std::shared_ptr<const PreparedPiece> prep = r.prepared;
    if (r.refine >= 0 || !prep || prep->key != preparedPieceKey(r.markings->digest, r.piece, r.mmPerUnit, r.h, r.options.preparation) || prep->piece != r.piece) {
        std::shared_ptr<PreparedPiece> own = std::make_shared<PreparedPiece>();
        if (r.refine >= 0) {
            own->piece = r.piece;
            own->refined = prepareReference(pieceMm, r.refine, r.options.preparation);
            own->nb = triangleNeighbours(own->refined.mesh);
        } else {
            preparePiece(*r.pieces, r.piece, r.markings->digest, r.mmPerUnit, r.h, r.options.preparation, *own);
        }
        prep = own;
    }
    if (r.layerOffsetMm != 0) prep = offsetReference(*prep, r.layerOffsetMm);
    cancelCheckpoint();
    const RefinedMesh& refined = prep->refined;
    const int levels = refined.levels;
    const TriMesh& mm = refined.mesh;
    RefinedMesh other;   // the controlled comparison's relaxation reference
    if (r.options.relaxOnOther) other = prepareReference(pieceMm, levels, r.options.relaxPreparation);
    ev.refineLevels = levels;
    ev.reference = mm;
    ev.preparation = toString(r.options.preparation);
    ev.relaxPreparation = toString(r.options.relaxOnOther ? r.options.relaxPreparation : r.options.preparation);
    ev.original = r.options.preparation == ReferencePreparation::Facets ? mm : refineMidpoint(pieceMm, levels).mesh;
    const TriMesh& judged = r.options.relaxOnOther ? other.mesh : mm;
    ev.referenceDigest = meshDigest(judged);
    ev.areaOriginal = ev.original.area();
    ev.areaReference = judged.area();
    const TriangleNeighbours& nb = prep->nb;
    const double toFile = 1.0 / r.mmPerUnit;
    char b[400];
    std::snprintf(b, sizeof b, "prepared reference: piece %d, %zu triangles, %s x%d -> %zu triangles (mm), digest %016llx", r.piece, pieceMm.tri.size(),
                  r.options.preparation == ReferencePreparation::Facets ? "midpoint-refined" : "Loop-subdivided to the limit", levels, mm.tri.size(),
                  static_cast<unsigned long long>(meshDigest(mm)));
    out.log.push_back(b);
    if (ev.relaxPreparation != "facets") {
        // the change from the imported surface, sampled at the reference's vertices both ways (the same split: as many each way)
        const MeshQuery qo(ev.original), qj(judged);
        std::vector<double> there, back;
        for (const Vec3& x : judged.pos) there.push_back(norm(qo.closest(x).p - x));
        for (const Vec3& x : ev.original.pos) back.push_back(norm(qj.closest(x).p - x));
        const Quantiles a = quantiles(there), c = quantiles(back);
        ev.changeMax = a.max;
        ev.changeP95 = a.p95;
        ev.changeBackMax = c.max;
        ev.changeBackP95 = c.p95;
        std::snprintf(b, sizeof b, "The %s reference differs from the imported surface by up to %.2f mm (95 %% within %.2f mm; back %.2f / %.2f mm); area %.0f mm2 against %.0f mm2. "
                      "Cells are judged against this reference; their distance from the imported surface is reported separately.",
                      r.options.relaxOnOther ? "relaxation's" : "prepared", std::max(a.max, c.max), std::max(a.p95, c.p95), c.max, c.p95, ev.areaReference, ev.areaOriginal);
        addFinding(review, 0, FindingLevel::Ok, FindingCause::Input, std::string("Smooth reference (Loop limit): ") + b);
        out.log.push_back(b);
    }

    // 1. this piece's markings, resolved, mapped to the reference and traced
    ResolvedMarkings rm;
    resolvePieceMarkings(r, refined, nb, review, rm);
    std::vector<const SurfaceMarking*> used;
    std::vector<MarkingTrace> traceMm;
    std::vector<std::vector<SurfacePoint>> localPoints;   // per used marking: its points on the reference
    for (const ResolvedMarking& u : rm.used) {
        used.push_back(u.marking);
        traceMm.push_back(u.trace);
        localPoints.push_back(u.local);
    }
    const MeaningChoice& choice = rm.choice;
    out.markingIds = rm.markingIds;
    out.previewIds = rm.previewIds;
    out.traces = rm.traces;
    if (r.options.preparation != ReferencePreparation::Facets) {
        for (const ResolvedMarking& u : rm.used) {
            if (!u.resolved) continue;
            // the same control points (the same refined triangle and weights) on the imported surface, and the same rule traced there
            MarkedEvaluation::MarkingShift ms;
            ms.id = u.marking->id;
            for (const SurfacePoint& lp : u.local) ms.controlMax = std::max(ms.controlMax, norm(surfacePosition(mm, lp) - surfacePosition(ev.original, lp)));
            const MarkingTrace onOriginal = traceMarking(ev.original, nb, u.local, u.marking->shape == MarkingShape::ClosedCurve);
            std::vector<Vec3> a, c;
            for (const TracePoint& x : u.trace.points) a.push_back(x.pos);
            for (const TracePoint& x : onOriginal.points) c.push_back(x.pos);
            auto oneWay = [](const std::vector<Vec3>& from, const std::vector<Vec3>& to) {
                double worst = 0;
                for (const Vec3& x : from) worst = std::max(worst, to.size() > 1 ? distanceToPolyline(x, to) : to.empty() ? 0.0 : norm(x - to[0]));
                return worst;
            };
            ms.traceHausdorff = std::max(oneWay(a, c), oneWay(c, a));
            ev.markingShift.push_back(ms);
        }
    }
    if (!ev.markingShift.empty()) {
        double control = 0, trace = 0;
        for (const MarkedEvaluation::MarkingShift& ms : ev.markingShift) {
            control = std::max(control, ms.controlMax);
            trace = std::max(trace, ms.traceHausdorff);
        }
        std::snprintf(b, sizeof b, "The markings keep their file triangles and weights; on the smooth reference their control points lie up to %.2f mm from the imported "
                      "surface's, and their traced paths up to %.2f mm from the paths traced there.", control, trace);
        addFinding(review, 0, FindingLevel::Ok, FindingCause::Input, b);
    }
    reviewMeaningChoice(rm, review);
    cancelCheckpoint();
    const bool poleRoute = choice.start >= 0 && used[choice.start]->shape == MarkingShape::Point;
    const bool finishPoint = choice.finish >= 0 && used[choice.finish]->shape == MarkingShape::Point;
    if (hasProblem(review)) {
        review.measures.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        review.recommendation = "Resolve the problems above, then build again.";
        return true;
    }
    const int startId = used[choice.start]->id, finishId = used[choice.finish]->id, pathId = used[choice.rowStart]->id;
    HiddenTest hidden;
    hiddenTest(*r.pieces, r.piece, r.mmPerUnit, hidden);
    const MarkingTrace& path = traceMm[choice.rowStart];
    std::vector<Vec3> pathLine;
    for (const TracePoint& p : path.points) pathLine.push_back(p.pos);
    Tail tail;
    tail.r = &r;
    tail.out = &out;
    tail.finishId = finishId;
    tail.pathLine = pathLine;
    tail.relaxTo = r.options.relaxOnOther ? &other.mesh : nullptr;
    // per reference triangle: knitted (true) or omitted, and omitted's kind; then per piece triangle by the larger part
    std::vector<MarkedRegion> refRegion(mm.tri.size(), MarkedRegion::Kept);
    std::vector<double> refKeptArea(mm.tri.size(), 0);
    auto omit = [&](const TriMesh& m, int t, MarkedRegion kind, int refTri) {
        const double a = triangleArea(m, t);
        review.measures.omittedAreaMm2 += a;
        const Vec3 c = (m.pos[m.tri[t][0]] + m.pos[m.tri[t][1]] + m.pos[m.tri[t][2]]) / 3.0;
        (hidden.inside(c) ? review.measures.omittedHiddenMm2 : review.measures.omittedVisibleMm2) += a;
        refRegion[refTri] = kind;
    };
    auto pieceRegions = [&](const TriMesh& m, const std::vector<int>& refOf, const std::vector<char>& keptTri, const std::vector<MarkedRegion>& kindOf) {
        out.region.assign(pieceMm.tri.size(), MarkedRegion::Kept);
        std::vector<double> largest(pieceMm.tri.size(), -1);
        for (int t = 0; t < static_cast<int>(m.tri.size()); ++t) {
            const int parent = refined.parent[refOf[t]];
            const double a = triangleArea(m, t);
            if (a <= largest[parent]) continue;
            largest[parent] = a;
            out.region[parent] = keptTri[t] ? MarkedRegion::Kept : kindOf[t];
        }
        out.regionHidden.assign(pieceMm.tri.size(), 0);
        for (size_t pt = 0; pt < pieceMm.tri.size(); ++pt)
            if (out.region[pt] != MarkedRegion::Kept) {
                const std::array<int, 3>& c = pieceMm.tri[pt];
                out.regionHidden[pt] = hidden.inside((pieceMm.pos[c[0]] + pieceMm.pos[c[1]] + pieceMm.pos[c[2]]) / 3.0) ? 1 : 0;
            }
        review.measures.omittedKnown = Known::Current;
    };
    std::vector<double> t;
    std::string fieldErr;
    int R = 0;
    TypedStitchMesh stored;

    if (!poleRoute) {
        // ------------------------------------------------------------------ band: start curve, finish curve
        ev.route = "band";
        // a start or finish lying on a border loop of the piece is that opening itself: no trim along it, nothing beyond it
        const std::vector<int>& startLoop = rm.used[choice.start].borderLoop;
        const std::vector<int>& finishLoop = rm.used[choice.finish].borderLoop;
        std::vector<MarkingTrace> cut;
        int startCurve = -1, finishCurve = -1;
        if (startLoop.empty()) {
            startCurve = static_cast<int>(cut.size());
            cut.push_back(traceMm[choice.start]);
        }
        if (finishLoop.empty()) {
            finishCurve = static_cast<int>(cut.size());
            cut.push_back(traceMm[choice.finish]);
        }
        TrimmedSurface trim;
        std::string trimErr;
        if (!splitAlongCurves(mm, cut, trim, &trimErr)) {
            addFinding(review, 0, FindingLevel::Problem, FindingCause::Marking, "The start and finish curves cannot trim the piece: " + trimErr + ".");
            review.recommendation = "Move the start or finish curve.";
            return true;
        }
        cancelCheckpoint();
        int band = -1;
        if (cut.size() == 2) {
            for (int g = 0; g < trim.regions; ++g)
                if (trim.touches[g][0] && trim.touches[g][1]) band = band < 0 ? g : -2;
        } else {
            // the region holding the opening's first edge, which must also border the drawn curve (if any)
            const std::vector<int>& loop = startLoop.empty() ? finishLoop : startLoop;
            const int g = regionOfEdge(trim, loop[0], loop[1]);
            band = g >= 0 && (cut.empty() || trim.touches[g][0]) ? g : -1;
        }
        if (band < 0) {
            addFinding(review, 0, FindingLevel::Problem, FindingCause::Marking,
                       "The start and finish curves must divide the piece into one band between them (each curve going once around the piece).");
            review.recommendation = "Draw each curve once around the piece, apart from each other.";
            return true;
        }
        std::vector<char> keptTri(trim.mesh.tri.size(), 0);
        std::vector<MarkedRegion> kindOf(trim.mesh.tri.size(), MarkedRegion::Kept);
        for (int tt = 0; tt < static_cast<int>(trim.mesh.tri.size()); ++tt) {
            const int g = trim.region[tt];
            keptTri[tt] = g == band;
            kindOf[tt] = g == band                                     ? MarkedRegion::Kept
                         : startCurve >= 0 && trim.touches[g][startCurve] ? MarkedRegion::OmittedStart
                         : choice.finishIsBoundary                        ? MarkedRegion::OutsideBoundary
                                                                          : MarkedRegion::OmittedFinish;
            if (g != band) omit(trim.mesh, tt, kindOf[tt], trim.parentTri[tt]);
        }
        pieceRegions(trim.mesh, trim.parentTri, keptTri, kindOf);
        const RegionMesh kept = regionMesh(trim, band);
        std::vector<int> castOn, bindOff;
        if (startLoop.empty()) castOn = kept.loops[startCurve];
        else
            for (int v : startLoop) castOn.push_back(kept.local[v]);
        if (finishLoop.empty()) bindOff = kept.loops[finishCurve];
        else
            for (int v : finishLoop) bindOff.push_back(kept.local[v]);
        const int chi = eulerCharacteristic(kept.mesh);
        if (chi != 0) {
            std::snprintf(b, sizeof b, "The surface between the start and finish curves is not one tube (Euler characteristic %d, a tube has 0).", chi);
            addFinding(review, 0, FindingLevel::Problem, FindingCause::Input, b);
            review.recommendation = "Method A needs one tube between the curves.";
            return true;
        }
        TimeFieldReport tfr;
        if (!geodesicDistance(kept.mesh, castOn, t, tfr, &fieldErr)) {
            addFinding(review, 0, FindingLevel::Problem, FindingCause::Input, "No time field: " + fieldErr);
            return true;
        }
        double T = 1e300;
        for (int v : bindOff) T = std::min(T, t[v]);
        R = std::max(1, static_cast<int>(std::lround(T / r.h)));
        std::vector<Vec3> loop0, loopR;
        for (int v : castOn) loop0.push_back(kept.mesh.pos[v] * toFile);
        for (int v : bindOff) loopR.push_back(kept.mesh.pos[v] * toFile);
        out.rings.push_back(loop0);
        for (int k = 1; k < R; ++k) {
            const RingLevel lv = extractLevel(kept.mesh, t, k * T / R);
            for (const Contour& c : lv.contours) {
                std::vector<Vec3> line;
                for (const RingPoint& p : c.points) line.push_back(p.p * toFile);
                if (c.closed && !line.empty()) line.push_back(line.front());
                out.rings.push_back(line);
            }
        }
        out.rings.push_back(loopR);
        const SurfacePathGuide guide = guideFromTrace(path, &kept.parentTri);
        const MeshQuery keptQuery(kept.mesh);
        const std::vector<MarkedPathPoint> pathField = pathOnField(path, kept.mesh, keptQuery, t);
        GuidedRowsRequest req;
        req.w = r.w;
        req.h = r.h;
        req.path = &guide;           // the boundary rings' origins (nearest point)
        req.pathField = &pathField;  // the interior rings' crossings and origins
        req.adjacentShaping = r.options.adjacentShaping;
        req.boundaryShaping = r.options.boundaryShaping;
        req.rowPhase = r.options.rowPhase;
        req.bindOffCount = r.options.bindOffCount;
        GuidedRowsProduct gp;
        std::string rowsErr;
        if (!buildGuidedRows(kept.mesh, t, castOn, bindOff, req, gp, &rowsErr)) {
            addFinding(review, 0, FindingLevel::Problem, FindingCause::Input, "Method A refused the band: " + rowsErr);
            return true;
        }
        for (const std::string& l : gp.log) out.log.push_back(l);
        review.measures.rows = gp.rows;
        review.measures.contourEvents = gp.contourEvents;
        review.measures.rowsKnown = review.measures.contourEventsKnown = review.measures.jumpKnown = Known::Current;
        {
            // the row field, kept for quick reviews of this layout's start and finish
            std::shared_ptr<MarkedField> f = std::make_shared<MarkedField>();
            f->key = markedFieldKey(rm, r.markings->digest, toString(r.options.preparation), r.piece, r.mmPerUnit, r.h);
            f->reference = prep;
            f->route = "band";
            f->mesh = kept.mesh;
            f->parentTri = kept.parentTri;
            f->t = t;
            f->query = std::make_unique<MeshQuery>(f->mesh);   // built on the field's own mesh, as the rows builder's lookup
            f->rows = gp.rows;
            f->rowHeight = gp.hPrime;
            f->contours = gp.contours;
            for (const RingSample& rs : gp.rings) f->contoursPerRing.push_back(rs.contours);
            out.field = f;
        }
        crossingFindings(rowCrossings(*out.field, path), "band", startId, pathId, toFile, review);
        knittedOverFindings(gp.knittedOver, startId, toFile, review);
        cancelCheckpoint();
        int maxJump = 0;
        for (size_t k = 1; k < gp.rings.size(); ++k) maxJump = std::max(maxJump, std::abs(gp.rings[k].n - gp.rings[k - 1].n));
        review.measures.maxStitchCountJump = maxJump;
        if (gp.status == Status::Fail || gp.mesh.mesh.faces() == 0) {
            review.measures.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            review.recommendation = "No stitch mesh: correct the marked rows above.";
            return true;
        }
        stored = gp.mesh;
        ev.knitRegion = kept.mesh;
        ev.regionLoops = {castOn, bindOff};
        tail.s = stored;
        ev.border.assign(stored.mesh.pos.size(), 0);
        for (const std::vector<int>& loop : polygonBorderLoops(stored.mesh))
            for (int v : loop) ev.border[v] = 1;
        tail.surfaceLoops = {castOn, bindOff};
        ev.relaxLoops = tail.surfaceLoops;
        if (r.options.anchors == RelaxAnchors::RowOriginsOnPath) {
            const MeshQuery q(kept.mesh);
            tail.guides.lines = {guideFromPath(path, kept.mesh, q)};
            tail.guides.lineOf.assign(stored.mesh.pos.size(), -1);
            tail.guides.pinned.assign(stored.mesh.pos.size(), 0);
            for (size_t k = 1; k + 1 < stored.rings.size(); ++k) tail.guides.lineOf[stored.rings[k][0]] = 0;
        } else if (r.options.anchors == RelaxAnchors::YarnEnds) {
            // the yarn's two ends stay where the layout put them (RelaxAnchors::YarnEnds)
            tail.guides.lineOf.assign(stored.mesh.pos.size(), -1);
            tail.guides.pinned.assign(stored.mesh.pos.size(), 0);
            tail.guides.pinned[stored.rings.front()[0]] = 1;
            tail.guides.pinned[stored.rings.back()[0]] = 1;
        }
        DeclaredBoundary d0, d1;
        d0.role = EdgeBoundary::CastOn;
        d1.role = EdgeBoundary::BindOff;
        for (int v : castOn) d0.loop.push_back(kept.mesh.pos[v]);
        for (int v : bindOff) d1.loop.push_back(kept.mesh.pos[v]);
        tail.declared = {d0, d1};
        ev.guides = tail.guides;
        // the walk runs on the band (its border is the two curves)
        TriMesh bandSurface = kept.mesh;
        tail.surface = &bandSurface;
        finishBuild(tail);
        ev.guides = tail.guides;
        review.measures.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::snprintf(b, sizeof b, "Built %d rows, %d cells (%d exceptional, %d without a construction).", gp.rows, out.stitchMesh.mesh.faces(),
                      review.measures.exceptionalCells, review.measures.unsupportedCells);
        review.recommendation = b;
        return true;
    }

    // ---------------------------------------------------------------------- pole: start point, finish curve or point
    ev.route = "pole";
    const std::vector<int>& finishOpening = rm.used[choice.finish].borderLoop;   // the finish is the piece's opening itself
    if (!finishOpening.empty() && eulerCharacteristic(mm) != 1) {   // a connected piece with chi 1 is a disc: one opening
        addFinding(review, finishId, FindingLevel::Problem, FindingCause::AlgorithmLimitation,
                   "A start point with an opening as the finish needs a piece with exactly one opening (a cup); draw a finish curve instead.");
        review.recommendation = "Draw a finish curve.";
        return true;
    }
    const SurfacePoint startPt = localPoints[choice.start][0];
    const std::string ring1 = markingSetting(*used[choice.start], "ring1");
    MarkedRowsRequest req;
    req.w = r.w;
    req.h = r.h;
    req.adjacentShaping = r.options.adjacentShaping;
    req.boundaryShaping = r.options.boundaryShaping;
    req.rowPhase = r.options.rowPhase;
    req.bindOffCount = r.options.bindOffCount;
    req.ring1 = ring1.empty() || ring1 == "auto" ? 0 : std::atoi(ring1.c_str());
    const TriMesh* rowsMesh = &mm;
    TriMesh keptMesh;
    std::vector<int> keptParent;
    if (finishPoint) {
        if (!fieldFromPoint(mm, startPt, t, &fieldErr)) {
            addFinding(review, 0, FindingLevel::Problem, FindingCause::Input, "No time field: " + fieldErr);
            return true;
        }
        ev.finishLevel = valueAt(mm, t, localPoints[choice.finish][0]);
        req.finishLevel = ev.finishLevel;
        keptMesh = clipBelowLevel(mm, t, ev.finishLevel, keptParent);
        // omitted: the reference above the level
        std::vector<int> aboveParent;
        std::vector<double> negated(t.size());
        for (size_t v = 0; v < t.size(); ++v) negated[v] = -t[v];
        const TriMesh above = clipBelowLevel(mm, negated, -ev.finishLevel, aboveParent);
        for (int tt = 0; tt < static_cast<int>(above.tri.size()); ++tt) omit(above, tt, MarkedRegion::OmittedFinish, aboveParent[tt]);
        cancelCheckpoint();
        TriMesh both = keptMesh;
        std::vector<int> refOf = keptParent;
        std::vector<char> keptTri(keptMesh.tri.size(), 1);
        std::vector<MarkedRegion> kindOf(keptMesh.tri.size(), MarkedRegion::Kept);
        const int base = static_cast<int>(both.pos.size());
        for (const Vec3& p : above.pos) both.pos.push_back(p);
        for (int tt = 0; tt < static_cast<int>(above.tri.size()); ++tt) {
            both.tri.push_back({above.tri[tt][0] + base, above.tri[tt][1] + base, above.tri[tt][2] + base});
            refOf.push_back(aboveParent[tt]);
            keptTri.push_back(0);
            kindOf.push_back(MarkedRegion::OmittedFinish);
        }
        pieceRegions(both, refOf, keptTri, kindOf);
        const MeshQuery q(mm);
        req.path = pathOnField(path, mm, q, t);
    } else if (!finishOpening.empty()) {
        // the opening as the finish: the whole piece is the disc round the start point, the opening its bind-off ring
        if (!fieldFromPoint(mm, startPt, t, &fieldErr)) {
            addFinding(review, 0, FindingLevel::Problem, FindingCause::Input, "No time field: " + fieldErr);
            return true;
        }
        if (r.options.openingRows == OpeningRows::OpeningIsALevel) {
            // t' = T t / (t + d), d the distance to the opening and T the mean of t along it: 0 at the start, T all along
            // the opening, so every row is a closed curve between the two and the last strip is as tall as the others
            std::vector<double> toOpening;
            TimeFieldReport rep;
            if (!geodesicDistance(mm, finishOpening, toOpening, rep, &fieldErr)) {
                addFinding(review, 0, FindingLevel::Problem, FindingCause::Input, "No time field: " + fieldErr);
                return true;
            }
            double sum = 0, nearest = 1e300, farthest = 0;
            for (int v : finishOpening) {
                sum += t[v];
                nearest = std::min(nearest, t[v]);
                farthest = std::max(farthest, t[v]);
            }
            const double level = sum / finishOpening.size();
            for (size_t v = 0; v < t.size(); ++v) t[v] = t[v] + toOpening[v] > 0 ? level * t[v] / (t[v] + toOpening[v]) : 0.0;
            char note[200];
            std::snprintf(note, sizeof note, "opening as a level: its distance from the start is %.1f to %.1f mm (mean %.1f); rows between the start and the opening", nearest,
                          farthest, level);
            out.log.push_back(note);
        }
        keptMesh = mm;   // the knit region is the whole piece (same vertex and triangle numbering)
        keptParent.resize(mm.tri.size());
        for (int tt = 0; tt < static_cast<int>(mm.tri.size()); ++tt) keptParent[tt] = tt;
        pieceRegions(mm, keptParent, std::vector<char>(mm.tri.size(), 1), std::vector<MarkedRegion>(mm.tri.size(), MarkedRegion::Kept));
        const MeshQuery q(keptMesh);
        req.bindOffLoop = finishOpening;
        req.path = pathOnField(path, keptMesh, q, t);
        rowsMesh = &keptMesh;
    } else {
        // a finish curve: the region it bounds around the start point
        TrimmedSurface trim;
        std::string trimErr;
        if (!splitAlongCurves(mm, {traceMm[choice.finish]}, trim, &trimErr)) {
            addFinding(review, finishId, FindingLevel::Problem, FindingCause::Marking, "The finish curve cannot trim the piece: " + trimErr + ".");
            return true;
        }
        const Vec3 sp = surfacePosition(mm, startPt);
        const MeshQuery tq(trim.mesh);
        const SurfacePoint inTrim = locate(tq, trim.mesh, sp);
        const int region = inTrim.tri >= 0 ? trim.region[inTrim.tri] : -1;
        std::vector<char> keptTri(trim.mesh.tri.size(), 0);
        std::vector<MarkedRegion> kindOf(trim.mesh.tri.size(), MarkedRegion::Kept);
        for (int tt = 0; tt < static_cast<int>(trim.mesh.tri.size()); ++tt) {
            keptTri[tt] = trim.region[tt] == region;
            kindOf[tt] = keptTri[tt] ? MarkedRegion::Kept : choice.finishIsBoundary ? MarkedRegion::OutsideBoundary : MarkedRegion::OmittedFinish;
            if (!keptTri[tt]) omit(trim.mesh, tt, kindOf[tt], trim.parentTri[tt]);
        }
        pieceRegions(trim.mesh, trim.parentTri, keptTri, kindOf);
        const RegionMesh kept = regionMesh(trim, region);
        keptMesh = kept.mesh;
        keptParent = kept.parentTri;
        if (eulerCharacteristic(keptMesh) != 1 || kept.loops.empty() || kept.loops[0].empty()) {
            std::snprintf(b, sizeof b, "The finish curve must enclose the start point in one disc (Euler characteristic %d, a disc has 1).", eulerCharacteristic(keptMesh));
            addFinding(review, finishId, FindingLevel::Problem, FindingCause::Marking, b);
            return true;
        }
        const MeshQuery kq(keptMesh);
        const SurfacePoint sk = locate(kq, keptMesh, sp);
        if (!fieldFromPoint(keptMesh, sk, t, &fieldErr)) {
            addFinding(review, 0, FindingLevel::Problem, FindingCause::Input, "No time field: " + fieldErr);
            return true;
        }
        req.bindOffLoop = kept.loops[0];
        req.path = pathOnField(path, keptMesh, kq, t);
        rowsMesh = &keptMesh;
    }
    req.pole = surfacePosition(mm, startPt);
    {
        const MeshQuery rq(*rowsMesh);
        const SurfacePoint ps = locate(rq, *rowsMesh, req.pole);
        req.poleBinding.tri = ps.tri;
        for (int k = 0; k < 3; ++k) req.poleBinding.bary[k] = ps.bary[k];
    }
    MarkedRowsProduct mp;
    std::string rowsErr;
    if (!buildMarkedRows(*rowsMesh, t, req, mp, &rowsErr)) {
        addFinding(review, 0, FindingLevel::Problem, FindingCause::Input, "Method A refused the piece: " + rowsErr);
        return true;
    }
    for (const std::string& l : mp.rows.log) out.log.push_back(l);
    R = mp.rows.rows;
    review.measures.rows = R;
    review.measures.contourEvents = mp.rows.contourEvents + mp.guideProblems;
    review.measures.rowsKnown = review.measures.contourEventsKnown = review.measures.jumpKnown = Known::Current;
    {
        // the row field, kept for quick reviews of this layout's start and finish
        std::shared_ptr<MarkedField> f = std::make_shared<MarkedField>();
        f->key = markedFieldKey(rm, r.markings->digest, toString(r.options.preparation), r.piece, r.mmPerUnit, r.h);
        f->reference = prep;
        f->route = "pole";
        f->mesh = *rowsMesh;
        if (rowsMesh != &mm) f->parentTri = keptParent;
        f->t = t;
        f->query = std::make_unique<MeshQuery>(f->mesh);
        f->rows = R;
        f->rowHeight = mp.rows.hPrime;
        f->loopFinish = !finishPoint;
        for (const RingSample& rs : mp.rows.rings) f->contoursPerRing.push_back(rs.contours);
        out.field = f;
    }
    cancelCheckpoint();
    out.rings.push_back({req.pole * toFile});
    for (int k = 1; k <= R; ++k) {
        if (finishPoint || k < R) {
            const RingLevel lv = extractLevel(*rowsMesh, t, k * mp.rows.hPrime);
            for (const Contour& c : lv.contours) {
                std::vector<Vec3> line;
                for (const RingPoint& p : c.points) line.push_back(p.p * toFile);
                if (c.closed && !line.empty()) line.push_back(line.front());
                out.rings.push_back(line);
            }
        } else {
            std::vector<Vec3> line;
            for (int v : req.bindOffLoop) line.push_back(rowsMesh->pos[v] * toFile);
            if (!line.empty()) line.push_back(line.front());
            out.rings.push_back(line);
        }
    }
    crossingFindings(rowCrossings(*out.field, path), "pole", startId, pathId, toFile, review);
    knittedOverFindings(mp.rows.knittedOver, startId, toFile, review);
    int maxJump = 0;
    for (size_t k = 2; k < mp.rows.rings.size(); ++k) maxJump = std::max(maxJump, std::abs(mp.rows.rings[k].n - mp.rows.rings[k - 1].n));
    review.measures.maxStitchCountJump = maxJump;
    if (req.ring1 == 0 && mp.ring1Auto > 0) {
        std::snprintf(b, sizeof b, "Pole start: first ring of %d stitches (automatic: round(L1 / w), at least 3, raised until the first strip's increases fit); set ring1 on the start point to choose.", mp.ring1Auto);
        addFinding(review, startId, FindingLevel::Ok, FindingCause::Marking, b);
    }
    if (mp.rows.status == Status::Fail || mp.rows.mesh.mesh.faces() == 0) {
        review.measures.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        review.recommendation = "No stitch mesh: correct the marked rows above.";
        return true;
    }
    stored = mp.rows.mesh;
    ev.knitRegion = keptMesh;
    ev.regionLoops = borderLoops(keptMesh);
    tail.s = stored;
    tail.pole = stored.rings[0][0];
    ev.border.assign(stored.mesh.pos.size(), 0);
    tail.guides.lineOf.assign(stored.mesh.pos.size(), -1);
    tail.guides.pinned.assign(stored.mesh.pos.size(), 0);
    if (r.options.anchors == RelaxAnchors::RowOriginsOnPath) tail.guides.pinned[tail.pole] = 1;
    DeclaredBoundary d1;
    d1.role = EdgeBoundary::BindOff;
    if (finishPoint) {
        // the walk on the whole reference; the bind-off ring slides on its level's contour
        const RingLevel fl = extractLevel(mm, t, ev.finishLevel);
        GuideLine finishLine;
        finishLine.closed = true;
        for (const RingPoint& p : fl.contours.front().points) {
            finishLine.points.push_back(p.p);
            finishLine.segmentTri.push_back(p.tri);
            d1.loop.push_back(p.p);
        }
        tail.guides.lines.push_back(finishLine);
        for (int v : stored.rings[R]) tail.guides.lineOf[v] = 0;
        tail.surface = &ev.reference;
        if (r.options.relaxOnOther && r.options.otherGuides != MarkedBuildOptions::OtherGuides::Carried) {
            // the same rules on the other preparation: the field from the start point, its level through the finish point, the
            // row-start path traced between the same control points
            std::vector<double> tOther;
            std::string e2;
            if (fieldFromPoint(other.mesh, startPt, tOther, &e2)) {
                // the level of the other field that the carried bind-off ring sits on, on average (the field's value at the
                // finish point itself differs: geodesic lengths change with the surface)
                const MeshQuery qm(mm);
                double levelOther = 0;
                for (int v : stored.rings[R]) levelOther += valueAt(other.mesh, tOther, locate(qm, mm, stored.mesh.pos[v]));
                levelOther /= std::max<size_t>(1, stored.rings[R].size());
                const double atFinishPoint = valueAt(other.mesh, tOther, localPoints[choice.finish][0]);
                if (r.options.otherGuides == MarkedBuildOptions::OtherGuides::FinishPoint) levelOther = atFinishPoint;
                std::vector<int> unusedParent;
                tail.relaxToRegion = clipBelowLevel(other.mesh, tOther, levelOther, unusedParent);
                const RingLevel flo = extractLevel(other.mesh, tOther, levelOther);
                GuideLine finishOther;
                finishOther.closed = true;
                if (!flo.contours.empty())
                    for (const RingPoint& p : flo.contours.front().points) {
                        finishOther.points.push_back(p.p);
                        finishOther.segmentTri.push_back(p.tri);
                    }
                tail.relaxToLines.push_back(finishOther);
                if (r.options.anchors == RelaxAnchors::RowOriginsOnPath) {
                    const MarkingTrace pathOther = traceMarking(other.mesh, nb, localPoints[choice.rowStart], false);
                    const MeshQuery qo(other.mesh);
                    tail.relaxToLines.push_back(guideFromPath(pathOther, other.mesh, qo));
                    for (const TracePoint& p : pathOther.points) tail.relaxToPath.push_back(p.pos);
                }
                std::snprintf(b, sizeof b, "controlled comparison: constraint lines derived again on the %s reference: finish contour at level %.3f mm "
                              "(the field there is %.3f mm at the finish point; %.3f mm on the layout's reference); the knit region clipped there",
                              toString(r.options.relaxPreparation), levelOther, atFinishPoint, ev.finishLevel);
                out.log.push_back(b);
            }
        }
    } else {
        for (int v : req.bindOffLoop) d1.loop.push_back(keptMesh.pos[v]);
        for (int v : stored.rings[R]) ev.border[v] = 1;
        tail.surfaceLoops = {req.bindOffLoop};
        ev.relaxLoops = tail.surfaceLoops;
        tail.surface = &keptMesh;
    }
    if (r.options.anchors == RelaxAnchors::RowOriginsOnPath) {
        const MeshQuery q(*tail.surface);
        tail.guides.lines.push_back(guideFromPath(path, *tail.surface, q));
        const int pathLine = static_cast<int>(tail.guides.lines.size()) - 1;
        for (int k = 1; k < R; ++k) tail.guides.lineOf[stored.rings[k][0]] = pathLine;
    } else if (r.options.anchors == RelaxAnchors::YarnEnds) {
        // the yarn's end, the first vertex of the bind-off ring, stays where the layout put it; the pole is free
        tail.guides.pinned[stored.rings[R][0]] = 1;
    }
    tail.declared = {d1};
    finishBuild(tail);
    ev.guides = tail.guides;
    review.measures.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::snprintf(b, sizeof b, "Built %d rows from a pole start, %d cells (%d exceptional, %d without a construction); pole moved %.2f mm.", R, out.stitchMesh.mesh.faces(),
                  review.measures.exceptionalCells, review.measures.unsupportedCells, ev.poleMoved);
    review.recommendation = b;
    return true;
}

} // namespace knit
