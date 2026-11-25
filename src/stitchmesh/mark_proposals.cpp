// Proposals for one piece (mark_proposals.h). The closed-piece rules and helpers are the lab's (stitchmesh_lab
// propose-marks), moved here unchanged; the open-tube, isolated and branching
// rules are new and use the same helpers.
#include "stitchmesh/mark_proposals.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <memory>
#include <set>
#include <utility>

#include "model/cancel.h"
#include "model/mesh_query.h"
#include "model/surface_path.h"
#include "stitchmesh/marked_build.h"
#include "stitchmesh/row_rings.h"
#include "stitchmesh/time_field.h"

namespace knit {

namespace {

constexpr double kNearEnd = 15;        // rule 1: the neighbourhood of a diameter end whose hidden share decides the start [mm]
constexpr int kCurvePoints = 12;       // control points of a level curve or a border-loop marking
constexpr int kWideFinishStitches = 8; // a finish ring below this many stitches is called narrow (the leg's 4-6 failed)

std::string text(const char* fmt, ...) {
    char b[1024];
    va_list args;
    va_start(args, fmt);
    std::vsnprintf(b, sizeof b, fmt, args);
    va_end(args);
    return b;
}

void advance(const ProposalRequest& r) {
    if (r.step) r.step->fetch_add(1);
    cancelCheckpoint();
}

// The piece the rules work on: in mm, refined to h / 2 on the request's preparation, with the closed pieces of the same
// component that can hide it (the lab's PieceField).
struct PieceField {
    const SurfacePiece* piece = nullptr;
    TriMesh pieceMm;                       // the piece in mm, unrefined (its vertex i is vertex i of the refined mesh)
    RefinedMesh refined;
    TriMesh mesh;                          // refined.mesh
    std::vector<TriMesh> otherMeshes;
    std::vector<Vec3> otherCentroid;
    std::vector<std::unique_ptr<MeshQuery>> others;
    double mmPerUnit = 1;
    bool hidden(const Vec3& p) const {
        for (const std::unique_ptr<MeshQuery>& q : others)
            if (q->inside(p)) return true;
        return false;
    }
};

bool closedMesh(const TriMesh& m) {
    std::map<std::pair<int, int>, int> use;
    for (const std::array<int, 3>& tr : m.tri)
        for (int k = 0; k < 3; ++k) use[{std::min(tr[k], tr[(k + 1) % 3]), std::max(tr[k], tr[(k + 1) % 3])}] += 1;
    for (const std::pair<const std::pair<int, int>, int>& e : use)
        if (e.second != 2) return false;
    return true;
}

int eulerCharacteristic(const TriMesh& m) {
    std::set<std::pair<int, int>> edges;
    for (const std::array<int, 3>& t : m.tri)
        for (int k = 0; k < 3; ++k) edges.insert({std::min(t[k], t[(k + 1) % 3]), std::max(t[k], t[(k + 1) % 3])});
    return static_cast<int>(m.pos.size()) - static_cast<int>(edges.size()) + static_cast<int>(m.tri.size());
}

void preparePieceField(const ProposalRequest& r, PieceField& f) {
    const SurfacePiece& piece = r.pieces->pieces[r.piece];
    f.piece = &piece;
    f.mmPerUnit = r.mmPerUnit;
    f.pieceMm = piece.mesh;
    f.pieceMm.transform(r.mmPerUnit, Vec3{0, 0, 0});
    f.refined = prepareReference(f.pieceMm, refinementLevels(f.pieceMm, 0.5 * r.h), r.preparation);
    f.mesh = f.refined.mesh;
    for (int i = 0; i < static_cast<int>(r.pieces->pieces.size()); ++i) {
        const SurfacePiece& o = r.pieces->pieces[i];
        if (i == r.piece || o.node != piece.node || !closedMesh(o.mesh)) continue;
        TriMesh om = o.mesh;
        om.transform(r.mmPerUnit, Vec3{0, 0, 0});
        f.otherMeshes.push_back(om);
        f.otherCentroid.push_back(om.centroid());
    }
    for (const TriMesh& om : f.otherMeshes) f.others.push_back(std::make_unique<MeshQuery>(om));
}

double pointValueOn(const TriMesh& m, const std::vector<double>& t, const SurfacePoint& p) {
    const std::array<int, 3>& tr = m.tri[p.tri];
    return p.bary[0] * t[tr[0]] + p.bary[1] * t[tr[1]] + p.bary[2] * t[tr[2]];
}

SurfacePoint vertexPt(const TriMesh& m, int v) {
    SurfacePoint p;
    for (int t = 0; t < static_cast<int>(m.tri.size()); ++t)
        for (int k = 0; k < 3; ++k)
            if (m.tri[t][k] == v) {
                p.tri = t;
                p.bary[0] = p.bary[1] = p.bary[2] = 0;
                p.bary[k] = 1;
                return p;
            }
    return p;
}

// A control point strictly inside triangle tri, 2 % of the way from p towards the triangle's centroid: a point on a vertex
// or an edge belongs to several triangles, and a marking path traced from it starts in the one it names.
SurfacePoint interior(const TriMesh& m, int tri, const Vec3& p) {
    const Vec3 c = (m.pos[m.tri[tri][0]] + m.pos[m.tri[tri][1]] + m.pos[m.tri[tri][2]]) / 3.0;
    const Vec3 q = 0.98 * p + 0.02 * c;
    const TrianglePoint tp = closestOnTriangle(q, m.pos[m.tri[tri][0]], m.pos[m.tri[tri][1]], m.pos[m.tri[tri][2]]);
    SurfacePoint s;
    s.tri = tri;
    for (int k = 0; k < 3; ++k) s.bary[k] = tp.bary[k];
    return s;
}

bool field(const TriMesh& m, const SurfacePoint& s, std::vector<double>& t) {
    const Vec3 x = surfacePosition(m, s);
    std::vector<int> seeds;
    std::vector<double> values;
    for (int k = 0; k < 3; ++k) {
        seeds.push_back(m.tri[s.tri][k]);
        values.push_back(norm(m.pos[m.tri[s.tri][k]] - x));
    }
    TimeFieldReport rep;
    std::string err;
    return geodesicDistance(m, seeds, t, rep, &err, &values);
}

// The finish level: the first level whose contour lies entirely inside another piece, plus half a row.
double finishLevel(const PieceField& f, const std::vector<double>& t, double h, double& fullyHidden) {
    const double tMax = *std::max_element(t.begin(), t.end());
    fullyHidden = -1;
    for (double L = 0.5 * h; L < tMax; L += 0.25) {
        const RingLevel lv = extractLevel(f.mesh, t, L);
        if (lv.contours.size() != 1 || !lv.contours.front().closed) continue;
        bool all = true;
        for (const RingPoint& p : lv.contours.front().points) all = all && f.hidden(p.p);
        if (all) {
            fullyHidden = L;
            break;
        }
    }
    return fullyHidden < 0 ? -1 : std::min(fullyHidden + 0.5 * h, tMax - h);
}

std::vector<std::set<int>> vertexNeighbours(const TriMesh& m) {
    std::vector<std::set<int>> nb(m.pos.size());
    for (const std::array<int, 3>& tr : m.tri)
        for (int k = 0; k < 3; ++k) {
            nb[tr[k]].insert(tr[(k + 1) % 3]);
            nb[tr[k]].insert(tr[(k + 2) % 3]);
        }
    return nb;
}

// Steepest descent of t from vertex v until a seed vertex, as the chain of vertices from v.
std::vector<int> descentChain(const TriMesh& m, const std::vector<double>& t, int v, const std::set<int>& seeds) {
    const std::vector<std::set<int>> nb = vertexNeighbours(m);
    std::vector<int> chain{v};
    for (int guard = 0; guard < 100000 && !seeds.count(v); ++guard) {
        int next = -1;
        double slope = 0;
        for (int w : nb[v]) {
            const double s = (t[v] - t[w]) / norm(m.pos[v] - m.pos[w]);
            if (s > slope) {
                slope = s;
                next = w;
            }
        }
        if (next < 0) break;
        v = next;
        chain.push_back(v);
    }
    return chain;
}

// A chain (reversed, from its seed end) as path control points: near every second chain edge's midpoint, inside one of
// its triangles, then the chain's end vertex.
void chainControlPoints(const TriMesh& m, std::vector<int> chain, std::vector<SurfacePoint>& out) {
    std::reverse(chain.begin(), chain.end());
    for (size_t i = 0; i + 1 < chain.size(); i += 2) {
        const int u = chain[i], v = chain[i + 1];
        int tri = -1;
        for (int t2 = 0; t2 < static_cast<int>(m.tri.size()) && tri < 0; ++t2) {
            int hits = 0;
            for (int k = 0; k < 3; ++k) hits += (m.tri[t2][k] == u || m.tri[t2][k] == v) ? 1 : 0;
            if (hits == 2) tri = t2;
        }
        out.push_back(interior(m, tri, 0.5 * (m.pos[u] + m.pos[v])));
    }
    const int last = chain.back();
    out.push_back(interior(m, vertexPt(m, last).tri, m.pos[last]));
}

// The row-start path: steepest descent of t from the given point of the level's contour to the start's triangle.
std::vector<SurfacePoint> descentFromContourPoint(const TriMesh& m, const std::vector<double>& t, const SurfacePoint& start, double level, const RingPoint& rp) {
    const int v = t[rp.a] >= level ? rp.a : rp.b;
    const std::set<int> seeds(m.tri[start.tri].begin(), m.tri[start.tri].end());
    std::vector<SurfacePoint> out{start};
    chainControlPoints(m, descentChain(m, t, v, seeds), out);
    return out;
}

// The lab's row-start path: from the finish contour's inner (or outer) point, facing the hiding piece, to the start.
std::vector<SurfacePoint> descentPath(const PieceField& f, const std::vector<double>& t, const SurfacePoint& start, double level, bool inner) {
    const TriMesh& m = f.mesh;
    const RingLevel lv = extractLevel(m, t, level);
    const Contour& c = lv.contours.front();
    Vec3 centroid{0, 0, 0};
    for (const RingPoint& p : c.points) centroid += p.p;
    centroid = centroid / double(c.points.size());
    // the hiding piece: the one containing most of the contour
    int best = -1, bestCount = 0;
    for (size_t i = 0; i < f.others.size(); ++i) {
        int n = 0;
        for (const RingPoint& p : c.points) n += f.others[i]->inside(p.p) ? 1 : 0;
        if (n > bestCount) {
            bestCount = n;
            best = static_cast<int>(i);
        }
    }
    const Vec3 axis = centroid - surfacePosition(m, start);
    const Vec3 ax = axis / norm(axis);
    Vec3 inward = (best >= 0 ? f.otherCentroid[best] : centroid) - centroid;
    inward = inward - dot(inward, ax) * ax;
    inward = inward / norm(inward);
    int pick = 0;
    double score = -1e300;
    for (size_t i = 0; i < c.points.size(); ++i) {
        const double s = dot(c.points[i].p - centroid, inward) * (inner ? 1 : -1);
        if (s > score) {
            score = s;
            pick = static_cast<int>(i);
        }
    }
    return descentFromContourPoint(m, t, start, level, c.points[pick]);
}

// A row-start path through the rows on one side: a control point on every mid-row level (k + 1/2) h', k = 0 .. R - 1, at
// the contour point furthest along `towards` (a side of a bend), or, without a direction, the point nearest the previous
// one, starting beside `startNear`. Every row k h' (k = 1 .. R - 1) lies between two consecutive points, so the path
// crosses each row once and stays on its side (a steepest descent drifts towards the shorter side of a bend).
std::vector<SurfacePoint> sidePath(const TriMesh& m, const std::vector<double>& t, int R, double hp, const Vec3* towards, const Vec3& startNear) {
    std::vector<SurfacePoint> out;
    Vec3 previous = startNear;
    // one point in the middle of every row; a piece of a single row (the rim strip of a cut ear, rim_cut.py) takes two, a
    // quarter and three quarters up, because a path needs two points
    const int points = R == 1 ? 2 : R;
    for (int k = 0; k < points; ++k) {
        const RingLevel lv = extractLevel(m, t, R == 1 ? (0.25 + 0.5 * k) * hp : (k + 0.5) * hp);
        if (lv.contours.size() != 1 || !lv.contours.front().closed) continue;
        const Contour& c = lv.contours.front();
        Vec3 centroid{0, 0, 0};
        for (const RingPoint& p : c.points) centroid += p.p;
        centroid = centroid / double(c.points.size());
        int pick = 0;
        double best = -1e300;
        for (int i = 0; i < static_cast<int>(c.points.size()); ++i) {
            const double score = towards ? dot(c.points[i].p - centroid, *towards) : -norm(c.points[i].p - previous);
            if (score > best && c.points[i].tri >= 0) {
                best = score;
                pick = i;
            }
        }
        previous = c.points[pick].p;
        out.push_back(interior(m, c.points[pick].tri, c.points[pick].p));
    }
    return out;
}

// A closed curve at a level: control points evenly spaced along the contour (empty when the level is not one contour).
std::vector<SurfacePoint> levelCurve(const TriMesh& m, const std::vector<double>& t, double level, int points) {
    std::vector<SurfacePoint> out;
    const RingLevel lv = extractLevel(m, t, level);
    if (lv.contours.size() != 1) return out;
    const Contour& c = lv.contours.front();
    const int n = static_cast<int>(c.points.size());
    for (int j = 0; j < points; ++j) {
        const RingPoint& p = c.points[(j * n) / points];
        out.push_back(interior(m, p.tri, p.p));
    }
    return out;
}

// The finish point: on the path where t crosses the level.
SurfacePoint onPathAtLevel(const TriMesh& m, const std::vector<double>& t, const std::vector<SurfacePoint>& path, double level) {
    for (size_t i = 0; i + 1 < path.size(); ++i) {
        const double a = pointValueOn(m, t, path[i]), b = pointValueOn(m, t, path[i + 1]);
        if (a < level && b >= level) {
            const double u = (level - a) / (b - a);
            const Vec3 x = (1 - u) * surfacePosition(m, path[i]) + u * surfacePosition(m, path[i + 1]);
            return closestSurfacePoint(m, x);
        }
    }
    return path.back();
}

double ringLength(const TriMesh& m, const std::vector<double>& t, double level) {
    const RingLevel lv = extractLevel(m, t, level);
    return lv.contours.size() == 1 ? lv.contours.front().length : 0.0;
}

// The pole ring: the smallest count the first strip's increases fit (adjacent shaping: at most n1 - 2 per strip).
int poleRing(const TriMesh& m, const std::vector<double>& t, double finish, double w, double h, int* automatic = nullptr, int* second = nullptr) {
    const int R = std::max(1, static_cast<int>(std::lround(finish / h)));
    const double hp = finish / R;
    const int auto1 = std::max(3, static_cast<int>(std::lround(ringLength(m, t, hp) / w)));
    const int n2 = static_cast<int>(std::lround(ringLength(m, t, 2 * hp) / w));
    int ring1 = auto1;
    while (n2 - ring1 > ring1 - 2) ++ring1;
    if (automatic) *automatic = auto1;
    if (second) *second = n2;
    return ring1;
}

// A point of the refined mesh as a point of the file (the refined split maps exactly back to the piece's triangle).
SurfacePoint toFile(const PieceField& f, const SurfacePoint& q) {
    const SurfacePoint lp = fromRefined(f.refined, q);
    SurfacePoint fp;
    fp.tri = f.piece->fileTriangle[lp.tri];
    for (int k = 0; k < 3; ++k) fp.bary[k] = lp.bary[k];
    return fp;
}

SurfaceMarking offered(const PieceField& f, int id, MarkingMeaning meaning, MarkingShape shape, const std::vector<SurfacePoint>& refinedPoints, int group,
                       const std::string& note, const std::string& reason) {
    SurfaceMarking mk;
    mk.id = id;
    mk.meaning = meaning;
    mk.shape = shape;
    for (const SurfacePoint& q : refinedPoints) mk.points.push_back(toFile(f, q));
    mk.origin = MarkingOrigin::Proposal;
    mk.proposal = ProposalState::Offered;
    mk.proposalId = id;
    mk.alternativeGroup = group;
    mk.note = note;
    mk.reason = reason;
    return mk;
}

// Rule 1: the end of the piece's geodesic diameter with the smaller hidden share nearby (a vertex of the refined mesh).
bool ruleOneStart(const PieceField& f, int& startVertex, std::string* err) {
    const TriMesh& m = f.mesh;
    std::vector<double> t, tB;
    TimeFieldReport tr;
    if (!geodesicDistance(m, {0}, t, tr, err)) return false;
    const int endA = static_cast<int>(std::max_element(t.begin(), t.end()) - t.begin());
    if (!geodesicDistance(m, {endA}, t, tr, err)) return false;
    const int endB = static_cast<int>(std::max_element(t.begin(), t.end()) - t.begin());
    if (!geodesicDistance(m, {endB}, tB, tr, err)) return false;
    int nA = 0, hA = 0, nB = 0, hB = 0;
    for (size_t v = 0; v < m.pos.size(); ++v) {
        if (t[v] < kNearEnd) {
            ++nA;
            hA += f.hidden(m.pos[v]);
        }
        if (tB[v] < kNearEnd) {
            ++nB;
            hB += f.hidden(m.pos[v]);
        }
    }
    startVertex = (nA ? double(hA) / nA : 0) <= (nB ? double(hB) / nB : 0) ? endA : endB;
    return true;
}

// Where an unconstrained relaxation of a pole layout settles the pole: one build (Relaxed) of the three markings as a
// labelled research preview, with nothing holding the pole. Position in mm on the refined mesh.
bool trialPole(const ProposalRequest& r, const PieceField& f, const std::vector<SurfaceMarking>& three, Vec3& settled, std::string* err) {
    MarkingSet probe;
    probe.digest = r.fileDigest;
    probe.preparation = toString(r.preparation);
    probe.markings = three;
    MarkedBuildRequest rq;
    rq.pieces = r.pieces;
    rq.piece = r.piece;
    rq.markings = &probe;
    rq.mmPerUnit = r.mmPerUnit;
    rq.w = r.w;
    rq.h = r.h;
    for (const SurfaceMarking& m : three) rq.previewOffered.push_back(m.id);
    rq.options.anchors = RelaxAnchors::None;   // the probe: where would the pole go if nothing held it
    rq.options.borderRings = BorderRings::SlideInRelaxation;
    rq.options.preparation = r.preparation;
    MarkedBuildProduct out;
    if (!buildFromMarkings(rq, out, err) || !out.evaluation.relaxed) {
        if (err && err->empty()) {
            *err = "the trial build did not relax";
            for (const MarkingFinding& x : out.review.findings)
                if (x.level == FindingLevel::Problem) *err += ": " + x.message;
        }
        return false;
    }
    settled = out.evaluation.after.mesh.pos[out.evaluation.after.rings[0][0]];
    (void) f;
    return true;
}

// The area outside [lo, hi] of t (per refined triangle by its centroid's value), split into hidden and visible [mm2].
void areaOutside(const PieceField& f, const std::vector<double>& t, double lo, double hi, double& visible, double& hidden) {
    visible = hidden = 0;
    const TriMesh& m = f.mesh;
    for (const std::array<int, 3>& tr : m.tri) {
        const double tc = (t[tr[0]] + t[tr[1]] + t[tr[2]]) / 3.0;
        if (tc >= lo && tc <= hi) continue;
        const Vec3 c = (m.pos[tr[0]] + m.pos[tr[1]] + m.pos[tr[2]]) / 3.0;
        const double a = 0.5 * norm(cross(m.pos[tr[1]] - m.pos[tr[0]], m.pos[tr[2]] - m.pos[tr[0]]));
        (f.hidden(c) ? hidden : visible) += a;
    }
}

double loopLength(const TriMesh& m, const std::vector<int>& loop) {
    double L = 0;
    for (size_t i = 0; i < loop.size(); ++i) L += norm(m.pos[loop[(i + 1) % loop.size()]] - m.pos[loop[i]]);
    return L;
}

Vec3 loopCentroid(const TriMesh& m, const std::vector<int>& loop) {
    Vec3 c{0, 0, 0};
    for (int v : loop) c += m.pos[v];
    return c / double(std::max<size_t>(1, loop.size()));
}

// A border loop of the piece (piece vertex ids in walking order) as closed-curve control points on the piece's own
// triangles: the midpoints of evenly spaced border edges, weights (0.5, 0.5, 0) on the triangle that owns the edge. These
// lie exactly on the border, so the mesher takes the loop itself as the cast-on or bind-off (pointsOnOneBorderLoop).
std::vector<SurfacePoint> borderLoopPoints(const PieceField& f, const std::vector<int>& loop) {
    const TriMesh& m = f.piece->mesh;
    std::vector<SurfacePoint> out;
    const int n = static_cast<int>(loop.size());
    const int count = std::min(kCurvePoints, n);
    for (int j = 0; j < count; ++j) {
        const int i = (j * n) / count;
        const int a = loop[i], b = loop[(i + 1) % n];
        for (int t = 0; t < static_cast<int>(m.tri.size()); ++t) {
            int ka = -1, kb = -1;
            for (int k = 0; k < 3; ++k) {
                if (m.tri[t][k] == a) ka = k;
                if (m.tri[t][k] == b) kb = k;
            }
            if (ka < 0 || kb < 0 || (ka + 1) % 3 != kb) continue;   // the triangle walking the border edge a -> b
            SurfacePoint p;
            p.tri = f.piece->fileTriangle[t];
            p.bary[0] = p.bary[1] = p.bary[2] = 0;
            p.bary[ka] = 0.5;
            p.bary[kb] = 0.5;
            out.push_back(p);
            break;
        }
    }
    return out;
}

// Markings made directly on file triangles (a border loop's points), not mapped from the refined mesh.
SurfaceMarking offeredOnFile(int id, MarkingMeaning meaning, MarkingShape shape, const std::vector<SurfacePoint>& filePoints, const std::string& note,
                             const std::string& reason) {
    SurfaceMarking mk;
    mk.id = id;
    mk.meaning = meaning;
    mk.shape = shape;
    mk.points = filePoints;
    mk.origin = MarkingOrigin::Proposal;
    mk.proposal = ProposalState::Offered;
    mk.proposalId = id;
    mk.note = note;
    mk.reason = reason;
    return mk;
}

// The refined mesh's border loop that contains the piece vertex v (vertex ids are shared: refinement appends midpoints).
std::vector<int> refinedLoopThrough(const TriMesh& refined, int v) {
    for (const std::vector<int>& loop : borderLoops(refined))
        for (int x : loop)
            if (x == v) return loop;
    return {};
}

std::vector<Vec3> fileUnits(const std::vector<Vec3>& mm, double mmPerUnit) {
    std::vector<Vec3> out;
    for (const Vec3& p : mm) out.push_back(p / mmPerUnit);
    return out;
}

std::vector<Vec3> contourSample(const Contour& c, int count) {
    std::vector<Vec3> out;
    const int n = static_cast<int>(c.points.size());
    for (int j = 0; j < std::min(count, n); ++j) out.push_back(c.points[(j * n) / std::min(count, n)].p);
    return out;
}

void addFinding(ProposalSet& out, FindingLevel level, FindingCause cause, const std::string& message, int row = -1, const std::vector<Vec3>& at = {}) {
    MarkingFinding f;
    f.level = level;
    f.cause = cause;
    f.message = message;
    f.row = row;
    f.at = at;
    out.findings.push_back(f);
}

std::string stitches(int n) { return text("%d stitch%s", n, n == 1 ? "" : "es"); }

// The known narrow-finish limitation, stated when a finish ring has fewer than kWideFinishStitches ("" otherwise).
std::string narrowFinishText(int finishStitches) {
    if (finishStitches >= kWideFinishStitches) return "";
    return text("The finish ring has only about %s: a bind-off this narrow may leave decreases without a place (a known limitation at narrow ends). ",
                stitches(finishStitches).c_str());
}

// -------------------------------------------------------------------------------------------------------------------
// open tube: two border loops

void proposeOpenTube(const ProposalRequest& r, const PieceField& f, const std::vector<std::vector<int>>& loops, ProposalSet& out) {
    const TriMesh& m = f.mesh;
    const double L[2] = {loopLength(f.pieceMm, loops[0]), loopLength(f.pieceMm, loops[1])};
    const bool sized = std::fabs(L[0] - L[1]) > r.w;   // the openings differ by more than a stitch: name them by size
    auto openingName = [&](int i) -> std::string {
        if (!sized) return i == 0 ? "opening 1" : "opening 2";
        return L[i] > L[1 - i] ? "the wider opening" : "the narrower opening";
    };
    int next = 1;
    for (int d = 0; d < 2; ++d) {
        const std::vector<int>& A = loops[d];
        const std::vector<int>& B = loops[1 - d];
        const std::vector<int> seedLoop = refinedLoopThrough(m, A[0]), endLoop = refinedLoopThrough(m, B[0]);
        std::vector<double> t;
        TimeFieldReport tr;
        std::string e;
        if (seedLoop.empty() || endLoop.empty() || !geodesicDistance(m, seedLoop, t, tr, &e)) continue;
        double T = 1e300;
        for (int v : endLoop) T = std::min(T, t[v]);
        const int R = std::max(1, static_cast<int>(std::lround(T / r.h)));
        const double hp = T / R;
        // every interior row one closed contour; the contours' centroids say whether the tube is bent
        const Vec3 cA = loopCentroid(m, seedLoop), cB = loopCentroid(m, endLoop);
        const Vec3 chord = cB - cA;
        double bend = 0;
        Vec3 bendDir{0, 0, 0};
        int splitRow = -1;
        double lastRowLength = L[d];   // the last interior row's length (the cast-on's for a one-row tube)
        for (int k = 1; k < R; ++k) {
            RingLevel lv = extractLevel(m, t, k * hp);
            dropContoursShorterThan(lv, kKnitOverStitches * r.w);   // a bump too small for a round is knitted over, as the layout does
            if (lv.contours.size() != 1 || !lv.contours.front().closed) {
                if (splitRow < 0) splitRow = k;
                continue;
            }
            lastRowLength = lv.contours.front().length;
            Vec3 c{0, 0, 0};
            for (const RingPoint& p : lv.contours.front().points) c += p.p;
            c = c / double(lv.contours.front().points.size());
            const Vec3 off = (c - cA) - (dot(c - cA, chord) / std::max(1e-12, dot(chord, chord))) * chord;
            if (norm(off) > bend) {
                bend = norm(off);
                bendDir = off / norm(off);
            }
        }
        if (splitRow >= 0) {
            addFinding(out, FindingLevel::Problem, FindingCause::AlgorithmLimitation,
                       text("Starting at %s, row %d splits into several pieces: this direction is not offered.", openingName(d).c_str(), splitRow), splitRow);
            continue;
        }
        const bool bent = bend > 0.5 * r.w;
        const int nA = static_cast<int>(std::lround(L[d] / r.w)), nB = static_cast<int>(std::lround(L[1 - d] / r.w));
        // the layout refuses a strip that changes by more than one stitch per stitch of its shorter ring (guided_rows.cpp):
        // that happens when the finish opening is far from a row of this start's field, so the direction is not offered
        const int nLast = static_cast<int>(std::lround(lastRowLength / r.w));
        if (std::abs(nLast - nB) > std::min(nLast, nB)) {
            addFinding(out, FindingLevel::Problem, FindingCause::AlgorithmLimitation,
                       text("Starting at %s, %s lies far from a row: the last row would go from about %d to %d stitches, more than a row can shape. "
                            "This direction is not offered; cut the piece so that the opening follows a row, or draw your own start and finish.",
                            openingName(d).c_str(), openingName(1 - d).c_str(), nLast, nB),
                       R - 1);
            continue;
        }
        const std::vector<SurfacePoint> startPts = borderLoopPoints(f, A), finishPts = borderLoopPoints(f, B);
        const int sides = bent ? 2 : 1;
        for (int side = 0; side < sides; ++side) {
            // the inner side of the bend lies towards its centre (against bendDir, the rows' offset from the chord), the
            // outer side away from it; a straight tube's path starts beside the start loop's first vertex
            const Vec3 towards = side == 0 ? -1.0 * bendDir : bendDir;
            const std::vector<SurfacePoint> path = sidePath(m, t, R, hp, bent ? &towards : nullptr, m.pos[A[0]]);
            std::vector<SurfacePoint> pathFile;
            for (const SurfacePoint& q : path) pathFile.push_back(toFile(f, q));
            const std::string sideWord = !bent ? "any" : side == 0 ? "inner" : "outer";
            ProposalSetup s;
            s.setup = next++;
            s.recipe = "tube/" + std::to_string(d) + "/" + sideWord;
            s.title = "Start at " + openingName(d) + ", finish at " + openingName(1 - d) + (bent ? (side == 0 ? " (row starts on the inside of the bend)" : " (row starts on the outside of the bend)") : "");
            s.purpose = "Knit the whole tube from " + openingName(d) + " to " + openingName(1 - d) + "; nothing is left out.";
            s.reason = text("The piece is an open tube: its two openings are the cast-on (about %s) and the bind-off (about %s), so no cap is needed. ", stitches(nA).c_str(),
                            stitches(nB).c_str()) +
                       (bent ? (side == 0 ? "The row starts follow the inside of the bend." : "The row starts follow the outside of the bend.")
                             : "The tube is nearly straight, so any side works for the row starts; this path begins at the opening's first vertex.");
            s.limitation = std::string(bent && side == 0 ? "The inside of the bend is the shorter side, where rows are closest together. " : "") +
                           "Row and stitch counts here are estimates from the piece's geometry; the layout build gives the actual rows and shaping, and only "
                           "the Knit step settles the yarn.";
            s.rowsAbout = R;
            s.firstRingStitches = nA;
            s.lastRingStitches = nB;
            s.markings = {offeredOnFile(1, MarkingMeaning::Start, MarkingShape::ClosedCurve, startPts, "cast-on: " + openingName(d), "a border loop of the piece"),
                          offeredOnFile(2, MarkingMeaning::Finish, MarkingShape::ClosedCurve, finishPts, "bind-off: " + openingName(1 - d), "a border loop of the piece"),
                          offeredOnFile(3, MarkingMeaning::RowStartPath, MarkingShape::Path, pathFile, "row starts (" + sideWord + " side)",
                                        "a point in the middle of every row, on that side, from the cast-on to the bind-off")};
            out.setups.push_back(s);
        }
    }
    out.shape = PieceShape::OpenTube;
    out.shapeText = "An open tube with two openings";
}

// -------------------------------------------------------------------------------------------------------------------
// closed pieces

// The pole-start setup texts and markings common to the hidden and isolated cases.
ProposalSetup poleSetup(const PieceField& f, const SurfacePoint& start, int ring1, const SurfaceMarking& finish, const std::vector<SurfacePoint>& path,
                        const std::string& pathNote) {
    ProposalSetup s;
    SurfaceMarking st = offered(f, 1, MarkingMeaning::Start, MarkingShape::Point, {start}, 0, "start point (pole start)", "the free end of the piece's longest geodesic path");
    st.settings = {{"ring1", std::to_string(ring1)}};
    s.markings = {st, finish, offered(f, 3, MarkingMeaning::RowStartPath, MarkingShape::Path, path, 0, pathNote, "steepest descent of the field from the finish to the start")};
    s.refinable = true;
    return s;
}

// -------------------------------------------------------------------------------------------------------------------
// one opening: a cup (an ear, a head, a sleeve with a closed end)

// The cup's start vertex: the one from which the opening is most nearly one row, i.e. whose distances to the opening's
// points vary least relative to their mean (sampled at up to kCupStartSamples points of the opening, one field each),
// among the vertices at least half as deep as the deepest (fromOpening: the distance from the whole opening). On a
// round cup it is the bottom; on a bent or lopsided one (a bunny's body with the neck on top) it moves to where the
// rows reach the opening evenly, instead of the farthest point, whose last row would have to absorb the slant.
int cupStart(const TriMesh& m, const std::vector<int>& opening, const std::vector<double>& fromOpening) {
    constexpr int kCupStartSamples = 24;
    const int nv = static_cast<int>(m.pos.size());
    const double deepest = *std::max_element(fromOpening.begin(), fromOpening.end());
    std::vector<double> lo(nv, 1e300), hi(nv, 0), sum(nv, 0);
    const int samples = std::min(kCupStartSamples, static_cast<int>(opening.size()));
    for (int k = 0; k < samples; ++k) {
        std::vector<double> t;
        TimeFieldReport tr;
        std::string e;
        if (!geodesicDistance(m, {opening[k * opening.size() / samples]}, t, tr, &e)) continue;
        for (int v = 0; v < nv; ++v) {
            lo[v] = std::min(lo[v], t[v]);
            hi[v] = std::max(hi[v], t[v]);
            sum[v] += t[v];
        }
    }
    int best = static_cast<int>(std::max_element(fromOpening.begin(), fromOpening.end()) - fromOpening.begin());
    double bestScore = 1e300;
    for (int v = 0; v < nv; ++v) {
        if (fromOpening[v] < 0.5 * deepest || sum[v] <= 0) continue;
        const double score = (hi[v] - lo[v]) / (sum[v] / samples);
        if (score < bestScore) {
            bestScore = score;
            best = v;
        }
    }
    return best;
}

// A pole start at the cup's bottom, the opening itself as the bind-off: the whole piece knitted, nothing left out. The
// opening must lie close to a row of the start's field, otherwise the
// last row would have to change by more than a row can shape (guided_rows.cpp) and nothing is offered. The start: cupStart.
void proposeCup(const ProposalRequest& r, const PieceField& f, const std::vector<int>& loop, ProposalSet& out) {
    const TriMesh& m = f.mesh;
    out.shape = PieceShape::Cup;
    out.shapeText = "A piece with one opening (a cup)";
    const std::vector<int> opening = refinedLoopThrough(m, loop[0]);
    std::vector<double> fromOpening;
    TimeFieldReport tr;
    std::string e;
    if (opening.empty() || !geodesicDistance(m, opening, fromOpening, tr, &e)) {
        addFinding(out, FindingLevel::Problem, FindingCause::Input, "No field from the opening: " + e);
        return;
    }
    const int deepest = cupStart(m, opening, fromOpening);
    const SurfacePoint start = interior(m, vertexPt(m, deepest).tri, m.pos[deepest]);
    std::vector<double> t;
    if (!field(m, start, t)) {
        addFinding(out, FindingLevel::Problem, FindingCause::Input, "No field from the cup's bottom.");
        return;
    }
    // the rows of the build (OpeningRows::OpeningIsALevel, marked_build.h): t' = T t / (t + d), the opening a level at the
    // mean of its distances from the start
    double level = 0;
    for (int v : opening) level += t[v];
    level /= opening.size();
    for (size_t v = 0; v < t.size(); ++v) t[v] = t[v] + fromOpening[v] > 0 ? level * t[v] / (t[v] + fromOpening[v]) : 0.0;
    double T = 1e300;
    for (int v : opening) T = std::min(T, t[v]);
    const int R = std::max(1, static_cast<int>(std::lround(T / r.h)));
    const double hp = T / R;
    const int nOpening = static_cast<int>(std::lround(loopLength(f.pieceMm, loop) / r.w));
    const int nLast = R >= 2 ? static_cast<int>(std::lround(ringLength(m, t, (R - 1) * hp) / r.w)) : nOpening;
    if (std::abs(nLast - nOpening) > std::min(nLast, nOpening)) {
        addFinding(out, FindingLevel::Problem, FindingCause::AlgorithmLimitation,
                   text("Starting at the cup's bottom, the opening lies far from a row: the last row would go from about %d to %d stitches, "
                        "more than a row can shape. Nothing is offered; cut the piece so that the opening follows a row, or draw your own start and finish.",
                        nLast, nOpening),
                   R - 1);
        return;
    }
    const int ring1 = poleRing(m, t, T, r.w, r.h);
    std::vector<SurfacePoint> path{start};
    for (const SurfacePoint& q : sidePath(m, t, R, hp, nullptr, surfacePosition(m, start))) path.push_back(q);
    const SurfaceMarking finish = offeredOnFile(2, MarkingMeaning::Finish, MarkingShape::ClosedCurve, borderLoopPoints(f, loop), "bind-off: the opening",
                                                "the piece's only border loop");
    ProposalSetup s = poleSetup(f, start, ring1, finish, path, "row starts");
    s.markings[2].reason = "a point in the middle of every row, on one side, from the start to the opening";
    s.recipe = "cup/pole";
    s.title = "Start at the bottom of the cup, finish at its opening";
    s.purpose = "Knit the whole piece from its bottom to the opening; nothing is left out.";
    s.reason = text("The piece has one opening, so it is knitted like a cup: a pole start at its bottom (where the rows reach the opening most evenly) and the opening "
                    "as the bind-off (about %s), about %d rows.",
                    stitches(nOpening).c_str(), R);
    s.limitation = "A start point uses the gathered pole start, an experimental construction. " + narrowFinishText(nOpening) +
                   "Row and stitch counts here are estimates from the piece's geometry; the layout build gives the actual rows and shaping.";
    s.rowsAbout = R;
    s.firstRingStitches = ring1;
    s.lastRingStitches = nOpening;
    out.setups.push_back(s);
}

void proposeClosed(const ProposalRequest& r, const PieceField& f, ProposalSet& out) {
    const TriMesh& m = f.mesh;
    int startVertex = -1;
    std::string e;
    if (!ruleOneStart(f, startVertex, &e)) {
        addFinding(out, FindingLevel::Problem, FindingCause::Input, "No field on this piece: " + e);
        return;
    }
    const SurfacePoint s1 = interior(m, vertexPt(m, startVertex).tri, m.pos[startVertex]);
    std::vector<double> t1;
    field(m, s1, t1);
    const double tMax = *std::max_element(t1.begin(), t1.end());
    double fullyHidden = 0;
    const double f1 = finishLevel(f, t1, r.h, fullyHidden);
    if (f1 >= 0) {
        // ----- hidden at one end (the arm, the leg)
        out.shape = PieceShape::ClosedHidden;
        out.shapeText = "A closed piece with one end inside a neighbouring piece";
        const int ring1 = poleRing(m, t1, f1, r.w, r.h);
        const int R = std::max(1, static_cast<int>(std::lround(f1 / r.h)));
        const double hp = f1 / R;
        const int nFinish = static_cast<int>(std::lround(ringLength(m, t1, f1) / r.w));
        double visible = 0, hidden = 0;
        areaOutside(f, t1, -1, f1, visible, hidden);
        const bool narrow = nFinish < kWideFinishStitches;
        const std::string narrowText = narrowFinishText(nFinish);
        const std::string poleText = "A start point uses the gathered pole start, an experimental construction. ";
        for (int side = 0; side < 2; ++side) {
            const bool inner = side == 0;
            const std::vector<SurfacePoint> p = descentPath(f, t1, s1, f1, inner);
            const SurfacePoint fin = onPathAtLevel(m, t1, p, f1);
            const SurfaceMarking finish = offered(f, 2, MarkingMeaning::Finish, MarkingShape::Point, {fin}, 0, "finish level inside the neighbouring piece",
                                                  "the first level entirely inside the neighbouring piece, plus half a row");
            ProposalSetup s = poleSetup(f, s1, ring1, finish, p, inner ? "row starts facing the neighbouring piece" : "row starts facing away from the neighbouring piece");
            s.recipe = inner ? "pole/inner" : "pole/outer";
            s.title = std::string("Start at the free end, finish inside the neighbouring piece") + (inner ? " (row starts on the hidden side)" : " (row starts on the visible side)");
            s.purpose = "Knit the piece from its free end; the end that lies inside the neighbouring piece is left out where it cannot be seen.";
            s.reason = text("The free end is the end of the piece's longest geodesic path furthest from the neighbouring piece. The finish is where the piece first lies "
                            "entirely inside the neighbouring piece, plus half a row: about %d rows. The row starts run %s.",
                            R, inner ? "on the side facing the neighbouring piece, where the helix seam is least visible" : "on the side facing away from the neighbouring piece");
            s.limitation = poleText + narrowText + text("Left out: %.0f mm2 hidden, %.0f mm2 visible (measured on the piece).", hidden, visible);
            s.rowsAbout = R;
            s.firstRingStitches = ring1;
            s.lastRingStitches = nFinish;
            s.omittedVisibleMm2 = visible;
            s.omittedHiddenMm2 = hidden;
            s.partial = visible > 1;
            out.setups.push_back(s);
        }
        // a band: the start curve 4 rows in (the cap inside it left out on purpose), the same finish as a curve
        const std::vector<SurfacePoint> startCurve = levelCurve(m, t1, 4 * hp, kCurvePoints), finishCurve = levelCurve(m, t1, f1, kCurvePoints);
        if (R > 6 && startCurve.size() == static_cast<size_t>(kCurvePoints) && finishCurve.size() == static_cast<size_t>(kCurvePoints)) {
            const std::vector<SurfacePoint> p = descentPath(f, t1, s1, f1, true);
            double bv = 0, bh = 0;
            areaOutside(f, t1, 4 * hp, f1, bv, bh);
            ProposalSetup s;
            s.recipe = "band";
            s.markings = {offered(f, 1, MarkingMeaning::Start, MarkingShape::ClosedCurve, startCurve, 0, "start curve 4 rows from the free end", "the field's level 4 rows from the free end"),
                          offered(f, 2, MarkingMeaning::Finish, MarkingShape::ClosedCurve, finishCurve, 0, "finish curve inside the neighbouring piece", "the same finish level as a curve"),
                          offered(f, 3, MarkingMeaning::RowStartPath, MarkingShape::Path, p, 0, "row starts facing the neighbouring piece", "steepest descent of the field")};
            s.title = "Cast on 4 rows from the free end, finish inside the neighbouring piece";
            s.purpose = "Knit a band without the experimental pole start: a cast-on curve instead of a start point.";
            s.reason = "A cast-on curve is a supported construction; the cap at the free end inside it is left out on purpose.";
            s.limitation = narrowText + text("The cap at the free end stays unknitted: %.0f mm2 visible and %.0f mm2 hidden are left out in all (measured on the piece).", bv, bh);
            s.rowsAbout = R - 4;
            s.firstRingStitches = static_cast<int>(std::lround(ringLength(m, t1, 4 * hp) / r.w));
            s.lastRingStitches = nFinish;
            s.omittedVisibleMm2 = bv;
            s.omittedHiddenMm2 = bh;
            s.partial = bv > 1;
            out.setups.push_back(s);
        }
        // a narrow hidden finish: a finish where the ring still has 8 stitches, and the visible surface it costs
        if (narrow) {
            double fAlt = -1;
            for (double Lv = f1; Lv > 2 * r.h; Lv -= 0.25)
                if (static_cast<int>(std::lround(ringLength(m, t1, Lv) / r.w)) >= kWideFinishStitches) {
                    fAlt = Lv;
                    break;
                }
            if (fAlt > 0) {
                const std::vector<SurfacePoint> p = descentPath(f, t1, s1, fAlt, true);
                const SurfacePoint fin = onPathAtLevel(m, t1, p, fAlt);
                const SurfaceMarking finish = offered(f, 2, MarkingMeaning::Finish, MarkingShape::Point, {fin}, 0, "finish where the ring has 8 stitches",
                                                      "the last level before the hidden finish whose ring has at least 8 stitches");
                const int nAlt = static_cast<int>(std::lround(ringLength(m, t1, fAlt) / r.w));
                double av = 0, ah = 0;
                areaOutside(f, t1, -1, fAlt, av, ah);
                ProposalSetup s = poleSetup(f, s1, poleRing(m, t1, fAlt, r.w, r.h), finish, p, "row starts facing the neighbouring piece");
                s.recipe = "pole/inner/wide";
                s.title = text("Start at the free end, finish earlier where the ring still has %s", stitches(nAlt).c_str());
                s.purpose = "Avoid the narrow bind-off inside the neighbouring piece.";
                s.reason = text("The hidden finish ring has only about %s; finishing %.0f mm earlier gives a ring of about %s.", stitches(nFinish).c_str(), f1 - fAlt,
                                stitches(nAlt).c_str());
                s.limitation = poleText + text("This costs visible surface: %.0f mm2 visible (%.0f mm2 more than the hidden finish) and %.0f mm2 hidden are left out (measured on the piece).", av,
                                               av - visible, ah);
                s.rowsAbout = std::max(1, static_cast<int>(std::lround(fAlt / r.h)));
                s.firstRingStitches = s.markings[0].settings.empty() ? 0 : std::atoi(s.markings[0].settings[0].second.c_str());
                s.lastRingStitches = nAlt;
                s.omittedVisibleMm2 = av;
                s.omittedHiddenMm2 = ah;
                s.partial = av > 1;
                out.setups.push_back(s);
            }
        }
        return;
    }
    // ----- nothing hides either end (the capsule): every setup leaves a visible cap
    out.shape = PieceShape::ClosedIsolated;
    out.shapeText = "A closed piece that nothing hides";
    const double finishAt = tMax - 4 * r.h;
    if (finishAt <= 4 * r.h) {
        addFinding(out, FindingLevel::Problem, FindingCause::AlgorithmLimitation, "This closed piece is too short for a start and a finish 4 rows from its ends; a closed finish is not available.");
        return;
    }
    const int R = std::max(1, static_cast<int>(std::lround(finishAt / r.h)));
    const double hp = finishAt / R;
    const RingLevel finishLevelContour = extractLevel(m, t1, finishAt);
    const std::vector<SurfacePoint> finishCurve = levelCurve(m, t1, finishAt, kCurvePoints), startCurve = levelCurve(m, t1, 4 * hp, kCurvePoints);
    if (finishLevelContour.contours.size() != 1 || finishCurve.size() != static_cast<size_t>(kCurvePoints)) {
        addFinding(out, FindingLevel::Problem, FindingCause::Input, "The level 4 rows before the far end is not one closed contour; no finish is proposed.");
        return;
    }
    const std::vector<SurfacePoint> p = descentFromContourPoint(m, t1, s1, finishAt, finishLevelContour.contours.front().points.front());
    const int ring1 = poleRing(m, t1, finishAt, r.w, r.h);
    const int nFinish = static_cast<int>(std::lround(ringLength(m, t1, finishAt) / r.w));
    {
        double visible = 0, hidden = 0;
        areaOutside(f, t1, -1, finishAt, visible, hidden);
        const SurfaceMarking finish = offered(f, 2, MarkingMeaning::Finish, MarkingShape::ClosedCurve, finishCurve, 0, "finish curve 4 rows before the far end", "the field's level 4 rows before the far end");
        ProposalSetup s = poleSetup(f, s1, ring1, finish, p, "row starts (any side: nothing hides this piece)");
        s.recipe = "pole/curve";
        s.title = "Start at one end, finish 4 rows before the other end";
        s.purpose = "Knit most of the piece from one end; a cap at the other end stays unknitted.";
        s.reason = "Nothing hides either end of this piece, and a closed finish is not available, so the finish is a curve 4 rows before the far end.";
        s.limitation = text("A start point uses the gathered pole start, an experimental construction. The cap beyond the finish, %.0f mm2, is visible and stays unknitted.", visible);
        s.rowsAbout = R;
        s.firstRingStitches = ring1;
        s.lastRingStitches = nFinish;
        s.omittedVisibleMm2 = visible;
        s.omittedHiddenMm2 = hidden;
        s.partial = true;
        out.setups.push_back(s);
    }
    if (startCurve.size() == static_cast<size_t>(kCurvePoints) && R > 6) {
        double visible = 0, hidden = 0;
        areaOutside(f, t1, 4 * hp, finishAt, visible, hidden);
        ProposalSetup s;
        s.recipe = "band";
        s.markings = {offered(f, 1, MarkingMeaning::Start, MarkingShape::ClosedCurve, startCurve, 0, "start curve 4 rows from one end", "the field's level 4 rows from the end"),
                      offered(f, 2, MarkingMeaning::Finish, MarkingShape::ClosedCurve, finishCurve, 0, "finish curve 4 rows before the other end", "the field's level 4 rows before the far end"),
                      offered(f, 3, MarkingMeaning::RowStartPath, MarkingShape::Path, p, 0, "row starts (any side)", "steepest descent of the field")};
        s.title = "Cast on 4 rows from one end, finish 4 rows before the other";
        s.purpose = "Knit the middle of the piece as a band, without the experimental pole start.";
        s.reason = "Cast-on and bind-off curves are supported constructions; both caps beyond them are left out on purpose.";
        s.limitation = text("Both caps stay unknitted: %.0f mm2 of visible surface (measured on the piece).", visible);
        s.rowsAbout = R - 4;
        s.firstRingStitches = static_cast<int>(std::lround(ringLength(m, t1, 4 * hp) / r.w));
        s.lastRingStitches = nFinish;
        s.omittedVisibleMm2 = visible;
        s.omittedHiddenMm2 = hidden;
        s.partial = true;
        out.setups.push_back(s);
    }
    addFinding(out, FindingLevel::Warning, FindingCause::AlgorithmLimitation,
               "Nothing hides either end of this piece, so every suggestion leaves a visible cap unknitted (its area is listed). A closed finish is not "
               "available; draw your own finish curve or boundary to choose what stays unknitted.");
}

// -------------------------------------------------------------------------------------------------------------------
// branching: rows split

void proposeBranching(const ProposalRequest& r, const PieceField& f, const std::vector<std::vector<int>>& loops, ProposalSet& out) {
    const TriMesh& m = f.mesh;
    out.shape = PieceShape::Branching;
    out.shapeText = text("A branching piece with %zu openings", loops.size());
    int widest = 0;
    for (int i = 1; i < static_cast<int>(loops.size()); ++i)
        if (loopLength(f.pieceMm, loops[i]) > loopLength(f.pieceMm, loops[widest])) widest = i;
    const std::vector<int> seedLoop = refinedLoopThrough(m, loops[widest][0]);
    std::vector<double> t;
    TimeFieldReport tr;
    std::string e;
    if (seedLoop.empty() || !geodesicDistance(m, seedLoop, t, tr, &e)) {
        addFinding(out, FindingLevel::Problem, FindingCause::Input, "No field from the widest opening: " + e);
        return;
    }
    const double tMax = *std::max_element(t.begin(), t.end());
    int split = -1;
    RingLevel at;
    for (int k = 1; k * r.h < tMax; ++k) {
        at = extractLevel(m, t, k * r.h);
        dropContoursShorterThan(at, kKnitOverStitches * r.w);   // a bump too small for a round is knitted over, as the layout does
        if (at.contours.size() != 1 || !at.contours.front().closed) {
            split = k;
            break;
        }
    }
    if (split < 0) {
        addFinding(out, FindingLevel::Problem, FindingCause::AlgorithmLimitation, "This piece has more openings than a tube; Method A needs one tube between a start and a finish.");
        return;
    }
    std::vector<Vec3> where;
    for (const Contour& c : at.contours)
        for (const Vec3& p : contourSample(c, 8)) where.push_back(p);
    addFinding(out, FindingLevel::Problem, FindingCause::AlgorithmLimitation,
               text("Rows split here (row %d from the widest opening falls apart into %zu pieces): this piece branches. Method A needs one tube between a start and a finish, and "
                    "splitting a piece is not available, so no complete setup is offered.",
                    split, at.contours.size()),
               split, fileUnits(where, r.mmPerUnit));
    const double finishAt = (split - 1) * r.h;
    if (split < 3) return;
    const RingLevel lv = extractLevel(m, t, finishAt);
    const std::vector<SurfacePoint> finishCurve = levelCurve(m, t, finishAt, kCurvePoints);
    if (lv.contours.size() != 1 || finishCurve.size() != static_cast<size_t>(kCurvePoints)) return;
    const std::set<int> seeds(seedLoop.begin(), seedLoop.end());
    const RingPoint& rp = lv.contours.front().points.front();
    std::vector<SurfacePoint> path;
    chainControlPoints(m, descentChain(m, t, t[rp.a] >= finishAt ? rp.a : rp.b, seeds), path);
    double visible = 0, hidden = 0;
    areaOutside(f, t, -1, finishAt, visible, hidden);
    ProposalSetup s;
    s.recipe = "branch/partial";
    s.markings = {offeredOnFile(1, MarkingMeaning::Start, MarkingShape::ClosedCurve, borderLoopPoints(f, loops[widest]), "cast-on: the widest opening", "a border loop of the piece"),
                  offered(f, 2, MarkingMeaning::Finish, MarkingShape::ClosedCurve, finishCurve, 0, "finish curve one row before the split", "the last level before the rows split"),
                  offered(f, 3, MarkingMeaning::RowStartPath, MarkingShape::Path, path, 0, "row starts", "steepest descent of the field")};
    s.title = "Partial: knit from the widest opening up to where the piece branches";
    s.purpose = "Knit the part before the branches; the branches stay unknitted.";
    s.reason = text("The rows split at row %d; the finish is one row before that.", split);
    s.limitation = text("Partial: %.0f mm2 of visible surface (the branches) stays unknitted. Splitting the piece into tubes is not available.", visible);
    s.rowsAbout = split - 1;
    s.firstRingStitches = static_cast<int>(std::lround(loopLength(f.pieceMm, loops[widest]) / r.w));
    s.lastRingStitches = static_cast<int>(std::lround(lv.contours.front().length / r.w));
    s.omittedVisibleMm2 = visible;
    s.omittedHiddenMm2 = hidden;
    s.partial = true;
    out.setups.push_back(s);
}

bool usable(const ProposalRequest& r, std::string* err) {
    if (!r.pieces || r.piece < 0 || r.piece >= static_cast<int>(r.pieces->pieces.size()) || r.mmPerUnit <= 0 || r.w <= 0 || r.h <= 0) {
        if (err) *err = "proposals: no such piece, or a unit or gauge that is not positive";
        return false;
    }
    return true;
}

// Numbers the setups 1, 2, ... and stamps their markings with the setup and the basis.
void stamp(ProposalSet& out) {
    for (size_t i = 0; i < out.setups.size(); ++i) {
        ProposalSetup& s = out.setups[i];
        s.setup = static_cast<int>(i) + 1;
        for (SurfaceMarking& m : s.markings) {
            m.setup = s.setup;
            m.basis = out.basis;
        }
    }
}

} // namespace

bool proposeMarkings(const ProposalRequest& r, ProposalSet& out, std::string* err) {
    out = ProposalSet{};
    if (!usable(r, err)) return false;
    const std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
    out.piece = r.piece;
    out.basis = proposalBasis(r.fileDigest, r.piece, toString(r.preparation), r.mmPerUnit, r.w, r.h);
    PieceField f;
    preparePieceField(r, f);
    advance(r);
    const std::vector<std::vector<int>> loops = borderLoops(f.pieceMm);
    const int chi = eulerCharacteristic(f.pieceMm);
    out.log.push_back(text("piece %d: %zu triangles, refined x%d -> %zu (mm, %.4g mm per unit); %zu border loops, Euler characteristic %d; %zu closed neighbouring pieces",
                           r.piece, f.pieceMm.tri.size(), f.refined.levels, f.mesh.tri.size(), r.mmPerUnit, loops.size(), chi, f.otherMeshes.size()));
    advance(r);
    if (loops.size() == 2 && chi == 0) {
        proposeOpenTube(r, f, loops, out);
    } else if (loops.size() == 1 && chi == 1) {
        proposeCup(r, f, loops[0], out);
    } else if (loops.empty() && chi == 2) {
        proposeClosed(r, f, out);
    } else if (loops.size() >= 3 || (!loops.empty() && chi < 0)) {
        proposeBranching(r, f, loops, out);
    } else {
        out.shape = PieceShape::Unsupported;
        out.shapeText = text("A piece with %zu openings and Euler characteristic %d", loops.size(), chi);
        addFinding(out, FindingLevel::Problem, FindingCause::AlgorithmLimitation,
                   "Suggestions need an open tube, a cup, a closed piece or a branching piece with openings; this piece's shape has no suggestion (a handle "
                   "or a hole in a closed piece is not available).");
    }
    advance(r);
    stamp(out);
    advance(r);
    out.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return true;
}

bool refineProposalStart(const ProposalRequest& r, const ProposalSetup& from, ProposalSetup& out, std::string* err) {
    if (!usable(r, err)) return false;
    if (!from.refinable || from.markings.size() != 3 || from.markings[0].meaning != MarkingMeaning::Start || from.markings[0].shape != MarkingShape::Point) {
        if (err) *err = "only a setup with a start point can be refined";
        return false;
    }
    PieceField f;
    preparePieceField(r, f);
    const TriMesh& m = f.mesh;
    std::vector<SurfaceMarking> three = from.markings;
    for (size_t i = 0; i < three.size(); ++i) {
        three[i].id = static_cast<int>(i) + 1;
        three[i].proposal = ProposalState::Offered;
    }
    Vec3 settled;
    if (!trialPole(r, f, three, settled, err)) return false;
    advance(r);
    // the start where the pole settled; the finish and path again by the setup's recipe
    int pieceIndex = -1;
    SurfacePoint local;
    pieceLocalPoint(*r.pieces, from.markings[0].points[0], pieceIndex, local);
    const Vec3 before = surfacePosition(m, toRefined(f.refined, local));
    const SurfacePoint s2 = closestSurfacePoint(m, settled);
    std::vector<double> t2;
    field(m, s2, t2);
    const double tMax = *std::max_element(t2.begin(), t2.end());
    const bool hiddenRecipe = from.recipe.rfind("pole/inner", 0) == 0 || from.recipe.rfind("pole/outer", 0) == 0;
    double finishAt = -1, fullyHidden = 0;
    if (hiddenRecipe) {
        finishAt = finishLevel(f, t2, r.h, fullyHidden);
        if (from.recipe == "pole/inner/wide" && finishAt > 0)
            for (double Lv = finishAt; Lv > 2 * r.h; Lv -= 0.25)
                if (static_cast<int>(std::lround(ringLength(m, t2, Lv) / r.w)) >= kWideFinishStitches) {
                    finishAt = Lv;
                    break;
                }
    } else {
        finishAt = tMax - 4 * r.h;
    }
    if (finishAt <= 2 * r.h) {
        if (err) *err = "no finish level in the refined start's field";
        return false;
    }
    std::vector<SurfacePoint> p;
    SurfaceMarking finish;
    if (hiddenRecipe) {
        p = descentPath(f, t2, s2, finishAt, from.recipe != "pole/outer");
        finish = offered(f, 2, MarkingMeaning::Finish, MarkingShape::Point, {onPathAtLevel(m, t2, p, finishAt)}, 0, from.markings[1].note, from.markings[1].reason);
    } else {
        const RingLevel lv = extractLevel(m, t2, finishAt);
        const std::vector<SurfacePoint> curve = levelCurve(m, t2, finishAt, kCurvePoints);
        if (lv.contours.size() != 1 || curve.size() != static_cast<size_t>(kCurvePoints)) {
            if (err) *err = "the refined start's finish level is not one closed contour";
            return false;
        }
        p = descentFromContourPoint(m, t2, s2, finishAt, lv.contours.front().points.front());
        finish = offered(f, 2, MarkingMeaning::Finish, MarkingShape::ClosedCurve, curve, 0, from.markings[1].note, from.markings[1].reason);
    }
    const int ring1 = poleRing(m, t2, finishAt, r.w, r.h);
    out = poleSetup(f, s2, ring1, finish, p, from.markings[2].note);
    const double moved = norm(settled - before);
    out.recipe = from.recipe + "/refined";
    out.title = text("Start where a trial relaxation settled the pole (%.1f mm from the suggested start)", moved);
    out.purpose = from.purpose;
    out.reason = text("A relaxation of the suggested layout with nothing holding the pole moved it %.1f mm; starting there avoids that drift. This is a deeper "
                      "calculation (one layout build with relaxation), not a verdict on the first suggestion.",
                      moved);
    double visible = 0, hidden = 0;
    areaOutside(f, t2, -1, finishAt, visible, hidden);
    out.rowsAbout = std::max(1, static_cast<int>(std::lround(finishAt / r.h)));
    out.firstRingStitches = ring1;
    out.lastRingStitches = static_cast<int>(std::lround(ringLength(m, t2, finishAt) / r.w));
    out.limitation = "A start point uses the gathered pole start, an experimental construction. " + narrowFinishText(out.lastRingStitches) +
                     text("Left out: %.0f mm2 hidden, %.0f mm2 visible (measured on the piece).", hidden, visible);
    out.omittedVisibleMm2 = visible;
    out.omittedHiddenMm2 = hidden;
    out.partial = visible > 1;
    const uint64_t basis = proposalBasis(r.fileDigest, r.piece, toString(r.preparation), r.mmPerUnit, r.w, r.h);
    for (SurfaceMarking& mk : out.markings) mk.basis = basis;
    advance(r);
    return true;
}

bool labProposalFile(const ProposalRequest& r, MarkingSet& set, std::string& report, std::string* err) {
    if (!usable(r, err)) return false;
    PieceField f;
    preparePieceField(r, f);
    const TriMesh& m = f.mesh;
    const SurfacePiece& piece = *f.piece;
    report = text("piece %d (%s, piece %d of its component), %zu triangles, refined x%d -> %zu (mm, %.4g mm per unit)\n", r.piece, r.componentName.c_str(), piece.ordinal, piece.mesh.tri.size(),
                  f.refined.levels, m.tri.size(), r.mmPerUnit);
    // rule 1: the visible end of the geodesic diameter
    int startVertex = -1;
    if (!ruleOneStart(f, startVertex, err)) return false;
    const SurfacePoint s1 = interior(m, vertexPt(m, startVertex).tri, m.pos[startVertex]);
    std::vector<double> t1;
    field(m, s1, t1);
    double hidden1 = 0;
    const double f1 = finishLevel(f, t1, r.h, hidden1);
    if (f1 < 0) {
        if (err) *err = "no level of the piece lies entirely inside another piece: no finish proposed";
        return false;
    }
    const std::vector<SurfacePoint> p1 = descentPath(f, t1, s1, f1, true);
    const SurfacePoint fin1 = onPathAtLevel(m, t1, p1, f1);
    report += text("rule 1 (start at the visible end of the geodesic diameter): start (%.2f, %.2f, %.2f) mm, first entirely hidden level %.2f, finish level %.2f mm\n",
                   surfacePosition(m, s1).x, surfacePosition(m, s1).y, surfacePosition(m, s1).z, hidden1, f1);
    // rule 2: the start where an unconstrained relaxation of rule 1's layout settles the pole (one build of the backend)
    Vec3 settled;
    std::vector<SurfaceMarking> three = {offered(f, 1, MarkingMeaning::Start, MarkingShape::Point, {s1}, 0, "", ""),
                                         offered(f, 2, MarkingMeaning::Finish, MarkingShape::Point, {fin1}, 0, "", ""),
                                         offered(f, 3, MarkingMeaning::RowStartPath, MarkingShape::Path, p1, 0, "", "")};
    if (!trialPole(r, f, three, settled, err)) return false;
    const SurfacePoint s2 = closestSurfacePoint(m, settled);
    std::vector<double> t2;
    field(m, s2, t2);
    double hidden2 = 0;
    const double f2 = finishLevel(f, t2, r.h, hidden2);
    const std::vector<SurfacePoint> p2 = descentPath(f, t2, s2, f2, true), p2o = descentPath(f, t2, s2, f2, false);
    const SurfacePoint fin2 = onPathAtLevel(m, t2, p2, f2);
    report += text("rule 2 (start where an unconstrained relaxation of rule 1's layout settles the pole, %.2f mm from rule 1's start): finish level %.2f mm\n",
                   norm(settled - surfacePosition(m, s1)), f2);
    int auto1 = 0, n2 = 0;
    const int ring1 = poleRing(m, t2, f2, r.w, r.h, &auto1, &n2);
    report += text("pole ring for rule 2: automatic %d, ring 2 has %d, so %d (the smallest whose first strip's increases fit n1 - 2)\n", auto1, n2, ring1);
    const int R2 = std::max(1, static_cast<int>(std::lround(f2 / r.h)));
    const double hp = f2 / R2;
    // a start curve for the band case: the level 4 rows from rule 2's start (the cap inside it is omitted on purpose)
    const std::vector<SurfacePoint> startCurve = levelCurve(m, t2, 4 * hp, kCurvePoints);
    const std::vector<SurfacePoint> finishCurve = levelCurve(m, t2, f2, kCurvePoints);
    set = MarkingSet{};
    set.source = r.source;
    set.digest = r.fileDigest;
    SurfaceMarking start2 = offered(f, 1, MarkingMeaning::Start, MarkingShape::Point, {s2}, 1, "start at the free end (pole start)",
                                    "rule 2: where an unconstrained relaxation of rule 1's layout settled the pole; a proposal, not a verdict on rule 1");
    start2.settings = {{"ring1", std::to_string(ring1)}};
    const SurfaceMarking start1 = offered(f, 4, MarkingMeaning::Start, MarkingShape::Point, {s1}, 1, "start at the free end's tip (pole start); pairs with path 5 and finish 6",
                                          "rule 1: the visible end of the piece's geodesic diameter");
    const SurfaceMarking startC = offered(f, 9, MarkingMeaning::Start, MarkingShape::ClosedCurve, startCurve, 1,
                                          "start curve 4 rows from the free end (band: the end cap is omitted on purpose); pairs with path 3 and finish curve 8",
                                          "the level 4 h' of rule 2's field");
    set.markings = {start2,
                    offered(f, 2, MarkingMeaning::Finish, MarkingShape::Point, {fin2}, 2, "finish level inside the piece that hides this end",
                            "the first level of rule 2's field entirely inside another piece, plus half a row"),
                    offered(f, 3, MarkingMeaning::RowStartPath, MarkingShape::Path, p2, 3, "row-start path, inner side (facing the hiding piece)",
                            "steepest descent of rule 2's field from the finish level's point facing the hiding piece"),
                    start1,
                    offered(f, 5, MarkingMeaning::RowStartPath, MarkingShape::Path, p1, 3, "row-start path, inner side, for start 4", "as path 3, in rule 1's field"),
                    offered(f, 6, MarkingMeaning::Finish, MarkingShape::Point, {fin1}, 2, "finish level inside the hiding piece, for start 4", "as finish 2, in rule 1's field"),
                    offered(f, 7, MarkingMeaning::RowStartPath, MarkingShape::Path, p2o, 3, "row-start path, outer side (facing away from the hiding piece)", "as path 3, from the point facing away"),
                    offered(f, 8, MarkingMeaning::Finish, MarkingShape::ClosedCurve, finishCurve, 2, "finish curve at the same level as finish 2",
                            "the contour of finish 2's level, 12 control points"),
                    startC};
    set.nextId = 10;
    return true;
}

} // namespace knit
