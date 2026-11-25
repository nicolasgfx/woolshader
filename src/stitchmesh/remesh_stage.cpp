#include "stitchmesh/remesh_stage.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <map>
#include <mutex>
#include <set>
#include <utility>

#include "rim_remesh.h"

namespace knit {

namespace {

constexpr int kAlgorithmVersion = 1;

// The bounds of GeometricQualification (its header documents them).
constexpr double kEdgeLow = 0.85, kEdgeHigh = 1.15;
constexpr double kMinAngleDeg = 45;
constexpr double kMaxSurfaceDistance = 0.10;
constexpr double kMaxFieldAlignmentDeg = 20;
constexpr double kMaxIrregularShare = 0.10;
constexpr double kMaxTriangleFraction = 0.20;
constexpr double kCollapsedArea = 1e-3;     // of coarse edge^2
constexpr double kExtremeEdgeLow = 0.5, kExtremeEdgeHigh = 2.0;  // local defects, / coarse edge
constexpr double kSharpFaceDeg = 20;
constexpr double kUnboundDistance = 0.5;    // coarse edges

// RIM keeps file-scope state: one remesh at a time in the process.
std::mutex g_remesh;

using Clock = std::chrono::steady_clock;
double msSince(Clock::time_point t) { return std::chrono::duration<double, std::milli>(Clock::now() - t).count(); }

void toArrays(const TriMesh& m, std::vector<std::array<double, 3>>& pos, std::vector<std::array<int, 3>>& tri) {
    pos.resize(m.pos.size());
    for (size_t i = 0; i < m.pos.size(); ++i) pos[i] = {m.pos[i].x, m.pos[i].y, m.pos[i].z};
    tri = m.tri;
}

Vec3 toVec(const std::array<double, 3>& a) { return {a[0], a[1], a[2]}; }

struct Variant {
    const char* name;
    int sweepsFactor;
    bool extrinsic;
};

std::string requestText(const RemeshRequest& r, const Variant& v) {
    char b[400];
    std::snprintf(b, sizeof b,
                  "    gauge = %.3f x %.3f mm\n    coarse edge = %.3f mm (factor %.2f)\n    variant = %s (sweeps per level %d, %s field)\n"
                  "    seed = %llu\n    axis guide = %s\n    attempt cap = %.0f s\n",
                  r.w, r.h, r.coarseEdge(), r.coarseFactor, v.name, r.sweepsPerLevel * v.sweepsFactor, v.extrinsic ? "extrinsic" : "intrinsic",
                  static_cast<unsigned long long>(r.seed), r.axisGuide.empty() ? "none (solved)" : "fixed per reference vertex", r.attemptCapS);
    return b;
}

// The interior edges of m as (origin, destination) with both faces present, each undirected edge once.
std::vector<std::pair<int, int>> interiorEdges(const PolygonMesh& m, const HalfEdgeMesh& he) {
    std::vector<std::pair<int, int>> out;
    for (int h = 0; h < static_cast<int>(he.next.size()); ++h) {
        if (he.twin[h] < 0 || he.twin[h] < h) continue;  // boundary, or the pair already listed from its twin
        out.push_back({he.origin[h], he.destination(h)});
    }
    (void) m;
    return out;
}

} // namespace

// The order of attempts: qualified first, then valid, then invalid; within a class the distortion, the surface fidelity,
// the irregular vertices, then the earlier attempt.
bool ranksBefore(const RemeshAttempt& a, int ia, const RemeshAttempt& b, int ib) {
    auto klass = [](const RemeshAttempt& x) {
        if (x.structural.status != Status::Pass) return 2;
        return x.quality.status == Status::Pass ? 0 : 1;
    };
    const int ka = klass(a), kb = klass(b);
    if (ka != kb) return ka < kb;
    if (a.quality.distortionScore != b.quality.distortionScore) return a.quality.distortionScore < b.quality.distortionScore;
    if (a.quality.fidelityScore != b.quality.fidelityScore) return a.quality.fidelityScore < b.quality.fidelityScore;
    if (a.quality.irregularVertices != b.quality.irregularVertices) return a.quality.irregularVertices < b.quality.irregularVertices;
    return ia < ib;
}

int bestAttempt(const std::vector<RemeshAttempt>& attempts) {
    int best = -1;
    for (size_t i = 0; i < attempts.size(); ++i)
        if (best < 0 || ranksBefore(attempts[i], static_cast<int>(i), attempts[best], best)) best = static_cast<int>(i);
    return best;
}

// --- the checks

std::string StructuralCheck::text() const {
    char b[256];
    std::snprintf(b, sizeof b, "structural %s: inverted %d, collapsed %d, unbound %d", toString(status), invertedFaces, collapsedFaces, unboundVertices);
    std::string t = b;
    if (!detail.empty()) t += " (" + detail + ")";
    return t + "\n  " + mesh.text();
}

StructuralCheck checkStructure(const PolygonMesh& m, const ReferenceSurface& s, const std::vector<VertexBinding>& binding, double coarseEdge) {
    StructuralCheck c;
    c.mesh = checkPolygonMesh(m);
    for (int f = 0; f < m.faces(); ++f) {
        if (m.faceArea(f) < kCollapsedArea * coarseEdge * coarseEdge) ++c.collapsedFaces;
        const MeshQuery::Hit hit = s.query->closest(m.faceCentroid(f));
        if (hit.tri >= 0 && dot(m.faceNormal(f), s.faceNormal[hit.tri]) < 0) ++c.invertedFaces;
    }
    for (const VertexBinding& b : binding)
        if (b.tri < 0 || b.distance > kUnboundDistance * coarseEdge) ++c.unboundVertices;
    const bool ok = c.mesh.status == Status::Pass && c.invertedFaces == 0 && c.collapsedFaces == 0 && c.unboundVertices == 0;
    c.status = ok ? Status::Pass : Status::Fail;
    if (!ok) {
        if (c.mesh.status != Status::Pass) c.detail = "mesh: " + c.mesh.detail;
        else if (c.invertedFaces) c.detail = "faces wound against the reference normal";
        else if (c.collapsedFaces) c.detail = "collapsed faces";
        else c.detail = "vertices off the reference surface";
    }
    return c;
}

std::string GeometricQualification::text() const {
    char b[1100];
    std::snprintf(b, sizeof b,
                  "qualification %s%s%s\n"
                  "  edge / coarse edge: min %.3f p5 %.3f p50 %.3f p95 %.3f max %.3f (%d edges between interior vertices; bounds p5 >= %.2f, p95 <= %.2f)\n"
                  "  edges touching a boundary / coarse edge: min %.3f p50 %.3f max %.3f (%d, reported, not bounded)\n"
                  "  min angle [deg]: min %.1f p5 %.1f p50 %.1f (bound p5 >= %.0f)\n"
                  "  surface distance / coarse edge: p50 %.4f p95 %.4f max %.4f (bound p95 <= %.2f)\n"
                  "  field alignment [deg]: p50 %.1f p95 %.1f max %.1f (bound p95 <= %.0f)\n"
                  "  irregular interior vertices %d of %d (bound <= %.0f %%), triangle fraction %.3f (bound <= %.2f)\n"
                  "  distortion %.4f, fidelity %.4f",
                  toString(status), detail.empty() ? "" : ": ", detail.c_str(), edgeLength.min, edgeLength.p5, edgeLength.p50, edgeLength.p95,
                  edgeLength.max, edgeLength.count, kEdgeLow, kEdgeHigh, boundaryEdgeLength.min, boundaryEdgeLength.p50, boundaryEdgeLength.max,
                  boundaryEdgeLength.count, minAngleDeg.min, minAngleDeg.p5, minAngleDeg.p50, kMinAngleDeg,
                  surfaceDistance.p50, surfaceDistance.p95, surfaceDistance.max, kMaxSurfaceDistance, fieldAlignmentDeg.p50,
                  fieldAlignmentDeg.p95, fieldAlignmentDeg.max, kMaxFieldAlignmentDeg, irregularVertices, interiorVertices,
                  kMaxIrregularShare * 100, triangleFraction, kMaxTriangleFraction, distortionScore, fidelityScore);
    std::string t = b;
    double worstEdge = 1, worstAngle = 180;
    for (const EdgeDefect& e : extremeEdges) worstEdge = std::fabs(e.length - 1) > std::fabs(worstEdge - 1) ? e.length : worstEdge;
    for (const FaceDefect& f : sharpFaces) worstAngle = std::min(worstAngle, f.minAngleDeg);
    std::snprintf(b, sizeof b, "\n  local defects: %zu edges outside [%.1f, %.1f] (worst %.2f), %zu faces with a corner under %.0f deg (worst %.1f)", extremeEdges.size(),
                  kExtremeEdgeLow, kExtremeEdgeHigh, worstEdge, sharpFaces.size(), kSharpFaceDeg, worstAngle);
    t += b;
    t += "\n  valence histogram:";
    for (size_t v = 0; v < valenceHistogram.size(); ++v)
        if (valenceHistogram[v]) t += " " + std::to_string(v) + ":" + std::to_string(valenceHistogram[v]);
    return t;
}

GeometricQualification qualifyGeometry(const PolygonMesh& m, const std::vector<Vec3>& axis, const std::vector<VertexBinding>& binding,
                                       double coarseEdgeAlong, double coarseEdgeAcross) {
    const double coarseEdge = 0.5 * (coarseEdgeAlong + coarseEdgeAcross);
    GeometricQualification q;
    HalfEdgeMesh he;
    std::string why;
    if (!buildHalfEdges(m, he, &why)) {
        q.status = Status::Fail;
        q.detail = "no half-edge structure: " + why;
        return q;
    }
    std::vector<double> edges, boundaryEdges, angles, distances, alignments;
    std::vector<int> valence(m.pos.size(), 0);
    for (const std::pair<int, int>& e : interiorEdges(m, he)) {
        const Vec3 d = m.pos[e.second] - m.pos[e.first];
        // the edge's own target: along the axis at its origin (within 45 degrees) or across it; the mean without an axis
        double target = coarseEdge;
        if (e.first < static_cast<int>(axis.size()) && norm(axis[e.first]) > 0 && norm(d) > 0)
            target = std::fabs(dot(normalized(d), normalized(axis[e.first]))) > 0.70710678 ? coarseEdgeAlong : coarseEdgeAcross;
        if (he.boundaryVertex[e.first] || he.boundaryVertex[e.second]) {
            boundaryEdges.push_back(norm(d) / target);
            continue;
        }
        edges.push_back(norm(d) / target);
        if (edges.back() < kExtremeEdgeLow || edges.back() > kExtremeEdgeHigh) q.extremeEdges.push_back({e.first, e.second, edges.back()});
        // the edge against the 2-RoSy axis at its origin: 0 along the axis or across it (either lattice axis), 45 at worst
        if (e.first < static_cast<int>(axis.size()) && norm(axis[e.first]) > 0 && norm(d) > 0) {
            const double c = std::clamp(std::fabs(dot(normalized(d), normalized(axis[e.first]))), 0.0, 1.0);
            const double deg = std::acos(c) * 180 / kPi;  // 0..90
            alignments.push_back(std::min(deg, 90 - deg));
        }
    }
    // the valence of every vertex: its undirected edges, boundary edges included
    std::set<std::pair<int, int>> undirected;
    for (int h = 0; h < static_cast<int>(he.next.size()); ++h) {
        const int a = he.origin[h], b = he.destination(h);
        undirected.insert({std::min(a, b), std::max(a, b)});
    }
    for (const std::pair<int, int>& e : undirected) {
        ++valence[e.first];
        ++valence[e.second];
    }
    for (size_t v = 0; v < m.pos.size(); ++v) {
        if (he.halfEdgeOfVertex[v] < 0 || he.boundaryVertex[v]) continue;
        ++q.interiorVertices;
        if (valence[v] >= static_cast<int>(q.valenceHistogram.size())) q.valenceHistogram.resize(valence[v] + 1, 0);
        ++q.valenceHistogram[valence[v]];
        if (valence[v] != 4) ++q.irregularVertices;
    }
    int triangles = 0;
    for (int f = 0; f < m.faces(); ++f) {
        const int s = m.sides(f);
        if (s == 3) ++triangles;
        double minAngle = 180;
        for (int k = 0; k < s; ++k) {
            const Vec3 a = m.pos[m.vertexAt(f, (k + s - 1) % s)] - m.pos[m.vertexAt(f, k)];
            const Vec3 b = m.pos[m.vertexAt(f, (k + 1) % s)] - m.pos[m.vertexAt(f, k)];
            if (!(norm(a) > 0 && norm(b) > 0)) continue;
            minAngle = std::min(minAngle, std::acos(std::clamp(dot(a, b) / (norm(a) * norm(b)), -1.0, 1.0)) * 180 / kPi);
        }
        angles.push_back(minAngle);
        if (minAngle < kSharpFaceDeg) q.sharpFaces.push_back({f, minAngle});
    }
    for (const VertexBinding& b : binding) distances.push_back(b.tri < 0 ? 1e9 : b.distance / coarseEdge);
    q.edgeLength = quantiles(edges);
    q.boundaryEdgeLength = quantiles(boundaryEdges);
    q.minAngleDeg = quantiles(angles);
    q.surfaceDistance = quantiles(distances);
    q.fieldAlignmentDeg = quantiles(alignments);
    q.triangleFraction = m.faces() ? static_cast<double>(triangles) / m.faces() : 0;
    q.distortionScore = std::fabs(q.edgeLength.p95 - 1) + std::fabs(q.edgeLength.p5 - 1) + (90 - q.minAngleDeg.p5) / 90;
    q.fidelityScore = q.surfaceDistance.p95;
    char b[160];
    if (edges.empty()) q.detail = "no edge between interior vertices";
    else if (q.edgeLength.p95 > kEdgeHigh || q.edgeLength.p5 < kEdgeLow) {
        std::snprintf(b, sizeof b, "edge length p5 %.3f / p95 %.3f outside [%.2f, %.2f]", q.edgeLength.p5, q.edgeLength.p95, kEdgeLow, kEdgeHigh);
        q.detail = b;
    } else if (q.minAngleDeg.p5 < kMinAngleDeg) {
        std::snprintf(b, sizeof b, "min angle p5 %.1f below %.0f degrees", q.minAngleDeg.p5, kMinAngleDeg);
        q.detail = b;
    } else if (q.surfaceDistance.p95 > kMaxSurfaceDistance) {
        std::snprintf(b, sizeof b, "surface distance p95 %.4f above %.2f coarse edges", q.surfaceDistance.p95, kMaxSurfaceDistance);
        q.detail = b;
    } else if (!alignments.empty() && q.fieldAlignmentDeg.p95 > kMaxFieldAlignmentDeg) {
        std::snprintf(b, sizeof b, "field alignment p95 %.1f above %.0f degrees", q.fieldAlignmentDeg.p95, kMaxFieldAlignmentDeg);
        q.detail = b;
    } else if (q.interiorVertices > 0 && q.irregularVertices > kMaxIrregularShare * q.interiorVertices) {
        std::snprintf(b, sizeof b, "%d irregular of %d interior vertices, above %.0f %%", q.irregularVertices, q.interiorVertices, kMaxIrregularShare * 100);
        q.detail = b;
    } else if (q.triangleFraction > kMaxTriangleFraction) {
        std::snprintf(b, sizeof b, "triangle fraction %.3f above %.2f", q.triangleFraction, kMaxTriangleFraction);
        q.detail = b;
    }
    q.status = q.detail.empty() ? Status::Pass : Status::Fail;
    return q;
}

// --- the stage

SizePreflight preflightSize(const ReferenceSurface& s, const RemeshRequest& r) {
    SizePreflight p;
    p.inputTriangles = static_cast<int>(s.mesh->tri.size());
    const double edge = r.coarseEdge();
    p.predictedCoarseFaces = edge > 0 ? static_cast<int>(s.area / (edge * edge)) : 0;
    char b[200];
    if (p.inputTriangles > r.limits.inputTriangles) {
        std::snprintf(b, sizeof b, "preflight: %d input triangles above the limit %d", p.inputTriangles, r.limits.inputTriangles);
        p.why = b;
        return p;
    }
    if (p.predictedCoarseFaces > r.limits.coarseFaces) {
        std::snprintf(b, sizeof b, "preflight: about %d coarse faces (area %.0f mm2 / edge %.2f mm squared) above the limit %d", p.predictedCoarseFaces,
                      s.area, edge, r.limits.coarseFaces);
        p.why = b;
        return p;
    }
    std::vector<std::array<double, 3>> pos;
    std::vector<std::array<int, 3>> tri;
    toArrays(*s.mesh, pos, tri);
    {
        const std::lock_guard<std::mutex> lock(g_remesh);
        p.predictedFinestVertices = knit_rim::predictedVertices(pos, tri, edge);
    }
    if (p.predictedFinestVertices > r.limits.finestVertices) {
        std::snprintf(b, sizeof b, "preflight: %d vertices after RIM's subdivision above the limit %d", p.predictedFinestVertices, r.limits.finestVertices);
        p.why = b;
        return p;
    }
    p.ok = true;
    return p;
}

bool remeshStage(const ReferenceSurface& s, const RemeshRequest& r, RemeshProduct& out, std::string* err) {
    out = RemeshProduct{};
    if (!s.mesh || !s.query) {
        if (err) *err = "remesh stage: no reference surface";
        return false;
    }
    const Clock::time_point start = Clock::now();
    StageReport& rep = out.report;
    rep.stage = "remesh";
    rep.algorithmVersion = kAlgorithmVersion;
    rep.inputDigest = s.digest;
    rep.seed = r.seed;
    rep.gaugeW = r.w;
    rep.gaugeH = r.h;
    if (!(r.w > 0) || !(r.h > 0) || !(r.coarseFactor > 0)) {
        rep.status = Status::Fail;
        rep.detail = "a positive gauge and coarse factor are required";
        rep.elapsedMs = msSince(start);
        return true;
    }
    out.preflight = preflightSize(s, r);
    rep.counts.push_back({"input triangles", out.preflight.inputTriangles});
    rep.counts.push_back({"predicted finest vertices", out.preflight.predictedFinestVertices});
    rep.counts.push_back({"predicted coarse faces", out.preflight.predictedCoarseFaces});
    if (!out.preflight.ok) {
        rep.status = Status::Fail;
        rep.detail = out.preflight.why;
        rep.elapsedMs = msSince(start);
        return true;
    }
    std::vector<std::array<double, 3>> pos;
    std::vector<std::array<int, 3>> tri;
    toArrays(*s.mesh, pos, tri);
    const Variant variants[3] = {{r.extrinsic ? "extrinsic-20" : "intrinsic-20", 1, r.extrinsic},
                                 {r.extrinsic ? "extrinsic-60" : "intrinsic-60", 3, r.extrinsic},
                                 {r.extrinsic ? "intrinsic-20" : "extrinsic-20", 1, !r.extrinsic}};
    std::vector<RemeshAttempt> attempts;
    const double coarseEdge = r.coarseEdge();
    for (int i = 0; i < std::min(r.maxAttempts, 3); ++i) {
        const Variant& v = variants[i];
        const Clock::time_point t0 = Clock::now();
        RemeshAttempt a;
        a.variant = v.name;
        a.report.stage = "remesh attempt";
        a.report.algorithmVersion = kAlgorithmVersion;
        a.report.inputDigest = s.digest;
        a.report.seed = r.seed;
        a.report.gaugeW = r.w;
        a.report.gaugeH = r.h;
        a.report.settings = requestText(r, v);
        knit_rim::Options o;
        o.spacing = r.coarseEdgeAlong();                       // along the field's axis q
        o.aspect = r.coarseEdgeAcross() / r.coarseEdgeAlong();  // across it (t = n x q): RIM's rectangular lattice, inactive at 1
        o.rosy = r.rosy;
        o.extrinsic = v.extrinsic && o.aspect == 1;               // the rectangular lattice needs the intrinsic field
        o.sweepsPerLevel = r.sweepsPerLevel * v.sweepsFactor;
        o.seed = r.seed;
        o.deadlineMs = r.attemptCapS > 0 ? r.attemptCapS * 1000 : 0;
        for (const Vec3& g : r.axisGuide) o.orientation.push_back({g.x, g.y, g.z});
        knit_rim::Result res;
        std::string why;
        bool ok = false;
        {
            const std::lock_guard<std::mutex> lock(g_remesh);
            ok = knit_rim::remesh(pos, tri, o, res, &why);
        }
        a.backendLog = res.log;
        a.report.counts.push_back({"finest vertices", res.inputVertices});
        a.report.counts.push_back({"ms hierarchy", res.msHierarchy});
        a.report.counts.push_back({"ms orientations", res.msOrientations});
        a.report.counts.push_back({"ms positions", res.msPositions});
        a.report.counts.push_back({"ms extraction", res.msExtraction});
        if (!ok) {
            a.structural.status = Status::Fail;
            a.structural.detail = why;
            a.quality.status = Status::Fail;
            a.quality.detail = why;
            a.report.status = Status::Fail;
            a.report.detail = res.deadlineHit ? why : "backend: " + why;
            a.report.elapsedMs = msSince(t0);
            attempts.push_back(std::move(a));
            continue;
        }
        PolygonMesh extracted;
        extracted.pos.reserve(res.pos.size());
        for (const std::array<double, 3>& p : res.pos) extracted.pos.push_back(toVec(p));
        for (const std::vector<int>& f : res.faces) extracted.addFace(f);
        a.normal.reserve(res.nrm.size());
        for (const std::array<double, 3>& n : res.nrm) a.normal.push_back(toVec(n));
        a.field.axis.reserve(res.axis.size());
        for (const std::array<double, 3>& q : res.axis) a.field.axis.push_back(toVec(q));
        for (const std::array<double, 3>& p : res.orientationSingularPoints) a.field.orientationSingularPoints.push_back(toVec(p));
        for (const std::array<double, 3>& p : res.positionSingularPoints) a.field.positionSingularPoints.push_back(toVec(p));
        a.mesh = normalizePolygons(extracted, a.field.axis, a.normal, a.normalization);
        a.binding = bindToSurface(s, a.mesh.pos);
        a.structural = checkStructure(a.mesh, s, a.binding, coarseEdge);
        a.quality = qualifyGeometry(a.mesh, a.field.axis, a.binding, r.coarseEdgeAlong(), r.coarseEdgeAcross());
        a.report.outputDigest = a.mesh.digest();
        a.report.counts.push_back({"faces", a.mesh.faces()});
        a.report.counts.push_back({"triangles", a.structural.mesh.triangles});
        a.report.counts.push_back({"quads", a.structural.mesh.quads});
        a.report.counts.push_back({"vertices", a.structural.mesh.vertices});
        a.report.counts.push_back({"irregular interior vertices", a.quality.irregularVertices});
        a.report.counts.push_back({"orientation singular points", static_cast<double>(a.field.orientationSingularPoints.size())});
        a.report.counts.push_back({"position singular points", static_cast<double>(a.field.positionSingularPoints.size())});
        a.report.counts.push_back({"large polygons split", a.normalization.largeSplits});
        a.report.counts.push_back({"pentagons split", a.normalization.pentagons});
        a.report.counts.push_back({"reflex quads split", a.normalization.reflexQuads});
        a.report.counts.push_back({"wu loop rounds", res.wuRounds});
        a.report.counts.push_back({"wu diagonals", res.wuDiagonals});
        a.report.counts.push_back({"wu loop splits", res.wuLoopSplits});
        if (a.structural.status != Status::Pass) {
            a.report.status = Status::Fail;
            a.report.detail = "structural: " + a.structural.detail;
        } else if (a.quality.status != Status::Pass) {
            a.report.status = Status::Unresolved;
            a.report.detail = "valid but not qualified: " + a.quality.detail;
        } else {
            a.report.status = Status::Pass;
            a.report.detail = "qualified";
        }
        a.report.elapsedMs = msSince(t0);
        const bool qualified = a.report.status == Status::Pass;
        attempts.push_back(std::move(a));
        if (qualified) break;  // no retry after a qualified attempt
    }
    // ranking
    std::vector<int> order(attempts.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = static_cast<int>(i);
    std::stable_sort(order.begin(), order.end(), [&](int x, int y) { return ranksBefore(attempts[x], x, attempts[y], y); });
    for (size_t k = 0; k < order.size(); ++k) {
        if (k == 0) out.kept = attempts[order[k]];
        else out.rejected.push_back(attempts[order[k]]);
    }
    rep.counts.push_back({"attempts", static_cast<double>(attempts.size())});
    rep.outputDigest = out.kept.mesh.digest();
    rep.settings = out.kept.report.settings;
    for (const std::pair<std::string, double>& c : out.kept.report.counts) rep.counts.push_back(c);
    if (attempts.empty()) {
        rep.status = Status::Fail;
        rep.detail = "no attempt ran";
    } else if (out.kept.structural.status != Status::Pass) {
        rep.status = Status::Fail;
        rep.detail = "no valid attempt; the best, " + out.kept.variant + ": " + out.kept.report.detail;
    } else if (out.kept.quality.status != Status::Pass) {
        rep.status = Status::Unresolved;
        rep.detail = "kept " + out.kept.variant + ", valid but not qualified: " + out.kept.quality.detail + "; local defects: " +
                     std::to_string(out.kept.quality.extremeEdges.size()) + " extreme edges, " + std::to_string(out.kept.quality.sharpFaces.size()) + " sharp faces";
    } else {
        rep.status = Status::Pass;
        rep.detail = "kept " + out.kept.variant + ", qualified";
    }
    rep.elapsedMs = msSince(start);
    return true;
}

std::vector<Vec3> globalAxisGuide(const std::string& name, const TriMesh& mesh) {
    if (name == "boundary") return boundaryAxisGuide(mesh);
    const bool ring = name.rfind("ring-", 0) == 0, along = name.rfind("along-", 0) == 0;
    if ((!ring && !along) || name.empty()) return {};
    const char axis = name.back();
    if (axis != 'x' && axis != 'y' && axis != 'z') return {};
    const Vec3 g = axis == 'x' ? Vec3{1, 0, 0} : axis == 'y' ? Vec3{0, 1, 0} : Vec3{0, 0, 1};
    std::vector<Vec3> guide(mesh.pos.size());
    for (size_t v = 0; v < mesh.pos.size(); ++v) {
        const Vec3 n = v < mesh.vnormal.size() ? mesh.vnormal[v] : Vec3{0, 0, 1};
        Vec3 d = ring ? cross(g, n) : g - n * dot(n, g);
        if (norm(d) < 1e-6) d = cross(n, Vec3{0.3, 0.5, 0.8});  // where the axis meets the normal: any tangent
        guide[v] = normalized(d);
    }
    return guide;
}

std::vector<Vec3> boundaryAxisGuide(const TriMesh& mesh) {
    // the boundary edges: used by one triangle
    std::map<std::pair<int, int>, int> uses;
    for (const std::array<int, 3>& t : mesh.tri)
        for (int k = 0; k < 3; ++k) {
            const int a = t[k], b = t[(k + 1) % 3];
            ++uses[{std::min(a, b), std::max(a, b)}];
        }
    std::vector<std::pair<int, int>> boundary;
    std::map<int, std::vector<int>> nextOf;   // the boundary's vertex adjacency
    for (const auto& [e, n] : uses)
        if (n == 1) {
            boundary.push_back(e);
            nextOf[e.first].push_back(e.second);
            nextOf[e.second].push_back(e.first);
        }
    if (boundary.empty()) return {};
    // the loops, by walking; per loop its 2-RoSy mean direction and its axis
    struct Loop {
        Vec3 direction{0, 0, 0};   // the dominant edge direction (unit), or zero when none
        Vec3 axis{0, 0, 1};        // the loop's plane normal (the ring's axis), from its area vector
        std::vector<int> vertices;
    };
    std::vector<Loop> loops;
    std::set<int> seen;
    for (const auto& [v0, nb] : nextOf) {
        if (seen.count(v0)) continue;
        Loop loop;
        int prev = -1, v = v0;
        while (v >= 0 && !seen.count(v)) {
            seen.insert(v);
            loop.vertices.push_back(v);
            const std::vector<int>& around = nextOf[v];
            int next = -1;
            for (int w : around)
                if (w != prev && !seen.count(w)) { next = w; break; }
            prev = v;
            v = next;
        }
        // the 2-RoSy mean: the edges' doubled angles in the loop's plane; the plane from the area vector
        Vec3 area{0, 0, 0}, centre{0, 0, 0};
        for (int q : loop.vertices) centre += mesh.pos[q];
        centre = centre * (1.0 / loop.vertices.size());
        for (size_t i = 0; i < loop.vertices.size(); ++i) {
            const Vec3& a = mesh.pos[loop.vertices[i]];
            const Vec3& b = mesh.pos[loop.vertices[(i + 1) % loop.vertices.size()]];
            area += cross(a - centre, b - centre);
        }
        loop.axis = norm(area) > 1e-12 ? normalized(area) : Vec3{0, 0, 1};
        const Vec3 u = normalized(std::fabs(loop.axis.x) < 0.9 ? cross(loop.axis, Vec3{1, 0, 0}) : cross(loop.axis, Vec3{0, 1, 0}));
        const Vec3 w = cross(loop.axis, u);
        double cx = 0, sx = 0, total = 0;
        for (size_t i = 0; i < loop.vertices.size(); ++i) {
            const Vec3 d = mesh.pos[loop.vertices[(i + 1) % loop.vertices.size()]] - mesh.pos[loop.vertices[i]];
            const double l = norm(d);
            if (l <= 1e-12) continue;
            const double ang = std::atan2(dot(d, w), dot(d, u));
            cx += l * std::cos(2 * ang);
            sx += l * std::sin(2 * ang);
            total += l;
        }
        if (total > 0 && std::hypot(cx, sx) > 0.5 * total) {
            const double ang = 0.5 * std::atan2(sx, cx);
            loop.direction = u * std::cos(ang) + w * std::sin(ang);
        }
        loops.push_back(loop);
    }
    // every vertex: its nearest loop's direction (or the ring round that loop's axis) in its tangent plane
    std::vector<Vec3> guide(mesh.pos.size());
    for (size_t v = 0; v < mesh.pos.size(); ++v) {
        const Vec3 n = v < mesh.vnormal.size() ? mesh.vnormal[v] : Vec3{0, 0, 1};
        double best = 1e300;
        const Loop* nearest = &loops[0];
        for (const Loop& loop : loops)
            for (int q : loop.vertices) {
                const double d = norm(mesh.pos[q] - mesh.pos[v]);
                if (d < best) { best = d; nearest = &loop; }
            }
        Vec3 d = norm(nearest->direction) > 0 ? nearest->direction : cross(nearest->axis, n);
        d = d - n * dot(d, n);
        if (norm(d) < 1e-6) d = cross(n, Vec3{0.3, 0.5, 0.8});
        guide[v] = normalized(d);
    }
    return guide;
}

} // namespace knit
