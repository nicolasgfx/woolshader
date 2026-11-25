// Inspections of finished fibres: the exported tubes, the reference and attachment checks, the
// measured angles and the exact text snapshots.
#include "fibres/fibre_checks.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>

#include "knitting/checks/checks.h"
#include "knitting/tube/param_text.h"
#include "model/cancel.h"
#include "model/hash.h"
#include "model/log.h"
#include "model/spatial_hash.h"

namespace knit {

namespace {

constexpr double kPolylineSagitta = 1e-3;  // [mm] the centreline polyline's bound
constexpr size_t kListedFibres = 100;      // fibre index lists keep this many

using Clock = std::chrono::steady_clock;
Vec3 minV(const Vec3& a, const Vec3& b) { return {std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z)}; }
Vec3 maxV(const Vec3& a, const Vec3& b) { return {std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)}; }

double pointSegment(const Vec3& p, const Vec3& a, const Vec3& b) {
    const Vec3 ab = b - a;
    const double l2 = dot(ab, ab);
    const double t = l2 > 0 ? std::clamp(dot(p - a, ab) / l2, 0.0, 1.0) : 0.0;
    return norm(p - (a + ab * t));
}

// The fibres' chunks in a parallel check, each with its own result.
int64_t fibreChunks(const FibreSet& fs, const ExecPolicy& exec) {
    const int64_t chunk = exec.chunk > 0 ? exec.chunk : 256;
    return (int64_t(fs.fibres.size()) + chunk - 1) / chunk;
}

void appendListed(std::vector<int>& to, const std::vector<int>& from) {
    for (int i : from)
        if (to.size() < kListedFibres) to.push_back(i);
}

} // namespace

// --- the centreline polyline

struct CentrePolyline {
    std::vector<Vec3> a, b;      // piece end points
    std::vector<double> u0, u1;  // curve parameters of the ends
    std::vector<int> yarn;
    std::vector<double> kappa;   // certified curvature bound of the piece's span and its neighbours
    SpatialHash3 hash{1.0};
    double eChord = 0;
};

namespace {

// Per span (in parallel): its pieces; then appended in span order, so the polyline (and its hash
// insertion order) is that of one pass over the spans.
CentrePolyline buildPolyline(const CurveSet& cs, const PatternParams& pp, const ExecPolicy& exec) {
    CentrePolyline pl;
    const Allowances al = Allowances::forExtent(cs.extent());
    struct SpanPieces {
        std::vector<Vec3> pts;  // m + 1 points, at the curve parameters k + q / m
        double kappaN = 0, eChord = 0;
    };
    for (size_t y = 0; y < cs.yarns.size(); ++y) {
        const YarnPath& yp = cs.yarns[y];
        const CurvatureResult cr = checkCurvature(yp, pp.radius(), al, exec);
        const int spans = yp.curve.numSpans();
        auto spanKappa = [&](int k) {
            return k >= 0 && k < static_cast<int>(cr.spanKappaHi.size()) ? std::min(cr.spanKappaHi[k], 1e6) : (k < 0 || k >= spans ? 0.0 : 1e6);
        };
        std::vector<SpanPieces> per(static_cast<size_t>(std::max(spans, 0)));
        parallelFor(spans, exec.chunk > 0 ? exec.chunk : 256, [&](int64_t kBegin, int64_t kEnd, int64_t) {
            cancelCheckpoint();
            for (int k = static_cast<int>(kBegin); k < static_cast<int>(kEnd); ++k) {
                SpanPieces& sp = per[k];
                const double hi = yp.lengthBounds(k, k + 1).hi;
                const double kappa = spanKappa(k);
                sp.kappaN = std::max({spanKappa(k - 1), kappa, spanKappa(k + 1)});
                const double cMax = kappa > 0 ? std::sqrt(8 * kPolylineSagitta / kappa) : hi;
                const int m = std::max(1, static_cast<int>(std::ceil(hi / std::max(cMax, 1e-6))));
                const double c = hi / m;
                sp.eChord = kappa * c * c / 8;
                for (int q = 0; q <= m; ++q) sp.pts.push_back(yp.curve.eval(k + double(q) / m));
            }
        }, exec);
        for (int k = 0; k < spans; ++k) {
            const SpanPieces& sp = per[k];
            pl.eChord = std::max(pl.eChord, sp.eChord);
            const int m = static_cast<int>(sp.pts.size()) - 1;
            for (int q = 1; q <= m; ++q) {
                const Vec3 &prev = sp.pts[q - 1], &p = sp.pts[q];
                pl.a.push_back(prev);
                pl.b.push_back(p);
                pl.u0.push_back(k + double(q - 1) / m);
                pl.u1.push_back(k + double(q) / m);
                pl.yarn.push_back(static_cast<int>(y));
                pl.kappa.push_back(sp.kappaN);
                pl.hash.insert(minV(prev, p), maxV(prev, p));
            }
        }
    }
    pl.hash.finalize();
    return pl;
}

// The part of the polyline a distance is measured to: everything except yarn exYarn on the
// parameter range [exU0, exU1] (pieces clipped exactly in their own linear parameter).
struct ExcludedArc {
    int yarn = -1;
    double u0 = 0, u1 = 0;
};

// Distance from p to piece i outside the excluded arc: the piece is clipped exactly in its own linear
// parameter; 1e300 when the excluded arc covers all of it.
double distanceToPiece(const CentrePolyline& pl, int i, const Vec3& p, const ExcludedArc& ex) {
    if (pl.yarn[i] != ex.yarn || pl.u1[i] <= ex.u0 || pl.u0[i] >= ex.u1) return pointSegment(p, pl.a[i], pl.b[i]);
    const double du = pl.u1[i] - pl.u0[i];
    const auto at = [&](double u) { return pl.a[i] + (pl.b[i] - pl.a[i]) * ((u - pl.u0[i]) / du); };
    double d = 1e300;
    if (pl.u0[i] < ex.u0) d = std::min(d, pointSegment(p, pl.a[i], at(ex.u0)));
    if (pl.u1[i] > ex.u1) d = std::min(d, pointSegment(p, at(ex.u1), pl.b[i]));
    return d;
}

// Distance from p to the polyline, without the excluded arc.
double polyDistance(const CentrePolyline& pl, const Vec3& p, double reach, const ExcludedArc& ex) {
    double best = 1e300;
    const Vec3 w{reach, reach, reach};
    pl.hash.query(p - w, p + w, [&](int i) { best = std::min(best, distanceToPiece(pl, i, p, ex)); });
    return best;
}

} // namespace

// --- the exported tubes

struct ExportTubeGeometry {
    const CurveSet* cs = nullptr;
    ExportTubeSpec spec;
    double r = 0.5;
    double scale = 1, rhoMax = 0;
    bool sharedSurface = false;
    std::vector<double> yarnLength;
    std::vector<uint32_t> yarnFirstVertex;
    CentrePolyline pl;
    // One slab: the two end planes (ring planes shared with the neighbours, oriented along the axis)
    // and the 2 x sides facet planes (outward normals), all in export units.
    struct Slab {
        Vec3 c0, n0, c1, n1;
        std::vector<std::pair<Vec3, Vec3>> facets;  // (point, outward unit normal)
        double tol = 0;
    };
    std::vector<std::vector<Slab>> slabs;     // [yarn][j]
    std::vector<std::pair<int, int>> slabId;  // hash entry -> (yarn, j)
    SpatialHash3 slabHash{1.0};               // slab bounding boxes (export units)
    // Clearance search: pieces with curvature above the fold limit 1 / r are "steep"; every other
    // piece's bulge is <= bulgeMax = chordFactor / r.
    double chordFactor = 0, bulgeMax = 0;
    std::vector<int> steep;

    int ringsOf(int y) const { return static_cast<int>((yarnFirstVertex[y + 1] - yarnFirstVertex[y]) / (spec.sides + 1)); }
    double slabDepth(int y, int j, const Vec3& Q) const;  // export units
    // A lower bound of the distance from Q to slab (y, j) when Q is outside it: the largest signed
    // outward plane distance of the convex slab (<= 0 if inside). Export units.
    double slabDistanceLower(int y, int j, const Vec3& Q) const;
    // Clearance restricted to the pieces that can bring it below `need`, so the decision c < need is
    // the full search's: a piece's term is below need only if its distance is below
    // need + r + e_c + its bulge. Pieces not `steep` lie within need + r + e_c + bulgeMax of p; the
    // steep ones are tested on their own, under the full search's box test.
    double clearance(const Vec3& p, double reach, double need, const ExcludedArc& ex) const;
};

namespace {

// Slab (j, j + 1) of a tube, from its exported vertices V(ring, i).
ExportTubeGeometry::Slab makeSlab(const std::vector<Vec3>& E, uint32_t base, int sides, int j) {
    ExportTubeGeometry::Slab sl;
    const int stride = sides + 1;
    const auto V = [&](int ring, int i) { return E[base + ring * stride + i]; };
    Vec3 c0{0, 0, 0}, c1{0, 0, 0};
    for (int i = 0; i < sides; ++i) {
        c0 += V(j, i);
        c1 += V(j + 1, i);
    }
    c0 = c0 * (1.0 / sides);
    c1 = c1 * (1.0 / sides);
    const Vec3 axis = normalized(c1 - c0), mid = (c0 + c1) * 0.5;
    // End planes = the rings' own planes (shared by adjacent slabs, so consecutive slabs tile the
    // tube without wedge gaps at bends), oriented along the axis.
    const auto ringNormal = [&](int ring) {
        Vec3 n{0, 0, 0};
        for (int i = 0; i < sides / 2; ++i)
            n += cross(V(ring, i) - V(ring, i + sides / 2), V(ring, (i + sides / 4) % sides) - V(ring, (i + sides / 4 + sides / 2) % sides));
        n = normalized(n);
        return dot(n, axis) < 0 ? -n : n;
    };
    // Between the two ring planes (a point ON the shared plane belongs to both slabs, so the end
    // planes carry no margin); the depth is the lateral one.
    sl.c0 = c0;
    sl.c1 = c1;
    sl.n0 = ringNormal(j);
    sl.n1 = ringNormal(j + 1);
    sl.tol = 1e-9 * (1 + norm(c0));
    for (int i = 0; i < sides; ++i) {
        const Vec3 A = V(j, i), B = V(j, i + 1), C = V(j + 1, i), D = V(j + 1, i + 1);
        for (const std::array<Vec3, 3>& tri : {std::array<Vec3, 3>{A, C, B}, std::array<Vec3, 3>{B, C, D}}) {
            Vec3 n = normalized(cross(tri[1] - tri[0], tri[2] - tri[0]));
            if (dot(n, (tri[0] + tri[1] + tri[2]) * (1.0 / 3) - mid) < 0) n = -n;
            sl.facets.push_back({tri[0], n});
        }
    }
    return sl;
}

} // namespace

double ExportTubeGeometry::slabDepth(int y, int j, const Vec3& Q) const {
    const Slab& sl = slabs[y][j];
    if (std::min(dot(Q - sl.c0, sl.n0), dot(sl.c1 - Q, sl.n1)) < -sl.tol) return -1e300;
    double depth = 1e300;
    for (const auto& [p, n] : sl.facets) depth = std::min(depth, -dot(Q - p, n));
    return depth;
}

double ExportTubeGeometry::slabDistanceLower(int y, int j, const Vec3& Q) const {
    const Slab& sl = slabs[y][j];
    double d = std::max(-dot(Q - sl.c0, sl.n0), dot(Q - sl.c1, sl.n1));
    for (const auto& [p, n] : sl.facets) d = std::max(d, dot(Q - p, n));
    return d;
}

double ExportTubeGeometry::clearance(const Vec3& p, double reach, double need, const ExcludedArc& ex) const {
    double best = 1e300;
    const double nearReach = std::min(reach, need + r + pl.eChord + bulgeMax + 1e-6);
    const Vec3 w{nearReach, nearReach, nearReach};
    // A piece's term: its distance minus the bulge chordFactor kappa_i (the outward chord sagitta of the
    // exported facets between rings on the inner side of a bend).
    pl.hash.query(p - w, p + w, [&](int i) { best = std::min(best, distanceToPiece(pl, i, p, ex) - chordFactor * pl.kappa[i]); });
    const Vec3 lo = p - Vec3{reach, reach, reach}, hi = p + Vec3{reach, reach, reach};
    for (int i : steep) {
        const Vec3 &a = pl.hash.lo(i), &b = pl.hash.hi(i);
        if (a.x > hi.x || b.x < lo.x || a.y > hi.y || b.y < lo.y || a.z > hi.z || b.z < lo.z) continue;
        best = std::min(best, distanceToPiece(pl, i, p, ex) - chordFactor * pl.kappa[i]);
    }
    return best - pl.eChord - r;
}

ExportTubeIndex::ExportTubeIndex(const CurveSet& cs, const PatternParams& pp, ExportTubeSpec spec, const ExecPolicy& exec,
                                 const SurfaceModel* shared, uint64_t sharedInputs) {
    auto g = std::make_unique<ExportTubeGeometry>();
    g->cs = &cs;
    g->spec = std::move(spec);
    const ExportTubeSpec& sp = g->spec;
    SurfaceParams surfaceParams;
    surfaceParams.compression = sp.compression;
    surfaceParams.maxCompression = sp.maxCompression;
    surfaceParams.ringSpacing = sp.modelRingSpacing;
    g->sharedSurface = shared && sharedInputs == surfaceBuildDigest(cs, pp, surfaceParams);
    SurfaceModel own;
    if (!g->sharedSurface) own = buildSurface(cs, pp, surfaceParams, exec);
    const SurfaceModel& surface = g->sharedSurface ? *shared : own;
    TubeParams tp;
    tp.ringSpacing = sp.ringSpacing;
    tp.sides = sp.sides;
    TubeMesh tm = buildTubeMesh(cs, surface, tp, exec);
    std::vector<Vec3> E(tm.pos.size());  // exported vertices (export units); not kept
    parallelFor(static_cast<int64_t>(E.size()), exec.chunk > 0 ? exec.chunk : 4096, [&](int64_t b, int64_t e, int64_t) {
        cancelCheckpoint();
        for (int64_t i = b; i < e; ++i) E[i] = sp.toExport(tm.pos[i]);
    }, exec);
    g->pl = buildPolyline(cs, pp, exec);
    g->r = pp.radius();
    g->chordFactor = sp.ringSpacing * sp.ringSpacing / 8;
    const double foldCurvature = 1.0 / g->r;
    g->bulgeMax = g->chordFactor * foldCurvature;
    for (size_t i = 0; i < g->pl.kappa.size(); ++i)
        if (g->pl.kappa[i] > foldCurvature) g->steep.push_back(static_cast<int>(i));
    for (const YarnSurface& ys : surface.yarns) {
        g->yarnLength.push_back(ys.length);
        for (const Ring& R : ys.rings) g->rhoMax = std::max({g->rhoMax, R.a, R.b});
    }
    // Export units per mm (the mapping is a similarity up to float rounding).
    g->scale = norm(sp.toExport(Vec3{1000, 0, 0}) - sp.toExport(Vec3{0, 0, 0})) / 1000;
    g->yarnFirstVertex = std::move(tm.yarnFirstVertex);
    g->slabs.resize(cs.yarns.size());
    g->slabHash = SpatialHash3(2 * sp.ringSpacing * g->scale);
    const int stride = sp.sides + 1;
    for (size_t y = 0; y < cs.yarns.size(); ++y) {
        const uint32_t base = g->yarnFirstVertex[y];
        const int rings = g->ringsOf(static_cast<int>(y));
        // The slabs in parallel (each from the exported vertices only), hashed in order.
        g->slabs[y].resize(static_cast<size_t>(std::max(rings - 1, 0)));
        parallelFor(std::max(rings - 1, 0), exec.chunk > 0 ? exec.chunk : 256, [&](int64_t b, int64_t e, int64_t) {
            cancelCheckpoint();
            for (int64_t j = b; j < e; ++j) g->slabs[y][j] = makeSlab(E, base, sp.sides, static_cast<int>(j));
        }, exec);
        for (int j = 0; j + 1 < rings; ++j) {
            Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
            for (int ring = j; ring <= j + 1; ++ring)
                for (int i = 0; i < sp.sides; ++i) {
                    lo = minV(lo, E[base + ring * stride + i]);
                    hi = maxV(hi, E[base + ring * stride + i]);
                }
            g->slabHash.insert(lo, hi);
            g->slabId.push_back({static_cast<int>(y), j});
        }
    }
    g->slabHash.finalize();
    g_ = std::move(g);
}

ExportTubeIndex::~ExportTubeIndex() = default;
double ExportTubeIndex::scale() const { return g_->scale; }
double ExportTubeIndex::rhoMax() const { return g_->rhoMax; }
bool ExportTubeIndex::sharedSurface() const { return g_->sharedSurface; }

double ExportTubeIndex::depth(int y, double s0, double s1, const Vec3& p, int* slabs) const {
    const ExportTubeGeometry& g = *g_;
    const int N = g.ringsOf(y) - 1;
    const double w = g.yarnLength[y] / N;
    const int j0 = std::max(0, static_cast<int>(std::floor(s0 / w))), j1 = std::min(N - 1, static_cast<int>(std::floor(s1 / w)));
    const Vec3 Q = g.spec.toExport(p);
    double best = -1e300;
    for (int j = j0; j <= j1; ++j) best = std::max(best, g.slabDepth(y, j, Q));
    if (slabs) *slabs = std::max(0, j1 - j0 + 1);
    return best <= -1e299 ? best : best / g.scale;
}

namespace {

double sampleStep(const Fibre& f) { return std::max(0.5 * f.radius0, 0.005); }

// The arc of the own tube that can hold the fibre: its length, the largest ring and a ring spacing
// around the anchor.
double ownReach(const ExportTubeIndex& tubes, const ExportTubeSpec& spec, const std::vector<Vec3>& pts) {
    double total = 0;
    for (size_t j = 1; j < pts.size(); ++j) total += norm(pts[j] - pts[j - 1]);
    return total + tubes.rhoMax() + spec.ringSpacing;
}

// The whole polyline, sampled from the root, leaves the own exported tube once and never re-enters.
bool leavesOwnTubeOnce(const ExportTubeIndex& tubes, const Fibre& f, const std::vector<Vec3>& pts, double W, double h) {
    const auto inside = [&](const Vec3& p) { return tubes.depth(f.yarnId, f.s - W, f.s + W, p) > 0; };
    int transitions = 0;
    bool prev = inside(pts[0]);
    const bool first = prev;
    for (size_t j = 0; j + 1 < pts.size(); ++j) {
        const double len = norm(pts[j + 1] - pts[j]);
        const int n = std::max(2, static_cast<int>(std::ceil(len / h)));
        for (int k = 1; k <= n; ++k) {
            const bool in = inside(lerp(pts[j], pts[j + 1], double(k) / n));
            transitions += in != prev;
            prev = in;
        }
    }
    return first && transitions == 1 && !prev;
}

// Every sample keeps its clearance. A near miss of the bound is decided on the exported slabs:
// the exact lower bound (convex: the largest signed plane distance), the own window excluded by arc.
bool staysClear(const ExportTubeGeometry& g, const Fibre& f, const std::vector<Vec3>& pts, double h) {
    const YarnPath& yp = g.cs->yarns[f.yarnId];
    const double ri = f.radius0, attachLen = f.attach;
    const double win = std::max(0.0, double(f.attachArc)) + g.r + ri;
    const ExcludedArc own{f.yarnId, yp.paramAtArcLength(std::max(0.0, f.s - win)), yp.paramAtArcLength(std::min(yp.length(), f.s + win))};
    const ExcludedArc nothing;
    const double need = ri + 0.5 * h + g.spec.eSer;
    const double reach = g.r + need + g.pl.eChord + 0.5;
    double l = -f.emergeL;
    for (size_t j = 0; j + 1 < pts.size(); ++j) {
        const double len = norm(pts[j + 1] - pts[j]);
        const int n = std::max(1, static_cast<int>(std::ceil(len / h)));
        for (int k = (j == 0 ? 0 : 1); k <= n; ++k) {
            const double t = double(k) / n;
            const Vec3 q = lerp(pts[j], pts[j + 1], t);
            const bool excluded = l + t * len <= attachLen;
            double c = g.clearance(q, reach, need, excluded ? own : nothing);
            if (c < need) {
                const Vec3 Q = g.spec.toExport(q);
                const double reachE = (need + 0.05) * g.scale;
                double cE = reachE;
                g.slabHash.query(Q - Vec3{reachE, reachE, reachE}, Q + Vec3{reachE, reachE, reachE}, [&](int id) {
                    const auto [y, slab] = g.slabId[id];
                    if (excluded && y == f.yarnId) {
                        const double w = g.yarnLength[y] / double(g.slabs[y].size());
                        if ((slab + 1) * w > f.s - win && slab * w < f.s + win) return;
                    }
                    cE = std::min(cE, g.slabDistanceLower(y, slab, Q));
                });
                c = std::max(c, cE / g.scale - g.spec.eSer);
            }
            if (c < need) return false;
        }
        l += len;
    }
    return true;
}

} // namespace

ExportTubeIndex::TubeVerdict ExportTubeIndex::check(const FibreSet& fs, const Fibre& f, double margin) const {
    TubeVerdict v;
    const std::vector<Vec3> pts = fibrePoints(fs, f);
    if (pts.size() < 2) return v;
    const double W = ownReach(*this, g_->spec, pts), h = sampleStep(f);
    v.rootDepth = depth(f.yarnId, f.s - W, f.s + W, pts[0], &v.slabs);
    v.buried = v.rootDepth >= f.radius0 + margin;
    v.singleExit = leavesOwnTubeOnce(*this, f, pts, W, h);
    v.clear = staysClear(*g_, f, pts, h);
    return v;
}

bool ExportTubeIndex::accepts(const FibreSet& fs, const Fibre& f, double margin) const {
    const std::vector<Vec3> pts = fibrePoints(fs, f);
    if (pts.size() < 2) return false;
    const double W = ownReach(*this, g_->spec, pts), h = sampleStep(f);
    return depth(f.yarnId, f.s - W, f.s + W, pts[0]) >= f.radius0 + margin && leavesOwnTubeOnce(*this, f, pts, W, h) && staysClear(*g_, f, pts, h);
}

// --- the reference check

FibreReferenceResult fibreReferenceCheck(const ExportTubeIndex& tubes, const FibreSet& fs, double hMin, const ExecPolicy& exec) {
    const auto t0 = Clock::now();
    const ExportTubeGeometry& g = tubes.geometry();
    const CentrePolyline& pl = g.pl;
    const double r = g.r, e = pl.eChord;
    FibreReferenceResult out;
    out.eChord = e;
    const int64_t chunk = exec.chunk > 0 ? exec.chunk : 256;
    std::vector<FibreReferenceResult> chunkResults(static_cast<size_t>(fibreChunks(fs, exec)));
    parallelFor(static_cast<int64_t>(fs.fibres.size()), chunk, [&](int64_t begin, int64_t end, int64_t ci) {
        FibreReferenceResult& res = chunkResults[static_cast<size_t>(ci)];
        cancelCheckpoint();  // per chunk: a cancelled job does not finish the diagnostics
        for (size_t fi = static_cast<size_t>(begin); fi < static_cast<size_t>(end); ++fi) {
            const Fibre& f = fs.fibres[fi];
            const std::vector<Vec3> pts = fibrePoints(fs, f);
            const YarnPath& yp = g.cs->yarns[f.yarnId];
            // As centreline balls the excluded cross-sections reach r + r_f further along the arc.
            const double win = f.attachArc + r + f.radius0;
            const ExcludedArc own{f.yarnId, yp.paramAtArcLength(std::max(0.0, f.s - win)), yp.paramAtArcLength(std::min(yp.length(), f.s + win))};
            const ExcludedArc nothing;
            const double R = r + f.radius0;
            bool failed = false, unresolved = false;
            double l = -f.emergeL;
            for (size_t j = 0; j + 1 < pts.size() && !failed; ++j) {
                const Vec3 A = pts[j], B = pts[j + 1];
                const double len = norm(B - A);
                ++res.segments;
                // Split at the end of the attachment exception so each part has one rule.
                std::vector<std::pair<double, double>> parts;
                const double tSplit = len > 0 ? (f.attach - l) / len : 2.0;
                if (tSplit > 0 && tSplit < 1) parts = {{0.0, tSplit}, {tSplit, 1.0}};
                else parts = {{0.0, 1.0}};
                for (const auto& [ta, tb] : parts) {
                    const bool attached = l + 0.5 * (ta + tb) * len <= f.attach;
                    std::vector<std::pair<double, double>> stack{{ta, tb}};
                    while (!stack.empty() && !failed) {
                        const auto [a, b] = stack.back();
                        stack.pop_back();
                        ++res.intervals;
                        const Vec3 m = A + (B - A) * (0.5 * (a + b));
                        const double h = 0.5 * (b - a) * len;
                        const double d = polyDistance(pl, m, R + e + h, attached ? own : nothing);
                        res.minMargin = std::min(res.minMargin, d - e - h - R);
                        if (d - e - h >= R) continue;
                        if (d + e < R) { failed = true; break; }
                        if (h <= hMin) { unresolved = true; continue; }
                        stack.push_back({a, 0.5 * (a + b)});
                        stack.push_back({0.5 * (a + b), b});
                    }
                }
                l += len;
            }
            if (failed) {
                ++res.fail;
                if (res.failing.size() < kListedFibres) res.failing.push_back(static_cast<int>(fi));
            } else if (unresolved) {
                ++res.unresolved;
                if (res.unresolvedFibres.size() < kListedFibres) res.unresolvedFibres.push_back(static_cast<int>(fi));
            } else {
                ++res.ok;
            }
        }
    }, exec);
    for (const FibreReferenceResult& p : chunkResults) {
        out.segments += p.segments;
        out.intervals += p.intervals;
        out.ok += p.ok;
        out.fail += p.fail;
        out.unresolved += p.unresolved;
        out.minMargin = std::min(out.minMargin, p.minMargin);
        appendListed(out.failing, p.failing);
        appendListed(out.unresolvedFibres, p.unresolvedFibres);
    }
    out.ms = msSince(t0);
    return out;
}

// --- the attachment check

AttachmentResult rootAttachmentCheck(const ExportTubeIndex& tubes, const FibreSet& fs, double margin, const ExecPolicy& exec) {
    const auto t0 = Clock::now();
    const double ringSpacing = tubes.geometry().spec.ringSpacing;
    AttachmentResult out;
    const int64_t chunk = exec.chunk > 0 ? exec.chunk : 256;
    std::vector<AttachmentResult> parts(static_cast<size_t>(fibreChunks(fs, exec)));
    parallelFor(static_cast<int64_t>(fs.fibres.size()), chunk, [&](int64_t begin, int64_t end, int64_t ci) {
        AttachmentResult& res = parts[static_cast<size_t>(ci)];
        cancelCheckpoint();  // per chunk: a cancelled job does not finish the diagnostics
        for (size_t fi = static_cast<size_t>(begin); fi < static_cast<size_t>(end); ++fi) {
            const Fibre& f = fs.fibres[fi];
            const std::vector<Vec3> pts = fibrePoints(fs, f);
            ++res.fibres;
            const double W = norm(pts[1] - pts[0]) + tubes.rhoMax() + ringSpacing;
            int slabs = 0;
            const double dRoot = tubes.depth(f.yarnId, f.s - W, f.s + W, pts[0], &slabs) * tubes.scale();  // export units
            res.maxSlabs = std::max(res.maxSlabs, slabs);
            res.minRootDepth = std::min(res.minRootDepth, dRoot);
            if (dRoot > margin) ++res.rootInside;
            else if (res.detached.size() < kListedFibres) res.detached.push_back(static_cast<int>(fi));
            if (tubes.depth(f.yarnId, f.s - W, f.s + W, pts.back()) <= 0) ++res.tipOutside;
            const ExportTubeIndex::TubeVerdict v = tubes.check(fs, f, margin / tubes.scale());
            res.buried += v.buried;
            res.singleExit += v.singleExit;
            res.clear += v.clear;
            res.maxSlabs = std::max(res.maxSlabs, v.slabs);
        }
    }, exec);
    for (const AttachmentResult& p : parts) {
        out.fibres += p.fibres;
        out.rootInside += p.rootInside;
        out.tipOutside += p.tipOutside;
        out.minRootDepth = std::min(out.minRootDepth, p.minRootDepth);
        appendListed(out.detached, p.detached);
        out.buried += p.buried;
        out.singleExit += p.singleExit;
        out.clear += p.clear;
        out.maxSlabs = std::max(out.maxSlabs, p.maxSlabs);
    }
    out.ms = msSince(t0);
    return out;
}

// --- the measured angles

namespace {

Vec3 tangentAtArc(const YarnPath& y, double s) {
    const double u = y.paramAtArcLength(std::clamp(s, 0.0, y.length()));
    return normalized(y.curve.d1(u));
}

double angleDeg(const Vec3& d, const Vec3& T) { return std::acos(std::clamp(std::fabs(dot(normalized(d), T)), 0.0, 1.0)) * 180 / kPi; }

void quantiles(std::vector<double> v, double out[3]) {
    if (v.empty()) return;
    std::sort(v.begin(), v.end());
    const auto q = [&](double p) { return v[std::min(v.size() - 1, static_cast<size_t>(p * (v.size() - 1) + 0.5))]; };
    out[0] = q(0.5);
    out[1] = q(0.95);
    out[2] = v.back();
}

} // namespace

FibreAngles measureFibreAngles(const CurveSet& yarns, const FibreSet& fs, const ExecPolicy& exec) {
    FibreAngles out;
    // Per fibre (independent), then collected in fibre order: the same values in the same order.
    std::vector<double> emF(fs.fibres.size(), -1), alF(fs.fibres.size(), -1);
    parallelFor(static_cast<int64_t>(fs.fibres.size()), 256, [&](int64_t b, int64_t e, int64_t) {
        for (int64_t fi = b; fi < e; ++fi) {
            const Fibre& f = fs.fibres[static_cast<size_t>(fi)];
            const std::vector<Vec3> pts = fibrePoints(fs, f);
            if (pts.size() < 3) continue;
            const int yi = f.parentYarn >= 0 ? f.parentYarn : f.yarnId;
            if (yi < 0 || yi >= static_cast<int>(yarns.yarns.size())) continue;
            const YarnPath& y = yarns.yarns[yi];
            emF[static_cast<size_t>(fi)] = angleDeg(pts[2] - pts[1], tangentAtArc(y, f.s));
            double visible = 0;
            for (size_t k = 1; k + 1 < pts.size(); ++k) visible += norm(pts[k + 1] - pts[k]);
            const double reach = visible + 1.0;
            double worst = 0;
            for (size_t k = 1; k + 1 < pts.size(); ++k) {
                const Vec3 m = (pts[k] + pts[k + 1]) * 0.5;
                double best = 1e300, nearestS = f.s;
                for (double s = f.s - reach; s <= f.s + reach; s += 0.02) {
                    const double sc = std::clamp(s, 0.0, y.length());
                    const double d = norm(y.positionAtParam(y.paramAtArcLength(sc)) - m);
                    if (d < best) {
                        best = d;
                        nearestS = sc;
                    }
                }
                worst = std::max(worst, angleDeg(pts[k + 1] - pts[k], tangentAtArc(y, nearestS)));
            }
            alF[static_cast<size_t>(fi)] = worst;
        }
    }, exec);
    std::vector<double> emergence[2], along[2];
    for (size_t fi = 0; fi < fs.fibres.size(); ++fi) {
        if (emF[fi] < 0) continue;  // skipped (angles are >= 0)
        const int p = static_cast<int>(fs.fibres[fi].population);
        emergence[p].push_back(emF[fi]);
        along[p].push_back(alF[fi]);
    }
    for (int p = 0; p < 2; ++p) {
        out.pop[p].fibres = static_cast<int>(emergence[p].size());
        quantiles(emergence[p], out.pop[p].emergence);
        quantiles(along[p], out.pop[p].along);
    }
    return out;
}

// --- snapshots

namespace {

void hexd(std::string& o, double v) {
    char b[64];
    std::snprintf(b, sizeof b, " %a", v);
    o += b;
}
void hexv(std::string& o, const Vec3& p) {
    hexd(o, p.x);
    hexd(o, p.y);
    hexd(o, p.z);
}
void num(std::string& o, long long v) { o += " " + std::to_string(v); }
void ints(std::string& o, const std::vector<int>& v) {
    o += " [";
    for (int i : v) num(o, i);
    o += " ]";
}

} // namespace

std::string reportSnapshot(const FibreReport& r) {
    std::string o = "report";
    for (int FibreReport::*const counter : fibreReportCounters()) num(o, r.*counter);
    hexd(o, r.totalLength);  // r.ms: timing, excluded
    o += " pop " + std::to_string(r.pop.size());
    for (const FibreReport& p : r.pop) o += " {" + reportSnapshot(p) + "}";
    return o;
}

namespace {

// One fibre's exact text: every deterministic field, its points last.
std::string fibreLine(const FibreSet& fs, const Fibre& f) {
    const std::span<const Vec3f> pts = fs.pointsOf(f);
    std::string o = "fibre";
    num(o, f.yarnId);
    num(o, f.parentYarn);
    num(o, static_cast<int>(f.population));
    hexd(o, f.s);
    hexd(o, f.theta);
    hexv(o, toVec3(f.anchor));
    hexv(o, toVec3(f.normal));
    num(o, f.truncated);
    num(o, 1);  // retired fields keep their place: stored digests stay valid
    num(o, 0);
    hexv(o, toVec3(pts[0]));
    hexd(o, f.targetLength);
    hexd(o, f.radius0);
    hexd(o, f.radius1);
    num(o, 1);
    hexd(o, f.emergeL);
    hexd(o, f.attach);
    hexd(o, f.attachArc);
    hexd(o, f.key);
    for (float v : {f.psi, f.lift, f.smoothLength, f.emergeAngle, f.helixActual, f.maxDev, f.minTurn, f.maxTurn}) hexd(o, v);
    o += " pts " + std::to_string(pts.size());
    for (const Vec3f& p : pts) hexv(o, toVec3(p));
    return o;
}

} // namespace

uint64_t fibreSnapshotChunkDigest(const FibreSet& fs, const ExecPolicy& exec) {
    const int64_t chunk = 4096, n = static_cast<int64_t>(fs.fibres.size());
    std::vector<uint64_t> part(static_cast<size_t>((n + chunk - 1) / chunk));
    parallelFor(n, chunk, [&](int64_t b, int64_t e, int64_t ci) {
        cancelCheckpoint();
        std::string t;
        for (int64_t i = b; i < e; ++i) t += fibreLine(fs, fs.fibres[i]) + "\n";
        part[ci] = fnv1a64(t);
    }, exec);
    std::string all = reportSnapshot(fs.rep) + "\n";
    char b[24];
    for (uint64_t d : part) {
        std::snprintf(b, sizeof b, "%016llx\n", static_cast<unsigned long long>(d));
        all += b;
    }
    return fnv1a64(all);
}

std::string fibreSnapshot(const FibreSet& fs) {
    std::string o = reportSnapshot(fs.rep) + "\n";
    for (const Fibre& f : fs.fibres) o += fibreLine(fs, f) + "\n";
    return o;
}

std::string checkSnapshot(const FibreReferenceResult& r) {
    std::string o = "reference";
    num(o, r.segments);
    num(o, r.intervals);
    num(o, r.ok);
    num(o, r.fail);
    num(o, r.unresolved);
    hexd(o, r.eChord);
    hexd(o, r.minMargin);
    ints(o, r.failing);
    ints(o, r.unresolvedFibres);
    return o;  // r.ms: timing, excluded
}

std::string checkSnapshot(const AttachmentResult& r) {
    std::string o = "attachment";
    for (int v : {r.fibres, r.rootInside, r.tipOutside, r.fibres, r.buried, r.singleExit, r.clear, r.maxSlabs}) num(o, v);
    hexd(o, r.minRootDepth);
    ints(o, r.detached);
    return o;  // r.ms: timing, excluded
}

} // namespace knit
