// Clearance (curve-set check).
//
// Polyline distances are exact segment-segment distances; the polylines are within
// eps of the curves, so true distances are within +-deltaD = 2 eps + 2 etaX.
//
// Same-yarn window (Schur's comparison theorem): with a curvature bound
// kappa <= 1/rho on the ENTIRE intervening arc, two points at arc distance
// ds <= pi*rho have chord >= 2 rho sin(ds / (2 rho)). Choosing
//   L = 2 rho asin(R / rho)   (R = r + c/2, threshold D + c)
// guarantees chord >= threshold for every ds in [L, pi*rho]:
//   ds_max <  L                     local, excluded (covered by the curvature check)
//   ds_max <= pi*rho                PASS by the theorem (any part with ds >= L)
//   ds_min >= L and ds_max > pi*rho distance check; closer than the threshold is a
//                                   genuine nonlocal approach
//   otherwise                       subdivide (L and pi*rho both inside the range)
// Window membership uses an enclosure of the arc length (chords below, Bezier
// control polygons above).
// The bound is local: for a pair of curve pieces [a0, a1] and [b0, b1] (a1 <= b0)
// every intervening arc lies in [a0, b1], so rho = 1 / max kappa_hi over the spans
// covering [a0, b1] (per-span certified bounds of the curvature check, range maximum
// by a sparse table). So one tight bend elsewhere on a long yarn does not widen the
// window of every pair, and a pair whose own arc contains the tight bend still gets
// it (tested: schur_local_tight_bend).
#include <algorithm>
#include <cmath>

#include "model/cancel.h"
#include "model/mesh.h"
#include "knitting/checks/checks.h"
#include "model/spatial_hash.h"

namespace knit {

namespace {

struct YarnWindow {
    double L = 0;           // local window length [mm]
    double rho = 0;         // certified curvature-radius lower bound [mm]
    double threshold = 0;   // D + c, or D in the tight-bend case
    bool tightBend = false; // D <= d < D + c is UNRESOLVED unless Schur gives D + c
    bool certifiable = true;
    // Schur chord lower bound for arc distance ds (valid for ds <= pi*rho).
    double chordLB(double ds) const { return 2 * rho * std::sin(std::min(ds, kPi * rho) / (2 * rho)); }
};

struct Chunk { int y; size_t s0, s1; Vec3 lo, hi; };

double boxDistance(const Chunk& a, const Chunk& b) {
    double dx = std::max({0.0, a.lo.x - b.hi.x, b.lo.x - a.hi.x});
    double dy = std::max({0.0, a.lo.y - b.hi.y, b.lo.y - a.hi.y});
    double dz = std::max({0.0, a.lo.z - b.hi.z, b.lo.z - a.hi.z});
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// Range maximum over per-span curvature bounds.
struct SparseMax {
    std::vector<std::vector<double>> t;
    void build(const std::vector<double>& v) {
        t.assign(1, v);
        for (size_t len = 2; len <= v.size(); len *= 2) {
            const std::vector<double>& prev = t.back();
            std::vector<double> next(v.size() - len + 1);
            for (size_t i = 0; i < next.size(); ++i) next[i] = std::max(prev[i], prev[i + len / 2]);
            t.push_back(std::move(next));
        }
    }
    double query(int i, int j) const {  // max over [i, j], inclusive
        int k = 0;
        while ((2 << k) <= j - i + 1) ++k;
        return std::max(t[k][i], t[k][j - (1 << k) + 1]);
    }
};

struct Ctx {
    const CurveSet& cs;
    const std::vector<Polyline>& polys;
    const PatternParams& pp;
    std::vector<YarnWindow> win;     // per yarn, from the global bound (fallback)
    std::vector<SparseMax> spanMax;  // per yarn, empty if no per-span bounds
    double delta = 0;  // largest pair allowance (reported)
    double etaX2 = 0;
    ClearanceResult res;
};

constexpr size_t kMaxIssues = 200000;  // issues kept (the status counts every one)

Status classify(double d, double thr, double delta) {
    if (d - delta >= thr) return Status::Pass;
    if (d + delta < thr) return Status::Fail;
    return Status::Unresolved;
}

// How an issue was judged: its threshold, the pair's allowance, the status and the reason.
struct IssueVerdict {
    double threshold, delta;
    Status status;
    const char* why;
};

// is holds the two closest points and their distance.
void record(const Ctx& c, ClearanceResult& out, ClearanceIssue is, const IssueVerdict& v) {
    is.threshold = v.threshold;
    is.status = v.status;
    is.reason = v.why;
    is.penetrating = v.status == Status::Fail && is.distance + v.delta < c.pp.yarnDiameter;
    if (out.issues.size() < kMaxIssues) out.issues.push_back(is);
    out.status = combine(out.status, v.status);
}

YarnWindow makeWindow(double rho, bool certifiable, double r, double D, double Dc) {
    YarnWindow w;
    const double Rc = 0.5 * Dc;
    w.rho = rho;
    // L is rounded up slightly so the theorem's inequality holds despite rounding.
    if (certifiable && rho >= Rc) {
        w.L = 2 * rho * std::asin(std::min(1.0, Rc / rho)) * (1 + 1e-9);
        w.threshold = Dc;
    } else if (certifiable && rho >= r) {
        w.L = 2 * rho * std::asin(std::min(1.0, r / rho)) * (1 + 1e-9);
        w.threshold = D;
        w.tightBend = true;
    } else {
        w.L = kPi * Rc;  // descriptive only; results are UNRESOLVED
        w.threshold = Dc;
        w.certifiable = false;
    }
    return w;
}

// Window for a same-yarn pair whose intervening arcs lie in [u0, u1]. With per-span
// bounds the pair is certifiable iff ITS arc has rho >= r (a fold elsewhere on the
// yarn does not matter to it); without them the yarn-global window applies.
YarnWindow windowFor(const Ctx& c, int y, double u0, double u1) {
    const YarnWindow& g = c.win[y];
    if (c.spanMax[y].t.empty()) return g;
    const int n = static_cast<int>(c.spanMax[y].t[0].size());
    int k0 = std::clamp(static_cast<int>(std::floor(u0)), 0, n - 1);
    int k1 = std::clamp(static_cast<int>(std::ceil(u1)) - 1, k0, n - 1);
    double kap = c.spanMax[y].query(k0, k1);
    double rho = kap > 0 ? 1.0 / kap : 1e300;
    return makeWindow(std::max(rho, g.rho), rho >= c.pp.radius(), c.pp.radius(), c.pp.yarnDiameter,
                      c.pp.yarnDiameter + c.pp.clearance);
}

// A piece [u0, u1] of a yarn with its chord and the chord's sampling bound.
struct ChordPiece {
    double u0, u1;
    Vec3 p0, p1;
    double eps;
};

// Two pieces of one yarn, a before b (a.u1 <= b.u0), and the arc distance range between their points.
struct SameYarnPair {
    int y;
    ChordPiece a, b;
    double dsMin, dsMax;
};

void checkSameYarnPair(const Ctx& c, ClearanceResult& out, const SameYarnPair& pair, int depth);

// Both L and pi*rho fall inside the pair's arc range: the four pairs of half pieces.
void subdivideSameYarnPair(const Ctx& c, ClearanceResult& out, const SameYarnPair& pair, int depth) {
    const ChordPiece& pa = pair.a;
    const ChordPiece& pb = pair.b;
    const YarnPath& yp = c.cs.yarns[pair.y];
    double am = 0.5 * (pa.u0 + pa.u1), bm = 0.5 * (pb.u0 + pb.u1);
    Vec3 pam = yp.curve.eval(am), pbm = yp.curve.eval(bm);
    const double A[3] = {pa.u0, am, pa.u1}, B[3] = {pb.u0, bm, pb.u1};
    const Vec3 PA[3] = {pa.p0, pam, pa.p1}, PB[3] = {pb.p0, pbm, pb.p1};
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) {
            // the shortest and the longest arc between the two pieces
            const double lo = yp.lengthBounds(A[i + 1], B[j]).lo, hi2 = yp.lengthBounds(A[i], B[j + 1]).hi;
            // A chord over half the interval has a quarter of the error bound.
            const ChordPiece qa{A[i], A[i + 1], PA[i], PA[i + 1], 0.25 * pa.eps};
            const ChordPiece qb{B[j], B[j + 1], PB[j], PB[j + 1], 0.25 * pb.eps};
            checkSameYarnPair(c, out, {pair.y, qa, qb, lo, hi2}, depth + 1);
        }
}

void checkSameYarnPair(const Ctx& c, ClearanceResult& out, const SameYarnPair& pair, int depth) {
    const ChordPiece& pa = pair.a;
    const ChordPiece& pb = pair.b;
    const int y = pair.y;
    const double dsMin = pair.dsMin, dsMax = pair.dsMax;
    const YarnWindow w = windowFor(c, y, pa.u0, pb.u1);
    // Allowance for this pair: each chord is within its own sampling bound.
    const double delta = pa.eps + pb.eps + c.etaX2;
    const double D = c.pp.yarnDiameter, Dc = D + c.pp.clearance;
    const ClosestSegmentPoints closestPair = closestSegmentPoints(pa.p0, pa.p1, pb.p0, pb.p1);
    const double d = closestPair.distance, s = closestPair.s, t = closestPair.t;
    ClearanceIssue at;
    at.yarnA = at.yarnB = y;
    at.uA = pa.u0 + s * (pa.u1 - pa.u0);
    at.uB = pb.u0 + t * (pb.u1 - pb.u0);
    at.pA = lerp(pa.p0, pa.p1, s);
    at.pB = lerp(pb.p0, pb.p1, t);
    at.distance = d;
    if (!w.certifiable) {
        if (dsMax >= w.L && d - delta < Dc) record(c, out, at, {Dc, delta, Status::Unresolved, "self: no curvature bound for window"});
        return;
    }
    if (dsMax < w.L) return;  // entirely local: covered by the curvature check
    const double piRho = kPi * w.rho;
    if (dsMax <= piRho) {
        // Every part of this pair with ds >= L has chord >= chordLB(max(L, dsMin)) >= threshold.
        ++out.schurCertified;
        if (w.tightBend && w.chordLB(std::max(w.L, dsMin)) < Dc && d - delta < Dc)
            record(c, out, at, {Dc, delta, Status::Unresolved, "self: tight bend, clearance window not certifiable"});
        return;
    }
    if (dsMin >= w.L) {
        out.minSelfNonlocal = std::min(out.minSelfNonlocal, d);
        Status st = classify(d, w.threshold, delta);
        if (st == Status::Pass && w.tightBend && d - delta < Dc) {
            record(c, out, at, {Dc, delta, Status::Unresolved, "self: tight bend, clearance window not certifiable"});
            return;
        }
        if (st != Status::Pass) record(c, out, at, {w.threshold, delta, st, "self: nonlocal approach"});
        return;
    }
    if (d - delta >= Dc) return;
    if (depth >= 12) {
        record(c, out, at, {Dc, delta, Status::Unresolved, "self: pair straddles local window"});
        return;
    }
    subdivideSameYarnPair(c, out, pair, depth);
}

// One yarn window per yarn from its global curvature bound, and the range maximum over the per-span
// bounds where the check has them.
void buildYarnWindows(Ctx& c, const std::vector<CurvatureResult>& curv) {
    const double r = c.pp.radius(), D = c.pp.yarnDiameter, Dc = D + c.pp.clearance;
    const int Y = static_cast<int>(c.cs.yarns.size());
    c.win.resize(Y);
    c.spanMax.resize(Y);
    for (int y = 0; y < Y; ++y) {
        const CurvatureResult& k = curv[y];
        const bool cert = k.regular == Status::Pass && k.noFold == Status::Pass;
        c.win[y] = makeWindow(k.rhoLo, cert, r, D, Dc);
        if (k.regular == Status::Pass && !k.spanKappaHi.empty() &&
            static_cast<int>(k.spanKappaHi.size()) == c.cs.yarns[y].curve.numSpans())
            c.spanMax[y].build(k.spanKappaHi);
    }
}

// Every polyline in chunks of 16 segments with their boxes.
std::vector<Chunk> segmentChunks(const std::vector<Polyline>& polys) {
    std::vector<Chunk> chunks;
    for (int y = 0; y < static_cast<int>(polys.size()); ++y) {
        cancelCheckpoint();
        const Polyline& pl = polys[y];
        for (size_t s0 = 0; s0 < pl.segments(); s0 += 16) {
            size_t s1 = std::min(pl.segments(), s0 + 16);
            Chunk ch{y, s0, s1, pl.p[s0], pl.p[s0]};
            for (size_t i = s0; i <= s1; ++i) {
                const Vec3& q = pl.p[i];
                ch.lo = {std::min(ch.lo.x, q.x), std::min(ch.lo.y, q.y), std::min(ch.lo.z, q.z)};
                ch.hi = {std::max(ch.hi.x, q.x), std::max(ch.hi.y, q.y), std::max(ch.hi.z, q.z)};
            }
            chunks.push_back(ch);
        }
    }
    return chunks;
}

// Segment i of A against segment j of B, of different yarns.
void checkInterYarnSegments(const Ctx& c, const Chunk& A, const Chunk& B, size_t i, size_t j, ClearanceResult& out) {
    const Polyline& PA = c.polys[A.y];
    const Polyline& PB = c.polys[B.y];
    const double Dc = c.pp.yarnDiameter + c.pp.clearance;
    const ClosestSegmentPoints closestPair = closestSegmentPoints(PA.p[i], PA.p[i + 1], PB.p[j], PB.p[j + 1]);
    const double d = closestPair.distance, s = closestPair.s, t = closestPair.t;
    out.minInterYarn = std::min(out.minInterYarn, d);
    const double delta = PA.segEps[i] + PB.segEps[j] + c.etaX2;
    const Status st = classify(d, Dc, delta);
    if (st == Status::Pass) return;
    ClearanceIssue at;
    at.yarnA = A.y;
    at.yarnB = B.y;
    at.uA = PA.u[i] + s * (PA.u[i + 1] - PA.u[i]);
    at.uB = PB.u[j] + t * (PB.u[j + 1] - PB.u[j]);
    at.pA = lerp(PA.p[i], PA.p[i + 1], s);
    at.pB = lerp(PB.p[j], PB.p[j + 1], t);
    at.distance = d;
    record(c, out, at, {Dc, delta, st, "inter-yarn"});
}

// Every segment pair of two chunks; within one chunk, each pair once.
void checkChunkPair(const Ctx& c, const Chunk& A, const Chunk& B, bool sameChunk, ClearanceResult& out) {
    const Polyline& PA = c.polys[A.y];
    for (size_t i = A.s0; i < A.s1; ++i) {
        size_t j0 = sameChunk ? i + 1 : B.s0;
        for (size_t j = j0; j < B.s1; ++j) {
            ++out.segmentPairs;
            if (A.y != B.y) {
                checkInterYarnSegments(c, A, B, i, j, out);
                continue;
            }
            const double dsMin = PA.cumLo[j] - PA.cumLo[i + 1];
            const double dsMax = PA.cumHi[j + 1] - PA.cumHi[i];
            const ChordPiece a{PA.u[i], PA.u[i + 1], PA.p[i], PA.p[i + 1], PA.segEps[i]};
            const ChordPiece b{PA.u[j], PA.u[j + 1], PA.p[j], PA.p[j + 1], PA.segEps[j]};
            checkSameYarnPair(c, out, {A.y, a, b, dsMin, dsMax}, 0);
        }
    }
}

// The candidate chunk pairs by rows ci, each once and in the order of one sorted pass: each block of
// rows fills its own result, and the blocks are merged in row order (issues in order up to the cap,
// status, counters, exact minima), so the result is that of one pass.
void checkChunkPairs(Ctx& c, const std::vector<Chunk>& chunks, const ExecPolicy& exec) {
    const double Dc = c.pp.yarnDiameter + c.pp.clearance;
    const double reach = Dc + c.delta;
    std::vector<Vec3> clo, chi;
    for (const Chunk& ch : chunks) { clo.push_back(ch.lo); chi.push_back(ch.hi); }
    const int nCh = static_cast<int>(chunks.size());
    double diag = 0;
    for (int i = 0; i < nCh; ++i) diag += norm(chi[i] - clo[i]);
    SpatialHash3 hash(std::max(diag / std::max(nCh, 1), 1e-9) + reach);
    for (int i = 0; i < nCh; ++i) hash.insert(clo[i], chi[i]);
    hash.finalize();
    const int64_t rowChunk = exec.chunk > 0 ? exec.chunk : 64;
    std::vector<ClearanceResult> parts(static_cast<size_t>((int64_t(nCh) + rowChunk - 1) / rowChunk));
    const Vec3 rv{reach, reach, reach};
    parallelFor(nCh, rowChunk, [&](int64_t rowBegin, int64_t rowEnd, int64_t pi) {
        ClearanceResult& out = parts[pi];
        std::vector<int> js;
        for (int64_t row = rowBegin; row < rowEnd; ++row) {
            if ((row & 255) == 0) cancelCheckpoint();
            const int ci = static_cast<int>(row);
            js.clear();
            hash.query(clo[ci] - rv, chi[ci] + rv, [&](int j) {
                if (j >= ci) js.push_back(j);
            });
            std::sort(js.begin(), js.end());
            for (int cj : js) {
                const Chunk& A = chunks[ci];
                const Chunk& B = chunks[cj];
                if (boxDistance(A, B) <= reach) checkChunkPair(c, A, B, ci == cj, out);
            }
        }
    }, exec);
    for (const ClearanceResult& p : parts) {
        for (const ClearanceIssue& is : p.issues)
            if (c.res.issues.size() < kMaxIssues) c.res.issues.push_back(is);
        c.res.status = combine(c.res.status, p.status);
        c.res.segmentPairs += p.segmentPairs;
        c.res.schurCertified += p.schurCertified;
        c.res.minInterYarn = std::min(c.res.minInterYarn, p.minInterYarn);
        c.res.minSelfNonlocal = std::min(c.res.minSelfNonlocal, p.minSelfNonlocal);
    }
}

// Regions: issues close along both yarns merged.
void groupIssuesIntoRegions(ClearanceResult& res) {
    std::vector<int> order(res.issues.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = static_cast<int>(i);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        const ClearanceIssue &x = res.issues[a], &y = res.issues[b];
        if (x.yarnA != y.yarnA) return x.yarnA < y.yarnA;
        if (x.yarnB != y.yarnB) return x.yarnB < y.yarnB;
        return x.uA < y.uA;
    });
    for (int idx : order) {
        const ClearanceIssue& is = res.issues[idx];
        ClearanceRegion* reg = nullptr;
        for (auto it = res.regions.rbegin(); it != res.regions.rend(); ++it) {
            if (it->yarnA != is.yarnA || it->yarnB != is.yarnB) break;
            if (is.uA - it->uA1 <= 0.5 && is.uB >= it->uB0 - 1.0 && is.uB <= it->uB1 + 1.0) { reg = &*it; break; }
        }
        if (!reg) {
            res.regions.push_back({is.yarnA, is.yarnB, is.uA, is.uA, is.uB, is.uB, is.distance, is.threshold,
                                   is.status, is.penetrating, idx});
            continue;
        }
        reg->uA0 = std::min(reg->uA0, is.uA); reg->uA1 = std::max(reg->uA1, is.uA);
        reg->uB0 = std::min(reg->uB0, is.uB); reg->uB1 = std::max(reg->uB1, is.uB);
        reg->status = combine(reg->status, is.status);
        reg->penetrating = reg->penetrating || is.penetrating;
        if (is.distance < reg->worstDistance) { reg->worstDistance = is.distance; reg->worstIssue = idx; }
    }
}

} // namespace

ClearanceResult checkClearance(const ClearanceInput& in, const PatternParams& pp, const Allowances& al, const ExecPolicy& exec) {
    Ctx c{in.cs, in.polys, pp, {}, {}, 0, 0, {}};
    c.etaX2 = 2 * al.etaX;
    double eps = 0;
    for (const Polyline& p : in.polys) eps = std::max(eps, p.epsBound);
    c.delta = 2 * eps + 2 * al.etaX;
    c.res.deltaD = c.delta;
    buildYarnWindows(c, in.curvature);
    checkChunkPairs(c, segmentChunks(in.polys), exec);
    groupIssuesIntoRegions(c.res);
    return c.res;
}

} // namespace knit
