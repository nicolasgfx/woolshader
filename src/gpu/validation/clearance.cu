// Gate (a): samplePath (src/model/yarn_path.cpp) and checkClearance
// (src/knitting/checks/clearance.cpp) on the device (checks.h).
//
// Sampling, per span (one thread): n = max(1, ceil(sqrt(M / (8 eps)))) capped at maxPerSpan,
// M = max |C''| at the span ends; then the vertices, u, chord / polygon length bounds and segEps at
// precomputed offsets; prefix sums (double) on the host in order; the arc-length estimate s of every
// vertex on the device (arcLengthAt: table + local G7/K15). The chord / polygon bounds use 1 -+ 1e-6
// (CPU 1 -+ 1e-12).
//
// Clearance has the structure of the CPU: chunks of 16 polyline segments; candidate chunk pairs (i <= j,
// box distance <= reach = D + c + delta) from a device grid; every segment pair of a candidate
// pair: inter-yarn classification against D + c, same-yarn sameYarnPair with its Schur /
// curvature windows and the subdivision (an explicit DFS stack in the CPU's recursion order,
// depth limit 12; a full stack would be UNRESOLVED "GPU work stack exhausted", never silent).
// One thread per (chunk pair, segment of the first chunk).
//
// Issues: every issue is counted per chunk row; the device stores them (all, or - if more than
// the buffer holds - exactly the rows up to the first 200000 in the CPU's row order, in a
// second pass); the host sorts them into the CPU's order (row ci, column cj, segment i,
// segment j, recursion order), keeps the first 200000 (kMaxIssues) and merges regions with
// the CPU's code. status / counters / minima cover every pair and issue.
//
// Float allowances (documented constants):
//   etaX       2e-7 Lref per position (CPU 1e-12 Lref): ~3 ulps of the largest coordinate
//   window L   factor (1 + 1e-6) (CPU 1 + 1e-9); window maths in float (rho clamped to 1e30)
//   lengthBounds: math.cuh kBoundRel 4e-6 (CPU 1e-12)
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstring>

#include "gpu/validation/checks.h"
#include "gpu/geometry/math.cuh"

namespace knit::gpu {

namespace {

constexpr float kChordRel = 1e-6f;

__global__ void countVerticesKernel(SpanList spans, float eps, int maxPerSpan, int* nOut, float* epsOut) {
    const CurveSetView& cs = spans.curves;
    const int g = threadItem();
    if (g >= spans.count) return;
    const CurveRecord fc = cs.curves[spans.yarn[g]];
    const Bezier z = localSpan(cs.ctrl + 3 * fc.ctrlOff, g - fc.spanOff).z;
    const f3 a = (z.b[2] - z.b[1] * 2.0f + z.b[0]) * 6.0f, b = (z.b[3] - z.b[2] * 2.0f + z.b[1]) * 6.0f;
    const float M = fmaxf_(norm(a), norm(b));
    int n = static_cast<int>(ceilf(sqrtf(M / (8.0f * eps))));
    n = n < 1 ? 1 : n;
    n = n < maxPerSpan ? n : maxPerSpan;
    nOut[g] = n;
    epsOut[g] = M / (8.0f * n * n);
}

struct PolylineOutput {
    float* p;
    double* u;
    double* segLo;
    double* segHi;
    float* segEps;
};

// Vertex `first[g] + j - 1` of span g (j = 1..n): the end of sub-piece j.
__global__ void fillVerticesKernel(SpanList spans, const int* nPer, const int* first, const float* spanEps, PolylineOutput out) {
    const CurveSetView& cs = spans.curves;
    const int g = threadItem();
    if (g >= spans.count) return;
    const CurveRecord fc = cs.curves[spans.yarn[g]];
    const int k = g - fc.spanOff;
    const LocalSpan s = localSpan(cs.ctrl + 3 * fc.ctrlOff, k);
    const int n = nPer[g];
    for (int j = 1; j <= n; ++j) {
        const float t0 = float(j - 1) / n, t1 = float(j) / n;
        const Bezier piece = bezSub(s.z, t0, t1);
        const int v = first[g] + j - 1;  // segment index == vertex index - 1 within the polyline buffer
        const f3 q = piece.b[3] + s.o;
        out.p[3 * (v + 1)] = q.x;
        out.p[3 * (v + 1) + 1] = q.y;
        out.p[3 * (v + 1) + 2] = q.z;
        out.u[v + 1] = k + double(t1);
        out.segLo[v] = double(norm(piece.b[3] - piece.b[0])) * (1.0 - kChordRel);
        out.segHi[v] = double(polygonLength(piece)) * (1.0 + kChordRel);
        out.segEps[v] = spanEps[g];
    }
}

__global__ void firstVertexKernel(CurveSetView cs, int nCurves, const int* vertexOff, float* p, double* u) {
    const int y = threadItem();
    if (y >= nCurves) return;
    const CurveRecord fc = cs.curves[y];
    const int v = vertexOff[y];
    u[v] = 0.0;
    if (fc.nCtrl < 4) return;
    const LocalSpan s = localSpan(cs.ctrl + 3 * fc.ctrlOff, 0);
    const f3 q = s.z.b[0] + s.o;
    p[3 * v] = q.x;
    p[3 * v + 1] = q.y;
    p[3 * v + 2] = q.z;
}

__global__ void arcKernel(CurveSetView cs, const int* vertYarn, int nVert, const double* u, double* s) {
    const int v = threadItem();
    if (v >= nVert) return;
    s[v] = arcLengthAt(devCurve(cs, vertYarn[v]), u[v]);
}

constexpr size_t kMaxIssues = 200000;  // checkClearance's cap (issues kept)
constexpr float kInfF = 1e30f;

struct Chunk {
    int y, s0, s1;
    float lo[3], hi[3];
};

struct YarnLayout {
    int vOff, nSpans;
    int stOff, stLevels;  // sparse table of the per-span kappa bounds (stLevels 0: none)
    float rhoG, LG, thrG;
    int certG, tightG;
};

struct Window {
    float L, rho, threshold;
    bool tightBend, certifiable;
    __device__ float chordLB(float ds) const { return 2 * rho * sinf(fminf_(ds, kPiF * rho) / (2 * rho)); }
};

__device__ Window makeWindow(float rho, bool certifiable, float r, float D, float Dc) {
    Window w;
    const float Rc = 0.5f * Dc;
    w.rho = rho;
    w.tightBend = false;
    w.certifiable = true;
    if (certifiable && rho >= Rc) {
        w.L = 2 * rho * asinf(fminf_(1.0f, Rc / rho)) * (1 + 1e-6f);
        w.threshold = Dc;
    } else if (certifiable && rho >= r) {
        w.L = 2 * rho * asinf(fminf_(1.0f, r / rho)) * (1 + 1e-6f);
        w.threshold = D;
        w.tightBend = true;
    } else {
        w.L = kPiF * Rc;
        w.threshold = Dc;
        w.certifiable = false;
    }
    return w;
}

struct ClearanceIssueRecord {
    int ci, cj, i, j, seq;
    int ya, yb;
    unsigned char status, pen, reason, pad;
    float d, thr;
    double uA, uB;
    float pA[3], pB[3];
};

struct ClearanceCounters {
    unsigned long long segmentPairs, schur;
    unsigned minInter, minSelf;  // float bits (non-negative)
    unsigned anyFail, anyUnres;
    unsigned long long issues;   // all issues (stored or not)
    unsigned long long stored;   // slots taken
};

struct ClearanceArgs {
    CurveSetView cs;
    const float* p;
    const double* u;
    const float* segEps;
    const double* cumLo;
    const double* cumHi;
    const YarnLayout* yi;
    const double* st;
    const Chunk* ch;
    float r, D, Dc, etaX2;
    ClearanceIssueRecord* issues;
    unsigned long long cap;
    int rowLimit;  // store issues of rows <= rowLimit (INT_MAX: all)
    unsigned long long* rowCount;
    ClearanceCounters* acc;
};

__device__ __forceinline__ f3 P(const ClearanceArgs& c, int y, int v) { return fromF(c.p + 3 * (c.yi[y].vOff + v)); }

__device__ Window windowFor(const ClearanceArgs& c, int y, double u0, double u1) {
    const YarnLayout& Y = c.yi[y];
    Window g;
    g.rho = Y.rhoG;
    g.L = Y.LG;
    g.threshold = Y.thrG;
    g.certifiable = Y.certG != 0;
    g.tightBend = Y.tightG != 0;
    if (Y.stLevels == 0) return g;
    const int n = Y.nSpans;
    int k0 = static_cast<int>(floor(u0));
    k0 = k0 < 0 ? 0 : (k0 > n - 1 ? n - 1 : k0);
    int k1 = static_cast<int>(ceil(u1)) - 1;
    k1 = k1 < k0 ? k0 : (k1 > n - 1 ? n - 1 : k1);
    int k = 0;
    while ((2 << k) <= k1 - k0 + 1) ++k;
    // level k of the table: entries [offset of level k, ...), level size n - 2^k + 1
    long long off = Y.stOff;
    for (int l = 0; l < k; ++l) off += n - (1 << l) + 1;
    const double kap = fmax(c.st[off + k0], c.st[off + k1 - (1 << k) + 1]);
    const double rho = kap > 0 ? 1.0 / kap : 1e300;
    const double rhoUse = fmax(rho, double(Y.rhoG));
    return makeWindow(float(fmin(rhoUse, double(kInfF))), rho >= c.r, c.r, c.D, c.Dc);
}

__device__ __forceinline__ unsigned char classify(float d, float thr, float delta) {
    if (d - delta >= thr) return 0;  // Pass
    if (d + delta < thr) return 1;   // Fail
    return 2;                        // Unresolved
}

// Where an issue was found: the chunk pair and the segment pair (with the pair's issue counter, the CPU's sort keys).
struct IssueSlot {
    int ci, cj, i, j;
};
// The closest points of an issue: yarn, curve parameter and position on each side, and their distance.
struct ContactPoints {
    int ya;
    double uA;
    f3 pA;
    int yb;
    double uB;
    f3 pB;
    float d;
};
// The verdict: the threshold, the allowance, the status (1 FAIL, 2 UNRESOLVED) and why (reasonText).
struct IssueVerdict {
    float thr, delta;
    unsigned char status, reason;
};

// Counts an issue and stores it (rows up to rowLimit); seq is the pair's issue counter.
__device__ void record(const ClearanceArgs& c, const IssueSlot& where, int& seq, const ContactPoints& points, const IssueVerdict& v) {
    const int s = seq++;
    atomicAdd(&c.acc->issues, 1ull);
    atomicAdd(&c.rowCount[where.ci], 1ull);
    if (v.status == 1) atomicOr(&c.acc->anyFail, 1u);
    if (v.status == 2) atomicOr(&c.acc->anyUnres, 1u);
    if (where.ci > c.rowLimit) return;
    const unsigned long long slot = atomicAdd(&c.acc->stored, 1ull);
    if (slot >= c.cap) return;  // counted in `stored`: the host grows the buffer and reruns
    ClearanceIssueRecord& o = c.issues[slot];
    o.ci = where.ci; o.cj = where.cj; o.i = where.i; o.j = where.j; o.seq = s;
    o.ya = points.ya; o.yb = points.yb;
    o.status = v.status;
    o.pen = v.status == 1 && points.d + v.delta < c.D;
    o.reason = v.reason;
    o.d = points.d; o.thr = v.thr; o.uA = points.uA; o.uB = points.uB;
    o.pA[0] = points.pA.x; o.pA[1] = points.pA.y; o.pA[2] = points.pA.z;
    o.pB[0] = points.pB.x; o.pB[1] = points.pB.y; o.pB[2] = points.pB.z;
}

struct Item {
    double a0, a1, b0, b1, dsMin, dsMax;
    f3 pa0, pa1, pb0, pb1;
    float epsA, epsB;
    int depth;
};
constexpr int kStack = 48;  // depth 12 with 4 children needs at most 3 * 12 + 1

// What a same-yarn pair adds to its thread's totals: pairs certified by the chord bound, and the smallest
// nonlocal distance.
struct SelfPairTotals {
    unsigned long long schur = 0;
    float minSelf = 3.0e38f;
};

__device__ SelfPairTotals sameYarnPair(const ClearanceArgs& c, const IssueSlot& where, int& seq, int y, const Item& top) {
    SelfPairTotals totals;
    Item stack[kStack];
    int n = 0;
    stack[n++] = top;
    const CurveView curve = devCurve(c.cs, y);
    const CurveRecord fc = c.cs.curves[y];
    const double* preLo = c.cs.preLo && fc.preOff >= 0 ? c.cs.preLo + fc.preOff : nullptr;
    const double* preHi = c.cs.preHi && fc.preOff >= 0 ? c.cs.preHi + fc.preOff : nullptr;
    while (n > 0) {
        const Item it = stack[--n];
        const Window w = windowFor(c, y, it.a0, it.b1);
        const float delta = it.epsA + it.epsB + c.etaX2;
        float s, t;
        const float d = segmentDistance(it.pa0, it.pa1, it.pb0, it.pb1, s, t);
        const double uA = it.a0 + s * (it.a1 - it.a0), uB = it.b0 + t * (it.b1 - it.b0);
        const f3 qA = lerp(it.pa0, it.pa1, s), qB = lerp(it.pb0, it.pb1, t);
        const ContactPoints points{y, uA, qA, y, uB, qB, d};
        if (!w.certifiable) {
            if (it.dsMax >= w.L && d - delta < c.Dc) record(c, where, seq, points, {c.Dc, delta, 2, 1});
            continue;
        }
        if (it.dsMax < w.L) continue;
        const float piRho = kPiF * w.rho;
        if (it.dsMax <= piRho) {
            ++totals.schur;
            if (w.tightBend && w.chordLB(fmaxf_(w.L, float(it.dsMin))) < c.Dc && d - delta < c.Dc)
                record(c, where, seq, points, {c.Dc, delta, 2, 2});
            continue;
        }
        if (it.dsMin >= w.L) {
            totals.minSelf = fminf_(totals.minSelf, d);
            const unsigned char st = classify(d, w.threshold, delta);
            if (st == 0 && w.tightBend && d - delta < c.Dc) {
                record(c, where, seq, points, {c.Dc, delta, 2, 2});
                continue;
            }
            if (st != 0) record(c, where, seq, points, {w.threshold, delta, st, 3});
            continue;
        }
        if (d - delta >= c.Dc) continue;
        if (it.depth >= 12) {
            record(c, where, seq, points, {c.Dc, delta, 2, 4});
            continue;
        }
        if (n + 4 > kStack) {
            record(c, where, seq, points, {c.Dc, delta, 2, 5});
            continue;
        }
        const double am = 0.5 * (it.a0 + it.a1), bm = 0.5 * (it.b0 + it.b1);
        const f3 pam = evalAt(curve, am), pbm = evalAt(curve, bm);
        const double A[3] = {it.a0, am, it.a1}, B[3] = {it.b0, bm, it.b1};
        const f3 PA[3] = {it.pa0, pam, it.pa1}, PB[3] = {it.pb0, pbm, it.pb1};
        // Children in reverse, so they are popped in the CPU's order (i outer, j inner).
        for (int q = 3; q >= 0; --q) {
            const int i = q >> 1, j = q & 1;
            const ArcLengthBounds inner = lengthBounds(curve, preLo, preHi, A[i + 1], B[j]);
            const ArcLengthBounds outer = lengthBounds(curve, preLo, preHi, A[i], B[j + 1]);
            Item ch;
            ch.a0 = A[i]; ch.a1 = A[i + 1]; ch.b0 = B[j]; ch.b1 = B[j + 1];
            ch.pa0 = PA[i]; ch.pa1 = PA[i + 1]; ch.pb0 = PB[j]; ch.pb1 = PB[j + 1];
            ch.dsMin = inner.lo;
            ch.dsMax = outer.hi;
            ch.epsA = 0.25f * it.epsA;
            ch.epsB = 0.25f * it.epsB;
            ch.depth = it.depth + 1;
            stack[n++] = ch;
        }
    }
    return totals;
}

struct ChunkList {
    const Chunk* items;
    int count;
};
struct PairOutput {
    int2* pairs;
    unsigned long long capacity;
    unsigned long long* count;
};

__global__ void pairKernel(ChunkList chunks, GridView g, float reach, PairOutput out) {
    const int i = threadItem();
    if (i >= chunks.count) return;
    const Chunk A = chunks.items[i];
    const f3 lo = fromF(A.lo) - f3{reach, reach, reach}, hi = fromF(A.hi) + f3{reach, reach, reach};
    for (GridCursor cursor(g, lo, hi); cursor.next();) {
        const int j = cursor.item();
        if (j < i) continue;
        const Chunk& B = chunks.items[j];
        const float dx = fmaxf_(fmaxf_(0.0f, A.lo[0] - B.hi[0]), B.lo[0] - A.hi[0]);
        const float dy = fmaxf_(fmaxf_(0.0f, A.lo[1] - B.hi[1]), B.lo[1] - A.hi[1]);
        const float dz = fmaxf_(fmaxf_(0.0f, A.lo[2] - B.hi[2]), B.lo[2] - A.hi[2]);
        if (sqrtf(dx * dx + dy * dy + dz * dz) > reach) continue;
        const unsigned long long slot = atomicAdd(out.count, 1ull);
        if (slot < out.capacity) out.pairs[slot] = make_int2(i, j);
    }
}

__global__ void workKernel(ClearanceArgs c, const int2* pairs, long long nPairs) {
    const long long t = threadItem64();
    if (t >= nPairs * 16) return;
    const int2 pr = pairs[t >> 4];
    const int ci = pr.x, cj = pr.y;
    const Chunk A = c.ch[ci], B = c.ch[cj];
    const int i = A.s0 + static_cast<int>(t & 15);
    if (i >= A.s1) return;
    unsigned long long segPairs = 0, schur = 0;
    float minInter = 3.0e38f, minSelf = 3.0e38f;
    const int ya = A.y, yb = B.y;
    const int va = c.yi[ya].vOff, vb = c.yi[yb].vOff;
    const f3 a0 = P(c, ya, i), a1 = P(c, ya, i + 1);
    const int j0 = ci == cj ? i + 1 : B.s0;
    for (int j = j0; j < B.s1; ++j) {
        ++segPairs;
        int seq = 0;
        const f3 b0 = P(c, yb, j), b1 = P(c, yb, j + 1);
        if (ya != yb) {
            float s, tt;
            const float d = segmentDistance(a0, a1, b0, b1, s, tt);
            minInter = fminf_(minInter, d);
            const float delta = c.segEps[va + i] + c.segEps[vb + j] + c.etaX2;
            const unsigned char st = classify(d, c.Dc, delta);
            if (st != 0) {
                const double uA = c.u[va + i] + s * (c.u[va + i + 1] - c.u[va + i]);
                const double uB = c.u[vb + j] + tt * (c.u[vb + j + 1] - c.u[vb + j]);
                record(c, {ci, cj, i, j}, seq, {ya, uA, lerp(a0, a1, s), yb, uB, lerp(b0, b1, tt), d}, {c.Dc, delta, st, 0});
            }
        } else {
            Item it;
            it.a0 = c.u[va + i];
            it.a1 = c.u[va + i + 1];
            it.b0 = c.u[va + j];
            it.b1 = c.u[va + j + 1];
            it.pa0 = a0; it.pa1 = a1; it.pb0 = b0; it.pb1 = b1;
            it.dsMin = c.cumLo[va + j] - c.cumLo[va + i + 1];
            it.dsMax = c.cumHi[va + j + 1] - c.cumHi[va + i];
            it.epsA = c.segEps[va + i];
            it.epsB = c.segEps[va + j];
            it.depth = 0;
            const SelfPairTotals self = sameYarnPair(c, {ci, cj, i, j}, seq, ya, it);
            schur += self.schur;
            minSelf = fminf_(minSelf, self.minSelf);
        }
    }
    atomicAdd(&c.acc->segmentPairs, segPairs);
    if (schur) atomicAdd(&c.acc->schur, schur);
    if (minInter < 3.0e38f) atomicMin(&c.acc->minInter, __float_as_uint(minInter));
    if (minSelf < 3.0e38f) atomicMin(&c.acc->minSelf, __float_as_uint(minSelf));
}

// Host makeWindow (double, then stored as float).
void hostWindow(double rho, bool certifiable, double r, double D, double Dc, YarnLayout& y) {
    const double Rc = 0.5 * Dc;
    y.tightG = 0;
    y.certG = 1;
    if (certifiable && rho >= Rc) {
        y.LG = float(2 * rho * std::asin(std::min(1.0, Rc / rho)) * (1 + 1e-6));
        y.thrG = float(Dc);
    } else if (certifiable && rho >= r) {
        y.LG = float(2 * rho * std::asin(std::min(1.0, r / rho)) * (1 + 1e-6));
        y.thrG = float(D);
        y.tightG = 1;
    } else {
        y.LG = float(kPi * Rc);
        y.thrG = float(Dc);
        y.certG = 0;
    }
    y.rhoG = float(std::min(rho, double(kInfF)));
}

const char* reasonText(int k) {
    static const char* t[] = {"inter-yarn", "self: no curvature bound for window", "self: tight bend, clearance window not certifiable",
                              "self: nonlocal approach", "self: pair straddles local window", "GPU work stack exhausted"};
    return k >= 0 && k < 6 ? t[k] : "?";
}

// The vertices of a sampling, downloaded.
struct SampledVertices {
    std::vector<float> p, segEps;
    std::vector<double> u, segLo, segHi;
};

// Per yarn: the global window and the sparse table of the per-span kappa bounds.
std::vector<YarnLayout> yarnLayouts(const ClearanceInput& in, double r, double D, double Dc, std::vector<double>& table) {
    const CurveSet& cs = *in.curves.host;
    const int Y = static_cast<int>(cs.yarns.size());
    std::vector<YarnLayout> layouts(Y);
    for (int y = 0; y < Y; ++y) {
        const CurvatureResult& k = in.curvature[y];
        YarnLayout& Yi = layouts[y];
        Yi.vOff = in.devicePolylines.vertexOff[y];
        Yi.nSpans = cs.yarns[y].curve.numSpans();
        hostWindow(k.rhoLo, k.regular == Status::Pass && k.noFold == Status::Pass, r, D, Dc, Yi);
        Yi.stOff = static_cast<int>(table.size());
        Yi.stLevels = 0;
        if (k.regular == Status::Pass && !k.spanKappaHi.empty() && static_cast<int>(k.spanKappaHi.size()) == Yi.nSpans) {
            const std::vector<double>& v = k.spanKappaHi;
            const size_t base = table.size();
            table.insert(table.end(), v.begin(), v.end());
            int levels = 1;
            size_t prev = base;
            for (size_t len = 2; len <= v.size(); len *= 2) {
                const size_t cnt = v.size() - len + 1, at = table.size();
                table.resize(at + cnt);
                for (size_t i = 0; i < cnt; ++i) table[at + i] = std::max(table[prev + i], table[prev + i + len / 2]);
                prev = at;
                ++levels;
            }
            Yi.stLevels = levels;
        }
    }
    return layouts;
}

// Chunks of 16 segments and their boxes (the vertices are floats: exact). `diagonals`: the boxes'
// summed diagonals, which size the grid.
std::vector<Chunk> makeChunks(const std::vector<Polyline>& polylines, BoxArrays& boxes, double& diagonals) {
    std::vector<Chunk> chunks;
    for (size_t y = 0; y < polylines.size(); ++y) {
        const Polyline& pl = polylines[y];
        const int segs = static_cast<int>(pl.segments());
        for (int s0 = 0; s0 < segs; s0 += 16) {
            const int s1 = std::min(segs, s0 + 16);
            Chunk c{static_cast<int>(y), s0, s1, {0, 0, 0}, {0, 0, 0}};
            Vec3 lo = pl.p[s0], hi = pl.p[s0];
            for (int i = s0; i <= s1; ++i) {
                const Vec3& q = pl.p[i];
                lo = {std::min(lo.x, q.x), std::min(lo.y, q.y), std::min(lo.z, q.z)};
                hi = {std::max(hi.x, q.x), std::max(hi.y, q.y), std::max(hi.z, q.z)};
            }
            const double L[3] = {lo.x, lo.y, lo.z}, H[3] = {hi.x, hi.y, hi.z};
            for (int a = 0; a < 3; ++a) {
                c.lo[a] = float(L[a]);
                c.hi[a] = float(H[a]);
                boxes.lo.push_back(c.lo[a]);
                boxes.hi.push_back(c.hi[a]);
            }
            diagonals += norm(hi - lo);
            chunks.push_back(c);
        }
    }
    return chunks;
}

// Candidate chunk pairs into `pairs` (grown and rerun on overflow); returns how many.
long long candidatePairs(const ChunkList& chunks, const GridView& grid, float reach, DeviceBuffer& pairs) {
    DeviceBuffer count;
    count.alloc(8);
    unsigned long long capacity = std::max<unsigned long long>(1024, 64ull * chunks.count), found = 0;
    for (;;) {
        pairs.alloc(sizeof(int2) * capacity);
        check(cudaMemset(count.data(), 0, 8), "GPU clearance");
        pairKernel<<<blocksFor(chunks.count), kThreadsPerBlock>>>(chunks, grid, reach, {static_cast<int2*>(pairs.data()), capacity, count.ulongs()});
        checkLaunch("clearance pair kernel");
        count.download(&found, 8);
        if (found <= capacity) return static_cast<long long>(found);
        capacity = found;
    }
}

struct IssueBuffers {
    DeviceBuffer issues, counters, rowCounts;
};

// The candidate chunk pairs on the device, and the chunks they index.
struct PairWork {
    const DeviceBuffer& pairs;  // int2 per pair
    long long count;
    int chunks;
};

// One pass of the work kernel, storing the issues of rows up to rowLimit.
ClearanceCounters runWork(ClearanceArgs& c, const PairWork& work, IssueBuffers& b, int rowLimit, unsigned long long capacity) {
    b.issues.alloc(sizeof(ClearanceIssueRecord) * std::max<unsigned long long>(capacity, 1));
    ClearanceCounters start{};
    start.minInter = start.minSelf = 0x7f7fffffu;  // FLT_MAX
    check(cudaMemcpy(b.counters.data(), &start, sizeof start, cudaMemcpyHostToDevice), "GPU clearance");
    check(cudaMemset(b.rowCounts.data(), 0, 8 * size_t(work.chunks)), "GPU clearance");
    c.issues = static_cast<ClearanceIssueRecord*>(b.issues.data());
    c.cap = capacity;
    c.rowLimit = rowLimit;
    const long long threads = work.count * 16;
    if (threads > 0) workKernel<<<static_cast<unsigned>((threads + 255) / 256), 256>>>(c, static_cast<int2*>(work.pairs.data()), work.count);
    checkLaunch("clearance work kernel");
    ClearanceCounters counters;
    b.counters.download(&counters, sizeof counters);
    return counters;
}

// Every pair's counters and the stored issues in the CPU's order. Up to 2^18 issues fit the first
// pass; above that, an exact second pass of the rows that hold the first kMaxIssues.
std::vector<ClearanceIssueRecord> collectIssues(ClearanceArgs& c, const PairWork& work, ClearanceCounters& counters) {
    const unsigned long long firstCapacity = 1ull << 18;
    IssueBuffers b;
    b.counters.alloc(sizeof(ClearanceCounters));
    b.rowCounts.alloc(8 * size_t(work.chunks));
    c.rowCount = b.rowCounts.ulongs();
    c.acc = static_cast<ClearanceCounters*>(b.counters.data());
    counters = runWork(c, work, b, INT_MAX, firstCapacity);
    unsigned long long stored = counters.stored;
    if (stored > firstCapacity) {
        std::vector<unsigned long long> rows(work.chunks);
        b.rowCounts.download(rows.data(), rows.size() * sizeof(rows[0]));
        unsigned long long sum = 0;
        int lastRow = work.chunks - 1;
        for (int i = 0; i < work.chunks; ++i) {
            sum += rows[i];
            if (sum >= kMaxIssues) { lastRow = i; break; }
        }
        unsigned long long need = 0;
        for (int i = 0; i <= lastRow; ++i) need += rows[i];
        stored = runWork(c, work, b, lastRow, need).stored;
        if (stored > need) throw DeviceError("clearance: issue capture inconsistent between passes");
    }
    std::vector<ClearanceIssueRecord> issues(static_cast<size_t>(stored));
    b.issues.download(issues.data(), issues.size() * sizeof(issues[0]));
    std::sort(issues.begin(), issues.end(), [](const ClearanceIssueRecord& x, const ClearanceIssueRecord& y) {
        if (x.ci != y.ci) return x.ci < y.ci;
        if (x.cj != y.cj) return x.cj < y.cj;
        if (x.i != y.i) return x.i < y.i;
        if (x.j != y.j) return x.j < y.j;
        return x.seq < y.seq;
    });
    if (issues.size() > kMaxIssues) issues.resize(kMaxIssues);
    return issues;
}

void addIssues(const std::vector<ClearanceIssueRecord>& issues, ClearanceResult& out) {
    out.issues.reserve(issues.size());
    for (const ClearanceIssueRecord& s : issues) {
        ClearanceIssue is;
        is.yarnA = s.ya;
        is.yarnB = s.yb;
        is.uA = s.uA;
        is.uB = s.uB;
        is.pA = {s.pA[0], s.pA[1], s.pA[2]};
        is.pB = {s.pB[0], s.pB[1], s.pB[2]};
        is.distance = s.d;
        is.threshold = s.thr;
        is.status = s.status == 1 ? Status::Fail : Status::Unresolved;
        is.penetrating = s.pen != 0;
        is.reason = reasonText(s.reason);
        out.issues.push_back(std::move(is));
    }
}

// Regions: the CPU's merge (clearance.cpp).
void mergeRegions(ClearanceResult& out) {
    std::vector<int> order(out.issues.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = static_cast<int>(i);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        const ClearanceIssue &x = out.issues[a], &y = out.issues[b];
        if (x.yarnA != y.yarnA) return x.yarnA < y.yarnA;
        if (x.yarnB != y.yarnB) return x.yarnB < y.yarnB;
        return x.uA < y.uA;
    });
    for (int idx : order) {
        const ClearanceIssue& is = out.issues[idx];
        ClearanceRegion* reg = nullptr;
        for (auto it = out.regions.rbegin(); it != out.regions.rend(); ++it) {
            if (it->yarnA != is.yarnA || it->yarnB != is.yarnB) break;
            if (is.uA - it->uA1 <= 0.5 && is.uB >= it->uB0 - 1.0 && is.uB <= it->uB1 + 1.0) { reg = &*it; break; }
        }
        if (!reg) {
            out.regions.push_back({is.yarnA, is.yarnB, is.uA, is.uA, is.uB, is.uB, is.distance, is.threshold, is.status, is.penetrating, idx});
            continue;
        }
        reg->uA0 = std::min(reg->uA0, is.uA); reg->uA1 = std::max(reg->uA1, is.uA);
        reg->uB0 = std::min(reg->uB0, is.uB); reg->uB1 = std::max(reg->uB1, is.uB);
        reg->status = combine(reg->status, is.status);
        reg->penetrating = reg->penetrating || is.penetrating;
        if (is.distance < reg->worstDistance) { reg->worstDistance = is.distance; reg->worstIssue = idx; }
    }
}

// Polyline y's vertices back on the host, with the prefix sums in segment order (double), epsBound.
void fillHostPolyline(int y, const SampledVertices& v, const DevicePolylines& device, Polyline& P) {
    const int v0 = device.vertexOff[y], v1 = device.vertexOff[y + 1], nv = v1 - v0;
    P.p.resize(nv);
    P.u.assign(v.u.begin() + v0, v.u.begin() + v1);
    for (int i = 0; i < nv; ++i) P.p[i] = {v.p[3 * (v0 + i)], v.p[3 * (v0 + i) + 1], v.p[3 * (v0 + i) + 2]};
    P.segLo.assign(v.segLo.begin() + v0, v.segLo.begin() + v1 - 1);
    P.segHi.assign(v.segHi.begin() + v0, v.segHi.begin() + v1 - 1);
    P.segEps.resize(nv - 1);
    for (int s = 0; s < nv - 1; ++s) P.segEps[s] = v.segEps[v0 + s];
    for (int s = 0; s < nv - 1; ++s) P.epsBound = std::max(P.epsBound, P.segEps[s]);
    P.cumLo.assign(nv, 0.0);
    P.cumHi.assign(nv, 0.0);
    for (int s = 0; s < nv - 1; ++s) {
        P.cumLo[s + 1] = P.cumLo[s] + P.segLo[s];
        P.cumHi[s + 1] = P.cumHi[s] + P.segHi[s];
    }
}

} // namespace

DevicePolylines samplePolylinesOnDevice(const UploadedCurves& d, double eps, int maxPerSpan, std::vector<Polyline>& out) {
    const CurveSet& cs = *d.host;
    eps = std::max(eps, 1e-9);
    const int Y = static_cast<int>(cs.yarns.size());
    const int nSpans = d.nSpans;
    DeviceBuffer spanYarn, nPer, spanEps;
    spanYarn.upload(yarnOfSpans(d));
    const SpanList spans{d.dev.view, spanYarn.ints(), nSpans};
    std::vector<int> n;
    if (nSpans > 0) {
        nPer.alloc(4 * size_t(nSpans));
        spanEps.alloc(4 * size_t(nSpans));
        countVerticesKernel<<<blocksFor(nSpans), kThreadsPerBlock>>>(spans, float(eps), maxPerSpan, nPer.ints(), spanEps.floats());
        checkLaunch("sampling count kernel");
        nPer.download(n, nSpans);
    }
    // Layout: polyline y has 1 + sum(n) vertices; segment arrays are indexed by (global vertex - 1)
    // within one buffer that also holds the first vertex of every polyline (its segment slot unused).
    DevicePolylines pl;
    pl.vertexOff.assign(Y + 1, 0);
    std::vector<int> first(static_cast<size_t>(nSpans));
    int total = 0;
    for (int y = 0; y < Y; ++y) {
        pl.vertexOff[y] = total;
        ++total;  // vertex 0
        for (int g = d.dev.spanOffset[y]; g < d.dev.spanOffset[y + 1]; ++g) {
            first[g] = total - 1;  // the segment slot of span g's first sub-piece is its end vertex - 1
            total += n[g];
        }
    }
    pl.vertexOff[Y] = total;
    DeviceBuffer firstB, offB;
    pl.p.alloc(12 * size_t(total));
    pl.u.alloc(8 * size_t(total));
    pl.segLo.alloc(8 * size_t(total));
    pl.segHi.alloc(8 * size_t(total));
    pl.segEps.alloc(4 * size_t(total));
    offB.upload(pl.vertexOff);
    if (nSpans > 0) {
        firstB.upload(first);
        fillVerticesKernel<<<blocksFor(nSpans), kThreadsPerBlock>>>(spans, nPer.ints(), firstB.ints(), spanEps.floats(),
                                                                    {pl.p.floats(), pl.u.doubles(), pl.segLo.doubles(), pl.segHi.doubles(), pl.segEps.floats()});
    }
    if (Y > 0) firstVertexKernel<<<blocksFor(Y), kThreadsPerBlock>>>(d.dev.view, Y, offB.ints(), pl.p.floats(), pl.u.doubles());
    SampledVertices v;
    checkLaunch("sampling fill kernel");
    pl.p.download(v.p, 3 * size_t(total));
    pl.u.download(v.u, total);
    pl.segLo.download(v.segLo, total);
    pl.segHi.download(v.segHi, total);
    pl.segEps.download(v.segEps, total);
    out.assign(Y, Polyline{});
    std::vector<double> cumLo, cumHi;  // per yarn: segments + 1
    std::vector<int> vertYarn(static_cast<size_t>(total));
    for (int y = 0; y < Y; ++y) {
        out[y].yarnId = cs.yarns[y].id;
        fillHostPolyline(y, v, pl, out[y]);
        for (int i = pl.vertexOff[y]; i < pl.vertexOff[y + 1]; ++i) vertYarn[i] = y;
        cumLo.insert(cumLo.end(), out[y].cumLo.begin(), out[y].cumLo.end());
        cumHi.insert(cumHi.end(), out[y].cumHi.begin(), out[y].cumHi.end());
    }
    DeviceBuffer vy;
    pl.cumLo.upload(cumLo);
    pl.cumHi.upload(cumHi);
    vy.upload(vertYarn);
    pl.s.alloc(8 * size_t(total));
    if (total > 0) arcKernel<<<blocksFor(total), kThreadsPerBlock>>>(d.dev.view, vy.ints(), total, pl.u.doubles(), pl.s.doubles());
    std::vector<double> s;
    checkLaunch("sampling arc-length kernel");
    pl.s.download(s, total);
    for (int y = 0; y < Y; ++y) out[y].s.assign(s.begin() + pl.vertexOff[y], s.begin() + pl.vertexOff[y + 1]);
    return pl;
}

ClearanceResult checkClearanceOnDevice(const ClearanceInput& in) {
    const UploadedCurves& d = in.curves;
    const DevicePolylines& dp = in.devicePolylines;
    const PatternParams& pp = in.pattern;
    ClearanceResult out;
    const double r = pp.radius(), D = pp.yarnDiameter, Dc = D + pp.clearance;
    const double etaX = 2e-7 * in.allowances.Lref;
    double eps = 0;
    for (const Polyline& p : in.polylines) eps = std::max(eps, p.epsBound);
    const double delta = 2 * eps + 2 * etaX;
    out.deltaD = delta;
    std::vector<double> kappaTable;
    const std::vector<YarnLayout> layouts = yarnLayouts(in, r, D, Dc, kappaTable);
    BoxArrays boxes;
    double diagonals = 0;
    const std::vector<Chunk> chunks = makeChunks(in.polylines, boxes, diagonals);
    const int nChunks = static_cast<int>(chunks.size());
    const double reach = Dc + delta;
    if (nChunks == 0) return out;
    boxes.cell = float(std::max(diagonals / nChunks, 1e-9) + reach);
    DeviceGrid grid;
    DeviceBuffer chunkB, layoutB, tableB, pairs;
    grid.upload(boxes);
    chunkB.upload(chunks.data(), chunks.size() * sizeof(chunks[0]));
    layoutB.upload(layouts.data(), layouts.size() * sizeof(layouts[0]));
    tableB.upload(kappaTable);
    const long long nPairs = candidatePairs({static_cast<Chunk*>(chunkB.data()), nChunks}, grid.view, float(reach), pairs);
    ClearanceArgs c{};
    c.cs = d.dev.view;
    c.p = dp.p.floats();
    c.u = dp.u.doubles();
    c.segEps = dp.segEps.floats();
    c.cumLo = dp.cumLo.doubles();
    c.cumHi = dp.cumHi.doubles();
    c.yi = static_cast<YarnLayout*>(layoutB.data());
    c.st = tableB.doubles();
    c.ch = static_cast<Chunk*>(chunkB.data());
    c.r = float(r);
    c.D = float(D);
    c.Dc = float(Dc);
    c.etaX2 = float(2 * etaX);
    ClearanceCounters counters;
    const std::vector<ClearanceIssueRecord> issues = collectIssues(c, {pairs, nPairs, nChunks}, counters);
    // Counters and status of every pair and issue (the first pass).
    out.segmentPairs = static_cast<long>(counters.segmentPairs);
    out.schurCertified = static_cast<long>(counters.schur);
    float minInter, minSelf;
    std::memcpy(&minInter, &counters.minInter, 4);
    std::memcpy(&minSelf, &counters.minSelf, 4);
    out.minInterYarn = counters.minInter == 0x7f7fffffu ? 1e300 : minInter;
    out.minSelfNonlocal = counters.minSelf == 0x7f7fffffu ? 1e300 : minSelf;
    out.status = counters.anyFail ? Status::Fail : counters.anyUnres ? Status::Unresolved : Status::Pass;
    addIssues(issues, out);
    mergeRegions(out);
    return out;
}

} // namespace knit::gpu
