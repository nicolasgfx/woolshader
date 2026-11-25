// Projected crossings (curve-set check).
//
// Front projection along N (drop z). For every pair of pieces:
//   * hull rejection on padded projected Bezier boxes;
//   * same piece: injective if its projected hodograph cone lies in an open
//     half-plane; adjacent pieces: they meet only at the join if the union cone
//     does (strict monotonicity along a direction);
//   * transversal pieces (line-direction cones disjoint) meet at most once
//     (chord-in-cone argument); Newton locates the candidate, a Krawczyk test on
//     the ACTUAL piecewise spline (all spans the box touches, clipped to the
//     physical domain) proves existence + uniqueness in the root box;
//   * depth sign is enclosed over the root box and contracted until resolved.
// Anything that cannot be resolved is reported as UNRESOLVED.
#include <algorithm>
#include <cmath>
#include <functional>

#include "knitting/checks/checks.h"

namespace knit {



namespace {

struct Box2 { Interval x, y; };
struct Box3 { Interval x, y, z; };

void accumulate(Box3& box, bool& first, const Vec3& v) {
    if (first) { box = {{v.x}, {v.y}, {v.z}}; first = false; return; }
    box.x = hull(box.x, Interval(v.x));
    box.y = hull(box.y, Interval(v.y));
    box.z = hull(box.z, Interval(v.z));
}

// The span pieces overlapping a global parameter interval [u0, u1] (may span several knots):
// span k, for k = first .. last, restricted to t in [max(0, u0 - k), min(1, u1 - k)].
struct PieceRange {
    int first, last;
    double u0, u1;
    double t0(int k) const { return std::max(0.0, u0 - k); }
    double t1(int k) const { return std::min(1.0, u1 - k); }
};

PieceRange pieceRange(const BSpline& s, double u0, double u1) {
    const int n = s.numSpans();
    u0 = std::max(0.0, u0);
    u1 = std::min(static_cast<double>(n), u1);
    const int k0 = std::min(n - 1, static_cast<int>(std::floor(u0)));
    int k1 = std::min(n - 1, static_cast<int>(std::floor(u1)));
    if (u1 == std::floor(u1) && k1 > k0) k1 = static_cast<int>(u1) - 1;
    return {k0, k1, u0, u1};
}

// Interval enclosures over [u0, u1]: hull of the restricted Bezier control values per component.
Box3 enclosePosition(const BSpline& s, double u0, double u1) {
    const PieceRange r = pieceRange(s, u0, u1);
    Box3 box;
    bool first = true;
    for (int k = r.first; k <= r.last; ++k) {
        const CubicBezier c = s.span(k).sub(r.t0(k), r.t1(k));
        for (const Vec3& v : c.b) accumulate(box, first, v);
    }
    return box;
}

Box3 encloseDerivative(const BSpline& s, double u0, double u1) {
    const PieceRange r = pieceRange(s, u0, u1);
    Box3 box;
    bool first = true;
    for (int k = r.first; k <= r.last; ++k) {
        const QuadBezier q = s.span(k).hodograph().sub(r.t0(k), r.t1(k));
        for (const Vec3& v : q.h) accumulate(box, first, v);
    }
    return box;
}

// All restricted hodograph control vectors over [u0, u1] (for cone tests).
void hodographControls(const BSpline& s, double u0, double u1, std::vector<Vec3>& out) {
    out.clear();
    const PieceRange r = pieceRange(s, u0, u1);
    for (int k = r.first; k <= r.last; ++k) {
        const QuadBezier q = s.span(k).hodograph().sub(r.t0(k), r.t1(k));
        for (const Vec3& v : q.h) out.push_back(v);
    }
}

struct Ctx {
    const CurveSet& cs;
    Allowances al;
    CrossingOptions opt;
    CrossingResult res;
    long budget = 4'000'000;
    double minNonTransversalWidth = 1e-4;  // parameter width below which tangency is unresolved
};

Box2 projBox(const Ctx& c, int y, double a0, double a1) {
    Box3 b = enclosePosition(c.cs.yarns[y].curve, a0, a1);
    double e = c.al.etaX;
    return {pad(b.x, e), pad(b.y, e)};
}

bool overlaps(const Box2& A, const Box2& B) { return A.x.overlaps(B.x) && A.y.overlaps(B.y); }

// Angular hull of projected vectors. Returns false if they do not fit in an open
// half-plane (or a vector is too short to have a direction).
bool angularHull(const std::vector<Vec3>& v, double minNorm, double& lo, double& hi) {
    std::vector<double> ang;
    ang.reserve(v.size());
    for (const Vec3& w : v) {
        double n = std::hypot(w.x, w.y);
        if (n <= minNorm) return false;
        ang.push_back(std::atan2(w.y, w.x));
    }
    std::sort(ang.begin(), ang.end());
    // Largest circular gap between consecutive angles.
    double bestGap = ang.front() + 2 * kPi - ang.back();
    size_t bestIdx = 0;  // gap precedes ang[bestIdx]
    for (size_t i = 1; i < ang.size(); ++i) {
        double g = ang[i] - ang[i - 1];
        if (g > bestGap) { bestGap = g; bestIdx = i; }
    }
    const double angTol = 1e-9;
    if (bestGap <= kPi + angTol) return false;  // spread >= pi: not in an open half-plane
    lo = ang[bestIdx];
    hi = bestIdx == 0 ? ang.back() : ang[bestIdx - 1] + 2 * kPi;
    lo -= angTol;
    hi += angTol;
    return true;
}

bool cone(const Ctx& c, int y, double a0, double a1, double& lo, double& hi) {
    std::vector<Vec3> h;
    hodographControls(c.cs.yarns[y].curve, a0, a1, h);
    return angularHull(h, 1e-9 * c.al.Lref, lo, hi);
}

bool unionCone(const Ctx& c, int y, double a0, double a1, double b0, double b1) {
    std::vector<Vec3> h, h2;
    hodographControls(c.cs.yarns[y].curve, a0, a1, h);
    hodographControls(c.cs.yarns[y].curve, b0, b1, h2);
    h.insert(h.end(), h2.begin(), h2.end());
    double lo, hi;
    return angularHull(h, 1e-9 * c.al.Lref, lo, hi);
}

// Line-direction cones [a0,a1] and [b0,b1] (each narrower than pi) disjoint mod pi?
bool disjointModPi(double a0, double a1, double b0, double b1) {
    for (int k = -4; k <= 4; ++k) {
        double s = k * kPi;
        if (a0 <= b1 + s && b0 + s <= a1) return false;
    }
    return true;
}

struct Krawczyk2 {
    bool ok = false;
    Interval Ka, Kb;
};

// Krawczyk operator for F(a,b) = A(a) - B(b) (projected), on X = Xa x Xb.
Krawczyk2 krawczyk(const Ctx& c, int ya, int yb, double a, double b, Interval Xa, Interval Xb) {
    const BSpline& A = c.cs.yarns[ya].curve;
    const BSpline& B = c.cs.yarns[yb].curve;
    Krawczyk2 out;
    Vec3 fa = A.eval(a), fb = B.eval(b);
    double e = 2 * c.al.etaX;
    Interval F0 = pad(Interval(fa.x) - Interval(fb.x), e);
    Interval F1 = pad(Interval(fa.y) - Interval(fb.y), e);
    Box3 dA = encloseDerivative(A, Xa.lo, Xa.hi);
    Box3 dB = encloseDerivative(B, Xb.lo, Xb.hi);
    auto padD = [&](Interval v) { return pad(v, c.al.rel * v.mag() + 1e-15 * c.al.Lref); };
    Interval J00 = padD(dA.x), J10 = padD(dA.y);
    Interval J01 = -padD(dB.x), J11 = -padD(dB.y);
    Vec3 pa = A.d1(a), pb = B.d1(b);
    double m00 = pa.x, m01 = -pb.x, m10 = pa.y, m11 = -pb.y;
    double det = m00 * m11 - m01 * m10;
    if (!(std::fabs(det) > 0)) return out;
    Interval Y00(m11 / det), Y01(-m01 / det), Y10(-m10 / det), Y11(m00 / det);
    Interval YF0 = Y00 * F0 + Y01 * F1;
    Interval YF1 = Y10 * F0 + Y11 * F1;
    Interval M00 = Interval(1.0) - (Y00 * J00 + Y01 * J10);
    Interval M01 = -(Y00 * J01 + Y01 * J11);
    Interval M10 = -(Y10 * J00 + Y11 * J10);
    Interval M11 = Interval(1.0) - (Y10 * J01 + Y11 * J11);
    Interval dXa = Xa - Interval(a), dXb = Xb - Interval(b);
    out.Ka = Interval(a) - YF0 + M00 * dXa + M01 * dXb;
    out.Kb = Interval(b) - YF1 + M10 * dXa + M11 * dXb;
    out.ok = out.Ka.strictlyInside(Xa) && out.Kb.strictlyInside(Xb);
    return out;
}

bool newton(const Ctx& c, int ya, int yb, double& a, double& b, Interval Ia, Interval Ib) {
    const BSpline& A = c.cs.yarns[ya].curve;
    const BSpline& B = c.cs.yarns[yb].curve;
    for (int it = 0; it < 50; ++it) {
        Vec3 F = A.eval(a) - B.eval(b);
        if (std::hypot(F.x, F.y) <= 1e-13 * c.al.Lref) return true;
        Vec3 da = A.d1(a), db = B.d1(b);
        double j00 = da.x, j01 = -db.x, j10 = da.y, j11 = -db.y;
        double det = j00 * j11 - j01 * j10;
        if (!(std::fabs(det) > 1e-300)) return false;
        double sa = -(j11 * F.x - j01 * F.y) / det;
        double sb = -(-j10 * F.x + j00 * F.y) / det;
        a = std::clamp(a + sa, Ia.lo, Ia.hi);
        b = std::clamp(b + sb, Ib.lo, Ib.hi);
    }
    Vec3 F = A.eval(a) - B.eval(b);
    return std::hypot(F.x, F.y) <= 1e-12 * c.al.Lref;
}

void addUnresolved(Ctx& c, int ya, double a0, double a1, int yb, double b0, double b1, const std::string& why) {
    // Merge with the previous item if contiguous (keeps reports short).
    if (!c.res.unresolved.empty()) {
        UnresolvedItem& p = c.res.unresolved.back();
        if (p.yarnA == ya && p.yarnB == yb && p.reason == why && std::fabs(p.uA1 - a0) < 1e-3 + (a1 - a0)) {
            p.uA0 = std::min(p.uA0, a0); p.uA1 = std::max(p.uA1, a1);
            p.uB0 = std::min(p.uB0, b0); p.uB1 = std::max(p.uB1, b1);
            return;
        }
    }
    c.res.unresolved.push_back({ya, yb, a0, a1, b0, b1, why});
}

Interval depthOver(const Ctx& c, int ya, int yb, Interval Xa, Interval Xb) {
    Box3 A = enclosePosition(c.cs.yarns[ya].curve, Xa.lo, Xa.hi);
    Box3 B = enclosePosition(c.cs.yarns[yb].curve, Xb.lo, Xb.hi);
    return pad(A.z - B.z, 2 * c.al.etaX);
}

// Verify a Newton candidate; on success append/merge a Crossing. Returns false if
// the root could not be certified (caller decides to subdivide or report).
bool verifyAndRecord(Ctx& c, int ya, int yb, double a, double b) {
    const double ua = c.cs.yarns[ya].uEnd(), ub = c.cs.yarns[yb].uEnd();
    for (double r : {1e-8, 1e-7, 1e-6, 1e-5, 1e-4}) {
        Interval Xa = intersect({a - r, a + r}, {0.0, ua});
        Interval Xb = intersect({b - r, b + r}, {0.0, ub});
        if (ya == yb && Xa.overlaps(Xb)) return false;  // diagonal: not a crossing
        Krawczyk2 k = krawczyk(c, ya, yb, a, b, Xa, Xb);
        if (!k.ok) continue;
        // Contract until the depth sign is resolved.
        Interval Ka = k.Ka, Kb = k.Kb;
        Interval depth = depthOver(c, ya, yb, Ka, Kb);
        for (int it = 0; it < 40 && depth.containsZero(); ++it) {
            Krawczyk2 k2 = krawczyk(c, ya, yb, a, b, Ka, Kb);
            Interval na = intersect(k2.Ka, Ka), nb = intersect(k2.Kb, Kb);
            if (na.width() >= Ka.width() && nb.width() >= Kb.width()) break;  // no progress
            Ka = na; Kb = nb;
            depth = depthOver(c, ya, yb, Ka, Kb);
        }
        Crossing x;
        x.yarnA = ya; x.yarnB = yb; x.boxA = Ka; x.boxB = Kb; x.uA = a; x.uB = b;
        x.pA = c.cs.yarns[ya].curve.eval(a);
        x.pB = c.cs.yarns[yb].curve.eval(b);
        x.depth = depth;
        if (depth.lo > 0) x.front = c.cs.yarns[ya].id;
        else if (depth.hi < 0) x.front = c.cs.yarns[yb].id;
        else { x.front = -1; x.status = Status::Unresolved; x.note = "depth sign unresolved"; }
        x.stitchA = c.cs.yarns[ya].stitchContaining(Ka.lo, Ka.hi);
        x.stitchB = c.cs.yarns[yb].stitchContaining(Kb.lo, Kb.hi);
        if (x.stitchA < 0 || x.stitchB < 0) {
            x.status = Status::Unresolved;
            x.note += x.note.empty() ? "on stitch boundary" : "; on stitch boundary";
        }
        // Deduplicate by root enclosure.
        for (Crossing& o : c.res.crossings) {
            if (o.yarnA != ya || o.yarnB != yb) continue;
            if (!o.boxA.overlaps(x.boxA) || !o.boxB.overlaps(x.boxB)) continue;
            Interval Ha = hull(o.boxA, x.boxA), Hb = hull(o.boxB, x.boxB);
            Ha = pad(Ha, 0.5 * Ha.width() + 1e-12); Hb = pad(Hb, 0.5 * Hb.width() + 1e-12);
            Ha = intersect(Ha, {0.0, ua}); Hb = intersect(Hb, {0.0, ub});
            Krawczyk2 kh = krawczyk(c, ya, yb, o.uA, o.uB, Ha, Hb);
            if (kh.ok) return true;  // same root
            o.status = Status::Unresolved;
            o.note += "; ambiguous duplicate";
            x.status = Status::Unresolved;
            x.note += "; ambiguous duplicate";
        }
        c.res.crossings.push_back(x);
        return true;
    }
    return false;
}

enum class Kind { Self, Adjacent, Distinct };

struct Piece { int y; double a0, a1; };

void recurse(Ctx& c, Piece A, Piece B, Kind kind, int depth) {
    if (--c.budget < 0) {
        if (c.budget == -1) addUnresolved(c, A.y, A.a0, A.a1, B.y, B.a0, B.a1, "work budget exhausted");
        return;
    }
    ++c.res.pairsVisited;
    const int maxDepth = c.opt.maxDepth;
    if (kind == Kind::Self) {
        double lo, hi;
        if (cone(c, A.y, A.a0, A.a1, lo, hi)) return;  // strictly monotone projection
        if (depth >= maxDepth || A.a1 - A.a0 < 1e-9) {
            addUnresolved(c, A.y, A.a0, A.a1, A.y, A.a0, A.a1, "projection of piece not certified injective");
            return;
        }
        double m = 0.5 * (A.a0 + A.a1);
        Piece L{A.y, A.a0, m}, R{A.y, m, A.a1};
        recurse(c, L, L, Kind::Self, depth + 1);
        recurse(c, R, R, Kind::Self, depth + 1);
        recurse(c, L, R, Kind::Adjacent, depth + 1);
        return;
    }
    if (kind == Kind::Adjacent) {
        if (unionCone(c, A.y, A.a0, A.a1, B.a0, B.a1)) return;  // meet only at the join
        if (depth >= maxDepth || A.a1 - A.a0 < 1e-9) {
            addUnresolved(c, A.y, A.a0, A.a1, B.y, B.a0, B.a1, "adjacent pieces not certified to meet only at join");
            return;
        }
        double ma = 0.5 * (A.a0 + A.a1), mb = 0.5 * (B.a0 + B.a1);
        Piece A1{A.y, A.a0, ma}, A2{A.y, ma, A.a1}, B1{B.y, B.a0, mb}, B2{B.y, mb, B.a1};
        recurse(c, A2, B1, Kind::Adjacent, depth + 1);
        recurse(c, A1, B1, Kind::Distinct, depth + 1);
        recurse(c, A1, B2, Kind::Distinct, depth + 1);
        recurse(c, A2, B2, Kind::Distinct, depth + 1);
        return;
    }
    // Distinct pieces.
    if (!overlaps(projBox(c, A.y, A.a0, A.a1), projBox(c, B.y, B.a0, B.a1))) return;
    double la, ha, lb, hb;
    bool transversal = cone(c, A.y, A.a0, A.a1, la, ha) && cone(c, B.y, B.a0, B.a1, lb, hb) &&
                       disjointModPi(la, ha, lb, hb);
    if (transversal) {
        // Expanded pieces: if still transversal, they contain at most one root, so a
        // certified root inside the expanded box settles this pair.
        double ma = 0.25 * (A.a1 - A.a0), mb = 0.25 * (B.a1 - B.a0);
        Interval Ea = intersect({A.a0 - ma, A.a1 + ma}, {0.0, c.cs.yarns[A.y].uEnd()});
        Interval Eb = intersect({B.a0 - mb, B.a1 + mb}, {0.0, c.cs.yarns[B.y].uEnd()});
        bool sameYarnOverlap = A.y == B.y && Ea.overlaps(Eb);
        double l2a, h2a, l2b, h2b;
        bool expTransversal = !sameYarnOverlap && cone(c, A.y, Ea.lo, Ea.hi, l2a, h2a) &&
                              cone(c, B.y, Eb.lo, Eb.hi, l2b, h2b) && disjointModPi(l2a, h2a, l2b, h2b);
        if (expTransversal) {
            double a = 0.5 * (A.a0 + A.a1), b = 0.5 * (B.a0 + B.a1);
            if (newton(c, A.y, B.y, a, b, Ea, Eb)) {
                if (verifyAndRecord(c, A.y, B.y, a, b)) return;
            }
        }
        if (depth >= maxDepth) {
            addUnresolved(c, A.y, A.a0, A.a1, B.y, B.a0, B.a1, "root not isolated at depth limit");
            return;
        }
    } else if (depth >= maxDepth || std::max(A.a1 - A.a0, B.a1 - B.a0) < c.minNonTransversalWidth) {
        addUnresolved(c, A.y, A.a0, A.a1, B.y, B.a0, B.a1, "near-tangent or overlapping projection");
        return;
    }
    double ma = 0.5 * (A.a0 + A.a1), mb = 0.5 * (B.a0 + B.a1);
    Piece A1{A.y, A.a0, ma}, A2{A.y, ma, A.a1}, B1{B.y, B.a0, mb}, B2{B.y, mb, B.a1};
    recurse(c, A1, B1, Kind::Distinct, depth + 1);
    recurse(c, A1, B2, Kind::Distinct, depth + 1);
    recurse(c, A2, B1, Kind::Distinct, depth + 1);
    recurse(c, A2, B2, Kind::Distinct, depth + 1);
}

} // namespace

CrossingResult findCrossings(const CurveSet& cs, const Allowances& al, const CrossingOptions& opt) {
    Ctx c{cs, al, opt, {}};
    const int Y = static_cast<int>(cs.yarns.size());
    // Span-level boxes for a cheap first rejection.
    std::vector<std::vector<Box2>> boxes(Y);
    for (int y = 0; y < Y; ++y)
        for (int k = 0; k < cs.yarns[y].curve.numSpans(); ++k) boxes[y].push_back(projBox(c, y, k, k + 1));
    for (int ya = 0; ya < Y; ++ya) {
        const int na = cs.yarns[ya].curve.numSpans();
        for (int yb = ya; yb < Y; ++yb) {
            const int nb = cs.yarns[yb].curve.numSpans();
            for (int ka = 0; ka < na; ++ka) {
                for (int kb = (ya == yb ? ka : 0); kb < nb; ++kb) {
                    Piece A{ya, double(ka), double(ka + 1)}, B{yb, double(kb), double(kb + 1)};
                    if (ya == yb && kb == ka) { recurse(c, A, A, Kind::Self, 0); continue; }
                    if (ya == yb && kb == ka + 1) { recurse(c, A, B, Kind::Adjacent, 0); continue; }
                    if (!overlaps(boxes[ya][ka], boxes[yb][kb])) continue;
                    recurse(c, A, B, Kind::Distinct, 0);
                }
            }
        }
    }
    return c.res;
}

} // namespace knit
