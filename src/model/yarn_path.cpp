#include "model/yarn_path.h"

#include <map>

#include <algorithm>
#include <cmath>

namespace knit {

namespace {
// QUADPACK qk15 nodes/weights.
constexpr double xgk[8] = {0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
                           0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
                           0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
                           0.207784955007898467600689403773245, 0.000000000000000000000000000000000};
constexpr double wgk[8] = {0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
                           0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
                           0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
                           0.204432940075298892414161999234649, 0.209482141084727828012999174891714};
constexpr double wg[4] = {0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
                          0.381830050505118944950369775488975, 0.417959183673469387755102040816327};

// One G7/K15 application on [a,b]; returns the Kronrod value and |K - G|.
double gk15(const BSpline& c, double a, double b, double& err) {
    double h = 0.5 * (b - a), m = 0.5 * (a + b);
    double fc = norm(c.d1(m));
    double K = fc * wgk[7], G = fc * wg[3];
    for (int j = 0; j < 7; ++j) {
        double f1 = norm(c.d1(m - h * xgk[j])), f2 = norm(c.d1(m + h * xgk[j]));
        K += wgk[j] * (f1 + f2);
        if (j % 2 == 1) G += wg[j / 2] * (f1 + f2);
    }
    err = std::fabs((K - G) * h);
    return K * h;
}

// A finished piece of the adaptive integration: it ends at u.
struct Leaf {
    double u, length, error;
};

// Bisects [a, b] until each piece meets its share of tol; appends the pieces in order.
void integrateAdaptive(const BSpline& c, double a, double b, double tol, int depth, std::vector<Leaf>& leaves) {
    double err;
    double v = gk15(c, a, b, err);
    if (err <= tol || depth >= 30) {
        leaves.push_back({b, v, err});
        return;
    }
    double m = 0.5 * (a + b);
    integrateAdaptive(c, a, m, 0.5 * tol, depth + 1, leaves);
    integrateAdaptive(c, m, b, 0.5 * tol, depth + 1, leaves);
}
} // namespace

void YarnPath::finalize(double tolS, const ExecPolicy& exec) {
    Tables t;
    t.tu.assign(1, 0.0);
    t.ts.assign(1, 0.0);
    int n = curve.numSpans();
    double errSum = 0, s = 0;
    // Per span (in parallel): its integration leaves and its whole-span bounds; then the sums
    // in span and leaf order, as one serial pass accumulates them.
    struct Span {
        std::vector<Leaf> leaves;
        ArcLengthBounds bounds;
    };
    std::vector<Span> spans(static_cast<size_t>(std::max(n, 0)));
    parallelFor(n, exec.chunk > 0 ? exec.chunk : 128, [&](int64_t kBegin, int64_t kEnd, int64_t) {
        for (int k = static_cast<int>(kBegin); k < static_cast<int>(kEnd); ++k) {
            Span& sp = spans[k];
            integrateAdaptive(curve, k, k + 1, tolS / std::max(1, n), 0, sp.leaves);
            sp.bounds = spanBounds(k, 0.0, 1.0);
        }
    }, exec);
    for (const Span& sp : spans)
        for (const Leaf& leaf : sp.leaves) {
            s += leaf.length;
            t.tu.push_back(leaf.u);
            t.ts.push_back(s);
            errSum += leaf.error;
        }
    t.length = s;
    t.lengthErr = errSum;
    // Prefix sums of whole-span bounds (Kahan-compensated).
    t.preLo.assign(n + 1, 0.0);
    t.preHi.assign(n + 1, 0.0);
    double cLo = 0, cHi = 0;
    for (int k = 0; k < n; ++k) {
        const double lo = spans[k].bounds.lo, hi = spans[k].bounds.hi;
        double yl = lo - cLo, tl = t.preLo[k] + yl;
        cLo = (tl - t.preLo[k]) - yl;
        t.preLo[k + 1] = tl;
        double yh = hi - cHi, th = t.preHi[k] + yh;
        cHi = (th - t.preHi[k]) - yh;
        t.preHi[k + 1] = th;
    }
    t.stitchesSorted = true;
    for (size_t i = 0; i < stitches.size(); ++i) {
        if (stitches[i].uBegin >= stitches[i].uEnd) t.stitchesSorted = false;
        if (i > 0 && stitches[i].uBegin < stitches[i - 1].uEnd) t.stitchesSorted = false;
    }
    t.sortedCount = stitches.size();
    tables_ = std::move(t);
}

ArcLengthBounds YarnPath::spanBounds(int k, double t0, double t1) const {
    CubicBezier piece = curve.span(k).sub(t0, t1);
    ArcLengthBounds b;
    b.hi = piece.polygonLength();
    Vec3 prev = piece.b[0];
    constexpr int chords = 8;
    for (int j = 1; j <= chords; ++j) {
        Vec3 p = piece.eval(double(j) / chords);
        b.lo += norm(p - prev);
        prev = p;
    }
    return b;
}

double YarnPath::arcLengthAt(double u) const {
    u = std::clamp(u, uBegin(), uEnd());
    const std::vector<double>& tu = tables_.tu;
    auto it = std::upper_bound(tu.begin(), tu.end(), u);
    size_t i = it == tu.begin() ? 0 : static_cast<size_t>(it - tu.begin()) - 1;
    if (i >= tu.size() - 1) return tables_.length;
    double err;
    return tables_.ts[i] + (u > tu[i] ? gk15(curve, tu[i], u, err) : 0.0);
}

double YarnPath::paramAtArcLength(double s) const {
    const std::vector<double>&tu = tables_.tu, &ts = tables_.ts;
    const double length = tables_.length;
    if (s <= 0) return uBegin();
    if (s >= length) return uEnd();
    auto it = std::upper_bound(ts.begin(), ts.end(), s);
    size_t i = static_cast<size_t>(it - ts.begin()) - 1;
    i = std::min(i, ts.size() - 2);
    double a = tu[i], b = tu[i + 1];
    double lo = a, hi = b;
    double u = a + (b - a) * (s - ts[i]) / std::max(1e-300, ts[i + 1] - ts[i]);
    for (int it2 = 0; it2 < 60; ++it2) {
        double err;
        double f = ts[i] + (u > a ? gk15(curve, a, u, err) : 0.0) - s;
        if (std::fabs(f) <= 1e-12 * std::max(1.0, length)) break;
        if (f > 0) hi = u; else lo = u;
        double sp = norm(curve.d1(u));
        double un = sp > 0 ? u - f / sp : 0.5 * (lo + hi);
        if (!(un > lo && un < hi)) un = 0.5 * (lo + hi);  // safeguard: stay in the bracket
        u = un;
    }
    return u;
}

ArcLengthBounds YarnPath::lengthBounds(double u0, double u1) const {
    ArcLengthBounds b;
    if (u1 <= u0) return b;
    u0 = std::max(u0, uBegin());
    u1 = std::min(u1, uEnd());
    const int n = curve.numSpans();
    int k0 = std::min(n - 1, static_cast<int>(std::floor(u0)));
    // Partial end spans exactly; whole spans in between from the prefix sums. The
    // difference of two prefix values is accurate to a few ulps of the PREFIX (not of
    // the difference): the allowance is relative to the accumulated magnitude.
    // Last span with a nonempty part of [u0, u1]: t1 = u1 - k1 in (0, 1].
    int k1 = std::clamp(static_cast<int>(std::ceil(u1)) - 1, k0, n - 1);
    if (k0 == k1) {
        b = spanBounds(k0, u0 - k0, u1 - k0);
        b.lo *= (1.0 - 1e-12);
        b.hi *= (1.0 + 1e-12);
        return b;
    }
    const ArcLengthBounds first = spanBounds(k0, u0 - k0, 1.0);
    b.lo += first.lo;
    b.hi += first.hi;
    const ArcLengthBounds last = spanBounds(k1, 0.0, u1 - k1);
    b.lo += last.lo;
    b.hi += last.hi;
    b.lo += tables_.preLo[k1] - tables_.preLo[k0 + 1];
    b.hi += tables_.preHi[k1] - tables_.preHi[k0 + 1];
    b.lo = b.lo * (1.0 - 1e-12) - 1e-13 * tables_.preLo[k1];
    b.hi = b.hi * (1.0 + 1e-12) + 1e-13 * tables_.preHi[k1];
    b.lo = std::max(0.0, b.lo);
    return b;
}

int YarnPath::stitchContaining(double u0, double u1) const {
    auto contains = [&](size_t i) {
        const StitchSpan& st = stitches[i];
        bool last = st.uEnd >= uEnd();
        return u0 >= st.uBegin && (u1 < st.uEnd || (last && u1 <= st.uEnd));
    };
    if (tables_.stitchesSorted && tables_.sortedCount == stitches.size()) {
        // Last stitch with uBegin <= u0; disjoint ranges, so it is the only candidate.
        auto it = std::upper_bound(stitches.begin(), stitches.end(), u0,
                                   [](double u, const StitchSpan& st) { return u < st.uBegin; });
        if (it == stitches.begin()) return -1;
        size_t i = static_cast<size_t>(it - stitches.begin()) - 1;
        return contains(i) ? static_cast<int>(i) : -1;
    }
    for (size_t i = 0; i < stitches.size(); ++i)
        if (contains(i)) return static_cast<int>(i);
    return -1;
}

Polyline samplePath(const YarnPath& y, double eps, int maxPerSpan, const ExecPolicy& exec) {
    Polyline pl;
    pl.yarnId = y.id;
    eps = std::max(eps, 1e-9);
    const int nSpans = y.curve.numSpans();
    pl.p.push_back(y.curve.eval(0.0));
    pl.u.push_back(0.0);
    // Per span (in parallel), its samples; appended in span order.
    struct Span {
        std::vector<Vec3> p;
        std::vector<double> u, segLo, segHi;
        double segEps = 0;
    };
    std::vector<Span> spans(static_cast<size_t>(std::max(nSpans, 0)));
    parallelFor(nSpans, exec.chunk > 0 ? exec.chunk : 256, [&](int64_t kBegin, int64_t kEnd, int64_t) {
        for (int k = static_cast<int>(kBegin); k < static_cast<int>(kEnd); ++k) {
            Span& sp = spans[k];
            CubicBezier bz = y.curve.span(k);
            double M = std::max(norm(bz.d2At0()), norm(bz.d2At1()));
            int n = std::max(1, static_cast<int>(std::ceil(std::sqrt(M / (8.0 * eps)))));
            n = std::min(n, maxPerSpan);
            sp.segEps = M / (8.0 * n * n);
            for (int j = 1; j <= n; ++j) {
                double t0 = double(j - 1) / n, t1 = double(j) / n;
                CubicBezier piece = bz.sub(t0, t1);
                sp.p.push_back(piece.b[3]);
                sp.u.push_back(k + t1);
                sp.segLo.push_back(norm(piece.b[3] - piece.b[0]) * (1.0 - 1e-12));
                sp.segHi.push_back(piece.polygonLength() * (1.0 + 1e-12));
            }
        }
    }, exec);
    for (const Span& sp : spans) {
        pl.epsBound = std::max(pl.epsBound, sp.segEps);
        pl.p.insert(pl.p.end(), sp.p.begin(), sp.p.end());
        pl.u.insert(pl.u.end(), sp.u.begin(), sp.u.end());
        pl.segLo.insert(pl.segLo.end(), sp.segLo.begin(), sp.segLo.end());
        pl.segHi.insert(pl.segHi.end(), sp.segHi.begin(), sp.segHi.end());
        pl.segEps.insert(pl.segEps.end(), sp.p.size(), sp.segEps);
    }
    pl.cumLo.assign(pl.segLo.size() + 1, 0.0);
    pl.cumHi.assign(pl.segHi.size() + 1, 0.0);
    for (size_t i = 0; i < pl.segLo.size(); ++i) {
        pl.cumLo[i + 1] = pl.cumLo[i] + pl.segLo[i];
        pl.cumHi[i + 1] = pl.cumHi[i] + pl.segHi[i];
    }
    pl.s.resize(pl.u.size());
    parallelFor(static_cast<int64_t>(pl.u.size()), exec.chunk > 0 ? exec.chunk : 4096, [&](int64_t b, int64_t e, int64_t) {
        for (int64_t i = b; i < e; ++i) pl.s[i] = y.arcLengthAt(pl.u[i]);
    }, exec);
    return pl;
}

double CurveSet::extent() const {
    double e = 0;
    for (const YarnPath& y : yarns)
        for (const Vec3& p : y.curve.ctrl) e = std::max({e, std::fabs(p.x), std::fabs(p.y), std::fabs(p.z)});
    return std::max(1.0, 2.0 * e);
}

void fillFlatTopology(CurveSet& cs) {
    struct Item { int course, wale, yarn, span; };
    std::vector<Item> items;
    for (size_t y = 0; y < cs.yarns.size(); ++y)
        for (size_t k = 0; k < cs.yarns[y].stitches.size(); ++k) {
            const StitchSpan& st = cs.yarns[y].stitches[k];
            items.push_back({st.course, st.wale, static_cast<int>(y), static_cast<int>(k)});
        }
    std::stable_sort(items.begin(), items.end(), [](const Item& a, const Item& b) {
        return a.course != b.course ? a.course < b.course : a.wale < b.wale;
    });
    cs.topo = KnitTopology{};
    std::map<std::pair<int, int>, int> at;  // (course, wale) -> id
    std::vector<double> ends;
    for (const YarnPath& y : cs.yarns) ends.push_back(y.uEnd());
    for (const Item& it : items) {
        StitchSpan& st = cs.yarns[it.yarn].stitches[it.span];
        StitchRecord r;
        r.id = st.id = static_cast<int>(cs.topo.stitches.size());
        r.yarn = it.yarn;
        r.uBegin = st.uBegin;
        r.uEnd = st.uEnd;
        r.turn = it.course;
        r.index = it.wale;
        auto p = at.find({it.course - 1, it.wale});
        if (p != at.end()) { r.parents[0] = p->second; r.nParents = 1; }
        at[{it.course, it.wale}] = r.id;
        cs.topo.stitches.push_back(r);
    }
    cs.topo.buildIndex(ends);
}

double YarnPath::thicknessAt(double u) const {
    if (thickness.empty()) return 1;
    const double c = std::clamp(u + 1, 0.0, static_cast<double>(thickness.size() - 1));  // control c sits at the knot u = c - 1
    const size_t k = static_cast<size_t>(c);
    const double f = c - static_cast<double>(k);
    return k + 1 < thickness.size() ? thickness[k] * (1 - f) + thickness[k + 1] * f : thickness[k];
}

double YarnPath::largestThickness() const {
    return thickness.empty() ? 1.0 : static_cast<double>(*std::max_element(thickness.begin(), thickness.end()));
}

} // namespace knit
