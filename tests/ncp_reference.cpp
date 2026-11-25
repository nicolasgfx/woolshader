// The CPU reference of the knitting problem (ncp_reference.h).
#include "ncp_reference.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>

#include "knitting/checks/checks.h"
#include "knitting/layout/uv_atlas.h"
#include "model/spatial_hash.h"

namespace knit {

namespace {

// Uniform cubic B-spline basis on a span, t in [0, 1] (settle.cpp basisWeights), and its first
// and second derivatives in t.
void basis(double t, double w[4], double d1[4], double d2[4]) {
    const double s = 1 - t;
    w[0] = s * s * s / 6;
    w[1] = (3 * t * t * t - 6 * t * t + 4) / 6;
    w[2] = (-3 * t * t * t + 3 * t * t + 3 * t + 1) / 6;
    w[3] = t * t * t / 6;
    d1[0] = -s * s / 2;
    d1[1] = (3 * t * t - 4 * t) / 2;
    d1[2] = (-3 * t * t + 2 * t + 1) / 2;
    d1[3] = t * t / 2;
    d2[0] = s;
    d2[1] = 3 * t - 2;
    d2[2] = -3 * t + 1;
    d2[3] = t;
}

uint64_t rowKey(RowFamily f, uint64_t a, uint64_t b) { return (uint64_t(f) << 60) ^ (a << 30) ^ b; }

void addEntry(Row& r, int var, const Vec3& g) {
    for (int k = 0; k < r.n; ++k)
        if (r.e[k].var == var) {
            r.e[k].g += g;
            return;
        }
    if (r.n >= kRowMaxEntries) throw std::runtime_error("knit problem: row with more than 10 variables");
    r.e[r.n++] = {var, g};
}

} // namespace

KnitProblem::KnitProblem(const CurveSet& cs, const Controls& P, const Controls& Prest, const ControlMask& freeCtrl,
                         const BaseSurface* surface, const SurfaceBinding* bind, const RelaxMoveParams& mp, const KnitProblemParams& kp)
    : cs_(cs), mp_(mp), kp_(kp), P0_(P), Xref_(P), free_(freeCtrl), surface_(surface) {
    (void)bind;
    // Normalization length of the contact and edge rows: the contact distance thr = D + c + margin.
    D_ = mp.thr > 0 ? mp.thr : 1.0;
    const int Y = static_cast<int>(P.size());
    varOf_.resize(Y);
    rest_.resize(Y);
    for (int y = 0; y < Y; ++y) {
        varOf_[y].assign(P[y].size(), -1);
        for (size_t c = 0; c < P[y].size(); ++c)
            if (freeCtrl[y][c]) {
                varOf_[y][c] = static_cast<int>(ctrlOf_.size());
                ctrlOf_.push_back({y, static_cast<int>(c)});
            }
        for (size_t c = 0; c + 1 < P[y].size(); ++c) rest_[y].push_back(norm(Prest[y][c + 1] - Prest[y][c]));
    }
    // Samples and active flags exactly as relaxFree builds them.
    double w4[4], d1[4], d2[4];
    for (int y = 0; y < Y; ++y) {
        const YarnPath& yp = cs.yarns[y];
        for (int k = 0; k < yp.curve.numSpans(); ++k) {
            bool act = false;
            for (int i = 0; i < 4; ++i) act = act || freeCtrl[y][k + i];
            int n = mp.samplesPerSpan;
            if (kp.maxSampleSpacing > 0)
                n = std::max(n, static_cast<int>(std::ceil((yp.arcLengthAt(k + 1.0) - yp.arcLengthAt(k)) / kp.maxSampleSpacing)));
            for (int j = 0; j < n; ++j) {
                Sample sm;
                sm.y = y;
                sm.k = k;
                sm.t = double(j) / n;
                basis(sm.t, w4, d1, d2);
                for (int i = 0; i < 4; ++i) sm.w[i] = w4[i];
                sm.s = yp.arcLengthAt(k + sm.t);
                samples_.push_back(sm);
                active_.push_back(act);
            }
        }
    }
}

Vec3 KnitProblem::evalSample(const Controls& P, const Sample& sm) const {
    const auto& c = P[sm.y];
    return c[sm.k] * sm.w[0] + c[sm.k + 1] * sm.w[1] + c[sm.k + 2] * sm.w[2] + c[sm.k + 3] * sm.w[3];
}

void KnitProblem::refresh(const Controls& P, SurfaceBinding* bind) {
    ++refreshes_;
    Xref_ = P;
    const double pad = mp_.thr + 2 * kp_.delta;
    std::vector<Vec3> pos(samples_.size());
    for (size_t i = 0; i < samples_.size(); ++i) pos[i] = evalSample(P, samples_[i]);
    auto segOk = [&](size_t i) { return i + 1 < samples_.size() && samples_[i + 1].y == samples_[i].y; };
    auto box = [&](size_t i, Vec3& lo, Vec3& hi) {
        lo = {std::min(pos[i].x, pos[i + 1].x), std::min(pos[i].y, pos[i + 1].y), std::min(pos[i].z, pos[i + 1].z)};
        hi = {std::max(pos[i].x, pos[i + 1].x), std::max(pos[i].y, pos[i + 1].y), std::max(pos[i].z, pos[i + 1].z)};
    };
    SpatialHash3 h(pad + 1.0);
    std::vector<size_t> segId;
    for (size_t i = 0; i < samples_.size(); ++i)
        if (segOk(i)) {
            Vec3 lo, hi;
            box(i, lo, hi);
            h.insert(lo, hi);
            segId.push_back(i);
        }
    h.finalize();
    cand_.clear();
    const Vec3 pv{pad, pad, pad};
    for (size_t i = 0; i < samples_.size(); ++i) {
        if (!segOk(i) || !active_[i]) continue;
        Vec3 lo, hi;
        box(i, lo, hi);
        h.query(lo - pv, hi + pv, [&](int id) {
            const size_t j = segId[id];
            if (j == i || (active_[j] && j < i)) return;  // each active pair once (relaxFree's rule)
            // Same yarn: the arc gap between the segments [s_i, s_i+1] and [s_j, s_j+1] (header:
            // CHANGE, stop point 1; adjacent segments have gap 0).
            if (samples_[j].y == samples_[i].y) {
                const double gap = j > i ? samples_[j].s - samples_[i + 1].s : samples_[i].s - samples_[j + 1].s;
                if (gap < mp_.window) return;
            }
            cand_.push_back({i, j});
        });
    }
    // Contact rows of the candidates by the contact model (KnitProblemParams::ContactModel); a
    // point-segment candidate gives two rows (point i - segment j, point j - segment i).
    {
        using CM = KnitProblemParams::ContactModel;
        std::vector<std::pair<size_t, size_t>> rowsC;
        std::vector<unsigned char> types;
        for (const auto& [i, j] : cand_) {
            int ty = kp_.contactModel == CM::Segment ? 0 : kp_.contactModel == CM::Point ? 1 : 2;
            if (kp_.contactModel == CM::Hybrid) {
                const Vec3 a = pos[i + 1] - pos[i], b = pos[j + 1] - pos[j];
                const double la = norm(a), lb = norm(b);
                ty = la > 0 && lb > 0 && std::fabs(dot(a, b)) > kp_.parallelCos * la * lb ? 2 : 0;
            }
            rowsC.push_back({i, j});
            types.push_back(static_cast<unsigned char>(ty));
            if (ty == 2) {
                rowsC.push_back({j, i});
                types.push_back(2);
            }
        }
        cand_.swap(rowsC);
        candType_.swap(types);
    }
    // Candidates touching each variable (block solvers).
    candOfVar_.assign(ctrlOf_.size(), {});
    for (size_t q = 0; q < cand_.size(); ++q) {
        int seen[kRowMaxEntries * 2];
        int ns = 0;
        for (size_t smp : {cand_[q].first, cand_[q].first + 1, cand_[q].second, cand_[q].second + 1}) {
            const Sample& s = samples_[smp];
            for (int m = 0; m < 4; ++m) {
                const int v = varOf_[s.y][s.k + m];
                if (v < 0 || std::find(seen, seen + ns, v) != seen + ns) continue;
                seen[ns++] = v;
                candOfVar_[v].push_back(static_cast<int>(q));
            }
        }
    }
    // Binding states continued to P (legacy refresh semantics: a state advances only on Owned).
    // These are the reference states the band rows continue from until the next refresh.
    bind_ = bind && bind->surface ? bind : nullptr;
    if (bind_)
        for (size_t v = 0; v < ctrlOf_.size(); ++v) {
            const auto [y, c] = ctrlOf_[v];
            ProjState next;
            const ProjResult r = bind->surface->continueFrom(bind->state[y][c], P[y][c], next, bind->cp);
            bind->result[y][c] = r;
            if (r == ProjResult::Owned) bind->state[y][c] = next;
        }
}

namespace {
// Signed distance from x to the foot, signed by the surface normal there, and its gradient (the
// unit direction from the foot; the normal itself at the foot).
void signedDistance(const Vec3& x, const Vec3& foot, const Vec3& n, double& z, Vec3& u) {
    const Vec3 d = x - foot;
    const double l = norm(d);
    const double s = dot(d, n) >= 0 ? 1.0 : -1.0;
    z = s * l;
    u = l > 1e-12 ? d * (s / l) : n;
}
} // namespace

bool KnitProblem::bandDistance(const Controls& P, int v, double& z, Vec3& u) const {
    const auto [y, c] = ctrlOf_[v];
    if (bind_) {
        ProjState next;
        if (bind_->surface->continueFrom(bind_->state[y][c], P[y][c], next, bind_->cp) != ProjResult::Owned) return false;
        signedDistance(P[y][c], next.foot, next.n, z, u);
        return true;
    }
    if (!surface_) return false;
    const SurfaceSample q = surface_->project(P[y][c]);
    signedDistance(P[y][c], q.p, q.n, z, u);
    return true;
}

double KnitProblem::motionSinceRefresh(const Controls& P) const {
    double m = 0;
    for (const auto& [y, c] : ctrlOf_) m = std::max(m, norm(P[y][c] - Xref_[y][c]));
    return m;
}

void KnitProblem::addSampleTerms(Row& r, const Sample& s, double f, const Vec3& n) const {
    for (int m = 0; m < 4; ++m) {
        const int v = varOf_[s.y][s.k + m];
        if (v >= 0 && s.w[m] != 0) addEntry(r, v, n * (f * s.w[m]));
    }
}

double KnitProblem::contactDistance(const Controls& P, size_t i, size_t j, int type, double& s, double& t) const {
    const Vec3 a0 = evalSample(P, samples_[i]), a1 = evalSample(P, samples_[i + 1]);
    const Vec3 b0 = evalSample(P, samples_[j]), b1 = evalSample(P, samples_[j + 1]);
    if (type == 1) {
        s = t = 0;
        return norm(a0 - b0);
    }
    const ClosestSegmentPoints closestPair = type == 2 ? closestSegmentPoints(a0, a0, b0, b1)  // s = 0
                                                : closestSegmentPoints(a0, a1, b0, b1);
    s = closestPair.s;
    t = closestPair.t;
    return closestPair.distance;
}

void KnitProblem::contactRows(const Controls& P, size_t i, size_t j, int type, Row& r, Row* target) const {
    const Vec3 a0 = evalSample(P, samples_[i]), a1 = evalSample(P, samples_[i + 1]);
    const Vec3 b0 = evalSample(P, samples_[j]), b1 = evalSample(P, samples_[j + 1]);
    double s, t;
    const double d = contactDistance(P, i, j, type, s, t);
    const Vec3 pa = lerp(a0, a1, s), pb = lerp(b0, b1, t);
    const Vec3 n = d > 1e-9 ? (pa - pb) / d : Vec3{0, 0, 1};
    r = Row{};
    r.family = RowFamily::Contact;
    r.value = (hardC() - d) / D_;
    r.key = rowKey(RowFamily::Contact, i, j) ^ (uint64_t(type) << 57);
    // g = (thr - n . (pa - pb)) / D with s, t, n frozen: dg/dx = -(1/D) n * (coefficient of x in pa - pb).
    const double f = -1.0 / D_;
    addSampleTerms(r, samples_[i], f * (1 - s), n);
    addSampleTerms(r, samples_[i + 1], f * s, n);
    addSampleTerms(r, samples_[j], -f * (1 - t), n);
    addSampleTerms(r, samples_[j + 1], -f * t, n);
    if (target) {
        *target = Row{};
        target->family = RowFamily::ContactTarget;
        target->key = rowKey(RowFamily::ContactTarget, i, j) ^ (uint64_t(type) << 57);
        if (targets() && d < mp_.thr) {  // residual max(0, thr - d) / D: same gradient as the hard row while active
            target->value = (mp_.thr - d) / D_;
            target->n = r.n;
            for (int e = 0; e < r.n; ++e) target->e[e] = r.e[e];
        }
    }
}

bool KnitProblem::curvatureRows(const Controls& P, int y, int k, int q, Row& r, Row* target, bool* degenerate) const {
    double w[4], b1[4], b2[4];
    basis(double(q) / (kp_.curvatureSamples - 1), w, b1, b2);
    Vec3 d1, d2;
    for (int i = 0; i < 4; ++i) {
        d1 += P[y][k + i] * b1[i];
        d2 += P[y][k + i] * b2[i];
    }
    const double L = norm(d1);
    *degenerate = L < kp_.kDegenerate;
    if (*degenerate) return false;
    const Vec3 c = cross(d1, d2);
    const double C = norm(c), L3 = L * L * L;
    const double kappa = C / L3;
    // dkappa/dd1 = (d2 x c) / (C L^3) - 3 C d1 / L^5; dkappa/dd2 = (c x d1) / (C L^3).
    Vec3 g1 = d1 * (-3 * C / (L3 * L * L)), g2;
    if (C > 1e-15) {
        g1 += cross(d2, c) / (C * L3);
        g2 = cross(c, d1) / (C * L3);
    }
    auto fill = [&](Row& row, RowFamily fam, double limit) {
        row = Row{};
        row.family = fam;
        row.value = kappa / limit - 1;
        row.key = rowKey(fam, uint64_t(y) << 20 | uint64_t(k), uint64_t(q));
        for (int i = 0; i < 4; ++i) {
            const int v = varOf_[y][k + i];
            if (v >= 0) addEntry(row, v, (g1 * b1[i] + g2 * b2[i]) / limit);
        }
    };
    fill(r, RowFamily::Curvature, hardK());
    if (target) {
        if (targets() && kappa > mp_.kappaMax) {
            fill(*target, RowFamily::CurvatureTarget, mp_.kappaMax);
        } else {
            *target = Row{};
            target->family = RowFamily::CurvatureTarget;
            target->key = rowKey(RowFamily::CurvatureTarget, uint64_t(y) << 20 | uint64_t(k), uint64_t(q));
        }
    }
    return true;
}

void KnitProblem::linearize(const Controls& P, std::vector<Row>& rows, double activeBelow) const {
    rows.clear();
    const int Y = static_cast<int>(P.size());
    // Edges.
    for (int y = 0; y < Y; ++y)
        for (size_t c = 0; c + 1 < P[y].size(); ++c) {
            const int va = varOf_[y][c], vb = varOf_[y][c + 1];
            if (va < 0 && vb < 0) continue;
            const Vec3 d = P[y][c + 1] - P[y][c];
            const double l = norm(d);
            Row r;
            r.family = RowFamily::Edge;
            r.value = (l - rest_[y][c]) / D_;
            r.key = rowKey(RowFamily::Edge, y, c);
            r.aux = l;
            if (l > 1e-12) {
                const Vec3 u = d / (l * D_);
                if (vb >= 0) addEntry(r, vb, u);
                if (va >= 0) addEntry(r, va, -u);
            }
            rows.push_back(r);
        }
    // Curvature: every span with a free support control.
    for (int y = 0; y < Y; ++y) {
        const int nSp = static_cast<int>(P[y].size()) - 3;
        for (int k = 0; k < nSp; ++k) {
            bool any = false;
            for (int i = 0; i < 4; ++i) any = any || varOf_[y][k + i] >= 0;
            if (!any) continue;
            for (int q = 0; q < kp_.curvatureSamples; ++q) {
                Row r, t;
                bool deg = false;
                if (!curvatureRows(P, y, k, q, r, &t, &deg)) continue;
                if (r.value > activeBelow) rows.push_back(r);
                if (t.n > 0) rows.push_back(t);
            }
        }
    }
    // Contacts.
    for (size_t q = 0; q < cand_.size(); ++q) {
        const auto [i, j] = cand_[q];
        Row r, t;
        contactRows(P, i, j, candType_[q], r, &t);
        if (r.value > activeBelow) rows.push_back(r);
        if (t.n > 0) rows.push_back(t);
    }
    // Band and cap.
    for (size_t v = 0; v < ctrlOf_.size(); ++v) {
        const auto [y, c] = ctrlOf_[v];
        anchorRows(P, static_cast<int>(v), rows);
        double z;
        Vec3 nrm;
        if (bandDistance(P, static_cast<int>(v), z, nrm)) {
            Row up, lo;
            up.family = RowFamily::BandUpper;
            up.value = z / mp_.bandHalf - 1;
            up.key = rowKey(RowFamily::BandUpper, y, c);
            addEntry(up, static_cast<int>(v), nrm / mp_.bandHalf);
            lo.family = RowFamily::BandLower;
            lo.value = -z / mp_.bandHalf - 1;
            lo.key = rowKey(RowFamily::BandLower, y, c);
            addEntry(lo, static_cast<int>(v), nrm * (-1.0 / mp_.bandHalf));
            if (up.value > activeBelow) rows.push_back(up);
            if (lo.value > activeBelow) rows.push_back(lo);
        }
        if (mp_.maxDisplacement > 0) {
            const Vec3 d = P[y][c] - P0_[y][c];
            const double cap2 = mp_.maxDisplacement * mp_.maxDisplacement;
            Row r;
            r.family = RowFamily::Cap;
            r.value = dot(d, d) / cap2 - 1;
            r.key = rowKey(RowFamily::Cap, y, c);
            addEntry(r, static_cast<int>(v), d * (2.0 / cap2));
            if (r.value > activeBelow) rows.push_back(r);
        }
    }
}

void KnitProblem::linearizeVar(const Controls& P, int v, std::vector<Row>& rows, double activeBelow) const {
    rows.clear();
    const auto [y, c] = ctrlOf_[v];
    const int nc = static_cast<int>(P[y].size());
    for (int e = c - 1; e <= c; ++e) {
        if (e < 0 || e + 1 >= nc) continue;
        const int va = varOf_[y][e], vb = varOf_[y][e + 1];
        const Vec3 d = P[y][e + 1] - P[y][e];
        const double l = norm(d);
        Row r;
        r.family = RowFamily::Edge;
        r.value = (l - rest_[y][e]) / D_;
        r.key = rowKey(RowFamily::Edge, y, e);
        r.aux = l;
        if (l > 1e-12) {
            const Vec3 u = d / (l * D_);
            if (vb >= 0) addEntry(r, vb, u);
            if (va >= 0) addEntry(r, va, -u);
        }
        rows.push_back(r);
    }
    for (int k = std::max(0, c - 3); k <= c && k + 3 < nc; ++k)
        for (int q = 0; q < kp_.curvatureSamples; ++q) {
            Row r, t;
            bool deg = false;
            if (!curvatureRows(P, y, k, q, r, &t, &deg)) continue;
            if (r.value > activeBelow) rows.push_back(r);
            if (t.n > 0) rows.push_back(t);
        }
    for (int q : candOfVar_[v]) {
        Row r, t;
        contactRows(P, cand_[q].first, cand_[q].second, candType_[q], r, &t);
        if (r.value > activeBelow) rows.push_back(r);
        if (t.n > 0) rows.push_back(t);
    }
    anchorRows(P, v, rows);
    double z;
    Vec3 nrm;
    if (bandDistance(P, v, z, nrm)) {
        Row up, lo;
        up.family = RowFamily::BandUpper;
        up.value = z / mp_.bandHalf - 1;
        up.key = rowKey(RowFamily::BandUpper, y, c);
        addEntry(up, v, nrm / mp_.bandHalf);
        lo.family = RowFamily::BandLower;
        lo.value = -z / mp_.bandHalf - 1;
        lo.key = rowKey(RowFamily::BandLower, y, c);
        addEntry(lo, v, nrm * (-1.0 / mp_.bandHalf));
        if (up.value > activeBelow) rows.push_back(up);
        if (lo.value > activeBelow) rows.push_back(lo);
    }
    if (mp_.maxDisplacement > 0) {
        const Vec3 d = P[y][c] - P0_[y][c];
        const double cap2 = mp_.maxDisplacement * mp_.maxDisplacement;
        Row r;
        r.family = RowFamily::Cap;
        r.value = dot(d, d) / cap2 - 1;
        r.key = rowKey(RowFamily::Cap, y, c);
        addEntry(r, v, d * (2.0 / cap2));
        if (r.value > activeBelow) rows.push_back(r);
    }
}

void KnitProblem::anchorRows(const Controls& P, int v, std::vector<Row>& rows) const {
    if (kp_.wA <= 0) return;
    const auto [y, c] = ctrlOf_[v];
    const Vec3 d = (P[y][c] - P0_[y][c]) / D_;
    const Vec3 unit[3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    for (int i = 0; i < 3; ++i) {
        Row r;
        r.family = RowFamily::Anchor;
        r.value = i == 0 ? d.x : i == 1 ? d.y : d.z;
        r.key = rowKey(RowFamily::Anchor, uint64_t(y) << 2 | uint64_t(i), uint64_t(c));
        addEntry(r, v, unit[i] / D_);
        rows.push_back(r);
    }
}

KnitProblemScore KnitProblem::scoreGeometry(const Controls& P, const SurfaceBinding* bind, double tol) const {
    KnitProblemScore sc;
    const int Y = static_cast<int>(P.size());
    auto viol = [&](RowFamily f, double g) {
        const int i = static_cast<int>(f);
        if (g > sc.maxViol[i]) sc.maxViol[i] = g;
        if (g > tol) ++sc.nViol[i];
    };
    long long edges = 0;
    double strain2 = 0;
    for (int y = 0; y < Y; ++y)
        for (size_t c = 0; c + 1 < P[y].size(); ++c) {
            if (varOf_[y][c] < 0 && varOf_[y][c + 1] < 0) continue;
            const double l = norm(P[y][c + 1] - P[y][c]);
            const double r = (l - rest_[y][c]) / D_;
            sc.objective += 0.5 * kp_.wL * r * r;
            sc.maxViol[0] = std::max(sc.maxViol[0], std::fabs(r));
            if (rest_[y][c] > 0) {
                const double st = std::fabs(l - rest_[y][c]) / rest_[y][c];
                sc.strainMax = std::max(sc.strainMax, st);
                strain2 += st * st;
                ++edges;
            }
        }
    sc.strainRms = edges ? std::sqrt(strain2 / edges) : 0;
    if (kp_.wA > 0)
        for (const auto& [y, c] : ctrlOf_) {
            const Vec3 d = (P[y][c] - P0_[y][c]) / D_;
            sc.objective += 0.5 * kp_.wA * dot(d, d);
        }
    for (int y = 0; y < Y; ++y) {
        const int nSp = static_cast<int>(P[y].size()) - 3;
        for (int k = 0; k < nSp; ++k) {
            bool any = false;
            for (int i = 0; i < 4; ++i) any = any || varOf_[y][k + i] >= 0;
            if (!any) continue;
            for (int q = 0; q < kp_.curvatureSamples; ++q) {
                Row r, t;
                bool deg = false;
                if (curvatureRows(P, y, k, q, r, &t, &deg)) {
                    viol(RowFamily::Curvature, r.value);
                    if (t.n > 0) {
                        viol(RowFamily::CurvatureTarget, t.value);
                        sc.objective += 0.5 * kp_.wT * t.value * t.value;
                    }
                } else if (deg) {
                    ++sc.degenerate;
                }
            }
        }
    }
    for (size_t q = 0; q < cand_.size(); ++q) {
        const auto [i, j] = cand_[q];
        double s, t;
        const double d = contactDistance(P, i, j, candType_[q], s, t);
        viol(RowFamily::Contact, (hardC() - d) / D_);
        if (targets() && d < mp_.thr) {
            const double r = (mp_.thr - d) / D_;
            viol(RowFamily::ContactTarget, r);
            sc.objective += 0.5 * kp_.wT * r * r;
        }
    }
    for (size_t v = 0; v < ctrlOf_.size(); ++v) {
        const auto [y, c] = ctrlOf_[v];
        Vec3 foot, n;
        bool on = false;
        if (bind && bind->surface) {
            ProjState next;
            if (bind->surface->continueFrom(bind->state[y][c], P[y][c], next, bind->cp) == ProjResult::Owned) {
                foot = next.foot;
                n = next.n;
                on = true;
            } else {
                ++sc.notOwned;
            }
        } else if (surface_) {
            const SurfaceSample q = surface_->project(P[y][c]);
            foot = q.p;
            n = q.n;
            on = true;
        }
        if (on) {
            double z;
            Vec3 u;
            signedDistance(P[y][c], foot, n, z, u);
            viol(RowFamily::BandUpper, z / mp_.bandHalf - 1);
            viol(RowFamily::BandLower, -z / mp_.bandHalf - 1);
        }
        if (mp_.maxDisplacement > 0) {
            const Vec3 d = P[y][c] - P0_[y][c];
            viol(RowFamily::Cap, dot(d, d) / (mp_.maxDisplacement * mp_.maxDisplacement) - 1);
        }
    }
    sc.coverage = motionSinceRefresh(P);
    sc.covered = sc.coverage <= kp_.delta * (1 + 1e-12);
    return sc;
}

bool KnitProblem::passThrough(Controls& P, int v, const Vec3& from, const Vec3& to, double near) const {
    const auto [vy, vc] = ctrlOf_[v];
    auto closest = [&](const Controls& P, size_t i, size_t j, Vec3& n) {
        const Vec3 a0 = evalSample(P, samples_[i]), a1 = evalSample(P, samples_[i + 1]);
        const Vec3 b0 = evalSample(P, samples_[j]), b1 = evalSample(P, samples_[j + 1]);
        const ClosestSegmentPoints closestPair = closestSegmentPoints(a0, a1, b0, b1);
        n = lerp(a0, a1, closestPair.s) - lerp(b0, b1, closestPair.t);
        return closestPair.distance;
    };
    bool through = false;
    for (int q : candOfVar_[v]) {
        const auto [i, j] = cand_[q];
        Vec3 n0, n1;
        P[vy][vc] = from;
        const double d0 = closest(P, i, j, n0);
        if (d0 >= near) continue;
        P[vy][vc] = to;
        const double d1 = closest(P, i, j, n1);
        if (d1 <= 1e-9 || (d0 > 1e-9 && dot(n0, n1) < 0)) {  // touching after the step, or sides switched
            through = true;
            break;
        }
    }
    P[vy][vc] = to;
    return through;
}

} // namespace knit
