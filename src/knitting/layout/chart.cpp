// The (a, b) chart of a sphere-like piece, from pole A to pole B (chart.h).
#include "knitting/layout/chart.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>

namespace knit {

Chart::Chart(const TriMesh& mesh, int rings) : mesh_(mesh), rings_(rings) {
    if (rings < 1) throw std::runtime_error("chart: rings must be >= 1");
    aMin_ = 1.0 / (rings + 1);
    aMax_ = double(rings) / (rings + 1);
    if (mesh_.uv.size() != mesh_.tri.size() || mesh_.tri.empty())
        throw std::runtime_error("chart: mesh has no per-corner chart coordinates");
    if (mesh_.vnormal.size() != mesh_.pos.size()) mesh_.computeNormals();
    const int T = static_cast<int>(mesh_.tri.size());
    // Poles: the vertices whose corners carry a = 0 / a = 1.
    for (int t = 0; t < T; ++t)
        for (int k = 0; k < 3; ++k) {
            if (mesh_.uv[t][k].x == 0.0) poleA_ = mesh_.tri[t][k];
            if (mesh_.uv[t][k].x == 1.0) poleB_ = mesh_.tri[t][k];
        }
    if (poleA_ < 0 || poleB_ < 0) throw std::runtime_error("chart: no pole corners (a = 0 / 1)");
    // The exported UVs are float32: snap the domain bounds to the stored ring values.
    double dA = 1e300, dB = 1e300, snapA = aMin_, snapB = aMax_;
    for (const auto& u : mesh_.uv)
        for (const Vec2& q : u) {
            if (std::fabs(q.x - aMin_) < dA) { dA = std::fabs(q.x - aMin_); snapA = q.x; }
            if (std::fabs(q.x - aMax_) < dB) { dB = std::fabs(q.x - aMax_); snapB = q.x; }
        }
    if (dA > 1e-6 || dB > 1e-6) throw std::runtime_error("chart: no corners on the first/last ring");
    aMin_ = snapA;
    aMax_ = snapB;

    // Per-triangle chart derivatives: solve [dp/da dp/db] * [da db]^T = dp on two edges.
    dpda_.resize(T);
    dpdb_.resize(T);
    b0_ = 1e300;
    b1_ = -1e300;
    for (int t = 0; t < T; ++t) {
        const std::array<int, 3>& f = mesh_.tri[t];
        const std::array<Vec2, 3>& u = mesh_.uv[t];
        Vec2 e1 = u[1] - u[0], e2 = u[2] - u[0];
        Vec3 p1 = mesh_.pos[f[1]] - mesh_.pos[f[0]], p2 = mesh_.pos[f[2]] - mesh_.pos[f[0]];
        double det = cross(e1, e2);
        if (det != 0) {
            dpda_[t] = (p1 * e2.y - p2 * e1.y) / det;
            dpdb_[t] = (p2 * e1.x - p1 * e2.x) / det;
        }
        for (int k = 0; k < 3; ++k) {
            b0_ = std::min(b0_, u[k].y);
            b1_ = std::max(b1_, u[k].y);
        }
    }
    buildTriangleGrid();
    buildTriangleHash();
}

// A 2D grid over the chart triangles in (a, b).
void Chart::buildTriangleGrid() {
    const int T = static_cast<int>(mesh_.tri.size());
    na_ = std::max(16, static_cast<int>(std::sqrt(double(T)) * 1.5));
    nb_ = na_;
    grid_.assign(size_t(na_) * nb_, {});
    const double sa = na_, sb = nb_ / (b1_ - b0_);
    for (int t = 0; t < T; ++t) {
        const std::array<Vec2, 3>& u = mesh_.uv[t];
        double alo = std::min({u[0].x, u[1].x, u[2].x}), ahi = std::max({u[0].x, u[1].x, u[2].x});
        double blo = std::min({u[0].y, u[1].y, u[2].y}), bhi = std::max({u[0].y, u[1].y, u[2].y});
        int i0 = std::clamp(static_cast<int>(std::floor(alo * sa)), 0, na_ - 1);
        int i1 = std::clamp(static_cast<int>(std::floor(ahi * sa)), 0, na_ - 1);
        int j0 = std::clamp(static_cast<int>(std::floor((blo - b0_) * sb)), 0, nb_ - 1);
        int j1 = std::clamp(static_cast<int>(std::floor((bhi - b0_) * sb)), 0, nb_ - 1);
        for (int i = i0; i <= i1; ++i)
            for (int j = j0; j <= j1; ++j) grid_[size_t(i) * nb_ + j].push_back(t);
    }
}

// A 3D hash over triangle boxes, cell ~ 2 mean edge lengths.
void Chart::buildTriangleHash() {
    const int T = static_cast<int>(mesh_.tri.size());
    double edge = 0;
    for (int t = 0; t < T; ++t) edge += norm(mesh_.pos[mesh_.tri[t][1]] - mesh_.pos[mesh_.tri[t][0]]);
    hashCell_ = std::max(1e-12, 2.0 * edge / T);
    hash_ = SpatialHash3(hashCell_);
    for (int t = 0; t < T; ++t) {
        const std::array<int, 3>& f = mesh_.tri[t];
        Vec3 lo = mesh_.pos[f[0]], hi = lo;
        for (int k = 1; k < 3; ++k) {
            const Vec3& q = mesh_.pos[f[k]];
            lo = {std::min(lo.x, q.x), std::min(lo.y, q.y), std::min(lo.z, q.z)};
            hi = {std::max(hi.x, q.x), std::max(hi.y, q.y), std::max(hi.z, q.z)};
        }
        hash_.insert(lo, hi);
    }
    hash_.finalize();
}

bool Chart::inTriangle(int t, double a, double b, double w[3], double tol) const {
    const std::array<Vec2, 3>& u = mesh_.uv[t];
    Vec2 e1 = u[1] - u[0], e2 = u[2] - u[0], q = Vec2{a, b} - u[0];
    double det = cross(e1, e2);
    if (det == 0) return false;
    w[1] = cross(q, e2) / det;
    w[2] = cross(e1, q) / det;
    w[0] = 1 - w[1] - w[2];
    return w[0] >= -tol && w[1] >= -tol && w[2] >= -tol;
}

ChartPoint Chart::at(int t, const double w[3]) const {
    ChartPoint s;
    const std::array<int, 3>& f = mesh_.tri[t];
    s.ok = true;
    s.tri = t;
    for (int k = 0; k < 3; ++k) s.bary[k] = w[k];
    s.p = mesh_.pos[f[0]] * w[0] + mesh_.pos[f[1]] * w[1] + mesh_.pos[f[2]] * w[2];
    s.n = normalized(mesh_.vnormal[f[0]] * w[0] + mesh_.vnormal[f[1]] * w[1] + mesh_.vnormal[f[2]] * w[2]);
    s.dpda = dpda_[t];
    s.dpdb = dpdb_[t];
    const std::array<Vec2, 3>& u = mesh_.uv[t];
    s.ab = u[0] * w[0] + u[1] * w[1] + u[2] * w[2];
    return s;
}

ChartPoint Chart::evalRaw(double a, double b) const {
    const double sb = nb_ / (b1_ - b0_);
    int i = std::clamp(static_cast<int>(std::floor(a * na_)), 0, na_ - 1);
    int j = std::clamp(static_cast<int>(std::floor((b - b0_) * sb)), 0, nb_ - 1);
    if (b < b0_ - 1e-12 || b > b1_ + 1e-12 || a < -1e-12 || a > 1 + 1e-12) return {};
    // The triangle containing (a, b), or on a shared edge the one with the largest
    // minimum weight (deterministic, and the seam side is decided by b itself).
    int best = -1;
    double bestMin = -1e300, wb[3] = {0, 0, 0};
    for (int t : grid_[size_t(i) * nb_ + j]) {
        double w[3];
        if (!inTriangle(t, a, b, w, 1e-12)) continue;
        double m = std::min({w[0], w[1], w[2]});
        if (m > bestMin) { bestMin = m; best = t; std::copy(w, w + 3, wb); }
    }
    if (best < 0) return {};
    return at(best, wb);
}

ChartPoint Chart::eval(double a, double b) const {
    double bw = b - std::floor(b);
    ChartPoint s = evalRaw(a, bw);
    if (!s.ok) s = evalRaw(a, bw + 1.0);
    return s;
}

ChartPoint Chart::closestPoint(const Vec3& p) const {
    ChartPoint best;
    double bestD = 1e300, bestW[3] = {0, 0, 0};
    int bestT = -1;
    for (double r = hashCell_; bestT < 0 || bestD > r; r *= 2) {
        Vec3 d{r, r, r};
        hash_.query(p - d, p + d, [&](int t) {
            const std::array<int, 3>& f = mesh_.tri[t];
            const TrianglePoint onTri = closestOnTriangle(p, mesh_.pos[f[0]], mesh_.pos[f[1]], mesh_.pos[f[2]]);
            double dd = norm(onTri.p - p);
            if (dd < bestD) { bestD = dd; bestT = t; std::copy(onTri.bary, onTri.bary + 3, bestW); }
        });
        if (r > 1e6 * hashCell_) break;
    }
    if (bestT < 0) return best;
    best = at(bestT, bestW);
    best.distance = bestD;
    return best;
}

double Chart::ringLength(double a) const {
    double L = 0;
    for (size_t t = 0; t < mesh_.tri.size(); ++t) {
        const std::array<int, 3>& f = mesh_.tri[t];
        const std::array<Vec2, 3>& u = mesh_.uv[t];
        Vec3 pts[3];
        int n = 0;
        for (int k = 0; k < 3; ++k) {
            int l = (k + 1) % 3;
            double ak = u[k].x - a, al = u[l].x - a;
            // Half-open rule: an edge is cut if its endpoints are on different sides,
            // with a == level counted as "above", so shared edges are counted once.
            if ((ak < 0) != (al < 0)) {
                double s = ak / (ak - al);
                pts[n++] = lerp(mesh_.pos[f[k]], mesh_.pos[f[l]], s);
            }
        }
        if (n == 2) L += norm(pts[1] - pts[0]);
    }
    return L;
}

Chart::RingAt Chart::longestRing() const {
    const int N = 200;
    int best = 0;
    std::vector<double> L(N + 1);
    for (int i = 0; i <= N; ++i) L[i] = ringLength(double(i) / N);
    for (int i = 1; i <= N; ++i)
        if (L[i] > L[best]) best = i;
    double lo = std::max(0, best - 1) / double(N), hi = std::min(N, best + 1) / double(N);
    const double g = 0.5 * (std::sqrt(5.0) - 1);
    for (int it = 0; it < 40; ++it) {
        double x1 = hi - g * (hi - lo), x2 = lo + g * (hi - lo);
        if (ringLength(x1) > ringLength(x2)) hi = x2; else lo = x1;
    }
    const double a = 0.5 * (lo + hi);
    return {a, ringLength(a)};
}

double Chart::meridianLength(double b, double a0, double a1, int n) const {
    double L = 0;
    if (a0 < aMin_) {
        L += norm(eval(aMin_, b).p - poleA());
        a0 = aMin_;
    }
    if (a1 > aMax_) {
        L += norm(eval(aMax_, b).p - poleB());
        a1 = aMax_;
    }
    Vec3 prev = eval(a0, b).p;
    for (int i = 1; i <= n; ++i) {
        Vec3 q = eval(a0 + (a1 - a0) * i / n, b).p;
        L += norm(q - prev);
        prev = q;
    }
    return L;
}

// Orientation in (a, b) and agreement with the outward 3D orientation.
void Chart::checkOrientation(ChartCheck& c) const {
    int pos2 = 0, neg2 = 0, pos3 = 0, neg3 = 0;
    for (int t = 0; t < static_cast<int>(mesh_.tri.size()); ++t) {
        const std::array<Vec2, 3>& u = mesh_.uv[t];
        double det = cross(u[1] - u[0], u[2] - u[0]);
        if (std::fabs(det) <= 1e-14) { ++c.degenerateTriangles; continue; }
        (det > 0 ? pos2 : neg2)++;
        double o = dot(cross(dpda_[t], dpdb_[t]), mesh_.faceNormal(t));
        (o > 0 ? pos3 : neg3)++;
    }
    c.flippedTriangles = std::min(pos2, neg2);
    c.flipped3D = std::min(pos3, neg3);
    c.orientation = (c.degenerateTriangles == 0 && c.flippedTriangles == 0 && c.flipped3D == 0) ? Status::Pass
                                                                                                 : Status::Fail;
}

// a per vertex (agreeing corners, no interior extremum, two poles), and b per vertex modulo 1.
void Chart::checkVertexValues(ChartCheck& c) const {
    const int T = static_cast<int>(mesh_.tri.size());
    const int V = static_cast<int>(mesh_.pos.size());
    std::vector<double> av(V, -1), bv(V, -1);
    std::vector<std::vector<int>> adj(V);
    for (int t = 0; t < T; ++t)
        for (int k = 0; k < 3; ++k) {
            int v = mesh_.tri[t][k];
            const Vec2& u = mesh_.uv[t][k];
            if (av[v] < 0) av[v] = u.x;
            else if (av[v] != u.x) ++c.aCornerMismatches;
            if (u.x >= aMin_ && u.x <= aMax_) {
                double bm = u.y - std::floor(u.y);
                if (bv[v] < 0) bv[v] = bm;
                else {
                    double d = std::fabs(bv[v] - bm);
                    d = std::min(d, 1 - d);
                    if (d > 1e-9) ++c.seamVertexMismatches;
                }
            }
            adj[v].push_back(mesh_.tri[t][(k + 1) % 3]);
            adj[v].push_back(mesh_.tri[t][(k + 2) % 3]);
        }
    for (int v = 0; v < V; ++v) {
        if (av[v] == 0.0 || av[v] == 1.0) {
            ++c.aPoles;
            continue;
        }
        bool lower = false, higher = false;
        for (int w : adj[v]) {
            lower = lower || av[w] < av[v];
            higher = higher || av[w] > av[v];
        }
        if (!lower || !higher) ++c.aExtrema;
    }
    c.monotoneA = (c.aCornerMismatches == 0 && c.aExtrema == 0 && c.aPoles == 2) ? Status::Pass : Status::Fail;
}

// Seam sweep: both sides of b = 0 == b = 1; a grid of interior points must all be found.
void Chart::checkSeam(ChartCheck& c) const {
    Vec3 lo, hi;
    mesh_.bounds(lo, hi);
    const double ext = std::max(1e-300, norm(hi - lo));
    const int N = 400;
    for (int i = 0; i <= N; ++i) {
        double a = aMin_ + (aMax_ - aMin_) * i / N;
        ChartPoint s0 = evalRaw(a, 0.0), s1 = evalRaw(a, 1.0);
        if (!s0.ok || !s1.ok) { ++c.domainGaps; continue; }
        ++c.seamSamples;
        c.seamMaxPositionGap = std::max(c.seamMaxPositionGap, norm(s0.p - s1.p));
        c.seamMaxNormalGap = std::max(c.seamMaxNormalGap, norm(s0.n - s1.n));
    }
    const double tolP = 1e-9 * ext, tolN = 1e-9;
    for (int i = 0; i <= 64; ++i)
        for (int j = 0; j < 64; ++j)
            if (!eval(aMin_ + (aMax_ - aMin_) * i / 64, (j + 0.5) / 64).ok) ++c.domainGaps;
    c.seam = (c.seamVertexMismatches == 0 && c.domainGaps == 0 && c.seamMaxPositionGap <= tolP && c.seamMaxNormalGap <= tolN)
                 ? Status::Pass
                 : Status::Fail;
}

ChartCheck Chart::check() const {
    ChartCheck c;
    checkOrientation(c);
    checkVertexValues(c);
    checkSeam(c);
    char buf[400];
    std::snprintf(buf, sizeof buf,
                  "orientation %s (degenerate %d, flipped %d, flipped vs outward %d); seam %s on a in [%.4f, %.4f] "
                  "(vertex mismatches %d, domain gaps %d, %d samples, max gap %.3g / normal %.3g); a %s (corner "
                  "mismatches %d, extrema %d, poles %d)",
                  toString(c.orientation), c.degenerateTriangles, c.flippedTriangles, c.flipped3D, toString(c.seam),
                  aMin_, aMax_, c.seamVertexMismatches, c.domainGaps, c.seamSamples, c.seamMaxPositionGap, c.seamMaxNormalGap,
                  toString(c.monotoneA), c.aCornerMismatches, c.aExtrema, c.aPoles);
    c.detail = buf;
    return c;
}

ChartSkew chartSkew(const TriMesh& charted) {
    const TriMesh& m = charted;
    std::vector<std::pair<double, double>> offByArea;
    double total = 0, weighted = 0;
    for (size_t t = 0; t < m.tri.size(); ++t) {
        const std::array<Vec2, 3>& uv = m.uv[t];
        if (uv[0].x <= 0 || uv[1].x <= 0 || uv[2].x <= 0 || uv[0].x >= 1 || uv[1].x >= 1 || uv[2].x >= 1) continue;
        const Vec3 e1 = m.pos[m.tri[t][1]] - m.pos[m.tri[t][0]], e2 = m.pos[m.tri[t][2]] - m.pos[m.tri[t][0]];
        const Vec3 n = cross(e1, e2);
        const double twiceArea = norm(n);
        if (twiceArea <= 0) continue;
        // grad f = (n x e2 df1 - n x e1 df2) / |n|^2 for the linear f with f0 = 0, f1 = df1, f2 = df2.
        const Vec3 ga = cross(n, e2) * (uv[1].x - uv[0].x) - cross(n, e1) * (uv[2].x - uv[0].x);
        const Vec3 gb = cross(n, e2) * (uv[1].y - uv[0].y) - cross(n, e1) * (uv[2].y - uv[0].y);
        if (norm(ga) <= 0 || norm(gb) <= 0) continue;
        const double off = std::asin(std::min(1.0, std::fabs(dot(ga, gb)) / (norm(ga) * norm(gb)))) * 180 / kPi;
        offByArea.push_back({off, twiceArea});
        total += twiceArea;
        weighted += off * twiceArea;
    }
    ChartSkew skew;
    if (offByArea.empty()) return skew;
    skew.mean = weighted / total;
    std::sort(offByArea.begin(), offByArea.end());
    double running = 0;
    skew.p95 = offByArea.back().first;
    for (const std::pair<double, double>& oa : offByArea) {
        running += oa.second;
        if (running >= 0.95 * total) {
            skew.p95 = oa.first;
            break;
        }
    }
    return skew;
}

} // namespace knit
