#include "knitting/layout/uv_atlas.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <numeric>
#include <stdexcept>

#include "model/cancel.h"

namespace knit {

Mat2 sqrtSpd(const Mat2& g) {
    // sqrt(G) = (G + s I) / t, s = sqrt(det G), t = sqrt(trace G + 2 s).
    const double s = std::sqrt(std::max(0.0, g.det())), t = std::sqrt(std::max(1e-300, g.a + g.d + 2 * s));
    return {(g.a + s) / t, g.b / t, g.c / t, (g.d + s) / t};
}

namespace {

double cross2(const Vec2& a, const Vec2& b) { return a.x * b.y - a.y * b.x; }
double dot2(const Vec2& a, const Vec2& b) { return a.x * b.x + a.y * b.y; }

// Area-weighted quantiles {min, p5, p50, p95, max}.
void quantiles(std::vector<std::pair<double, double>> v, double out[5]) {
    if (v.empty()) return;
    std::sort(v.begin(), v.end());
    double total = 0;
    for (const auto& x : v) total += x.second;
    out[0] = v.front().first;
    out[4] = v.back().first;
    const double qs[3] = {0.05, 0.5, 0.95};
    for (int k = 0; k < 3; ++k) {
        double acc = 0;
        out[k + 1] = v.back().first;
        for (const auto& x : v) {
            acc += x.second;
            if (acc >= qs[k] * total) {
                out[k + 1] = x.first;
                break;
            }
        }
    }
}

// Interiors of two 2D triangles overlap by more than tol (separating axis on the 6 edge normals).
bool overlap2(const std::array<Vec2, 3>& a, const std::array<Vec2, 3>& b, double tol) {
    for (const auto* t : {&a, &b})
        for (int k = 0; k < 3; ++k) {
            const Vec2 e = (*t)[(k + 1) % 3] - (*t)[k];
            const Vec2 n{-e.y, e.x};
            const double l = std::sqrt(dot2(n, n));
            if (l < 1e-300) continue;
            double a0 = 1e300, a1 = -1e300, b0 = 1e300, b1 = -1e300;
            for (int q = 0; q < 3; ++q) {
                a0 = std::min(a0, dot2(a[q], n) / l);
                a1 = std::max(a1, dot2(a[q], n) / l);
                b0 = std::min(b0, dot2(b[q], n) / l);
                b1 = std::max(b1, dot2(b[q], n) / l);
            }
            if (std::max(b0 - a1, a0 - b1) > -tol) return false;  // separated or only touching
        }
    return true;
}

// The uv of vertex v's corner in triangle t.
Vec2 uvAt(const TriMesh& m, int t, int v) {
    for (int k = 0; k < 3; ++k)
        if (m.tri[t][k] == v) return m.uv[t][k];
    return Vec2{};
}

bool sameUv(const Vec2& a, const Vec2& b, double tol) { return std::fabs(a.x - b.x) <= tol && std::fabs(a.y - b.y) <= tol; }

} // namespace

// What the constructor's stages share. The corners and their box are those of the chart in work.
struct UvAtlas::Construction {
    Vec3 up;
    double uvTol, lo, hi;
    std::map<std::pair<int, int>, std::vector<int>> edgeTris;  // undirected edge -> its triangles
    std::vector<double> detU, area3;  // per triangle: the uv determinant, the 3D area
    std::vector<char> degenerate;
    std::vector<std::array<Vec2, 3>> corners;  // uv' per triangle of the chart
    Vec2 cornerLo, cornerHi;
    long long steps = 0;

    // Long loops poll for cancellation every so many steps.
    void tick(int every = 0xffff) {
        if ((++steps & every) == 0) cancelCheckpoint();
    }
    bool continuousAcross(const TriMesh& m, const std::pair<int, int>& edge, int t0, int t1) const {
        return sameUv(uvAt(m, t0, edge.first), uvAt(m, t1, edge.first), uvTol) && sameUv(uvAt(m, t0, edge.second), uvAt(m, t1, edge.second), uvTol);
    }
};

// Triangles joined across edges whose two corner uvs agree on both sides.
void UvAtlas::joinCharts(Construction& b) {
    const TriMesh& m = m_;
    const int T = static_cast<int>(m.tri.size());
    for (int t = 0; t < T; ++t)
        for (int k = 0; k < 3; ++k) {
            b.tick();
            const int v0 = m.tri[t][k], v1 = m.tri[t][(k + 1) % 3];
            b.edgeTris[{std::min(v0, v1), std::max(v0, v1)}].push_back(t);
        }
    std::vector<int> par(T);
    std::iota(par.begin(), par.end(), 0);
    auto find = [&](int x) {
        while (par[x] != x) x = par[x] = par[par[x]];
        return x;
    };
    for (const auto& [e, ts] : b.edgeTris) {
        b.tick();
        if (ts.size() == 2 && b.continuousAcross(m, e, ts[0], ts[1])) par[find(ts[0])] = find(ts[1]);
    }
    std::map<int, int> rootChart;
    for (int t = 0; t < T; ++t) {
        b.tick();
        auto it = rootChart.find(find(t));
        if (it == rootChart.end()) it = rootChart.insert({find(t), static_cast<int>(charts_.size())}).first, charts_.emplace_back();
        chartOf_[t] = it->second;
        charts_[it->second].tris.push_back(t);
    }
}

// Per triangle: the uv determinant, the 3D area and the Jacobian J = E U^-1, U = [u1 u2] (columns).
void UvAtlas::measureTriangles(Construction& b) {
    const TriMesh& m = m_;
    const int T = static_cast<int>(m.tri.size());
    b.detU.resize(T);
    b.area3.resize(T);
    b.degenerate.assign(T, 0);
    for (int t = 0; t < T; ++t) {
        b.tick();
        const Vec3 e1 = m.pos[m.tri[t][1]] - m.pos[m.tri[t][0]], e2 = m.pos[m.tri[t][2]] - m.pos[m.tri[t][0]];
        const Vec2 u1 = m.uv[t][1] - m.uv[t][0], u2 = m.uv[t][2] - m.uv[t][0];
        b.detU[t] = cross2(u1, u2);
        b.area3[t] = 0.5 * norm(cross(e1, e2));
        if (b.detU[t] != 0) {
            const double k = 1.0 / b.detU[t];
            const double i00 = u2.y * k, i01 = -u2.x * k, i10 = -u1.y * k, i11 = u1.x * k;
            const Vec3 Ju = e1 * i00 + e2 * i10, Jv = e1 * i01 + e2 * i11;
            J_[t] = {Ju.x, Ju.y, Ju.z, Jv.x, Jv.y, Jv.z};
        }
    }
}

// The chart's majority orientation, its metric M and centre over the majority-oriented,
// non-degenerate triangles, which become usable. false: no usable triangle (UNRESOLVED).
bool UvAtlas::chartMetric(int c, Construction& b) {
    const TriMesh& m = m_;
    UvChart& ch = charts_[c];
    double meanDet = 0;
    for (int t : ch.tris) meanDet += std::fabs(b.detU[t]);
    meanDet /= std::max<size_t>(1, ch.tris.size());
    double pos = 0, neg = 0;
    for (int t : ch.tris) {
        b.tick();
        ch.area += b.area3[t];
        const Vec3 Ju{J_[t][0], J_[t][1], J_[t][2]}, Jv{J_[t][3], J_[t][4], J_[t][5]};
        if (std::fabs(b.detU[t]) < 1e-12 * meanDet || norm(cross(Ju, Jv)) <= 1e-12 * norm(Ju) * norm(Jv) || b.area3[t] <= 0) {
            b.degenerate[t] = 1;
            ch.degenerateArea += b.area3[t];
            continue;
        }
        (b.detU[t] > 0 ? pos : neg) += b.area3[t];
    }
    const double sign = pos >= neg ? 1.0 : -1.0;
    ch.mirrored = sign < 0;
    Mat2 G{0, 0, 0, 0};
    double wsum = 0;
    Vec2 cen{0, 0};
    for (int t : ch.tris) {
        b.tick();
        if (b.degenerate[t]) continue;
        if (b.detU[t] * sign < 0) {
            ch.foldedArea += b.area3[t];
            continue;
        }
        const Vec3 Ju{J_[t][0], J_[t][1], J_[t][2]}, Jv{J_[t][3], J_[t][4], J_[t][5]};
        const double guu = dot(Ju, Ju), guv = dot(Ju, Jv), gvv = dot(Jv, Jv);
        G = {G.a + b.area3[t] * guu, G.b + b.area3[t] * guv, G.c + b.area3[t] * guv, G.d + b.area3[t] * gvv};
        wsum += b.area3[t];
        const Vec2 uc = (m.uv[t][0] + m.uv[t][1] + m.uv[t][2]) * (1.0 / 3);
        cen = cen + uc * b.area3[t];
        usable_[t] = 1;
    }
    if (wsum <= 0) {
        ch.status = Status::Unresolved;
        ch.detail = "no usable triangle";
        return false;
    }
    G = {G.a / wsum, G.b / wsum, G.c / wsum, G.d / wsum};
    ch.M = sqrtSpd(G);
    if (ch.mirrored) ch.M = Mat2{-1, 0, 0, 1} * ch.M;
    ch.centre = cen * (1.0 / wsum);
    return true;
}

// The chart's corners in uv' and a grid over them: a cell about two triangle sizes wide, each
// triangle in every cell its box meets.
void UvAtlas::buildChartGrid(int c, Construction& b) {
    const UvChart& ch = charts_[c];
    b.corners.resize(ch.tris.size());
    Vec2 blo{1e300, 1e300}, bhi{-1e300, -1e300};
    double uvArea = 0;
    for (size_t q = 0; q < ch.tris.size(); ++q) {
        b.tick();
        b.corners[q] = uvp(ch.tris[q]);
        for (const Vec2& v : b.corners[q]) {
            blo = {std::min(blo.x, v.x), std::min(blo.y, v.y)};
            bhi = {std::max(bhi.x, v.x), std::max(bhi.y, v.y)};
        }
        uvArea += 0.5 * std::fabs(cross2(b.corners[q][1] - b.corners[q][0], b.corners[q][2] - b.corners[q][0]));
    }
    b.cornerLo = blo;
    b.cornerHi = bhi;
    Grid& g = grid_[c];
    g.lo = blo;
    g.cell = std::max(1e-9, 2 * std::sqrt(uvArea / std::max<size_t>(1, ch.tris.size())));
    g.nx = std::max(1, static_cast<int>(std::ceil((bhi.x - blo.x) / g.cell)) + 1);
    g.ny = std::max(1, static_cast<int>(std::ceil((bhi.y - blo.y) / g.cell)) + 1);
    g.cells.assign(static_cast<size_t>(g.nx) * g.ny, {});
    for (size_t q = 0; q < ch.tris.size(); ++q) {
        b.tick();
        double ax = 1e300, bx = -1e300, ay = 1e300, by = -1e300;
        for (const Vec2& v : b.corners[q]) ax = std::min(ax, v.x), bx = std::max(bx, v.x), ay = std::min(ay, v.y), by = std::max(by, v.y);
        const int x0 = std::clamp(static_cast<int>(std::floor((ax - g.lo.x) / g.cell)), 0, g.nx - 1);
        const int x1 = std::clamp(static_cast<int>(std::floor((bx - g.lo.x) / g.cell)), 0, g.nx - 1);
        const int y0 = std::clamp(static_cast<int>(std::floor((ay - g.lo.y) / g.cell)), 0, g.ny - 1);
        const int y1 = std::clamp(static_cast<int>(std::floor((by - g.lo.y) / g.cell)), 0, g.ny - 1);
        for (int y = y0; y <= y1; ++y)
            for (int x = x0; x <= x1; ++x) g.cells[static_cast<size_t>(y) * g.nx + x].push_back(static_cast<int>(q));
    }
}

// Overlaps (non-injective uv'): non-adjacent usable pairs with intersecting interiors become
// unusable; a chart with folded + overlap area > 1 % is UNRESOLVED.
void UvAtlas::chartValidity(int c, Construction& b) {
    const TriMesh& m = m_;
    UvChart& ch = charts_[c];
    const Grid& g = grid_[c];
    const std::vector<std::array<Vec2, 3>>& P = b.corners;
    const double tol = 1e-9 * std::max(1e-12, std::max(b.cornerHi.x - b.cornerLo.x, b.cornerHi.y - b.cornerLo.y));
    std::vector<char> overlapped(ch.tris.size(), 0);
    for (const std::vector<int>& cell : g.cells)
        for (size_t i = 0; i < cell.size(); ++i) {
            b.tick(0xfff);
            for (size_t j = i + 1; j < cell.size(); ++j) {
                const int qa = cell[i], qb = cell[j];
                const int ta = ch.tris[qa], tb = ch.tris[qb];
                if (!usable_[ta] || !usable_[tb] || (overlapped[qa] && overlapped[qb])) continue;
                bool adjacent = false;
                for (int u = 0; u < 3 && !adjacent; ++u)
                    for (int v = 0; v < 3; ++v) adjacent = adjacent || m.tri[ta][u] == m.tri[tb][v];
                if (adjacent || !overlap2(P[qa], P[qb], tol)) continue;
                overlapped[qa] = overlapped[qb] = 1;
            }
        }
    for (size_t q = 0; q < ch.tris.size(); ++q)
        if (overlapped[q] && usable_[ch.tris[q]]) {
            usable_[ch.tris[q]] = 0;
            ch.overlapArea += b.area3[ch.tris[q]];
        }
    if (ch.foldedArea + ch.overlapArea > 0.01 * ch.area) ch.status = Status::Unresolved;
}

// The singular values of K_t = J_t M^-1 per usable triangle and their area-weighted quantiles.
void UvAtlas::chartDistortion(int c, Construction& b) {
    UvChart& ch = charts_[c];
    std::vector<std::pair<double, double>> q1, q2, qr, qa;
    double outside = 0, usableArea = 0;
    for (int t : ch.tris) {
        b.tick();
        if (!usable_[t]) continue;
        const double area = b.area3[t];
        const SingularValues sv = singularValues(t);
        s1_[t] = sv.s1;
        s2_[t] = sv.s2;
        q1.push_back({s1_[t], area});
        q2.push_back({s2_[t], area});
        qr.push_back({s1_[t] / std::max(1e-300, s2_[t]), area});
        qa.push_back({s1_[t] * s2_[t], area});
        usableArea += area;
        if (s1_[t] < b.lo || s1_[t] > b.hi || s2_[t] < b.lo || s2_[t] > b.hi) outside += area;
    }
    quantiles(q1, ch.distortion.s1);
    quantiles(q2, ch.distortion.s2);
    quantiles(qr, ch.distortion.ratio);
    quantiles(qa, ch.distortion.areaScale);
    ch.distortion.outsideFraction = usableArea > 0 ? outside / usableArea : 0;
}

// World up pulled back into uv' per usable triangle (d = M J^+ u_t, J^+ = (J^T J)^-1 J^T), its
// area-weighted mean direction, or the uv' v axis where that mean is too short to mean anything.
void UvAtlas::chartOrientation(int c, Construction& b) {
    const TriMesh& m = m_;
    UvChart& ch = charts_[c];
    double usableArea = 0;
    Vec2 upSum{0, 0};
    for (int t : ch.tris) {
        b.tick();
        if (!usable_[t]) continue;
        usableArea += b.area3[t];
        const Vec3 Ju{J_[t][0], J_[t][1], J_[t][2]}, Jv{J_[t][3], J_[t][4], J_[t][5]};
        const Vec3 n = normalized(cross(m.pos[m.tri[t][1]] - m.pos[m.tri[t][0]], m.pos[m.tri[t][2]] - m.pos[m.tri[t][0]]));
        const Vec3 u = b.up - n * dot(b.up, n);
        if (norm(u) < 1e-6) continue;
        const Mat2 JtJ{dot(Ju, Ju), dot(Ju, Jv), dot(Jv, Ju), dot(Jv, Jv)};
        const Vec2 duv = JtJ.inverse() * Vec2{dot(Ju, u), dot(Jv, u)};
        const Vec2 d = ch.M * duv;
        const double l = std::sqrt(dot2(d, d));
        if (l > 0) upSum = upSum + d * (b.area3[t] / l);
    }
    const Vec2 mean = upSum * (1.0 / std::max(1e-300, usableArea));
    ch.upMeanLength = std::sqrt(dot2(mean, mean));
    if (ch.upMeanLength >= 0.2) ch.upDir = mean * (1.0 / ch.upMeanLength);
    else {
        ch.upDir = {0, 1};
        ch.upFromWorld = false;
    }
}

// The boundary of the usable region, in uv', and its part on open edges of the mesh.
void UvAtlas::collectBoundaries(int c, Construction& b) {
    const TriMesh& m = m_;
    const UvChart& ch = charts_[c];
    for (int t : ch.tris) {
        b.tick();
        if (!usable_[t]) continue;
        const std::array<Vec2, 3> pt = uvp(t);
        for (int k = 0; k < 3; ++k) {
            const int v0 = m.tri[t][k], v1 = m.tri[t][(k + 1) % 3];
            const std::pair<int, int> edge{std::min(v0, v1), std::max(v0, v1)};
            const std::vector<int>& ts = b.edgeTris[edge];
            // Inner only if the neighbour is usable, in this chart, AND uv-continuous across this
            // edge: a seam can run inside one chart (a cut cylinder is one chart joined the long way
            // round), and it bounds the uv' region there too.
            bool inner = false;
            if (ts.size() == 2) {
                const int o = ts[0] == t ? ts[1] : ts[0];
                inner = chartOf_[o] == c && usable_[o] && b.continuousAcross(m, edge, t, o);
            }
            if (!inner) boundary_[c].push_back({pt[k], pt[(k + 1) % 3]});
            if (ts.size() == 1) openBoundary_[c].push_back({pt[k], pt[(k + 1) % 3]});
        }
    }
}

UvAtlas::UvAtlas(const TriMesh& m, const Vec3& up, double uvTol, double lo, double hi) : m_(m) {
    const int T = static_cast<int>(m.tri.size());
    if (m.uv.size() != m.tri.size()) throw std::invalid_argument("UvAtlas: mesh without per-corner uv");
    chartOf_.assign(T, -1);
    usable_.assign(T, 0);
    J_.resize(T);
    s1_.assign(T, 0);
    s2_.assign(T, 0);
    Construction b{up, uvTol, lo, hi};
    joinCharts(b);
    measureTriangles(b);
    boundary_.resize(charts_.size());
    openBoundary_.resize(charts_.size());
    grid_.resize(charts_.size());
    for (int c = 0; c < static_cast<int>(charts_.size()); ++c) {
        if (!chartMetric(c, b)) continue;
        buildChartGrid(c, b);
        chartValidity(c, b);
        chartDistortion(c, b);
        chartOrientation(c, b);
        collectBoundaries(c, b);
        UvChart& ch = charts_[c];
        char buf[256];
        std::snprintf(buf, sizeof buf, "%zu tris, area %.0f mm2 (degenerate %.1f, folded %.1f, overlap %.1f)%s; s1 p50 %.3f max %.3f, s2 min %.3f p50 %.3f; outside %.1f %%; up %s (|mean| %.2f)",
                      ch.tris.size(), ch.area, ch.degenerateArea, ch.foldedArea, ch.overlapArea, ch.mirrored ? ", mirrored" : "",
                      ch.distortion.s1[2], ch.distortion.s1[4], ch.distortion.s2[0], ch.distortion.s2[2],
                      100 * ch.distortion.outsideFraction, ch.upFromWorld ? "world" : "uv v (fallback)", ch.upMeanLength);
        ch.detail = buf;
    }
}

std::array<Vec2, 3> UvAtlas::uvp(int t) const {
    const UvChart& ch = charts_[chartOf_[t]];
    return {ch.M * (m_.uv[t][0] - ch.centre), ch.M * (m_.uv[t][1] - ch.centre), ch.M * (m_.uv[t][2] - ch.centre)};
}

Vec3 UvAtlas::pushForward(int t, const Vec2& d) const {
    const Vec2 duv = charts_[chartOf_[t]].M.inverse() * d;
    const std::array<double, 6>& J = J_[t];
    return Vec3{J[0], J[1], J[2]} * duv.x + Vec3{J[3], J[4], J[5]} * duv.y;
}

UvAtlas::SingularValues UvAtlas::singularValues(int t) const {
    const Vec3 a = pushForward(t, {1, 0}), b = pushForward(t, {0, 1});
    const double p = dot(a, a), q = dot(a, b), r = dot(b, b);
    const double tr = 0.5 * (p + r), disc = std::sqrt(std::max(0.0, 0.25 * (p - r) * (p - r) + q * q));
    return {std::sqrt(tr + disc), std::sqrt(std::max(0.0, tr - disc))};
}

std::vector<int> UvAtlas::locate(int chart, const Vec2& p, double bary[3]) const {
    std::vector<int> hits;
    const Grid& g = grid_[chart];
    const int x = static_cast<int>(std::floor((p.x - g.lo.x) / g.cell)), y = static_cast<int>(std::floor((p.y - g.lo.y) / g.cell));
    if (x < 0 || y < 0 || x >= g.nx || y >= g.ny) return hits;
    const UvChart& ch = charts_[chart];
    for (int q : g.cells[static_cast<size_t>(y) * g.nx + x]) {
        const int t = ch.tris[q];
        if (!usable_[t]) continue;
        const std::array<Vec2, 3> P = uvp(t);
        const double A = cross2(P[1] - P[0], P[2] - P[0]);
        if (A == 0) continue;
        const double b1 = cross2(p - P[0], P[2] - P[0]) / A, b2 = cross2(P[1] - P[0], p - P[0]) / A, b0 = 1 - b1 - b2;
        const double eps = -1e-12;
        if (b0 >= eps && b1 >= eps && b2 >= eps) {
            if (hits.empty()) {
                bary[0] = std::max(0.0, b0);
                bary[1] = std::max(0.0, b1);
                bary[2] = std::max(0.0, b2);
                const double s = bary[0] + bary[1] + bary[2];
                for (int k = 0; k < 3; ++k) bary[k] /= s;
            }
            hits.push_back(t);
        }
    }
    // A point on a shared edge or vertex is in several adjacent triangles: one location.
    if (hits.size() > 1) {
        bool allAdjacent = true;
        for (size_t i = 1; i < hits.size() && allAdjacent; ++i) {
            bool shares = false;
            for (int u = 0; u < 3; ++u)
                for (int v = 0; v < 3; ++v) shares = shares || m_.tri[hits[0]][u] == m_.tri[hits[i]][v];
            allAdjacent = shares;
        }
        if (allAdjacent) hits.resize(1);
    }
    return hits;
}

double UvAtlas::minS2In(int chart, const Vec2& lo, const Vec2& hi) const {
    const Grid& g = grid_[chart];
    const UvChart& ch = charts_[chart];
    const int x0 = std::clamp(static_cast<int>(std::floor((lo.x - g.lo.x) / g.cell)), 0, g.nx - 1);
    const int x1 = std::clamp(static_cast<int>(std::floor((hi.x - g.lo.x) / g.cell)), 0, g.nx - 1);
    const int y0 = std::clamp(static_cast<int>(std::floor((lo.y - g.lo.y) / g.cell)), 0, g.ny - 1);
    const int y1 = std::clamp(static_cast<int>(std::floor((hi.y - g.lo.y) / g.cell)), 0, g.ny - 1);
    double s = 1e300;
    for (int y = y0; y <= y1; ++y)
        for (int x = x0; x <= x1; ++x)
            for (int q : g.cells[static_cast<size_t>(y) * g.nx + x])
                if (usable_[ch.tris[q]]) s = std::min(s, s2_[ch.tris[q]]);
    return s;
}

} // namespace knit
