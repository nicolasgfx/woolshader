// How much of the base mesh the rendered yarn covers, and where the holes are (coverage.h).
#include "knitting/knit/coverage.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <queue>
#include <unordered_map>

#include "model/cancel.h"
#include "model/spatial_hash.h"

namespace knit {

namespace {

uint64_t mix64(uint64_t x) {
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdull;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ull;
    x ^= x >> 33;
    return x;
}
double hash01(uint64_t seed, uint64_t a, uint64_t b) {
    return static_cast<double>(mix64(seed ^ mix64(a * 0x9e3779b97f4a7c15ull + b)) >> 11) * (1.0 / 9007199254740992.0);
}

struct Sample {
    int tri = -1;
    double b[3] = {0, 0, 0};
    Vec3 p, n;
};

// Segment q0 -> q1 against triangle (a, b, c) (Moller-Trumbore, both sides; t in [0, 1]).
bool segmentHitsTriangle(const Vec3& q0, const Vec3& q1, const Vec3& a, const Vec3& b, const Vec3& c) {
    const Vec3 d = q1 - q0, e1 = b - a, e2 = c - a;
    const Vec3 pv = cross(d, e2);
    const double det = dot(e1, pv);
    if (std::fabs(det) < 1e-18) return false;
    const double inv = 1.0 / det;
    const Vec3 tv = q0 - a;
    const double u = dot(tv, pv) * inv;
    if (u < 0 || u > 1) return false;
    const Vec3 qv = cross(tv, e1);
    const double v = dot(d, qv) * inv;
    if (v < 0 || u + v > 1) return false;
    const double t = dot(e2, qv) * inv;
    return t >= 0 && t <= 1;
}

// The tube's triangles in a spatial hash, for the samples' normal segments.
class TubeHash {
public:
    explicit TubeHash(const TubeMesh& m) : m_(m) {
        const std::vector<Vec3>& pos = m.pos;
        for (size_t k = 0; k + 2 < m.tri.size(); k += 3) {
            const Vec3 &a = pos[m.tri[k]], &b = pos[m.tri[k + 1]], &c = pos[m.tri[k + 2]];
            hash_.insert({std::min({a.x, b.x, c.x}), std::min({a.y, b.y, c.y}), std::min({a.z, b.z, c.z})},
                         {std::max({a.x, b.x, c.x}), std::max({a.y, b.y, c.y}), std::max({a.z, b.z, c.z})});
        }
        hash_.finalize();
    }
    // Every tube triangle the segment meets, into met (cleared first).
    void trianglesMet(const Vec3& q0, const Vec3& q1, std::vector<int>& met) const {
        met.clear();
        hash_.query(minCorner(q0, q1), maxCorner(q0, q1), [&](int id) {
            if (segmentMeets(id, q0, q1)) met.push_back(id);
        });
    }

private:
    static Vec3 minCorner(const Vec3& a, const Vec3& b) { return {std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z)}; }
    static Vec3 maxCorner(const Vec3& a, const Vec3& b) { return {std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)}; }
    bool segmentMeets(int triangle, const Vec3& q0, const Vec3& q1) const {
        const size_t k = static_cast<size_t>(triangle) * 3;
        return segmentHitsTriangle(q0, q1, m_.pos[m_.tri[k]], m_.pos[m_.tri[k + 1]], m_.pos[m_.tri[k + 2]]);
    }
    const TubeMesh& m_;
    SpatialHash3 hash_{0.5};
};

// Triangle adjacency of a welded mesh (edges by vertex pair; non-manifold edges join all).
std::vector<std::vector<int>> triangleAdjacency(const TriMesh& m) {
    std::unordered_map<uint64_t, std::vector<int>> edges;
    edges.reserve(m.tri.size() * 2);
    for (size_t t = 0; t < m.tri.size(); ++t)
        for (int e = 0; e < 3; ++e) {
            const uint32_t a = static_cast<uint32_t>(m.tri[t][e]), b = static_cast<uint32_t>(m.tri[t][(e + 1) % 3]);
            edges[(uint64_t(std::min(a, b)) << 32) | std::max(a, b)].push_back(static_cast<int>(t));
        }
    std::vector<std::vector<int>> adj(m.tri.size());
    for (const auto& [k, ts] : edges)
        for (int a : ts)
            for (int b : ts)
                if (a != b) adj[a].push_back(b);
    for (auto& v : adj) {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    }
    return adj;
}

// Bounded A* over the triangle graph (centroid distances): is f within G of t?
class Ownership {
public:
    Ownership(const std::vector<std::vector<int>>& adj, const std::vector<Vec3>& centroid, const std::vector<Vec3>& faceN, double G, double minCos)
        : adj_(adj), c_(centroid), n_(faceN), G_(G), minCos_(minCos), g_(adj.size(), 0), stamp_(adj.size(), 0) {}
    bool owned(int t, const Vec3& sampleN, int f) {
        if (f < 0) return false;
        if (dot(n_[f], sampleN) <= minCos_) return false;
        if (t == f) return true;
        if (norm(c_[f] - c_[t]) > G_) return false;
        const uint64_t key = (uint64_t(uint32_t(t)) << 32) | uint32_t(f);
        if (auto it = cache_.find(key); it != cache_.end()) return it->second;
        const bool r = search(t, f);
        if (cache_.size() > 4000000) cache_.clear();
        cache_[key] = r;
        return r;
    }

private:
    bool search(int t, int f) {
        if (++cur_ == 0) {
            std::fill(stamp_.begin(), stamp_.end(), 0);
            cur_ = 1;
        }
        using Item = std::pair<double, int>;  // (g + h, triangle)
        std::priority_queue<Item, std::vector<Item>, std::greater<Item>> open;
        stamp_[t] = cur_;
        g_[t] = 0;
        open.push({norm(c_[f] - c_[t]), t});
        long expansions = 0;
        while (!open.empty() && ++expansions < 200000) {
            const auto [fh, u] = open.top();
            open.pop();
            if (u == f) return true;
            if (fh > G_) return false;
            const double gu = g_[u];
            if (fh - norm(c_[f] - c_[u]) > gu + 1e-12) continue;  // stale entry
            for (int v : adj_[u]) {
                const double gv = gu + norm(c_[v] - c_[u]);
                if (gv > G_) continue;
                if (stamp_[v] == cur_ && g_[v] <= gv) continue;
                stamp_[v] = cur_;
                g_[v] = gv;
                open.push({gv + norm(c_[f] - c_[v]), v});
            }
        }
        return false;
    }
    const std::vector<std::vector<int>>& adj_;
    const std::vector<Vec3>& c_;
    const std::vector<Vec3>& n_;
    double G_, minCos_;
    std::vector<double> g_;
    std::vector<uint32_t> stamp_;
    uint32_t cur_ = 0;
    std::unordered_map<uint64_t, char> cache_;
};

struct SampleNeighbour {
    int sample;
    double distance;
};

// The sample set's scales: the spacing, the neighbour radius of a region, the area per sample.
struct SampleGeometry {
    double spacing, radius, sampleArea;
};

// Fixed-radius neighbour lookup over the samples (sorted cell keys).
class PointGrid {
public:
    PointGrid(const std::vector<Sample>& s, double cell) : s_(s), cell_(cell) {
        keys_.resize(s.size());
        for (size_t i = 0; i < s.size(); ++i) keys_[i] = {key(cellOf(s[i].p)), static_cast<int>(i)};
        std::sort(keys_.begin(), keys_.end());
    }
    // The samples within r of sample i, into out (cleared first).
    void near(int i, double r, std::vector<SampleNeighbour>& out) const {
        out.clear();
        const std::array<int64_t, 3> c = cellOf(s_[i].p);
        for (int dz = -1; dz <= 1; ++dz)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dx = -1; dx <= 1; ++dx) {
                    const uint64_t k = key({c[0] + dx, c[1] + dy, c[2] + dz});
                    auto it = std::lower_bound(keys_.begin(), keys_.end(), std::make_pair(k, -1));
                    for (; it != keys_.end() && it->first == k; ++it) {
                        const int j = it->second;
                        if (j == i) continue;
                        const double d = norm(s_[j].p - s_[i].p);
                        if (d <= r) out.push_back({j, d});
                    }
                }
    }

private:
    std::array<int64_t, 3> cellOf(const Vec3& p) const {
        return {static_cast<int64_t>(std::floor(p.x / cell_)), static_cast<int64_t>(std::floor(p.y / cell_)),
                static_cast<int64_t>(std::floor(p.z / cell_))};
    }
    static uint64_t key(const std::array<int64_t, 3>& c) {
        const uint64_t m = (1ull << 21) - 1;
        return ((uint64_t(c[0]) & m) << 42) | ((uint64_t(c[1]) & m) << 21) | (uint64_t(c[2]) & m);
    }
    const std::vector<Sample>& s_;
    double cell_;
    std::vector<std::pair<uint64_t, int>> keys_;
};

// Is one of these tube triangles' vertices bound to a base triangle that owns the sample?
bool anyEligible(const std::vector<int>& triangles, const TubeMesh& tube, const std::vector<int>& feet, const Sample& q, Ownership& ow) {
    for (int id : triangles)
        for (int c = 0; c < 3; ++c) {
            const uint32_t v = tube.tri[static_cast<size_t>(id) * 3 + c];
            if (v < feet.size() && ow.owned(q.tri, q.n, feet[v])) return true;
        }
    return false;
}

// Connected regions of the flagged samples: per sample its region (-1: none) and, once
// spreadRegionDistances has run, its distance to the nearest sample outside the region.
struct RegionLabels {
    std::vector<int> label;
    std::vector<double> dist;
    std::vector<int> size;
    std::vector<Vec3> lo, hi;
};

RegionLabels labelRegions(const std::vector<Sample>& s, const std::vector<char>& flag, const PointGrid& grid, double radius) {
    RegionLabels out;
    const size_t n = s.size();
    out.label.assign(n, -1);
    out.dist.assign(n, 1e300);
    std::vector<SampleNeighbour> around;
    for (size_t i = 0; i < n; ++i) {
        if (!flag[i] || out.label[i] >= 0) continue;
        const int id = static_cast<int>(out.size.size());
        out.size.push_back(0);
        out.lo.push_back(s[i].p);
        out.hi.push_back(s[i].p);
        std::vector<int> stack{static_cast<int>(i)};
        out.label[i] = id;
        while (!stack.empty()) {
            const int u = stack.back();
            stack.pop_back();
            ++out.size[id];
            out.lo[id] = {std::min(out.lo[id].x, s[u].p.x), std::min(out.lo[id].y, s[u].p.y), std::min(out.lo[id].z, s[u].p.z)};
            out.hi[id] = {std::max(out.hi[id].x, s[u].p.x), std::max(out.hi[id].y, s[u].p.y), std::max(out.hi[id].z, s[u].p.z)};
            grid.near(u, radius, around);
            for (const SampleNeighbour& nb : around) {
                const int v = nb.sample;
                if (dot(s[v].n, s[u].n) <= 0.3) continue;
                if (!flag[v]) {
                    if (nb.distance < out.dist[u]) out.dist[u] = nb.distance;
                } else if (out.label[v] < 0) {
                    out.label[v] = id;
                    stack.push_back(v);
                }
            }
        }
        cancelCheckpoint();
    }
    return out;
}

// Dijkstra from the samples next to the outside, over the flagged samples.
void spreadRegionDistances(const std::vector<Sample>& s, const std::vector<char>& flag, const PointGrid& grid, double radius, std::vector<double>& dist) {
    using Item = std::pair<double, int>;
    std::priority_queue<Item, std::vector<Item>, std::greater<Item>> pq;
    for (size_t i = 0; i < s.size(); ++i)
        if (flag[i] && dist[i] < 1e300) pq.push({dist[i], static_cast<int>(i)});
    std::vector<SampleNeighbour> around;
    while (!pq.empty()) {
        const auto [d, u] = pq.top();
        pq.pop();
        if (d > dist[u]) continue;
        grid.near(u, radius, around);
        for (const SampleNeighbour& nb : around) {
            const int v = nb.sample;
            if (!flag[v] || dot(s[v].n, s[u].n) <= 0.3) continue;
            if (d + nb.distance < dist[v]) {
                dist[v] = d + nb.distance;
                pq.push({dist[v], v});
            }
        }
    }
}

// Connected regions of the flagged samples and their widths (header comment).
CoverageRegions regions(const std::vector<Sample>& s, const std::vector<char>& flag, const PointGrid& grid, const SampleGeometry& g) {
    CoverageRegions out;
    const size_t n = s.size();
    RegionLabels labels = labelRegions(s, flag, grid, g.radius);
    spreadRegionDistances(s, flag, grid, g.radius, labels.dist);
    const std::vector<double>& dist = labels.dist;
    std::vector<double> maxD(labels.size.size(), 0);
    std::vector<int> deepest(labels.size.size(), -1);
    for (size_t i = 0; i < n; ++i)
        if (labels.label[i] >= 0) {
            const int r = labels.label[i];
            const double d = dist[i] >= 1e300 ? 0 : dist[i];
            if (deepest[r] < 0 || d > (dist[deepest[r]] >= 1e300 ? 0 : dist[deepest[r]])) deepest[r] = static_cast<int>(i);
            maxD[r] = std::max(maxD[r], dist[i]);
        }
    out.count = static_cast<int>(labels.size.size());
    for (size_t r = 0; r < labels.size.size(); ++r) {
        const double area = labels.size[r] * g.sampleArea;
        // No sample outside the region next to it (the whole component): the box diagonal.
        const double width = maxD[r] >= 1e300 ? norm(labels.hi[r] - labels.lo[r]) : std::max(g.spacing, 2 * maxD[r] - g.spacing);
        out.totalArea += area;
        out.maxArea = std::max(out.maxArea, area);
        out.maxWidth = std::max(out.maxWidth, width);
        CoverageRegion item;
        item.area = area;
        item.width = width;
        item.lo = labels.lo[r];
        item.hi = labels.hi[r];
        if (deepest[r] >= 0) {
            const Sample& c = s[deepest[r]];
            item.centreTri = c.tri;
            for (int k = 0; k < 3; ++k) item.centreBary[k] = c.b[k];
            item.centre = c.p;
        }
        out.items.push_back(item);
    }
    std::stable_sort(out.items.begin(), out.items.end(), [](const CoverageRegion& a, const CoverageRegion& b) { return a.area > b.area; });
    return out;
}

// Inside CoverageInput's measured region (strictly).
bool insideRegion(const CoverageInput& in, const Vec3& p) {
    return p.x > in.regionLo.x && p.x < in.regionHi.x && p.y > in.regionLo.y && p.y < in.regionHi.y && p.z > in.regionLo.z &&
           p.z < in.regionHi.z;
}

// Per triangle of the base mesh: centroid, face normal and area.
struct TriangleMeasures {
    std::vector<Vec3> centroid, faceN;
    std::vector<double> area;
    double total = 0;
};

TriangleMeasures measureTriangles(const TriMesh& m) {
    TriangleMeasures out;
    const size_t nt = m.tri.size();
    out.centroid.resize(nt);
    out.faceN.resize(nt);
    out.area.resize(nt);
    for (size_t t = 0; t < nt; ++t) {
        const Vec3 &a = m.pos[m.tri[t][0]], &b = m.pos[m.tri[t][1]], &c = m.pos[m.tri[t][2]];
        out.centroid[t] = (a + b + c) * (1.0 / 3);
        const Vec3 cr = cross(b - a, c - a);
        out.area[t] = 0.5 * norm(cr);
        out.faceN[t] = out.area[t] > 0 ? cr * (0.5 / out.area[t]) : Vec3{0, 0, 1};
        out.total += out.area[t];
    }
    return out;
}

// Stratified samples per triangle: an R2 sequence from a hashed offset, kept inside the measured region;
// generated counts them before that test.
std::vector<Sample> stratifiedSamples(const CoverageInput& in, const CoverageParams& p, const TriangleMeasures& tris, double spacing, long long& generated) {
    const TriMesh& m = *in.mesh;
    std::vector<Sample> smp;
    generated = 0;
    for (size_t t = 0; t < m.tri.size(); ++t) {
        const double expected = tris.area[t] / (spacing * spacing);
        const long k = static_cast<long>(std::floor(expected + hash01(p.seed, t, 0)));
        const double o1 = hash01(p.seed, t, 1), o2 = hash01(p.seed, t, 2);
        const std::array<int, 3>& tv = m.tri[t];
        for (long j = 0; j < k; ++j) {
            double r1 = o1 + 0.7548776662466927 * (j + 1), r2 = o2 + 0.5698402909980532 * (j + 1);
            r1 -= std::floor(r1);
            r2 -= std::floor(r2);
            const double su = std::sqrt(r1);
            Sample q;
            q.tri = static_cast<int>(t);
            q.b[0] = 1 - su;
            q.b[1] = su * (1 - r2);
            q.b[2] = su * r2;
            q.p = m.pos[tv[0]] * q.b[0] + m.pos[tv[1]] * q.b[1] + m.pos[tv[2]] * q.b[2];
            Vec3 n = tris.faceN[t];
            if (m.vnormal.size() == m.pos.size()) {
                const Vec3 vn = m.vnormal[tv[0]] * q.b[0] + m.vnormal[tv[1]] * q.b[1] + m.vnormal[tv[2]] * q.b[2];
                if (norm(vn) > 1e-9) n = normalized(vn);
            }
            q.n = n;
            ++generated;
            if (!insideRegion(in, q.p)) continue;
            smp.push_back(q);
        }
    }
    return smp;
}

// Per sample, what its normal segment met: eligible yarn (covered), or only ineligible yarn (unowned).
struct SampleTraces {
    std::vector<char> covered, unowned;
};

// The segment of every sample against the tube, in parallel.
SampleTraces traceSamples(const CoverageInput& in, const CoverageParams& p, const std::vector<Sample>& smp, const TriangleMeasures& tris, const ExecPolicy& exec) {
    const std::vector<std::vector<int>> adj = triangleAdjacency(*in.mesh);
    const TubeHash own(*in.own.tube);
    SampleTraces out;
    out.covered.assign(smp.size(), 0);
    out.unowned.assign(smp.size(), 0);
    const std::vector<int>& feet = in.own.vertexFoot;
    const TubeMesh& tm = *in.own.tube;
    const int64_t chunk = 8192;
    parallelFor(static_cast<int64_t>(smp.size()), chunk,
                [&](int64_t b, int64_t e, int64_t) {
                    cancelCheckpoint();
                    Ownership ow(adj, tris.centroid, tris.faceN, p.ownershipRadius, p.minNormalCos);
                    std::vector<int> met;
                    for (int64_t i = b; i < e; ++i) {
                        const Sample& q = smp[i];
                        const Vec3 q0 = q.p - q.n * p.reach, q1 = q.p + q.n * p.reach;
                        own.trianglesMet(q0, q1, met);
                        out.covered[i] = !met.empty() && (feet.empty() || anyEligible(met, tm, feet, q, ow));
                        out.unowned[i] = !out.covered[i] && !met.empty();
                    }
                },
                exec);
    return out;
}

// The holes among the samples (the uncovered ones), counted into the report.
std::vector<char> countHoles(const SampleTraces& traces, CoverageReport& rep) {
    const size_t n = traces.covered.size();
    std::vector<char> holes(n);
    for (size_t i = 0; i < n; ++i) {
        holes[i] = !traces.covered[i];
        rep.uncovered += holes[i];
        rep.unownedHits += traces.unowned[i];
    }
    rep.uncoveredFrac = rep.uncovered / double(rep.samples);
    return holes;
}

} // namespace

std::string CoverageReport::summary() const {
    char b[640];
    std::snprintf(b, sizeof b,
                  "coverage: %lld samples at %.3f mm on %.0f mm2; uncovered %.2f %% (%d regions, max %.1f mm2, max width %.2f mm) "
                  "(%.0f ms)",
                  samples, spacing, area, 100 * uncoveredFrac, holes.count, holes.maxArea, holes.maxWidth, ms);
    return b;
}

CoverageReport measureCoverage(const CoverageInput& in, const CoverageParams& p, const ExecPolicy& exec) {
    const auto t0 = std::chrono::steady_clock::now();
    CoverageReport rep;
    const TriangleMeasures tris = measureTriangles(*in.mesh);
    rep.area = tris.total;
    const double s = std::max(p.spacing, p.maxSamples > 0 ? std::sqrt(tris.total / double(p.maxSamples)) : 0.0);
    rep.spacing = s;
    long long generated = 0;
    const std::vector<Sample> smp = stratifiedSamples(in, p, tris, s, generated);
    rep.samples = static_cast<long long>(smp.size());
    rep.sampleArea = generated ? tris.total / double(generated) : 0;  // the domain's samples have the mesh's density
    if (smp.empty()) return rep;
    cancelCheckpoint();
    const SampleTraces traces = traceSamples(in, p, smp, tris, exec);
    const std::vector<char> holes = countHoles(traces, rep);
    const SampleGeometry g{s, p.neighbourFactor * s, rep.sampleArea};
    const PointGrid grid(smp, g.radius);
    rep.holes = regions(smp, holes, grid, g);
    rep.ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    return rep;
}

std::vector<int> tubeVertexFeet(const TubeMesh& tube, const CurveSet& cs, const std::vector<std::vector<ProjState>>& states) {
    std::vector<int> feet(tube.pos.size(), -1);
    for (size_t v = 0; v < tube.pos.size(); ++v) {
        const int y = v < tube.yarn.size() ? tube.yarn[v] : -1;
        if (y < 0 || y >= static_cast<int>(cs.yarns.size()) || y >= static_cast<int>(states.size())) continue;
        const YarnPath& yp = cs.yarns[y];
        const std::vector<ProjState>& st = states[y];
        if (st.empty()) continue;
        const double u = yp.paramAtArcLength(std::clamp(tube.s[v], 0.0, yp.length()));
        const int c = std::clamp(static_cast<int>(std::lround(u + 1)), 0, static_cast<int>(st.size()) - 1);
        feet[v] = st[c].tri;
    }
    return feet;
}

// --- a flat swatch (knit_check --flat-dmax)

InteriorWindow interiorWindow(const CurveSet& cs) {
    std::vector<std::vector<Vec3>> centre(cs.courses, std::vector<Vec3>(cs.wales));
    for (const YarnPath& y : cs.yarns)
        for (const StitchSpan& sp : y.stitches) {
            Vec3 c{};
            for (int i = 0; i < 32; ++i) c += y.curve.eval(sp.uBegin + (i + 0.5) / 32 * (sp.uEnd - sp.uBegin));
            centre[sp.course][sp.wale] = c / 32.0;
        }
    InteriorWindow w{1e300, 1e300, -1e300, -1e300};
    double sw = 0, sc = 0;
    int nw = 0, nc = 0;
    for (int r = 1; r + 1 < cs.courses; ++r)
        for (int i = 1; i + 1 < cs.wales; ++i) {
            w.x0 = std::min(w.x0, centre[r][i].x); w.x1 = std::max(w.x1, centre[r][i].x);
            w.y0 = std::min(w.y0, centre[r][i].y); w.y1 = std::max(w.y1, centre[r][i].y);
            if (i + 2 < cs.wales) { sw += std::fabs(centre[r][i + 1].x - centre[r][i].x); ++nw; }
            if (r + 2 < cs.courses) { sc += std::fabs(centre[r + 1][i].y - centre[r][i].y); ++nc; }
        }
    double hw = nw ? 0.5 * sw / nw : 0, hc = nc ? 0.5 * sc / nc : 0;
    w.x0 -= hw; w.x1 += hw; w.y0 -= hc; w.y1 += hc;
    return w;
}

double meshOpenArea(const TubeMesh& m, const InteriorWindow& w, int res) {
    if (w.x1 <= w.x0 || w.y1 <= w.y0) return 0;
    const int nx = res, ny = std::max(2, static_cast<int>(std::round(res * (w.y1 - w.y0) / (w.x1 - w.x0))));
    const double dx = (w.x1 - w.x0) / nx, dy = (w.y1 - w.y0) / ny;
    std::vector<char> cov(size_t(nx) * ny, 0);
    auto edge = [](double ax, double ay, double bx, double by, double px, double py) {
        return (bx - ax) * (py - ay) - (by - ay) * (px - ax);
    };
    for (size_t k = 0; k + 2 < m.tri.size(); k += 3) {
        const Vec3 &a = m.pos[m.tri[k]], &b = m.pos[m.tri[k + 1]], &c = m.pos[m.tri[k + 2]];
        int ix0 = std::max(0, int((std::min({a.x, b.x, c.x}) - w.x0) / dx));
        int ix1 = std::min(nx - 1, int((std::max({a.x, b.x, c.x}) - w.x0) / dx));
        int iy0 = std::max(0, int((std::min({a.y, b.y, c.y}) - w.y0) / dy));
        int iy1 = std::min(ny - 1, int((std::max({a.y, b.y, c.y}) - w.y0) / dy));
        for (int iy = iy0; iy <= iy1; ++iy)
            for (int ix = ix0; ix <= ix1; ++ix) {
                double px = w.x0 + (ix + 0.5) * dx, py = w.y0 + (iy + 0.5) * dy;
                double e0 = edge(a.x, a.y, b.x, b.y, px, py), e1 = edge(b.x, b.y, c.x, c.y, px, py),
                       e2 = edge(c.x, c.y, a.x, a.y, px, py);
                if ((e0 >= 0 && e1 >= 0 && e2 >= 0) || (e0 <= 0 && e1 <= 0 && e2 <= 0)) cov[size_t(iy) * nx + ix] = 1;
            }
    }
    size_t open = 0;
    for (char c : cov) open += c == 0;
    return double(open) / cov.size();
}

} // namespace knit
