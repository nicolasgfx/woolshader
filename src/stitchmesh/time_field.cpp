#include "stitchmesh/time_field.h"

#include <chrono>
#include <cmath>
#include <limits>
#include <queue>

#include "model/cancel.h"

namespace knit {

namespace {

// The value at C from the planar front through A (ta) and B (tb) (Kimmel and Sethian 1998), or infinity when the front
// is not causal there. In the triangle's plane, with AB along x and C at (cx, cy), cy > 0: the front's unit direction g
// satisfies g . (B - A) = tb - ta, so gx = (tb - ta) / |AB|, and points into the triangle, gy = sqrt(1 - gx^2). Then
// t(C) = ta + g . (C - A). Causal when the characteristic through C, followed back along -g, crosses AB between A and B.
double planarFront(const Vec3& A, double ta, const Vec3& B, double tb, const Vec3& C) {
    const double inf = std::numeric_limits<double>::infinity();
    const Vec3 ab = B - A;
    const double c = norm(ab);
    if (!(c > 0)) return inf;
    const Vec3 ex = ab / c;
    const Vec3 ac = C - A;
    const double cx = dot(ac, ex);
    const double cy = norm(ac - ex * cx);   // C in the frame of AB: (cx, cy), cy >= 0
    if (!(cy > 0)) return inf;
    const double gx = (tb - ta) / c;
    if (std::fabs(gx) >= 1) return inf;     // |tb - ta| >= |AB|: the front runs along the edge, no planar update
    const double gy = std::sqrt(1.0 - gx * gx);
    const double foot = cx - gx * (cy / gy);   // where the characteristic through C meets the line of AB
    if (foot < 0 || foot > c) return inf;
    return ta + gx * cx + gy * cy;
}

struct Item {
    double t;
    int v;
    bool operator>(const Item& o) const { return t > o.t || (t == o.t && v > o.v); }
};

} // namespace

bool geodesicDistance(const TriMesh& m, const std::vector<int>& sources, std::vector<double>& t, TimeFieldReport& report,
                      std::string* err, const std::vector<double>* seedValues) {
    const auto start = std::chrono::steady_clock::now();
    report = TimeFieldReport{};
    const int n = (int) m.pos.size();
    if (sources.empty()) {
        if (err) *err = "time field: no source vertex";
        return false;
    }
    std::vector<std::vector<int>> trianglesOf((size_t) n);
    for (int f = 0; f < (int) m.tri.size(); ++f)
        for (int k = 0; k < 3; ++k) trianglesOf[(size_t) m.tri[f][k]].push_back(f);
    const double inf = std::numeric_limits<double>::infinity();
    t.assign((size_t) n, inf);
    std::vector<char> accepted((size_t) n, 0), fromTriangle((size_t) n, 0);
    std::priority_queue<Item, std::vector<Item>, std::greater<Item>> heap;
    if (seedValues && seedValues->size() != sources.size()) {
        if (err) *err = "time field: one seed value per source is needed";
        return false;
    }
    for (size_t i = 0; i < sources.size(); ++i) {
        const int s = sources[i];
        if (s < 0 || s >= n) {
            if (err) *err = "time field: source vertex " + std::to_string(s) + " is not in the mesh";
            return false;
        }
        const double v = seedValues ? (*seedValues)[i] : 0.0;
        t[(size_t) s] = std::min(t[(size_t) s], v);
        heap.push({t[(size_t) s], s});
    }
    while (!heap.empty()) {
        const Item top = heap.top();
        heap.pop();
        if (accepted[(size_t) top.v] || top.t > t[(size_t) top.v]) continue;
        accepted[(size_t) top.v] = 1;
        ++report.accepted;
        if ((report.accepted & 0xfff) == 0) cancelCheckpoint();
        if (top.t > 0) (fromTriangle[(size_t) top.v] ? report.triangleUpdates : report.edgeUpdates) += 1;
        for (int f : trianglesOf[(size_t) top.v]) {
            for (int k = 0; k < 3; ++k) {
                const int c = m.tri[f][k];
                if (accepted[(size_t) c]) continue;
                const int a = m.tri[f][(k + 1) % 3], b = m.tri[f][(k + 2) % 3];
                double best = t[(size_t) c];
                bool viaTriangle = false;
                // edge updates from every accepted corner
                if (accepted[(size_t) a]) best = std::min(best, t[(size_t) a] + norm(m.pos[c] - m.pos[a]));
                if (accepted[(size_t) b]) best = std::min(best, t[(size_t) b] + norm(m.pos[c] - m.pos[b]));
                if (accepted[(size_t) a] && accepted[(size_t) b]) {
                    const double front = planarFront(m.pos[a], t[(size_t) a], m.pos[b], t[(size_t) b], m.pos[c]);
                    if (front < best) {
                        best = front;
                        viaTriangle = true;
                    }
                }
                if (best < t[(size_t) c]) {
                    t[(size_t) c] = best;
                    fromTriangle[(size_t) c] = viaTriangle ? 1 : 0;
                    heap.push({best, c});
                }
            }
        }
    }
    report.ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
    return true;
}

} // namespace knit
