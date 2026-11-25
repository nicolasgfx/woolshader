// The border ring's corners placed on their boundary loop for the smallest side-to-loop distance (border_placement.h).
#include "stitchmesh/border_placement.h"

#include <algorithm>
#include <cmath>

namespace knit {

namespace {

// Points of the closed loop at equal arc steps of about `step` (the loop's length over a whole number of steps), in
// walking order from its first point. spacing: the step used.
std::vector<Vec3> pointsAlong(const std::vector<Vec3>& loop, double step, double& spacing) {
    const size_t n = loop.size();
    double length = 0;
    for (size_t i = 0; i < n; ++i) length += norm(loop[(i + 1) % n] - loop[i]);
    const int count = std::max(3, static_cast<int>(std::round(length / step)));
    spacing = length / count;
    std::vector<Vec3> out;
    double walked = 0;   // arc length at loop[i]
    size_t i = 0;
    for (int k = 0; k < count; ++k) {
        const double target = k * spacing;
        while (i + 1 < n && walked + norm(loop[i + 1] - loop[i]) < target) {
            walked += norm(loop[i + 1] - loop[i]);
            ++i;
        }
        const Vec3 a = loop[i], b = loop[(i + 1) % n];
        const double segment = norm(b - a);
        out.push_back(segment > 0 ? a + (b - a) * std::clamp((target - walked) / segment, 0.0, 1.0) : a);
    }
    return out;
}

double pointToSegment(const Vec3& p, const Vec3& a, const Vec3& b) {
    const Vec3 d = b - a;
    const double dd = dot(d, d);
    const double t = dd > 0 ? std::clamp(dot(p - a, d) / dd, 0.0, 1.0) : 0.0;
    return norm(p - (a + d * t));
}

// The side from point i over k points: the larger of (its arc's points to the chord) and (the chord's points, a quarter
// millimetre apart, to the arc's points).
double sideDistance(const std::vector<Vec3>& pts, int i, int k) {
    const int n = static_cast<int>(pts.size());
    const Vec3 a = pts[i], b = pts[(i + k) % n];
    double worst = 0;
    for (int j = 1; j < k; ++j) worst = std::max(worst, pointToSegment(pts[(i + j) % n], a, b));
    const int parts = std::max(1, static_cast<int>(std::ceil(norm(b - a) / 0.25)));
    for (int m = 1; m < parts; ++m) {
        const Vec3 c = a + (b - a) * (static_cast<double>(m) / parts);
        double nearest = 1e300;
        for (int j = 0; j <= k; ++j) nearest = std::min(nearest, norm(c - pts[(i + j) % n]));
        worst = std::max(worst, nearest);
    }
    return worst;
}

int nearestPoint(const std::vector<Vec3>& pts, const Vec3& p) {
    int best = 0;
    double bestDistance = 1e300;
    for (size_t i = 0; i < pts.size(); ++i) {
        const double d = norm(pts[i] - p);
        if (d < bestDistance) {
            bestDistance = d;
            best = static_cast<int>(i);
        }
    }
    return best;
}

} // namespace

BorderRingPlacement placeBorderRing(const std::vector<Vec3>& loop, const std::vector<Vec3>& ring, const BorderRingRequest& request) {
    BorderRingPlacement out;
    const int count = static_cast<int>(ring.size());
    if (count < 3 || loop.size() < 3) return out;
    double spacing = request.step;
    std::vector<Vec3> pts = pointsAlong(loop, request.step, spacing);
    int n = static_cast<int>(pts.size());
    if (n < 2 * count) return out;
    // the ring's corners as points of the loop; the loop is walked the way the ring runs
    std::vector<int> at(count);
    for (int k = 0; k < count; ++k) at[k] = nearestPoint(pts, ring[k]);
    // the direction: the sum of the shortest signed steps from each corner to the next (a relaxed ring can have two
    // corners in the wrong order, so the steps need not all point one way)
    long long turn = 0;
    for (int k = 0; k < count; ++k) {
        int step = ((at[(k + 1) % count] - at[k]) % n + n) % n;
        if (step > n / 2) step -= n;
        turn += step;
    }
    if (turn < 0) {
        std::reverse(pts.begin(), pts.end());
        for (int& i : at) i = n - 1 - i;
    }
    long long forward = 0;
    for (int k = 0; k < count; ++k) forward += ((at[(k + 1) % count] - at[k]) % n + n) % n;
    out.inOrder = forward == n;
    for (int k = 0; k < count && out.inOrder; ++k) {
        const int steps = ((at[(k + 1) % count] - at[k]) % n + n) % n;
        out.distanceBefore = std::max(out.distanceBefore, sideDistance(pts, at[k], steps == 0 ? n : steps));
    }
    const double mean = n * spacing / count;   // the loop's length over the corners
    const int reach = static_cast<int>(std::floor(request.anchorFreedom / spacing));
    for (double spread : request.spreads) {
        const int kMax = std::min(n - 1, static_cast<int>(std::ceil((1 + spread) * mean / spacing * 1.6)));
        // distance[i][k - 1] of the side from point i over k points; 1e300 where its chord is shorter or longer than the
        // spread allows
        std::vector<std::vector<double>> side(n, std::vector<double>(kMax, 1e300));
        for (int i = 0; i < n; ++i)
            for (int k = 1; k <= kMax; ++k) {
                const double chord = norm(pts[(i + k) % n] - pts[i]);
                if (chord < (1 - spread) * mean || chord > (1 + spread) * mean) continue;
                side[i][k - 1] = sideDistance(pts, i, k);
            }
        // the first corner where it is, then one point further each way in turn, as far as the request allows; a place
        // further away is taken only when it is closer to the loop
        bool placedHere = false;
        for (int offset = 0; offset <= 2 * reach; ++offset) {
            const int shift = offset % 2 == 1 ? (offset + 1) / 2 : -(offset / 2);
            const int anchor = ((at[0] + shift) % n + n) % n;
            std::vector<std::vector<double>> best(count + 1, std::vector<double>(n + 1, 1e300));
            std::vector<std::vector<int>> from(count + 1, std::vector<int>(n + 1, -1));
            best[0][0] = 0;
            for (int c = 0; c < count; ++c)
                for (int s = 0; s < n; ++s) {
                    if (best[c][s] > 1e299) continue;
                    for (int k = 1; k <= kMax && s + k <= n; ++k) {
                        const double d = std::max(best[c][s], side[(anchor + s) % n][k - 1]);
                        if (d < best[c + 1][s + k]) {
                            best[c + 1][s + k] = d;
                            from[c + 1][s + k] = s;
                        }
                    }
                }
            if (best[count][n] > 1e299) continue;
            if (out.placed && best[count][n] >= out.distance - 1e-9) continue;   // an earlier spread or a nearer first corner is at least as close
            std::vector<int> steps(count);
            int s = n;
            for (int c = count; c > 0; --c) {
                s = from[c][s];
                steps[c - 1] = s;
            }
            out.placed = true;
            placedHere = true;
            out.distance = best[count][n];
            out.spread = spread;
            out.firstCornerMoved = shift * spacing;
            out.corners.assign(count, Vec3{});
            for (int k = 0; k < count; ++k) out.corners[k] = pts[(anchor + steps[k]) % n];
            out.shortestSide = 1e300;
            out.longestSide = 0;
            for (int k = 0; k < count; ++k) {
                const double length = norm(out.corners[(k + 1) % count] - out.corners[k]);
                out.shortestSide = std::min(out.shortestSide, length);
                out.longestSide = std::max(out.longestSide, length);
            }
        }
        (void) placedHere;
        if (out.placed && out.distance <= request.limit) break;   // the tightest spread that meets the limit
    }
    return out;
}

} // namespace knit
