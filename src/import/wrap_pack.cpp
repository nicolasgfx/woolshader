// Stage 10: the flattened charts laid out in the unit square. Each chart (an island of uv) is turned so that its smallest
// enclosing rectangle stands with the long side along v, then the rectangles, each with a gutter, are placed corner to corner
// by a greedy packer that keeps the layout square, and one scale maps the layout into [0, 1].
#include "import/wrap_stages.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace knit {

namespace {

constexpr double kGutter = 0.002;            // the gutter on each side of an island, as a share of sqrt(sum of rectangle areas)
constexpr double kEqualAreas = 1e-12;        // rectangle areas this close (relative) are equal: the smaller angle wins
constexpr double kEqualSides = 1e-6;         // rectangle sides this close (relative) are equal: the smaller turn wins

// The faces of one chart, ascending.
struct Island {
    std::vector<int> faces;
};

// The charts of the flagged faces: faces joined across edges that are no seam. Numbered by their lowest face.
std::vector<Island> islandsOf(const FaceGraph& graph, const std::vector<char>& seam, const std::vector<char>& faces) {
    std::vector<Island> islands;
    std::vector<char> reached(faces.size(), 0);
    for (int seed = 0; seed < graph.faceCount(); ++seed) {
        if (!faces[seed] || reached[seed]) continue;
        Island island;
        reached[seed] = 1;
        island.faces.push_back(seed);
        for (size_t next = 0; next < island.faces.size(); ++next) {
            const int f = island.faces[next];
            for (int k = graph.start[f]; k < graph.start[f + 1]; ++k) {
                const FaceNeighbour& n = graph.neighbour[k];
                if (seam[n.edge] || !faces[n.face] || reached[n.face]) continue;
                reached[n.face] = 1;
                island.faces.push_back(n.face);
            }
        }
        std::sort(island.faces.begin(), island.faces.end());
        islands.push_back(island);
    }
    return islands;
}

bool pointComesFirst(const Vec2& a, const Vec2& b) { return a.x != b.x ? a.x < b.x : a.y < b.y; }
bool samePoint(const Vec2& a, const Vec2& b) { return a.x == b.x && a.y == b.y; }

// Twice the signed area of the triangle (o, a, b): positive when a to b turns counter-clockwise about o.
double turn(const Vec2& o, const Vec2& a, const Vec2& b) { return cross(a - o, b - o); }

// The convex hull, counter-clockwise, without collinear points: Andrew's monotone chain (Andrew 1979). Fewer than three points
// when all points lie on one line.
std::vector<Vec2> convexHull(std::vector<Vec2> points) {
    std::sort(points.begin(), points.end(), pointComesFirst);
    points.erase(std::unique(points.begin(), points.end(), samePoint), points.end());
    if (points.size() < 3) return points;
    std::vector<Vec2> hull(points.size() * 2);
    size_t size = 0;
    for (size_t i = 0; i < points.size(); ++i) {  // the lower chain, left to right
        while (size >= 2 && turn(hull[size - 2], hull[size - 1], points[i]) <= 0) --size;
        hull[size++] = points[i];
    }
    const size_t lower = size + 1;
    for (size_t i = points.size() - 1; i-- > 0;) {  // the upper chain, right to left
        while (size >= lower && turn(hull[size - 2], hull[size - 1], points[i]) <= 0) --size;
        hull[size++] = points[i];
    }
    hull.resize(size - 1);  // the last point repeats the first
    return hull;
}

// An enclosing rectangle: its first side runs along (cos angle, sin angle) with length `width`, the second along
// (-sin angle, cos angle) with length `height`.
struct Rectangle {
    double angle = 0;  // in [0, pi / 2)
    double width = 0, height = 0;
    Vec2 centre;
};

// The angle of a direction, reduced to [0, pi / 2): a rectangle is the same after a quarter turn.
double quarterTurnAngle(const Vec2& direction) {
    const double quarter = 0.5 * kPi;
    double angle = std::fmod(std::atan2(direction.y, direction.x), quarter);
    if (angle < 0) angle += quarter;
    if (angle >= quarter) angle = 0;
    return angle;
}

// The rectangle at `angle` that just encloses the points.
Rectangle rectangleAt(const std::vector<Vec2>& points, double angle) {
    const Vec2 first(std::cos(angle), std::sin(angle));
    const Vec2 second(-first.y, first.x);
    double lowA = 1e300, highA = -1e300, lowB = 1e300, highB = -1e300;
    for (const Vec2& p : points) {
        lowA = std::min(lowA, dot(p, first));
        highA = std::max(highA, dot(p, first));
        lowB = std::min(lowB, dot(p, second));
        highB = std::max(highB, dot(p, second));
    }
    Rectangle r;
    r.angle = angle;
    r.width = highA - lowA;
    r.height = highB - lowB;
    r.centre = first * (0.5 * (lowA + highA)) + second * (0.5 * (lowB + highB));
    return r;
}

// The index after i on a closed polygon of n points.
size_t nextOf(size_t i, size_t n) { return i + 1 == n ? 0 : i + 1; }

// The smallest-area enclosing rectangle of a convex polygon. One of its sides lies on a hull edge (Freeman and Shapira 1975), and
// rotating calipers find, for each edge in turn, the extreme points along the edge, against it and across it (Toussaint 1983):
// each caliper only moves forward, so the whole search is linear in the hull's size.
Rectangle smallestRectangle(const std::vector<Vec2>& hull) {
    if (hull.empty()) return Rectangle{};
    if (hull.size() == 1) return rectangleAt(hull, 0.0);
    if (hull.size() == 2) return rectangleAt(hull, quarterTurnAngle(hull[1] - hull[0]));
    const size_t n = hull.size();
    size_t along = 0, across = 0, behind = 0;  // the calipers: farthest along the edge, from it, and back along it
    double bestArea = -1, bestAngle = 0;
    for (size_t i = 0; i < n; ++i) {
        const Vec2 edge = hull[nextOf(i, n)] - hull[i];
        const double length = norm(edge);
        if (!(length > 0)) continue;
        const Vec2 u = edge * (1.0 / length);
        const Vec2 v(-u.y, u.x);  // into the hull: it lies on the left of its counter-clockwise edges
        if (i == 0) {
            for (size_t k = 0; k < n; ++k) {
                if (dot(hull[k], u) > dot(hull[along], u)) along = k;
                if (dot(hull[k], v) > dot(hull[across], v)) across = k;
                if (dot(hull[k], u) < dot(hull[behind], u)) behind = k;
            }
        } else {
            for (size_t step = 0; step < n && dot(hull[nextOf(along, n)], u) > dot(hull[along], u); ++step) along = nextOf(along, n);
            for (size_t step = 0; step < n && dot(hull[nextOf(across, n)], v) > dot(hull[across], v); ++step) across = nextOf(across, n);
            for (size_t step = 0; step < n && dot(hull[nextOf(behind, n)], u) < dot(hull[behind], u); ++step) behind = nextOf(behind, n);
        }
        const double area = (dot(hull[along], u) - dot(hull[behind], u)) * (dot(hull[across], v) - dot(hull[i], v));
        const double angle = quarterTurnAngle(u);
        const bool equal = bestArea >= 0 && std::fabs(area - bestArea) <= kEqualAreas * std::max(area, bestArea);
        if (bestArea < 0 || (!equal && area < bestArea) || (equal && angle < bestAngle)) {
            bestArea = area;
            bestAngle = angle;
        }
    }
    return rectangleAt(hull, bestAngle);
}

// The turn that brings the rectangle's longer side onto the v axis. Each side can be brought up by two turns half a turn apart;
// the smaller one is taken, and with equal sides the smaller over both sides.
double turnToStandUp(const Rectangle& r) {
    const double raiseWidth = 0.5 * kPi - r.angle;  // the first side, at `angle`, onto v
    const double raiseHeight = -r.angle;           // the second side, at angle + pi / 2, onto v
    const bool equal = std::fabs(r.width - r.height) <= kEqualSides * std::max(r.width, r.height);
    if (equal) return std::fabs(raiseWidth) < std::fabs(raiseHeight) ? raiseWidth : raiseHeight;
    return r.width > r.height ? raiseWidth : raiseHeight;
}

Vec2 rotated(const Vec2& p, double cosine, double sine) { return {cosine * p.x - sine * p.y, sine * p.x + cosine * p.y}; }

// An island after its turn: uv' = R(turn)(uv - centre), whose bounding box is [low, high].
struct TurnedIsland {
    double turn = 0;
    Vec2 centre;
    Vec2 low, high;
    double rectangleArea = 0;
};

TurnedIsland turnIsland(const Island& island, const QuadUv& uv) {
    std::vector<Vec2> points;
    points.reserve(island.faces.size() * 4);
    for (int f : island.faces)
        for (const Vec2& p : uv[f]) points.push_back(p);
    const Rectangle rectangle = smallestRectangle(convexHull(points));
    TurnedIsland turned;
    turned.turn = turnToStandUp(rectangle);
    turned.centre = rectangle.centre;
    turned.rectangleArea = rectangle.width * rectangle.height;
    const double cosine = std::cos(turned.turn), sine = std::sin(turned.turn);
    turned.low = {1e300, 1e300};
    turned.high = {-1e300, -1e300};
    for (const Vec2& p : points) {
        const Vec2 q = rotated(p - turned.centre, cosine, sine);
        turned.low = {std::min(turned.low.x, q.x), std::min(turned.low.y, q.y)};
        turned.high = {std::max(turned.high.x, q.x), std::max(turned.high.y, q.y)};
    }
    return turned;
}

// A box of the layout, the gutter included: [u0, u1] x [v0, v1].
struct LayoutBox {
    double u0 = 0, v0 = 0, u1 = 0, v1 = 0;
};

bool boxesOverlap(const LayoutBox& a, const LayoutBox& b) { return a.u0 < b.u1 && b.u0 < a.u1 && a.v0 < b.v1 && b.v0 < a.v1; }

struct BoxToPlace {
    int island = 0;
    double width = 0, height = 0;
};

bool placedEarlier(const BoxToPlace& a, const BoxToPlace& b) {
    if (a.height != b.height) return a.height > b.height;
    if (a.width != b.width) return a.width > b.width;
    return a.island < b.island;
}

// Places the boxes, tallest first, each at the free position that keeps the layout's larger extent smallest. The candidate
// points are the origin and three corners of every box placed, (u1, v0), (u0, v1) and (u1, v1); a box may sit with any of its four
// corners on a candidate. Equal extents keep the earlier candidate. Returns each island's box.
std::vector<LayoutBox> layOut(std::vector<BoxToPlace> boxes) {
    std::sort(boxes.begin(), boxes.end(), placedEarlier);
    std::vector<LayoutBox> boxOfIsland(boxes.size());
    std::vector<LayoutBox> placed;
    double extentU = 0, extentV = 0;
    for (const BoxToPlace& box : boxes) {
        std::vector<Vec2> candidates = {{0.0, 0.0}};
        for (const LayoutBox& p : placed) {
            candidates.push_back({p.u1, p.v0});
            candidates.push_back({p.u0, p.v1});
            candidates.push_back({p.u1, p.v1});
        }
        bool found = false;
        LayoutBox best;
        double bestExtent = 0;
        for (const Vec2& c : candidates) {
            const LayoutBox corners[4] = {{c.x, c.y, c.x + box.width, c.y + box.height},
                                          {c.x - box.width, c.y, c.x, c.y + box.height},
                                          {c.x, c.y - box.height, c.x + box.width, c.y},
                                          {c.x - box.width, c.y - box.height, c.x, c.y}};
            for (const LayoutBox& option : corners) {
                if (option.u0 < 0 || option.v0 < 0) continue;
                const double extent = std::max(std::max(extentU, option.u1), std::max(extentV, option.v1));
                if (found && extent >= bestExtent) continue;
                bool free = true;
                for (const LayoutBox& p : placed)
                    if (boxesOverlap(option, p)) {
                        free = false;
                        break;
                    }
                if (!free) continue;
                found = true;
                best = option;
                bestExtent = extent;
            }
        }
        // Always found: the origin when nothing is placed, else (u1, v0) of the box reaching farthest along u.
        placed.push_back(best);
        boxOfIsland[box.island] = best;
        extentU = std::max(extentU, best.u1);
        extentV = std::max(extentV, best.v1);
    }
    return boxOfIsland;
}

} // namespace

void packCharts(const FaceGraph& graph, const std::vector<char>& seam, const std::vector<char>& faces, QuadUv& uv) {
    const std::vector<Island> islands = islandsOf(graph, seam, faces);
    if (islands.empty()) return;
    std::vector<TurnedIsland> turned;
    turned.reserve(islands.size());
    double areaSum = 0;
    for (const Island& island : islands) {
        turned.push_back(turnIsland(island, uv));
        areaSum += turned.back().rectangleArea;
    }
    const double gutter = kGutter * std::sqrt(areaSum);

    std::vector<BoxToPlace> boxes(islands.size());
    for (size_t i = 0; i < islands.size(); ++i) {
        boxes[i].island = static_cast<int>(i);
        boxes[i].width = turned[i].high.x - turned[i].low.x + 2.0 * gutter;
        boxes[i].height = turned[i].high.y - turned[i].low.y + 2.0 * gutter;
    }
    const std::vector<LayoutBox> placed = layOut(boxes);
    double extent = 0;
    for (const LayoutBox& box : placed) extent = std::max({extent, box.u1, box.v1});
    const double scale = extent > 0 ? 1.0 / extent : 1.0;

    // uv' = R(turn)(uv - centre) * scale + t, with t putting the island's turned box inside its gutter.
    for (size_t i = 0; i < islands.size(); ++i) {
        const TurnedIsland& island = turned[i];
        const double cosine = std::cos(island.turn), sine = std::sin(island.turn);
        const Vec2 t((placed[i].u0 + gutter - island.low.x) * scale, (placed[i].v0 + gutter - island.low.y) * scale);
        for (int f : islands[i].faces)
            for (Vec2& p : uv[f]) {
                const Vec2 q = rotated(p - island.centre, cosine, sine);
                p = {std::clamp(q.x * scale + t.x, 0.0, 1.0), std::clamp(q.y * scale + t.y, 0.0, 1.0)};
            }
    }
}

} // namespace knit
