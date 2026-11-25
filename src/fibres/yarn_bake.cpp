// The bake of one yarn's fibres into maps on its tube; see yarn_bake.h. In reading order: the fibres back
// in nominal space, a grid over the chart, the rays, the texels, the height field's AO, the atlas.
#include "fibres/yarn_bake.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>

#include "knitting/tube/yarn_surface.h"
#include "model/cancel.h"
#include "model/colour.h"
#include "model/yarn_path.h"

namespace knit {
namespace {

constexpr double kGridCell = 0.1;       // [mm] the chart grid's cell along s and around (as arc length)
constexpr double kRayStart = 1.02;      // x r: the rays start just outside the nominal yarn
constexpr double kFootprintMargin = 1.2;// x radius: a segment's grid footprint, for the chord's dip towards the axis
constexpr double kAoSmoothFibres = 2;   // fibre diameters: the box the heights are smoothed over before the horizon scan
constexpr double kEnvelopeFibres = 3;   // fibre diameters: the box the fibre envelope (the local highest hit) is taken over
constexpr int kAoMaxSteps = 16;         // samples per direction: a reach of more texels is sampled with a stride (cost)
constexpr int kAoDirections = 8;
constexpr double kAxisMargin = 1e-6;    // [mm] a fibre point this close to the axis covers every angle
constexpr double kGapFloor = 0.5;       // x r: where a ray that met no fibre and no core ends (about the ply centres)

// A fibre segment in nominal space: a round linear tube from a to b whose radius runs from ra to rb.
struct NominalSegment {
    Vec3 a, b;
    float ra = 0, rb = 0;
};

// The yarn's fibres on its chart: the segments and a grid over (s, theta) listing which segments can be
// met by a radial ray at a cell (CSR: cellStart[cell] .. cellStart[cell + 1] into items).
struct ChartGrid {
    std::vector<NominalSegment> segments;
    int cellsAlong = 0, cellsAround = 0;
    std::vector<uint32_t> cellStart, items;
};

// The texel's rays summed: the hits' tangent-space normals, heights and counts, before the maps.
struct TexelSum {
    Vec3 normal;        // sum of unit tangent-space normals
    double height = 0;  // sum of rho / r
    int hits = 0;       // rays that met a fibre
    int coreHits = 0;   // rays that met a ply's core through a gap
};

// What one radial ray met: a fibre, else a ply's core, else nothing (the floor, between the plies).
struct RayHit {
    bool fibre = false, core = false;
    double rho = 0;     // [mm] distance from the axis
    Vec3 normal;        // unit, nominal space: the fibre's or the core's normal; at the floor the nearest ply's flank
    Vec3 direction;     // the fibre's unit direction in nominal space
};

// The plies along one yarn: per texel row i along s, the axis of every ply at the row's centre
// (axes[i * plies + ply], fibres/yarn_body.h bodyPlyAxes). Each ply's core is a round cylinder of radius
// rCore about its axis (none when rCore <= 0); a ray that reaches the floor takes its nearest ply's flank.
struct PlyCores {
    int plies = 0;
    double rCore = 0;  // [mm] BodyLayout::rCore
    std::vector<PlyAxis> axes;
};

double wrapAngle(double theta) {
    theta = std::fmod(theta, 2 * kPi);
    return theta < 0 ? theta + 2 * kPi : theta;
}

// The inverse of a symmetric 2 x 2 (m11, m12, m22).
std::array<double, 3> inverse2(const std::array<double, 3>& m) {
    const double det = m[0] * m[2] - m[1] * m[1];
    return {m[2] / det, -m[1] / det, m[0] / det};
}

Vec2 apply2(const std::array<double, 3>& m, Vec2 v) { return {m[0] * v.x + m[1] * v.y, m[1] * v.x + m[2] * v.y}; }

// --- the fibres back in nominal space

// A body point X at arc length s back to the nominal section: (x, y) = M(s)^-1 [n1 . (X - c), n2 . (X - c)].
Vec3 toNominal(const CurveSet& cs, const YarnSurface& ys, const SectionTable& section, const Vec3& X, double s) {
    const Ring g = ringAt(cs, ys, s);
    const Vec3 d = X - g.c;
    const Vec2 xy = apply2(inverse2(section.at(s)), {dot(d, g.n1), dot(d, g.n2)});
    return {xy.x, xy.y, s};
}

std::vector<NominalSegment> nominalSegments(const CurveSet& cs, const YarnSurface& ys, const BodyStrands& body, int yarn) {
    const SectionTable section = sectionTable(ys);
    std::vector<NominalSegment> out;
    std::vector<Vec3> nominal;
    for (size_t k = 0; k < body.strands(); ++k) {
        if (body.yarnOf[k] != yarn) continue;
        const uint64_t first = body.first[k], end = body.first[k + 1];
        nominal.clear();
        for (uint64_t v = first; v < end; ++v) nominal.push_back(toNominal(cs, ys, section, body.points[v], body.s[v]));
        for (uint64_t v = first; v + 1 < end; ++v) {
            NominalSegment seg;
            seg.a = nominal[v - first];
            seg.b = nominal[v + 1 - first];
            seg.ra = static_cast<float>(body.radius[v]);
            seg.rb = static_cast<float>(body.radius[v + 1]);
            out.push_back(seg);
        }
    }
    return out;
}

// --- the grid over the chart

// The cells a segment can be met in: s within its ends plus the radius; theta within the shorter arc
// between its ends widened by the angle its radius spans at its smaller distance from the axis. A
// segment at the axis covers every angle.
struct CellRange {
    int s0 = 0, s1 = -1;
    int t0 = 0, t1 = -1;  // t1 may exceed cellsAround - 1: wrap when visiting
};

CellRange footprint(const NominalSegment& seg, const ChartGrid& grid, double length) {
    const double rad = kFootprintMargin * std::max(seg.ra, seg.rb);
    CellRange c;
    const double sLo = std::clamp(std::min(seg.a.z, seg.b.z) - rad, 0.0, length), sHi = std::clamp(std::max(seg.a.z, seg.b.z) + rad, 0.0, length);
    c.s0 = std::min(grid.cellsAlong - 1, static_cast<int>(sLo / kGridCell));
    c.s1 = std::min(grid.cellsAlong - 1, static_cast<int>(sHi / kGridCell));
    const double rhoA = std::hypot(seg.a.x, seg.a.y), rhoB = std::hypot(seg.b.x, seg.b.y), rhoMin = std::min(rhoA, rhoB);
    if (rhoMin <= rad + kAxisMargin) {
        c.t0 = 0;
        c.t1 = grid.cellsAround - 1;
        return c;
    }
    const double thA = std::atan2(seg.a.y, seg.a.x), thB = std::atan2(seg.b.y, seg.b.x);
    double span = thB - thA;
    while (span > kPi) span -= 2 * kPi;
    while (span < -kPi) span += 2 * kPi;
    const double widen = std::asin(std::min(1.0, rad / rhoMin));
    const double lo = std::min(thA, thA + span) - widen, hi = std::max(thA, thA + span) + widen;
    const double cell = 2 * kPi / grid.cellsAround;
    c.t0 = static_cast<int>(std::floor(lo / cell));
    c.t1 = static_cast<int>(std::floor(hi / cell));
    if (c.t1 - c.t0 >= grid.cellsAround) c.t1 = c.t0 + grid.cellsAround - 1;
    return c;
}

int cellIndex(const ChartGrid& grid, int sCell, int tCell) {
    const int t = ((tCell % grid.cellsAround) + grid.cellsAround) % grid.cellsAround;
    return sCell * grid.cellsAround + t;
}

ChartGrid buildChartGrid(std::vector<NominalSegment> segments, double length, double r) {
    ChartGrid grid;
    grid.segments = std::move(segments);
    grid.cellsAlong = std::max(1, static_cast<int>(std::ceil(length / kGridCell)));
    grid.cellsAround = std::max(8, static_cast<int>(std::lround(2 * kPi * r / kGridCell)));
    const size_t cells = static_cast<size_t>(grid.cellsAlong) * grid.cellsAround;
    std::vector<uint32_t> counts(cells, 0);
    std::vector<CellRange> ranges(grid.segments.size());
    for (size_t k = 0; k < grid.segments.size(); ++k) {
        ranges[k] = footprint(grid.segments[k], grid, length);
        for (int s = ranges[k].s0; s <= ranges[k].s1; ++s)
            for (int t = ranges[k].t0; t <= ranges[k].t1; ++t) ++counts[cellIndex(grid, s, t)];
    }
    grid.cellStart.assign(cells + 1, 0);
    for (size_t c = 0; c < cells; ++c) grid.cellStart[c + 1] = grid.cellStart[c] + counts[c];
    grid.items.resize(grid.cellStart.back());
    std::vector<uint32_t> fill(grid.cellStart.begin(), grid.cellStart.end() - 1);
    for (size_t k = 0; k < grid.segments.size(); ++k)
        for (int s = ranges[k].s0; s <= ranges[k].s1; ++s)
            for (int t = ranges[k].t0; t <= ranges[k].t1; ++t) grid.items[fill[cellIndex(grid, s, t)]++] = static_cast<uint32_t>(k);
    return grid;
}

// --- the rays

// The ray o + t d against the round linear segment: the cylinder about its axis at the radius of the
// closest axis parameter (one refinement from the mean radius), then the end spheres. The nearest t >= 0.
bool hitSegment(const Vec3& o, const Vec3& d, const NominalSegment& seg, double& tHit, Vec3& normal) {
    const Vec3 ab = seg.b - seg.a;
    const double L2 = dot(ab, ab);
    if (L2 <= 0) return false;
    const Vec3 m = o - seg.a;
    const Vec3 w = m - ab * (dot(m, ab) / L2), v = d - ab * (dot(d, ab) / L2);
    double radius = 0.5 * (seg.ra + seg.rb);
    bool found = false;
    double best = 0;
    for (int pass = 0; pass < 2; ++pass) {
        const double A = dot(v, v), B = 2 * dot(w, v), C = dot(w, w) - radius * radius;
        if (A <= 0) break;
        const double disc = B * B - 4 * A * C;
        if (disc < 0) break;
        const double t = (-B - std::sqrt(disc)) / (2 * A);
        if (t < 0) break;
        const double u = dot(m + d * t, ab) / L2;
        if (u < 0 || u > 1) break;
        radius = seg.ra + (seg.rb - seg.ra) * u;
        found = true;
        best = t;
    }
    if (found) {
        const double u = std::clamp(dot(m + d * best, ab) / L2, 0.0, 1.0);
        normal = normalized((o + d * best) - (seg.a + ab * u));
        tHit = best;
        return true;
    }
    // The end spheres (the axis parameter fell outside the segment).
    for (int end = 0; end < 2; ++end) {
        const Vec3 c = end == 0 ? seg.a : seg.b;
        const double rr = end == 0 ? seg.ra : seg.rb;
        const Vec3 oc = o - c;
        const double B = 2 * dot(oc, d), C = dot(oc, oc) - rr * rr;
        const double disc = B * B - 4 * C;
        if (disc < 0) continue;
        const double t = (-B - std::sqrt(disc)) / 2;
        if (t < 0 || (found && t >= tHit)) continue;
        found = true;
        tHit = t;
        normal = normalized(o + d * t - c);
    }
    return found;
}

// The radial ray o + t d of the section plane (z = 0 relative to the row) against a ply's core: the round
// cylinder of radius rc about the axis through centre c along unit tangent a. A point P of the plane lies on
// it where |P - c|^2 - ((P - c) . a)^2 = rc^2, a quadratic in t. The nearest t >= 0, and the outward normal
// there: the part of P - c across the axis.
bool hitCore(Vec2 o, Vec2 d, const PlyAxis& axis, double rc, double& tHit, Vec3& normal) {
    if (rc <= 0) return false;
    const Vec2 w{o.x - axis.centre.x, o.y - axis.centre.y};
    const Vec2 a{axis.tangent.x, axis.tangent.y};  // (P - c) has no z: only the tangent's section part counts
    const double wa = dot(w, a), da = dot(d, a);
    const double A = dot(d, d) - da * da, B = 2 * (dot(w, d) - wa * da), C = dot(w, w) - wa * wa - rc * rc;
    const double disc = B * B - 4 * A * C;
    if (A <= 0 || disc < 0) return false;
    const double t = (-B - std::sqrt(disc)) / (2 * A);
    if (t < 0) return false;
    const Vec3 p{w.x + d.x * t, w.y + d.y * t, 0};
    tHit = t;
    normal = normalized(p - axis.tangent * dot(p, axis.tangent));
    return true;
}

// The floor point q of the section plane, between the plies: the normal of the nearest ply's flank there, the part
// of q - c across that ply's axis (the wall of the valley q lies in). On the yarn's axis of a single ply: e_rho.
Vec3 floorNormal(Vec2 q, const PlyAxis* axes, int plies, const Vec3& eRho) {
    int nearest = -1;
    double best = 0;
    for (int k = 0; k < plies; ++k) {
        const double d = std::hypot(q.x - axes[k].centre.x, q.y - axes[k].centre.y);
        if (nearest < 0 || d < best) {
            nearest = k;
            best = d;
        }
    }
    if (nearest < 0) return eRho;
    const Vec3 p{q.x - axes[nearest].centre.x, q.y - axes[nearest].centre.y, 0};
    const Vec3 across = p - axes[nearest].tangent * dot(p, axes[nearest].tangent);
    return norm(across) > 1e-9 ? normalized(across) : eRho;
}

// The radial ray at (s, theta): from kRayStart r towards the axis, the nearest fibre of its cell; through a
// gap, the nearest ply core of its row (`axes`: the row's plies); else the floor, with its ply's flank.
RayHit castRadialRay(const ChartGrid& grid, const PlyCores& cores, const PlyAxis* axes, double s, double theta, double r) {
    const Vec3 eRho{std::cos(theta), std::sin(theta), 0};
    const Vec3 o = eRho * (kRayStart * r) + Vec3{0, 0, s}, d = -eRho;
    const int sCell = std::min(grid.cellsAlong - 1, static_cast<int>(s / kGridCell));
    const int tCell = static_cast<int>(theta / (2 * kPi) * grid.cellsAround);
    const int cell = cellIndex(grid, sCell, tCell);
    RayHit hit;
    double tBest = (kRayStart - kGapFloor) * r;  // no fibre before the gap floor: a gap
    for (uint32_t k = grid.cellStart[cell]; k < grid.cellStart[cell + 1]; ++k) {
        const NominalSegment& seg = grid.segments[grid.items[k]];
        double t = 0;
        Vec3 n;
        if (!hitSegment(o, d, seg, t, n) || t >= tBest) continue;
        tBest = t;
        hit.fibre = true;
        hit.normal = n;
        hit.direction = normalized(seg.b - seg.a);
    }
    if (!hit.fibre) {
        const Vec2 o2{o.x, o.y}, d2{d.x, d.y};
        for (int k = 0; k < cores.plies; ++k) {
            double t = 0;
            Vec3 n;
            if (!hitCore(o2, d2, axes[k], cores.rCore, t, n) || t >= tBest) continue;
            tBest = t;
            hit.core = true;
            hit.normal = n;
        }
        if (!hit.core) hit.normal = floorNormal({o2.x + d2.x * tBest, o2.y + d2.y * tBest}, axes, cores.plies, eRho);
    }
    hit.rho = kRayStart * r - tBest;
    return hit;
}

// --- the texels

// A nominal normal in the tube's tangent space at (s, theta) on the compressed section: the normal
// transforms by M^-T (M symmetric: M^-1), the frame is T = e_s, B = M e_theta normalised, N = M^-1 e_rho
// normalised (orthogonal: e_theta . e_rho = 0 survives the map).
Vec3 tangentSpaceNormal(const Vec3& n, double theta, const std::array<double, 3>& M) {
    const std::array<double, 3> Minv = inverse2(M);
    const Vec2 eRho{std::cos(theta), std::sin(theta)}, eTheta{-std::sin(theta), std::cos(theta)};
    const Vec2 nxy = apply2(Minv, {n.x, n.y});
    Vec2 B = apply2(M, eTheta), N = apply2(Minv, eRho);
    B = B * (1 / norm(B));
    N = N * (1 / norm(N));
    return normalized(Vec3{n.z, dot(nxy, B), dot(nxy, N)});
}

// Every texel's rays summed, rows of s in parallel. Returns the per-yarn fibre direction sums (cos 2 psi,
// sin 2 psi) and the hit count through the arguments.
void castTexelRays(const ChartGrid& grid, const PlyCores& cores, const SectionTable& section, const BakeParams& p, double r,
                   const BakedYarnMaps& m, std::vector<TexelSum>& sums, double& cos2Sum, double& sin2Sum, long long& hits,
                   const ExecPolicy& exec) {
    const int n = p.samples;
    std::vector<double> rowCos(static_cast<size_t>(m.along), 0), rowSin(static_cast<size_t>(m.along), 0);
    std::vector<long long> rowHits(static_cast<size_t>(m.along), 0);
    parallelFor(m.along, exec.chunk > 0 ? exec.chunk : 64, [&](int64_t begin, int64_t end, int64_t) {
        cancelCheckpoint();  // a chunk is well under a second: a cancel never waits for the whole yarn
        for (int i = static_cast<int>(begin); i < static_cast<int>(end); ++i) {
            const double sCentre = (i + 0.5) * m.texel;
            const std::array<double, 3> M = section.at(sCentre);
            const PlyAxis* rowAxes = &cores.axes[static_cast<size_t>(i) * cores.plies];
            for (int j = 0; j < m.around; ++j) {
                TexelSum& sum = sums[static_cast<size_t>(i) * m.around + j];
                for (int a = 0; a < n; ++a)
                    for (int b = 0; b < n; ++b) {
                        const double s = (i + (a + 0.5) / n) * m.texel;
                        const double theta = 2 * kPi * (j + (b + 0.5) / n) / m.around;
                        const RayHit hit = castRadialRay(grid, cores, rowAxes, s, theta, r);
                        sum.height += hit.rho / r;
                        sum.normal += tangentSpaceNormal(hit.normal, theta, M);
                        if (hit.core) ++sum.coreHits;
                        if (!hit.fibre) continue;
                        ++sum.hits;
                        const Vec3 eTheta{-std::sin(theta), std::cos(theta), 0};
                        const double psi = std::atan2(dot(hit.direction, eTheta), hit.direction.z);
                        rowCos[static_cast<size_t>(i)] += std::cos(2 * psi);
                        rowSin[static_cast<size_t>(i)] += std::sin(2 * psi);
                        ++rowHits[static_cast<size_t>(i)];
                    }
            }
        }
    }, exec);
    for (int i = 0; i < m.along; ++i) {
        cos2Sum += rowCos[static_cast<size_t>(i)];
        sin2Sum += rowSin[static_cast<size_t>(i)];
        hits += rowHits[static_cast<size_t>(i)];
    }
}

// --- the height field's ambient occlusion

// The heights [mm] over a (2 half + 1)^2 box of texels, their mean (largest = false) or their largest: theta
// wraps, s is clamped at the yarn's ends.
std::vector<float> boxFiltered(const std::vector<float>& heightMm, const BakedYarnMaps& m, int half, bool largest) {
    if (half <= 0) return heightMm;
    std::vector<float> rows(heightMm.size()), out(heightMm.size());
    for (int i = 0; i < m.along; ++i)  // around first (theta wraps)
        for (int j = 0; j < m.around; ++j) {
            double sum = 0, top = -1e30;
            for (int k = -half; k <= half; ++k) {
                const double h = heightMm[static_cast<size_t>(i) * m.around + ((j + k) % m.around + m.around) % m.around];
                sum += h;
                top = std::max(top, h);
            }
            rows[static_cast<size_t>(i) * m.around + j] = static_cast<float>(largest ? top : sum / (2 * half + 1));
        }
    for (int i = 0; i < m.along; ++i)  // then along s
        for (int j = 0; j < m.around; ++j) {
            double sum = 0, top = -1e30;
            for (int k = -half; k <= half; ++k) {
                const double h = rows[static_cast<size_t>(std::clamp(i + k, 0, m.along - 1)) * m.around + j];
                sum += h;
                top = std::max(top, h);
            }
            out[static_cast<size_t>(i) * m.around + j] = static_cast<float>(largest ? top : sum / (2 * half + 1));
        }
    return out;
}

// Where the horizon scan looks, in texels: from `first` to `last` by `stride` (at most kAoMaxSteps samples).
struct AoReach {
    int first = 1, last = 1, stride = 1;
};

// The scan of a body: from just outside the smoothing box to the ply's radius, the scale of the valleys
// between the plies and of the gaps down to the cores.
AoReach aoReachOf(const BakedYarnMaps& m, int smoothHalf, double plyRadius) {
    AoReach a;
    a.first = smoothHalf + 1;
    a.last = std::max(a.first, static_cast<int>(std::ceil(plyRadius / m.texel)));
    a.stride = std::max(1, (a.last - a.first + kAoMaxSteps - 1) / kAoMaxSteps);
    return a;
}

// The horizon scan of the ply maps (ply_profile.cpp heightFieldAO), on the smoothed heights [mm]: in each
// direction the steepest elevation to a texel within the reach; AO = 1 - mean sin(elevation). theta
// wraps, s stops at the yarn's ends.
double horizonAO(const std::vector<float>& heightMm, const BakedYarnMaps& m, const AoReach& reach, int i, int j) {
    constexpr int dx[kAoDirections] = {1, 1, 0, -1, -1, -1, 0, 1}, dy[kAoDirections] = {0, 1, 1, 1, 0, -1, -1, -1};
    const double h0 = heightMm[static_cast<size_t>(i) * m.around + j];
    double occlusion = 0;
    for (int d = 0; d < kAoDirections; ++d) {
        double horizon = 0;
        for (int k = reach.first; k <= reach.last; k += reach.stride) {
            const int ii = i + k * dx[d];
            if (ii < 0 || ii >= m.along) break;
            const int jj = ((j + k * dy[d]) % m.around + m.around) % m.around;
            const double dist = k * m.texel * std::hypot(static_cast<double>(dx[d]), static_cast<double>(dy[d]));
            horizon = std::max(horizon, std::atan((heightMm[static_cast<size_t>(ii) * m.around + jj] - h0) / dist));
        }
        occlusion += std::sin(horizon);
    }
    return 1 - occlusion / kAoDirections;
}

// --- the atlas

uint8_t byteOf(double unit) { return static_cast<uint8_t>(std::clamp(std::lround(unit * 255), 0L, 255L)); }

void putRgb(std::vector<uint8_t>& map, const BakedYarnMaps& m, int x, int row, uint8_t r, uint8_t g, uint8_t b) {
    const size_t at = 3 * (static_cast<size_t>(row) * m.width + x);
    map[at] = r;
    map[at + 1] = g;
    map[at + 2] = b;
}

void putGrey(std::vector<uint8_t>& map, const BakedYarnMaps& m, int x, int row, double unit) {
    const uint8_t v = byteOf(unit);
    putRgb(map, m, x, row, v, v, v);
}

} // namespace

int atlasStrips(double length, const BakeParams& p) {
    const int along = std::max(1, static_cast<int>(std::ceil(length / p.texel)));
    return (along + p.atlasWidth - 1) / p.atlasWidth;
}

bool planAtlas(double length, double r, const BakeParams& p, BakedYarnMaps& out, std::string* err) {
    out = BakedYarnMaps{};
    if (!(p.texel > 0) || p.samples < 1 || p.atlasWidth < 1 || !(r > 0) || !(length > 0)) {
        if (err) *err = "yarn bake: texel, samples and atlas width must be positive on a yarn of positive length";
        return false;
    }
    out.texel = p.texel;
    out.around = std::max(8, static_cast<int>(std::lround(2 * kPi * r / p.texel)));
    out.along = std::max(1, static_cast<int>(std::ceil(length / p.texel)));
    out.width = p.atlasWidth;
    out.stripMm = p.texel * p.atlasWidth;
    out.strips = atlasStrips(length, p);
    out.height = out.strips * out.around;
    if (out.height > p.atlasWidth) {
        char line[200];
        std::snprintf(line, sizeof line, "yarn bake: %d strips of %d texels around (%d rows) exceed the atlas width %d: a coarser texel",
                      out.strips, out.around, out.height, p.atlasWidth);
        if (err) *err = line;
        return false;
    }
    return true;
}

void relabelAtlasUv(TubeMesh& tube, const BakeParams& p) {
    const double stripMm = p.texel * p.atlasWidth;
    for (size_t y = 0; y + 1 < tube.yarnFirstVertex.size(); ++y) {
        const uint32_t first = tube.yarnFirstVertex[y], end = tube.yarnFirstVertex[y + 1];
        double length = 0;
        for (uint32_t v = first; v < end; ++v) length = std::max(length, tube.s[v]);
        const int strips = atlasStrips(length, p);
        for (uint32_t v = first; v < end; ++v) {
            const double u = tube.uv[2 * v];
            const int strip = std::clamp(static_cast<int>(std::lround(tube.s[v] / stripMm - u)), 0, strips - 1);
            tube.uv[2 * v + 1] = static_cast<float>((strip + tube.theta[v] / (2 * kPi)) / strips);
        }
    }
}

void atlasPixel(const BakedYarnMaps& m, int i, int j, int& x, int& row) {
    const int strip = i / m.width;
    x = i % m.width;
    row = m.height - 1 - (strip * m.around + j);
}

bool chartTexel(const BakedYarnMaps& m, int x, int row, int& i, int& j) {
    const int fromBottom = m.height - 1 - row;
    const int strip = fromBottom / m.around;
    j = fromBottom % m.around;
    i = strip * m.width + x;
    return i < m.along;
}

std::vector<uint8_t> cavityMapOf(const std::vector<uint8_t>& ao, double cavity) {
    std::vector<uint8_t> out(ao.size());
    for (size_t k = 0; k < ao.size(); ++k) out[k] = byteOf(srgbEncode(1 - cavity * (1 - ao[k] / 255.0)));
    return out;
}

std::vector<uint8_t> roughnessMapOf(const std::vector<uint8_t>& spread, double baseRoughness) {
    std::vector<uint8_t> out(spread.size());
    for (size_t k = 0; k < spread.size(); ++k)
        out[k] = byteOf(std::clamp(std::sqrt(baseRoughness * baseRoughness + spread[k] / 255.0), baseRoughness, 1.0));
    return out;
}

bool bakeYarnMaps(const CurveSet& cs, const YarnSurface& ys, const BodyStrands& body, const YarnBodyParams& bodyParams, uint64_t partSeed,
                  int yarn, const BakeParams& p, BakedYarnMaps& out, std::string* err, const ExecPolicy& exec) {
    const auto t0 = std::chrono::steady_clock::now();
    const double r = ys.radius;
    if (!planAtlas(ys.length, r, p, out, err)) return false;
    std::vector<NominalSegment> segments = nominalSegments(cs, ys, body, yarn);
    if (segments.empty()) {
        if (err) *err = "yarn bake: the body holds no fibre of yarn " + std::to_string(yarn);
        return false;
    }
    const ChartGrid grid = buildChartGrid(std::move(segments), ys.length, r);
    const SectionTable section = sectionTable(ys);
    const BodyLayout layout = bodyLayout(bodyParams, r);
    PlyCores cores;  // the rows' ply axes: bodyPlyAxes gives one ply's rows, interleaved here per row
    cores.plies = bodyParams.plies;
    cores.rCore = layout.rCore;
    cores.axes.resize(static_cast<size_t>(out.along) * cores.plies);
    for (int k = 0; k < cores.plies; ++k) {
        const std::vector<PlyAxis> axes = bodyPlyAxes(bodyParams, r, partSeed, yarn, k, out.along * p.texel, p.texel);
        for (int i = 0; i < out.along; ++i) cores.axes[static_cast<size_t>(i) * cores.plies + k] = axes[static_cast<size_t>(i)];
    }
    std::vector<TexelSum> sums(static_cast<size_t>(out.along) * out.around);
    double cos2 = 0, sin2 = 0;
    long long hits = 0;
    castTexelRays(grid, cores, section, p, r, out, sums, cos2, sin2, hits, exec);
    const double rays = static_cast<double>(out.along) * out.around * p.samples * p.samples;
    out.coverage = hits / rays;
    long long coreHits = 0;
    for (const TexelSum& sum : sums) coreHits += sum.coreHits;
    out.cores = coreHits / rays;
    out.fibreAngleDeg = hits > 0 ? 0.5 * std::atan2(sin2, cos2) * 180 / kPi : 0;
    // The height field in mm (a gap counts as the floor), then the maps.
    const double samples2 = static_cast<double>(p.samples) * p.samples;
    std::vector<float> heightMm(sums.size());
    for (size_t k = 0; k < sums.size(); ++k) heightMm[k] = static_cast<float>(r * sums[k].height / samples2);
    const int smoothHalf = static_cast<int>(std::lround(0.5 * kAoSmoothFibres * 2 * bodyParams.fibreRadius / p.texel));
    const std::vector<float> smoothMm = boxFiltered(heightMm, out, smoothHalf, false);
    const AoReach reach = aoReachOf(out, smoothHalf, layout.rp);
    // The shell's transmittance: a hit below the fibre envelope (the highest hit in a box of kEnvelopeFibres fibre
    // diameters, averaged over the same box) gets its ambient light through the fibres above it, exp(-depth / l)
    // with l = pi d / (4 phi) the mean free path between fibres of diameter d at the shell's packing phi. A gap's
    // core lies deep under the envelope and reads dark, as the fibres around it shade it.
    const int envelopeHalf = static_cast<int>(std::lround(0.5 * kEnvelopeFibres * 2 * bodyParams.fibreRadius / p.texel));
    const std::vector<float> envelopeMm = boxFiltered(boxFiltered(heightMm, out, envelopeHalf, true), out, envelopeHalf, false);
    const double freePath = layout.phiShell > 0 ? kPi * 2 * bodyParams.fibreRadius / (4 * layout.phiShell) : 0;
    const size_t pixels = static_cast<size_t>(out.width) * out.height;
    out.normal.assign(3 * pixels, 0);
    out.specular.assign(3 * pixels, 0);
    out.ao.assign(3 * pixels, 255);
    out.spread.assign(3 * pixels, 255);
    out.heightMap.assign(3 * pixels, 0);
    for (size_t k = 0; k < pixels; ++k) {  // the unused tail: flat, unoccluded, no gloss, rough
        out.normal[3 * k] = out.normal[3 * k + 1] = 128;
        out.normal[3 * k + 2] = 255;
    }
    for (int i = 0; i < out.along; ++i)
        for (int j = 0; j < out.around; ++j) {
            const TexelSum& sum = sums[static_cast<size_t>(i) * out.around + j];
            int x = 0, row = 0;
            atlasPixel(out, i, j, x, row);
            const size_t k = static_cast<size_t>(i) * out.around + j;
            const double depth = std::max(0.0, static_cast<double>(envelopeMm[k] - heightMm[k]));
            const double transmittance = freePath > 0 ? std::exp(-depth / freePath) : 1.0;
            putGrey(out.ao, out, x, row, horizonAO(smoothMm, out, reach, i, j) * transmittance);
            const Vec3 n = normalized(sum.normal);
            putRgb(out.normal, out, x, row, byteOf(0.5 * (n.x + 1)), byteOf(0.5 * (n.y + 1)), byteOf(0.5 * (n.z + 1)));
            putGrey(out.specular, out, x, row, sum.hits / samples2);
            // Toksvig 2005: the mean normal's length over the samples says how far they spread.
            const double meanLength = std::max(1e-3, norm(sum.normal) / samples2);
            putGrey(out.spread, out, x, row, std::min(1.0, (1 - meanLength) / meanLength));
            putGrey(out.heightMap, out, x, row, heightMm[static_cast<size_t>(i) * out.around + j] / r);
        }
    out.ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    return true;
}

} // namespace knit
