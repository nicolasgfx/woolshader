// The fur field: cached candidate roots and one disc per fibre (fur_field.h has the method). The expressions
// keep the prototype's operation order (feature/furry-teddy fur.cpp sampleFollicles), so the frozen baseline's
// fibres come out bit for bit; reordering a product changes the last bits and moves every fibre.
#include "fibres/fur_field.h"

#include <algorithm>
#include <cmath>
#include <unordered_map>

#include "model/cancel.h"
#include "model/hash.h"
#include "model/noise.h"

namespace knit {

namespace {

constexpr double kDeg = kPi / 180;
constexpr double kHiddenTestLift = 1e-3;  // [mm] the hidden test's point above the surface
const Vec3 kComb{0, 0, -1};               // the combing direction (library, Z up): the pile leans downward

// Hash channels (fur_field.h). Never repurposed: a changed meaning moves every fibre.
enum RootChannel {
    kCountChannel = 0, kPlacementChannel1 = 1, kPlacementChannel2 = 2, kCoverageChannel = 3, kLengthChannel = 4,
    kLengthChannel2 = 5, kLiftChannel = 6, kCombChannel = 7, kTangentChannel = 8, kKeyChannel = 9, kPlaneChannel = 10,
    kRadiusChannel = 11, kThicknessChannel = 15, kLengthShareChannel = 16, kBendChannel = 17, kDensityChannel = 18
};
// The irregularities' own seeds: clump centres, islands, the patch noise.
constexpr uint64_t kClumpSalt = 0x636c755fULL, kIslandSalt = 0x69736c5fULL, kPatchSalt = 0x7061745fULL;  // "clu_", "isl_", "pat_"
constexpr int kPatchOctaves = 3;
constexpr double kWearRadius = 5;  // [mm] a place curved tighter than this radius counts as fully exposed (wear e = 1)

constexpr uint64_t kFibresSalt = 0x746f705fULL;  // "top_": the fine fibres' own roots

// The start direction at a root with face normal n: the comb on the tangent plane (when the comb is along n,
// a tangent turned by fallbackTurn), turned by `turn` about n, lifted by `lift` [rad] above the tangent plane.
Vec3 furDirection(const Vec3& n, const Vec3& comb, double fallbackTurn, double turn, double lift) {
    Vec3 ct = comb - n * dot(comb, n);
    if (norm(ct) < 1e-9) {
        const Vec3 t1 = normalized(cross(n, std::fabs(n.x) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0}));
        ct = rotateAbout(t1, n, fallbackTurn);
    }
    ct = rotateAbout(normalized(ct), n, turn);
    return normalized(ct * std::cos(lift) + n * std::sin(lift));
}

// The direction toward a loop's centre: in the plane across d, at angle `phase` from e1 (the normal's part
// across d, away from the skin), flipped so that the loop curls away from the skin (c . n >= 0).
Vec3 discCentreDirection(const Vec3& n, const Vec3& d, double phase) {
    Vec3 e1 = n - d * dot(n, d);
    if (norm(e1) < 1e-9) e1 = cross(d, std::fabs(d.x) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0});
    e1 = normalized(e1);
    const Vec3 e2 = cross(d, e1);
    Vec3 c = e1 * std::cos(phase) + e2 * std::sin(phase);
    if (dot(c, n) < 0) c = c * -1.0;
    return c;
}

// A loop fibre's disc at a root (the draws: the cache seed's channels of this candidate).
FurDisc loopDisc(const FurRoot& r, uint64_t seed, const FurLook& look) {
    const FurLoops& l = look.loops;
    const FurVariation& var = look.variation;
    const int t = r.tri, i = r.index;
    const auto draw = [&](int channel) { return hash01(seed, t, i, channel); };
    FurDisc f;
    f.root = r.position;
    f.normal = r.normal;
    f.tri = t;
    f.index = i;
    f.budgetKey = r.budgetKey;
    f.layer = FurLayer::Loops;
    // Direction: the comb on the tangent plane, turned about the normal, lifted.
    const double v = std::clamp(var.direction, 0.0, 1.0);
    const double combJitterDeg = 180 * v, liftDeg = 90 - 50 * v, liftJitterDeg = 40 * v;
    const double turn = (2 * draw(kCombChannel) - 1) * combJitterDeg * kDeg;
    const double lift = (liftDeg + (2 * draw(kLiftChannel) - 1) * liftJitterDeg) * kDeg;
    f.d = furDirection(f.normal, kComb, 2 * kPi * draw(kTangentChannel), turn, lift);
    // Thickness: the root diameter +- its share; the tip keeps the taper.
    const double r0 = 0.5 * l.rootDiameter, r1 = r0 * l.taper;
    const double thick = var.thickness > 0 ? 1 + std::clamp(var.thickness, 0.0, 1.0) * (2 * draw(kThicknessChannel) - 1) : 1.0;
    f.rootRadius = r0 * thick;
    f.tipRadius = r1 * thick;
    const double lengthShare = 1 - std::clamp(var.length, 0.0, 1.0) * draw(kLengthShareChannel);
    if (!(l.curl > 0)) {
        f.straightLength = l.straightLength * lengthShare;
        return f;
    }
    // The loop: its plane turned about d at random, its radius +- half its share, its rise +- its share.
    const double phase = 2 * kPi * draw(kPlaneChannel);
    f.radius = l.curl * (1 + 0.5 * std::clamp(var.curl, 0.0, 1.0) * (2 * draw(kRadiusChannel) - 1));
    f.rise = l.bend * (1 + std::clamp(var.bend, 0.0, 1.0) * (2 * draw(kBendChannel) - 1));
    f.turns = l.turns * lengthShare;
    f.c = discCentreDirection(f.normal, f.d, phase);
    return f;
}

// A fine fibre's disc (radius 0) at a root.
FurDisc fibreDisc(const FurRoot& r, uint64_t seed, const FurLook& look) {
    const FurFibres& fb = look.fibres;
    const int t = r.tri, i = r.index;
    const auto draw = [&](int channel) { return hash01(seed, t, i, channel); };
    FurDisc f;
    f.root = r.position;
    f.normal = r.normal;
    f.tri = t;
    f.index = i;
    f.budgetKey = r.budgetKey;
    f.layer = FurLayer::Fibres;
    const double v = std::clamp(fb.direction, 0.0, 1.0);
    const double combJitterDeg = 180 * v, liftJitterDeg = 40 * v;
    const double turn = (2 * draw(kCombChannel) - 1) * combJitterDeg * kDeg;
    const double lift = (fb.liftDeg + (2 * draw(kLiftChannel) - 1) * liftJitterDeg) * kDeg;
    f.d = furDirection(f.normal, kComb, 2 * kPi * draw(kTangentChannel), turn, lift);
    // Length: log-normal about the median (Box-Muller, clamped at 3 sigma), then [median / 2, 2 median].
    const double z = hashNormal(seed, t, i, kLengthChannel, kLengthChannel2, 3.0);
    f.straightLength = std::clamp(fb.lengthMedian * std::exp(fb.lengthSpread * z), 0.5 * fb.lengthMedian, 2 * fb.lengthMedian);
    f.rootRadius = 0.5 * fb.rootDiameter;
    f.tipRadius = f.rootRadius * fb.taper;
    return f;
}

// --- the irregularities

// What the shared fields say at one root.
struct FurLocal {
    double patch = 0, wear = 0;
    IslandWeight island;
    int clump = -1;
    double clumpWeight = 0;
};

int64_t cellKey(const Vec3& p, double cell) {
    const int64_t x = static_cast<int64_t>(std::floor(p.x / cell)), y = static_cast<int64_t>(std::floor(p.y / cell)),
                  z = static_cast<int64_t>(std::floor(p.z / cell));
    return (x & 0x1FFFFF) | (y & 0x1FFFFF) << 21 | (z & 0x1FFFFF) << 42;
}

constexpr int64_t kMaxClumpGridCells = int64_t(1) << 21;  // 8 MB of cell starts per part at most

// The clumps' grid (FurClumpGrid): counted per cell, then filled in clump order.
FurClumpGrid buildClumpGrid(const std::vector<FurClump>& clumps, double reach) {
    FurClumpGrid g;
    g.reach = reach;
    if (clumps.empty()) return g;
    Vec3 lo = clumps[0].centre, hi = lo;
    for (const FurClump& c : clumps) {
        lo = {std::min(lo.x, c.centre.x), std::min(lo.y, c.centre.y), std::min(lo.z, c.centre.z)};
        hi = {std::max(hi.x, c.centre.x), std::max(hi.y, c.centre.y), std::max(hi.z, c.centre.z)};
    }
    const Vec3 extent = hi - lo;
    g.cell = reach;
    // Enlarge the cell until the box (plus one cell on each side) fits the cap.
    const auto cellsAlong = [&](double length) { return static_cast<int64_t>(std::floor(length / g.cell)) + 3; };
    while (cellsAlong(extent.x) * cellsAlong(extent.y) * cellsAlong(extent.z) > kMaxClumpGridCells) g.cell *= 1.25;
    g.origin = lo - Vec3{g.cell, g.cell, g.cell};
    g.nx = static_cast<int>(cellsAlong(extent.x));
    g.ny = static_cast<int>(cellsAlong(extent.y));
    g.nz = static_cast<int>(cellsAlong(extent.z));
    std::vector<int32_t> cellOf(clumps.size());
    g.start.assign(static_cast<size_t>(g.nx) * g.ny * g.nz + 1, 0);
    for (size_t i = 0; i < clumps.size(); ++i) {
        const Vec3 q = (clumps[i].centre - g.origin) / g.cell;
        cellOf[i] = static_cast<int32_t>(int(q.x) + g.nx * (int(q.y) + g.ny * int(q.z)));
        ++g.start[static_cast<size_t>(cellOf[i]) + 1];
    }
    for (size_t i = 1; i < g.start.size(); ++i) g.start[i] += g.start[i - 1];
    g.clump.resize(clumps.size());
    std::vector<int32_t> fill(g.start.begin(), g.start.end() - 1);
    for (size_t i = 0; i < clumps.size(); ++i) g.clump[static_cast<size_t>(fill[static_cast<size_t>(cellOf[i])]++)] = static_cast<int32_t>(i);
    return g;
}

// The nearest clump within reach on the same side (n . n_k > 0.5), and its weight strength (1 - dist / reach)^2.
void nearestClump(const FurPartModifiers& m, const Vec3& p, const Vec3& n, double strength, FurLocal& local) {
    const FurClumpGrid& g = m.clumpGrid;
    if (g.clump.empty()) return;
    const Vec3 q = (p - g.origin) / g.cell;
    const int cx = static_cast<int>(std::floor(q.x)), cy = static_cast<int>(std::floor(q.y)), cz = static_cast<int>(std::floor(q.z));
    double best = g.reach;
    for (int z = std::max(0, cz - 1); z <= std::min(g.nz - 1, cz + 1); ++z)
        for (int y = std::max(0, cy - 1); y <= std::min(g.ny - 1, cy + 1); ++y)
            for (int x = std::max(0, cx - 1); x <= std::min(g.nx - 1, cx + 1); ++x) {
                const size_t cell = static_cast<size_t>(x) + static_cast<size_t>(g.nx) * (static_cast<size_t>(y) + static_cast<size_t>(g.ny) * z);
                for (int32_t k = g.start[cell]; k < g.start[cell + 1]; ++k) {
                    const int i = g.clump[static_cast<size_t>(k)];
                    const FurClump& c = m.clumps[static_cast<size_t>(i)];
                    const double d = norm(p - c.centre);
                    if (d < best && dot(n, c.normal) > 0.5) {
                        best = d;
                        local.clump = i;
                    }
                }
            }
    if (local.clump >= 0) local.clumpWeight = strength * (1 - best / g.reach) * (1 - best / g.reach);
}

// The vertex fields at a root (patch, wear, islands), interpolated by its barycentric weights; no clump yet.
FurLocal vertexFieldsAt(const FurPartModifiers& m, const FurRoot& r) {
    FurLocal l;
    const std::array<int, 3>& v = m.tri[static_cast<size_t>(r.tri)];
    const double b[3] = {r.bary.x, r.bary.y, r.bary.z};
    for (int k = 0; k < 3; ++k) {
        const size_t i = static_cast<size_t>(v[k]);
        l.patch += b[k] * m.patch[i];
        l.wear += b[k] * m.wear[i];
        l.island.fuzzy += static_cast<float>(b[k] * m.island[i].fuzzy);
        l.island.matted += static_cast<float>(b[k] * m.island[i].matted);
        l.island.worn += static_cast<float>(b[k] * m.island[i].worn);
    }
    return l;
}

// The local density over the layer's density (thin, fluffy for the fine fibres, patches, wear on exposed places).
double densityFactor(const FurLocal& l, const FurLook& look, FurLayer layer) {
    double f = (1 + std::clamp(look.variation.density, 0.0, 1.0) * l.patch) * (1 - l.island.worn);
    if (layer == FurLayer::Fibres) f *= 1 + l.island.fuzzy;
    if (l.wear > 0) f *= 1 - 0.3 * look.irregular.wear * l.wear;
    return std::max(0.0, f);
}

// The frame (d, c) rotated rigidly so that d becomes `to` (the minimal rotation; antiparallel: about the normal).
void turnFrame(FurDisc& f, const Vec3& to) {
    const Vec3 axis = cross(f.d, to);
    const double s = norm(axis), c = dot(f.d, to);
    if (s < 1e-12 && c > 0) return;
    const Vec3 k = s < 1e-12 ? f.normal : axis / s;
    const double angle = s < 1e-12 ? kPi : std::atan2(s, c);
    f.d = normalized(rotateAbout(f.d, k, angle));
    if (f.radius > 0) f.c = normalized(rotateAbout(f.c, k, angle));
}

// d's elevation above the tangent plane scaled by `factor` (clamped to [0, 90] deg); the frame turns with it.
void scaleLift(FurDisc& f, double factor) {
    const double up = std::clamp(dot(f.d, f.normal), -1.0, 1.0);
    const Vec3 tangent = f.d - f.normal * up;
    if (norm(tangent) < 1e-9 || factor == 1) return;
    const double lift = std::clamp(std::asin(up) * factor, 0.0, kPi / 2);
    turnFrame(f, normalized(normalized(tangent) * std::cos(lift) + f.normal * std::sin(lift)));
}

// The composition: clump -> island shape effects -> patches -> wear; the segment rules follow.
void applyModifiers(FurDisc& f, const FurLocal& l, const FurPartModifiers& m, const FurLook& look) {
    const bool loop = f.layer == FurLayer::Loops;
    const double matted = l.island.matted, fluffy = l.island.fuzzy;
    // Clump: lean toward the nub's direction by its weight (matted islands clump tighter), take its factors.
    if (l.clump >= 0) {
        const FurClump& c = m.clumps[static_cast<size_t>(l.clump)];
        const double w = std::min(1.0, l.clumpWeight + (loop ? 0.5 * matted : 0.0));
        const double angle = std::acos(std::clamp(dot(f.d, c.d), -1.0, 1.0));
        if (w > 0 && angle > 1e-9) {
            Vec3 axis = cross(f.d, c.d);
            axis = norm(axis) < 1e-12 ? f.normal : normalized(axis);
            const Vec3 to = normalized(rotateAbout(f.d, axis, w * angle));
            turnFrame(f, to);
        }
        if (loop && f.radius > 0) {
            f.radius *= 1 + w * (c.radiusFactor - 1);
            f.turns *= 1 + w * (c.turnsFactor - 1);
        }
    }
    double lift = 1, length = 1;
    if (loop) {
        lift *= 1 - 0.6 * matted;
        length *= (1 - 0.3 * matted) * (1 - 0.25 * fluffy);
        if (f.radius > 0) f.radius *= 1 + 0.3 * fluffy;
        lift *= 1 + look.irregular.patchLift * l.patch;
        length *= 1 + look.irregular.patchLength * l.patch;
        if (l.wear > 0) {
            length *= 1 - 0.4 * look.irregular.wear * l.wear;
            lift *= 1 - 0.5 * look.irregular.wear * l.wear;
        } else {
            length *= 1 + 0.2 * look.irregular.wear * -l.wear;
        }
    } else {
        lift *= 1 - 0.7 * matted;
        if (l.wear > 0) length *= 1 - 0.4 * look.irregular.wear * l.wear;
    }
    scaleLift(f, lift);
    if (f.radius > 0) {
        f.turns = std::max(0.05, f.turns * length);
        f.radius = std::max(0.05, f.radius);
    } else {
        f.straightLength = std::max(0.05, f.straightLength * length);
    }
}

// The discs of one cache's active candidates at `density`, in the cache's order, appended to `out`.
void addLayerDiscs(const FurRootCache& cache, double density, const FurLook& look, const FurCoverage& coverage, const FurPartModifiers* mods,
                   const ExecPolicy& exec, float keyBelow, std::vector<FurDisc>& out, FurLayerCounts& counts) {
    const bool modified = mods && mods->active;
    // The candidates a sampler at `density` draws: the first n_t(density) of each triangle (fur_field.h).
    std::vector<int> active(cache.triangleArea.size(), 0);
    for (size_t t = 0; t < active.size(); ++t)
        active[t] = static_cast<int>(std::floor(density * cache.triangleArea[t] + cache.triangleOffset[t]));
    // Two passes, so no disc is written twice: which candidates stay (active, inside the coverage), then each kept
    // one's disc straight into its slot (the order stays the cache's).
    const size_t n = cache.roots.size();
    std::vector<char> state(n, 0);  // 0 inactive, 1 kept, 2 masked, 3 thinned by the local density
    // Each candidate's nearest clump, found once in the first pass (the costliest lookup) and read in the second.
    std::vector<int32_t> clumpOf(modified ? n : 0, -1);
    std::vector<double> clumpWeightOf(modified ? n : 0, 0.0);
    const double layerDensity = cache.layer == FurLayer::Loops ? look.loops.density : look.fibres.density;  // `density` is its bound
    parallelFor(static_cast<int64_t>(n), 4096, [&](int64_t b, int64_t e, int64_t) {
        cancelCheckpoint();
        for (int64_t q = b; q < e; ++q) {
            const FurRoot& r = cache.roots[static_cast<size_t>(q)];
            if (r.index >= active[static_cast<size_t>(r.tri)] || !(r.budgetKey < keyBelow)) continue;
            const bool masked = coverage && hash01(cache.seed, r.tri, r.index, kCoverageChannel) >= std::clamp(coverage(r.tri, r.bary), 0.0, 1.0);
            if (masked) {
                state[static_cast<size_t>(q)] = 2;
                continue;
            }
            // The local density: kept with probability D_local / D_bound (its own channel; never the cache's capacity).
            if (modified) {
                FurLocal local = vertexFieldsAt(*mods, r);
                nearestClump(*mods, r.position, r.normal, look.irregular.clumpStrength, local);
                clumpOf[static_cast<size_t>(q)] = local.clump;
                clumpWeightOf[static_cast<size_t>(q)] = local.clumpWeight;
                if (hash01(cache.seed, r.tri, r.index, kDensityChannel) >= layerDensity * densityFactor(local, look, cache.layer) / density) {
                    state[static_cast<size_t>(q)] = 3;
                    continue;
                }
            }
            state[static_cast<size_t>(q)] = 1;
        }
    }, exec);
    std::vector<int32_t> slot(n, -1);
    int64_t kept = 0, masked = 0, thinned = 0, sampled = 0;
    for (size_t q = 0; q < n; ++q) {
        masked += state[q] == 2;
        thinned += state[q] == 3;
        if (state[q] == 1) slot[q] = static_cast<int32_t>(kept++);
    }
    const size_t base = out.size();
    out.resize(base + static_cast<size_t>(kept));
    parallelFor(static_cast<int64_t>(n), 4096, [&](int64_t b, int64_t e, int64_t) {
        cancelCheckpoint();
        for (int64_t q = b; q < e; ++q) {
            if (slot[static_cast<size_t>(q)] < 0) continue;
            const FurRoot& r = cache.roots[static_cast<size_t>(q)];
            FurDisc& d = out[base + static_cast<size_t>(slot[static_cast<size_t>(q)])];
            d = cache.layer == FurLayer::Loops ? loopDisc(r, cache.seed, look) : fibreDisc(r, cache.seed, look);
            if (modified) {
                FurLocal local = vertexFieldsAt(*mods, r);
                local.clump = clumpOf[static_cast<size_t>(q)];
                local.clumpWeight = clumpWeightOf[static_cast<size_t>(q)];
                applyModifiers(d, local, *mods, look);
            }
        }
    }, exec);
    // Counted at the density asked for, not the capacity: the hidden ones are the active candidates the cache dropped.
    for (int a : active) sampled += std::max(0, a);
    if (keyBelow < 1) sampled = static_cast<int64_t>(std::llround(double(sampled) * keyBelow));  // the share asked for
    counts.sampled += sampled;
    counts.masked += masked;
    counts.thinned += thinned;
    counts.kept += kept;
    counts.hidden += sampled - kept - masked - thinned;
}

} // namespace

FurRootCache sampleFurRoots(const TriMesh& m, FurLayer layer, double capacity, uint64_t seed, const FurHiddenTest& hidden,
                            const ExecPolicy& exec) {
    FurRootCache c;
    c.layer = layer;
    c.seed = seed;
    c.capacity = capacity;
    const size_t nt = m.tri.size();
    c.triangleArea.assign(nt, 0);
    c.triangleOffset.assign(nt, 0);
    // The candidates: n_t = floor(D A_t + u_t) per triangle, in triangle order.
    std::vector<FurRoot> cand;
    for (size_t t = 0; t < nt; ++t) {
        const Vec3 &a = m.pos[m.tri[t][0]], &b = m.pos[m.tri[t][1]], &cc = m.pos[m.tri[t][2]];
        const double area = 0.5 * norm(cross(b - a, cc - a));
        if (!(area > 0)) continue;
        c.triangleArea[t] = area;
        c.triangleOffset[t] = hash01(seed, int(t), 0, kCountChannel);
        const int n = static_cast<int>(std::floor(capacity * area + c.triangleOffset[t]));
        for (int i = 0; i < n; ++i) {
            FurRoot r;
            r.tri = int(t);
            r.index = i;
            cand.push_back(r);
        }
    }
    // Placement, normal and the hidden test of every candidate (independent: parallel).
    std::vector<char> visible(cand.size(), 0);
    parallelFor(static_cast<int64_t>(cand.size()), 2048, [&](int64_t b, int64_t e, int64_t) {
        cancelCheckpoint();
        for (int64_t q = b; q < e; ++q) {
            FurRoot& r = cand[static_cast<size_t>(q)];
            const int t = r.tri, i = r.index;
            const Vec3 &A = m.pos[m.tri[t][0]], &B = m.pos[m.tri[t][1]], &C = m.pos[m.tri[t][2]];
            const double s = std::sqrt(hash01(seed, t, i, kPlacementChannel1)), u = hash01(seed, t, i, kPlacementChannel2);
            r.bary = {1 - s, s * (1 - u), s * u};
            r.position = A * r.bary.x + B * r.bary.y + C * r.bary.z;
            const Vec3 faceCross = cross(B - A, C - A);
            r.normal = faceCross / norm(faceCross);
            r.budgetKey = float(hash01(seed, t, i, kKeyChannel));
            visible[static_cast<size_t>(q)] = !(hidden && hidden(r.position + r.normal * kHiddenTestLift));
        }
    }, exec);
    c.sampled = static_cast<int64_t>(cand.size());
    c.roots.reserve(cand.size());
    for (size_t q = 0; q < cand.size(); ++q)
        if (visible[q]) c.roots.push_back(cand[q]);
    c.hidden = c.sampled - static_cast<int64_t>(c.roots.size());
    return c;
}

uint64_t furPartSeed(uint64_t lookSeed, FurLayer layer, uint64_t partNameHash) {
    const uint64_t layerSeed = layer == FurLayer::Fibres ? splitmix64(lookSeed ^ kFibresSalt) : lookSeed;
    return partNameHash ? splitmix64(layerSeed ^ splitmix64(partNameHash)) : layerSeed;
}

double furDensityBound(const FurLook& look, FurLayer layer) {
    const double density = layer == FurLayer::Loops ? look.loops.density : look.fibres.density;
    const FurIslands& is = look.islands;
    const bool fluffyIslands = layer == FurLayer::Fibres && is.rate > 0 && is.fluffy > 0 && is.strength > 0;
    return density * (1 + std::clamp(look.variation.density, 0.0, 1.0)) * (fluffyIslands ? 1 + std::clamp(is.strength, 0.0, 1.0) : 1.0);
}

FurPartModifiers buildFurModifiers(const TriMesh& m, const FurLook& look, uint64_t partNameHash) {
    FurPartModifiers out;
    const FurIrregularity& ir = look.irregular;
    const bool islands = look.islands.rate > 0 && look.islands.strength > 0 && look.islands.matted + look.islands.fluffy + look.islands.thin > 0;
    out.active = look.variation.density > 0 || ir.patchLength > 0 || ir.patchLift > 0 || ir.wear > 0 || ir.clumpStrength > 0 || islands;
    if (!out.active) return out;
    const size_t nv = m.pos.size();
    out.tri = m.tri;
    out.patch.assign(nv, 0);
    out.wear.assign(nv, 0);
    out.island.assign(nv, IslandWeight{});
    const uint64_t seed = splitmix64(look.seed ^ splitmix64(partNameHash));
    // Patches: fractal noise at the feature size.
    if (look.variation.density > 0 || ir.patchLength > 0 || ir.patchLift > 0) {
        const double scale = std::max(0.5, look.variation.densityScale);
        for (size_t i = 0; i < nv; ++i) out.patch[i] = fractalNoise(m.pos[i] * (1 / scale), seed ^ kPatchSalt, kPatchOctaves);
    }
    // Islands: blobs on the surface, weights per vertex.
    if (islands) {
        SurfaceIslands si;
        si.ratePerCm2 = look.islands.rate;
        si.size = look.islands.size;
        si.strength = std::clamp(look.islands.strength, 0.0, 1.0);
        si.fluffy = look.islands.fluffy;
        si.matted = look.islands.matted;
        si.thin = look.islands.thin;
        const std::vector<IslandBlob> blobs = surfaceIslandBlobs(m, si, seed ^ kIslandSalt);
        for (size_t i = 0; i < nv; ++i) out.island[i] = islandWeightAt(blobs, m.pos[i]);
    }
    // Wear: the mean curvature per vertex from the umbrella vector, kappa = 2 n . (p - mean(neighbours)) / edge^2 (1/R on
    // a sphere; convex > 0), as e = kappa x kWearRadius clamped to [-1, 1], smoothed twice over the neighbours.
    if (ir.wear > 0) {
        std::vector<std::vector<int>> nb(nv);
        for (const std::array<int, 3>& t : m.tri)
            for (int k = 0; k < 3; ++k) {
                nb[static_cast<size_t>(t[k])].push_back(t[(k + 1) % 3]);
                nb[static_cast<size_t>(t[k])].push_back(t[(k + 2) % 3]);
            }
        std::vector<Vec3> normal(nv);
        for (const std::array<int, 3>& t : m.tri) {
            const Vec3 fn = cross(m.pos[t[1]] - m.pos[t[0]], m.pos[t[2]] - m.pos[t[0]]);
            for (int k = 0; k < 3; ++k) normal[static_cast<size_t>(t[k])] += fn;
        }
        for (size_t i = 0; i < nv; ++i) {
            if (nb[i].empty() || norm(normal[i]) < 1e-12) continue;
            Vec3 mean;
            double edge = 0;
            for (int j : nb[i]) {
                mean += m.pos[static_cast<size_t>(j)] / double(nb[i].size());
                edge += norm(m.pos[static_cast<size_t>(j)] - m.pos[i]) / double(nb[i].size());
            }
            const double kappa = edge > 0 ? 2 * dot(normalized(normal[i]), m.pos[i] - mean) / (edge * edge) : 0.0;
            out.wear[i] = std::clamp(kappa * kWearRadius, -1.0, 1.0);
        }
        for (int round = 0; round < 2; ++round) {
            std::vector<double> next = out.wear;
            for (size_t i = 0; i < nv; ++i) {
                if (nb[i].empty()) continue;
                double sum = out.wear[i];
                for (int j : nb[i]) sum += out.wear[static_cast<size_t>(j)];
                next[i] = sum / double(nb[i].size() + 1);
            }
            out.wear = std::move(next);
        }
    }
    // Clumps: Poisson-disc centres (area-weighted candidates, kept in order when 0.8 size from every kept one).
    if (ir.clumpStrength > 0) {
        const double size = std::max(0.5, ir.clumpSize), minGap = 0.8 * size;
        const uint64_t cs = seed ^ kClumpSalt;
        std::unordered_map<int64_t, std::vector<int>> kept;  // cell of minGap -> its clumps
        for (size_t t = 0; t < m.tri.size(); ++t) {
            const Vec3 &A = m.pos[m.tri[t][0]], &B = m.pos[m.tri[t][1]], &C = m.pos[m.tri[t][2]];
            const Vec3 crossed = cross(B - A, C - A);
            const double area = 0.5 * norm(crossed);
            if (!(area > 0)) continue;
            const int n = static_cast<int>(std::floor(2 / (size * size) * area + hash01(cs, int(t), 0, 0)));
            for (int i = 0; i < n; ++i) {
                const auto draw = [&](int channel) { return hash01(cs, int(t), i + 1, channel); };
                const double sq = std::sqrt(draw(1)), u = draw(2);
                FurClump c;
                c.centre = A * (1 - sq) + B * (sq * (1 - u)) + C * (sq * u);
                c.normal = crossed / (2 * area);
                bool free = true;
                for (int dx = -1; dx <= 1 && free; ++dx)
                    for (int dy = -1; dy <= 1 && free; ++dy)
                        for (int dz = -1; dz <= 1 && free; ++dz) {
                            const int64_t key = cellKey(c.centre + Vec3{dx * minGap, dy * minGap, dz * minGap}, minGap);
                            const auto cellIt = kept.find(key);
                            if (cellIt == kept.end()) continue;
                            for (int j : cellIt->second) free &= norm(out.clumps[static_cast<size_t>(j)].centre - c.centre) >= minGap;
                        }
                if (!free) continue;
                c.d = furDirection(c.normal, kComb, 2 * kPi * draw(7), 2 * kPi * draw(3), (20 + 40 * draw(4)) * kDeg);
                c.radiusFactor = 1 + 0.3 * (2 * draw(5) - 1);
                c.turnsFactor = 1 + 0.2 * (2 * draw(6) - 1);
                kept[cellKey(c.centre, minGap)].push_back(static_cast<int>(out.clumps.size()));
                out.clumps.push_back(c);
            }
        }
        out.clumpGrid = buildClumpGrid(out.clumps, size);
    }
    return out;
}

FurField buildFurField(const FurFieldRequest& r, const ExecPolicy& exec) {
    FurField f;
    const FurLook& look = r.look;
    addLayerDiscs(r.loops, furDensityBound(look, FurLayer::Loops), look, r.coverage, r.modifiers, exec, r.keyBelow, f.discs, f.loops);
    if (look.fibres.enabled && look.fibres.density > 0 && !r.fibres.roots.empty())
        addLayerDiscs(r.fibres, furDensityBound(look, FurLayer::Fibres), look, r.coverage, r.modifiers, exec, r.keyBelow, f.discs, f.fibres);
    return f;
}

} // namespace knit
