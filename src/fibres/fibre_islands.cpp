#include "fibres/fibre_islands.h"

#include <algorithm>
#include <cmath>

#include "model/choices.h"
#include "model/hash.h"

namespace knit {

namespace {

// Hash channels of an island's draws. A channel never changes meaning: that would move every island.
enum IslandChannel { kPlaceChannel = 0, kKindChannel = 1, kSizeChannel = 2, kSatelliteChannel = 3, kSurfacePlaceChannel = 9 };

constexpr int kSatellites = 3;
constexpr double kSatelliteReach = 0.8;  // of the island's radius: how far out a satellite's centre sits

// The kind is drawn by the shares; strength scales it.
IslandWeight peakOf(const FibreIslands& islands, double draw) {
    const float strength = static_cast<float>(islands.strength / 100);
    const double x = draw * (islands.fuzzy + islands.matted + islands.worn);
    if (x < islands.fuzzy) return {strength, 0, 0};
    if (x < islands.fuzzy + islands.matted) return {0, strength, 0};
    return {0, 0, strength};
}

std::vector<IslandBlob> placeIslands(const IslandRequest& request) {
    const FibreIslands& islands = request.islands;
    const double chance = islands.rate * request.stitchArea / 100;  // 100 mm2 to the cm2
    std::vector<IslandBlob> blobs;
    for (const StitchRecord& st : request.stitches.stitches) {
        const auto draw = [&](int channel) { return hash01(request.seed, st.turn, st.index, channel); };
        if (draw(kPlaceChannel) >= chance) continue;
        const double radius = 0.5 * islands.size * (0.7 + 0.6 * draw(kSizeChannel));
        addIsland(blobs, st.origin, st.U, st.V, radius, peakOf(islands, draw(kKindChannel)), request.seed, st.turn, st.index);
    }
    return blobs;
}

} // namespace

void addIsland(std::vector<IslandBlob>& blobs, const Vec3& centre, const Vec3& U, const Vec3& V, double radius, IslandWeight peak, uint64_t seed,
               int a, int b) {
    blobs.push_back({centre, radius, peak});
    for (int k = 0; k < kSatellites; ++k) {
        const double angle = 2 * kPi * hash01(seed, a, b, kSatelliteChannel + 2 * k), share = hash01(seed, a, b, kSatelliteChannel + 2 * k + 1);
        const Vec3 reach = (U * std::cos(angle) + V * std::sin(angle)) * (kSatelliteReach * radius);
        blobs.push_back({centre + reach, radius * (0.5 + 0.3 * share), peak});
    }
}

IslandWeight islandWeightAt(const std::vector<IslandBlob>& blobs, const Vec3& p) {
    IslandWeight w;
    for (const IslandBlob& b : blobs) {
        const double d = norm(p - b.centre);
        if (d >= b.radius) continue;
        const double x = 1 - d / b.radius;
        const float fade = static_cast<float>(x * x * (3 - 2 * x));
        w.fuzzy = std::max(w.fuzzy, b.peak.fuzzy * fade);
        w.matted = std::max(w.matted, b.peak.matted * fade);
        w.worn = std::max(w.worn, b.peak.worn * fade);
    }
    return w;
}

std::vector<IslandBlob> surfaceIslandBlobs(const TriMesh& m, const SurfaceIslands& islands, uint64_t seed) {
    std::vector<IslandBlob> blobs;
    const double shares = islands.fluffy + islands.matted + islands.thin;
    if (!(islands.ratePerCm2 > 0) || !(shares > 0) || !(islands.strength > 0)) return blobs;
    const float strength = static_cast<float>(islands.strength);
    for (size_t t = 0; t < m.tri.size(); ++t) {
        const Vec3 &A = m.pos[m.tri[t][0]], &B = m.pos[m.tri[t][1]], &C = m.pos[m.tri[t][2]];
        const Vec3 crossed = cross(B - A, C - A);
        const double area = 0.5 * norm(crossed);
        if (!(area > 0)) continue;
        const int n = static_cast<int>(std::floor(islands.ratePerCm2 * area / 100 + hash01(seed, int(t), 0, kPlaceChannel)));
        for (int i = 0; i < n; ++i) {
            const auto draw = [&](int channel) { return hash01(seed, int(t), i + 1, channel); };
            const double s = std::sqrt(draw(kPlaceChannel)), u = draw(kSurfacePlaceChannel);
            const Vec3 centre = A * (1 - s) + B * (s * (1 - u)) + C * (s * u);
            const Vec3 normal = crossed / (2 * area), U = normalized(B - A), V = cross(normal, U);
            const double x = draw(kKindChannel) * shares;
            const IslandWeight peak = x < islands.fluffy ? IslandWeight{strength, 0, 0} : x < islands.fluffy + islands.matted ? IslandWeight{0, strength, 0}
                                                                                                                                : IslandWeight{0, 0, strength};
            addIsland(blobs, centre, U, V, 0.5 * islands.size * (0.7 + 0.6 * draw(kSizeChannel)), peak, seed, int(t), i + 1);
        }
    }
    return blobs;
}

std::vector<std::vector<IslandWeight>> islandWeights(const IslandRequest& request) {
    const std::vector<IslandBlob> blobs = placeIslands(request);
    std::vector<std::vector<IslandWeight>> out;
    out.reserve(request.strands.yarns.size());
    for (const YarnPath& y : request.strands.yarns) {
        std::vector<IslandWeight> weights;
        weights.reserve(y.curve.ctrl.size());
        for (const Vec3& p : y.curve.ctrl) weights.push_back(islandWeightAt(blobs, p));
        out.push_back(std::move(weights));
    }
    return out;
}

IslandWeight islandAt(const std::vector<IslandWeight>& strand, double u) {
    const double c = std::clamp(u + 1, 0.0, static_cast<double>(strand.size() - 1));  // control c sits at the knot u = c - 1
    const size_t k = static_cast<size_t>(c), next = std::min(k + 1, strand.size() - 1);
    const float f = static_cast<float>(c - static_cast<double>(k));
    const IslandWeight& a = strand[k];
    const IslandWeight& b = strand[next];
    return {a.fuzzy * (1 - f) + b.fuzzy * f, a.matted * (1 - f) + b.matted * f, a.worn * (1 - f) + b.worn * f};
}

} // namespace knit
