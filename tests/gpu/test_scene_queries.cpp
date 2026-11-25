// GPU port tests: scene queries (gpu/surface/scene_queries.cu) against the CPU loop they replaced.
#include <algorithm>
#include <cstdint>
#include <map>
#include <random>
#include <vector>

#include "gpu/device.h"
#include "gpu/surfaces.h"
#include "gpu_common.h"
#include "model/cancel.h"
#include "model/mesh.h"
#include "model/spatial_hash.h"

using namespace knit;

namespace {

// The written reference: the CPU loop of the cross-part clearance, as it ran before the port. Every
// pair of different groups the hash returns is measured; the nearest counts only within the reach.
gpu::SceneClearance cpuSceneClearance(const std::vector<gpu::ClearanceSegment>& segments, double clearance) {
    gpu::SceneClearance out;
    double largest = 0;
    for (const gpu::ClearanceSegment& g : segments) largest = std::max(largest, g.radius);
    const double reach = 2 * largest + clearance;
    SpatialHash3 hash(std::max(reach, 0.5));
    for (const gpu::ClearanceSegment& g : segments) {
        hash.insert({std::min(g.a.x, g.b.x), std::min(g.a.y, g.b.y), std::min(g.a.z, g.b.z)},
                    {std::max(g.a.x, g.b.x), std::max(g.a.y, g.b.y), std::max(g.a.z, g.b.z)});
    }
    std::map<std::pair<int, int>, gpu::ClosePartPair> close;
    const Vec3 pad{reach, reach, reach};
    for (size_t i = 0; i < segments.size(); ++i) {
        const gpu::ClearanceSegment& a = segments[i];
        hash.query(hash.lo(static_cast<int>(i)) - pad, hash.hi(static_cast<int>(i)) + pad, [&](int j) {
            if (j <= static_cast<int>(i)) return;
            const gpu::ClearanceSegment& b = segments[j];
            if (b.group == a.group) return;
            const double d = closestSegmentPoints(a.a, a.b, b.a, b.b).distance;
            if (d < reach) out.minDistance = std::min(out.minDistance, d);
            const double threshold = a.radius + b.radius + clearance;
            if (d >= threshold) return;
            const std::pair<int, int> key{std::min(a.part, b.part), std::max(a.part, b.part)};
            gpu::ClosePartPair& pair = close.try_emplace(key, gpu::ClosePartPair{key.first, key.second, 0, 1e300, threshold}).first->second;
            ++pair.segments;
            pair.minDistance = std::min(pair.minDistance, d);
        });
    }
    for (const auto& [key, pair] : close) out.pairs.push_back(pair);
    return out;
}

// Short polylines of several parts in one 40 mm box, in groups: parts 0 and 1 are one group.
std::vector<gpu::ClearanceSegment> randomSegments(int count) {
    std::mt19937 random(5);
    std::uniform_real_distribution<double> position(0, 40), step(-1, 1);
    std::vector<gpu::ClearanceSegment> segments;
    for (int i = 0; i < count; ++i) {
        const int part = i % 7;
        const Vec3 a{position(random), position(random), position(random)};
        segments.push_back({a, a + Vec3{step(random), step(random), step(random)}, 0.5 + 0.05 * part, part, part <= 1 ? 0 : part});
    }
    return segments;
}

} // namespace

TEST(gpu_scene_clearance_gpu_vs_cpu) {
    if (!gpu::available(nullptr)) SKIP("no usable CUDA GPU");
    const std::vector<gpu::ClearanceSegment> segments = randomSegments(6000);
    const double clearance = 0.1;
    const gpu::SceneClearance cpu = cpuSceneClearance(segments, clearance);
    gpu::SceneClearance device;
    std::string err;
    CHECK(gpu::checkSceneClearance({segments, clearance}, device, &err));
    std::printf("    %zu part pairs, nearest %.6f mm\n", cpu.pairs.size(), cpu.minDistance);
    CHECK(cpu.pairs.size() >= 10);
    CHECK(device.pairs.size() == cpu.pairs.size());
    CHECK_NEAR(device.minDistance, cpu.minDistance, 1e-12);
    for (size_t k = 0; k < cpu.pairs.size() && k < device.pairs.size(); ++k) {
        const gpu::ClosePartPair &c = cpu.pairs[k], &d = device.pairs[k];
        CHECK(c.partA == d.partA && c.partB == d.partB && c.segments == d.segments);
        CHECK_NEAR(c.minDistance, d.minDistance, 1e-12);  // the device contracts to fused multiply-add: a few ulps
        CHECK(c.threshold == d.threshold);
        CHECK(!(c.partA == 0 && c.partB == 1));  // one group: never compared
    }
}

TEST(gpu_scene_clearance_empty_and_far) {
    gpu::SceneClearance out;
    std::string err;
    CHECK(gpu::checkSceneClearance({{}, 0.1}, out, &err) && out.pairs.empty());
    if (!gpu::available(nullptr)) SKIP("no usable CUDA GPU");
    const std::vector<gpu::ClearanceSegment> far = {{{0, 0, 0}, {1, 0, 0}, 0.5, 0, 0}, {{0, 100, 0}, {1, 100, 0}, 0.5, 1, 1}};
    CHECK(gpu::checkSceneClearance({far, 0.1}, out, &err) && out.pairs.empty() && out.minDistance == 1e300);
}

// A cancelled caller stops before the first batch ends: the device call reports Cancelled, as every
// long loop does.
TEST(gpu_scene_clearance_cancels) {
    if (!gpu::available(nullptr)) SKIP("no usable CUDA GPU");
    const std::vector<gpu::ClearanceSegment> segments = randomSegments(2000);
    CancelToken token;
    token.cancel();
    ThreadCancelScope scope(&token);
    gpu::SceneClearance out;
    std::string err;
    bool stopped = false;
    try {
        gpu::checkSceneClearance({segments, 0.1}, out, &err);
    } catch (const Cancelled&) {
        stopped = true;
    }
    CHECK(stopped);
}
