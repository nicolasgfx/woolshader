// GPU backend: queries over a whole scene's parts (gpu/surfaces.h checkSceneClearance).
//
// Cross-part clearance: every pair of sampled segments of different groups closer than
// R_a + R_b + c. Positions and distances are double, as on the CPU: the same pairs, with distances
// within a few ulps (the device contracts to fused multiply-add); float appears only in the grid's
// boxes, rounded outward. Variable-length list in three steps: count
// per segment, offsets on the host, fill. Batches paced below the display watchdog.
#include <algorithm>
#include <cmath>
#include <map>
#include <vector>

#include "gpu/geometry/device_geometry.h"
#include "gpu/geometry/math.cuh"
#include "gpu/device.h"
#include "gpu/surfaces.h"
#include "model/cancel.h"

namespace knit::gpu {

namespace {

// One close pair: the two parts, part A the smaller, and their segments' distance.
struct ClosePairRecord {
    int partA, partB;
    double distance;
};

struct ClearanceView {
    const ClearanceSegment* segments;
    GridView grid;
    double clearance, reach;  // reach: the largest threshold of any pair
    float queryPad;           // reach, rounded up
};

// What a pass writes for its batch. pairs null: only the count and the nearest distance.
struct ClearanceOutput {
    int* counts;
    double* nearest;
    const long long* offsets;
    ClosePairRecord* pairs;
};

__host__ __device__ __forceinline__ d3 toD3(const Vec3& v) { return {v.x, v.y, v.z}; }

__global__ void closePairsKernel(ClearanceView v, int first, int count, ClearanceOutput out) {
    const int item = threadItem();
    if (item >= count) return;
    const int i = first + item;
    const ClearanceSegment& a = v.segments[i];
    const f3 pad{v.queryPad, v.queryPad, v.queryPad};
    const f3 lo = fromF(v.grid.boxLo + 3 * i) - pad, hi = fromF(v.grid.boxHi + 3 * i) + pad;
    int close = 0;
    double nearest = 1e300;
    for (GridCursor cursor(v.grid, lo, hi); cursor.next();) {
        const int j = cursor.item();
        if (j <= i) continue;
        const ClearanceSegment& b = v.segments[j];
        if (b.group == a.group) continue;
        double s, t;
        const double d = segmentDistance(toD3(a.a), toD3(a.b), toD3(b.a), toD3(b.b), s, t);
        if (d < v.reach) nearest = d < nearest ? d : nearest;
        if (d >= a.radius + b.radius + v.clearance) continue;
        if (out.pairs) out.pairs[out.offsets[i] + close] = {a.part < b.part ? a.part : b.part, a.part < b.part ? b.part : a.part, d};
        ++close;
    }
    out.counts[i] = close;
    out.nearest[i] = nearest;
}

// Every segment once, in batches.
void runClosePairs(const ClearanceView& view, int segments, const ClearanceOutput& out) {
    LaunchPacer pace(1ll << 14, 256, 1ll << 22);
    for (int first = 0, n = 0; first < segments; first += n) {
        cancelCheckpoint();
        n = static_cast<int>(std::min<long long>(pace.next(), segments - first));
        const auto started = Clock::now();
        closePairsKernel<<<blocksFor(n), kThreadsPerBlock>>>(view, first, n, out);
        checkLaunch("scene clearance");
        waitForKernels("scene clearance");
        pace.measured(n, elapsedMs(started));
    }
}

// The segments' boxes, rounded outward to float.
BoxArrays segmentBoxes(std::span<const ClearanceSegment> segments, float cell) {
    BoxArrays boxes;
    boxes.cell = cell;
    boxes.lo.reserve(3 * segments.size());
    boxes.hi.reserve(3 * segments.size());
    const auto below = [](double v) { const float f = static_cast<float>(v); return double(f) > v ? std::nextafter(f, -3.0e38f) : f; };
    const auto above = [](double v) { const float f = static_cast<float>(v); return double(f) < v ? std::nextafter(f, 3.0e38f) : f; };
    for (const ClearanceSegment& g : segments) {
        boxes.lo.insert(boxes.lo.end(), {below(std::min(g.a.x, g.b.x)), below(std::min(g.a.y, g.b.y)), below(std::min(g.a.z, g.b.z))});
        boxes.hi.insert(boxes.hi.end(), {above(std::max(g.a.x, g.b.x)), above(std::max(g.a.y, g.b.y)), above(std::max(g.a.z, g.b.z))});
    }
    return boxes;
}

} // namespace

bool checkSceneClearance(const SceneClearanceRequest& request, SceneClearance& out, std::string* err) try {
    out = SceneClearance{};
    const std::span<const ClearanceSegment> segments = request.segments;
    if (segments.empty()) return true;
    if (!available(err)) return false;
    Reservation lease(Reservation::exclusiveBytes());  // runs alone on the device (runtime.h)
    std::map<int, double> radiusOf;
    double largest = 0;
    for (const ClearanceSegment& g : segments) {
        radiusOf[g.part] = g.radius;
        largest = std::max(largest, g.radius);
    }
    const double reach = 2 * largest + request.clearance;
    const int count = static_cast<int>(segments.size());

    DeviceBuffer records, counts, nearest, offsets, pairs;
    DeviceGrid grid;
    records.upload(segments.data(), segments.size_bytes());
    grid.upload(segmentBoxes(segments, static_cast<float>(std::max(reach, 0.5))));
    counts.alloc(4ull * count);
    nearest.alloc(8ull * count);
    const ClearanceView view{static_cast<const ClearanceSegment*>(records.data()), grid.view, request.clearance, reach,
                             std::nextafter(static_cast<float>(reach), 3.0e38f)};
    runClosePairs(view, count, {counts.ints(), nearest.doubles(), nullptr, nullptr});

    std::vector<int> closePerSegment;
    std::vector<double> nearestPerSegment;
    counts.download(closePerSegment, count);
    nearest.download(nearestPerSegment, count);
    for (double d : nearestPerSegment) out.minDistance = std::min(out.minDistance, d);
    std::vector<long long> starts(count);
    long long total = 0;
    for (int i = 0; i < count; ++i) {
        starts[i] = total;
        total += closePerSegment[i];
    }
    if (total == 0) return true;
    offsets.upload(starts);
    pairs.alloc(static_cast<size_t>(total) * sizeof(ClosePairRecord));
    runClosePairs(view, count, {counts.ints(), nearest.doubles(), offsets.longs(), static_cast<ClosePairRecord*>(pairs.data())});

    std::vector<ClosePairRecord> found(static_cast<size_t>(total));
    pairs.download(found.data(), found.size() * sizeof(ClosePairRecord));
    std::map<std::pair<int, int>, ClosePartPair> byParts;
    for (const ClosePairRecord& r : found) {
        ClosePartPair& pair = byParts[{r.partA, r.partB}];
        if (pair.segments == 0) {
            pair = {r.partA, r.partB, 0, 1e300, radiusOf[r.partA] + radiusOf[r.partB] + request.clearance};
        }
        ++pair.segments;
        pair.minDistance = std::min(pair.minDistance, r.distance);
    }
    for (const auto& [parts, pair] : byParts) out.pairs.push_back(pair);
    return true;
} catch (const DeviceError& e) {
    return reportFailure(e, err);
}

} // namespace knit::gpu
