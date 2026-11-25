// Guards of the GPU simplification: what later refactors must not change without a test saying so.
//   - The solver's identity text is hashed into every knit cache key; a different byte re-knits
//     every scene.
//   - The uniform grid (fibres, surfaces, the solver's broad phase, the mesh queries) is built on
//     the device by a sort; a bucket's order decides ties between queries, so it is pinned against
//     plain loops.
//   - A mesh shared by the pieces of a scene is on the device once.
#include <algorithm>
#include <climits>
#include <cmath>
#include <random>
#include <string>
#include <vector>

#include "gpu/geometry/device_geometry.h"
#include "gpu_common.h"
#include "model/mesh_query.h"
#include "gpu/device.h"
#include "gpu/knitting.h"

using namespace knit;

TEST(gpu_solver_identity_text) {
    CHECK(gpu::ncpSolverIdentity() ==
          "ncp-gpu v5 (GEOMETRY 8.21; device-controlled iteration, matrix-free PCG, compact contact rows; early stops: progress 0.01, infeasible window 3 at rho 1000, "
          "step budget 900; settle: one round, progress 0.05, contact tolViol 0.005; repair: one certified move)");
}

namespace {

gpu::BoxArrays randomBoxes(int n, float extent, float maxSize, float cell, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> at(-extent, extent), size(0.0f, maxSize);
    gpu::BoxArrays b;
    b.cell = cell;
    for (int i = 0; i < n; ++i)
        for (int a = 0; a < 3; ++a) {
            const float lo = at(rng);
            b.lo.push_back(lo);
            b.hi.push_back(lo + size(rng));
        }
    return b;
}

struct BuiltGrid {
    float cell = 1;
    long long lo[3] = {0, 0, 0}, hi[3] = {-1, -1, -1};
    std::vector<unsigned> cellStart, entries;
};

// The grid written as plain loops: each cell's items in increasing id order. Items with an inverted
// box stay out.
BuiltGrid referenceGrid(const gpu::BoxArrays& b, float cell) {
    BuiltGrid g;
    g.cell = cell;
    const float inv = 1.0f / cell;
    const auto cellOf = [inv](float v) { return static_cast<long long>(std::floor(v * inv)); };
    const auto included = [&b](int i) { return b.lo[3 * i] <= b.hi[3 * i]; };
    long long lowest[3] = {LLONG_MAX, LLONG_MAX, LLONG_MAX}, highest[3] = {LLONG_MIN, LLONG_MIN, LLONG_MIN};
    for (int i = 0; i < b.count(); ++i)
        for (int a = 0; a < 3 && included(i); ++a) {
            lowest[a] = std::min(lowest[a], cellOf(b.lo[3 * i + a]));
            highest[a] = std::max(highest[a], cellOf(b.hi[3 * i + a]));
        }
    if (lowest[0] > highest[0]) return g;
    long long n[3];
    for (int a = 0; a < 3; ++a) {
        g.lo[a] = lowest[a];
        g.hi[a] = highest[a];
        n[a] = highest[a] - lowest[a] + 1;
    }
    std::vector<std::vector<unsigned>> buckets(static_cast<size_t>(n[0] * n[1] * n[2]));
    for (int i = 0; i < b.count(); ++i) {
        if (!included(i)) continue;
        const long long first[3] = {cellOf(b.lo[3 * i]), cellOf(b.lo[3 * i + 1]), cellOf(b.lo[3 * i + 2])};
        const long long last[3] = {cellOf(b.hi[3 * i]), cellOf(b.hi[3 * i + 1]), cellOf(b.hi[3 * i + 2])};
        for (long long x = first[0]; x <= last[0]; ++x)
            for (long long y = first[1]; y <= last[1]; ++y)
                for (long long z = first[2]; z <= last[2]; ++z) {
                    const unsigned flags = (x == first[0] ? 1u : 0u) | (y == first[1] ? 2u : 0u) | (z == first[2] ? 4u : 0u);
                    const long long at = ((x - lowest[0]) * n[1] + (y - lowest[1])) * n[2] + (z - lowest[2]);
                    buckets[static_cast<size_t>(at)].push_back(static_cast<unsigned>(i) | flags << 29);
                }
    }
    g.cellStart.push_back(0);
    for (const std::vector<unsigned>& bucket : buckets) {
        g.entries.insert(g.entries.end(), bucket.begin(), bucket.end());
        g.cellStart.push_back(static_cast<unsigned>(g.entries.size()));
    }
    return g;
}

gpu::GridView hostView(const BuiltGrid& g, const gpu::BoxArrays& b) {
    gpu::GridView view;
    view.cell = g.cell;
    for (int a = 0; a < 3; ++a) {
        view.lo[a] = g.lo[a];
        view.hi[a] = g.hi[a];
    }
    view.cellStart = g.cellStart.data();
    view.entries = g.entries.data();
    view.boxLo = b.lo.data();
    view.boxHi = b.hi.data();
    view.nItems = b.count();
    return view;
}

BuiltGrid downloaded(const gpu::DeviceGrid& g) {
    BuiltGrid a;
    a.cell = g.view.cell;
    if (g.view.nItems == 0) return a;
    long long cells = 1;
    for (int k = 0; k < 3; ++k) {
        a.lo[k] = g.view.lo[k];
        a.hi[k] = g.view.hi[k];
        cells *= g.view.hi[k] - g.view.lo[k] + 1;
    }
    g.cellStart.download(a.cellStart, static_cast<size_t>(cells) + 1);
    g.entries.download(a.entries, a.cellStart.back());
    return a;
}

void checkSame(const BuiltGrid& a, const BuiltGrid& b) {
    CHECK(a.cell == b.cell);
    for (int k = 0; k < 3; ++k) CHECK(a.lo[k] == b.lo[k] && a.hi[k] == b.hi[k]);
    CHECK(a.cellStart == b.cellStart);
    CHECK(a.entries == b.entries);
}

} // namespace

// The device's sort gives the buckets of plain loops, in order; items with an inverted box stay out.
TEST(gpu_grid_matches_reference) {
    REQUIRE_GPU();
    gpu::BoxArrays same;
    same.cell = 3.0f;
    for (int i = 0; i < 64; ++i) {
        same.lo.insert(same.lo.end(), {-1.5f, 0.25f, 3.0f});
        same.hi.insert(same.hi.end(), {-1.0f, 0.75f, 3.5f});
    }
    gpu::BoxArrays leftOut = randomBoxes(1000, 10.0f, 2.0f, 3.0f, 13);
    for (int i = 0; i < leftOut.count(); i += 7) {
        leftOut.lo[3 * i] = 1;
        leftOut.hi[3 * i] = 0;
    }
    gpu::BoxArrays nothing = leftOut;
    for (int i = 0; i < nothing.count(); ++i) {
        nothing.lo[3 * i] = 1;
        nothing.hi[3 * i] = 0;
    }
    for (const gpu::BoxArrays& b : {randomBoxes(5000, 50.0f, 2.0f, 3.0f, 11), randomBoxes(2000, 1.0f, 0.5f, 3.0f, 12), same, leftOut, gpu::BoxArrays{}}) {
        gpu::DeviceGrid g;
        g.upload(b);
        checkSame(downloaded(g), referenceGrid(b, b.cell));
        CHECK(g.view.nItems == b.count());
    }
    gpu::DeviceGrid empty;
    empty.upload(nothing);
    CHECK(empty.view.nItems == 0);
    CHECK(gpu::BoxArrays{}.gridBytes() < 1024);
}

// A cell too small for the limit grows, and every item stays. A rebuild in the same buffers with
// smaller boxes gives the same grid as a new one.
TEST(gpu_grid_cell_grows_and_rebuilds) {
    REQUIRE_GPU();
    const gpu::BoxArrays b = randomBoxes(5000, 50.0f, 2.0f, 0.05f, 11);
    gpu::DevicePrimitives prims;
    gpu::GridScratch scratch;
    gpu::GridLimits limits;
    limits.maxCells = 4096;
    gpu::DeviceGrid g;
    g.boxLo.upload(b.lo);
    g.boxHi.upload(b.hi);
    g.build(b.count(), b.cell, prims, scratch, limits);
    long long cells = 1;
    for (int k = 0; k < 3; ++k) cells *= g.view.hi[k] - g.view.lo[k] + 1;
    std::printf("    cell 0.05 grew to %.3f: %lld cells\n", g.view.cell, cells);
    CHECK(g.view.cell > 0.05f);
    CHECK(cells <= limits.maxCells);
    checkSame(downloaded(g), referenceGrid(b, g.view.cell));
    const gpu::BoxArrays smaller = randomBoxes(300, 8.0f, 1.0f, 2.0f, 14);
    g.boxLo.upload(smaller.lo);
    g.boxHi.upload(smaller.hi);
    g.build(smaller.count(), smaller.cell, prims, scratch);
    checkSame(downloaded(g), referenceGrid(smaller, smaller.cell));
}

// More cell registrations than arrays can index fail before anything is allocated for them.
TEST(gpu_grid_capacity_fails_explicitly) {
    REQUIRE_GPU();
    const gpu::BoxArrays b = randomBoxes(1000, 10.0f, 2.0f, 1.0f, 15);
    gpu::DevicePrimitives prims;
    gpu::GridScratch scratch;
    gpu::GridLimits limits;
    limits.maxEntries = 100;
    gpu::DeviceGrid g;
    g.boxLo.upload(b.lo);
    g.boxHi.upload(b.hi);
    bool failed = false;
    try {
        g.build(b.count(), b.cell, prims, scratch, limits);
    } catch (const gpu::DeviceError& e) {
        failed = std::string(e.what()).find("capacity exceeded") != std::string::npos;
    }
    CHECK(failed);
}

// BoxArrays::gridBytes sizes a run's reservation before anything is on the device: it must cover
// what building the grid allocates at its peak, and not be far above it.
TEST(gpu_grid_bytes_cover_the_build) {
    REQUIRE_GPU();
    const gpu::BoxArrays b = randomBoxes(200000, 30.0f, 1.0f, 2.0f, 16);
    gpu::resetPeak();
    const size_t before = gpu::memStats().current;
    size_t used = 0;
    {
        gpu::DeviceGrid g;
        g.upload(b);
        used = gpu::memStats().peak - before;
    }
    std::printf("    estimate %.1f MB, peak %.1f MB\n", b.gridBytes() / 1048576.0, used / 1048576.0);
    CHECK(b.gridBytes() >= used);
    CHECK(b.gridBytes() < 3 * used);
}

// GridCursor, the query every kernel walks, against brute force: every overlapping item once, and
// nothing else. The same order on every run: a bucket's order decides ties between queries.
TEST(gpu_grid_cursor_matches_brute_force) {
    const gpu::BoxArrays boxes = randomBoxes(3000, 20.0f, 3.0f, 2.5f, 21);
    const BuiltGrid built = referenceGrid(boxes, boxes.cell);
    const gpu::GridView grid = hostView(built, boxes);
    std::mt19937 rng(22);
    std::uniform_real_distribution<float> at(-25.0f, 25.0f), size(0.0f, 8.0f);
    long long hits = 0;
    for (int q = 0; q < 2000; ++q) {
        const float lo[3] = {at(rng), at(rng), at(rng)}, hi[3] = {lo[0] + size(rng), lo[1] + size(rng), lo[2] + size(rng)};
        std::vector<int> expected;
        for (int i = 0; i < boxes.count(); ++i) {
            const float* a = boxes.lo.data() + 3 * i;
            const float* b = boxes.hi.data() + 3 * i;
            if (a[0] <= hi[0] && b[0] >= lo[0] && a[1] <= hi[1] && b[1] >= lo[1] && a[2] <= hi[2] && b[2] >= lo[2]) expected.push_back(i);
        }
        const std::vector<int> visited = gpu::itemsOverlapping(grid, lo, hi);
        CHECK(visited == gpu::itemsOverlapping(grid, lo, hi));
        std::vector<int> sorted = visited;
        std::sort(sorted.begin(), sorted.end());
        CHECK(sorted == expected);
        hits += static_cast<long long>(visited.size());
    }
    std::printf("    2000 queries, %lld items visited\n", hits);
    CHECK(hits > 1000);
    // No items, and a query entirely outside the grid.
    const float far[3] = {1000, 1000, 1000}, farther[3] = {1001, 1001, 1001};
    CHECK(gpu::itemsOverlapping(gpu::GridView{}, far, farther).empty());
    CHECK(gpu::itemsOverlapping(grid, far, farther).empty());
}

// One mesh, whoever asks for it: the same object on the device, counted once. Its other use (with
// the walk's adjacency) is a second entry, until releaseResident().
TEST(gpu_resident_mesh_uploaded_once) {
    REQUIRE_GPU();
    gpu::releaseResident();
    TriMesh m;
    m.pos = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    m.tri = {{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}};
    m.computeNormals();
    const MeshQuery query(m);
    const std::shared_ptr<gpu::DeviceMesh> first = gpu::residentMesh(query);
    const size_t once = gpu::memStats().resident;
    CHECK(once > 0);
    CHECK(gpu::residentMesh(query) == first);
    CHECK(gpu::memStats().resident == once);
    const MeshQuery sameContent(m);
    CHECK(gpu::residentMesh(sameContent) == first);
    const std::shared_ptr<gpu::DeviceMesh> walked = gpu::residentMesh(m, gpu::MeshWalk::Adjacency);
    CHECK(walked != first);
    CHECK(walked->adjacency.nbr != nullptr);
    CHECK(gpu::residentMesh(m, gpu::MeshWalk::Adjacency) == walked);
    CHECK(gpu::memStats().resident == once + walked->bytes);
    gpu::releaseResident();
    CHECK(gpu::memStats().resident == 0);
}
