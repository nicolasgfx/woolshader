#include "model/spatial_hash.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace knit {

namespace {

constexpr int kFlagShift = 29;
constexpr std::uint32_t kIdMask = (1u << kFlagShift) - 1;
constexpr long long kDenseCells = 1 << 22;

// A cell's key: 21 bits per axis, offset so that cells in +-1e6 stay positive.
constexpr std::int64_t kKeyOffset = 1 << 20;
constexpr std::int64_t kKeyMask = 0x1FFFFF;
std::int64_t cellKey(std::int64_t i, std::int64_t j, std::int64_t k) {
    return ((i + kKeyOffset) & kKeyMask) | (((j + kKeyOffset) & kKeyMask) << 21) | (((k + kKeyOffset) & kKeyMask) << 42);
}
std::int64_t keyI(std::int64_t key) { return (key & kKeyMask) - kKeyOffset; }
std::int64_t keyJ(std::int64_t key) { return ((key >> 21) & kKeyMask) - kKeyOffset; }
std::int64_t keyK(std::int64_t key) { return ((key >> 42) & kKeyMask) - kKeyOffset; }

} // namespace

std::int64_t SpatialHash3::cellOf(double v) const { return static_cast<std::int64_t>(std::floor(v / cell_)); }

size_t SpatialHash3::denseIndex(std::int64_t i, std::int64_t j, std::int64_t k) const {
    const std::int64_t nj = occHi_[1] - occLo_[1] + 1, nk = occHi_[2] - occLo_[2] + 1;
    return static_cast<size_t>(((i - occLo_[0]) * nj + (j - occLo_[1])) * nk + (k - occLo_[2]));
}

const SpatialHash3::Bucket* SpatialHash3::bucketAt(std::int64_t i, std::int64_t j, std::int64_t k) const {
    if (!dense_.empty()) {
        const std::int32_t b = dense_[denseIndex(i, j, k)];
        return b >= 0 ? &denseBuckets_[b] : nullptr;
    }
    const auto it = cells_.find(cellKey(i, j, k));
    return it != cells_.end() ? &it->second : nullptr;
}

int SpatialHash3::insert(const Vec3& lo, const Vec3& hi) {
    if (!dense_.empty()) throw std::logic_error("SpatialHash3::insert after finalize");
    const int id = static_cast<int>(lo_.size());
    lo_.push_back(lo);
    hi_.push_back(hi);
    const std::int64_t i0 = cellOf(lo.x), j0 = cellOf(lo.y), k0 = cellOf(lo.z);
    const std::int64_t i1 = cellOf(hi.x), j1 = cellOf(hi.y), k1 = cellOf(hi.z);
    for (std::int64_t i = i0; i <= i1; ++i)
        for (std::int64_t j = j0; j <= j1; ++j)
            for (std::int64_t k = k0; k <= k1; ++k) {
                const std::uint32_t lowest = (i == i0 ? 1u : 0u) | (j == j0 ? 2u : 0u) | (k == k0 ? 4u : 0u);
                cells_[cellKey(i, j, k)].push_back(static_cast<std::uint32_t>(id) | lowest << kFlagShift);
            }
    occLo_[0] = std::min(occLo_[0], i0);
    occLo_[1] = std::min(occLo_[1], j0);
    occLo_[2] = std::min(occLo_[2], k0);
    occHi_[0] = std::max(occHi_[0], i1);
    occHi_[1] = std::max(occHi_[1], j1);
    occHi_[2] = std::max(occHi_[2], k1);
    return id;
}

void SpatialHash3::finalize() {
    if (cells_.empty()) return;
    long long n = 1;
    for (int a = 0; a < 3; ++a) n *= occHi_[a] - occLo_[a] + 1;
    if (n > kDenseCells) return;
    dense_.assign(static_cast<size_t>(n), -1);
    for (auto& [key, bucket] : cells_) {
        dense_[denseIndex(keyI(key), keyJ(key), keyK(key))] = static_cast<std::int32_t>(denseBuckets_.size());
        denseBuckets_.push_back(std::move(bucket));
    }
    cells_.clear();
}

void SpatialHash3::query(const Vec3& lo, const Vec3& hi, const std::function<void(int)>& visit) const {
    if (lo_.empty()) return;
    const std::int64_t qi = cellOf(lo.x), qj = cellOf(lo.y), qk = cellOf(lo.z);
    const std::int64_t i0 = std::max(qi, occLo_[0]), i1 = std::min(cellOf(hi.x), occHi_[0]);
    const std::int64_t j0 = std::max(qj, occLo_[1]), j1 = std::min(cellOf(hi.y), occHi_[1]);
    const std::int64_t k0 = std::max(qk, occLo_[2]), k1 = std::min(cellOf(hi.z), occHi_[2]);
    for (std::int64_t i = i0; i <= i1; ++i)
        for (std::int64_t j = j0; j <= j1; ++j)
            for (std::int64_t k = k0; k <= k1; ++k) {
                const Bucket* bucket = bucketAt(i, j, k);
                if (!bucket) continue;
                const std::uint32_t atQueryLow = (i == qi ? 1u : 0u) | (j == qj ? 2u : 0u) | (k == qk ? 4u : 0u);
                for (std::uint32_t e : *bucket) {
                    if (((e >> kFlagShift) | atQueryLow) != 7u) continue;  // not the item's first shared cell
                    const int id = static_cast<int>(e & kIdMask);
                    const Vec3 &a = lo_[id], &b = hi_[id];
                    if (a.x > hi.x || b.x < lo.x || a.y > hi.y || b.y < lo.y || a.z > hi.z || b.z < lo.z) continue;
                    visit(id);
                }
            }
}

} // namespace knit
