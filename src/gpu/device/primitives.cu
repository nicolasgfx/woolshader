// GPU backend: the data-parallel steps (primitives.h). Each CUB call runs twice, first to ask how much
// scratch it needs and then with that scratch; the overloads differ only in their types.
#include "gpu/device/primitives.h"

#include <algorithm>
#include <string>

#include <cub/device/device_radix_sort.cuh>
#include <cub/device/device_scan.cuh>
#include <cub/device/device_select.cuh>

namespace knit::gpu {

unsigned long long DevicePrimitives::countsToOffsets(unsigned long long* counts, unsigned long long* offsets, long long n) {
    size_t bytes = 0;
    cub::DeviceScan::ExclusiveSum(nullptr, bytes, counts, offsets, static_cast<int>(n + 1));
    scratch_.reserve(bytes);
    bytes = scratch_.bytes();
    check(cub::DeviceScan::ExclusiveSum(scratch_.data(), bytes, counts, offsets, static_cast<int>(n + 1)), "scan");
    unsigned long long total = 0;
    check(cudaMemcpy(&total, offsets + n, sizeof total, cudaMemcpyDeviceToHost), "scan total");
    return total;
}

long long DevicePrimitives::countsToOffsets(long long* counts, long long* offsets, long long n) {
    size_t bytes = 0;
    cub::DeviceScan::ExclusiveSum(nullptr, bytes, counts, offsets, static_cast<int>(n + 1));
    scratch_.reserve(bytes);
    bytes = scratch_.bytes();
    check(cub::DeviceScan::ExclusiveSum(scratch_.data(), bytes, counts, offsets, static_cast<int>(n + 1)), "scan");
    long long total = 0;
    check(cudaMemcpy(&total, offsets + n, sizeof total, cudaMemcpyDeviceToHost), "scan total");
    return total;
}

int DevicePrimitives::countsToOffsets(int* counts, int* offsets, long long n) {
    size_t bytes = 0;
    cub::DeviceScan::ExclusiveSum(nullptr, bytes, counts, offsets, static_cast<int>(n + 1));
    scratch_.reserve(bytes);
    bytes = scratch_.bytes();
    check(cub::DeviceScan::ExclusiveSum(scratch_.data(), bytes, counts, offsets, static_cast<int>(n + 1)), "scan");
    int total = 0;
    check(cudaMemcpy(&total, offsets + n, sizeof total, cudaMemcpyDeviceToHost), "scan total");
    return total;
}

unsigned DevicePrimitives::countsToOffsets(unsigned* counts, unsigned* offsets, long long n) {
    size_t bytes = 0;
    cub::DeviceScan::ExclusiveSum(nullptr, bytes, counts, offsets, static_cast<int>(n + 1));
    scratch_.reserve(bytes);
    bytes = scratch_.bytes();
    check(cub::DeviceScan::ExclusiveSum(scratch_.data(), bytes, counts, offsets, static_cast<int>(n + 1)), "scan");
    unsigned total = 0;
    check(cudaMemcpy(&total, offsets + n, sizeof total, cudaMemcpyDeviceToHost), "scan total");
    return total;
}

void DevicePrimitives::sortKeys(const KeySort& s) {
    size_t bytes = 0;
    cub::DeviceRadixSort::SortKeys(nullptr, bytes, s.keys, s.sorted, static_cast<int>(s.n), 0, s.endBit);
    scratch_.reserve(bytes);
    bytes = scratch_.bytes();
    check(cub::DeviceRadixSort::SortKeys(scratch_.data(), bytes, s.keys, s.sorted, static_cast<int>(s.n), 0, s.endBit), "sort");
}

void DevicePrimitives::sortPairs(const KeySort& s, int* values, int* sortedValues) {
    size_t bytes = 0;
    cub::DeviceRadixSort::SortPairs(nullptr, bytes, s.keys, s.sorted, values, sortedValues, static_cast<int>(s.n), 0, s.endBit);
    scratch_.reserve(bytes);
    bytes = scratch_.bytes();
    check(cub::DeviceRadixSort::SortPairs(scratch_.data(), bytes, s.keys, s.sorted, values, sortedValues, static_cast<int>(s.n), 0, s.endBit), "sort");
}

void DevicePrimitives::sortPairs(const KeySort& s, double* values, double* sortedValues) {
    size_t bytes = 0;
    cub::DeviceRadixSort::SortPairs(nullptr, bytes, s.keys, s.sorted, values, sortedValues, static_cast<int>(s.n), 0, s.endBit);
    scratch_.reserve(bytes);
    bytes = scratch_.bytes();
    check(cub::DeviceRadixSort::SortPairs(scratch_.data(), bytes, s.keys, s.sorted, values, sortedValues, static_cast<int>(s.n), 0, s.endBit), "sort");
}

long long DevicePrimitives::selectFlagged(int* in, unsigned char* flags, int* out, long long n) {
    selected_.reserve(sizeof(long long));
    size_t bytes = 0;
    cub::DeviceSelect::Flagged(nullptr, bytes, in, flags, out, selected_.longs(), static_cast<int>(n));
    scratch_.reserve(bytes);
    bytes = scratch_.bytes();
    check(cub::DeviceSelect::Flagged(scratch_.data(), bytes, in, flags, out, selected_.longs(), static_cast<int>(n)), "select");
    long long total = 0;
    selected_.download(&total, sizeof total);
    return total;
}

long long DevicePrimitives::selectFlagged(unsigned long long* in, unsigned char* flags, unsigned long long* out, long long n) {
    selected_.reserve(sizeof(long long));
    size_t bytes = 0;
    cub::DeviceSelect::Flagged(nullptr, bytes, in, flags, out, selected_.longs(), static_cast<int>(n));
    scratch_.reserve(bytes);
    bytes = scratch_.bytes();
    check(cub::DeviceSelect::Flagged(scratch_.data(), bytes, in, flags, out, selected_.longs(), static_cast<int>(n)), "select");
    long long total = 0;
    selected_.download(&total, sizeof total);
    return total;
}

void DevicePrimitives::reserveScratch(long long maxItems) { scratch_.reserve(largestScratch(maxItems)); }

// What DeviceBuffer::reserve allocates for it: 25 % more.
size_t DevicePrimitives::scratchBytesFor(long long maxItems) {
    const size_t most = largestScratch(maxItems);
    return std::max<size_t>(most + most / 4, 64);
}

// The largest scratch any on-stream step needs for maxItems items.
size_t DevicePrimitives::largestScratch(long long maxItems) {
    const int n = static_cast<int>(maxItems + 1);
    size_t most = 16, bytes = 0;
    unsigned long long* u64 = nullptr;
    int* i32 = nullptr;
    unsigned* u32 = nullptr;
    double* f64 = nullptr;
    cub::DeviceScan::ExclusiveSum(nullptr, bytes, u64, u64, n);
    most = std::max(most, bytes);
    cub::DeviceScan::ExclusiveSum(nullptr, bytes, i32, i32, n);
    most = std::max(most, bytes);
    cub::DeviceScan::ExclusiveSum(nullptr, bytes, u32, u32, n);
    most = std::max(most, bytes);
    cub::DeviceRadixSort::SortKeys(nullptr, bytes, u64, u64, n, 0, 64);
    most = std::max(most, bytes);
    cub::DoubleBuffer<unsigned long long> keys(u64, u64);
    cub::DoubleBuffer<double> values(f64, f64);
    cub::DeviceRadixSort::SortPairs(nullptr, bytes, keys, values, n, 0, 64);
    most = std::max(most, bytes);
    return most;
}

void DevicePrimitives::reserveScratch(const ScratchPlan& plan) { scratch_.reserve(largestScratch(plan)); }

// The largest scratch of the plan's steps, each at its own list length.
size_t DevicePrimitives::largestScratch(const ScratchPlan& plan) {
    size_t most = 16, bytes = 0;
    unsigned long long* u64 = nullptr;
    int* i32 = nullptr;
    unsigned* u32 = nullptr;
    double* f64 = nullptr;
    const int scan = static_cast<int>(plan.scanItems + 1);
    cub::DeviceScan::ExclusiveSum(nullptr, bytes, u64, u64, scan);
    most = std::max(most, bytes);
    cub::DeviceScan::ExclusiveSum(nullptr, bytes, i32, i32, scan);
    most = std::max(most, bytes);
    cub::DeviceScan::ExclusiveSum(nullptr, bytes, u32, u32, scan);
    most = std::max(most, bytes);
    cub::DeviceRadixSort::SortKeys(nullptr, bytes, u64, u64, static_cast<int>(plan.keySortItems + 1), 0, 64);
    most = std::max(most, bytes);
    cub::DoubleBuffer<unsigned long long> keys(u64, u64);
    cub::DoubleBuffer<double> values(f64, f64);
    cub::DeviceRadixSort::SortPairs(nullptr, bytes, keys, values, static_cast<int>(plan.pairSortItems + 1), 0, 64);
    most = std::max(most, bytes);
    return most;
}

void DevicePrimitives::needScratch(size_t bytes) const {
    if (bytes > scratch_.bytes()) throw DeviceError("primitives: the reserved scratch (" + std::to_string(scratch_.bytes()) + " bytes) is smaller than " + std::to_string(bytes));
}

void DevicePrimitives::countsToOffsetsOnStream(const unsigned long long* counts, unsigned long long* offsets, long long n, cudaStream_t stream) {
    size_t bytes = 0;
    cub::DeviceScan::ExclusiveSum(nullptr, bytes, counts, offsets, static_cast<int>(n + 1), stream);
    needScratch(bytes);
    bytes = scratch_.bytes();
    check(cub::DeviceScan::ExclusiveSum(scratch_.data(), bytes, counts, offsets, static_cast<int>(n + 1), stream), "scan");
}

void DevicePrimitives::countsToOffsetsOnStream(const int* counts, int* offsets, long long n, cudaStream_t stream) {
    size_t bytes = 0;
    cub::DeviceScan::ExclusiveSum(nullptr, bytes, counts, offsets, static_cast<int>(n + 1), stream);
    needScratch(bytes);
    bytes = scratch_.bytes();
    check(cub::DeviceScan::ExclusiveSum(scratch_.data(), bytes, counts, offsets, static_cast<int>(n + 1), stream), "scan");
}

void DevicePrimitives::countsToOffsetsOnStream(const unsigned* counts, unsigned* offsets, long long n, cudaStream_t stream) {
    size_t bytes = 0;
    cub::DeviceScan::ExclusiveSum(nullptr, bytes, counts, offsets, static_cast<int>(n + 1), stream);
    needScratch(bytes);
    bytes = scratch_.bytes();
    check(cub::DeviceScan::ExclusiveSum(scratch_.data(), bytes, counts, offsets, static_cast<int>(n + 1), stream), "scan");
}

void DevicePrimitives::sortKeysOnStream(const KeySort& s, cudaStream_t stream) {
    size_t bytes = 0;
    cub::DeviceRadixSort::SortKeys(nullptr, bytes, s.keys, s.sorted, static_cast<int>(s.n), 0, s.endBit, stream);
    needScratch(bytes);
    bytes = scratch_.bytes();
    check(cub::DeviceRadixSort::SortKeys(scratch_.data(), bytes, s.keys, s.sorted, static_cast<int>(s.n), 0, s.endBit, stream), "sort");
}

bool DevicePrimitives::sortPairsAlternatingOnStream(const PairSort& s, cudaStream_t stream) {
    cub::DoubleBuffer<unsigned long long> keys(s.keys, s.keysAlt);
    cub::DoubleBuffer<double> values(s.values, s.valuesAlt);
    size_t bytes = 0;
    cub::DeviceRadixSort::SortPairs(nullptr, bytes, keys, values, static_cast<int>(s.n), 0, s.endBit, stream);
    needScratch(bytes);
    bytes = scratch_.bytes();
    check(cub::DeviceRadixSort::SortPairs(scratch_.data(), bytes, keys, values, static_cast<int>(s.n), 0, s.endBit, stream), "sort");
    return keys.Current() == s.keysAlt;
}

} // namespace knit::gpu
