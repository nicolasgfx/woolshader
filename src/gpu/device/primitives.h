#pragma once
// GPU backend: the data-parallel steps every list on the device is built from. A list of unknown
// length takes three steps: a kernel counts per item, countsToOffsets turns the counts into offsets,
// a kernel fills. CUB's scans and sorts need scratch memory; here it comes from a tracked
// DeviceBuffer, so it counts in the budget and under the cap (CUB's own single-call overloads would
// take it from a memory resource the budget cannot see).
#include <cstddef>

#include <cuda_runtime.h>

#include "gpu/device/runtime.h"

namespace knit::gpu {

// Sums use a fixed grid of blocks, so a rerun adds the same numbers in the same order: no floating
// point atomics.
constexpr int kSumBlocks = 256, kSumThreads = 256;

// Keys with double values, sorted through a double buffer: the second arrays are where CUB writes
// alternately, so the scratch stays small (sortPairsAlternatingOnStream).
struct PairSort {
    unsigned long long* keys;
    unsigned long long* keysAlt;
    double* values;
    double* valuesAlt;
    long long n;
    int endBit = 64;
};

struct KeySort {
    unsigned long long* keys;
    unsigned long long* sorted;
    long long n;
    int endBit = 64;  // keys are sorted on bits [0, endBit)
};

// The longest list each kind of on-stream step will be given, when the kinds see lists of very different
// lengths (ncp-gpu's live contact storage scans a long row list but sorts only a short table): the
// scratch is the largest any one of them needs, not every kind at the longest length.
struct ScratchPlan {
    long long scanItems = 0;      // countsToOffsetsOnStream (every element type)
    long long keySortItems = 0;   // sortKeysOnStream
    long long pairSortItems = 0;  // sortPairsAlternatingOnStream
};

class DevicePrimitives {
public:
    // offsets[i] = counts[0] + ... + counts[i - 1], for i in [0, n]; returns offsets[n]. Both arrays
    // hold n + 1 entries.
    unsigned long long countsToOffsets(unsigned long long* counts, unsigned long long* offsets, long long n);
    long long countsToOffsets(long long* counts, long long* offsets, long long n);
    int countsToOffsets(int* counts, int* offsets, long long n);
    unsigned countsToOffsets(unsigned* counts, unsigned* offsets, long long n);

    void sortKeys(const KeySort& sort);
    void sortPairs(const KeySort& sort, int* values, int* sortedValues);
    void sortPairs(const KeySort& sort, double* values, double* sortedValues);

    // Copies the flagged items to `out` in order; returns how many.
    long long selectFlagged(int* in, unsigned char* flags, int* out, long long n);
    long long selectFlagged(unsigned long long* in, unsigned char* flags, unsigned long long* out, long long n);

    // The same steps for an iteration captured into a graph (ncp-gpu): enqueued on `stream`, with no
    // host read and no allocation. reserveScratch reserves, once and beforehand, the scratch for the
    // largest number of items they will be given.
    void reserveScratch(long long maxItems);
    void reserveScratch(const ScratchPlan& plan);
    // The bytes reserveScratch(maxItems) would hold (no allocation): what a capacity plan counts.
    static size_t scratchBytesFor(long long maxItems);
    void countsToOffsetsOnStream(const unsigned long long* counts, unsigned long long* offsets, long long n, cudaStream_t stream);
    void countsToOffsetsOnStream(const int* counts, int* offsets, long long n, cudaStream_t stream);
    void countsToOffsetsOnStream(const unsigned* counts, unsigned* offsets, long long n, cudaStream_t stream);
    void sortKeysOnStream(const KeySort& sort, cudaStream_t stream);
    // Returns true when the sorted pairs are in keysAlt / valuesAlt (decided when enqueued).
    bool sortPairsAlternatingOnStream(const PairSort& sort, cudaStream_t stream);

    size_t scratchBytes() const { return scratch_.bytes(); }  // the memory report

private:
    static size_t largestScratch(long long maxItems);
    static size_t largestScratch(const ScratchPlan& plan);
    void needScratch(size_t bytes) const;  // the reserved scratch is too small: a DeviceError
    DeviceBuffer scratch_, selected_;
};

} // namespace knit::gpu
