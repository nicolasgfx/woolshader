#pragma once
// What every kernel starts with. One thread handles one item: a launch of n items runs
// blocksFor(n) blocks of kThreadsPerBlock threads, and the threads past the last item return at once.
//     const int i = threadItem();
//     if (i >= n) return;
#include <cuda_runtime.h>

namespace knit::gpu {

constexpr int kThreadsPerBlock = 128;
inline unsigned blocksFor(long long items, int threads = kThreadsPerBlock) {
    return static_cast<unsigned>((items + threads - 1) / threads);
}

__device__ __forceinline__ int threadItem() { return blockIdx.x * blockDim.x + threadIdx.x; }
__device__ __forceinline__ long long threadItem64() { return static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x; }

// The abort flag (runtime.h abortFlagDevice), polled every 1024 steps of a long per-thread loop.
__device__ __forceinline__ bool abortRequested(const volatile int* flag, unsigned step) { return (step & 1023u) == 0 && *flag != 0; }

// Non-negative doubles order like their bit patterns: max / min as integer atomics, so the result
// does not depend on the order the threads arrive in.
__device__ __forceinline__ void atomicMaxPositive(double* a, double v) {
    if (!(v > 0)) return;
    atomicMax(reinterpret_cast<unsigned long long*>(a), static_cast<unsigned long long>(__double_as_longlong(v)));
}
__device__ __forceinline__ void atomicMinPositive(double* a, double v) {
    if (v < 0) v = 0;
    atomicMin(reinterpret_cast<unsigned long long*>(a), static_cast<unsigned long long>(__double_as_longlong(v)));
}

} // namespace knit::gpu
