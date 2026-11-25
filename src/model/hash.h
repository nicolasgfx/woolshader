#pragma once
// FNV-1a 64: the identities of cached results. Stored keys depend on the starting value, so both
// values in use stay: kFnvBasis (fnv1a64: knit digests, cache checksums, mesh digests) and the
// standard basis (fnv1a: asset content hashes, stage keys). splitmix64 and hash01: the deterministic
// draws of the fibre generators, fur and hair, each a pure function of its identity.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>

#include "model/vec.h"

namespace knit {

inline constexpr uint64_t kFnvBasis = 1469598103934665603ull;
inline constexpr uint64_t kFnvStandardBasis = 0xcbf29ce484222325ull;

inline uint64_t fnv1a64Bytes(const void* data, size_t bytes, uint64_t h = kFnvBasis) {
    const unsigned char* p = static_cast<const unsigned char*>(data);
    for (size_t i = 0; i < bytes; ++i) h = (h ^ p[i]) * 1099511628211ull;
    return h;
}

inline uint64_t fnv1a64(const std::string& text, uint64_t h = kFnvBasis) { return fnv1a64Bytes(text.data(), text.size(), h); }

inline uint64_t fnv1a(const void* data, size_t bytes, uint64_t h = kFnvStandardBasis) { return fnv1a64Bytes(data, bytes, h); }

inline uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ull;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
    return x ^ (x >> 31);
}

// Uniform [0, 1) from (seed, a, b, channel). The GPU kernels hold the same construction.
inline double hash01(uint64_t seed, int a, int b, int channel) {
    uint64_t h = splitmix64(seed);
    h = splitmix64(h ^ (uint64_t(uint32_t(a)) << 32 | uint32_t(b)));
    h = splitmix64(h ^ uint64_t(channel));
    return double(h >> 11) * (1.0 / 9007199254740992.0);
}

// A standard normal by Box-Muller from two channels, clamped to +-limit.
inline double hashNormal(uint64_t seed, int a, int b, int channel0, int channel1, double limit) {
    const double u1 = 1.0 - hash01(seed, a, b, channel0), u2 = hash01(seed, a, b, channel1);
    return std::clamp(std::sqrt(-2 * std::log(u1)) * std::cos(2 * kPi * u2), -limit, limit);
}

} // namespace knit
