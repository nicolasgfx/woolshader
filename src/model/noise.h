#pragma once
// Deterministic 3D gradient noise and its fractal sum, seeded by a hash: the spatially correlated fields
// of the Fur workflow's irregularities (patches of density and length that neighbouring fibres share).
//
// gradientNoise: Perlin's improved noise ("Improving Noise", SIGGRAPH 2002) with each lattice corner's
// gradient drawn from a hash of (seed, corner) instead of a permutation table: the quintic fade
// 6t^5 - 15t^4 + 10t^3 and the dot product of each corner's unit gradient with the offset to it, scaled by
// 2 / sqrt(3) (unit gradients bound the raw sum by sqrt(3) / 2), so its value lies in [-1, 1]. It is 0 at
// every lattice point and continuous with a continuous gradient.
//
// fractalNoise: octaves of gradientNoise at doubling frequency and halving amplitude (fBm), divided by the
// sum of the amplitudes, so it lies in [-1, 1] too. Its feature size is about the unit length: scale the
// point (p / size) for a field of that size.
#include <cstdint>

#include "model/vec.h"

namespace knit {

double gradientNoise(const Vec3& p, uint64_t seed);
double fractalNoise(const Vec3& p, uint64_t seed, int octaves);

} // namespace knit
