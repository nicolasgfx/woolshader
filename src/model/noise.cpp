// Seeded gradient noise and its fractal sum (noise.h has the method and its source).
#include "model/noise.h"

#include <algorithm>
#include <cmath>

#include "model/hash.h"

namespace knit {

namespace {

// Perlin's quintic fade: 0 at 0, 1 at 1, first and second derivatives 0 at both ends.
double fade(double t) { return t * t * t * (t * (t * 6 - 15) + 10); }

double lerp(double a, double b, double t) { return a + (b - a) * t; }

// The unit gradient at lattice corner (i, j, k): a uniform direction on the sphere from two hashes.
Vec3 cornerGradient(int i, int j, int k, uint64_t seed) {
    const uint64_t cell = splitmix64(seed ^ splitmix64(uint64_t(uint32_t(k))));
    const double z = 2 * hash01(cell, i, j, 0) - 1, phi = 2 * kPi * hash01(cell, i, j, 1);
    const double r = std::sqrt(std::max(0.0, 1 - z * z));
    return {r * std::cos(phi), r * std::sin(phi), z};
}

} // namespace

double gradientNoise(const Vec3& p, uint64_t seed) {
    const double fx = std::floor(p.x), fy = std::floor(p.y), fz = std::floor(p.z);
    const int i = static_cast<int>(fx), j = static_cast<int>(fy), k = static_cast<int>(fz);
    const Vec3 f{p.x - fx, p.y - fy, p.z - fz};
    // The eight corners' contributions, blended along x, then y, then z.
    double c[2][2][2];
    for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b)
            for (int e = 0; e < 2; ++e) c[a][b][e] = dot(cornerGradient(i + a, j + b, k + e, seed), f - Vec3{double(a), double(b), double(e)});
    const double u = fade(f.x), v = fade(f.y), w = fade(f.z);
    const double x00 = lerp(c[0][0][0], c[1][0][0], u), x10 = lerp(c[0][1][0], c[1][1][0], u);
    const double x01 = lerp(c[0][0][1], c[1][0][1], u), x11 = lerp(c[0][1][1], c[1][1][1], u);
    // With unit gradients the sum reaches at most sqrt(3) / 2 in magnitude; scaling by 2 / sqrt(3) bounds it by 1.
    return lerp(lerp(x00, x10, v), lerp(x01, x11, v), w) * (2 / std::sqrt(3.0));
}

double fractalNoise(const Vec3& p, uint64_t seed, int octaves) {
    double sum = 0, amplitude = 1, total = 0, frequency = 1;
    for (int o = 0; o < octaves; ++o) {
        sum += amplitude * gradientNoise(p * frequency, splitmix64(seed + uint64_t(o)));
        total += amplitude;
        amplitude *= 0.5;
        frequency *= 2;
    }
    return total > 0 ? sum / total : 0;
}

} // namespace knit
