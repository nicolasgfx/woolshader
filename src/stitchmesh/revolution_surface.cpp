#include "stitchmesh/revolution_surface.h"

#include <algorithm>
#include <cmath>

namespace knit {

namespace {

// splitmix64 (Steele, Lea, Flood 2014): a small, fixed, platform-independent generator for variant J.
struct SplitMix64 {
    uint64_t state = 0;
    uint64_t next() {
        uint64_t z = (state += 0x9E3779B97F4A7C15ull);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
        return z ^ (z >> 31);
    }
    double uniform(double lo, double hi) {   // [lo, hi)
        return lo + (hi - lo) * (double) (next() >> 11) * (1.0 / 9007199254740992.0);
    }
};

double interpolate(const std::vector<double>& xs, const std::vector<double>& ys, double x) {
    if (x <= xs.front()) return ys.front();
    if (x >= xs.back()) return ys.back();
    const size_t i = (size_t) (std::upper_bound(xs.begin(), xs.end(), x) - xs.begin());
    const double f = (x - xs[i - 1]) / (xs[i] - xs[i - 1]);
    return ys[i - 1] + f * (ys[i] - ys[i - 1]);
}

} // namespace

double RevolutionProfile::radius(double z) const {
    if (z <= z0) return r0;
    if (z >= z1) return r1;
    const double t = (z - z0) / (z1 - z0);
    return r1 + (r0 - r1) * (1.0 + std::cos(kPi * t)) / 2.0;
}

double RevolutionProfile::slope(double z) const {
    if (z <= z0 || z >= z1) return 0;
    return -(r0 - r1) * kPi / (2.0 * (z1 - z0)) * std::sin(kPi * (z - z0) / (z1 - z0));
}

double RevolutionProfile::secondDerivative(double z) const {
    if (z <= z0 || z >= z1) return 0;
    const double L = z1 - z0;
    return -(r0 - r1) * kPi * kPi / (2.0 * L * L) * std::cos(kPi * (z - z0) / L);
}

void MeridianTable::build(const RevolutionProfile& p, int samples) {
    z.assign((size_t) samples, 0);
    s.assign((size_t) samples, 0);
    double prev = std::sqrt(1.0 + p.slope(0) * p.slope(0));
    for (int i = 1; i < samples; ++i) {
        z[i] = p.z2 * i / (samples - 1);
        const double cur = std::sqrt(1.0 + p.slope(z[i]) * p.slope(z[i]));
        s[i] = s[i - 1] + 0.5 * (prev + cur) * (z[i] - z[i - 1]);
        prev = cur;
    }
}

double MeridianTable::arcAt(double zq) const { return interpolate(z, s, zq); }
double MeridianTable::heightAt(double sq) const { return interpolate(s, z, sq); }

PrincipalCurvatures principalCurvatures(const RevolutionProfile& p, double z) {
    const double d1 = p.slope(z), d2 = p.secondDerivative(z);
    const double g = std::sqrt(1.0 + d1 * d1);
    PrincipalCurvatures k;
    k.meridian = -d2 / (g * g * g);
    k.parallel = 1.0 / (p.radius(z) * g);
    return k;
}

Vec3 revolutionNormal(const RevolutionProfile& p, double z, double phi) {
    return normalized(Vec3{std::cos(phi), std::sin(phi), -p.slope(z)});
}

bool makeRevolutionMesh(const RevolutionProfile& p, const RevolutionMeshParams& params, RevolutionMesh& out, std::string* err) {
    out = RevolutionMesh{};
    if (params.segments < 3 || !(params.ringSpacing > 0) || params.tableSamples < 2) {
        if (err) *err = "revolution mesh: needs at least 3 segments, a positive ring spacing and a meridian table";
        return false;
    }
    MeridianTable table;
    table.build(p, params.tableSamples);
    const double T = table.length();
    const int bands = std::max(1, (int) std::lround(T / params.ringSpacing));
    const int M = params.segments;
    out.rings = bands + 1;
    out.meridianLength = T;
    SplitMix64 rng{params.jitterSeed};
    for (int i = 0; i <= bands; ++i) {
        for (int j = 0; j < M; ++j) {
            double s = T * i / bands;
            double phi = (2.0 * kPi * j + kPi * i) / M;
            if (params.jitter > 0 && i > 0 && i < bands) {
                const double ds = rng.uniform(-params.jitter, params.jitter);
                const double dt = rng.uniform(-params.jitter, params.jitter);
                s += ds;
                phi += dt / p.radius(table.heightAt(s));
            }
            const double z = table.heightAt(s);
            const double r = p.radius(z);
            out.mesh.pos.push_back(Vec3{r * std::cos(phi), r * std::sin(phi), z});
            out.arc.push_back(s);
            out.height.push_back(z);
            out.phi.push_back(phi);
        }
    }
    // Band i: lower ring i, upper ring i + 1, whose vertex j sits half a segment past lower vertex j.
    // (L_j, L_j+1, U_j) and (U_j, L_j+1, U_j+1): counter-clockwise seen from outside (outward normals).
    for (int i = 0; i < bands; ++i)
        for (int j = 0; j < M; ++j) {
            const int lj = i * M + j, lj1 = i * M + (j + 1) % M;
            const int uj = (i + 1) * M + j, uj1 = (i + 1) * M + (j + 1) % M;
            out.mesh.tri.push_back({lj, lj1, uj});
            out.mesh.tri.push_back({uj, lj1, uj1});
        }
    for (int j = 0; j < M; ++j) {
        out.bottomLoop.push_back(j);
        out.topLoop.push_back(bands * M + j);
    }
    out.mesh.computeNormals();
    // exact area: 2 pi * integral of r(z) sqrt(1 + r'^2) dz over the table's samples
    double area = 0;
    for (size_t k = 1; k < table.z.size(); ++k) area += 0.5 * (p.radius(table.z[k]) + p.radius(table.z[k - 1])) * (table.s[k] - table.s[k - 1]);
    out.exactArea = 2.0 * kPi * area;
    return true;
}

} // namespace knit
