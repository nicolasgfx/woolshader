#pragma once
// The simple reference surface of the stockinette stitch-mesh experiment (recipe
// tube-taper-v1): an open surface of revolution about +z, its exact profile and meridian arc length, and its
// triangle mesh. The exact quantities (arc length per vertex, normals, curvatures) are what the analysis of
// checkpoint 1 measures the computed fields against.
#include <cstdint>
#include <string>
#include <vector>

#include "model/mesh.h"

namespace knit {

// r(z): r0 up to z0, a cosine taper to r1 between z0 and z1 (C1: slope 0 at both ends), r1 up to z2. All in mm.
struct RevolutionProfile {
    double r0 = 30, r1 = 20, z0 = 60, z1 = 120, z2 = 150;
    double radius(double z) const;
    double slope(double z) const;              // dr/dz
    double secondDerivative(double z) const;   // d2r/dz2
};

// The meridian arc length s(z) from z = 0, tabulated by the trapezoid rule on sqrt(1 + r'^2), and its inverse.
struct MeridianTable {
    std::vector<double> z, s;
    void build(const RevolutionProfile& p, int samples);
    double length() const { return s.empty() ? 0 : s.back(); }
    double arcAt(double zq) const;       // linear interpolation in the table
    double heightAt(double sq) const;
};

// Principal curvatures at height z [1/mm], positive where the surface bends away from its outward normal: along the
// meridian (-r'' / (1 + r'^2)^1.5) and along the parallel (1 / (r sqrt(1 + r'^2))).
struct PrincipalCurvatures {
    double meridian = 0, parallel = 0;
};
PrincipalCurvatures principalCurvatures(const RevolutionProfile& p, double z);

// The outward unit normal at height z and azimuth phi.
Vec3 revolutionNormal(const RevolutionProfile& p, double z, double phi);

struct RevolutionMeshParams {
    int segments = 192;          // vertices per ring
    double ringSpacing = 1.0;    // meridian arc length between rings [mm] (rounded to fill the meridian exactly)
    int tableSamples = 600001;   // the meridian table
    // Variant J: every vertex off the two boundary rings moves by a uniform random offset in [-jitter, jitter] along
    // the meridian and along the parallel (each in mm), then lies on the exact surface again. jitter 0: no variant.
    double jitter = 0;
    uint64_t jitterSeed = 0;
};

struct RevolutionMesh {
    TriMesh mesh;                     // ring i, vertex j at angle (2 pi j + i pi) / segments: every ring turned half a
                                      // segment against the one below, so the mesh's edges follow no row direction
    std::vector<double> arc, height, phi;   // per vertex: exact meridian arc length from z = 0, z, azimuth
    std::vector<int> bottomLoop, topLoop;   // the boundary rings' vertices in +phi order (z = 0, z = z2)
    int rings = 0;                    // ring count (bands + 1)
    double meridianLength = 0;        // T [mm]
    double exactArea = 0;             // 2 pi * integral of r ds [mm2]
};

// False only for unusable parameters (fewer than 3 segments, no spacing).
bool makeRevolutionMesh(const RevolutionProfile& p, const RevolutionMeshParams& params, RevolutionMesh& out, std::string* err);

} // namespace knit
