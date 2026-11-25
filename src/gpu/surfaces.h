#pragma once
// The yarn surface on the device: the surface model, gate (d), the tube mesh, the plies, and the
// clearance between the parts of a scene (the CPU twins are knitting/tube/yarn_surface.h and plies.h).
// What the GPU results are: the existing algorithms and settings evaluated in float with the
// documented scheduling changes. Their verdicts are GPU diagnostic verdicts, not the CPU's certified
// PASS / FAIL / UNRESOLVED; deviations from the CPU are measured by the tests, never filtered or replaced
// by CPU results.
#include <span>
#include <string>
#include <vector>

#include "knitting/tube/plies.h"
#include "knitting/tube/yarn_surface.h"

namespace knit::gpu {

// buildSurface (yarn_surface.h) on the device, downloaded: rings in parallel, the rotation-minimizing
// frames in ordered chunks, with SurfaceParams::compression the compressed sections (the same device
// model buildTubeMesh meshes; gpu_surface_compression_gpu_vs_cpu).
struct SurfaceGpuStats {
    long long rings = 0;
    double msBuild = 0, msIndex = 0, msCheck = 0, msDownload = 0;
    std::string line() const;
};
struct SurfaceRequest {
    const CurveSet& curves;
    const PatternParams& pattern;
    const SurfaceParams& surface;
};
bool buildSurface(const SurfaceRequest& request, SurfaceModel& out, SurfaceGpuStats* st, std::string* err);
// Gate (d) on the device: buildSurface + checkSurface (folds, sampled collisions), the model
// kept on the device. A GPU diagnostic result.
bool checkSurface(const SurfaceRequest& request, const SurfaceCheckParams& p, SurfaceCheck& out, SurfaceGpuStats* st, std::string* err);

// buildTubeMesh of buildSurface on the device: the export's tubes.
bool buildTubeMesh(const SurfaceRequest& request, const TubeParams& tp, TubeMesh& m, SurfaceGpuStats* st, std::string* err);

// Cross-part clearance over a whole scene. A flat record per sampled segment, library mm; radius: its
// part's final-surface support radius; segments of one group are not compared (the parts of one
// component overlap by design, and a part never meets itself).
struct ClearanceSegment {
    Vec3 a, b;
    double radius;
    int part, group;
};
struct ClosePartPair {
    int partA, partB;  // partA < partB
    long long segments;
    double minDistance, threshold;
};
struct SceneClearanceRequest {
    std::span<const ClearanceSegment> segments;
    double clearance;
};
struct SceneClearance {
    std::vector<ClosePartPair> pairs;  // sorted by (partA, partB)
    double minDistance = 1e300;        // the nearest pair of different groups within the largest threshold
};
// Every pair of segments closer than R_a + R_b + c, in double.
bool checkSceneClearance(const SceneClearanceRequest& request, SceneClearance& out, std::string* err);

// buildPlies (plies.h) on the device: samples, the arc-length tables and the containment
// sampling; the per-ply tridiagonal spline fit stays on the host (a sequential recurrence).
struct PliesGpuStats {
    int plies = 0;
    double msSamples = 0, msSplinesHost = 0, msFinalize = 0, msContainment = 0;
    std::string line() const;
};
struct PlyRequest {
    const CurveSet& curves;
    const SurfaceModel& surface;
    const PlyGeometry& geometry;
    double radius;
};
bool buildPlies(const PlyRequest& request, PlySet& out, PliesGpuStats* st, std::string* err);

} // namespace knit::gpu
