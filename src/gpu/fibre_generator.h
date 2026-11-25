#pragma once
// The fibre generator on the device (the CPU twin is fibres/fibre_generator.h).
// What the GPU results are: the existing algorithms and settings evaluated in float with the
// documented scheduling changes. Their verdicts are GPU diagnostic verdicts, not the CPU's certified
// PASS / FAIL / UNRESOLVED; deviations from the CPU are measured by the tests, never filtered or replaced
// by CPU results.
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "fibres/fibre.h"
#include "fibres/fibre_record.h"
#include "model/mesh_query.h"

namespace knit::gpu {

// Flat replacements of FibreFrame's callbacks (fibre.h): the GPU evaluates them itself.
struct GpuFibreInputs;  // the flattened host inputs of a run (fibres.cu)
struct FibreFrameGpu {
    // Front / exposure normal (FibreFrame::outward): the closest point of this mesh (library
    // mm) with its interpolated vertex normals (Chart / MeshSurface convention). nullptr: the
    // flat patch (+z, unbounded exposure).
    const TriMesh* outward = nullptr;
    double exposureLength = 0;
    double backDepth = 0;  // FibreFrame::backDepth: the fabric's back lies this far behind the mesh
    // FibreFrame::hidden: inside any of these meshes (MeshQuery pseudo-normal convention), in
    // the common frame q = p * hiddenScale + hiddenOffset of the library point p.
    std::vector<const MeshQuery*> hidden;
    double hiddenScale = 1;
    Vec3 hiddenOffset;
    // FibreFrame::quantize: the export's rounding (library -> scene -> float32 -> library).
    bool quantize = false;
    double quantSigma = 1;
    Vec3 quantCentre;
    uint64_t pieceSeed = 0;
    // Plies (FibreFrame::plies / parentCurves / parentSurface / candidates).
    const PlySet* plies = nullptr;
    const CurveSet* parentCurves = nullptr;
    const SurfaceModel* parentSurface = nullptr;
    std::vector<int> candidates;
    const std::vector<std::vector<IslandWeight>>* islands = nullptr;  // FibreFrame::islands; null or empty: none
    std::function<void()> checkpoint;  // between batches (cancel.h: may throw Cancelled)
    // After each batch: the covered fraction of the candidates in [0, 1], monotonic (kernel B's
    // survivors of a kernel-A batch fill that batch's share of the candidate range).
    std::function<void(double)> progress;
    // Budget order (fibre_budget.h): only the candidates whose float key lies in [keyLo, keyHi);
    // the others are not candidates of the run (FibreReport and the outcomes leave them out).
    float keyLo = 0.0f, keyHi = 2.0f;
    // Set: the flattened inputs are kept here by the first call and reused by the later ones (the
    // caller guarantees the same curves, surface model and meshes, e.g. one part's budget windows).
    std::shared_ptr<GpuFibreInputs>* inputs = nullptr;
};

// Builds the GPU frame from a CPU FibreFrame's non-callback fields (setupPieceFibres).
FibreFrameGpu frameFrom(const FibreFrame& fr);

struct FibreRunReport {
    long long candidates = 0, exposed = 0, accepted = 0;
    long long rootCheckSkipped = 0;  // accepted without the serialized root check (omitted on the GPU)
    double prepMs = 0;               // host: flat arrays, grids
    double waitMs = 0;               // device admission (budget shared with other jobs)
    double uploadMs = 0, residentMs = 0, kernelAMs = 0, selectMs = 0, kernelBMs = 0, compactMs = 0, downloadMs = 0;
    double mergeMs = 0;              // host: records -> FibreSet, report
    double totalMs = 0;              // wall time of generateFibres
    int batchesA = 0, batchesB = 0;
    double slowestLaunchMs = 0;  // the longest single kernel launch
    double inputMB = 0, residentMB = 0, reservedMB = 0, peakMB = 0;
    std::string line() const;        // one log line
};

// The fibre generator on the device: both populations, envelope tubes and explicit plies (parent
// coordinates, density per parent), random overrides, candidate budgets. The serialized root check
// is not applied (FibreRunReport::rootCheckSkipped). `outcomes`: optional per-candidate trace
// (FibreOutcome, as FibreFrame::outcomes). Fibres carry Fibre::candidate.
struct FibreRequest {
    const CurveSet& curves;
    const SurfaceModel& surface;
    const PatternParams& pattern;
    const FibreParams& fibres;
    const FibreFrameGpu& frame;
};
bool generateFibres(const FibreRequest& request, FibreSet& out, FibreRunReport* rep, std::vector<uint8_t>* outcomes, std::string* err);

// The fabric-back test's closest-point search started from a triangle (meshClosestFrom, surface.cuh) next to the
// plain search (meshClosest), for a test: one record per query point and starting triangle.
struct MeshClosestPair {
    int tri[2];       // [0] meshClosest, [1] meshClosestFrom
    float dist[2];
    float bary[2][3];
    Vec3f p[2];
};
bool compareMeshClosest(const TriMesh& mesh, const std::vector<Vec3>& points, const std::vector<int>& startTris, std::vector<MeshClosestPair>& out,
                        std::string* err);

} // namespace knit::gpu
