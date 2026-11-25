#pragma once
// Fibre islands: patches of the yarn where the fuzz differs. Each is a sphere with a few satellites
// around a stitch, so its edge is irregular. Weights are stored per control of the supporting strands.
#include <cstdint>
#include <vector>

#include "fibres/fibre_record.h"
#include "model/choices.h"
#include "model/mesh.h"
#include "model/yarn_path.h"

namespace knit {

struct IslandRequest {
    const CurveSet& strands;       // the yarn or its plies
    const KnitTopology& stitches;  // of the yarn
    const FibreIslands& islands;
    double stitchArea;             // [mm2] turns the rate per cm2 into a chance per stitch
    uint64_t seed;
};

// Per strand, per control.
std::vector<std::vector<IslandWeight>> islandWeights(const IslandRequest& request);

// Interpolated between a strand's controls. gpu/fibres/fibre_candidate.cuh holds the twin.
IslandWeight islandAt(const std::vector<IslandWeight>& strand, double u);

// --- the placement both the wool islands and the Fur workflow's islands use

// One sphere of an island, and what the fibres become at its centre (the three weight slots: wool fuzzy, matted,
// worn; Fur fluffy, matted, thin).
struct IslandBlob {
    Vec3 centre;
    double radius = 0;
    IslandWeight peak;
};
// An island at `centre` with its tangent frame (U, V): a sphere of `radius` and three satellites around it (their
// directions and sizes drawn from hash channels 3-8 of (seed, a, b)), so its edge is irregular.
void addIsland(std::vector<IslandBlob>& blobs, const Vec3& centre, const Vec3& U, const Vec3& V, double radius, IslandWeight peak, uint64_t seed,
               int a, int b);
// The strongest of each kind at p: full at a sphere's centre, fading smoothly (smoothstep) to nothing at its edge.
IslandWeight islandWeightAt(const std::vector<IslandBlob>& blobs, const Vec3& p);

// The Fur workflow's islands on a mesh surface: centres area-weighted at `ratePerCm2` (n_t = floor(rate A_t / 100
// + u_t), the root sampler's rule), each of diameter `size` x (0.7 .. 1.3), its kind drawn by the shares (slot fuzzy:
// fluffy, matted, worn: thin) and its peak `strength` (0..1). All-zero shares: no islands.
struct SurfaceIslands {
    double ratePerCm2 = 0, size = 15, strength = 0.6;
    double fluffy = 1, matted = 1, thin = 1;
};
std::vector<IslandBlob> surfaceIslandBlobs(const TriMesh& m, const SurfaceIslands& islands, uint64_t seed);

} // namespace knit
