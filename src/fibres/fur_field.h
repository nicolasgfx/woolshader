#pragma once
// The fur field of one part: one disc per fibre, built from cached roots. Parameters -> discs (cheap,
// deterministic, no mesh queries) -> curves (fur_loops.h growLoops: settle and stuck test).
//
//   1. sampleFurRoots: the candidate roots of one layer, area-weighted on the part's mesh at a capacity
//      density D_cache. Triangle t holds n_t(D) = floor(D A_t + u_t) candidates (u_t = hash(seed, t, 0, 0)),
//      so E[n_t] = D A_t exactly; candidate i sits at the barycentric weights (1 - s, s (1 - u), s u),
//      s = sqrt(hash(seed, t, i, 1)), u = hash(seed, t, i, 2) (the square-root warp: uniform over the
//      triangle). A candidate inside another solid (tested 1e-3 mm above the surface) is dropped there.
//      The cache depends on the mesh, the solids, the capacity, the layer and the seed only.
//   2. buildFurField: a disc per active candidate. At density D <= D_cache the candidates i < n_t(D) of each
//      triangle are active: the same set a sampler at D draws, whatever the capacity, so cached surplus
//      candidates never change the result. A coverage mask keeps a candidate with probability c (its own hash
//      channel). Each kept candidate draws its shape from its own hash channels.
//
// A disc is the control and preview of one fibre: the circle of radius r through the root, centred at
// root + r c, in the plane (d, c). It is an approximate envelope of the curve fur_loops.h grows (the
// curvature ramp, the climb across the disc and the settle move the curve off it), never a constraint.
// A fine fibre is a disc of radius 0: a straight needle along d.
//
// Identity: (part, layer, triangle, candidate index, seed); the seed holds the part and the layer. Every
// curve input stays double precision (only a display buffer may use floats), so a field reproduces the
// prototype's fibres bit for bit (handover rev. 3, section 3).
//
// Hash channels (fur_field.cpp; a channel never changes meaning, that would move every fibre): 0 the count
// offset, 1 / 2 the placement, 3 coverage, 4 / 5 the fine fibres' length, 6 lift, 7 comb turn, 8 the tangent
// when the comb is along the normal, 9 budget key, 10 loop plane, 11 radius, 15 thickness, 16 length share,
// 17 bend. 12-14 belonged to the prototype's helix and tufts and stay unused.
#include <array>
#include <cstdint>
#include <functional>
#include <vector>

#include "fibres/fibre_islands.h"
#include "fibres/fur_look.h"
#include "model/mesh.h"
#include "model/parallel.h"
#include "model/vec.h"

namespace knit {

enum class FurLayer : uint8_t { Loops, Fibres };

// One candidate root (visible: not inside another solid).
struct FurRoot {
    int32_t tri = -1;                 // triangle of the part's mesh
    int32_t index = 0;                // candidate number on that triangle
    Vec3 bary;                        // barycentric weights on the triangle
    Vec3 position, normal;            // surface point, face normal (outward, unit)
    float budgetKey = 0;              // channel 9: the strand budget keeps the smallest keys
};

// The candidates of one layer on one part up to a capacity density, the hidden ones dropped. Immutable once
// built; a field reads it, never copies it.
struct FurRootCache {
    FurLayer layer = FurLayer::Loops;
    uint64_t seed = 0;                // the part's seed of this layer
    double capacity = 0;              // [1/mm^2] the density the candidates were drawn at
    std::vector<double> triangleArea, triangleOffset;  // per triangle: n_t(D) = floor(D area + offset)
    std::vector<FurRoot> roots;       // in triangle order, then candidate number
    int64_t sampled = 0;              // candidates at capacity, hidden ones included
    int64_t hidden = 0;
};

// A root inside another solid (null: none).
using FurHiddenTest = std::function<bool(const Vec3& p)>;
// The coverage in [0, 1] at a point of triangle `tri` with barycentric weights `bary` (null: 1 everywhere).
using FurCoverage = std::function<double(int tri, const Vec3& bary)>;

FurRootCache sampleFurRoots(const TriMesh& m, FurLayer layer, double capacity, uint64_t seed, const FurHiddenTest& hidden,
                            const ExecPolicy& exec);

// One fibre's control: its root, its frame and its shape.
struct FurDisc {
    Vec3 root;                        // surface point (no edit ever moves it)
    Vec3 normal;                      // face normal at the root (outward, unit)
    Vec3 d;                           // start direction (unit)
    Vec3 c;                           // toward the disc centre: unit, across d, c . normal >= 0 (radius > 0)
    double radius = 0;                // loop radius [mm]; 0: a straight fibre
    double turns = 0;                 // the loop's nominal length multiple, length variation folded in
    double rise = 1;                  // the curvature reaches 1 / radius after (pi radius / 2) rise of arc
    double rootRadius = 0, tipRadius = 0;  // [mm]
    double straightLength = 0;        // [mm] a straight fibre's length (radius 0)
    int32_t tri = -1, index = 0;      // the root's identity in its cache
    float budgetKey = 0;
    FurLayer layer = FurLayer::Loops;
};

// One part's discs, both layers: the loops in their cache's order, then the fine fibres.
struct FurField {
    std::vector<FurDisc> discs;
    FurLayerCounts loops, fibres;     // sampled, hidden, masked and kept
};

// The seed of one part's layer: the look's seed mixed with the part's name hash (the loops), and the fine
// fibres' own salt "top_" on top.
uint64_t furPartSeed(uint64_t lookSeed, FurLayer layer, uint64_t partNameHash);

// --- the irregularities (fields, coefficients, composition order)

// A nub: loops near it lean toward its start direction d and take its radius and turns factors.
struct FurClump {
    Vec3 centre, normal, d;
    double radiusFactor = 1, turnsFactor = 1;
};

// The clumps of one part in a dense grid over their bounding box (one more cell on each side): a root looks at
// its own cell and the 26 around it, which holds every clump within reach because a cell is at least the reach.
// Cells are enlarged when the box would need more than kMaxClumpGridCells (a large part, small clumps): still
// correct, a few more clumps per cell. Cell (x, y, z) holds clump[start[i]] .. clump[start[i + 1] - 1],
// i = x + nx (y + ny z), in clump order.
struct FurClumpGrid {
    Vec3 origin;                               // the grid's lowest corner [mm]
    double cell = 0, reach = 0;                // [mm]; cell >= reach
    int nx = 0, ny = 0, nz = 0;
    std::vector<int32_t> start, clump;
};

// The fields of one part that neighbouring fibres share, built once per part from its mesh (scene mm): per mesh
// vertex the patch value f, the island weights and the wear e (a root interpolates them by its barycentric weights),
// and the clumps with their grid for the nearest one. Empty when every irregularity is off.
struct FurPartModifiers {
    bool active = false;                       // any irregularity on (else the fields are not built)
    std::vector<std::array<int, 3>> tri;       // the part's triangles (vertex indices)
    std::vector<double> patch, wear;           // per vertex, in [-1, 1]
    std::vector<IslandWeight> island;          // per vertex: fluffy (slot fuzzy), matted, thin (slot worn)
    std::vector<FurClump> clumps;
    FurClumpGrid clumpGrid;
};
FurPartModifiers buildFurModifiers(const TriMesh& m, const FurLook& look, uint64_t partNameHash);

// The highest local density a layer can reach under the look's modifiers [1/mm^2]: its candidates are active at this
// density, and each is kept with probability D_local / D_bound (hash channel 18). All off: the layer's density.
double furDensityBound(const FurLook& look, FurLayer layer);

// What one part's field is built from.
struct FurFieldRequest {
    const FurRootCache& loops;
    const FurRootCache& fibres;                // an empty cache when the layer is off
    const FurLook& look;
    FurCoverage coverage;                      // null: the whole surface
    const FurPartModifiers* modifiers = nullptr;  // null or inactive: none
    float keyBelow = 2.0f;                     // < 1: only the budget keys below it (the disc preview's draft); the
                                               // counts are then those of that share
};

// The discs of the active candidates of both caches.
FurField buildFurField(const FurFieldRequest& request, const ExecPolicy& exec);

} // namespace knit
