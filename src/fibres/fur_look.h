#pragma once
// The Fur workflow's look in the generator's own units (mm, fibres per mm^2, degrees where named): what
// fur_field.h turns into one disc per fibre and fur_loops.h into curves. scene/stage_params.cpp maps the
// pile.* settings onto it (the settings' cm^2 and m/s^2 become mm^2 and a gravity scale there).
//
// The look has two layers on every part, each with its own roots (fur_field.h FurRootCache):
//   loops:  thick fibres, each curling round once in its own plane (a disc through its root, turned about
//           its start direction at random), so the pile reads as small arches; the approved look of the
//           furry teddy (frozen baseline gate/top15, 2026-11-19);
//   fibres: fine straight fibres standing over the loops (the silhouette's halo).
// Variation is drawn per fibre, independently (uniform shares, fur_field.cpp lists the hash channels).
#include <cstdint>

#include "fibres/fur.h"

namespace knit {

struct FurLoops {
    double density = 4;               // roots per mm^2
    double turns = 1;                 // the loop's nominal length, in multiples of one loop (fur_loops.h)
    double straightLength = 4;        // [mm] each fibre's length at curl 0 (no loop)
    double rootDiameter = 0.6;        // [mm]
    double taper = 0.8;               // tip diameter / root diameter
    double bend = 1;                  // the curvature reaches 1 / curl after this many quarter turns' arc
    double curl = 0.9;                // loop radius [mm]; 0: straight fibres
};

// Shares in [0, 1] (densityScale in mm). Each fibre: length x (1 - length u); diameter x (1 +- thickness);
// rise x (1 +- bend); radius x (1 +- curl / 2); lift 90 - 50 v +- 40 v deg and comb turn +- 180 v deg for
// direction v. density / densityScale: spatial patches (not used before the irregularities milestone).
struct FurVariation {
    double density = 0, densityScale = 10;
    double length = 0.5, thickness = 0.33, direction = 1, bend = 0.2, curl = 1;
};

// The fine straight fibres: length log-normal about the median, clamped to [median / 2, 2 median]; lift
// liftDeg +- 40 v deg and comb turn +- 180 v deg for direction v; no thickness variation.
struct FurFibres {
    bool enabled = true;
    double density = 1.5;             // roots per mm^2
    double lengthMedian = 4;          // [mm]
    double lengthSpread = 0.3;        // sigma of ln(length)
    double rootDiameter = 0.06;       // [mm]
    double taper = 0.3;
    double liftDeg = 85;
    double direction = 1;
};

// The irregularities: fields neighbouring fibres share.
// The patches' feature size is FurVariation::densityScale. All off at 0.
struct FurIrregularity {
    double clumpSize = 2, clumpStrength = 0;  // [mm] reach of a nub, [0, 1]
    double patchLength = 0, patchLift = 0;    // [0, 1]: the patch field on the loops' turns and lift
    double wear = 0;                          // [0, 1]: exposed, convex places flatter, shorter, sparser
};

// Islands of matted, fluffy or thin pile (fibres/fibre_islands.h places them, as the wool islands).
struct FurIslands {
    double rate = 0;                          // per cm^2
    double size = 15;                         // [mm] diameter
    double strength = 0.6;                    // [0, 1]
    double matted = 1, fluffy = 1, thin = 1;  // relative shares of the kinds (all 0: none)
};

// The settle every fibre of both layers takes (fur.h settleHair).
struct FurSettle {
    StrandSettleParams strand;
    double gravityFullLength = 30;    // [mm] (0: every fibre the full gravity)
    double collisionOffset = 0.02;    // [mm] beyond each vertex's radius
};

struct FurLook {
    FurLoops loops;
    FurVariation variation;
    FurFibres fibres;
    FurIrregularity irregular;
    FurIslands islands;
    FurSettle settle;
    uint64_t seed = 0x6675725fULL;    // "fur_" (mixed with pile.seed and each part's name)
};

// The fibres of one layer through the generator (fur_field.h, fur_loops.h): candidates on the mesh, and
// where each one left. Every count stays separate per layer.
struct FurLayerCounts {
    int64_t sampled = 0;              // candidate roots at the layer's density
    int64_t hidden = 0;               // inside another solid (dropped)
    int64_t masked = 0;               // outside the coverage mask (dropped)
    int64_t thinned = 0;              // left out by the local density (patches, islands, wear)
    int64_t kept = 0;                 // discs: visible candidates before the settle
    int64_t clamped = 0;              // loops whose length the minimum-length rule raised (fur_loops.h)
    int64_t overBudget = 0;           // dropped by the segment budget
    int64_t stuck = 0;                // rejected after the settle (fur.h hairInsideBody)
    int64_t grown = 0;                // exported
    int64_t segments = 0;             // visible segments of the grown fibres
};

struct FurReport {
    FurLayerCounts loops, fibres;
    double maxResidual = 0;           // largest last-step move over the grown fibres [mm]
    double msSample = 0, msField = 0, msSettle = 0;
};

} // namespace knit
