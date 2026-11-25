#pragma once
// The RIM remesher (knit_rim, third_party/rim) behind a plain-array interface, C++14: stage 1 of the stitch-mesh
// experiment's Method B calls it from src/stitchmesh/remesh_stage.cpp and nothing else (CHANGES.md: upstream, licences).
// The pipeline of RIM's batch mode on a surface: hierarchy (RIM's own normals, boundary handling and seeded random
// start), orientation sweeps then position sweeps from level L-2 down to 0 (RIM's Optimizer loop without its thread:
// `sweepsPerLevel` per level instead of its fixed 100), then meshExtraction2D with the flags its batch mode sets.
// Wu 2018 section 4 changes (knit_rim_hooks.h): 2-RoSy orientation compatibility and the post-processing inside the
// extraction loop. Not thread-safe (RIM's extraction keeps file-scope state): the caller serializes.
#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace knit_rim {

struct Options {
    double spacing = 0;            // target edge length [mm] along q (RIM's mScale)
    // The lattice's t = n x q spacing relative to q's (rim.cmake): h / w. 1 here (M5 qualifies the rectangular
    // lattice); it needs the intrinsic 2-RoSy field.
    double aspect = 1;
    int rosy = 2;                  // 2 (Wu 2018 section 4) or 4 (RIM)
    bool extrinsic = false;        // Wu 2018's models use the intrinsic field (section 9)
    int sweepsPerLevel = 20;
    bool randomization = true;     // RIM's seeded neighbour shuffle (its iteration counters and the row index)
    uint64_t seed = 0;             // build()'s random start: the frames and lattice origins (knit_rim_seed, rim.cmake)
    bool wuLoop = true;            // Wu 2018 section 4 post-processing inside the extraction loop (knit_wu_loop)
    // Per input vertex, a tangent direction fixing the orientation field (empty: solved). Coarser levels
    // take their first child's direction projected into their tangent plane; no orientation sweeps.
    std::vector<std::array<double, 3>> orientation;
    // Input subdivision before the hierarchy: when the longest edge exceeds half the spacing or twice the mean
    // edge, edges are split down to min(spacing / 2, 2 x mean). predictedVertices() tells the count beforehand.
    bool subdivide = true;
    // A time cap [ms] on the whole call (0: none), checked after the hierarchy, after every sweep and before the
    // extraction; the extraction itself runs to its end. A hit returns false with Result::deadlineHit set.
    double deadlineMs = 0;
    // The orientation field only (the stockinette experiment's checkpoint 1: input analysis before any mesh): the
    // hierarchy and the orientation sweeps (or the fixed orientation), then fieldQ and the orientation singularities;
    // no position sweeps and no extraction (pos, faces, fieldO stay empty).
    bool fieldsOnly = false;
};

struct Result {
    std::vector<std::array<double, 3>> pos, nrm;   // extracted vertices (lattice positions), unit normals
    std::vector<std::array<double, 3>> axis;       // per extracted vertex: the field q of the nearest finest vertex, in its tangent plane
    std::vector<std::vector<int>> faces;           // polygons as extracted (3 or more sides)
    std::vector<int> levelSizes;                   // vertices per hierarchy level (finest first)
    std::vector<std::array<double, 3>> fieldQ, fieldO;  // the finest level's solved fields, per vertex after subdivision
    std::vector<std::array<double, 3>> orientationSingularPoints, positionSingularPoints;  // on the input surface
    int inputVertices = 0;                         // the finest level: the input after subdivision
    int wuRounds = 0, wuDiagonals = 0, wuLoopSplits = 0;  // Wu section 4 operations inside the extraction loop
    int genusIn = 0, genusOut = 0;
    double msHierarchy = 0, msOrientations = 0, msPositions = 0, msExtraction = 0;  // stage times
    bool deadlineHit = false;                      // stopped by Options::deadlineMs; `log` ends with the stage
    std::string log;                               // RIM's console output
};

// How many vertices the input has after RIM's subdivision rule at this spacing (the finest level): the rule run
// on a copy, before anything of the hierarchy is built. The input's own count when no edge needs a split.
int predictedVertices(const std::vector<std::array<double, 3>>& pos, const std::vector<std::array<int, 3>>& tri, double spacing);

bool remesh(const std::vector<std::array<double, 3>>& pos, const std::vector<std::array<int, 3>>& tri, const Options& o, Result& out,
            std::string* err);

} // namespace knit_rim
