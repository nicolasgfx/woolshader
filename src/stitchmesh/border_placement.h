#pragma once
// Where the border vertices of a stitch mesh sit on the boundary loop of their knit region.
//
// The border of a stitch mesh is a polygon whose corners lie on the boundary loop. The coverage check
// (stitchmesh/coverage.h) asks that polygon and loop stay within a quarter of a row height of each other. A side that
// spans a bend of the loop cuts it off: on the teddy head one side across a bend stands 0.88 mm off, while the same 25
// corners, one of them on the bend, stay within 0.26 mm. The row layout spaces the corners evenly from the row's first
// vertex and the relaxation slides them by its spring forces; neither looks at the loop's bends.
//
// placeBorderRing chooses the corners for the fit: the same number of corners, in the same order around the loop, the
// first one (the row's first vertex, an end of the yarn) where it is, every side's length within a share of the mean
// side, and the largest distance between a side and its arc of the loop as small as it can be. The search is exact on
// the loop's points taken every `step`: dynamic programming over (corners placed, point reached),
//   best[c + 1][s + k] = min over k of max(best[c][s], distance of the side from point s over k points),
// with the distance of a side the larger of its arc's distance to the chord and the chord's distance to the arc (the
// two directions of the check's Hausdorff distance, side by side).
#include <vector>

#include "model/vec.h"

namespace knit {

struct BorderRingRequest {
    double step = 0.25;                       // spacing of the loop's points the corners may take [mm]
    std::vector<double> spreads = {0.1, 0.2, 0.3};   // side lengths within this share of the mean side, tried in this order
    double limit = 0.85;                      // the first spread whose placement stays within this distance is taken [mm]
    // How far along the loop the first corner may be moved [mm]. 0: it stays (a held end of the yarn). The ring of a
    // relaxation that lets it turn as one has no place of its own within a stitch: half a stitch width lets a corner
    // reach the nearest bend.
    double anchorFreedom = 0;
};

struct BorderRingPlacement {
    bool placed = false;          // false: fewer than three corners, a loop shorter than its corners, or no placement in any spread
    std::vector<Vec3> corners;    // per ring vertex, on the loop (placed only)
    bool inOrder = true;          // the ring as given goes once around the loop (false: two corners have changed places)
    double distanceBefore = 0;    // largest side-to-arc distance of the ring as given, when in order [mm]
    double distance = 0;          // ... of the placement [mm]
    double spread = 0;            // the spread used
    double firstCornerMoved = 0;  // along the loop, in the ring's direction [mm]
    double shortestSide = 0, longestSide = 0;   // of the placement [mm]
};

// loop: a closed polyline (its last point joins the first). ring: the ring's vertices in ring order, on or near the loop;
// ring[0] keeps its place on the loop, or stays within request.anchorFreedom of it.
BorderRingPlacement placeBorderRing(const std::vector<Vec3>& loop, const std::vector<Vec3>& ring, const BorderRingRequest& request);

} // namespace knit
