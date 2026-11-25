#pragma once
// Coverage of a knitted surface: how much of the base mesh the rendered yarn covers, and where the
// holes are. The overlap mode's gap patches go into the holes (scene/knit/gap_patches.cpp), and
// knit_check --flat-dmax measures flat swatches with it.
//
// Domain (fixed): area-uniform stratified samples on the component's base mesh. The samples never
// depend on the knit, so dropped stitches count as uncovered.
//
// A sample p (triangle t, normal n) is covered iff a rendered tube triangle of ELIGIBLE yarn meets
// the segment p +- n H. A tube triangle is eligible iff one of its vertices' foot triangles f (the
// base triangle that yarn point is bound to, SurfaceBinding) satisfies
//   - geodesic distance t -> f <= G (bounded search over the triangle adjacency, centroid edge
//     lengths: an upper bound of the surface distance between the centroids), and
//   - n_f . n > minNormalCos (face normals),
// so yarn of the opposite sheet of a thin ear and of a nearby fold of the same surface never
// counts. Without feet (a flat swatch) every hit counts.
//
// Regions: uncovered samples joined to neighbours within 1.6 spacing in 3D whose normals agree
// (n . n' > 0.3). Width of a region = 2 max_s d(s) - spacing, d = graph distance from the region's
// samples to the nearest sample outside it (a strip k samples wide reads about k spacings); it is an
// estimate at the sample spacing, reported with the spacing.
#include <cstdint>
#include <string>
#include <vector>

#include "knitting/layout/base_surface.h"
#include "knitting/tube/yarn_surface.h"
#include "model/mesh.h"
#include "model/parallel.h"
#include "model/yarn_path.h"

namespace knit {

struct CoverageParams {
    double spacing = 0.25;           // target sample spacing [mm]
    long long maxSamples = 400000;   // the spacing grows to keep the count below this
    double reach = 2.0;              // H: half length of the normal segment [mm]
    double ownershipRadius = 11.0;   // G [mm]
    double minNormalCos = 0.0;
    uint64_t seed = 0x5eed5eedull;
    double neighbourFactor = 1.6;    // region adjacency radius in spacings
};

// The component's rendered yarn: its triangles and, per vertex, the base triangle it is bound to
// (-1: unbound, never eligible).
struct CoverageTube {
    const TubeMesh* tube = nullptr;
    std::vector<int> vertexFoot;     // empty: every hit counts (a flat swatch)
};

struct CoverageInput {
    const TriMesh* mesh = nullptr;   // the component's base mesh (normals computed)
    CoverageTube own;
    // The measured region: the samples are drawn on the whole mesh (their density does not depend on
    // it), then those outside this box are dropped. Default: everything. Tests measure a strip.
    Vec3 regionLo{-1e300, -1e300, -1e300}, regionHi{1e300, 1e300, 1e300};
};

// One connected region: its area and width as above, its bounding box and its centre (the sample
// farthest inside it: the largest distance to an unflagged sample; seeds the gap patches).
struct CoverageRegion {
    double area = 0, width = 0;
    Vec3 lo, hi;
    int centreTri = -1;
    double centreBary[3] = {0, 0, 0};
    Vec3 centre;
};

struct CoverageRegions {
    int count = 0;
    double maxArea = 0, maxWidth = 0;        // [mm^2], [mm]
    double totalArea = 0;
    std::vector<CoverageRegion> items;       // per region (largest first)
};

struct CoverageReport {
    long long samples = 0;
    double spacing = 0, area = 0, sampleArea = 0;  // spacing used, mesh area, area per sample
    long long uncovered = 0, unownedHits = 0;      // unownedHits: samples hit only by ineligible yarn
    double uncoveredFrac = 0;
    CoverageRegions holes;
    double ms = 0;
    std::string summary() const;
};

CoverageReport measureCoverage(const CoverageInput& in, const CoverageParams& p, const ExecPolicy& exec = {});

// Per tube vertex, the base triangle of the nearest control's binding state: tube vertex arc
// length s -> spline parameter -> the control with the nearest Greville abscissa (uniform cubic
// B-spline: control c at u = c - 1, clamped). states[y][c].tri as bound by placement / settle;
// -1 where a yarn has no states.
std::vector<int> tubeVertexFeet(const TubeMesh& tube, const CurveSet& cs, const std::vector<std::vector<ProjState>>& states);

// --- a flat swatch (knit_check --flat-dmax)

// The window spanned by the interior stitch centres (wale 1..n-2, course 1..R-2), widened by half the
// mean spacing.
struct InteriorWindow { double x0 = 0, y0 = 0, x1 = 0, y1 = 0; };
InteriorWindow interiorWindow(const CurveSet& cs);

// Fraction of the window not covered by the front projection of the tube mesh's triangles.
double meshOpenArea(const TubeMesh& m, const InteriorWindow& w, int resolution = 500);

} // namespace knit
