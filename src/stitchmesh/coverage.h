#pragma once
// Coverage of the reference by a typed stitch mesh, reported apart from cell
// quality: a missing ring or a shortened patch shows up here even when every cell is good.
//   covered     every reference sample (vertices and triangle centroids, each weighted by its share of the area) has an
//               output face within coverDistance; the uncovered area fraction, outside a band of boundaryBand at the
//               reference boundaries, must stay within uncoveredLimit
//   twice       the segment through the sample along the reference normal, coverDistance each way, crosses the output
//               surface at most once (a fold or an overlap crosses it twice)
//   area        the output area projected onto the reference normal / the reference area, within [areaLow, areaHigh]
//   boundaries  every output boundary loop matches one reference loop, Hausdorff distance within hausdorffLimit
// Output faces are ear-clipped in their own plane (earClipFace) for the distance and ray tests.
#include <string>
#include <vector>

#include "model/mesh.h"
#include "model/status.h"
#include "stitchmesh/polygon_mesh.h"

namespace knit {

struct CoverageOptions {
    double coverDistance = 0.55;     // 0.1 w
    double boundaryBand = 2.75;      // 0.5 h
    double hausdorffLimit = 1.375;   // 0.25 h
    double uncoveredLimit = 0.005;
    double areaLow = 0.98, areaHigh = 1.02;
};

struct CoverageSample {
    Vec3 p;
    double weight = 0;     // area [mm^2]
    bool vertex = false;   // a reference vertex (else a triangle centroid)
    bool excluded = false; // inside the boundary band
    double distance = 0;   // to the output surface
    int face = -1;         // the covering output face (-1: none within reach)
    int crossings = 0;     // of the normal segment with the output surface
};

struct LoopCorrespondence {
    int outputLoop = -1, referenceLoop = -1;
    int outputVertices = 0;
    double outputToReference = 0, referenceToOutput = 0, hausdorff = 0;
};

struct CoverageReport {
    Status status = Status::Unresolved;
    int samples = 0, excluded = 0, uncovered = 0, coveredTwice = 0;
    double referenceArea = 0, consideredArea = 0, uncoveredArea = 0, uncoveredFraction = 0;
    // false when every sample lies in the boundary band (a strip narrower than one row): the uncovered share is not
    // measurable then; uncoveredFraction is 0 and the status at best Unresolved, never Pass
    bool uncoveredMeasured = false;
    double outputArea = 0, projectedArea = 0, projectedRatio = 0;
    int outputLoops = 0, referenceLoops = 0;
    std::vector<LoopCorrespondence> loops;
    double worstHausdorff = 0;
    std::vector<CoverageSample> perSample;
    std::string detail;
    std::string text() const;
};

// refLoops: the reference boundary loops as vertex ids in walking order.
CoverageReport measureCoverage(const TriMesh& reference, const std::vector<std::vector<int>>& refLoops, const PolygonMesh& out, const CoverageOptions& o);

// The border loops of a polygon mesh, each as vertex ids in walking order.
std::vector<std::vector<int>> polygonBorderLoops(const PolygonMesh& m);

} // namespace knit
