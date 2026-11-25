#pragma once
// Gate (a) on the device: the twins of knitting/checks (checkCurvature, samplePath + checkClearance,
// checkLocalSignature) in batches over one uploaded curve set, in single precision:
//   curvature.cu   checkCurvature per span
//   clearance.cu   samplePath and checkClearance
//   signature.cu   checkLocalSignature, CLOSE drawstring passes
// Results are GPU diagnostic verdicts (float32 evaluation of the CPU algorithms with float allowances),
// never the CPU's certified PASS / FAIL / UNRESOLVED. The port (port/knit_numerics.cpp) composes them.
#include <vector>

#include "gpu/device/runtime.h"
#include "gpu/geometry/device_geometry.h"
#include "knitting/checks/checks.h"

namespace knit::gpu {

// The sampled polylines on the device. Polyline y occupies global vertices [vertexOff[y],
// vertexOff[y + 1]); segment s of polyline y (vertices s, s + 1) is stored at global index
// vertexOff[y] + s in the segment arrays; cumLo and cumHi hold segments + 1 entries per polyline
// from vertexOff[y].
struct DevicePolylines {
    DeviceBuffer p;             // float, 3 per vertex
    DeviceBuffer u, s;          // double per vertex
    DeviceBuffer segLo, segHi;  // double per segment slot
    DeviceBuffer segEps;        // float per segment slot
    DeviceBuffer cumLo, cumHi;  // double
    std::vector<int> vertexOff;
};

std::vector<CurvatureResult> checkCurvatureOnDevice(const UploadedCurves& d, double radius, const Allowances& al);
DevicePolylines samplePolylinesOnDevice(const UploadedCurves& d, double eps, int maxPerSpan, std::vector<Polyline>& polylines);
struct ClearanceInput {
    const UploadedCurves& curves;
    const DevicePolylines& devicePolylines;
    const std::vector<Polyline>& polylines;
    const std::vector<CurvatureResult>& curvature;
    const PatternParams& pattern;
    const Allowances& allowances;
};
ClearanceResult checkClearanceOnDevice(const ClearanceInput& in);
LocalSignatureResult checkSignatureOnDevice(const UploadedCurves& d);

} // namespace knit::gpu
