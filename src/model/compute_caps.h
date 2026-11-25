#pragma once
// What this build and machine compute with: a usable CUDA device, and whether the knitting solver (ncp)
// can run on it; the reasons otherwise. Detected by gpu::detectComputeCaps (gpu/knitting.h),
// read by the settings profiles (computeNotes), the jobs and the viewer. And the build configuration that
// bit-identity is scoped to.
#include <string>

namespace knit {

struct ComputeCaps {
    bool device = false, ncp = false;
    std::string deviceName, whyNoDevice, whyNoNcp;
};

// The build configuration that bit-identity is scoped to: compiler version, floating-point mode, OpenMP,
// NDEBUG and the architecture macros of the build. The knit cache's key and the fibre digests carry it.
std::string buildConfiguration();

} // namespace knit
