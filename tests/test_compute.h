#pragma once
// What the scene / pipeline tests compute with: scenes compute on the GPU (the tool has no CPU
// scene path); the direct knitMesh calls fall back to the CPU test double without a device.
// Tests whose subject is a CPU implementation the tool keeps (exposed placement, the reference
// preview, the flat swatch's gate (a)) pin the CPU explicitly (CMakeLists.txt, labels
// retained_cpu and cpu_reference).
#include <stdexcept>
#include <memory>
#include <string>

#include "cpu_numerics.h"
#include "gpu/device.h"
#include "gpu/knitting.h"

namespace knit_test {

// A usable device (the GPU halves of the direct-call tests run only then).
inline bool useGpu() {
    static const bool on = knit::gpu::available();
    return on;
}

// The numerics of a direct knitMesh call: the device's when usable, else the CPU double (one
// device object per thread, kept for the process).
inline knit::KnitNumerics& testNumerics() {
    if (!useGpu()) return cpuNumerics();
    thread_local std::unique_ptr<knit::KnitNumerics> num;
    if (!num) {
        std::string e;
        num = knit::gpu::makeKnitNumerics(&e);
        if (!num) throw std::runtime_error("GPU knit numerics: " + e);
    }
    return *num;
}

} // namespace knit_test
