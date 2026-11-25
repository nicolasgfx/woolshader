// The build configuration: what bit-identity and the caches' keys are scoped to.
#include "model/compute_caps.h"

namespace knit {

std::string buildConfiguration() {
    std::string o;
#if defined(_MSC_FULL_VER)
    o += "MSVC " + std::to_string(_MSC_FULL_VER);
#elif defined(__clang__)
    o += "clang " __clang_version__;
#elif defined(__GNUC__)
    o += "gcc " __VERSION__;
#else
    o += "unknown-compiler";
#endif
#if defined(_M_FP_FAST)
    o += " fp:fast";
#elif defined(_M_FP_STRICT)
    o += " fp:strict";
#elif defined(_M_FP_PRECISE)
    o += " fp:precise";
#endif
#if defined(_OPENMP)
    o += " openmp " + std::to_string(_OPENMP);
#else
    o += " no-openmp";
#endif
#if defined(NDEBUG)
    o += " release";
#else
    o += " debug";
#endif
#if defined(_M_X64) || defined(__x86_64__)
    o += " x64";
#endif
#if defined(__AVX512F__)
    o += " avx512";
#elif defined(__AVX2__)
    o += " avx2";
#elif defined(__AVX__)
    o += " avx";
#endif
    return o;
}

} // namespace knit
