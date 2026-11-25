#pragma once
// Colour conventions every module shares: the sRGB transfer (IEC 61966-2-1) and the default yarn colour.
// Materials and averages work in linear light; image files and the display are sRGB encoded.
#include <cmath>

#include "model/vec.h"

namespace knit {

inline double srgbEncode(double linear) {
    return linear <= 0.0031308 ? 12.92 * linear : 1.055 * std::pow(linear, 1.0 / 2.4) - 0.055;
}
inline double srgbDecode(double encoded) {
    return encoded <= 0.04045 ? encoded / 12.92 : std::pow((encoded + 0.055) / 1.055, 2.4);
}

// Linear. The colour of undyed wool, and of a surface that has no colour of its own.
inline constexpr Vec3 kDefaultYarnColour{0.80, 0.78, 0.72};

} // namespace knit
