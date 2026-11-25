#pragma once
// Images: the writers of the shading maps and base colours (PNG through stb_image_write, EXR by hand), their
// readers, decoding of any base-colour file through stb_image, and reading a colour at a uv.
//   PNG  RGB8, deflate-compressed by stb_image_write (a yarn's baked maps are hundreds of MB uncompressed);
//   EXR  scanline, NO_COMPRESSION, channels B G R as 32-bit FLOAT.
// Pixels are row-major with row 0 at the TOP (the photon_path_tracer loader reads top-row-first and
// its sampler flips v).
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "model/vec.h"

namespace knit {

// One decoded base-colour image: sRGB bytes, rows top first, as stored in the file.
struct ColourImage {
    std::string path;
    int width = 0, height = 0;
    std::vector<uint8_t> rgba;
};

bool writePngRgb8(const std::string& path, int width, int height, const std::vector<uint8_t>& rgb);
bool writeExrRgbF32(const std::string& path, int width, int height, const std::vector<float>& rgb);
bool readPngRgb8(const std::string& path, int& width, int& height, std::vector<uint8_t>& rgb);
// Reads exactly the files above: false for anything else.
bool readExrRgbF32(const std::string& path, int& width, int& height, std::vector<float>& rgb);

// PNG, JPEG, TGA or BMP. The size first, so a caller can refuse before decoding; `why` carries the
// decoder's reason. An embedded image is refused beyond maxSide on a side.
bool readImageSize(const std::string& path, int& width, int& height, std::string* why);
bool decodeImageFile(const std::string& path, ColourImage& out, std::string* why);
// One grey byte per texel, rows top first; colours are converted by the decoder.
bool decodeGreyImageFile(const std::string& path, int& width, int& height, std::vector<uint8_t>& value, std::string* why);
bool decodeImageBytes(const unsigned char* bytes, size_t size, int maxSide, ColourImage& out, std::string* why);
// At most maxSide on a side, by a box filter over a power of two.
ColourImage shrunkToSide(const ColourImage& image, int maxSide);

// The nearest texel (repeat wrap) at uv in the OBJ convention (v up; row 0 is the image's top), as
// 0xRRGGBB of its sRGB bytes.
uint32_t nearestTexel(const ColourImage& t, Vec2 uvVUp);
// The bilinear colour at uv (same convention; texel centres at (i + 1/2) / size), interpolated in
// linear light.
Vec3 bilinearTexel(const ColourImage& t, Vec2 uvVUp);
// The mean linear colour at a triangle's interior points, by nearest texel.
Vec3 triangleTextureColour(const ColourImage& t, const std::array<Vec2, 3>& uv);
// One sRGB byte -> linear value, and back.
double srgbByteToLinear(uint32_t byte);
uint32_t linearToSrgbByte(double linear);

} // namespace knit
