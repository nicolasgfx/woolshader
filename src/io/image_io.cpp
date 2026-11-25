// Images: PNG and EXR writers, base-colour decoding through stb_image, colours at a uv (image_io.h).
#include "io/image_io.h"

#include <algorithm>
#include <cmath>
#include <cstring>

#include "io/files.h"
#include "model/colour.h"
#include "stb_image.h"
#include "stb_image_write.h"

namespace knit {

namespace {

void le32(std::vector<uint8_t>& o, uint32_t v) {
    for (int s = 0; s < 32; s += 8) o.push_back(uint8_t(v >> s));
}
void le64(std::vector<uint8_t>& o, uint64_t v) {
    for (int s = 0; s < 64; s += 8) o.push_back(uint8_t(v >> s));
}
uint32_t rd32le(const uint8_t* p) { return p[0] | (uint32_t(p[1]) << 8) | (uint32_t(p[2]) << 16) | (uint32_t(p[3]) << 24); }

bool writeBytes(const std::string& path, const std::vector<uint8_t>& b) {
    return replaceFile(path, std::string(b.begin(), b.end()), nullptr);
}

// stb_image_write's sink: the PNG's bytes into the vector passed as the context, so the file goes through
// writeBytes like every other image (stb's own file writing opens narrow paths).
void appendPngBytes(void* context, void* data, int size) {
    std::vector<uint8_t>* out = static_cast<std::vector<uint8_t>*>(context);
    const uint8_t* bytes = static_cast<const uint8_t*>(data);
    out->insert(out->end(), bytes, bytes + size);
}

void attr(std::vector<uint8_t>& o, const char* name, const char* type, const std::vector<uint8_t>& v) {
    o.insert(o.end(), name, name + std::strlen(name) + 1);
    o.insert(o.end(), type, type + std::strlen(type) + 1);
    le32(o, uint32_t(v.size()));
    o.insert(o.end(), v.begin(), v.end());
}

} // namespace

bool writePngRgb8(const std::string& path, int w, int h, const std::vector<uint8_t>& rgb) {
    if (w <= 0 || h <= 0 || rgb.size() != size_t(w) * h * 3) return false;
    std::vector<uint8_t> o;
    o.reserve(rgb.size() / 3);
    if (!stbi_write_png_to_func(appendPngBytes, &o, w, h, 3, rgb.data(), 3 * w)) return false;
    return writeBytes(path, o);
}

bool readPngRgb8(const std::string& path, int& w, int& h, std::vector<uint8_t>& rgb) {
    int channels = 0;
    unsigned char* px = stbi_load(path.c_str(), &w, &h, &channels, 3);
    if (!px) return false;
    rgb.assign(px, px + size_t(w) * h * 3);
    stbi_image_free(px);
    return true;
}

bool writeExrRgbF32(const std::string& path, int w, int h, const std::vector<float>& rgb) {
    if (w <= 0 || h <= 0 || rgb.size() != size_t(w) * h * 3) return false;
    std::vector<uint8_t> o = {0x76, 0x2f, 0x31, 0x01, 2, 0, 0, 0};
    std::vector<uint8_t> ch;
    for (const char* c : {"B", "G", "R"}) {  // alphabetical
        ch.insert(ch.end(), c, c + 2);
        le32(ch, 2);                            // FLOAT
        ch.insert(ch.end(), {0, 0, 0, 0});      // pLinear, reserved
        le32(ch, 1);
        le32(ch, 1);
    }
    ch.push_back(0);
    attr(o, "channels", "chlist", ch);
    attr(o, "compression", "compression", {0});
    std::vector<uint8_t> box;
    for (int v : {0, 0, w - 1, h - 1}) le32(box, uint32_t(v));
    attr(o, "dataWindow", "box2i", box);
    attr(o, "displayWindow", "box2i", box);
    attr(o, "lineOrder", "lineOrder", {0});
    std::vector<uint8_t> one, v2;
    float f1 = 1.0f, f0 = 0.0f;
    uint32_t u1, u0;
    std::memcpy(&u1, &f1, 4);
    std::memcpy(&u0, &f0, 4);
    le32(one, u1);
    le32(v2, u0);
    le32(v2, u0);
    attr(o, "pixelAspectRatio", "float", one);
    attr(o, "screenWindowCenter", "v2f", v2);
    attr(o, "screenWindowWidth", "float", one);
    o.push_back(0);
    const size_t lineBytes = size_t(w) * 3 * 4, table = o.size();
    for (int y = 0; y < h; ++y) le64(o, uint64_t(table + size_t(h) * 8 + size_t(y) * (8 + lineBytes)));
    for (int y = 0; y < h; ++y) {
        le32(o, uint32_t(y));
        le32(o, uint32_t(lineBytes));
        for (int c : {2, 1, 0})  // B, G, R
            for (int x = 0; x < w; ++x) {
                uint32_t u;
                std::memcpy(&u, &rgb[(size_t(y) * w + x) * 3 + c], 4);
                le32(o, u);
            }
    }
    return writeBytes(path, o);
}

bool readExrRgbF32(const std::string& path, int& w, int& h, std::vector<float>& rgb) {
    std::string text;
    if (!readWholeFile(path, text, nullptr)) return false;
    const std::vector<uint8_t> b(text.begin(), text.end());
    if (b.size() < 8 || rd32le(b.data()) != 20000630u) return false;
    size_t p = 8;
    int x0 = 0, y0 = 0, x1 = -1, y1 = -1;
    while (p < b.size() && b[p] != 0) {
        const std::string name(reinterpret_cast<const char*>(&b[p]));
        p += name.size() + 1;
        const std::string type(reinterpret_cast<const char*>(&b[p]));
        p += type.size() + 1;
        const uint32_t n = rd32le(&b[p]);
        p += 4;
        if (name == "dataWindow") {
            x0 = int(rd32le(&b[p]));
            y0 = int(rd32le(&b[p + 4]));
            x1 = int(rd32le(&b[p + 8]));
            y1 = int(rd32le(&b[p + 12]));
        }
        if (name == "compression" && b[p] != 0) return false;
        p += n;
    }
    w = x1 - x0 + 1;
    h = y1 - y0 + 1;
    if (w <= 0 || h <= 0) return false;
    p += 1 + size_t(h) * 8;
    rgb.assign(size_t(w) * h * 3, 0.f);
    for (int y = 0; y < h; ++y) {
        if (p + 8 + size_t(w) * 12 > b.size()) return false;
        const int yy = int(rd32le(&b[p])) - y0;
        p += 8;
        for (int c : {2, 1, 0})
            for (int x = 0; x < w; ++x, p += 4) {
                const uint32_t u = rd32le(&b[p]);
                std::memcpy(&rgb[(size_t(yy) * w + x) * 3 + c], &u, 4);
            }
    }
    return true;
}

bool readImageSize(const std::string& path, int& width, int& height, std::string* why) {
    int channels = 0;
    if (stbi_info(path.c_str(), &width, &height, &channels)) return true;
    return fail(why, stbi_failure_reason());
}

namespace {

// The stb result as an image, freed either way.
bool takePixels(unsigned char* px, int width, int height, ColourImage& out, std::string* why) {
    if (!px) return fail(why, stbi_failure_reason());
    out.width = width;
    out.height = height;
    out.rgba.assign(px, px + size_t(width) * height * 4);
    stbi_image_free(px);
    return true;
}

} // namespace

bool decodeImageFile(const std::string& path, ColourImage& out, std::string* why) {
    int w = 0, h = 0, channels = 0;
    unsigned char* px = stbi_load(path.c_str(), &w, &h, &channels, 4);
    out.path = path;
    return takePixels(px, w, h, out, why);
}

bool decodeGreyImageFile(const std::string& path, int& width, int& height, std::vector<uint8_t>& value, std::string* why) {
    int channels = 0;
    unsigned char* px = stbi_load(path.c_str(), &width, &height, &channels, 1);
    if (!px) return fail(why, stbi_failure_reason() ? stbi_failure_reason() : "?");
    value.assign(px, px + size_t(width) * size_t(height));
    stbi_image_free(px);
    return true;
}

bool decodeImageBytes(const unsigned char* bytes, size_t size, int maxSide, ColourImage& out, std::string* why) {
    int w = 0, h = 0, channels = 0;
    if (!stbi_info_from_memory(bytes, static_cast<int>(size), &w, &h, &channels)) return fail(why, stbi_failure_reason());
    if (w > maxSide || h > maxSide) return fail(why, "larger than the texture limit");
    unsigned char* px = stbi_load_from_memory(bytes, static_cast<int>(size), &w, &h, &channels, 4);
    return takePixels(px, w, h, out, why);
}

ColourImage shrunkToSide(const ColourImage& image, int maxSide) {
    int f = 1;
    while (image.width / f > maxSide || image.height / f > maxSide) f *= 2;
    if (f == 1) return image;
    ColourImage t;
    t.path = image.path;
    t.width = std::max(1, image.width / f);
    t.height = std::max(1, image.height / f);
    t.rgba.resize(size_t(t.width) * t.height * 4);
    for (int y = 0; y < t.height; ++y)
        for (int x = 0; x < t.width; ++x)
            for (int c = 0; c < 4; ++c) {
                unsigned sum = 0;
                for (int j = 0; j < f; ++j)
                    for (int i = 0; i < f; ++i) sum += image.rgba[(size_t(y * f + j) * image.width + size_t(x * f + i)) * 4 + c];
                t.rgba[(size_t(y) * t.width + x) * 4 + c] = uint8_t(sum / unsigned(f * f));
            }
    return t;
}

double srgbByteToLinear(uint32_t b) { return srgbDecode(b / 255.0); }

namespace {

// A table: a bilinear lookup converts 12 bytes, and every stitch makes 16 lookups.
const std::array<double, 256>& srgbByteToLinearTable() {
    static const std::array<double, 256> table = [] {
        std::array<double, 256> t{};
        for (uint32_t b = 0; b < 256; ++b) t[b] = srgbByteToLinear(b);
        return t;
    }();
    return table;
}

} // namespace

uint32_t linearToSrgbByte(double c) {
    return static_cast<uint32_t>(std::lround(srgbEncode(std::clamp(c, 0.0, 1.0)) * 255));
}

uint32_t nearestTexel(const ColourImage& t, Vec2 uv) {
    // Repeat wrap, v up (row 0 is the top of the image), nearest texel: as the importer's sample().
    const double u = uv.x - std::floor(uv.x), v = uv.y - std::floor(uv.y);
    const int x = std::clamp(int(u * t.width), 0, t.width - 1), y = std::clamp(int((1.0 - v) * t.height), 0, t.height - 1);
    const uint8_t* p = &t.rgba[(size_t(y) * t.width + x) * 4];
    return uint32_t(p[0]) << 16 | uint32_t(p[1]) << 8 | uint32_t(p[2]);
}

Vec3 bilinearTexel(const ColourImage& t, Vec2 uv) {
    const std::array<double, 256>& linear = srgbByteToLinearTable();
    const double fromLeft = (uv.x - std::floor(uv.x)) * t.width - 0.5, fromTop = (std::ceil(uv.y) - uv.y) * t.height - 0.5;
    const int x0 = static_cast<int>(std::floor(fromLeft)), y0 = static_cast<int>(std::floor(fromTop));
    const double fx = fromLeft - x0, fy = fromTop - y0;
    auto texel = [&](int x, int y) {
        const uint8_t* p = &t.rgba[(size_t((y + t.height) % t.height) * t.width + (x + t.width) % t.width) * 4];
        return Vec3{linear[p[0]], linear[p[1]], linear[p[2]]};
    };
    const Vec3 top = texel(x0, y0) * (1 - fx) + texel(x0 + 1, y0) * fx, bottom = texel(x0, y0 + 1) * (1 - fx) + texel(x0 + 1, y0 + 1) * fx;
    return top * (1 - fy) + bottom * fy;
}

Vec3 triangleTextureColour(const ColourImage& t, const std::array<Vec2, 3>& uv) {
    // Interior sample points (barycentric): the centroid and three points towards the corners.
    static const double B[4][3] = {{1 / 3.0, 1 / 3.0, 1 / 3.0}, {2 / 3.0, 1 / 6.0, 1 / 6.0}, {1 / 6.0, 2 / 3.0, 1 / 6.0}, {1 / 6.0, 1 / 6.0, 2 / 3.0}};
    Vec3 c{0, 0, 0};
    for (const auto& b : B) {
        const uint32_t texel = nearestTexel(t, {b[0] * uv[0].x + b[1] * uv[1].x + b[2] * uv[2].x, b[0] * uv[0].y + b[1] * uv[1].y + b[2] * uv[2].y});
        c = c + Vec3{srgbByteToLinear(texel >> 16), srgbByteToLinear((texel >> 8) & 255), srgbByteToLinear(texel & 255)};
    }
    return c * 0.25;
}

} // namespace knit
