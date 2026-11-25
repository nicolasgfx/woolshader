// The ply shading maps (ply_maps.h): written with the constants a material needs, and kept in their
// content-addressed cache.
#include "scene/appearance/ply_maps.h"
#include "io/cache_paths.h"
#include "model/choices.h"
#include "model/colour.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <vector>

#include "io/image_io.h"
#include "io/publish.h"
#include "model/cancel.h"
#include "model/hash.h"

namespace knit {

PlyMapSet plyMapSetOf(const PlyMapSpec& spec) {
    PlyMapSet set;
    set.prof = buildPlyShadingProfile(spec.g, 0.5 * spec.yarnD, spec.pp);
    set.uvRepeat = spec.g.pitch;
    set.bumpScale = (set.prof.rhoMax - set.prof.rhoMin) / spec.sigma;
    set.height = "ply_height.exr";
    set.depth = "ply_depth.exr";
    set.normal = "ply_normal.png";
    set.normalFlipGreen = "ply_normal_flipg.png";
    set.ao = "ply_ao.png";
    set.cavity = "ply_cavity.png";
    return set;
}

bool writePlyMaps(const std::string& dir, const PlyMapSpec& spec, PlyMapSet& out, std::string* err) {
    if (spec.g.count <= 1) {
        if (err) *err = "one ply: no maps (the plain yarn tube)";
        return false;
    }
    std::filesystem::create_directories(dir);
    out = plyMapSetOf(spec);
    const int W = spec.size, H = spec.size;
    const PlyMaps m = makePlyMaps(out.prof, out.uvRepeat, W, H);
    const size_t n = size_t(W) * H;
    std::vector<float> height(3 * n), depth(3 * n);
    std::vector<uint8_t> nrm(3 * n), flip(3 * n), ao(3 * n), cav(3 * n);
    auto q8 = [](double x) { return uint8_t(std::clamp(std::lround(x * 255), 0L, 255L)); };
    for (size_t i = 0; i < n; ++i) {
        height[3 * i] = height[3 * i + 1] = height[3 * i + 2] = m.H[i];
        depth[3 * i] = depth[3 * i + 1] = depth[3 * i + 2] = 1 - m.H[i];
        for (int c = 0; c < 3; ++c) nrm[3 * i + c] = q8(0.5 * (m.normal[3 * i + c] + 1));
        flip[3 * i] = nrm[3 * i];
        flip[3 * i + 1] = q8(0.5 * (1 - m.normal[3 * i + 1]));
        flip[3 * i + 2] = nrm[3 * i + 2];
        ao[3 * i] = ao[3 * i + 1] = ao[3 * i + 2] = q8(m.ao[i]);
        cav[3 * i] = cav[3 * i + 1] = cav[3 * i + 2] = q8(srgbEncode(1 - spec.cavity * (1 - m.ao[i])));
    }
    const bool ok = writeExrRgbF32(dir + "/" + out.height, W, H, height) && writeExrRgbF32(dir + "/" + out.depth, W, H, depth) &&
                    writePngRgb8(dir + "/" + out.normal, W, H, nrm) &&
                    writePngRgb8(dir + "/" + out.normalFlipGreen, W, H, flip) && writePngRgb8(dir + "/" + out.ao, W, H, ao) &&
                    writePngRgb8(dir + "/" + out.cavity, W, H, cav);
    if (!ok && err) *err = "cannot write the maps in " + dir;
    // The constants a material needs.
    writePlyMapsJson(dir, spec, out);
    return ok;
}

void writePlyMapsJson(const std::string& dir, const PlyMapSpec& spec, const PlyMapSet& out) {
    char b[512];
    std::snprintf(b, sizeof b,
                  "{\n  \"material\": {\"uv_repeat_mm\": %.6g, \"bump_scale\": %.8g, \"sigma\": %.6g, \"cavity_k\": %.4g},\n"
                  "  \"profile\": {\"rho_min_mm\": %.6g, \"rho_max_mm\": %.6g, \"yarn_d_mm\": %.6g, \"plies\": %d}\n}\n",
                  out.uvRepeat, out.bumpScale, spec.sigma, spec.cavity, out.prof.rhoMin, out.prof.rhoMax, spec.yarnD, spec.g.count);
    std::ofstream(dir + "/ply_maps.json") << b;
}

namespace fs = std::filesystem;

namespace {

// Bumped whenever writePlyMaps' output changes for the same spec.
constexpr int kPlyMapGenerator = 1;

// The cached set is written for sigma = 1 (its ply_maps.json gives bump_scale in mm).
PlyMapSpec cached(const PlyMapSpec& spec) {
    PlyMapSpec c = spec;
    c.sigma = 1;
    return c;
}

bool reread(const std::string& dir, const PlyMapSet& s, int size, std::string* err) {
    int w = 0, h = 0;
    for (const std::string& f : {s.height, s.depth}) {
        std::vector<float> px;
        if (!readExrRgbF32(dir + "/" + f, w, h, px) || w != size || h != size) {
            if (err) *err = "map set: " + f + " does not read back as " + std::to_string(size) + "^2";
            return false;
        }
    }
    for (const std::string& f : {s.normal, s.normalFlipGreen, s.ao, s.cavity}) {
        std::vector<uint8_t> px;
        if (!readPngRgb8(dir + "/" + f, w, h, px) || w != size || h != size) {
            if (err) *err = "map set: " + f + " does not read back as " + std::to_string(size) + "^2";
            return false;
        }
    }
    return true;
}

} // namespace

std::string plyMapSetKey(const PlyMapSpec& spec) {
    const PlyGeometry& g = spec.g;
    const PlyProfileParams& p = spec.pp;
    char b[640];
    std::snprintf(b, sizeof b,
                  "plymaps v%d plies=%d pitch=%.17g hand=%d phase=%.17g fill=%.17g model=%s blend_tilt=%.17g min_blend=%.17g "
                  "gap_floor=%.17g ao=%s ao_dirs=%d table=%d yarn_d=%.17g size=%d cavity=%.17g",
                  kPlyMapGenerator, g.count, g.pitch, g.handedness, g.phase, g.fill, p.model == PlyProfileModel::Swept ? "swept" : "circles",
                  p.blendTiltDeg, p.minBlend, p.gapFloor, "height-field",  // the one AO since the cross-section variant went: the key text stays
                  p.aoDirections, p.tableSamples, spec.yarnD, spec.size, spec.cavity);
    return b;
}

std::string plyMapSetId(const PlyMapSpec& spec) {
    const std::string k = plyMapSetKey(spec);
    char b[24];
    std::snprintf(b, sizeof b, "%016llx", static_cast<unsigned long long>(fnv1a(k.data(), k.size())));
    return b;
}

std::string plyMapCacheRoot() { return outputCacheDir("plymaps"); }

bool ensurePlyMapSet(const std::string& root, const PlyMapSpec& spec, PlyMapSetResult& out, std::string* err) {
    const auto t0 = std::chrono::steady_clock::now();
    out = PlyMapSetResult{};
    const PlyMapSpec c = cached(spec);
    const std::string key = plyMapSetKey(c);
    out.dir = root + "/" + plyMapSetId(c);
    recoverPublished(out.dir);
    out.set = plyMapSetOf(c);
    if (isPublished(out.dir, key)) {
        out.hit = true;
    } else {
        std::error_code ec;
        fs::create_directories(root, ec);
        Publisher pub;
        if (!pub.begin(out.dir, err, PublishPolicy::CacheEntry)) return false;
        PlyMapSet written;
        if (!writePlyMaps(pub.staging(), c, written, err)) return false;
        cancelCheckpoint();
        if (!reread(pub.staging(), written, c.size, err)) return false;
        cancelCheckpoint();
        if (!pub.commit(key, err)) return false;
        out.set = std::move(written);
    }
    out.ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    return true;
}

std::vector<PlyMapSpec> defaultPlyMapSpecs(double yarnD, const PlyGeometry& g) {
    std::vector<PlyMapSpec> v;
    for (PlyProfileModel m : {PlyProfileModel::Swept, PlyProfileModel::Circles})
        for (int n = 2; n <= 4; ++n) {
            PlyMapSpec s;
            s.g = g;
            s.g.count = n;
            s.pp.model = m;
            s.yarnD = yarnD;
            v.push_back(s);
        }
    return v;
}

} // namespace knit
