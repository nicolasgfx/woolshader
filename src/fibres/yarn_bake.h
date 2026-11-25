#pragma once
// One yarn's fibre surface baked into maps on its tube (yarn_body.render = maps):
// height, normal, specular share, ambient occlusion and normal spread, the geometry of the fibres on the
// plain single-ply tube; a package derives the cavity and roughness maps its material shows from them
// (cavityMapOf, roughnessMapOf), so the body's fibres need not be drawn and a material edit never bakes.
//
// The body (yarn_body.h) builds every fibre in a straight nominal yarn (x, y, z), z the arc length, and
// bends it with the tube's own rings: X = c(s) + [n1 n2] M(s) (x, y). The tube mesh's texture
// coordinates are u ~ s and v = theta / 2 pi with theta the angle from n1. So the scan of a yarn is a
// grid of texels (s, theta) in nominal space: the fibre points go back to (x, y, s) by M^-1, and at every
// texel a radial ray from outside the yarn towards its axis finds the first fibre. The hit's height,
// normal and direction become the maps. A ray that meets no fibre looks through a gap between fibres:
// it meets the ply's core (the round body of the ply inside its shell, where the core proxy lies), whose
// height and normal it takes, so a gap reads as deep and round, not as the smooth tube. A ray that meets
// no core either (between the plies) stops at the ply centres with the normal of its nearest ply's flank,
// so the valley between two plies reads as one. The ambient occlusion is the product of two shadings: the
// horizon scan at the ply's scale on the height field smoothed over two fibre diameters (the valleys between
// the plies), and the shell's transmittance down to the hit, exp(-depth / mean free path) below the local
// fibre envelope (a fibre under others, and a gap's core, read dark; the top fibres stay lit).
//
// A yarn's maps form one atlas: the chart's s axis is cut into strips of atlasWidth texels (the
// renderer's largest texture side), stacked in v. The tube's uv is relabelled to match
// (collect_scene.cpp relabelAtlasUv): u = the position in the strip, v = (strip + theta / 2 pi) / strips.
#include <cstdint>
#include <string>
#include <vector>

#include "fibres/yarn_body.h"
#include "model/parallel.h"
#include "model/vec.h"

namespace knit {

struct CurveSet;
struct TubeMesh;
struct YarnSurface;

// The baked tube: finer than the export tube (fibres/fibre.h kExportRing 0.4 mm, 8 sides), which
// facets under a normal map at the swatch's 1.9 mm yarn.
constexpr double kBakeRing = 0.2;  // [mm] between tube rings
constexpr int kBakeSides = 16;

struct BakeParams {
    double texel = 0.040;        // [mm] along the yarn and around it (around = round(2 pi r / texel))
    int samples = 3;             // samples x samples stratified rays per texel
    int atlasWidth = 4096;       // texels per strip (the renderer caps a texture side there); tests use less
};

// One yarn's atlas, every map width x height RGB8 with row 0 at the top (v = 1), as io/image_io.h writes.
struct BakedYarnMaps {
    int width = 0, height = 0;     // the atlas: atlasWidth x (strips * around)
    int strips = 0, around = 0;    // strips along s, texels around the yarn
    int along = 0;                 // texels used along s (ceil(length / texel)); the last strip's tail is unused
    double texel = 0, stripMm = 0; // [mm]: texel size, and texel * atlasWidth of yarn per strip
    std::vector<uint8_t> normal;   // tangent space (x along +s, y along +theta, z out), (n + 1) / 2
    std::vector<uint8_t> specular; // grey: the share of rays that met a fibre (gaps reflect nothing)
    std::vector<uint8_t> ao;       // grey, linear: the horizon scan's ambient occlusion times the shell's transmittance
    std::vector<uint8_t> spread;   // grey, linear: Toksvig's spread (1 - |n|) / |n| of the texel's rays, at most 1
    std::vector<uint8_t> heightMap;// grey: the hit's radius over r (a gap: the floor near the ply centres)
    double fibreAngleDeg = 0;      // mean fibre direction against +s on the chart (-90, 90], from the hits
    double coverage = 0;           // share of all rays that met a fibre
    double cores = 0;              // share of all rays that met a ply's core through a gap
    double ms = 0;
};

// The atlas of a yarn of `length` [mm] and radius r, its maps still empty: around, along, strips, width,
// height. Refused when strips * around would exceed the atlas width (a longer yarn, a finer texel).
[[nodiscard]] bool planAtlas(double length, double r, const BakeParams& p, BakedYarnMaps& out, std::string* err);
// The strips planAtlas gives a yarn of `length`: ceil(ceil(length / texel) / atlasWidth).
int atlasStrips(double length, const BakeParams& p);

// Chart texel (i along s, j around) <-> atlas pixel (x, row); row 0 is the top. strip = i / atlasWidth,
// x = i % atlasWidth, row = height - 1 - (strip * around + j). chartTexel is false on the unused tail.
void atlasPixel(const BakedYarnMaps& m, int i, int j, int& x, int& row);
bool chartTexel(const BakedYarnMaps& m, int x, int row, int& i, int& j);

// The tube's uv relabelled for the atlas. The tube was built with uvRepeat = texel * atlasWidth and
// uvWrap = 1, so u runs 0..1 along each strip and a ring is doubled where a strip ends (no triangle
// straddles two strips). Here strip = round(s / stripMm - u) and v = (strip + theta / 2 pi) / strips of
// the yarn (its length: the largest s of its vertices).
void relabelAtlasUv(TubeMesh& tube, const BakeParams& p);

// The material maps a package shows, from the baked ones (RGB8 in, RGB8 out, the same size): the cavity map
// 1 - cavity (1 - AO), sRGB-encoded grey, which the material's albedo (the dye) multiplies; the roughness map
// sqrt(a0^2 + spread) clamped to [a0, 1], linear grey (Toksvig 2005: a0 the smooth fibre's roughness).
std::vector<uint8_t> cavityMapOf(const std::vector<uint8_t>& ao, double cavity);
std::vector<uint8_t> roughnessMapOf(const std::vector<uint8_t>& spread, double baseRoughness);

// The maps of yarn `yarn` (its index in the surface model and BodyStrands::yarnOf) from the body built
// on that surface with `bodyParams` and `partSeed` (they place the plies whose cores the gaps show).
// Refused when strips * around would exceed the atlas width (a longer yarn, a finer texel) or the body
// holds no fibre of the yarn. exec: texel rows in parallel.
[[nodiscard]] bool bakeYarnMaps(const CurveSet& cs, const YarnSurface& ys, const BodyStrands& body, const YarnBodyParams& bodyParams,
                                uint64_t partSeed, int yarn, const BakeParams& p, BakedYarnMaps& out, std::string* err,
                                const ExecPolicy& exec = {});

} // namespace knit
