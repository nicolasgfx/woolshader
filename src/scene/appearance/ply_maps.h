#pragma once
// The shading maps of the plies (knit_check --ply-map-sets): their files, the
// constants a material needs (ply_maps.json) and their content-addressed cache. Attaching a set to a
// package is the export's (export/collect_scene.h).
#include <optional>
#include <string>
#include <vector>

#include "knitting/tube/ply_profile.h"

namespace knit {

struct PlyMapSpec {
    PlyGeometry g;
    PlyProfileParams pp;
    double yarnD = 1.12;    // reference diameter of the map set [mm]
    int size = 1024;        // width = height [texels]
    double sigma = 68.103;  // library mm per scene unit (bump_scale is in scene metres)
    double cavity = 0.5;    // k in the albedo multiplier c = 1 - k (1 - ao)
};

// Files (names relative to the output directory) and the constants a material needs.
//   height  ply_height.exr     float, linear: H = (rho - rhoMin) / (rhoMax - rhoMin), R = G = B
//   depth   ply_depth.exr      float, 1 - H: for displacement that keeps the crowns on the tube
//                              (photon_path_tracer shape displacement: units_to_object = -bumpScale)
//   normal  ply_normal.png     8-bit, linear: (n + 1) / 2, T = +u, B = +v (green up, "OpenGL")
//   normalFlipGreen            the same with green inverted (convention fixture only)
//   ao      ply_ao.png         8-bit, linear: cross-section AO (R = G = B)
//   cavity  ply_cavity.png     8-bit, sRGB-encoded c = 1 - k (1 - ao): decodes to c
struct PlyMapSet {
    std::string height, depth, normal, normalFlipGreen, ao, cavity;
    double uvRepeat = 0;   // = pitch: the tube's u period [mm]
    double bumpScale = 0;  // (rhoMax - rhoMin) / sigma [scene units]
    PlyShadingProfile prof;
};

// The spec's file names, profile and constants, nothing written.
PlyMapSet plyMapSetOf(const PlyMapSpec& spec);
// Writes the maps and <dir>/ply_maps.json, the constants a material needs.
bool writePlyMaps(const std::string& dir, const PlyMapSpec& spec, PlyMapSet& out, std::string* err);
// <dir>/ply_maps.json alone: uv repeat, bump scale at spec.sigma, cavity k and the profile. The one writer
// of that file, for the cache and for a package.
void writePlyMapsJson(const std::string& dir, const PlyMapSpec& spec, const PlyMapSet& set);

// The cache, one directory per key under out/cache/plymaps/<id>:
//   * the key is the canonical text of everything the maps depend on: ply geometry (count, pitch,
//     fill, handedness, phase), the profile parameters, the map diameter, size, cavity k and a
//     generator version. NOT sigma: bump_scale is computed per scene when a set is attached, so one
//     set serves every scene at the same D / pitch;
//   * a set is generated into a staging directory, every map re-read (decode, dimensions) and then
//     published (publish.h). A cache hit is a verified publication of the key (isPublished); staging,
//     unmarked or corrupt directories never are and are regenerated;
//   * cancellation (cancel.h) aborts before publication: the previous state is kept.
std::string plyMapSetKey(const PlyMapSpec& spec);  // canonical text
std::string plyMapSetId(const PlyMapSpec& spec);   // 16 hex digits of FNV-1a 64 of the key
std::string plyMapCacheRoot();                     // <repo>/out/cache/plymaps

struct PlyMapSetResult {
    std::string dir;   // the published set
    bool hit = false;  // found verified (not generated)
    PlyMapSet set;     // file names, uvRepeat and the profile; bumpScale for sigma = 1 (mm)
    double ms = 0;
};
bool ensurePlyMapSet(const std::string& root, const PlyMapSpec& spec, PlyMapSetResult& out, std::string* err);

// The default sets at map diameter D: {2, 3, 4} plies x {Swept, Circles}, the given ply
// geometry otherwise.
std::vector<PlyMapSpec> defaultPlyMapSpecs(double yarnD, const PlyGeometry& g);

} // namespace knit
