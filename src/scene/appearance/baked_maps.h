#pragma once
// The maps baked from a yarn body's fibres (fibres/yarn_bake.h), one atlas per yarn, in a content-addressed
// cache, and the material maps a package shows from them. The bake is a pure function of the knit's curves,
// the body's parameters and the texel, and it takes minutes, so it is baked once per identity and attached
// to every export and render preview from the cache (the ply map set, ply_maps.h, is the precedent).
//
// What the entry holds is geometry only: height, normal, specular share, ambient occlusion and Toksvig's
// normal spread. The material values (the dye, the cavity strength, the roughness floor, the highlight) are
// applied when a package attaches a yarn (attachBakedYarn), so a slider drag never bakes again.
#include <functional>
#include <string>
#include <vector>

#include "fibres/yarn_bake.h"
#include "fibres/yarn_body.h"

namespace knit {

struct CurveSet;
struct PatternParams;
struct SceneSettings;
struct SurfaceModel;

// The key's inputs. Not inputs: the core proxy (the bake never sees it), the material values, the tube mesh.
struct BakedMapSpec {
    uint64_t surfaceDigest = 0;  // surfaceBuildDigest (knitting/tube/param_text.h): the curves at the rendered diameter, compression
    YarnBodyParams body;         // the fibres' structure and seed
    uint64_t partSeed = 0;       // fnv1a64 of the part's name, as the export seeds the body
    double texel = 0.040;        // [mm]
    int samples = 3;             // rays per texel side
    int atlasWidth = 4096;       // the renderer's largest texture side
};

// One yarn's atlas as planned and measured by the bake (bake.txt in the entry).
struct BakedYarnPlan {
    int along = 0, around = 0, strips = 0, width = 0, height = 0;
    double fibreAngleDeg = 0, coverage = 0;
};

// A published entry: out/cache/bakes/<16 hex>/y%03d_{height,normal,specular,ao,spread}.png and bake.txt.
struct BakedMapEntry {
    std::string dir;
    bool hit = false;                 // found verified, not baked now
    std::vector<BakedYarnPlan> yarns;
    double ms = 0;
};

// The spec of a knitted part's maps from the settings and the part's rendered strands (at the rendered diameter
// and the surface the export meshes): the one place the export and the viewport compute it. False with the
// reason when the body's settings are refused.
[[nodiscard]] bool bakedMapSpecOf(const SceneSettings& s, const std::string& partName, const CurveSet& strands, const PatternParams& pp,
                                  BakedMapSpec& out, std::string* err);
std::string bakedMapKey(const BakedMapSpec& spec);  // canonical text
std::string bakedMapId(const BakedMapSpec& spec);   // 16 hex digits of FNV-1a 64 of the key
std::string bakeCacheRoot();                        // <repo>/out/cache/bakes
std::string bakedMapFile(int yarn, const char* map); // "y%03d_<map>.png"
// The newest entries kept after a bake publishes (each a few hundred MB: the one cache here that grows by the GB).
constexpr int kBakeCacheEntries = 8;

// Whether a bake in flight may be cancelled by the thread's token. FinishOnce: the render preview's exporter
// is cancelled whenever its package key changes; a bake it started is cache-bound work and finishes once.
enum class BakeInterrupt : int { Cancellable, FinishOnce };
// Progress of a bake, a yarn at a time; an empty function is silent.
struct BakeProgress {
    std::function<void(int yarnsDone, int yarns)> report;
};

// The entry of `spec`: a hit returns the published entry; a miss builds the body on `surface` (the same surface
// model the export meshes the tube from), bakes every yarn (at most 8 at once), writes bake.txt, rereads every
// map, publishes, then drops entries beyond kBakeCacheEntries. False with the reason; a cancelled bake leaves
// no entry.
[[nodiscard]] bool ensureBakedMapSet(const std::string& root, const BakedMapSpec& spec, const CurveSet& cs, const SurfaceModel& surface,
                                     BakeInterrupt interrupt, const BakeProgress& progress, BakedMapEntry& out, std::string* err);

// The published entry of `spec`, if the cache holds it: nothing is baked. False without one.
bool loadBakedMapSet(const std::string& root, const BakedMapSpec& spec, BakedMapEntry& out);

// The atlas a tube of `yarnLength` and radius r at `spec`'s texel would get is the plan's: the maps fit the tube.
bool bakedPlanMatches(const BakedYarnPlan& plan, double yarnLength, double r, const BakedMapSpec& spec);

// One yarn's cavity map from an entry (RGB8 grey, rows top first, the atlas's size): the viewport's texture.
[[nodiscard]] bool bakedCavityMap(const BakedMapEntry& entry, int yarn, double cavity, int& width, int& height, std::vector<uint8_t>& rgb,
                                  std::string* err);

// The five files of one yarn as written into a package, relative to it.
struct BakedYarnFiles {
    std::string cavity, normal, specular, roughness, height;
};
// One yarn's maps into a package: the cavity map 1 - cavity (1 - AO) and the roughness map sqrt(a0^2 + spread)
// derived from the entry, the others copied, as <stem>_{cavity,normal,specular,roughness,height}.png under
// packageDir (`stem` relative to it, e.g. "textures/baked_wool_white_y003").
[[nodiscard]] bool attachBakedYarn(const BakedMapEntry& entry, int yarn, const std::string& packageDir, const std::string& stem,
                                   double cavity, double baseRoughness, BakedYarnFiles& out, std::string* err);

} // namespace knit
