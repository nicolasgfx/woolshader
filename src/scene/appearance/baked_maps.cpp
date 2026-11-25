// The baked yarn maps' cache; see baked_maps.h. In reading order: the key, the entry's files, the bake on a miss
// (publish.h protocol), the cap, then what a package takes from an entry.
#include "scene/appearance/baked_maps.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <sstream>

#include "io/cache_paths.h"
#include "io/files.h"
#include "io/image_io.h"
#include "io/publish.h"
#include "knitting/tube/param_text.h"
#include "knitting/tube/yarn_surface.h"
#include "model/cancel.h"
#include "model/hash.h"
#include "model/parallel.h"
#include "model/yarn_path.h"
#include "scene/stage_params.h"

namespace knit {

namespace fs = std::filesystem;

namespace {

// Bumped whenever the maps change for the same spec (the bake's rays, maps or atlas layout).
constexpr int kBakeGenerator = 2;  // 2: gaps show the ply cores; the AO is scanned at the ply's scale, times the shell's transmittance
constexpr int kBakeYarnsAtOnce = 8;  // a yarn's bake holds ~200 MB (segments, grid, texel sums, maps)
const char* const kMapNames[5] = {"height", "normal", "specular", "ao", "spread"};
constexpr const char* kPlanFile = "bake.txt";

// --- the entry's files

bool writePlans(const std::string& dir, const std::vector<BakedYarnPlan>& plans) {
    std::ofstream f(dir + "/" + kPlanFile);
    f << "bake v" << kBakeGenerator << " yarns " << plans.size() << "\n";
    for (size_t y = 0; y < plans.size(); ++y) {
        const BakedYarnPlan& p = plans[y];
        char line[200];
        std::snprintf(line, sizeof line, "yarn %zu along %d around %d strips %d width %d height %d angle_deg %.6f coverage %.6f\n", y,
                      p.along, p.around, p.strips, p.width, p.height, p.fibreAngleDeg, p.coverage);
        f << line;
    }
    return static_cast<bool>(f);
}

// The next word of a plan line is this label (each value in bake.txt follows its label).
bool readLabel(std::istream& in, const char* label) {
    std::string word;
    return in >> word && word == label;
}

bool readPlans(const std::string& dir, std::vector<BakedYarnPlan>& plans, std::string* err) {
    std::ifstream f(dir + "/" + kPlanFile);
    std::string line;
    int version = 0;
    size_t yarns = 0;
    std::istringstream header(std::getline(f, line) ? line : std::string());
    char v = 0;
    const bool headerRead = readLabel(header, "bake") && header >> v >> version && v == 'v' && readLabel(header, "yarns") && header >> yarns;
    if (!headerRead || version != kBakeGenerator) return fail(err, "baked maps: " + dir + "/" + kPlanFile + " is not this generator's");
    plans.assign(yarns, BakedYarnPlan{});
    for (size_t y = 0; y < yarns; ++y) {
        BakedYarnPlan& p = plans[y];
        size_t index = 0;
        std::istringstream row(std::getline(f, line) ? line : std::string());
        const bool rowRead = readLabel(row, "yarn") && row >> index && readLabel(row, "along") && row >> p.along && readLabel(row, "around") &&
                             row >> p.around && readLabel(row, "strips") && row >> p.strips && readLabel(row, "width") && row >> p.width &&
                             readLabel(row, "height") && row >> p.height && readLabel(row, "angle_deg") && row >> p.fibreAngleDeg &&
                             readLabel(row, "coverage") && row >> p.coverage;
        if (!rowRead || index != y) return fail(err, "baked maps: " + dir + "/" + kPlanFile + ": yarn " + std::to_string(y) + " unreadable");
    }
    return true;
}

// Every map of every yarn decodes with the plan's size.
bool reread(const std::string& dir, const std::vector<BakedYarnPlan>& plans, std::string* err) {
    for (size_t y = 0; y < plans.size(); ++y)
        for (const char* map : kMapNames) {
            int w = 0, h = 0;
            std::vector<uint8_t> px;
            const std::string file = bakedMapFile(static_cast<int>(y), map);
            if (!readPngRgb8(dir + "/" + file, w, h, px) || w != plans[y].width || h != plans[y].height)
                return fail(err, "baked maps: " + file + " does not read back as " + std::to_string(plans[y].width) + " x " +
                                     std::to_string(plans[y].height));
        }
    return true;
}

bool writeYarnMaps(const std::string& dir, int yarn, const BakedYarnMaps& m) {
    const std::vector<uint8_t>* maps[5] = {&m.heightMap, &m.normal, &m.specular, &m.ao, &m.spread};
    for (int k = 0; k < 5; ++k)
        if (!writePngRgb8(dir + "/" + bakedMapFile(yarn, kMapNames[k]), m.width, m.height, *maps[k])) return false;
    return true;
}

BakedYarnPlan planOf(const BakedYarnMaps& m) {
    BakedYarnPlan p;
    p.along = m.along;
    p.around = m.around;
    p.strips = m.strips;
    p.width = m.width;
    p.height = m.height;
    p.fibreAngleDeg = m.fibreAngleDeg;
    p.coverage = m.coverage;
    return p;
}

BakeParams bakeParamsOf(const BakedMapSpec& spec) {
    BakeParams p;
    p.texel = spec.texel;
    p.samples = spec.samples;
    p.atlasWidth = spec.atlasWidth;
    return p;
}

// --- the cap: the newest kBakeCacheEntries published entries stay (by their marker's time), `keep` always

void dropOldEntries(const std::string& root, const std::string& keep) {
    std::error_code ec;
    std::vector<std::pair<fs::file_time_type, fs::path>> entries;
    for (const fs::directory_entry& e : fs::directory_iterator(root, ec)) {
        if (!e.is_directory(ec)) continue;
        const fs::path marker = e.path() / kPublishMarker;
        if (!fs::exists(marker, ec) || fs::equivalent(e.path(), fs::path(keep), ec)) continue;
        entries.push_back({fs::last_write_time(marker, ec), e.path()});
    }
    std::sort(entries.begin(), entries.end(), [](const std::pair<fs::file_time_type, fs::path>& a, const std::pair<fs::file_time_type, fs::path>& b) {
        return a.first > b.first;  // newest first
    });
    for (size_t k = static_cast<size_t>(kBakeCacheEntries) - 1; k < entries.size(); ++k) fs::remove_all(entries[k].second, ec);
}

// --- the bake of every yarn into the staging directory

bool bakeAllYarns(const std::string& staging, const BakedMapSpec& spec, const CurveSet& cs, const SurfaceModel& surface,
                  const BakeProgress& progress, std::vector<BakedYarnPlan>& plans, std::string* err) {
    YarnBody body;
    if (!buildYarnBody(cs, surface, spec.body, spec.partSeed, body, err, ExecPolicy{true})) return false;
    const int yarns = static_cast<int>(surface.yarns.size());
    plans.assign(static_cast<size_t>(yarns), BakedYarnPlan{});
    std::vector<std::string> failures(static_cast<size_t>(yarns));
    std::atomic<int> done{0};
    std::mutex reporting;
    const BakeParams bp = bakeParamsOf(spec);
    ExecPolicy policy;
    policy.parallel = true;
    policy.maxConcurrency = kBakeYarnsAtOnce;
    parallelFor(yarns, 1, [&](int64_t begin, int64_t end, int64_t) {
        for (int y = static_cast<int>(begin); y < static_cast<int>(end); ++y) {
            BakedYarnMaps maps;
            std::string why;
            if (!bakeYarnMaps(cs, surface.yarns[static_cast<size_t>(y)], body.strands, spec.body, spec.partSeed, y, bp, maps, &why) || !writeYarnMaps(staging, y, maps)) {
                failures[static_cast<size_t>(y)] = why.empty() ? "cannot write the maps of yarn " + std::to_string(y) : why;
                continue;
            }
            plans[static_cast<size_t>(y)] = planOf(maps);
            const int n = ++done;
            if (progress.report) {
                const std::lock_guard<std::mutex> lock(reporting);
                progress.report(n, yarns);
            }
        }
    }, policy);
    for (const std::string& f : failures)
        if (!f.empty()) return fail(err, f);
    return true;
}

} // namespace

std::string bakedMapFile(int yarn, const char* map) {
    char b[48];
    std::snprintf(b, sizeof b, "y%03d_%s.png", yarn, map);
    return b;
}

bool bakedMapSpecOf(const SceneSettings& s, const std::string& partName, const CurveSet& strands, const PatternParams& pp, BakedMapSpec& out,
                    std::string* err) {
    out = BakedMapSpec{};
    if (!yarnBodyParams(s, out.body, err)) return false;
    out.surfaceDigest = surfaceBuildDigest(strands, pp, yarnSurfaceParams(s));
    out.partSeed = fnv1a64(partName);
    const BakeParams bp = bakeParams(s);
    out.texel = bp.texel;
    out.samples = bp.samples;
    out.atlasWidth = bp.atlasWidth;
    return true;
}

std::string bakedMapKey(const BakedMapSpec& spec) {
    const YarnBodyParams& b = spec.body;
    char t[512];
    std::snprintf(t, sizeof t,
                  "bakes v%d surface=%016llx body=plies %d fibres %d interior %d radius %.17g shell %.17g fibre_angle %.17g ply_angle %.17g "
                  "irregularity %.17g seed %llu part_seed=%016llx texel=%.17g samples=%d atlas=%d",
                  kBakeGenerator, static_cast<unsigned long long>(spec.surfaceDigest), b.plies, b.fibres, b.interiorFibres, b.fibreRadius,
                  b.shellDepth, b.fibreAngleDeg, b.plyAngleDeg, b.irregularity, static_cast<unsigned long long>(b.seed),
                  static_cast<unsigned long long>(spec.partSeed), spec.texel, spec.samples, spec.atlasWidth);
    return t;
}

std::string bakedMapId(const BakedMapSpec& spec) {
    const std::string k = bakedMapKey(spec);
    char b[24];
    std::snprintf(b, sizeof b, "%016llx", static_cast<unsigned long long>(fnv1a(k.data(), k.size())));
    return b;
}

std::string bakeCacheRoot() { return outputCacheDir("bakes"); }

bool loadBakedMapSet(const std::string& root, const BakedMapSpec& spec, BakedMapEntry& out) {
    out = BakedMapEntry{};
    out.dir = root + "/" + bakedMapId(spec);
    recoverPublished(out.dir);
    if (!isPublished(out.dir, bakedMapKey(spec))) return false;
    std::string why;
    if (!readPlans(out.dir, out.yarns, &why)) return false;
    out.hit = true;
    return true;
}

bool ensureBakedMapSet(const std::string& root, const BakedMapSpec& spec, const CurveSet& cs, const SurfaceModel& surface,
                       BakeInterrupt interrupt, const BakeProgress& progress, BakedMapEntry& out, std::string* err) {
    const auto t0 = std::chrono::steady_clock::now();
    if (loadBakedMapSet(root, spec, out)) {
        out.ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        return true;
    }
    const std::string key = bakedMapKey(spec);
    {
        std::error_code ec;
        fs::create_directories(root, ec);
        Publisher pub;
        if (!pub.begin(out.dir, err, PublishPolicy::CacheEntry)) return false;
        // FinishOnce: no token on this thread from the bake to the publish, so a preview key change (which cancels
        // the exporter's token) waits for the entry instead of throwing minutes away; the cancel acts right after.
        const ThreadCancelScope quiet(interrupt == BakeInterrupt::FinishOnce ? nullptr : threadCancelToken());
        if (!bakeAllYarns(pub.staging(), spec, cs, surface, progress, out.yarns, err)) return false;
        cancelCheckpoint();
        if (!writePlans(pub.staging(), out.yarns)) return fail(err, "baked maps: cannot write " + std::string(kPlanFile));
        if (!reread(pub.staging(), out.yarns, err)) return false;
        cancelCheckpoint();
        if (isPublished(out.dir, key)) {
            pub.abort();  // another thread baked the same entry meanwhile: theirs stands
            if (!readPlans(out.dir, out.yarns, err)) return false;
            out.hit = true;
        } else if (!pub.commit(key, err)) {
            return false;
        }
        dropOldEntries(root, out.dir);
    }
    out.ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    return true;
}

bool bakedCavityMap(const BakedMapEntry& entry, int yarn, double cavity, int& width, int& height, std::vector<uint8_t>& rgb, std::string* err) {
    if (yarn < 0 || yarn >= static_cast<int>(entry.yarns.size())) return fail(err, "baked maps: no yarn " + std::to_string(yarn) + " in " + entry.dir);
    std::vector<uint8_t> ao;
    if (!readPngRgb8(entry.dir + "/" + bakedMapFile(yarn, "ao"), width, height, ao))
        return fail(err, "baked maps: the ao map of yarn " + std::to_string(yarn) + " does not read back from " + entry.dir);
    rgb = cavityMapOf(ao, cavity);
    return true;
}

bool bakedPlanMatches(const BakedYarnPlan& plan, double yarnLength, double r, const BakedMapSpec& spec) {
    BakedYarnMaps m;
    if (!planAtlas(yarnLength, r, bakeParamsOf(spec), m, nullptr)) return false;
    return m.along == plan.along && m.around == plan.around && m.strips == plan.strips && m.width == plan.width && m.height == plan.height;
}

bool attachBakedYarn(const BakedMapEntry& entry, int yarn, const std::string& packageDir, const std::string& stem, double cavity,
                     double baseRoughness, BakedYarnFiles& out, std::string* err) {
    if (yarn < 0 || yarn >= static_cast<int>(entry.yarns.size())) return fail(err, "baked maps: no yarn " + std::to_string(yarn) + " in " + entry.dir);
    std::error_code ec;
    fs::create_directories(fs::path(packageDir + "/" + stem).parent_path(), ec);
    const BakedYarnPlan& plan = entry.yarns[static_cast<size_t>(yarn)];
    int w = 0, h = 0;
    std::vector<uint8_t> ao, spread;
    if (!readPngRgb8(entry.dir + "/" + bakedMapFile(yarn, "ao"), w, h, ao) || w != plan.width || h != plan.height ||
        !readPngRgb8(entry.dir + "/" + bakedMapFile(yarn, "spread"), w, h, spread) || w != plan.width || h != plan.height)
        return fail(err, "baked maps: the ao or spread map of yarn " + std::to_string(yarn) + " does not read back from " + entry.dir);
    out.cavity = stem + "_cavity.png";
    out.roughness = stem + "_roughness.png";
    out.normal = stem + "_normal.png";
    out.specular = stem + "_specular.png";
    out.height = stem + "_height.png";
    if (!writePngRgb8(packageDir + "/" + out.cavity, w, h, cavityMapOf(ao, cavity)) ||
        !writePngRgb8(packageDir + "/" + out.roughness, w, h, roughnessMapOf(spread, baseRoughness)))
        return fail(err, "baked maps: cannot write the material maps of yarn " + std::to_string(yarn));
    const std::pair<const char*, const std::string*> copies[3] = {{"normal", &out.normal}, {"specular", &out.specular}, {"height", &out.height}};
    for (const auto& [map, target] : copies) {
        fs::copy_file(entry.dir + "/" + bakedMapFile(yarn, map), packageDir + "/" + *target, fs::copy_options::overwrite_existing, ec);
        if (ec) return fail(err, "baked maps: cannot copy " + bakedMapFile(yarn, map) + ": " + ec.message());
    }
    return true;
}

} // namespace knit
