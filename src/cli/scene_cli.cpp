// knit_check --scene <id> (S6): run a registered scene from the scene settings.
//
//   knit_check --scene teddy [--settings <file>] [--set key=value]... [--write-settings <file>]
//              [--export-ppt <dir>] [--export-format blend|ppt|obj] [--fiber-check] [--closeup N] [--part name] [-v]
//              (--fiber-check: the fibre diagnostics, FibreRunOptions::fibreCheck; without it the light checks, as the viewer)
//              [--stage-profile] [--fibre-digest] [--threads N]
//              [--save-knit <dir>] [--load-knit <dir>]   spiral pieces (the teddy): the per-piece knit cache, knit_cache.h;
//                                                          a scene knitted on its mesh (imports): <dir>/knit.bin (saved_knit.h)
//                                                          with <dir>/knit.key, loaded only for the same knitting settings
//              [--save-fibres <dir>] [--load-fibres <dir>]   (wool fibres per part, fibre_files.h: loaded instead of grown)
//              [--taper-points dense|existing|onset|refined] [--taper-tolerance f] [--taper-report <file.json>]
//                                                          (developer: where the late taper's radii are carried,
//                                                          StrandTaper; the report measures the written strands)
//              [--body-file-segments N]  (developer: the most segments per yarn-body strand file, 2.5 M)
//              [--fibres-full]                             (= --set fibres.generation=full: every candidate)
//              [--fibre-support-mb MB]                     (bounded fibre-support residency in budget mode)
//              [--time-cap <s>]                            hard cap: after s seconds the run stops with exit code 4
//              [--part-cap <s>]                            each part's time cap (the run's cap is knit.time_cap_s)
//              [--stitch-spread]                           how far plain stitches stray from their mean shape
//   knit_check --scene import --set import.identity=<16 hex> --set import.source=<file> ...
//              a prepared import (out/cache/imports, or KNIT_IMPORT_CACHE), as the viewer knits it
//   knit_check --scene list            the registry
//   knit_check --scene <id> --cache-build [--cache dir] [--settings f] [--set k=v]...
//              computes what the scene viewer shows for the viewer's defaults (the preview
//              profile, scene_builder.h) or the given settings, and stores it in the
//              persistent scene cache (scene_cache.h), so the viewer opens it at once.
//
// The settings file is "knit-settings v1" (src/scene/scene_settings.h); --set applies
// after it, in order. Run-control flags (export, checks, parts) are not settings.
// Every scene runs one path, the viewer's: knitScene knits it (the teddy, knitted, as spiral
// pieces), growSceneFibres grows its fibres and exportScene writes them with the knit. With
// knit.enabled = false the teddy is its GLB mesh, as every scene.
//
// Exit code: 0 every selected part has its result (the report gives each verdict, FAIL included);
// 1 incomplete: a selected part lacks its knit (partIncomplete: a time cap, the device memory budget,
// a spiral piece that could not be knitted), the export is still written; 2 an error (refused
// settings or flags, a failed run or export); 4 the --time-cap watchdog ended the run.
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <thread>
#include <sstream>
#include <string>

#include "knitting/stitch/hand_tension.h"
#include "scene/fibres/fibre_files.h"
#include "model/parallel.h"
#include "io/glb_load.h"
#include "jobs/scene_builder.h"
#include "scene/source/sewn_pieces.h"
#include "jobs/scene_cache.h"
#include "scene/stage_params.h"
#include "settings/scene_registry.h"
#include "settings/scene_settings.h"
#include "settings/settings_profiles.h"
#include "settings/yarn_look.h"
#include "scene/run_log.h"
#include "scene/scene.h"
#include "settings/settings_mapping.h"
#include "scene/source/scene_parts.h"
#include "scene/knit/knit_stage.h"
#include "scene/knit/knit_report.h"
#include "scene/knit/saved_knit.h"
#include "scene/result_keys.h"
#include "scene/fibres/fibre_stage.h"
#include "io/files.h"
#include "export/export.h"

using namespace knit;

namespace {

// Flags of retired routes: refused with what replaces them, never silently ignored. Empty for any
// other argument.
std::string retiredFlagMessage(const std::string& a) {
    const std::string route = " was a flag of the teddy's own pieces route, retired (S2): the teddy's fibres and export are the "
                              "common ones (growSceneFibres, exportScene)";
    if (a == "--piece") return "--piece N" + route + "; one piece: --part <piece name>, e.g. --part piece_bear_002_00";
    if (a == "--fibre-trace") return a + " (per-piece candidate outcome sidecars)" + route;
    if (a == "--fibre-domain-dump" || a == "--fibre-domain-box") return a + " (exposure domain dumps)" + route;
    if (a == "--no-root-check") return a + " (CPU fibres without the serialized root check)" + route;
    if (a == "--fiber-samples") return a + " (the gap-ray check)" + route;
    if (a == "--surface-diag") return a + " (the compression fit variants and their causes) was retired: one fit remains";
    if (a == "--fiber-profile") return a + " (per-stage fibre timers) was retired: --stage-profile keeps the wall times";
    if (a == "--coverage-report") return a + " (per-part verdicts and coverage as JSON lines) was retired: the report and report.json keep the verdicts";
    if (a == "--mesh-route" || a == "--xatlas-charts") return a + " (the teddy on its own mesh) was retired: the experiment ended, the teddy knits as spiral pieces";
    if (a == "--jobs") return a + " was retired: parts knit one at a time, because the device solves them one at a time";
    if (a == "--fibre-regions") return a + " (the CPU generator's hidden-test cells) was retired with the visible-sides placement";
    if (a == "--dump-relax" || a == "--dump-filter" || a == "--dump-exit" || a == "--dump-trace")
        return a + " (relax-call dumps) was retired: nothing writes or reads those files any more";
    if (a == "--no-attachment-check")
        return a + route + "; the common export runs the attachment and reference checks only with --fiber-check";
    if (a == "--viewer-checks" || a == "--no-gpu-audit" || a == "--no-fibre-angles")
        return a + " was retired: the light checks are the default, as the viewer's; --fiber-check runs the CPU audit, the angles and the checks";
    return "";
}

// What the command line asks for.
struct SceneCommandOptions {
    std::string id, settingsFile, writeFile, exportDir, exportFormat = "ppt";
    std::vector<std::pair<std::string, std::string>> sets;
    int closeup = -1;
    std::string onlyPart;
    bool fibreCheck = false, verbose = false, noSelf = false, loadOnly = false, cacheBuild = false;
    std::string cacheRoot = sceneCacheRoot();
    bool stageProfile = false, fibreDigest = false;  // profiling (run_log.h), digests
    bool stitchSpread = false;
    int threads = 0;
    std::string saveKnit, loadKnit, saveFibres, loadFibres;
    std::string taperPoints, taperReport;  // developer: the late taper's points, and its report (fibre_strands.h)
    double taperTolerance = 0.05;
    long long bodyFileSegments = 2'500'000;  // developer: ExportOptions::bodyFileSegments
    bool viewerSettings = false;
    double partCapS = 0;  // --part-cap: each part's time cap [s] (the scene cap is knit.time_cap_s)
    double fibreSupportMB = 3000;
};

constexpr int kContinue = -1;  // a stage's result: go on (any other value is the exit code)

int parseSceneOptions(int argc, char** argv, SceneCommandOptions& o) {
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (const std::string retired = retiredFlagMessage(a); !retired.empty()) {
            std::fprintf(stderr, "%s\n", retired.c_str());
            return 2;
        }
        if (a == "--scene" && i + 1 < argc) o.id = argv[++i];
        else if (a == "--settings" && i + 1 < argc) o.settingsFile = argv[++i];
        else if (a == "--write-settings" && i + 1 < argc) o.writeFile = argv[++i];
        else if (a == "--set" && i + 1 < argc) {
            const std::string kv = argv[++i];
            const size_t eq = kv.find('=');
            if (eq == std::string::npos) { std::fprintf(stderr, "--set %s: expected key=value\n", kv.c_str()); return 2; }
            o.sets.push_back({kv.substr(0, eq), kv.substr(eq + 1)});
        }
        else if (a == "--export-ppt" && i + 1 < argc) o.exportDir = argv[++i];
        else if (a == "--export-format" && i + 1 < argc) o.exportFormat = argv[++i];
        else if (a == "--time-cap" && i + 1 < argc) {
            // A hard cap on the whole run (knit, fibres, export): a run over it is a FAIL, never a
            // result. The watchdog ends the process; nothing partial is reported as done.
            const double cap = std::atof(argv[++i]);
            if (!(cap > 0)) { std::fprintf(stderr, "--time-cap: expected seconds > 0\n"); return 2; }
            std::thread([cap] {
                std::this_thread::sleep_for(std::chrono::duration<double>(cap));
                std::fprintf(stderr, "time cap of %.0f s reached: run stopped (FAIL)\n", cap);
                std::fflush(stdout);
                std::fflush(stderr);
                std::_Exit(4);
            }).detach();
        }
        else if (a == "--closeup" && i + 1 < argc) o.closeup = std::atoi(argv[++i]);
        else if (a == "--part" && i + 1 < argc) o.onlyPart = argv[++i];         // the parts (spiral: pieces) whose name contains it
        else if (a == "--fiber-check") o.fibreCheck = true;
        else if (a == "--stage-profile") o.stageProfile = true;
        else if (a == "--fibre-digest") o.fibreDigest = true;
        else if (a == "--stitch-spread") o.stitchSpread = true;
        else if (a == "--threads" && i + 1 < argc) o.threads = std::atoi(argv[++i]);
        else if (a == "--save-knit" && i + 1 < argc) o.saveKnit = argv[++i];
        else if (a == "--load-knit" && i + 1 < argc) o.loadKnit = argv[++i];
        else if (a == "--save-fibres" && i + 1 < argc) o.saveFibres = argv[++i];
        else if (a == "--load-fibres" && i + 1 < argc) o.loadFibres = argv[++i];
        else if (a == "--taper-points" && i + 1 < argc) {
            o.taperPoints = argv[++i];
            TaperPoints p;
            if (!parseTaperPoints(o.taperPoints, p)) { std::fprintf(stderr, "--taper-points: dense, existing, onset or refined\n"); return 2; }
        }
        else if (a == "--taper-tolerance" && i + 1 < argc) o.taperTolerance = std::atof(argv[++i]);
        else if (a == "--body-file-segments" && i + 1 < argc) o.bodyFileSegments = std::atoll(argv[++i]);  // developer: yarn-body file split
        else if (a == "--taper-report" && i + 1 < argc) o.taperReport = argv[++i];
        else if (a == "--viewer-defaults") o.viewerSettings = true;          // start from the viewer's settings (Draft, this machine's devices)
        else if (a == "--fibres-full") o.sets.push_back({"fibres.generation", "full"});  // rev. 2 5c: every candidate
        else if (a == "--fibre-support-mb" && i + 1 < argc) o.fibreSupportMB = std::atof(argv[++i]);  // bounded support residency
        else if (a == "--no-self-intersection") o.noSelf = true;  // GLB load: skip the triangle-pair test
        else if (a == "--load-only") o.loadOnly = true;            // GLB: load and mesh checks only
        else if (a == "--cache-build") o.cacheBuild = true;
        else if (a == "--part-cap" && i + 1 < argc) o.partCapS = std::atof(argv[++i]);
        else if (a == "--cache" && i + 1 < argc) o.cacheRoot = argv[++i];
        else if (a == "-v" || a == "--verbose") o.verbose = true;
        else { std::fprintf(stderr, "unknown argument %s\n", a.c_str()); return 2; }
    }
    // The relax dumps carry the command line as their identity (the run's identity.txt has the rest).
    return kContinue;
}

// The scene's entry and its settings: the defaults, the settings file, --set; the refusals that need
// them. kContinue, or the exit code.
int resolveScene(const SceneCommandOptions& o, const SceneEntry*& e, SceneSettings& s) {
    // An import resolves its entry from the import.* settings (below); its defaults are any tiled
    // mesh's (D 1.10, the fallback chain, knitted).
    SceneEntry importDefaults;
    importDefaults.id = "import";
    importDefaults.kind = SceneEntry::Kind::Mesh;
    e = o.id == "import" ? &importDefaults : findScene(o.id);
    if (!e) { std::fprintf(stderr, "unknown scene %s (knit_check --scene list)\n", o.id.c_str()); return 2; }
    s.scene = o.id;
    applySceneDefaults(s, *e);  // GLB scenes: D 1.10 (the largest flat PASS at 5.5 x 3.4); no scene sets a workflow
    if (o.cacheBuild || o.viewerSettings) {
        s = viewerDefaults(*e);
        s.scene = o.id;
    }
    const SceneSettings defaults = s;
    std::string err;
    if (!o.settingsFile.empty()) {
        std::ifstream f(o.settingsFile);
        std::stringstream ss;
        ss << f.rdbuf();
        std::vector<std::string> retired;
        if (!f || !readSettings(ss.str(), s, &err, &retired)) { std::fprintf(stderr, "%s: %s\n", o.settingsFile.c_str(), err.c_str()); return 2; }
        for (const std::string& k : retired) std::fprintf(stderr, "%s: %s is no longer used (ignored)\n", o.settingsFile.c_str(), k.c_str());
        if (s.scene != o.id) { std::fprintf(stderr, "%s is for scene %s, not %s\n", o.settingsFile.c_str(), s.scene.c_str(), o.id.c_str()); return 2; }
    }
    for (const auto& [k, v] : o.sets)
        if (!setSetting(s, k, v, &err)) { std::fprintf(stderr, "--set: %s\n", err.c_str()); return 2; }
    if (o.id == "import") {
        e = sceneEntryFor(s);
        if (!e) { std::fprintf(stderr, "--scene import: set import.identity (16 hex digits) and import.source of a prepared import\n"); return 2; }
    }
    if (s.stitchSize == "yarn")
        for (const auto& [k, v] : o.sets)
            if (k == "strand.wale" || k == "strand.course")
                std::fprintf(stderr, "note: --set %s=%s is replaced by the stitch size from the yarn (strand.stitch_size = yarn): %.2f x %.2f mm\n",
                             k.c_str(), v.c_str(), s.wale, s.course);
    if (const std::string c = yarnLookConflict(s); !c.empty()) std::fprintf(stderr, "warning: %s\n", c.c_str());
    if (!o.writeFile.empty()) {
        std::ofstream f(o.writeFile);
        f << writeSettings(s);
        if (!f) { std::fprintf(stderr, "cannot write %s\n", o.writeFile.c_str()); return 2; }
    }
    const std::vector<std::string> changed = diffSettings(defaults, s);
    std::printf("scene %s (%s): %zu settings differ from the defaults%s", e->id.c_str(), e->label.c_str(), changed.size(),
                changed.empty() ? "\n" : ":");
    for (size_t k = 0; k < changed.size(); ++k) std::printf(" %s%s", changed[k].c_str(), k + 1 < changed.size() ? "," : "\n");
    if (furWorkflow(s) && !pileSettingsValid(s, &err)) {
        std::fprintf(stderr, "refused (workflow = fur): %s\n", err.c_str());
        return 2;
    }
    return kContinue;
}

// --cache-build: what the viewer shows, computed and stored in the scene cache.
int buildSceneCache(const SceneCommandOptions& o, const SceneSettings& s) {
    std::string err;
    KnitRequest job;
    job.settings = std::make_shared<const SceneSettings>(s);
    job.cacheRoot = o.cacheRoot;
    auto b = knitObject(job, JobControl{}, &err);
    if (!b) { std::fprintf(stderr, "%s\n", err.c_str()); return 2; }
    int pass = 0, fail = 0, unres = 0, plain = 0;
    for (const PartRow& r : b->parts) {
        if (!r.knitted) ++plain;
        else if (r.verdict == Status::Pass) ++pass;
        else if (r.verdict == Status::Fail) ++fail;
        else ++unres;
    }
    std::printf("%s: %zu parts (%d PASS, %d FAIL, %d UNRESOLVED, %d plain), %zu display meshes; %s in %.1f s -> %s\n", o.id.c_str(),
                b->parts.size(), pass, fail, unres, plain, b->meshes.size(), b->cached ? "already cached, read" : "computed and stored",
                b->ms / 1000, o.cacheRoot.c_str());
    if (b->report.find("scene cache not written") != std::string::npos) { std::fprintf(stderr, "%s", b->report.c_str()); return 2; }
    return 0;
}

// The scene's source as it loads, with its mesh checks. kContinue, or the exit code.
int reportLoadedScene(const SceneCommandOptions& o, const SceneEntry* e, const SceneSettings& s, bool spiral) {
    std::string err;
    if (spiral) {
        std::printf("scene %s: knitted as spiral pieces (%s)\n", o.id.c_str(), e->piecesManifest.c_str());
        if (o.loadOnly) {  // the pieces and kept parts as the scene's primitives, with their mesh checks
            Scene t;
            SpiralSource src;
            if (!loadPiecesScene(*e, s, src, t, &err)) { std::fprintf(stderr, "%s\n", err.c_str()); return 2; }
            std::printf("loaded %s: %zu primitives (%zu knitted pieces); bounds [%.1f %.1f %.1f] .. [%.1f %.1f %.1f] mm\n",
                        t.source.path.c_str(), t.source.prims.size(), t.parts.size(), t.source.lo.x, t.source.lo.y, t.source.lo.z, t.source.hi.x,
                        t.source.hi.y, t.source.hi.z);
            for (const SourcePrimitive& p : t.source.prims) std::printf("  [%s / %s]\n%s", p.mesh.c_str(), p.material.c_str(), p.check.summary().c_str());
            return 0;
        }
    } else if (e->kind == SceneEntry::Kind::Glb) {
        const std::string path = findProjectFile(e->path);
        if (path.empty()) { std::fprintf(stderr, "scene %s: %s not found\n", o.id.c_str(), e->path.c_str()); return 2; }
        SourceScene g;
        SourceOptions go;
        go.selfIntersection = !o.noSelf;
        if (!loadGlb(path, e->glb, g, &err, go)) { std::fprintf(stderr, "%s\n", err.c_str()); return 2; }
        std::printf("loaded %s (FNV-1a %016llx): %zu primitives; bounds [%.1f %.1f %.1f] .. [%.1f %.1f %.1f] mm\n", path.c_str(),
                    static_cast<unsigned long long>(g.fileHash), g.prims.size(), g.lo.x, g.lo.y, g.lo.z, g.hi.x, g.hi.y, g.hi.z);
        for (const std::string& l : g.log) std::printf("  %s\n", l.c_str());
        for (const SourcePrimitive& p : g.prims) std::printf("  [%s / %s]\n%s", p.mesh.c_str(), p.material.c_str(), p.check.summary().c_str());
    } else {
        std::printf("import %s (%s)\n", s.importSource.c_str(), e->path.c_str());
    }
        return o.loadOnly ? 0 : kContinue;
}

// --stitch-spread: each knitted part's spread, what settle moved, and the raised loops left standing.
void printStitchSpread(const Scene& t, const SceneSettings& s) {
    for (const ScenePart& p : t.parts) {
        if (!p.knit) continue;
        const StitchSpread sp = measureStitchSpread(p.knit->t.curves);
        const std::vector<int>& defects = p.knit->settle.defectsPerRound;
        std::printf("stitch spread %s: %d of %zu stitches, shape %.4f mm, height %.4f mm; settle %d rounds, first round %d stitches, max move %.3f mm\n",
                    p.name.c_str(), sp.stitches, p.knit->t.curves.topo.stitches.size(), sp.shape, sp.height, p.knit->settle.rounds,
                    defects.empty() ? 0 : defects.front(), p.knit->settle.maxMove);
        const CurveSet& curves = p.knit->t.curves;
        const RaisedLoops loops = measureRaisedLoops(curves, drawStitchTensions(curves.topo, handRequest(s, p.name)));
        if (loops.asked > 0) std::printf("raised loops %s: %d of %d standing\n", p.name.c_str(), loops.standing, loops.asked);
    }
}

// Knits the scene, reports it, and grows or loads its fibres and exports them as asked.
// One TaperReport as a JSON object (lengths and radii in metres, errors as fractions).
std::string taperReportJson(const TaperReport& r) {
    std::ostringstream j;
    j.precision(10);
    j << "{\"fibres\": " << r.fibres << ", \"tapered\": " << r.tapered << ", \"fallback_no_tail\": " << r.fallbackNoTail
      << ", \"fallback_rule\": {";
    for (int k = 1; k < kStrandFaultCount; ++k)
        j << (k > 1 ? ", " : "") << "\"" << strandFaultName(StrandFault(k)) << "\": " << r.fallbackRule[k];
    j << "}, \"onset_exact\": " << r.onsetExact << ", \"onset_shifted\": " << r.onsetShifted
      << ", \"onset_impossible\": " << r.onsetImpossible << ", \"nonconforming\": " << r.nonconforming
      << ", \"nonconforming_first\": [";
    for (size_t i = 0; i < r.nonconformingFirst.size(); ++i) j << (i ? ", " : "") << r.nonconformingFirst[i];
    const double areaError = r.tailAreaExact > 0 ? std::fabs(r.tailAreaLinear - r.tailAreaExact) / r.tailAreaExact : 0;
    j << "], \"original_segments\": " << r.originalSegments << ", \"final_segments\": " << r.finalSegments
      << ", \"max_error_over_r0\": " << r.maxErrorOverR0 << ", \"max_error_over_local_r_sampled\": " << r.maxErrorOverLocalR
      << ", \"max_early_thinning_over_r0\": " << r.maxEarlyThinningOverR0 << ", \"tail_area_error\": " << areaError
      << ", \"min_radius_m\": " << r.minRadius << ", \"max_radius_m\": " << r.maxRadius
      << ", \"max_stored_radius_ratio\": " << r.maxStoredRatio << "}";
    return j.str();
}

// knit_check --taper-report: the late taper of every part's exported fibres, measured on the strands as
// written (measureTaper), per part and in total. The fibre indices of nonconforming strands are the
// part's.
bool writeTaperReport(const std::string& path, const Scene& t, const SceneSettings& s, const SceneFibres& fibres,
                      const SceneCommandOptions& o, std::string* err) {
    StrandTaper taper{s.tailTaper, s.tailTaperMinDiameters, s.tailTaperTip, TaperPoints::Dense, o.taperTolerance};
    if (!o.taperPoints.empty()) parseTaperPoints(o.taperPoints, taper.points);
    std::ostringstream j;
    j << "{\"points\": \"" << taperPointsName(taper.points) << "\", \"tolerance\": " << taper.tolerance
      << ", \"tail\": " << taper.tail << ", \"min_diameters\": " << taper.minDiameters << ", \"tip\": " << taper.tip
      << ", \"parts\": {";
    TaperReport total;
    for (size_t i = 0; i < fibres.parts.size(); ++i) {
        const PartFibres& pf = fibres.parts[i];
        const TaperReport r = partTaperReport(pf, t.parts[static_cast<size_t>(pf.part)].centre, taper);
        addTaperReport(total, r);
        j << (i ? ", " : "") << "\"" << pf.name << "\": " << taperReportJson(r);
    }
    j << "}, \"total\": " << taperReportJson(total) << "}\n";
    std::ofstream f(path, std::ios::binary);
    f << j.str();
    if (!f) {
        if (err) *err = "cannot write " + path;
        return false;
    }
    return true;
}

// A scene knitted on its mesh (an import, a plane) kept as one file: knit.bin (saved_knit.h) and
// knit.key, the geometry key of the settings it was knitted with. A loaded knit is the one that
// was reviewed: knitting is not run again (the solver is not deterministic), and gate (a) is
// verified on the stored curves (verifyKnit).
const char* kKnitKeyFile = "knit.key";

bool saveTiledKnit(const std::string& dir, const SceneSettings& s, const Scene& t, std::string* err) {
    std::error_code ec;
    std::filesystem::create_directories(dir, ec);
    std::ofstream bin(dir + "/knit.bin", std::ios::binary), key(dir + "/" + kKnitKeyFile, std::ios::binary);
    bin << serializeKnit(t);
    key << std::hex << geometryKey(s) << "\n";
    if (!bin || !key) {
        if (err) *err = "--save-knit: cannot write " + dir;
        return false;
    }
    return true;
}

bool loadTiledKnit(const std::string& dir, const SceneEntry& e, const SceneSettings& s, Scene& t, std::string* err) {
    std::ifstream key(dir + "/" + kKnitKeyFile, std::ios::binary), bin(dir + "/knit.bin", std::ios::binary);
    uint64_t stored = 0;
    if (!(key >> std::hex >> stored) || !bin) {
        if (err) *err = "--load-knit: " + dir + " has no knit.bin and knit.key";
        return false;
    }
    if (stored != geometryKey(s)) {
        if (err) *err = "--load-knit: " + dir + " was knitted with other knitting settings (a Load, Knit or Surface setting differs)";
        return false;
    }
    std::stringstream bytes;
    bytes << bin.rdbuf();
    return restoreKnit(bytes.str(), e, s, t, err) && verifyKnit(t, s, err);
}

int knitAndExport(const SceneCommandOptions& o, const SceneEntry* e, const SceneSettings& s) {
    std::string err;
    if (o.threads > 0) setThreadBudget(o.threads);
    Scene t;
    KnitRunOptions ro;
    ro.verbose = o.verbose;
    ro.source.selfIntersection = !o.noSelf;
    ro.stageProfile = o.stageProfile;
    ro.capS = s.timeCapS;
    ro.partCapS = o.partCapS;
    ro.only = o.onlyPart;
    // Spiral pieces: the per-piece knit cache (knitScene); a scene knitted on its mesh: one knit file here.
    const bool tiledFile = knitWorkflow(s) && !knitsSewnPieces(*e, s);
    if (!tiledFile) {
        ro.saveKnitCache = o.saveKnit;
        ro.loadKnitCache = o.loadKnit;
    }
    const StageClock spRun;
    // Fur: nothing is knitted, the pile grows on the mesh (as the viewer).
    bool loaded = false;
    if (tiledFile && !o.loadKnit.empty()) {
        loaded = loadTiledKnit(o.loadKnit, *e, s, t, &err);
        if (loaded) std::printf("knit: loaded from %s (not knitted again; gate (a) verified on the stored curves)\n", o.loadKnit.c_str());
    } else {
        loaded = knitWorkflow(s) ? knitScene(*e, s, ro, t, &err) : loadBaseParts(*e, s, ro.source, t, &err);
    }
    if (!loaded) { std::fprintf(stderr, "%s\n", err.c_str()); return 2; }
    if (tiledFile && !o.saveKnit.empty()) {
        if (!saveTiledKnit(o.saveKnit, s, t, &err)) { std::fprintf(stderr, "%s\n", err.c_str()); return 2; }
        std::printf("knit: saved to %s\n", o.saveKnit.c_str());
    }
    std::printf("%s", formatKnitReport(t, s, o.verbose).c_str());
    if (o.stageProfile) std::printf("profile: %s\n", spRun.line("knit (all parts)").c_str());
    if (o.stitchSpread) printStitchSpread(t, s);
    // Incomplete (exit 1, after the export) is not a failed check: a selected part lacks its knit.
    std::vector<std::string> incomplete;
    for (const ScenePart& r : t.parts)
        if (partIncomplete(r, o.onlyPart)) incomplete.push_back(r.name + ": " + r.note);
    if (!incomplete.empty()) {
        std::printf("incomplete: %zu selected part(s) not knitted (exit code 1)\n", incomplete.size());
        for (const std::string& l : incomplete) std::printf("  %s\n", l.c_str());
    }
    const bool exportFibres = !o.exportDir.empty() && (s.fibres || furWorkflow(s));  // Fur: the pile is the export's content (as the viewer)
    // The fibres, grown as the viewer grows them (with knit_check's diagnostics) or loaded.
    SceneFibres fibres;
    if (!o.loadFibres.empty()) {
        if (!loadSceneFibres(o.loadFibres, t, s, fibres, &err)) { std::fprintf(stderr, "fibres: %s\n", err.c_str()); return 2; }
        std::printf("%s", fibres.log.c_str());
    } else if (exportFibres || !o.saveFibres.empty()) {
        {
            FibreRunOptions fo;
            fo.fibreCheck = o.fibreCheck;
            fo.stageProfile = o.stageProfile;
            fo.fibreDigest = o.fibreDigest;
            fo.fibreSupportMB = o.fibreSupportMB;
            const StageClock spFibres;
            std::optional<SceneFibres> grown = growSceneFibres(*e, t, s, fo, &err);
            if (!grown) { std::fprintf(stderr, "fibres: %s\n", err.c_str()); return 2; }
            fibres = std::move(*grown);
            if (o.stageProfile) std::printf("profile: %s\n", spFibres.line("fibres (all parts)").c_str());
        }
        if (!o.saveFibres.empty() && !saveSceneFibres(o.saveFibres, t, s, fibres, &err)) { std::fprintf(stderr, "fibres: %s\n", err.c_str()); return 2; }
    }
    if (!o.exportDir.empty()) {
        ExportOptions eo;
        eo.format = o.exportFormat;
        eo.folder = o.exportDir;
        eo.stageProfile = o.stageProfile;
        eo.closeupPart = o.closeup;
        eo.taperPoints = o.taperPoints;
        eo.taperTolerance = o.taperTolerance;
        eo.bodyFileSegments = o.bodyFileSegments;
        std::string log;
        const StageClock spExport;
        const bool ok = exportScene({*e, t, s, exportFibres ? &fibres : nullptr}, eo, &log, &err);
        std::printf("%s", log.c_str());
        if (o.stageProfile) std::printf("profile: %s\n", spExport.line("export (all parts)").c_str());
        if (!ok) { std::fprintf(stderr, "export: %s\n", err.c_str()); return 2; }
        if (!o.taperReport.empty() && exportFibres && !writeTaperReport(o.taperReport, t, s, fibres, o, &err)) {
            std::fprintf(stderr, "taper report: %s\n", err.c_str());
            return 2;
        }
    }
    return incomplete.empty() ? 0 : 1;
}

} // namespace

int sceneCommand(int argc, char** argv) {
    SceneCommandOptions o;
    if (const int exitCode = parseSceneOptions(argc, argv, o); exitCode != kContinue) return exitCode;
    if (o.id == "list" || o.id.empty()) {
        for (const SceneEntry& e : sceneRegistry())
            std::printf("%-17s %-38s %s%s\n", e.id.c_str(), e.label.c_str(), e.path.c_str(),
                        findProjectFile(e.path).empty() ? "  (missing)" : "");
        return o.id.empty() ? 2 : 0;
    }
    const SceneEntry* e = nullptr;
    SceneSettings s;
    if (const int exitCode = resolveScene(o, e, s); exitCode != kContinue) return exitCode;
    if (o.cacheBuild) return buildSceneCache(o, s);
    // Every scene runs the common path below (knitScene: a knitted scene of spiral pieces, the
    // teddy, is knitted as spirals there; growSceneFibres grows the fibres, exportScene writes them).
    const bool spiral = knit::knitsSewnPieces(*e, s);
    if (furWorkflow(s) && (!o.saveKnit.empty() || !o.loadKnit.empty())) {
        std::fprintf(stderr, "--save-knit / --load-knit: the Fur workflow (workflow = fur) has no knit\n");
        return 2;
    }
    const ExportFormat* format = findExportFormat(o.exportFormat);
    if (!format || !format->available) {
        std::fprintf(stderr, "--export-format %s: unknown or unavailable\n", o.exportFormat.c_str());
        return 2;
    }
    if (o.closeup >= 0 && o.exportDir.empty()) {
        std::fprintf(stderr, "--closeup N: the close-up cameras are part of the export (--export-ppt dir)\n");
        return 2;
    }
    const StageClock reportClock;
    if (const int exitCode = reportLoadedScene(o, e, s, spiral); exitCode != kContinue) return exitCode;
    if (o.stageProfile) std::printf("profile: %s\n", reportClock.line("source report").c_str());
    return knitAndExport(o, e, s);
}
