// knit_check: command-line validation of stockinette patches.
//
//   knit_check --all                       validate every preset, exit 1 on unexpected verdict
//   knit_check --preset <name> [-v]        validate one preset (verbose lists issues)
//   knit_check --flat-dmax [...]           the largest D whose flat 8 x 8 swatch passes gate (a)
//   knit_check --ply-map-sets [...]        S6: the default texture map sets in the cache (ply_map_sets.h)
//   knit_check --retaper <in> <out> <tail> [<min diameters> <tip>]
//                                          a photon_path_tracer package's fibres with a late taper, the
//                                          knit frozen (fibres.tail_taper; fibre_strands.h taperStrandSet)
//   Overrides: --w <mm> --h <mm> --D <mm> --c <mm> --A <xD> --eps <mm> --head <x> --leg <x> --n <count>
//              --pattern <name> --trans <neckY neckZ legY> --rib-layer <xD> --rib-spread <w> --garter-wave <xD>
//              --local  (also the per-stitch signature with crossing attribution, as on a surface)
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <vector>

#include <filesystem>
#include <fstream>
#include <sstream>

#include "fibres/fibre_strands.h"
#include "knitting/checks/checks.h"
#include "knitting/knit/coverage.h"
#include "knitting/tube/yarn_surface.h"
#include "knitting/knit/fixtures.h"
#include "knitting/stitch/stockinette.h"
#include "jobs/crash_report.h"
#include "settings/scene_registry.h"

using namespace knit;

namespace {

ValidationReport run(const PatternParams& pp, CurveSet* keep = nullptr) {
    CurveSet cs = generateStockinette(pp);
    ValidationReport r = validate(cs, pp);
    if (keep) *keep = std::move(cs);
    return r;
}

// The largest D whose flat 8 x 8 swatch passes gate (a) at the given
// gauge (scan upward in dStep from dStart; every verdict is printed, a non-monotone sequence is
// visible), with the front open area of the round tube at that D (coverage.h meshOpenArea) and the
// largest connected hole width (coverage.h regions on the interior window).
struct FlatRow { double D = 0; Status verdict = Status::Unresolved; double open = -1, holeWidth = -1; };
FlatRow flatRow(const PatternParams& pp, bool measure) {
    FlatRow r;
    r.D = pp.yarnDiameter;
    const CurveSet cs = generateStockinette(pp);
    r.verdict = validate(cs, pp).verdict;
    if (!measure || r.verdict != Status::Pass) return r;
    const SurfaceParams sp{};
    const SurfaceModel sm = buildSurface(cs, pp, sp);
    TubeParams tp;
    tp.ringSpacing = 0.05;
    tp.sides = 16;
    const TubeMesh tube = buildTubeMesh(cs, sm, tp);
    const InteriorWindow w = interiorWindow(cs);
    r.open = meshOpenArea(tube, w, 600);
    TriMesh wm;
    wm.pos = {{w.x0, w.y0, 0}, {w.x1, w.y0, 0}, {w.x1, w.y1, 0}, {w.x0, w.y1, 0}};
    wm.tri = {{0, 1, 2}, {0, 2, 3}};
    wm.computeNormals();
    CoverageInput in;
    in.mesh = &wm;
    in.own.tube = &tube;
    CoverageParams cp;
    cp.spacing = 0.04;
    cp.reach = pp.centrelineDepth() + pp.yarnDiameter + 1;
    r.holeWidth = measureCoverage(in, cp, ExecPolicy{true}).holes.maxWidth;
    return r;
}

void flatDmax(PatternParams pp, double dStart, double dStep) {
    pp.wales = pp.courses = 8;
    std::printf("flat D_max at %.3f x %.3f mm (c %.3f, A %.2f D), 8 x 8 swatch; D from %.3f step %.3f\n", pp.waleSpacing, pp.courseSpacing,
                pp.clearance, pp.depthAmplitude, dStart, dStep);
    double dMax = -1;
    int failsInRow = 0;
    for (double D = dStart; failsInRow < 3 && D < 5; D += dStep) {
        pp.yarnDiameter = D;
        const FlatRow r = flatRow(pp, false);
        std::printf("  D %.3f: %s\n", D, toString(r.verdict));
        if (r.verdict == Status::Pass) {
            if (failsInRow == 0) dMax = D;
            else std::printf("  (non-monotone: PASS again at D %.3f after a non-PASS)\n", D);
        } else {
            ++failsInRow;
        }
    }
    if (dMax < 0) {
        std::printf("D_max: none (no PASS from D %.3f)\n", dStart);
        return;
    }
    pp.yarnDiameter = dMax;
    const FlatRow r = flatRow(pp, true);
    std::printf("D_max %.3f mm (D/w %.3f): front open area %.4f, max connected hole width %.3f mm (%.3f w)\n", dMax, dMax / pp.waleSpacing,
                r.open, r.holeWidth, r.holeWidth / pp.waleSpacing);
}

// Flags of removed features: refused with what became of them. Empty for any other argument.
std::string retiredFlagMessage(const std::string& a) {
    if (a == "--contract" || a == "--k" || a == "--core-c") return a + " (the compress contact) was retired: yarns keep their clearance";
    if (a == "--relax-cell" || a == "--relaxed" || a == "--variant" || a == "--tau" || a == "--cell-time")
        return a + " (the relaxed stitch cell) was retired: stitches keep the fixed template";
    return "";
}

} // namespace

int shapingDemo(int argc, char** argv);    // src/cli/pieces_cli.cpp
int strandDemo(int argc, char** argv);     // src/cli/pieces_cli.cpp
int plyMapSetsCli(int argc, char** argv);  // src/cli/pieces_cli.cpp
int sceneCommand(int argc, char** argv);   // src/cli/scene_cli.cpp (S6)

// knit_check --pieces <manifest> [--gauge <piece file> <stitches>] [the --scene options]: the manifest
// as a scene (knit::manifestScene), then the one scene command.
int piecesCommand(int argc, char** argv) {
    std::vector<std::string> args{argv[0]};
    std::string manifest, gaugeFile;
    double gaugeStitches = 0;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--pieces" && i + 1 < argc) manifest = argv[++i];
        else if (a == "--gauge" && i + 2 < argc) {
            gaugeFile = argv[++i];
            gaugeStitches = std::atof(argv[++i]);
        } else args.push_back(a);
    }
    std::string err;
    const knit::SceneEntry* e = manifest.empty() ? nullptr : knit::manifestScene(manifest, gaugeFile, gaugeStitches, &err);
    if (!e) {
        std::fprintf(stderr, "knit_check --pieces: %s\n", err.empty() ? "a manifest file is required" : err.c_str());
        return 2;
    }
    args.push_back("--scene");
    args.push_back(e->id);
    std::vector<char*> sceneArgv;
    for (std::string& a : args) sceneArgv.push_back(a.data());
    return sceneCommand(static_cast<int>(sceneArgv.size()), sceneArgv.data());
}

// knit_check --retaper: a copy of a package with every strand file late-tapered (taperStrandSet) and the
// renderer's whole-fibre taper switched off (tip_scale 1, as an export with fibres.tail_taper writes it),
// so a taper is compared on one frozen knit instead of a new, different one. The exported radii are
// the untapered export's; a strand whose taper would break a strand rule is copied as it was.
int retaperCli(int argc, char** argv) {
    namespace fs = std::filesystem;
    int at = 0;
    for (int i = 1; i < argc; ++i)
        if (std::string(argv[i]) == "--retaper") at = i;
    if (at + 3 >= argc) {
        std::fprintf(stderr, "usage: knit_check --retaper <package> <out> <tail> [<min diameters> <tip>]\n");
        return 2;
    }
    const fs::path in = argv[at + 1], out = argv[at + 2];
    knit::StrandTaper taper;
    taper.tail = std::atof(argv[at + 3]);
    if (at + 5 < argc) {
        taper.minDiameters = std::atof(argv[at + 4]);
        taper.tip = std::atof(argv[at + 5]);
    }
    if (!(taper.tail > 0 && taper.tail <= 1) || !(taper.tip > 0 && taper.tip <= 1) || taper.minDiameters < 0) {
        std::fprintf(stderr, "--retaper: tail above 0 up to 1, tip above 0 up to 1, min diameters 0 or more\n");
        return 2;
    }
    std::error_code ec;
    if (!fs::exists(in / "scene.json") || fs::exists(out)) {
        std::fprintf(stderr, "--retaper: %s must be a package and %s must not exist\n", in.string().c_str(), out.string().c_str());
        return 2;
    }
    fs::copy(in, out, fs::copy_options::recursive, ec);
    if (ec) {
        std::fprintf(stderr, "--retaper: cannot copy %s: %s\n", in.string().c_str(), ec.message().c_str());
        return 2;
    }
    long long tapered = 0, kept = 0;
    for (const fs::directory_entry& e : fs::directory_iterator(out / "strands")) {
        if (e.path().extension() != ".pptstrands") continue;
        knit::StrandSetF s;
        std::string err;
        if (!knit::readStrandsFile(e.path().string(), s, &err)) {
            std::fprintf(stderr, "--retaper: %s\n", err.c_str());
            return 2;
        }
        const knit::TaperedStrands t = knit::taperStrandSet(s, taper);
        if (!knit::writeStrandsFile(e.path().string(), t.strands, knit::StrandRules{}, &err)) {
            std::fprintf(stderr, "--retaper: %s\n", err.c_str());
            return 2;
        }
        tapered += t.tapered;
        kept += t.untapered;
    }
    // The renderer's tip_scale (format_ppt.cpp writeInstances): 1 under a late taper.
    std::ifstream f(out / "scene.json");
    std::stringstream text;
    text << f.rdbuf();
    f.close();
    std::string scene = text.str();
    const std::string from = "\"tip_scale\": 0.2", to = "\"tip_scale\": 1";
    int shapes = 0;
    for (size_t pos = scene.find(from); pos != std::string::npos; pos = scene.find(from, pos + to.size())) {
        scene.replace(pos, from.size(), to);
        ++shapes;
    }
    std::ofstream(out / "scene.json", std::ios::binary) << scene;
    std::printf("retaper: tail %.3f, at least %.1f diameters, tip %.2f: %lld fibres tapered, %lld kept their taper; "
                "%d fibre shapes set to tip_scale 1 -> %s\n", taper.tail, taper.minDiameters, taper.tip, tapered, kept, shapes,
                out.string().c_str());
    return 0;
}

int main(int argc, char** argv) {
    // Unbuffered when redirected (seamless plan M0; MSVC has no line buffering): a run killed at its
    // time cap keeps its log.
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    knit::installCrashReport("");  // stderr
    for (int i = 1; i < argc; ++i)
        if (std::string(argv[i]) == "--ply-map-sets") return plyMapSetsCli(argc, argv);
    for (int i = 1; i < argc; ++i)
        if (std::string(argv[i]) == "--scene") return sceneCommand(argc, argv);
    for (int i = 1; i < argc; ++i)
        if (std::string(argv[i]) == "--pieces") return piecesCommand(argc, argv);
    for (int i = 1; i < argc; ++i)
        if (std::string(argv[i]) == "--shaping-demo") return shapingDemo(argc, argv);
    for (int i = 1; i < argc; ++i)
        if (std::string(argv[i]) == "--strand-demo") return strandDemo(argc, argv);
    for (int i = 1; i < argc; ++i)
        if (std::string(argv[i]) == "--retaper") return retaperCli(argc, argv);
    for (int i = 1; i < argc; ++i)
        if (std::string(argv[i]) == "--relax-frames") {
            std::fprintf(stderr, "--relax-frames (centrelines of relax-call dumps) was retired with the dumps\n");
            return 2;
        }
    std::string preset = "expanded";
    bool all = false, verbose = false, flatD = false, local = false;
    double dStart = 0.5, dStep = 0.02;
    PatternParams pp = findPreset(preset)->params;
    bool overridden = false;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (const std::string retired = retiredFlagMessage(a); !retired.empty()) {
            std::fprintf(stderr, "%s\n", retired.c_str());
            return 2;
        }
        auto num = [&]() { return i + 1 < argc ? std::atof(argv[++i]) : 0.0; };
        if (a == "--all") all = true;
        else if (a == "--flat-dmax") flatD = true;
        else if (a == "--local") local = true;  // also the per-stitch signature (as on a surface), with its failures
        else if (a == "--d-start") dStart = num();
        else if (a == "--d-step") dStep = num();
        else if (a == "-v" || a == "--verbose") verbose = true;
        else if (a == "--preset" && i + 1 < argc) {
            preset = argv[++i];
            const Preset* p = findPreset(preset);
            if (!p) { std::fprintf(stderr, "unknown preset %s\n", preset.c_str()); return 2; }
            pp = p->params;
        }
        else if (a == "--w") { pp.waleSpacing = num(); overridden = true; }
        else if (a == "--h") { pp.courseSpacing = num(); overridden = true; }
        else if (a == "--D") { pp.yarnDiameter = num(); overridden = true; }
        else if (a == "--c") { pp.clearance = num(); overridden = true; }
        else if (a == "--A") { pp.depthAmplitude = num(); overridden = true; }
        else if (a == "--eps") { pp.sampleTolerance = num(); overridden = true; }
        else if (a == "--head") { pp.headHalfWidth = num(); overridden = true; }
        else if (a == "--leg") { pp.legHalfWidth = num(); overridden = true; }
        else if (a == "--n") { pp.wales = pp.courses = static_cast<int>(num()); overridden = true; }
        else if (a == "--trans" && i + 3 < argc) {  // S6 garter transition: neckY neckZ legY
            pp.transNeckY = num();
            pp.transNeckZ = num();
            pp.transLegY = num();
            overridden = true;
        }
        else if (a == "--garter-wave") { pp.garterWave = num(); overridden = true; }  // garter ridges [D]
        else if (a == "--rib-spread") { pp.ribPairSpread = num(); overridden = true; }  // rib 2x2 pair spread [w]
        else if (a == "--rib-layer") { pp.ribLayerOffset = num(); overridden = true; }  // rib layers [D]
        else if (a == "--pattern" && i + 1 < argc) {  // S6: stockinette | reverse-stockinette | garter | rib1x1 | rib2x2
            const std::optional<StitchPattern> pattern = parsePattern(argv[++i]);
            if (!pattern) { std::fprintf(stderr, "unknown pattern %s\n", argv[i]); return 2; }
            pp.pattern = *pattern;
            overridden = true;
        }
        else { std::fprintf(stderr, "unknown argument %s\n", a.c_str()); return 2; }
    }
    if (flatD) { flatDmax(pp, dStart, dStep); return 0; }
    if (all) {
        int bad = 0;
        for (const Preset& p : presets()) {
            CurveSet cs;
            ValidationReport r = run(p.params, &cs);
            bool asExpected = r.verdict == p.expected;
            bad += !asExpected;
            std::printf("=== preset '%s' (expected %s) -> %s %s\n    %s\n", p.name.c_str(), toString(p.expected),
                        toString(r.verdict), asExpected ? "[as expected]" : "[UNEXPECTED]", p.description.c_str());
            std::printf("%s\n", formatReport(r, cs, p.params, verbose).c_str());
        }
        return bad ? 1 : 0;
    }
    CurveSet cs;
    ValidationReport r = run(pp, &cs);
    std::printf("=== %s%s\n%s", preset.c_str(), overridden ? " (with overrides)" : "",
                formatReport(r, cs, pp, verbose).c_str());
    if (local) {
        // The flat patch's stitches checked as a surface knit checks them: each in its own frame, with
        // crossing attribution (the window of the scene knits).
        const LocalSignatureResult ls = checkLocalSignature(cs, pp, ExecPolicy{true}, ValidateOptions{}.attributionWindow);
        std::printf("  local signature: %d PASS, %d FAIL, %d UNRESOLVED; %s\n", ls.pass, ls.fail, ls.unresolved, ls.detail.c_str());
        for (const LocalPairResult& lp : ls.pairs)
            if (lp.status != Status::Pass)
                std::printf("    stitch %d (parent %d): %s, %d crossings, %d unexpected: %s\n", lp.child, lp.parent, toString(lp.status),
                            lp.crossings, lp.unexpected, lp.detail.c_str());
    }
    return r.verdict == Status::Pass ? 0 : 1;
}
