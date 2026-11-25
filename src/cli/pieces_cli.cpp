// knit_check's research commands around the pieces' yarn: --strand-demo (one fibred strand as a
// photon_path_tracer package) and --shaping-demo (INC / DEC / cone patches). knit_check --pieces, the pieces with their own
// fibres and export, was retired in S2 (the teddy scene knits them: src/scene/knit_spiral_part.h).
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <atomic>
#include <memory>
#include <mutex>
#include <thread>
#include <string>
#include <vector>

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#include <psapi.h>
#pragma comment(lib, "psapi.lib")
#endif

#include "model/choices.h"
#include "model/coordinates.h"
#include "knitting/layout/chart.h"
#include "knitting/knit/fixtures.h"
#include "fibres/fibre_generator.h"
#include "knitting/tube/plies.h"
#include "knitting/layout/stitch_layout.h"
#include "knitting/knit/piece_knit.h"
#include "knitting/knit/pieces.h"
#include "knitting/relax/settle.h"
#include "knitting/tube/yarn_surface.h"
#include "scene/appearance/ply_maps.h"
#include "export/collect_scene.h"

#include "fibres/fibre_strands.h"
#include "scene/knit/knit_report.h"

using namespace knit;


// Consumes argv[i] (and its value) if it is one of PlyModeOptions' options.
bool parsePlyOption(PlyModeOptions& o, int& i, int argc, char** argv) {
    const std::string a = argv[i];
    const bool hasValue = i + 1 < argc;
    if (!hasValue) return false;
    const std::string v = argv[i + 1];
    if (a == "--ply-mode") o.mode = v == "geometry" ? PlyModeOptions::Mode::Geometry : v == "none" ? PlyModeOptions::Mode::None : PlyModeOptions::Mode::Texture;
    else if (a == "--ply-shading") o.shading = v == "bump" ? PlyModeOptions::Shading::Bump : v == "normal-flipg" ? PlyModeOptions::Shading::NormalFlipGreen : PlyModeOptions::Shading::Normal;
    else if (a == "--ply-cavity") o.cavity = v == "on";
    else if (a == "--ply-cavity-k") o.spec.cavity = std::atof(v.c_str());
    else if (a == "--ply-map-d") o.spec.yarnD = std::atof(v.c_str());
    else if (a == "--ply-map-size") o.spec.size = std::max(8, std::atoi(v.c_str()));
    else return false;
    ++i;
    return true;
}

// knit_check --strand-demo <dir> [--plies N] [--ply-pitch mm] [--ply-hand Z|S] [--ply-fill f]
//            [--ply-mode texture|geometry|none] [S5 texture options, ply_maps.h]
//            [--fiber-random f] [--fiber-density d] [--fiber-density-per-ply] [--fiber-material yarn|fibre]
//            [--yarn-d mm] [--length mm] [--bend deg] [--ply-sides N] [--ply-ring-spacing mm]
// S4 close-ups: one straight (plied) strand with both fibre populations, generated and
// checked as the pieces export does (shared setup, float32 quantization, serialized root
// check), exported as a photon_path_tracer scene with fixed "closeup" / "medium" cameras.
// S5: --bend replaces the straight strand by an arc turning by that angle (the curved
// fixture of the matched texture-vs-geometry renders); --tube-sides N sets the yarn tube's
// sides (default kExportSides), --ply-sides / --ply-ring-spacing the explicit ply tubes' (default
// kExportPlySides / kExportPlyRing); --sky s scales the uniform sky; --no-ground omits the ground quad.
int strandDemo(int argc, char** argv) {
    std::string dir;
    PlyGeometry plyGeom;
    PlyModeOptions plyOpt;
    FibreParams fp = pieceFibres();
    bool fibreMaterial = false, ground = true;
    double yarnD = 1.16, length = 12.0, bendDeg = 0, sky = 1.0;
    int tubeSides = 0;  // --tube-sides: sides of the (non-ply) yarn tube (0: kExportSides)
    int plySides = 0;   // --ply-sides: sides of the explicit ply tubes (0: kExportPlySides)
    double plyRing = 0; // --ply-ring-spacing: their ring spacing [mm] (0: kExportPlyRing)
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto num = [&]() { return i + 1 < argc ? std::atof(argv[++i]) : 0.0; };
        if (parsePlyOption(plyOpt, i, argc, argv)) continue;
        if (a == "--strand-demo" && i + 1 < argc) dir = argv[++i];
        else if (a == "--bend") bendDeg = num();
        else if (a == "--sky") sky = num();
        else if (a == "--no-ground") ground = false;
        else if (a == "--tube-sides") tubeSides = static_cast<int>(num());
        else if (a == "--ply-sides") plySides = static_cast<int>(num());
        else if (a == "--ply-ring-spacing") plyRing = num();
        else if (a == "--light") {
            std::fprintf(stderr, "--light was retired: the export writes no area lights; --sky scales the uniform sky\n");
            return 2;
        }
        else if (a == "--plies") plyGeom.count = std::max(1, static_cast<int>(num()));
        else if (a == "--ply-pitch") plyGeom.pitch = num();
        else if (a == "--ply-hand" && i + 1 < argc) plyGeom.handedness = std::string(argv[++i]) == "S" ? -1 : 1;
        else if (a == "--ply-fill") plyGeom.fill = num();
        else if (a == "--fiber-random") fp.randomFraction = num();
        else if (a == "--fiber-density") fp.density = num();
        else if (a == "--fiber-density-per-ply") fp.densityPerParent = false;
        else if (a == "--fiber-material" && i + 1 < argc) fibreMaterial = std::string(argv[++i]) == "fibre";
        else if (a == "--yarn-d") yarnD = num();
        else if (a == "--length") length = num();
        else { std::fprintf(stderr, "unknown argument %s\n", a.c_str()); return 2; }
    }
    if (dir.empty()) { std::fprintf(stderr, "--strand-demo needs an output directory\n"); return 2; }
    std::filesystem::create_directories(dir + "/meshes");
    StrandDemo d = makeStrandDemo(length, yarnD);
    if (bendDeg != 0) {
        // Control points moved onto an arc of radius length / bend through the origin in
        // the x-y plane (the fibre front +z is kept); the B-spline follows it closely.
        const double Rb = length / (bendDeg * kPi / 180);
        YarnPath& y = d.curves.yarns[0];
        for (Vec3& c : y.curve.ctrl) c = {Rb * std::sin(c.x / Rb), Rb * (1 - std::cos(c.x / Rb)), 0};
        y.finalize();
    }
    const PatternParams& pp = d.pp;
    SurfaceParams spY;
    spY.compression = false;
    const SurfaceModel parentSm = buildSurface(d.curves, pp, spY);
    const bool plied = plyOpt.mode == PlyModeOptions::Mode::Geometry && plyGeom.count > 1;
    std::unique_ptr<PlySet> plies;
    if (plied) plies = std::make_unique<PlySet>(buildPlies(d.curves, parentSm, plyGeom, pp.radius()));
    const CurveSet& fcs = plied ? plies->curves : d.curves;
    PatternParams fpp = pp;
    if (plied) fpp.yarnDiameter = 2 * plies->layout.rp;
    ExportScene exp(dir);
    exp.fibreLook.hairMaterials = fibreMaterial;
    TubeStyle style;
    if (tubeSides >= 3) style.sides = tubeSides;
    const double ring = !plied ? style.ring : plyRing > 0 ? plyRing : kExportPlyRing;
    const int sides = !plied ? style.sides : plySides >= 3 ? plySides : kExportPlySides;
    SceneMapping map;
    map.sigma = 1000;  // mm per scene unit: the package in true metres, as the bears' and the tiled scenes'
    double kappaMax = 0;
    for (const YarnPath& y : fcs.yarns) {
        const CurvatureResult cr = checkCurvature(y, fpp.radius(), Allowances::forExtent(fcs.extent()));
        kappaMax = std::max(kappaMax, std::min(cr.kappaHi, 1e3));
    }
    FibreFrame fr;
    FibreSetup fset;
    fset.tubeRadius = fpp.radius();
    fset.sides = sides;
    fset.ring = ring;
    fset.kappaMax = kappaMax;
    if (plied) {
        fset.plies = plies.get();
        fset.parentCurves = &d.curves;
        fset.parentSurface = &parentSm;
    }
    setupPieceFibres(fp, fr, pp, fset);
    fr.fabric = [](const Vec3& p) { return FabricPoint{{p.x, p.y, 0}, {0, 0, 1}}; };
    const double maxCoord = (0.5 * length + 2 * fp.lengthMax) / map.sigma;
    const double eSer = 2 * maxCoord * 1.1920929e-7 * map.sigma;
    fp.serializedMargin = eSer;
    fp.exportChordFactor = ring * ring / 8;
    SurfaceParams sp;
    sp.compression = false;
    const SurfaceModel sm = buildSurface(fcs, fpp, sp);
    const FibreSet fs = growYarnFibres(fcs, sm, fpp, fp, &fr);
    const FibreReport& r = fs.rep;
    plyOpt.spec.g = plyGeom;
    std::printf("strand %.1f mm (bend %.0f deg), D %.2f mm, %d plies (%s); density %.1f/mm (%s); random fraction %.2f\n", length,
                bendDeg, yarnD, plyGeom.count,
                plied ? "geometry" : plyOpt.textured() ? "texture" : "plain tube", fp.density,
                fp.densityPerParent ? "per mm of yarn" : "per mm of strand", fp.randomFraction);
    for (int pi = 0; pi < 2 && r.pop.size() == 2; ++pi)
        std::printf("  %s: %d candidates -> %d exposed -> %d accepted (%d truncated)\n", pi ? "random" : "directed",
                    r.pop[pi].candidates, r.pop[pi].exposed, r.pop[pi].accepted, r.pop[pi].truncated);
    std::string err;
    const std::optional<PlyShading> shading = setupPlyTextures(plyOpt, exp.folder(), map.sigma, &err);
    if (!shading) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 2;
    }
    exp.plyShading = *shading;
    style.uvRepeat = shading->uvRepeat;
    style.uvWrap = shading->uvWrap;
    const std::string tubeMat = plyOpt.textured() ? ExportScene::yarnMaterial("wool_white") : "wool_white";
    style.ring = ring;
    style.sides = sides;
    const std::optional<ExportTube> tube = buildExportTube(fcs, fpp, map, style, &err);
    if (!tube || !exp.addMeshPart("strand", tubeMat, tube->mesh, &err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 2;
    }
    if (!exp.addStrands("strand_fibres_oriented", exp.orientedFibreMaterial("wool_white"),
                        fibresToStrands(fs, map, int(FibrePopulation::Directed)), &err) ||
        !exp.addStrands("strand_fibres_random", exp.randomFibreMaterial("wool_white"),
                        fibresToStrands(fs, map, int(FibrePopulation::Random)), &err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 2;
    }
    // Cameras from the fibre front (+z), tilted towards -y (the export's up axis is fixed).
    const Vec3 view = normalized(Vec3{0, -1, 1.5});
    exp.extraCameras.push_back({"closeup", map.toExport(view * 25.0), map.toExport({}), 20.0});
    exp.extraCameras.push_back({"medium", map.toExport(view * 60.0), map.toExport({}), 20.0});
    exp.extraCameras.push_back({"detail", map.toExport(view * 10.0), map.toExport({}), 20.0});
    exp.skyIntensity = sky;
    exp.groundPlane = ground;
    exp.sortParts();
    if (!writeSceneJson(exp, &err)) { std::fprintf(stderr, "%s\n", err.c_str()); return 2; }
    std::printf("wrote %s/scene.json (cameras closeup, medium)\n", dir.c_str());
    return 0;
}

// knit_check --shaping-demo inc|dec|cone [-v] [--shift x] [--depth x] [--legdepth x]
int shapingDemo(int argc, char** argv) {
    std::string kind = "inc";
    bool verbose = false;
    ShapingParams sh;
    GaugeParams gauge;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--shaping-demo" && i + 1 < argc) kind = argv[++i];
        else if (a == "--w" && i + 1 < argc) gauge.w = std::atof(argv[++i]);
        else if (a == "--h" && i + 1 < argc) gauge.h = std::atof(argv[++i]);
        else if (a == "-v") verbose = true;
        else if (a == "--shift" && i + 1 < argc) sh.decShift = std::atof(argv[++i]);
        else if (a == "--depth" && i + 1 < argc) sh.decDepth = std::atof(argv[++i]);
        else if (a == "--legdepth" && i + 1 < argc) sh.decLegDepth = std::atof(argv[++i]);
        else if (a == "--blend" && i + 2 < argc) { sh.headBlend0 = std::atof(argv[++i]); sh.headBlend1 = std::atof(argv[++i]); }
        else if (a == "--decblend" && i + 2 < argc) { sh.decBlend0 = std::atof(argv[++i]); sh.decBlend1 = std::atof(argv[++i]); }
        else if (a == "--decramp" && i + 2 < argc) { sh.decDepth0 = std::atof(argv[++i]); sh.decDepth1 = std::atof(argv[++i]); }
    }
    ShapingDemo demo;
    if (!makeShapingDemo(kind, sh, demo, gauge)) { std::fprintf(stderr, "unknown demo %s\n", kind.c_str()); return 2; }
    const CurveSet& cs = demo.curves;
    const PatternParams& pp = demo.pp;
    ValidationReport r = validate(cs, pp, {.generatorChecks = false});  // shaped: no generator periodicity
    std::printf("%s", formatReport(r, cs, pp, verbose).c_str());
    if (verbose) {
        const KnitTopology& topo = cs.topo;
        for (const CurvatureResult& c : r.curvature)
            if (c.noFold != Status::Pass) {
                int y = 0;
                while (y < static_cast<int>(cs.yarns.size()) && cs.yarns[y].id != c.yarnId) ++y;
                int t = topo.stitchAt(y, c.uRhoEst);
                std::printf("    fold: yarn %d u %.3f rho_est %.4f mm, stitch %d (%s, turn %d, index %d), local u %.2f\n",
                            c.yarnId, c.uRhoEst, c.rhoEst, t, t >= 0 ? toString(topo.stitches[t].op) : "-",
                            t >= 0 ? topo.stitches[t].turn : -1, t >= 0 ? topo.stitches[t].index : -1,
                            t >= 0 ? c.uRhoEst - topo.stitches[t].uBegin : 0.0);
            }
    }
    if (verbose && cs.kind == CurveKind::Surface) {
        printDefects(cs, r, pp);
        int shown = 0;
        for (const LocalPairResult& lp : r.localSignature.pairs) {
            if (lp.status == Status::Pass || ++shown > 3) continue;
            std::printf("    around stitch %d:\n", lp.child);
            for (int k = std::max(0, lp.child - 3); k <= std::min(static_cast<int>(demo.lay.nodes.size()) - 1, lp.child + 3); ++k) {
                const LayoutNode& nd = demo.lay.nodes[k];
                const StitchRecord& st = demo.prog.stitches[k];
                double dPrev = k > 0 ? norm(nd.p - demo.lay.nodes[k - 1].p) : 0;
                std::printf("      %d %s slot %d turn %d idx %d parents %d,%d  w %.2f h %.2f  dprev %.2f  c.v %.2f\n", k,
                            toString(st.op), st.incSlot, st.turn, st.index, st.parents[0], st.parents[1], nd.w, nd.h,
                            dPrev, dot(nd.c, nd.v));
            }
        }
    }
    if (cs.kind == CurveKind::Surface)
        for (const LocalPairResult& lp : r.localSignature.pairs)
            if (lp.status != Status::Pass && verbose)
                std::printf("    stitch %d (%s): %s, %d crossings, %d unexpected\n", lp.child,
                            toString(cs.topo.stitches[lp.child].op), lp.detail.c_str(), lp.crossings, lp.unexpected);
    if (cs.kind == CurveKind::FlatPatch && verbose)
        for (const PairSignature& ps : r.signature.pairs)
            if (ps.status != Status::Pass)
                std::printf("    pair child %d (%s) parent %d: %s\n", ps.child, toString(cs.topo.stitches[ps.child].op),
                            ps.parent, ps.detail.c_str());
    return r.verdict == Status::Pass ? 0 : 1;
}

// knit_check --ply-map-sets [--yarn-d mm] [--ply-pitch mm] [--ply-fill f] [--size px]
//            [--cache dir]: ensures the default sets, prints hit / generated per set.
int plyMapSetsCli(int argc, char** argv) {
    double yarnD = 1.10;
    PlyGeometry g;
    int size = PlyMapSpec{}.size;
    std::string root = plyMapCacheRoot();
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto num = [&]() { return i + 1 < argc ? std::atof(argv[++i]) : 0.0; };
        if (a == "--ply-map-sets") continue;
        else if (a == "--yarn-d") yarnD = num();
        else if (a == "--ply-pitch") g.pitch = num();
        else if (a == "--ply-fill") g.fill = num();
        else if (a == "--size") size = std::max(8, static_cast<int>(num()));
        else if (a == "--cache" && i + 1 < argc) root = argv[++i];
        else {
            std::fprintf(stderr, "--ply-map-sets: unknown option %s\n", a.c_str());
            return 2;
        }
    }
    int rc = 0;
    for (PlyMapSpec s : defaultPlyMapSpecs(yarnD, g)) {
        s.size = size;
        PlyMapSetResult r;
        std::string err;
        if (!ensurePlyMapSet(root, s, r, &err)) {
            std::fprintf(stderr, "%d plies %s: %s\n", s.g.count, s.pp.model == PlyProfileModel::Swept ? "swept" : "circles", err.c_str());
            rc = 1;
            continue;
        }
        std::printf("%d plies %-7s D %.2f pitch %.2f: %s %s (%.0f ms)\n", s.g.count, s.pp.model == PlyProfileModel::Swept ? "swept" : "circles",
                    s.yarnD, s.g.pitch, r.hit ? "hit" : "generated", r.dir.c_str(), r.ms);
    }
    return rc;
}
