// The Fibres step's reference samples; see reference_samples.h.
#include "scene/fibres/reference_samples.h"

#include <algorithm>
#include <cmath>

#include "model/choices.h"
#include "model/mesh_query.h"
#include "fibres/fibre_generator.h"
#include "fibres/fur.h"
#include "fibres/fur_field.h"
#include "fibres/fur_loops.h"
#include "knitting/checks/checks.h"
#include "knitting/stitch/stockinette.h"
#include "scene/appearance/rendered_part.h"
#include "scene/fibres/fur_stage.h"
#include "scene/stage_params.h"
#include "settings/settings_mapping.h"

namespace knit {

namespace {

TriMesh uvSphere(double R, int lat, int lon) {
    TriMesh m;
    for (int i = 0; i <= lat; ++i)
        for (int j = 0; j < lon; ++j) {
            const double th = kPi * i / lat, ph = 2 * kPi * j / lon;
            m.pos.push_back({R * std::sin(th) * std::cos(ph), R * std::sin(th) * std::sin(ph), R * std::cos(th)});
        }
    for (int i = 0; i < lat; ++i)
        for (int j = 0; j < lon; ++j) {
            const int a = i * lon + j, b = i * lon + (j + 1) % lon, c = a + lon, e = b + lon;
            if (i > 0) m.tri.push_back({a, c, b});          // outward (counter-clockwise from outside)
            if (i + 1 < lat) m.tri.push_back({b, c, e});
        }
    m.computeNormals();
    return m;
}

} // namespace

std::optional<WoolSample> growWoolSample(const SceneSettings& s, std::string* err) {
    WoolSample out;
    const std::optional<StitchPattern> pattern = parsePattern(s.pattern);
    if (!pattern) {
        if (err) *err = "unknown pattern " + s.pattern;
        return std::nullopt;
    }
    PatternParams pp = flatPatternParams(s);
    applyPatternShape(pp, s, *pattern);  // the pattern's relaxed shape
    const StitchSize st{pp.waleSpacing, pp.courseSpacing};
    // A square patch: 5 wales and as many courses as fill the same length, less one for the loops'
    // overhang beyond the last course (7 at 5.5 x 3.4 mm: 27.5 x 27.4 mm).
    pp.wales = 5;
    pp.courses = std::max(1, int(std::lround(pp.wales * st.w / st.h - 1)));
    const CurveSet cs = generateStockinette(pp);
    const Allowances al = Allowances::forExtent(cs.extent());
    double kappa = 0;
    for (const YarnPath& y : cs.yarns) kappa = std::max(kappa, std::min(checkCurvature(y, pp.radius(), al).kappaHi, 1e3));
    const std::optional<RenderedPart> drawn = renderedStrands(cs, pp, kappa, s, err);
    if (!drawn) {
        if (err) *err = "the reference patch: " + *err;
        return std::nullopt;
    }
    const RenderedPart& rp = *drawn;
    const bool plied = rp.plies != nullptr;
    SurfaceParams sp;
    sp.compression = rp.fs.compression && !plied;
    sp.maxCompression = rp.fs.maxCompression;
    out.surface = buildSurface(*rp.strands, rp.strandPp, sp, ExecPolicy{true});
    const SurfaceModel& sm = out.surface;
    out.plyTubes = plied;
    // The twisted-ply appearance (single tube only).
    PlyProfileModel model = PlyProfileModel::Swept;
    if (!plied && parseTextureProfile(s.textureProfile, out.reliefPlies, model, nullptr) && out.reliefPlies > 1) {
        PlyGeometry g;
        g.count = out.reliefPlies;
        g.pitch = s.plyPitch;
        g.fill = s.plyFill;
        g.handedness = s.plyHandedness;
        PlyProfileParams p;
        p.model = model;
        out.relief = std::make_unique<PlyShadingProfile>(buildPlyShadingProfile(g, rp.ppR.radius(), p));
    }
    // Fibres: the model's generator and settings on the patch.
    FibreParams fp;
    if (!fibreParamsFromSettings(s, fp, err)) return std::nullopt;
    FibreFrame fr;
    FibreSetup fset;
    fset.tubeRadius = rp.strandPp.radius();
    fset.sides = rp.sides;
    fset.ring = rp.ring;
    fset.kappaMax = rp.kappaMax;
    fset.strands = rp.strands;
    fset.stitches = &cs.topo;
    if (plied) {
        fset.plies = rp.plies.get();
        fset.parentCurves = rp.yarn;
        fset.parentSurface = rp.parentSm.get();
    }
    setupPieceFibres(fp, fr, pp, fset);
    fr.exec.parallel = true;
    const FibreSet all = growYarnFibres(*rp.strands, sm, rp.strandPp, fp, &fr);
    // The central square, 60 % of the patch's side (no fibres from the patch's cut ends).
    Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    for (const YarnSurface& ys : sm.yarns)
        for (const Ring& R : ys.rings) {
            lo = {std::min(lo.x, R.c.x), std::min(lo.y, R.c.y), 0};
            hi = {std::max(hi.x, R.c.x), std::max(hi.y, R.c.y), 0};
        }
    const double half = 0.3 * std::min(hi.x - lo.x, hi.y - lo.y);
    const Vec3 mid = (lo + hi) * 0.5;
    for (const Fibre& f : all.fibres) {
        const Vec3 a = toVec3(all.pointsOf(f)[1]);
        if (std::abs(a.x - mid.x) > half || std::abs(a.y - mid.y) > half) continue;
        (f.population == FibrePopulation::Random ? out.stray : out.along).addFibre(f, all.pointsOf(f));
    }
    out.candidates = all.rep.candidates;
    out.square = 2 * half;
    return out;
}

FurSample growFurSample(const SceneSettings& s) {
    FurSample out;
    FurLook look = furLookFromSettings(s);
    const double R = std::max(10.0, 3 * std::max(look.fibres.lengthMedian, look.loops.curl > 0 ? 2 * look.loops.curl : look.loops.straightLength));
    out.radius = R;
    out.sphere = uvSphere(R, 48, 96);
    const TriMesh& sphere = out.sphere;
    const MeshQuery body(sphere);
    // The cap: the upper half and a band below the equator (fibres drape over it).
    std::vector<int> capTris;
    TriMesh cap;
    cap.pos = sphere.pos;
    for (size_t t = 0; t < sphere.tri.size(); ++t) {
        const Vec3 c = (sphere.pos[sphere.tri[t][0]] + sphere.pos[sphere.tri[t][1]] + sphere.pos[sphere.tri[t][2]]) / 3.0;
        if (c.z > -0.3 * R) {
            cap.tri.push_back(sphere.tri[t]);
            capTris.push_back(int(t));
        }
    }
    cap.computeNormals();
    // At most ~4000 fibres of both layers: the preview thins a dense coat (stated in the note).
    const double expected = (look.loops.density + (look.fibres.enabled ? look.fibres.density : 0)) * cap.area();
    double thin = 1;
    if (expected > 4000) {
        thin = 4000 / expected;
        look.loops.density *= thin;
        look.fibres.density *= thin;
    }
    out.thinning = thin;
    const ExecPolicy exec{true};
    const FurRootCache loops = sampleFurRoots(cap, FurLayer::Loops, furDensityBound(look, FurLayer::Loops), furPartSeed(look.seed, FurLayer::Loops, 0), nullptr, exec);
    const FurRootCache fibres = sampleFurRoots(cap, FurLayer::Fibres, furDensityBound(look, FurLayer::Fibres), furPartSeed(look.seed, FurLayer::Fibres, 0), nullptr, exec);
    const FurPartModifiers mods = buildFurModifiers(cap, look, 0);
    const FurField field = buildFurField({loops, fibres, look, nullptr, &mods}, exec);
    FurReport rep;
    const std::vector<Hair> hairs = growLoops(field, body, look.settle, exec, &rep);
    for (const Hair& h : hairs) addHairFibre(out.hairs, h, {});
    out.sampled = field.loops.sampled + field.fibres.sampled;
    out.stuck = rep.loops.stuck + rep.fibres.stuck;
    return out;
}

} // namespace knit
