// GPU port tests: the alongYarn fibre generator (gpu::generateFibres) on a small deterministic
// patch. (1) Properties of every GPU fibre that hold by construction: finite geometry, the root
// inside its own yarn's tube, lengths and radii within their draw bounds, the population and the
// target length equal to the candidate's draw (drawFibre: a pure function of seed, yarn and
// candidate), unique candidate identities. (2) The CPU generator on the same input and frame, both
// under the light rules. Classification: both generators trace every
// candidate's outcome and population (FibreFrame::outcomes, in the same order); a FLIP is a
// candidate whose outcome or population differs between them. The traces carry no margin to the
// deciding threshold, so no flip can be attributed to a float32 boundary: every flip fails the
// test (printed with both outcomes). Geometry: matched fibres agree within a float32 error
// budget. Not byte equality: the device evaluates in float32. (3) The fabric-back test's closest-point search
// started from a triangle (meshClosestFrom) returns the plain search's result bit for bit, ties included.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <map>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "gpu_common.h"
#include "scene/stage_params.h"
#include "settings/file_presets.h"
#include "settings/settings_mapping.h"
#include "knitting/checks/checks.h"
#include "fibres/fibre_generator.h"
#include "fibres/fibre_shape.h"
#include "fibres/fibre_islands.h"
#include "knitting/knit/fixtures.h"
#include "knitting/stitch/stockinette.h"
#include "knitting/tube/yarn_surface.h"
#include "model/mesh.h"
#include "gpu/fibre_generator.h"

using namespace knit;

namespace {

struct FibrePatch {
    PatternParams pp = findPreset("tightest")->params;
    CurveSet cs = generateStockinette(pp);
    SurfaceModel sm;
    FibreParams fp = pieceFibres();
    std::vector<std::vector<Vec3>> centre;  // dense centreline samples per yarn
    double maxChord = 0;                    // largest distance between consecutive samples [mm]
    FibrePatch() {
        sm = buildSurface(cs, pp, SurfaceParams{});
        fp.density = 30;  // candidates per mm of yarn: a few thousand on this patch
        for (const YarnPath& y : cs.yarns) {
            std::vector<Vec3> c;
            const double u1 = y.uEnd();
            for (double u = 0; u <= u1; u += 0.002) {
                c.push_back(y.curve.eval(u));
                if (c.size() > 1) maxChord = std::max(maxChord, norm(c.back() - c[c.size() - 2]));
            }
            centre.push_back(std::move(c));
        }
    }
    // Distance from p to yarn y's centreline (dense samples: overestimates by at most maxChord / 2).
    double toCentre(int y, const Vec3& p) const {
        double d = 1e300;
        for (const Vec3& c : centre[static_cast<size_t>(y)]) d = std::min(d, norm(p - c));
        return d;
    }
    int yarnIndex(int yarnId) const {
        for (size_t k = 0; k < cs.yarns.size(); ++k)
            if (cs.yarns[k].id == yarnId) return static_cast<int>(k);
        return -1;
    }
};

bool finite(const Vec3f& v) { return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z); }

} // namespace

namespace {

// The checks below on the patch with the parameters fp.
void agreementOn(const FibrePatch& P, FibreParams fp, const char* what) {
    std::printf("  %s\n", what);
    fp.lightChecks = true;  // the rules the CPU generator runs
    FibreFrame cpuFrame;  // no callbacks: the flat patch
    if (fp.islands.inUse()) cpuFrame.islands = islandWeights({P.cs, P.cs.topo, fp.islands, fp.stitchArea, 5});
    std::vector<uint8_t> cpuOut, gpuOut;
    cpuFrame.outcomes = &cpuOut;
    const FibreSet cpu = growYarnFibres(P.cs, P.sm, P.pp, fp, &cpuFrame);
    FibreSet gpu;
    gpu::FibreRunReport rep;
    std::string err;
    CHECK(gpu::generateFibres({P.cs, P.sm, P.pp, fp, gpu::frameFrom(cpuFrame)}, gpu, &rep, &gpuOut, &err));
    std::printf("    patch: %zu yarns; GPU %lld candidates -> %lld exposed -> %zu fibres; CPU %zu fibres\n", P.cs.yarns.size(),
                rep.candidates, rep.exposed, gpu.fibres.size(), cpu.fibres.size());
    CHECK(gpu.fibres.size() > 200 && cpu.fibres.size() > 200);

    // (1) Every GPU fibre, by construction.
    const double r = P.pp.radius();
    const double sampleErr = P.maxChord / 2;
    std::map<std::pair<int, int>, const Fibre*> byId;
    int random = 0;
    for (const Fibre& f : gpu.fibres) {
        CHECK(f.pointCount >= 3 && f.candidate >= 0);
        for (const Vec3f& q : gpu.pointsOf(f)) CHECK(finite(q));
        CHECK(finite(f.anchor) && finite(f.normal) && std::isfinite(f.radius0) && std::isfinite(f.radius1));
        const int y = P.yarnIndex(f.yarnId);
        CHECK(y >= 0 && f.parentYarn == f.yarnId);  // no plies: the yarn supports itself
        CHECK(byId.emplace(std::make_pair(f.yarnId, f.candidate), &f).second);  // candidate identity unique
        // The candidate's draws (seed = fp.seed without a piece seed).
        const IslandWeight island = cpuFrame.islands.empty() ? IslandWeight{} : islandAt(cpuFrame.islands[y], P.cs.yarns[y].paramAtArcLength(f.s));
        const FibrePopulation pop = drawPopulation(fp, fp.seed, f.yarnId, f.candidate, island);
        CHECK(f.population == pop);
        random += f.population == FibrePopulation::Random;
        const FibreParams v = populationParams(fp, pop);
        const double target = drawFibre(fp, fp.seed, f.yarnId, f.candidate, island).length;
        CHECK(std::fabs(f.targetLength - target) <= 1e-5 * target);  // float32 copy of the draw
        const double gain = 1 + kFuzzyLength * island.fuzzy, loss = 1 - kMattedLength * island.matted;
        CHECK(f.targetLength >= v.lengthMin * loss * (1 - 1e-6) && f.targetLength <= v.lengthMax * gain * (1 + 1e-6));
        // Visible length: the polyline is within maxCurveDev of a curve of the target length, so
        // it cannot exceed the target by more than twice that (both chord ends). A curled stray's
        // span is solved on 16 chords over up to 3 turns: within 10 %.
        CHECK(visibleLength(gpu, f) > 0 && visibleLength(gpu, f) <= (fp.randomCurlRadius > 0 && pop == FibrePopulation::Random ? 1.1 * f.targetLength : f.targetLength + 2 * v.maxCurveDev));
        // Radius: thickness / 2 * exp(sigma z), z clamped to +-2.
        const double r0 = v.thickness / 2;
        CHECK(f.radius0 >= r0 * std::exp(-2 * v.thicknessSigma) * (1 - 1e-5) && f.radius0 <= r0 * std::exp(2 * v.thicknessSigma) * (1 + 1e-5));
        // Root placement: the buried first vertex lies inside the own tube (nominal radius r;
        // compression only shrinks the tube).
        CHECK(P.toCentre(y, toVec3(gpu.pointsOf(f)[0])) <= r + sampleErr);
    }
    std::printf("    GPU fibres: %d random, %zu directed; properties hold for all\n", random, gpu.fibres.size() - random);
    CHECK(random > 0 && random < static_cast<int>(gpu.fibres.size()));
    // No point from the curve start on, and no segment midpoint, behind the patch's back (the lowest point
    // of the tubes; the generators' capsules reach a little below the rings), on either generator.
    double back = 1e300;
    for (const YarnSurface& ys : P.sm.yarns)
        for (const Ring& R : ys.rings) back = std::min(back, R.c.z - std::max(R.a, R.b));
    const FibreSet* generated[] = {&cpu, &gpu};
    for (const FibreSet* set : generated)
        for (const Fibre& f : set->fibres) {
            const std::span<const Vec3f> pts = set->pointsOf(f);
            double low = pts[1].z;
            for (size_t k = 2; k < pts.size(); ++k) low = std::min({low, double(pts[k].z), 0.5 * (double(pts[k - 1].z) + pts[k].z)});
            CHECK(low >= back - 0.01);
        }

    // (2) CPU vs GPU on the behaviour they share.
    // Classification of every candidate (outcome in bits 0..6, population in bit 7).
    CHECK(cpuOut.size() == gpuOut.size() && static_cast<long long>(gpuOut.size()) == rep.candidates);
    int flips = 0;
    for (size_t i = 0; i < std::min(cpuOut.size(), gpuOut.size()); ++i)
        if (cpuOut[i] != gpuOut[i] && ++flips <= 10)
            std::printf("    FLIP candidate #%zu: CPU %s%s, GPU %s%s\n", i, toString(static_cast<FibreOutcome>(cpuOut[i] & 0x7f)),
                        (cpuOut[i] & 0x80) ? " (random)" : "", toString(static_cast<FibreOutcome>(gpuOut[i] & 0x7f)),
                        (gpuOut[i] & 0x80) ? " (random)" : "");
    // The accepted fibres matched by (strand, candidate).
    int matched = 0, cpuOnly = 0, segDiff = 0;
    double maxAnchor = 0, maxPoint = 0, maxLen = 0;
    for (const Fibre& c : cpu.fibres) {
        const auto it = byId.find({c.yarnId, c.candidate});
        if (it == byId.end()) {
            ++cpuOnly;
            continue;
        }
        const Fibre& g = *it->second;
        ++matched;
        CHECK(g.population == c.population);
        maxAnchor = std::max(maxAnchor, norm(toVec3(g.anchor) - toVec3(c.anchor)));
        // Every point: the packed pool puts each fibre's points in its place, at any stride.
        if (g.pointCount == c.pointCount)
            for (int k = 0; k < g.pointCount; ++k)
                maxPoint = std::max(maxPoint, norm(toVec3(gpu.pointsOf(g)[k]) - toVec3(cpu.pointsOf(c)[k])));
        else
            ++segDiff;
        maxLen = std::max(maxLen, std::fabs(visibleLength(gpu, g) - visibleLength(cpu, c)));
    }
    const int gpuOnly = static_cast<int>(gpu.fibres.size()) - matched;
    std::printf("    CPU vs GPU: %zu candidates, %d flips; %d accepted fibres matched, %d CPU only, %d GPU only, %d with another "
                "segment count; max |anchor| %.2e, |point| %.2e, |visible length| %.2e mm\n",
                gpuOut.size(), flips, matched, cpuOnly, gpuOnly, segDiff, maxAnchor, maxPoint, maxLen);
    CHECK(flips == 0 && cpuOnly == 0 && gpuOnly == 0 && segDiff == 0);
    // Geometry tolerance (an error budget, not the observed error). The device works in float32
    // on coordinates below 20 mm on this patch: 2^-24 relative is ~1.2e-6 mm per rounding. E* and
    // the tip come from short iterations (Newton on the ring, the bow evaluation) whose error grows
    // to at most a few hundred roundings: 1e-3 mm (0.2 % of the yarn radius, 1/40 of the thinnest
    // fibre's root radius) bounds them with margin.
    CHECK(maxAnchor <= 1e-3 && maxPoint <= 1e-3 && maxLen <= 1e-3);
}

} // namespace

TEST(gpu_fibres_properties_and_cpu_agreement) {
    REQUIRE_GPU();
    const FibrePatch P;
    agreementOn(P, P.fp, "pieceFibres (the former defaults)");
    // The viewer's default wool preset (presets/fibres, the lowest order; fibre_presets.h), through
    // the viewer's conversion.
    CHECK(!presets(PresetGroup::Fibres).empty());
    if (presets(PresetGroup::Fibres).empty()) return;
    const FilePreset& def = presets(PresetGroup::Fibres).front();
    SceneSettings s;
    applyPreset(PresetGroup::Fibres, s, def);
    FibreParams fp;
    std::string err;
    CHECK(fibreParamsFromSettings(s, fp, &err));
    fp.density = P.fp.density;
    CHECK(fp.maxSegments == def.values.fibreMaxSegments);
    agreementOn(P, fp, def.label.c_str());
    // Islands of all three kinds, read from one table by both generators.
    FibreParams islands = P.fp;
    islands.islands = {1.5, 10, 100, 1, 1, 1};
    agreementOn(P, islands, "fibre islands");
    // Fluffy wool: the long curled strays on both generators. Mohair: fuzz up to 12 mm and strays up to 30 mm.
    for (const FilePreset& p : presets(PresetGroup::Fibres))
        if (p.id == "fluffy" || p.id == "mohair") {
            SceneSettings t;
            applyPreset(PresetGroup::Fibres, t, p);
            FibreParams soft;
            CHECK(fibreParamsFromSettings(t, soft, &err));
            soft.density = P.fp.density;
            agreementOn(P, soft, p.label.c_str());
        }
    // Crimped directed fibres (fibres.crimp_*): short, fine and finely chorded, as a crimp needs.
    SceneSettings c;
    c.directed.lengthMin = 0.5;
    c.directed.lengthMedian = 0.8;
    c.directed.lengthMax = 1.5;
    c.directed.thickness = 0.025;
    c.fibreSegmentLength = 0.12;
    c.fibreMaxSegments = 24;
    c.crimpAmplitude = 0.08;
    c.crimpWavelength = 1.5;
    // The strays of this case keep the profile it was written with (any direction, rise up to 60 deg):
    // the case is about the crimp. Strays along the yarn are compared in every preset case above.
    c.randomSpreadDeg = 180;
    c.random.maxAngleDeg = 60;
    FibreParams crimped;
    CHECK(fibreParamsFromSettings(c, crimped, &err));
    crimped.density = P.fp.density;
    agreementOn(P, crimped, "crimped directed fibres");
    // Standard wool hugging the yarn along the ply helix (fibres.directed.hug_height, follow_plies).
    for (const FilePreset& p : presets(PresetGroup::Fibres))
        if (p.id == "standard") {
            SceneSettings t;
            applyPreset(PresetGroup::Fibres, t, p);
            t.directedHugHeight = 0.15;
            t.directedFollowPlies = true;
            FibreParams hug;
            CHECK(fibreParamsFromSettings(t, hug, &err));
            CHECK(hug.hugHeight == 0.15 && hug.plyPitch == t.plyPitch);
            hug.density = P.fp.density;
            agreementOn(P, hug, "Standard wool, hugging along the ply helix");
        }
    // Strays along the yarn (fibres.random.spread_deg), curled, long-tailed and a larger share.
    for (const FilePreset& p : presets(PresetGroup::Fibres))
        if (p.id == "soft") {
            SceneSettings t;
            applyPreset(PresetGroup::Fibres, t, p);
            t.randomFraction = 0.3;
            t.randomSpreadDeg = 30;
            t.random.maxAngleDeg = 30;
            t.random.lengthMax = 24;
            t.random.lengthSpread = 0.85;
            FibreParams spread;
            CHECK(fibreParamsFromSettings(t, spread, &err));
            CHECK(spread.randomSpreadDeg == 30);
            spread.density = P.fp.density;
            agreementOn(P, spread, "Soft wool, strays within 30 deg of the yarn");
        }
}

namespace {

FibreSet gpuPatchFibres(const FibrePatch& P, const FibreParams& fp, std::vector<uint8_t>& outcomes) {
    FibreSet out;
    gpu::FibreRunReport rep;
    std::string err;
    CHECK(gpu::generateFibres({P.cs, P.sm, P.pp, fp, gpu::frameFrom(FibreFrame{})}, out, &rep, &outcomes, &err));
    return out;
}

int countOutcome(const std::vector<uint8_t>& outcomes, FibreOutcome o) {
    int n = 0;
    for (uint8_t code : outcomes) n += (code & 0x7f) == static_cast<uint8_t>(o);  // bit 7: the population
    return n;
}

} // namespace

// The full checks repeat the anchor's tests at E*, where the root segment leaves the own strand (EmergeRejected).
// E* lies within about 1e-4 mm of the anchor, so a fixed configuration rarely lands between the two; this test
// finds one. The light rules record the anchor's normal, the full rules the normal at E*; on the flat patch the
// front value is n.z. A candidate accepted by both whose front value is lower at E* is rejected at the emergence by
// a threshold between its two values: frontCos is read only by these two tests, and candidates are independent.
TEST(gpu_fibres_full_checks_reject_at_the_emergence) {
    REQUIRE_GPU();
    const FibrePatch P;
    FibreParams light = P.fp, full = P.fp;
    light.lightChecks = true;
    full.lightChecks = false;
    std::vector<uint8_t> outcomes;
    std::map<std::pair<int, int>, float> anchorFront;
    for (const Fibre& f : gpuPatchFibres(P, light, outcomes).fibres) anchorFront[{f.strand, f.candidate}] = f.normal.z;
    float largestDrop = 0, threshold = 0;
    std::pair<int, int> chosen{-1, -1};
    int bothAccepted = 0;
    for (const Fibre& f : gpuPatchFibres(P, full, outcomes).fibres) {
        const auto it = anchorFront.find({f.strand, f.candidate});
        if (it == anchorFront.end()) continue;
        ++bothAccepted;
        if (it->second - f.normal.z > largestDrop) {
            largestDrop = it->second - f.normal.z;
            threshold = 0.5f * (it->second + f.normal.z);
            chosen = it->first;
        }
    }
    std::printf("  %d candidates accepted by both rules; largest front drop from the anchor to E* %.3g (strand %d, candidate %d)\n",
                bothAccepted, largestDrop, chosen.first, chosen.second);
    CHECK(largestDrop > 1e-6f);  // far above the last bits two kernels may differ in
    FibreParams raised = full;
    raised.frontCos = threshold;
    const FibreSet rejected = gpuPatchFibres(P, raised, outcomes);
    const int atEmergence = countOutcome(outcomes, FibreOutcome::EmergeRejected);
    std::printf("  front threshold %.7f: %d rejected at the emergence\n", threshold, atEmergence);
    CHECK(atEmergence >= 1);
    for (const Fibre& f : rejected.fibres) CHECK(!(f.strand == chosen.first && f.candidate == chosen.second));
}

namespace {

// A grid of nx x ny quads over [x0, x0 + w] x [y0, y0 + h], two triangles each, at the height z(x, y).
struct SheetGrid {
    double x0, y0, w, h;
    int nx, ny;
};
void addSheet(TriMesh& m, const SheetGrid& g, double (*z)(double, double)) {
    const int first = static_cast<int>(m.pos.size());
    for (int j = 0; j <= g.ny; ++j)
        for (int i = 0; i <= g.nx; ++i) {
            const double x = g.x0 + g.w * i / g.nx, y = g.y0 + g.h * j / g.ny;
            m.pos.push_back({x, y, z(x, y)});
        }
    const auto at = [&](int i, int j) { return first + j * (g.nx + 1) + i; };
    for (int j = 0; j < g.ny; ++j)
        for (int i = 0; i < g.nx; ++i) {
            m.tri.push_back({at(i, j), at(i + 1, j), at(i + 1, j + 1)});
            m.tri.push_back({at(i, j), at(i + 1, j + 1), at(i, j + 1)});
        }
}
double wavyHeight(double x, double y) { return 1.5 * std::sin(x / 7) * std::cos(y / 5); }
double foldHeight(double x, double) { return std::fabs(x); }  // a 90 deg V fold along the y axis
double zeroHeight(double, double) { return 0; }
double sixHeight(double, double) { return 6; }

TriMesh sphereMesh(double r, int nu, int nv) {
    TriMesh m;
    for (int j = 0; j <= nv; ++j)
        for (int i = 0; i < nu; ++i) {
            const double th = kPi * j / nv, ph = 2 * kPi * i / nu;
            m.pos.push_back({r * std::sin(th) * std::cos(ph), r * std::sin(th) * std::sin(ph), r * std::cos(th)});
        }
    for (int j = 0; j < nv; ++j)
        for (int i = 0; i < nu; ++i) {
            const int a = j * nu + i, b = j * nu + (i + 1) % nu, c = (j + 1) * nu + (i + 1) % nu, d = (j + 1) * nu + i;
            m.tri.push_back({a, d, c});
            m.tri.push_back({a, c, b});
        }
    return m;
}

// Stray-like polylines: from a point on or near the surface, rising up to `reach`, curling, then falling back
// through the surface (or past its edge) and behind it. Consecutive samples, as the fabric-back test asks them.
std::vector<Vec3> strayPolylines(Vec3 lo, Vec3 hi, double reach, int curves, int samples, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(0, 1);
    std::vector<Vec3> pts;
    for (int c = 0; c < curves; ++c) {
        const Vec3 p0{lo.x + (hi.x - lo.x) * u(rng), lo.y + (hi.y - lo.y) * u(rng), lo.z + (hi.z - lo.z) * u(rng)};
        const double az = 2 * kPi * u(rng);
        const Vec3 v{reach * std::cos(az), reach * std::sin(az), reach * (0.3 + u(rng))};
        const double fall = 2.5 * reach * u(rng), curl = 0.05 * reach * u(rng);
        for (int k = 0; k < samples; ++k) {
            const double t = double(k) / (samples - 1);
            pts.push_back(p0 + v * t + Vec3{curl * std::cos(25 * t), curl * std::sin(25 * t), -fall * t * t});
        }
    }
    return pts;
}

// Both searches on every point, twice: started from scattered triangles (mostly far from the answer, on the other
// sheet or the other side), then each from the previous point's closest triangle as along a fibre. Mismatches of
// any field, bit for bit.
int closestMismatches(const char* name, const TriMesh& mesh, const std::vector<Vec3>& pts) {
    const int nTri = static_cast<int>(mesh.tri.size());
    std::vector<int> start(pts.size());
    for (size_t k = 0; k < pts.size(); ++k) start[k] = static_cast<int>((k * 7919u) % static_cast<unsigned>(nTri));
    int bad = 0;
    for (int pass = 0; pass < 2; ++pass) {
        std::vector<gpu::MeshClosestPair> r;
        std::string err;
        if (!gpu::compareMeshClosest(mesh, pts, start, r, &err)) {
            std::printf("  %s: %s\n", name, err.c_str());
            return -1;
        }
        int passBad = 0;
        for (size_t k = 0; k < r.size(); ++k) {
            const gpu::MeshClosestPair& c = r[k];
            const bool same = c.tri[0] == c.tri[1] && std::memcmp(&c.dist[0], &c.dist[1], sizeof(float)) == 0 &&
                              std::memcmp(c.bary[0], c.bary[1], sizeof c.bary[0]) == 0 && std::memcmp(&c.p[0], &c.p[1], sizeof c.p[0]) == 0;
            if (!same && passBad++ < 3)
                std::printf("  %s pass %d point %zu: tri %d / %d, dist %.9g / %.9g\n", name, pass, k, c.tri[0], c.tri[1], c.dist[0], c.dist[1]);
        }
        std::printf("  %s pass %d: %zu points, %d mismatches\n", name, pass, r.size(), passBad);
        bad += passBad;
        for (size_t k = 1; k < pts.size(); ++k) start[k] = r[k - 1].tri[0];
    }
    return bad;
}

} // namespace

namespace {

// The Mohair profile's fibres (presets/fibres/mohair.preset) at the patch's density, drape on or off.
FibreParams mohairFibres(const FibrePatch& P, bool drape) {
    FibreParams fp;
    for (const FilePreset& p : presets(PresetGroup::Fibres))
        if (p.id == "mohair") {
            SceneSettings s;
            applyPreset(PresetGroup::Fibres, s, p);
            std::string err;
            CHECK(fibreParamsFromSettings(s, fp, &err));
        }
    CHECK(fp.maxSegments > 3);  // found
    fp.density = P.fp.density;
    fp.lightChecks = true;
    fp.drape = drape;
    return fp;
}

int truncatedDirected(const FibreSet& set) {
    int n = 0;
    for (const Fibre& f : set.fibres) n += f.truncated && f.population == FibrePopulation::Directed;
    return n;
}

// The deepest any accepted fibre's vertex or segment midpoint after the curve's start lies inside a strand, as a
// share of the yarn radius (the deep-contact limit is 0.2).
double deepestContact(const FibrePatch& P, const FibreSet& set) {
    const SurfaceIndex index(P.cs, P.sm);
    double deepest = 0;
    for (const Fibre& f : set.fibres) {
        const std::span<const Vec3f> pts = set.pointsOf(f);
        for (size_t k = 2; k < pts.size(); ++k)
            for (const Vec3& q : {0.5 * (toVec3(pts[k - 1]) + toVec3(pts[k])), toVec3(pts[k])}) {
                const SurfaceIndex::Hit h = index.exactInside(q, 0.0, {}, {});
                if (h.inside) deepest = std::max(deepest, h.depth / P.pp.radius());
            }
    }
    return deepest;
}

// Turns at the interior vertices of an accepted fibre [deg].
double largestTurn(const FibreSet& set, const Fibre& f) {
    const std::span<const Vec3f> p = set.pointsOf(f);
    double largest = 0;
    for (size_t k = 1; k + 1 < p.size(); ++k) {
        const Vec3 a = normalized(toVec3(p[k]) - toVec3(p[k - 1])), b = normalized(toVec3(p[k + 1]) - toVec3(p[k]));
        largest = std::max(largest, std::acos(std::clamp(dot(a, b), -1.0, 1.0)) * 180 / kPi);
    }
    return largest;
}

bool sameFibres(const FibreSet& a, const FibreSet& b) {
    if (a.fibres.size() != b.fibres.size()) return false;
    for (size_t i = 0; i < a.fibres.size(); ++i) {
        const std::span<const Vec3f> pa = a.pointsOf(a.fibres[i]), pb = b.pointsOf(b.fibres[i]);
        if (a.fibres[i].candidate != b.fibres[i].candidate || pa.size() != pb.size()) return false;
        if (pa.size() > 0 && std::memcmp(pa.data(), pb.data(), pa.size() * sizeof(Vec3f)) != 0) return false;
    }
    return true;
}

} // namespace

TEST(gpu_fibre_drape_lifts_instead_of_cutting) {
    // Mohair's fuzz (6-12 mm, lifted to 45 deg) reaches the neighbouring stitches. With the drape (fibres.directed.drape)
    // a fibre that would sink into one rests on it instead of being cut there, on both generators: fewer truncated
    // fuzz fibres, none deeper than the deep-contact limit, no draped fibre longer than its drawn length allows (the tip
    // gives back what the lift added), the turn limit kept, and the two generators still agree.
    REQUIRE_GPU();
    const FibrePatch P;
    const FibreParams off = mohairFibres(P, false), on = mohairFibres(P, true);
    FibreFrame frame;
    const FibreSet cpuOff = growYarnFibres(P.cs, P.sm, P.pp, off, &frame), cpuOn = growYarnFibres(P.cs, P.sm, P.pp, on, &frame);
    std::vector<uint8_t> outcomes;
    const FibreSet gpuOff = gpuPatchFibres(P, off, outcomes), gpuOn = gpuPatchFibres(P, on, outcomes);
    std::printf("  truncated fuzz fibres: CPU %d -> %d, GPU %d -> %d (of %zu / %zu accepted with the drape)\n", truncatedDirected(cpuOff),
                truncatedDirected(cpuOn), truncatedDirected(gpuOff), truncatedDirected(gpuOn), cpuOn.fibres.size(), gpuOn.fibres.size());
    CHECK(truncatedDirected(cpuOn) < truncatedDirected(cpuOff) && truncatedDirected(gpuOn) < truncatedDirected(gpuOff));
    const double deepCpu = deepestContact(P, cpuOn), deepGpu = deepestContact(P, gpuOn);
    std::printf("  deepest contact with the drape: CPU %.4f r, GPU %.4f r (limit 0.2 r)\n", deepCpu, deepGpu);
    CHECK(deepCpu <= 0.2 + 1e-9 && deepGpu <= 0.2 + 1e-4);  // GPU: float32 points checked in double
    int draped = 0;
    double worstExcess = -1e300;
    for (const FibreSet* set : {&cpuOn, &gpuOn})
        for (const Fibre& f : set->fibres) {
            if (f.population != FibrePopulation::Directed) continue;
            worstExcess = std::max(worstExcess, visibleLength(*set, f) - (f.targetLength + 2 * on.maxCurveDev));
            CHECK(largestTurn(*set, f) <= on.maxTurnDeg + 1e-3);
        }
    // The draped fibres: accepted uncut with the drape where the same candidate was cut without it.
    std::map<std::pair<int, int>, const Fibre*> cutWithout;
    for (const Fibre& f : cpuOff.fibres)
        if (f.truncated) cutWithout[{f.yarnId, f.candidate}] = &f;
    for (const Fibre& f : cpuOn.fibres) draped += !f.truncated && cutWithout.count({f.yarnId, f.candidate});
    std::printf("  CPU: %d fibres draped; largest visible length over target + 2 maxCurveDev: %.4f mm\n", draped, worstExcess);
    CHECK(draped > 0 && worstExcess <= 0);
    agreementOn(P, on, "Mohair with the drape");
}

TEST(gpu_fibre_drape_full_checks_unchanged) {
    // The drape is a rule of the light checks. Under the full checks (another collision path, no CPU twin) the key
    // changes nothing: the same outcomes and the same fibres, bit for bit.
    REQUIRE_GPU();
    const FibrePatch P;
    FibreParams off = mohairFibres(P, false), on = mohairFibres(P, true);
    off.lightChecks = on.lightChecks = false;
    std::vector<uint8_t> outOff, outOn;
    const FibreSet a = gpuPatchFibres(P, off, outOff), b = gpuPatchFibres(P, on, outOn);
    std::printf("  full checks: %zu fibres without, %zu with the drape key\n", a.fibres.size(), b.fibres.size());
    CHECK(!a.fibres.empty() && outOff == outOn && sameFibres(a, b));
}

TEST(gpu_mesh_closest_from_matches_full) {
    // The fabric-back test (keptInFront) starts each closest-point search at the previous sample's triangle and
    // visits only the cells within that distance (meshClosestFrom). It must return meshClosest's MeshHit bit for
    // bit, ties included: on folds and on parallel sheets, equally distant points of different triangles give
    // heights of opposite sign, and the tie is left to meshClosest.
    REQUIRE_GPU();
    TriMesh sheet;  // (a) an open, gently wavy swatch: strays to 50 mm, past the open edges, behind the back
    addSheet(sheet, {0, 0, 60, 40, 60, 40}, wavyHeight);
    sheet.computeNormals();
    const TriMesh sphere = [] { TriMesh m = sphereMesh(20, 48, 24); m.computeNormals(); return m; }();  // (b) closed
    TriMesh fold;  // (c) a V fold: points on the bisector plane x = 0 are equally far from both halves
    addSheet(fold, {-20, 0, 40, 30, 40, 30}, foldHeight);
    fold.computeNormals();
    TriMesh sheets;  // (c) two parallel sheets 6 mm apart: the mid-plane z = 3 is equally far from both
    addSheet(sheets, {0, 0, 30, 30, 30, 30}, zeroHeight);
    addSheet(sheets, {0, 0, 30, 30, 30, 30}, sixHeight);
    sheets.computeNormals();

    std::vector<Vec3> foldPts, sheetsPts;
    std::mt19937 rng(11);
    std::uniform_real_distribution<double> u(0, 1);
    for (int k = 0; k < 4000; ++k) {
        const double y = 30 * u(rng), z = 25 * u(rng);
        foldPts.push_back({0, y, z});                        // on the bisector: exact ties by symmetry
        foldPts.push_back({1e-4 * (u(rng) - 0.5), y, z});    // just off it
        sheetsPts.push_back({30 * u(rng), 30 * u(rng), 3});  // on the mid-plane
        sheetsPts.push_back({30 * u(rng), 30 * u(rng), 3 + 1e-4 * (u(rng) - 0.5)});
    }
    const std::vector<Vec3> foldStrays = strayPolylines({-15, 2, 0}, {15, 28, 15}, 30, 60, 40, 5);
    foldPts.insert(foldPts.end(), foldStrays.begin(), foldStrays.end());
    const std::vector<Vec3> sheetsStrays = strayPolylines({2, 2, 0}, {28, 28, 6}, 20, 60, 40, 6);
    sheetsPts.insert(sheetsPts.end(), sheetsStrays.begin(), sheetsStrays.end());

    const int a = closestMismatches("sheet", sheet, strayPolylines({0, 0, -1}, {60, 40, 1}, 50, 400, 50, 1));
    const int b = closestMismatches("sphere", sphere, strayPolylines({-20, -20, -20}, {20, 20, 20}, 50, 300, 50, 2));
    const int c = closestMismatches("fold", fold, foldPts);
    const int d = closestMismatches("parallel sheets", sheets, sheetsPts);
    CHECK(a == 0 && b == 0 && c == 0 && d == 0);
}
