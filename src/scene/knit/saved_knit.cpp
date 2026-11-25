// Writes a knitted scene to bytes and restores it, refusing anything that does not match the scene.
#include "scene/knit/saved_knit.h"
#include "scene/stage_params.h"

#include <cstring>

#include "model/cancel.h"
#include "model/hash.h"
#include "model/parallel.h"
#include "scene/knit/knit_cache.h"
#include "knitting/relax/settle.h"
#include "gpu/knitting.h"
#include "scene/scene.h"
#include "settings/settings_mapping.h"
#include "scene/source/scene_parts.h"

namespace knit {

namespace {

const char kMagic[] = "tiled-knit v1\n";

void putChart(ByteWriter& w, const ChartTiling& c) {
    w.i64(c.chart); w.i64(static_cast<int>(c.status)); w.f64(c.baseAngleDeg); w.i64(c.upFromWorld);
    w.i64(c.candidates); w.i64(c.kept); w.i64(c.outside); w.i64(c.ambiguous); w.i64(c.boundary); w.i64(c.insideOther); w.i64(c.shortRun);
    w.i64(c.runEnd); w.i64(c.runs); w.i64(c.droppedRuns); w.f64(c.coverage);
}
ChartTiling getChart(ByteReader& r) {
    ChartTiling c;
    c.chart = int(r.i64()); c.status = static_cast<Status>(r.i64()); c.baseAngleDeg = r.f64(); c.upFromWorld = r.i64() != 0;
    c.candidates = int(r.i64()); c.kept = int(r.i64()); c.outside = int(r.i64()); c.ambiguous = int(r.i64()); c.boundary = int(r.i64());
    c.insideOther = int(r.i64()); c.shortRun = int(r.i64()); c.runEnd = int(r.i64()); c.runs = int(r.i64()); c.droppedRuns = int(r.i64());
    c.coverage = r.f64();
    return c;
}

// v4: the gap patches' seeds (made again on restore, before the parts are compared).
void writeGapPatches(ByteWriter& w, const Scene& t) {
    int64_t patches = 0;
    for (const ScenePart& r : t.parts) patches += r.gapPatch ? 1 : 0;
    w.i64(patches);
    for (const ScenePart& r : t.parts) {
        if (!r.gapPatch) continue;
        w.string64(r.name); w.i64(r.prim); w.i64(r.component);
        w.i64(r.patch.tri); w.f64(r.patch.bary[0]); w.f64(r.patch.bary[1]); w.f64(r.patch.bary[2]);
        w.vec3(r.patch.wale); w.f64(r.patch.radius);
    }
}

// The gap patches, made again as addGapPatches made them (from their component's first part).
bool readGapPatches(ByteReader& r, Scene& out, std::string& why) {
    const size_t patches = r.count(8);
    for (size_t k = 0; k < patches && r.ok(); ++k) {
        ScenePart pr;
        pr.name = r.string64();
        pr.prim = int(r.i64());
        pr.component = int(r.i64());
        pr.gapPatch = true;
        pr.patch.tri = int(r.i64());
        for (double& b : pr.patch.bary) b = r.f64();
        pr.patch.wale = r.vec3();
        pr.patch.radius = r.f64();
        const ScenePart* first = nullptr;
        for (const ScenePart& p : out.parts)
            if (!p.gapPatch && p.prim == pr.prim && p.component == pr.component) {
                first = &p;
                break;
            }
        if (!first) {
            why = "gap patch " + pr.name + ": its component is not in the scene";
            return false;
        }
        pr.componentName = first->componentName;
        pr.centre = first->centre;
        pr.mesh = first->mesh;
        out.parts.push_back(std::move(pr));
    }
    return true;
}

void writePart(ByteWriter& w, const ScenePart& r) {
    w.string64(r.name); w.i64(r.prim); w.i64(r.component); w.vec3(r.centre); w.u64(meshDigest(*r.mesh));
    w.f64(r.yarnD); w.string64(r.note); w.f64(r.ms); w.string64(r.profile);
    w.i64(r.knit ? 1 : 0);
    if (!r.knit) return;
    const PartKnit& k = *r.knit;
    writePatternParams(w, k.t.pp);
    writeCurves(w, k.t.curves);
    writeStitches(w, k.t.prog.stitches);
    w.i64(int64_t(k.t.charts.size()));
    for (const ChartTiling& c : k.t.charts) putChart(w, c);
    w.i64(static_cast<int>(k.t.lattice.status)); w.i64(k.t.lattice.stitches); w.i64(k.t.lattice.parentless); w.i64(k.t.lattice.uncaught);
    w.string64(k.t.lattice.detail);
    w.i64(static_cast<int>(k.t.layout.status)); w.i64(k.t.layout.cells); w.i64(k.t.layout.outOfBounds); w.i64(k.t.layout.skewed);
    w.i64(k.t.layout.inverted); w.string64(k.t.layout.detail);
    w.string64(k.t.detail);
    w.i64(static_cast<int>(k.boundary.b1)); w.i64(static_cast<int>(k.boundary.b2)); w.f64(k.boundary.minDistance); w.i64(k.boundary.b1Fail);
    w.i64(k.boundary.b1Unresolved); w.i64(k.boundary.b2Samples);
    for (int u : k.boundary.b2Unresolved) w.i64(u);
    w.string64(k.boundary.detail);
    w.i64(static_cast<int>(k.gateDStatus)); w.i64(k.projectionUnresolved); w.i64(static_cast<int>(k.verdict)); w.string64(k.detail);
    w.i64(static_cast<int>(k.after.verdict));
    w.u64(fnv1a64(validationSnapshot(k.after)));
    w.i64(k.gateAOnCpu ? 1 : 0);
    w.f64(k.msTile); w.f64(k.msValidate); w.f64(k.msSettle); w.f64(k.msBoundary); w.f64(k.msGateD);
}

// Part `index` of the scene from the record: its identity must match the scene's. The atlas (UV tiling
// only) and the mesh surface are rebuilt from the part's mesh (deterministic).
bool readPart(ByteReader& r, int64_t version, size_t index, Scene& out, std::string& why) {
    ScenePart& p = out.parts[index];
    const std::string name = r.string64();
    const int prim = int(r.i64()), comp = int(r.i64());
    const Vec3 centre = r.vec3();
    const uint64_t md = r.u64();
    if (!r.ok() || name != p.name || prim != p.prim || comp != p.component) {
        why = "part " + std::to_string(index) + " differs (" + name + ")";
        return false;
    }
    if (md != meshDigest(*p.mesh) || centre.x != p.centre.x || centre.y != p.centre.y || centre.z != p.centre.z) {
        why = "the mesh of part " + name + " differs";
        return false;
    }
    p.yarnD = r.f64(); p.note = r.string64(); p.ms = r.f64(); p.profile = r.string64();
    if (r.i64() == 0) return true;
    auto k = std::make_unique<PartKnit>();
    k->t.pp = readPatternParams(r);
    readCurves(r, k->t.curves);
    k->t.prog.stitches = readStitches(r);
    if (!r.ok()) {
        why = "truncated curves of part " + name;
        return false;
    }
    k->t.charts.resize(r.count(8));
    for (ChartTiling& c : k->t.charts) c = getChart(r);
    k->t.lattice.status = static_cast<Status>(r.i64()); k->t.lattice.stitches = int(r.i64()); k->t.lattice.parentless = int(r.i64());
    k->t.lattice.uncaught = int(r.i64()); k->t.lattice.detail = r.string64();
    k->t.layout.status = static_cast<Status>(r.i64()); k->t.layout.cells = int(r.i64()); k->t.layout.outOfBounds = int(r.i64());
    k->t.layout.skewed = int(r.i64()); k->t.layout.inverted = int(r.i64()); k->t.layout.detail = r.string64();
    k->t.detail = r.string64();
    k->boundary.b1 = static_cast<Status>(r.i64()); k->boundary.b2 = static_cast<Status>(r.i64()); k->boundary.minDistance = r.f64();
    k->boundary.b1Fail = int(r.i64()); k->boundary.b1Unresolved = int(r.i64()); k->boundary.b2Samples = static_cast<long>(r.i64());
    for (int& u : k->boundary.b2Unresolved) u = int(r.i64());
    k->boundary.detail = r.string64();
    k->gateDStatus = static_cast<Status>(r.i64()); k->projectionUnresolved = int(r.i64()); k->verdict = static_cast<Status>(r.i64());
    k->detail = r.string64();
    k->after.verdict = static_cast<Status>(r.i64());
    out.afterDigests[index] = r.u64();
    k->gateAOnCpu = version >= 3 && r.i64() != 0;
    k->msTile = r.f64(); k->msValidate = r.f64(); k->msSettle = r.f64(); k->msBoundary = r.f64(); k->msGateD = r.f64();
    if (!r.ok()) {
        why = "truncated part " + name;
        return false;
    }
    if (p.layout == KnitLayout::UvTiling) k->atlas = std::make_unique<UvAtlas>(*p.mesh, Vec3{0, 0, 1});
    k->surface = std::make_unique<MeshSurface>(*p.mesh);
    p.knit = std::move(k);
    return true;
}

void writeCrossReport(ByteWriter& w, const CrossComponentReport& cross) {
    w.i64(static_cast<int>(cross.status)); w.i64(cross.pairs); w.i64(cross.closeSegments); w.f64(cross.minDistance);
    w.i64(int64_t(cross.lines.size()));
    for (const std::string& l : cross.lines) w.string64(l);
}

void readCrossReport(ByteReader& r, CrossComponentReport& cross) {
    cross.status = static_cast<Status>(r.i64()); cross.pairs = int(r.i64()); cross.closeSegments = static_cast<long>(r.i64());
    cross.minDistance = r.f64();
    cross.lines.resize(r.count(8));
    for (std::string& l : cross.lines) l = r.string64();
}

} // namespace

std::string serializeKnit(const Scene& t) {
    ByteWriter w;
    w.bytes(kMagic, std::strlen(kMagic));
    w.i64(kSavedKnitVersion);
    writeGapPatches(w, t);
    w.i64(int64_t(t.parts.size()));
    for (const ScenePart& r : t.parts) writePart(w, r);
    writeCrossReport(w, t.cross);
    w.u64(fnv1a64(w.written()));  // the checksum of the whole record
    return w.written();
}

bool restoreKnit(const std::string& bytes, const SceneEntry& e, const SceneSettings& s, Scene& out, std::string* err) {
    auto fail = [&](const std::string& m) {
        if (err) *err = "knit cache: " + m;
        return false;
    };
    const size_t ml = std::strlen(kMagic);
    if (bytes.compare(0, ml, kMagic) != 0) return fail("not a tiled knit");
    uint64_t sum = 0;
    if (bytes.size() < ml + sizeof sum) return fail("truncated");
    std::memcpy(&sum, bytes.data() + bytes.size() - sizeof sum, sizeof sum);
    if (sum != fnv1a64(bytes.substr(0, bytes.size() - sizeof sum))) return fail("damaged (checksum)");
    if (!loadPartsToKnit(e, s, SourceOptions{}, out, err)) return false;
    ByteReader r(bytes, ml);
    const int64_t version = r.i64();
    if (version != kSavedKnitVersion && version != 3 && version != 2) return fail("another version");
    std::string why;
    if (version >= 4 && !readGapPatches(r, out, why)) return fail(why);
    const size_t n = r.count(8);
    if (!r.ok() || n != out.parts.size()) return fail("the part count differs from the scene's");
    out.restored = true;
    out.afterDigests.assign(n, 0);
    for (size_t i = 0; i < n; ++i) {
        cancelCheckpoint();
        if (!readPart(r, version, i, out, why)) return fail(why);
    }
    for (ScenePart& p : out.parts) {
        const HandRequest hand = handRequest(s, p.name);
        if (p.knit && hand.inUse()) setYarnThickness(p.knit->t.curves, drawStitchTensions(p.knit->t.curves.topo, hand), hand);
    }
    readCrossReport(r, out.cross);
    r.u64();  // the checksum (compared above)
    if (!r.ok() || r.at() != bytes.size()) return fail("truncated or trailing data");
    out.verified = false;
    return true;
}

bool verifyKnit(Scene& t, const SceneSettings& s, std::string* err) {
    if (!t.restored || t.verified) return true;
    std::string e;
    std::unique_ptr<KnitNumerics> numerics = gpu::makeKnitNumerics(&e);
    if (!numerics) {
        if (err) *err = "verifying a cached knit needs a usable CUDA GPU: " + e;
        return false;
    }
    KnitNumerics& num = *numerics;
    for (size_t i = 0; i < t.parts.size(); ++i) {
        cancelCheckpoint();
        ScenePart& p = t.parts[i];
        if (!p.knit) continue;
        // Gate (a) as the knitting ran it: with the generator checks, or without them after a settle
        // or repair round (settle.cpp); by the knitting's numerics, or by the CPU validator for a
        // spiral piece (its audit). The stored digest must equal one of the two variants (the
        // curves are the stored ones either way).
        bool same = false;
        for (const bool gen : {true, false}) {
            const ValidateOptions vo{.generatorChecks = gen, .exec = ExecPolicy{true}};
            ValidationReport r = p.knit->gateAOnCpu ? validate(p.knit->t.curves, p.knit->t.pp, vo)
                                                    : num.validate(p.knit->t.curves, p.knit->t.pp, vo);
            if (fnv1a64(validationSnapshot(r)) == t.afterDigests[i]) {
                p.knit->after = std::move(r);
                same = true;
                break;
            }
            cancelCheckpoint();
        }
        if (!same) {
            if (err) *err = "gate (a) on the cached curves of " + p.name + " differs from the knitting run";
            return false;
        }
    }
    t.verified = true;
    return true;
}

} // namespace knit
