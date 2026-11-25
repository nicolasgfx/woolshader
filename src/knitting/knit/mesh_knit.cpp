// Knitting one mesh component end to end, stage by stage, and its verdict (mesh_knit.h has the order).
#include "knitting/knit/mesh_knit.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <memory>

#include "model/cancel.h"
#include "model/log.h"

namespace knit {

namespace {

using Clock = std::chrono::steady_clock;

// The verdict of a knitted part: gate (a), lattice, layout, projection, B1, B2, and gate (d) when
// it counts.
Status partVerdict(const PartKnit& out, bool gateDInVerdict) {
    Status v = out.after.verdict;
    v = combine(v, out.t.lattice.status);
    v = combine(v, out.t.layout.status);
    v = combine(v, out.projectionUnresolved ? Status::Unresolved : Status::Pass);
    v = combine(v, out.boundary.b1);
    v = combine(v, out.boundary.b2);
    if (gateDInVerdict) v = combine(v, out.gateDStatus);
    return v;
}

} // namespace

bool knitMesh(const TriMesh& mesh, const MeshKnitParams& p, KnitNumerics& numerics, PartKnit& out, std::string* err) {
    out = PartKnit{};
    FinalSurfaceSpec fs = p.surface;
    fs.D = p.tile.pp.yarnDiameter;
    fs.c = p.tile.pp.clearance;
    if (std::string why; !fs.valid(&why)) {
        if (err) *err = why;
        return false;
    }
    TileParams tp = p.tile;
    if (tp.supportRadius <= 0) tp.supportRadius = fs.supportRadius();
    auto t0 = Clock::now();
    out.atlas = std::make_unique<UvAtlas>(mesh, p.tile.up);
    out.surface = numerics.meshSurface(mesh);
    if (!tileMesh({*out.atlas, *out.surface}, tp, numerics, out.t, err)) return false;
    out.msTile = msSince(t0);
    cancelCheckpoint();
    const PatternParams pp = out.t.pp;  // gate (a), settle
    HandResult hand;
    if (p.hand.inUse()) hand = knitByHand(out.t.curves, p.hand, numerics, p.exec);
    t0 = Clock::now();
    out.before = numerics.validate(out.t.curves, pp, ValidateOptions{.generatorChecks = true, .exec = p.exec});
    out.after = out.before;
    out.msValidate = msSince(t0);
    cancelCheckpoint();
    t0 = Clock::now();
    const ValidateStats v0 = validateStats();
    SettleParams sp = p.sp;
    sp.exec = p.exec;
    if (p.settle) out.settle = settle({out.t.curves, pp, out.surface.get(), &out.t.bind, p.hand.inUse() ? &hand.rest : nullptr}, sp, numerics, &out.after);
    out.msSettle = msSince(t0);
    out.vSettle = validateStats().since(v0);
    cancelCheckpoint();
    out.projectionUnresolved = out.t.bind.unresolved();
    t0 = Clock::now();
    out.boundary = numerics.checkBoundary(out.t.curves, out.t.bind, *out.surface, tp.supportRadius);
    out.msBoundary = msSince(t0);
    cancelCheckpoint();
    char b[400];
    std::snprintf(b, sizeof b, "%s | gate a %s | lattice %s | layout %s | projection %d unresolved | %s", out.t.detail.c_str(),
                  toString(out.after.verdict), toString(out.t.lattice.status), toString(out.t.layout.status),
                  out.projectionUnresolved, out.boundary.detail.c_str());
    out.detail = b;
    runSurfaceGate(out, p, numerics);
    return true;
}

void runSurfaceGate(PartKnit& out, const MeshKnitParams& p, KnitNumerics& numerics) {
    FinalSurfaceSpec fs = p.surface;
    const PatternParams pp = out.t.pp;
    fs.D = pp.yarnDiameter;
    fs.c = pp.clearance;
    out.gateD.clear();
    out.compressionClamped = -1;
    out.msGateDBuild.clear();
    out.msGateDCheck.clear();
    out.gateDStatus = Status::Pass;
    if (p.surfaceGate) {
        // The sampled surface check on the RENDERED surface (FinalSurfaceSpec).
        const auto t0 = Clock::now();
        const PatternParams ppR = fs.rendered(pp);
        for (int pass = 0; pass < (p.surfaceHalf ? 2 : 1); ++pass) {
            SurfaceParams sp;
            sp.compression = fs.compression;
            sp.maxCompression = fs.maxCompression;
            sp.ringSpacing = pass == 0 ? 0.025 : 0.0125;
            const auto tb = Clock::now();
            SurfaceCheck sc = numerics.checkRenderedSurface(out.t.curves, ppR, sp);
            if (fs.compression) out.compressionClamped = std::max(out.compressionClamped, sc.compressionClamped);
            out.msGateDBuild.push_back(0);  // model and check are one device call
            out.msGateDCheck.push_back(msSince(tb));
            out.gateD.push_back(std::move(sc));
            out.gateDStatus = combine(out.gateDStatus, combine(out.gateD.back().fold, out.gateD.back().collisions));
            cancelCheckpoint();
        }
        out.msGateD = msSince(t0);
    }
    out.verdict = partVerdict(out, p.gateDInVerdict);
}

} // namespace knit
