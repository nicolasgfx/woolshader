#pragma once
// Knitting one mesh component end to end: atlas -> lattice -> placement -> gate (a) -> settle (with
// the binding) -> gate (a) again -> boundary checks -> gate (d) -> verdict. Each stage is a call to
// the file that owns it; this file holds their order and the verdict.
#include <memory>
#include <string>
#include <vector>

#include "knitting/layout/base_surface.h"
#include "knitting/layout/chart_lattice.h"
#include "knitting/checks/checks.h"
#include "knitting/numerics.h"
#include "knitting/relax/settle.h"
#include "knitting/tube/yarn_surface.h"
#include "knitting/layout/uv_atlas.h"
#include "model/mesh.h"

namespace knit {

struct PartKnit {
    std::unique_ptr<UvAtlas> atlas;
    std::unique_ptr<MeshSurface> surface;
    LatticeKnit t;
    ValidationReport before, after;
    // `after` was computed by the CPU validator (the spiral pieces' audit of their GPU curves), not by
    // the knitting's own numerics: a cache restore verifies it with the same validator.
    bool gateAOnCpu = false;
    SettleReport settle;
    BoundaryCheck boundary;
    // Gate (d) on the rendered surface (FinalSurfaceSpec), when requested: one pass at 0.025 mm
    // rings, or two (and the half step 0.0125 mm).
    std::vector<SurfaceCheck> gateD;
    Status gateDStatus = Status::Pass;
    // Rings that needed more than SurfaceParams::maxCompression (-1: not measured).
    int compressionClamped = -1;
    int projectionUnresolved = 0;   // controls whose last continuation was not Owned
    Status verdict = Status::Unresolved;  // gate (a) + lattice + layout + projection + B1 + B2
    std::string detail;
    double msTile = 0, msValidate = 0, msSettle = 0, msBoundary = 0, msGateD = 0;
    // Profiling (--stage-profile): validate() inside settle; gate (d) per pass, model build and check.
    ValidateStats vSettle;
    std::vector<double> msGateDBuild, msGateDCheck;
};

struct MeshKnitParams {
    TileParams tile;
    bool settle = true;
    SettleParams sp;
    HandRequest hand;
    // The final surface (D and c come from tile.pp): inflation / compression / plies; it sets
    // R_s of the footprint inset and of B1, and the surface gate (d) checks it.
    FinalSurfaceSpec surface;
    bool surfaceGate = false, surfaceHalf = false;
    // false: gate (d) is reported but not part of the verdict (a compressed yarn is drawn wider
    // than its knit diameter and overlaps its neighbours by design).
    bool gateDInVerdict = true;
    ExecPolicy exec;  // gate (a), settle's validations (identical results)
};

bool knitMesh(const TriMesh& mesh, const MeshKnitParams& p, KnitNumerics& numerics, PartKnit& out, std::string* err);

// Gate (d) on a knitted part (when p.surfaceGate) and its verdict. knitMesh calls it last; a
// caller trying several diameters may knit with surfaceGate off and call it once on the kept
// attempt (the diameter choice reads gate (a) only, so the kept result is the same).
void runSurfaceGate(PartKnit& out, const MeshKnitParams& p, KnitNumerics& numerics);

} // namespace knit
