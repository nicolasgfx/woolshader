#pragma once
// The one port between knitting and the compute behind it. knit_knitting owns this interface; knit_gpu
// implements it (GpuKnitNumerics) and the tests implement a CPU double (tests/cpu_numerics.h). The
// orchestration (rounds, free sets, rollback, repair stages, verdicts) is host code on this side;
// every numerical batch goes through here. Errors are thrown (std::runtime_error, or PartFailure
// below); Cancelled propagates.
//
// This header includes only src/model and names every knitting type by forward declaration, so any
// part of knitting may include it without depending on another part.
#include <memory>
#include <stdexcept>
#include <vector>

#include "model/mesh.h"
#include "model/parallel.h"
#include "model/params.h"
#include "model/yarn_path.h"

namespace knit {

// A failure that ends one part, not the run: the part is reported as not knitted.
struct PartFailure : std::runtime_error {
    using std::runtime_error::runtime_error;
};

class BaseSurface;
class Chart;
class MeshSurface;
class UvAtlas;
struct Allowances;
struct BoundaryCheck;
struct CurvatureResult;
struct LatticeChartQuery;
struct LatticeNodeClass;
struct LatticeSolids;
struct SurfaceBinding;
struct SurfaceCheck;
struct SurfaceParams;
struct ValidateOptions;
struct ValidationReport;

using Controls = std::vector<std::vector<Vec3>>;
using ControlMask = std::vector<std::vector<char>>;

// Every yarn's control points, in yarn order.
inline Controls controlsOf(const CurveSet& cs) {
    Controls P;
    P.reserve(cs.yarns.size());
    for (const YarnPath& y : cs.yarns) P.push_back(y.curve.ctrl);
    return P;
}

// Who is asking for a relaxation, for the solver's hard thresholds and the --stage-profile lines.
struct RelaxCall {
    const char* phase = "";    // "settle" | "repair"; "": outside both, the move's thresholds are hard
    int round = -1, stage = -1;
    PatternParams validatePP;  // the validation that judges this move: settle's pp, repair's current D
};

struct RelaxMoveParams {
    double thr = 1.15;        // contact distance [mm]
    double kappaMax = 1.0;    // bend limit [1/mm]
    double window = 2.0;      // same-yarn arc window for contacts [mm]
    int iterations = 120, refresh = 10, samplesPerSpan = 4;
    double lengthStiffness = 0.3;
    double bandHalf = 2.0;
    double maxDisplacement = 0;  // > 0: no free control moves farther than this from its start
    ExecPolicy exec;             // the refresh (positions, candidates, band) in parallel; identical
    RelaxCall call;
};

// One relaxation: the projections on the free controls of P (bend limiter, soft edge lengths with
// rest lengths from Prest, contacts, band, displacement cap).
struct RelaxRequest {
    Controls& P;
    const Controls& Prest;
    const ControlMask& freeCtrl;
    const CurveSet& cs;
    const BaseSurface* surface;  // the band's; null: no band
    const RelaxMoveParams& mp;
    SurfaceBinding* bind;        // the band by ownership-preserving projection (base_surface.h); null: by surface
};

class KnitNumerics {
public:
    virtual ~KnitNumerics() = default;
    virtual const char* name() const = 0;
    // Gate (a): curvature, sampling, local signature / crossings, clearance (validate in checks.h);
    // checkCurvature of one yarn.
    virtual ValidationReport validate(const CurveSet& cs, const PatternParams& pp, const ValidateOptions& opt) = 0;
    virtual CurvatureResult curvature(const YarnPath& y, double radius, const Allowances& al, const ExecPolicy& exec) = 0;
    // Settle and repair: moves the controls of P; stores P into cs with every yarn's arc-length
    // tables (YarnPath::finalize).
    virtual void relax(const RelaxRequest& request) = 0;
    virtual void store(CurveSet& cs, const Controls& P, const ExecPolicy& exec) = 0;
    // Surface queries: the mesh a part is placed on; a piece's batched chart queries (PieceSurface);
    // a chart's candidate lattice nodes in one batch.
    virtual std::unique_ptr<MeshSurface> meshSurface(const TriMesh& mesh) = 0;
    virtual std::unique_ptr<BaseSurface> pieceQuery(const Chart& chart) = 0;
    virtual void classifyLattice(const UvAtlas& atlas, const LatticeChartQuery& q, const LatticeSolids& solids,
                                 std::vector<LatticeNodeClass>& out) = 0;
    // Gates on finished yarn: gate (d) on the rendered surface; the boundary checks B1 and B2.
    virtual SurfaceCheck checkRenderedSurface(const CurveSet& cs, const PatternParams& ppR, const SurfaceParams& sp) = 0;
    virtual BoundaryCheck checkBoundary(const CurveSet& cs, const SurfaceBinding& bind, const MeshSurface& surf, double Rs) = 0;
};

} // namespace knit
