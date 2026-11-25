#pragma once
// The exact input of one relax call as a file, for diagnostics: what GpuKnitNumerics::relax hands to
// solveNcp, written before the solve (NcpNumericsOptions::captureDir) and read back by replayNcpCall
// (knitting.h), so one large call can be measured and compared without running its whole pipeline.
// Version 1 holds the curve set's controls and arc-length tables, the controls, rest controls and free
// mask, the band's triangle mesh (a MeshSurface, or no band), the move parameters, and the derived
// problem and stopping parameters, followed by an FNV-1a 64 digest of everything before it.
// A binding or another surface type is refused: no capture needs them yet.
#include <string>

#include "knitting/numerics.h"
#include "knitting/relax/knit_problem.h"
#include "model/mesh.h"
#include "model/yarn_path.h"

namespace knit::gpu {

struct NcpCallRecord {
    CurveSet curves;  // per yarn: the controls and the arc-length tables (all the solver reads)
    Controls P, Prest;
    ControlMask freeCtrl;
    bool hasMesh = false;  // the band's surface is a MeshSurface of `mesh`
    TriMesh mesh;          // positions, triangles, vertex normals
    RelaxMoveParams move;
    std::string phase;  // move.call.phase (the record owns the text)
    KnitProblemParams problem;
    NcpSolveParams solve;
};

// false with `err` when the request's surface or binding is not supported, or the file cannot be written.
[[nodiscard]] bool writeNcpCall(const std::string& path, const RelaxRequest& request, const KnitProblemParams& problem, const NcpSolveParams& solve,
                                std::string* err);
// false with `err` for a missing, truncated, foreign or corrupted file (digest mismatch).
[[nodiscard]] bool readNcpCall(const std::string& path, NcpCallRecord& out, std::string* err);

} // namespace knit::gpu
