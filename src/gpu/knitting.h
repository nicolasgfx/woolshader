#pragma once
// The knitting port (knitting/numerics.h) on the device, and the solver behind KnitNumerics::relax.
// What the GPU results are: the existing algorithms and settings evaluated in float with the
// documented scheduling changes. Their verdicts are GPU diagnostic verdicts, not the CPU's certified
// PASS / FAIL / UNRESOLVED; deviations from the CPU are measured by the tests, never filtered or replaced
// by CPU results.
#include <memory>
#include <string>
#include <vector>

#include "knitting/numerics.h"
#include "knitting/relax/knit_problem.h"
#include "model/compute_caps.h"

namespace knit::gpu {

// The contact storage the environment selects (ContactStorage, knitting/relax/knit_problem.h):
// KNIT_NCP_CONTACT_STORAGE=live gives LiveRows, anything else (or nothing) CompactRows. The app's runtime switch for
// the experimental storage; the default is unchanged.
ContactStorage contactStorageFromEnvironment();
// The curvature rows the environment selects (CurvatureRows, knit_problem.h): KNIT_NCP_CURVATURE_ROWS=span-max gives
// SamplesAndSpanMaximum, anything else (or nothing) Samples.
CurvatureRows curvatureRowsFromEnvironment();

// The knitting solver: KnitNumerics::relax is solveNcp on the problem
// of the call (hard thresholds at the call's validation parameters, RelaxMoveParams::call, plus
// margins; the move's thresholds as targets). Each solve holds the whole device
// (Reservation::exclusiveBytes). The tuned constants are in knit_numerics.cpp.
struct NcpNumericsOptions {
    ContactStorage contactStorage = contactStorageFromEnvironment();
    CurvatureRows curvatureRows = curvatureRowsFromEnvironment();
    SlowCall slowCall = SlowCall::StopWhenSlow;   // knit_problem.h; not in the environment
    // Device memory cap of a solve [MiB]: the smaller of this and solveCapBytes (KNIT_GPU_BUDGET_MB, else
    // the tracked memory plus the device's free memory at the call minus a reserve); 0: solveCapBytes
    // alone. A solve that would exceed it is refused (DeviceBudgetExceeded): with 8 parts knitting side by
    // side the dense man head (330k unknowns) paged the device and ran ~30x slower.
    double memoryCapMB = 0;
    // Diagnostics: when set, every relax call first writes its exact solver input to
    // <captureDir>/call_<n>.ncpcall (n = 1, 2, ... per numerics object), which replayNcpCall reads.
    std::string captureDir;
};
// The port (numerics.h) on the device: gate (a), the relaxation of the explicit problem, arc-length
// tables, the surfaces, the lattice classification, gate (d) and the boundary checks. Its validations
// are GPU diagnostic verdicts (ValidationReport::gpuDiagnostic). null with `err` without a device.
// One object per run: its placement statistics take no lock.
// The relaxation needs ncpAvailable(); without it a relax() call fails with its reason.
std::unique_ptr<KnitNumerics> makeKnitNumerics(std::string* err, const NcpNumericsOptions& o = {});
// Two lines: the relax calls, and the batches (timings, counts).
std::string knitNumericsStats(const KnitNumerics& n);
// One line: the placement's surface queries so far.
std::string knitPlacementStats(const KnitNumerics& n);
// The ncp solver's identity in the cache keys (teddy knit cache, scene stage keys): the revision
// and the defaults that change its output. Bump the revision with every other change of behaviour.
// The memory cap is not part of it: a refused part is never cached (KnitStop::DeviceBudget). The contact
// storage and the curvature rows are: ncpSolverIdentity() names the ones the environment selects
// (contactStorageFromEnvironment, curvatureRowsFromEnvironment), the overload the given ones; CompactRows with
// Samples keeps the v5 text.
std::string ncpSolverIdentity();
std::string ncpSolverIdentity(ContactStorage storage, CurvatureRows curvature);
// One line per relax call so far (phase, round, size, LM steps, PCG iterations, ms, status): the stage profile.
std::vector<std::string> ncpNumericsCalls(const KnitNumerics& n);

// One captured relax call (NcpNumericsOptions::captureDir) run again on its own: the whole solve, or
// only its starting sizes (no contact buffers: how large the call is, and what it would allocate).
// Frozen: two linearizations at the call's start held fixed (NcpDevice::frozenLinearizations), written
// to `frozenOut` for a comparison of layouts or with a reference.
enum class NcpReplayMode { Solve, CountOnly, Frozen };
struct NcpReplayRequest {
    std::string callFile;
    NcpReplayMode mode = NcpReplayMode::Solve;
    ContactStorage contactStorage = ContactStorage::CompactRows;  // not part of the call file: chosen per replay
    CurvatureRows curvatureRows = CurvatureRows::Samples;         // not part of the call file: chosen per replay
    double memoryCapMB = 0;   // as NcpNumericsOptions::memoryCapMB
    std::string controlsOut;  // Solve: the final controls as raw doubles (x y z per control, yarn order); "": not written
    std::string frozenOut;    // Frozen: the passes' merits, gradients, diagonal blocks and products (raw doubles)
    double frozenRho = 10;    // Frozen: the penalty of both passes
    long long maxLmSteps = 0; // Solve: > 0 replaces the call's step budget (how two solves part, step by step)
};
struct NcpReplayResult {
    bool refused = false;      // the cap refused the solve (DeviceBudgetExceeded); `status` says why
    std::string status;        // the solver's status, or the refusal
    std::string callLine;      // phase, size, LM steps, PCG iterations, ms (as ncpNumericsCalls)
    std::string statsLine;     // NcpStats::line
    std::string memoryLine;    // the memory report (buffer groups, lists, capacities, peaks)
    unsigned long long controlsDigest = 0;  // FNV-1a 64 of the final controls' bits (Solve)
    double capMB = 0, ms = 0;
};
// false with `err` when the file cannot be read or the device fails; a refusal by the cap is a result
// (refused, true).
[[nodiscard]] bool replayNcpCall(const NcpReplayRequest& request, NcpReplayResult& result, std::string* err);

// The knitting solver can run: a usable device (the solver needs no library beyond the CUDA runtime).
// The reason otherwise.
[[nodiscard]] bool ncpAvailable(std::string* why = nullptr);
// available (device.h), ncpAvailable and deviceSummary in one value; the callers cache it
// (jobs/jobs.h).
ComputeCaps detectComputeCaps();

} // namespace knit::gpu
