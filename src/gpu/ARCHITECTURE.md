# GPU backend: file tree and architecture

The CUDA backend: validation, the knitting solver, surfaces, placement and fibres. This file is the
way into the folder: what each file owns, how data moves through them, and what must stay true when
you change them. The rest of the program uses four public headers, the only headers at this folder's root,
and the compiler holds it: every other header, in the subfolders, is private (`CMakeLists.txt`; the mechanism is described in src/ARCHITECTURE.md).

## File tree

```text
src/gpu/
  ARCHITECTURE.md        this file
  CMakeLists.txt         the component knit_gpu: sources, the public and private headers, CUDA flags;
                         depends on model, knitting and fibres

  the public headers: one function per capability
  device.h               available, deviceSummary, releaseResident; DeviceError and
                         DeviceBudgetExceeded; residency: PriorityLease, ResidencyTenant,
                         memoryAvailability (the live preview).
                         CUDA-free, no project include
  knitting.h             makeKnitNumerics (the KnitNumerics port on the device), the solver's options,
                         statistics and cache identity (ncpSolverIdentity), ncpAvailable, detectComputeCaps
  surfaces.h             the surface model, gate (d), the tube mesh, the plies, the cross-part clearance
  fibre_generator.h      generateFibres and its frame (the CPU twin is fibres/fibre_generator.h)

  device/                what any CUDA program needs; no knit types
  residency.cpp          the priority lease of a GPU job over the tenants (the live preview's renderer
                         process), the DXGI budget reading; its own lock, never the allocator's
  runtime.h, .cu         failures and the checks, the device, the memory budget, DeviceBuffer, leases,
                         shared residents, waiting with cancel, launch pacing
  primitives.h, .cu      scan, sort and select, on tracked scratch memory
  kernel.cuh             what every kernel starts with

  geometry/              yarns and meshes on the device
  device_geometry.h, .cu views, host arrays, owners, flattening, the grid, resident meshes, the uploaded
                         curve set and its spans
  math.cuh               f3 and d3, Bezier spans, arc length, distances, hashing, grid queries
  surface.cuh            rings, surface and mesh queries, the continuation walk
  arc_length.cu          YarnPath::finalize: the arc-length tables of a curve set

  validation/            gate (a)
  checks.h               the batches of gate (a) on one uploaded curve set
  curvature.cu           checkCurvature per span
  clearance.cu           samplePath and checkClearance
  signature.cu           checkLocalSignature

  solver/                the knitting problem
  ncp.h                  solveNcp, NcpDevice: the problem resident on the device
  ncp_control.cu         the iteration as one CUDA graph: conditional WHILE / IF nodes, the controller
                         kernels that decide (augmented Lagrangian over Levenberg-Marquardt steps), solveNcp
  ncp_device.cuh         what the parts share: rows, parameters, counts, the control record, NcpDevice::State
  ncp_setup.cu           one relax call on the device: samples, general rows, capacities, upload, download;
                         the single operations of the tests
  ncp_evaluate.cu        row values, merit, multipliers; the span-maximum curvature row (CurvatureRows, experimental)
  ncp_live.cu            ContactStorage::LiveRows (experimental): the near cutoff, the live list's renumbering, the
                         multipliers' move to their new slots, the key carry-over at a refresh
  ncp_refresh.cu         in run order: sample points, motion, broad phase, contact rows, active rows
  ncp_step.cu            assembly (gradient, diagonal blocks), trial, swept test, pinning, KKT
  ncp_pcg.cu             matrix-free conjugate gradients with a 3x3 block-Jacobi preconditioner

  surface/               surfaces and placement
  placement.h            the placement's surface queries, the lattice and the boundary checks
  surface_build.cu       in run order: the surface model, its index, compression, gate (d), the tube mesh
  plies.cu               plies from a yarn's surface
  piece_surface.cu       PieceSurface queries in batches
  mesh_surface.cu        MeshSurface queries in batches, boundary checks B1 and B2
  lattice.cu             the tiled lattice's node classification
  scene_queries.cu       queries over a whole scene's parts: cross-part clearance

  fibres/
  fibres.h               the flat parameters and the records of the generator
  fibre_candidate.cuh    one candidate on the device
  fibres.cu              the kernels, the batched run, generateFibres

  port/                  the knitting port at the top
  knit_numerics.cpp      GpuKnitNumerics, KnitNumerics on the device: gate (a), relax, store, surfaces,
                         lattice, gate (d), boundary checks; the tuned constants; ncpAvailable;
                         replayNcpCall (one captured relax call solved or sized again)
  ncp_call_file.h, .cpp  the exact input of one relax call as a file (NcpNumericsOptions::captureDir)
```

**Reading order:** `device.h`, `knitting.h`, `surfaces.h`, `fibre_generator.h`, `device/runtime.h`, `geometry/device_geometry.h`, `geometry/math.cuh`,
then `surface/plies.cu` as the shortest complete example of the whole recipe below. After that any
capability can be read on its own; the reading paths below follow one request each.

## The common recipe

Well-built GPU components share one shape, whichever library they use:

1. **A runtime owns the device and its memory.** Everything else asks it for memory (Kokkos memory
   spaces, RMM memory resources).
2. **Buffers own; plain views travel.** A kernel receives a struct of device pointers by value and
   owns nothing.
3. **A kernel is a pure function over views,** one thread per item, writing only that item's slot.
4. **Variable-length lists are built in three steps:** count per item, turn the counts into offsets
   with a scan, fill. No list is ever guessed and grown inside a kernel.
5. **The host decides; the data stays.** Control flow lives on the host and reads small summaries.
   Whole arrays cross only at the start and the end. The exception is the knitting solver's
   iteration: there the device decides too (below), so no read on the host stands between two steps.
6. **One failure type, every call checked.** A CUDA call's status goes through `check()` and becomes a
   `DeviceError` (libcudf's `cuda_error`); the public entry points turn it into `false` and a reason.

## What this project adds

| Need | Where |
|---|---|
| Several jobs share one device, so N jobs never hold N uncontrolled uploads | `Reservation` (a FIFO lease on the budget) |
| An overcommitted GPU pages to system memory and runs 10-30x slower on Windows | `AllocationCap`, `solveCapBytes` |
| A launch longer than the display watchdog (2 s) resets the device | `LaunchPacer`: batches sized from the last batch's measured time |
| Cancel must reach a kernel that runs for seconds | the abort flag in host-mapped memory, polled by long loops |
| Reruns should agree, and floating-point atomics depend on the order | fixed-grid sum reductions, integer atomics on bit patterns, radix sorts |
| Meshes are shared by the pieces of a scene | `residentMesh`, kept across the jobs of one scene and counted in the budget; the job queue calls `releaseResident` when a job starts on another scene |
| A yarn is 52 m long, and float would quantize arc length to 4 um | the precision rules below |
| Verdicts must stay honest | GPU results are diagnostic float32 verdicts; a failure never falls back to the CPU |
| A failure anywhere must leave nothing behind | RAII owners and leases; `gpu_injected_failures_release_everything` fails each entry point at its checks |

## Layers

An arrow means "is used by". The includes run one way, and `architecture_gpu_layers` holds them.

```text
device ──> geometry ──┬──> solver ─────┐
                      ├──> validation ─┤
                      ├──> surface ────┼──> port
                      └──> fibres
```

| Folder | Includes folders | Includes public headers |
|---|---|---|
| `device/` | device | `device.h` |
| `geometry/` | device, geometry | `device.h` |
| `solver/` | device, geometry, solver | `device.h` |
| `validation/` | device, geometry, validation | `device.h` |
| `surface/` | device, geometry, surface | `device.h`, `surfaces.h` |
| `fibres/` | device, geometry, fibres | `device.h`, `fibre_generator.h` |
| `port/` | device, geometry, solver, validation, surface | `device.h`, `knitting.h`, `surfaces.h` |

- `device/` and `device.h`, which everything includes, include no knit header except `model/cancel.h`.
- `port/` composes the capabilities; none of them knows another, and none calls the port.
- Everything outside this folder includes only the four public headers; the others are not on its include
  path. Each user includes the one it needs, so a change to the fibre generator's interface does not
  recompile the knitting stages.
- `knit_knitting` and `knit_fibres` never name this folder. Knitting reaches the device through one port, `KnitNumerics`
  (`src/knitting/numerics.h`); `GpuKnitNumerics` implements it, and each worker owns one object, so its
  placement statistics take no lock.

## The data

Every concept exists in up to three forms, which is the common recipe's second step plus a mirror:

| Form | Example | Lives on | Role |
|---|---|---|---|
| View | `CurveSetView`, `MeshView`, `GridView` | host and device | device pointers and sizes, passed to kernels by value |
| Arrays | `CurveArrays`, `MeshArrays`, `BoxArrays` | host | the bytes the device receives; a run sizes its reservation from them before it uploads |
| Owner | `DeviceCurves`, `DeviceMesh`, `DeviceGrid` | host | owns the `DeviceBuffer`s behind one view |

**Tables per control.** A value that varies along a yarn slowly sits in a flat array parallel to `ctrl`,
at the same `ctrlOff`, not in `RingRecord`: rings are rebuilt per consumer and are about 50 times as
many. Thickness is read once, by `ringKernel`, and island weights by both fibre kernels. Neither exists
while unused. A metre of yarn holds 3.2 MB of rings, 3.4 KB of thickness and 10 KB of weights.

**Precision.** 3-D geometry, frames, distances and local integrals are float. One-dimensional yarn
coordinates that grow with the yarn (the arc length s, the spline parameter u, arc-length tables) are
double, and u is split into a span and a float local t before any geometric evaluation. The solver's
rows, products and reductions are double, because its decisions compare merits to about 1e-12; its
conjugate gradients iterate in double.

**Results are compiled code.** nvcc contracts multiplies and adds into fused ones (no fast math, but
`--fmad=true`), and its choices follow the code around an expression: moving float arithmetic between
functions, or removing a branch that never runs, can change results in the last bit. A refactoring of
kernel code is checked field by field against the code before it, not by tolerance.

## One request, end to end

Gate (d), `checkSurface` in `surface/surface_build.cu`, shows every step:

```text
checkSurface
├─ available                      no device: the reason, never a CPU stand-in
├─ Reservation                    wait for the turn on the device
├─ DeviceSurfaceModel             the model, top to bottom in the file:
│  ├─ planRings                   curves to CurveArrays to DeviceCurves, ring and yarn records uploaded
│  ├─ computeRings                ringKernel: rings in parallel
│  ├─ propagateFrames             chunkMapKernel, chunkStartKernel, chunkPropagateKernel: frames in ordered chunks
│  └─ compressAxes                with compression: a capsule index, compressReqKernel, compressApplyKernel
├─ DeviceCapsuleIndex             capsuleKernel and a DeviceGrid: the surface index
├─ checkFolds, checkCollisions    foldKernel, collisionKernel: per-item results
├─ reduce kernels                 one key per check: the worst item, by integer atomics
└─ download                       a few numbers; the host computes the verdict
```

The solver follows the same shape with the loop on the device: `ncp_control.cu` captures the
iteration once into a CUDA graph whose WHILE / IF nodes are driven by one-thread controller kernels
(`cudaGraphSetConditional`); the host launches it, waits, and reads the control record. Sizes that change
(contact rows, active rows) live in `NcpCounts` on the device; kernels over them run a fixed grid and loop
to the device count; capacities are sized once from the starting state, and an overflow ends the launch
with an explicit status instead of truncating. A contact row stores only its pair, type, multiplier,
values and terms; controls, keys and gradients are derived where they are needed (`ncp_device.cuh`). When
a refresh outgrows the contact rows, `solveNcp` grows the capacity within the memory cap and launches the
graph again with `NcpControl::resume`, continuing the same iteration (the host's only decision).

## Reading paths

**Gate (a), one validation.** `GpuKnitNumerics::validate` (`port/knit_numerics.cpp`): a `Reservation`;
`uploadCurves` (geometry) flattens and uploads the controls, tables and prefix sums once. Then the checks of
`validation/`:
- `checkCurvatureOnDevice`: `curvKernel` per span, `waitForKernels` with the abort flag, the per-span results down.
- `samplePolylinesOnDevice`: `countVerticesKernel`, offsets, `fillVerticesKernel`, `firstVertexKernel`, `arcKernel`;
  the polylines down for the report, kept on the device for the clearance.
- `checkSignatureOnDevice`: the host builds the items; `transformKernel`, `findCrossingsKernel` (paced, with
  `waitForKernels`), a sort, `gatherKernel`, `evaluateSignatureKernel`; one outcome per item down. The host
  then applies the crossing attribution (CPU).
- `checkClearanceOnDevice`: a grid of chunks, `pairKernel`, `workKernel` (`record`, `sameYarnPair`), a
  second pass if the issues overflow; counters and issues down, sorted in the CPU's order. The host
  summarizes and decides the verdict.

**One relax call.** `GpuKnitNumerics::relax`: the problem's constants, a `Reservation` and an
`AllocationCap`; `solveNcp` (`solver/ncp_control.cu`):
- the `NcpDevice` constructor (`ncp_setup.cu`) uploads samples and general rows, fixes the broad phase's
  grid, counts the starting state's rows (`startCountsKernel`: the only counts read before the iteration)
  and sizes the contact capacity to fit the memory cap (`chooseCapacity`);
- `SolveGraph::build` captures the enqueue sequences of `ncp_refresh.cu`, `ncp_evaluate.cu`,
  `ncp_step.cu` and `ncp_pcg.cu` into the bodies of the graph's conditional nodes, around the controller
  kernels (`innerHeadKernel`, `afterAssembleKernel`, `trialDecideKernel`, `attemptEndKernel`,
  `outerDecideKernel`, ...), which hold the stopping rules and the early stops (`earlyStopRule`);
- one `cudaGraphLaunch`; `waitForStream` (cancel raises the abort flag, which `innerHeadKernel` reads); on
  a contact-row overflow `regrowContactRows` and a resumed launch (`LaunchKind::Resume`);
- the control record and `downloadResult`: the status, the statistics, the controls back.

**One surface model with gate (d).** The tree above.

**One placement batch.** `GpuKnitNumerics::pieceQuery` → `makePieceQuery` (`surface/piece_surface.cu`)
uploads the chart once per object; each batch uploads its points, runs `atKernel` or `projectKernel`
(one thread per query, the mesh's `DeviceGrid`) and downloads the results; `PlacementStats` counts them.

**One fibre run.** `generateFibres` (`fibres/fibres.cu`): `checkRequest`; `flattenInputs` (or the kept
inputs), `flatParameters`, `listStrands`; a `Reservation` sized from the inputs and batches; `runFibres`
uploads (`DeviceFibreInputs`, resident meshes) and synchronizes. `FibreBatches::run` then works on its own
non-blocking stream (`TimedStream`), checking for cancel before each batch:
- `placeCandidates`: `placeCandidatesKernel` (kernel A: anchor, front, hidden, ply side, exposure,
  thinning), a CUB select of the survivors;
- `growBatch`: `growFibresKernel` (kernel B: `growFibre`, the device twin of `growCandidate`: the same
  stages under the same names, in its order, on a `GrowingFibre`; finished by `finishLight` or `finishFull`),
  a scan, `packPointsKernel`, the compaction;
- each section ends with `TimedStream::finish` (an event wait); `LaunchPacer` sizes the next batch.

The records and points come down; `mergeFibres` builds the `FibreSet` and the report on the host.

## What must stay true

| Rule | Held by |
|---|---|
| `ncpSolverIdentity()` never changes by accident: it is hashed into every knit cache key | `gpu_solver_identity_text` |
| A grid's buckets list items in increasing id order: the order decides ties | `gpu_grid_matches_reference` |
| Every GPU result is compared with the CPU's where a CPU version exists | `gpu_*_gpu_vs_cpu`, `gpu_ncp_gpu_rows_match_cpu_*`, `gpu_fibres_properties_and_cpu_agreement` |
| The includes run one way, the public headers included as the layers say, `device.h` free of knit types | `architecture_gpu_layers` |
| This file names every file here | `architecture_gpu_doc_complete` |
| No templates and no generic lambdas: write the concrete function, or a plain loop such as `GridCursor` | `architecture_gpu_no_templates` |
| Every CUDA call is checked, handled on purpose, or cleanup that never throws | `architecture_gpu_cuda_calls_checked` |
| The solver's iteration reads nothing on the host: one graph launch per relax call, every decision in a controller kernel | `gpu_ncp_gpu_*` (results) |
| A device failure is a `DeviceError` below the public headers, and `false` with a reason at them; it frees what it took and releases its lease | `gpu_device_failure_is_false_and_err`, `gpu_injected_failures_release_everything` |
| The kernels' outcomes are the CPU twin's `FibreOutcome` | the type itself (`fibres/fibres.h` `outcomeCode`) |
| A kernel or a function with more than about five parameters takes a struct of what they are | review |

## Where to change what

| To change | Touch |
|---|---|
| A solver constant | `port/knit_numerics.cpp` or `solver/ncp_control.cu`. If it changes results, bump the revision in `ncpSolverIdentity`. |
| A row of the knitting problem | `tests/ncp_reference.*` and its device twin in `solver/ncp_evaluate.cu`; `gpu_ncp_gpu_rows_match_cpu_*` must agree |
| A fibre rule | `src/fibres/fibre_generator.cpp` and `fibres/fibres.cu` together |
| A new check on the surface | a kernel and a reduction in `surface/surface_build.cu`, with a CPU agreement test |
| A gate (a) check | `knitting/checks/` and its twin in `validation/`, with a CPU agreement test |
| A new kernel | see below |

## Adding a kernel

A kernel handles one item per thread, and its launch names how many:

```cpp
// The sample's position on the device, one thread per sample.
__global__ void samplePointsKernel(SpanList spans, const double* u, float* points) {
    const int g = threadItem();          // this thread's item
    if (g >= spans.count) return;        // the last block has spare threads
    ...
}

samplePointsKernel<<<blocksFor(n), kThreadsPerBlock>>>(spans, u.doubles(), points.floats());
checkLaunch("sample points");            // a launch that failed reports here
```

1. Read views, write only your item's slot. Five parameters at most: group the rest in a struct.
2. Reduce with integer atomics or partial sums on the fixed `kSumBlocks` grid, never with floating-point atomics.
3. Size the buffers before the run's `Reservation`. A list of unknown length is count, offsets, fill.
4. A launch that can run for seconds goes through a `LaunchPacer` and `waitForKernels`.
5. Every CUDA call around it goes through `check()`.
6. Add a test against the CPU version where one exists.

## Literature

- Blelloch 1990, *Vector Models for Data-Parallel Computing*: scan as a primitive, the count, offsets,
  fill idiom.
- Ericson 2004, *Real-Time Collision Detection*: segment and triangle distances (5.1.5, 5.1.9).
- Madsen, Nielsen, Tingleff 2004, *Methods for Non-Linear Least Squares Problems*: the damping update.
- CCCL (Thrust, CUB), Kokkos and RMM: the primitives and the separation of memory from execution;
  libcudf's error utilities: one failure type at the boundary.
