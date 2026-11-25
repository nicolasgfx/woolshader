# Knitting: file tree and architecture

How a triangle mesh becomes yarn: which stitches exist, where they sit, how the yarn runs through
them, and whether the result is really knitting. This file is the way into the folder. The kernels
are in `src/gpu`.

## From a stitch to a yarn tube

A **stitch** is a loop of yarn pulled through the loop of the row below. One stitch is a template of
16 control points of a cubic B-spline, and a row of stitches is one continuous curve. A column of
stitches is a *wale*, a row is a *course*, and the gauge `w` by `h` is the spacing of the two in mm.

Six steps turn a surface into yarn:

1. **The program** says which stitches exist and who is whose parent: knit (K), increase (INC),
   decrease (DEC). It holds no geometry.
2. **The layout** gives each stitch a node on the surface and a frame around it.
3. **The placement** maps the template through each frame and projects the controls onto the surface.
   The result is a `CurveSet`: the centrelines of the yarns.
4. **The check** decides whether those curves are valid knitting.
5. **The relaxation** moves the controls of the stitches that failed, then checks again.
6. **The tube** puts rings around each centreline. That surface is what the renderer draws and what
   the fibres grow on.

There are two layouts. **Mesh charts** cover each UV chart with a regular lattice of stitches, one open
yarn per run: imports and every preset but the teddy. **Spiral pieces** knit a sphere-like piece in the
round from pole to pole, as one yarn with increases and decreases: the teddy's sewn pieces, whose (a, b)
chart is modelled. Imports knitted their closed toy pieces this way on a computed chart until 2026-11-11.

## File tree

```text
src/knitting/
  ARCHITECTURE.md       this file
  CMakeLists.txt        the component knit_knitting (src/ARCHITECTURE.md): every header public; depends on
                        model only
  numerics.h            KnitNumerics, the port to the compute device

  stitch/               what a stitch is: no surface, no scene
    stockinette.h, .cpp       the 16-point template and the flat patch generator
    expected_signature.h      the crossings a correct stitch must show
    stitch_program.h, .cpp    turns, K / INC / DEC, the program and its check
    hand_tension.h, .cpp      a knitter's hand: every stitch's tension, raised loops, the yarn's thickness

  checks/               gate (a): is this centreline valid knitting?
    interval.h          interval arithmetic for the certified bounds
    checks.h            every result type and check function; flatPatternReference (its result type
                        PatternReference is in src/model)
    curvature.cpp       certified speed and curvature bounds
    crossings.cpp       the projected crossings, each proved by a Krawczyk test
    clearance.cpp       how close the yarns come
    signature.cpp       global and local crossing signatures, crossing attribution
    validate.cpp        runs the checks: verdict, report, its exact text, flat swatches

  layout/               where stitches go on a surface
    base_surface.h, .cpp      BaseSurface, cone, plane, MeshSurface, the ownership-preserving projection
    chart.h, .cpp             the (a, b) chart of a sphere-like piece, and how far its columns are from perpendicular
    uv_atlas.h, .cpp          UV charts with a metric, and which triangles are usable
    stitch_layout.h, .cpp     nodes, frames and yarn placement: placeSpiral, placeRows
    chart_lattice.h, .cpp     the stitch lattice on UV charts: tileMesh

  relax/                moving the yarn until it is valid
    knit_problem.h, .cpp      the constrained problem's types, shared with the GPU solver
    settle.h, .cpp            settle and the certified repair

  tube/                 the yarn as a solid
    compression_fit.h         the cross-section fit, compiled by host and device
    yarn_surface.h, .cpp      rings, tube mesh, SurfaceIndex, gate (d), FinalSurfaceSpec
    plies.h, .cpp             explicit twisted plies
    ply_profile.h, .cpp       plies as shading maps
    param_text.h, .cpp        parameters and curves as exact text: the inputs of cache keys

  knit/                 the use cases, end to end
    mesh_knit.h, .cpp         knit one mesh component; PartKnit, the common result of a part
    piece_knit.h, .cpp        knit one piece as a spiral
    pieces.h, .cpp            the pieces manifest, gauge, loading a piece: the one file that reads files
    coverage.h, .cpp          how much of the surface the yarn covers; a flat swatch's open area
    fixtures.h, .cpp          validation presets, shaping and strand demos
```

**Reading order:** this file, `stitch/stockinette.h`, `stitch/stitch_program.h`,
`layout/stitch_layout.h`, `knit/piece_knit.cpp`, `checks/checks.h`, `relax/settle.h`, then
`layout/chart_lattice.h`, `knit/mesh_knit.cpp`, `tube/yarn_surface.h`.

## Layers

The includes run one way. An arrow means "may include".

```text
                 ┌──> checks ──┐
knit ──> relax ──┤             ├──> stitch ──> model
                 ├──> layout ──┘                 ^
                 └──> tube ──────────────────────┘

numerics.h sits at the root. Any layer may include it. It includes only src/model.
```

A folder may also include what its targets may, so `knit` may include every folder here and `relax`
everything but `knit`. Nothing here includes anything outside `src/model` and `src/knitting`.

- `stitch` knows no surface, so everything else may depend on it.
- `checks` and `layout` do not see each other. A check judges curves, however they were placed.
- `tube` includes only `model`. That is why the fibres can use the yarn's surface without knowing how
  a stitch is laid out.
- `relax` moves curves on a surface and judges them, so it sits above `checks` and `layout`.
- `knit` holds the order of the stages and the verdict.
- `numerics.h` names every knitting type by forward declaration, so including it adds no edge
  between the folders. That is what lets `layout` and `relax` call the device.

## How a part is knitted

The scene layer picks the layout and tries the diameters. Each attempt runs the road below. Both
layouts take the same steps in the same order: placement, the knitter's hand, gate (a), settle with the
hand's rest shape.

### On mesh charts

`scene/knit/mesh_part.cpp`, `knitMeshPart`, tries each diameter until gate (a) passes.

```text
knitMesh                            knit/mesh_knit.cpp
├─ FinalSurfaceSpec::valid          the rendered surface; its support radius sets the footprint inset and B1
├─ UvAtlas                          layout/uv_atlas.cpp: joinCharts, measureTriangles, then per chart chartMetric,
│                                   buildChartGrid, chartValidity, chartDistortion, chartOrientation, collectBoundaries
├─ numerics.meshSurface             a MeshSurface with the device's batched queries
├─ tileMesh                         layout/chart_lattice.cpp
│    buildRuns                      per chart: candidate nodes, numerics.classifyLattice, runs of kept nodes
│                                   (a gap patch: walkPatch from a seed instead)
│    buildStitchProgram             parents, children and pattern face from the lattice coordinates
│    layOutNodes                    nodes and frames on the surface (computeFrames)
│    placeYarn                      placeRows, the SurfaceBinding, numerics.store
│    checkLatticeProgram, checkLayout, recordChartCoverage
├─ knitByHand                       stitch/hand_tension.cpp, when the hand is on
├─ numerics.validate                gate (a)
├─ settle                           relax/settle.cpp, keeping the binding, resting on the hand's rest shape
├─ numerics.checkBoundary           B1, B2
└─ runSurfaceGate                   gate (d) on the rendered surface, then partVerdict
```

`knitMeshPart` runs `knitMesh` with the surface gate off and calls `runSurfaceGate` once, on the kept
attempt: the diameter choice reads gate (a) only.

### As spiral pieces

The chart comes first, from the scene: a teddy piece carries its modelled chart as uv. `knitSpiralPart`
(`scene/knit/spiral_part.cpp`) then runs:

```text
chartPiece                          the chart's check and skew, in the piece's log
knit cache, or knitWithDiameterFallback, which per diameter runs knitAttempt:
  knitScenePiece                    knit/piece_knit.cpp
  ├─ prepareScenePiece              Chart of the piece in memory, numerics.pieceQuery, PieceSurface, the pole tests
  ├─ the poles                      start where the end pole hides in another piece; the end gathered where it shows,
  │                                 threaded inside another piece
  ├─ knitPiece
  │    measureRows                  stitch counts from row lengths: scheduleCounts, buildProgram
  │    knittingDirection            which way round b runs, so that c x v points outwards
  │    planTurns                    where each turn runs in b
  │    placeTurns                   turnsOnChart, helicalLayout (and the virtual row), placeSpiral, numerics.store
  │    closeLastTurn                CAST at the start pole, CLOSE (threaded or gathered) at the end pole
  └─ knitByHand                     stitch/hand_tension.cpp, when the hand is on
  checkProgram, checkLayout
  numerics.validate                 gate (a), float32
  settle                            relax/settle.cpp, resting on the hand's rest shape
  certifiedRepair                   relax/settle.cpp, if gate (a) still fails
then once, on the kept attempt:
  validate                          the CPU audit of the GPU curves, double precision (not for a cached piece)
  savePieceToCache                  the knit cache entry
  gate (d)                          on the rendered surface: gpu::checkSurface, called directly, not through the port
```

`spiralKnit` copies the kept attempt into a `PartKnit`, so everything downstream sees one kind of part.

### The tube

`buildSurface` makes the rings: `sampleRings`, `transportFrames`, then with compression
`collectContactRequirements`, `smoothRequirements` (with smoothing), `fitSections`. `buildTubeMesh` turns
them into triangles, `buildPlies` into plies. `measureCoverage` then reports how much of the mesh the tubes
hide; the overlap mode's gap patches go into its holes.

## The gates

A gate is a check the result must pass before anything downstream may trust it. Every check returns
`Status`, in three values:

- **PASS**: the claim holds after subtracting every allowance.
- **FAIL**: the claim is broken after adding every allowance.
- **UNRESOLVED**: the result lies within the allowances, or a precondition of the method was not met.

FAIL wins over UNRESOLVED, which wins over PASS. UNRESOLVED never counts as a pass.

**Gate (a) asks: are these curves knitting?** It looks at centrelines only.

- The yarn bends no tighter than its radius (`curvature.cpp`, certified per span).
- The crossings are the crossings of a correct stitch (`crossings.cpp`, `signature.cpp`): the yarn
  goes through the loop of its parent, and not past it. Flat patches compare one signature;
  surfaces compare each stitch in its own frame, the *local signature*.
- Different yarns, and distant parts of one yarn, stay at least `D + c` apart (`clearance.cpp`), with
  `D` the yarn diameter and `c` the clearance.

**Gate (d) asks: is the surface we draw sound?** `checkSurface` looks at the rendered tube: no fold
where the yarn bends inside its own cross-section, no point of a tube inside another tube. It is
sampled, so it is evidence and not a proof. The letters (a) and (d) are the names everyone uses; there
is no gate (b) or (c).

**B1 and B2 ask: does settled yarn stay on its part?** Mesh charts only. B1 proves that the yarn stays
away from the mesh's open boundary. B2 samples the yarn and checks that each projection still belongs
to the part (`BoundaryCheck`, `layout/base_surface.h`).

**Certified, or a diagnostic.** The checks in `checks/` give certified bounds in double. The GPU runs
the same gate in float32, and its report says so (`gpuDiagnostic`). The CPU code is the written
specification of the GPU twins, and it runs in production: the spiral audit, crossing attribution and
the flat swatches behind `flatPatternReference`.

## The one port

`KnitNumerics` in `numerics.h` is everything knitting asks of the compute device. Settle's free set,
repair's start diameter and verdicts stay on this side, as host code that runs each step once (no
loop decides on a result read back); each relax call's iteration runs and is decided on the device. `knit_knitting` owns the interface
because it cannot see `knit_gpu`; `src/gpu` implements it.

It has four groups of methods: gate (a) (`validate`, `curvature`), settle and repair (`relax`, `store`),
surface queries (`meshSurface`, `pieceQuery`, `classifyLattice`) and the gates on finished yarn
(`checkRenderedSurface`, `checkBoundary`). `PartFailure` ends one part and not the run.

- `GpuKnitNumerics` (`gpu/port/knit_numerics.cpp`) is production. Each worker owns one object, and
  there is no CPU fallback. One device call bypasses the port: the spiral pieces' gate (d)
  (`scene/knit/spiral_part.cpp`, `checkFinalSurface`) calls `gpu::checkSurface` directly.
- `CpuNumerics` (`tests/cpu_numerics.*`) is the same port in double, for the tests. `ForwardingNumerics`,
  in the same file, is a base for tests that override one method.

`BaseSurface` is the only other virtual interface, because the GPU's subclasses override the batched
queries. Everything else is a concrete function.

## What must stay true

| Rule | Held by |
|---|---|
| The component depends on `model` only and no other component gains or loses an edge to it (`src/ARCHITECTURE.md`). | the compiler (forwarding headers), `architecture_module_dependencies`, `cmake/component_report.py` |
| The includes follow the layers. `numerics.h` includes only `src/model`. | `architecture_knitting_layers`; outside this folder only `src/model` is on the include path, so the compiler refuses anything else (`CMakeLists.txt` here) |
| Only `knit/pieces.cpp` touches files. | `architecture_knitting_no_file_io` |
| No templates and no generic lambdas. | `architecture_no_new_templates`, `architecture_no_new_generic_lambdas` |
| A saved knit stays valid: key text and report text do not change by accident. Retired fields print as constants. | `knit_params_text_pinned`, `knit_stage_key_pinned`, `io_pin_validation_snapshot`, the `static_assert`s in `tube/param_text.cpp` |
| The file tree above names every file of the folder. | `architecture_knitting_doc_complete` |
| A spiral's placement and hand, the coverage the gap patches read and the ply-map key stay bit for bit. | `spiral_placement_pinned`, `gpu_spiral_placement_pinned`, `coverage_regions_pinned`, `ply_map_key_pinned` |
| A chart's skew is the angle between its rows and columns off 90 degrees. | `chart_skew_measures_sheared_columns` |
| A local signature agrees with the flat one; INC and DEC show the signatures `stitch/expected_signature.h` describes. | `local_signature_matches_flat`, `flat_inc_signature`, `flat_dec_signature` |
| Settle keeps the relaxed state whole; a repair move gate (a) rejects leaves curves and binding as they were. | `settle_keeps_the_relaxed_state_whole`, `repair_rejected_attempt_keeps_projection_state` |
| The repair makes one move, kept only when it passes at the target and is certified; its start diameter comes from the report. | `repair_one_certified_move`, `repair_start_diameter_from_the_report` |
| A surface that cannot be judged stays UNRESOLVED, never PASS. | `reference_check_unresolved_when_compressed` |
| The port gives one knit on both sides. | `knit_mesh_port_cpu_gpu` |
| The solver's rows match their written reference and their derivatives. | `ncp_gpu_rows_match_cpu_cylinder`, `ncp_jacobian_fd_edge` and its siblings |
| The fit on the CPU and on the device agree. | `surface_compression_gpu_vs_cpu`, `fold_radius_cpu_gpu_agree` |
| Parallel runs equal serial runs. | `validate_parallel_identical`, `gate_d_parallel_identical`, `settle_repair_parallel_identical`, `tiled_parallel_identical` |
| A device failure on one part does not stop the run. | `tiled_classify_failure_ends_the_part` |
| A long loop can be cancelled. | `tiled_validation_cancel` |

## Where to change what

| To change | Touch |
|---|---|
| A stitch's shape | `stitchTemplate` in `stitch/stockinette.cpp`, then `expected_signature.h` if the crossings differ. |
| How uneven a hand is | `stitch/hand_tension.cpp`, a `HandKnit` field in `model/params.h`, one row of `settingFields`. Hash channels keep their meaning. |
| A stitch pattern | `StitchPattern` in `model/params.h`, `flatStitchFaces`, `buildStitchProgram`. `flatPatternReference` tells you whether it fits at a gauge. |
| Which stitches a chart keeps | `LatticeRunBuilder` in `layout/chart_lattice.cpp`, and `gpu/surface/lattice.cu`, which classifies the nodes. |
| A check | The file in `checks/`, its GPU twin in `gpu/validation/`, and a case comparing them. |
| A row of the solver | `relax/knit_problem.h`, `tests/ncp_reference.*`, `gpu/solver/ncp_*`. |
| The cross-section | `tube/compression_fit.h`, which holds the host and device twins, and `gpu/surface/surface_build.cu`. |
| A `PatternParams` or `SurfaceParams` field | `tube/param_text.cpp`, or the build fails. Keep retired fields as constants, or every saved knit goes stale. |
| A new device batch | `numerics.h`, `GpuKnitNumerics`, `CpuNumerics`. |

## Literature

- TopoKnit (arXiv:2101.04560): how a knit stitch is drawn through its parent, the source of
  `expected_signature.h`.
- Yuksel, Kaldor, James, Marschner 2012, *Stitch Meshes for Modeling Knitted Clothing with Yarn-Level
  Detail*: the overall approach of yarn from a stitch layout; its §6.1 is the 1D relaxation of positions
  along a course, `turnFractions`.
- Kaldor et al. (YarnCloth): yarns as B-spline tubes.
- Krawczyk test: a crossing exists and is unique, `crossings.cpp`, with the intervals of `interval.h`.
- Schur's comparison theorem: the same-yarn window of `clearance.cpp`.
- Wang, Jüttler, Zheng, Liu 2008, *Computation of Rotation Minimizing Frames*: `transportFrames`.
- Zhao, Luan, Bala 2016, *Fitting Procedural Yarn Models for Realistic Cloth Rendering* (§3.2,
  Algorithm 1), and Wu, Yuksel 2017, *Real-time Fiber-level Cloth Rendering* (§3.1): plies around the
  yarn centre, `tube/plies.h`.
- Kaldor et al. 2010, Leaf et al. 2018, Sperl et al. 2021: yarn-level methods do not require
  overlap-free round tubes, the reason for `compress_to_fit`.
