# Scene: file tree and architecture

`src/scene` takes a scene from its source through Knit and Fibres to what the export writes. This
file is the way into the folder: what each file owns, how a scene moves through it, and what must
stay true when you change it.

## File tree

The folders follow the pipeline. A `*` marks a header private to `src/scene`: only this folder's files
(and its tests) can include it. The compiler holds that: `CMakeLists.txt` lists the public and private
headers, and other folders see only the public ones (src/ARCHITECTURE.md).

```text
src/scene/
  ARCHITECTURE.md             this file
  CMakeLists.txt              the component knit_scene: its public and private headers; depends on model,
                              knitting, fibres, io and settings, privately on gpu, import and stitchmesh
  scene.h, .cpp               Scene, ScenePart, KnitLayout, KnitStop, questions
  run_log.h, .cpp             StageClock (wall, CPU and validate time per stage), process memory
  result_keys.h, .cpp         stage keys: what makes a result current (they hash the solver's identity)
  stage_params.h, .cpp        settings to knitting's, fibres' and the yarn body's parameters: hand,
                              final surface, filling depth, fibres, fur, yarn body

  source/                     from a scene entry to parts, nothing knitted
    scene_parts.h, .cpp       load a GLB or import and its solids, one part per component or chart
    guided_source.h, .cpp     the guided route: the source file's pieces, the marked one knitted, the others plain
    sewn_pieces.h, .cpp       the teddy: manifest, gauge, origin, its solids, one part per piece
    scene_solids.h*, .cpp     closed components and their inside tests (Scene::solids, made at load)

  appearance/                 how a part looks: derived, never stored in a knit
    rendered_part.h, .cpp     the strands as drawn: the yarn or its plies
    stitch_colours.h, .cpp    one colour per stitch, the scene palette (made by the first reader)
    ply_maps.h, .cpp          ply shading maps, ply_maps.json and their cache (the export attaches a set
                              to a package: attachPlySet)
    baked_maps.h, .cpp        the maps baked from a yarn body's fibres (the bake lives in the fibres component),
                              one atlas per yarn, in a content-addressed cache keyed by the knit's surface, the
                              body and the texel; a package derives its cavity and roughness maps from an entry
                              (the material values never bake). Modelled on ply_maps
    filling.h, .cpp           the filling behind a knitted part: its surface moved behind the yarn
    yarn_tube.h, .cpp         a yarn's tube mesh on the device: the one mesher for display, export, coverage;
                              and its surface model (the yarn body)

  knit/                       the Knit stage
    knit_stage.h, .cpp        knitScene: the stage list, time caps, cross-part clearance
    knit_run.h*               what every part of a run shares, and the two operations a part knits by
    mesh_part.cpp             one part on its UV charts, trying each diameter (declared in knit_run.h)
    guided_part.cpp           the guided part: stitch mesh, yarn, settle, stitch records (declared in knit_run.h)
    spiral_part.h*, .cpp      one spiral part (a sewn piece, an import's sphere), stage by stage
    spiral_cache.h*, .cpp     the per-piece knit cache files
    gap_patches.h*, .cpp      coverage of a component's parts at a fixed tube, patches for what is left open
    knit_report.h, .cpp       the text report of a knitted scene
    saved_knit.h, .cpp        a knit stored with its cache entry, restored and verified
    knit_cache.h, .cpp        .knitcache: one knitted piece with its attempts, and the curve encodings
                              knit.bin shares

  fibres/                     the Fibres stage
    fibre_stage.h, .cpp       growSceneFibres, wool in budget order: the stage list
    scene_fibres.h, .cpp      what the stage makes: PartFibres, SceneFibres, and what the export reads of them
    part_fibre_job.h*, .cpp   PartFibreJob: prepare, generate, finish one part
    fur_stage.h, .cpp         fur on the meshes, with its masks
    fibre_files.h, .cpp       saved fibres
    reference_samples.h, .cpp the Fibres step's preview samples: a flat wool patch, a fur sphere
```

**Reading order:** `scene.h`, `source/scene_parts.h`, `knit/knit_stage.cpp`, `knit/knit_run.h`,
`knit/mesh_part.cpp`, `fibres/scene_fibres.h`, `fibres/fibre_stage.cpp`, `appearance/rendered_part.h`. What the user chose is
the settings component, `src/settings/ARCHITECTURE.md`.

## The pipeline

A plain description comes first, stages compute from it, and each result flows one way.

```text
SceneEntry + SceneSettings
    │  loadPartsToKnit, or loadBaseParts when the Knit step is skipped
    ▼
Scene: source and parts, nothing knitted
    │  knitScene: every part in order, then gap patches, clearance
    ▼
KnitResult                  what the user reviews (jobs/results.h)
    │  growSceneFibres: every object's candidates in a budget order of its own, or fur on the meshes
    ▼
FibreResult
    │  exportScene (src/export)
    ▼
the package
```

`knitScene` and `growWoolFibres` are the two functions to read first. Each is a list of named stages
in the order they run, the way a renderer's scene update is. The source loads first, and its solids with
it (`loadMeshScene` or `loadPiecesScene`), then the parts. A part knits by one of two operations, declared
together in `knit_run.h`, and `knitPart` in `knit_stage.cpp` is the one place that chooses: `knitMeshPart`
on UV charts, `knitSpiralPart` for a spiral part, `knitGuidedPart` for the guided route's one piece. Which one is the part's own `ScenePart::layout`: `sewn_pieces.cpp`
sets it for the teddy's pieces; every part of an import is on charts (`collectParts`). A spiral part's mesh carries
its chart as uv, so nothing reads a manifest by part index. A spiral part's stitch colours come
from the source primitive's uv through the source triangle, never from that chart.

## Layers

What a folder may include (`architecture_scene_layers` refuses anything else; it does not require an edge
to be used, so removing an include never fails it):

| Folder | May include `scene/...` from |
|---|---|
| root (`scene.h`, `run_log`, `result_keys`, `stage_params`) | root |
| `source/` | root, `source/` |
| `appearance/` | root, `appearance/` |
| `knit/` | root, `source/`, `appearance/`, `knit/` |
| `fibres/` | root, `source/`, `appearance/`, `fibres/` |

The edges the code has today:

| Edge | Made by |
|---|---|
| source -> root | `scene.h`; `stage_params.h` (`scene_parts.cpp`, `sewn_pieces.cpp`) |
| appearance -> root | `scene.h`; `stage_params.h` (`rendered_part.cpp`) |
| knit -> root | `scene.h`, `stage_params.h`, `run_log.h` |
| knit -> source | `knit_run.h`, `spiral_cache.h`, `knit_stage.cpp`, `saved_knit.cpp` |
| knit -> appearance | only `gap_patches.cpp` -> `yarn_tube.h` (coverage measures a tube) |
| fibres -> root | `scene.h`, `stage_params.h`, `run_log.h` |
| fibres -> source | `fibre_stage.cpp`, `fur_stage.cpp`, `part_fibre_job.cpp` (`scene_solids.h`) |
| fibres -> appearance | `part_fibre_job.h`, `reference_samples.cpp` -> `rendered_part.h`; `fur_stage.cpp` -> `stitch_colours.h` |

- The settings are a component below this one (`src/settings`): they know no scene. The functions that
  turn settings into knitting's and fibres' types are here, in `stage_params`.
- `source/` and `appearance/` do not include each other; `knit/` and `fibres/` do not either. Fibres read
  a knit through `Scene`.

## Design

**Records and stages.** `Scene` and `ScenePart` are plain records, and stages are free functions over
them. A class appears only where an invariant needs guarding: `PartFibreJob` keeps a part's support
and its surface parameters in agreement, `SceneSolids` keeps meshes alive behind their queries.

**Results are objects.** The solver is not deterministic, so knitting the same settings again need
not give the reviewed result. The fibre and export jobs take `KnitResult` and `FibreResult`
(`jobs/results.h`), never settings, and read them only.

**Keys and computation are separate.** `result_keys.*` says what makes a result current. No
stage computes a key, and there are no dirty flags: a key either matches or it does not.

**A failed part is a result.** `ScenePart::stop` says why a part that should have a knit has none. A
run fails only for an error of the run, such as a device failure, which is never replaced by a CPU
stand-in.

**One part at a time.** The solves take the device one after another, and parts knitting side by
side only added their device memory to the running solve. So `knitParts` is a plain loop on the
caller's thread, with one set of numerics. Parallelism lives inside a part (`ExecPolicy`).

**What may change after it is made.** A `Scene` is filled by the load and the Knit stage and then read.
Two things change later, each with its lock: `verifyKnit` fills a restored part's `knit->after` once,
under `KnitResult::m` (jobs), and the stitch colours are made by their first reader. The solids are made
when the source loads (`loadMeshScene`, `loadPiecesScene`) from `source` alone, which never changes after,
so they need no lock.

**The stitch colours' contract** (`stitch_colours.h`). Inputs: every part's curves and stitch topology and
the source textures. Owner: the `Scene`, for its lifetime. Made by the first `sceneStitchColours`, which
always comes after `knitScene`, `restoreKnit` or `loadBaseParts` returned (no stage calls it), so the parts
are complete and nothing invalidates it. Lazy because colourwork is an Export setting: switching it on must
not knit again, and knit-only runs and fur do not need the colours. The function-local lock stays because
the process-wide last knit and last base scene (`jobs/job_steps.cpp`) let two builders in one process reach
one scene; it costs one uncontended lock per call.

**Progress, not dispatch.** `KnitRunOptions::onLoaded` and `onPart`, and `FibreRunOptions::onPart`, are the
only callbacks. They run on the caller's thread, in order: `onLoaded` each time the part list is set (after
loading, and again after gap patches add parts), `onPart` when a part starts and ends. The builder
publishes partial scenes from them while the stage runs on its thread; a polled record would need a
second thread and a lock over the parts, and a step-by-step interface would copy the stage order into
every caller.

**What a knit depends on.** Only its key. Gap coverage measures a fixed tube (0.6 mm rings, 6 sides,
`gap_patches.cpp`), not the tessellation (`strand.tube_*`, Fibres rows), so the tessellation never changes
which gap patches exist.

**Static storage, by kind.** Identity constants hashed into keys and formats are `constexpr`
(`kGeneratorVersion`, `kFurVersion`, `kRouting`, the files' `kMagic`, `kRetiredInputs`,
`kPlyMapGenerator`). Locks: the colour memo's and the pieces' GLB colours (`piecesObjectColours`). One cached
value lives for the process: those GLB colours, keyed by scene id (`sewn_pieces.cpp`).

**Deliberately absent.** No dependency graph: the pipeline has one fixed order, and nobody wires it.
No node classes, observers, plugin registry or entity system.

## GPU-ready data

A scene holds each concept in one flat form, and only `src/gpu` builds device views from it, behind its
four public headers (`gpu/device.h`, `knitting.h`, `surfaces.h`, `fibre_generator.h`).

| Data | Flat form | Read by |
|---|---|---|
| Yarn thickness | one `float` per control on `YarnPath`, derived from the stitches, so a saved knit keeps its bytes | the ring builders |
| Fibre islands | one `IslandWeight` per control of the strands, on `FibreFrame` | both fibre generators |
| Fibres | `PartFibres::fs`: `Fibre` records and a packed `Vec3f` pool, shared with the kernels | device, viewer, export |
| Stitch colours | one `uint32_t` per stitch, per part | viewer, export, fibre lines |
| Meshes | `shared_ptr<const TriMesh>` on a part; on the device `residentMesh`, once per content | every stage |
| Inside tests | one `SceneSolids` per scene (`Scene::solids`), made when the source loads | Knit, Fibres, fur |
| Clearance | one `gpu::ClearanceSegment` per sampled segment | `checkSceneClearance` |

Nothing keeps a second copy for convenience. The export makes strands from fibres with `partStrands`
when it writes them.

**When a loop moves to the device.** All three must hold:

1. It takes at least 10 % of a Draft run on a quiet device. Timings beside another GPU job mean
   nothing: the same man run took 60 s alone and 402 s shared.
2. Existing structures carry it (`DeviceGrid`, `residentMesh`, `DevicePrimitives`), in about 150
   lines of kernel code.
3. Its result equals the CPU code's in a test, or for colours a shown render accepts the difference.

The CPU version then moves to `tests/` as the written reference. Gap coverage stays on the CPU by
rule 3: it decides which patches exist, so a difference changes the knit, and fused multiply-add on
the device moves a double distance by a few ulps.

## What must stay true

| Rule | Held by |
|---|---|
| The settings table's order and stage text never change by accident: stored results are keyed by them. | `knit_stage_key_pinned`, `fibre_stage_text_pinned`, `io_pin_stage_keys` |
| A retired key keeps loading. If it was in a stage's text, its line stays there (`kRetiredInputs`). | `scene_settings_retired_keys`, `scene_settings_retired_knit_keys`, `scene_settings_retired_hair_keys` |
| Saved knits and fibres keep their bytes. | `io_pin_*`, `tiled_knit_cache_roundtrip`, `saved_fibre_key_check` |
| Parts knit in scene order, one at a time. | `knit_parts_in_scene_order` |
| Time caps and cancellation reach every part. | `deadline_marks_part_unresolved_and_rest_not_started`, `teddy_spiral_pieces_time_cap` |
| The device finds the CPU's close pairs. | `gpu_scene_clearance_gpu_vs_cpu` |
| A knit depends only on its key: the tessellation does not change the gap patches. | `gap_coverage_ignores_tube_settings` |
| Every load makes the scene's solids. | `scene_solids_made_at_load` |
| A second part list keeps the partial scene's textures. | `knit_stand_ins_survive_the_second_part_list` |
| The export's format writers see no scene type. | `architecture_export_door` |
| Includes follow the layers. | `architecture_scene_layers` |
| Private headers stay private; only the declared folders are included. | the compiler (`CMakeLists.txt`) |
| This file names every file, and no file that is gone. | `architecture_scene_doc_complete` |

## Where to change what

| To change | Touch |
|---|---|
| A setting | one row in `makeTable`, `src/settings/scene_settings.cpp`. Its place in the table is part of the stage keys. |
| Retire a setting | delete its row, add it to `kRetiredKeys` with its last default; if it was in a default stage text, also to `kRetiredInputs`. |
| What a knit depends on | a Knit row, or `stageInputs` in `result_keys.cpp`. |
| A stage of the Knit step | `knitScene` in `knit/knit_stage.cpp`. Keep `cancelCheckpoint` in its loops. |
| One part's knitting | `knit/mesh_part.cpp` or `knit/spiral_part.cpp`. |
| A fibre stage | `growWoolFibres` in `fibres/fibre_stage.cpp`; one part's steps are in `fibres/part_fibre_job.cpp`. |
| A new preset group | `src/settings/file_presets.*`: a `PresetGroup` value, its keys, its messages, its rules. |
| A new scene | `src/settings/scene_registry.cpp`, with its model's credit. |

## Literature

- Pharr, Jakob, Humphreys, *Physically Based Rendering* (4th ed.): a plain scene description, then one
  creation step.
- OpenUSD Hydra: scene data apart from its render backends.
- Gaffer: a node's `hash()` and `compute()` are separate.
- Ousterhout, *A Philosophy of Software Design*: deep modules, one entry per stage.
