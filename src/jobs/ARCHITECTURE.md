# Jobs: file tree and architecture

What the viewer and `knit_check` start: four jobs (Knit, the strand mesh rebuild, Fibres, Export), the queue that
runs them on a worker thread, the results the user reviews, and the display data the viewer draws. This file is the
way into the folder: what each file owns, how one job runs from request to published scene, and what must stay true
when you change it. The viewer's screens are in `src/viewer`; what the jobs compute is in `src/scene` (and below it
`src/knitting`, `src/fibres`, `src/gpu`); the files written by an export are in `src/export`.

## File tree

```text
src/jobs/
  ARCHITECTURE.md        this file
  CMakeLists.txt         the component knit_jobs (src/ARCHITECTURE.md): 9 public headers, 2 private;
                         public on model, settings, export, import, io; private on knitting, fibres, gpu, scene

  Level 0: what the viewer reads (public)
  results.h, .cpp        KnitResult and FibreResult, the results the user reviews; whether they are current
                         and whether they may be exported (exportEligibility)
  built_scene.h          BuiltScene: what the viewer draws after a job, as CPU arrays

  Level 1: the jobs' interface, the display builders, the cache
  jobs.h                 the four requests, JobControl, JobReporter, the four job functions, computeCaps (public)
  display_scene.h, .cpp  how a BuiltScene is built: source meshes, yarn tubes, fibre lines and their chunks,
                         bounds, the texture preview (private)
  scene_cache.h, .cpp    the persistent cache of Knit results, so a scene opens without knitting (public)

  Level 2: what the jobs share (private)
  job_steps.h, .cpp      which scene a job runs on, reporting, usable knits, ids, the last knit in memory,
                         the partial scenes of the Knit job, the progress of the Fibres job

  Level 3: one file per job
  knit_job.cpp           knitObject: the scene cache, the knit with partial scenes, the display, the cache write
  strand_mesh_job.cpp    rebuildStrandMesh: the shown knit's tubes at a new strand tessellation, or under their
                         yarn maps (StrandMeshRequest::bakeMaps bakes the maps the cache lacks first)
  fibre_job.cpp          growFibres: the fibres of the reviewed knit (or hair on the meshes) and their lines
  export_job.cpp         exportResults: the export rule on the job's snapshot, then the export's door

  Level 4: the queue (public)
  scene_builder.h, .cpp  SceneBuilder: one worker thread, latest wins, one published scene per finished job

  Beside the queue: the viewer's other workers (public)
  previews.h, .cpp       the Fibres step's reference preview and the flat-swatch references, each on its own
                         worker; the ply cross-section sketch; the hair coverage overlay
  fur_preview.h, .cpp    the Fur step's disc preview on its own worker: the loaded scene and each part's root caches
                         kept, the fields and display buffers per edit, drafts while a slider moves
  import_session.h, .cpp the Import step's worker: read and analyse a file, then prepare and store it
  marking_session.h, .cpp the Mark step's two lanes. Deep: read the imported file, its pieces and their display; Method A
                         on a marked piece (stitchmesh, private) as display data and a review; proposals for a piece and rule
                         2's refinement; cancellable (a newer request cancels the running one). Quick: the quick review while
                         drawing on its own worker, against the kept row field; published only for the newest request.
                         Every finished build is published per piece and knitting component (Snapshot::builds), also a failed
                         one: its polygons and what the checks measured as model data (StitchMeshData, StitchMeshQuality), or why
                         there is no mesh; the Stitch Mesh view's display meshes and legend colours
  render_preview.h, .cpp the live render preview's worker: one renderer process kept warm, its package
                         export, and its residency under the GPU jobs' priority lease
  crash_report.h, .cpp   what the program could not catch, written down before the process ends
```

**Reading order:**

1. `jobs.h`: the four jobs, what each takes and what it returns.
2. `results.h`: what a job hands on to the next one, and the export rule.
3. `knit_job.cpp` from `knitObject` upwards: one job from its request to its result scene.
4. `scene_builder.h`: how the viewer runs the jobs and reads their results.
5. `built_scene.h`, then `display_scene.h`: what the viewer draws, and how it is made.

## The common recipe

The jobs follow the way interactive programs run long work in the background (threaded resource
loading): an **active object** owns one worker thread; callers post a request and read an immutable result.

- **Requests in, immutable results out.** A request holds a settings snapshot and the reviewed results it consumes.
  A finished job publishes one `std::shared_ptr<const BuiltScene>`; nothing changes it afterwards.
- **Latest wins.** The newest pending request replaces an older one, and a newer request cancels a running Knit,
  strand mesh or Fibres job. An Export job finishes from its snapshot unless `cancel()` stops it.
- **Cooperative cancellation.** The worker installs a `CancelToken` for the thread (`model/cancel.h`); the long
  loops in scene, knitting and fibres call `cancelCheckpoint()`, which throws `Cancelled`. A cancelled or failed job
  publishes nothing, so the previous result stays.
- **One coherent read.** The viewer polls `SceneBuilder::snapshot()` once per frame; the fields of a `JobSnapshot`
  are taken under one lock and agree with each other.
- **Use cases.** Each job is one file with one entry function: a request in, a result scene out (Clean
  Architecture's interactor, without its interfaces). The queue is a scheduler and knows nothing of what a job does.

What this program adds on top:

- **Reviewed results flow downstream.** Fibres take a `KnitResult`, the export a `KnitResult` and a `FibreResult`;
  nothing knits again from settings keys, because the NCP solver is not deterministic.
- **Partial scenes.** A long Knit job shows its stand-ins first and each part as it finishes. A partial scene never
  carries a result.
- **The scene cache.** A Knit job whose result is stored (`scene_cache.h`, keyed by the settings and assets that
  determine it) is read, not computed; the stored knit (`knit.bin`) makes it a computational knit again.

## Layers

```text
level 4  scene_builder ................................ jobs
level 3  knit_job  strand_mesh_job  fibre_job  export_job ... job_steps
         render_preview ............................... job_steps (the preview's package), jobs
level 2  job_steps .................................... jobs, display_scene, scene_cache
         previews, fur_preview, import_session ........ built_scene, display_scene
         marking_session .............................. built_scene
level 1  jobs ......................................... built_scene, results
         display_scene, scene_cache ................... built_scene
level 0  built_scene .................................. results
         results, crash_report ........................ (no jobs header)
```

A file includes jobs headers of lower levels only (and its own header). The jobs never see the queue: a job is a
function that tests and `knit_check` call directly, and the queue only calls it.

## One job, end to end: the Knit job

1. The viewer fills a `KnitRequest` (the settings snapshot, the cache root) and calls `SceneBuilder::startKnit`.
2. `queue` numbers it (its revision), replaces a pending request and cancels a running job that is not an Export.
3. The worker (`loop`) takes it, installs its `CancelToken`, releases the GPU meshes of another scene
   (`gpu::releaseResident`), and calls `runJob`, which calls `knitObject` (`knit_job.cpp`).
4. `knitObject`: `jobScene` finds the scene entry (or refuses); with a cache root, `restoreCachedKnit` reads a
   stored result and its `knit.bin`. Otherwise:
   - `newResultScene` stamps the result with the job's kind, revision and settings;
   - `rememberedOrNewKnit` takes the last knit in memory when its Knit-stage key matches, or runs `knitScene` with
     `PartialKnitScenes` (stand-ins first, then each part's coarse tube);
   - `knitDisplay` makes the rows and the tubes; `adoptKnit` and `summarize` fill the `KnitResult`;
   - `storeInSceneCache` writes the entry.
5. The worker assigns the result under its lock (`latest_`, `complete_`); the viewer sees it in its next snapshot.

The strand mesh rebuild copies the shown scene and replaces its yarn tubes; the Fibres job grows fibres on the
reviewed knit (`growSceneFibres`) and draws their lines; the Export job checks `exportEligibility` on its own
snapshot and calls the export's door (`exportScene`).

## The one interface and the one global

- **`JobReporter`** (`jobs.h`) is where a running job reports its progress and partial scenes. The queue's own
  reporter puts them into its snapshot; tests record them. It is the folder's one virtual interface.
- **The last knit and the last base scene** (`job_steps.cpp`) are kept per process, so a job that changes only later
  settings (plies, materials, fibres) shares the knit instead of knitting again. It is the folder's one piece of
  hidden state; a knit with a refused or stopped part is never kept.

## What must stay true

| Rule | Held by |
|---|---|
| A cancelled or failed job publishes nothing; the previous result stays | `s6_builder_*` |
| A partial scene never carries a result | `s6_builder_progressive_publish`, `workflow_partial_base_inside_fabric` |
| Downstream jobs consume the reviewed results; a changed knitting setting refuses them | `usableKnit`, `workflow_unified_workflow_rules` |
| One export rule for the UI and the job | `exportEligibility`, `workflow_export_eligibility_rows`, `workflow_export_job_enforces_eligibility` |
| A cache hit is a verified publication of its key | `scene_cache.cpp`, `s6_scene_cache_roundtrip`, `s6_scene_cache_key` |
| A knit with a refused or stopped part is neither remembered nor cached | `refusedParts` |
| A null texture preview keeps the shown texture; an empty one (the map set failed) replaces it | `texturePreview` (`display_scene.h`) |
| The jobs include only lower layers; a public header includes no private one | `architecture_jobs_layers` |
| This file names every file of the folder | `architecture_jobs_doc_complete` |

## The scene cache

`scene_cache.h` stores a Knit job's `BuiltScene` (`scene.bin`) and its computational knit (`knit.bin`) under a key
of the scene, its Knit and Surface stage keys, the asset hash, the strand tessellation and two versions:

- **`kFormat`** is the layout of `scene.bin`. Change it with any field of `BuiltScene`, `DisplayMesh`,
  `DisplayLines`, `DisplayTexture` or `PartRow` that is stored; old entries then miss once.
- **`kGenerator`** is bumped when the Knit job draws something else for the same settings and assets.

The file is written with `io/files.h`'s `ByteWriter` (little-endian, counts as i64) and published through
`io/publish.h`, so a hit is a verified publication and a corrupt or partial entry is a miss. Material and display
settings are not in the key: the viewer applies them when drawing.

## Where to change what

| To | Change |
|---|---|
| Add a job | its request and function in `jobs.h`, its file at level 3, a branch in `runJob`, `jobKind` and `jobSettings` |
| Change what the viewer draws | `built_scene.h` (and the scene cache format, `kFormat`), `display_scene.cpp` |
| Change the Knit job's display for the same settings | `display_scene.cpp`, and bump `kGenerator` in `scene_cache.cpp` |
| Change the export rule | `exportEligibility` in `results.cpp` |
| Change what a job reports | the stage texts in the job's file; the viewer classifies them (`viewer/presentation.cpp`, `progressView`) |

## Known limits

- The five workers (`SceneBuilder`, `FibrePreviewWorker`, `FlatReferences`, `ImportSession`,
  `RenderPreviewWorker`) are written separately: their rules differ (partial scenes; a 150 ms pause in the edits;
  a result per key; request numbers retired at once; a renderer process and its residency), and a shared worker
  would need callbacks or templates.
- The knit's part threads report outside `PartialKnitScenes`'s lock, so two partial scenes can reach the reporter
  out of order; the next one, or the result, replaces them.
- Two "most common yarn diameter" rules exist on purpose (`mostCommonD` over the rows, `mostCommonKnitD` over the
  knit's parts): they count different parts.

## Literature

- R. G. Lavender, D. C. Schmidt, "Active Object" (Pattern Languages of Program Design 2, 1996); H. Sutter, "Prefer
  Using Active Objects Instead of Naked Threads" (Dr. Dobb's, 2010).
- S. Parent, "Better Code: Concurrency" (2016): immutable values shared by pointer instead of shared mutable state.
- R. C. Martin, *Clean Architecture* (2017), chapters 20 and 22: use cases with request and response models;
  chapters 13 and 14: the component principles.
