# Viewer: file tree and architecture

`knit_viewer`: the window, the four workflow steps and the 3D view. This file is the way into the
folder: what each file owns, how a frame runs, and what must stay true when you change them. The
viewer only displays and edits: all geometry, checks and exports come from the libraries, through
the jobs in `src/jobs`; it edits settings with `src/settings`.

## File tree

```text
src/viewer/
  ARCHITECTURE.md        this file
  CMakeLists.txt         the component (the state library knit_viewer_state) and the program knit_viewer;
                         every header private; depends on model, io, settings, export and jobs

  The program            GLFW, the ImGui backends and the native dialogs live here
  main.cpp               process, window, fonts on the GPU, frame(), the window's answers, the shot, the bench
  native_dialogs.h, .cpp file and folder dialogs, the file manager

  The state              no ImGui, no GL, no GLFW
  workflow.h             SceneUi and its groups, every action, the steps' words
  workflow.cpp           what CHANGES SceneUi: open, choose, start, cancel, poll
  presentation.cpp       what READS it: stage words, progress, export texts, the Import step's open file
  marking_flow.h, .cpp   the Mark step: MarkingFlow (pieces, the set, undo, the overlay, suggestions and their preview,
                         quick-review requests and latency, the layout baseline and its automatic refresh, every piece's
                         last build as it arrived, the Stitch Mesh view's switches) and what changes it, Export Stitch Mesh
  view_content.h, .cpp   what the viewport shows: its content, the stages and Display rows, the comparison,
                         the scenes on the GPU, the counts, picking
  viewer_script.h, .cpp  the command line, its settings and object, the --do steps, when a shot is due
  viewer_prefs.h, .cpp   viewer.ini and the ImGui layout file
  camera.h, .cpp         Camera, CameraDrag, CameraRig, the presets, the ray through a pixel
  about.h, .cpp          what the About window says: author, licences, third-party notices, credits, Blender, privacy
  about_texts.h          the licence texts and the build's commit, compiled in (cmake/embed_texts.cmake,
                         cmake/build_identity.cmake write their definitions into the build folder)

  The interface          ImGui, no GL, no GLFW
  interface.h            drawSceneUi, the look, the widgets, the step panels
  interface.cpp          strip, view toolbar and its popovers, viewport labels and status line, status bar,
                         Report and Log
  widgets.cpp            colours, scale, font atlas, glyphs, rows, lists
  panel_import.cpp       the Import step: the object so far, or the open file in five sections and its footer
  panel_mark.cpp         the Mark step: pieces, suggestions (cards), the quick check while drawing, the layout's findings and
                         measures with their availability, markings, preview settings, its footer
  panel_knit.cpp         the Knit step (construction tiles, fit checks, character tiles, detail nodes), its footer, the parts table
  panel_fibres.cpp       the Fibres step and its footer (the Fur step's footer too)
  panel_fur.cpp          the Fur step: coverage, loops and their variation, fine fibres, appearance (workflow = fur)
  panel_yarn_fibres.cpp  the Fibres step's Yarn fibres section (a yarn made of fibres: its look, structure, shading, maps)
  panel_export.cpp       the Export step and its footer
  about_window.cpp       the About window (status bar, About): five tabs, links, licence texts to copy

  The 3D view            GL, no ImGui
  gl_functions.h         the OpenGL entry points above version 1.1
  renderer.h, .cpp       one frame of the view, whole or split for the comparison; everything on the GPU
  mesh_lighting.h, .cpp  Flat and Studio lighting, the shaders, the SSAO kernel
  fibre_lines.h, .cpp    fibres as lines: chunks, culling, thinning
  fur_discs.h, .cpp      the Fur step's disc preview: one triangle cut to a circle per loop, marks and needles as lines
  overlay_lines.h, .cpp  the Mark step's markings, rings and stitch outlines over the meshes (client arrays)
```

**Reading order:** `workflow.h`, `presentation.cpp`, `view_content.h`, `frame()` in `main.cpp`, `interface.h`, one
panel, `camera.h`, `renderer.h`, `mesh_lighting.h`.

## The common recipe

Viewers built on an immediate-mode UI share one shape:

1. **One loop.** Poll events, start the UI frame, build the UI, draw the 3D view, draw the UI, swap.
2. **The UI is code.** It runs every frame, top to bottom. No widget objects to keep in sync.
3. **One-way flow.** Input calls an action, the action changes the state, the next frame draws the
   state.
4. **A humble view.** Logic lives where a test can reach it; the drawing code only says what the
   logic says.
5. **Input as data.** The camera takes one frame of mouse and needs no window system; what the
   interface asks of the window is a request, answered after the frame.

That is enough for a model viewer. Ours also runs knits that take minutes on a background thread,
so most of `workflow.cpp` is about a result arriving while the interface goes on, and most of
`presentation.cpp` is about saying honestly how that result relates to the settings on screen.

## Layers

An arrow means "includes". There are no cycles, so each file can be read knowing only the files to
its right.

```text
main.cpp ──┬──> native_dialogs.h
           ├──> renderer.h ──┬──> mesh_lighting.h ─────────────┐
           │                 ├──> fibre_lines.h ──> camera.h,  ├──> gl_functions.h
           │                 │                      jobs (display types)
           │                 └──> view_content.h ──┐
           ├──> interface.h ───────────────────────┤
           ├──> viewer_script.h ───────────────────┤
           │                                       v
           └──────────────────────────────> workflow.h ──> camera.h, viewer_prefs.h, jobs, settings, export
```

| Layer | Files | May use | Held by |
|---|---|---|---|
| The state | `workflow.*`, `marking_flow.*`, `presentation.cpp`, `view_content.*`, `viewer_script.*`, `viewer_prefs.*`, `camera.*`, `about.*`, `about_texts.h` | model, io, settings, export, jobs | the target `knit_viewer_state` has no ImGui or GLFW on its include path; `architecture_viewer_layers` (GL, ImGui, components) |
| The interface | `interface.*`, `widgets.cpp`, `panel_*.cpp`, `about_window.cpp` | ImGui, the state; model, settings, export, jobs | `architecture_viewer_layers` |
| The 3D view | `renderer.*`, `mesh_lighting.*`, `fibre_lines.*`, `fur_discs.*`, `overlay_lines.*`, `gl_functions.h` | GL, GLFW's proc address, the state; model, io, jobs | `architecture_viewer_layers` |
| The program | `main.cpp`, `native_dialogs.*` | everything of the viewer; model, jobs | `architecture_viewer_layers` |

"Uses a component" means a direct `#include "<component>/..."`. The import's vocabulary
(`ImportOptions`, `ImportAnalysis`) reaches the viewer through `jobs/import_session.h`, whose header
exposes it: the viewer calls nothing of the import's store or pipeline itself. Keeping the state free
of ImGui is what makes it testable: `knit_tests_ui` drives all of it without a window. In
`presentation.cpp` and `view_content.cpp` no function takes a mutable `SceneUi&`, so they can only
read. The panels do not call each other: what they share is in `workflow.h` or in the widgets.

## The data

`SceneUi` is plain data. The groups say who writes what:

| Group | Holds | Written by |
|---|---|---|
| `s` | the settings being edited | panels (bound fields and the settings' edit functions), `switchScene`, `openExample`, `chooseWorkflow`, `runExportFix`, `knitAtFinalQuality`, `applyExportDetail`, `takeFibreResult`, `importPrepared`, `loadMaskImage`, `reloadMaskImage` |
| `builder`, `cacheRoot` | the job queue and the scene cache's folder | `runViewer` in `main.cpp` |
| `stage` | the step on show | the strip, actions |
| `results` | each stage's last completed result and attempt; the other workflow's, parked | `pollSceneUi`, `resetResults`, `chooseWorkflow`, `takeFibreResult` |
| `job` | the one running job and its progress snapshot | `jobStarting`, `cancelJob`, `pollSceneUi` |
| `import` | the import session and the open file's choices | `openMeshFile`, `changeImport`, `prepareKnittingCharts`, `presetImportUnits`, `confirmImport`, `closeImport`, the choice actions (`chooseImportUnit` ... `chooseImportUp`), node selection (`toggleNode`, `selectAllNodes`, `selectNoNodes`), focus (`focusNode`, `pickNode`), `startImport`, `adoptSettings`, `openExample`, `pollSceneUi`; the Import panel only for plain option flags (the engine, the repairs, `revealFocus` once its row is shown) |
| `exportTo` | destination and format | Export panel (format, name), `chooseExportFolder`, `resetDestination`, `switchScene` |
| `view` | what the viewport shows and how; none of it is a setting | toolbar, `setUpWindow`, actions; `frame()` the camera's preset (`cameraPreset`) |
| `coverage` | the fur coverage overlay | `pollSceneUi`, `resetResults`, `chooseWorkflow` |
| `workflowChosen` | whether Knit or Fur is chosen (the draft's unchosen state; with the Fur switch off, `knit::furEnabled`, always Knit) | `chooseWorkflow`, `openExample`, `openMeshFile`, `importPrepared`, `loadSettings` (scripts) |
| `previews` | what follows the settings without a job | `pollSceneUi` |
| `panels` | what the panels keep that is no setting (the pop-out windows open, About among them), the mask image check, the ply sketches | panels, the status bar, `applyExportDetail`, `checkMaskImage` |
| `layout` | sizes of the last frame | `drawSceneUi` |
| `camera` | requests to the camera | toolbar, `resetResults`, `pollSceneUi`; `CameraRig::apply` clears them |
| `platform` | what the window provides, and the one request to it | `setUpWindow` (the flags); panels (the request); `answerPlatformRequest` clears it |
| `prefs` | viewer.ini | `setUpWindow`; `openMeshFile`, the status bar's scale |

A `BuiltScene` is immutable once published (`src/jobs/scene_builder.h`), so `results` holds
`shared_ptr`s and nothing copies a mesh. A stage's result and its attempt are separate: a cancelled
knit records the attempt and keeps the last result on screen. `results.partial` is the running knit's
scene so far. It is a preview, never a result, and offers no Fibres or Export.

**Threads.** Only the main thread touches `SceneUi` and GL. Four workers publish under their own
locks, and the viewer reads them through `snapshot()`, `latest()` and `get()`: the builder's job
thread (one immutable `BuiltScene` per job, from a settings snapshot), the import session's worker,
`FibrePreviewWorker` and `FlatReferences`.

## One frame

`frame()` in `main.cpp` is the whole story:

```text
glfwPollEvents            a dropped file -> openMeshFile
pollSceneUi               job's end, completed result, progress (one snapshot of the builder), the
                          import's phase, the import's charts check, previews
ScriptedRun::beforeFrame  confirms a scripted import, runs the next --do step
beginImGuiFrame           scale, fonts on the GPU, the backends' NewFrame
CameraDrag::update        the mouse, with the layout of the last frame; a click -> pickNode
updateFrameStats          F3
drawSceneUi               the interface: a click calls one action, or sets a PlatformRequest
answerPlatformRequest     the dialog the interface asked for, then its action (openMeshFile,
                          loadMaskImage, chooseExportFolder) or the file manager
viewContent               what the viewport shows, now that the interface and the window have acted
CameraRig                 this view's camera, then the requests; ScriptedRun::frameShot
Renderer::draw            meshes under the chosen lighting, then fibres
frameStatsLabel           F3
ImGui::Render             the interface over the 3D view
Renderer::releaseUnused
shot, or swap             the shot reads the back buffer, so it comes before the swap
```

### One job at a time

`startKnit`, `startFibres`, `startExport` and `startStrandMesh` take a snapshot of the settings
(`jobStarting`) and hand a request to the builder. While a job runs, the panels say why they wait and
nothing restarts by itself. A cancelled or failed job publishes nothing, so the previous result stays;
its attempt carries the message.

Each step has one footer below its scrolling settings: the state, then one action, always with its
room. The status bar shows progress and Cancel only when the running job's own step is not on show.

### Honest states

The strip's words come from `importState`, `knitState`, `fibreState` and `exportState`. Each state has
a glyph and a word, so colour is never the only cue:

| Status | Means |
|---|---|
| `NotStarted`, `Optional`, `Blocked` | nothing yet; skipped or off; waiting for an earlier step |
| `Running`, `Cancelling`, `Cancelled`, `Failed` | the job |
| `Ready`, `Diagnostic` | current; or exportable with failed parts, and marked so |
| `NeedsUpdate`, `ViewOnly` | the settings moved on; or the result cannot be computed from |

A result is current when its key equals the key of the current settings (`knit::knitOutdated`,
`knit::fibresOutdated` in `jobs/results.h`). The words are tested as text in `ui_text_snapshot`.

### Choices and rules

A control bound to a plain field writes it directly, as `dragRow(&s.depthAmplitude)` does for a
setting. A choice that carries a rule calls a named action, so the rule is in one place and the tests
reach it without ImGui: the import's unit and height (`chooseImportUnit`, `setCustomImportUnit`,
`chooseSizeByHeight`, `setImportHeight`, `chooseImportUp`; a size that cannot be used is a unit of 0, and
`sizeProblem` says why), the workflow (`chooseWorkflow`, workflow.cpp: Knit or Fur; under Fur `nodeKeptAsMesh`
and `importRequest` give every component to the fur while the ticks wait in the draft, results are parked per
workflow and source; a scripted import follows the same rule), Knit draft and Knit final (`knitAtDraftQuality`, `knitAtFinalQuality`;
their notes and tooltips from the settings with that profile applied, `knitActionView`), Review Import
(`reviewImport`: the same file prepared again; `ImportFlow::keepKnitChoices` makes `importPrepared` keep
every Knit setting instead of switching scenes), and the shown tubes
following a new tessellation (`shownTubesDiffer`, then `startStrandMesh`). The Gauge list is one
choice over three keys (`knit::gaugeChoices`, `currentGaugeChoice`, `applyGaugeChoice`: a preset, a
stitch sized for the yarn at a tightness, or Custom). In yarn mode the stitch
size follows its inputs where they change (`knit::applyStitchSize` after the gauge rows, the clearance
and a scene switch), never once per frame. The stitch size is the stockinette stitch; the gauge line
under it gives a rib's or garter's own (`knit::patternStitch`).

### Requests to the window

What the interface asks of the window is data, as `CameraRequests` is for the camera: a panel sets
`ui.platform.request` (`OpenMesh`, `OpenMaskImage`, `ChooseExportFolder`, `ShowExportFolder`), and
`answerPlatformRequest` in `main.cpp` runs the native dialog after the interface and gives its answer
to the named action. Without dialogs (`ui.platform.dialogs`) the controls say how instead. A test sets
the flags and answers the request itself.

## The Import step

Stage 1 is the object so far and the ways to open another; stage 2 is the open file: read, then its
sections (what to make, components, real-world size, surface preparation, colour and compatibility),
then prepared (`importStartView`, `importFileView`, `importFooter`). The footer's words and whether its
action is enabled come from `importFooterView`, the sizes from `importSize`. A new file resets the draft
(`openMeshFile`); a failed or cancelled preparation keeps it. The
import session (`jobs/import_session.h`) reads, prepares and stores on its own thread; `pollSceneUi`
integrates each new phase once (`importAnalysed`: suggestions, or for an import named in the settings
`ImportSession::storedImportMatches`; `importPrepared`: the new scene). A knitted import is always
rebuilt (`applySurface`); one that kept the file's mesh because its knitting was skipped is rebuilt
before it knits (`refreshImportCharts`, the Knit step's "Rebuild the surface").

**Picking a component.** A click is a press and release within 4 pixels (`CameraDrag::takeClick`).
`pickNode` casts the ray through that pixel (`viewRay`), takes it back into the file's coordinates
through the inverse of the preview's transform (`rayInFile`, next to `addImportPreview` in
`view_content.cpp`), and tests every triangle (Moller-Trumbore): up to 2 million, once per click, so
nothing is indexed or kept. The nearest hit's component is focused, never ticked: the renderer outlines it
(`previewNodeIsFocused`), and `revealFocus` scrolls the form to the list, then the list to its row (only
once the list is on screen: ImGui lays out nothing in a hidden child). The background clears the focus.

## The Export step's forecast

Each detail level (`kExportDetails`: Full, Reduced, Light) sets the yarn tessellation and the fibre
segment budget. `exportDetailForecast` says what an export at that level would hold: the tube
triangles from the knit's tubed yarn length, the fibre segments from the result's count per object
(scaled to another budget), and the renderer's GPU memory (`knit::rendererMemoryGiB`, a fit in
`export/export.h`). `exportDetailMemoryWarning` warns above `knit::kRendererMemoryCapGiB`. It lives
here, next to the levels: they are the Export step's own choice, and the forecast changes with them.

## The 3D view

`Renderer::draw` clears, loads the camera and turns each shown mesh into a `MeshLook`: a colour, a
texture and whether its normals are trusted. Both lighting paths take the same looks.

**Studio** is four fixed passes in `StudioLighting::draw`:

```text
drawNormalDepth    view-space normal and distance of every pixel, off screen
computeOcclusion   SSAO at two radii: the yarn's grooves, and contact between parts
blurOcclusion      4 x 4, depth-aware
drawLit            per-pixel light: a key from the upper left, a fill, a sky / ground hemisphere
```

**Flat** is the fixed-function headlight. It draws from the same buffers, also serves wireframe, and
is the fallback where Studio cannot run. The window asks for GL 2.1, which has buffer objects, so
the loader goes on past a missing framebuffer function.

**Compare with source** draws two passes into one viewport (`ViewportFrame` from `viewportFrame`): the
source on the left, the stage on the right. Both take the whole viewport's camera and projection; only
the scissor (`ViewColumns`) differs, so the halves are two halves of one view and an orbit crosses the
divider without a seam. The framebuffer is cleared once (the halves never overlap); Studio's off-screen
passes cover the whole viewport each time, so the occlusion at the divider sees its real neighbours, and
only its lit pass is clipped. The source's buffers upload on first use, as for the Source view.

**Fibres** are lines, not tubes, and there are tens of millions of segments. `FibreLines` keeps them
in buffers and draws them by `knit::FibreChunk`: a chunk outside the frustum is skipped, and while the
camera moves every visible chunk draws the same share of its fibres. Their order inside a chunk is
random, so a prefix is an even thinning.

### One lifetime rule

Whatever is on the GPU belongs to a scene in `MeshOwners`, the scenes a view can show, and `Renderer`
keeps that scene alive while it holds its buffers. So a mesh's address cannot be reused under a stale
buffer. Scenes stop being listed when they are replaced; `releaseUnused` then frees their meshes, and
any texture no kept scene lists. A texture can be shared: the coverage overlay uses the hair's. A
parked result of the other fibre profile is not listed: it stays in memory, not on the GPU.

## Scripted runs

`viewer_script.*` is a run with nobody at the keyboard, and it has no GL, so `ui_script_steps` can test
it. `loadSettings` and `openObject` take the command line's scene, settings and object. `--do` runs its
steps when nothing else is running, `--shot` writes the window once everything has finished, `--bench`
times frames first. A scripted import comes in two steps on purpose: `presetImportUnits` before the
first poll, which skips the unit suggestion, and `confirmImport` once the file is analysed, with the
components `--plain-nodes` leaves out.

## What must stay true

| Rule | Held by |
|---|---|
| The state compiles without ImGui or GLFW. | the target `knit_viewer_state` |
| The viewer reaches the scene, knitting, fibres, the import's pipeline and the GPU only through the jobs: none of them is a dependency. | the compiler (`CMakeLists.txt`), `architecture_module_dependencies` |
| The interface has no GL, the 3D view no ImGui, `presentation.cpp` and `view_content.cpp` only read, and each layer includes only its components. | `architecture_viewer_layers` |
| This file names every file in the folder, and only those. | `architecture_viewer_doc_complete` |
| Every word on screen is deliberate. | `ui_text_snapshot`, against `tests/data/ui_text_snapshot.txt` |
| No two widgets share an id, and ids do not change: the tests and `imgui.ini` read them. | the `ui_*_inventory` and `ui_panel_unique_ids_*` cases |
| Every action can be reached at the window sizes we support. | the `ui_*_footer_reachable` cases |
| A view shows only what a kept scene holds, both halves of a comparison included. | `ui_shown_meshes_cover_view`, `ui_compare_with_source` |
| A stage with nothing to draw is disabled with its reason; an open file disables all three and keeps the chosen stage; fibres over a partial build are never shown. | `ui_view_mode_availability`, `ui_view_mode_open_file`, `ui_view_mode_partial_rebuild` |
| A Display row that does not apply is hidden and keeps its value; Unknitted parts is judged on the knit's scene, not on the view. | `ui_display_rows`, `ui_view_popovers` |
| Both halves of a comparison share one camera and the whole viewport's projection; only the scissor differs. | `Renderer::draw`; the pixel check in the toolbar handoff's report |
| The toolbar stays inside its window and each popover inside the viewport, at 1.5 x on a narrow window. | `ui_toolbar_narrow` |
| A press over the interface is no camera drag; a Shift-drag locks one axis. | `ui_camera_drag_lock` |
| A click is a press and release within 4 pixels: an orbit never focuses a component. | `ui_import_drag_is_no_click` |
| A click focuses the nearest component under the cursor, in either up axis, and reveals its row; ticks are the checkboxes' only. | `ui_import_pick_focuses`, `ui_import_pick_reveals_row`, `ui_import_focus_and_tick_are_separate`, `ui_view_ray_matches_the_camera` |
| Fur gives every component to the fur, the ticks wait in the draft, and switching back selects the rebuilt surface. | `ui_import_workflow_switch_selects_source` |
| The Fur step opens on the disc preview, which follows every edit. | `ui_fur_disc_preview` |
| A size that cannot be used disables the import with a reason; Set height and back keeps the unit. | `ui_import_invalid_size`, `ui_import_custom_survives_height` |
| The import's rules hold where a choice is made, not per frame; a scripted import follows them. | `ui_import_unit_choices`, `ui_import_skip_rules`, `ui_import_rules_without_a_frame`, `ui_import_scripted_skip_keeps_charts` |
| The window's dialogs are requests, answered after the interface. | `ui_platform_requests` |
| The mask image check reads the file's size and write time where the warning is drawn, and hashes on a change only. | `ui_fibre_mask_image_changed` |
| The command line and the `--do` steps behave as specified. | `ui_viewer_options`, `ui_script_steps` |
| The draw path does no scene-scale work: facts are computed when a result arrives. | the "Viewer draw paths" section of `STYLE.md` at the repository root |
| A result is never read from a half-edited settings object: a job takes a snapshot. | `jobStarting` |

## Where to change what

| To change | Touch |
|---|---|
| A setting's control | its panel: one `dragRow`. A new setting also needs `settings/scene_settings.*`. |
| What a state says | `presentation.cpp`, then re-pin `ui_text_snapshot` on purpose. |
| What a button does | an action in `workflow.cpp`; the panel only calls it. |
| An Import step choice with a rule | an action next to `chooseImportUnit`; the panel calls it. |
| A dialog | a `PlatformRequest` value, its button, its case in `answerPlatformRequest`, the action for its answer. |
| A new job | `jobStarting`, `integrateCompletedResult`, the step's footer. |
| What the viewport shows | one function under `viewContent` (`view_content.cpp`). A new scene in it also goes in `meshOwners`. |
| The look | `applySceneTheme` and the colours in `interface.h`. |
| A new widget | `widgets.cpp`, declared in `interface.h`, with `testLabel` if ImGui reports no label. |
| The lighting | `mesh_lighting.cpp`: the shaders are the strings at its top. |
| A GL function | one line in `gl_functions.h`, one in `loadGlFunctions`. |
| A command line flag | `parseViewerOptions`, `ui_viewer_options`. |
| A `--do` step | `applyScriptStep`, `ui_script_steps`. |

## Literature

- Muratori 2005, *Immediate-Mode Graphical User Interfaces*, and the Dear ImGui wiki, *About the IMGUI
  paradigm*: the panels.
- Dear ImGui's `example_glfw_opengl2/main.cpp`: the frame's order.
- Czaplicki, *The Elm Architecture*, and the Elm guide's *Commands and Subscriptions*: the one-way flow
  from actions to state to view; requests the runtime performs (`PlatformRequest`).
- Feathers 2002, *The Humble Dialog Box*, and Fowler, *Humble Object*: why `presentation.cpp` and the
  actions exist.
- Filament's viewer library: a GUI apart from the viewer, and an automation engine for batch runs.
- Lakos 1996, *Large-Scale C++ Software Design*: levelization, and what a header costs its includers.
- Martin 2002, *Agile Software Development*: acyclic dependencies, common closure, common reuse.
- Parnas 1972, *On the Criteria To Be Used in Decomposing Systems into Modules*; Ousterhout 2018,
  *A Philosophy of Software Design*: what the import session hides from the viewer.
- Mittring 2007, *Finding Next Gen*: screen-space ambient occlusion. J. de Vries, *LearnOpenGL*, "SSAO":
  the normal-oriented hemisphere and the kernel in `mesh_lighting.cpp`.
- Gribb, Hartmann 2001, *Fast Extraction of Viewing Frustum Planes*: `chunkInFrustum`.
- Moller, Trumbore 1997, *Fast, Minimum Storage Ray-Triangle Intersection*: `rayMeetsTriangle`.
