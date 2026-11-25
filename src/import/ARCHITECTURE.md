# Import: file tree and architecture

From a mesh file somebody made for something else to an object this program can knit. This file is
the way into the folder: what each file owns, how data moves through them, and what must stay true
when you change them. The Import step's screen is in `src/viewer`, its thread in `src/jobs`, and the
reader in `src/io`.

## File tree

```text
src/import/
  ARCHITECTURE.md        this file
  CMakeLists.txt         the component knit_import (src/ARCHITECTURE.md): the analysis, the preparation, the
                         pipeline and the store public, the rest private; depends on model and io;
                         knit_flatten inside (wrap_flatten.cpp, wrap_quad_split.cpp and the SLIM solver)

  What a file is, and what is made of it (public)
  import_analysis.h, .cpp
                         after Open: the sanity check and the nodes
  prepare_import.h, .cpp the options, the repairs, library space, the identity
  import_pipeline.h, .cpp
                         one preparation top to bottom: prepareAndStore, the rebuild's way in and out
  import_store.h, .cpp   the cache entry: mesh.kmesh, plain.kmesh, textures, and the way back as a scene
  repair_rules.h, .cpp   the rules the repairs apply and the analysis counts by (private)

  The rebuild, stage by stage (private)
  wrap_surface.h, .cpp   the recipe and the whole run; the three files are written here
  wrap_stages.h          stages 1 to 11 in run order with their types: the table of contents
  wrap_grid.cpp          1 to 3: voxels, the seal, the blurred field
  wrap_nets.cpp          4, 5: surface nets and Taubin smoothing
  wrap_remesh.cpp        6: the voxel remesh through third_party/vdb, one outward winding, the sorted edges
  wrap_quad_split.cpp    how every stage splits a quad into two triangles
  wrap_tidy.cpp          7: fragments dropped, one closed manifold kept
  wrap_charts.cpp        8: a few wide charts
  wrap_flatten.cpp       9: each chart flattened by SLIM from a Tutte start, on a mollified rest shape
  wrap_pack.cpp          10: charts turned to stand along v and packed
  wrap_bake.cpp          11: the colours baked by rays, the closest surface where a ray misses
```

**Reading order:**

1. `io/mesh_read.h`: what a file holds.
2. `import_analysis.h`: what the Import step tells the user after Open.
3. `prepare_import.h`: what the user chooses, and what a preparation makes of the file.
4. `import_pipeline.cpp` from `prepareAndStore` upwards: the order of a preparation.
5. `import_store.h`: what is kept, and how it comes back as a scene.

The `wrap_*` files are only needed to change the rebuild: then `wrap_surface.cpp` first, with
`wrap_stages.h` as the table of contents.

## The common recipe

Assimp, three.js's loaders, the import pipelines of game engines and asset databases share one shape:

1. **A reader fills one neutral result, and knows nothing else.** Assimp's importer only fills an
   `aiScene`, Unreal Interchange's translators only fill a container of neutral nodes; ours only fills
   an `ImportedMesh`.
2. **Nodes name parts of flat data.** An `aiNode` refers to meshes by index, and a tinyobjloader
   `shape_t` is a name with its faces over file-wide vertex arrays. Here `nodeOfTriangle` over one
   indexed mesh. OBJ and PLY have no transforms and no hierarchy, so the tree is one level deep.
3. **Passes run afterwards, one by one.** Assimp's post-processing steps (`JoinIdenticalVertices`,
   `FindDegenerates`, `FixInfacingNormals`) run in a fixed order; here `finishNodes`, then
   `analyzeImport`, and in a preparation the repairs in their order.
4. **The user's choices live beside the data, keyed by name.** Godot's per-node `import/skip_import`,
   and a USD stage's population mask, a set of paths that limits what is loaded. Here
   `ImportOptions::plainNodes`.
5. **Source plus choices identify the product.** Unity's asset database keys an artifact by the source,
   the importer's settings and its version; change any and it is another product. Here
   `importIdentity`, which is why a preview can stay cheap and a stored result can be trusted.

What we decline from them: step classes chosen by flags, plugin registries, post-import scripts and
`.meta` files beside the user's file. Our steps are functions called in order, and nothing is written
next to the source.

We add three things, and each has a section below: a file without names gets its loose parts as nodes;
the nodes left out stay in the object beside a surface that was rebuilt; and the rebuild can run in
another program, which reads files only.

## Layers

What each file (header and source, by name) may include from this folder. The includes run one way, and
`architecture_import_layers` holds exactly this table; removing an include stays legal.

```text
public    prepare_import    repair_rules
          import_analysis   prepare_import, repair_rules
          import_store      prepare_import
          import_pipeline   prepare_import, import_store, wrap_surface

private   repair_rules      nothing
          wrap_stages       nothing
          wrap_grid, wrap_nets, wrap_remesh, wrap_quad_split, wrap_tidy, wrap_charts, wrap_flatten, wrap_pack, wrap_bake
                            wrap_stages
          wrap_surface      wrap_stages, prepare_import
```

- `prepare_import.h` and `import_analysis.h` include no other file of the folder; `prepare_import.cpp`
  and `import_analysis.cpp` both include `repair_rules.h`, so the analysis counts a defect by the rule
  the repair applies. The analysis includes the preparation's header, never the other way.
- Nothing here includes `knitting/`: the import makes a surface to knit and never asks the knit
  whether it will.
- Third-party code is named per file: OpenVDB only behind `third_party/vdb/level_set_remesh.h`, which only
  `wrap_remesh.cpp` includes; SLIM and its Eigen only in `wrap_flatten.cpp`.
- A public header includes no private one, so jobs and scene see neither the rules nor the engine.
- `wrap_flatten.cpp` and `wrap_quad_split.cpp` are their own library, `knit_flatten`, because the flatten
  builds on the solver's pinned Eigen; `quadSplit` lives there so that the flatten and the stages above it
  share one definition. No Eigen type appears in a header. The solver's eight files are vendored MPL-2.0 code
  in `third_party/slim/src`, with our changes in `third_party/slim/CHANGES.md`.
- Only `import_pipeline.cpp` calls `storeImport` in the program (tests store entries too). The rebuild
  makes a surface and knows nothing of the cache, which is what keeps a preparation readable in one place.
- This folder knows no viewer and no session. `jobs/import_session.*` runs it on its thread, and tests
  call `prepareAndStore` directly. Some stages split their loops over the executor (`model/parallel.h`)
  with results that do not depend on the number of threads.

## The data

### What a file holds

`ImportedMesh` (`io/mesh_read.h`) is one indexed mesh in file units and file axes:

```text
pos        [ v0 v1 v2 v3 v4 v5 v6 v7 ... ]
tri        [ t0 t1 t2 t3 ... ]                  three indices into pos
uv         [ t0 t1 t2 t3 ... ]                  three per triangle: the file's corners, or empty
triMaterial[ 0  0  1  1  ... ]                  index into materials
nodeOfTriangle [ 0 0 0 1 ... ]                  index into nodes

nodes      [ {"bow", triangles, area, inTheRound}, {"brown_bear_eyes", ...}, ... ]
```

A **node** is what the Import step calls a component. The word differs on purpose: "component"
already means a connected component in `ComponentCheck`, `ScenePart` and `knit.coverage`.

The reader records names by the rule of three.js's `OBJLoader`: `o` and `g` are one rule, and a name
without faces is no node, so a stray `g default` makes no row. Faces under one name are one node
wherever they stand in the file.

`finishNodes` is the pass after reading. A file naming fewer than two nodes gets its connected pieces
as nodes, called `Part 1`, `Part 2`, on the welded mesh, because a PLY repeats a corner for every
face; a single named node that is one piece keeps its name. Then every node's triangles and area are
counted. `bears.obj` has 9 named objects and 95 connected pieces, and the 9 are what a person means
by "the eyes, the nose, the bow".

`analyzeImport` counts the file's defects with the repairs' own rules (`repair_rules.h`): coincident
vertices, zero-area, duplicate and unreferenced elements, boundary, non-manifold and inconsistent
edges. The counts describe the file as read, so they can differ from what the repairs report after
welding. The repairs act only on the file's own mesh, that is when the knit is skipped (below).

### The two surfaces

A knitted import is always rebuilt (`Surface::Wrapped`): the ticked nodes become one closed surface at
stitch scale with its own charts. The file's own mesh (`Surface::KeepFile`, `import.surface = file`)
is never a choice of the dialog. It serves an import whose knitting is skipped (hair grows on the
file's surface), the unticked nodes, and scripted `knit_check` runs that knit a file's own uv (the
swatch's open plane, the test fixtures). Until 2026-11-11 the dialog offered it for knitting too, with
generated charts (xatlas) and closed toy pieces knitted in the round.

### What the user chose

`ImportOptions` is everything that decides the result. `plainNodes` is a set of **names**, so it stays
true when a file changes, and a name the file no longer has matches nothing. `Part 3` is a position
and not a name, so after a nameless source file changes it can mean another piece.

`importIdentity` hashes the importer's version, a frozen token (it named xatlas), `importOptionsText`, the source
file's bytes, every dependency's bytes and the names of missing dependencies. Options write a line
only when they differ from what older imports did: a selection adds `plain_nodes`, a rebuilt surface
its `wrap_*` lines with `kWrapVersion`. So every older keep-file identity holds; the line `charts = keep`
stays in every identity for that reason, though charts are no longer a choice. A rebuilt surface names no
engine, because the import has one (since 2026-11-12): imports stored with a `wrap_engine` line are
rebuilt once. The constants live in
`prepare_import.cpp`, beside the two functions that read them.

In the settings the selection is `import.plain_nodes`, with `|` between names, and the reader turns a
`|` in a name into `_`. It is no stage input, because `import.identity` already depends on it. If it
were one, adding the row would have moved every import's keys.

### The cache entry

```text
out/cache/imports/<identity>/
  mesh.kmesh    the surface to knit: library mm, Z up, knit uv, colour uv, vertex colours, materials
  plain.kmesh   the nodes left plain, in the same form and the same space. None: no file
  textures/     copies of the base colour images of both, named <index>_<file name>
  import.txt    source path, options, what the preparation changed
```

The entry is published as one directory (`io/publish.h`), and is never deleted as a cache: saved
settings refer to it. `loadImportedScene` turns it into a `SourceScene`: one `Class::Knit` primitive
per material of `mesh.kmesh`, then one `Class::Plain` primitive per material of `plain.kmesh`.
From there on an import is like any other scene, and plain primitives are what the viewer and all
three export formats already handle, with their textures.

## How a file becomes an object to knit

```text
stage 1  a file is dropped or opened
         ImportSession::open queues the file; readAndAnalyse runs on the session's thread:
         ├─ readMeshFile                  limits first: 200 MiB of mesh data, 2 M triangles
         ├─ finishNodes
         ├─ the node limit                over 200 components: refused, with its message
         ├─ analyzeImport                 what is wrong, and what can be repaired
         ├─ loadImportTextures            the base colour images, under the same limits
         └─ makePreview                   one display mesh per node and material

stage 2  the Import tab shows the object on the live 3D view
         ticks and clicks edit ImportOptions::plainNodes; the choices below them edit the rest
         a click: viewRay -> rayInFile -> nodeUnderPixel -> nodeOfTriangle -> pickNode, toggleNode

         Import
         prepareAndStore                  import_pipeline.cpp, on the session's thread
         ├─ a stored rebuild answers at once: its identity fixes its result
         ├─ splitByNodes                  only when something is left plain
         ├─ prepareSurface                of the ticked nodes
         │    knitted: rebuildSurface, then readRebuiltSurface
         │    knitting skipped (the file's mesh): prepareImport
         ├─ preparePlainNodes             the rest as the file has it; both moved to the whole object's place
         └─ storeImport                   one published entry

         the Knit step opens, as for any object
```

The repairs of `prepareImport` run in this order, because each relies on the one before: merge
coincident vertices, drop what became degenerate, drop duplicates, compact, orient.

**Failures.** A step that can fail returns `[[nodiscard]] bool`, takes its inputs first, then what it
fills, then `std::string* err` for a message the user reads, then the progress callback last
(`readMeshFile`, `prepareImport`, `prepareAndStore`, `wrapFileToFolder`). A cancelled job throws
`Cancelled` from a checkpoint; the session catches it, and `RebuildFolder` removes the temporary folder
however the job ends.

### The selection

Two things make a partly knitted object work.

**Everything stands where the whole object stands.** `toLibrarySpace` centres a mesh over the origin
with its lowest point on the floor. Done for each part alone, the parts would drift apart, so
`preparePlainNodes` takes the offset from the whole file, the one the preview uses, and moves both
meshes to it. A rebuilt surface keeps the file's frame, so one offset lines it
up with what was left plain. With nothing left plain no code runs here, and the output is the same,
byte for byte.

**The rebuild sees only the ticked nodes.** `splitByNodes` hands it the knitted side as a mesh of its
own, so nothing in the rebuild knows about selections.

**Size decides whether a plain detail shows.** The yarn stands off the surface by a fixed length, whatever
the object's size, and a rebuilt surface seals what it does not hold. Measured on `bears.obj` with the
eyes, noses, mouths and bow plain: at 300 mm they all show, and at 140 mm, the unit the import suggests, the
eyes and noses are nearly buried. The bow stands proud at both.

`RebuildFolder` lives in `prepareAndStore` and outlives `storeImport`, because the rebuilt surface's
texture is copied from it when the entry is stored.

### The rebuild

The rebuild writes three files into that folder (`kRebuiltSurfaceName`): `wrapped.obj` in mm, Y up,
with its charts as uv, `wrapped.mtl`, and `wrapped_basecolor.png`. `readRebuiltSurface` takes them back
through the ordinary reader and repairs, so a rebuilt import is an ordinary prepared import.

It is this application's own engine and, since 2026-11-12, the only one; settings that still name an
engine or a recipe (`import.wrap_engine`, `import.wrap_recipe`) load and are ignored. `wrapFileToFolder` is
the whole run, and `wrap_stages.h` lists the stages in this order:

```text
1-3   wrap_grid.cpp          triangles to voxels, sealed, flooded from outside, blurred
4-5   wrap_nets.cpp          surface nets, Taubin smoothing
6     wrap_remesh.cpp        the voxel remesh through OpenVDB, then one outward winding per part
7     wrap_tidy.cpp          fragments dropped; a large second part or an open edge stops the import
8     wrap_charts.cpp        charts grown, opened into discs, judged, split or slit
9     wrap_flatten.cpp       per chart: a mollified rest shape, a Tutte start map, SLIM; one common scale
10    wrap_pack.cpp          charts turned so their smallest rectangle stands along v, laid out by corners
11    wrap_bake.cpp          colours by rays from outside, the closest surface as the fallback
```

Stages from 6 on compute in doubles, and no result depends on the number of threads.
Each `wrap_*` file states its stage's method in its first lines, with the papers it follows where there are any. Only
the base colour is baked, as stored.

## What must stay true

| Rule | Held by |
|---|---|
| A file with every node ticked gives the same identity, keys and stored bytes as before selections existed. | `import_identity_unchanged`, `import_stage_keys_unchanged`, `io_pin_import_identity`, `import_selection_identity_and_keys` |
| Two different selections are two entries. | `import_selection_identity_and_keys` |
| Knitted and plain parts stand where they stood in the file. | `import_selection_keep_file`, `import_selection_built_in` |
| A plain node keeps its texture, and reaches every export. | `import_selection_keep_file`, `export_selection_keeps_the_plain_node` |
| Every triangle of a file lands on exactly one side of a split. | `import_split_by_nodes` |
| More than 200 components or 200 MiB of mesh data is refused on load, with a message. | `import_session_refuses_too_many_nodes`, `import_oversize_rejected_before_parse` |
| A cancelled job leaves no folder behind. | `import_wrap_session_cancel` |
| A failed rebuild never falls back to the file's mesh by itself; the user may press "Keep the file's own surface" (Knit workflow; the Mark step and the guided stitch mesh work on it). | `ui_import_failed_wrap_keeps_the_choice`, `ui_import_failed_wrap_keep_file_on_request` |
| A knitted import is rebuilt; the dialog offers no other surface. | `ui_import_default_rebuilds`, `ui_import_scripted_skip_keeps_charts`, `ui_import_rules_without_a_frame` |
| Old settings with `import.charts` still load. | `import_settings_old_file` |
| Old settings with `import.wrap_engine` or `import.wrap_recipe` still load and import as without them. | `import_retired_engine_keys_load` |
| Each rebuild stage holds its properties: consistent winding, discs, no flipped or overlapping uv, uv in the unit square, the same bytes at any thread count. | `tests_wrap_properties.cpp`, `tests_wrap_flatten.cpp`, fixtures in `tests/data/wrap_fixtures` |
| The file format of `mesh.kmesh` does not change by accident. | `io_pin_*` in `tests_io_pins.cpp` |
| The includes follow the layers above. | `architecture_import_layers` |
| This file names every file in the folder, and only those. | `architecture_import_doc_complete` |

## Where to change what

| To change | Touch |
|---|---|
| What a preparation does, or its order | `prepareAndStore` and the functions just above it |
| A limit | `ImportLimits`, the text in `viewer/panel_import.cpp` |
| What counts as a node | the OBJ reader's `nodeOfFace`, or `finishNodes` |
| What the analysis reports | `analyzeImport` and its `add*Issues`; a rule it shares with a repair, in `repair_rules.cpp` |
| An option | `ImportOptions`, `importOptionsText` with a line only when it is not the default, `importOptionsFrom` and `applyTo` in `jobs/import_session.cpp`, a row in `settings/scene_settings.cpp` |
| What is stored | `writeKmesh`, `readKmesh`, then re-pin `io_pin_*` on purpose |
| A repair | `prepareImport`, an issue in `addRepairIssues`, a line in `viewer/panel_import.cpp` |
| A rebuild stage | its `wrap_*` file, its property tests, and `kWrapVersion` in `prepare_import.cpp` |
| The look of an unticked component | `kPlainNodeGrey` in `viewer/renderer.cpp` |

## Known limits

- Three ray-triangle tests stay apart, because their contracts differ: `wrap_bake.cpp` (the bake,
  with a running limit), the knit's coverage (a segment) and the viewer's click (a ray from the
  camera).

## Literature

- Assimp, *The Asset-Importer-Lib Documentation*: a neutral `aiScene`, nodes by index, readers apart
  from post-processing steps.
- three.js, `OBJLoader.js`: `o` and `g` as one rule, and names without faces skipped.
- tinyobjloader, `tiny_obj_loader.h`: file-wide attributes under named shapes.
- Godot, *Import configuration*, and `resource_importer_scene.cpp`: per-node import options,
  `import/skip_import`.
- Unity, *The Asset Database* and *AssetImporter*: artifacts keyed by the source, the importer's
  settings and its version.
- Unreal Engine, *Interchange framework*: translators, pipelines and factories.
- OpenUSD, `UsdStagePopulationMask`: a set of paths that limits what a stage populates.
- Möller, Trumbore 1997, *Fast, Minimum Storage Ray-Triangle Intersection*: the click, and `wrap_bake`.
- Taubin 1995, *A Signal Processing Approach to Fair Surface Design*: `smoothTaubin`.
- Rabinovich et al. 2017, *Scalable Locally Injective Mappings*: `flattenFaces`.
- Sharp, Crane 2020, *A Laplacian for Nonmanifold Triangle Meshes*: the mollified rest shape.
- Floater 2003, *Mean value coordinates*; Tutte 1963, *How to draw a graph*: the start map.
- Kahan 2014, *Miscalculating Area and Angles of a Needle-like Triangle*: triangle areas from lengths.
- Andrew 1979; Freeman, Shapira 1975; Toussaint 1983: hull and smallest rectangle in `packCharts`.
- Ericson 2005, *Real-Time Collision Detection* 5.1.5: `TriangleBvh::closestPoint`.
- Thürmer, Wüthrich 1998, *Computing vertex normals from polygonal facets*: `vertexNormalsOf`.
- Museth 2013, *VDB*, and the OpenVDB documentation: the level set remesh.
