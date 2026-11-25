# Export: file tree and architecture

Reviewed results become files: a photon_path_tracer package, a Blender file or a Wavefront OBJ. This
file is the way into the folder: what each file owns, how data moves through them, how to add a
format, and what must stay true when you change them. The folder computes no knit, no fibres and no
colour. Those are results it is given, and it writes them. Its tubes come from the scene's one mesher
(`scene/appearance/yarn_tube.h`), so the folder never runs the GPU.

## File tree

```text
src/export/
  ARCHITECTURE.md        this file
  CMakeLists.txt         the component knit_export: export.h is its only public header; depends on io, and
                         privately on model, fibres, settings and scene

  The door
  export.h, .cpp         the format table, validateExportName, rendererMemoryGiB, exportScene, report.json

  The neutral scene
  export_scene.h, .cpp   ExportScene: records in memory, parts in files; materials resolved; the
                         part mesh writer and reader; PlyShading; the format writers' declarations

  From reviewed results to the neutral scene
  collect_scene.h        buildExportTube (a scene tube in Y-up metres), PlyModeOptions and attaching a
                         cached ply map set to the package (attachPlySet, setupPlyTextures),
                         PartMaterials, collectScene
  collect_scene.cpp      the stages, in the order the output lists them
  collect_materials.cpp  every material, and which one each part takes

  The formats: each reads an ExportScene and writes files
  format_ppt.cpp         scene.json, the one place the renderer's scene format is written
  format_obj.cpp         .obj and .mtl
  format_blend.cpp       Blender run on blender/ppt_package_to_blend.py
```

**Reading order:** `export.h`, `exportScene` at the end of `export.cpp`, `export_scene.h`,
`format_obj.cpp` (the shortest full format), `collectScene` at the end of `collect_scene.cpp`, then
`format_ppt.cpp`.

## The common recipe

Good exporters look alike, whatever they export:

1. **One door.** A small API: list the formats, ask about them, export. (Assimp's `Exporter`.)
2. **Formats are a table.** Everything said about a format is one row. (Assimp's
   `ExportFormatEntry`, FFmpeg's muxers.)
3. **One neutral scene.** What is exported, independent of any format. N sources and M formats cost
   N + M, not N x M. (`aiScene`, pandoc's AST, LLVM IR.)
4. **One writer per format.** A function from the neutral scene to files, which knows nothing of how
   the scene was made. (trimesh, Open3D.)
5. **Collect apart from write.** The domain is read in one place only.
6. **Atomic output at the edge.** A failed writer leaves nothing behind.
7. **Say what was lost.** Each format states what it cannot carry.

That is enough for most programs. Ours adds needs of its own, and most of the rules below exist for them:

- **Geometry does not fit in memory.** The teddy at Fine with 4 plies is 349,579,840 triangles, about
  10 GB as floats. So the neutral scene holds records, and each part goes to a file as it is made.
- **Blender is another process.** It reads a package on disk, so the package is both a format and the
  input of another format.
- **Bytes are pinned.** Material order, part order, number formats and text or binary open modes are
  all part of the output. The same inputs give the same files.
- **Exports land in user folders.** A folder appears complete or not at all (`io/publish.h`), and a
  cancelled export leaves the previous one.
- **Verdicts and credit travel with every format.** `report.json` and `CREDITS.txt` go with every
  export, and a preset's `.blend` embeds the credit.
- **No part disappears.** What has no yarn is written as its source mesh, by `addPlainMesh` alone, for
  presets and imports alike; a knitted part has its filling behind the yarn (`filling.enabled`).

The door is a header of free functions, not a class, because nothing lives between two calls.
There is no format base class: as everywhere in this project, dispatch is a `switch` over an
`enum class`, so the whole of `writeFormat` fits on a screen. A format is a row, a `case` and a function.

## Layers

An arrow means "includes". The includes run one way, so each file can be read knowing only the files
to its right.

```text
export.cpp ──> collect_scene.h ──> export_scene.h
format_*.cpp ─────────────────────────┘

everything else in src/ ──> export.h                            (the only public header: the compiler holds it)
src/cli/pieces_cli.cpp  ──> collect_scene.h                     (the one named exception: src/cli/CMakeLists.txt)
```

- `export.h` includes only `io/publish.h` and forward-declares the scene types, so `viewer/workflow.h`
  stays light.
- `export_scene.*` and `format_*.cpp` never see `Scene`, `SceneEntry`, `SceneFibres` or the GPU. A format
  cannot read the domain. Of the settings they see only `MaterialSpec` (`settings/material_spec.h`), and of
  the scene nothing: `PlyShading`, what a package needs to shade its tubes with ply maps, is the neutral
  scene's own record (`export_scene.h`), so a new setting or a change in knitting does not recompile the
  formats. `architecture_export_door` refuses any `scene/` include there.
- `knit_check --strand-demo` builds a one-strand package with no scene, no report and no publication,
  so it goes straight to `ExportScene`. A facade does not forbid the subsystem to those who need it:
  `knit_check` is the one program with this folder's private headers (`PRIVATE_ACCESS export` in
  `src/cli/CMakeLists.txt`), and `architecture_export_door` names the same exception.
- Nothing in `src/viewer` or `src/jobs` compares a format's id. They read its row.

## The data

### The neutral scene

`ExportScene` is what every format reads. It has two kinds of content:

```text
in memory                                on disk, under folder()
  projectName, cameras, sky, ground        meshes/<name>.ply         one per mesh part
  baseMaterials, fibreLook, plyShading     strands/<name>.pptstrands one per strand part
  placements: part -> object, group        textures/colour_<n>.png
  bounds, one per part and the whole       textures/...              ply maps
```

`add*` writes a part's file at once and keeps a `PartFile` (name, path, material) and its bounds. A
format walks `meshParts()` and `strandParts()` and reads each part back with `readPartMesh` or
`readStrandsFile`. It holds one part at a time, however large the scene. `sortParts` is the last step of
collecting, so the same inputs give the same files whichever order the parts came in. No format sorts.

Everything in the scene is in Y-up metres. The OBJ alone turns back into an import's own units, by
`ExportScene::importFrame`. Unset, its numbers pass untouched, so a package and its OBJ hold the same values.

### A part file

`writePartMesh` and `readPartMesh` sit side by side and are the only code that knows the layout: a text
header, then binary little-endian floats and lists.

```text
vertex   x y z nx ny nz [u v]     stride 6, or 8 with uvs; u = s / uvRepeat along the yarn, v = theta / 2 pi
face     uchar 3, uint a b c
```

Tubes are triangle meshes because the renderer has no curve primitive for them. Fibres are `.pptstrands`
(`fibres/strands_file.h`): round linear curves with a radius per vertex.

### A material, from choice to resolved form

A material passes through three hands, and each owns one question:

```text
PartMaterials        which entry does this part take?       adds baseMaterials, says the id
ExportScene          which entries are referenced?          materials(): resolved ExportMaterial
format_ppt.cpp       how does the renderer spell it?        writeMaterial, schema 12
```

`materials()` lists the base entries a part references, then the ply-mapped yarns, then the fibre
materials, then the ground. An entry left out can still colour a variant or a fibre, which is why the
table holds more than the file does.

A fibre material is resolved late, because it depends on its yarn: `fibreMaterialSpec` takes the yarn's
colour through the hair maths. The renderer's hair colour is that of a dense mass, so one thin fibre
of the yarn's colour is almost clear (Chiang et al. 2016, Eq. 9). Fuzz takes the geometric mean of the
dense mass and one fibre; a dense coat stays close to the yarn.

`MaterialKind` says which layers a material has: `Surface`, `Hair`, `PlyMappedYarn`, a yarn colour
under the ply maps, or `BakedYarn`, one yarn's tube under the maps baked from its fibres
(`fibres/yarn_bake.h`; `ExportScene::bakedMaps` holds their paths and the renderer's constants). Spell the kind, not `spec.model`: a ply-mapped yarn can have a fibre model and is
still written as a surface. A baked yarn is written in one of two materials (`BakedMaps::fabric`, from
`yarn_body.map_material`), both with the highlight of the roughness and specular maps (`writeBakedHighlight`, the
renderer's `specular` block): the yarn's fabric (sheen, the highlight under it, translucency; the same `writeSheenKeys`
as a surface) under the cavity and normal maps, or a plastic with the highlight alone, opaque.

### The scene graph

`Placement` puts a part under its source object, in a child group by fibre type. `format_ppt.cpp` builds
one root over every object, so the package imports as one object:

```text
project
├─ object (one bear)          pivot at its bottom centre: the editor turns a bear about its feet
│  ├─ tube parts
│  ├─ fibres_oriented
│  └─ fibres_random
└─ Ground
```

## How an export runs

`exportScene` in `export.cpp` is the whole story in about 40 lines.

```text
exportScene
│
├─ exportRefusal         before anything is written: format, close-up part
├─ recoverPublished, Publisher::begin(folder, policy)
├─ partsFolder           staging for the package, else staging/.package
│
├─ collectScene          ExportScene filled, sorted, strands read back
│    addMaterials, addPlainPrimitives, addPart per part,
│    addCloseupCameras, addFibres, sortParts, readStrandsBack
│
├─ frameScene            refused if nothing can be framed: every format
├─ writeSidecars         report.json, CREDITS.txt into partsFolder (Blender reads the credit)
│
├─ writeFormat           switch (format.kind)
│    Ppt    writeSceneJson
│    Obj    writeObjFile
│    Blend  writeSceneJson, then writeBlendFile
│
├─ partsFolder != staging: writeSidecars into staging, remove partsFolder
└─ Publisher::commit
```

The caller supplies a `PublishPolicy`: `IntoFolder` for the viewer, whose user chose a folder that may
hold other files, and `Strict` for `knit_check`. Eligibility stays with the caller (`exportEligibility`,
`jobs/results.h`), because the command line and the tests have no `KnitResult` to ask.

### Collect

`collectScene` runs stages over one `Collector`, the state they share. The order is not free:
`colour_<n>.png`, `tex<n>_` and the material table follow the order of the calls.

| Stage | What it adds |
|---|---|
| `addMaterials` | every material entry, then the ply map set, its files and the tube's uvs |
| `addPlainPrimitives` | surfaces that were never for knitting, as their source meshes |
| `addPart` | the part's tubes and filling; not knitted, its surface as a plain mesh |
| `addCloseupCameras` | `closeup_<part>` and `medium_<part>`, for the article's renders |
| `addFibres` | each part's fibres below the budget's key: one part per population and colour |
| `readStrandsBack` | every strand file read back, before anything is published |

`TubeStyle` is built in three steps that must all reach the tube: ring and sides from the settings,
compression from the final surface, uvs from the map set. A knitted part with colourwork becomes one
mesh per stitch colour: `addTubeGroups` splits it by the stitch that owns each triangle.

## The formats

A row of `exportFormats()` holds everything said about a format: names, what it holds, what it leaves
out, its coordinates, the extensions of the files that take the export's name, whether it can carry
fibres, and why it may be unavailable. The Export step draws these rows as they are.

**Package (`ppt`).** The parts are already in the output, so `writeSceneJson` adds scene.json. It writes
one `PhotonScene`, v2, or v5 with strand parts because the renderer requires v5 for root `strands[]`.
It adds a default front camera, a grey ground quad under the parts, and a uniform analytical sky.
Nothing is left out.

**OBJ (`obj`).** `writeObjFile` writes the `.mtl`, then one object per mesh and one per strand part, with
fibres as polylines. v and vt keep separate running numbers, because only meshes with uvs write vt. It
leaves out the fibres' radii, which OBJ lines cannot hold (each object's comment gives the mean), the
ply maps and the ground.

**Blender (`blend`).** Blender writes the file itself: `writeBlendFile` runs it headless on the script.
So the format starts from a package folder, written first by `writeSceneJson`, and not from an
`ExportScene`. Its tests hand it folders they wrote themselves.

## Adding a format

A new format needs six small edits, and none in `viewer/`, `jobs/` or `cli/`:

1. `ExportFormatKind`: one value (`export.h`).
2. `exportFormats()`: one row (`export.cpp`).
3. `writeFormat`: one `case` (`export.cpp`).
4. `export_scene.h`: one declaration.
5. `format_x.cpp`: the writer.
6. `src/CMakeLists.txt`: one line.

The writer is a walk over the parts. `readPartMesh` and `readStrandsFile` hand it one part at a time,
and `materials()` gives it every colour:

```cpp
// format_stats.cpp: a text file of the scene's parts.
bool writeStatsFile(const ExportScene& scene, const FormatTarget& target, std::string* err) {
    std::ofstream f(target.folder + "/" + target.baseName + ".txt");
    for (const PartFile& part : scene.meshParts()) {
        const std::optional<StoredMesh> mesh = readPartMesh(scene.pathOf(part), err);
        if (!mesh) return false;
        f << part.name << " " << part.material << " " << mesh->tri.size() / 3 << " triangles\n";
    }
    for (const PartFile& part : scene.strandParts()) {
        StrandSetF strands;
        if (!readStrandsFile(scene.pathOf(part), strands, err)) return false;
        f << part.name << " " << part.material << " " << strands.strands() << " strands\n";
    }
    return f ? true : fail(err, "cannot write the stats file");
}
```

MSVC's `/W4` does not warn about a missing `case`, so `export_every_format_writes` holds step 3: every
row must reach a writer and leave `report.json`.

If the format needs what only `format_ppt.cpp` computes today, the ground quad, say, move that into
`export_scene.cpp` when the second caller arrives, not before. `frameScene` went that way.

## Checks

A check sits where its answer is first known, and each runs once.

- **`exportEligibility`** (`jobs/results.h`, a rule over the reviewed results, which this folder never
  sees) decides before a job starts: Ready, Diagnostic or Blocked. The Export step and the job call the
  same function on the same inputs, so the button and the refusal never disagree. It has no refusal of
  its own to report: every outcome is a verdict; `ExportOptions::status` carries it into report.json.
- **`exportRefusal`** is what no caller may skip: a close-up of a part the scene does not have is refused
  before anything is written.
- **`frameScene`** runs for every format, because a scene with no finite bounds would otherwise publish an
  OBJ of infinite vertices.
- **`readStrandsBack`** reads every strand file before anything is published.
- **The strand rules** (`checkStrands`) refuse a strand part the renderer could not load.

## What must stay true

| Rule | Held by |
|---|---|
| The program includes `export.h` only; a format never sees the domain; no format id outside the table. | `architecture_export_door` |
| Every format in the table reaches a writer and leaves `report.json`. | `export_every_format_writes` |
| No written byte changes by accident. A moved pin is a failure, never a re-pin. | `io_pin_scene_json`, `io_pin_obj_files`, `io_pin_report_json` |
| `scene.json`'s material entry and schema version are written in one function. | `writeMaterial`, `ppt_material_schema12` |
| A preset's `.blend` holds its credit: `CREDITS.txt` is in the parts folder before Blender runs. | `blender_export_keeps_the_credit` |
| A preset writes its unknitted parts as source meshes, its patches as moved copies, and a filling per knitted part. | `preset_export_writes_its_details` |
| The fibres the export writes are the given result, below its key. | `workflow_plane_knit_fibres_export` |
| `sortParts` and the strand read-back end collecting, for every format. | `collectScene`, `export_every_format_writes` |
| Stage order sets file names and the material table. | `io_pin_scene_json`, the gate of a change |
| Without an import the OBJ's numbers pass untouched: no renormalised normal, unit exactly 1. | `io_pin_obj_files` |
| `.mtl` is sorted by id, a duplicate id takes the last entry, an unknown id a light grey. | `io_pin_obj_files` |
| `Kd` keeps `scene.json`'s rounding, by `colourAsWritten`, on the material's kind. | `io_pin_obj_files` |
| A mesh named `_ground` stays out of the OBJ, like the ground quad. | `io_pin_obj_files` |
| `part_fibre_shapes` goes on wool fibre parts only, by `isWoolFibreMaterial`, not by hair material. | `ppt_export_fibre_groups`, `io_pin_scene_json` |
| JSON, `.mtl` and credits open in text mode, `.ply` and `.obj` in binary: CRLF on Windows is in the bytes. | `io_pin_report_json`, `io_pin_obj_files` |
| A report stays valid JSON whatever the settings hold. | `export_report_escapes_the_settings` |

## Where to change what

| To change | Touch |
|---|---|
| What a format says in the interface | its row in `exportFormats()` |
| Whether an export may run | `exportEligibility` (`jobs/results.h`); a refusal no caller may skip, `exportRefusal` |
| Which material a part takes | `PartMaterials` in `collect_materials.cpp` |
| A fibre's colour | `fibreMaterialSpec` and the hair colour functions in `export_scene.cpp` |
| The renderer's material schema | `writeMaterial` in `format_ppt.cpp`, with the renderer's own migration |
| What a tube looks like | `tubeStyleOf` and `buildExportTube`; `TubeStyle` is `scene/appearance/yarn_tube.h`; the look's settings live in `settings` |
| How parts are grouped | `Placement` in `collect_scene.cpp`, `buildGroups` in `format_ppt.cpp` |
| The `.ply` layout | `writePartMesh` and `readPartMesh` together; `ppt_package_to_blend.py` reads it too |
| What `report.json` holds | `writeReport` in `export.cpp` |
| A format's bytes on purpose | Bump its version where it has one, re-pin, and prove the re-pin by undoing the change and hashing back |

## Literature

- Gamma et al. 1994, *Design Patterns*: Facade (one door, the subsystem does not know it) and Builder
  (one construction process, several representations): `export.h` and `collectScene`.
- Fowler 2002, *Patterns of Enterprise Application Architecture*: Two Step View, domain to logical form,
  logical form to output.
- Lattner, Adve 2004, *LLVM: A Compilation Framework*: one intermediate form between N sources and M targets.
- Pharr, Jakob, Humphreys, *Physically Based Rendering* (4th ed.): a scene streamed into a target, not built in
  memory first.
- Assimp's `Exporter`, trimesh's `exchange`, Open3D's `io`: a table and one function per format.
- Chiang, Bitterli, Tappan, Burley 2016, *A Practical and Controllable Hair and Fur Model*: the absorption
  behind `hairChannelOfFibre`.
