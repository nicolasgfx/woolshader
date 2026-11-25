# IO: file tree and architecture

Bytes and files: what every format stands on, the publication protocol, and the formats several
modules share, where the caches live, and the Blender process of the `.blend` export. This file is the way into the folder:
what each file owns, how data moves through them, and what must stay true when you change them. A
format with one user lives with that user (see *Every format*), so this folder knows no settings, no
scenes, no knitting and no GPU.

## File tree

```text
src/io/
  ARCHITECTURE.md       this file
  CMakeLists.txt        the component knit_io: every header public, depends on model only; cgltf and stb
                        private; the one place KNIT_SOURCE_DIR is defined

  What every format stands on
  files.h, .cpp         ByteWriter, ByteReader, whole-file read, atomic replace, file hash, project paths,
                        sourceTreeDir (the only reader of KNIT_SOURCE_DIR)
  cache_paths.h, .cpp   out/cache/<name>, the import store's folder names (identity in hex), the Mark step's autosave root

  Getting files into place
  publish.h, .cpp       Publisher: a folder appears complete or not at all

  Formats several modules share
  image_io.h, .cpp      ColourImage, PNG and EXR, decoding through stb_image, colours at a uv
  mesh_read.h, .cpp     foreign meshes: OBJ with MTL, PLY, their textures
  marking_file.h, .cpp  "knit-markings v4" (reads v1 to v3): surface markings as text (the Mark step and stitchmesh_lab)
  stitch_mesh_obj.h, .cpp  a built stitch mesh as OBJ (Export Stitch Mesh: polygons kept, file units, row and kind
                        groups) and a polygon OBJ read back (tests, checks)
  glb_load.h, .cpp      GLB into library space, and SourceScene, what every scene is loaded into
  glb_selection.h       GlbSelection: what a scene takes from its GLB file, apart from the reader

  Other programs
  blender_process.h, .cpp  finding Blender and running the .blend export's script, with a time cap and cancel
  child_process.h, .cpp  a program in the background: its log, its stdin, terminate and the confirmed exit
                        (Blender for the .blend export, and the live preview's `ppt serve`)
  preview_pipe.h, .cpp  the live preview's named pipe: the renderer's STATUS and FRAME messages
```

**Reading order:** `files.h`, `publish.h`, `image_io.h`, `mesh_read.h`, `glb_load.h`.

## The common recipe

Clean file code has the same shape in most systems, whatever the format:

1. **Bytes.** Moving bytes between memory and disk is one small layer with no format in it.
2. **Codec.** A format turns a record into bytes and back, in memory, with no disk and no policy.
3. **Store.** A place decides what is kept under which name, using the codec and the layer below.
4. **Replace, never overwrite.** New content is written beside the old and renamed over it, so a
   reader never meets half a file.
5. **Refuse at the boundary.** A reader returns a complete record or a reason. Input is hostile:
   every count is checked against what is left before anything is allocated.

That is enough for most programs. Ours adds four needs, and most of the rules below exist for them:

- **Identities are frozen.** Stored keys hash exact text, so an innocent change of order or number
  format makes every saved result stale.
- **Doubles round-trip bit for bit.** A saved knit must verify against its digest.
- **Packages are byte-identical.** The same inputs must give the same files, whatever order parts
  finish in.
- **Exports land in user folders.** Nothing may be deleted or replaced that we did not create.

## Layers

An arrow means "includes". The includes run one way, so each file can be read knowing only the
files to its right.

```text
glb_load.h ──> glb_selection.h
     └──────────┐
mesh_read.h ────┴──> image_io.h

publish.h, files.h, cache_paths.h, blender_process.h, child_process.h, preview_pipe.h    include nothing of this folder
```

Every `.cpp` includes `files.h` for its bytes, and `glb_load.cpp` also takes the limits of `mesh_read.h`.

- `src/io` includes only `model` and itself. The compiler holds that: no other folder is on its include
  path (`CMakeLists.txt`). Nothing here knits or runs the GPU, and nothing knows `SceneSettings` or
  `SceneEntry`.
- The knit cache stores a `PieceKnit`, so it lives with its user, `scene/knit/knit_cache`.
- `loadGlb` takes a `GlbSelection`, the ten fields it reads, so a registry entry holds one and this
  folder never sees the entry. The modules above decide what a scene is. A preset's mesh selection is
  `plainNodes`, by glTF node name as an import's `import.plain_nodes`, decided before the material lists;
  `listGlbNodes` lists the nodes and their selection for the Import step from the file's JSON alone.
- `SourceScene` lives in `glb_load.h` because the GLB reader was first. A prepared import is loaded
  into the same type (`import/import_store.h`), so the rest of the program sees one kind of source.

## The data

### A binary record

The knit and fibre formats share one shape. The magic is written twice, so a truncated file is
refused by its end as well as by a count (`.kmesh` opens with a magic and a version and checks its
arrays instead; `tiled-knit` closes with a checksum of the whole record):

```text
"knit-cache v2\n"   the magic: which format, which version
fields              written in order by ByteWriter, read in the same order by ByteReader
"knit-cache v2\n"   the end marker
```

Widths are named in the code (`u8`, `i32`, `u32`, `i64`, `u64`, `f64`) because formats differ: `.kmesh`
writes 32-bit indices and string lengths, the knit formats 64-bit integers everywhere. There is
deliberately no generic `write(value)`: a change of type would silently change the file.

`ByteReader` makes a reader short. Past the end every read gives zero and `ok()` stays false, so a
reader states its fields top to bottom and asks once at the end. `count(bytesPerItem)` reads a count
and refuses one the remaining bytes cannot hold, so a corrupt file never allocates beyond its size.

### Every format

| File | Format | Lives in | Why there |
|---|---|---|---|
| `.knitcache` | `knit-cache v2`, binary | `scene/knit/knit_cache` | one user, the scene's knit stage; it stores knitting's types |
| `knit.bin` | `tiled-knit v1`, record v4; reads v2 to v4 | `scene/knit/saved_knit` | one user; shares `writeCurves` |
| `.fibrecache` | `fibre-cache v6`, binary; older refused | `scene/fibres/fibre_files` | one user |
| `mesh.kmesh` | `KMESH1`, version 1, binary | `import/import_store` | one user, with what it stores |
| settings | `knit-settings v1`, text | `settings/scene_settings` | beside the struct it serialises |
| `.preset` | `knit-preset v1`, text | `settings/file_presets` | shares the key = value grammar |
| `scene.json`, `.ply` | PhotonScene v2 or v5, schema 12 | `export/format_ppt`, `export/export_scene` | one user |
| `.pptstrands` | `PPTSTRND`, version 1, binary | `fibres/strands_file` | below this folder: fibres write it |
| PNG, EXR | stored blocks, no compression | `io/image_io` | used everywhere |
| OBJ, MTL, PLY, GLB | foreign, read only | `io/mesh_read`, `io/glb_load` | the outside world |

The rule: a format used by one module lives in that module, on top of `files.h`. This folder holds
the primitives and the formats several modules use.

## How a mesh arrives

```text
readMeshFile                      limits first: file size before opening
├─ readObjFile                    statement by statement, counts checked as they grow
│    readObjVertex, readObjFace, readMtl
└─ readPlyFile
     parsePlyHeader
     checkPlyLimits               counts against limits, payload against file size, before allocating
     findPlyColumns               which property holds what
     readPlyElements              one PlyValues source: ascii and binary are one path from here

import/import_analysis            finishNodes, analyzeImport: what the file is, before any choice
import/prepare_import             prepareImport: repairs, charts, library space
import/import_store               storeImport: a publication; loadImportedScene: into a SourceScene
```

`ImportedMesh` is what the file holds, in file units and file axes, with the names of its parts as
`nodes` over one mesh (`src/import/ARCHITECTURE.md`). Nothing is repaired while reading.
A fatal problem returns false with a message that says what is wrong, where, and what to do.

A GLB takes a shorter road: `loadGlb` opens the file in a scope that frees it, walks the nodes
(`readNode`, `readPrimitive`), then classifies, scales and welds what it copied.

## How a folder is published

Writers never open their final paths. `Publisher` makes the whole folder appear at once:

```text
begin(target)            may the target be replaced? (PublishPolicy)
                         a staging folder beside it, <target>.staging-<16 hex>, with an owner token
  the writer fills       pub.staging()
commit(key)
  writeMarker            key, then size and checksum of every file: written LAST
  swapIntoTarget         checked again; the old folder gets a token, moves to .old-<16 hex>,
                         staging becomes the target, the old one is deleted
abort()                  deletes staging; the earlier output stays byte-identical
```

`PublishPolicy::IntoFolder` is the export into a folder the user chose. That folder is never replaced:
staging goes inside it, and `commitIntoFolder` moves each file to its place, the marker last. A file
of the same name goes to `<staging>/.replaced` first, so a failed move puts everything back.

`recoverPublished` runs at startup. It acts only on siblings with the exact name pattern and this
target's owner token, so a user's `export.old-x` folder is left alone. A cache hit is a target whose
marker holds the requested key and whose checksums verify; staging folders never are.

## Checks

A file is refused rather than guessed at. Refusals are ordinary returns (`bool` with a message), never
exceptions. After a refusal the record is not to be used: readers clear it first and fill it as they go.

- **Limits before allocation.** `ImportLimits` bounds a foreign file: size before opening, counts from
  a PLY header before anything is reserved, per statement while reading OBJ.
- **Counts against the rest.** `ByteReader::count` and `string32`.
- **Whole records.** The end marker, or an exact end position, or `.kmesh`'s check of indices and array sizes.
- **Stale keys.** `staleKeyField` names the first input that differs. A knit is never reused on a guess.
- **Checksums.** Publication checks size and FNV-1a of every file; `tiled-knit` checks its own record.

## What must stay true

| Rule | Held by |
|---|---|
| `src/io` includes only `model` and itself. | the compiler (`CMakeLists.txt`), `architecture_module_dependencies` |
| No stored format changes by accident: any change of a written byte is a failure, never a re-pin. | `io_pin_knit_formats`, `io_pin_prepared_import`, `io_pin_scene_json` |
| Stage keys and import identities hash exact text: field order, `%.17g`, the options' lines. | `io_pin_stage_keys`, `io_pin_import_identity`, `io_pin_settings_text` |
| What a reader reports and refuses stays the same for the same file. | `io_pin_mesh_read` |
| The two FNV bases stay apart: `fnv1a64` for digests, `fnv1a` and `fnv1aFile` for assets and stage keys. | the pins above, `model/hash.h` |
| A corrupt count never allocates. | `io_byte_reader_refuses` |
| `validationSnapshot` cuts each piece at 511 characters: stored digests hold that cut. | `io_pin_validation_snapshot` |
| Publication never touches what it did not create, and survives a crash at any step. | `knit_tests_publish` |
| A read past the end of a `cgltf` file cannot happen: it is freed by scope. | `glb_required_extension_is_named` |

## Where to change what

| To change | Touch |
|---|---|
| A field of `.knitcache` | `writeKnitCache` and `readKnitCache` together, the same line in each, and a new magic: old files are refused. `kKnitCacheVersion` is different: raise it when knitting's output changes. |
| A shared curve field | `writeCurves` and `readCurves`. It reaches `knit.bin` too: raise `kSavedKnitVersion` in `saved_knit.h`. |
| A new binary format | `ByteWriter` and `ByteReader`, in the module that uses it. Magic first and last, `replaceFile` or a `Publisher`, a pin in `tests/tests_io_pins.cpp`. |
| A new foreign mesh format | `mesh_read.cpp`: a reader that fills `ImportedMesh`, and the dispatch in `readMeshFile`. |
| A new image use | `image_io`. Pixels are row-major with row 0 at the top. |
| What may replace a folder | `replaceable` in `publish.cpp`, with a case in `tests/tests_publish.cpp`. |
| A format's bytes on purpose | Bump its version, then re-pin. Prove the re-pin by undoing the change and hashing back to the old value. |

## Literature

- Cockburn, *Hexagonal architecture*: dependencies point at the domain, never out of it.
- Bernhardt, *Functional core, imperative shell*: codecs as pure functions, the disk at the edge.
- King, *Parse, don't validate*: a reader returns a complete record or a refusal.
- Kleppmann, *Designing Data-Intensive Applications*, ch. 4: one writer version, listed readers.
- Pillai et al. 2014, *All File Systems Are Not Created Equal*: why write-then-rename, and its limits.
