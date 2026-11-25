# src: the components and their dependency graph

Each folder of `src/` is one component: a library with a declared interface and declared dependencies.
A change inside a component recompiles only that component. A change to one of its public headers
recompiles only the components that include that header. This file is the way into `src/`: what each
component is, the graph between them, and how the build makes sure the code follows it. Each folder's
own `ARCHITECTURE.md` describes the inside of that folder.

## The components

| folder | library | what it is |
|---|---|---|
| `model/` | `knit_model` | what a knit is: vectors, meshes, yarn curves, stitch records, threads and cancellation. Standard C++ only. |
| `knitting/` | `knit_knitting` | how a part is knitted and checked: layouts, checks, the yarn surface, settling, one subfolder per stage. |
| `stitchmesh/` | `knit_stitchmesh` | the stockinette stitch-mesh experiment before yarn: the reference surface, Method A's time field, rings and guided rows, Method B's W18 pipeline (RIM and HiGHS, private), the typed mesh with its independent check, the virtual shift, coverage, the construction fixtures, the Y12 mesh relaxation, the marking records of the imported arm, the first yarn embedding with its checks, and the adapter of that yarn to production curves. |
| `fibres/` | `knit_fibres` | fibres on knitted yarn, fur on meshes, the strand files the renderer reads. |
| `gpu/` | `knit_gpu` | the CUDA backend: validation, the knitting solver, surfaces, fibres. |
| `io/` | `knit_io` | files, the formats several components share, where the caches live, the Blender process of the `.blend` export. |
| `settings/` | `knit_settings` | what the user chose: the settings document, presets, the scene table, Draft and Final. |
| `import/` | `knit_import` | a user's mesh made knittable: analysed, repaired, rebuilt as a closed surface with new charts. |
| `scene/` | `knit_scene` | a scene's parts and its stages: Knit and Fibres, their caches and keys. |
| `export/` | `knit_export` | the reviewed results written as a package, a Blender file or an OBJ. |
| `jobs/` | `knit_jobs` | what the viewer and `knit_check` start: the four jobs, the queue that runs them, the reviewed results and the export rule, the display data. |
| `viewer/` | `knit_viewer_state`, `knit_viewer` | the window, the four workflow steps and the 3D view. |
| `cli/` | `knit_check`, `stitchmesh_lab` (with `stitchmesh_arm.cpp`, the imported-arm commands, `stitchmesh_smooth.cpp`, the smooth-reference commands, and `stitchmesh_yarn.cpp`, the yarn commands, and `stitchmesh_settle.cpp`, the yarn-settle commands) | the command line: programs, no library. |

## The dependency graph

A component depends on another in one of two ways:
- **public (P):** one of its public headers includes the other, so its users see the other component too;
- **private (p):** only its sources and private headers include the other, and its users do not see it.

| folder | public (P) | private (p) | headers public / private |
|---|---|---|---|
| `model` | none | none | 25 / 0 |
| `knitting` | model | none | 24 / 0 |
| `stitchmesh` | model, knitting | (RIM, pcg32, SLIM's Eigen, HiGHS) | 30 / 0 |
| `fibres` | model | knitting | 15 / 0 |
| `gpu` | model, knitting, fibres | (CUDA) | 4 / 13 |
| `io` | model | (cgltf, stb) | 12 / 0 |
| `settings` | model, io | fibres | 8 / 0 |
| `import` | model, io | (OpenVDB, SLIM's Eigen) | 4 / 3 |
| `scene` | model, knitting, fibres, io, settings | gpu, import, stitchmesh | 22 / 6 |
| `export` | io | model, fibres, settings, scene | 1 / 2 |
| `jobs` | model, settings, export, import, io | knitting, fibres, gpu, scene, stitchmesh | 11 / 2 |
| `viewer` | none | model, io, settings, export, jobs | 0 / 16 |
| `cli` | none | model, knitting, gpu, stitchmesh, import, fibres, io, settings, scene, jobs; export's private headers | none |

The same graph in levels. A component depends only on components of lower levels (in brackets: what it
uses besides model), so there is no cycle:

```text
level 7   viewer (io, settings, export, jobs)             cli (knitting, gpu, stitchmesh, import, fibres, io, settings, scene, export, jobs)
level 6   jobs (knitting, fibres, gpu, io, settings, import, scene, export, stitchmesh)
level 5   export (fibres, io, settings, scene)
level 4   scene (knitting, fibres, gpu, io, settings, import, stitchmesh)
level 3   gpu (knitting, fibres)     settings (io, fibres)
level 2   fibres (knitting)          import (io)          stitchmesh (knitting)
level 1   knitting                   io
level 0   model
```

For a change it means:
- **The viewer** sees jobs, settings, export, io and model only; the import's vocabulary reaches it through
  `jobs/import_session.h`, and it calls nothing of the import's store or pipeline. A change to a header of
  knitting, fibres, the GPU or the scene recompiles none of its files.
- **The export's format writers** see only `MaterialSpec` (`settings/material_spec.h`) of the settings
  and nothing of the scene (`architecture_export_door`). A new setting or a change in knitting or the
  scene does not recompile them.
- **A private header** is compiled only by its own component and that component's test program.

**Named exceptions**, each written in a `CMakeLists.txt`:
- `knit_check`'s strand demo builds a one-strand package by hand, so the program has the export's private
  headers (`PRIVATE_ACCESS export`, `cli/CMakeLists.txt`).
- `knit_tests_workflow` is the integration test program: it inspects jobs', scene's and export's internals
  (`knit_integration_test_access`, `tests/CMakeLists.txt`).

## How the build ensures it

**Each component declares itself.** Its folder's `CMakeLists.txt` calls `knit_component()`
(`cmake/components.cmake`) once, with the component's sources, its public and private headers, and its
public and private dependencies. Only the component folders have one; their subfolders do not. Settings
that belong to one component sit there too: the gpu's CUDA flags, import's Eigen sub-library.

**The compiler refuses what the graph does not allow.** `src/` is on no include path. For every header,
`knit_component()` writes a two-line forwarding header into the build folder:

```text
build/components/<folder>/public/<folder>/<path>.h    for each public header
build/components/<folder>/private/<folder>/<path>.h   for each private header
```

A forwarding header holds only `#include "<the real header's absolute path>"`. The library gets the
public root as a PUBLIC include directory and the private root as a PRIVATE one, and links its
dependencies PUBLIC or PRIVATE as declared; CMake hands the public roots of PUBLIC dependencies on to the
users. Includes keep their spelling (`#include "knitting/checks/checks.h"`), but resolve only where the
graph allows them. An include of an undeclared component, or of another component's private header,
stops the build:

```text
src\model\log.cpp(2): fatal error C1083: Cannot open include file: 'knitting/checks/checks.h'
```

In the editor, opening an include lands on the forwarding header first; the real file is one step on.
Edit the real file: the forwarding headers are written again at every configure.

**The configure step checks the declarations** (`knit_components_finish`):
- every header of a folder is in exactly one of its two lists, and every listed header exists;
- a dependency is a component declared before (`src/CMakeLists.txt` adds them in the order of the
  layers), so the graph cannot have a cycle;
- the forwarding headers of removed headers and of removed components are deleted;
- no header shadows another component's header: a quoted include looks in the including file's folder
  first, so a file at `scene/fibres/fibres/x.h` would win over the forwarding header for `fibres/x.h`;
- it writes `build/components.json`, the graph and the header lists, for the checks and the tools below.

**What the compiler cannot see is tested** (`knit_tests_architecture`):
- `architecture_module_dependencies`: every component a file includes directly is a declared dependency
  of its folder (read from `components.json`). A PUBLIC dependency hands its own dependencies' headers
  on, so the compiler alone would accept them. The test also refuses `../` and absolute include
  spellings, which would get around the forwarding headers.
- `architecture_knitting_layers`, `architecture_scene_layers`, `architecture_gpu_layers`,
  `architecture_viewer_layers`, `architecture_import_layers`, `architecture_jobs_layers`: the order of the
  subfolders or files inside a component; the viewer's test also says which components each of its layers may
  include.
- `architecture_export_door`: the format writers stay away from the domain, and the viewer and the jobs
  name no format.
- The `*_doc_complete` cases: a folder's `ARCHITECTURE.md` names every file of the folder.

**Every public header compiles on its own.** The target `knit_header_checks` compiles one translation unit
per public header, as a user sees it. That catches a public header that misses an include, includes a
private header, or needs a dependency its component declared private.

## The test programs

Each test program belongs to one component and compiles like that component's sources: with its private
headers and dependencies (`knit_component_private_access`). A change to one component recompiles and
relinks only the test programs that depend on it.

| program | component |
|---|---|
| `knit_tests` | knitting |
| `knit_tests_fp` | fibres |
| `knit_tests_gpu` | gpu |
| `knit_tests_s6` | scene and settings |
| `knit_tests_export` | export |
| `knit_tests_import` | import |
| `knit_tests_ui` | viewer |
| `knit_tests_workflow` | the integration program (see the exceptions above) |
| `knit_tests_publish` | the publication protocol (io and the export's door) |
| `knit_tests_architecture` | the checks above |

A case that moves to another program keeps its CTest name (`knit_register_cases`). Helpers several
programs share are in `tests/strand_fixtures.h` and `tests/scene_fixtures.h`.

## Changing the graph

| To | Do |
|---|---|
| Add a source or header | list it in its folder's `CMakeLists.txt`, as public or private. Configure fails until it is listed. |
| Let other components use a header | move it to the public list. Everything its includes reach in this component must be public too; `knit_header_checks` says so if not. |
| Add a dependency | add it to the folder's PUBLIC or PRIVATE list: PUBLIC only when a public header includes it. Then update the table above and the folder's `ARCHITECTURE.md`. |
| Add a component | a folder with a `CMakeLists.txt` calling `knit_component()`, added to `src/CMakeLists.txt` above its users. Then this file and the folder's `ARCHITECTURE.md`. |
| A program, no library | `knit_program()`, as `cli/CMakeLists.txt` does. |

## Measuring it

`cmake/component_report.py <build folder>` reads `components.json` and the header dependencies Ninja
recorded, and prints which translation units each component's public and private headers reach, and
how many recompile when a given header changes; `--verify` compares that with a real build.
