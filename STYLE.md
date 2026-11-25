# Code style

The rules for code in this repository. Prefer clear, direct code and one canonical implementation.

## Educational clarity, for every new or rewritten line

Code is **readable by an engineering student** who knows C++, the geometry and the knitting terms but not this codebase. Opening one file tells them what each type represents, what each function does, and in which order the algorithm runs, without chasing indirection. Clarity beats cleverness; the only exception is a measured performance need, and it needs a recorded number (in the commit message).

- **Names state the domain meaning.** Types are nouns from the geometry or the product (`StitchProgram`, `CurveSet`, `KnitResult`, `FibreResult`); functions are verbs saying what they do (`placeSpiral`, `checkLocalSignature`, `applyStitchSize`). No names that describe nothing (`Manager`, `Helper`, `Util`, `Impl`, `Data`, `Info`, `Handler`, `process`, `doWork`). Abbreviate only standard domain terms (`uv`, `ncp`, `pcg`, `glb`, `ao`, `D` for the yarn diameter, `w`/`h` for wale and course spacing); single letters only inside a formula that cites its source.
- **A class interface is a short list of capabilities.** Each public method does the one thing its name says. No boolean parameter that switches behaviour: write two named functions or take an `enum class`. More than about five parameters become a named struct. Prefer return values to output parameters.
- **The algorithm reads top to bottom, in the order it happens.** Dispatch is a `switch` over an `enum class`, not a virtual hierarchy, function pointer, callback registry or template. Data flows through arguments and return values, not hidden state.
- **No templates in project code.** Do not write class, function or variable templates, specializations, generic lambdas (`auto` parameters), CRTP, SFINAE, concepts, `if constexpr` type switches, type erasure or metaprogramming. Write the concrete function for the concrete type: two short concrete functions are clearer than one generic one. Standard-library containers and vocabulary types are fine (`std::vector`, `std::string`, `std::optional`, `std::span`, `std::unique_ptr`, `std::shared_ptr` for published immutable results). Existing project templates may be called, never extended or copied.
- **No new macros** beyond the existing ones: the test harness (`TEST`, `CHECK`, `REQUIRE_GPU`, `SKIP`), the host/device qualifier `KNIT_CFIT_HD`, and the platform defines (`NOMINMAX`, `WIN32_LEAN_AND_MEAN`, the third-party implementation switches). A constant is `constexpr`; a function is a function.
- **Explicit over implicit.** `explicit` single-argument constructors; `auto` only when the type is spelled on the same line; no operator overloading beyond the math vector and matrix types (`Vec2`, `Vec3`, `Mat2`); allocation, ownership and synchronization visible at the call site.
- **Each file says what it is for** in its first comment lines. A formula cites its source (paper, equation). Comments explain why, never restate what.
- **Generalize only for a second real caller.** One concrete use is written concretely, with a comment naming the general concept when there is one.
- **Before finishing**, re-read the new code as that student: can each function be explained in one sentence from its name and body? If not, rename or split it.

This applies to code you write or rewrite; untouched code is not refactored only to comply. It is enforced, not just stated: `architecture_no_new_templates` and `architecture_no_new_generic_lambdas` (`tests/tests_architecture.cpp`, CTest `s6_architecture_*`) pin today's counts in `src/` and `third_party/` (not the `shim/` folders, which stand in for third-party template APIs), which may only fall. A change that adds one fails by name and lists the places; a change that removes one lowers the pinned number in the same commit. An unavoidable new template is a decision recorded with its reason and number, never a silent raise of the constant.

## Project invariants

- **Reviewed results flow downstream.** Fibres consume a `KnitResult`, the export job a `KnitResult` and a `FibreResult` (`src/jobs/results.h`); nothing knits again implicitly from settings keys (the NCP solver is not deterministic). A cancelled or failed job publishes nothing; the previous result stays.
- **The builder is GL-free.** `SceneBuilder` publishes one immutable `BuiltScene`, assigned under the queue's lock; only the viewer's main thread uploads to GL. Long loops poll `model/cancel.h` checkpoints.
- **One settings table.** Every scene setting is one row of `settingFields()` (`src/settings/scene_settings.cpp`) with its stage; reader, writer, diff and UI bindings use it. Unknown keys are refused, so a removed key goes into an explicit retired-key list that is accepted and ignored on read: saved settings files stay loadable.
- **Verdicts are honest.** Checks return `Status::Pass`, `Fail` or `Unresolved` (`src/knitting/checks/checks.h`); `Unresolved` never counts as a pass, and a gate is never weakened to make a run pass.
- **One place per output format.** The scene.json material entry and its schema version are written in one function (`src/export/format_ppt.cpp`); exported lengths are true metres; the scene graph is root > objects > fibre groups + Ground. Persisted-field changes cover writer and reader in the same commit.
- **Caches are keyed by identity.** Scene, knit and ply-map caches derive their key from the settings and inputs that determine the content; prepared imports (`out/cache/imports`) are referenced from saved settings and are never treated as disposable.
- **Tests report honestly.** A case that cannot run for a missing fixture calls `SKIP(reason)` (exit 77 when every case run skipped, CTest `SKIP_RETURN_CODE`); it never returns early as a pass.

## Architecture, when changing boundaries or ownership

- Base redesigns on a concrete defect: duplicate ownership/logic, an illegal dependency, an unnecessary API, or measured cost. Metrics diagnose problems; they are not quotas.
- Keep one owner per fact and cohesive responsibilities. Prefer deleting duplication over adding layers; avoid speculative abstractions and extension points.
- Public APIs expose capabilities needed by real callers, with explicit ownership and lifetime. Move misplaced responsibilities to their owner rather than merely hiding methods.
- Enforce changed invariants with the earliest practical mechanism: compiler, build/link boundary, then tests. Verify a new check detects its intended failure. Documentation alone is not enforcement.
- Resolve ownership/type/boundary decisions before propagating compiler fixes. Verify a cache, synchronization mechanism, or compatibility path is unnecessary before deleting it.
- Split a file only when that materially improves understanding; merge files only on demonstrated redundancy (a forwarding wrapper, a tiny header/source pair with one caller). Third-party implementation units (`cgltf_impl.cpp`, `stb_impl.cpp`) keep their boundary.

## C++/CUDA conventions

- Names: `camelCase` functions and locals, trailing `_` private members, `PascalCase` types and enum values, `kPascalCase` constants, `UPPER_CASE` macros. Namespaces `knit`, `knit::gpu`, `viewer`. Use positive booleans and standard domain abbreviations.
- K&R braces, four spaces, roughly 120 columns (existing longer lines are not re-wrapped for that alone). Unbraced bodies only for a single statement on the same line. Use single blank lines and no decorative banners; section labels may use `// --- name`.
- Includes: own header, standard library, project, separated by blank lines. Prefer `#pragma once`. Header implementations are limited to inline/constexpr/trivial functions.
- Prefer guard clauses and flat control flow. Extract at natural responsibility boundaries, not line-count targets. Structs hold simple public data; classes enforce invariants. Prefer composition.
- Comments explain non-obvious math, thresholds, tradeoffs, ownership, invariants, or traps. Keep them short; no restating code, density quotas, or migration history in headers. Use lowercase `//` comments except proper names/acronyms.
- Prefer `constexpr`, `enum class`, const inputs, const non-mutating methods, and `[[nodiscard]]` for failure returns. Use `auto` only when the type is spelled on the same line.
- Use explicit RAII resource wrappers, move-only ownership, and non-owning views (`std::span<const T>` where supported). Keep allocation and synchronization visible. GPU code keeps device buffers in owning wrappers and checks every CUDA call at the API boundary with the function `check()`; explain non-obvious launch dimensions. Shared host/device helpers use `__host__ __device__` (or `KNIT_CFIT_HD` in `compression_fit.h`).
- Prefer failure returns (`bool` with an `std::string* err` message, or a `Status`) over exceptions for control flow; the builder's job boundary reports what a stage throws. Below the GPU's public headers (`gpu/device.h`, `knitting.h`, `surfaces.h`, `fibre_generator.h`), a device failure is a `DeviceError`, because it is not control flow; those entry points turn it into `bool` and `err`. Warnings and errors go to stderr. Avoid `std::endl`, `using namespace std`, and hot-path strings/heap churn.
- Prefer disjoint thread ownership (`ExecPolicy` chunks); explain atomic ordering. Parallel and serial runs give identical results where a test says so. Profile before significant performance optimization and verify its effect.

## Viewer draw paths, only when touched

For changes to `drawSceneUi`, the panels in `src/viewer/panel_*.cpp`, the per-frame loop in `src/viewer/main.cpp`, or the GL uploads in `src/viewer/renderer.cpp`, apply the viewer draw-path checks below.

### Viewer draw-path checks

The viewer redraws continuously. A loop over scene-scale data (stitches, yarn control points,
fibres, triangles) inside `drawSceneUi` or the per-frame loop of `src/viewer/main.cpp` runs 60
times a second for a panel nobody is looking at. Apply these checks to changed draw functions and
their dependencies.

**Ownership first:**

- [ ] Is the value a fact about the built scene (counts, verdicts, part lists, bounds)? Then the
      builder computes it once when it publishes the `BuiltScene` (`src/jobs/scene_builder.h`),
      and the panel reads it. The panel never derives it from geometry per frame.

**Per-frame cost:**

- [ ] Loops over `DisplayMesh` vertices, curve control points or fibre strands? Only in the GL
      upload when a new `BuiltScene` is shown, never per frame.
- [ ] Builds a `set`/`map`/`vector`/combo string from scene state per frame? Compute it when the
      shown scene or the settings revision changes, and keep it next to that state.
- [ ] Does a `find` or linear search run inside a loop? Hoist a lookup map.
- [ ] Is the scan above a `CollapsingHeader` or tab `if`? A collapsed section must not cost
      scene-scale work.
- [ ] Is a GL object (buffer, texture) cached in a map keyed by a pointer? It needs explicit
      release when its scene is no longer shown (`Renderer::releaseUnused`: what is on the GPU
      belongs to a scene in `MeshOwners`); a stale key frees nothing.

**Two traps:**

- [ ] A request to the builder takes an immutable settings snapshot. The panel never reads a
      half-edited settings object from another thread, and the builder never reads UI state.
- [ ] Partial scenes (`BuiltScene::partial`) are previews, never results: a panel that offers
      Fibres or Export must check that the shown scene has its result object.
