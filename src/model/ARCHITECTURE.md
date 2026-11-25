# Model: file tree and architecture

What a knit is: vectors, meshes, yarn curves and stitch records, what a knit asks for, plus the threads,
cancellation and identities every module shares. This folder is the bottom of the program. It includes
nothing outside itself, reads no files and knows nothing about parts or scenes. All lengths are
millimetres.

## File tree

The folder is flat. The four groups are its layers (below).

```text
src/model/
  ARCHITECTURE.md       this file
  CMakeLists.txt        the component knit_model: every header public, no dependencies

  foundation: no knitting in here
    vec.h               Vec2, Vec3, kPi, componentMin/Max: the values every file needs
    coordinates.h       where a point is: SceneMapping (library mm to scene units and the Y-up
                        export), UpAxis and toFileCoordinates (back to a source file's units and axes)
    colour.h            the sRGB transfer, the default yarn colour
    status.h            Status: Pass, Fail, Unresolved
    hash.h              FNV-1a for cache identities; splitmix64, hash01, hashNormal for repeatable draws
    noise.h, .cpp       seeded 3D gradient noise and its fractal sum in [-1, 1] (the Fur workflow's patches)
    env.h               reading environment variables
    log.h, .cpp         appendFormat, msSince, secondsSince: the stages' log text and timings
    compute_caps.h, .cpp
                        ComputeCaps: what this machine computes with (a CUDA device, the solver's
                        libraries), detected by gpu and read by settings, jobs and the viewer; and
                        buildConfiguration, what bit-identity is scoped to (cache keys, fibre digests)
    cancel.h            CancelToken, ThreadCancelScope, cancelCheckpoint
    parallel.h, .cpp    ExecPolicy, parallelFor, the process-wide CPU budget

  vocabulary: what a knit asks for
    params.h            StitchPattern, PatternParams, PatternReference (a flat swatch's verdict),
                        FailCategory and failCategoryText
    choices.h           HandKnit, FibreIslands, PlyProfileModel: choices a settings document holds
                        and knitting or fibres carry out

  curves: a knit as yarns and stitches
    bspline.h, .cpp     QuadBezier, CubicBezier (+ BezierHalves), BSpline (+ SpanPosition)
    topology.h, .cpp    StitchRecord, TerminalRecord, KnitTopology
    yarn_path.h, .cpp   StitchSpan, YarnPath (+ ArcLengthBounds), Polyline, samplePath, CurveSet

  surfaces: triangle meshes and how to search them
    spatial_hash.h, .cpp  SpatialHash3: candidate boxes, never a decision
    mesh.h, .cpp          TriMesh, closestOnTriangle, closestSegmentPoints, meshDigest
    mesh_query.h, .cpp    MeshQuery: closest point and inside test on a closed mesh
    mesh_check.h, .cpp    weldVertices, checkMesh, triangleIntersection: is a mesh fit to knit on?

  markings: what a user draws on an imported surface for the stitch mesher (the Mark step)
    surface_marking.h, .cpp  MarkingMeaning, SurfacePoint, SurfaceMarking, MarkingSet: control points on the
                        file's triangles, meanings, proposals and their setups; names, shape rules, construction roles
    marking_review.h    ReviewDepth, MarkingFinding (with its cause, row and places), LayoutMeasures (each measure Known or
                        not), MarkingReview; the proposals as data: PieceShape, ProposalSetup, ProposalSet
    surface_path.h, .cpp  TriangleNeighbours, TracePoint, MarkingTrace: the path between control points (a plane
                        section walked over the triangles, edge path as the fallback), closestSurfacePoint; a marking on
                        one border loop (an opening) and its trace along the rim
    surface_pieces.h, .cpp  SurfacePiece(s): an imported file's connected pieces per component, welded; file to
                        piece triangles; pieceWeldTolerance, surfaceDigest; proposal setups: basis, role conflicts,
                        acceptance, the preview's ephemeral set, merging
    stitch_mesh_data.h, .cpp  StitchMeshData (a built stitch mesh's polygons in file units, per-face kind, row and
                        findings), StitchMeshQuality (what its checks measured), stitchMeshUse (the only draw and save
                        gate), stitchMeshFindings and stitchMeshHeadline (generation and quality in words)
```

**Reading order:**

| Step | Header | What it answers |
|---|---|---|
| 1 | `vec.h` | What are the basic values, and in which unit? Millimetres, in double. |
| 2 | `coordinates.h` | Where is a point: in the library, in the scene, in an exported file, in a source file? |
| 3 | `bspline.h` | What is a yarn's shape? A uniform cubic B-spline, every span an exact cubic Bezier. |
| 4 | `yarn_path.h` | What is one yarn? A curve with arc-length tables and enclosures. And a `CurveSet`: the yarns of one part. |
| 5 | `topology.h` | Which stitch is drawn through which? |
| 6 | `params.h`, `choices.h` | What does a knit ask for, and which choices of the settings do knitting and fibres carry out? |
| 7 | `mesh.h`, `mesh_query.h`, `mesh_check.h` | What is a surface, what is closest to a point on it, and is it fit to knit on? |
| 8 | `status.h`, `cancel.h`, `parallel.h` | How does a long check report, stop and share the CPU? |
| 9 | `hash.h`, `log.h`, `env.h`, `colour.h`, `compute_caps.h` | Identities and draws, log text, environment, colour, the machine. |
| 10 | `surface_marking.h`, `surface_path.h`, `surface_pieces.h`, `marking_review.h` | What does a user tell the stitch mesher on a surface, and what does the mesher answer? |

## Layers

What a file may include from this folder (`architecture_model_layers` in `tests/tests_architecture.cpp` refuses
anything else; it does not require an edge to be used, so removing an include never fails it):

| Layer | Files | May include `model/` headers from |
|---|---|---|
| foundation | `vec`, `coordinates`, `colour`, `status`, `hash`, `env`, `log`, `compute_caps`, `cancel`, `parallel`, `noise` | foundation |
| vocabulary | `params`, `choices` | foundation, vocabulary |
| curves | `bspline`, `topology`, `yarn_path` | foundation, curves |
| surfaces | `spatial_hash`, `mesh`, `mesh_query`, `mesh_check` | foundation, surfaces |
| markings | `surface_marking`, `marking_review`, `surface_path`, `surface_pieces` | foundation, surfaces, markings |

Curves and surfaces never include each other, and nothing in the folder includes the vocabulary.

The edges the code has today (an arrow means "includes"):

```text
foundation   coordinates.h, colour.h, hash.h -> vec.h          parallel.cpp -> cancel.h
             noise.h -> vec.h                                   noise.cpp -> hash.h
             vec.h, status.h, cancel.h, env.h, log.h, compute_caps.h, parallel.h: nothing from here
vocabulary   params.h -> status.h                               choices.h: nothing from here
curves       bspline.h, topology.h -> vec.h                     yarn_path.h -> bspline.h, topology.h, parallel.h
surfaces     spatial_hash.h, mesh.h -> vec.h                    mesh_query.h -> mesh.h, spatial_hash.h
             mesh_check.h -> mesh.h, status.h                   mesh.cpp -> hash.h
             mesh_check.cpp -> spatial_hash.h
```

`parallel.cpp` includes `cancel.h` because a range's chunks keep the caller's token on any thread.
`yarn_path.h` includes `parallel.h` because `finalize` and `samplePath` take an `ExecPolicy`.

**What an edit costs.** `vec.h` reaches about 180 translation units: every file needs a `Vec3`. So what only
some files need lives in a header of its own: the coordinate conventions (`coordinates.h`, used by import,
export, settings and the writers, not by the solver), the settings' choices (`choices.h`), and the search
over a mesh (`mesh_query.h`, so a `TriMesh` user does not see the spatial hash).

## What belongs here

`model` is the bottom level: it includes no other component, reads no files, has no GPU qualifiers, and owns
no scene orchestration, file format or device operation (`architecture_model_is_plain`). Something lives here
when one of these holds:

1. **Shared values:** at least two components must agree on it (`Vec3`, `TriMesh`, `CurveSet`,
   `StitchRecord`, `Status`).
2. **Operations that describe such a value**, even with one caller: a member, or a free function whose result
   is a property of the value (the closest point on a triangle, a mesh's digest, an arc-length enclosure, a
   mesh's validity). One caller does not move a basic operation out; a generator never comes in
   (`generateStockinette` is knitting's).
3. **Settings vocabulary:** a record that the settings document holds and a component below settings carries
   out (`HandKnit`, `FibreIslands`, `PlyProfileModel`, `StitchPattern`). It cannot live with the component
   that carries it out: settings may not see knitting, and the viewer, which reads settings' headers, may see
   neither knitting nor fibres. Settings depending on knitting and fibres instead would point a stable
   component, which scene, export, jobs and the viewer depend on, at changing ones. Each such type's comment
   names its implementer.
4. **Runtime records and infrastructure every component uses, depending on nothing:** `status`, `cancel`,
   `parallel`, `log`, `env`, `hash`, `colour`, `ComputeCaps`. `ComputeCaps` is detected by gpu and read by
   settings, jobs and the viewer; the only components below both gpu and settings are model and fibres, and
   it has nothing to do with fibres.

A type that only one component uses, and that is none of 2 to 4, lives in that component: `QuadMesh` is
import's. Kept by this rule although one component calls them: `meshDigest` (a mesh's identity, bit for bit,
read by the scene's stored knits), `hashNormal` (the normal draw beside `hash01`, which knitting and fibres
use), `SceneMapping::toExportF`, `triangleIntersection` (inside `checkMesh`; public for its test).

## The data, as it is held

```text
CurveSet                       the knit of one part
  yarns: YarnPath[]            one continuous strand each
    curve: BSpline             control points; u in [0, spans], knots at the integers
    tables: Tables             u to arc length s, and bounds on s (made by finalize)
    stitches: StitchSpan[]     each stitch's range of u on this yarn, and its id in topo
  topo: KnitTopology           the graph, shared by all yarns
    stitches: StitchRecord[]   one loop: its yarn, range of u, parents, children, face, local frame
    terminals: TerminalRecord[]  cast-on rings and closing drawstrings
```

A yarn has two parameters, kept apart: `u`, the spline parameter and material coordinate, and `s`, the
arc length in mm. Arc length is an estimate; where membership in a window matters, `lengthBounds`
encloses the true length (`ArcLengthBounds`: chords below, control polygons above). A `Polyline` is a yarn
sampled with a proven chord-error bound.

A yarn's `Tables` are valid from `finalize` on; they describe the curve as it was then. Moving the controls
calls for `finalize` again (or for `setTables`, when the GPU has computed them). `stitchesSorted` and
`sortedCount` remember whether the spans were sorted and disjoint at that time, so `stitchContaining` may
search them; spans edited later fall back to the linear scan.

A stitch's range of u is held twice, on purpose: `StitchRecord` in the graph, `StitchSpan` in the yarn's own
list, linked by id (owner decision 2026-11-06). The knit cache stores both, and `curvesDigest` hashes the spans.

A `CurveSet` always carries its topology. Layouts fill it as they place stitches; a set built by hand
calls `fillFlatTopology`.

A `TriMesh` shares positions per vertex and stores uv per triangle corner, so a seam is an explicit
discontinuity. A `MeshQuery` borrows its mesh, which must outlive it, and never changes after its
constructor, so concurrent queries are safe. `SpatialHash3` only lists candidates, and the caller still
runs its own exact test.

## Where this folder sits

Includes run one way: `scene` includes `knitting`, which includes `model`, never the reverse. The
compiler holds it: the library `knit_model` declares no dependency, and its include path holds only this
folder's headers (`cmake/components.cmake`), so an include of another folder does not compile.

| Module | The question it answers | Its unit |
|---|---|---|
| `model` | What is a knit? | one yarn set, one mesh |
| `knitting` | How is it made, and is it valid? | one part |
| `scene` | What does the whole object hold, and in which order do the stages run? | all parts |

The values are written by several modules. This folder only defines them.

| Value | Who makes it | Who reads it |
|---|---|---|
| `TriMesh` | `io` and `import` read and prepare it | everyone below the viewer |
| `CurveSet`, `KnitTopology` | `knitting` places the stitches; `gpu` rewrites controls and tables and builds plies; `scene` restores it from a `.knitcache` | knitting, fibres, gpu, scene, export |
| `Polyline` | `knitting` and `gpu` sample it | knitting's checks |
| `PatternParams` | settings (`settings_mapping`) and knitting's fixtures | knitting, fibres, gpu, scene, export |
| `PatternReference` | knitting (`flatPatternReference`) | jobs, the viewer |
| `HandKnit`, `FibreIslands`, `PlyProfileModel` | the settings document | knitting, fibres; shown by the viewer |
| `ComputeCaps` | gpu (`detectComputeCaps`) | settings, jobs, the viewer |

## The common recipe

Renderers and geometry tools organise their bottom layer the same way:

1. Foundation below geometry below domain (USD `base/`).
2. One header per concern; what few files need is not in what every file includes (USD `gf`, `vt`, `work`).
3. Data are plain structs; algorithms are free functions or live elsewhere (libigl).
4. Accelerators are apart from the geometry they search (Embree):
   `mesh_query.h` apart from `mesh.h`.
5. File readers and generators are not in the data layer.
6. A set of curves is arrays plus offsets (`StrandSetF`).
7. Knit structure is a graph of stitches, separate from its curves (AutoKnit, Stitch Meshes).

`CurveSet` differs from rule 6, see below.

## What this project adds

| Need | Where |
|---|---|
| Verdicts that can say "not decided"; `Unresolved` never counts as a pass | `status.h` |
| Stopping any long loop from another thread | `cancelCheckpoint` in `cancel.h` |
| One CPU budget over nested parallel ranges, no deadlock | `parallelFor` in `parallel.h` |
| Results that do not depend on the number of threads | `ExecPolicy`, and sums made serially after parallel spans |
| Arc lengths with rigorous enclosures | `YarnPath::Tables`, `lengthBounds` |
| Polylines with a proven chord error | `samplePath` |
| Stable identities over saved work | the two FNV bases in `hash.h`, `meshDigest`, `Vec3`'s raw bytes |
| Stitch patterns on a lattice | `stitchFaceIsPurl` in `params.h` |
| Choices the viewer shows and knitting or fibres carry out | `choices.h` |

## CPU and GPU

No kernel reads a class from this folder: a `std::vector` cannot go to the device. `src/gpu` flattens
these values into the `gpu::*Arrays` records, one array per field with offset tables, and those arrays are
what the kernels see. So the layout here cannot change kernel speed, only the host copy in between. The
first pass recorded that copy as small (0.26 % of validation time; part of the 7.5 % fibres spend preparing
on the host); those are historical numbers, to be measured again before they justify a change. This is why
`CurveSet` stays a list of yarns and not one flat array.

Coordinates stay double here. The GPU works in float where it can; the NCP solver and the enclosures
need double. The GPU keeps its own copies of a few routines here (segment distance, closest point on a
triangle, arc-length bounds, `hash01`), each citing the function it follows.

## Do not change these

| What | Why |
|---|---|
| `Vec3` is three doubles | `curvesDigest`, `meshDigest` and the knit cache hash its raw bytes; a `static_assert` holds it |
| Both FNV starting values in `hash.h` | stored keys depend on each; merging them changes every key |
| The order of `.knitcache` fields | saved knits must keep loading; a retired field keeps its slot |
| `Tables::stitchesSorted` | `.knitcache` stores it and the GPU's arc-length tables write it; spans that are unsorted, empty or edited after `finalize` still need the linear scan |
| `CubicBezier::split` next to `sub` | `sub` gives a result some ulps away, which would change output |

## Further reading

- [src/fibres/ARCHITECTURE.md](../fibres/ARCHITECTURE.md): the next layer of the same idea, for fibres.
- Kaldor, Marschner, James 2008, *Simulating Knitted Cloth at the Yarn Level*: a yarn as a cubic B-spline.
- Yuksel, Kaldor, James, Marschner 2012, *Stitch Meshes for Modeling Knitted Clothing with Yarn-level Detail*.
- Narayanan et al. 2018, *Automatic Machine Knitting of 3D Meshes* (AutoKnit): a knit as a stitch graph.
- Baerentzen, Aanaes 2005, *Signed Distance Computation Using the Angle Weighted Pseudonormal*: `MeshQuery`.
- Ericson 2005, *Real-Time Collision Detection*, 5.1.5 and 5.1.9: `closestOnTriangle`, `closestSegmentPoints`.
- Lakos 1996, *Large-Scale C++ Software Design*: levelization, and what a header costs its includers.
- Martin 2017, *Clean Architecture*, part IV: the stable-dependencies and common-closure principles.
