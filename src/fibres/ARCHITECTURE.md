# Fibres: file tree and architecture

Short fibres on knitted yarn, fur on meshes, and the strand files the renderer reads. This file is the
way into the folder: what each file owns, how data moves through them, and what must stay true when you
change them. Scenes grow their wool fibres on the GPU, with kernels in the gpu component; the generator
here is the same rules written to be read, and the specification the kernels follow.

## File tree

```text
src/fibres/
  ARCHITECTURE.md       this file
  CMakeLists.txt        the component knit_fibres (src/ARCHITECTURE.md): every header public; depends on
                        model, and privately on knitting (the headers name the yarn surface through
                        pointers only)

  Wool fibres on knitted yarn
  fibre_record.h        Fibre, Vec3f, IslandWeight: plain data compiled by host and GPU
  fibre.h, .cpp         the vocabulary: FibreParams and the preset, FibreOutcome, FibreReport
                        (addCandidate), FibreSet, FibreFrame, setupPieceFibres, fibreParamsText
  fibre_shape.h, .cpp   one fibre as maths: draws, DirectedCurve, RandomCurve, SmoothCurve. No tubes
  fibre_generator.h, .cpp
                        the rules, one candidate at a time: candidateSource, growYarnFibres
  fibre_islands.h, .cpp patches of different fuzz: weights per control of the strands; the island blobs and
                        their weight at a point (addIsland, islandWeightAt), shared with the fur
  fibre_budget.h, .cpp  BudgetSchedule: which key windows to grow until the segment budget's threshold
                        is known
  fibre_checks.h, .cpp  everything that inspects finished fibres; none of it changes them

  Strands for the renderer
  strands_file.h, .cpp  the .pptstrands format: StrandSetF, the renderer's rules (StrandRules,
                        checkStrands), reading and writing, keepStrands; the margins generators keep on
                        those rules (kMinSegmentRadii, kTaperStep, taperedRadius)
  fibre_strands.h, .cpp a FibreSet as strands: fibresToStrands, and the late taper (StrandTaper,
                        fibres.tail_taper*) with its report

  Fur on meshes (the Fur workflow, no knitting)
  fur_look.h            FurLook: the pile's loops, their variation, the fine fibres and the settle in the
                        generator's units; FurLayerCounts, FurReport
  fur_field.h, .cpp     the cached candidate roots of one part and layer (sampleFurRoots, FurRootCache) and one
                        disc per fibre (buildFurField, FurDisc): the control and preview of each fibre; the
                        per-part islands, patches, clumps and wear (buildFurModifiers) and the density bound
  fur_loops.h, .cpp     a disc's rest-shape curve (loop or straight, equal chords, the renderer's joint rule for
                        loops) and growLoops: settle and stuck test, counted per layer
  fur.h, .cpp           the static settle (settleHair, Follow-The-Leader) and the stuck test (hairInsideBody)

  The yarn body of structural fibres (a prototype)
  yarn_body.h, .cpp     1 to 3 plies of long fibres in a shell on the yarn's envelope, built on the CPU:
                        buildYarnBody, measureOverlap. Behind yarn_body.enabled (off by default; the Knit
                        step's Construction "Made of fibres", the blue swatch preset). It runs on the CPU,
                        not a second fibre path: a device port is the owner's decision
  yarn_bake.h, .cpp     the body's fibres baked into per-yarn maps on the tube (albedo, normal, specular,
                        roughness, height): radial rays on the (s, theta) chart in nominal space, the
                        atlas layout. A pilot behind yarn_body.bake (the yarn body's design note, section 6)
```

**Reading order:** `fibre_record.h`, `fibre.h`, `fibre_shape.h`, `fibre_generator.cpp`, then `fibre_budget.h`,
`fibre_strands.h`, `fibre_checks.h`, then for fur `fur_look.h`, `fur_field.h`, `fur_loops.h`, `fur.h`. The kernels' twin is `src/gpu/fibres/fibre_candidate.cuh` and
`fibres.cu`; the scene's use of all of this is `src/scene/fibres/` (`fibre_stage.cpp`, `part_fibre_job.cpp`).

## The common recipe

Procedural yarn, real-time cloth and groom tools share one shape:

1. **Enumerate identities.** A fibre is (seed, strand, index), so it needs no stored state.
2. **Draw from a hash.** Every random value is a pure function of the identity.
3. **Evaluate a curve in local coordinates.** Here: Eq. 11 of Zhao et al. 2016, in yarn coordinates.
4. **Map it to the world.** A fibre is a polyline of equal-length chords.
5. **Thin by rank.** A lower key threshold keeps fewer fibres and never reshuffles the rest.

That is enough for a renderer that invents fibres on the fly. We export real geometry next to real
tubes, so each fibre also has to fit there: a buried root and an emergence point, exposure and cuts
against the tubes, the renderer's strand rules, an exact budget, and the reason every refused candidate
was refused. Most of the stages under *One candidate* below exist for that reason.

## Layers

An arrow means "includes". The includes run one way, so each file can be read knowing only the files
to its right. `architecture_fibres_layers` holds this, for the headers and the `.cpp` files.

```text
fibre_budget.h ──> fibre_generator.h ──┐
fibre_shape.h ─────────────────────────┤
fibre_checks.h ────────────────────────┼──> fibre.h ──> fibre_record.h <── fibre_islands.h
fibre_strands.h ───────────────────────┘        │
        └──────────────────────────────────> strands_file.h <── fur.cpp

.cpp only: fibre.cpp ──> fibre_islands.h, strands_file.h
           fibre_generator.cpp ──> fibre_islands.h, fibre_shape.h, strands_file.h
leaves:    fibre_record.h, strands_file.h, fur.h, yarn_body.h include nothing of this folder
fur:       fur_loops.h ──> fur_field.h ──> fur_look.h ──> fur.h;   fur_field.h ──> fibre_islands.h;   fur.cpp, fur_loops.cpp ──> strands_file.h
           yarn_bake.h ──> yarn_body.h
```

- `fibre_record.h` includes nothing at all: the CUDA compiler reads it.
- No header includes knitting. `fibre.h` names the yarn's surface (`SurfaceModel`, `YarnSurface`, `PlySet`)
  through pointers; the `.cpp` files include what they query (`fibre.cpp`, `fibre_generator.cpp`,
  `fibre_shape.cpp`, `fibre_checks.cpp`, `yarn_body.cpp`), so the users of this folder do not compile knitting.
- Fur does not include `fibre.h`. It works on a `TriMesh` and produces `Hair` values; the scene copies those
  into a `FibreSet`, so the export sees one kind of result (`src/scene/fibres/fur_stage.cpp`).
- This folder never includes `src/gpu`: `gpu/fibre_generator.h` includes `fibre.h`, and the kernels'
  `fibres.h` includes `fibre_record.h`. Only `src/model` and `src/knitting` are on this folder's include
  path (`CMakeLists.txt`), so the compiler holds it.

## The data

### A candidate

A candidate is a fibre that might exist, identified by (strand, index). Strand `y` holds
`FibreFrame::candidates[y]` or `floor(density L_y)` candidates. `candidateSource` lists them with no
geometry, and `CandidateSource::key` is the candidate's rank: hash channel 9, known before anything is
computed. Every draw is a hash of the identity too, so a candidate gives the same fibre whichever others
run with it.

Each candidate ends in one `FibreOutcome`: `Accepted`, or the reason it was refused.
`FibreReport::addCandidate` counts it, in the total and in its population; both generators count through
it.

### Records and the point pool

A `FibreSet` holds one `Fibre` per accepted fibre and one packed array of points. The record finds its
points by offset:

```text
fibres   [ Fibre 0 ][ Fibre 1 ][ Fibre 2 ] ...
             │ firstPoint = 0, pointCount = 4
             │          │ firstPoint = 4, pointCount = 3
             v          v
points   [ r  c  v  v ][ r  c  v ][ r  c  v  v  v  v ] ...

r  the buried root, inside the tube
c  the curve start; the fibre leaves the tube on the segment r-c, at Fibre::anchor
v  the visible vertices, up to the tip
```

The kernel writes this layout and the host keeps it, so nothing is converted in between. The records are
136 bytes; points are float32, as the kernels compute them. Most fibres have about 4 points and a few have
50, which is why the points are not stored in the record. The layout is the one `StrandSetF` writes to disk.

## How a scene grows its fibres

Scenes grow wool fibres on the GPU only; `growWoolFibres` in the scene's fibre stage drives it.

```text
growWoolFibres
│
├─ every part: PartFibreJob::prepare   setupPieceFibres: root depth, exposure length, seed, plies, islands
│                                      candidateSource: identities and keys, no geometry
│
├─ generateInBudgetOrder               fibres.generation = budget, the default (full: every candidate
│    │                                 in one call, generateEveryCandidate)
│    │  per object: candidatesPerBucket for every part, then a BudgetSchedule
│    │  while not finished: nextWindow ─ evaluateWindow ─ addWindow
│    │                                    │
│    │                                    └─ PartFibreJob::generate ─ gpu::generateFibres
│    │                                         placeCandidatesKernel   anchors: which survive
│    │                                         growFibresKernel        grow the survivors
│    │                                         packPointsKernel        records + one pool
│    │                                         mergeFibres             FibreReport::addCandidate
│    └─ finish: every part's fibres in candidate order, the threshold q
│
├─ every part: PartFibreJob::finish   fibre_checks.*, when asked for (--fiber-check)
│
└─ applyObjectBudgets                 the export keeps float(key) < q
```

The Fibres step's reference preview and the strand demo take a shorter road: `growYarnFibres` on the CPU,
with no budget.

### One candidate

`growCandidate` in `fibre_generator.cpp` is the whole story in about 35 lines. `placeAnchor` runs first
(where on the strand, which population; drawn before any test), then these tests in order; the first
refusal ends the candidate:

| Step | Question | Refused as |
|---|---|---|
| `facesOutward` | Does the anchor face out of the fabric? | `BackFacing` |
| `onOuterPlySide` | On plies: is it on the yarn's outer side? | `PlyInnerEarly` |
| `anchorExposure` | Is it outside every tube, and is the ray to the outside clear? | `OriginInside`, `Occluded` |
| `keptByThinning` | Does it survive the artistic thinning? | `Probability` |
| `drawFibre`, `chooseSize` | Length, radius, number of segments: too few? | `TooShortToBend` |
| `setUpCurve` | The curve; on plies in the parent yarn's coordinates (`toParentYarn`) | `PlyInnerLate` |
| the root's depth | Can the root be buried this deep inside the strand? | `RootDepth` |
| `fitPolyline` | The length's fit: equal chords, `placeRoot`, `findEmergence` | `YarnEnd` |
| what `fitPolyline` found | Root inside the strand's arc, first chord steep enough, an emergence? | `RootRange`, `NoEmergenceElev`, `NoEmergence`, `RootNotInside` |
| `meetsStrandRules` | Are the segments legal for the renderer (`kMinSegmentRadii`, `kTaperStep`)? | `RuleRejected` |
| `measureTurns` | Is every vertex's turn within the bounds? | `TurnTooSharp`, `NoBend` |
| `SmoothCurve::deviationFrom` | Is the polyline close to its curve? | `ApproxError` |
| `cutAtDeepContact` | Cut before anything goes deeper than 0.2 r into a strand, or reject. With `FibreParams::drape` (directed only, light checks only: `finishFull` has no such rule) a fibre first tries to rest on the strand (`drapeOverStrands`: lifted along the strand's normal, never longer than before, all checks again). | `CollisionShort`, `Collision` |
| `cutBehindFabric` | Cut before anything goes behind the fabric's back, or reject. | `BehindFabric` |
| `recordAccepted`, `toFibre` | Count it, write the record. | |

Everything up to `keptByThinning` happens at the anchor. On the GPU that is `placeCandidatesKernel`, where
most candidates end; it also refuses `HiddenAnchor` (inside another solid) and `OccludedSelf`, which the
CPU's preview does not test. The rest is `growFibre`: the same stages under the same names, in the same order,
finished by `finishLight` (these rules) or `finishFull` (`fibres.checks = full`, no CPU twin).

### The budget

The export holds at most `fibres.segment_budget` segments per object and keeps the fibres with the lowest
keys. The threshold `q` depends on how many segments the accepted fibres have, which is only known after
growing them. `BudgetSchedule` avoids growing all of them:

1. Count every candidate's key into buckets (`candidatesPerBucket`), with no geometry.
2. `nextWindow`: a pilot of about 2 % of the candidates, then the estimated remainder x 1.1, then at
   least doubling, always whole buckets.
3. `addWindow`: the window's accepted segments, summed in bucket order, until the sum would cross the
   budget.
4. At the crossing bucket its fibres are sorted and `budgetThreshold`'s scan finds `q` exactly.

The fibres below `q` are therefore the ones a full run would keep, and the export is byte-identical;
candidates far above `q` are never grown. The schedule is data: the scene evaluates its windows, holds
or releases each part's device support, polls cancellation and writes the log. One window is pending at a
time, and the schedule refuses a window added twice or a part that evaluated other candidates.

### On the GPU

`gpu::generateFibres` runs a key window in batches:

```text
placeCandidatesKernel   anchor tests per candidate; survivors flagged
Flagged                 compact the survivors' seeds
growFibresKernel        grow each survivor into slot t: its record, and points at slotPoints + t * stride
ExclusiveSum            each record's offset into the pool
packPointsKernel        copy the points into the pool, and set firstPoint
Flagged                 compact the accepted records
```

`stride = 2 + max(maxSegments, max(2, minSegments))`: room for the root, the curve start and the segments
of the run's own settings.

## Writing the strands

`fibresToStrands` writes each fibre's polyline as it was generated, with `vertexRadius` (`taperedRadius`
from the root to the tip). With `fibres.tail_taper > 0` the export uses `fibresToTaperedStrands` instead: a
fibre keeps its root radius along its body and narrows over its tail by a cosine. The taper belongs to
writing the fibre, not to the fibre: the generator, the cache and the digests see the untapered fibre.

Where the tail's radii are carried is a `TaperPoints` mode. `Dense`, the export's, adds a vertex wherever
the radius falls by 1.35, inside the renderer's 50 % step rule. The sparse modes (`Existing`, `Onset`,
`Refined`; `knit_check --taper-points`, a developer switch) keep the generated vertices and let one segment
step from r0 to tip x r0. They are checked with `taperRules`, never with the default rules, and
`measureTaper` grades them on the stored float32 strands. A fibre whose tapered strand would break a rule
keeps its generated strand.

`writeStrandsFile` refuses a set `checkStrands` refuses; `readStrandsFile` refuses what the renderer's
reader refuses.

## Checks

`fibre_checks.*` answers a different question from the generator: is the result right against the
geometry actually exported? It shares no rule code with the generator, so a bug in one does not hide in
the other.

- `ExportTubeIndex`: the exported polygonal tubes, indexed once and queried per fibre.
- `fibreReferenceCheck`, `rootAttachmentCheck`: `--fiber-check`.
- `ExportTubeIndex::accepts`: the root check the GPU leaves out. `knit_check --fiber-check` audits it and
  reports the count; it is not applied.
- `measureFibreAngles`: reported angles, never a bound (`--fiber-check`).
- `fibreSnapshot`, `fibreSnapshotChunkDigest`: exact text and digests, for `--fibre-digest`, saved fibres and
  the regression tests.

## Fur

Fur reads like a groom tool's modifier stack: `sampleFollicles` (area-weighted roots on a mesh, coverage,
hidden roots dropped), `groomHair` (the rest shape), then a static settle per hair (Dynamic
Follow-The-Leader with a shape constraint, and the body pushed out along its face normals). `growFur`
applies them to every follicle; a hair still inside the body is counted, never kept.

## The yarn body

`yarn_body.*` builds each yarn as plies of long fibres in a shell, bent onto the surface model's rings. It
runs behind `yarn_body.enabled` (the Knit step's Construction "Made of fibres"), on the CPU until the owner
decides on a device port; the export writes its fibres in place of the tubes, in strand files of at most
2.5 M segments. Its offsets are 64 bit (a body may pass 2^32 points) and each ply zone is a cancel point.

`yarn_bake.*` bakes that body into maps on the yarn's plain tube (a pilot behind `yarn_body.bake`): every
fibre point goes back to the body's nominal yarn by the inverse section map, radial rays on the tube's
(s, theta) chart find the first fibre, and the hits become albedo (dye under the height field's AO), a
tangent-space normal, the specular share and Toksvig's roughness, laid out as one atlas per yarn. The
export (`export/collect_scene.cpp`) meshes the tube with the matching uv and writes a `plastic` material
per yarn. README section 6 has the design and the measured outcome.

## What must stay true

| Rule | Held by |
|---|---|
| The CPU and GPU make the same decisions, in the same order. A rule changes in `fibre_generator.cpp` and the kernels together. | `gpu_fibres_properties_and_cpu_agreement` (0 flips) |
| Hash channels never change meaning. A changed channel moves every fibre. | snapshots in `tests/data` (`fibers_snapshots`) |
| `Fibre`'s layout matches on both compilers. | `static_assert` in `fibre_record.h` |
| `FibreOutcome` and the kernels' codes agree, value by value. | `static_assert` in the kernels' `fibres.cu` |
| The cache key text never changes by accident. `fibreParamsText` prints retired settings as constants, so old caches stay valid. | `fibre_params_text_pinned`, `sizeof(FibreParams)` in `fibre.cpp` |
| Values derived from a fibre are computed in double (`vertexRadius`, `fibreLength`). In float the export's bytes change. | `--fibre-digest`, the sha1 of the exported strand files |
| The budget keeps the full run's fibres below q, with the same windows as before it became data. | `budget_sched_*`, `budget_sched_windows_pinned`, `fibre_budget_matches_full_per_backend` |
| Fur stays deterministic. | `fur_strands_pinned` |
| No fibre goes deeper than 0.2 r into a strand. | `fibres_cut_before_deep_contact` |
| No fibre point from the curve start on lies behind the fabric's back: on a flat patch the lowest point of the tubes, on a part `FibreFrame::backDepth` behind the base surface. | `fibres_stay_in_front`, `gpu_fibres_properties_and_cpu_agreement` |
| The late taper: dense by default; the sparse modes are checked with `taperRules`. | the `late_taper_*` cases |
| A zero island weight changes nothing. | `fibre_islands_zero_weight_unchanged`, `fibre_islands_cover_what_they_ask` |
| The includes follow the layers above, and this file names every file of the folder. | `architecture_fibres_layers`, `architecture_fibres_doc_complete` |

## Where to change what

| To change | Touch |
|---|---|
| A fibre's shape | `fibre_shape.*` and its device twin `gpu/fibres/fibre_candidate.cuh` |
| A rule | `growCandidate` and the stage of the same name, or `finishLight`, in `gpu/fibres/fibres.cu`. A new reason also needs `FibreOutcome` (the kernels return it too) and `FibreReport::count`. |
| A `FibreParams` field | `fibre.h`, `fibreParamsText` in `fibre.cpp`, the flat parameters in `gpu/fibres/fibres.h` and `fibres.cu`. A setting also needs `settings/scene_settings.*` and `scene/stage_params.cpp`. |
| A field of the record | `fibre_record.h`, the kernels, and `fibreLine` in `fibre_checks.cpp`. Saved fibres have their own format in `scene/fibres/fibre_files.cpp`. |
| The budget | `fibre_budget.*`; the scene's `generateInBudgetOrder` drives it. Its result must stay exact: the same `q` as evaluating everything. |
| What the export writes for a fibre | `fibre_strands.*`; the file format itself in `strands_file.*`, beside the renderer's reader. |
| A groom step | `fur.*`, then re-pin `fur_strands_pinned` on purpose. |

## Literature

- Zhao, Luan, Bala 2016, *Fitting Procedural Yarn Models for Realistic Cloth Rendering*: yarn
  coordinates, flyaways as a separate population.
- Wu, Yuksel 2017, *Real-time Fiber-level Cloth Rendering*: fibres generated on the GPU, thinning as
  level of detail.
- Luan, Zhao, Bala 2017, *Fiber-Level On-the-Fly Procedural Textiles*: a fibre as a function of its
  identity.
- Müller et al. 2012 and Han, Harada 2012: Follow-The-Leader with a shape constraint, fur's settle.
- Trunz et al. 2024, *Neural Inverse Procedural Modeling of Knitting Yarns from Images*: plies swept in
  the helix's local frame (`yarn_body`); Zhao et al.'s migration (Eq. 3) within a shell.
