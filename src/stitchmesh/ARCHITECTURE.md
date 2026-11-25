# src/stitchmesh: the stockinette stitch-mesh experiment

The stages of the stockinette stitch-mesh experiment that come before any yarn: the reference surface of the simple case, Method A's row-time field, rings and guided rows, Method B's W18 pipeline (ported from the recovery branch dbf5f41: remesh, labels, directions, subdivision), the typed stitch mesh both methods produce with its independent check, the shift-path's virtual helical connectivity, coverage, and the construction fixtures with their negative controls. Checkpoint 3 adds the mesh relaxation of Y12 section 6.1. The imported arm adds the marking records, Method A driven by them (a pole start, a finish level, a row-start path on the surface), the pole-start face kind and its fixtures, and the relaxation's marking constraints. The first-yarn checkpoint adds yarn embedded in the typed mesh and its checks; the closeout adds the adapter that hands that yarn to the existing knitting solver, and the guided Knit route (scene/knit/guided_part.cpp) adds the settle itself (moved here from the lab, so the lab and the app run one implementation) and the stitch records production's later stages read.

- One component, `knit_stitchmesh`, at level 2.
- Public dependencies: `model` and `knitting`. `yarn_settle.h` takes the knitting solver (`KnitNumerics`) and returns its validation reports; the adapter interpolates with `knitting/tube/plies.h`.
- The remesher RIM (`third_party/rim`, library `knit_rim`, with pcg32 and SLIM's Eigen) and the integer solver HiGHS are private. No RIM, Eigen or HiGHS type appears in a header.
- There is no file IO: the lab (`src/cli/stitchmesh_lab.cpp`) reads and writes every file.

## File tree

```text
src/stitchmesh/
  ARCHITECTURE.md          this file
  CMakeLists.txt           the component declaration
  revolution_surface.h/.cpp  the simple reference surface: the exact profile, meridian table, normals, curvatures, and its triangle mesh with variant J
  time_field.h/.cpp        Method A: geodesic distance from source vertices by fast marching (Kimmel and Sethian 1998)
  row_rings.h/.cpp         Method A: the level sets of the time field as oriented contours with bindings; the crossing with a cut half-plane
  orientation_field.h/.cpp Method B: RIM's 2-RoSy orientation field, solved or fixed by an axis guide, without positions or extraction
  guided_rows.h/.cpp       Method A: ring vertices by arc fraction from the cut guide, the slot lookup, shaping placement (logged), faces, shift-path;
                           buildMarkedRows: the same from markings (pole start, finish level, the row-start path's crossings), the adjacent-shaping variant
  reference_markings.h/.cpp  knit::refmark: the imported arm's reference-bound records of one prepared component (start, finish, growth, row-start path, cut, attachment) bound by
                           (triangle, barycentrics) to a reference digest; read, write, validate, densify a path, checked rebinding
  polygon_mesh.h/.cpp      ported (dbf5f41): polygon faces, half-edges, the polygon check, W18 sec. 4 normalization
  reference_surface.h/.cpp ported: the reference triangles prepared for closest-point queries and vertex bindings
  stage_report.h/.cpp      ported: identities, counts, time and verdict of a stage
  quantiles.h/.cpp         ported (model/quantiles at dbf5f41): the five-number summary the remesh qualification reports
  remesh_stage.h/.cpp      Method B, ported: RIM's field-aligned coarse mesh with its structural check and qualification (W18 sec. 4)
  labels.h/.cpp            Method B, ported: course / wale per coarse edge by HiGHS, with the repairs (W18 sec. 5)
  directions.h/.cpp        Method B, ported: the wale direction per row by HiGHS, anchored by the growth guide (W18 sec. 6)
  stitch_mesh.h/.cpp       Method B, ported: the subdivision into K, kp, d12k and special quads, strands, loops, its own validator (W18 sec. 7)
  stitch_fixtures.h/.cpp   ported: dbf5f41's hand-authored coarse fixtures (the special-quad neighbourhood uses fixture 2)
  typed_mesh.h/.cpp        the typed mesh of both methods (roles, boundary roles, strips) and its independent check: connectivity, geometry, cells
  virtual_shift.h/.cpp     the shift-path's virtual helical connectivity (W19 sec. 3.4) and its validator, which traces the virtual polygons
  coverage.h/.cpp          coverage of the reference: uncovered area, double coverage, projected area, boundary correspondence
  construction_fixtures.h/.cpp  the CP2 fixtures and negative controls, and the arm's (pole start, adjacent shaping), shared by the lab and the tests
  relax_mesh.h/.cpp        CP3: Y12 sec. 6.1 mesh relaxation (backward Euler, filtered CG; stretch, shear, wale strut; sub-quads) on the reference by a directional walk;
                           optional marking constraints (pinned vertices, guide lines) and the smooth filter normals of solver note R1
  border_placement.h/.cpp  the border rings placed on their boundary loops (BorderRings::PlacedFirst / PlacedOnLoop), held for the relaxation
  surface_trim.h/.cpp      markings: a piece split along closed marking curves into regions; the band between start and finish as its own mesh
  marked_build.h/.cpp      markings: Method A on a marked piece, one backend for the app and the lab: the prepared reference, the band route
                           (start curve; a closed start or finish on a border loop is that opening itself, untrimmed) and the experimental
                           pole route (start point; finish curve or level), offered proposals only as a labelled research preview,
                           relaxation with the markings' constraints, coverage, the review and the lab's evaluation; the shared marking
                           resolution, the row field it keeps (MarkedField) and the row-crossing count and findings the quick review reuses;
                           cancellation checkpoints between phases
  prepared_reference.h/.cpp  markings: the piece midpoint-refined (facets) or Loop-subdivided to the limit (the smooth-reference variant), one split with an exact map back to the file's triangles; clipping below a field level
  mark_proposals.h/.cpp    markings: proposals for one piece from its own geometry, one backend for the app and the lab: compatible
                           setups (start, termination, row-start path from one field) for an open tube (its border loops), a cup
                           (one opening: a pole where the opening is most nearly one row, the opening as the bind-off), a closed
                           piece hidden at one end, a closed piece nothing hides (visible caps disclosed) and a branching piece
                           (refused at the split, partial setup); rule 2's trial-relaxation refinement on request; the lab's
                           nine-marking file
  quick_review.h/.cpp      markings: the quick review while drawing, without a build: the shared resolution and checks, a path or
                           curve crossing itself, and the rows checked against a kept field only when its key matches (otherwise
                           "row checks pending"); it marks which measures it knows
  yarn_embed.h/.cpp        the first yarn experiment: canonical templates (K, kp, d12k, gathered pole start) mapped into the virtual faces
                           (shared ports, mean value coordinates, normal offsets), bind-off caps, the binding thread, the pole-round join; the
                           strands; the yarn checks (ports, continuity, crossing order per face in the reference frame and in the yarn's own
                           best-fit plane, interlocks, yarn ends, clearance, bends); tubes
  yarn_curves.h/.cpp       the yarn-settle adapter: the embedded strand as production's uniform cubic B-spline (equal-arc points, natural
                           ends) with each piece's parameter interval, the measured conversion error, and the embedding re-evaluated on moved
                           curves (the closeout's settle)
  yarn_settle.h/.cpp       the settle: conversion, raw/converted/settled checks, radius continuation on KnitNumerics::relax with
                           the certified gate after every call (closeout settings); one implementation for stitchmesh_lab and the app;
                           a failed crossing order is refused (lab), returned unsettled without a relax call (app) or, on the
                           lab's explicit opt-in, settled and marked; a broken strand or port is always refused
  guided_records.h/.cpp    production stitch records derived from the stitch mesh and its yarn: one per loop, parents from the heads,
                           frames; the counts cross-checked against the yarn check (the local signature does not apply to them)
```

**Reading order:** `revolution_surface.h` (what is measured against what), `time_field.h`, `row_rings.h`, `orientation_field.h`, then `typed_mesh.h` (the shared record and check), `guided_rows.h` (Method A), `stitch_mesh.h` (Method B's subdivision), `virtual_shift.h`, `coverage.h`, `construction_fixtures.h`, `relax_mesh.h`, then `surface_trim.h`, `prepared_reference.h` and `marked_build.h` (the
Mark step's input), and `reference_markings.h` (the imported arm's reference-bound records); last `yarn_embed.h` (yarn on the typed mesh) and `yarn_curves.h` (that yarn as production curves for the existing solver), `yarn_settle.h` (the settle on that solver) and `guided_records.h` (the stitch records of the settled yarn).

## Rules

- Exact quantities (arc length, normals, curvatures) come from `revolution_surface` and are used only to measure the computed fields. The computed fields never read them.
- A level with anything but one closed contour is a contour event. It is reported, not repaired (`row_rings.h`).
- `computeOrientationField` refuses an input that RIM would subdivide, so the field is always per input vertex.
- The typed mesh's check re-derives kinds, interfaces, strands and the count invariant from the faces and roles; it never reads the builder's bookkeeping. Connectivity, geometry and cell dimensions are three separate verdicts; cell bounds are reported, not a CP2 verdict.
- The virtual-shift validator rebuilds adjacency from the re-read polygons and traces the yarn; it never applies the shift rule itself.
- The relaxation moves vertices only, along the reference surface (border vertices along their reference loop); it never retypes a face or changes the topology. It keeps the order of the border vertices around their loop, gives every stored border edge its spring, turns no clearly convex face corner over and takes back a step that raises the forces. The row layout turns each ring for upright wales (`RowPhase`), and the marked build places the border rings on their loops before the relaxation and holds them. Its solver faces are given by the caller (Method A: the virtual polygons).
- Markings are data: a set names one component and one reference digest, and is refused on any other until it is rebound by the checked rebinding. The mesher reads them from the file; a row-start path never cuts the surface and an attachment region never makes a join.
