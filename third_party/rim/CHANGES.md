# RIM in Wool Shader

RIM is the field-aligned remesher behind the stitch-mesh experiment's Method B (`src/stitchmesh`): Gao, Jakob,
Tarini, Panozzo 2017, *Robust Hex-Dominant Mesh Generation using Field-Guided Polyhedral Agglomeration*.
CMake downloads it at a pinned commit (`cmake/dependencies.cmake`); none of its source is kept in this repository.

| | |
|---|---|
| upstream | https://github.com/gaoxifeng/robust_hex_dominant_meshing, commit `83dd574d5e87a4491d0cf6d5970a4117a8ab1845` |
| licence | MIT, Copyright (c) 2025 Xifeng Gao (its `LICENSE`, shown in the viewer's About window) |
| compiled | `hierarchy`, `orientations`, `positions`, `extract`, `dedge`, `subdivide`, `bvh` (surface path only) |

Three of those units and their headers (`dedge`, `subdivide`, `bvh`, with `aabb.h`, `bvh.h`, `dedge.h`,
`subdivide.h`) come from **Instant Meshes** (Jakob, Panozzo, Tarini, Sorkine-Hornung 2015) and keep its BSD
licence; RIM's tree does not include that licence's text, so it is here: `INSTANT_MESHES_LICENSE.txt` (from
https://github.com/wjakob/instant-meshes, `LICENSE.txt`), also shown in the About window. `hierarchy` includes
`tri_tri_intersection.h`, Tomas Möller's triangle-triangle test (1997; the "opttritri" variant, in double precision)
as RIM ships it, without a notice; its published code is MIT: `MOLLER_TRITRI_LICENSE.txt`, also in the About window. RIM's TetGen (AGPL) is
not compiled: the volume path is removed by the patch below.

The files in this folder are Wool Shader's own: the build patch, the driver, the hooks and serial stand-ins.

- `rim.cmake` — the patch applied to the downloaded sources at configure time. Every change, as its header lists:
  - TetGen removed (the include and the body of the volume path's `tet_meshing`);
  - an Eigen 3.4 type fix (the per-face centres are `Vector4f` = `Vector4d`, not Eigen's float);
  - `system("PAUSE")` on three inconsistency paths becomes an exception (headless: the remesh fails loudly);
  - the seed: `build()`'s random start reads `knit_rim_seed`; 0 keeps pcg32's default state, RIM's own start;
  - the 2-RoSy boundary: boundary vertices smoothed and snapped along or across the boundary;
  - two hook calls into `knit_rim_hooks.h` (written from Wu et al. 2018 section 4): the 2-RoSy orientation
    compatibility in `smoothOrientationsTri` and the post-processing inside `meshExtraction2D`'s collapse / split
    loop;
  - the rectangular lattice (an aspect ratio, inactive at 1);
  - three crash fixes: `simple_polygon` (an empty boundary; a discarded vertex on the boundary) and
    `remove_doublets2D` (a shared-edge chain without two end vertices).
- `rim_remesh.h`, `rim_remesh.cpp` — the driver behind a plain-array interface. It follows RIM's batch mode
  (`batch.cpp`: the order of hierarchy, orientation and position sweeps, extraction, subdivision); `retag()`
  follows the first loop of RIM's `init_edge_tagging2D` (`extract.cpp`) and the subdivision uses RIM's batch-mode
  rule. Those parts are adapted from RIM and are covered by its MIT licence.
- `knit_rim_hooks.h` — the Wu 2018 section 4 hooks the patched sources call.
- `prelude.h`, `shim/` — the prelude every unit is compiled with (namespace `rim`, SLIM's Eigen), and serial
  stand-ins for the parts of TBB and nanogui RIM includes.
