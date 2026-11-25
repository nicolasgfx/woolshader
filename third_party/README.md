# Third-party code

CMake downloads each dependency at a pinned revision (`cmake/dependencies.cmake`), and this folder holds the
project's own integration with them. The one vendored source is `slim/src`: eight MPL-2.0 files of the SLIM solver
(`slim/README.md` names their origin and hashes, `slim/CHANGES.md` every change made to them).
Offline: `-DFETCHCONTENT_SOURCE_DIR_<NAME>=<path>`.

| dependency | used for | licence | revision |
|---|---|---|---|
| [cgltf](https://github.com/jkuhlmann/cgltf) | GLB reading | MIT | `52c23814` (v1.14) |
| [stb](https://github.com/nothings/stb) | image decoding, PNG writing | MIT / public domain | `2c980bb5` |
| [OpenVDB](https://github.com/AcademySoftwareFoundation/openvdb) | the import's voxel remesh | Apache-2.0 | v13.0.0 |
| [oneTBB](https://github.com/uxlfoundation/oneTBB) | OpenVDB only | Apache-2.0 | v2022.3.0 |
| SLIM solver (Rabinovich et al. 2017), eight files from [Blender's `intern/slim`](https://github.com/blender/blender/tree/v5.1.1/intern/slim/intern), which derives them from libigl | the import's flatten | MPL-2.0 | v5.1.1, vendored in `slim/src`, upstream SHA-256 in `slim/README.md` |
| [Eigen](https://gitlab.com/libeigen/eigen), a pinned commit, fetched as `slim_eigen` | SLIM and the import's flatten, under the namespace name `SlimEigen` | MPL-2.0 | `8a1083e9bf41` |
| [GLFW](https://github.com/glfw/glfw) | viewer window | zlib | 3.3.9 |
| [Dear ImGui](https://github.com/ocornut/imgui) | viewer interface | MIT | v1.91.8 |
| [RIM](https://github.com/gaoxifeng/robust_hex_dominant_meshing) (robust_hex_dominant_meshing), patched at configure time; its `dedge`, `subdivide`, `bvh` units are [Instant Meshes](https://github.com/wjakob/instant-meshes) code (BSD, text in `rim/INSTANT_MESHES_LICENSE.txt`) | the stitch-mesh experiment's field-aligned remesher (Method B) | MIT (RIM), BSD (Instant Meshes parts) | `83dd574d`, `rim/CHANGES.md` |
| [pcg32](https://github.com/wjakob/pcg32) | RIM's random numbers | Apache-2.0 | `0ef13e68` |
| [HiGHS](https://github.com/ERGO-Code/HiGHS), with its bundled filereaderlp (MIT) and pdqsort (zlib) | the stitch-mesh experiment's integer programs (Method B's labels and directions) | MIT | v1.11.0 (`364c83a5`) |

The CUDA toolkit is installed separately (see "Build and run" in the top-level README).

The viewer's About window (status bar, About) shows each of these with its licence text, compiled in from the fetched
sources (`src/viewer/about.h`, `cmake/embed_texts.cmake`). A new dependency needs a notice there: `ui_about_licence_texts`
fails for a `FetchContent_Declare` without one.

Blender is not a build dependency. The viewer starts an installed Blender as a separate program only for the
Blender file (the `.blend` export), and works without it. The one script it runs is `blender/ppt_package_to_blend.py`,
under GPL-3.0-or-later (`blender/README.md`).

## Local parts

- `cgltf_impl.cpp`, `stb_impl.cpp`: the single implementation units of the header-only libraries.
- `vdb/`: the one call into OpenVDB, behind plain arrays, so no OpenVDB or TBB header reaches the project.
- `slim/src/`: the eight MPL-2.0 files of the SLIM solver, each with its own copyright and licence header. Every change made
  to them stays under the MPL-2.0 and is listed in `slim/CHANGES.md`; the About window shows that list with the licence text.
- `slim/shim/`: our own stand-ins (MIT) for the three small headers the SLIM files include: an assert, pi and a serial loop.
  Only the eight files in `slim/src` are third-party code.
- `rim/`: the build patch for the downloaded RIM sources (`rim.cmake`), the driver behind plain arrays (`rim_remesh.*`,
  partly adapted from RIM's batch mode), the Wu et al. 2018 hooks and serial stand-ins for TBB and nanogui (`shim/`).
  `rim/CHANGES.md` lists every change and the licences involved; no RIM source is kept here.
