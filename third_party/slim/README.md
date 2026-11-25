# SLIM solver (MPL-2.0)

`src/` holds eight files of the SLIM solver (Rabinovich, Poranne, Panozzo, Sorkine-Hornung 2017, *Scalable Locally Injective
Mappings*). They were taken from Blender 5.1.1, `intern/slim/intern`
(https://github.com/blender/blender/tree/v5.1.1/intern/slim/intern), which derives them from libigl. They are not Wool
Shader's own code: each file keeps its own copyright and MPL-2.0 header. The Mozilla Public License 2.0 applies to these
files and to every change made to them; the rest of Wool Shader is MIT. The licence text is in `LICENSE` beside this file
(also shown in the viewer's About window).

| file | SHA-256 of the file as it was vendored, before the changes in CHANGES.md |
|---|---|
| slim.cpp | ab3bba6d4f84fd693ec93441ab81b00a538095055a3f0b2108da7ed5195224bd |
| slim.h | 299faa9325534756f09ce769bf667efe0af3b5092dd40dbfb729c21564e7c55d |
| flip_avoiding_line_search.cpp | 35a2af73db3e9b50be59f49b5359044e1c9158ae230cf1e394e354c22ff456c4 |
| flip_avoiding_line_search.h | c95b3cd801cfea6eeb08be756dae74b0e470814a20248403853b3e7c4258aab0 |
| doublearea.cpp | b3e0434c658ef6a3825f3da34ea22f8d23486c48e3d227f1ebee37df1d98867b |
| doublearea.h | 0067400560075358c7d808894f67221022ab51ddf7c1102dbe43d7f42271aeb7 |
| edge_lengths.cpp | 1ac7c0282958b359b8f0fc61f3bfdaf73f594b404e5450778586e790e645ebb9 |
| edge_lengths.h | a0329ab340b8fbd8265d3580e301b80d3e7a1d5876425e889c273117f329a927 |

`CHANGES.md` lists every change made here. `shim/` holds our own stand-ins (MIT) for the three small headers the files
include: `slim_assert.h` (an assert that checks nothing), `slim_math_base.h` (pi, through `<cmath>`) and `slim_parallel.hh`
(a serial loop in the namespace `slim_support`); the files need nothing else besides Eigen.

How they build (`cmake/dependencies.cmake`): the library `knit_slim` compiles `slim.cpp`, whose headers include the other
three `.cpp` files, with `shim/` on its private include path. It uses a pinned Eigen commit, fetched as `slim_eigen` and
compiled under the namespace name `SlimEigen` (`Eigen=SlimEigen`), so it cannot meet the Eigen 3.4 of the other libraries
at link time. Only `knit_slim` and the import's `knit_flatten` see that Eigen.
