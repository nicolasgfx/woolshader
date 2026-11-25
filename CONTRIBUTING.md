# Contributing

## Build and test

The build targets Windows with Visual Studio 2022, CMake, Ninja and the CUDA toolkit. Builds live outside
the checkout.

```powershell
.\build.bat              # configure if needed and build, in ..\<checkout>-build\release
.\build.bat test         # build and run CTest
.\build.bat run          # build and start the viewer
.\run.bat                # the same as build.bat run; arguments go to the viewer, e.g. --scene wool-swatch
```

`CMakePresets.json` gives the same with plain CMake:

```powershell
cmake --preset release
cmake --build --preset release
ctest --preset release
```

A build folder belongs to one checkout: `build.bat` refuses a folder configured for another one. Set
`KNIT_BUILD_DIR` to choose a folder, `KNIT_CUDA_ROOT` for another CUDA Toolkit and `KNIT_CMAKE_ARGS` for extra
configure options of a new build folder. Long runs go through `tests\capped_run.ps1`, which enforces a hard
time cap.

Tests that need the separate renderer (`photon_path_tracer`) skip without it; set `KNIT_PPT_EXE` to a `ppt.exe`
before configuring to run them. Tests that need an example model the repository does not ship (the teddy bears)
skip as well. CTest switches the unreleased Fur workflow on for every case (`KNIT_FUR=1`); the viewer starts with
it off.

## Where code goes

[src/ARCHITECTURE.md](src/ARCHITECTURE.md) is the way into the code and names the owner of each concern; each
component folder has its own `ARCHITECTURE.md`. Modules include only the ones below them, and
`architecture_module_dependencies` checks it. [STYLE.md](STYLE.md) holds the code style: concrete, readable
code, one implementation of each thing.

Every change comes with a test of its behaviour. A test that cannot run for a missing fixture calls
`SKIP(reason)`; it never returns early as a pass.
