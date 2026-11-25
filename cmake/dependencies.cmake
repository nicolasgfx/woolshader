# The pinned third-party dependencies (third_party/README.md). Offline: -DFETCHCONTENT_SOURCE_DIR_<NAME>=<path>.
# S6: cgltf (header-only glTF / GLB reader), pinned to the v1.14 commit. SOURCE_SUBDIR
# points at a directory without a CMakeLists.txt, so nothing is added as a subproject.
# Offline: -DFETCHCONTENT_SOURCE_DIR_CGLTF=<path>
include(FetchContent)
FetchContent_Declare(cgltf
    GIT_REPOSITORY https://github.com/jkuhlmann/cgltf.git
    GIT_TAG        52c23814dbb60c2fbb54750ddf41d342d432a498
    SOURCE_SUBDIR  no-cmake)
FetchContent_MakeAvailable(cgltf)

# Mesh import (src/io/mesh_read.cpp): stb_image (MIT / public domain) decodes base-colour textures. Pinned commit;
# SOURCE_SUBDIR as for cgltf (no subproject). Offline: -DFETCHCONTENT_SOURCE_DIR_STB=<path>
FetchContent_Declare(stb
    GIT_REPOSITORY https://github.com/nothings/stb.git
    GIT_TAG        2c980bb59875b0d32144a71867fbdebb2f77cd20
    SOURCE_SUBDIR  no-cmake)
FetchContent_MakeAvailable(stb)

# S6: export, texture maps, scene settings / registry / pipeline (scene/). Shared by
# knit_check, the viewer and the tests; model, knitting and fibres stay standard-library only.
# The built-in import's voxel remesh (third_party/vdb): OpenVDB 13.0.0 (Apache-2.0) on oneTBB v2022.3.0 (Apache-2.0). Static,
# with no file formats: Blosc, zlib, OpenEXR and Boost stay out. Only knit_vdb sees
# their headers. Offline as above.
set(TBB_TEST OFF CACHE BOOL "" FORCE)
set(TBB_STRICT OFF CACHE BOOL "" FORCE)
set(TBBMALLOC_BUILD OFF CACHE BOOL "" FORCE)
set(TBBMALLOC_PROXY_BUILD OFF CACHE BOOL "" FORCE)
set(USE_BLOSC OFF CACHE BOOL "" FORCE)
set(USE_ZLIB OFF CACHE BOOL "" FORCE)
set(USE_IMATH_HALF OFF CACHE BOOL "" FORCE)
set(USE_NANOVDB OFF CACHE BOOL "" FORCE)
set(USE_CCACHE OFF CACHE BOOL "" FORCE)
set(OPENVDB_BUILD_BINARIES OFF CACHE BOOL "" FORCE)
set(OPENVDB_CORE_SHARED OFF CACHE BOOL "" FORCE)
set(OPENVDB_CORE_STATIC ON CACHE BOOL "" FORCE)
set(OPENVDB_USE_DELAYED_LOADING OFF CACHE BOOL "" FORCE)
set(OPENVDB_ENABLE_UNINSTALL OFF CACHE BOOL "" FORCE)
set(OPENVDB_INSTALL_CMAKE_MODULES OFF CACHE BOOL "" FORCE)
set(OPENVDB_FUTURE_DEPRECATION OFF CACHE BOOL "" FORCE)
set(DISABLE_DEPENDENCY_VERSION_CHECKS ON CACHE BOOL "" FORCE)
# A static oneTBB loads no DLL; stated explicitly, oneTBB no longer warns that the check is off by default.
set(TBB_VERIFY_DEPENDENCY_SIGNATURE OFF CACHE BOOL "" FORCE)
# Only OpenVDB's programs and unit tests (both off) link an allocator; "Auto" warns that Jemalloc is missing.
set(CONCURRENT_MALLOC None CACHE STRING "" FORCE)
FetchContent_Declare(tbb
    GIT_REPOSITORY https://github.com/uxlfoundation/oneTBB.git
    GIT_TAG        v2022.3.0
    GIT_SHALLOW    TRUE
    OVERRIDE_FIND_PACKAGE)
FetchContent_Declare(openvdb
    GIT_REPOSITORY https://github.com/AcademySoftwareFoundation/openvdb.git
    GIT_TAG        v13.0.0
    GIT_SHALLOW    TRUE)
set(BUILD_SHARED_LIBS OFF)
set(CMAKE_MSVC_RUNTIME_LIBRARY MultiThreadedDLL)  # left unset, OpenVDB's static library imposes /MT on all who link it
# oneTBB and OpenVDB warn unconditionally about three of the choices above: a static oneTBB, the dynamic CRT under a
# static OpenVDB, and no Blosc. While they configure, exactly those three warnings are dropped (the message() override
# below), and so are deprecation notices (OpenVDB's: it will later need CMake 4). Every other message passes through.
function(message)
    set(text "")
    foreach(i RANGE 1 ${ARGC})  # each argument read on its own keeps its semicolons
        if(i LESS ARGC)
            string(APPEND text "${ARGV${i}}")
        endif()
    endforeach()
    if(KNIT_QUIET_DEPENDENCY_ADVISORIES AND ARGV0 STREQUAL "WARNING" AND (
            text MATCHES "^You are building oneTBB as a static library\\." OR
            text MATCHES "^Unexpected value for the Windows CRT.*use the dyanmic CRT \\(/MD\\) with a static configuration" OR
            text MATCHES "^Blosc support is disabled\\."))
        return()
    endif()
    if(ARGV0 MATCHES "^(FATAL_ERROR|SEND_ERROR|WARNING|AUTHOR_WARNING|DEPRECATION|NOTICE|STATUS|VERBOSE|DEBUG|TRACE|CHECK_START|CHECK_PASS|CHECK_FAIL|CONFIGURE_LOG)$")
        _message(${ARGV0} "${text}")
    else()
        _message("${ARGV0}${text}")
    endif()
endfunction()
set(KNIT_QUIET_DEPENDENCY_ADVISORIES ON)
set(CMAKE_WARN_DEPRECATED OFF)
FetchContent_MakeAvailable(tbb openvdb)
unset(CMAKE_WARN_DEPRECATED)
unset(KNIT_QUIET_DEPENDENCY_ADVISORIES)
unset(CMAKE_MSVC_RUNTIME_LIBRARY)
unset(BUILD_SHARED_LIBS)
add_library(knit_vdb STATIC third_party/vdb/level_set_remesh.cpp)
target_include_directories(knit_vdb PUBLIC third_party/vdb)
target_link_libraries(knit_vdb PRIVATE openvdb_static)
target_compile_definitions(knit_vdb PRIVATE NOMINMAX __TBB_NO_IMPLICIT_LINKAGE=1)
if(MSVC)
    target_compile_options(knit_vdb PRIVATE /bigobj)
endif()


# The import's flatten: the SLIM solver (Rabinovich et al. 2017, derived from libigl). Its eight files are MPL-2.0 (their own
# headers) and are vendored in third_party/slim/src, with our changes listed in third_party/slim/CHANGES.md; changes to them
# stay under the MPL. The three small headers they include are our own stand-ins (third_party/slim/shim).
#
# They build on a pinned Eigen commit (MPL-2.0). Eigen 3.4 stays for the other libraries, so this one is compiled under its own
# namespace name, and nothing of it can meet them at link time. Only knit_slim and knit_flatten (src/import/CMakeLists.txt) see it.
FetchContent_Declare(slim_eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/8a1083e9bf41b91fdea6546681f806154efdc25a/eigen-8a1083e9bf41b91fdea6546681f806154efdc25a.tar.gz
    URL_HASH SHA256=cc28a84fdec496c6777596350ea805519bf10f717d21044ae6ba3dd562183a26
    SOURCE_SUBDIR only_the_headers_are_used)
FetchContent_MakeAvailable(slim_eigen)
add_library(knit_slim STATIC third_party/slim/src/slim.cpp)
target_include_directories(knit_slim PRIVATE third_party/slim/shim)
target_include_directories(knit_slim SYSTEM PUBLIC third_party/slim/src ${slim_eigen_SOURCE_DIR})
target_compile_definitions(knit_slim PRIVATE _USE_MATH_DEFINES NOMINMAX PUBLIC Eigen=SlimEigen)
if(MSVC)
    target_compile_options(knit_slim PRIVATE /W0 /bigobj)
endif()

# The stitch-mesh experiment (src/stitchmesh; ported from knit/stitch-meshing-v2-recovery dbf5f41): RIM, the remesher of Wu et al. 2018 (Gao, Jakob, Tarini, Panozzo 2017; MIT,
# (c) Xifeng Gao) at a pinned commit, surface path only, patched (third_party/rim/rim.cmake: TetGen removed, a type fix, the
# seed, two hooks into third_party/rim/knit_rim_hooks.h). Its sources are compiled inside namespace rim (one generated unit
# per source, third_party/rim/prelude.h) against SLIM's Eigen under its namespace name, so one Eigen is in the program;
# pcg32 (its generator, Apache-2.0) at the commit RIM's submodule names; TBB and nanogui are serial stand-ins
# (third_party/rim/shim), so a remesh is deterministic. C++14, as upstream. Only knit_rim sees their headers. Offline as above.
# Three of the compiled units (dedge, subdivide, bvh) are Instant Meshes code under its BSD licence, whose text RIM's tree
# lacks: third_party/rim/INSTANT_MESHES_LICENSE.txt (third_party/rim/CHANGES.md).
FetchContent_Declare(rim
    GIT_REPOSITORY https://github.com/gaoxifeng/robust_hex_dominant_meshing.git
    GIT_TAG        83dd574d5e87a4491d0cf6d5970a4117a8ab1845
    SOURCE_SUBDIR  no-cmake
    # KNIT_RIM_PATCH: the patch's revision (a new revision runs the patch again on an existing source tree).
    PATCH_COMMAND  ${CMAKE_COMMAND} -DKNIT_RIM_PATCH=6 -P ${CMAKE_CURRENT_SOURCE_DIR}/third_party/rim/rim.cmake)
FetchContent_Declare(pcg32
    GIT_REPOSITORY https://github.com/wjakob/pcg32.git
    GIT_TAG        0ef13e68ca0be5506e1cfc0db76831e6f916e9e9
    SOURCE_SUBDIR  no-cmake)
FetchContent_MakeAvailable(rim pcg32)
set(KNIT_RIM_UNITS)
foreach(unit hierarchy orientations positions extract dedge subdivide bvh)
    file(CONFIGURE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/knit_rim/unit_${unit}.cpp
         CONTENT "// Generated (cmake/dependencies.cmake): RIM's ${unit}.cpp inside namespace rim (third_party/rim/prelude.h).\n#include \"prelude.h\"\nnamespace rim {\n#include \"${unit}.cpp\"\n}\n")
    list(APPEND KNIT_RIM_UNITS ${CMAKE_CURRENT_BINARY_DIR}/knit_rim/unit_${unit}.cpp)
endforeach()
add_library(knit_rim STATIC third_party/rim/rim_remesh.cpp ${KNIT_RIM_UNITS})
set_target_properties(knit_rim PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)
target_include_directories(knit_rim PUBLIC third_party/rim)
target_include_directories(knit_rim PRIVATE third_party/rim/shim ${rim_SOURCE_DIR}/src)
target_include_directories(knit_rim SYSTEM PRIVATE ${slim_eigen_SOURCE_DIR} ${pcg32_SOURCE_DIR})
target_compile_definitions(knit_rim PRIVATE _USE_MATH_DEFINES NOMINMAX _CRT_SECURE_NO_WARNINGS Eigen=SlimEigen)
if(MSVC)
    # RIM needs the permissive mode: the project-wide /permissive- comes off this target (D9025 "overriding" otherwise).
    get_target_property(KNIT_RIM_OPTIONS knit_rim COMPILE_OPTIONS)
    list(REMOVE_ITEM KNIT_RIM_OPTIONS "$<$<COMPILE_LANGUAGE:CXX>:/permissive->")
    set_target_properties(knit_rim PROPERTIES COMPILE_OPTIONS "${KNIT_RIM_OPTIONS}")
    target_compile_options(knit_rim PRIVATE /W0 /bigobj /permissive)
endif()

# The integer solver of Method B's label and direction programs (W18 sections 5, 6; ported from dbf5f41): HiGHS (MIT)
# v1.11.0, static, without its executable, tests or zlib; its BUILD_TESTING is scoped to it. It compiles its bundled
# filereaderlp (MIT, Michael Feldmeier; its notice is in the About window) and pdqsort (zlib). Only stitchmesh links it,
# privately. Offline: -DFETCHCONTENT_SOURCE_DIR_HIGHS=<path>.
set(BUILD_TESTING OFF)
set(BUILD_CXX_EXE OFF)
set(ZLIB OFF)
set(FAST_BUILD ON)
FetchContent_Declare(highs
    GIT_REPOSITORY https://github.com/ERGO-Code/HiGHS.git
    GIT_TAG        364c83a51e44ba6c27def9c8fc1a49b1daf5ad5c
    SYSTEM)
FetchContent_MakeAvailable(highs)
unset(BUILD_TESTING)

# Viewer (optional). Offline: -DFETCHCONTENT_SOURCE_DIR_GLFW=<path> -DFETCHCONTENT_SOURCE_DIR_IMGUI=<path>
if(KNIT_BUILD_VIEWER)
    include(FetchContent)
    set(GLFW_BUILD_DOCS     OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_TESTS    OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(GLFW_INSTALL        OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(glfw
        GIT_REPOSITORY https://github.com/glfw/glfw.git
        GIT_TAG        3.3.9
        GIT_SHALLOW    TRUE)
    FetchContent_Declare(imgui
        GIT_REPOSITORY https://github.com/ocornut/imgui.git
        GIT_TAG        v1.91.8
        GIT_SHALLOW    TRUE)
    FetchContent_MakeAvailable(glfw imgui)

    find_package(OpenGL REQUIRED)
endif()
