cmake_minimum_required(VERSION 3.24)
project(knit_proto LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(KNIT_BUILD_VIEWER "Build the GLFW/ImGui viewer" ON)
option(KNIT_BUILD_TESTS  "Build the geometry tests"      ON)

if(MSVC)
    # C++ only (nvcc receives its host flags through -Xcompiler, see knit_gpu).
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:/permissive-> $<$<COMPILE_LANGUAGE:CXX>:/utf-8>)
    # Strict IEEE semantics: the interval/allowance logic relies on round-to-nearest
    # and on the compiler not reassociating floating-point expressions.
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:/fp:precise>)
endif()

# Warning level for this project's own targets only: a global /W4 also reaches the
# fetched GLFW / ImGui targets, which set their own level (D9025 "overriding /W4"). C++ only
# (the CUDA sources of knit_gpu take their host flags through -Xcompiler).
function(knit_warnings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE $<$<COMPILE_LANGUAGE:CXX>:/W4> $<$<COMPILE_LANGUAGE:CXX>:/wd4100>)
    endif()
endfunction()

# Executables and their DLLs sit at the top of the build tree, whatever directory defines them.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

enable_language(CUDA)
find_package(CUDAToolkit REQUIRED)

include(cmake/dependencies.cmake)
include(cmake/components.cmake)
add_subdirectory(src)
if(KNIT_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
knit_components_finish()
