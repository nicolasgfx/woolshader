# Files, the formats several components share, the caches' places and the Blender process of the .blend export
# (ARCHITECTURE.md). Every header is public. cgltf and stb are private: only glb_load.cpp and
# image_io.cpp include them. KNIT_SOURCE_DIR reaches sourceTreeDir() alone.
knit_component(io
    SOURCES
        blender_process.cpp cache_paths.cpp child_process.cpp files.cpp glb_load.cpp image_io.cpp marking_file.cpp mesh_read.cpp preview_pipe.cpp publish.cpp stitch_mesh_obj.cpp
        ${CMAKE_SOURCE_DIR}/third_party/cgltf_impl.cpp ${CMAKE_SOURCE_DIR}/third_party/stb_impl.cpp
    PUBLIC_HEADERS
        blender_process.h cache_paths.h child_process.h files.h glb_load.h glb_selection.h image_io.h marking_file.h mesh_read.h preview_pipe.h publish.h stitch_mesh_obj.h
    PUBLIC_DEPS model
)
# cgltf and stb as system includes: their implementations compile without warnings (/external:W0) while
# the units themselves stay at the project's level.
target_include_directories(knit_io SYSTEM PRIVATE ${cgltf_SOURCE_DIR} ${stb_SOURCE_DIR})
target_compile_definitions(knit_io PRIVATE KNIT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
