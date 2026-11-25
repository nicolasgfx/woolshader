# The viewer (ARCHITECTURE.md): its state, a library without ImGui or GLFW on its include path, so the
# compiler holds that layer, and the program with the panels and the 3D view. A leaf: every header is
# private, and the dependencies too (import's types come through jobs). The program sees the state's private headers and dependencies
# (knit_component_private_access). GL 2.1 compatibility profile + ImGui OpenGL2 backend; no GL loader.
# The About window's texts (about_texts.h): the licence files, word for word, written into a C++ file when the build is configured
# (cmake/embed_texts.cmake), and the commit and CUDA version, written at every build (cmake/build_identity.cmake).
include(${CMAKE_SOURCE_DIR}/cmake/embed_texts.cmake)
set(KNIT_ABOUT_DIR ${CMAKE_BINARY_DIR}/generated/viewer)
knit_embed_texts(OUTPUT ${KNIT_ABOUT_DIR}/about_texts.cpp INCLUDE viewer/about_texts.h TEXTS
    kWoolShaderLicenceText ${CMAKE_SOURCE_DIR}/LICENSE
    kImguiLicenceText      ${imgui_SOURCE_DIR}/LICENSE.txt
    kGlfwLicenceText       ${glfw_SOURCE_DIR}/LICENSE.md
    kCgltfLicenceText      ${cgltf_SOURCE_DIR}/LICENSE
    kStbLicenceText        ${stb_SOURCE_DIR}/LICENSE
    kOpenVdbLicenceText    ${openvdb_SOURCE_DIR}/LICENSE
    kTbbLicenceText        ${tbb_SOURCE_DIR}/LICENSE.txt
    kEigenLicenceNote      ${slim_eigen_SOURCE_DIR}/COPYING.README
    kMplLicenceText        ${slim_eigen_SOURCE_DIR}/COPYING.MPL2
    kRimLicenceText        ${rim_SOURCE_DIR}/LICENSE
    kHighsLicenceText      ${highs_SOURCE_DIR}/LICENSE.txt
    kFilereaderlpLicenceText ${highs_SOURCE_DIR}/extern/filereaderlp/LICENSE
    kInstantMeshesLicenceText ${CMAKE_SOURCE_DIR}/third_party/rim/INSTANT_MESHES_LICENSE.txt
    kMollerTriTriLicenceText ${CMAKE_SOURCE_DIR}/third_party/rim/MOLLER_TRITRI_LICENSE.txt
    kPcg32Notice           ${pcg32_SOURCE_DIR}/pcg32.h
    kSlimChangesText       ${CMAKE_SOURCE_DIR}/third_party/slim/CHANGES.md)
add_custom_target(knit_build_identity
    COMMAND ${CMAKE_COMMAND} -DSOURCE_DIR=${CMAKE_SOURCE_DIR} -DOUTPUT=${KNIT_ABOUT_DIR}/build_identity.cpp
            -DINCLUDE=viewer/about_texts.h -DCUDA_VERSION=${CUDAToolkit_VERSION} -P ${CMAKE_SOURCE_DIR}/cmake/build_identity.cmake
    BYPRODUCTS ${KNIT_ABOUT_DIR}/build_identity.cpp
    COMMENT "The commit for the About window"
    VERBATIM)

knit_component(viewer TARGET knit_viewer_state
    SOURCES
        about.cpp camera.cpp marking_flow.cpp presentation.cpp view_content.cpp viewer_prefs.cpp viewer_script.cpp workflow.cpp
        ${KNIT_ABOUT_DIR}/about_texts.cpp ${KNIT_ABOUT_DIR}/build_identity.cpp
    PRIVATE_HEADERS
        about.h about_texts.h camera.h fibre_lines.h fur_discs.h gl_functions.h interface.h marking_flow.h mesh_lighting.h native_dialogs.h overlay_lines.h renderer.h
        viewer_prefs.h view_content.h viewer_script.h workflow.h
    PRIVATE_DEPS model io settings export jobs
)
add_dependencies(knit_viewer_state knit_build_identity)
# The Render preview's renderer: a photon_path_tracer ppt.exe with `ppt serve`, a separate program. Its default is
# KNIT_PPT_EXE from the environment at configure time ("" without it: no Render preview until one is named);
# KNIT_PPT_EXE in the environment also overrides it at run time.
set(KNIT_PREVIEW_PPT_EXE "$ENV{KNIT_PPT_EXE}" CACHE FILEPATH "ppt.exe with `serve` for the Render preview")
target_compile_definitions(knit_viewer_state PRIVATE KNIT_DEFAULT_PPT_EXE="${KNIT_PREVIEW_PPT_EXE}")

add_library(imgui_lib STATIC
    ${imgui_SOURCE_DIR}/imgui.cpp
    ${imgui_SOURCE_DIR}/imgui_draw.cpp
    ${imgui_SOURCE_DIR}/imgui_tables.cpp
    ${imgui_SOURCE_DIR}/imgui_widgets.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl2.cpp)
target_include_directories(imgui_lib PUBLIC ${imgui_SOURCE_DIR} ${imgui_SOURCE_DIR}/backends)
target_link_libraries(imgui_lib PUBLIC glfw)
if(MSVC)
    target_compile_options(imgui_lib PRIVATE /W0)
endif()

add_executable(knit_viewer
    main.cpp interface.cpp about_window.cpp widgets.cpp panel_import.cpp panel_mark.cpp panel_knit.cpp panel_fibres.cpp panel_fur.cpp panel_yarn_fibres.cpp
    panel_export.cpp
    native_dialogs.cpp fibre_lines.cpp fur_discs.cpp mesh_lighting.cpp overlay_lines.cpp renderer.cpp)
knit_component_private_access(knit_viewer viewer)
target_link_libraries(knit_viewer PRIVATE imgui_lib glfw OpenGL::GL)
if(WIN32)
    target_link_libraries(knit_viewer PRIVATE ole32 shell32 uuid)  # file / folder dialogs (native_dialogs.cpp)
endif()
set_target_properties(knit_viewer PROPERTIES CXX_SCAN_FOR_MODULES OFF)
knit_warnings(knit_viewer)
# The UI font (widgets.cpp), next to the executable.
file(COPY ${imgui_SOURCE_DIR}/misc/fonts/Roboto-Medium.ttf DESTINATION ${CMAKE_BINARY_DIR}/fonts)
target_compile_definitions(knit_viewer PRIVATE KNIT_FONT_DIR="${CMAKE_BINARY_DIR}/fonts")
