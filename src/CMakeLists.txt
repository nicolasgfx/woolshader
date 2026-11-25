# The components (cmake/components.cmake; the graph and its rules: ARCHITECTURE.md here), lowest first: each folder declares itself,
# its public and private headers and its dependencies in its CMakeLists.txt. Headers sit next to their
# sources; includes name the component (`#include "knitting/relax/settle.h"`).
add_subdirectory(model)
add_subdirectory(knitting)
add_subdirectory(stitchmesh)
add_subdirectory(fibres)
add_subdirectory(gpu)
add_subdirectory(io)
add_subdirectory(settings)
add_subdirectory(import)
add_subdirectory(scene)
add_subdirectory(export)
add_subdirectory(jobs)
add_subdirectory(cli)

# Viewer (optional).
if(KNIT_BUILD_VIEWER)
    add_subdirectory(viewer)
endif()
