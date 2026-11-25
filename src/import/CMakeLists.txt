# Turning a user's mesh into a knittable one (ARCHITECTURE.md).
# Public: the analysis, the preparation, the pipeline and the store. The rebuild's stages (wrap_*) are private.

# The wrap's flattening, on the solver's pinned Eigen (cmake/dependencies.cmake): a sub-library of this
# component because it compiles with that Eigen. It sees this component's headers and calls nothing back;
# no Eigen type appears in a header.
add_library(knit_flatten STATIC wrap_flatten.cpp wrap_quad_split.cpp)
target_include_directories(knit_flatten PRIVATE ${KNIT_COMPONENT_ROOT}/import/public ${KNIT_COMPONENT_ROOT}/import/private)
target_link_libraries(knit_flatten PRIVATE knit_model knit_io knit_slim)  # import's public dependencies
set_target_properties(knit_flatten PROPERTIES CXX_SCAN_FOR_MODULES OFF)
knit_warnings(knit_flatten)

knit_component(import
    SOURCES
        import_analysis.cpp import_pipeline.cpp import_store.cpp prepare_import.cpp repair_rules.cpp
        wrap_bake.cpp wrap_charts.cpp wrap_grid.cpp wrap_nets.cpp wrap_pack.cpp
        wrap_remesh.cpp wrap_surface.cpp wrap_tidy.cpp
    PUBLIC_HEADERS
        import_analysis.h import_pipeline.h import_store.h prepare_import.h
    PRIVATE_HEADERS
        repair_rules.h wrap_stages.h wrap_surface.h
    PUBLIC_DEPS model io
    PRIVATE_LIBS knit_flatten knit_vdb
)
