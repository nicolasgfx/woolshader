# The stockinette stitch-mesh experiment's stages before relaxation (ARCHITECTURE.md): the reference surface, Method A's
# time field, rings and guided rows, Method B's W18 pipeline (ported from dbf5f41: remesh, labels, directions, subdivision),
# the shared typed mesh with its independent check, the virtual shift, coverage and the construction fixtures. Every
# header is public; it depends on model and knitting publicly (yarn_settle.h takes the knitting solver, KnitNumerics, and returns its validation reports). The remesher (knit_rim) and the integer solver (HiGHS) stay private: no
# RIM, Eigen or HiGHS header reaches a user of this component.
knit_component(stitchmesh
    SOURCES
        border_placement.cpp construction_fixtures.cpp coverage.cpp directions.cpp guided_rows.cpp labels.cpp mark_proposals.cpp marked_build.cpp orientation_field.cpp polygon_mesh.cpp prepared_reference.cpp quick_review.cpp reference_markings.cpp
        quantiles.cpp reference_surface.cpp relax_mesh.cpp remesh_stage.cpp revolution_surface.cpp row_rings.cpp stage_report.cpp stitch_fixtures.cpp
        stitch_mesh.cpp surface_trim.cpp time_field.cpp typed_mesh.cpp virtual_shift.cpp guided_records.cpp yarn_curves.cpp yarn_embed.cpp yarn_settle.cpp
    PUBLIC_HEADERS
        border_placement.h construction_fixtures.h coverage.h directions.h guided_rows.h labels.h mark_proposals.h marked_build.h orientation_field.h polygon_mesh.h prepared_reference.h quantiles.h quick_review.h reference_markings.h
        reference_surface.h relax_mesh.h remesh_stage.h revolution_surface.h row_rings.h stage_report.h stitch_fixtures.h stitch_mesh.h
        surface_trim.h time_field.h typed_mesh.h virtual_shift.h guided_records.h yarn_curves.h yarn_embed.h yarn_settle.h
    PUBLIC_DEPS model knitting
    PRIVATE_LIBS knit_rim highs
)
