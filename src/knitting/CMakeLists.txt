# How a part is knitted and checked, one folder per stage (ARCHITECTURE.md). Every header is public: gpu,
# fibres, scene, the command line and the tests use the stage headers directly.
knit_component(knitting
    SOURCES
        checks/clearance.cpp checks/crossings.cpp checks/curvature.cpp checks/signature.cpp checks/validate.cpp
        knit/coverage.cpp knit/fixtures.cpp knit/mesh_knit.cpp knit/piece_knit.cpp
        knit/pieces.cpp
        layout/base_surface.cpp layout/chart.cpp layout/chart_lattice.cpp layout/stitch_layout.cpp
        layout/uv_atlas.cpp
        relax/knit_problem.cpp relax/settle.cpp
        stitch/hand_tension.cpp stitch/stitch_program.cpp stitch/stockinette.cpp
        tube/param_text.cpp tube/plies.cpp tube/ply_profile.cpp tube/yarn_surface.cpp
    PUBLIC_HEADERS
        numerics.h
        checks/checks.h checks/interval.h
        knit/coverage.h knit/fixtures.h knit/mesh_knit.h knit/piece_knit.h knit/pieces.h
        layout/base_surface.h layout/chart.h layout/chart_lattice.h layout/stitch_layout.h
        layout/uv_atlas.h
        relax/knit_problem.h relax/settle.h
        stitch/expected_signature.h stitch/hand_tension.h stitch/stitch_program.h stitch/stockinette.h
        tube/compression_fit.h tube/param_text.h tube/plies.h tube/ply_profile.h tube/yarn_surface.h
    PUBLIC_DEPS model
)
