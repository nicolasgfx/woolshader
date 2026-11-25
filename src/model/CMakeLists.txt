# What a knit is: plain values, standard C++ only (ARCHITECTURE.md). Every header is public.
knit_component(model
    SOURCES
        bspline.cpp compute_caps.cpp log.cpp mesh.cpp mesh_check.cpp mesh_query.cpp noise.cpp parallel.cpp spatial_hash.cpp stitch_mesh_data.cpp surface_marking.cpp
        surface_path.cpp surface_pieces.cpp topology.cpp yarn_path.cpp
    PUBLIC_HEADERS
        bspline.h cancel.h choices.h colour.h compute_caps.h coordinates.h env.h hash.h log.h mesh.h mesh_check.h
        marking_review.h mesh_query.h noise.h parallel.h params.h spatial_hash.h status.h stitch_mesh_data.h surface_marking.h surface_path.h
        surface_pieces.h topology.h vec.h yarn_path.h
)
