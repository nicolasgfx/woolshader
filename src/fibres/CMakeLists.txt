# Fibres and fur, grown from the yarn surface and from meshes (ARCHITECTURE.md). The public headers name the yarn
# surface's types through pointers only, so knitting is a private dependency.
knit_component(fibres
    SOURCES
        fibre.cpp fibre_budget.cpp fibre_checks.cpp fibre_generator.cpp fibre_islands.cpp fibre_shape.cpp fibre_strands.cpp fur.cpp fur_field.cpp fur_loops.cpp
        strands_file.cpp yarn_bake.cpp yarn_body.cpp
    PUBLIC_HEADERS
        fibre.h fibre_budget.h fibre_checks.h fibre_generator.h fibre_islands.h fibre_record.h fibre_shape.h fibre_strands.h fur.h fur_field.h fur_look.h fur_loops.h
        strands_file.h yarn_bake.h yarn_body.h
    PUBLIC_DEPS model
    PRIVATE_DEPS knitting
)
