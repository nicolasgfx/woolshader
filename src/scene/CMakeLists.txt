# A scene from its source through Knit and Fibres to what the export reads (ARCHITECTURE.md). Public:
# the records, the stage entry points, and what export, jobs and knit_check read; a part's own knitting
# and fibre jobs, the solids and the caches' files are private. The GPU, the
# import store and the stitch mesh (the guided route) are private dependencies.
knit_component(scene
    SOURCES
        result_keys.cpp run_log.cpp scene.cpp stage_params.cpp
        appearance/baked_maps.cpp appearance/filling.cpp appearance/ply_maps.cpp
        appearance/rendered_part.cpp appearance/stitch_colours.cpp appearance/yarn_tube.cpp
        fibres/fibre_files.cpp fibres/fibre_stage.cpp fibres/fur_stage.cpp fibres/part_fibre_job.cpp fibres/scene_fibres.cpp
        fibres/reference_samples.cpp
        knit/gap_patches.cpp knit/guided_part.cpp knit/knit_cache.cpp knit/knit_report.cpp knit/knit_stage.cpp knit/mesh_part.cpp
        knit/saved_knit.cpp knit/spiral_cache.cpp knit/spiral_part.cpp
        source/guided_source.cpp source/scene_parts.cpp source/scene_solids.cpp source/sewn_pieces.cpp
    PUBLIC_HEADERS
        result_keys.h run_log.h scene.h stage_params.h
        appearance/baked_maps.h appearance/filling.h appearance/ply_maps.h appearance/rendered_part.h appearance/stitch_colours.h appearance/yarn_tube.h
        fibres/fibre_files.h fibres/fibre_stage.h fibres/fur_stage.h fibres/reference_samples.h fibres/scene_fibres.h
        knit/knit_cache.h knit/knit_report.h knit/knit_stage.h knit/saved_knit.h
        source/guided_source.h source/scene_parts.h source/sewn_pieces.h
    PRIVATE_HEADERS
        fibres/part_fibre_job.h
        knit/gap_patches.h knit/knit_run.h knit/spiral_cache.h knit/spiral_part.h
        source/scene_solids.h
    PUBLIC_DEPS model knitting fibres io settings
    PRIVATE_DEPS gpu import stitchmesh
)
