# What the viewer and knit_check start: the four jobs (jobs.h), the queue that runs them, the reviewed results
# and the display data (built_scene.h). Public: the jobs and the queue, the results and the questions about
# them, the display data, the previews, the import session, the scene cache, the crash report. The display
# builders and the jobs' shared steps are private; so are scene, io, gpu, knitting and fibres, which the jobs
# run and read but never hand on.
knit_component(jobs
    SOURCES
        crash_report.cpp display_scene.cpp export_job.cpp fibre_job.cpp import_session.cpp job_steps.cpp knit_job.cpp marking_session.cpp
        fur_preview.cpp previews.cpp render_preview.cpp results.cpp scene_builder.cpp scene_cache.cpp strand_mesh_job.cpp
    PUBLIC_HEADERS
        built_scene.h crash_report.h fur_preview.h import_session.h jobs.h marking_session.h previews.h render_preview.h results.h scene_builder.h
        scene_cache.h
    PRIVATE_HEADERS
        display_scene.h job_steps.h
    PUBLIC_DEPS model settings export import io
    PRIVATE_DEPS knitting fibres gpu scene stitchmesh
)
