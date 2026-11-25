# Writing reviewed results: the format table, exportScene and one writer per format
# (ARCHITECTURE.md). The door, export.h, is the only public header; it needs io's publication
# policy. Everything the formats stand on (the neutral scene, the collection from a scene) is private.
knit_component(export
    SOURCES
        collect_materials.cpp collect_scene.cpp export.cpp export_scene.cpp format_blend.cpp format_obj.cpp
        format_ppt.cpp
    PUBLIC_HEADERS
        export.h
    PRIVATE_HEADERS
        collect_scene.h export_scene.h
    PUBLIC_DEPS io
    PRIVATE_DEPS model fibres settings scene
)
