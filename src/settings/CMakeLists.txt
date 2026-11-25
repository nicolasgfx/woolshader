# What the user chose: the settings document, its field table and text form, the presets, the scene
# registry, Draft and Final, the yarn look and the gauge fit, in the settings' own terms
# (ARCHITECTURE.md). Every header is public. The fibre model's defaults bound the fibre settings
# (fibreSettingsValid), so fibres is a private dependency; the mapping to knitting's and fibres' types
# is scene/stage_params.h.
knit_component(settings
    SOURCES
        file_presets.cpp scene_registry.cpp scene_settings.cpp settings_mapping.cpp settings_profiles.cpp
        thickness_presets.cpp yarn_look.cpp
    PUBLIC_HEADERS
        file_presets.h material_spec.h scene_registry.h scene_settings.h settings_mapping.h settings_profiles.h
        thickness_presets.h yarn_look.h
    PUBLIC_DEPS model io
    PRIVATE_DEPS fibres
)
