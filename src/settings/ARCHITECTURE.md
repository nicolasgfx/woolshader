# Settings: file tree and architecture

What the user chose, in the settings' own terms: the settings document and its text form, the preset
files, the table of scenes, the Draft and Final profiles, the yarn's look and the gauge fit. This file is
the way into the folder. It knows no scene and no knitting: `src/scene` turns these values into the
stages' parameters (`scene/stage_params.h`) and keys (`scene/result_keys.h`).

## File tree

```text
src/settings/
  ARCHITECTURE.md             this file
  CMakeLists.txt              the component knit_settings: every header public; depends on model and io,
                              and privately on fibres (the fibre model's defaults bound its settings)
  scene_settings.h, .cpp      SceneSettings, the field table (settingFields), settings text, diff
  material_spec.h             MaterialSpec: a renderer material as the user sets it (the export's formats
                              read this header alone)
  scene_registry.h, .cpp      the presets table, workflow rules, credits
  settings_mapping.h, .cpp    settings to numbers in their own terms: scene defaults, diameters to try,
                              gauge fit, the scene's metric mapping, the fibre settings' rules
  settings_profiles.h, .cpp   Draft and Final, opening settings, what this machine cannot compute
  file_presets.h, .cpp        the .preset format and its groups: `PresetGroup`, one list
  thickness_presets.h, .cpp   yarn and stitch presets, fit notes and warnings
  yarn_look.h, .cpp           the yarn's look mapped onto its settings keys
```

**Reading order:** `scene_settings.h`, `scene_registry.h`, `settings_mapping.h`, then the others as a
reference.

## Layers

An arrow means "includes".

```text
settings_mapping.h ───┐
settings_profiles.h ──┼──> scene_registry.h ──┐
thickness_presets.h ──┘                       ├──> scene_settings.h ──> material_spec.h
file_presets.h, yarn_look.h ──────────────────┘
```

The three on the left include `scene_settings.h` directly too.

- Only `src/model`, `src/io` and `src/fibres` are on this folder's include path, and fibres only for the
  fibre model's defaults in `settings_mapping.cpp` (`fibreSettingsValid`); the compiler refuses anything
  else (`CMakeLists.txt`, `src/ARCHITECTURE.md`).
- The viewer edits settings through these headers alone, so a change inside `src/scene` does not
  recompile the panels.

## What must stay true

| Rule | Held by |
|---|---|
| One settings table: every setting is one row of `settingFields()` with its stage. Reader, writer, diff and the viewer's bindings use it. | `scene_settings_roundtrip`, `scene_settings_refusals` |
| Unknown keys are refused; a removed key is retired (`kRetiredKeys`) and still loads. | `scene_settings_retired_keys`, `scene_settings_retired_knit_keys`, `scene_settings_retired_hair_keys` |
| A preset file is checked at load, with the same rules the stages apply (`fibreSettingsValid`, `tubeSettingsValid`). | `hand_presets_load_with_an_off`, `fibre_presets_assign_all_values` |
| The table's order is part of the stage keys (`scene/result_keys.cpp`). | `knit_stage_key_pinned`, `io_pin_stage_keys` |

## Where to change what

| To change | Touch |
|---|---|
| A setting | one row in `makeTable`, `scene_settings.cpp`. Its place in the table is part of the stage keys. |
| Retire a setting | delete its row, add it to `kRetiredKeys` with its last default; if it was in a default stage text, also to `kRetiredInputs` (`scene/result_keys.cpp`). |
| A new preset group | `file_presets.*`: a `PresetGroup` value, its keys, its messages, its rules. |
| A new scene | `scene_registry.cpp`, with its model's credit. |
