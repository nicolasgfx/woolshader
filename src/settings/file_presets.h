#pragma once
// Presets the user edits by hand (presets/<kind>/*.preset in the repository): the file format and its
// groups. One file per preset:
//
//   knit-preset v1
//   # comments
//   id = standard              (unique within the kind; the settings files name it)
//   label = Standard wool      (what the viewer shows)
//   order = 1                  (position in the list; the lowest is the viewer's default)
//   note = ...                 (one line: what it is for)
//   strays = loose             (Fibres presets only: the stray preset it comes with)
//   fibres.random_fraction = 0.2   (settings keys, scene_settings.h, with their values)
//
// A preset sets exactly the keys its kind owns, and its values are validated. The viewer reads the
// files at start (edits apply after a restart); a file that fails is left out and its error reported.
//
//   Fibres (presets/fibres): the two populations together, the fibres along the yarn ("directed",
//     oriented) and the stray ones ("random"), with the stray share and the segments per fibre. The
//     density is not part of them. The lowest `order` is the viewer's default for new objects. The
//     label (fibres.preset) only names the values: it is not part of any result key. Editing a value
//     makes the settings Custom. Each names its stray preset (strays = <id>), applied with it.
//   Strays (presets/strays): how far the stray fibres deviate from the yarn, in direction and in
//     length: fibres.random.spread_deg, max_angle_deg and length_spread. The Fibres group does not own
//     these three, so choosing a stray preset keeps the wool preset's name, and the reverse. No label
//     is stored: the matching preset is found by value.
//   Density (presets/density): fibres.density alone, within the Density control's 1..1000. No label is
//     stored: the matching preset is found by value, Custom when none matches.
//   Strand mesh (presets/strand_mesh): the yarn tube's sides and ring spacing and the plies' (explicit
//     plies only), within tubeSettingsValid's bounds.
//   Hand knit (presets/hand_knit): how unevenly the stitches are knitted (knit.hand.*, but its seed).
//   Fibre islands (presets/fibre_islands): the patches of different fuzz (fibres.islands.*).
//   In these two an all-zero preset is the Off of the list: there is no switch.
//   Yarn fibres (presets/yarn_fibres): the look of a yarn made of fibres (Knit step, Construction "Made of
//     fibres", "Look"), for both renderings: the body's structure (yarn_body.plies, fibres,
//     fibre_diameter_um, shell_depth, fibre_angle_deg, ply_angle_deg, irregularity), the baked maps' material
//     (yarn_body.map_cavity, map_roughness, map_specular, map_anisotropy, map_material) and the real fibres' hair shading
//     (material.yarn_white.fibre_*). No label is stored: the matching look is found by value. The interior
//     fibres, the core proxy, the seed, the texel and the rendering stay yours.
#include <string>
#include <utility>
#include <vector>

#include "settings/scene_settings.h"

namespace knit {

// A line of a "key = value" text file: settings and presets share this grammar.
struct KeyValueLine {
    int line = 0;
    std::string key, value;
};
// Blank lines and # comments are skipped; the first other line must equal versionLine.
bool parseKeyValueText(const std::string& text, const std::string& versionLine, std::vector<KeyValueLine>& out, std::string* err);
std::string trimmed(const std::string& text);  // no blanks, tabs or \r at either end

struct PresetFile {
    std::string path, id, label, note;
    std::string strays;  // the `strays = <id>` line ("" without one)
    int order = 0;
    std::vector<std::pair<std::string, std::string>> values;  // settings keys, in file order
};

// Parses one file's text; false with the reason (and the line) otherwise.
bool parsePresetFile(const std::string& text, PresetFile& out, std::string* err);
// The folder of a kind (presets/<kind>, from the working directory or the source tree); "" if absent.
std::string presetFolder(const std::string& kind);
// Every *.preset file of the folder, parsed, sorted by order then id; failures go to `errors`
// ("<file>: <reason>"). Duplicate ids: the later file is refused.
std::vector<PresetFile> readPresetFiles(const std::string& folder, std::vector<std::string>* errors);

// The folders of presets.
enum class PresetGroup { Fibres, Density, StrandMesh, HandKnit, FibreIslands, Strays, YarnFibres };

// One preset file: what the viewer shows, and the settings it sets.
struct FilePreset {
    std::string id, label, note, path;  // id: unique in its group; path: the file it came from
    SceneSettings values;               // every key of presetKeys as the file sets it
    std::string strays;                 // Fibres: the id of its stray preset (applied with it); else ""
};

// The settings keys a preset of the group sets: all of them, each exactly once.
const std::vector<std::string>& presetKeys(PresetGroup group);
// One preset from a parsed file; false with the reason.
bool presetFromFile(PresetGroup group, const PresetFile& f, FilePreset& out, std::string* err);
// The group's presets, read once per process (sorted by order); the files left out.
const std::vector<FilePreset>& presets(PresetGroup group);
const std::vector<std::string>& presetErrors(PresetGroup group);
// The group's preset of that id; nullptr if no file has it.
const FilePreset* presetWithId(PresetGroup group, const std::string& id);
// The preset whose values the settings hold (nullptr: Custom).
const FilePreset* matchingFilePreset(PresetGroup group, const SceneSettings& s);
// Sets every key the group owns to the preset's values; a Fibres preset then applies its stray preset.
void applyPreset(PresetGroup group, SceneSettings& s, const FilePreset& p);

} // namespace knit
