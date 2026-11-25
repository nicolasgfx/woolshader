// The yarn's construction in the Knit step, mapped onto its settings keys.
#include "settings/yarn_look.h"

#include <algorithm>

#include "model/choices.h"
#include "settings/settings_mapping.h"

namespace knit {

std::string textureProfileId(int appearancePlies, PlyProfileModel style) {
    if (appearancePlies <= 1) return "none";
    return std::to_string(appearancePlies) + (style == PlyProfileModel::Circles ? "-circles" : "-swept");
}

YarnBodyRender yarnBodyRenderOf(const SceneSettings& s) { return s.yarnBody.render == "maps" ? YarnBodyRender::Maps : YarnBodyRender::Fibres; }

bool yarnOfFibres(YarnConstruction c) { return c == YarnConstruction::FibreStrands || c == YarnConstruction::FibreMaps; }

YarnConstruction yarnConstructionOf(const SceneSettings& s) {
    if (s.yarnBody.enabled) return yarnBodyRenderOf(s) == YarnBodyRender::Maps ? YarnConstruction::FibreMaps : YarnConstruction::FibreStrands;
    return s.plies > 1 ? YarnConstruction::PlyTubes : YarnConstruction::SingleTube;
}

YarnLook yarnLookOf(const SceneSettings& s, const YarnLookMemory& mem) {
    YarnLook l;
    l.construction = yarnConstructionOf(s);
    const bool plyTubes = l.construction == YarnConstruction::PlyTubes;
    l.tubePlies = plyTubes ? s.plies : std::clamp(mem.tubePlies, 2, 4);
    // The single tube's appearance: the settings' profile (a yarn of fibres keeps it, unused), or with
    // ply tubes the remembered one.
    int texPlies = 1;
    PlyProfileModel model = PlyProfileModel::Swept;
    if (!parseTextureProfile(plyTubes ? mem.textureProfile : s.textureProfile, texPlies, model, nullptr)) texPlies = 1;
    l.appearancePlies = texPlies > 1 ? texPlies : 0;
    l.style = model;
    return l;
}

void applyYarnLook(SceneSettings& s, const YarnLook& look, YarnLookMemory& mem) {
    const YarnConstruction was = yarnConstructionOf(s);
    const YarnConstruction next = look.construction;
    // Ply tubes clear the single tube's maps and flattening: they are remembered for the way back. A yarn
    // of fibres keeps the profile in the settings (unused while it is on), so switching to it and back
    // changes no knitting key and the knit stays current.
    if (was != YarnConstruction::PlyTubes && next == YarnConstruction::PlyTubes) {
        mem.textureProfile = s.textureProfile;
        mem.flatten = s.flatten;
        mem.flattenInflate = s.flattenInflate;
    }
    if (was == YarnConstruction::PlyTubes && next != YarnConstruction::PlyTubes) {
        mem.tubePlies = s.plies;
        s.flatten = mem.flatten;
        s.flattenInflate = mem.flattenInflate;
    }
    switch (next) {
    case YarnConstruction::PlyTubes:
        s.yarnBody.enabled = false;
        s.plies = std::clamp(look.tubePlies, 2, 4);
        mem.tubePlies = s.plies;
        s.textureProfile = "none";
        s.flatten = false;
        return;
    case YarnConstruction::FibreStrands:
    case YarnConstruction::FibreMaps:
        s.yarnBody.enabled = true;
        s.yarnBody.render = next == YarnConstruction::FibreMaps ? "maps" : "fibres";
        s.plies = 1;
        s.textureProfile = textureProfileId(look.appearancePlies, look.style);
        return;
    case YarnConstruction::SingleTube:
        s.yarnBody.enabled = false;
        s.plies = 1;
        s.textureProfile = textureProfileId(look.appearancePlies, look.style);
        return;
    }
}

std::string yarnLookConflict(const SceneSettings& s) {
    if (s.yarnBody.enabled && s.plies > 1)
        return "yarn_body.enabled with strand.plies = " + std::to_string(s.plies) +
               " is refused by the export: the yarn body makes its own plies (yarn_body.plies) on the single drawn yarn";
    // A yarn of fibres keeps the single tube's profile for the way back: not a conflict.
    if (s.textureProfile == "none" || s.yarnBody.enabled || s.plies <= 1) return "";
    return "strand.texture_profile = " + s.textureProfile + " is ignored: with strand.plies = " + std::to_string(s.plies) +
           " the plies are separate tubes (real geometry), and the twisted-ply maps apply to a single tube only";
}

} // namespace knit
