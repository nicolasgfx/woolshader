#pragma once
// The yarn's construction in the Knit step: how each knitted yarn is built for the export and the render
// preview, as the user chooses it on four tiles, mapped onto the settings keys strand.plies /
// strand.texture_profile / strand.flatten / yarn_body.enabled / yarn_body.render.
//   Single tube (the main workflow): one validated tube per yarn. Its twisted-ply appearance is a
//     set of surface maps (height, normal, AO, cavity) of 2..4 plies baked around that tube
//     (ply_profile.h); "off" is the plain tube. Contact flattening may apply.
//   Ply tubes: 2..4 real ply tubes inside the drawn yarn (plies.h). The maps do not apply (they would
//     depict the same plies again) and neither does flattening (the plies are round).
//   Fibre strands: the yarn body (fibres/yarn_body.h), 1..3 plies of long fibres on the single tube's
//     surface model, in place of the tube in every export. The maps do not apply (the fibres are the
//     plies; the export attaches none), but strand.texture_profile keeps its value, a knitting key: so
//     switching to fibres and back knits nothing again. Flattening applies (the fibres follow the
//     flattened surface). The viewport shows the tube.
//   Fibre maps: the same yarn body, rendered as one tube per yarn under the maps baked from its fibres
//     (fibres/yarn_bake.h): yarn_body.render = maps. Everything else as Fibre strands.
// Switching to ply tubes and back restores the single tube's appearance and flattening settings
// (YarnLookMemory), so a look is never lost by trying another one. The yarn body's own settings
// (yarn_body.*, Export-stage keys) are edited in the Fibres step (viewer/panel_yarn_fibres.cpp).
#include <string>

#include "model/choices.h"
#include "settings/scene_settings.h"

namespace knit {

enum class YarnConstruction { SingleTube, PlyTubes, FibreStrands, FibreMaps };
// How a yarn made of fibres is rendered (yarn_body.render): its fibres as strands, or one tube per yarn under
// the maps baked from them (fibres/yarn_bake.h). Both are equal choices of the same yarn body.
enum class YarnBodyRender { Fibres, Maps };
YarnBodyRender yarnBodyRenderOf(const SceneSettings& s);
// The yarn is made of fibres, rendered either way.
bool yarnOfFibres(YarnConstruction c);

struct YarnLook {
    YarnConstruction construction = YarnConstruction::SingleTube;
    int tubePlies = 3;               // ply tubes: 2..4
    int appearancePlies = 3;         // single tube: 0 (off) or 2..4 plies in the maps
    PlyProfileModel style = PlyProfileModel::Swept;
    bool operator==(const YarnLook&) const = default;
};

// The single tube's settings while another construction is chosen (and the tube count while ply tubes
// are not), restored when switching back.
struct YarnLookMemory {
    std::string textureProfile = "3-swept";
    bool flatten = false;
    double flattenInflate = 0.15;
    int tubePlies = 3;
};

// How the settings build the yarn: a yarn of fibres when yarn_body.enabled (strands or maps by
// yarn_body.render), ply tubes when strand.plies > 1.
YarnConstruction yarnConstructionOf(const SceneSettings& s);
// The look the settings describe; the values of the constructions not in use come from the memory.
YarnLook yarnLookOf(const SceneSettings& s, const YarnLookMemory& mem);
// Writes the look into the settings (and the memory: the values of the mode left).
void applyYarnLook(SceneSettings& s, const YarnLook& look, YarnLookMemory& mem);
// "none" or "<n>-swept" / "<n>-circles".
std::string textureProfileId(int appearancePlies, PlyProfileModel style);
// A settings combination that drops a key (ply tubes with a texture profile) or that the export refuses
// (the yarn body with ply tubes): what and why, for the log; "" when there is none.
std::string yarnLookConflict(const SceneSettings& s);

} // namespace knit
