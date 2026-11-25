#pragma once
// From settings to the parameters of knitting and fibres: a part's hand, its final surface, the
// filling depth, and the fibre and fur parameters. The settings-only mapping is
// settings/settings_mapping.h.
#include <string>

#include "fibres/fibre.h"
#include "fibres/fur_look.h"
#include "fibres/yarn_bake.h"
#include "fibres/yarn_body.h"
#include "knitting/stitch/hand_tension.h"
#include "knitting/tube/yarn_surface.h"
#include "model/params.h"
#include "settings/scene_settings.h"

namespace knit {

// The fibre parameters of the settings (pieceFibres() plus the settings' fields; random-population
// overrides only where they differ from the directed values). False with a reason for inconsistent
// settings (fibreSettingsValid).
bool fibreParamsFromSettings(const SceneSettings& s, FibreParams& fp, std::string* err);
// The Fur workflow's look of the settings (pile.*): cm^2 to mm^2, gravity [m/s^2] to the settle's scale. The
// bounds are pileSettingsValid's (settings_mapping.h), checked by the Fur stage.
FurLook furLookFromSettings(const SceneSettings& s);
// A part's hand, from the settings and the part's name.
HandRequest handRequest(const SceneSettings& s, const std::string& partName);
// The final surface of a part knitted at diameter D.
FinalSurfaceSpec finalSurface(const SceneSettings& s, double D);
// The yarn body's parameters (yarn_body.*). Refused: explicit ply tubes (strand.plies > 1), whose round
// plies keep their own path. The geometric checks need the yarn's radius: yarnBodyParamsValid.
bool yarnBodyParams(const SceneSettings& s, YarnBodyParams& p, std::string* err);
// The bake's parameters (yarn_bake.*): the texel from the settings (the material values are applied when a
// package attaches the maps, scene/appearance/baked_maps.h).
BakeParams bakeParams(const SceneSettings& s);
// The surface model the export's tubes and the yarn body follow: the model's rings at the export's spacing,
// flattened (compressed to the flattening limit, 0.025 mm rings) where the settings say so.
SurfaceParams yarnSurfaceParams(const SceneSettings& s);
// How far below the source surface (along -N) the fabric's inner extent lies, plus the filling
// margin: A D + r_R + c + filling.margin (r_R of the final surface at D). The knit is centred on the
// surface (controls at +-A D), so anything drawn behind the yarn sits at least this deep: the
// export's filling (scene/appearance/filling.h) and the viewer's stand-ins of parts still being knitted.
double fillingOffset(const SceneSettings& s, double D);

// The flat swatch at the settings' gauge: the knit (contact) diameter k D (yarnFit), w, h, A and c; the
// pattern and the patch size are the caller's.
PatternParams flatPatternParams(const SceneSettings& s);
// The flat reference verdict of pattern p at the settings' gauge (flatPatternReference,
// knitting/checks/checks.h, at flatPatternParams).
PatternReference flatPatternReference(const SceneSettings& s, StitchPattern p);

} // namespace knit
