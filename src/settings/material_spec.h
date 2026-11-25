#pragma once
// A renderer material as the user sets it (SceneSettings holds several), apart from the settings
// document, so the export's formats compile against the material alone.
#include "model/colour.h"
#include "model/vec.h"

namespace knit {

// The yarn's sheen takes the yarn's own colour (the wool look, 2026-11-01). Yarn materials only:
// MaterialSpec's own default also reaches materials that are no yarn.
constexpr double kYarnSheenTint = 1.0;

// Renderer material (photon_path_tracer layer types; export/format_ppt.cpp writes them).
struct MaterialSpec {
    enum class Model { DiffuseSheen, Fibre } model = Model::DiffuseSheen;  // Fibre: strands only (hair BSDF)
    Vec3 color = kDefaultYarnColour;
    // Sheen layer weight (0: no sheen layer).
    double sheen = 1.0, sheenTint = 0.5, sheenRoughness = 0.5;
    bool clearcoat = false;                                                 // dielectric coat
    double clearcoatWeight = 0.3, clearcoatIor = 1.5, clearcoatRoughness = 0.3;
    double fibreLongRoughness = 0.3, fibreAzimRoughness = 0.3, fibreIor = 1.55, fibreCuticleDeg = 2.0;
    // A fabric's far side (the renderer's TranslucencyParams): the share of the body's light that
    // leaves through the far side, scattered, instead of being reflected; one budget. 0: opaque.
    // Written for a fabric (sheen > 0) material only; a hair material has its own transmission.
    double translucency = 0;
    bool operator==(const MaterialSpec&) const = default;
};

} // namespace knit
