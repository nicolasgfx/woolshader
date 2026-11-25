#pragma once
// A part's yarn as a tube mesh, built on the device: the one place a drawn tube is meshed, for the
// viewer's display, the export and the gap coverage. Each caller turns the mesh into its own form.
#include <string>

#include "fibres/fibre.h"
#include "knitting/tube/yarn_surface.h"
#include "model/params.h"
#include "model/yarn_path.h"

namespace knit {

// How a yarn tube is meshed. The defaults are the export's tube (fibres/fibre.h, where the fibres'
// root depth uses the same polygons).
struct TubeStyle {
    double modelRing = kExportModelRing;  // surface model rings [mm]
    double ring = kExportRing;            // tube rings [mm]
    int sides = kExportSides;
    bool compression = false;             // elliptical rings
    double maxCompression = kDefaultMaxCompression;  // their flattening limit (SurfaceParams)
    double uvRepeat = 2.5;                // [mm] of yarn per u
    int uvWrap = 0;                       // TubeParams::uvWrap
    double uvPhaseJitter = 0, uvPitchJitter = 0;  // TubeParams::uvPhaseJitter, uvPitchJitter
};
// The surface model's parameters of a style (rings, compression, flattening limit): what both builders build
// from, and what surfaceBuildDigest (knitting/tube/param_text.h) keys a bake by.
SurfaceParams surfaceParamsOf(const TubeStyle& style);
// buildTubeMesh of buildSurface on the device; false with the device's reason.
bool buildYarnTube(const CurveSet& cs, const PatternParams& pp, const TubeStyle& style, TubeMesh& out, std::string* err);
// The surface model buildYarnTube meshes for the same style (its rings, compression and thickness), on
// the device, downloaded: what the yarn body's fibres follow (fibres/yarn_body.h).
bool buildYarnSurface(const CurveSet& cs, const PatternParams& pp, const TubeStyle& style, SurfaceModel& out, std::string* err);
bool buildYarnSurface(const CurveSet& cs, const PatternParams& pp, const SurfaceParams& sp, SurfaceModel& out, std::string* err);

} // namespace knit
