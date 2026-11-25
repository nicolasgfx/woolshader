#pragma once
// The Mark step's overlay in the 3D view (GL, no ImGui): the markings' paths, rings and stitch outlines as lines and
// the control points as dots, over the shaded meshes with the matrices currently loaded. What lies in front of the
// surface is drawn bright; what lies behind it is drawn faintly, so a marking going round a piece stays findable.
// Client-side arrays from the OverlayGeometry (tens of thousands of vertices at most): no buffer to own or release.
#include "viewer/marking_flow.h"

namespace viewer {

void drawOverlay(const OverlayGeometry& g);

} // namespace viewer
