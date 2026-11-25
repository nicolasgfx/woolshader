#pragma once
// The fibres of a fur field (fur_field.h): each disc's rest-shape curve, then its static settle (fur.h).
//
// The curve of a loop (radius r > 0), from the root along d: the tangent turns in the plane (d, c) with a
// curvature rising linearly from 0 to 1/r over the arc (pi r / 2) rise, then constant; across the plane,
// along w = d x c, the curve climbs kDiscClimbDiameters root diameters per turn, so a loop past one turn
// coils beside itself instead of retracing. Its requested visible length is
//     L = max(4 kMinSegmentRadii r0, turns (2 pi r + (pi r / 4) rise)),
// one nominal loop (a full circle plus the quarter-turn ramp's shortfall) times `turns`; the minimum keeps
// a short, thick fibre within the strand rules (it counts as "clamped"). L is the polyline's length, not a
// guaranteed turning angle: the ramp and the equal-chord vertices change the actual shape (measure it).
// A straight fibre (radius 0) runs along d for straightLength.
//
// Vertices: a buried root one strand-rule segment (kMinSegmentRadii r0) below the surface along -normal, the
// root, then n visible segments of equal chord L / n, each vertex on a fine curve (kCurlStepsPerSegment
// samples per segment) exactly one chord from the previous one. n (furSegmentCount) follows the strand
// rules: about 0.5 mm per segment, at most an eighth of a loop's circle, and every segment at least
// kDiscMinSegmentRadii root radii for a loop (the renderer's joint rule, segment >= one diameter, with a
// margin; the exporter checks kDiscMinSegOverR) or kMinSegmentRadii for a straight fibre, then more segments
// while the thin tip needs them. On a tight, thick loop the joint rule wins: fewer segments than eight per
// turn (the baseline's tightest loops hold 3-4).
#include <vector>

#include "fibres/fur.h"
#include "fibres/fur_field.h"
#include "fibres/fur_look.h"
#include "model/mesh_query.h"
#include "model/parallel.h"

namespace knit {

constexpr double kDiscClimbDiameters = 1.1;   // a loop's climb across its disc per turn [root diameters]
constexpr double kDiscMinSegmentRadii = 2.3;  // a loop's shortest segment [root radii]: the renderer's 2 r, with a margin
constexpr double kDiscMinSegOverR = 2.2;      // the exporter's rule for the Fur workflow (strands_file.h StrandRules)
constexpr double kFurSegmentTarget = 0.5;     // [mm] a segment's length where no rule asks for less

// The requested visible length [mm] (see above).
double furFibreLength(const FurDisc& d);
// Whether the minimum-length rule raised a loop's length.
bool furLengthClamped(const FurDisc& d);
// The visible segments (the polyline has one more, the buried root).
int furSegmentCount(const FurDisc& d);
// The rest shape: buried root, root, then the visible vertices.
Hair furRestShape(const FurDisc& d);

// Every disc's fibre settled (in parallel; the result does not depend on the thread count), in the field's
// order; the stuck ones dropped. rep: each layer's stuck, grown and segment counts, the largest residual
// and the settle time are added.
std::vector<Hair> growLoops(const FurField& field, const MeshQuery& body, const FurSettle& settle, const ExecPolicy& exec,
                            FurReport* rep);

} // namespace knit
