#pragma once
// The static settle of a fibre grown on a mesh, and the test that it does not stick in the body. The Fur
// workflow's fibres (fur_loops.h) are shaped in their rest pose first, then settled here once:
//
//   settleHair: a static shape under gravity, bend stiffness and body collision, solved as a fixed number of
//   damped time steps of
//     - Dynamic Follow-The-Leader (Mueller, Kim, Chentanez, "Fast Simulation of Inextensible Hair and Fur",
//       VRIPHYS 2012): predicted positions are projected one segment at a time onto the exact rest length
//       from the root outward; the velocity update subtracts a share s_damp of the next vertex's correction
//       (their eq. for v_i);
//     - a local shape constraint after Han and Harada ("Real-time Hair Simulation with Efficient Hair Style
//       Preservation", VRIPHYS 2012), simplified: each vertex is pulled toward its rest direction relative
//       to the previous segment, the rest direction being carried by the minimal rotation between the
//       previous segment's rest and current directions (no twist frame), with stiffness k_i from the root's
//       to the tip's value;
//     - collision with the body: a vertex closer than its radius plus an offset to the body
//       (MeshQuery::closest, the LOCAL pseudo-normal side - never a global inside test, which is undefined on
//       open meshes) is moved to the surface point plus that distance along the face normal.
//   The two root vertices are fixed.
//
//   hairInsideBody: the final pass. A fibre with a vertex inside the body, or a centreline segment crossing
//   its surface (MeshQuery::segmentCrosses, exact Moeller-Trumbore), is "stuck": rejected and counted, never
//   kept silently. A thick fibre's surface dipping into the body is not detected (its centreline stays out).
//
// Gravity tied to length is a styling rule, not a calibrated one: a fibre of visible length L settles under
// gravity * min(1, (L / gravityFullLength)^3), motivated by a cantilever's tip angle under its own weight
// (q L^3 / 6EI); off at 0. Limits (stated, not hidden): a plausible static groom, not a mechanical model -
// the droop depends on the stiffness and the step count together; fibres do not collide with each other;
// there is no animation. Units: library mm; gravity 9810 mm/s^2 times StrandSettleParams::gravity.
#include <vector>

#include "model/mesh_query.h"
#include "model/vec.h"

namespace knit {

// The static settle of one fibre (see above).
struct StrandSettleParams {
    int steps = 60;
    double gravity = 1.0;             // gravity scale (1: 9810 mm/s^2)
    double timeStep = 1.0 / 240;      // [s]
    double damping = 0.9;             // DFTL s_damp
    double velocityKeep = 0.8;        // global damping per step
    double stiffness = 0.6;           // shape constraint at the root, in [0, 1]
    double tipStiffness = 0.2;        // ... at the tip (linear in between)
};

// A grown fibre: points[0] is the buried root, points[1] the surface point, the rest visible.
struct Hair {
    int tri = -1, index = 0;          // triangle and candidate number of its root (with the seed: the identity)
    std::vector<Vec3> points;
    double rootRadius = 0, tipRadius = 0;
    Vec3 normal;
    float key = 0;                    // budget key
    double residual = 0;              // largest vertex move of the last settle step [mm]
};

// The static settle of one fibre, in place; gravity scaled for a fibre of visible length gravityLength
// (gravityFullLength, 0: off), the body kept collisionOffset beyond each vertex's radius. False: stuck in
// the body (hairInsideBody).
bool settleHair(Hair& h, const StrandSettleParams& settle, double gravityFullLength, double collisionOffset, const MeshQuery& body,
                double gravityLength);
// Follow-The-Leader from points[2] on: each vertex at its rest length L[i] from its parent, pushed out
// of the body to radii[i] + offset (a few rounds). d (optional): each vertex's total correction.
void followTheLeader(std::vector<Vec3>& x, const std::vector<double>& L, const std::vector<double>& radii, const MeshQuery& body,
                     double offset, std::vector<Vec3>* d = nullptr);
// Whether a vertex from points[2] on lies inside the body, or a segment from points[1] on crosses its
// surface (the first one past kRootCrossingSkip from its root, which lies on the surface).
bool hairInsideBody(const Hair& h, const MeshQuery& body);
constexpr double kRootCrossingSkip = 1e-4;  // [mm]
// Visible length: the polyline from points[1] on [mm].
double visibleLength(const Hair& h);
// v rotated by angle a [rad] about the unit axis k (Rodrigues' formula).
Vec3 rotateAbout(const Vec3& v, const Vec3& k, double a);

} // namespace knit
