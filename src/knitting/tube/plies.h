#pragma once
// Explicit twisted plies inside the yarn.
//
// The yarn's knitted centreline and its validated envelope (radius r) stay the
// authority; n plies are built INSIDE that envelope, following Zhao, Luan & Bala 2016
// (§3.2, Algorithm 1: plies twisted around the yarn centre) and Wu & Yuksel 2017 (§3.1,
// ply phase offsets 2 pi j / n).
//
// Packing (n >= 2): circles of radius r_p0 = r sin(pi/n) / (1 + sin(pi/n)) touching each
// other around a centre at distance R = r - r_p0, inside the circle of radius r. The
// built plies have radius r_p = fill r_p0 (fill < 1: a margin; whether a ply folds is
// decided by the certified curvature check, rho_lo >= r_p). n = 1 is the plain yarn
// (R = 0, r_p = r).
//
// Ply j's centre, in the surface model's yarn coordinates (yarnCoordinates, the
// rotation-minimizing ring frame, so the plies carry only the authored twist):
//   C_j(s) = X(s, theta_j(s), R / r),  theta_j(s) = phase + h 2 pi s / pitch + 2 pi j / n,
// sampled every min(pitch / samplesPerTwist, maxSampleStep) and interpolated by a uniform
// cubic B-spline.
//
// Containment: a point within r_p of the fitted ply centre lies within r_p + R + e of
// the yarn centreline, e = the measured (sampled) distance of the fitted ply centre from
// the circle of radius R around the yarn centreline. With e < c / 2, plies of different
// yarn portions stay apart because the yarn envelopes are certified by gates (a) / (d).
// A yarn compressed to fit (scene_settings.h yarnFit) is the exception: r is its drawn radius,
// past the validated k D / 2, so round plies of neighbouring yarns may intersect at contacts.
#include <vector>

#include "model/parallel.h"
#include "knitting/tube/yarn_surface.h"

namespace knit {

struct PlyGeometry {
    int count = 3;              // 1 = the plain yarn tube
    double pitch = 2.5;         // [mm per full twist]
    int handedness = +1;        // +1 Z twist, -1 S twist
    double phase = 0.0;         // [rad]
    double fill = 0.95;         // r_p = fill * r_p0
    int samplesPerTwist = 24;
    double maxSampleStep = 0.1;  // [mm] the samples also resolve the yarn's own bending
    bool operator==(const PlyGeometry&) const = default;
};

struct PlyLayout {
    double rp0 = 0, rp = 0, R = 0;
};
PlyLayout plyLayout(const PlyGeometry& g, double r);

struct PlySet {
    CurveSet curves;            // one YarnPath per ply (kind Surface)
    std::vector<int> yarnOf, plyOf;
    PlyLayout layout;
    double containment = 0;     // e: max over dense samples of (distance of the ply centre to the
                                // yarn centreline) - R [mm] (sampled, both signs: |.|)
    double minCentreDistance = 0;  // min distance between plies of one yarn at equal s [mm]
    std::vector<double> parentScale;  // per ply: parent arc length per unit of ply parameter u
    std::vector<double> parentOffset; // per ply: parent arc length at u = 0 (the plies start r / 2 inside the yarn)
    double ms = 0;
};

// Plies of every yarn of `cs` from its surface model `sm` (round or compressed rings:
// the centre follows the ring ellipse scaled by R / r). For count == 1 the result's
// curves are copies of cs's yarns.
// samplePlyCentres, then fitPly per ply. exec: samples and plies in parallel, appended in ply order
// (identical set).
PlySet buildPlies(const CurveSet& cs, const SurfaceModel& sm, const PlyGeometry& g, double r, const ExecPolicy& exec = {});

// (parent yarn, parent arc length) of ply `ply` at its curve parameter u (FibreFrame::plies).
std::pair<int, double> plyParent(const PlySet& ps, int ply, double u);

// Fibre candidates per ply (FibreFrame::candidates) for `density` per mm of the
// PARENT yarn: N_y = floor(density L_y) for parent yarn y (L_y its arc length in
// `parent`), split over its plies in proportion to their arc lengths by largest
// remainder (ties: lower ply index), so the counts of a yarn's plies sum to N_y exactly
// and pitch, ply count and fill do not change the total.
std::vector<int> plyCandidateBudget(const PlySet& ps, const CurveSet& parent, double density);

// Interpolating uniform cubic B-spline through points P[0..m] at u = 0..m (end
// conditions: zero second difference of the control polygon), m >= 1.
BSpline interpolateBSpline(const std::vector<Vec3>& P);

} // namespace knit
