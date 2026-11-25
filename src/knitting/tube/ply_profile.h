#pragma once
// The plies of plies.h as shading maps on the single yarn tube.
//
// Coordinates: a straight yarn of radius r along +s with the ring frame of the surface
// model (n1, n2); theta is the ring angle of buildTubeMesh (v = theta / 2 pi). Ply j's centre
// is at theta_j(s) = phase + h 2 pi s / P + 2 pi j / n (plies.h), so with the twist-relative
// angle
//   phi = theta - h 2 pi s / P - phase
// ply j sits at phi = 2 pi j / n and, on a straight yarn, the outer surface is a function of
// phi alone: rho(theta, s) = rho(phi).
//
// Three objects are kept apart:
//  (1) the circle layout: n circles of radius r_p at distance R in the PARENT plane (closed
//      form). It is not the cross-section of the plies of plies.h;
//  (2) the swept-ply model: the union of balls of radius r_p centred on the ply helices of plies.h
//      (round tubes perpendicular to each ply's own tangent; a union of balls because the
//      helix curvature radius exceeds r_p). Its parent-plane cut is the union over t of the
//      discs the balls at C(t) cut, so the outer radial profile of ply 0 is
//        rho_0(alpha) = max_t [ R cos(alpha - w t) + sqrt(r_p^2 - t^2 - R^2 sin^2(alpha - w t)) ],
//      w = h 2 pi / P, over the t where the root is real (sweptPlyHit);
//  (3) the tessellated ply mesh, built from the plies' curves.
//
// The shading profile (PlyShadingProfile) is an explicitly defined C1 approximation: it equals
// model (1) or (2) on each ply's flank and replaces the junction between neighbouring plies
// (a gap, a crease, or a cusp) by a bridge; see buildPlyShadingProfile.
#include <vector>

#include "model/choices.h"
#include "model/parallel.h"
#include "knitting/tube/plies.h"
#include "knitting/tube/yarn_surface.h"

namespace knit {

// Outermost intersection of the radial ray at angle phi with the plies (rho > 0), and its
// derivative d rho / d phi; hit = false where the ray misses every ply.
struct RadialHit {
    bool hit = false;
    double rho = 0, slope = 0;
};

// (2), ply 0 only (centre at phi = 0), w = h 2 pi / P [rad/mm]. Symmetric in alpha and in h.
RadialHit sweptPlyHit(const PlyLayout& l, double w, double alpha);

struct PlyProfileParams {
    PlyProfileModel model = PlyProfileModel::Swept;
    // The bridge starts where the flank's in-plane normal is tilted this far from radial
    // (atan(|d rho / d phi| / rho)), or earlier (minBlend).
    double blendTiltDeg = 60.0;
    double minBlend = 0.03;   // minimum bridge half-width [fraction of the half period pi / n]
    // Valley floor d_gap: the bridge never goes below this radius [mm]. < 0: the default,
    // R cos(pi / n), the tangency point of touching (fill 1) circles.
    double gapFloor = -1.0;
    // The maps' AO is the cosine-weighted visibility of the shading profile as the 3D helical height
    // field rho(theta - w s) (heightFieldAO): the grooves run at the helix angle, so a cross-section
    // extruded along the yarn axis would miss most of the occlusion.
    int aoDirections = 256;   // AO quadrature (per point)
    int tableSamples = 4096;  // rho table per half period (rhoMin, rhoMax; the AO marches surfaceTable)
};

// rho(phi) for phi in [0, pi / n] (ply 0's half period; the profile is even and
// 2 pi / n periodic):
//   [0, aL]          the model's flank (exact, model (1) or (2));
//   [aL, aL + Lh]    cubic Hermite from (rho(aL), rho'(aL)) to (floor, 0);
//   [aL + Lh, pi/n]  the constant floor;
// floor = max(d_gap, rho(aL) - |rho'(aL)| b / 2), b = pi / n - aL (the parabola that
// matches the flank and meets its mirror image with zero slope), and
// Lh = min(b, 3 (rho(aL) - floor) / |rho'(aL)|) (Fritsch-Carlson: the Hermite piece is then
// monotone, so the bridge stays >= floor). C1 at aL, aL + Lh, pi / n and phi = 0.
struct PlyShadingProfile {
    PlyGeometry g;
    PlyLayout layout;
    PlyProfileParams p;
    double r = 0;           // yarn radius [mm]
    double w = 0;           // twist rate h 2 pi / P [rad/mm]
    double half = 0;        // pi / n
    double edge = 0;        // end of ply 0's ray-hit region in [0, pi / n] (= pi / n: crease, no gap)
    bool gap = false;       // a radial ray misses every ply near phi = pi / n
    double aL = 0, Lh = 0;  // bridge start (angle from the ply centre) and Hermite length [rad]
    double p0 = 0, m0 = 0;  // flank value and slope at aL
    double floorRho = 0, gapFloor = 0;
    double rhoMin = 0, rhoMax = 0;
    std::vector<double> table;  // rho at half * i / (tableSamples - 1)
    // The model's own surface (no bridge) where a radial ray hits a ply, else the bridge: the
    // AO's points and occluders (AO needs no C1 profile, so it keeps the real crease).
    std::vector<double> surfaceTable;
    std::vector<double> aoTable;  // AO at half * i / 256

    // rho(phi) [mm] and d rho / d phi; n == 1: rho = r, slope 0.
    double eval(double phi, double* slope = nullptr) const;
    // rho and slope of the model's own surface where a radial ray hits a ply (exact), else of
    // the bridge; fastSurface interpolates surfaceTable.
    double surface(double phi, double* slope = nullptr) const;
    double fastSurface(double phi) const;
    double ao(double phi) const;    // linear interpolation of aoTable; n == 1: 1
    // Twist-relative angle of a yarn-surface point (plies.h convention).
    double phiAt(double s, double theta) const;
};

PlyShadingProfile buildPlyShadingProfile(const PlyGeometry& g, double r, const PlyProfileParams& p = {});

// 3D AO at phi: cosine-weighted visibility over the hemisphere of the height-field normal,
// against the helical height field r <= rho(theta - w s - phase) of the model's own surface
// (PlyShadingProfile::surface: the swept or circle union with its real crease; the bridge
// only where a ray misses every ply), exact for that surface on a straight yarn; `rays`
// stratified directions, deterministic.
double heightFieldAO(const PlyShadingProfile& prof, double phi, int rays);

// ---------------------------------------------------------------- maps ----
// Map domain: u in [0, 1) <-> s in [0, uvRepeat), v in [0, 1) <-> theta in [0, 2 pi), the
// tube's uvRepeat = P. Tangent space: T along +u (+s), B along +v (+theta), N outward.
struct PlyMapSample {
    double rho = 0, slope = 0;  // rho(phi), d rho / d phi
    double H = 0;               // (rho - rhoMin) / (rhoMax - rhoMin)
    Vec3 nTangent;              // normalize(-d rho/ds, -(1/rho) d rho/d theta, 1)
    double ao = 1;
};
PlyMapSample plyMapSample(const PlyShadingProfile& prof, double uvRepeat, double u, double v);

// Texel grids of the maps, row-major, row y = 0 at the TOP: texel (x, y) is the sample at
// u = (x + 0.5) / width, v = 1 - (y + 0.5) / height (the renderer's sampler flips v).
struct PlyMaps {
    int width = 0, height = 0;
    std::vector<float> H, ao;   // one per texel
    std::vector<float> normal;  // three per texel (tangent space, unit)
};
PlyMaps makePlyMaps(const PlyShadingProfile& prof, double uvRepeat, int width, int height);

} // namespace knit
