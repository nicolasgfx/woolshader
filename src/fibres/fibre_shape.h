#pragma once
// One fibre as pure mathematics, with no tube anywhere: the draws of a candidate, the two kinds of
// smooth curve, and the arc-length table that turns a curve into a polyline of a chosen length.
// This is the part the papers describe; fibre_generator.cpp adds what explicit tubes require.
//
// Draws are a pure function of (seed, yarn, i), so a fibre can be regenerated from its identity
// alone and the distributions can be tested before any rejection. They are the specification the
// GPU kernels implement (gpu/fibres/fibre_candidate.cuh): the hash channels may not change.
#include <array>
#include <cstdint>
#include <vector>

#include "fibres/fibre.h"

namespace knit {

struct FibreDraw {
    FibrePopulation population = FibrePopulation::Directed;
    double length = 0;        // target visible length [mm]
    double radiusFactor = 1;  // exp(thicknessSigma z)
    double key = 0;           // uniform [0, 1): thinning rank (the segment budget)
    // Directed.
    int sgn = 1;              // along +s or -s
    double psi = 0;           // nominal helix angle [rad] (signed by handedness)
    double lift = 0;          // nominal lift angle [rad]
    double bowTheta = 0;      // [rad]
    double bowRadial = 0;     // [mm]
    double crimpPhase = 0;    // start of the crimp wave [rad]; channel 28, drawn only when crimped
    double crimpPlane = 0;    // the wave's plane from the yarn surface towards its normal [rad]; channel 29
    // Random.
    double azimuth = 0;       // start tangent around the normal from the projected strand tangent [rad]
    double elev = 0;          // start tangent above the tangent plane [rad]
    double bendAngle = 0;     // bow direction around the start tangent [rad]
    double bendAmp = 0;       // signed bow amplitude a
    // Curled strays (radius 0: not curled), capped by the fibre's length.
    double curlRadius = 0;    // [mm]
    double curlTurns = 0;     // turns along the fibre, at most maxSegments / kCurlChordsPerTurn
    double curlPhase = 0, curlDrift = 0;  // start angle, phase of the irregularity [rad]
    int curlSense = 1;
};

FibreDraw drawFibre(const FibreParams& fp, uint64_t seed, int yarn, int i, const IslandWeight& island = {});
// Hash channel 17 alone, drawn before any rejection: Random iff u17 < randomFraction, so raising
// the fraction only turns Directed candidates Random. An island raises it where it lies.
FibrePopulation drawPopulation(const FibreParams& fp, uint64_t seed, int yarn, int i, const IslandWeight& island = {});
// The target length: with lengthSpread > 0 a truncated log-normal sampled by inverse CDF (no mass at
// the bounds), else the clamped exponential.
double fibreTargetLength(const FibreParams& fp, uint64_t seed, int yarn, int i);
// The parameters a population draws from (Random: its overrides applied).
FibreParams populationParams(const FibreParams& fp, FibrePopulation p);
bool hasRandomOverrides(const FibreParams& fp);
double normalCdf(double x);
double normalQuantile(double p);  // Acklam's rational approximation, refined by two Newton steps

// The curve's yarn coordinates at u and their derivatives.
struct YarnCoords {
    double s = 0, theta = 0, lambda = 0, ds = 0, dTheta = 0, dLambda = 0;
};

// A directed fibre, u in [0, 1] from the anchor (lambda = lambda0) to the tip (Zhao Eq. 11, linear
// in (s, theta, lambda), plus a half-sine bow):
//   s(u)      = sE + sgn S u
//   theta(u)  = thetaE + sgn (S p / r) u + bowTheta sin(pi u)
//   lambda(u) = lambda0 + (S q / r) u + bowLambda sin(pi u)
// X(u) = yarnCoordinates(s, theta, lambda). p and q are NOMINAL slopes; the physical angles follow
// from dX. A hugging fibre (hug = H > 0) rises with the same slope at the anchor and levels off
// under H, so its tip lies along the yarn (bowLambda is 0 then):
//   lambda(u) = lambda0 + H (1 - exp(-k u)),  k = S q / (r H)
// On the ply helix (rampWinding) the winding ramps in over kWindRampLength along the yarn: a steep
// helix from the anchor would cut back into the round tube on the first chord, and a sudden one
// turns too sharply at a vertex.
//   d theta / du = sgn (S p / r) (1 - exp(-S u / kWindRampLength))^2 + bow
constexpr double kWindRampLength = 0.8;  // [mm] along the yarn; with 0.45 mm chords every turn stays below 35 deg
struct DirectedCurve {
    const CurveSet* cs = nullptr;
    const YarnSurface* ys = nullptr;
    double sE = 0, thetaE = 0;
    int sgn = 1;
    double S = 0, p = 0, q = 0, bowTheta = 0, bowLambda = 0, r = 1;
    double hug = 0;      // FibreParams::hugHeight; 0: the linear rise
    bool rampWinding = false;  // FibreParams::plyPitch > 0
    double lambda0 = 1;  // 1: on the yarn surface; < 1: a ply inside it
    // Crimp: + A f(u) sin(2 pi (S / wavelength) u + phase), A in mm, split between theta (A cos(plane) / r)
    // and lambda (A sin(plane) / r); f a smoothstep from 0 at the start to 1 half a wave later, so
    // X'(0) and the root are the plain curve's. crimpAmplitude 0: no crimp.
    double crimpAmplitude = 0, crimpWavelength = 1, crimpPhase = 0, crimpPlane = 0;
    YarnCoords coords(double u) const;
    Vec3 X(double u) const;
    Vec3 dX(double u) const;  // analytic
};

// A random flyaway, u in [0, 1] from the anchor A:
//   X(u) = A + S (d0 u + a w u^2),  X'(0) = S d0,
// d0 the unit start tangent (azimuth / elevation in the anchor frame n, t1, t2 with t1 the strand
// tangent projected onto the tangent plane), w a unit vector perpendicular to d0 (angle bendAngle
// around it), a = bendAmp. A curled stray adds a helix around that axis, in the plane e1, e2
// perpendicular to d0:
//   + rho f(u) g(u) (e1 cos th(u) + e2 sin th(u)),
//   th(u) = phase + sense t(u) + 0.7 sin(0.6 t(u) + drift),  t(u) = 2 pi turns u,
//   g(u) = 1 + 0.3 cos(0.37 t(u) + drift),
// the sines keeping the coils from looking like a spring, f a smoothstep from 0 at the anchor to 1
// half a turn later (X'(0) stays S d0, so root and emergence are those of the plain curve).
constexpr int kCurlChordsPerTurn = 16;   // polyline chords per turn: vertex turns stay below 35 deg
constexpr double kCurlMaxChord = 2.0;    // [mm] longest chord of a curled stray
struct RandomCurve {
    Vec3 A, d0, w;
    double a = 0, S = 0;
    Vec3 e1, e2;
    double rho = 0, turns = 0, phase = 0, drift = 0, rampEnd = 1;
    int sense = 1;
    Vec3 X(double u) const;
    Vec3 dX(double u) const;
    Vec3 curl(double u) const;   // the helix offset at u
    Vec3 dCurl(double u) const;  // its derivative
};
// The curve of a draw at anchor A with surface normal n and strand tangent T (S = 0).
RandomCurve randomCurve(const FibreDraw& d, const Vec3& A, const Vec3& n, const Vec3& T);

constexpr int kTableChords = 32;
// One fibre's curve and its arc-length table: fits the curve's span S to a target length and
// samples it at equal arc length.
struct SmoothCurve {
    bool random = false;
    DirectedCurve directed;
    RandomCurve stray;
    double bowRadial = 0;   // the draw's outward bow [mm], capped anew for every span
    double maxSpan = 0;     // the solve's bound on S
    double lastSpan = -1;   // the last solved span: the next solve starts there
    std::array<Vec3, kTableChords + 1> table;
    std::array<double, kTableChords + 1> arc{};
    Vec3 point(double u) const { return random ? stray.X(u) : directed.X(u); }
    double length() const { return arc[kTableChords]; }
    // Builds the table for span S over `chords` chords (kTableChords divides by it); returns the length.
    double setSpan(double S, int chords);
    // Secant steps on S, bracketed by [0, maxSpan], then the final table. False: no span within the
    // bound reaches the target, or the secant does not converge.
    bool solveLength(double target);
    // nv points at equal arc length after the start, the last exactly the curve's end.
    std::vector<Vec3> equalArcVertices(int nv) const;
    // Largest distance of a table point to the polyline.
    double deviationFrom(const std::vector<Vec3>& polyline) const;
};

} // namespace knit
