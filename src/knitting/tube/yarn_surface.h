#pragma once
// Yarn surfaces (tubes): the model of rings along each yarn, display meshes, a spatial index, gate
// (d), and the description of the final rendered surface. Compression is a contact-driven
// approximation.
//
// The authoritative surface model is sampled at a FIXED fine arc-length spacing
// per yarn (independent of display tessellation):
//   * frames: rotation-minimizing (double reflection, Wang et al. 2008), one
//     continuous frame along the whole yarn, so the angle reference (and therefore
//     the UV phase) is consistent across spline spans;
//   * cross-section: an ellipse per ring (semi-axes a, b; a-axis at angle phi in
//     the frame), circular radius r = D/2 unless compressed.
// Display meshes resample this model; fibres anchor on it.
//
// Compression (geometric approximation, NOT mechanics): every ring takes its
// requirement from the actual centreline distance d to nearby tube portions; where
// d < 2r (an overlap of the nominal circles) BOTH sides are flattened along the
// contact direction so that their supports meet (h = d/2 each), never further, so
// contact is kept rather than replaced by a gap. The other axis keeps r at most, so a
// squashed yarn loses section area; no equilibrium or pressure is computed. Gate (d)
// (checkSurface) tests the result independently.
#include <cstdint>
#include <string>
#include <vector>

#include "model/parallel.h"
#include "model/params.h"
#include "model/spatial_hash.h"
#include "model/status.h"
#include "model/yarn_path.h"

namespace knit {

struct Ring {
    double s = 0, u = 0;
    Vec3 c, T, n1, n2;       // centre, unit tangent, RMF normal / binormal
    double a = 0, b = 0;     // semi-axes [mm]
    double phi = 0;          // direction of the a-axis in the (n1, n2) plane [rad]
};

struct YarnSurface {
    int yarnId = 0;
    double radius = 0;       // nominal r
    double length = 0;
    double ds = 0;           // ring spacing
    std::vector<Ring> rings;
};

// The compression's cross-section fit (compression_fit.h): the support is
// bounded towards every contact direction (up to kMaxReq) and along the principal normal (no fold
// of the swept cross-section); the free axis keeps r, so a squashed yarn loses section area.
// A scene's flattening limit is strand.max_flattening (FinalSurfaceSpec::maxCompression).
constexpr double kDefaultMaxCompression = 0.4;
constexpr double kMinMaxCompression = 0.05, kMaxMaxCompression = 0.6;  // the setting's range

struct SurfaceParams {
    double ringSpacing = 0.025;  // fixed internal spacing of the model [mm]
    bool compression = true;
    double maxCompression = kDefaultMaxCompression;  // b >= (1 - maxCompression) r
    double maxExpansion = 0.5;   // the own-yarn window of the search, pi (1 + maxExpansion) r; FinalSurfaceSpec::supportRadius
    double smoothing = 0.05;     // sigma of the max-filter that widens the flattening along s [mm]; 0 = off
    bool operator==(const SurfaceParams&) const = default;
};

struct CompressionReport {
    int clamped = 0;              // needed more than maxCompression
    int dropped = 0;              // CPU: rings whose full requirement set gave up a binding h
    double minAxisRatio = 1;      // min b / r
    double maxAreaChange = 0;     // max |a b / r^2 - 1|
};

struct SurfaceModel {
    std::vector<YarnSurface> yarns;
    CompressionReport compression;
    double ms = 0;
};

// sampleRings, transportFrames, then with compression collectContactRequirements,
// smoothRequirements, fitSections. exec: per-ring evaluation and the compression requirement in
// parallel (the rotation-minimizing frame stays sequential); identical model.
SurfaceModel buildSurface(const CurveSet& cs, const PatternParams& pp, const SurfaceParams& sp, const ExecPolicy& exec = {});

// Interpolated ring at arc length s (centre and tangent evaluated on the spline).
Ring ringAt(const CurveSet& cs, const YarnSurface& ys, double s);
// Support of the cross-section along a unit direction in the ring plane.
double ringSupport(const Ring& r, const Vec3& dir);

struct YarnSurfacePoint {
    Vec3 p, n;     // position and outward normal of the final surface
    Ring ring;
};
// theta is the ellipse parameter (angle measured from n1).
YarnSurfacePoint evalSurface(const CurveSet& cs, const YarnSurface& ys, double s, double theta);

// Yarn coordinates (s, theta, lambda): X = c(s) + lambda * E(s, theta), with E the ring
// ellipse point of ringAt(s) at theta (evalSurface's convention; lambda = 1 is the
// surface), and the analytic partial derivatives. Within one
// model-ring interval the ring frame, axes and phi are interpolated linearly, so the
// derivative along s is that of ringAt's interpolation (kinks at model rings). s must lie
// in [0, length]: ringAt clamps outside, so the caller rejects such points.
struct YarnCoordPoint {
    Vec3 X, dS, dTheta, dLambda;
};
YarnCoordPoint yarnCoordinates(const CurveSet& cs, const YarnSurface& ys, double s, double theta, double lambda);
// yarnCoordinates(...).X alone, by the same positional arithmetic (no second
// derivative, no frame derivatives): bit-identical to it.
Vec3 yarnPosition(const CurveSet& cs, const YarnSurface& ys, double s, double theta, double lambda);

// Display mesh: rings every `ringSpacing` (resampled from the model), `sides` around,
// seam duplicated so UV v runs 0..1 without wrapping.
struct TubeParams {
    double ringSpacing = 0.1;
    int sides = 16;
    double uvRepeat = 2.5;   // u = s / uvRepeat (physical arc length)
    // > 0 keeps u below uvWrap + 1 (float32 u loses texel precision on long yarns: 2^-9 at
    // u = 27000). Where u would reach uvWrap, the ring is emitted twice, at the same positions
    // and normals, the copy with u lowered by uvWrap (an integer: maps of period 1 in u stay
    // continuous); the strip continues from the copy. 0: u = s / uvRepeat throughout.
    int uvWrap = 0;
    // Per-yarn variation of the ply map (strand.ply_phase_jitter, ply_pitch_jitter): each yarn's u gets its own
    // phase in [0, uvPhaseJitter) and its own pitch factor in [1 - uvPitchJitter, 1 + uvPitchJitter], from a
    // hash of its id, so neighbouring yarns' ridges stop lining up. 0, 0: u = s / uvRepeat exactly.
    double uvPhaseJitter = 0, uvPitchJitter = 0;
    bool operator==(const TubeParams&) const = default;
};

// One yarn's ply-map variation (TubeParams::uvPhaseJitter, uvPitchJitter), from a hash of its index:
// u = (s / uvRepeat) / pitchFactor + phase. Without jitter 1 and 0, which leave u bit for bit.
// The CPU tube (buildTubeMesh) and the device's (gpu::buildTubeMesh) both take it from here.
struct YarnUv {
    double pitchFactor = 1, phase = 0;
};
YarnUv yarnUvOf(const TubeParams& tp, int yarnIndex);

struct TubeMesh {
    std::vector<Vec3> pos, nrm;
    std::vector<float> uv;            // 2 per vertex: (s / uvRepeat, theta / 2 pi)
    std::vector<double> s, t;         // arc length and t = s / L of the ring
    std::vector<double> theta;        // ellipse parameter
    std::vector<int> yarn;
    std::vector<float> axisRatio;     // b / r of the ring (compression)
    std::vector<uint32_t> tri;
    std::vector<uint32_t> yarnFirstVertex;  // per yarn, plus the end
    double ms = 0;
};

// exec: rings in parallel at precomputed offsets (identical mesh).
TubeMesh buildTubeMesh(const CurveSet& cs, const SurfaceModel& sm, const TubeParams& tp, const ExecPolicy& exec = {});

// ------------------------------------------------------------------------------
// Spatial index over the final surfaces.
//   * conservative enclosure: capsules (~0.1 mm) around the centreline with
//     radius max(a, b) + chord sag + margin; a miss is certain, a hit only
//     "may intersect";
//   * exact test: the ring ellipse of the final cross-section (inflated by a given
//     radius) in the plane through the query point.
// Exclusions are explicit arc-length ranges on one yarn, never whole yarns.
class SurfaceIndex {
public:
    SurfaceIndex(const CurveSet& cs, const SurfaceModel& sm);

    struct Exclude {
        int yarn = -1;
        double s0 = 0, s1 = -1;  // empty when s1 < s0
        bool covers(int y, double s) const { return y == yarn && s >= s0 && s <= s1; }
    };
    bool pointInEnclosure(const Vec3& q, double inflate, const Exclude& ex) const;
    bool segmentHitsEnclosure(const Vec3& a, const Vec3& b, double inflate, const Exclude& ex) const;

    struct Hit {
        bool inside = false;
        int yarn = -1;
        double s = 0;
        double depth = 0;  // approximate penetration depth of q [mm]
        Vec3 normal;       // the strand's outward normal in the ring plane at q (the ellipse's gradient)
    };
    // Exact inside test; `only` (if yarn >= 0) restricts the test to one yarn range,
    // `ex` removes a range. Returns the deepest hit.
    Hit exactInside(const Vec3& q, double inflate, const Exclude& only, const Exclude& ex) const;
    struct NearbyAxis {
        int yarn;
        double s;  // arc length of the closest axis point
    };
    // Every capsule axis within `radius` of q that is not excluded, into out (cleared first).
    void nearbyAxes(const Vec3& q, double radius, const Exclude& ex, std::vector<NearbyAxis>& out) const;
    double zTop() const { return zTop_; }
    double zBottom() const { return zBottom_; }  // the lowest point of the capsules

private:
    struct Capsule { Vec3 a, b; double R, s0, s1; int yarn; };
    const CurveSet& cs_;
    const SurfaceModel& sm_;
    std::vector<Capsule> caps_;
    SpatialHash3 hash_;  // capsule boxes (3D: correct and fast on curved pieces too)
    double zTop_ = 0, zBottom_ = 0;
    // q inside capsule c's strand at its nearest ring: the hit's yarn, ring arc length, depth and normal.
    bool ringInside(const Capsule& c, const Vec3& q, double inflate, Hit& hit) const;
};

// ------------------------------------------------------------------------------
// Independent check of the complete final surface.
//   * local folds: for a cross-section in the normal plane swept along the
//     centreline, the sweep is locally injective where support(kappa_hat) * kappa
//     < 1 (kappa_hat: principal normal). Evaluated on every model ring (variation of
//     the ellipse along s is ignored);
//   * collisions: surface points sampled every `ringStride` model rings and
//     `angles` around are tested exactly against the other tubes and against
//     distant portions of the same yarn (arc distance > pi * extent).
// This is a sampled check (stated resolution), not a certificate.
struct SurfaceCheckParams {
    int ringStride = 2;           // every 2nd model ring (0.05 mm)
    int angles = 24;
    double depthTolerance = 1e-3; // [mm] (1 um; the polyline sampling bound is 5 um)
};

struct SurfaceCheck {
    // folds
    double minFoldRatio = 1e300;  // min rho / support over rings (< 1: fold)
    int foldRings = 0;
    int worstFoldYarn = -1;
    double worstFoldS = 0;
    Status fold = Status::Pass;
    // collisions
    long pointsTested = 0;
    long penetratingPoints = 0;
    double maxDepth = 0;
    int worstYarnA = -1, worstYarnB = -1;
    double worstSA = 0, worstSB = 0;
    Vec3 worstPoint;
    Status collisions = Status::Pass;
    // Rings clamped at maxCompression (-1: not measured; the device's gate (d) fills it, the CPU gate
    // reads the model's CompressionReport).
    int compressionClamped = -1;
    std::string resolution;
    double ms = 0;
};

// checkFolds, then checkCollisions. exec: the sampled rings in parallel, merged in ring order
// (identical result).
SurfaceCheck checkSurface(const CurveSet& cs, const SurfaceModel& sm, const SurfaceCheckParams& p = {}, const ExecPolicy& exec = {});

// ------------------------------------------------------------------------------
// ONE description of the final (rendered) yarn surface, passed to every surface consumer: gate (d),
// fibre setup and collision tests, the export tubes, the cross-component report, the viewer's
// stand-in inset and the boundary check B1.
//
// Gate (a) validates the CENTRELINE at D (and c). The rendered tube may be inflated to
// D_r = D (1 + e) with compression, so contacts with a centre distance d < D_r are flattened
// (validated centrelines never overlap at D: compression alone is a no-op, as measured). Gate (d)
// then checks that rendered surface. Unsupported: flattening set by hand with explicit plies
// (strand.flatten; a yarn compressed to fit draws its round plies around the full yarn, past the
// validated envelope); inflation without compression (overlapping tubes). The round-tube fibre
// reference check (fibreReferenceCheck) is not applicable to compressed tubes: UNRESOLVED, never
// PASS.
struct FinalSurfaceSpec {
    double D = 1.0, c = 0.1;      // centreline validation (gate a)
    double inflate = 0;           // e: rendered D_r = D (1 + e)
    bool compression = false;     // flattening of the rendered surface at contacts
    int plies = 1;                // > 1: explicit round plies inside the drawn yarn
    bool compressedYarn = false;  // knitted compressed (yarnFit) and drawn at its full diameter
    double maxCompression = kDefaultMaxCompression;  // the flattening limit (strand.max_flattening)

    double renderedDiameter() const { return D * (1 + inflate); }
    // How far the rendered yarn reaches from its centreline: D_r / 2, or at most
    // (1 + SurfaceParams::maxExpansion) D_r / 2 with compression (a bound on the ellipse's long
    // axis, which the fit keeps at r: stored footprints depend on this value); explicit plies
    // stay inside the round yarn D_r / 2.
    double supportRadius() const;
    // pp with the rendered diameter (geometry consumers); gate (a) keeps the centreline pp.
    PatternParams rendered(const PatternParams& centreline) const;
    bool referenceCheckApplicable() const { return !compression; }
    bool valid(std::string* why) const;
};

} // namespace knit
