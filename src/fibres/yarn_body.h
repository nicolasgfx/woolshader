#pragma once
// The yarn body made of structural fibres (prototype): 1 to 3 plies per
// yarn, each a bundle of long fibres in a shell of finite, irregular depth, optionally with interior
// fibres and a core proxy. A TEMPORARY CPU PROTOTYPE behind yarn_body.enabled (off by default) for a
// visual experiment; it is not the project's GPU path and is not meant to stay as one.
//
// Construction: every fibre is built in a straight nominal yarn (x, y, z), z the yarn's
// arc length and (x, y) its round section of radius r, then bent onto the envelope by
//   X(z, x, y) = c(z) + [n1 n2](z) M(z) (x, y),   M = Rot(phi) diag(a / r, b / r) Rot(-phi),
// the rings of the surface model (ringAt: the rotation-minimizing frame, thickness and compression).
// Plies orbit the axis at the rate omega = tan(ply angle) / R (zero for one ply); each ply carries the
// helix's own frame (u = radial, v = t x u) and its fibres rotate about the ply axis at the rate the
// fibre angle asks for, measured against a non-twisting ply frame (the helix frame's torsion is taken
// off). Rates are integrated on a fixed arc-length grid, so phases are continuous and independent of
// where the curves are sampled. Everything is a function of (part, yarn, ply, slot, piece) seeds and
// z: a rebuild is bit-identical, and nothing restarts at a stitch boundary.
//
// The section map M is the rings' ellipses smoothed along the yarn (sigma 0.1 mm): the compression fit
// turns an ellipse abruptly where its binding contact changes, which fibres would show as hairpins.
// Guarantees: every fibre centreline lies in the smoothed ellipse; with the fibre's own radius it
// protrudes by at most r_max (1 - its shortest semi-axis / r) where the section is compressed (measured:
// maxProtrusion). Against each true ring the centreline's excess is measured, not bounded (maxEnvelope).
// NOT guaranteed: fibres do not interpenetrate (measured in windows: BodyStats overlap), and yarns
// overlap wherever their drawn tubes do.
#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "model/mesh.h"
#include "model/parallel.h"
#include "model/vec.h"

namespace knit {

struct CurveSet;
struct Ring;
struct SurfaceModel;
struct YarnSurface;

// The section map of one yarn: M = Rot(phi) diag(a, b) / r Rot(-phi) of every model ring as its
// symmetric entries (m11, m12, m22) in the ring's (n1, n2), smoothed along s by a Gaussian of 0.1 mm.
// The compression fit turns a ring's ellipse abruptly where its binding contact changes (up to 0.9 rad
// between neighbouring rings on the swatch): a closed tube hides that, a fibre would show it as a
// hairpin. A mean of symmetric positive definite matrices is one again; the fibres then follow the
// smoothed ellipse, and their excess over each true ring is measured (maxEnvelope). The bake
// (yarn_bake.h) inverts the same map to put a fibre point back on the tube's chart.
struct SectionTable {
    double ds = 1;
    std::vector<std::array<double, 3>> m;
    std::array<double, 3> at(double s) const;
};
SectionTable sectionTable(const YarnSurface& ys);
// The smoothed section map applied to nominal (x, y): X = c + [n1 n2] M (x, y).
Vec3 sectionOffset(const Ring& g, const std::array<double, 3>& m, double x, double y);

enum class BodyCore { None, Opaque, Translucent };

struct YarnBodyParams {
    int plies = 1;                // 1 to 3
    int fibres = 160;             // shell fibres per yarn, split evenly across its plies
    int interiorFibres = 0;       // fibres inside the shell, per yarn (the dense reference); 0 = none
    double fibreRadius = 0.012;   // [mm]
    double shellDepth = 0.35;     // fraction of the ply's centreline radius; 1 = filled
    double fibreAngleDeg = 20;    // surface fibres against the ply axis
    double plyAngleDeg = 20;      // ply centre against the yarn axis (unused for one ply)
    double irregularity = 1;      // 0 = regular fibres (README: what stays random), 2 = most
    BodyCore core = BodyCore::None;
    uint64_t seed = 1;
};

// The radii a ply's fibres are placed in [mm].
struct BodyLayout {
    double r = 0;         // nominal yarn radius
    double R = 0, rp = 0; // ply centre distance and ply radius (plies.h plyLayout)
    double rMax = 0;      // largest fibre radius
    double crimp = 0;     // crimp amplitude bound
    double rOut = 0;      // outermost centreline radius around the ply centre: rp - rMax - crimp
    double rIn = 0;       // innermost shell centreline radius: rOut (1 - shellDepth)
    double rCore = 0;     // core proxy radius: rIn - rMax
    double phiShell = 0;  // fibre area / shell annulus area, per ply
    double phiPly = 0;    // fibre area / whole ply area (Trunz et al.'s definition)
};
BodyLayout bodyLayout(const YarnBodyParams& p, double yarnRadius);

// What a refused yarn body is refused for, so a caller can name the setting to change.
enum class BodyRefusal { None, Plies, FibreCount, FibreTooThin, FibreTooLarge, ShellDepth, Angles, Irregularity, ShellFull,
                         InteriorFull, NoCoreRoom };

// Refused with the reason (and the largest count or depth that would be accepted): plies outside 1-3,
// fibres too large for the ply, more fibres than the shell or interior can hold, a core with no room.
BodyRefusal yarnBodyRefusal(const YarnBodyParams& p, double yarnRadius, std::string* err);
inline bool yarnBodyParamsValid(const YarnBodyParams& p, double yarnRadius, std::string* err) {
    return yarnBodyRefusal(p, yarnRadius, err) == BodyRefusal::None;
}

// The fibres as polylines in the curves' space [mm]: the shell's strands first (yarn, ply, slot,
// piece order), then the interior's, so the shell is the same with or without interior fibres.
struct BodyStrands {
    std::vector<Vec3> points;
    std::vector<double> radius;          // per point
    std::vector<double> s;               // per point: the yarn arc length its section was taken at
    std::vector<uint64_t> first{0};      // per strand, plus the end (64 bit: a body may pass 2^32 points)
    std::vector<int> yarnOf;             // per strand
    size_t strands() const { return yarnOf.size(); }
};
// Appends from's strands after to's: their point offsets shifted by to's point count (first.back()).
void appendStrands(BodyStrands& to, const BodyStrands& from);

// Interpenetration in sample windows, by segment-to-segment distance.
struct OverlapWindow {
    int yarn = 0;
    double s0 = 0, s1 = 0;
};
struct OverlapStats {
    double fibreLength = 0;                     // [mm] in the windows
    long long pairsWithinYarn = 0, pairsBetweenYarns = 0;  // penetrating segment pairs
    long long segments = 0, segmentsInvolved = 0;
    double maxDepth = 0;                        // max (r_i + r_j - d) / (r_i + r_j)
};
OverlapStats measureOverlap(const BodyStrands& strands, const std::vector<OverlapWindow>& windows);

struct BodyStats {
    BodyLayout layout;
    long long fibres = 0, pieces = 0, droppedPieces = 0, segments = 0;
    long long softTurnSegments = 0;  // at the shortest allowed length, still above the 6 degree target
    double maxTurnDeg = 0, p99TurnDeg = 0;
    double yarnLength = 0, fibreLength = 0;  // [mm]
    double maxProtrusion = 0;   // [mm] upper bound over vertices of radius - distance to the (smoothed) envelope
    double maxEnvelope = 0;     // max centreline radius in its true ring ellipse (1 = on it; above 1: the smoothing's excess)
    OverlapStats overlap;
    double ms = 0;
    std::string line() const;
};

struct YarnBody {
    BodyStrands strands;
    TriMesh core;   // the core proxy's closed tubes (empty without one)
    BodyStats stats;
};

// A ply's axis in nominal coordinates (x, y in the yarn's section, z along it) at yarn arc length z:
// its centre and unit tangent. For the checks that measure fibre angles against the ply (tests, README).
struct PlyAxis {
    Vec2 centre;
    Vec3 tangent;
};
PlyAxis bodyPlyAxis(const YarnBodyParams& p, double yarnRadius, uint64_t partSeed, int yarn, int ply, double z);
// The same axis at the centre of every step along a yarn of `length`: z_i = (i + 0.5) step for
// i < ceil(length / step), the twist integrated once (the bake's ply cores, yarn_bake.h).
std::vector<PlyAxis> bodyPlyAxes(const YarnBodyParams& p, double yarnRadius, uint64_t partSeed, int yarn, int ply, double length,
                                 double step);

// The body of every yarn of `yarns` on its surface model `surface` (built from the same curves, at the
// rendered diameter). partSeed: the part's identity (its name's hash). False with the reason when the
// parameters are refused or a fibre cannot be placed (never fewer fibres than asked).
bool buildYarnBody(const CurveSet& yarns, const SurfaceModel& surface, const YarnBodyParams& p, uint64_t partSeed,
                   YarnBody& out, std::string* err, const ExecPolicy& exec = {});

} // namespace knit
