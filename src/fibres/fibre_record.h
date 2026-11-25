#pragma once
// One accepted fibre, as the GPU kernel writes it and the host keeps it: plain data compiled by both
// sides, so nothing is converted in between. Its points are not here: a run's fibres share one
// packed pool (FibreSet::points), the layout of StrandSetF, of OptiX curve buffers and of USD
// curves. Most fibres have about 4 points and a few have 50.
#include <cstddef>
#include <cstdint>

namespace knit {

struct Vec3f {
    float x = 0, y = 0, z = 0;
};

enum class FibrePopulation : uint8_t { Directed = 0, Random = 1 };

// How much of each island kind a point of the yarn has, 0 .. 1, strength included. All zero: the even
// fuzz, and every formula that reads it reduces to its plain form.
struct IslandWeight {
    float fuzzy = 0, matted = 0, worn = 0;
};
// What a full island does. Fuzzy: more strays, longer fibres. Matted: every fibre a stray lying flat,
// and shorter. Worn: fewer fibres, so the yarn shows.
constexpr float kFuzzyStrays = 0.35f, kFuzzyLength = 1.5f;
constexpr float kMattedFlat = 0.6f, kMattedLength = 0.3f;

struct Fibre {
    uint32_t firstPoint = 0;           // into FibreSet::points
    int pointCount = 0;                // [0] the buried root, [1] the curve start, then the visible vertices
    int strand = 0;                    // position of the supporting strand in SurfaceModel::yarns
    int candidate = -1;                // index on it: with strand, the fibre's stable identity
    int yarnId = 0;                    // supporting strand (a ply, or the yarn itself)
    int parentYarn = 0;                // the yarn (plyParent on plies; else yarnId)
    FibrePopulation population = FibrePopulation::Directed;
    bool truncated = false;
    double s = 0;                      // the candidate's anchor: arc length along the strand
    double key = 0;                    // thinning rank
    float theta = 0;                   // and ellipse parameter
    Vec3f anchor, normal;              // the emergence point; the surface normal at the candidate's anchor
    float emergeL = 0;                 // polyline arc length from the root to the emergence
    float radius0 = 0, radius1 = 0;    // root and tip radius
    float targetLength = 0;            // drawn visible length [mm]
    float attach = -1, attachArc = -1; // the fibre's attachment window (fibre / yarn arc) [mm]
    // Reported only (deg, mm).
    float psi = 0, lift = 0;           // nominal draws (Random: azimuth, start elevation)
    float smoothLength = 0;            // visible length of the smooth curve
    float emergeAngle = 0;             // first chord to the surface
    float helixActual = 0;             // physical angle of the smooth curve to +-T at the anchor
    float maxDev = 0;                  // smooth curve to polyline
    float minTurn = 0, maxTurn = 0;    // world turn at interior vertices
    int segments() const { return pointCount - 1; }
};

// Device and host must agree on this layout: both compilers evaluate these.
static_assert(sizeof(Vec3f) == 12 && sizeof(Fibre) == 136, "the fibre record's layout is shared with the GPU kernels");
static_assert(offsetof(Fibre, s) == 32 && offsetof(Fibre, theta) == 48 && offsetof(Fibre, maxTurn) == 128, "the fibre record's layout");

} // namespace knit
