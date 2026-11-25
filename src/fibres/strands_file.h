#pragma once
// photon_path_tracer strand files (.pptstrands v1), and the renderer's rules every written strand keeps.
//
// Layout (little-endian; photon_path_tracer/src/authoring/strands_file.cpp):
//   char[8] "PPTSTRND" | u32 version = 1 | u32 basis = 0 (round linear)
//   u64 V (vertices) | u64 N (strands) | u64 reserved = 0          -- 40 bytes
//   f32 xyz[V] | f32 radius[V] | u32 strand_first[N + 1]
// File size is exactly 40 + 16 V + 4 (N + 1). Positions and radii are in the
// scene's units (meters, Y-up for the renderer); radius is a radius, not a width.
//
// checkStrands validates the float32 values that are actually written, against
// the renderer's reader (strand_set_error) and the admissibility bounds of its
// round-linear primitive (the renderer's strand-format decision):
//   - every strand has >= 2 vertices; strand_first runs 0 .. V, monotone;
//   - positions and radii finite, radius >= 1e-6 (object space);
//   - segment length within [minSegOverR, maxSegOverR] times the smaller end
//     radius (the renderer needs >= 2 r and length / (2 r) <= 100);
//   - adjacent vertex radii differ by at most maxRadiusStep (50 %);
//   - discrete curvature radius (circumradius of consecutive vertex triples)
//     >= the smallest radius of the triple;
//   - precision envelope (a renderer WARNING, not a refusal): radius >=
//     precisionK * max |coordinate| of the vertex.
#include <cstdint>
#include <string>
#include <vector>
#include "model/coordinates.h"

namespace knit {

struct StrandSetF {
    std::vector<float> xyz;              // 3 V
    std::vector<float> radius;           // V
    std::vector<uint32_t> first{0};      // N + 1, first[0] = 0, first[N] = V
    size_t vertices() const { return radius.size(); }
    size_t strands() const { return first.empty() ? 0 : first.size() - 1; }
    // Append one strand (positions and radii already in output units).
    void add(const std::vector<float>& pts, const std::vector<float>& radii);
};

struct StrandRules {
    double minRadius = 1e-6;
    double minSegOverR = 4.0;     // stricter than the renderer's 2 r
    double maxSegOverR = 200.0;   // length / (2 r) <= 100
    double maxRadiusStep = 0.5;   // max(r_a, r_b) / min(r_a, r_b) <= 1 + step
    double minSegOverRMax = 0;    // segment >= this x the LARGER end radius (0: not checked; taperRules)
    double precisionK = 2.4e-6;   // renderer envelope (warning only)
};

// How the generators keep these rules with a margin: every segment at least kMinSegmentRadii radii (the
// rule's 4 r), and adjacent radii in a ratio of at most kTaperStep (the rule allows 1.5).
constexpr double kMinSegmentRadii = 4.2;
constexpr double kTaperStep = 1.45;
// The radius at vertex k of a strand of `segments` segments: geometric from r0 towards r1, each step at most
// kTaperStep, so a strand with too few segments for r1 ends thicker than r1.
double taperedRadius(double r0, double r1, int k, int segments);

// The first rule a strand set breaks, so a caller can count refusals by reason.
enum class StrandFault { None, Structure, NonFinite, RadiusFloor, SegmentShort, SegmentLong, RadiusStep, Curvature };
constexpr int kStrandFaultCount = 8;
const char* strandFaultName(StrandFault f);

struct StrandCheck {
    bool ok = true;
    std::string error;            // first violation (with the strand / vertex index)
    StrandFault fault = StrandFault::None;
    double minSegOverR = 1e300, maxSegOverR = 0, maxRadiusRatio = 1, minCurvOverR = 1e300;
    double minRadius = 1e300, maxRadius = 0, maxCoord = 0;
    int precisionWarnings = 0;    // vertices outside the precision envelope
    double worstPrecision = 0;    // max over vertices of precisionK |p|inf / r (< 1 is inside)
};

StrandCheck checkStrands(const StrandSetF& s, const StrandRules& rules = {});

// Refuses (returns false, writes nothing) if checkStrands fails or the set is empty.
bool writeStrandsFile(const std::string& path, const StrandSetF& s, const StrandRules& rules, std::string* err,
                      StrandCheck* check = nullptr);

// Independent reader with the renderer's refusal rules (magic, version, basis,
// plausible counts, exact file size, no trailing bytes, CSR table, >= 2 vertices,
// finite values, radius >= 1e-6).
bool readStrandsFile(const std::string& path, StrandSetF& out, std::string* err);

// The strands of `s` (one key per strand) with key < q, copied whole.
StrandSetF keepStrands(const StrandSetF& s, const std::vector<float>& keys, float q);

} // namespace knit
