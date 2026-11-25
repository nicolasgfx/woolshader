#pragma once
// A FibreSet written as renderer strands (strands_file.h): the generated polylines with their radii
// (fibresToStrands), or with a late taper, which keeps a fibre's thickness along its body and narrows only
// its tail. The taper belongs to how a fibre is written, not to the fibre: the generator, the cache and
// the digests see the untapered fibre.
#include <string>
#include <vector>

#include "fibres/fibre.h"
#include "fibres/strands_file.h"

namespace knit {

// Fibres as strands, in scene units: positions through map.toExport, radii over map.sigma, the
// polyline exactly as generated with Fibre::vertexRadius. `population` >= 0 (FibrePopulation): only
// the fibres of that population, in fibre order.
StrandSetF fibresToStrands(const FibreSet& fs, const SceneMapping& map, int population = -1);

// A late taper (fibres.tail_taper*): the fibre keeps its root radius along its body and narrows
// only over its tail, by a cosine in arc length from the root radius to tip x root radius, so
// there is no shoulder and no point. The tail is the larger of `tail` x the visible length and
// `minDiameters` root diameters (never more than the visible length): a round-linear strand needs
// every segment at least 2 r long (the renderer's strand-format decision), and a short
// tail on a thick fibre has no room for the steps of a smooth narrowing. tail = 0 is no late taper.
//
// Where the tail's radii are carried (a developer switch, knit_check --taper-points; Dense is the export's):
//   Dense     a vertex wherever the radius has fallen by 1.35 (the renderer's 50 % step rule)
//   Existing  the profile evaluated at the generated vertices only: the body may thin before x0
//   Onset     Existing plus a vertex at the taper start x0, so the body keeps r0 up to it
//   Refined   Onset plus vertices until |r - r_linear| <= tolerance x r0 on every interval
// The sparse modes allow one radius step from r0 to tip x r0 within a segment (taperRules).
enum class TaperPoints { Dense, Existing, Onset, Refined };
bool parseTaperPoints(const std::string& text, TaperPoints& out);  // "dense", "existing", "onset", "refined"
const char* taperPointsName(TaperPoints p);

struct StrandTaper {
    double tail = 0;           // fraction of the visible length
    double minDiameters = 8;   // the shortest tail, in root diameters
    double tip = 0.1;          // the tip radius over the root radius
    TaperPoints points = TaperPoints::Dense;
    double tolerance = 0.05;   // Refined: the largest |r - r_linear| allowed, over r0
};

// The rules a tapered strand set is written with. Dense: StrandRules{}. The sparse modes: one step
// from r0 to tip x r0 in a segment, every segment at least 4 x its larger radius: the renderer's strand
// oracle reads internal joints as surface at 2 r_max with a 3.33x step and agrees at 4 r_max
// (photon_path_tracer dee6851a1, --strands-oracle); it also keeps |dr| < L. The step bound has a float32
// margin: two radii rounded independently change their ratio by at most about 2u (u = 2^-24), so
// the bound is (1 / tip)(1 + 4u) - 1, never below the default 0.5 (fallback strands step by 1.45).
StrandRules taperRules(const StrandTaper& taper);

// The strands of fibresToStrands with a late taper. The path is the generated polyline: new
// vertices lie on its chords. A fibre whose tapered strand would break a StrandRules bound keeps
// fibresToStrands' strand and is counted in `untapered`.
struct TaperedStrands {
    StrandSetF strands;
    long long tapered = 0, untapered = 0;
};
TaperedStrands fibresToTaperedStrands(const FibreSet& fs, const SceneMapping& map, int population, const StrandTaper& taper);
// The segments the export writes for the fibres with key below keyBelow (both populations), with the
// late taper if taper.tail > 0, counted exactly without keeping the strands: the taper inserts about
// 4 vertices per fibre at a tip of 0.3 and 8 at 0.1 (radius steps of 1.35).
long long taperedSegmentCount(const FibreSet& fs, const SceneMapping& map, const StrandTaper& taper, float keyBelow);
// The same late taper on strands already written (knit_check --retaper): one frozen export, tapered
// without knitting again. Each strand's root radius is its first vertex's; the emergence, which a
// strand file does not hold, is taken at its second vertex (the curve start), so the visible length
// is at most the root segment longer than the export's.
TaperedStrands taperStrandSet(const StrandSetF& s, const StrandTaper& taper);

// The late taper measured on the strands the export writes (fibres with key below keyBelow, both
// populations), from their stored float32 positions and radii (knit_check --taper-report). The
// error is |r(s) - r_linear(s)| between consecutive written vertices, its maximum found exactly
// (largestLinearError in the .cpp), against the profile of lateTaperRadius.
struct TaperReport {
    long long fibres = 0, tapered = 0;
    long long fallbackNoTail = 0;                      // no visible length or no tail: the generated strand
    long long fallbackRule[kStrandFaultCount] = {};    // the generated strand, by the rule the taper broke
    long long onsetExact = 0, onsetShifted = 0, onsetImpossible = 0;  // Onset, Refined (strands_file.cpp)
    long long nonconforming = 0;                       // Refined: max error above tolerance x r0
    std::vector<long long> nonconformingFirst;         // their fibre indices, at most 20
    long long originalSegments = 0, finalSegments = 0; // as generated, as written
    double maxErrorOverR0 = 0;
    double maxErrorOverLocalR = 0;                     // sampled (64 points and the extrema per interval)
    double maxEarlyThinningOverR0 = 0;                 // max of r0 - r_linear before x0
    double tailAreaExact = 0, tailAreaLinear = 0;      // sum over tails of the integral of r^2 (Simpson)
    double minRadius = 1e300, maxRadius = 0, maxStoredRatio = 1;
};
TaperReport measureTaper(const FibreSet& fs, const SceneMapping& map, const StrandTaper& taper, float keyBelow);
void addTaperReport(TaperReport& total, const TaperReport& part);

} // namespace knit
