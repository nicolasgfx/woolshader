#pragma once
// What the Fibres stage makes (fibre_stage.h, fur_stage.h, or saved fibres, fibre_files.h): one
// PartFibres per part, in the packed form the device, the viewer and the export all read, and what the
// export and the display read from it: the strands as written, the segments kept by the budget.
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "fibres/fibre.h"
#include "fibres/fibre_strands.h"
#include "scene/scene.h"

namespace knit {

// One part's fibres, in the one packed form the device, the viewer and the export all read.
struct PartFibres {
    int part = -1;                  // index into Scene::parts
    std::string name;               // the part's name: the export's "<name>_fibres_oriented" / "_random" (fur: "_fibres")
    FibreSet fs;                    // library mm relative to the part's centre, before the segment budget's key threshold
    // Fur: per fibre the texture colour at its root (0xRRGGBB, reduced to the scene's palette,
    // stitch_colours.h); empty: wool, or a part without a texture. Whatever material.fur_colour_source
    // says: the source is an Export setting, shown and written without a new fur job.
    std::vector<uint32_t> fibreColour;
    // The segment budget (fibres.segment_budget) holds per object (sceneObjectOf: a bear, an import):
    // the export writes this part's fibres with key below keyBelow, the threshold of its object's
    // fibres together (applyObjectBudgets), whole fibres, never shortened. A second bear leaves the
    // first one's fibres as they were.
    float keyBelow = 2.0f;
    float coveredBelow = 2.0f;      // every candidate with key below it was evaluated (2: all of them)
};

// One population as the export writes it: scene units (Y-up metres), with each strand's budget key.
// With a late taper (taper.tail > 0) the counts say how many fibres took it (StrandTaper).
struct PartStrands {
    StrandSetF strands;
    std::vector<float> keys;
    long long tapered = 0, untapered = 0;
};
PartStrands partStrands(const PartFibres& fibres, const Vec3& centreMm, FibrePopulation population, const StrandTaper& taper);
// One part's fibre segments below its keyBelow: as grown, and as the export writes them (the late
// taper's vertices included), without building the strands.
struct KeptSegments {
    long long kept = 0, exported = 0;
};
KeptSegments keptSegments(const PartFibres& fibres, const Vec3& centreMm, const StrandTaper& taper);
// The late taper of one part's exported fibres, measured (measureTaper, knit_check --taper-report).
TaperReport partTaperReport(const PartFibres& fibres, const Vec3& centreMm, const StrandTaper& taper);

// Adds every fibre's (key, segments), what a segment budget's threshold is decided on.
void addKeySegments(const FibreSet& fibres, std::vector<std::pair<float, int>>& out);
long long segmentCount(const FibreSet& fibres);
// Sets every part's keyBelow: the threshold of the segment budget over its object's parts together
// (budgetThreshold). Returns one log line per object.
std::string applyObjectBudgets(const Scene& scene, long long budget, std::vector<PartFibres>& parts);

// The scene's fibres exactly as the export makes them (every part's candidates in one budget order, or
// full generation; the export's model ring and checks), kept so the export writes the fibres the user
// reviewed.
struct SceneFibres {
    std::vector<PartFibres> parts;
    bool budgetMode = false;       // wool: only the candidates below the budget's key were evaluated
    std::string log;               // the generation's log lines (as the export prints them)
};

} // namespace knit
