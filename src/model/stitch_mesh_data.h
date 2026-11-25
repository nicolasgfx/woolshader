#pragma once
// A built stitch mesh as the app shows and saves it, and what its checks measured, as plain data: the viewer and the
// file writer see the polygons without seeing the mesher (stitchmesh is a private dependency of jobs).
//
// Generation and quality are two facts. A mesh is generated when the build produced polygon faces that can be drawn
// (stitchMeshUse: faces, finite coordinates, valid indices); its quality is what the existing checks measured (fit to the
// piece, coverage, boundary, relaxation, structure). A mesh with findings is still shown and saved; only a mesh that
// cannot be drawn is refused, with the reason. The reference viewer (Wu et al., stitchMeshing, convertStitchMesh2Rend and
// exportResult) shows and writes whatever was generated; this keeps that and reports the checks next to it.
#include <cstdint>
#include <string>
#include <vector>

#include "model/vec.h"

namespace knit {

// A face's stitch, as the mesher typed it (the order of stitchmesh's StitchFaceKind).
enum class StitchKind : uint8_t { Regular, Increase, Decrease, SpecialQuad, Unsupported, PoleStart };

// What the checks found on one face (bits of StitchMeshData::finding).
constexpr uint8_t kFaceOverFit = 1;          // farther from the piece than the fit target
constexpr uint8_t kFaceInverted = 2;         // wound against the piece's surface
constexpr uint8_t kFaceDegenerate = 4;       // of almost no area
constexpr uint8_t kFaceUntriangulable = 8;   // its polygon cannot be ear-clipped in its own plane
constexpr uint8_t kFaceUnsupported = 16;     // a cell without a knitting construction
constexpr uint8_t kFaceStructural = kFaceInverted | kFaceDegenerate | kFaceUntriangulable;

// The polygon faces, in the imported file's units and axes (exactly the coordinates the import read), not lifted off
// the surface for display.
struct StitchMeshData {
    std::vector<Vec3> pos;              // vertices
    std::vector<int> faceStart;         // faces + 1 entries: face f has the corners corner[faceStart[f]] .. corner[faceStart[f + 1] - 1]
    std::vector<int> corner;            // vertex indices, counter-clockwise seen from outside
    std::vector<StitchKind> kind;       // per face
    std::vector<uint8_t> finding;       // per face: kFace* bits
    std::vector<float> deviationMm;     // per face: centroid and edge midpoints to the piece [mm]
    std::vector<int> row;               // per face: its row (strip), counted from the start; -1 when unknown
    int faces() const { return faceStart.empty() ? 0 : static_cast<int>(faceStart.size()) - 1; }
    int sides(int f) const { return faceStart[f + 1] - faceStart[f]; }
};

// The build's own measurements, copied (nothing is measured again here). Limits are the ones the checks used.
struct StitchMeshQuality {
    bool relaxed = false;               // the build relaxed the mesh (ReviewDepth::Relaxed); false: an unrelaxed layout
    int faces = 0, rows = 0;
    int regular = 0, increase = 0, decrease = 0, special = 0, unsupported = 0, poleStart = 0;
    std::vector<int> facesPerRow;       // per row from the start: rings that grow and shrink
    // fit to the piece
    double fitTargetMm = 0, maxDeviationMm = 0;
    int facesOverFit = 0;
    // structure
    int inverted = 0, degenerate = 0, untriangulable = 0, intersectingPairs = 0;
    bool connectivityPass = false;
    std::string connectivityDetail;     // the failing or open connectivity checks ("" when it passes)
    // coverage and boundary
    double uncoveredFraction = 0, uncoveredLimit = 0;
    bool uncoveredMeasured = true;      // false: the piece lies within half a row of its border everywhere (not measurable)
    int coveredTwice = 0;
    double areaRatio = 0, areaLow = 0, areaHigh = 0;
    double boundaryMm = 0, boundaryLimitMm = 0;
    bool coveragePass = false;
    // relaxation
    bool relaxRan = false;
    std::string relaxStop;              // converged | stalled | step cap | time cap | diverged ("" when it did not run)
    std::string relaxNote;              // why it did not run
    bool clean = false;                 // the layout passes its checks (LayoutMeasures::checksPass, unchanged)
};

// Whether the mesh can be drawn and saved; reason says why not.
struct StitchMeshUse {
    bool drawable = false;
    std::string reason;
};

// The only gate for drawing and saving: faces exist, every coordinate is finite, every face has at least three corners
// and every index names a vertex. Fit, coverage or relaxation findings never make a mesh undrawable.
StitchMeshUse stitchMeshUse(const StitchMeshData& mesh);

// What the layout checks found, as lines for people, most serious first (structure, connectivity, fit, coverage,
// boundary); empty when the mesh passes them.
std::vector<std::string> stitchMeshFindings(const StitchMeshQuality& q);

// The relaxation's outcome when it did not converge or did not run ("" otherwise): shown next to the findings, but not one
// of the layout checks (LayoutMeasures::checksPass does not count it), so a mesh that passes its checks still reads clean.
std::string stitchMeshRelaxationNote(const StitchMeshQuality& q);

// One line: "Mesh generated · 563 faces · 25 rows · 117 faces exceed the 0.55 mm fit target", "... · passes its checks",
// or "No stitch mesh: <reason>".
std::string stitchMeshHeadline(const StitchMeshUse& use, const StitchMeshQuality& q);

// The axis-aligned bounds of the vertices (lo = hi = 0 without vertices).
void stitchMeshBounds(const StitchMeshData& mesh, Vec3& lo, Vec3& hi);

} // namespace knit
