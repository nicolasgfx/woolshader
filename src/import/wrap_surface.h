#pragma once
// The import's rebuild from start to end: a source file becomes three files in a folder (wrapped.obj, wrapped.mtl,
// wrapped_basecolor.png), and the application reads them back through the ordinary reader (readRebuiltSurface).
//
//   the file's triangles -> voxels, sealed -> blurred field -> surface nets -> Taubin -> remesh -> tidy -> charts, flattened and packed
//   -> colours baked -> wrapped.obj, wrapped.mtl, wrapped_basecolor.png
#include <cstdint>
#include <string>
#include <vector>

#include "import/prepare_import.h"
#include "import/wrap_stages.h"

namespace knit {

constexpr double kMaxWrapVoxels = 64000000; // about 400 mm cubed at 1 mm
constexpr const char* kRebuiltSurfaceName = "wrapped";  // the rebuild writes wrapped.obj, wrapped.mtl, wrapped_basecolor.png

// The recipe: one set of values, the lengths scaled by the stitch width.
struct WrapRecipe {
    double voxelMm = 1.0;
    double closingMm = 8.0;
    double offsetMm = 1.0;
    double remeshMm = 1.5;
    double rayMm = 5.0;
    int blurPasses = 14;
    int taubinIterations = 25;
    ChartSettings charts;
};

// The values for a stitch width, with the voxel enlarged until the grid stays within kMaxWrapVoxels. The blur keeps its width in mm.
WrapRecipe wrapRecipeFor(double stitchMm, const std::vector<TriangleCorners>& soup);

// The closed exterior of the soup, as tidy quads with numbered edges. The closing widens while a large part stays separate. Lines for
// the import's change list go to `changes`.
[[nodiscard]] bool wrapSoup(const std::vector<TriangleCorners>& soup, const WrapRecipe& recipe, QuadMesh& surface,
                            std::vector<std::string>& changes, std::string* err, const ImportProgress& progress);

// A few wide charts on the surface, flattened and packed: each corner's uv. False with the flattening's message if a chart cannot
// be flattened.
[[nodiscard]] bool chartSurface(const FaceGraph& graph, const ChartSettings& settings, QuadUv& uv, std::vector<std::string>& changes,
                                std::string* err, const ImportProgress& progress);

// The surface as an OBJ the application's reader takes back: millimetres, Y up, one uv per quad corner, each quad as the two
// triangles of quadSplit (the ones the flatten measured), numbers in their shortest exact form.
[[nodiscard]] bool writeSurfaceObj(const std::string& path, const QuadMesh& mesh, const QuadUv& uv);

// Runs everything and writes wrapped.obj, wrapped.mtl and wrapped_basecolor.png into `folder`.
[[nodiscard]] bool wrapFileToFolder(const ImportedMesh& file, const std::vector<ColourImage>& textures, const ImportOptions& options,
                                    const std::string& folder, std::vector<std::string>& changes, std::string* err,
                                    const ImportProgress& progress);

} // namespace knit
