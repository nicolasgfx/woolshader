// The built-in import's stages, one after the other. Each stage has its own file; here they meet, and the result is written as
// the three files the application reads back (wrapped.obj, wrapped.mtl, wrapped_basecolor.png).
#include "import/wrap_surface.h"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>

#include "import/wrap_stages.h"
#include "io/cache_paths.h"
#include "io/image_io.h"
#include "io/files.h"
#include "model/cancel.h"

namespace knit {

namespace {

constexpr double kReferenceStitchMm = 5.5;  // the stitch width the recipe's millimetres were found at
constexpr int kClosingRetries = 2;
constexpr double kClosingGrowth = 1.5;
constexpr double kSeparatePartShare = 0.01;  // a second part above this share of the area is no fragment
constexpr int kBleedRounds = 4;

void addChange(std::vector<std::string>& changes, const char* format, double a, double b = 0, double c = 0, double d = 0, double e = 0, double f = 0, double g = 0,
               double h = 0, double i = 0) {
    char line[288];
    std::snprintf(line, sizeof line, format, a, b, c, d, e, f, g, h, i);
    changes.push_back(line);
}

// One wrap at one closing width.
bool wrapOnce(const std::vector<TriangleCorners>& soup, const WrapRecipe& recipe, double closingMm, QuadMesh& surface, std::string* err) {
    const int closing = static_cast<int>(std::nearbyint(closingMm / recipe.voxelMm));
    const int offset = static_cast<int>(std::nearbyint(recipe.offsetMm / recipe.voxelMm));
    const GridPlacement grid = placeGrid(soup, recipe.voxelMm, closing + offset + recipe.blurPasses + 4);
    const SolidGrid marked = markTriangles(soup, grid);
    const SolidGrid blocked = dilateByBall(marked, closing);
    const SolidGrid air = floodFromGridBorder(blocked);
    const SolidGrid solid = dilateByBall(invertedGrid(dilateByBall(air, closing)), offset);
    const FieldGrid field = blurredSolidField(solid, recipe.blurPasses);
    QuadMesh nets = extractSurfaceNets(field, grid);
    smoothTaubin(nets, recipe.taubinIterations);
    if (!remeshToEvenQuads(nets, recipe.remeshMm, surface, err)) return false;
    surface.edge = sortedEdges(surface.quad);
    return true;
}

uint8_t byteOf(float value) {
    if (value <= 0.0f) return 0;
    if (value > 1.0f - 0.5f / 255.0f) return 255;
    return static_cast<uint8_t>(255.0f * value + 0.5f);
}

// Rows run from the image's top in a file, from v = 0 in the bake.
bool writeBaseColourPng(const std::string& path, const std::vector<float>& rgb, int side) {
    std::vector<uint8_t> bytes(static_cast<size_t>(side) * side * 3);
    for (int y = 0; y < side; ++y)
        for (int x = 0; x < side * 3; ++x) bytes[static_cast<size_t>(y) * side * 3 + x] = byteOf(rgb[static_cast<size_t>(side - 1 - y) * side * 3 + x]);
    return writePngRgb8(path, side, side, bytes);
}

// A number in its shortest form that reads back to the same double, independent of the locale.
void appendNumber(std::string& text, double value) {
    char digits[32];
    const std::to_chars_result written = std::to_chars(digits, digits + sizeof digits, value);
    text.append(digits, written.ptr);
}

// The flagged faces' charts flattened, then packed.
bool flattenAndPack(const FaceGraph& graph, const std::vector<char>& seam, const std::vector<char>& faces, QuadUv& uv, FlattenReport& report,
                    std::string* err) {
    if (!flattenFaces(graph, seam, faces, uv, report, err)) return false;
    packCharts(graph, seam, faces, uv);
    return true;
}

bool writeMaterial(const std::string& path, const Vec3& meanColour, const std::string& sourceName) {
    char kd[128];
    std::snprintf(kd, sizeof kd, "Kd %.6f %.6f %.6f\n", meanColour.x, meanColour.y, meanColour.z);
    std::ofstream file(path);
    file << "# converted from " << sourceName << ", by the built-in import\nnewmtl " << kRebuiltSurfaceName << "\n" << kd << "map_Kd "
         << kRebuiltSurfaceName << "_basecolor.png\n";
    file.close();
    return !file.fail();
}

} // namespace

WrapRecipe wrapRecipeFor(double stitchMm, const std::vector<TriangleCorners>& soup) {
    WrapRecipe recipe;
    const double scale = stitchMm / kReferenceStitchMm;
    recipe.voxelMm *= scale;
    recipe.closingMm *= scale;
    recipe.offsetMm *= scale;
    recipe.remeshMm *= scale;
    recipe.rayMm *= scale;

    Vec3 low(1e300, 1e300, 1e300), high(-1e300, -1e300, -1e300);
    for (const TriangleCorners& t : soup)
        for (const Vec3& p : t) {
            low = {std::min(low.x, p.x), std::min(low.y, p.y), std::min(low.z, p.z)};
            high = {std::max(high.x, p.x), std::max(high.y, p.y), std::max(high.z, p.z)};
        }
    const Vec3 extent = high - low;
    const double padded = 2.0 * (recipe.closingMm + recipe.offsetMm + (recipe.blurPasses + 4) * recipe.voxelMm);
    const double wanted = recipe.voxelMm;
    double voxel = wanted;
    while ((extent.x + padded) / voxel * ((extent.y + padded) / voxel) * ((extent.z + padded) / voxel) > kMaxWrapVoxels) voxel *= 1.05;
    const double ratio = wanted / voxel;
    recipe.blurPasses = static_cast<int>(std::nearbyint(recipe.blurPasses * (ratio * ratio)));
    recipe.voxelMm = voxel;
    return recipe;
}

bool wrapSoup(const std::vector<TriangleCorners>& soup, const WrapRecipe& recipe, QuadMesh& surface, std::vector<std::string>& changes,
              std::string* err, const ImportProgress& progress) {
    double closingMm = recipe.closingMm;
    for (int attempt = 0;; ++attempt) {
        if (progress) progress(attempt == 0 ? "wrapping the surface" : "wrapping the surface again, closing wider", 0.02);
        if (!wrapOnce(soup, recipe, closingMm, surface, err)) return false;
        if (surface.quad.empty()) return fail(err, "nothing is left of the object after the blur: it is too small for this stitch width");
        if (attempt == kClosingRetries || secondPartShare(surface) <= kSeparatePartShare) break;
        addChange(changes, "separate parts at a closing of %.3g mm: closing again, wider", closingMm);
        closingMm *= kClosingGrowth;
    }
    int genus = 0;
    if (!dropFragments(surface, genus, err)) return false;
    addChange(changes, "wrapped at a voxel of %.3g mm: %.0f quads, one closed manifold of genus %.0f", recipe.voxelMm, double(surface.quad.size()), double(genus));
    return true;
}

bool chartSurface(const FaceGraph& graph, const ChartSettings& settings, QuadUv& uv, std::vector<std::string>& changes, std::string* err,
                  const ImportProgress& progress) {
    if (progress) progress("cutting and flattening charts", 0.4);
    const QuadMesh& mesh = *graph.mesh;
    ChartRound round;
    round.label = growCharts(graph, spreadSeeds(graph, settings.charts));
    round.seam = cutToDiscs(graph, round.label);
    round.slit.assign(round.seam.size(), 0);
    uv.assign(mesh.quad.size(), {});
    std::vector<char> faces(mesh.quad.size(), 1);
    FlattenReport report;  // of the latest flattening
    if (!flattenAndPack(graph, round.seam, faces, uv, report, err)) return false;
    for (int r = 0; r < settings.maxRounds; ++r) {
        cancelCheckpoint();
        const ChartQuality quality = measureCharts(mesh, round.label, uv);
        const ChartStep step = decideNextRound(graph, quality, settings, round);
        if (step.decision == ChartDecision::Finished) break;
        if (!flattenAndPack(graph, round.seam, step.facesToFlatten, uv, report, err)) return false;
    }
    // all charts again: one scale, packed; this flattening's map is the one written, so its report is the one shown
    if (!flattenAndPack(graph, round.seam, faces, uv, report, err)) return false;

    const ChartQuality quality = measureCharts(mesh, round.label, uv);
    std::vector<double> areaOfChart(*std::max_element(round.label.begin(), round.label.end()) + 1, 0.0);
    for (size_t f = 0; f < mesh.quad.size(); ++f) areaOfChart[round.label[f]] += graph.area[f];
    double smallest = 1e300, largest = 0;
    for (double area : areaOfChart)
        if (area > 0) {
            smallest = std::min(smallest, area);
            largest = std::max(largest, area);
        }
    const double worstUnusable = *std::max_element(quality.unusable.begin(), quality.unusable.end());
    addChange(changes, "%.0f charts, %.0f slits, %.1f %% outside the stretch bounds, worst fold or overlap %.1f %%, chart area %.0f .. %.0f mm^2, %.0f seam edges, %.0f charts stopped early, %.0f at the solver's iteration cap",
              double(quality.chart.size()), double(round.slitCount), 100.0 * quality.totalOutside, 100.0 * worstUnusable, smallest, largest,
              double(std::count(round.seam.begin(), round.seam.end(), char(1))), double(report.stoppedEarly), double(report.reachedCap));
    return true;
}

bool writeSurfaceObj(const std::string& path, const QuadMesh& mesh, const QuadUv& uv) {
    std::string text = std::string("# wrapped by the built-in import\nmtllib ") + kRebuiltSurfaceName + ".mtl\no " + kRebuiltSurfaceName + "\n";
    text.reserve(text.size() + mesh.pos.size() * 64 + mesh.quad.size() * 200);
    for (const Vec3& p : mesh.pos) {  // millimetres, Y up: (x, y, z) -> (x, z, -y)
        text += "v ";
        appendNumber(text, p.x);
        text += ' ';
        appendNumber(text, p.z);
        text += ' ';
        appendNumber(text, -p.y);
        text += '\n';
    }
    for (const std::array<Vec2, 4>& corners : uv)
        for (const Vec2& c : corners) {
            text += "vt ";
            appendNumber(text, c.x);
            text += ' ';
            appendNumber(text, c.y);
            text += '\n';
        }
    text += std::string("usemtl ") + kRebuiltSurfaceName + "\n";
    // Each quad as the two triangles of quadSplit on this very surface, the split the flatten, the chart measure and the bake
    // used: a reader triangulating a quad on its own may take the other diagonal, which can fold in uv.
    for (size_t f = 0; f < mesh.quad.size(); ++f) {
        const std::array<int, 4>& q = mesh.quad[f];
        const std::array<Vec3, 4> corners = {mesh.pos[q[0]], mesh.pos[q[1]], mesh.pos[q[2]], mesh.pos[q[3]]};
        for (const std::array<int, 3>& triangle : quadSplit(corners)) {
            text += 'f';
            for (int k : triangle) {
                text += ' ';
                text += std::to_string(q[k] + 1);
                text += '/';
                text += std::to_string(4 * f + k + 1);
            }
            text += '\n';
        }
    }
    std::ofstream file(path, std::ios::binary);
    file.write(text.data(), static_cast<std::streamsize>(text.size()));
    file.close();
    return !file.fail();
}

bool wrapFileToFolder(const ImportedMesh& file, const std::vector<ColourImage>& textures, const ImportOptions& options, const std::string& folder,
                      std::vector<std::string>& changes, std::string* err, const ImportProgress& progress) {
    if (!(options.unitMm > 0) || !std::isfinite(options.unitMm)) return fail(err, "the unit size must be a positive number of millimetres");
    if (!(options.wrapStitchMm > 0) || !std::isfinite(options.wrapStitchMm)) return fail(err, "the stitch width must be a positive number of millimetres");
    if (file.tri.empty()) return fail(err, "the file has no faces");

    if (progress) progress("reading the file", 0);
    const WrapSource source = wrapSourceOf(file, options.unitMm, options.up);
    const WrapRecipe recipe = wrapRecipeFor(options.wrapStitchMm, source.soup);
    QuadMesh surface;
    if (!wrapSoup(source.soup, recipe, surface, changes, err, progress)) return false;
    const FaceGraph graph = buildFaceGraph(surface);
    QuadUv uv;
    if (!chartSurface(graph, recipe.charts, uv, changes, err, progress)) return false;

    if (progress) progress("baking colours", 0.92);
    const BakedColours colours = bakeSurfaceColours(graph, uv, source, file, textures, recipe.rayMm);
    addChange(changes, "colours baked: %.0f texels, %.1f %% took the closest surface", double(colours.texels),
              colours.texels ? 100.0 * colours.closestTexels / colours.texels : 0.0);

    if (progress) progress("writing the surface", 0.98);
    const std::string stem = (std::filesystem::path(folder) / kRebuiltSurfaceName).string();
    if (!writeSurfaceObj(stem + ".obj", surface, uv)) return fail(err, "cannot write " + stem + ".obj");
    if (!writeBaseColourPng(stem + "_basecolor.png", bledColours(colours, kBleedRounds), colours.side)) return fail(err, "cannot write the base colour texture");
    if (!writeMaterial(stem + ".mtl", meanLinearColourOf(colours), std::filesystem::path(file.path).filename().string())) return fail(err, "cannot write " + stem + ".mtl");
    return true;
}

} // namespace knit
