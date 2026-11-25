#pragma once
// The built-in import's stages, 1 to 11, in the order wrap_surface.cpp runs them: the table of contents of the engine. Each
// stage has its own .cpp and its own tests (tests/tests_wrap.cpp, tests/tests_wrap_properties.cpp).
//
//   1-3  wrap_grid.cpp     the file's triangles -> voxels, sealed -> the blurred field
//   4-5  wrap_nets.cpp     surface nets, Taubin smoothing
//   6    wrap_remesh.cpp   the voxel remesh through OpenVDB, edges and winding
//   7    wrap_tidy.cpp     fragments dropped, one closed manifold kept
//   8    wrap_charts.cpp   a few wide charts
//   9    wrap_flatten.cpp  each chart flattened by SLIM (third_party/slim, MPL-2.0)
//   10   wrap_pack.cpp     charts turned and packed
//   11   wrap_bake.cpp     the colours baked by rays
#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "io/image_io.h"
#include "io/mesh_read.h"
#include "model/coordinates.h"
#include "model/mesh.h"

namespace knit {

// --- 1 to 3 (wrap_grid.cpp): from a triangle soup to the smooth field whose zero level is the outer surface. Each triangle is
// marked into a voxel grid, the marks are thickened by the closing so gaps seal, air is flooded in from the grid's border,
// whatever the air did not reach is the solid, and a blur makes the solid's staircase into a smooth field.

// The grid's voxels in C order: the x index runs slowest. Voxel (0, 0, 0) lies at `origin`.
struct GridPlacement {
    Vec3 origin;
    double voxelMm = 1;
    std::array<int, 3> size = {0, 0, 0};
    int64_t voxelCount() const { return int64_t(size[0]) * size[1] * size[2]; }
};

// The index of voxel (x, y, z) in a grid of `size`, in C order.
inline int64_t voxelIndex(const std::array<int, 3>& size, int x, int y, int z) {
    return (int64_t(x) * size[1] + y) * size[2] + z;
}

// One byte per voxel, 0 or 1.
struct SolidGrid {
    std::array<int, 3> size = {0, 0, 0};
    std::vector<uint8_t> voxels;
};

struct FieldGrid {
    std::array<int, 3> size = {0, 0, 0};
    std::vector<float> values;
};

using TriangleCorners = std::array<Vec3, 3>;

// The file's triangles in millimetres, Z up: a Y-up file turned by an exact quarter turn, (x, y, z) -> (x, -z, y). A quad of the
// file is split by quadSplit, then a triangle whose three vertices repeat an earlier triangle's is dropped. `fileTriangle` says
// which of the file's triangles each came from, and `vertex` its corners there, for reading its uv and colours.
struct WrapSource {
    std::vector<TriangleCorners> soup;
    std::vector<int> fileTriangle;
    std::vector<std::array<int, 3>> vertex;  // the file's vertex numbers
    std::vector<std::array<Vec2, 3>> uv;     // empty if the file has none
    std::vector<int> material;
};

WrapSource wrapSourceOf(const ImportedMesh& file, double unitMm, UpAxis up);

// A grid around the soup with `paddingVoxels` to spare on every side.
GridPlacement placeGrid(const std::vector<TriangleCorners>& soup, double voxelMm, int paddingVoxels);

// The voxels the triangles pass through: points on each triangle at most half a voxel apart, each rounded to its voxel.
SolidGrid markTriangles(const std::vector<TriangleCorners>& soup, const GridPlacement& grid);

// Dilation by a ball of `radiusVoxels`, wrapping around at the border. The grid's padding keeps that harmless.
SolidGrid dilateByBall(const SolidGrid& mask, int radiusVoxels);

// The voxels reachable from the grid's six faces through voxels that are not blocked.
SolidGrid floodFromGridBorder(const SolidGrid& blocked);

SolidGrid invertedGrid(const SolidGrid& mask);

// The solid as +1 and -1, blurred by a [1 2 1] / 4 filter along each axis `passes` times, wrapping at the border, with the
// six faces of the grid then set to -1 so that the surface closes inside it.
FieldGrid blurredSolidField(const SolidGrid& solid, int passes);

// --- 4, 5 (wrap_nets.cpp): the field's zero level as quads, then smoothed

// A mesh of quads, as the surface stages hand it on. Corners run counter-clockwise seen from outside. TriMesh holds triangles
// only, and the chart stage works on quads.
struct QuadMesh {
    std::vector<Vec3> pos;
    std::vector<std::array<int, 4>> quad;
    std::vector<std::array<int, 2>> edge;  // (lower vertex, higher vertex). Empty until numbered: seams refer to these numbers
};

// Surface nets: one vertex per grid cell the surface crosses, at the mean of its edge crossings, and one quad per crossed grid
// edge. Vertices come in the cells' C order, quads by axis and then in C order. Quads face out of the field's positive side.
QuadMesh extractSurfaceNets(const FieldGrid& field, const GridPlacement& grid);

// Taubin's low-pass filter over the mesh edges: a shrinking step and an inflating one per iteration, which removes the voxel
// terraces and keeps the volume (Taubin 1995).
void smoothTaubin(QuadMesh& mesh, int iterations);

// --- 6 (wrap_remesh.cpp): even quads through a level set (OpenVDB), then edges and a consistent outward winding

// The surface through an OpenVDB level set of voxel `voxelMm` and back to quads at the zero level, OpenVDB's documented defaults
// otherwise. Each connected part's quads are first made to agree with their neighbours, then turned so that the part's signed
// volume is positive: corners run counter-clockwise seen from outside. False with a message if the level set fails, gives
// triangles, or a part's volume is too small to tell its outside.
[[nodiscard]] bool remeshToEvenQuads(const QuadMesh& surface, double voxelMm, QuadMesh& remeshed, std::string* err);

// The two triangles of a quad, as corner positions 0 to 3, split on its shorter diagonal: (0, 1, 2) and (0, 2, 3), or (0, 1, 3)
// and (1, 2, 3). Equal diagonals take 0-2. A split whose triangle is degenerate or faces against the quad's Newell normal is
// not used when the other split is valid.
std::array<std::array<int, 3>, 2> quadSplit(const std::array<Vec3, 4>& corners);

// The mesh's unique edges as (lower vertex, higher vertex), sorted: the same list on every machine and at any thread count.
std::vector<std::array<int, 2>> sortedEdges(const std::vector<std::array<int, 4>>& quads);

// --- 7 (wrap_tidy.cpp): small loose fragments dropped; a large separate part or an open or branching edge stops the import

// The share of the area in the second largest connected part, 0 for a single part. The closing is widened while this is large.
double secondPartShare(const QuadMesh& mesh);

// Needs mesh.edge. Surviving vertices, faces and edges keep their order. What is left is one closed manifold, whose genus is
// returned.
[[nodiscard]] bool dropFragments(QuadMesh& mesh, int& genus, std::string* err);

// --- 8 (wrap_charts.cpp): a few wide charts on the quads of the closed surface. Charts grow from spread seeds over a graph of
// faces, each is opened into a disc, flattened by the caller, and judged. A chart that folds or overlaps is split, one that is
// only stretched is slit from a seam, and only the changed charts are flattened again. Chart placement decides how the knit
// looks, so ties are broken the same way every time: the lower index first.

constexpr double kCreaseCost = 4.0;      // chart borders prefer creases: crossing one is this much dearer
constexpr double kTurnCost = 8.0;        // a chart growing away from its mean direction: this much dearer

struct ChartSettings {
    int charts = 12;
    int maxCharts = 19;
    int maxSlits = 6;
    int maxRounds = 30;
};

struct FaceNeighbour {
    int face = 0;
    double cost = 0;
    int edge = 0;
};

// Which faces touch across which edge, with the cost of walking from face to face: distance, made dearer across a crease.
struct FaceGraph {
    const QuadMesh* mesh = nullptr;
    std::vector<double> area;
    std::vector<Vec3> normal;
    std::vector<int> start;                 // per face, where its neighbours begin; each face's list is in edge order
    std::vector<FaceNeighbour> neighbour;
    std::vector<std::array<int, 2>> facesOfEdge;
    std::vector<std::array<int, 4>> edgeOfCorner;  // per face, the edge from each corner to the next
    std::vector<int> edgeStart;             // per vertex, where its edges begin in edgeOfVertex
    std::vector<int> edgeOfVertex;
    std::vector<double> edgeLength;
    int faceCount() const { return static_cast<int>(start.size()) - 1; }
};

FaceGraph buildFaceGraph(const QuadMesh& mesh);

// `count` faces spread over the surface, each the farthest from those before it. The first is the farthest from face 0.
std::vector<int> spreadSeeds(const FaceGraph& graph, int count);

// Geodesic patches around the seeds, relaxed so that each seed sits in the middle of its patch. Returns each face's chart.
std::vector<int> growCharts(const FaceGraph& graph, std::vector<int> seeds);

// The seams, one flag per edge: every chart border, plus cuts that open each chart into a disc.
std::vector<char> cutToDiscs(const FaceGraph& graph, const std::vector<int>& label);

// One corner's uv per face corner, four per quad, in the mesh's order.
using QuadUv = std::vector<std::array<Vec2, 4>>;

// How each chart looks to knitting (src/knitting/layout/uv_atlas.h), measured on triangles.
struct ChartQuality {
    double totalOutside = 0;             // share of the surface outside the stretch bounds
    std::vector<int> chart;              // the chart numbers in use, ascending
    std::vector<double> stretched;       // per entry of `chart`: share outside the bounds
    std::vector<double> unusable;        // per entry of `chart`: share folded or overlapping
    std::vector<double> worst;           // per triangle: over 1 means outside the bounds, infinite if it cannot be measured
    std::vector<int> faceOfTriangle;
};

ChartQuality measureCharts(const QuadMesh& mesh, const std::vector<int>& label, const QuadUv& uv);

// The state of the loop between two flattenings.
struct ChartRound {
    std::vector<int> label;
    std::vector<char> seam;
    std::vector<char> slit;
    int slitCount = 0;
};

enum class ChartDecision { FlattenChanged, Finished };

struct ChartStep {
    ChartDecision decision = ChartDecision::Finished;
    std::vector<char> facesToFlatten;  // with FlattenChanged, one flag per face
};

// One round's decision after measuring: split the worst broken chart, or stop, or slit the most stretched chart.
ChartStep decideNextRound(const FaceGraph& graph, const ChartQuality& quality, const ChartSettings& settings, ChartRound& round);

// --- 9 (wrap_flatten.cpp): faces split along the seams into charts, each flattened with least distortion by SLIM (Rabinovich et
// al. 2017; the solver is third_party/slim, MPL-2.0). A chart must be a disc: connected, every edge used once or twice, one
// boundary cycle. Its rest shape is made safe by intrinsic mollification (Sharp and Crane 2020), its start map is Tutte's with
// mean value weights on a circle (Floater 2003) scaled to the chart's area, and SLIM (Rabinovich et al. 2017) runs from there
// until the energy stops falling.

// Why a chart's solver stopped. Converged says only that the energy stopped falling; whether the map is good enough is judged
// by measureCharts and the chart rules.
enum class SolveStop {
    Converged,        // one iteration lowered the energy by less than the tolerance share of it
    IterationCap,     // the iteration cap was reached first
    SolverFailed,     // a global solve failed: the chart keeps its last flip-free map
    NonFiniteEnergy,  // the map's energy is not finite (a start map degenerate in floating point): the last finite map is kept
};

struct ChartSolve {
    int iterations = 0;
    SolveStop stop = SolveStop::Converged;
    double energy = 0;  // the symmetric Dirichlet energy per unit rest area of the kept map (4 at an isometry)
};

struct FlattenReport {
    int charts = 0;
    int stoppedEarly = 0;  // charts stopped by SolverFailed or NonFiniteEnergy
    int reachedCap = 0;    // charts stopped by IterationCap
    std::vector<ChartSolve> solves;  // per chart, in chart order
};

// Flattens the flagged faces and leaves the uv of the others alone. All charts of a call share one scale, about their own
// centres. False with a message naming the chart if a chart is not a disc, or has no positive length or area.
[[nodiscard]] bool flattenFaces(const FaceGraph& graph, const std::vector<char>& seam, const std::vector<char>& faces, QuadUv& uv,
                                FlattenReport& report, std::string* err);

// For tests only: while set, every chart's solver is handed a system it cannot solve, so the failure path can be checked.
void breakFlattenSolverForTests(bool broken);

// For tests only: the solver's iteration cap while positive; 0 restores the built-in cap.
void capFlattenIterationsForTests(int cap);

// --- 10 (wrap_pack.cpp): each flattened chart turned so its smallest bounding rectangle stands with the long side along v, then
// all laid out in the unit square at one scale

// Packs the flagged faces' charts. A chart is what the seams and the flags leave connected. The uv of other faces stays.
void packCharts(const FaceGraph& graph, const std::vector<char>& seam, const std::vector<char>& faces, QuadUv& uv);

// --- 11 (wrap_bake.cpp): the source's colours baked into the charts' texture. Each texel on the surface shoots a ray from outside
// along the surface's normal, and takes the colour where it first meets the source, or at the closest point of the source if it
// meets nothing there (a sealed gap). Colours are as stored in the file, sRGB values, nothing invented.

constexpr int kBakedTextureSide = 1024;  // about 0.5 mm a texel, a tenth of a stitch

// Rows run with v: row 0 lies at v = 0.
struct BakedColours {
    int side = kBakedTextureSide;
    std::vector<float> rgb;     // sRGB values, three per texel
    std::vector<char> filled;   // per texel: lies on a chart
    long long closestTexels = 0;  // those that took the closest surface
    long long texels = 0;
};

// A bounding volume hierarchy over a triangle soup: median split on the longest axis of the triangle centres, leaves of 4. It
// answers the first hit along a ray (Moller and Trumbore 1997) and the closest point (Ericson 2005, 5.1.5). It keeps a reference
// to `soup`, which must outlive it.
class TriangleBvh {
public:
    explicit TriangleBvh(const TriMesh& soup);

    struct Hit {
        int tri = -1;
        double t = 0, b1 = 0, b2 = 0;  // the point is corner 0 + b1 (corner 1 - corner 0) + b2 (corner 2 - corner 0)
    };

    struct Nearest {
        int tri = -1;  // -1: an empty soup or a query that is not finite
        double distSq = 0, b1 = 0, b2 = 0;  // the point is corner 0 + b1 (corner 1 - corner 0) + b2 (corner 2 - corner 0)
    };

    // The nearest triangle met along origin + t direction, 0 <= t <= maxT.
    Hit firstHit(const Vec3& origin, const Vec3& direction, double maxT) const;

    // The closest point of the soup to `q`. Equal distances go to the lower triangle index.
    Nearest closestPoint(const Vec3& q) const;

private:
    struct Node {
        Vec3 lo, hi;
        int first = 0, count = 0;  // a leaf: order_[first .. first + count)
        int right = 0;             // an inner node: children are this node + 1 and `right`
    };

    int build(int begin, int end);
    void intersect(int t, const Vec3& origin, const Vec3& direction, double& limit, Hit& best) const;
    void closestOnTriangleOfSoup(int t, const Vec3& q, Nearest& best) const;

    const TriMesh& soup_;
    std::vector<int> order_;
    std::vector<Vec3> centre_;
    std::vector<Node> nodes_;
};

// Unit vertex normals: each face's normal weighted by its corner angle at the vertex (Thuermer and Wuethrich 1998). A vertex whose
// weights sum to zero takes the normal of its largest face.
std::vector<Vec3> vertexNormalsOf(const FaceGraph& graph);

BakedColours bakeSurfaceColours(const FaceGraph& graph, const QuadUv& uv, const WrapSource& source, const ImportedMesh& file,
                                const std::vector<ColourImage>& textures, double rayMm);

// Grows the colours past the chart borders, so that no lookup lands on empty texture. Returns the new sRGB values.
std::vector<float> bledColours(const BakedColours& colours, int rounds);

// The mean linear colour of the filled texels.
Vec3 meanLinearColourOf(const BakedColours& colours);

} // namespace knit
