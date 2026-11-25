// The properties of the built-in import's stages (import/wrap_stages.h), part of knit_tests_import: quadSplit, the file's
// triangles, the level set remesh and its winding, the edge list, the packer, the closest point, the vertex normals, the surface
// file, and the whole wrap of the fixtures in tests/data/wrap_fixtures. Small generated meshes, on the CPU. The flatten has its
// own file, tests_wrap_flatten.cpp.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <set>
#include <string>
#include <vector>

#include "harness.h"
#include "import/import_pipeline.h"
#include "import/wrap_stages.h"
#include "import/wrap_surface.h"
#include "io/mesh_read.h"
#include "knitting/layout/uv_atlas.h"
#include "model/mesh.h"
#include "model/parallel.h"

using namespace knit;
namespace fs = std::filesystem;

namespace {

std::string freshDir(const std::string& name) {
    const fs::path p = fs::temp_directory_path() / "knit_wrap_properties" / name;
    std::error_code ec;
    fs::remove_all(p, ec);
    fs::create_directories(p, ec);
    return p.string();
}

std::string readBytes(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    std::stringstream bytes;
    bytes << file.rdbuf();
    return bytes.str();
}

bool sameBits(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }

using Split = std::array<std::array<int, 3>, 2>;
const Split kAlongFirst = {{{0, 1, 2}, {0, 2, 3}}};
const Split kAlongSecond = {{{0, 1, 3}, {1, 2, 3}}};

std::array<Vec3, 4> flatQuad(const Vec2& a, const Vec2& b, const Vec2& c, const Vec2& d) {
    return {Vec3{a.x, a.y, 0}, Vec3{b.x, b.y, 0}, Vec3{c.x, c.y, 0}, Vec3{d.x, d.y, 0}};
}

// A cube of side `side` from the origin, its quads counter-clockwise seen from outside, its edges numbered.
QuadMesh cubeMesh(double side) {
    QuadMesh mesh;
    mesh.pos = {{0, 0, 0}, {side, 0, 0}, {side, side, 0}, {0, side, 0}, {0, 0, side}, {side, 0, side}, {side, side, side}, {0, side, side}};
    mesh.quad = {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {3, 7, 6, 2}, {0, 4, 7, 3}, {1, 2, 6, 5}};
    mesh.edge = sortedEdges(mesh.quad);
    return mesh;
}

// Six times the signed volume of a closed quad mesh, each quad as (0, 1, 2) and (0, 2, 3).
double sixTimesVolume(const QuadMesh& mesh) {
    double sum = 0;
    for (const std::array<int, 4>& q : mesh.quad)
        sum += dot(mesh.pos[q[0]], cross(mesh.pos[q[1]], mesh.pos[q[2]])) + dot(mesh.pos[q[0]], cross(mesh.pos[q[2]], mesh.pos[q[3]]));
    return sum;
}

// Every edge used by exactly two quads, once in each direction.
bool everyEdgeTwiceOpposite(const QuadMesh& mesh) {
    std::map<std::pair<int, int>, int> directed;
    for (const std::array<int, 4>& q : mesh.quad)
        for (int k = 0; k < 4; ++k) ++directed[{q[k], q[(k + 1) % 4]}];
    for (const std::pair<const std::pair<int, int>, int>& entry : directed) {
        if (entry.second != 1) return false;
        const std::map<std::pair<int, int>, int>::const_iterator reverse = directed.find({entry.first.second, entry.first.first});
        if (reverse == directed.end() || reverse->second != 1) return false;
    }
    return true;
}

// The number of connected parts of the quads, joined through shared vertices.
int partCount(const QuadMesh& mesh) {
    std::vector<int> parent(mesh.pos.size());
    std::iota(parent.begin(), parent.end(), 0);
    std::vector<char> used(mesh.pos.size(), 0);
    for (const std::array<int, 4>& q : mesh.quad)
        for (int k = 0; k < 4; ++k) {
            used[q[k]] = 1;
            int a = q[0], b = q[k];
            while (parent[a] != a) a = parent[a];
            while (parent[b] != b) b = parent[b];
            parent[std::max(a, b)] = std::min(a, b);
        }
    int parts = 0;
    for (size_t v = 0; v < parent.size(); ++v) parts += used[v] && parent[v] == static_cast<int>(v);
    return parts;
}

// A face graph of faces with no neighbours: what packCharts reads of separate quads.
FaceGraph unlinkedFaces(const QuadMesh& mesh) {
    FaceGraph graph;
    graph.mesh = &mesh;
    graph.start.assign(mesh.quad.size() + 1, 0);
    return graph;
}

double signedArea(const Vec2& a, const Vec2& b, const Vec2& c) { return 0.5 * cross(b - a, c - a); }

// Whether two triangles overlap with some area. Touching edges and shared corners do not count: only edges that cross properly,
// or a corner strictly inside the other triangle.
bool pointStrictlyInside(const Vec2& p, const std::array<Vec2, 3>& t) {
    const double s0 = cross(t[1] - t[0], p - t[0]), s1 = cross(t[2] - t[1], p - t[1]), s2 = cross(t[0] - t[2], p - t[2]);
    return (s0 > 0 && s1 > 0 && s2 > 0) || (s0 < 0 && s1 < 0 && s2 < 0);
}

bool segmentsCross(const Vec2& a, const Vec2& b, const Vec2& c, const Vec2& d) {
    const double d1 = cross(b - a, c - a), d2 = cross(b - a, d - a), d3 = cross(d - c, a - c), d4 = cross(d - c, b - c);
    return ((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0));
}

bool trianglesOverlap(const std::array<Vec2, 3>& a, const std::array<Vec2, 3>& b) {
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            if (segmentsCross(a[i], a[(i + 1) % 3], b[j], b[(j + 1) % 3])) return true;
    for (int i = 0; i < 3; ++i)
        if (pointStrictlyInside(a[i], b) || pointStrictlyInside(b[i], a)) return true;
    return false;
}

// Whether triangles of two different islands overlap. Triangles are binned by their boxes in a grid over [0, 1]^2, and only
// pairs sharing a cell are tested.
bool islandsOverlap(const std::vector<std::array<Vec2, 3>>& triangles, const std::vector<int>& islandOfTriangle) {
    constexpr int kCells = 64;
    std::vector<std::vector<int>> cell(kCells * kCells);
    for (size_t t = 0; t < triangles.size(); ++t) {
        double lowU = 1e300, lowV = 1e300, highU = -1e300, highV = -1e300;
        for (const Vec2& p : triangles[t]) {
            lowU = std::min(lowU, p.x), highU = std::max(highU, p.x), lowV = std::min(lowV, p.y), highV = std::max(highV, p.y);
        }
        const int u0 = std::clamp(static_cast<int>(lowU * kCells), 0, kCells - 1), u1 = std::clamp(static_cast<int>(highU * kCells), 0, kCells - 1);
        const int v0 = std::clamp(static_cast<int>(lowV * kCells), 0, kCells - 1), v1 = std::clamp(static_cast<int>(highV * kCells), 0, kCells - 1);
        for (int u = u0; u <= u1; ++u)
            for (int v = v0; v <= v1; ++v) cell[u * kCells + v].push_back(static_cast<int>(t));
    }
    for (const std::vector<int>& members : cell)
        for (size_t i = 0; i < members.size(); ++i)
            for (size_t j = i + 1; j < members.size(); ++j)
                if (islandOfTriangle[members[i]] != islandOfTriangle[members[j]] && trianglesOverlap(triangles[members[i]], triangles[members[j]]))
                    return true;
    return false;
}

// Separate quads whose uv are rectangles `size` turned by `angle` about `centre`: the packer's islands of known shape.
struct RectangleIslands {
    QuadMesh mesh;
    QuadUv uv;
};

RectangleIslands rectangleIslands(const std::vector<Vec2>& size, const std::vector<double>& angle, const std::vector<Vec2>& centre) {
    RectangleIslands r;
    for (size_t i = 0; i < size.size(); ++i) {
        const int first = static_cast<int>(r.mesh.pos.size());
        r.mesh.pos.insert(r.mesh.pos.end(), {{0, 0, 0}, {size[i].x, 0, 0}, {size[i].x, size[i].y, 0}, {0, size[i].y, 0}});
        r.mesh.quad.push_back({first, first + 1, first + 2, first + 3});
        const double c = std::cos(angle[i]), s = std::sin(angle[i]);
        std::array<Vec2, 4> corners;
        const Vec2 local[4] = {{-0.5 * size[i].x, -0.5 * size[i].y}, {0.5 * size[i].x, -0.5 * size[i].y}, {0.5 * size[i].x, 0.5 * size[i].y},
                               {-0.5 * size[i].x, 0.5 * size[i].y}};
        for (int k = 0; k < 4; ++k) corners[k] = {centre[i].x + c * local[k].x - s * local[k].y, centre[i].y + s * local[k].x + c * local[k].y};
        r.uv.push_back(corners);
    }
    return r;
}

const std::vector<Vec2> kFiveSizes = {{3, 1}, {1, 2.5}, {2, 0.5}, {4, 3}, {1.5, 1.2}};

RectangleIslands fiveRectangles() {
    return rectangleIslands(kFiveSizes, {0.3, 1.2, -0.7, 2.5, 0.0}, {{5, 5}, {-3, 2}, {0, 0}, {10, -4}, {1, 1}});
}

void packAll(const RectangleIslands& r, QuadUv& uv) {
    const FaceGraph graph = unlinkedFaces(r.mesh);
    const std::vector<char> seam;
    const std::vector<char> faces(r.mesh.quad.size(), 1);
    uv = r.uv;
    packCharts(graph, seam, faces, uv);
}

// One random triangle soup in a 10 mm box, seeded.
TriMesh randomSoup(int triangles, uint32_t seed) {
    std::mt19937 random(seed);
    std::uniform_real_distribution<double> coordinate(0.0, 10.0), offset(-1.0, 1.0);
    TriMesh soup;
    for (int t = 0; t < triangles; ++t) {
        const Vec3 centre{coordinate(random), coordinate(random), coordinate(random)};
        for (int k = 0; k < 3; ++k) soup.pos.push_back(centre + Vec3{offset(random), offset(random), offset(random)});
        soup.tri.push_back({3 * t, 3 * t + 1, 3 * t + 2});
    }
    return soup;
}

// The closest point by trying every triangle, equal distances to the lower index.
TriangleBvh::Nearest closestByBruteForce(const TriMesh& soup, const Vec3& q) {
    TriangleBvh::Nearest best;
    for (int t = 0; t < static_cast<int>(soup.tri.size()); ++t) {
        const Vec3& a = soup.pos[soup.tri[t][0]];
        const Vec3& b = soup.pos[soup.tri[t][1]];
        const Vec3& c = soup.pos[soup.tri[t][2]];
        const TrianglePoint point = closestOnTriangle(q, a, b, c);
        const Vec3 p = a + (b - a) * point.bary[1] + (c - a) * point.bary[2];
        const double distSq = dot(q - p, q - p);
        if (best.tri < 0 || distSq < best.distSq) best = {t, distSq, point.bary[1], point.bary[2]};
    }
    return best;
}

// --- the whole wrap

// The wrapped surface as written: Z-up triangles (each quad as quadSplit's two) with their corners' uv, read back through the
// application's reader.
struct WrittenSurface {
    std::vector<Vec3> pos;
    std::vector<std::array<int, 3>> tri;
    std::vector<std::array<Vec2, 3>> uv;
};

WrittenSurface readWrittenSurface(const std::string& path) {
    ImportedMesh file;
    std::string err;
    if (!readMeshFile(path, ImportLimits{}, file, &err)) throw std::runtime_error("cannot read " + path + ": " + err);
    CHECK(file.quadFirstTriangle.empty() && !file.tri.empty() && file.uv.size() == file.tri.size());
    WrittenSurface s;
    for (const Vec3& p : file.pos) s.pos.push_back({p.x, -p.z, p.y});  // Y up back to Z up, exactly
    s.tri = file.tri;
    s.uv = file.uv;
    return s;
}

// The islands of a written surface: triangles joined across an edge whose two ends carry the same uv in both.
std::vector<int> islandsOfWritten(const WrittenSurface& s) {
    std::vector<int> parent(s.tri.size());
    std::iota(parent.begin(), parent.end(), 0);
    std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> sidesOfEdge;  // (lo, hi) -> (triangle, corner)
    for (int t = 0; t < static_cast<int>(s.tri.size()); ++t)
        for (int k = 0; k < 3; ++k) {
            const int a = s.tri[t][k], b = s.tri[t][(k + 1) % 3];
            sidesOfEdge[{std::min(a, b), std::max(a, b)}].push_back({t, k});
        }
    for (const std::pair<const std::pair<int, int>, std::vector<std::pair<int, int>>>& entry : sidesOfEdge) {
        if (entry.second.size() != 2) continue;
        const int f = entry.second[0].first, k = entry.second[0].second, g = entry.second[1].first, j = entry.second[1].second;
        // f runs a -> b at corners k, k + 1; g runs b -> a at corners j, j + 1
        const Vec2 fa = s.uv[f][k], fb = s.uv[f][(k + 1) % 3], gb = s.uv[g][j], ga = s.uv[g][(j + 1) % 3];
        if (!(fa.x == ga.x && fa.y == ga.y && fb.x == gb.x && fb.y == gb.y)) continue;
        int x = f, y = g;
        while (parent[x] != x) x = parent[x];
        while (parent[y] != y) y = parent[y];
        parent[std::max(x, y)] = std::min(x, y);
    }
    std::vector<int> island(s.tri.size());
    for (size_t t = 0; t < island.size(); ++t) {
        int x = static_cast<int>(t);
        while (parent[x] != x) x = parent[x];
        island[t] = x;
    }
    return island;
}

ImportOptions fixtureOptions() {
    ImportOptions o;
    o.unitMm = 1;
    o.up = UpAxis::Z;
    o.surface = ImportOptions::Surface::Wrapped;
    return o;
}

void wrapFixtureInto(const std::string& fixture, const std::string& folder) {
    ImportedMesh file;
    std::string err;
    CHECK(readMeshFile(std::string(KNIT_TEST_DATA) + "/wrap_fixtures/" + fixture + "/fixture.obj", ImportLimits{}, file, &err));
    std::vector<std::string> warnings;
    const std::vector<ColourImage> textures = loadImportTextures(file, ImportLimits{}, &warnings);
    std::vector<std::string> changes;
    if (!wrapFileToFolder(file, textures, fixtureOptions(), folder, changes, &err, {})) throw std::runtime_error(fixture + ": " + err);
}

constexpr const char* kChildFolderVariable = "KNIT_WRAP_PROPERTIES_CHILD_OUT";
const char* const kWrittenFiles[3] = {"wrapped.obj", "wrapped.mtl", "wrapped_basecolor.png"};

// Wraps a fixture at 16 threads and checks its charts. Then the same case runs again in a child process at 1 thread (the
// executor's budget is fixed once it starts), and the three written files must be equal byte for byte.
void checkWholeWrap(const std::string& testName, const std::string& fixture) {
    char* childFolder = nullptr;
    size_t length = 0;
    if (_dupenv_s(&childFolder, &length, kChildFolderVariable) == 0 && childFolder) {
        const std::string folder = childFolder;
        std::free(childFolder);
        setThreadBudget(1);
        wrapFixtureInto(fixture, folder);
        return;
    }
    setThreadBudget(16);
    const std::string many = freshDir(fixture + "_16"), one = freshDir(fixture + "_1");
    wrapFixtureInto(fixture, many);

    const WrittenSurface s = readWrittenSurface(many + "/wrapped.obj");
    for (const std::array<Vec2, 3>& flat : s.uv) {  // the triangles as the reader gives them
        for (const Vec2& p : flat) CHECK(p.x >= 0 && p.x <= 1 && p.y >= 0 && p.y <= 1);
        CHECK(signedArea(flat[0], flat[1], flat[2]) > 0);  // no flipped uv triangle
    }
    CHECK(!islandsOverlap(s.uv, islandsOfWritten(s)));

    char* program = nullptr;
    CHECK(_get_pgmptr(&program) == 0 && program);
    CHECK(_putenv_s(kChildFolderVariable, one.c_str()) == 0);
    const std::string command = std::string("\"\"") + program + "\" " + testName + " > \"" + one + "/child.log\" 2>&1\"";
    const int status = std::system(command.c_str());
    _putenv_s(kChildFolderVariable, "");
    CHECK(status == 0);
    for (const char* name : kWrittenFiles) {
        const std::string a = readBytes(many + "/" + name), b = readBytes(one + "/" + name);
        CHECK(!a.empty() && a == b);
    }
}

} // namespace

// --- quadSplit

TEST(wrap_quad_split_planar) {
    CHECK(quadSplit(flatQuad({0, 0}, {1, 0}, {1, 1}, {0, 1})) == kAlongFirst);
    // an elongated rhombus: the short diagonal, whichever corners it joins
    CHECK(quadSplit(flatQuad({-2, 0}, {0, -1}, {2, 0}, {0, 1})) == kAlongSecond);
    CHECK(quadSplit(flatQuad({0, -1}, {2, 0}, {0, 1}, {-2, 0})) == kAlongFirst);
}

TEST(wrap_quad_split_folded) {
    // A dart: corner 3 lies inside the triangle of the others. The diagonal 0-2 is shorter, but its triangle (0, 2, 3) faces
    // against the quad's normal, so the split along 1-3 is taken.
    const std::array<Vec3, 4> dart = flatQuad({-1, 0}, {0, -3}, {1, 0}, {0, -0.5});
    CHECK(norm(dart[2] - dart[0]) < norm(dart[3] - dart[1]));
    CHECK(quadSplit(dart) == kAlongSecond);
    // folded out of its plane along 0-2: the valley split faces against the normal
    std::array<Vec3, 4> folded = dart;
    folded[3].z = 0.4;
    CHECK(quadSplit(folded) == kAlongSecond);
    // nothing valid: four corners on one point give the first split
    CHECK(quadSplit(flatQuad({1, 1}, {1, 1}, {1, 1}, {1, 1})) == kAlongFirst);
}

TEST(wrap_quad_split_uv_follow_corners) {
    // One quad written with the long diagonal 0-2 first: the source splits it along 1-3, and each corner keeps its own uv.
    ImportedMesh file;
    file.pos = {{-2, 0, 0}, {0, -1, 0}, {2, 0, 0}, {0, 1, 0}};
    file.tri = {{0, 1, 2}, {0, 2, 3}};
    file.quadFirstTriangle = {0};
    const Vec2 uvOf[4] = {{0.1, 0.2}, {0.3, 0.4}, {0.5, 0.6}, {0.7, 0.8}};
    file.uv = {{uvOf[0], uvOf[1], uvOf[2]}, {uvOf[0], uvOf[2], uvOf[3]}};
    file.triMaterial = {0, 0};
    const WrapSource source = wrapSourceOf(file, 1.0, UpAxis::Z);
    CHECK(source.soup.size() == 2 && source.uv.size() == 2);
    for (size_t t = 0; t < 2; ++t)
        for (int k = 0; k < 3; ++k) {
            const int v = source.vertex[t][k];
            CHECK(source.uv[t][k].x == uvOf[v].x && source.uv[t][k].y == uvOf[v].y);
            CHECK(source.soup[t][k] == file.pos[v]);
        }
}

// --- the file's triangles

TEST(wrap_source_soup_y_up_turn) {
    ImportedMesh file;
    file.pos = {{1, 2, 3}, {-4.25, 5, 0.125}, {7, -8, 9}};
    file.tri = {{0, 1, 2}};
    const WrapSource source = wrapSourceOf(file, 2.5, UpAxis::Y);
    CHECK(source.soup.size() == 1 && source.uv.empty() && source.material == std::vector<int>{-1});
    for (int k = 0; k < 3; ++k) {
        const Vec3 mm = file.pos[k] * 2.5;
        CHECK(source.soup[0][k] == (Vec3{mm.x, -mm.z, mm.y}));
    }
    const WrapSource zUp = wrapSourceOf(file, 2.5, UpAxis::Z);
    CHECK(zUp.soup[0][1] == file.pos[1] * 2.5);
}

TEST(wrap_source_soup_duplicate_dropped) {
    ImportedMesh file;
    file.pos = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    file.tri = {{0, 1, 2}, {2, 0, 1}, {0, 1, 3}, {1, 0, 2}};  // 1 and 3 repeat triangle 0's vertices
    file.triMaterial = {4, 5, 6, 7};
    const WrapSource source = wrapSourceOf(file, 1.0, UpAxis::Z);
    CHECK(source.fileTriangle == (std::vector<int>{0, 2}));
    CHECK(source.material == (std::vector<int>{4, 6}));
    CHECK(source.vertex[1] == (std::array<int, 3>{0, 1, 3}));
}

TEST(wrap_source_soup_quads_resplit) {
    // two quads: one written along its long diagonal (re-split), one along its short one (kept)
    ImportedMesh file;
    file.pos = {{-2, 0, 0}, {0, -1, 0}, {2, 0, 0}, {0, 1, 0}, {10, -1, 0}, {12, 0, 0}, {10, 1, 0}, {8, 0, 0}};
    file.tri = {{0, 1, 2}, {0, 2, 3}, {4, 5, 6}, {4, 6, 7}};
    file.quadFirstTriangle = {0, 2};
    file.triMaterial = {1, 1, 2, 2};
    const WrapSource source = wrapSourceOf(file, 1.0, UpAxis::Z);
    CHECK(source.soup.size() == 4);
    CHECK(source.vertex[0] == (std::array<int, 3>{0, 1, 3}) && source.vertex[1] == (std::array<int, 3>{1, 2, 3}));
    CHECK(source.vertex[2] == (std::array<int, 3>{4, 5, 6}) && source.vertex[3] == (std::array<int, 3>{4, 6, 7}));
    CHECK(source.fileTriangle == (std::vector<int>{0, 1, 2, 3}));
    CHECK(source.material == (std::vector<int>{1, 1, 2, 2}));
}

// --- the remesh and its winding

TEST(wrap_remesh_cube) {
    const QuadMesh cube = cubeMesh(20.0);
    QuadMesh remeshed;
    std::string err;
    CHECK(remeshToEvenQuads(cube, 1.5, remeshed, &err));
    CHECK(remeshed.quad.size() > 100 && remeshed.edge.empty());
    CHECK(partCount(remeshed) == 1);
    CHECK(everyEdgeTwiceOpposite(remeshed));
    CHECK_NEAR(sixTimesVolume(remeshed) / 6.0, 8000.0, 0.15 * 8000.0);
}

TEST(wrap_remesh_inside_out) {
    const QuadMesh cube = cubeMesh(20.0);
    QuadMesh insideOut = cube;
    for (std::array<int, 4>& q : insideOut.quad) std::swap(q[1], q[3]);
    QuadMesh outward, fromInsideOut;
    std::string err;
    CHECK(remeshToEvenQuads(cube, 1.5, outward, &err));
    CHECK(remeshToEvenQuads(insideOut, 1.5, fromInsideOut, &err));
    CHECK(sixTimesVolume(fromInsideOut) > 0 && everyEdgeTwiceOpposite(fromInsideOut));
    CHECK(fromInsideOut.quad == outward.quad && fromInsideOut.pos == outward.pos);
}

// --- the chart rules

// One round of the chart rules on a remeshed cube grown into two charts, chart 0 judged `stretchedShare` outside the bounds
// (its triangles twice past them), chart 1 inside.
struct RuleRound {
    ChartRound round;
    ChartStep step;
    size_t chartsBefore = 0, chartsAfter = 0;
};
RuleRound oneRuleRound(double stretchedShare) {
    QuadMesh surface;
    std::string err;
    CHECK(remeshToEvenQuads(cubeMesh(20.0), 1.5, surface, &err));
    surface.edge = sortedEdges(surface.quad);
    const FaceGraph graph = buildFaceGraph(surface);
    RuleRound r;
    r.round.label = growCharts(graph, spreadSeeds(graph, 2));
    r.round.seam = cutToDiscs(graph, r.round.label);
    r.round.slit.assign(r.round.seam.size(), 0);
    ChartQuality quality;
    quality.chart = {0, 1};
    quality.stretched = {stretchedShare, 0.0};
    quality.unusable = {0.0, 0.0};
    quality.totalOutside = 0.5 * stretchedShare;
    for (size_t f = 0; f < surface.quad.size(); ++f)
        for (int k = 0; k < 2; ++k) {
            quality.faceOfTriangle.push_back(static_cast<int>(f));
            quality.worst.push_back(r.round.label[f] == 0 ? 2.0 : 0.5);
        }
    r.chartsBefore = std::set<int>(r.round.label.begin(), r.round.label.end()).size();
    r.step = decideNextRound(graph, quality, ChartSettings{}, r.round);
    r.chartsAfter = std::set<int>(r.round.label.begin(), r.round.label.end()).size();
    return r;
}

// A slit relieves a spike; a chart mostly outside the bounds (the tank's gun barrel, a tube closed at one end) does not lie flat
// with any slit, so it is split, as a folded chart is. A chart only partly outside is still slit.
TEST(wrap_chart_rules_split_a_mostly_stretched_chart) {
    const RuleRound mostly = oneRuleRound(0.9);
    CHECK(mostly.chartsBefore == 2);
    CHECK(mostly.step.decision == ChartDecision::FlattenChanged);
    CHECK(mostly.chartsAfter == 3);
    CHECK(mostly.round.slitCount == 0);

    const RuleRound partly = oneRuleRound(0.3);
    CHECK(partly.step.decision == ChartDecision::FlattenChanged);
    CHECK(partly.chartsAfter == 2);
    CHECK(partly.round.slitCount == 1);
}

// --- the edge list

TEST(wrap_sorted_edges) {
    const QuadMesh cube = cubeMesh(1.0);
    const std::vector<std::array<int, 2>> edges = sortedEdges(cube.quad);
    CHECK(edges.size() == 12);
    for (size_t e = 0; e < edges.size(); ++e) {
        CHECK(edges[e][0] < edges[e][1]);
        if (e > 0) CHECK(edges[e - 1] < edges[e]);
    }
    // the same list for the quads in another order, each starting at another corner
    std::vector<std::array<int, 4>> shuffled = cube.quad;
    std::mt19937 random(7);
    std::shuffle(shuffled.begin(), shuffled.end(), random);
    for (size_t f = 0; f < shuffled.size(); ++f) std::rotate(shuffled[f].begin(), shuffled[f].begin() + f % 4, shuffled[f].end());
    CHECK(sortedEdges(shuffled) == edges);
}

// --- the packer

TEST(wrap_pack_rectangles_stand_up) {
    const RectangleIslands r = fiveRectangles();
    QuadUv uv;
    packAll(r, uv);
    const double degree = kPi / 180.0;
    for (size_t f = 0; f < uv.size(); ++f) {
        const Vec2 first = uv[f][1] - uv[f][0], second = uv[f][2] - uv[f][1];
        const Vec2 longer = norm(first) >= norm(second) ? first : second;
        const Vec2 shorter = norm(first) >= norm(second) ? second : first;
        CHECK(std::acos(std::min(1.0, std::fabs(longer.y) / norm(longer))) < 0.1 * degree);
        const Vec2& size = kFiveSizes[f];
        CHECK_NEAR(norm(longer) / norm(shorter), std::max(size.x, size.y) / std::min(size.x, size.y), 1e-9);
    }
}

TEST(wrap_pack_islands_do_not_overlap) {
    // the six faces of a cube, each its own island (every edge a seam), with uneven random quads as their uv
    const QuadMesh cube = cubeMesh(1.0);
    const FaceGraph graph = buildFaceGraph(cube);
    const std::vector<char> seam(cube.edge.size(), 1);
    const std::vector<char> faces(cube.quad.size(), 1);
    QuadUv uv(cube.quad.size());
    std::mt19937 random(11);
    std::uniform_real_distribution<double> jitter(-0.2, 0.2), size(0.5, 4.0), place(-5.0, 5.0);
    for (std::array<Vec2, 4>& corners : uv) {
        const double w = size(random), h = size(random), cu = place(random), cv = place(random);
        corners = {Vec2{cu + jitter(random), cv + jitter(random)}, Vec2{cu + w + jitter(random), cv + jitter(random)},
                   Vec2{cu + w + jitter(random), cv + h + jitter(random)}, Vec2{cu + jitter(random), cv + h + jitter(random)}};
    }
    packCharts(graph, seam, faces, uv);
    std::vector<std::array<Vec2, 3>> triangles;
    std::vector<int> islandOfTriangle;
    for (size_t f = 0; f < uv.size(); ++f)
        for (const std::array<int, 3>& t : kAlongFirst) {
            triangles.push_back({uv[f][t[0]], uv[f][t[1]], uv[f][t[2]]});
            islandOfTriangle.push_back(static_cast<int>(f));
        }
    for (size_t i = 0; i < triangles.size(); ++i)
        for (size_t j = i + 1; j < triangles.size(); ++j)
            if (islandOfTriangle[i] != islandOfTriangle[j]) CHECK(!trianglesOverlap(triangles[i], triangles[j]));
    CHECK(!islandsOverlap(triangles, islandOfTriangle));
}

TEST(wrap_pack_uv_in_unit_square) {
    QuadUv uv;
    packAll(fiveRectangles(), uv);
    double highest = 0;
    for (const std::array<Vec2, 4>& corners : uv)
        for (const Vec2& p : corners) {
            CHECK(p.x >= 0 && p.x <= 1 && p.y >= 0 && p.y <= 1);
            highest = std::max({highest, p.x, p.y});
        }
    CHECK(highest > 0.9);  // the layout fills the square along its larger extent
}

TEST(wrap_pack_gutter_holds) {
    const RectangleIslands r = fiveRectangles();
    QuadUv uv;
    packAll(r, uv);
    struct Box {
        double u0 = 1e300, v0 = 1e300, u1 = -1e300, v1 = -1e300;
    };
    std::vector<Box> boxes(uv.size());
    double areaSum = 0;
    for (size_t f = 0; f < uv.size(); ++f) {
        for (const Vec2& p : uv[f]) {
            boxes[f].u0 = std::min(boxes[f].u0, p.x), boxes[f].u1 = std::max(boxes[f].u1, p.x);
            boxes[f].v0 = std::min(boxes[f].v0, p.y), boxes[f].v1 = std::max(boxes[f].v1, p.y);
        }
        areaSum += norm(uv[f][1] - uv[f][0]) * norm(uv[f][2] - uv[f][1]);
    }
    const double gutter = 0.002 * std::sqrt(areaSum);  // kGutter, at the final scale
    const double slack = 1e-12;
    for (size_t a = 0; a < boxes.size(); ++a) {
        CHECK(boxes[a].u0 >= gutter - slack && boxes[a].v0 >= gutter - slack);
        CHECK(boxes[a].u1 <= 1 - gutter + slack && boxes[a].v1 <= 1 - gutter + slack);
        for (size_t b = a + 1; b < boxes.size(); ++b) {
            const double gapU = std::max(boxes[b].u0 - boxes[a].u1, boxes[a].u0 - boxes[b].u1);
            const double gapV = std::max(boxes[b].v0 - boxes[a].v1, boxes[a].v0 - boxes[b].v1);
            CHECK(std::max(gapU, gapV) >= 2.0 * gutter - slack);
        }
    }
}

TEST(wrap_pack_unflagged_untouched) {
    const RectangleIslands r = fiveRectangles();
    const FaceGraph graph = unlinkedFaces(r.mesh);
    const std::vector<char> seam;
    const std::vector<char> faces = {1, 0, 1, 0, 1};
    QuadUv uv = r.uv;
    packCharts(graph, seam, faces, uv);
    for (size_t f = 0; f < uv.size(); ++f)
        for (int k = 0; k < 4; ++k) {
            const bool same = sameBits(uv[f][k].x, r.uv[f][k].x) && sameBits(uv[f][k].y, r.uv[f][k].y);
            CHECK(same == !faces[f]);
        }
}

// --- the closest point

TEST(wrap_closest_point_regions) {
    TriMesh soup;
    soup.pos = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
    soup.tri = {{0, 1, 2}};
    const TriangleBvh bvh(soup);
    struct Expected {
        Vec3 q;
        double distSq, b1, b2;
    };
    const Expected cases[] = {
        {{-1, -1, 1}, 3.0, 0.0, 0.0},     // the region of corner 0
        {{2, -1, 0}, 2.0, 1.0, 0.0},      // corner 1
        {{0.5, -1, 0}, 1.0, 0.5, 0.0},    // edge 0-1
        {{1, 1, 0}, 0.5, 0.5, 0.5},       // edge 1-2
        {{0.25, 0.25, 2}, 4.0, 0.25, 0.25},  // the face
    };
    for (const Expected& e : cases) {
        const TriangleBvh::Nearest n = bvh.closestPoint(e.q);
        CHECK(n.tri == 0);
        CHECK_NEAR(n.distSq, e.distSq, 1e-15);
        CHECK_NEAR(n.b1, e.b1, 1e-15);
        CHECK_NEAR(n.b2, e.b2, 1e-15);
    }
    // a triangle of zero area: the closest point of its edges
    TriMesh line;
    line.pos = {{0, 0, 0}, {1, 0, 0}, {2, 0, 0}};
    line.tri = {{0, 1, 2}};
    const TriangleBvh lineBvh(line);
    const TriangleBvh::Nearest n = lineBvh.closestPoint({1.5, 1, 0});
    CHECK(n.tri == 0);
    CHECK_NEAR(n.distSq, 1.0, 1e-15);
    CHECK_NEAR(n.b1 * 1.0 + n.b2 * 2.0, 1.5, 1e-15);
}

TEST(wrap_closest_point_tie_lower_index) {
    // far triangles around two equidistant ones, at 7 (above the query) and 13 (below): the lower index wins either way round
    for (int round = 0; round < 2; ++round) {
        TriMesh soup = randomSoup(20, 3);
        for (Vec3& p : soup.pos) p = p + Vec3{100, 100, 100};
        const double above = round == 0 ? 1.0 : -1.0;
        soup.pos[soup.tri[7][0]] = {0, 0, above}, soup.pos[soup.tri[7][1]] = {1, 0, above}, soup.pos[soup.tri[7][2]] = {0, 1, above};
        soup.pos[soup.tri[13][0]] = {0, 0, -above}, soup.pos[soup.tri[13][1]] = {1, 0, -above}, soup.pos[soup.tri[13][2]] = {0, 1, -above};
        const TriangleBvh bvh(soup);
        const TriangleBvh::Nearest n = bvh.closestPoint({0.2, 0.2, 0});
        CHECK(n.tri == 7);
        CHECK(n.distSq == 1.0);
    }
}

TEST(wrap_closest_point_brute_force) {
    const TriMesh soup = randomSoup(400, 5);
    const TriangleBvh bvh(soup);
    std::mt19937 random(9);
    std::uniform_real_distribution<double> coordinate(-3.0, 13.0);
    for (int i = 0; i < 2000; ++i) {
        const Vec3 q{coordinate(random), coordinate(random), coordinate(random)};
        const TriangleBvh::Nearest fast = bvh.closestPoint(q);
        const TriangleBvh::Nearest slow = closestByBruteForce(soup, q);
        CHECK(fast.tri == slow.tri);
        CHECK_NEAR(fast.distSq, slow.distSq, 1e-12);
    }
}

TEST(wrap_closest_point_empty_and_not_finite) {
    const TriMesh empty;
    CHECK(TriangleBvh(empty).closestPoint({0, 0, 0}).tri == -1);
    const TriMesh soup = randomSoup(10, 1);
    const TriangleBvh bvh(soup);
    CHECK(bvh.closestPoint({std::numeric_limits<double>::quiet_NaN(), 0, 0}).tri == -1);
    CHECK(bvh.closestPoint({0, std::numeric_limits<double>::infinity(), 0}).tri == -1);
    CHECK(bvh.closestPoint({1, 2, 3}).tri >= 0);
}

// --- vertex normals

TEST(wrap_vertex_normals_cube) {
    const QuadMesh cube = cubeMesh(2.0);
    const FaceGraph graph = buildFaceGraph(cube);
    const std::vector<Vec3> normals = vertexNormalsOf(graph);
    CHECK(normals.size() == 8);
    const double third = 1.0 / std::sqrt(3.0);
    for (size_t v = 0; v < 8; ++v) {
        const Vec3 expected{cube.pos[v].x > 0 ? third : -third, cube.pos[v].y > 0 ? third : -third, cube.pos[v].z > 0 ? third : -third};
        CHECK_NEAR(norm(normals[v] - expected), 0.0, 1e-12);
    }
}

// --- the surface file

// The file holds every position bit for bit (Y up), one uv per quad corner bit for bit, and each quad as quadSplit's two
// triangles in quad order, with the uv of their corners.
TEST(wrap_surface_file_roundtrip) {
    QuadMesh mesh = cubeMesh(1.0);
    for (size_t v = 0; v < mesh.pos.size(); ++v) mesh.pos[v] = mesh.pos[v] * (1.0 / 3.0) + Vec3{0.1 * v, -1e-7 * v, 12345.678901234567};
    QuadUv uv(mesh.quad.size());
    for (size_t f = 0; f < uv.size(); ++f)
        for (int k = 0; k < 4; ++k) uv[f][k] = {(f * 4 + k) / 23.0, 1.0 - std::sqrt(0.5 + 0.01 * k)};
    const std::string path = freshDir("surface_file") + "/wrapped.obj";
    CHECK(writeSurfaceObj(path, mesh, uv));
    CHECK(readBytes(path).rfind("# wrapped by the built-in import\nmtllib wrapped.mtl\no wrapped\n", 0) == 0);

    ImportedMesh file;
    std::string err;
    CHECK(readMeshFile(path, ImportLimits{}, file, &err));
    CHECK(file.pos.size() == mesh.pos.size());
    for (size_t v = 0; v < mesh.pos.size(); ++v) {  // Y up, every bit kept
        CHECK(sameBits(file.pos[v].x, mesh.pos[v].x) && sameBits(file.pos[v].y, mesh.pos[v].z) && sameBits(file.pos[v].z, -mesh.pos[v].y));
    }
    CHECK(file.quadFirstTriangle.empty() && file.tri.size() == 2 * mesh.quad.size());
    for (size_t f = 0; f < mesh.quad.size(); ++f) {
        const std::array<int, 4>& q = mesh.quad[f];
        const std::array<std::array<int, 3>, 2> split = quadSplit({mesh.pos[q[0]], mesh.pos[q[1]], mesh.pos[q[2]], mesh.pos[q[3]]});
        for (int i = 0; i < 2; ++i) {
            const size_t t = 2 * f + i;
            for (int k = 0; k < 3; ++k) {
                CHECK(file.tri[t][k] == q[split[i][k]]);
                CHECK(sameBits(file.uv[t][k].x, uv[f][split[i][k]].x) && sameBits(file.uv[t][k].y, uv[f][split[i][k]].y));
            }
        }
    }
}

namespace {

// A cube of side 20 mm with its corner 6 pulled out along the diagonal: the three faces at that corner become kites whose
// shorter diagonal avoids corner 6, so quadSplit cuts them along 1-3 (corner 6 is corner 2 of each); the other three stay
// squares and are cut along 0-2. Every face is its own uv chart, laid out apart from the others. Face 1's uv is a dart whose
// reflex corner is corner 1: cut along 1-3 both triangles keep their orientation, cut along 0-2 the triangle (0, 1, 2) folds.
struct KiteCube {
    QuadMesh mesh;
    QuadUv uv;
};

KiteCube kiteCube() {
    KiteCube c;
    c.mesh = cubeMesh(20.0);
    c.mesh.pos[6] = c.mesh.pos[6] * 1.3;
    c.uv.resize(c.mesh.quad.size());
    for (size_t f = 0; f < c.mesh.quad.size(); ++f) {
        const Vec2 offset = {30.0 * f, 0.0};
        c.uv[f] = {offset + Vec2{0, 0}, offset + Vec2{20, 0}, offset + Vec2{20, 20}, offset + Vec2{0, 20}};
    }
    c.uv[1] = {Vec2{30, 0}, Vec2{40, 8}, Vec2{50, 0}, Vec2{40, 20}};
    return c;
}

// The index of the vertex of `pos` nearest to `p`.
int nearestVertex(const std::vector<Vec3>& pos, const Vec3& p) {
    int best = 0;
    for (size_t v = 1; v < pos.size(); ++v)
        if (norm(pos[v] - p) < norm(pos[best] - p)) best = static_cast<int>(v);
    return best;
}

} // namespace

// What the knit consumes: the surface written, then read back as the application reads a rebuilt surface (the OBJ reader and
// the repairs of prepareImport). Every triangle it gets is one of quadSplit's, with the same vertices, winding and uv, each
// once; and the knit's own chart check (UvAtlas) finds no unresolved chart, the dart face included.
TEST(wrap_surface_file_keeps_the_flattened_triangles) {
    const KiteCube cube = kiteCube();
    std::map<std::array<int, 3>, std::array<Vec2, 3>> expected;  // vertex triple, rotated to start at its lowest vertex -> uv
    int splitsAlong13 = 0;
    for (size_t f = 0; f < cube.mesh.quad.size(); ++f) {
        const std::array<int, 4>& q = cube.mesh.quad[f];
        const std::array<std::array<int, 3>, 2> split = quadSplit({cube.mesh.pos[q[0]], cube.mesh.pos[q[1]], cube.mesh.pos[q[2]], cube.mesh.pos[q[3]]});
        splitsAlong13 += split[0] == std::array<int, 3>{0, 1, 3} ? 1 : 0;
        for (const std::array<int, 3>& t : split) {
            int first = 0;
            for (int k = 1; k < 3; ++k)
                if (q[t[k]] < q[t[first]]) first = k;
            std::array<int, 3> vertices;
            std::array<Vec2, 3> corners;
            for (int k = 0; k < 3; ++k) {
                vertices[k] = q[t[(first + k) % 3]];
                corners[k] = cube.uv[f][t[(first + k) % 3]];
            }
            expected[vertices] = corners;
        }
    }
    CHECK(splitsAlong13 == 3);  // the three kites; the dart is one of them

    const std::string folder = freshDir("surface_file_triangles");
    CHECK(writeSurfaceObj(folder + "/wrapped.obj", cube.mesh, cube.uv));
    PreparedMesh prepared;
    std::string err;
    CHECK(readRebuiltSurface(folder, ImportOptions{}, {}, prepared, &err, {}));
    const TriMesh& read = prepared.mesh;
    CHECK(read.tri.size() == expected.size() && read.uv.size() == read.tri.size());

    std::vector<Vec3> library(cube.mesh.pos.size());  // the cube where the repairs put it: Z up, plus their offset
    for (size_t v = 0; v < library.size(); ++v) library[v] = cube.mesh.pos[v] + prepared.offset;
    std::map<std::array<int, 3>, int> seen;
    for (size_t t = 0; t < read.tri.size(); ++t) {
        std::array<int, 3> vertices;
        for (int k = 0; k < 3; ++k) {
            const Vec3& p = read.pos[read.tri[t][k]];
            vertices[k] = nearestVertex(library, p);
            CHECK(norm(library[vertices[k]] - p) < 1e-9);
        }
        int first = 0;
        for (int k = 1; k < 3; ++k)
            if (vertices[k] < vertices[first]) first = k;
        const std::array<int, 3> rotated = {vertices[first], vertices[(first + 1) % 3], vertices[(first + 2) % 3]};
        const std::map<std::array<int, 3>, std::array<Vec2, 3>>::const_iterator it = expected.find(rotated);
        CHECK(it != expected.end());  // same vertices in the same winding
        if (it == expected.end()) continue;
        ++seen[rotated];
        for (int k = 0; k < 3; ++k) {
            const Vec2& got = read.uv[t][(first + k) % 3];
            CHECK(sameBits(got.x, it->second[k].x) && sameBits(got.y, it->second[k].y));
        }
    }
    CHECK(seen.size() == expected.size());
    for (const std::pair<const std::array<int, 3>, int>& entry : seen) CHECK(entry.second == 1);

    const UvAtlas atlas(read);
    CHECK(atlas.charts().size() == cube.mesh.quad.size());
    for (const UvChart& chart : atlas.charts()) {
        CHECK(chart.status != Status::Unresolved);
        CHECK(chart.foldedArea == 0);
    }
}

// --- the whole wrap

TEST(wrap_whole_wrap_box) { checkWholeWrap("wrap_whole_wrap_box", "box"); }
TEST(wrap_whole_wrap_boxes_apart) { checkWholeWrap("wrap_whole_wrap_boxes_apart", "boxes_apart"); }
TEST(wrap_whole_wrap_touching_boxes) { checkWholeWrap("wrap_whole_wrap_touching_boxes", "touching_boxes"); }
TEST(wrap_whole_wrap_torus) { checkWholeWrap("wrap_whole_wrap_torus", "torus"); }
