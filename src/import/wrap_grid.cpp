// The voxel stages of the built-in import: the file's triangles in millimetres, then a sealed solid and its blurred field. Voxel
// coordinates round half to even, and the grid's border wraps around in the dilation and the blur.
#include "import/wrap_stages.h"

#include <algorithm>
#include <cmath>

#include "model/cancel.h"
#include "model/coordinates.h"
#include "model/parallel.h"

namespace knit {

namespace {

constexpr double kSampleSpacingVoxels = 0.5;

// The smallest non-negative value congruent to `value` modulo `size`: an index that wraps around the grid.
int wrapped(int value, int size) {
    const int rest = value % size;
    return rest < 0 ? rest + size : rest;
}

double longestEdge(const TriangleCorners& t) {
    return std::max({norm(t[0] - t[1]), norm(t[1] - t[2]), norm(t[2] - t[0])});
}

// Rounds halves to even: the default floating-point rounding mode.
int nearestVoxel(double coordinate) { return static_cast<int>(std::nearbyint(coordinate)); }

// All integer offsets inside a ball, as (x, y, z).
std::vector<std::array<int, 3>> ballOffsets(int radius) {
    std::vector<std::array<int, 3>> offsets;
    for (int x = -radius; x <= radius; ++x)
        for (int y = -radius; y <= radius; ++y)
            for (int z = -radius; z <= radius; ++z)
                if (x * x + y * y + z * z <= radius * radius) offsets.push_back({x, y, z});
    return offsets;
}

// One blur pass along one axis. `step` is the stride of that axis in the array, `length` its size.
void blurAlongAxis(const std::array<int, 3>& size, int axis, const std::vector<float>& in, std::vector<float>& out) {
    parallelFor(size[0], 1,
                [&](int64_t begin, int64_t end, int64_t) {
                    for (int x = static_cast<int>(begin); x < end; ++x)
                        for (int y = 0; y < size[1]; ++y)
                            for (int z = 0; z < size[2]; ++z) {
                                std::array<int, 3> below = {x, y, z};
                                std::array<int, 3> above = {x, y, z};
                                below[axis] = wrapped(below[axis] - 1, size[axis]);
                                above[axis] = wrapped(above[axis] + 1, size[axis]);
                                const float a = in[voxelIndex(size, below[0], below[1], below[2])];
                                const float b = in[voxelIndex(size, x, y, z)];
                                const float c = in[voxelIndex(size, above[0], above[1], above[2])];
                                out[voxelIndex(size, x, y, z)] = 0.25f * a + 0.5f * b + 0.25f * c;
                            }
                },
                ExecPolicy{.parallel = true});
}

// A file position in millimetres, Z up. The Y-up turn (x, y, z) -> (x, -z, y) is exact: it only moves and negates numbers.
Vec3 millimetresZUp(const Vec3& p, double unitMm, UpAxis up) {
    const Vec3 mm = p * unitMm;
    return up == UpAxis::Y ? Vec3{mm.x, -mm.z, mm.y} : mm;
}

// One triangle of the source before duplicates are dropped: the file triangle it stands for, and which of the file's corners
// (and their uv) it takes.
struct SourceTriangle {
    int fileTriangle = 0;
    std::array<int, 3> vertex = {0, 0, 0};
    std::array<Vec2, 3> uv;
};

// A triangle and its vertex numbers in ascending order, for finding repeats.
struct SortedTriple {
    std::array<int, 3> vertex = {0, 0, 0};
    int triangle = 0;
};

bool tripleComesFirst(const SortedTriple& a, const SortedTriple& b) {
    if (a.vertex != b.vertex) return a.vertex < b.vertex;
    return a.triangle < b.triangle;
}

// The file's triangles in file order, each quad of the file split again by quadSplit on its own four corners. The quad's two
// triangles take the numbers of the file's two.
std::vector<SourceTriangle> trianglesOfFile(const ImportedMesh& file, const std::vector<Vec3>& positionMm) {
    const bool hasUv = !file.uv.empty();
    std::vector<char> startsQuad(file.tri.size(), 0);
    for (int first : file.quadFirstTriangle) startsQuad[first] = 1;
    std::vector<SourceTriangle> out;
    out.reserve(file.tri.size());
    for (int t = 0; t < static_cast<int>(file.tri.size()); ++t) {
        if (!startsQuad[t]) {
            SourceTriangle triangle;
            triangle.fileTriangle = t;
            triangle.vertex = file.tri[t];
            if (hasUv) triangle.uv = file.uv[t];
            out.push_back(triangle);
            continue;
        }
        // corners t[first][0..2] and t[first + 1][2]
        const std::array<int, 4> vertex = {file.tri[t][0], file.tri[t][1], file.tri[t][2], file.tri[t + 1][2]};
        std::array<Vec2, 4> uv;
        if (hasUv) uv = {file.uv[t][0], file.uv[t][1], file.uv[t][2], file.uv[t + 1][2]};
        const std::array<Vec3, 4> corners = {positionMm[vertex[0]], positionMm[vertex[1]], positionMm[vertex[2]], positionMm[vertex[3]]};
        const std::array<std::array<int, 3>, 2> split = quadSplit(corners);
        for (int half = 0; half < 2; ++half) {
            SourceTriangle triangle;
            triangle.fileTriangle = t + half;
            for (int k = 0; k < 3; ++k) {
                triangle.vertex[k] = vertex[split[half][k]];
                triangle.uv[k] = uv[split[half][k]];
            }
            out.push_back(triangle);
        }
        ++t;  // the quad's second triangle is done
    }
    return out;
}

// Per triangle: whether an earlier triangle has the same three vertices, in any order.
std::vector<char> repeatedTriangles(const std::vector<SourceTriangle>& triangles) {
    std::vector<SortedTriple> triples(triangles.size());
    for (size_t t = 0; t < triangles.size(); ++t) {
        triples[t].vertex = triangles[t].vertex;
        std::sort(triples[t].vertex.begin(), triples[t].vertex.end());
        triples[t].triangle = static_cast<int>(t);
    }
    std::sort(triples.begin(), triples.end(), tripleComesFirst);
    std::vector<char> repeated(triangles.size(), 0);
    for (size_t i = 1; i < triples.size(); ++i)
        if (triples[i].vertex == triples[i - 1].vertex) repeated[triples[i].triangle] = 1;
    return repeated;
}

} // namespace

WrapSource wrapSourceOf(const ImportedMesh& file, double unitMm, UpAxis up) {
    std::vector<Vec3> positionMm(file.pos.size());
    for (size_t v = 0; v < file.pos.size(); ++v) positionMm[v] = millimetresZUp(file.pos[v], unitMm, up);
    const std::vector<SourceTriangle> triangles = trianglesOfFile(file, positionMm);
    const std::vector<char> repeated = repeatedTriangles(triangles);

    WrapSource source;
    for (size_t t = 0; t < triangles.size(); ++t) {
        if (repeated[t]) continue;
        const SourceTriangle& triangle = triangles[t];
        source.soup.push_back({positionMm[triangle.vertex[0]], positionMm[triangle.vertex[1]], positionMm[triangle.vertex[2]]});
        source.fileTriangle.push_back(triangle.fileTriangle);
        source.vertex.push_back(triangle.vertex);
        if (!file.uv.empty()) source.uv.push_back(triangle.uv);
        const bool hasMaterial = triangle.fileTriangle < static_cast<int>(file.triMaterial.size());
        source.material.push_back(hasMaterial ? file.triMaterial[triangle.fileTriangle] : -1);
    }
    return source;
}

GridPlacement placeGrid(const std::vector<TriangleCorners>& soup, double voxelMm, int paddingVoxels) {
    Vec3 low(1e300, 1e300, 1e300), high(-1e300, -1e300, -1e300);
    for (const TriangleCorners& t : soup)
        for (const Vec3& p : t) {
            low = {std::min(low.x, p.x), std::min(low.y, p.y), std::min(low.z, p.z)};
            high = {std::max(high.x, p.x), std::max(high.y, p.y), std::max(high.z, p.z)};
        }
    GridPlacement grid;
    grid.voxelMm = voxelMm;
    const double pad = paddingVoxels * voxelMm;
    grid.origin = {low.x - pad, low.y - pad, low.z - pad};
    grid.size = {static_cast<int>(std::ceil((high.x + pad - grid.origin.x) / voxelMm)) + 1,
                 static_cast<int>(std::ceil((high.y + pad - grid.origin.y) / voxelMm)) + 1,
                 static_cast<int>(std::ceil((high.z + pad - grid.origin.z) / voxelMm)) + 1};
    return grid;
}

SolidGrid markTriangles(const std::vector<TriangleCorners>& soup, const GridPlacement& grid) {
    SolidGrid surface;
    surface.size = grid.size;
    surface.voxels.assign(static_cast<size_t>(grid.voxelCount()), 0);
    const double spacing = kSampleSpacingVoxels * grid.voxelMm;
    for (const TriangleCorners& t : soup) {
        cancelCheckpoint();
        const int steps = std::max(1, static_cast<int>(std::ceil(longestEdge(t) / spacing)));
        for (int i = 0; i <= steps; ++i)
            for (int j = 0; i + j <= steps; ++j) {
                const double a = static_cast<double>(i) / steps;
                const double b = static_cast<double>(j) / steps;
                const double c = 1.0 - a - b;
                const Vec3 p(t[0].x * c + t[1].x * a + t[2].x * b, t[0].y * c + t[1].y * a + t[2].y * b,
                             t[0].z * c + t[1].z * a + t[2].z * b);
                const int x = nearestVoxel((p.x - grid.origin.x) / grid.voxelMm);
                const int y = nearestVoxel((p.y - grid.origin.y) / grid.voxelMm);
                const int z = nearestVoxel((p.z - grid.origin.z) / grid.voxelMm);
                surface.voxels[voxelIndex(grid.size, x, y, z)] = 1;
            }
    }
    return surface;
}

SolidGrid dilateByBall(const SolidGrid& mask, int radiusVoxels) {
    if (radiusVoxels == 0) return mask;
    const std::array<int, 3>& size = mask.size;
    const std::vector<std::array<int, 3>> offsets = ballOffsets(radiusVoxels);

    // Most rows are empty, and skipping them gives the same union much sooner.
    std::vector<uint8_t> rowHasVoxels(static_cast<size_t>(size[0]) * size[1], 0);
    for (int x = 0; x < size[0]; ++x)
        for (int y = 0; y < size[1]; ++y) {
            const uint8_t* row = mask.voxels.data() + voxelIndex(size, x, y, 0);
            rowHasVoxels[static_cast<size_t>(x) * size[1] + y] = std::any_of(row, row + size[2], [](uint8_t v) { return v != 0; });
        }

    SolidGrid out;
    out.size = size;
    out.voxels.assign(mask.voxels.size(), 0);
    // A chunk owns its x planes of the result, so no two threads write the same voxel.
    parallelFor(size[0], 1,
                [&](int64_t begin, int64_t end, int64_t) {
                    for (int x = static_cast<int>(begin); x < end; ++x) {
                        cancelCheckpoint();
                        for (const std::array<int, 3>& offset : offsets) {
                            const int sourceX = wrapped(x - offset[0], size[0]);
                            for (int y = 0; y < size[1]; ++y) {
                                const int sourceY = wrapped(y - offset[1], size[1]);
                                if (!rowHasVoxels[static_cast<size_t>(sourceX) * size[1] + sourceY]) continue;
                                const uint8_t* source = mask.voxels.data() + voxelIndex(size, sourceX, sourceY, 0);
                                uint8_t* target = out.voxels.data() + voxelIndex(size, x, y, 0);
                                const int shift = wrapped(offset[2], size[2]);
                                for (int z = 0; z < size[2] - shift; ++z) target[z + shift] |= source[z];
                                for (int z = size[2] - shift; z < size[2]; ++z) target[z + shift - size[2]] |= source[z];
                            }
                        }
                    }
                },
                ExecPolicy{.parallel = true});
    return out;
}

SolidGrid floodFromGridBorder(const SolidGrid& blocked) {
    const std::array<int, 3>& size = blocked.size;
    SolidGrid reached;
    reached.size = size;
    reached.voxels.assign(blocked.voxels.size(), 0);
    std::vector<std::array<int, 3>> frontier;
    const auto visit = [&](int x, int y, int z) {
        if (x < 0 || y < 0 || z < 0 || x >= size[0] || y >= size[1] || z >= size[2]) return;
        const int64_t i = voxelIndex(size, x, y, z);
        if (blocked.voxels[i] || reached.voxels[i]) return;
        reached.voxels[i] = 1;
        frontier.push_back({x, y, z});
    };
    for (int x = 0; x < size[0]; ++x)
        for (int y = 0; y < size[1]; ++y)
            for (int z = 0; z < size[2]; ++z)
                if (x == 0 || y == 0 || z == 0 || x == size[0] - 1 || y == size[1] - 1 || z == size[2] - 1) visit(x, y, z);
    // A flood fill: any order of filling reaches the same voxels.
    while (!frontier.empty()) {
        const std::array<int, 3> v = frontier.back();
        frontier.pop_back();
        visit(v[0] - 1, v[1], v[2]);
        visit(v[0] + 1, v[1], v[2]);
        visit(v[0], v[1] - 1, v[2]);
        visit(v[0], v[1] + 1, v[2]);
        visit(v[0], v[1], v[2] - 1);
        visit(v[0], v[1], v[2] + 1);
    }
    return reached;
}

SolidGrid invertedGrid(const SolidGrid& mask) {
    SolidGrid out = mask;
    for (uint8_t& v : out.voxels) v = v ? 0 : 1;
    return out;
}

FieldGrid blurredSolidField(const SolidGrid& solid, int passes) {
    FieldGrid field;
    field.size = solid.size;
    field.values.resize(solid.voxels.size());
    for (size_t i = 0; i < solid.voxels.size(); ++i) field.values[i] = solid.voxels[i] ? 1.0f : -1.0f;
    std::vector<float> next(field.values.size());
    for (int pass = 0; pass < passes; ++pass)
        for (int axis = 0; axis < 3; ++axis) {
            cancelCheckpoint();
            blurAlongAxis(field.size, axis, field.values, next);
            field.values.swap(next);
        }
    const std::array<int, 3>& size = field.size;
    for (int x = 0; x < size[0]; ++x)
        for (int y = 0; y < size[1]; ++y)
            for (int z = 0; z < size[2]; ++z)
                if (x == 0 || y == 0 || z == 0 || x == size[0] - 1 || y == size[1] - 1 || z == size[2] - 1)
                    field.values[voxelIndex(size, x, y, z)] = -1.0f;
    return field;
}

} // namespace knit
