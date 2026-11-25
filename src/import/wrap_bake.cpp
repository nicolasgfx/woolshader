// The colour bake. Each surface triangle is sampled on a fine lattice and the last sample to land on a texel owns it. A texel
// takes the colour where a ray from just outside first meets the source, or the closest point of the source if it meets nothing.
#include "import/wrap_stages.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <map>
#include <numeric>
#include <queue>
#include <utility>

#include "io/image_io.h"
#include "model/cancel.h"
#include "model/colour.h"
#include "model/mesh.h"
#include "model/parallel.h"

namespace knit {

namespace {

constexpr int kLatticeSteps = 12;
constexpr int kBvhLeafTriangles = 4;
const Vec3 kColourWithoutMaterial{0.5, 0.5, 0.5};  // sRGB values

// A material's albedo, in sRGB values as stored: an image, one colour, or the grey of no material.
class SourceColours {
public:
    SourceColours(const WrapSource& source, const ImportedMesh& file, const std::vector<ColourImage>& textures)
        : source_(source), file_(file) {
        std::map<std::string, const ColourImage*> byPath;
        for (const ColourImage& t : textures) byPath[t.path] = &t;
        imageOf_.assign(file.materials.size(), nullptr);
        for (size_t m = 0; m < file.materials.size(); ++m) {
            const auto it = byPath.find(file.materials[m].texture);
            if (it != byPath.end()) imageOf_[m] = it->second;
        }
    }

    // At weights (w0, w1, w2) of the soup's triangle t. Stored vertex colours come before the material, as in the application's reader.
    Vec3 at(int t, double w0, double w1, double w2) const {
        if (!file_.vcolor.empty()) {
            const std::array<int, 3>& v = source_.vertex[t];
            return file_.vcolor[v[0]] * w0 + file_.vcolor[v[1]] * w1 + file_.vcolor[v[2]] * w2;
        }
        const int m = source_.material[t];
        if (m < 0 || m >= static_cast<int>(file_.materials.size())) return kColourWithoutMaterial;
        if (imageOf_[m] && !source_.uv.empty()) {
            const std::array<Vec2, 3>& q = source_.uv[t];
            return texel(*imageOf_[m], w0 * q[0].x + w1 * q[1].x + w2 * q[2].x, w0 * q[0].y + w1 * q[1].y + w2 * q[2].y);
        }
        if (file_.materials[m].hasKd) {
            const Vec3& kd = file_.materials[m].kd;
            return {srgbEncode(std::clamp(kd.x, 0.0, 1.0)), srgbEncode(std::clamp(kd.y, 0.0, 1.0)), srgbEncode(std::clamp(kd.z, 0.0, 1.0))};
        }
        return kColourWithoutMaterial;
    }

private:
    // The nearest texel, repeating, with v counted from the image's bottom row. Rows are stored top first.
    static Vec3 texel(const ColourImage& image, double u, double v) {
        const double fractionU = u - std::floor(u), fractionV = v - std::floor(v);
        const int x = static_cast<int>(std::floor(fractionU * image.width)) % image.width;
        const int fromBottom = static_cast<int>(std::floor(fractionV * image.height)) % image.height;
        const uint8_t* p = &image.rgba[(static_cast<size_t>(image.height - 1 - fromBottom) * image.width + x) * 4];
        return {p[0] / 255.0, p[1] / 255.0, p[2] / 255.0};
    }

    const WrapSource& source_;
    const ImportedMesh& file_;
    std::vector<const ColourImage*> imageOf_;
};

struct SurfaceTexels {
    std::vector<int> index;  // texel numbers, ascending
    std::vector<Vec3> position, direction;
};

// Where each texel lies on the surface. The last sample to land on a texel owns it.
SurfaceTexels texelsOnSurface(const QuadMesh& mesh, const QuadUv& uv, const std::vector<Vec3>& normals, int side) {
    const size_t texels = static_cast<size_t>(side) * side;
    std::vector<Vec3> position(texels), direction(texels);
    std::vector<char> filled(texels, 0);
    for (size_t f = 0; f < mesh.quad.size(); ++f) {
        if ((f & 0xfff) == 0) cancelCheckpoint();
        const std::array<int, 4>& q = mesh.quad[f];
        const std::array<Vec3, 4> corners = {mesh.pos[q[0]], mesh.pos[q[1]], mesh.pos[q[2]], mesh.pos[q[3]]};
        for (const std::array<int, 3>& t : quadSplit(corners)) {
            const Vec3 p[3] = {corners[t[0]], corners[t[1]], corners[t[2]]};
            const Vec3 m[3] = {normals[q[t[0]]], normals[q[t[1]]], normals[q[t[2]]]};
            const Vec2 c[3] = {uv[f][t[0]], uv[f][t[1]], uv[f][t[2]]};
            for (int i = 0; i <= kLatticeSteps; ++i)
                for (int j = 0; i + j <= kLatticeSteps; ++j) {
                    const double w0 = static_cast<double>(i) / kLatticeSteps;
                    const double w1 = static_cast<double>(j) / kLatticeSteps;
                    const double w2 = static_cast<double>(kLatticeSteps - i - j) / kLatticeSteps;
                    const int x = std::clamp(static_cast<int>((w0 * c[0].x + w1 * c[1].x + w2 * c[2].x) * side), 0, side - 1);
                    const int y = std::clamp(static_cast<int>((w0 * c[0].y + w1 * c[1].y + w2 * c[2].y) * side), 0, side - 1);
                    const size_t at = static_cast<size_t>(y) * side + x;
                    position[at] = p[0] * w0 + p[1] * w1 + p[2] * w2;
                    direction[at] = m[0] * w0 + m[1] * w1 + m[2] * w2;
                    filled[at] = 1;
                }
        }
    }
    SurfaceTexels out;
    for (size_t i = 0; i < texels; ++i)
        if (filled[i]) {
            out.index.push_back(static_cast<int>(i));
            out.position.push_back(position[i]);
            out.direction.push_back(direction[i]);
        }
    return out;
}

// The squared distance from q to an axis-aligned box, 0 inside it.
double squaredDistanceToBox(const Vec3& q, const Vec3& lo, const Vec3& hi) {
    const double dx = std::max({lo.x - q.x, 0.0, q.x - hi.x});
    const double dy = std::max({lo.y - q.y, 0.0, q.y - hi.y});
    const double dz = std::max({lo.z - q.z, 0.0, q.z - hi.z});
    return dx * dx + dy * dy + dz * dz;
}

// Whether the ray o + t d, 0 <= t <= limit, passes through the box (slab test; `inverse` holds 1 / d per axis).
bool meetsBox(const Vec3& lo, const Vec3& hi, const Vec3& o, const Vec3& inverse, double limit) {
    double near = 0, far = limit;
    const double tx1 = (lo.x - o.x) * inverse.x, tx2 = (hi.x - o.x) * inverse.x;
    near = std::max(near, std::min(tx1, tx2)), far = std::min(far, std::max(tx1, tx2));
    const double ty1 = (lo.y - o.y) * inverse.y, ty2 = (hi.y - o.y) * inverse.y;
    near = std::max(near, std::min(ty1, ty2)), far = std::min(far, std::max(ty1, ty2));
    const double tz1 = (lo.z - o.z) * inverse.z, tz2 = (hi.z - o.z) * inverse.z;
    near = std::max(near, std::min(tz1, tz2)), far = std::min(far, std::max(tz1, tz2));
    return near <= far;
}

// The parameter s in [0, 1] of the point a + s (b - a) closest to q (Ericson 2005, 5.1.2). A zero-length segment gives 0.
double closestOnSegment(const Vec3& q, const Vec3& a, const Vec3& b) {
    const Vec3 ab = b - a;
    const double lengthSq = dot(ab, ab);
    if (!(lengthSq > 0)) return 0.0;
    return std::clamp(dot(q - a, ab) / lengthSq, 0.0, 1.0);
}

double squaredLength(const Vec3& v) { return dot(v, v); }

using BoxQueue = std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>>;

} // namespace

// --- the hierarchy

TriangleBvh::TriangleBvh(const TriMesh& soup) : soup_(soup) {
    order_.resize(soup.tri.size());
    std::iota(order_.begin(), order_.end(), 0);
    centre_.resize(soup.tri.size());
    for (size_t t = 0; t < soup.tri.size(); ++t) centre_[t] = (soup.pos[soup.tri[t][0]] + soup.pos[soup.tri[t][1]] + soup.pos[soup.tri[t][2]]) / 3.0;
    nodes_.reserve(soup.tri.size());
    build(0, static_cast<int>(soup.tri.size()));
}

int TriangleBvh::build(int begin, int end) {
    const int index = static_cast<int>(nodes_.size());
    nodes_.emplace_back();
    Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300}, centreLo = lo, centreHi = hi;
    for (int i = begin; i < end; ++i) {
        for (int v : soup_.tri[order_[i]]) {
            lo = componentMin(lo, soup_.pos[v]);
            hi = componentMax(hi, soup_.pos[v]);
        }
        centreLo = componentMin(centreLo, centre_[order_[i]]);
        centreHi = componentMax(centreHi, centre_[order_[i]]);
    }
    nodes_[index].lo = lo;
    nodes_[index].hi = hi;
    if (end - begin <= kBvhLeafTriangles) {
        nodes_[index].first = begin;
        nodes_[index].count = end - begin;
        return index;
    }
    const Vec3 size = centreHi - centreLo;
    const int axis = size.x >= size.y && size.x >= size.z ? 0 : size.y >= size.z ? 1 : 2;
    const int middle = (begin + end) / 2;
    std::nth_element(order_.begin() + begin, order_.begin() + middle, order_.begin() + end, [&](int a, int b) {
        const Vec3& p = centre_[a];
        const Vec3& q = centre_[b];
        const double x = axis == 0 ? p.x : axis == 1 ? p.y : p.z, y = axis == 0 ? q.x : axis == 1 ? q.y : q.z;
        return x != y ? x < y : a < b;
    });
    build(begin, middle);
    const int right = build(middle, end);
    nodes_[index].right = right;
    return index;
}

TriangleBvh::Hit TriangleBvh::firstHit(const Vec3& origin, const Vec3& direction, double maxT) const {
    Hit best;
    double limit = maxT;
    const Vec3 inverse{1.0 / direction.x, 1.0 / direction.y, 1.0 / direction.z};
    int stack[64];
    int top = 0;
    stack[top++] = 0;
    while (top > 0) {
        const Node& n = nodes_[stack[--top]];
        if (!meetsBox(n.lo, n.hi, origin, inverse, limit)) continue;
        if (n.count == 0) {
            stack[top++] = n.right;
            stack[top++] = static_cast<int>(&n - nodes_.data()) + 1;
            continue;
        }
        for (int i = n.first; i < n.first + n.count; ++i) intersect(order_[i], origin, direction, limit, best);
    }
    return best;
}

// Moller and Trumbore 1997: the ray's parameter and the barycentric weights from one determinant.
void TriangleBvh::intersect(int t, const Vec3& o, const Vec3& d, double& limit, Hit& best) const {
    const Vec3& a = soup_.pos[soup_.tri[t][0]];
    const Vec3 e1 = soup_.pos[soup_.tri[t][1]] - a, e2 = soup_.pos[soup_.tri[t][2]] - a;
    const Vec3 p = cross(d, e2);
    const double det = dot(e1, p);
    if (std::abs(det) < 1e-14) return;
    const double inverse = 1.0 / det;
    const Vec3 s = o - a;
    const double b1 = dot(s, p) * inverse;
    if (b1 < 0 || b1 > 1) return;
    const Vec3 q = cross(s, e1);
    const double b2 = dot(d, q) * inverse;
    if (b2 < 0 || b1 + b2 > 1) return;
    const double along = dot(e2, q) * inverse;
    if (along < 0 || along > limit) return;
    limit = along;
    best = {t, along, b1, b2};
}

// Triangle t's closest point to q, kept if it is nearer than `best`, or as near with a lower index. A triangle with area takes
// Ericson's region test (2005, 5.1.5: closestOnTriangle in model/mesh.h); one of zero area the closest point of its three edges.
void TriangleBvh::closestOnTriangleOfSoup(int t, const Vec3& q, Nearest& best) const {
    const Vec3& a = soup_.pos[soup_.tri[t][0]];
    const Vec3& b = soup_.pos[soup_.tri[t][1]];
    const Vec3& c = soup_.pos[soup_.tri[t][2]];
    double b1 = 0, b2 = 0;
    if (squaredLength(cross(b - a, c - a)) > 0) {
        const TrianglePoint point = closestOnTriangle(q, a, b, c);
        b1 = point.bary[1];
        b2 = point.bary[2];
    } else {
        const double alongAb = closestOnSegment(q, a, b), alongBc = closestOnSegment(q, b, c), alongCa = closestOnSegment(q, c, a);
        const double onAb = squaredLength(q - (a + (b - a) * alongAb));
        const double onBc = squaredLength(q - (b + (c - b) * alongBc));
        const double onCa = squaredLength(q - (c + (a - c) * alongCa));
        if (onAb <= onBc && onAb <= onCa) {
            b1 = alongAb;
        } else if (onBc <= onCa) {
            b1 = 1.0 - alongBc;
            b2 = alongBc;
        } else {
            b2 = 1.0 - alongCa;
        }
    }
    const double distSq = squaredLength(q - (a + (b - a) * b1 + (c - a) * b2));
    if (best.tri < 0 || distSq < best.distSq || (distSq == best.distSq && t < best.tri)) best = {t, distSq, b1, b2};
}

// Best first: boxes come off a min-heap by their distance to q, and a box farther than the best point so far is not opened. A
// box exactly as far is still opened, so that an equally near triangle of lower index is found.
TriangleBvh::Nearest TriangleBvh::closestPoint(const Vec3& q) const {
    Nearest best;
    if (soup_.tri.empty() || !std::isfinite(q.x) || !std::isfinite(q.y) || !std::isfinite(q.z)) return best;
    BoxQueue queue;
    queue.push({squaredDistanceToBox(q, nodes_[0].lo, nodes_[0].hi), 0});
    while (!queue.empty()) {
        const std::pair<double, int> box = queue.top();
        queue.pop();
        if (best.tri >= 0 && box.first > best.distSq) break;  // every box left is at least as far
        const Node& n = nodes_[box.second];
        if (n.count > 0) {
            for (int i = n.first; i < n.first + n.count; ++i) closestOnTriangleOfSoup(order_[i], q, best);
            continue;
        }
        const int left = box.second + 1;
        queue.push({squaredDistanceToBox(q, nodes_[left].lo, nodes_[left].hi), left});
        queue.push({squaredDistanceToBox(q, nodes_[n.right].lo, nodes_[n.right].hi), n.right});
    }
    return best;
}

// --- normals and the bake

std::vector<Vec3> vertexNormalsOf(const FaceGraph& graph) {
    const QuadMesh& mesh = *graph.mesh;
    std::vector<Vec3> sum(mesh.pos.size());
    std::vector<int> largestFace(mesh.pos.size(), -1);
    for (int f = 0; f < graph.faceCount(); ++f) {
        const std::array<int, 4>& q = mesh.quad[f];
        for (int k = 0; k < 4; ++k) {
            const int v = q[k];
            if (largestFace[v] < 0 || graph.area[f] > graph.area[largestFace[v]]) largestFace[v] = f;
            const Vec3 toPrev = mesh.pos[q[(k + 3) % 4]] - mesh.pos[v];
            const Vec3 toNext = mesh.pos[q[(k + 1) % 4]] - mesh.pos[v];
            const double prevLength = norm(toPrev), nextLength = norm(toNext);
            if (!(prevLength > 0) || !(nextLength > 0)) continue;  // a zero-length edge adds nothing
            const double angle = std::acos(std::clamp(dot(toPrev / prevLength, toNext / nextLength), -1.0, 1.0));
            sum[v] += graph.normal[f] * angle;
        }
    }
    std::vector<Vec3> normals(mesh.pos.size());
    for (size_t v = 0; v < mesh.pos.size(); ++v) {
        const double length = norm(sum[v]);
        if (length > 0) normals[v] = sum[v] / length;
        else if (largestFace[v] >= 0) normals[v] = graph.normal[largestFace[v]];
    }
    return normals;
}

BakedColours bakeSurfaceColours(const FaceGraph& graph, const QuadUv& uv, const WrapSource& source, const ImportedMesh& file,
                                const std::vector<ColourImage>& textures, double rayMm) {
    const QuadMesh& mesh = *graph.mesh;
    BakedColours baked;
    const int side = baked.side;
    const SurfaceTexels texels = texelsOnSurface(mesh, uv, vertexNormalsOf(graph), side);

    TriMesh soup;
    soup.pos.reserve(source.soup.size() * 3);
    soup.tri.reserve(source.soup.size());
    for (const TriangleCorners& t : source.soup) {
        const int first = static_cast<int>(soup.pos.size());
        soup.pos.insert(soup.pos.end(), t.begin(), t.end());
        soup.tri.push_back({first, first + 1, first + 2});
    }
    const TriangleBvh bvh(soup);
    const SourceColours colours(source, file, textures);

    baked.rgb.assign(static_cast<size_t>(side) * side * 3, 0.0f);
    baked.filled.assign(static_cast<size_t>(side) * side, 0);
    baked.texels = static_cast<long long>(texels.index.size());
    std::atomic<long long> missed{0};
    ExecPolicy parallel;
    parallel.parallel = true;
    parallelFor(baked.texels, 4096,
                [&](int64_t begin, int64_t end, int64_t) {
                    cancelCheckpoint();
                    for (int64_t k = begin; k < end; ++k) {
                        const Vec3& p = texels.position[k];
                        const Vec3 d = texels.direction[k] / std::max(norm(texels.direction[k]), 1e-9);
                        const TriangleBvh::Hit hit = bvh.firstHit(p + d * rayMm, -d, 2.0 * rayMm);
                        Vec3 colour = kColourWithoutMaterial;
                        if (hit.tri >= 0) {
                            colour = colours.at(hit.tri, 1.0 - hit.b1 - hit.b2, hit.b1, hit.b2);
                        } else {  // a sealed gap with nothing behind it: the closest point of the source
                            missed.fetch_add(1, std::memory_order_relaxed);
                            const TriangleBvh::Nearest nearest = bvh.closestPoint(p);
                            if (nearest.tri >= 0) colour = colours.at(nearest.tri, 1.0 - nearest.b1 - nearest.b2, nearest.b1, nearest.b2);
                        }
                        const size_t at = static_cast<size_t>(texels.index[k]);
                        baked.rgb[at * 3] = static_cast<float>(colour.x);
                        baked.rgb[at * 3 + 1] = static_cast<float>(colour.y);
                        baked.rgb[at * 3 + 2] = static_cast<float>(colour.z);
                        baked.filled[at] = 1;
                    }
                },
                parallel);
    baked.closestTexels = missed;
    return baked;
}

std::vector<float> bledColours(const BakedColours& colours, int rounds) {
    const int side = colours.side;
    std::vector<float> rgb = colours.rgb;
    std::vector<char> filled = colours.filled;
    const int steps[4][2] = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
    for (int round = 0; round < rounds; ++round) {
        std::vector<float> next = rgb;
        std::vector<char> nextFilled = filled;
        for (int y = 0; y < side; ++y)
            for (int x = 0; x < side; ++x) {
                const size_t at = static_cast<size_t>(y) * side + x;
                if (filled[at]) continue;
                float total[3] = {0.0f, 0.0f, 0.0f};
                float count = 0.0f;
                for (const auto& step : steps) {  // the texture wraps around its edges
                    const int nx = (x + step[0] + side) % side, ny = (y + step[1] + side) % side;
                    const size_t from = static_cast<size_t>(ny) * side + nx;
                    if (!filled[from]) continue;
                    for (int c = 0; c < 3; ++c) total[c] += rgb[from * 3 + c];
                    count += 1.0f;
                }
                if (count == 0.0f) continue;
                for (int c = 0; c < 3; ++c) next[at * 3 + c] = total[c] / count;
                nextFilled[at] = 1;
            }
        rgb.swap(next);
        filled.swap(nextFilled);
    }
    return rgb;
}

Vec3 meanLinearColourOf(const BakedColours& colours) {
    Vec3 sum;
    long long count = 0;
    for (size_t i = 0; i < colours.filled.size(); ++i) {
        if (!colours.filled[i]) continue;
        sum += Vec3{srgbDecode(colours.rgb[i * 3]), srgbDecode(colours.rgb[i * 3 + 1]), srgbDecode(colours.rgb[i * 3 + 2])};
        ++count;
    }
    return count ? sum / static_cast<double>(count) : Vec3{};
}

} // namespace knit
