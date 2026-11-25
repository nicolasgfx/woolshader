// Stage 6: the smoothed surface rebuilt as even quads through a level set (third_party/vdb), then given one consistent outward
// winding. Also the sorted list of edges the later stages share (quadSplit is in wrap_quad_split.cpp).
#include "import/wrap_stages.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "io/files.h"
#include "level_set_remesh.h"

namespace knit {

namespace {

constexpr double kFlatComponent = 1e-9;        // a part whose |volume| is below this share of its box volume has no outside

// One side of a quad, found under its undirected key (lower vertex, higher vertex).
struct QuadSide {
    int low = 0, high = 0;
    int quad = 0;
    bool runsUp = false;  // in the quad's own corner order, the side runs from `low` to `high`
};

bool sideComesFirst(const QuadSide& a, const QuadSide& b) {
    if (a.low != b.low) return a.low < b.low;
    if (a.high != b.high) return a.high < b.high;
    return a.quad < b.quad;
}

// For each quad side shared by exactly two quads, the other quad and whether the two sides run the same way in the input. Sides
// used by one quad or by more than two link nothing: such edges are not manifold, and stage 7 refuses them.
struct QuadLink {
    int quad = 0;
    int other = 0;
    bool sameDirection = false;
};

std::vector<QuadLink> linksBetweenQuads(const std::vector<std::array<int, 4>>& quads) {
    std::vector<QuadSide> sides;
    sides.reserve(quads.size() * 4);
    for (int q = 0; q < static_cast<int>(quads.size()); ++q)
        for (int k = 0; k < 4; ++k) {
            const int a = quads[q][k], b = quads[q][(k + 1) % 4];
            sides.push_back({std::min(a, b), std::max(a, b), q, a < b});
        }
    std::sort(sides.begin(), sides.end(), sideComesFirst);
    std::vector<QuadLink> links;
    size_t run = 0;
    while (run < sides.size()) {
        size_t end = run + 1;
        while (end < sides.size() && sides[end].low == sides[run].low && sides[end].high == sides[run].high) ++end;
        if (end - run == 2 && sides[run].quad != sides[run + 1].quad) {
            const bool same = sides[run].runsUp == sides[run + 1].runsUp;
            links.push_back({sides[run].quad, sides[run + 1].quad, same});
            links.push_back({sides[run + 1].quad, sides[run].quad, same});
        }
        run = end;
    }
    std::sort(links.begin(), links.end(), [](const QuadLink& a, const QuadLink& b) { return a.quad != b.quad ? a.quad < b.quad : a.other < b.other; });
    return links;
}

void reverseQuad(std::array<int, 4>& q) { std::swap(q[1], q[3]); }  // (a, b, c, d) -> (a, d, c, b)

// Six times the signed volume enclosed by the component's quads, each as the triangles (0, 1, 2) and (0, 2, 3), measured from
// `centre`. On a closed surface the result does not depend on the centre; near the part it loses less to rounding.
double sixTimesVolume(const QuadMesh& mesh, const std::vector<int>& members, const Vec3& centre) {
    double sum = 0;
    for (int q : members) {
        const std::array<int, 4>& c = mesh.quad[q];
        const Vec3 p0 = mesh.pos[c[0]] - centre, p1 = mesh.pos[c[1]] - centre, p2 = mesh.pos[c[2]] - centre, p3 = mesh.pos[c[3]] - centre;
        sum += dot(p0, cross(p1, p2)) + dot(p0, cross(p2, p3));
    }
    return sum;
}

// Gives every connected part of the quads one winding, outward. Breadth first from each part's lowest quad, a neighbour whose
// shared side runs the same way as the current quad's is reversed; meeting a quad again the wrong way means the part is not
// orientable. Then a part with a negative signed volume is reversed whole.
bool orientOutward(QuadMesh& mesh, std::string* err) {
    const int quadCount = static_cast<int>(mesh.quad.size());
    const std::vector<QuadLink> links = linksBetweenQuads(mesh.quad);
    std::vector<int> firstLink(quadCount + 1, 0);
    for (const QuadLink& link : links) ++firstLink[link.quad + 1];
    for (int q = 0; q < quadCount; ++q) firstLink[q + 1] += firstLink[q];

    std::vector<int> component(quadCount, -1);
    std::vector<char> reversed(quadCount, 0);
    std::vector<std::vector<int>> members;
    for (int seed = 0; seed < quadCount; ++seed) {
        if (component[seed] >= 0) continue;
        const int part = static_cast<int>(members.size());
        members.emplace_back();
        component[seed] = part;
        members[part].push_back(seed);
        for (size_t next = 0; next < members[part].size(); ++next) {
            const int q = members[part][next];
            for (int i = firstLink[q]; i < firstLink[q + 1]; ++i) {
                const QuadLink& link = links[i];
                // The other quad agrees when, after both reversals, the shared side runs opposite ways.
                const bool otherReversed = link.sameDirection != (reversed[q] != 0);
                if (component[link.other] < 0) {
                    component[link.other] = part;
                    reversed[link.other] = otherReversed;
                    members[part].push_back(link.other);
                } else if ((reversed[link.other] != 0) != otherReversed) {
                    return fail(err, "the remeshed surface is not orientable: quad " + std::to_string(link.other) + " cannot agree with its neighbours");
                }
            }
        }
    }
    for (int q = 0; q < quadCount; ++q)
        if (reversed[q]) reverseQuad(mesh.quad[q]);

    for (size_t part = 0; part < members.size(); ++part) {
        Vec3 low(1e300, 1e300, 1e300), high(-1e300, -1e300, -1e300);
        for (int q : members[part])
            for (int v : mesh.quad[q]) {
                low = componentMin(low, mesh.pos[v]);
                high = componentMax(high, mesh.pos[v]);
            }
        const Vec3 size = high - low;
        const double volume = sixTimesVolume(mesh, members[part], (low + high) * 0.5) / 6.0;
        if (!(std::fabs(volume) > kFlatComponent * size.x * size.y * size.z))
            return fail(err, "cannot tell the outside of a remeshed part of " + std::to_string(members[part].size()) + " quads: it encloses no volume");
        if (volume < 0)
            for (int q : members[part]) reverseQuad(mesh.quad[q]);
    }
    return true;
}

} // namespace

std::vector<std::array<int, 2>> sortedEdges(const std::vector<std::array<int, 4>>& quads) {
    std::vector<std::array<int, 2>> edges;
    edges.reserve(quads.size() * 4);
    for (const std::array<int, 4>& q : quads)
        for (int k = 0; k < 4; ++k) edges.push_back({std::min(q[k], q[(k + 1) % 4]), std::max(q[k], q[(k + 1) % 4])});
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    return edges;
}

bool remeshToEvenQuads(const QuadMesh& surface, double voxelMm, QuadMesh& remeshed, std::string* err) {
    remeshed = QuadMesh{};
    LevelSetInput input;
    input.voxelSize = static_cast<float>(voxelMm);
    input.positions.reserve(surface.pos.size());
    for (const Vec3& p : surface.pos) input.positions.push_back({static_cast<float>(p.x), static_cast<float>(p.y), static_cast<float>(p.z)});
    input.triangles.reserve(surface.quad.size() * 2);
    for (const std::array<int, 4>& q : surface.quad) {
        const std::array<Vec3, 4> corners = {surface.pos[q[0]], surface.pos[q[1]], surface.pos[q[2]], surface.pos[q[3]]};
        for (const std::array<int, 3>& t : quadSplit(corners))
            input.triangles.push_back({static_cast<uint32_t>(q[t[0]]), static_cast<uint32_t>(q[t[1]]), static_cast<uint32_t>(q[t[2]])});
    }

    LevelSetOutput output;
    if (!remeshThroughLevelSet(input, output, err)) return false;
    if (!output.triangles.empty())
        return fail(err, "the level set remesh gave " + std::to_string(output.triangles.size()) + " triangles besides its quads");

    remeshed.pos.reserve(output.positions.size());
    for (const std::array<float, 3>& p : output.positions) remeshed.pos.push_back({p[0], p[1], p[2]});
    remeshed.quad.reserve(output.quads.size());
    for (const std::array<uint32_t, 4>& q : output.quads)
        remeshed.quad.push_back({static_cast<int>(q[0]), static_cast<int>(q[1]), static_cast<int>(q[2]), static_cast<int>(q[3])});
    if (!orientOutward(remeshed, err)) {
        remeshed = QuadMesh{};
        return false;
    }
    return true;
}

} // namespace knit
