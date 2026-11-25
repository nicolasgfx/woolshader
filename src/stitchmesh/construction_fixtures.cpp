#include "stitchmesh/construction_fixtures.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>

#include "model/mesh_query.h"
#include "stitchmesh/guided_rows.h"
#include "stitchmesh/stitch_fixtures.h"
#include "stitchmesh/stitch_mesh.h"

namespace knit {

namespace {

constexpr double kTau = 6.283185307179586;

TriMesh cylinder(double r, double z0, double z1, double step) {
    TriMesh t;
    const int around = 192, rings = std::max(2, static_cast<int>(std::ceil((z1 - z0) / step)) + 1);
    for (int i = 0; i < rings; ++i) {
        const double z = z0 + (z1 - z0) * i / (rings - 1);
        for (int j = 0; j < around; ++j) t.pos.push_back({std::cos(kTau * j / around) * r, std::sin(kTau * j / around) * r, z});
    }
    for (int i = 0; i + 1 < rings; ++i)
        for (int j = 0; j < around; ++j) {
            const int a = i * around + j, b = i * around + (j + 1) % around, c = (i + 1) * around + (j + 1) % around, d = (i + 1) * around + j;
            t.tri.push_back({a, b, c});
            t.tri.push_back({a, c, d});
        }
    t.computeNormals();
    return t;
}

void bindAll(TypedStitchMesh& s, const TriMesh& reference) {
    const MeshQuery q(reference);
    s.binding.clear();
    for (const Vec3& p : s.mesh.pos) {
        const MeshQuery::Hit h = q.closest(p);
        VertexBinding b;
        b.tri = h.tri;
        b.distance = h.dist;
        if (h.tri >= 0) {
            const TrianglePoint tp = closestOnTriangle(p, reference.pos[reference.tri[h.tri][0]], reference.pos[reference.tri[h.tri][1]],
                                                       reference.pos[reference.tri[h.tri][2]]);
            for (int k = 0; k < 3; ++k) b.bary[k] = tp.bary[k];
        }
        s.binding.push_back(b);
    }
}

std::vector<Vec3> ringLoop(const TypedStitchMesh& s, int k) {
    std::vector<Vec3> l;
    for (int v : s.rings[k]) l.push_back(s.mesh.pos[v]);
    return l;
}

FixtureCase tube(const std::string& name, const std::string& entry, const std::string& purpose, const std::vector<int>& counts,
                 const std::vector<std::vector<int>>& shaping) {
    FixtureCase c;
    c.name = name;
    c.entry = entry;
    c.purpose = purpose;
    std::string err;
    ringTubeFixture(counts, shaping, 7.0, 5.5, c.mesh, c.reference, c.declared, &err);
    return c;
}

// The special-quad neighbourhood: the recovery branch's coarse fixture 2 (4 around, the middle row short, two row-end
// triangles) through the ported subdivision; the reference is the coarse mesh's own facets.
FixtureCase specialQuadCase() {
    FixtureCase c;
    c.name = "special-quad-neighbourhood";
    c.entry = "special quad (3.8)";
    c.purpose = "dbf5f41 SpecialQuads coarse fixture subdivided by buildStitchMesh: two special quads where the short middle coarse row ends; "
                "the yarn is traced through them";
    const CoarseFixture fx = stitchFixture(StitchFixture::SpecialQuads);
    StitchMeshProduct p;
    std::string err;
    buildStitchMesh(fx.mesh, fx.top, p, &err);
    TriMesh ref;
    ref.pos = fx.mesh.pos;
    for (int f = 0; f < fx.mesh.faces(); ++f)
        for (int k = 1; k + 1 < fx.mesh.sides(f); ++k) ref.tri.push_back({fx.mesh.vertexAt(f, 0), fx.mesh.vertexAt(f, k), fx.mesh.vertexAt(f, k + 1)});
    ref.computeNormals();
    c.reference = ref;
    // the coarse border loops: z = 0 (cast-on) and the top ring (bind-off)
    double zTop = 0;
    for (const Vec3& q : fx.mesh.pos) zTop = std::max(zTop, q.z);
    DeclaredBoundary lo, hi;
    lo.id = "z-low";
    lo.role = EdgeBoundary::CastOn;
    hi.id = "z-high";
    hi.role = EdgeBoundary::BindOff;
    for (const Vec3& q : fx.mesh.pos) {
        if (q.z < 1e-9) lo.loop.push_back(q);
        if (q.z > zTop - 1e-9) hi.loop.push_back(q);
    }
    c.declared = {lo, hi};
    c.mesh = typedFromStitchMesh(p.stitch, c.declared, 1e-6);
    bindAll(c.mesh, ref);
    return c;
}

} // namespace

bool ringTubeFixture(const std::vector<int>& counts, const std::vector<std::vector<int>>& shaping, double radius, double rowHeight,
                     TypedStitchMesh& mesh, TriMesh& reference, std::vector<DeclaredBoundary>& declared, std::string* err, bool boundaryShaping) {
    mesh = TypedStitchMesh{};
    const int R = static_cast<int>(counts.size()) - 1;
    if (R < 1 || static_cast<int>(shaping.size()) != R) {
        if (err) *err = "ringTubeFixture: counts and shaping disagree";
        return false;
    }
    for (int k = 0; k <= R; ++k) {
        std::vector<int> ring;
        for (int j = 0; j < counts[k]; ++j) {
            ring.push_back(static_cast<int>(mesh.mesh.pos.size()));
            mesh.mesh.pos.push_back({std::cos(kTau * j / counts[k]) * radius, std::sin(kTau * j / counts[k]) * radius, rowHeight * k});
        }
        mesh.rings.push_back(ring);
    }
    for (int k = 0; k < R; ++k) {
        std::vector<int> slots = shaping[k];
        std::sort(slots.begin(), slots.end());
        if (!assembleStrip(k, mesh.rings[k], mesh.rings[k + 1], slots, k == 0 ? EdgeBoundary::CastOn : EdgeBoundary::Interior,
                           k == R - 1 ? EdgeBoundary::BindOff : EdgeBoundary::Interior, mesh, nullptr, err, boundaryShaping))
            return false;
        mesh.shiftPath.push_back(ShiftPathEdge{k, mesh.rings[k][0], mesh.rings[k + 1][0]});
    }
    reference = cylinder(radius, -rowHeight, rowHeight * (R + 1), rowHeight / 4);
    bindAll(mesh, reference);
    declared.clear();
    DeclaredBoundary lo, hi;
    lo.id = "ring-0";
    lo.role = EdgeBoundary::CastOn;
    lo.loop = ringLoop(mesh, 0);
    hi.id = "ring-" + std::to_string(R);
    hi.role = EdgeBoundary::BindOff;
    hi.loop = ringLoop(mesh, R);
    declared = {lo, hi};
    return true;
}

TriMesh cylinderReference(double r, double z0, double z1, double step) { return cylinder(r, z0, z1, step); }

TriMesh revolutionThroughRings(const std::vector<double>& rho, const std::vector<double>& z, int around, int sub) {
    TriMesh m;
    std::vector<double> rr, zz;
    for (size_t k = 0; k + 1 < rho.size(); ++k)
        for (int s = 0; s < sub; ++s) {
            const double f = static_cast<double>(s) / sub;
            rr.push_back(rho[k] * (1 - f) + rho[k + 1] * f);
            zz.push_back(z[k] * (1 - f) + z[k + 1] * f);
        }
    rr.push_back(rho.back());
    zz.push_back(z.back());
    for (size_t k = 0; k < rr.size(); ++k)
        for (int j = 0; j < around; ++j) m.pos.push_back({rr[k] * std::cos(kTau * j / around), rr[k] * std::sin(kTau * j / around), zz[k]});
    auto id = [&](int k, int j) { return k * around + ((j % around) + around) % around; };
    for (int k = 0; k + 1 < static_cast<int>(rr.size()); ++k)
        for (int j = 0; j < around; ++j) {
            m.tri.push_back({id(k, j), id(k, j + 1), id(k + 1, j + 1)});
            m.tri.push_back({id(k, j), id(k + 1, j + 1), id(k + 1, j)});
        }
    m.computeNormals();
    return m;
}

namespace {

// A flat disc of radius r in z = 0 (a polar grid, 96 around, rings `step` apart), counter-clockwise from +z.
TriMesh disc(double r, double step) {
    TriMesh t;
    const int around = 96, rings = std::max(2, static_cast<int>(std::ceil(r / step)));
    t.pos.push_back({0, 0, 0});
    for (int k = 1; k <= rings; ++k)
        for (int j = 0; j < around; ++j) t.pos.push_back({std::cos(kTau * j / around) * r * k / rings, std::sin(kTau * j / around) * r * k / rings, 0});
    auto id = [&](int k, int j) { return k == 0 ? 0 : 1 + (k - 1) * around + ((j % around) + around) % around; };
    for (int j = 0; j < around; ++j) t.tri.push_back({0, id(1, j), id(1, j + 1)});
    for (int k = 1; k < rings; ++k)
        for (int j = 0; j < around; ++j) {
            t.tri.push_back({id(k, j), id(k + 1, j), id(k + 1, j + 1)});
            t.tri.push_back({id(k, j), id(k + 1, j + 1), id(k, j + 1)});
        }
    t.computeNormals();
    return t;
}

} // namespace

bool poleDiscFixture(const std::vector<int>& counts, const std::vector<std::vector<int>>& shaping, double rowHeight, TypedStitchMesh& mesh, TriMesh& reference,
                     std::vector<DeclaredBoundary>& declared, std::string* err) {
    mesh = TypedStitchMesh{};
    const int K = static_cast<int>(counts.size());
    if (K < 2 || static_cast<int>(shaping.size()) != K - 1) {
        if (err) *err = "poleDiscFixture: counts and shaping disagree";
        return false;
    }
    mesh.mesh.pos.push_back({0, 0, 0});
    mesh.rings.push_back({0});
    for (int k = 1; k <= K; ++k) {
        std::vector<int> ring;
        const int n = counts[k - 1];
        for (int j = 0; j < n; ++j) {
            ring.push_back(static_cast<int>(mesh.mesh.pos.size()));
            mesh.mesh.pos.push_back({std::cos(-kTau * j / n) * rowHeight * k, std::sin(-kTau * j / n) * rowHeight * k, 0});
        }
        mesh.rings.push_back(ring);
    }
    const std::vector<int>& r1 = mesh.rings[1];
    const int n1 = static_cast<int>(r1.size());
    for (int j = 0; j < n1; ++j) {
        mesh.mesh.addFace({0, r1[(j + 1) % n1], r1[j]});
        mesh.kind.push_back(StitchFaceKind::PoleStart);
        mesh.reason.push_back("");
        mesh.role.push_back({StitchEdgeRole::Wale, StitchEdgeRole::Top, StitchEdgeRole::Wale});
        mesh.boundary.push_back({EdgeBoundary::Interior, K == 1 ? EdgeBoundary::BindOff : EdgeBoundary::Interior, EdgeBoundary::Interior});
        mesh.strip.push_back(0);
        mesh.slot.push_back(j);
    }
    for (int k = 1; k < K; ++k) {
        std::vector<int> slots = shaping[k - 1];
        std::sort(slots.begin(), slots.end());
        if (!assembleStrip(k, mesh.rings[k], mesh.rings[k + 1], slots, EdgeBoundary::Interior, k == K - 1 ? EdgeBoundary::BindOff : EdgeBoundary::Interior, mesh, nullptr,
                           err))
            return false;
        mesh.shiftPath.push_back(ShiftPathEdge{k, mesh.rings[k][0], mesh.rings[k + 1][0]});
    }
    reference = disc(rowHeight * (K + 1), rowHeight / 4);
    bindAll(mesh, reference);
    DeclaredBoundary hi;
    hi.id = "ring-" + std::to_string(K);
    hi.role = EdgeBoundary::BindOff;
    hi.loop = ringLoop(mesh, K);
    declared = {hi};
    return true;
}

namespace {

FixtureCase poleDisc(const std::string& name, const std::string& entry, const std::string& purpose, const std::vector<int>& counts,
                     const std::vector<std::vector<int>>& shaping) {
    FixtureCase c;
    c.name = name;
    c.entry = entry;
    c.purpose = purpose;
    std::string err;
    poleDiscFixture(counts, shaping, 3.4, c.mesh, c.reference, c.declared, &err);
    c.runShift = true;
    c.expectShift = Status::Pass;
    return c;
}

} // namespace

std::vector<FixtureCase> armConstructionFixtureCases() {
    std::vector<FixtureCase> all;
    // positive: the gathered pole start, then a strip with two increases, then the bind-off
    all.push_back(poleDisc("pole-start-5-7-7", "pole start (proposed)",
                           "pole + ring of 5 (5 pole-start triangles), 5 -> 7 with two kp at slots 1 and 3, 7 -> 7 to the bind-off; the pole round "
                           "closed, the helix from strip 1",
                           {5, 7, 7}, {{1, 3}, {}}));
    // positive: adjacent kp pentagons (the adjacent-shaping variant's neighbourhood), on a tube and next to a pole start
    {
        FixtureCase c = tube("adjacent-kp-6-9", "adjacent shaping (variant)",
                             "ring 6 -> 9 with kp pentagons at slots 1, 2, 3 side by side (sharing wale edges), between K neighbours", {6, 6, 9, 9}, {{}, {1, 2, 3}, {}});
        c.runShift = true;
        c.expectShift = Status::Pass;
        all.push_back(c);
    }
    all.push_back(poleDisc("pole-start-5-8-adjacent", "pole start + adjacent shaping", "pole + ring of 5, 5 -> 8 with kp at slots 1, 2, 3 side by side, 8 -> 8",
                           {5, 8, 8}, {{1, 2, 3}, {}}));
    // negative controls
    {
        FixtureCase c = poleDisc("neg-pole-face-missing", "pole start", "one pole-start triangle removed: a hole at the pole", {5, 7, 7}, {{1, 3}, {}});
        TypedStitchMesh& s = c.mesh;
        TypedStitchMesh t;
        t.mesh.pos = s.mesh.pos;
        t.binding = s.binding;
        t.rings = s.rings;
        t.shiftPath = s.shiftPath;
        for (int f = 1; f < s.mesh.faces(); ++f) {
            std::vector<int> corners;
            for (int k = 0; k < s.mesh.sides(f); ++k) corners.push_back(s.mesh.vertexAt(f, k));
            t.mesh.addFace(corners);
            t.kind.push_back(s.kind[f]);
            t.role.push_back(s.role[f]);
            t.boundary.push_back(s.boundary[f]);
            t.reason.push_back(s.reason[f]);
            t.strip.push_back(s.strip[f]);
            t.slot.push_back(s.slot[f]);
        }
        s = t;
        c.negative = true;
        c.expectConnectivity = Status::Fail;
        c.expectShift = Status::Fail;
        all.push_back(c);
    }
    {
        FixtureCase c = poleDisc("neg-pole-face-typed-k", "pole start", "a pole-start triangle stored as K (kind against roles)", {5, 7, 7}, {{1, 3}, {}});
        c.mesh.kind[2] = StitchFaceKind::Regular;
        c.negative = true;
        c.expectConnectivity = Status::Fail;
        c.expectShift = Status::Fail;
        all.push_back(c);
    }
    {
        FixtureCase c = poleDisc("neg-pole-round-has-path-edge", "pole start", "a shift-path edge in the pole round (pole to ring 1 vertex 0)", {5, 7, 7}, {{1, 3}, {}});
        c.mesh.shiftPath.insert(c.mesh.shiftPath.begin(), ShiftPathEdge{0, 0, c.mesh.rings[1][0]});
        c.negative = true;
        c.expectConnectivity = Status::Pass;
        c.expectShift = Status::Fail;
        all.push_back(c);
    }
    {
        FixtureCase c = tube("neg-adjacent-kp-beside-path", "adjacent shaping", "ring 6 -> 9 with kp at slots 0, 1, 2: slot 0 is beside the path, which the lookup refuses",
                             {6, 6, 9, 9}, {{}, {0, 1, 2}, {}});
        c.runShift = true;
        c.negative = true;
        c.expectConnectivity = Status::Unresolved;
        c.expectShift = Status::Pass;
        all.push_back(c);
    }
    // boundary shaping: kp in the cast-on strip, d12k in the bind-off strip
    auto boundaryTube = [](const std::string& name, const std::string& purpose, const std::vector<int>& counts, const std::vector<std::vector<int>>& shaping,
                           bool variant) {
        FixtureCase c;
        c.name = name;
        c.entry = "boundary shaping (variant)";
        c.purpose = purpose;
        std::string err;
        ringTubeFixture(counts, shaping, 7.0, 5.5, c.mesh, c.reference, c.declared, &err, variant);
        c.runShift = true;
        c.expectShift = Status::Pass;
        return c;
    };
    all.push_back(boundaryTube("kp-cast-on-8-9", "ring 8 -> 9 in the cast-on strip, the kp pentagon at slot 4 with its bottom edge a cast-on edge", {8, 9, 9}, {{4}, {}},
                               true));
    all.push_back(boundaryTube("d12k-bind-off-9-8", "ring 9 -> 8 in the bind-off strip, the d12k pentagon at slot 4 with its top edge a bind-off edge", {9, 9, 8},
                               {{}, {4}}, true));
    {
        FixtureCase c = boundaryTube("neg-d12k-bind-off-without-variant", "the same d12k in the bind-off strip with the variant off: the CP2 rule types it Unsupported",
                                     {9, 9, 8}, {{}, {4}}, false);
        c.negative = true;
        c.expectConnectivity = Status::Unresolved;
        all.push_back(c);
    }
    {
        FixtureCase c = boundaryTube("neg-bind-off-on-d12k-bottom", "the bind-off role put on the d12k's two bottom edges (interior edges shared with strip 0)", {9, 9, 8},
                                     {{}, {4}}, true);
        for (int f = 0; f < c.mesh.mesh.faces(); ++f)
            if (c.mesh.kind[f] == StitchFaceKind::Decrease)
                for (size_t k = 0; k < c.mesh.role[f].size(); ++k)
                    if (c.mesh.role[f][k] == StitchEdgeRole::Bottom) c.mesh.boundary[f][k] = EdgeBoundary::BindOff;
        c.negative = true;
        c.expectConnectivity = Status::Fail;
        all.push_back(c);
    }
    return all;
}

std::vector<std::vector<int>> triangleBoundaryLoops(const TriMesh& m) {
    std::map<std::pair<int, int>, int> uses;
    for (const auto& t : m.tri)
        for (int k = 0; k < 3; ++k) uses[{std::min(t[k], t[(k + 1) % 3]), std::max(t[k], t[(k + 1) % 3])}] += 1;
    std::map<int, int> nextOf;
    for (const auto& t : m.tri)
        for (int k = 0; k < 3; ++k) {
            const int a = t[k], b = t[(k + 1) % 3];
            if (uses[{std::min(a, b), std::max(a, b)}] == 1) nextOf[a] = b;
        }
    std::vector<std::vector<int>> loops;
    std::set<int> seen;
    for (const auto& [v0, unused] : nextOf) {
        (void) unused;
        if (seen.count(v0)) continue;
        std::vector<int> loop;
        for (int v = v0; !seen.count(v) && nextOf.count(v); v = nextOf[v]) {
            seen.insert(v);
            loop.push_back(v);
        }
        loops.push_back(loop);
    }
    return loops;
}

std::vector<FixtureCase> constructionFixtureCases() {
    std::vector<FixtureCase> all;
    const std::vector<int> k8 = {8, 8, 8, 8};
    const std::vector<std::vector<int>> none = {{}, {}, {}};
    // positive fixtures
    {
        FixtureCase c = tube("k-open-tube-8x3", "K, cast-on, bind-off, shift-path (3.1, 3.5, 3.6, 3.7)",
                             "8 around, 3 strips: every bottom-ring edge CastOn, every top-ring edge BindOff, the vertex-0 wale edges as the shift-path; "
                             "virtual helix validated by tracing the virtual polygons",
                             k8, none);
        c.runShift = true;
        c.expectShift = Status::Pass;
        all.push_back(c);
    }
    {
        FixtureCase c = tube("kp-ring-8-9", "INC kp (3.3)", "ring 8 -> 9 in the middle strip, the kp pentagon at slot 4 (half a turn from the path)",
                             {8, 8, 9, 9}, {{}, {4}, {}});
        c.runShift = true;
        c.expectShift = Status::Pass;
        all.push_back(c);
    }
    {
        FixtureCase c = tube("d12k-ring-9-8", "DEC d12k (3.4)", "ring 9 -> 8 in the middle strip, the d12k pentagon at slot 4", {9, 9, 8, 8}, {{}, {4}, {}});
        c.runShift = true;
        c.expectShift = Status::Pass;
        all.push_back(c);
    }
    all.push_back(specialQuadCase());
    // negative controls: connectivity
    {
        FixtureCase c = tube("neg-cast-on-on-interior-ring", "cast-on (3.5)", "one bottom edge of strip 1 (an interior ring) marked CastOn", k8, none);
        for (int f = 0; f < c.mesh.mesh.faces(); ++f)
            if (c.mesh.strip[f] == 1 && c.mesh.slot[f] == 3) c.mesh.boundary[f][0] = EdgeBoundary::CastOn;
        c.negative = true;
        c.expectConnectivity = Status::Fail;
        all.push_back(c);
    }
    {
        FixtureCase c = tube("neg-bind-off-with-face-above", "bind-off (3.6)", "one top edge of strip 1 (a face above it) marked BindOff", k8, none);
        for (int f = 0; f < c.mesh.mesh.faces(); ++f)
            if (c.mesh.strip[f] == 1 && c.mesh.slot[f] == 3) c.mesh.boundary[f][2] = EdgeBoundary::BindOff;
        c.negative = true;
        c.expectConnectivity = Status::Fail;
        all.push_back(c);
    }
    {
        FixtureCase c = tube("neg-first-strip-dec", "cast-on (3.5) / lookup (6)", "a DEC in the first strip: the lookup refuses shaping-at-boundary-strip",
                             {9, 8, 8, 8}, {{4}, {}, {}});
        c.negative = true;
        c.expectConnectivity = Status::Unresolved;
        all.push_back(c);
    }
    {
        FixtureCase c = tube("neg-shaping-beside-path", "shift-path (3.7) / lookup (6)", "the kp at slot 0, beside the shift-path: refused", {8, 8, 9, 9},
                             {{}, {0}, {}});
        c.negative = true;
        c.expectConnectivity = Status::Unresolved;
        all.push_back(c);
    }
    {
        FixtureCase c = tube("neg-kp-roles-out-of-order", "INC kp (3.3)", "the kp face stored as INC with roles (B, T, W, T, W): the roles do not support it",
                             {8, 8, 9, 9}, {{}, {4}, {}});
        for (int f = 0; f < c.mesh.mesh.faces(); ++f)
            if (c.mesh.kind[f] == StitchFaceKind::Increase)
                c.mesh.role[f] = {StitchEdgeRole::Bottom, StitchEdgeRole::Top, StitchEdgeRole::Wale, StitchEdgeRole::Top, StitchEdgeRole::Wale};
        c.negative = true;
        c.expectConnectivity = Status::Fail;
        all.push_back(c);
    }
    {
        FixtureCase c = tube("neg-face-flipped", "interfaces (7)", "one K face of strip 1 with B and T exchanged: both-top and both-bottom edges, yarn against its row",
                             k8, none);
        for (int f = 0; f < c.mesh.mesh.faces(); ++f)
            if (c.mesh.strip[f] == 1 && c.mesh.slot[f] == 3) {
                c.mesh.role[f][0] = StitchEdgeRole::Top;
                c.mesh.role[f][2] = StitchEdgeRole::Bottom;
            }
        c.negative = true;
        c.expectConnectivity = Status::Unresolved;
        all.push_back(c);
    }
    {
        FixtureCase c = tube("neg-missing-face", "count invariant (7)", "one face of the d12k fixture's last strip removed: the strip's count no longer meets its rings",
                             {9, 9, 8, 8}, {{}, {4}, {}});
        TypedStitchMesh& s = c.mesh;
        const int drop = s.mesh.faces() - 3;
        TypedStitchMesh t = s;
        t.mesh.faceStart = {0};
        t.mesh.corner.clear();
        t.kind.clear();
        t.role.clear();
        t.boundary.clear();
        t.reason.clear();
        t.strip.clear();
        t.slot.clear();
        for (int f = 0; f < s.mesh.faces(); ++f) {
            if (f == drop) continue;
            std::vector<int> corners;
            for (int k = 0; k < s.mesh.sides(f); ++k) corners.push_back(s.mesh.vertexAt(f, k));
            t.mesh.addFace(corners);
            t.kind.push_back(s.kind[f]);
            t.role.push_back(s.role[f]);
            t.boundary.push_back(s.boundary[f]);
            t.reason.push_back(s.reason[f]);
            t.strip.push_back(s.strip[f]);
            t.slot.push_back(s.slot[f]);
        }
        c.mesh = t;
        c.negative = true;
        c.expectConnectivity = Status::Fail;
        all.push_back(c);
    }
    // negative controls: the virtual shift
    {
        FixtureCase c = tube("neg-shift-not-applied", "shift-path (3.7)", "the stored closed rings read as the virtual mesh: three closed strands, no helix", k8, none);
        c.runShift = true;
        c.reading = ShiftReading::None;
        c.negative = true;
        c.expectShift = Status::Fail;
        all.push_back(c);
    }
    {
        FixtureCase c = tube("neg-shift-read-upward", "shift-path (3.7)", "slot 0 read one ring higher: the helix would start at the bind-off", k8, none);
        c.runShift = true;
        c.reading = ShiftReading::PathSideUp;
        c.negative = true;
        c.expectShift = Status::Fail;
        all.push_back(c);
    }
    {
        FixtureCase c = tube("neg-shift-path-skips-ring", "shift-path (3.7)", "the path has no edge in strip 1", k8, none);
        c.mesh.shiftPath.erase(c.mesh.shiftPath.begin() + 1);
        c.runShift = true;
        c.negative = true;
        c.expectShift = Status::Fail;
        all.push_back(c);
    }
    {
        FixtureCase c = tube("neg-shift-path-two-edges-in-ring", "shift-path (3.7)", "the path has a second wale edge in strip 1 (vertex 4)", k8, none);
        c.mesh.shiftPath.push_back(ShiftPathEdge{1, c.mesh.rings[1][4], c.mesh.rings[2][4]});
        c.runShift = true;
        c.negative = true;
        c.expectShift = Status::Fail;
        all.push_back(c);
    }
    // negative control: the special quad
    {
        FixtureCase c = specialQuadCase();
        c.name = "neg-special-quad-wale-edges-apart";
        c.purpose = "one special quad stored as such with roles (B, W, T, W): its wale edges no longer meet, the roles do not support the kind";
        for (int f = 0; f < c.mesh.mesh.faces(); ++f)
            if (c.mesh.kind[f] == StitchFaceKind::SpecialQuad) {
                // rotate so the course edges alternate with the wale edges
                std::vector<StitchEdgeRole>& r = c.mesh.role[f];
                std::vector<StitchEdgeRole> alt(4, StitchEdgeRole::Wale);
                alt[0] = StitchEdgeRole::Bottom;
                alt[2] = StitchEdgeRole::Top;
                r = alt;
                break;
            }
        c.negative = true;
        c.expectConnectivity = Status::Fail;
        all.push_back(c);
    }
    return all;
}

} // namespace knit
