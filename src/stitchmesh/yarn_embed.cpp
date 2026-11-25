// The first yarn-generation experiment (yarn_embed.h).
#include "stitchmesh/yarn_embed.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <map>
#include <memory>
#include <set>
#include <unordered_map>

#include "model/mesh_query.h"

namespace knit {

const char* toString(KnotPlacement p) { return p == KnotPlacement::ReferenceProjection ? "reference" : "face"; }

const char* toString(YarnPieceKind k) {
    switch (k) {
    case YarnPieceKind::Head: return "head";
    case YarnPieceKind::Row: return "row";
    case YarnPieceKind::Cap: return "cap";
    case YarnPieceKind::Binding: return "binding";
    case YarnPieceKind::PoleJoin: return "pole_join";
    }
    return "row";
}

bool YarnPort::operator<(const YarnPort& o) const {
    if (kind != o.kind) return kind < o.kind;
    if (a != o.a) return a < o.a;
    if (b != o.b) return b < o.b;
    return index < o.index;
}

namespace {

constexpr double kYarnPi = 3.14159265358979323846;

// ------------------------------------------------------------------------------------------------ canonical templates
// Canonical face coordinates (u, v, z): u along the yarn's traversal (W- to W+), v from the bottom course edge(s) to the
// top, z the normal offset in yarn diameters. Every template's crossings are checked (the canonical
// prototype; and in 3D on every embedded face by checkYarn).
enum class EdgeType : char { B, T, WP, WM };
struct PortSpec {
    int edge = 0, index = 0;
};
struct PieceSpec {
    YarnPieceKind kind = YarnPieceKind::Row;
    int index = 0;
    PortSpec start, end;
    std::vector<std::array<double, 3>> knots;   // interior controls (u, v, z)
};
struct Template {
    std::string name;
    std::vector<std::array<double, 2>> domain;   // canonical corners, in the face's corner order (rotated)
    std::vector<EdgeType> edges;                 // edge k: corner k to corner k + 1
    std::vector<PieceSpec> pieces;
    bool pole = false;                           // wale ports at poleWalePort
    // required crossings: per head, along the row yarn and along the head; for two heads, the later over the earlier
    std::vector<std::string> headAlongNewer, headAlongOlder;
    std::string laterOverEarlier;
};

Template templateK() {
    Template t;
    t.name = "K";
    t.domain = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    t.edges = {EdgeType::B, EdgeType::WP, EdgeType::T, EdgeType::WM};
    t.pieces.push_back({YarnPieceKind::Head, 0, {0, 0}, {0, 1},
                        {{0.25, 0.18, 0.25}, {0.205, 0.40, 0.5}, {0.22, 0.63, 0}, {0.36, 0.78, -0.5}, {0.5, 0.82, -0.5}, {0.64, 0.78, -0.5}, {0.78, 0.63, 0},
                         {0.795, 0.40, 0.5}, {0.75, 0.18, 0.25}}});
    t.pieces.push_back({YarnPieceKind::Row, 0, {3, 0}, {2, 0}, {{0.12, 0.37, -0.35}, {0.23, 0.42, -0.5}, {0.34, 0.56, 0}, {0.37, 0.77, 0.5}, {0.32, 0.92, 0.1}}});
    t.pieces.push_back({YarnPieceKind::Row, 1, {2, 1}, {1, 0}, {{0.68, 0.92, 0.1}, {0.63, 0.77, 0.5}, {0.66, 0.56, 0}, {0.77, 0.42, -0.5}, {0.88, 0.37, -0.35}}});
    t.headAlongNewer = {"BFFB"};   // expected_signature.h kNewerAlongNewer
    t.headAlongOlder = {"BFFB"};   // kNewerAlongOlder
    return t;
}

// The calibrated K (the flat calibration's "calibrated" set): depth
// +-0.65 D (crossing gap >= 1.2 D > D + c); head sides at u 0.5 -/+ 0.385 and legs at 0.5 -/+ 0.12 over the depth flips, so a
// leg and the head side it swaps depth with are >= D + c apart in plane (production's proportions: 0.125 / 0.395 w); the
// leg's flip between its two crossings stretched over v 0.35 - 0.9. Used with the calibrated ports (YarnSettings).
Template templateKCalibrated() {
    Template t = templateK();
    t.name = "K";
    const std::vector<std::array<double, 3>> head = {{0.21, 0.15, 0.6}, {0.135, 0.34, 0.65}, {0.115, 0.56, 0.15}, {0.16, 0.77, -0.45}, {0.30, 0.90, -0.65}, {0.5, 0.94, -0.65},
                                                    {0.70, 0.90, -0.65}, {0.84, 0.77, -0.45}, {0.885, 0.56, 0.15}, {0.865, 0.34, 0.65}, {0.79, 0.15, 0.6}};
    const std::vector<std::array<double, 3>> row0 = {{0.13, 0.35, -0.55}, {0.26, 0.38, -0.6}, {0.355, 0.50, -0.3}, {0.38, 0.66, 0.2}, {0.375, 0.80, 0.6}, {0.345, 0.91, 0.6}};
    std::vector<std::array<double, 3>> row1;
    for (auto it = row0.rbegin(); it != row0.rend(); ++it) row1.push_back({1 - (*it)[0], (*it)[1], (*it)[2]});
    t.pieces[0].knots = head;
    t.pieces[1].knots = row0;
    t.pieces[2].knots = row1;
    return t;
}

// kp (W18 section 8): one loop of the row below, a knit slot (T1) and a purl slot (T2). Canonical trapezoid, B 1 wide, T 2.
Template templateKp() {
    Template t;
    t.name = "kp";
    t.domain = {{0.5, 0}, {1.5, 0}, {2, 1}, {1, 1}, {0, 1}};
    t.edges = {EdgeType::B, EdgeType::WP, EdgeType::T, EdgeType::T, EdgeType::WM};   // edge 2 = T2 (purl), edge 3 = T1 (knit)
    t.pieces.push_back({YarnPieceKind::Head, 0, {0, 0}, {0, 1},
                        {{0.7, 0.18, 0.25}, {0.50, 0.40, 0.5}, {0.42, 0.62, 0}, {0.42, 0.76, -0.5}, {0.62, 0.84, -0.5}, {1.0, 0.87, 0}, {1.38, 0.84, 0.5},
                         {1.58, 0.76, 0.5}, {1.58, 0.62, 0}, {1.50, 0.40, -0.5}, {1.3, 0.18, -0.25}}});
    t.pieces.push_back({YarnPieceKind::Row, 0, {4, 0}, {3, 0}, {{0.42, 0.38, -0.35}, {0.52, 0.44, -0.5}, {0.56, 0.58, 0}, {0.44, 0.82, 0.5}, {0.32, 0.93, 0.1}}});
    t.pieces.push_back({YarnPieceKind::Row, 1, {3, 1}, {2, 0},
                        {{0.68, 0.93, 0.1}, {0.64, 0.83, 0.5}, {0.78, 0.62, 0}, {1.0, 0.58, 0}, {1.22, 0.62, 0}, {1.36, 0.83, -0.5}, {1.32, 0.93, -0.1}}});
    t.pieces.push_back({YarnPieceKind::Row, 2, {2, 1}, {1, 0}, {{1.68, 0.93, -0.1}, {1.58, 0.82, -0.5}, {1.44, 0.58, 0}, {1.48, 0.44, 0.5}, {1.58, 0.38, 0.35}}});
    t.headAlongNewer = {"BFFBBF"};   // knit slot B,F,F then purl slot B,B,F (dbf5f41 kIncSlot0 + kIncSlot1Purl)
    t.headAlongOlder = {"BFFBBF"};
    return t;
}

// d12k (W18 section 8, k2tog): two loops of the row below, the later (B2) in front, one stitch drawn through both.
Template templateD12k() {
    Template t;
    t.name = "d12k";
    t.domain = {{0, 0}, {1, 0}, {2, 0}, {1.5, 1}, {0.5, 1}};
    t.edges = {EdgeType::B, EdgeType::B, EdgeType::WP, EdgeType::T, EdgeType::WM};
    t.pieces.push_back({YarnPieceKind::Head, 0, {0, 0}, {0, 1},
                        {{0.28, 0.16, 0.25}, {0.36, 0.38, 0.5}, {0.52, 0.60, 0}, {0.74, 0.72, -0.7}, {1.0, 0.76, -0.7}, {1.26, 0.72, -0.7}, {1.36, 0.58, -0.3},
                         {1.30, 0.44, 0.5}, {1.06, 0.28, -0.5}, {0.78, 0.14, -0.25}}});
    t.pieces.push_back({YarnPieceKind::Head, 1, {1, 0}, {1, 1},
                        {{1.22, 0.14, 0.25}, {0.94, 0.28, 0.5}, {0.70, 0.44, 0.5}, {0.64, 0.58, -0.3}, {0.74, 0.80, -0.2}, {1.0, 0.84, -0.2}, {1.26, 0.80, -0.2},
                         {1.48, 0.60, 0}, {1.64, 0.38, 0.5}, {1.72, 0.16, 0.25}}});
    t.pieces.push_back({YarnPieceKind::Row, 0, {4, 0}, {3, 0},
                        {{0.26, 0.37, -0.35}, {0.40, 0.42, -0.5}, {0.62, 0.50, -0.5}, {0.72, 0.62, 0}, {0.76, 0.80, 0.6}, {0.79, 0.93, 0.1}}});
    t.pieces.push_back({YarnPieceKind::Row, 1, {3, 1}, {2, 0},
                        {{1.21, 0.93, 0.1}, {1.24, 0.80, 0.6}, {1.28, 0.62, 0}, {1.38, 0.50, -0.5}, {1.60, 0.42, -0.5}, {1.74, 0.37, -0.35}}});
    t.headAlongNewer = {"BFFB", "BFFB"};   // a plain pair with each parent (expected_signature.h DEC)
    t.headAlongOlder = {"BFFB", "BFFB"};
    t.laterOverEarlier = "FF";            // the later parent in front where the two heads cross
    return t;
}

// The calibrated kp and d12k: the same principles where their constructions allow. Depth x 1.3 (crossing gaps from 1.0 to 1.3 D);
// knots next to a calibrated port at its side of the surface (row ends: behind at W-, in front at a knit slot's top). kp's
// head, both B edges of kp and d12k, kp's purl slot T2 and kp's W+ keep the frozen ports (their yarn changes depth next to
// the port; see the edge rule in embedYarn). d12k's two parent heads, 0.2 D apart in the baseline where they overlap: the
// earlier head's top lower and deeper (v 0.68, -1.1 D), the later one's higher and shallower (v 0.86, -0.3 D), so they are about
// D apart there while the legs still cross both in front.
Template scaledDepth(Template t, double k) {
    for (PieceSpec& p : t.pieces)
        for (auto& kn : p.knots) kn[2] *= k;
    return t;
}

Template templateKpCalibrated() {
    Template t = scaledDepth(templateKp(), 1.3);
    t.pieces[1].knots.back()[2] = 0.6;    // row 0 into T1 (knit slot, calibrated port in front)
    t.pieces[2].knots.front()[2] = 0.6;   // row 1 out of T1
    return t;
}

Template templateD12kCalibrated() {
    Template t = scaledDepth(templateD12k(), 1.3);
    // earlier head (B1 below, behind): top lower and deeper
    t.pieces[0].knots = {{0.28, 0.16, 0.33}, {0.36, 0.38, 0.65}, {0.50, 0.56, 0}, {0.72, 0.66, -1.1}, {1.0, 0.68, -1.1}, {1.26, 0.66, -1.1}, {1.36, 0.54, 0.15},
                         {1.30, 0.44, 0.65}, {1.06, 0.28, -0.65}, {0.78, 0.14, -0.33}};
    // later head (in front of the earlier one where they cross): top higher and shallower
    t.pieces[1].knots = {{1.22, 0.14, 0.33}, {0.94, 0.28, 0.65}, {0.70, 0.44, 0.65}, {0.64, 0.60, -0.3}, {0.74, 0.82, -0.3}, {1.0, 0.86, -0.3}, {1.26, 0.82, -0.3},
                         {1.48, 0.62, 0}, {1.64, 0.38, 0.65}, {1.72, 0.16, 0.33}};
    t.pieces[2].knots.back()[2] = 0.6;    // row 0 into T (calibrated port in front)
    t.pieces[3].knots.front()[2] = 0.6;   // row 1 out of T
    return t;
}

// Gathered pole start: a triangle (pole, ring 1 j + 1, ring 1 j), roles (W+, T, W-); no head (the
// first round's cast-on loops are gathered at the pole). Proposed.
Template templatePole() {
    Template t;
    t.name = "pole_start";
    t.domain = {{0.5, 0}, {1, 1}, {0, 1}};
    t.edges = {EdgeType::WP, EdgeType::T, EdgeType::WM};
    t.pole = true;
    t.pieces.push_back({YarnPieceKind::Row, 0, {2, 0}, {1, 0}, {{0.27, 0.72, 0}, {0.29, 0.90, 0}}});
    t.pieces.push_back({YarnPieceKind::Row, 1, {1, 1}, {0, 0}, {{0.71, 0.90, 0}, {0.73, 0.72, 0}}});
    return t;
}

StitchEdgeRole roleOf(EdgeType e) { return e == EdgeType::B ? StitchEdgeRole::Bottom : e == EdgeType::T ? StitchEdgeRole::Top : StitchEdgeRole::Wale; }

// The rotation r with roles[(k + r) % n] matching the template's edge k for every k; -1 when none.
int rotationFor(const Template& t, const std::vector<StitchEdgeRole>& roles) {
    const int n = static_cast<int>(roles.size());
    if (n != static_cast<int>(t.edges.size())) return -1;
    for (int r = 0; r < n; ++r) {
        bool ok = true;
        for (int k = 0; k < n && ok; ++k) ok = roles[(k + r) % n] == roleOf(t.edges[k]);
        if (ok) return r;
    }
    return -1;
}

// Mean value coordinates of an interior point of a convex polygon (Floater 2003).
std::vector<double> meanValue(const std::vector<std::array<double, 2>>& poly, double x, double y) {
    const int n = static_cast<int>(poly.size());
    std::vector<double> r(n), w(n), ang(n);
    for (int i = 0; i < n; ++i) r[i] = std::hypot(poly[i][0] - x, poly[i][1] - y);
    for (int i = 0; i < n; ++i) {
        const int j = (i + 1) % n;
        const double ax = poly[i][0] - x, ay = poly[i][1] - y, bx = poly[j][0] - x, by = poly[j][1] - y;
        ang[i] = std::atan2(ax * by - ay * bx, ax * bx + ay * by);
    }
    double sum = 0;
    for (int i = 0; i < n; ++i) {
        w[i] = (std::tan(0.5 * ang[(i + n - 1) % n]) + std::tan(0.5 * ang[i])) / std::max(1e-12, r[i]);
        sum += w[i];
    }
    for (double& x2 : w) x2 /= sum;
    return w;
}

Vec3 unit(const Vec3& v) {
    const double l = norm(v);
    return l > 0 ? v / l : v;
}

struct Surface {
    const TriMesh* m = nullptr;
    std::unique_ptr<MeshQuery> q;
    // the closest point and the vertex normals interpolated there
    void project(const Vec3& x, Vec3& p, Vec3& n) const {
        const MeshQuery::Hit h = q->closest(x);
        if (h.tri < 0) {
            p = x;
            n = Vec3{0, 0, 1};
            return;
        }
        const auto& t = m->tri[h.tri];
        const TrianglePoint tp = closestOnTriangle(x, m->pos[t[0]], m->pos[t[1]], m->pos[t[2]]);
        p = tp.p;
        Vec3 nn{0, 0, 0};
        for (int k = 0; k < 3; ++k) nn = nn + m->vnormal[t[k]] * tp.bary[k];
        n = unit(nn);
        if (norm(nn) == 0) n = m->faceNormal(h.tri);
    }
};

// Cubic Hermite through the controls with the given tangents; samplesPerSpan per span, the last point once. spanScaled: the
// tangents are directions, scaled on each span by its chord (the same direction on both sides of a control: G1).
std::vector<Vec3> sampleHermite(const std::vector<Vec3>& P, const std::vector<Vec3>& T, int per, bool spanScaled) {
    std::vector<Vec3> out;
    for (size_t i = 0; i + 1 < P.size(); ++i) {
        const double L = norm(P[i + 1] - P[i]);
        const Vec3 t0 = spanScaled ? unit(T[i]) * L : T[i], t1 = spanScaled ? unit(T[i + 1]) * L : T[i + 1];
        for (int k = 0; k < per; ++k) {
            const double s = static_cast<double>(k) / per;
            const double h00 = 2 * s * s * s - 3 * s * s + 1, h10 = s * s * s - 2 * s * s + s, h01 = -2 * s * s * s + 3 * s * s, h11 = s * s * s - s * s;
            out.push_back(P[i] * h00 + t0 * h10 + P[i + 1] * h01 + t1 * h11);
        }
    }
    out.push_back(P.back());
    return out;
}

// Port tangents given at the ends (unit, or zero: take the chord); Catmull-Rom inside; the end tangents scaled by their chord.
void finishPiece(YarnPiece& p, Vec3 t0, Vec3 t1, int per, bool spanScaled) {   // by value: callers pass the piece's own tangents
    const size_t n = p.controls.size();
    p.tangents.assign(n, Vec3{0, 0, 0});
    for (size_t i = 1; i + 1 < n; ++i) p.tangents[i] = (p.controls[i + 1] - p.controls[i - 1]) * 0.5;
    const double c0 = norm(p.controls[1] - p.controls[0]), c1 = norm(p.controls[n - 1] - p.controls[n - 2]);
    p.tangents[0] = norm(t0) > 0 ? unit(t0) * c0 : p.controls[1] - p.controls[0];
    p.tangents[n - 1] = norm(t1) > 0 ? unit(t1) * c1 : p.controls[n - 1] - p.controls[n - 2];
    p.samples = sampleHermite(p.controls, p.tangents, per, spanScaled);
}

struct PortGeom {
    YarnPort id;
    Vec3 p, t, n;   // point (on the reference, or z D off it), unit tangent in yarn direction, normal
    double z = 0;   // the port's offset along the normal [D] (calibrated ports; else 0)
};

struct FaceFrame {
    std::vector<int> c;                 // corners, rotated to the template
    std::vector<EdgeType> edges;
    std::vector<EdgeBoundary> boundary; // per rotated edge
    const Template* t = nullptr;
};

} // namespace

bool embedYarn(const TypedStitchMesh& s, const std::vector<std::vector<int>>& faces, const TriMesh& reference, double w, double h, const YarnSettings& st,
               YarnEmbedding& out, std::string* err) {
    out = YarnEmbedding{};
    const int nF = static_cast<int>(faces.size());
    if (nF == 0 || nF != s.mesh.faces() || s.role.size() != faces.size()) {
        if (err) *err = "embedYarn: no faces, or the virtual faces and the typed mesh disagree";
        return false;
    }
    TriMesh ref = reference;
    if (ref.vnormal.size() != ref.pos.size()) ref.computeNormals();
    Surface surf;
    surf.m = &ref;
    surf.q = std::make_unique<MeshQuery>(ref);
    // "mixed": the calibrated K and kp with the baseline d12k (the calibrated d12k loses a crossing pair while settling);
    // every edge of a d12k face then keeps the frozen ports its baseline knots were made for
    const bool mixed = st.templateSet == "mixed";
    const bool calibrated = st.templateSet == "calibrated" || mixed;
    if (!calibrated && st.templateSet != "baseline") {
        if (err) *err = "embedYarn: unknown template set " + st.templateSet;
        return false;
    }
    const Template tk = calibrated ? templateKCalibrated() : templateK(), tkp = calibrated ? templateKpCalibrated() : templateKp(),
                   td = calibrated && !mixed ? templateD12kCalibrated() : templateD12k(), tp = templatePole();
    const std::vector<Vec3>& X = s.mesh.pos;
    const double D = st.diameter;
    char b[300];
    int uniqueEnd = 0;
    out.faceTemplate.assign(nF, "");

    std::vector<FaceFrame> frames(nF);
    for (int f = 0; f < nF; ++f) {
        const Template* t = nullptr;
        switch (s.kind[f]) {
        case StitchFaceKind::Regular: t = &tk; break;
        case StitchFaceKind::Increase: t = &tkp; break;
        case StitchFaceKind::Decrease: t = &td; break;
        case StitchFaceKind::PoleStart: t = &tp; break;
        default: break;
        }
        const int r = t ? rotationFor(*t, s.role[f]) : -1;
        if (f == st.skipFace) continue;   // negative control: no yarn here
        if (!t || r < 0 || faces[f].size() != s.role[f].size()) {
            out.facesWithoutTemplate.push_back(f);
            continue;
        }
        FaceFrame& fr = frames[f];
        fr.t = t;
        const int n = static_cast<int>(faces[f].size());
        for (int k = 0; k < n; ++k) {
            fr.c.push_back(faces[f][(k + r) % n]);
            fr.edges.push_back(t->edges[k]);
            fr.boundary.push_back(s.boundary.size() == faces.size() && s.boundary[f].size() == static_cast<size_t>(n) ? s.boundary[f][(k + r) % n] : EdgeBoundary::Interior);
        }
        out.faceTemplate[f] = t->name;
    }

    // calibrated ports: an edge is calibrated unless one of its faces rules it out; both faces of an edge see the same answer.
    // Ruled out: any edge of a pole face; kp's purl slot T2 and its W+ (the purl's yarn meets the next knit slot there); the B
    // edges of kp and d12k (their heads' ends change depth next to the port).
    std::map<std::pair<int, int>, bool> allRegular;
    if (calibrated)
        for (int f = 0; f < nF; ++f) {
            const FaceFrame& fr = frames[f];
            if (!fr.t) continue;
            const int n = static_cast<int>(fr.c.size());
            for (int k = 0; k < n; ++k) {
                const int a = fr.c[k], bb = fr.c[(k + 1) % n];
                if (a == bb) continue;
                const std::pair<int, int> key{std::min(a, bb), std::max(a, bb)};
                const EdgeType et = fr.edges[k];
                bool ok = fr.t == &tk || fr.t == &tkp || (fr.t == &td && !mixed);
                if (fr.t == &tkp && et == EdgeType::WP) ok = false;
                if ((fr.t == &tkp || fr.t == &td) && et == EdgeType::B) ok = false;
                if (fr.t == &tkp && et == EdgeType::T && k == 2) ok = false;   // kp's edge 2 is T2, the purl slot (templateKp)
                const auto it = allRegular.find(key);
                allRegular[key] = (it == allRegular.end() ? true : it->second) && ok;
            }
        }
    auto calibratedEdge = [&](int a, int bb) {
        const auto it = allRegular.find({std::min(a, bb), std::max(a, bb)});
        return calibrated && it != allRegular.end() && it->second;
    };
    const double lean = st.calibratedLean * kYarnPi / 180.0;

    // The reference normal at every mesh vertex, and a point of an edge (a, b) at parameter prm with its normal, by the
    // placement rule: on the reference, or on the edge with the two corner normals blended. Both faces of an edge call
    // this with the same vertices in the same order, so they get the same point.
    const bool onReference = st.placement == KnotPlacement::ReferenceProjection;
    std::vector<Vec3> cornerNormal(X.size());
    for (size_t v = 0; v < X.size(); ++v) {
        Vec3 unused;
        surf.project(X[v], unused, cornerNormal[v]);
    }
    auto place = [&](const Vec3& point, int a, int bb, double prm, Vec3& p, Vec3& nn) {
        const Vec3 blended = cornerNormal[a] * (1 - prm) + cornerNormal[bb] * prm;
        if (onReference || norm(blended) == 0) {
            surf.project(point, p, nn);
            return;
        }
        p = point;
        nn = unit(blended);
    };
    // a port of a face: its point, shared tangent and normal
    auto port = [&](const FaceFrame& fr, const PortSpec& ps, PortGeom& g) {
        const int n = static_cast<int>(fr.c.size());
        const int a = fr.c[ps.edge], bb = fr.c[(ps.edge + 1) % n];
        const EdgeType et = fr.edges[ps.edge];
        if (a == bb) {   // a collapsed corner (the start step's W-)
            g.id = YarnPort{2, a, -1, 0};
            surf.project(X[a], g.p, g.n);
            g.t = Vec3{0, 0, 0};
            g.z = 0;
            return;
        }
        Vec3 point;
        if (et == EdgeType::B || et == EdgeType::T) {
            const int wSide = et == EdgeType::B ? a : bb, other = et == EdgeType::B ? bb : a;
            const double prm = ps.index == 0 ? 0.5 - st.courseOffset : 0.5 + st.courseOffset;
            point = X[wSide] + (X[other] - X[wSide]) * prm;
            g.id = YarnPort{0, wSide, other, ps.index};
            place(point, wSide, other, prm, g.p, g.n);
            const Vec3 up = unit(cross(g.n, X[other] - X[wSide]));
            g.t = ps.index == 0 ? up : up * -1.0;
            g.z = 0;
            if (calibratedEdge(wSide, other)) {
                // in front of the surface, the direction leaning towards the W- side (up at port 0, down at port 1)
                Vec3 e = X[other] - X[wSide];
                e = unit(e - g.n * dot(e, g.n));
                g.t = unit(g.t * std::cos(lean) - e * std::sin(lean));
                g.z = st.calibratedCourseDepth;
                g.p = g.p + g.n * (g.z * st.diameter);
            }
        } else {
            const int bottom = et == EdgeType::WP ? a : bb, top = et == EdgeType::WP ? bb : a;
            const double prm = fr.t->pole ? st.poleWalePort : st.walePort;
            point = X[bottom] + (X[top] - X[bottom]) * prm;
            g.id = YarnPort{1, bottom, top, 0};
            place(point, bottom, top, prm, g.p, g.n);
            g.t = unit(cross(X[top] - X[bottom], g.n));
            g.z = 0;
            if (calibratedEdge(bottom, top)) {   // behind the surface, as the sinker between two knit stitches
                g.z = st.calibratedWaleDepth;
                g.p = g.p + g.n * (g.z * st.diameter);
            }
        }
    };
    auto interior = [&](const FaceFrame& fr, const std::array<double, 3>& k, double& z) {
        const std::vector<double> lam = meanValue(fr.t->domain, k[0], k[1]);
        Vec3 foot{0, 0, 0}, blended{0, 0, 0};
        for (size_t i = 0; i < lam.size(); ++i) {
            foot = foot + X[fr.c[i]] * lam[i];
            blended = blended + cornerNormal[fr.c[i]] * lam[i];
        }
        Vec3 p = foot, nn = unit(blended);
        if (onReference || norm(blended) == 0) surf.project(foot, p, nn);
        z = k[2];
        return p + nn * (k[2] * D);
    };

    // the faces' pieces
    for (int f = 0; f < nF; ++f) {
        const FaceFrame& fr = frames[f];
        if (!fr.t) continue;
        for (const PieceSpec& ps : fr.t->pieces) {
            if (ps.kind == YarnPieceKind::Head && fr.boundary[ps.start.edge] == EdgeBoundary::CastOn) continue;   // cast-on: no loop below (Y12 4.3)
            YarnPiece piece;
            piece.face = f;
            piece.kind = ps.kind;
            piece.index = ps.index;
            piece.templ = fr.t->name;
            PortGeom g0, g1;
            port(fr, ps.start, g0);
            port(fr, ps.end, g1);
            piece.start = g0.id;
            piece.end = g1.id;
            piece.controls.push_back(g0.p);
            piece.depth.push_back(g0.z);
            for (const auto& k0 : ps.knots) {
                double z = 0;
                std::array<double, 3> k = k0;
                if (f == st.flipHeadsOfFace && ps.kind == YarnPieceKind::Head) k[2] = -k[2];   // negative control
                piece.controls.push_back(interior(fr, k, z));
                piece.depth.push_back(z);
            }
            piece.controls.push_back(g1.p);
            piece.depth.push_back(g1.z);
            finishPiece(piece, g0.t, g1.t, st.samplesPerSpan, st.spanTangents);
            out.pieces.push_back(piece);
        }
        // bind-off caps: the last row's loop closes above its top edge (proposed)
        for (size_t e = 0; e < fr.edges.size(); ++e) {
            if (fr.edges[e] != EdgeType::T || fr.boundary[e] != EdgeBoundary::BindOff) continue;
            PortGeom g0, g1;
            port(fr, PortSpec{static_cast<int>(e), 0}, g0);
            port(fr, PortSpec{static_cast<int>(e), 1}, g1);
            YarnPiece cap;
            cap.face = f;
            cap.kind = YarnPieceKind::Cap;
            cap.index = static_cast<int>(e);
            cap.templ = std::string(fr.t->name) + "+bind-off cap";
            cap.start = g0.id;
            cap.end = g1.id;
            // the loop widens above the edge (sides 0.35 w from its middle at half its height) so the binding thread can pass through it;
            // its first side behind, its second in front (the binding goes in front of the first and behind the second)
            const Vec3 up = unit(g0.t - g1.t);
            const Vec3 mid = (g0.p + g1.p) * 0.5;
            const Vec3 along = unit(g1.p - g0.p);
            Vec3 pm, nm;
            surf.project(mid, pm, nm);
            const double ch = st.capHeight * h;
            cap.controls = {g0.p, mid - along * (0.35 * w) + up * (0.5 * ch) - nm * (0.4 * D), mid + up * ch, mid + along * (0.35 * w) + up * (0.5 * ch) + nm * (0.4 * D), g1.p};
            cap.depth = {g0.z, -0.4, 0, 0.4, g1.z};
            finishPiece(cap, g0.t, g1.t, st.samplesPerSpan, st.spanTangents);
            out.pieces.push_back(cap);
            ++out.caps;
        }
    }

    // the gathered pole round opens at the path vertex and continues into strip 1 (proposed)
    int startFace = -1, startCorner = -1;
    for (int f = 0; f < nF; ++f) {
        const FaceFrame& fr = frames[f];
        if (!fr.t || fr.t->pole) continue;
        for (size_t e = 0; e < fr.edges.size(); ++e)
            if (fr.edges[e] == EdgeType::WM && fr.c[e] == fr.c[(e + 1) % fr.c.size()]) {
                startFace = f;
                startCorner = fr.c[e];
            }
    }
    if (st.joinPoleRound && startCorner >= 0) {
        int exitPiece = -1, enterPiece = -1, pole = -1;
        for (size_t i = 0; i < out.pieces.size(); ++i) {
            const YarnPiece& p = out.pieces[i];
            if (p.kind != YarnPieceKind::Row || p.templ != "pole_start") continue;
            if (p.end.kind == 1 && p.end.b == startCorner) {
                exitPiece = static_cast<int>(i);
                pole = p.end.a;
            }
            if (p.start.kind == 1 && p.start.b == startCorner) enterPiece = static_cast<int>(i);
        }
        if (exitPiece >= 0 && enterPiece >= 0) {
            YarnPiece& enter = out.pieces[enterPiece];
            // the round's yarn end (the gathered tail) lower on the same edge, so it does not coincide with the continuing yarn
            Vec3 tail, tn;
            surf.project(X[pole] + (X[startCorner] - X[pole]) * st.poleTailPort, tail, tn);
            enter.start = YarnPort{3, uniqueEnd++, -1, 0};
            enter.controls[0] = tail;
            finishPiece(enter, enter.tangents[0], enter.tangents.back(), st.samplesPerSpan, st.spanTangents);
            const YarnPiece& exit = out.pieces[exitPiece];
            YarnPiece join;
            join.face = exit.face;
            join.kind = YarnPieceKind::PoleJoin;
            join.templ = "pole_join";
            join.start = exit.end;
            join.end = YarnPort{2, startCorner, -1, 0};
            Vec3 c, cn;
            surf.project(X[startCorner], c, cn);
            const Vec3 x = exit.controls.back();
            join.controls = {x, (x + c) * 0.5, c};
            join.depth = {0, 0, 0};
            const Vec3 tc = unit(c - x);
            finishPiece(join, exit.tangents.back(), tc, st.samplesPerSpan, st.spanTangents);
            out.pieces.push_back(join);
            // the start step's row yarn leaves the corner in the join's direction
            for (YarnPiece& p : out.pieces)
                if (p.face == startFace && p.kind == YarnPieceKind::Row && p.start.kind == 2) finishPiece(p, tc, p.tangents.back(), st.samplesPerSpan, st.spanTangents);
            out.poleJoined = true;
            out.log.push_back("pole round joined at the path vertex: the gathered round's yarn end on its join edge, the round continuing into strip 1");
        }
    }

    // where pieces start and end
    std::map<YarnPort, std::vector<int>> starts, ends;
    for (size_t i = 0; i < out.pieces.size(); ++i) {
        starts[out.pieces[i].start].push_back(static_cast<int>(i));
        ends[out.pieces[i].end].push_back(static_cast<int>(i));
    }

    // the binding thread: from the helix's end through every cap once, in yarn order (proposed)
    if (st.bindingThread && out.caps > 0) {
        std::vector<int> openEnds;
        for (const auto& [p, list] : ends)
            if (p.kind == 1 && !starts.count(p)) openEnds.push_back(list[0]);
        if (openEnds.size() == 1) {
            const YarnPiece& last = out.pieces[openEnds[0]];
            std::vector<int> caps;
            for (size_t i = 0; i < out.pieces.size(); ++i)
                if (out.pieces[i].kind == YarnPieceKind::Cap) caps.push_back(static_cast<int>(i));
            auto order = [&](int i) {
                const int f = out.pieces[i].face;
                return std::make_pair(s.strip.empty() ? 0 : s.strip[f], s.slot.empty() ? f : s.slot[f]);
            };
            std::sort(caps.begin(), caps.end(), [&](int x, int y) { return order(x) < order(y); });
            YarnPiece bind;
            bind.face = -1;
            bind.kind = YarnPieceKind::Binding;
            bind.templ = "binding";
            bind.start = last.end;
            bind.end = YarnPort{3, uniqueEnd++, -1, 0};
            bind.controls.push_back(last.controls.back());
            bind.depth.push_back(0);
            // through each loop at half its height: in front of its first side, behind its second (knots inside the loop, so each side is
            // crossed once, between knots); between two loops, half way between their ports (on the surface) at the same height
            Vec3 along{0, 0, 0}, lastP1{0, 0, 0}, lastUp{0, 0, 0};
            const double lift = 0.5 * st.capHeight * h;
            for (size_t j = 0; j < caps.size(); ++j) {
                const YarnPiece& cap = out.pieces[caps[j]];
                const Vec3 p0 = cap.controls.front(), p1 = cap.controls.back();
                const Vec3 up = unit(cap.tangents.front() - cap.tangents.back());
                along = unit(p1 - p0);
                Vec3 pc, nc;
                surf.project((p0 + p1) * 0.5, pc, nc);
                const Vec3 mid = (p0 + p1) * 0.5 + up * lift;
                if (j == 0) bind.controls.push_back(p0 - along * (0.3 * w) + up * lift);
                else bind.controls.push_back((lastP1 + p0) * 0.5 + unit(lastUp + up) * lift);
                bind.depth.push_back(0);
                bind.controls.push_back(mid - along * (0.15 * w) + nc * (0.6 * D));
                bind.depth.push_back(0.6);
                bind.controls.push_back(mid + along * (0.15 * w) - nc * (0.6 * D));
                bind.depth.push_back(-0.6);
                lastP1 = p1;
                lastUp = up;
            }
            // the yarn's end: 0.6 w past the last loop, but never as far as the first loop, which on a closed bind-off ring
            // is the next one along (at most 0.4 of the way to its first port, so the end cannot thread it a second time)
            const double tail = std::min(0.6 * w, 0.4 * norm(out.pieces[caps.front()].controls.front() - lastP1));
            // one closing knot when the tail is short: two knots closer than 0.3 w bend the yarn's held end tighter than D / 2
            if (tail >= 0.6 * w) {
                bind.controls.push_back(lastP1 + along * (0.5 * tail) + lastUp * lift);
                bind.depth.push_back(0);
            }
            bind.controls.push_back(lastP1 + along * tail + lastUp * lift);
            bind.depth.push_back(0);
            finishPiece(bind, last.tangents.back(), Vec3{0, 0, 0}, st.samplesPerSpan, st.spanTangents);
            out.pieces.push_back(bind);
            starts[bind.start].push_back(static_cast<int>(out.pieces.size()) - 1);
            ends[bind.end].push_back(static_cast<int>(out.pieces.size()) - 1);
            std::snprintf(b, sizeof b, "binding thread through %zu caps, from the helix's end", caps.size());
            out.log.push_back(b);
        } else {
            std::snprintf(b, sizeof b, "no binding thread: %zu open row ends (one expected)", openEnds.size());
            out.log.push_back(b);
        }
    }

    // strands: follow each piece's end port to the piece that starts there
    const int nP = static_cast<int>(out.pieces.size());
    std::vector<int> next(nP, -1), prev(nP, -1);
    for (int i = 0; i < nP; ++i) {
        const auto it = starts.find(out.pieces[i].end);
        if (out.pieces[i].end.kind != 3 && it != starts.end() && it->second.size() == 1 && ends[out.pieces[i].end].size() == 1) {
            next[i] = it->second[0];
            prev[it->second[0]] = i;
        }
    }
    std::vector<char> used(nP, 0);
    auto addStrand = [&](int first, bool closed) {
        YarnStrand sd;
        sd.closed = closed;
        for (int i = first; i >= 0 && !used[i]; i = next[i]) {
            used[i] = 1;
            sd.pieces.push_back(i);
            const YarnPiece& p = out.pieces[i];
            for (size_t k = sd.samples.empty() ? 0 : 1; k < p.samples.size(); ++k) {
                sd.samples.push_back(p.samples[k]);
                sd.sampleFace.push_back(p.face);
                sd.samplePiece.push_back(i);
            }
        }
        sd.arc.assign(sd.samples.size(), 0);
        for (size_t k = 1; k < sd.samples.size(); ++k) sd.arc[k] = sd.arc[k - 1] + norm(sd.samples[k] - sd.samples[k - 1]);
        out.strands.push_back(sd);
    };
    for (int i = 0; i < nP; ++i)
        if (prev[i] < 0 && !used[i]) addStrand(i, false);
    for (int i = 0; i < nP; ++i)
        if (!used[i]) addStrand(i, true);
    std::snprintf(b, sizeof b, "%d pieces, %zu strands, %zu faces without a template, %d caps%s", nP, out.strands.size(), out.facesWithoutTemplate.size(), out.caps,
                  out.poleJoined ? ", pole round joined" : "");
    out.log.push_back(b);
    return true;
}

// ------------------------------------------------------------------------------------------------------------- checks
namespace {

struct Crossing {
    double alongA = 0, alongB = 0;
    char letter = 'B';   // A's depth relative to B
    double dz = 0;       // |height difference| [mm]
    double u = 0, v = 0; // where in the plane [mm]
};

// Projected crossings of two sampled curves in a plane (origin o, basis e1, e2, normal n).
std::vector<Crossing> planeCrossings(const std::vector<Vec3>& A, const std::vector<Vec3>& B, const Vec3& o, const Vec3& e1, const Vec3& e2, const Vec3& n) {
    std::vector<Crossing> out;
    auto pr = [&](const Vec3& x, double& u, double& v, double& z) {
        const Vec3 d = x - o;
        u = dot(d, e1);
        v = dot(d, e2);
        z = dot(d, n);
    };
    std::vector<std::array<double, 3>> a(A.size()), bb(B.size());
    for (size_t i = 0; i < A.size(); ++i) pr(A[i], a[i][0], a[i][1], a[i][2]);
    for (size_t i = 0; i < B.size(); ++i) pr(B[i], bb[i][0], bb[i][1], bb[i][2]);
    for (size_t i = 0; i + 1 < a.size(); ++i)
        for (size_t j = 0; j + 1 < bb.size(); ++j) {
            const double px = a[i][0], py = a[i][1], qx = a[i + 1][0], qy = a[i + 1][1];
            const double rx = bb[j][0], ry = bb[j][1], sx = bb[j + 1][0], sy = bb[j + 1][1];
            const double d = (qx - px) * (sy - ry) - (qy - py) * (sx - rx);
            if (std::fabs(d) < 1e-14) continue;
            const double t = ((rx - px) * (sy - ry) - (ry - py) * (sx - rx)) / d;
            const double u = ((rx - px) * (qy - py) - (ry - py) * (qx - px)) / d;
            if (t < 0 || t >= 1 || u < 0 || u >= 1) continue;
            const double za = a[i][2] + t * (a[i + 1][2] - a[i][2]), zb = bb[j][2] + u * (bb[j + 1][2] - bb[j][2]);
            out.push_back(Crossing{i + t, j + u, za > zb ? 'F' : 'B', std::fabs(za - zb), px + t * (qx - px), py + t * (qy - py)});
        }
    return out;
}

// How far a point of the plane lies outside a polygon given by its corners in order [same units]: 0 inside (even-odd
// rule), else the distance to the nearest edge.
double outsidePolygon(const std::vector<std::array<double, 2>>& polygon, double x, double y) {
    const size_t n = polygon.size();
    bool inside = false;
    double nearest = 1e300;
    for (size_t i = 0; i < n; ++i) {
        const double ax = polygon[i][0], ay = polygon[i][1], bx = polygon[(i + 1) % n][0], by = polygon[(i + 1) % n][1];
        if ((ay > y) != (by > y) && x < ax + (y - ay) * (bx - ax) / (by - ay)) inside = !inside;
        const double ex = bx - ax, ey = by - ay, l2 = ex * ex + ey * ey;
        const double t = l2 > 0 ? std::max(0.0, std::min(1.0, ((x - ax) * ex + (y - ay) * ey) / l2)) : 0.0;
        nearest = std::min(nearest, std::hypot(x - (ax + t * ex), y - (ay + t * ey)));
    }
    return inside ? 0.0 : nearest;
}

std::string lettersAlongA(std::vector<Crossing> c) {
    std::sort(c.begin(), c.end(), [](const Crossing& x, const Crossing& y) { return x.alongA < y.alongA; });
    std::string s;
    for (const Crossing& x : c) s += x.letter;
    return s;
}
std::string lettersAlongB(std::vector<Crossing> c) {
    std::sort(c.begin(), c.end(), [](const Crossing& x, const Crossing& y) { return x.alongB < y.alongB; });
    std::string s;
    for (const Crossing& x : c) s += x.letter;
    return s;
}

// The normal of the plane that fits points best (least squares): the eigenvector of their covariance with the smallest
// eigenvalue, by Jacobi rotations of the 3 x 3 symmetric matrix.
Vec3 bestFitNormal(const std::vector<Vec3>& pts, Vec3& centroid) {
    centroid = Vec3{0, 0, 0};
    for (const Vec3& p : pts) centroid = centroid + p;
    centroid = centroid / static_cast<double>(std::max<size_t>(1, pts.size()));
    double A[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (const Vec3& p : pts) {
        const double d[3] = {p.x - centroid.x, p.y - centroid.y, p.z - centroid.z};
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) A[i][j] += d[i] * d[j];
    }
    double V[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    for (int sweep = 0; sweep < 50; ++sweep) {
        const double off = A[0][1] * A[0][1] + A[0][2] * A[0][2] + A[1][2] * A[1][2];
        if (off < 1e-30) break;
        for (int p = 0; p < 2; ++p)
            for (int q = p + 1; q < 3; ++q) {
                if (std::fabs(A[p][q]) < 1e-300) continue;
                const double theta = 0.5 * (A[q][q] - A[p][p]) / A[p][q];
                const double t = (theta >= 0 ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta * theta + 1));
                const double cs = 1 / std::sqrt(t * t + 1), sn = t * cs;
                for (int k = 0; k < 3; ++k) {   // A <- J^T A J
                    const double akp = A[k][p], akq = A[k][q];
                    A[k][p] = cs * akp - sn * akq;
                    A[k][q] = sn * akp + cs * akq;
                }
                for (int k = 0; k < 3; ++k) {
                    const double apk = A[p][k], aqk = A[q][k];
                    A[p][k] = cs * apk - sn * aqk;
                    A[q][k] = sn * apk + cs * aqk;
                }
                for (int k = 0; k < 3; ++k) {
                    const double vkp = V[k][p], vkq = V[k][q];
                    V[k][p] = cs * vkp - sn * vkq;
                    V[k][q] = sn * vkp + cs * vkq;
                }
            }
    }
    int m = 0;
    for (int i = 1; i < 3; ++i)
        if (A[i][i] < A[m][m]) m = i;
    return unit(Vec3{V[0][m], V[1][m], V[2][m]});
}

// How often a curve passes through a loop: the loop (closed by the straight segment from its last sample to its first) is
// spanned by a fan from its centroid; every curve segment through a fan triangle counts once (Moller-Trumbore).
int threadings(const std::vector<Vec3>& loop, const std::vector<Vec3>& curve) {
    if (loop.size() < 3) return 0;
    Vec3 c{0, 0, 0};
    for (const Vec3& p : loop) c = c + p;
    c = c / static_cast<double>(loop.size());
    int count = 0;
    for (size_t i = 0; i < loop.size(); ++i) {
        const Vec3 a = loop[i], b = loop[(i + 1) % loop.size()];
        const Vec3 e1 = a - c, e2 = b - c;
        for (size_t k = 0; k + 1 < curve.size(); ++k) {
            const Vec3 o = curve[k], d = curve[k + 1] - curve[k];
            const Vec3 pv = cross(d, e2);
            const double det = dot(e1, pv);
            if (std::fabs(det) < 1e-14) continue;
            const Vec3 tv = o - c;
            const double u = dot(tv, pv) / det;
            if (u < 0 || u > 1) continue;
            const Vec3 qv = cross(tv, e1);
            const double v = dot(d, qv) / det;
            if (v < 0 || u + v > 1) continue;
            const double t = dot(e2, qv) / det;
            if (t >= 0 && t < 1) ++count;
        }
    }
    return count;
}

double angleDeg(const Vec3& a, const Vec3& b) {
    const double la = norm(a), lb = norm(b);
    if (la == 0 || lb == 0) return 0;
    return std::acos(std::max(-1.0, std::min(1.0, dot(a, b) / (la * lb)))) * 180.0 / kYarnPi;
}

} // namespace

YarnCheck checkYarn(const TypedStitchMesh& s, const std::vector<std::vector<int>>& faces, const TriMesh& reference, const YarnEmbedding& e, const YarnSettings& st) {
    YarnCheck c;
    const double D = st.diameter;
    TriMesh ref = reference;
    if (ref.vnormal.size() != ref.pos.size()) ref.computeNormals();
    const MeshQuery q(ref);
    const std::vector<Vec3>& X = s.mesh.pos;
    char b[400];
    c.pieces = static_cast<int>(e.pieces.size());
    c.facesWithoutTemplate = static_cast<int>(e.facesWithoutTemplate.size());
    for (int f : e.facesWithoutTemplate) {
        std::snprintf(b, sizeof b, "face %d (%s, strip %d slot %d): no yarn template for this construction", f, toString(s.kind[f]), s.strip.empty() ? -1 : s.strip[f],
                      s.slot.empty() ? -1 : s.slot[f]);
        c.defects.push_back(b);
    }
    // ports
    std::map<YarnPort, int> in, outn;
    for (const YarnPiece& p : e.pieces) {
        outn[p.start] += 1;
        in[p.end] += 1;
    }
    std::set<YarnPort> all;
    for (const auto& [p, n] : in) all.insert(p);
    for (const auto& [p, n] : outn) all.insert(p);
    for (const YarnPort& p : all) {
        const int a = in.count(p) ? in[p] : 0, o = outn.count(p) ? outn[p] : 0;
        if (a == 1 && o == 1) ++c.portsSharedOk;
        else if (a > 1 || o > 1) {
            ++c.portsOverused;
            std::snprintf(b, sizeof b, "port (kind %d, %d, %d, %d) used by %d pieces arriving and %d leaving", p.kind, p.a, p.b, p.index, a, o);
            c.defects.push_back(b);
        } else if (p.kind == 3 || (p.kind == 2 && a == 0)) ++c.yarnEnds;   // a declared yarn end, or the helix's start corner
        else {
            ++c.portsUnmatched;
            std::snprintf(b, sizeof b, "port (kind %s, vertices %d-%d, index %d) has %s only: a yarn end where the construction has none", p.kind == 0 ? "course" : p.kind == 1 ? "wale" : "corner",
                          p.a, p.b, p.index, a ? "a piece arriving" : "a piece leaving");
            c.defects.push_back(b);
        }
    }
    // strands, continuity
    c.strands = static_cast<int>(e.strands.size());
    for (const YarnStrand& sd : e.strands) {
        (sd.closed ? c.closedStrands : c.openStrands) += 1;
        if (!sd.arc.empty()) c.totalLength += sd.arc.back();
        for (size_t k = 0; k + 1 < sd.pieces.size(); ++k) {
            const YarnPiece& a = e.pieces[sd.pieces[k]];
            const YarnPiece& bb = e.pieces[sd.pieces[k + 1]];
            c.maxPortGap = std::max(c.maxPortGap, norm(a.controls.back() - bb.controls.front()));
            const double jump = angleDeg(a.tangents.back(), bb.tangents.front());
            c.maxTangentJump = std::max(c.maxTangentJump, jump);
            if (jump > 1.0) {
                std::snprintf(b, sizeof b, "direction jumps %.1f deg from a %s piece (face %d) to a %s piece (face %d)", jump, toString(a.kind), a.face, toString(bb.kind), bb.face);
                c.defects.push_back(b);
            }
        }
    }
    for (const YarnPiece& p : e.pieces)
        for (const YarnPort* port : {&p.start, &p.end})
            if (port->kind == 0 || port->kind == 1) {
                const Vec3 x = port == &p.start ? p.controls.front() : p.controls.back();
                c.maxPortOffSurface = std::max(c.maxPortOffSurface, norm(q.closest(x).p - x));
            }
    // per face: crossings in its own plane
    std::vector<std::vector<int>> piecesOf(faces.size());
    int binding = -1;
    for (size_t i = 0; i < e.pieces.size(); ++i) {
        if (e.pieces[i].face >= 0 && e.pieces[i].face < static_cast<int>(faces.size())) piecesOf[e.pieces[i].face].push_back(static_cast<int>(i));
        if (e.pieces[i].kind == YarnPieceKind::Binding) binding = static_cast<int>(i);
    }
    std::vector<int> faceIndex(faces.size(), -1);
    for (size_t f = 0; f < faces.size(); ++f) {
        if (e.faceTemplate[f].empty()) continue;
        YarnFaceCheck fc;
        fc.face = static_cast<int>(f);
        fc.templ = e.faceTemplate[f];
        // the frame: centroid of the distinct corners on the reference, the reference normal there
        std::vector<int> distinct;
        for (int v : faces[f])
            if (std::find(distinct.begin(), distinct.end(), v) == distinct.end()) distinct.push_back(v);
        Vec3 o{0, 0, 0};
        for (int v : distinct) o = o + X[v];
        o = o / static_cast<double>(distinct.size());
        const MeshQuery::Hit h = q.closest(o);
        Vec3 n = h.tri >= 0 ? ref.faceNormal(h.tri) : Vec3{0, 0, 1};
        if (h.tri >= 0) {
            const auto& t = ref.tri[h.tri];
            const TrianglePoint tp = closestOnTriangle(o, ref.pos[t[0]], ref.pos[t[1]], ref.pos[t[2]]);
            Vec3 nn{0, 0, 0};
            for (int k = 0; k < 3; ++k) nn = nn + ref.vnormal[t[k]] * tp.bary[k];
            if (norm(nn) > 0) n = unit(nn);
            o = tp.p;
        }
        Vec3 e1 = X[distinct[1 % distinct.size()]] - X[distinct[0]];
        e1 = unit(e1 - n * dot(e1, n));
        const Vec3 e2 = cross(n, e1);
        std::vector<int> heads(2, -1), caps;
        std::vector<Vec3> row;
        std::vector<std::pair<int, int>> rows;   // (index, piece)
        for (int i : piecesOf[f]) {
            const YarnPiece& p = e.pieces[i];
            if (p.kind == YarnPieceKind::Head) heads[p.index] = i;
            else if (p.kind == YarnPieceKind::Row) rows.push_back({p.index, i});
            else if (p.kind == YarnPieceKind::Cap) caps.push_back(i);
        }
        std::sort(rows.begin(), rows.end());
        for (const auto& [k, i] : rows)
            for (size_t m = row.empty() ? 0 : 1; m < e.pieces[i].samples.size(); ++m) row.push_back(e.pieces[i].samples[m]);
        // the template's required crossings per head, judged in a frame (origin, basis, normal): heads against the row, the
        // decrease's later head in front of the earlier, no crossings the construction does not have
        const std::string tn = fc.templ;
        const std::vector<std::string> needNew = tn == "K" ? std::vector<std::string>{"BFFB"} : tn == "kp" ? std::vector<std::string>{"BFFBBF"}
                                               : tn == "d12k" ? std::vector<std::string>{"BFFB", "BFFB"} : std::vector<std::string>{};
        std::string expected;
        for (size_t k = 0; k < needNew.size(); ++k)
            if (heads[k] >= 0) expected += "head" + std::to_string(k) + ":" + needNew[k] + "/" + needNew[k] + " ";
        if (tn == "d12k" && heads[0] >= 0 && heads[1] >= 0) expected += "h1/h0:FF ";
        struct FrameResult {
            bool ok = true;
            std::string found;
            int interlocks = 0, unexpected = 0;
            double minDepthGap = 1e300;
        };
        // `record`: the head-row crossings one by one, for the diagnosis (YarnFaceCheck::crossings); null: not kept
        auto judge = [&](const Vec3& fo, const Vec3& f1, const Vec3& f2, const Vec3& fn, std::vector<YarnCrossing>* record) {
            FrameResult r;
            for (size_t k = 0; k < needNew.size(); ++k) {
                if (heads[k] < 0) continue;   // a cast-on bottom edge: no loop below
                // per row piece, in order along the row (the pieces are not joined: a straight segment between two of them, from
                // one top port to the next, is no yarn and must not cross anything)
                std::vector<Crossing> cr;
                double along = 0;
                const double headSpans = static_cast<double>(std::max<size_t>(2, e.pieces[heads[k]].samples.size()) - 1);
                for (const auto& [ri, pi] : rows) {
                    const double rowSpans = static_cast<double>(std::max<size_t>(2, e.pieces[pi].samples.size()) - 1);
                    for (Crossing x : planeCrossings(e.pieces[pi].samples, e.pieces[heads[k]].samples, fo, f1, f2, fn)) {
                        if (record) record->push_back(YarnCrossing{static_cast<int>(k), ri, x.alongA / rowSpans, x.alongB / headSpans, x.letter, x.dz / D, x.u, x.v});
                        x.alongA += along;
                        cr.push_back(x);
                    }
                    along += static_cast<double>(e.pieces[pi].samples.size());
                }
                const std::string alongNew = lettersAlongA(cr), alongOld = lettersAlongB(cr);
                r.found += "head" + std::to_string(k) + ":" + alongNew + "/" + alongOld + " ";
                if (alongNew == needNew[k] && alongOld == needNew[k]) ++r.interlocks;
                else r.ok = false;
                for (const Crossing& x : cr) r.minDepthGap = std::min(r.minDepthGap, x.dz / D);
            }
            if (tn == "d12k" && heads[0] >= 0 && heads[1] >= 0) {
                const std::vector<Crossing> cr = planeCrossings(e.pieces[heads[1]].samples, e.pieces[heads[0]].samples, fo, f1, f2, fn);
                const std::string l = lettersAlongA(cr);
                r.found += "h1/h0:" + l + " ";
                if (l != "FF") r.ok = false;
                for (const Crossing& x : cr) r.minDepthGap = std::min(r.minDepthGap, x.dz / D);
            }
            // crossings the construction does not have: row with itself (pieces apart), heads with caps
            for (size_t a = 0; a < rows.size(); ++a)
                for (size_t bb = a + 2; bb < rows.size(); ++bb)
                    r.unexpected += static_cast<int>(planeCrossings(e.pieces[rows[a].second].samples, e.pieces[rows[bb].second].samples, fo, f1, f2, fn).size());
            for (int cp : caps)
                for (int hd : heads)
                    if (hd >= 0) r.unexpected += static_cast<int>(planeCrossings(e.pieces[hd].samples, e.pieces[cp].samples, fo, f1, f2, fn).size());
            if (r.unexpected > 0) r.ok = false;
            return r;
        };
        // the frame of the construction: the reference at the face
        const FrameResult atReference = judge(o, e1, e2, n, &fc.crossings);
        std::string found = atReference.found;
        // diagnosis: the face's shape in that frame, and how far its yarn leaves the face's polygon
        {
            std::vector<std::array<double, 2>> polygon;
            for (int v : distinct) {
                const Vec3 d = X[v] - o;
                polygon.push_back({dot(d, e1), dot(d, e2)});
                fc.nonPlanar = std::max(fc.nonPlanar, std::fabs(dot(d, n)));
            }
            Vec3 bottomMid{0, 0, 0}, topMid{0, 0, 0}, bottomDir{0, 0, 0};
            int bottoms = 0, tops = 0;
            fc.waleMin = 1e300;
            const std::vector<int>& corners = faces[f];
            for (size_t k = 0; k < corners.size() && k < s.role[f].size(); ++k) {
                const Vec3 a = X[corners[k]], bb = X[corners[(k + 1) % corners.size()]];
                const double length = norm(bb - a);
                if (s.role[f][k] == StitchEdgeRole::Bottom) {
                    fc.bottomLength += length;
                    bottomMid = bottomMid + (a + bb) * 0.5;
                    bottomDir = bottomDir + (bb - a);
                    ++bottoms;
                } else if (s.role[f][k] == StitchEdgeRole::Top) {
                    fc.topLength += length;
                    topMid = topMid + (a + bb) * 0.5;
                    ++tops;
                } else if (length > 0) {   // a wale edge (a collapsed one, the helix's start step, has no length)
                    fc.waleMin = std::min(fc.waleMin, length);
                    fc.waleMax = std::max(fc.waleMax, length);
                }
            }
            if (fc.waleMin > 1e299) fc.waleMin = 0;
            if (bottoms > 0 && tops > 0) {
                Vec3 up = topMid / static_cast<double>(tops) - bottomMid / static_cast<double>(bottoms);
                up = up - n * dot(up, n);
                bottomDir = bottomDir - n * dot(bottomDir, n);
                if (norm(up) > 0 && norm(bottomDir) > 0) fc.shear = std::fabs(90.0 - angleDeg(up, bottomDir));
            }
            auto outside = [&](const std::vector<Vec3>& samples) {
                double worst = 0;
                for (const Vec3& p : samples) worst = std::max(worst, outsidePolygon(polygon, dot(p - o, e1), dot(p - o, e2)));
                return worst;
            };
            for (int hd : heads)
                if (hd >= 0) fc.headOutside = std::max(fc.headOutside, outside(e.pieces[hd].samples));
            fc.rowOutside = outside(row);
        }
        fc.crossingsOk = atReference.ok;
        fc.unexpectedCrossings = atReference.unexpected;
        fc.minDepthGap = atReference.minDepthGap;
        c.interlocksExpected += static_cast<int>(std::count_if(heads.begin(), heads.end(), [](int x) { return x >= 0; }));
        c.interlocksVerified += atReference.interlocks;
        // the frame of the yarn as it is: the plane fitting the face's head and row samples best, its normal on the reference's side
        {
            std::vector<Vec3> pts = row;
            for (int hd : heads)
                if (hd >= 0) pts.insert(pts.end(), e.pieces[hd].samples.begin(), e.pieces[hd].samples.end());
            Vec3 oc;
            Vec3 nc = bestFitNormal(pts, oc);
            if (dot(nc, n) < 0) nc = nc * -1.0;
            Vec3 c1 = e1 - nc * dot(e1, nc);
            c1 = norm(c1) > 0 ? unit(c1) : e1;
            const FrameResult own = judge(oc, c1, cross(nc, c1), nc, nullptr);
            fc.ownPlaneOk = own.ok;
            fc.foundOwnPlane = own.found;
            fc.ownPlaneTilt = std::acos(std::max(-1.0, std::min(1.0, dot(nc, n)))) * 180.0 / kYarnPi;
            c.interlocksOwnPlane += own.interlocks;
        }
        // the binding thread through each cap: exactly once through the loop the cap and its edge span (no frame)
        for (size_t k = 0; k < caps.size(); ++k) {
            if (binding < 0) continue;
            const int n2 = threadings(e.pieces[caps[k]].samples, e.pieces[binding].samples);
            expected += "cap:through1 ";
            found += "cap:through" + std::to_string(n2) + " ";
            if (n2 != 1) {
                fc.crossingsOk = false;
                fc.ownPlaneOk = false;
            }
        }
        fc.expected = expected;
        fc.found = found;
        ++c.facesChecked;
        c.facesCrossingsOk += fc.crossingsOk;
        c.facesOwnPlaneOk += fc.ownPlaneOk;
        c.unexpectedCrossings += fc.unexpectedCrossings;
        if (fc.minDepthGap < 1e299) c.minDepthGap = std::min(c.minDepthGap, fc.minDepthGap);
        if (!fc.crossingsOk) {
            std::snprintf(b, sizeof b, "face %zu (%s, strip %d slot %d): crossings %s expected %s (%d unexpected)", f, tn.c_str(), s.strip.empty() ? -1 : s.strip[f],
                          s.slot.empty() ? -1 : s.slot[f], found.c_str(), expected.c_str(), fc.unexpectedCrossings);
            c.defects.push_back(b);
        }
        faceIndex[f] = static_cast<int>(c.faces.size());
        c.faces.push_back(fc);
    }
    // clearance between yarn more than 3 D apart along the yarn, and bends
    struct Sample {
        int strand;
        double arc;
        Vec3 p;
        int face;
        int piece;
    };
    std::vector<Sample> all2;
    for (size_t si = 0; si < e.strands.size(); ++si)
        for (size_t k = 0; k < e.strands[si].samples.size(); ++k)
            all2.push_back({static_cast<int>(si), e.strands[si].arc[k], e.strands[si].samples[k], e.strands[si].sampleFace[k], e.strands[si].samplePiece[k]});
    std::map<std::pair<int, int>, YarnCheck::Contact> pairs;
    c.samples = static_cast<int>(all2.size());
    std::unordered_map<long long, std::vector<int>> grid;
    auto key = [&](const Vec3& p) {
        const long long ix = static_cast<long long>(std::floor(p.x / D)), iy = static_cast<long long>(std::floor(p.y / D)), iz = static_cast<long long>(std::floor(p.z / D));
        return (ix * 73856093LL) ^ (iy * 19349663LL) ^ (iz * 83492791LL);
    };
    for (size_t i = 0; i < all2.size(); ++i) grid[key(all2[i].p)].push_back(static_cast<int>(i));
    for (size_t i = 0; i < all2.size(); ++i) {
        const Sample& a = all2[i];
        double best = 1e300;
        int bestPiece = -1;
        const long long ix = static_cast<long long>(std::floor(a.p.x / D)), iy = static_cast<long long>(std::floor(a.p.y / D)), iz = static_cast<long long>(std::floor(a.p.z / D));
        for (long long dx = -1; dx <= 1; ++dx)
            for (long long dy = -1; dy <= 1; ++dy)
                for (long long dz = -1; dz <= 1; ++dz) {
                    const long long k = ((ix + dx) * 73856093LL) ^ ((iy + dy) * 19349663LL) ^ ((iz + dz) * 83492791LL);
                    const auto it = grid.find(k);
                    if (it == grid.end()) continue;
                    for (int j : it->second) {
                        const Sample& o2 = all2[j];
                        if (o2.strand == a.strand && std::fabs(o2.arc - a.arc) < 3 * D) continue;
                        const double d = norm(o2.p - a.p);
                        if (d < best) {
                            best = d;
                            bestPiece = o2.piece;
                        }
                    }
                }
        c.minClearance = std::min(c.minClearance, best);
        if (best < 0.8 * D && bestPiece >= 0) {
            const std::pair<int, int> key2{std::min(a.piece, bestPiece), std::max(a.piece, bestPiece)};
            YarnCheck::Contact& ct = pairs[key2];
            if (ct.samples == 0 || best < ct.distance) ct.distance = best;
            ct.pieceA = key2.first;
            ct.pieceB = key2.second;
            ++ct.samples;
        }
        const bool contact = best < D;
        c.contactSamples += contact;
        c.overlapSamples += best < 0.8 * D;
        c.severeSamples += best < 0.5 * D;
        if (a.face >= 0 && faceIndex[a.face] >= 0) {
            YarnFaceCheck& fc = c.faces[faceIndex[a.face]];
            fc.minClearance = std::min(fc.minClearance, best);
            fc.contactSamples += contact;
        }
    }
    for (const auto& [k2, ct] : pairs) c.contacts.push_back(ct);
    std::sort(c.contacts.begin(), c.contacts.end(), [](const YarnCheck::Contact& x, const YarnCheck::Contact& y) { return x.distance < y.distance; });
    for (const YarnStrand& sd : e.strands)
        for (size_t k = 1; k + 1 < sd.samples.size(); ++k) {
            const Vec3 a = sd.samples[k] - sd.samples[k - 1], bb = sd.samples[k + 1] - sd.samples[k];
            const double la = norm(a), lb = norm(bb), lc = norm(a + bb), cr = norm(cross(a, bb));
            if (la < 1e-9 || lb < 1e-9) continue;
            const double r = cr > 1e-15 ? la * lb * lc / (2 * cr) : 1e300;
            c.minRadius = std::min(c.minRadius, r);
            const int kind = static_cast<int>(e.pieces[sd.samplePiece[k]].kind);
            ++c.samplesByKind[kind];
            if (r < 0.5 * D) {
                ++c.kinkSamples;
                ++c.kinkByKind[kind];
            }
            const int f = sd.sampleFace[k];
            if (f >= 0 && faceIndex[f] >= 0) c.faces[faceIndex[f]].minRadius = std::min(c.faces[faceIndex[f]].minRadius, r);
        }
    if (c.openStrands != 1 || c.closedStrands != 0) {
        std::snprintf(b, sizeof b, "%d open and %d closed strands; one open strand (two yarn ends) expected", c.openStrands, c.closedStrands);
        c.defects.insert(c.defects.begin(), b);
    }
    const bool construction = c.portsUnmatched == 0 && c.portsOverused == 0 && c.facesWithoutTemplate == 0 && c.facesCrossingsOk == c.facesChecked && c.openStrands == 1 &&
                              c.closedStrands == 0 && c.maxPortGap < 1e-6 && c.maxTangentJump < 1.0;
    c.status = construction ? Status::Pass : Status::Fail;
    return c;
}

std::string YarnCheck::text() const {
    char b[1000];
    std::snprintf(b, sizeof b,
                  "yarn construction %s: %d pieces, %d strands (%d open, %d closed), %d yarn ends; ports: %d shared, %d unmatched, %d over-used; continuity at ports: gap %.2g mm, "
                  "tangent %.2g deg; ports off the reference %.2g mm\n"
                  "crossings: %d of %d faces as their construction requires (in the yarn's own plane: %d), %d faces without a template, interlocks %d of %d, %d "
                  "unexpected crossings, smallest depth gap %.2f D\n"
                  "centrelines: length %.0f mm, %d samples; nearest other yarn min %.3f mm; closer than D %d (%.1f %%), than 0.8 D %d, than 0.5 D %d; tightest bend %.3f mm, "
                  "%d samples bend tighter than D / 2\n",
                  toString(status), pieces, strands, openStrands, closedStrands, yarnEnds, portsSharedOk, portsUnmatched, portsOverused, maxPortGap, maxTangentJump,
                  maxPortOffSurface, facesCrossingsOk, facesChecked, facesOwnPlaneOk, facesWithoutTemplate, interlocksVerified, interlocksExpected, unexpectedCrossings,
                  minDepthGap < 1e299 ? minDepthGap : 0.0, totalLength, samples, minClearance < 1e299 ? minClearance : 0.0, contactSamples,
                  samples ? 100.0 * contactSamples / samples : 0.0, overlapSamples, severeSamples, minRadius < 1e299 ? minRadius : 0.0, kinkSamples);
    return b;
}

TriMesh yarnTubes(const YarnEmbedding& e, double r, int sides) {
    TriMesh m;
    for (const YarnStrand& sd : e.strands) {
        const std::vector<Vec3>& P = sd.samples;
        if (P.size() < 2) continue;
        // rotation-minimising frames by the double reflection method (Wang et al. 2008)
        std::vector<Vec3> T(P.size()), U(P.size());
        for (size_t i = 0; i < P.size(); ++i) T[i] = unit(i + 1 < P.size() ? P[i + 1] - P[i] : P[i] - P[i - 1]);
        Vec3 a = std::fabs(T[0].x) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
        U[0] = unit(a - T[0] * dot(a, T[0]));
        for (size_t i = 0; i + 1 < P.size(); ++i) {
            const Vec3 v1 = P[i + 1] - P[i];
            const double c1 = dot(v1, v1);
            if (c1 < 1e-24) {
                U[i + 1] = U[i];
                continue;
            }
            const Vec3 rL = U[i] - v1 * (2.0 / c1 * dot(v1, U[i]));
            const Vec3 tL = T[i] - v1 * (2.0 / c1 * dot(v1, T[i]));
            const Vec3 v2 = T[i + 1] - tL;
            const double c2 = dot(v2, v2);
            U[i + 1] = c2 < 1e-24 ? rL : rL - v2 * (2.0 / c2 * dot(v2, rL));
        }
        const int base = static_cast<int>(m.pos.size());
        for (size_t i = 0; i < P.size(); ++i) {
            const Vec3 V = cross(T[i], U[i]);
            for (int k = 0; k < sides; ++k) {
                const double ang = 2 * kYarnPi * k / sides;
                m.pos.push_back(P[i] + (U[i] * std::cos(ang) + V * std::sin(ang)) * r);
            }
        }
        for (size_t i = 0; i + 1 < P.size(); ++i)
            for (int k = 0; k < sides; ++k) {
                const int a0 = base + static_cast<int>(i) * sides + k, a1 = base + static_cast<int>(i) * sides + (k + 1) % sides;
                const int b0 = a0 + sides, b1 = a1 + sides;
                m.tri.push_back({a0, a1, b1});
                m.tri.push_back({a0, b1, b0});
            }
    }
    m.computeNormals();
    return m;
}

std::vector<Vec3> sampleYarnHermite(const std::vector<Vec3>& knots, const std::vector<Vec3>& tangents, int per, bool spanScaled) {
    return sampleHermite(knots, tangents, per, spanScaled);
}


} // namespace knit
