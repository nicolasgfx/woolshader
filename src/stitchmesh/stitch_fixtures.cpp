// The hand-authored coarse inputs of stage 4 (stitch_fixtures.h), written as the fixtures document writes them.
#include "stitchmesh/stitch_fixtures.h"

#include <algorithm>
#include <cmath>

namespace knit {

namespace {

constexpr double kTau = 6.283185307179586;

// The roles of one face's corner edges as top flags: b bottom, t top, w wale, in corner order.
std::vector<char> flags(const char* roles) {
    std::vector<char> t;
    for (const char* c = roles; *c; ++c) t.push_back(*c == 'b' ? 0 : *c == 't' ? 1 : -1);
    return t;
}

// Rings of vertices around the z axis, ring i with sizes[i] vertices evenly spaced at z = rowHeight i on radius 10; the
// first id of each ring. Between a ring of n and one of n + 1 (an increase or decrease) the quads shear by one stitch
// over the row; concentrating the extra vertex next to the triangle instead (a half-width column) was measured worse
// for the yarn's local signature (10 failing stitches against 1 or 2), so the even spacing stays.
std::vector<int> rings(PolygonMesh& m, const std::vector<int>& sizes, double rowHeight) {
    std::vector<int> first;
    for (size_t i = 0; i < sizes.size(); ++i) {
        first.push_back(static_cast<int>(m.pos.size()));
        for (int j = 0; j < sizes[i]; ++j) m.pos.push_back({std::cos(kTau * j / sizes[i]) * 10, std::sin(kTau * j / sizes[i]) * 10, rowHeight * i});
    }
    return first;
}

// A full row of quads between two rings of n vertices, growing up ("bwtw") or, flipped, down ("twbw").
void fullRow(CoarseFixture& f, int lower, int upper, int n, bool flipped) {
    for (int j = 0; j < n; ++j) {
        f.mesh.addFace({lower + j, lower + (j + 1) % n, upper + (j + 1) % n, upper + j});
        f.top.push_back(flags(flipped ? "twbw" : "bwtw"));
    }
}

void face(CoarseFixture& f, const std::vector<int>& corners, const char* roles) {
    f.mesh.addFace(corners);
    f.top.push_back(flags(roles));
}

CoarseFixture closedRows(bool reverseFirstRow, bool reverseSecondRow, int n, double rowHeight) {
    CoarseFixture f;
    const std::vector<int> r = rings(f.mesh, {n, n, n}, rowHeight);
    fullRow(f, r[0], r[1], n, reverseFirstRow);
    fullRow(f, r[1], r[2], n, reverseSecondRow);
    return f;
}

// Fixtures 5d and 5e: two rows of n rectangles (w = 10 by rowHeight) on the sheet z = 0, counter-clockwise seen from +z;
// the shared line y = rowHeight is both top (the upper row flipped) or both bottom (the lower row flipped).
CoarseFixture mismatchFlat(bool flipLower, int n, double rowHeight) {
    CoarseFixture f;
    const double w = 10;
    for (int row = 0; row < 3; ++row)
        for (int j = 0; j <= n; ++j) f.mesh.pos.push_back({w * j, rowHeight * row, 0});
    auto v = [&](int row, int j) { return row * (n + 1) + j; };
    for (int j = 0; j < n; ++j) face(f, {v(0, j), v(0, j + 1), v(1, j + 1), v(1, j)}, flipLower ? "twbw" : "bwtw");
    for (int j = 0; j < n; ++j) face(f, {v(1, j), v(1, j + 1), v(2, j + 1), v(2, j)}, flipLower ? "bwtw" : "twbw");
    return f;
}

// Fixture 2: rings 0 and 3 of n vertices, rings 1 and 2 of n - 1 plus the pinch vertex a they share (the document: n = 4,
// ring 1 = 4 5 6, a = 7, ring 2 = 8 9 10, ring 3 = 11..14).
CoarseFixture specialQuads(int n, double rowHeight) {
    CoarseFixture f;
    // every ring at the n angles of ring 0, rings 1 and 2 without the last angle, where a sits at mid height
    auto at = [&](int j, double z) { return Vec3{std::cos(kTau * j / n) * 10, std::sin(kTau * j / n) * 10, z}; };
    for (int j = 0; j < n; ++j) f.mesh.pos.push_back(at(j, 0));
    for (int j = 0; j < n - 1; ++j) f.mesh.pos.push_back(at(j, rowHeight));
    const int a = 2 * n - 1;
    f.mesh.pos.push_back(at(n - 1, 1.5 * rowHeight));
    for (int j = 0; j < n - 1; ++j) f.mesh.pos.push_back(at(j, 2 * rowHeight));
    for (int j = 0; j < n; ++j) f.mesh.pos.push_back(at(j, 3 * rowHeight));
    auto v = [&](int ring, int j) { return ring == 0 ? j : ring == 1 ? n + j : ring == 2 ? 2 * n + j : 3 * n - 1 + j; };
    for (int j = 0; j < n - 2; ++j) face(f, {v(0, j), v(0, j + 1), v(1, j + 1), v(1, j)}, "bwtw");
    face(f, {v(0, n - 2), v(0, n - 1), a, v(1, n - 2)}, "bwtw");
    face(f, {v(0, n - 1), v(0, 0), v(1, 0), a}, "bwtw");
    for (int j = 0; j < n - 2; ++j) face(f, {v(1, j), v(1, j + 1), v(2, j + 1), v(2, j)}, "bwtw");
    face(f, {v(1, n - 2), a, v(2, n - 2)}, "btw");   // T+: bottom, top, wale
    face(f, {a, v(1, 0), v(2, 0)}, "bwt");           // T-: bottom, wale, top
    for (int j = 0; j < n - 2; ++j) face(f, {v(2, j), v(2, j + 1), v(3, j + 1), v(3, j)}, "bwtw");
    face(f, {v(2, n - 2), a, v(3, n - 1), v(3, n - 2)}, "bwtw");
    face(f, {a, v(2, 0), v(3, 0), v(3, n - 1)}, "bwtw");
    return f;
}

// Fixture 3: the middle row grows from n to n + 1 faces around; the triangle's apex (the first vertex of ring 1) is at
// the bottom (the document: n = 3, rings 0..2, 3..5, 6..9, 10..13).
CoarseFixture increase(int n, double rowHeight) {
    CoarseFixture f;
    const std::vector<int> r = rings(f.mesh, {n, n, n + 1, n + 1}, rowHeight);
    fullRow(f, r[0], r[1], n, false);
    for (int j = 0; j < n; ++j) face(f, {r[1] + j, r[1] + (j + 1) % n, r[2] + j + 1, r[2] + j}, "bwtw");
    face(f, {r[1], r[2], r[2] + n}, "wtw");   // T: wale, top, wale
    fullRow(f, r[2], r[3], n + 1, false);
    return f;
}

// Fixture 4: the middle row shrinks from n + 1 to n; the triangle's apex (the first vertex of ring 2) is at the top (the
// document: n = 3, rings 0..3, 4..7, 8..10, 11..13).
CoarseFixture decrease(int n, double rowHeight) {
    CoarseFixture f;
    const std::vector<int> r = rings(f.mesh, {n + 1, n + 1, n, n}, rowHeight);
    fullRow(f, r[0], r[1], n + 1, false);
    for (int j = 0; j < n; ++j) face(f, {r[1] + j, r[1] + j + 1, r[2] + (j + 1) % n, r[2] + j}, "bwtw");
    face(f, {r[1] + n, r[1], r[2]}, "bww");   // T: bottom, wale, wale
    fullRow(f, r[2], r[3], n, false);
    return f;
}

// Fixture 5b: arms a (y = 0) and b (y = 20) with the fold line s (y = 10) between them, the fold ending short of the
// sheet's right border so the turn triangle is the notch between the arms' ends.
CoarseFixture hairpin() {
    CoarseFixture f;
    for (int row = 0; row < 3; ++row)
        for (int i = 0; i < 3; ++i) f.mesh.pos.push_back({10.0 * i, 10.0 * row, 0});   // a = 0 1 2, s = 3 4 5, b = 6 7 8
    f.mesh.pos[5].x = 15;
    face(f, {0, 1, 4, 3}, "bwtw");   // A0
    face(f, {1, 2, 5, 4}, "bwtw");   // A1
    face(f, {2, 8, 5}, "bww");       // T: bottom a2 b2, wale b2 s2, wale s2 a2
    face(f, {8, 7, 4, 5}, "bwtw");   // B1
    face(f, {7, 6, 3, 4}, "bwtw");   // B0
    return f;
}

// Fixtures 3 and 4 on a flat sheet: `columns` stitch columns of width 10 and rows of height `rowHeight`; the middle row
// gains (or loses) one column through a full-width triangle at column `at`, every quad a rectangle: the shaping with no
// shear anywhere. Rows are open (yarn ends at both sides).
CoarseFixture shapingFlat(bool increase, int columns, int at, double rowHeight) {
    CoarseFixture f;
    const double w = 10;
    // ring sizes (vertex counts): the narrow rings have columns + 1 vertices, the wide ones columns + 2
    const int narrow = columns + 1, wide = columns + 2;
    auto ring = [&](int vertices, double z, int extraAt) {   // extraAt: the index of the extra vertex, or -1
        const int first = static_cast<int>(f.mesh.pos.size());
        for (int j = 0; j < vertices; ++j) {
            const double x = (extraAt >= 0 && j > extraAt) ? (j - 1) * w : j * w;
            f.mesh.pos.push_back({extraAt >= 0 && j == extraAt ? at * w : x, 0, z});
        }
        return first;
    };
    // the wide ring: vertices 0..at at x = j w, the extra one at index at + 1 at x = (at + 1) w, the rest shifted by one
    auto wideRing = [&](double z) {
        const int first = static_cast<int>(f.mesh.pos.size());
        for (int j = 0; j < wide; ++j) f.mesh.pos.push_back({j * w, 0, z});
        return first;
    };
    auto narrowRing = [&](double z) { return ring(narrow, z, -1); };
    const int r0 = increase ? narrowRing(0) : wideRing(0);
    const int r1 = increase ? narrowRing(rowHeight) : wideRing(rowHeight);
    const int r2 = increase ? wideRing(2 * rowHeight) : narrowRing(2 * rowHeight);
    const int r3 = increase ? wideRing(3 * rowHeight) : narrowRing(3 * rowHeight);
    const int n0 = increase ? columns : columns + 1, n2 = increase ? columns + 1 : columns;
    for (int j = 0; j < n0; ++j) face(f, {r0 + j, r0 + j + 1, r1 + j + 1, r1 + j}, "bwtw");
    if (increase) {
        // row 1: quads left of the triangle, the triangle (apex v(1, at) at the bottom, top edge v(2, at) v(2, at + 1)),
        // quads right of it shifted by one column above
        for (int j = 0; j < at; ++j) face(f, {r1 + j, r1 + j + 1, r2 + j + 1, r2 + j}, "bwtw");
        face(f, {r1 + at, r2 + at + 1, r2 + at}, "wtw");   // apex, top right, top left: the top edge runs right to left as the quads' do
        for (int j = at; j < columns; ++j) face(f, {r1 + j, r1 + j + 1, r2 + j + 2, r2 + j + 1}, "bwtw");
    } else {
        // row 1: the triangle's apex v(2, at) at the top, bottom edge v(1, at) v(1, at + 1)
        for (int j = 0; j < at; ++j) face(f, {r1 + j, r1 + j + 1, r2 + j + 1, r2 + j}, "bwtw");
        face(f, {r1 + at, r1 + at + 1, r2 + at}, "bww");
        for (int j = at; j < columns; ++j) face(f, {r1 + j + 1, r1 + j + 2, r2 + j + 1, r2 + j}, "bwtw");
    }
    for (int j = 0; j < n2; ++j) face(f, {r2 + j, r2 + j + 1, r3 + j + 1, r3 + j}, "bwtw");
    return f;
}

// The cut-out of the frozen sphere (stitch_fixtures.h CastOnSideBySide), its numbers as stitchmesh_lab wrote them
// (cutout.obj, cutout_roles.txt): vertex 2 is the pinch where the special quad's course edges meet.
CoarseFixture castOnSideBySide() {
    CoarseFixture f;
    f.mesh.pos = {{9.711385, -13.872085, -25.177562}, {19.671862, -13.275999, -18.765396}, {17.604689, -18.546280, -16.176872},
                  {7.963108, -22.141353, -19.148981}, {16.575765, -24.115910, -7.980436},  {6.985583, -27.707319, -10.125869},
                  {23.752201, -18.071315, -5.358448}};
    face(f, {0, 1, 2}, "wbt");
    face(f, {2, 3, 0}, "twb");
    face(f, {2, 4, 5, 3}, "wtwb");
    face(f, {2, 1, 6, 4}, "bwtw");
    // its own facets as the reference surface (the frozen sphere's yarn is placed on its coarse mesh's triangles)
    for (const Vec3& p : f.mesh.pos) f.surface.pos.push_back(p);
    f.surface.tri = {{0, 1, 2}, {2, 3, 0}, {2, 4, 5}, {2, 5, 3}, {2, 1, 6}, {2, 6, 4}};
    f.surface.computeNormals();
    return f;
}

CoarseFixture castOnKpReceivers() {
    CoarseFixture f;
    f.mesh.pos = {{-11.587963, -10.470695, -26.015205}, {-0.984311, -13.396062, -27.163754}, {-2.811451, -22.141953, -20.570973},
                  {-13.546843, -18.971028, -19.427501}, {-20.388896, -5.852829, -21.629792},  {-22.414675, -13.099833, -15.674262},
                  {-3.910817, -27.915998, -11.105376},  {-14.671774, -24.574525, -10.051303}, {-23.720960, -17.843171, -6.159942}};
    face(f, {0, 1, 2, 3}, "twbw");
    face(f, {4, 0, 3}, "tww");
    face(f, {3, 5, 4}, "btw");
    face(f, {3, 2, 6, 7}, "bwtw");
    face(f, {5, 3, 7, 8}, "bwtw");
    for (const Vec3& p : f.mesh.pos) f.surface.pos.push_back(p);
    f.surface.tri = {{0, 1, 2}, {0, 2, 3}, {4, 0, 3}, {3, 5, 4}, {3, 2, 6}, {3, 6, 7}, {5, 3, 7}, {5, 7, 8}};
    f.surface.computeNormals();
    return f;
}

CoarseFixture castOnPoleKp() {
    CoarseFixture f;
    f.mesh.pos = {{-16.245946, 10.697224, -23.268563}, {-6.816252, 10.146975, -27.769254}, {-8.929530, -0.395916, -28.991967},
                  {-17.932135, 2.204485, -24.397470},  {-5.555183, 19.327012, -22.632444}, {-23.642169, 14.397062, -12.458380},
                  {-15.425155, 19.759728, -17.051561}, {-4.758364, 23.826177, -18.078012}, {-24.499018, 3.946672, -17.437234}};
    face(f, {0, 1, 2}, "twb");
    face(f, {2, 3, 0}, "wtb");
    face(f, {4, 1, 0}, "wbt");
    face(f, {0, 5, 6}, "wtw");
    face(f, {4, 0, 6, 7}, "bwtw");
    face(f, {0, 3, 8, 5}, "bwtw");
    for (const Vec3& p : f.mesh.pos) f.surface.pos.push_back(p);
    f.surface.tri = {{0, 1, 2}, {2, 3, 0}, {4, 1, 0}, {0, 5, 6}, {4, 0, 6}, {4, 6, 7}, {0, 3, 8}, {0, 8, 5}};
    f.surface.computeNormals();
    return f;
}

} // namespace

const char* stitchFixtureName(StitchFixture which) {
    switch (which) {
    case StitchFixture::ClosedRows: return "closed-rows";
    case StitchFixture::SpecialQuads: return "special-quads";
    case StitchFixture::Increase: return "kp";
    case StitchFixture::Decrease: return "d12k";
    case StitchFixture::MismatchBetweenRows: return "mismatch-between-rows";
    case StitchFixture::MismatchSplitting: return "mismatch-splitting";
    case StitchFixture::MismatchFlatJoining: return "mismatch-flat-joining";
    case StitchFixture::MismatchFlatSplitting: return "mismatch-flat-splitting";
    case StitchFixture::MismatchWithinRow: return "mismatch-within-row";
    case StitchFixture::IncreaseFlat: return "kp-flat";
    case StitchFixture::CastOnSideBySide: return "cast-on-side-by-side";
    case StitchFixture::CastOnKpReceivers: return "cast-on-kp-receivers";
    case StitchFixture::CastOnPoleKp: return "cast-on-pole-kp";
    default: return "d12k-flat";
    }
}

bool parseStitchFixture(const std::string& name, StitchFixture& out) {
    const StitchFixture all[14] = {StitchFixture::ClosedRows, StitchFixture::SpecialQuads, StitchFixture::Increase, StitchFixture::Decrease,
                                   StitchFixture::MismatchBetweenRows, StitchFixture::MismatchSplitting, StitchFixture::MismatchFlatJoining,
                                   StitchFixture::MismatchFlatSplitting, StitchFixture::MismatchWithinRow, StitchFixture::IncreaseFlat, StitchFixture::DecreaseFlat,
                                   StitchFixture::CastOnSideBySide, StitchFixture::CastOnKpReceivers,
                                   StitchFixture::CastOnPoleKp};
    for (StitchFixture f : all)
        if (name == stitchFixtureName(f)) {
            out = f;
            return true;
        }
    return false;
}

// The cylinder of radius 10 the ring fixtures lie on, from zLo to zHi: 144 segments round, rings every millimetre,
// triangles counter-clockwise seen from outside.
TriMesh cylinderSurface(double zLo, double zHi) {
    TriMesh t;
    const int around = 144, rings = std::max(2, static_cast<int>(std::ceil(zHi - zLo)) + 1);
    for (int i = 0; i < rings; ++i) {
        const double z = zLo + (zHi - zLo) * i / (rings - 1);
        for (int j = 0; j < around; ++j) t.pos.push_back({std::cos(kTau * j / around) * 10, std::sin(kTau * j / around) * 10, z});
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

// The plane of a flat fixture: its vertices' bounding box in the plane of its first face, grown by `margin`, as a grid
// of 2 mm cells wound like that face.
TriMesh planeSurface(const PolygonMesh& m, double margin) {
    TriMesh t;
    const Vec3 n = m.faceNormal(0);
    Vec3 lo = m.pos[0], hi = m.pos[0];
    for (const Vec3& p : m.pos) { lo = componentMin(lo, p); hi = componentMax(hi, p); }
    const Vec3 c = (lo + hi) * 0.5;
    // the two in-plane axes: the longest bounding-box axis projected into the plane, and the normal's cross with it
    Vec3 u = hi - lo;
    u = u - n * dot(u, n);
    if (norm(u) < 1e-9) u = std::fabs(n.x) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
    u = normalized(u - n * dot(u, n));
    const Vec3 v = cross(n, u);
    double eu = 0, ev = 0;
    for (const Vec3& p : m.pos) { eu = std::max(eu, std::fabs(dot(p - c, u))); ev = std::max(ev, std::fabs(dot(p - c, v))); }
    eu += margin;
    ev += margin;
    const int nu = std::max(2, static_cast<int>(std::ceil(eu))), nv = std::max(2, static_cast<int>(std::ceil(ev)));
    for (int i = 0; i <= nv; ++i)
        for (int j = 0; j <= nu; ++j) t.pos.push_back(c + u * (-eu + 2 * eu * j / nu) + v * (-ev + 2 * ev * i / nv));
    for (int i = 0; i < nv; ++i)
        for (int j = 0; j < nu; ++j) {
            const int a = i * (nu + 1) + j, b = a + 1, cc = a + nu + 2, d = a + nu + 1;
            t.tri.push_back({a, b, cc});
            t.tri.push_back({a, cc, d});
        }
    t.computeNormals();
    return t;
}

CoarseFixture stitchFixture(StitchFixture which, int around, double radius) {
    const bool fine = around > 0;
    const double rowHeight = fine ? 5.0 : 10.0;
    CoarseFixture f;
    switch (which) {
    case StitchFixture::ClosedRows: f = closedRows(false, false, fine ? around : 3, rowHeight); break;
    case StitchFixture::SpecialQuads: f = specialQuads(fine ? around : 4, rowHeight); break;
    case StitchFixture::Increase: f = increase(fine ? around : 3, rowHeight); break;
    case StitchFixture::Decrease: f = decrease(fine ? around : 3, rowHeight); break;
    case StitchFixture::MismatchBetweenRows: f = closedRows(false, true, fine ? around : 3, rowHeight); break;
    case StitchFixture::MismatchSplitting: f = closedRows(true, false, fine ? around : 3, rowHeight); break;
    case StitchFixture::MismatchFlatJoining: f = mismatchFlat(false, fine ? around : 3, rowHeight); break;
    case StitchFixture::MismatchFlatSplitting: f = mismatchFlat(true, fine ? around : 3, rowHeight); break;
    case StitchFixture::MismatchWithinRow: f = hairpin(); break;
    case StitchFixture::IncreaseFlat: f = shapingFlat(true, fine ? around : 6, (fine ? around : 6) / 2, rowHeight); break;
    case StitchFixture::CastOnSideBySide: f = castOnSideBySide(); f.name = stitchFixtureName(which); return f;
    case StitchFixture::CastOnKpReceivers: f = castOnKpReceivers(); f.name = stitchFixtureName(which); return f;
    case StitchFixture::CastOnPoleKp: f = castOnPoleKp(); f.name = stitchFixtureName(which); return f;
    default: f = shapingFlat(false, fine ? around : 6, (fine ? around : 6) / 2, rowHeight); break;
    }
    f.name = stitchFixtureName(which);
    const bool onCylinder = which != StitchFixture::MismatchWithinRow && which != StitchFixture::IncreaseFlat && which != StitchFixture::DecreaseFlat &&
                            which != StitchFixture::MismatchFlatJoining && which != StitchFixture::MismatchFlatSplitting;
    if (onCylinder) {
        double zLo = f.mesh.pos[0].z, zHi = zLo;
        for (const Vec3& p : f.mesh.pos) { zLo = std::min(zLo, p.z); zHi = std::max(zHi, p.z); }
        f.surface = cylinderSurface(zLo - rowHeight, zHi + rowHeight);
        if (radius != 10.0) {
            const double k = radius / 10.0;
            for (Vec3& p : f.mesh.pos) p = {p.x * k, p.y * k, p.z};
            for (Vec3& p : f.surface.pos) p = {p.x * k, p.y * k, p.z};
            f.surface.computeNormals();
        }
    } else f.surface = planeSurface(f.mesh, rowHeight);
    return f;
}

} // namespace knit
