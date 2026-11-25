// The yarn body of structural fibres (src/fibres/yarn_body.h): its consequential invariants on small
// hand-made yarns. `knit_tests_yarn_body --list`.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <utility>
#include <vector>

#include "fibres/strands_file.h"
#include "fibres/yarn_bake.h"
#include "fibres/yarn_body.h"
#include "harness.h"
#include "model/cancel.h"
#include "knitting/tube/param_text.h"
#include "model/colour.h"
#include "knitting/tube/yarn_surface.h"
#include "strand_fixtures.h"

using namespace knit;

namespace {

constexpr double kD = 1.5;          // the swatch's drawn diameter [mm]
constexpr uint64_t kPartSeed = 7;

PatternParams yarnParams() {
    PatternParams p;
    p.yarnDiameter = kD;
    p.clearance = 0.1;
    return p;
}

// A straight yarn along x, `length` mm, controls every 0.5 mm.
YarnPath straightYarn(double length, double y = 0, int id = 0) {
    std::vector<Vec3> c;
    for (int i = -1; i <= static_cast<int>(std::lround(length / 0.5)) + 1; ++i) c.push_back({0.5 * i, y, 0});
    return makeYarn(c, id);
}

// A yarn with an inflection (as rmf_frames), scaled to a knitting yarn's bends (radius about 2 mm);
// tight: loop-like bends of about 1 mm radius, where the inner fibres turn on 0.3 mm.
YarnPath bentYarn(bool tight = false) {
    std::vector<Vec3> c;
    for (int i = -1; i <= 40; ++i) {
        const double x = 0.5 * i;
        c.push_back(tight ? Vec3{x, 1.2 * std::sin(0.9 * x), 0.15 * x} : Vec3{x, 2.5 * std::sin(0.45 * x), 0.15 * x});
    }
    return makeYarn(c, 0);
}

SurfaceModel roundSurface(const CurveSet& cs, double ringSpacing = 0.025) {
    SurfaceParams sp;
    sp.compression = false;
    sp.ringSpacing = ringSpacing;
    return buildSurface(cs, yarnParams(), sp);
}

YarnBody bodyOf(const CurveSet& cs, const SurfaceModel& sm, const YarnBodyParams& p, uint64_t seed = kPartSeed) {
    YarnBody b;
    std::string err;
    if (!buildYarnBody(cs, sm, p, seed, b, &err, ExecPolicy{true})) std::printf("    %s\n", err.c_str());
    CHECK(err.empty());
    return b;
}

StrandSetF asFloat(const BodyStrands& b) {
    StrandSetF f;
    f.xyz.reserve(3 * b.points.size());
    for (size_t v = 0; v < b.points.size(); ++v) {
        f.xyz.push_back(static_cast<float>(b.points[v].x));
        f.xyz.push_back(static_cast<float>(b.points[v].y));
        f.xyz.push_back(static_cast<float>(b.points[v].z));
        f.radius.push_back(static_cast<float>(b.radius[v]));
    }
    f.first.clear();
    for (const uint64_t v : b.first) f.first.push_back(static_cast<uint32_t>(v));  // the test bodies are small
    return f;
}

// The ellipse radius of a point in its ring: hypot(along e_a / a, along e_b / b); 1 = on the ellipse.
double normalizedRadius(const Ring& g, const Vec3& p) {
    const double cp = std::cos(g.phi), sp = std::sin(g.phi);
    const Vec3 ea = g.n1 * cp + g.n2 * sp, eb = g.n2 * cp - g.n1 * sp;
    const Vec3 d = p - g.c;
    return std::hypot(dot(d, ea) / g.a, dot(d, eb) / g.b);
}

} // namespace

// Every radius finite, positive and at most r_max; every strand obeys the renderer's hard rules
// (segment >= 4 x its larger radius, <= 200 x the smaller, radius step <= 50 %, curvature), on tight
// loop-like bends whose section is flattened to 0.6 and turned (a compressed yarn at its contacts).
TEST(yarn_body_radii_and_rules) {
    const CurveSet cs = makeSet({bentYarn(true)});
    SurfaceModel sm = roundSurface(cs);
    for (Ring& g : sm.yarns[0].rings) {
        g.b *= 0.6 + 0.4 * std::fabs(std::cos(0.7 * g.s));
        g.phi = 0.5 * g.s;
    }
    for (int plies = 1; plies <= 3; ++plies)
        for (double iota : {0.0, 1.0}) {
            YarnBodyParams p;
            p.plies = plies;
            p.irregularity = iota;
            p.fibres = 90;
            const YarnBody b = bodyOf(cs, sm, p);
            std::printf("    plies %d iota %.0f: %s\n", plies, iota, b.stats.line().c_str());
            CHECK(b.strands.strands() > 0 && b.stats.segments > 0);
            for (double r : b.strands.radius) CHECK(std::isfinite(r) && r > 0 && r <= b.stats.layout.rMax * (1 + 1e-12));
            StrandRules rules;
            rules.minSegOverRMax = 4;
            const StrandCheck check = checkStrands(asFloat(b.strands), rules);
            if (!check.ok) std::printf("    %s\n", check.error.c_str());
            CHECK(check.ok);
        }
}

// The angle of the surface fibres against their ply's axis is the asked fibre angle, measured against
// the ply's tangent (not the yarn's): a straight yarn, regular fibres in a thin shell. A fibre point
// belongs to the ply's section at its nominal z (the plane normal to the ply there), which lies off the
// yarn's section by the axial offset; the tangent is taken there.
TEST(yarn_body_fibre_angle_straight) {
    const CurveSet cs = makeSet({straightYarn(30)});
    const SurfaceModel sm = roundSurface(cs);
    const YarnSurface& ys = sm.yarns[0];
    const Ring g0 = ringAt(cs, ys, 0);
    for (int plies = 1; plies <= 3; ++plies)
        for (double plyAngle : {0.0, 20.0, 35.0}) {
            if (plies == 1 && plyAngle > 0) continue;
            YarnBodyParams p;
            p.plies = plies;
            p.plyAngleDeg = plyAngle;
            p.irregularity = 0;
            p.shellDepth = 0.05;
            p.fibres = 10 * plies;
            const YarnBody b = bodyOf(cs, sm, p);
            const int perPly = p.fibres / plies;
            // Nominal coordinates of a point of the straight yarn: (n1, n2) components and arc length.
            const auto nominal = [&](const Vec3& q) {
                const Vec3 d = q - g0.c;
                return Vec3{dot(d, g0.n1), dot(d, g0.n2), dot(d, g0.T)};
            };
            // The z whose ply section (the plane through the ply centre normal to its tangent) holds F.
            const auto sectionZ = [&](int ply, const Vec3& F) {
                double z = F.z;
                for (int it = 0; it < 30; ++it) {
                    const PlyAxis a = bodyPlyAxis(p, ys.radius, kPartSeed, 0, ply, z);
                    const double f = dot(Vec3{F.x - a.centre.x, F.y - a.centre.y, F.z - z}, a.tangent);
                    z += f / a.tangent.z;
                    if (std::fabs(f) < 1e-12) break;
                }
                return z;
            };
            double sum = 0, lo = 90, hi = 0;
            int count = 0;
            for (size_t i = 0; i < b.strands.strands(); ++i) {
                const int ply = static_cast<int>(i) / perPly;
                for (uint64_t v = b.strands.first[i]; v + 1 < b.strands.first[i + 1]; ++v) {
                    const Vec3 F0 = nominal(b.strands.points[v]), F1 = nominal(b.strands.points[v + 1]);
                    const Vec3 mid = (F0 + F1) * 0.5;
                    if (mid.z < 2 || mid.z > ys.length - 2) continue;  // away from the yarn's ends
                    const Vec3 t = bodyPlyAxis(p, ys.radius, kPartSeed, 0, ply, sectionZ(ply, mid)).tangent;
                    const Vec3 d = F1 - F0;
                    const double angle = std::acos(std::clamp(dot(d, t) / norm(d), -1.0, 1.0)) * 180 / kPi;
                    sum += angle;
                    lo = std::min(lo, angle);
                    hi = std::max(hi, angle);
                    ++count;
                }
            }
            const double mean = sum / std::max(count, 1);
            std::printf("    plies %d, ply angle %2.0f: fibre vs ply %.2f deg (range %.1f to %.1f, %d segments; asked %.0f at the shell's "
                        "outer radius)\n", plies, plyAngle, mean, lo, hi, count, p.fibreAngleDeg);
            CHECK(count > 50 && std::fabs(mean - p.fibreAngleDeg) <= 1.0);
        }
}

// Through bends: no frame flip and no jump. Every segment stays short (the strand rules), the ply axis
// moves continuously along the yarn, and the turn per segment stays near its 6 degree target.
TEST(yarn_body_continuous_through_bends) {
    const CurveSet cs = makeSet({bentYarn()});
    const SurfaceModel sm = roundSurface(cs);
    for (int plies = 1; plies <= 3; ++plies) {
        YarnBodyParams p;
        p.plies = plies;
        p.fibres = 60;
        const YarnBody b = bodyOf(cs, sm, p);
        std::printf("    plies %d: soft-turn segments %lld of %lld, p99 turn %.1f deg, max %.1f\n", plies, b.stats.softTurnSegments,
                    b.stats.segments, b.stats.p99TurnDeg, b.stats.maxTurnDeg);
        CHECK(b.stats.p99TurnDeg <= 6.5);
        // The ply's phase: continuous on a step far below the integration grid.
        const double R = b.stats.layout.R, omega = plies > 1 ? std::tan(20 * kPi / 180) / R : 0;
        for (double z = 0; z + 0.01 < sm.yarns[0].length; z += 0.01) {
            const Vec2 a = bodyPlyAxis(p, sm.yarns[0].radius, kPartSeed, 0, 0, z).centre;
            const Vec2 c = bodyPlyAxis(p, sm.yarns[0].radius, kPartSeed, 0, 0, z + 0.01).centre;
            CHECK(norm(Vec2{c.x - a.x, c.y - a.y}) <= 1.2 * R * omega * 0.01 + 1e-12);
        }
    }
}

// The same inputs give the same fibres, bit for bit; another part, yarn or ply gives other fibres.
TEST(yarn_body_deterministic_unique) {
    const CurveSet cs = makeSet({straightYarn(20, 0, 0), straightYarn(20, 5, 1)});
    const SurfaceModel sm = roundSurface(cs);
    YarnBodyParams p;
    p.plies = 2;
    p.fibres = 40;
    const YarnBody a = bodyOf(cs, sm, p), b = bodyOf(cs, sm, p), other = bodyOf(cs, sm, p, kPartSeed + 1);
    CHECK(a.strands.points.size() == b.strands.points.size());
    for (size_t v = 0; v < a.strands.points.size(); ++v)
        CHECK(a.strands.points[v].x == b.strands.points[v].x && a.strands.points[v].y == b.strands.points[v].y &&
              a.strands.points[v].z == b.strands.points[v].z && a.strands.radius[v] == b.strands.radius[v]);
    CHECK(a.strands.points[0].y != other.strands.points[0].y || a.strands.points[0].z != other.strands.points[0].z);
    // The two yarns are the same straight yarn 5 mm apart: their first fibres differ in their sections.
    size_t firstOfYarn1 = 0;
    while (a.strands.yarnOf[firstOfYarn1] == 0) ++firstOfYarn1;
    const Vec3 p0 = a.strands.points[a.strands.first[0]], p1 = a.strands.points[a.strands.first[firstOfYarn1]];
    CHECK(std::fabs(p0.z - p1.z) > 1e-6 || std::fabs(p0.y - (p1.y - 5)) > 1e-6);
    // Ply 0 and ply 1 of one yarn: different start depths.
    std::vector<double> depth0, depth1;
    for (size_t i = 0; i < a.strands.strands() && a.strands.yarnOf[i] == 0; ++i) {
        const Vec3 q = a.strands.points[a.strands.first[i]];
        (i < 20 ? depth0 : depth1).push_back(std::hypot(q.y, q.z));
    }
    std::sort(depth0.begin(), depth0.end());
    std::sort(depth1.begin(), depth1.end());
    CHECK(depth0 != depth1);
}

// Every centreline lies inside its ring's ellipse, through thick and thin spots and an elliptical
// (compressed) section; the fibres reach the local surface everywhere (the envelope is followed, not a
// fixed cylinder); the measured protrusion stays within r_max (1 - min(a, b) / r).
TEST(yarn_body_follows_envelope) {
    YarnPath y = straightYarn(24);
    y.thickness.assign(y.curve.ctrl.size(), 1.0f);
    for (size_t c = 0; c < y.thickness.size(); ++c) y.thickness[c] = static_cast<float>(1.0 + 0.4 * std::sin(0.35 * static_cast<double>(c)));
    const CurveSet cs = makeSet({y});
    SurfaceModel sm = roundSurface(cs);
    const double r = sm.yarns[0].radius;
    double minAxis = 1e300;
    for (size_t k = 0; k < sm.yarns[0].rings.size(); ++k) {
        Ring& g = sm.yarns[0].rings[k];
        if (g.s > 12) {  // the second half: flattened to 0.7 and turned
            g.b *= 0.7;
            g.phi = 0.4 * g.s;
        }
        minAxis = std::min(minAxis, std::min(g.a, g.b));
    }
    for (int plies = 1; plies <= 3; ++plies) {
        YarnBodyParams p;
        p.plies = plies;
        p.fibres = 120;
        const YarnBody b = bodyOf(cs, sm, p);
        const BodyLayout& l = b.stats.layout;
        std::vector<double> reach(24, 0.0);
        double worst = 0, worstAway = 0;
        for (size_t v = 0; v < b.strands.points.size(); ++v) {
            const Ring g = ringAt(cs, sm.yarns[0], b.strands.s[v]);
            const double t = normalizedRadius(g, b.strands.points[v]);
            worst = std::max(worst, t);
            // The smoothed section (sigma 0.1 mm) differs from the rings only near the abrupt step at s = 12.
            if (std::fabs(b.strands.s[v] - 12) > 0.4) worstAway = std::max(worstAway, t);
            const size_t bucket = std::min<size_t>(23, static_cast<size_t>(b.strands.s[v]));
            reach[bucket] = std::max(reach[bucket], t);
        }
        std::printf("    plies %d: centreline in its ring's ellipse at most %.4f (%.4f away from the step at s = 12; bound %.4f), "
                    "max protrusion %.4f mm (bound %.4f)\n", plies, worst, worstAway, (r - l.rMax) / r, b.stats.maxProtrusion,
                    l.rMax * std::max(0.0, 1 - minAxis / r));
        CHECK(std::fabs(b.stats.maxEnvelope - worst) < 1e-12);    // the generator reports what is there
        // Away from the abrupt step every centreline lies in its true ring's ellipse (the smoothed map of an
        // ellipse that turns 0.4 rad / mm is a little rounder, so the radius reserve is not kept exactly);
        // at the step of a 30 % flattening the excess stays small.
        CHECK(worstAway <= 1.0);
        CHECK(worst <= 1.15);
        CHECK(b.stats.maxProtrusion <= l.rMax * std::max(0.0, 1 - minAxis / r) + 1e-12);
        for (size_t k = 1; k + 1 < reach.size(); ++k) CHECK(reach[k] >= 0.85 * (r - l.rMax) / r);
    }
}

// Regular fibres (irregularity 0), as README defines them: every radius r_f, one piece per slot, no
// crimp in the layout, nothing dropped.
TEST(yarn_body_regular_definition) {
    const CurveSet cs = makeSet({bentYarn()});
    const SurfaceModel sm = roundSurface(cs);
    for (int plies = 1; plies <= 3; ++plies) {
        YarnBodyParams p;
        p.plies = plies;
        p.irregularity = 0;
        p.fibres = 60;
        const YarnBody b = bodyOf(cs, sm, p);
        CHECK(b.stats.layout.rMax == p.fibreRadius && b.stats.layout.crimp == 0);
        CHECK(static_cast<long long>(b.strands.strands()) == p.fibres && b.stats.pieces == p.fibres && b.stats.droppedPieces == 0);
        for (double rr : b.strands.radius) CHECK(rr == p.fibreRadius);
    }
}

// Interior fibres come from their own stream: the shell's strands are the same with or without them.
TEST(yarn_body_interior_keeps_shell) {
    const CurveSet cs = makeSet({bentYarn()});
    const SurfaceModel sm = roundSurface(cs);
    YarnBodyParams p;
    p.plies = 2;
    p.fibres = 60;
    p.shellDepth = 0.4;
    const YarnBody shell = bodyOf(cs, sm, p);
    p.interiorFibres = 80;
    const YarnBody both = bodyOf(cs, sm, p);
    CHECK(both.strands.strands() > shell.strands.strands());
    for (size_t v = 0; v < shell.strands.points.size(); ++v)
        CHECK(shell.strands.points[v].x == both.strands.points[v].x && shell.strands.points[v].y == both.strands.points[v].y &&
              shell.strands.points[v].z == both.strands.points[v].z && shell.strands.radius[v] == both.strands.radius[v]);
}

// Refusals are predictable and say what would be accepted; a full shell is either placed completely or
// refused with the count placed, never built with fewer fibres.
TEST(yarn_body_validation) {
    const double r = 0.5 * kD;
    std::string err;
    YarnBodyParams p;
    for (int plies : {0, 4}) {
        p.plies = plies;
        CHECK(!yarnBodyParamsValid(p, r, &err));
    }
    p = YarnBodyParams{};
    p.shellDepth = 1;
    p.core = BodyCore::Opaque;
    CHECK(!yarnBodyParamsValid(p, r, &err) && err.find("no room for a core") != std::string::npos);
    p = YarnBodyParams{};
    p.fibres = 100000;
    CHECK(!yarnBodyParamsValid(p, r, &err) && err.find("at most") != std::string::npos);
    p = YarnBodyParams{};
    p.fibreRadius = 0.2;
    CHECK(!yarnBodyParamsValid(p, r, &err) && err.find("too large") != std::string::npos);
    p = YarnBodyParams{};
    p.interiorFibres = 10;
    p.shellDepth = 1;
    CHECK(!yarnBodyParamsValid(p, r, &err));
    // At the capacity limit: all placed or refused by count.
    p = YarnBodyParams{};
    p.plies = 3;
    p.irregularity = 0;
    p.fibres = 100000;
    CHECK(!yarnBodyParamsValid(p, r, &err));
    const int perYarn = std::stoi(err.substr(err.find(", ") + 2));  // "(at most N per ply, M per yarn)"
    p.fibres = perYarn;
    CHECK(yarnBodyParamsValid(p, r, &err));
    const CurveSet cs = makeSet({straightYarn(10)});
    const SurfaceModel sm = roundSurface(cs);
    YarnBody b;
    err.clear();
    const bool built = buildYarnBody(cs, sm, p, kPartSeed, b, &err, ExecPolicy{true});
    std::printf("    %d fibres per yarn at the limit: %s\n", perYarn, built ? "placed" : err.c_str());
    CHECK(built ? static_cast<int>(b.strands.strands()) == perYarn : err.find("placed") != std::string::npos);
}

// The overlap measure finds two fibres whose segments cross between far-apart vertices, and nothing
// between separated ones.
TEST(yarn_body_overlap_measure) {
    const auto crossing = [](double gap) {
        BodyStrands b;
        b.points = {{-1, 0, 0}, {1, 0, 0}, {0, -1, gap}, {0, 1, gap}};
        b.radius = {0.012, 0.012, 0.012, 0.012};
        b.s = {4, 6, 4, 6};
        b.first = {0, 2, 4};
        b.yarnOf = {0, 0};
        return measureOverlap(b, {{0, 0, 10}});
    };
    const OverlapStats touching = crossing(0.01), apart = crossing(0.05);
    CHECK(touching.pairsWithinYarn == 1 && touching.segmentsInvolved == 2 && touching.maxDepth > 0.5);
    CHECK(apart.pairsWithinYarn == 0 && apart.segmentsInvolved == 0);
}

// The fibres do not depend on how densely the envelope was sampled: at model rings 0.025 and 0.05 mm
// every vertex of one lies on the other's polyline, within the chord tolerance (0.25 r_f) and a micron.
TEST(yarn_body_sampling_independent) {
    const CurveSet cs = makeSet({bentYarn()});
    const SurfaceModel fine = roundSurface(cs, 0.025), coarse = roundSurface(cs, 0.05);
    YarnBodyParams p;
    p.plies = 3;
    p.fibres = 45;
    const YarnBody a = bodyOf(cs, fine, p), b = bodyOf(cs, coarse, p);
    CHECK(a.strands.strands() == b.strands.strands());
    double worst = 0;
    for (size_t i = 0; i < a.strands.strands(); ++i)
        for (uint64_t v = a.strands.first[i]; v < a.strands.first[i + 1]; ++v) {
            double best = 1e300;
            for (uint64_t w = b.strands.first[i]; w + 1 < b.strands.first[i + 1]; ++w)
                best = std::min(best, closestSegmentPoints(a.strands.points[v], a.strands.points[v], b.strands.points[w],
                                                           b.strands.points[w + 1]).distance);
            worst = std::max(worst, best);
        }
    std::printf("    largest distance of a vertex from the other sampling's polyline: %.2f um\n", worst * 1000);
    CHECK(worst <= 0.25 * p.fibreRadius + 0.001);
}

// A body may pass 2^32 points (400 fibres per yarn on a long yarn): the merge's offsets are 64 bit. A
// body table whose end lies past 2^32 (the points themselves are not allocated: the offset is the
// table's end) takes a zone's strands at the right offsets, not wrapped.
TEST(yarn_body_offsets_64bit) {
    BodyStrands to;
    to.first = {0, 5'000'000'000ull};
    to.yarnOf = {0};
    BodyStrands zone;
    zone.points.resize(5);
    zone.radius.resize(5);
    zone.s.resize(5);
    zone.first = {0, 2, 5};
    zone.yarnOf = {1, 1};
    appendStrands(to, zone);
    CHECK(to.strands() == 3 && to.first.size() == 4);
    CHECK(to.first[2] == 5'000'000'002ull && to.first[3] == 5'000'000'005ull);
}

// A cancelled export stops in the yarn body's build (a checkpoint per ply zone), serial and parallel.
TEST(yarn_body_cancel) {
    const CurveSet cs = makeSet({straightYarn(20, 0, 0), straightYarn(20, 5, 1)});
    const SurfaceModel sm = roundSurface(cs);
    YarnBodyParams p;
    p.plies = 2;
    p.fibres = 40;
    for (bool parallel : {false, true}) {
        CancelToken token;
        token.cancel();
        const ThreadCancelScope scope(&token);
        YarnBody b;
        std::string err;
        bool cancelled = false;
        try {
            (void)buildYarnBody(cs, sm, p, kPartSeed, b, &err, ExecPolicy{parallel});
        } catch (const Cancelled&) {
            cancelled = true;
        }
        CHECK(cancelled);
    }
}

// --- the bake (fibres/yarn_bake.h)

namespace {

// The decoded tangent-space normal and a grey map's value of an atlas pixel.
Vec3 bakedNormal(const BakedYarnMaps& m, int x, int row) {
    const size_t at = 3 * (static_cast<size_t>(row) * m.width + x);
    return {m.normal[at] / 127.5 - 1, m.normal[at + 1] / 127.5 - 1, m.normal[at + 2] / 127.5 - 1};
}
double bakedGrey(const std::vector<uint8_t>& map, const BakedYarnMaps& m, int x, int row) {
    return map[3 * (static_cast<size_t>(row) * m.width + x)] / 255.0;
}

} // namespace

// A straight yarn with regular fibres along it in a shell. Every decoded normal is a unit vector; the
// normals point out of the tube and have no component along the yarn; the hits' heights lie in the shell;
// the fibre direction reads as along s; and the sign convention is plyMapSample's: on the -theta flank
// of a fibre crest the tangent-space normal's y is negative, on the +theta flank positive.
TEST(yarn_bake_straight_regular) {
    const CurveSet cs = makeSet({straightYarn(10)});
    const SurfaceModel sm = roundSurface(cs);
    YarnBodyParams p;
    p.plies = 1;
    p.irregularity = 0;
    p.fibreAngleDeg = 0;
    p.fibres = 120;
    p.shellDepth = 0.5;
    const YarnBody b = bodyOf(cs, sm, p);
    const BodyLayout& l = b.stats.layout;
    BakeParams bp;
    bp.texel = 0.005;
    bp.samples = 2;
    bp.atlasWidth = 2048;
    BakedYarnMaps m;
    std::string err;
    const bool baked = bakeYarnMaps(cs, sm.yarns[0], b.strands, p, kPartSeed, 0, bp, m, &err, ExecPolicy{true});
    if (!baked) std::printf("    %s\n", err.c_str());
    CHECK(baked);
    std::printf("    %d x %d texels, coverage %.3f, fibre angle %.2f deg, %.0f ms\n", m.along, m.around, m.coverage, m.fibreAngleDeg, m.ms);
    CHECK(m.coverage > 0.3 && m.coverage <= 1.0);
    CHECK(std::fabs(m.fibreAngleDeg) < 2.0);
    const double hLo = (l.rIn - l.rMax) / l.r - 0.01, hHi = (l.rOut + l.rMax) / l.r + 0.01;
    Vec3 mean;
    int texels = 0, hitTexels = 0, alongYarn = 0, unit = 0, inShell = 0;
    for (int i = 0; i < m.along; ++i)
        for (int j = 0; j < m.around; ++j) {
            int x = 0, row = 0;
            atlasPixel(m, i, j, x, row);
            const Vec3 n = bakedNormal(m, x, row);
            mean += n;
            ++texels;
            if (norm(n) > 0.95 && norm(n) < 1.05) ++unit;
            if (bakedGrey(m.specular, m, x, row) < 0.999) continue;  // not every ray met a fibre
            ++hitTexels;
            if (std::fabs(n.x) < 0.1) ++alongYarn;
            const double h = bakedGrey(m.heightMap, m, x, row);
            if (h >= hLo && h <= hHi) ++inShell;
        }
    mean = mean / texels;
    std::printf("    hit texels %d of %d; mean normal (%.3f, %.3f, %.3f)\n", hitTexels, texels, mean.x, mean.y, mean.z);
    CHECK(unit == texels);
    CHECK(hitTexels > 0 && alongYarn >= 0.99 * hitTexels && inShell >= 0.99 * hitTexels);
    CHECK(mean.z >= 0.7 && std::fabs(mean.y) < 0.1 && std::fabs(mean.x) < 0.05);
    // The sign convention, on six columns: a run of consecutive fully-hit texels in theta is a fibre's
    // visible arc (or several overlapping ones); along +theta its normal's y turns from negative to
    // positive. Every run of at least five texels must start negative and end positive.
    int runs = 0, signedRuns = 0;
    for (int i = m.along / 7; i < m.along; i += m.along / 7) {
        int j = 0;
        while (j < m.around) {
            int x = 0, row = 0;
            atlasPixel(m, i, j, x, row);
            if (bakedGrey(m.specular, m, x, row) < 0.999) {
                ++j;
                continue;
            }
            int j1 = j;
            while (j1 + 1 < m.around) {
                atlasPixel(m, i, j1 + 1, x, row);
                if (bakedGrey(m.specular, m, x, row) < 0.999) break;
                ++j1;
            }
            if (j1 - j + 1 >= 5) {
                ++runs;
                int x0 = 0, r0 = 0, x1 = 0, r1 = 0;
                atlasPixel(m, i, j, x0, r0);
                atlasPixel(m, i, j1, x1, r1);
                const double yFirst = bakedNormal(m, x0, r0).y, yLast = bakedNormal(m, x1, r1).y;
                if (yFirst < 0 && yLast > 0) ++signedRuns;
                else if (runs - signedRuns <= 3)
                    std::printf("    column %d run %d..%d: n.y %.3f .. %.3f\n", i, j, j1, yFirst, yLast);
            }
            j = j1 + 1;
        }
    }
    std::printf("    runs of 5+ fibre texels: %d, of which %d turn from -y to +y\n", runs, signedRuns);
    CHECK(runs >= 20 && signedRuns == runs);
}

// The gaps of a plied yarn (the blue swatch's body, thinned to this yarn, at the production texel): a ray
// through a gap meets the ply's core or, between the plies, the floor with its ply's flank, so a gap reads
// as deep and round, never as the smooth tube, and it is shadowed. The defect this guards (Figure 23 b of
// the article, 2026-11-22): a third of the texels kept the tube's flat normal, and the AO was brighter in
// the gaps (0.74) than on the fibres (0.67), so the tube rendered plain. Checks: almost no texel keeps the
// flat normal; the gaps' AO is clearly below the outer fibres' (the upper half of the fibre texels by
// height), and so is the inner fibres' (the shell's transmittance); cores are met; a gap texel's normal
// leans to the side of its nearest ply's centre it lies on (n.y < 0 before the centre in theta, > 0 after).
TEST(yarn_bake_gaps_shadowed) {
    const CurveSet cs = makeSet({straightYarn(10)});
    const SurfaceModel sm = roundSurface(cs);
    YarnBodyParams p;
    p.plies = 3;
    p.fibres = 180;  // the swatch's 400 on its 1.9 mm yarn; this 1.5 mm yarn's shell holds 195
    p.fibreRadius = 0.018;
    p.shellDepth = 0.6;
    p.fibreAngleDeg = 20;
    p.plyAngleDeg = 20;
    p.irregularity = 1;
    p.core = BodyCore::Opaque;
    const YarnBody b = bodyOf(cs, sm, p);
    BakeParams bp;  // the production texel and rays
    BakedYarnMaps m;
    std::string err;
    const bool baked = bakeYarnMaps(cs, sm.yarns[0], b.strands, p, kPartSeed, 0, bp, m, &err, ExecPolicy{true});
    if (!baked) std::printf("    %s\n", err.c_str());
    CHECK(baked);
    if (!baked) return;
    const double r = sm.yarns[0].radius;
    int texels = 0, flat = 0, gaps = 0, cores = 0, leaningTexels = 0, leaning = 0;
    double gapAo = 0;
    std::vector<std::pair<double, double>> fibreTexels;  // (height, AO) of the texels every ray of which met a fibre
    for (int i = 0; i < m.along; ++i) {
        const double s = (i + 0.5) * m.texel;
        double plyTurn[3];
        for (int k = 0; k < p.plies; ++k) {
            const Vec2 c = bodyPlyAxis(p, r, kPartSeed, 0, k, s).centre;
            plyTurn[k] = std::atan2(c.y, c.x) / (2 * kPi);
        }
        for (int j = 0; j < m.around; ++j) {
            int x = 0, row = 0;
            atlasPixel(m, i, j, x, row);
            const Vec3 n = bakedNormal(m, x, row);
            const double coverage = bakedGrey(m.specular, m, x, row), ao = bakedGrey(m.ao, m, x, row), h = bakedGrey(m.heightMap, m, x, row);
            ++texels;
            if (std::fabs(n.x) < 0.02 && std::fabs(n.y) < 0.02 && coverage < 0.5) ++flat;
            if (coverage >= 0.999) {
                fibreTexels.push_back({h, ao});
            } else if (coverage <= 0.001) {
                ++gaps;
                gapAo += ao;
                if (h > 0.52) ++cores;  // a core, not the floor at 0.5 r
                double offset = 1;      // the texel's turn from the nearest ply centre, in (-1/2, 1/2]
                for (int k = 0; k < p.plies; ++k) {
                    double d = (j + 0.5) / m.around - plyTurn[k];
                    d -= std::floor(d + 0.5);
                    if (std::fabs(d) < std::fabs(offset)) offset = d;
                }
                if (std::fabs(offset) > 0.5 / m.around) {  // not at the crest itself
                    ++leaningTexels;
                    if ((offset < 0) == (n.y < 0)) ++leaning;
                }
            }
        }
    }
    gapAo /= std::max(1, gaps);
    std::sort(fibreTexels.begin(), fibreTexels.end());  // by height
    const size_t fibres = fibreTexels.size(), half = fibres / 2;
    double innerAo = 0, outerAo = 0;
    for (size_t k = 0; k < fibres; ++k) (k < half ? innerAo : outerAo) += fibreTexels[k].second;
    innerAo /= static_cast<double>(std::max(size_t{1}, half));
    outerAo /= static_cast<double>(std::max(size_t{1}, fibres - half));
    std::printf("    %d x %d texels: coverage %.3f, cores %.3f; flat normals %.1f %%; AO gaps %.3f (%d), inner fibres %.3f, outer fibres "
                "%.3f (%zu); core texels %d; gap texels leaning to their ply's side %d of %d\n",
                m.along, m.around, m.coverage, m.cores, 100.0 * flat / texels, gapAo, gaps, innerAo, outerAo, fibres, cores, leaning,
                leaningTexels);
    CHECK(flat < 0.02 * texels);
    CHECK(gaps > 0.05 * texels && fibres > 0.05 * texels);
    CHECK(gapAo < outerAo - 0.05 && innerAo < outerAo - 0.05);
    CHECK(cores > 0.02 * texels && leaning >= 0.9 * leaningTexels);
}

// The atlas layout: strips of atlasWidth texels along s stacked in v; chart texel <-> pixel is a
// bijection on the used texels; a chart taller than the atlas is refused; atlasStrips agrees.
TEST(yarn_bake_atlas_layout) {
    BakeParams p;
    p.texel = 0.0625;
    p.atlasWidth = 256;
    BakedYarnMaps m;
    std::string err;
    CHECK(planAtlas(30, 0.75, p, m, &err));
    std::printf("    along %d around %d strips %d: %d x %d\n", m.along, m.around, m.strips, m.width, m.height);
    CHECK(m.along == 480 && m.around == 75 && m.strips == 2 && m.width == 256 && m.height == 150);
    std::vector<int> seen(static_cast<size_t>(m.width) * m.height, 0);
    bool roundTrip = true;
    for (int i = 0; i < m.along; ++i)
        for (int j = 0; j < m.around; ++j) {
            int x = 0, row = 0, i2 = -1, j2 = -1;
            atlasPixel(m, i, j, x, row);
            if (x < 0 || x >= m.width || row < 0 || row >= m.height) {
                roundTrip = false;
                continue;
            }
            ++seen[static_cast<size_t>(row) * m.width + x];
            roundTrip = roundTrip && chartTexel(m, x, row, i2, j2) && i2 == i && j2 == j;
        }
    CHECK(roundTrip);
    int once = 0, unused = 0, tailRefused = 0;
    for (int row = 0; row < m.height; ++row)
        for (int x = 0; x < m.width; ++x) {
            const int k = seen[static_cast<size_t>(row) * m.width + x];
            if (k == 1) ++once;
            if (k == 0) {
                ++unused;
                int i2 = 0, j2 = 0;
                if (!chartTexel(m, x, row, i2, j2)) ++tailRefused;
            }
        }
    CHECK(once == m.along * m.around && unused == m.width * m.height - m.along * m.around && tailRefused == unused);
    CHECK(atlasStrips(30, p) == 2 && atlasStrips(16, p) == 1 && atlasStrips(16.0625, p) == 2);
    p.atlasWidth = 64;
    CHECK(!planAtlas(30, 0.75, p, m, &err) && err.find("exceed") != std::string::npos);
}

// The tube built for the atlas (uvRepeat = one strip of yarn, uvWrap = 1) and relabelled: every vertex's
// (u, v) names the chart texel at its (s, theta) on its own yarn's atlas, and no triangle straddles two
// strips. Two yarns of different lengths, so the strips differ per yarn.
TEST(yarn_bake_uv_matches_tube) {
    const CurveSet cs = makeSet({straightYarn(30), straightYarn(20, 3, 1)});
    const SurfaceModel sm = roundSurface(cs, 0.1);
    BakeParams p;
    p.texel = 0.0625;
    p.atlasWidth = 256;
    TubeParams tp;
    tp.ringSpacing = 0.2;
    tp.sides = 16;
    tp.uvRepeat = p.texel * p.atlasWidth;
    tp.uvWrap = 1;
    TubeMesh tube = buildTubeMesh(cs, sm, tp);
    relabelAtlasUv(tube, p);
    CHECK(tube.yarnFirstVertex.size() == 3);
    int wrong = 0, vertices = 0;
    for (size_t y = 0; y + 1 < tube.yarnFirstVertex.size(); ++y) {
        const uint32_t first = tube.yarnFirstVertex[y], end = tube.yarnFirstVertex[y + 1];
        double length = 0;
        for (uint32_t v = first; v < end; ++v) length = std::max(length, tube.s[v]);
        BakedYarnMaps m;
        std::string err;
        CHECK(planAtlas(length, sm.yarns[y].radius, p, m, &err));
        CHECK(m.strips == (y == 0 ? 2 : 2) && atlasStrips(length, p) == m.strips);
        for (uint32_t v = first; v < end; ++v) {
            ++vertices;
            const double u = tube.uv[2 * v], vv = tube.uv[2 * v + 1], turn = tube.theta[v] / (2 * kPi);
            const int strip = static_cast<int>(std::lround(vv * m.strips - turn));
            const double iReal = strip * m.width + u * m.width, iWant = tube.s[v] / p.texel;
            // The pixel the renderer samples at (u, v): its chart texel. The seam vertex (theta = 2 pi) sits
            // on the boundary row between its strip and the next: its texel is read a hair inside the strip.
            // Likewise the last ring (s = length) sits on the edge of the first unused tail texel.
            const double uInside = tube.s[v] >= length - 1e-9 ? u - 0.5 / m.width : u;
            const int x = std::clamp(static_cast<int>(uInside * m.width), 0, m.width - 1);
            const double vInside = turn > 0.999 ? vv - 0.5 / m.height : vv;
            const int row = std::clamp(m.height - 1 - static_cast<int>(vInside * m.height), 0, m.height - 1);
            int i = -1, j = -1;
            const bool used = chartTexel(m, x, row, i, j);
            const double jWant = turn * m.around;
            const double jDiff = std::min(std::fabs(j - jWant), m.around - std::fabs(j - jWant));
            if (strip < 0 || strip >= m.strips || std::fabs(iReal - iWant) > 1e-3 || !used || std::fabs(i - iWant) > 1.001 || jDiff > 1.001) {
                if (wrong < 5)
                    std::printf("    yarn %zu vertex %u: s %.4f theta %.4f -> u %.5f v %.5f strip %d i %.2f (want %.2f) texel (%d, %d) want j %.2f\n",
                                y, v, tube.s[v], tube.theta[v], u, vv, strip, iReal, iWant, i, j, jWant);
                ++wrong;
            }
        }
    }
    std::printf("    %d vertices checked, %d wrong\n", vertices, wrong);
    CHECK(wrong == 0);
    int straddling = 0;
    for (size_t t = 0; t + 2 < tube.tri.size(); t += 3) {
        int strips[3];
        for (int k = 0; k < 3; ++k) {
            const uint32_t v = tube.tri[t + k];
            const double length = tube.s[tube.yarnFirstVertex[tube.yarn[v] + 1] - 1];
            strips[k] = static_cast<int>(std::lround(tube.uv[2 * v + 1] * atlasStrips(length, p) - tube.theta[v] / (2 * kPi)));
        }
        if (strips[0] != strips[1] || strips[1] != strips[2]) ++straddling;
    }
    CHECK(straddling == 0);
}

// The material maps a package derives from the baked ones: the cavity map darkens where the ambient occlusion
// is low (sRGB-encoded), the roughness map starts at the smooth fibre's roughness and grows with Toksvig's
// spread, never above 1; unoccluded, unspread texels keep the material's own values.
TEST(yarn_bake_material_maps) {
    const std::vector<uint8_t> ao = {0, 0, 0, 128, 128, 128, 255, 255, 255};
    const std::vector<uint8_t> cavity = cavityMapOf(ao, 0.6);
    CHECK(cavity.size() == ao.size());
    const double dark = srgbDecode(cavity[0] / 255.0), mid = srgbDecode(cavity[3] / 255.0), open = srgbDecode(cavity[6] / 255.0);
    std::printf("    cavity 0.6: ao 0 -> %.3f, ao 0.5 -> %.3f, ao 1 -> %.3f\n", dark, mid, open);
    CHECK(std::fabs(dark - 0.4) < 0.01 && std::fabs(mid - (1 - 0.6 * (1 - 128 / 255.0))) < 0.01 && std::fabs(open - 1.0) < 0.005);
    const std::vector<uint8_t> spread = {0, 0, 0, 64, 64, 64, 255, 255, 255};
    const std::vector<uint8_t> rough = roughnessMapOf(spread, 0.3);
    const double smooth = rough[0] / 255.0, some = rough[3] / 255.0, full = rough[6] / 255.0;
    std::printf("    roughness floor 0.3: spread 0 -> %.3f, 0.25 -> %.3f, 1 -> %.3f\n", smooth, some, full);
    CHECK(std::fabs(smooth - 0.3) < 0.005 && std::fabs(some - std::sqrt(0.09 + 64 / 255.0)) < 0.005 && full == 1.0);
}
