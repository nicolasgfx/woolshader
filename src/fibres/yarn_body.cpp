// The yarn body made of structural fibres (prototype); see yarn_body.h.
#include "fibres/yarn_body.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <unordered_map>

#include "knitting/tube/plies.h"
#include "knitting/tube/yarn_surface.h"
#include "model/cancel.h"
#include "model/hash.h"
#include "model/yarn_path.h"

namespace knit {
namespace {

// Engineering choices; paper values are named where they come from one.
constexpr double kGridStep = 0.05;          // [mm] the rate-integration grid along a yarn
constexpr double kSoftDeviation = 0.25;     // x r_f: midpoint deviation target (the renderer's strand-format decision)
constexpr double kSoftTurnDeg = 6;          // turn per segment: highlight facets below the hair lobe's roughness
constexpr double kMinSegOverR = 4.4;        // hard rule 4 x the larger end radius (strands_file.h taperRules), 10 % margin for float32 metres
constexpr double kMaxSegOverR = 150;        // hard rule is 200 x the smaller radius; margin for the merged last segment
constexpr double kMaxRadiusStep = 1.45;     // neighbouring radii: hard rule 1.5 (strands_file.h taperRules), with margin
constexpr double kPlyFill = 0.95;           // r_p = 0.95 r_p0, as strand.ply_fill
constexpr double kMinSpacingOverR = 2.1;    // dart throwing: start centres >= 2.1 r_f apart (an even spread, no guarantee)
constexpr int kPlacementAttempts = 128;
constexpr double kCapacityFraction = 0.35;  // well below random sequential adsorption's jamming (about 0.547): placeable
constexpr double kYarnNoiseLength = 10, kGroupNoiseLength = 4, kRadiusNoiseLength = 20;  // [mm]
constexpr int kGroupSize = 8;               // fibres that share a group's twist and depth noise
constexpr double kSurfaceShare = 0.3;       // fibres that barely migrate (irregular body)
constexpr double kStapleMedian = 60, kStapleSigma = 0.3;  // [mm] wool staple (worsted: 60-100 mm)
constexpr double kEndZone = 1.5;            // [mm] a piece's end dives into the shell and tapers over this
constexpr double kEndRadius = 0.5;          // the end's radius factor
constexpr double kCrimpWavelength = 1.5;    // [mm] fine merino, about 6.7 crimps / cm (AWTA)
constexpr double kMigrationMinAngleDeg = 5; // the migration period of a fibre with no twist
constexpr double kCoreRing = 0.3;           // [mm]
constexpr int kCoreSides = 8;
// Thinnest fibre: 10 um. A piece's end (kEndRadius) at the most irregular radius stays above the strand
// files' 1e-6 m radius floor (strands_file.h).
constexpr double kMinFibreRadius = 0.005;  // [mm]
constexpr double kOverlapWindow = 10;       // [mm] of yarn per overlap window
constexpr double kSectionSmoothing = 0.1;   // [mm] sigma of the section map's smoothing along the yarn (SectionTable)

// The draws of one fibre (hash01 channels): each a pure function of (ply seed, slot, piece).
enum Draw { kStartAngle, kStartDepth, kAttempt, kMigrationKind, kMigrationAmplitude, kMigrationPhase, kMigrationRatio,
            kRadiusNormal0, kRadiusNormal1, kCrimpPhase, kCrimpDirection, kCrimpLength, kStapleOffset, kStapleNormal0,
            kStapleNormal1, kPlyPhase };

double degToRad(double d) { return d * kPi / 180; }
double lerp(double a, double b, double t) { return a + (b - a) * t; }
double smoothstep(double t) {
    t = std::clamp(t, 0.0, 1.0);
    return t * t * (3 - 2 * t);
}

// Smooth value noise in [-1, 1] along x: hashed values on the integers, joined by a cubic
// smoothstep. A pure function of (seed, x): the same wherever and however densely it is sampled.
double smoothNoise(uint64_t seed, double x) {
    const double f = std::floor(x);
    const long long i = static_cast<long long>(f);
    const double v0 = 2 * hash01(seed, static_cast<int>(i & 0x3fffffff), static_cast<int>(i >> 30), 0) - 1;
    const double v1 = 2 * hash01(seed, static_cast<int>((i + 1) & 0x3fffffff), static_cast<int>((i + 1) >> 30), 0) - 1;
    return v0 + (v1 - v0) * smoothstep(x - f);
}

// The seed of a yarn: the part, the body's seed, the yarn's index (its row: stable for a frozen knit).
uint64_t yarnSeedOf(uint64_t partSeed, uint64_t bodySeed, int yarn) {
    const uint64_t base = splitmix64(partSeed ^ splitmix64(bodySeed));
    return splitmix64(base ^ (0x9E3779B97F4A7C15ull * static_cast<uint64_t>(yarn + 1)));
}

// A linear table on the fixed grid z_i = i * kGridStep (cumulative trapezoid of a rate).
struct GridTable {
    std::vector<double> v;
    double at(double z) const {
        const double x = std::clamp(z / kGridStep, 0.0, static_cast<double>(v.size() - 1));
        const size_t i = std::min(static_cast<size_t>(x), v.size() - 2);
        return v[i] + (v[i + 1] - v[i]) * (x - static_cast<double>(i));
    }
};
size_t gridPoints(double length) { return static_cast<size_t>(std::ceil(length / kGridStep)) + 2; }

} // namespace

// The section table (yarn_body.h): the rings' maps smoothed by kSectionSmoothing.
std::array<double, 3> SectionTable::at(double s) const {
    const double x = std::clamp(s / ds, 0.0, static_cast<double>(m.size() - 1));
    const size_t i = std::min(static_cast<size_t>(x), m.size() - 2);
    const double t = x - static_cast<double>(i);
    return {lerp(m[i][0], m[i + 1][0], t), lerp(m[i][1], m[i + 1][1], t), lerp(m[i][2], m[i + 1][2], t)};
}

SectionTable sectionTable(const YarnSurface& ys) {
    const size_t n = ys.rings.size();
    std::vector<std::array<double, 3>> raw(n);
    for (size_t k = 0; k < n; ++k) {
        const Ring& g = ys.rings[k];
        const double cp = std::cos(g.phi), sp = std::sin(g.phi), A = g.a / ys.radius, B = g.b / ys.radius;
        raw[k] = {A * cp * cp + B * sp * sp, (A - B) * cp * sp, A * sp * sp + B * cp * cp};
    }
    SectionTable t;
    t.ds = ys.ds;
    t.m.resize(std::max<size_t>(n, 2), raw.empty() ? std::array<double, 3>{1, 0, 1} : raw.back());
    const double sigma = kSectionSmoothing / ys.ds;  // in rings
    const int half = static_cast<int>(std::ceil(3 * sigma));
    for (size_t k = 0; k < n; ++k) {
        std::array<double, 3> sum{0, 0, 0};
        double weights = 0;
        for (int j = -half; j <= half; ++j) {
            const long long i = static_cast<long long>(k) + j;
            if (i < 0 || i >= static_cast<long long>(n)) continue;
            const double w = sigma > 0 ? std::exp(-0.5 * (j / sigma) * (j / sigma)) : (j == 0 ? 1.0 : 0.0);
            for (int e = 0; e < 3; ++e) sum[e] += w * raw[static_cast<size_t>(i)][e];
            weights += w;
        }
        for (int e = 0; e < 3; ++e) t.m[k][e] = sum[e] / weights;
    }
    return t;
}

Vec3 sectionOffset(const Ring& g, const std::array<double, 3>& m, double x, double y) {
    return g.n1 * (m[0] * x + m[1] * y) + g.n2 * (m[1] * x + m[2] * y);
}

namespace {

// The smaller eigenvalue of a symmetric 2 x 2 (m11, m12, m22): the section's shortest semi-axis over r.
double smallerAxis(const std::array<double, 3>& m) {
    const double mean = 0.5 * (m[0] + m[2]), half = 0.5 * (m[0] - m[2]);
    return mean - std::sqrt(half * half + m[1] * m[1]);
}

// The ply orbit along one yarn: omega(z) = omega0 (1 + 0.1 iota n(z)), integrated to theta(z).
struct YarnTwist {
    double R = 0, omega0 = 0, iota = 0;
    uint64_t seed = 0;
    GridTable theta;
    double omega(double z) const { return omega0 * (1 + 0.1 * iota * smoothNoise(seed, z / kYarnNoiseLength)); }
};

YarnTwist yarnTwist(const BodyLayout& l, const YarnBodyParams& p, uint64_t yarnSeed, double length) {
    YarnTwist t;
    t.R = l.R;
    t.iota = p.irregularity;
    t.seed = splitmix64(yarnSeed ^ 0x7177u);
    // One ply does not orbit (R = 0); the rate is the helix's: tan(gamma) = R omega.
    t.omega0 = p.plies > 1 && l.R > 0 ? std::tan(degToRad(p.plyAngleDeg)) / l.R : 0;  // handedness +1 (Z ply)
    t.theta.v.resize(gridPoints(length));
    t.theta.v[0] = 2 * kPi * hash01(yarnSeed, 0, 0, kPlyPhase);
    for (size_t i = 1; i < t.theta.v.size(); ++i) {
        const double z = static_cast<double>(i) * kGridStep;
        t.theta.v[i] = t.theta.v[i - 1] + 0.5 * kGridStep * (t.omega(z - kGridStep) + t.omega(z));
    }
    return t;
}

// The ply's local frame in nominal space at z: C its centre (x, y, 0),
// u radial, v = t x u, tau the helix frame's twist per ply arc length, cosGamma = dz / ds_ply.
struct PlyFrame {
    Vec3 C, u, v;
    double tau = 0, cosGamma = 1;
};

PlyFrame plyFrame(const YarnTwist& tw, int ply, int plies, double z) {
    const double w = tw.omega(z);
    const double th = tw.theta.at(z) + 2 * kPi * ply / plies;
    const Vec3 er{std::cos(th), std::sin(th), 0}, et{-std::sin(th), std::cos(th), 0}, ez{0, 0, 1};
    const double g = std::atan(tw.R * w);
    PlyFrame f;
    f.C = er * tw.R;
    f.u = er;
    f.v = et * std::cos(g) - ez * std::sin(g);
    f.cosGamma = std::cos(g);
    f.tau = w * f.cosGamma * f.cosGamma;
    return f;
}

// One fibre slot: its start in the section and its own draws.
struct Slot {
    int index = 0;
    double psi0 = 0, m0 = 0;      // start angle about the ply axis, start depth (0 outer .. 1 inner)
    double amplitude = 0.5, phase = 0, period = 0;   // migration
    uint64_t groupSeed = 0;
    double radiusFactor = 1;      // exp(0.2 iota g)
    double crimpPhase = 0, crimpDirection = 0, crimpLength = kCrimpWavelength;
};

// Where a ply's fibres start (z = 0): dart throwing in the annulus [rIn, rOut] (rIn = 0: the disc),
// centres at least kMinSpacingOverR r_f apart, at most kPlacementAttempts tries each. False when a
// fibre finds no place (never fewer fibres than asked).
bool placeStarts(uint64_t seed, int count, double rIn, double rOut, double rf, std::vector<Vec2>& out) {
    out.clear();
    const double minD2 = (kMinSpacingOverR * rf) * (kMinSpacingOverR * rf);
    for (int k = 0; k < count; ++k) {
        bool placed = false;
        for (int a = 0; a < kPlacementAttempts && !placed; ++a) {
            const double rr = std::sqrt(lerp(rIn * rIn, rOut * rOut, hash01(seed, k, a, kStartDepth)));
            const double an = 2 * kPi * hash01(seed, k, a, kStartAngle);
            const Vec2 q{rr * std::cos(an), rr * std::sin(an)};
            placed = std::none_of(out.begin(), out.end(), [&](const Vec2& o) {
                const double dx = o.x - q.x, dy = o.y - q.y;
                return dx * dx + dy * dy < minD2;
            });
            if (placed) out.push_back(q);
        }
        if (!placed) return false;
    }
    return true;
}

// The slots of one ply zone from their starts; iota scales every irregular draw:
// at iota = 0 the amplitude is 0.5, the period ratio 1.5 and the radius r_f for every fibre.
std::vector<Slot> makeSlots(uint64_t seed, const std::vector<Vec2>& starts, double rIn, double rOut, double period, double iota) {
    std::vector<Slot> slots(starts.size());
    const double irregular = std::min(iota, 1.0);
    for (size_t k = 0; k < starts.size(); ++k) {
        Slot& s = slots[k];
        const int i = static_cast<int>(k);
        s.index = i;
        s.psi0 = std::atan2(starts[k].y, starts[k].x);
        const double rho = std::hypot(starts[k].x, starts[k].y);
        s.m0 = rOut > rIn ? std::clamp((rOut - rho) / (rOut - rIn), 0.0, 1.0) : 0;
        const bool surface = hash01(seed, i, 0, kMigrationKind) < kSurfaceShare;
        const double u = hash01(seed, i, 0, kMigrationAmplitude);
        s.amplitude = lerp(0.5, surface ? 0.2 * u : lerp(0.4, 1.0, u), irregular);                // iota
        s.phase = 2 * kPi * hash01(seed, i, 0, kMigrationPhase);
        s.period = period / lerp(1.5, lerp(1.1, 1.9, hash01(seed, i, 0, kMigrationRatio)), irregular);  // iota; Zhao Table 2: s 1.1-1.9
        s.radiusFactor = std::exp(0.2 * iota * hashNormal(seed, i, 0, kRadiusNormal0, kRadiusNormal1, 2.0));  // iota
        s.crimpPhase = 2 * kPi * hash01(seed, i, 0, kCrimpPhase);
        s.crimpDirection = 2 * kPi * hash01(seed, i, 0, kCrimpDirection);
        s.crimpLength = kCrimpWavelength * (1 + 0.2 * irregular * (2 * hash01(seed, i, 0, kCrimpLength) - 1));  // iota
    }
    // Groups of neighbours around the ply share a smooth noise on twist rate and depth.
    std::vector<int> order(slots.size());
    for (size_t k = 0; k < order.size(); ++k) order[k] = static_cast<int>(k);
    std::sort(order.begin(), order.end(), [&](int a, int b) { return slots[a].psi0 < slots[b].psi0; });
    for (size_t rank = 0; rank < order.size(); ++rank)
        slots[order[rank]].groupSeed = splitmix64(seed ^ (0x51ull + rank / kGroupSize));
    return slots;
}

// A slot's sequence of pieces along [begin, end]: continuous at iota = 0, else log-normal staple lengths
// from a random offset. Each piece knows whether its ends are the yarn's (cut) or its own (dive).
struct Piece {
    double z0 = 0, z1 = 0;
    bool ownStart = false, ownEnd = false;
};

std::vector<Piece> slotPieces(uint64_t seed, int slot, double begin, double end, double iota) {
    if (iota <= 0) return {{begin, end, false, false}};
    std::vector<Piece> pieces;
    double z = begin - kStapleMedian * hash01(seed, slot, 0, kStapleOffset);
    for (int k = 1; z < end; ++k) {
        const double len = kStapleMedian * std::exp(kStapleSigma * hashNormal(seed, slot, k, kStapleNormal0, kStapleNormal1, 3.0));
        const Piece pc{std::max(z, begin), std::min(z + len, end), z > begin, z + len < end};
        if (pc.z1 > pc.z0) pieces.push_back(pc);
        z += len;
    }
    return pieces;
}

// What one ply zone (shell or interior) needs to evaluate its fibres.
struct ZoneContext {
    const CurveSet* cs = nullptr;
    const YarnSurface* ys = nullptr;
    const YarnTwist* twist = nullptr;
    const SectionTable* section = nullptr;
    const BodyLayout* layout = nullptr;
    const YarnBodyParams* p = nullptr;
    int ply = 0;
    double trim = 0;                // plied fibres end r_p inside the yarn's ends: their axial offset stays on the yarn
    double rOuter = 0, depth = 1;   // centreline radius at m = 0, and the depth span (rho = rOuter (1 - depth m))
    double omegaFibre = 0;          // target rotation about the ply axis per ply arc length
};

struct FibrePoint {
    Vec3 X;
    double radius = 0, s = 0, protrusion = 0;
    double envelope = 0;   // the centreline's radius in its true ring ellipse (1 = on it)
};

// A point's radius in the true ellipse of ring g: hypot(along e_a / a, along e_b / b).
double ellipseRadius(const Ring& g, const Vec3& X) {
    const double cp = std::cos(g.phi), sp = std::sin(g.phi);
    const Vec3 d = X - g.c, ea = g.n1 * cp + g.n2 * sp, eb = g.n2 * cp - g.n1 * sp;
    return std::hypot(dot(d, ea) / g.a, dot(d, eb) / g.b);
}

// A fibre of slot `s` at yarn arc length z, with psi(z) from its table and the piece's end factor.
FibrePoint fibreAt(const ZoneContext& c, const Slot& s, const GridTable& psi, const Piece& pc, double z) {
    const double iota = c.p->irregularity;
    const PlyFrame f = plyFrame(*c.twist, c.ply, c.p->plies, z);
    // Migration (Zhao et al. Eq. 3 within the zone) and the group's depth noise; a piece's own ends dive.
    double m = s.m0 + 0.5 * s.amplitude * (std::cos(2 * kPi * z / s.period + s.phase) - std::cos(s.phase)) +
               0.15 * iota * (smoothNoise(s.groupSeed, z / kGroupNoiseLength) - smoothNoise(s.groupSeed, 0));
    m = std::clamp(m, 0.0, 1.0);
    double end = 1;
    if (pc.ownStart) end = std::min(end, smoothstep((z - pc.z0) / kEndZone));
    if (pc.ownEnd) end = std::min(end, smoothstep((pc.z1 - z) / kEndZone));
    m = 1 - (1 - m) * end;
    const double rho = c.rOuter * (1 - c.depth * m);
    const double ps = psi.at(z);
    const double crimp = iota * c.p->fibreRadius * std::sin(2 * kPi * z / s.crimpLength + s.crimpPhase);  // iota
    const Vec3 Q = f.C + (f.u * std::cos(ps) + f.v * std::sin(ps)) * rho +
                   (f.u * std::cos(s.crimpDirection) + f.v * std::sin(s.crimpDirection)) * crimp;
    FibrePoint out;
    out.s = std::clamp(z + Q.z, 0.0, c.ys->length);
    const Ring g = ringAt(*c.cs, *c.ys, out.s);
    const std::array<double, 3> M = c.section->at(out.s);
    out.X = g.c + sectionOffset(g, M, Q.x, Q.y);
    out.radius = c.p->fibreRadius * s.radiusFactor * (1 + 0.1 * iota * smoothNoise(s.groupSeed ^ 0xabcu, z / kRadiusNoiseLength)) *  // iota
                 lerp(kEndRadius, 1.0, end);
    // The centreline at normalized radius t lies in t E (E the smoothed ellipse), (1 - t) r_min(E) inside E.
    const double t = std::hypot(Q.x, Q.y) / c.layout->r;
    out.protrusion = out.radius - (1 - t) * smallerAxis(M) * c.layout->r;
    out.envelope = ellipseRadius(g, out.X);
    return out;
}

// psi(z) of a slot over the whole yarn: d psi / dz = (Omega (1 + 0.15 iota n_g) - tau) / cos(gamma),
// the target rotation about the ply axis less the helix frame's own twist.
GridTable slotPsi(const ZoneContext& c, const Slot& s) {
    GridTable t;
    t.v.resize(gridPoints(c.ys->length));
    const auto rate = [&](double z) {
        const PlyFrame f = plyFrame(*c.twist, c.ply, c.p->plies, z);
        const double omega = c.omegaFibre * (1 + 0.15 * c.p->irregularity * smoothNoise(s.groupSeed ^ 0x5eedu, z / kGroupNoiseLength));
        return (omega - f.tau) / f.cosGamma;
    };
    t.v[0] = s.psi0;
    double previous = rate(0);
    for (size_t i = 1; i < t.v.size(); ++i) {
        const double next = rate(static_cast<double>(i) * kGridStep);
        t.v[i] = t.v[i - 1] + 0.5 * kGridStep * (previous + next);
        previous = next;
    }
    return t;
}

double angleDeg(const Vec3& a, const Vec3& b) {
    const double la = norm(a), lb = norm(b);
    if (la <= 0 || lb <= 0) return 0;
    return std::acos(std::clamp(dot(a, b) / (la * lb), -1.0, 1.0)) * 180 / kPi;
}

// The result of one zone (one ply's shell or interior along one yarn), merged in order afterwards.
struct ZoneResult {
    BodyStrands strands;
    long long pieces = 0, dropped = 0, softTurn = 0;
    std::array<long long, 1800> turns{};   // 0.1 degree bins
    double fibreLength = 0, maxProtrusion = 0, maxEnvelope = 0;
};

// The next vertex of a fibre from `from` at z. First the shortest admissible step:
// grown from the proposal h until the chord reaches the hard minimum (4.4 x the larger radius) or the
// piece ends. Then halved while a soft target is missed (midpoint deviation, turn) or the chord is above
// the hard maximum, as long as the half step is still admissible; a step that cannot be halved further
// is kept and counted (soft). The accepted vertex meets the hard minimum unless it ends the piece.
struct Step {
    FibrePoint e;
    double h = 0;
    bool soft = false;
};

Step nextStep(const ZoneContext& c, const Slot& s, const GridTable& psi, const Piece& pc, const FibrePoint& from, const Vec3& prevDir,
              double z, double h) {
    const double rest = pc.z1 - z;
    const auto minLength = [&](const FibrePoint& e) { return kMinSegOverR * std::max(e.radius, from.radius); };
    Step st;
    st.h = std::min(h, rest);
    st.e = fibreAt(c, s, psi, pc, z + st.h);
    for (int it = 0; it < 40 && norm(st.e.X - from.X) < minLength(st.e) && st.h < rest; ++it) {
        const double chord = norm(st.e.X - from.X);
        st.h = std::min(st.h * std::max(1.25, 1.05 * minLength(st.e) / std::max(chord, 1e-9)), rest);
        st.e = fibreAt(c, s, psi, pc, z + st.h);
    }
    for (int it = 0; it < 40; ++it) {
        const Vec3 d = st.e.X - from.X;
        const FibrePoint mid = fibreAt(c, s, psi, pc, z + 0.5 * st.h);
        const bool missed = norm(d) > kMaxSegOverR * std::min(st.e.radius, from.radius) ||
                            norm(mid.X - (st.e.X + from.X) * 0.5) > kSoftDeviation * c.p->fibreRadius ||
                            (norm(prevDir) > 0 && angleDeg(prevDir, d) > kSoftTurnDeg) ||
                            std::max(st.e.radius, from.radius) > kMaxRadiusStep * std::min(st.e.radius, from.radius);
        if (!missed) return st;
        if (norm(mid.X - from.X) < minLength(mid)) {  // the half step would break the hard minimum: keep, count
            st.soft = true;
            return st;
        }
        st.h *= 0.5;
        st.e = mid;
    }
    st.soft = true;
    return st;
}

// One piece as a polyline: nextStep from start to end. A segment still below the hard minimum (a piece's
// last, or a step that could not grow) is merged into its neighbour; a piece too short for two vertices
// is dropped.
void samplePiece(const ZoneContext& c, const Slot& s, const GridTable& psi, const Piece& pc, int yarn, ZoneResult& out) {
    std::vector<FibrePoint> pts{fibreAt(c, s, psi, pc, pc.z0)};
    double z = pc.z0, h = 0.1;
    Vec3 prevDir{0, 0, 0};
    while (z < pc.z1 - 1e-12) {
        const bool inEnd = (pc.ownStart && z - pc.z0 < 2 * kEndZone) || (pc.ownEnd && pc.z1 - z < 2 * kEndZone);
        const Step st = nextStep(c, s, psi, pc, pts.back(), prevDir, z, std::min(h, inEnd ? kEndZone / 4 : 1.0));
        out.softTurn += st.soft;
        const Vec3 d = st.e.X - pts.back().X;
        if (norm(prevDir) > 0) ++out.turns[std::min<size_t>(static_cast<size_t>(angleDeg(prevDir, d) * 10), out.turns.size() - 1)];
        prevDir = d;
        pts.push_back(st.e);
        z += st.h;
        h = 1.5 * st.h;
    }
    const auto tooShort = [&](size_t a, size_t b) {
        return norm(pts[b].X - pts[a].X) < kMinSegOverR * std::max(pts[a].radius, pts[b].radius);
    };
    for (size_t i = 1; i + 1 < pts.size();)  // an inner short segment: its end vertex goes
        if (tooShort(i - 1, i)) pts.erase(pts.begin() + static_cast<std::ptrdiff_t>(i));
        else ++i;
    while (pts.size() >= 3 && tooShort(pts.size() - 2, pts.size() - 1)) pts.erase(pts.end() - 2);
    if (pts.size() < 2 || tooShort(0, 1)) {
        ++out.dropped;
        return;
    }
    ++out.pieces;
    // Neighbouring radii within kMaxRadiusStep (the strand rules allow 50 %): only ever narrows a vertex.
    for (size_t i = 1; i < pts.size(); ++i) pts[i].radius = std::min(pts[i].radius, kMaxRadiusStep * pts[i - 1].radius);
    for (size_t i = pts.size() - 1; i-- > 0;) pts[i].radius = std::min(pts[i].radius, kMaxRadiusStep * pts[i + 1].radius);
    BodyStrands& st = out.strands;
    for (const FibrePoint& q : pts) {
        st.points.push_back(q.X);
        st.radius.push_back(q.radius);
        st.s.push_back(q.s);
        out.maxProtrusion = std::max(out.maxProtrusion, q.protrusion);
        out.maxEnvelope = std::max(out.maxEnvelope, q.envelope);
    }
    for (size_t i = 1; i < pts.size(); ++i) out.fibreLength += norm(pts[i].X - pts[i - 1].X);
    st.first.push_back(st.points.size());
    st.yarnOf.push_back(yarn);
}

// Every slot of one zone along one yarn.
void buildZone(const ZoneContext& c, const std::vector<Slot>& slots, uint64_t seed, int yarn, ZoneResult& out) {
    for (const Slot& s : slots) {
        const GridTable psi = slotPsi(c, s);
        for (const Piece& pc : slotPieces(seed, s.index, c.trim, c.ys->length - c.trim, c.p->irregularity)) samplePiece(c, s, psi, pc, yarn, out);
    }
}

// The core proxy of one ply: a closed tube at rCore (nominal), mapped like the fibres.
void addCoreTube(const ZoneContext& c, double rCore, TriMesh& mesh) {
    const double span = c.ys->length - 2 * c.trim;
    const int rings = std::max(2, static_cast<int>(std::ceil(span / kCoreRing)) + 1);
    const int base = static_cast<int>(mesh.pos.size());
    const auto at = [&](double z, double rr, double alpha) {
        const PlyFrame f = plyFrame(*c.twist, c.ply, c.p->plies, z);
        const Vec3 Q = f.C + (f.u * std::cos(alpha) + f.v * std::sin(alpha)) * rr;
        const double s = std::clamp(z + Q.z, 0.0, c.ys->length);
        return ringAt(*c.cs, *c.ys, s).c + sectionOffset(ringAt(*c.cs, *c.ys, s), c.section->at(s), Q.x, Q.y);
    };
    for (int i = 0; i < rings; ++i) {
        const double z = c.trim + span * i / (rings - 1);
        for (int k = 0; k < kCoreSides; ++k) mesh.pos.push_back(at(z, rCore, 2 * kPi * k / kCoreSides));
    }
    const int capStart = static_cast<int>(mesh.pos.size());
    mesh.pos.push_back(at(c.trim, 0, 0));
    mesh.pos.push_back(at(c.ys->length - c.trim, 0, 0));
    const auto vertex = [&](int ring, int k) { return base + ring * kCoreSides + (k % kCoreSides); };
    for (int i = 0; i + 1 < rings; ++i)
        for (int k = 0; k < kCoreSides; ++k) {
            // Counter-clockwise seen from outside: (u, v, t) is right-handed, alpha turns from u to v.
            mesh.tri.push_back({vertex(i, k), vertex(i + 1, k + 1), vertex(i + 1, k)});
            mesh.tri.push_back({vertex(i, k), vertex(i, k + 1), vertex(i + 1, k + 1)});
        }
    for (int k = 0; k < kCoreSides; ++k) {
        mesh.tri.push_back({capStart, vertex(0, k + 1), vertex(0, k)});
        mesh.tri.push_back({capStart + 1, vertex(rings - 1, k), vertex(rings - 1, k + 1)});
    }
}

} // namespace

BodyLayout bodyLayout(const YarnBodyParams& p, double yarnRadius) {
    BodyLayout l;
    l.r = yarnRadius;
    PlyGeometry g;
    g.count = p.plies;
    g.fill = kPlyFill;
    const PlyLayout pl = plyLayout(g, yarnRadius);
    l.R = pl.R;
    l.rp = pl.rp;
    const double iota = p.irregularity;
    l.rMax = p.fibreRadius * std::exp(0.4 * iota) * (1 + 0.1 * iota);
    l.crimp = iota * p.fibreRadius;
    l.rOut = l.rp - l.rMax - l.crimp;
    l.rIn = l.rOut * (1 - p.shellDepth);
    l.rCore = l.rIn - l.rMax;
    const double fibreArea = kPi * p.fibreRadius * p.fibreRadius * p.fibres / std::max(1, p.plies);
    const double shell = kPi * (l.rOut * l.rOut - l.rIn * l.rIn);
    l.phiShell = shell > 0 ? fibreArea / shell : 0;
    l.phiPly = fibreArea / (kPi * l.rp * l.rp);
    return l;
}

void appendStrands(BodyStrands& to, const BodyStrands& from) {
    const uint64_t offset = to.first.back();
    to.points.insert(to.points.end(), from.points.begin(), from.points.end());
    to.radius.insert(to.radius.end(), from.radius.begin(), from.radius.end());
    to.s.insert(to.s.end(), from.s.begin(), from.s.end());
    for (size_t i = 1; i < from.first.size(); ++i) to.first.push_back(offset + from.first[i]);
    to.yarnOf.insert(to.yarnOf.end(), from.yarnOf.begin(), from.yarnOf.end());
}

BodyRefusal yarnBodyRefusal(const YarnBodyParams& p, double yarnRadius, std::string* err) {
    const auto refuse = [&](BodyRefusal what, const std::string& why) {
        if (err) *err = "yarn body: " + why;
        return what;
    };
    if (p.plies < 1 || p.plies > 3) return refuse(BodyRefusal::Plies, "plies must be 1, 2 or 3");
    if (p.fibres < 1 || p.interiorFibres < 0) return refuse(BodyRefusal::FibreCount, "at least one shell fibre per yarn, no negative interior count");
    if (!(p.fibreRadius > 0) || !(yarnRadius > 0)) return refuse(BodyRefusal::FibreTooThin, "fibre and yarn radius must be positive");
    if (p.fibreRadius < kMinFibreRadius) return refuse(BodyRefusal::FibreTooThin, "fibres thinner than 10 um are below the strand files' radius floor");
    if (!(p.shellDepth > 0 && p.shellDepth <= 1)) return refuse(BodyRefusal::ShellDepth, "the shell depth lies in (0, 1]");
    if (!(p.fibreAngleDeg >= 0 && p.fibreAngleDeg <= 60) || !(p.plyAngleDeg >= 0 && p.plyAngleDeg <= 60))
        return refuse(BodyRefusal::Angles, "fibre and ply angles lie in [0, 60] degrees");
    if (!(p.irregularity >= 0 && p.irregularity <= 2)) return refuse(BodyRefusal::Irregularity, "the irregularity lies in [0, 2]");
    const BodyLayout l = bodyLayout(p, yarnRadius);
    char line[256];
    if (l.rOut < 0.5 * l.rp) {
        std::snprintf(line, sizeof line, "fibres of radius %.4f mm (up to %.4f with crimp and irregularity) are too large for a ply of radius %.3f mm",
                      p.fibreRadius, l.rMax + l.crimp, l.rp);
        return refuse(BodyRefusal::FibreTooLarge, line);
    }
    // Capacity: the start centres' area against each fibre's spacing disc, below jamming.
    const double perFibre = kPi * 0.25 * (kMinSpacingOverR * p.fibreRadius) * (kMinSpacingOverR * p.fibreRadius);
    const int perPly = (p.fibres + p.plies - 1) / p.plies;
    const int shellCap = static_cast<int>(kCapacityFraction * kPi * (l.rOut * l.rOut - l.rIn * l.rIn) / perFibre);
    if (perPly > shellCap) {
        std::snprintf(line, sizeof line, "%d fibres per ply do not fit a shell of depth %.2f (at most %d per ply, %d per yarn)",
                      perPly, p.shellDepth, shellCap, shellCap * p.plies);
        return refuse(BodyRefusal::ShellFull, line);
    }
    if (p.interiorFibres > 0) {
        const double rInt = l.rIn - 2 * l.rMax;
        const int interiorCap = rInt > 0 ? static_cast<int>(kCapacityFraction * kPi * rInt * rInt / perFibre) : 0;
        const int interiorPerPly = (p.interiorFibres + p.plies - 1) / p.plies;
        if (p.shellDepth >= 1 || interiorPerPly > interiorCap) {
            std::snprintf(line, sizeof line, "%d interior fibres per ply do not fit inside a shell of depth %.2f (at most %d per ply)",
                          interiorPerPly, p.shellDepth, interiorCap);
            return refuse(BodyRefusal::InteriorFull, line);
        }
    }
    if (p.core != BodyCore::None && l.rCore < 0.1 * l.rp) {
        const double maxDepth = 1 - (0.1 * l.rp + l.rMax) / l.rOut;
        std::snprintf(line, sizeof line, "no room for a core at shell depth %.2f (a core needs a depth of at most %.2f)", p.shellDepth, maxDepth);
        return refuse(BodyRefusal::NoCoreRoom, line);
    }
    return BodyRefusal::None;
}

PlyAxis bodyPlyAxis(const YarnBodyParams& p, double yarnRadius, uint64_t partSeed, int yarn, int ply, double z) {
    const YarnTwist tw = yarnTwist(bodyLayout(p, yarnRadius), p, yarnSeedOf(partSeed, p.seed, yarn), z + 2 * kGridStep);
    const PlyFrame f = plyFrame(tw, ply, p.plies, z);
    return {{f.C.x, f.C.y}, cross(f.u, f.v)};
}

std::vector<PlyAxis> bodyPlyAxes(const YarnBodyParams& p, double yarnRadius, uint64_t partSeed, int yarn, int ply, double length,
                                 double step) {
    const int steps = std::max(1, static_cast<int>(std::ceil(length / step)));
    const YarnTwist tw = yarnTwist(bodyLayout(p, yarnRadius), p, yarnSeedOf(partSeed, p.seed, yarn), steps * step + 2 * kGridStep);
    std::vector<PlyAxis> out(static_cast<size_t>(steps));
    for (int i = 0; i < steps; ++i) {
        const PlyFrame f = plyFrame(tw, ply, p.plies, (i + 0.5) * step);
        out[static_cast<size_t>(i)] = {{f.C.x, f.C.y}, cross(f.u, f.v)};
    }
    return out;
}

OverlapStats measureOverlap(const BodyStrands& b, const std::vector<OverlapWindow>& windows) {
    OverlapStats out;
    struct Segment {
        Vec3 a, b;
        double ra, rb;
        int strand, yarn;
        bool inWindow;
    };
    for (const OverlapWindow& w : windows) {
        std::vector<Segment> segs;
        Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
        for (size_t i = 0; i < b.strands(); ++i) {
            if (b.yarnOf[i] != w.yarn) continue;
            for (uint64_t v = b.first[i]; v + 1 < b.first[i + 1]; ++v) {
                if (b.s[v] < w.s0 || b.s[v + 1] > w.s1) continue;
                segs.push_back({b.points[v], b.points[v + 1], b.radius[v], b.radius[v + 1], static_cast<int>(i), w.yarn, true});
                for (const Vec3& p : {b.points[v], b.points[v + 1]}) {
                    lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
                    hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
                }
                out.fibreLength += norm(b.points[v + 1] - b.points[v]);
            }
        }
        if (segs.empty()) continue;
        const auto inside = [&](const Vec3& p) {
            return p.x >= lo.x && p.y >= lo.y && p.z >= lo.z && p.x <= hi.x && p.y <= hi.y && p.z <= hi.z;
        };
        for (size_t i = 0; i < b.strands(); ++i) {
            if (b.yarnOf[i] == w.yarn) continue;
            for (uint64_t v = b.first[i]; v + 1 < b.first[i + 1]; ++v)
                if (inside(b.points[v]) || inside(b.points[v + 1]))
                    segs.push_back({b.points[v], b.points[v + 1], b.radius[v], b.radius[v + 1], static_cast<int>(i), b.yarnOf[i], false});
        }
        double cell = 1e-3;
        for (const Segment& s : segs) cell = std::max(cell, norm(s.b - s.a) + s.ra + s.rb);
        const auto key = [&](long long x, long long y, long long z) { return (x * 73856093LL) ^ (y * 19349663LL) ^ (z * 83492791LL); };
        std::unordered_map<long long, std::vector<int>> grid;
        for (size_t k = 0; k < segs.size(); ++k) {
            const Vec3 m = (segs[k].a + segs[k].b) * 0.5;
            grid[key(std::llround(std::floor(m.x / cell)), std::llround(std::floor(m.y / cell)), std::llround(std::floor(m.z / cell)))].push_back(static_cast<int>(k));
        }
        std::vector<char> involved(segs.size(), 0);
        for (size_t k = 0; k < segs.size(); ++k) {
            const Segment& A = segs[k];
            if (!A.inWindow) continue;
            ++out.segments;
            const Vec3 m = (A.a + A.b) * 0.5;
            const long long cx = std::llround(std::floor(m.x / cell)), cy = std::llround(std::floor(m.y / cell)), cz = std::llround(std::floor(m.z / cell));
            for (long long dx = -1; dx <= 1; ++dx)
                for (long long dy = -1; dy <= 1; ++dy)
                    for (long long dz = -1; dz <= 1; ++dz) {
                        const auto it = grid.find(key(cx + dx, cy + dy, cz + dz));
                        if (it == grid.end()) continue;
                        for (int j : it->second) {
                            const Segment& B = segs[static_cast<size_t>(j)];
                            if (B.strand == A.strand || (B.inWindow && j <= static_cast<int>(k))) continue;
                            const ClosestSegmentPoints closestPair = closestSegmentPoints(A.a, A.b, B.a, B.b);
                            const double d = closestPair.distance;
                            const double reach = lerp(A.ra, A.rb, closestPair.s) + lerp(B.ra, B.rb, closestPair.t);
                            if (d >= reach) continue;
                            (A.yarn == B.yarn ? out.pairsWithinYarn : out.pairsBetweenYarns) += 1;
                            involved[k] = 1;
                            involved[static_cast<size_t>(j)] = 1;
                            out.maxDepth = std::max(out.maxDepth, (reach - d) / reach);
                        }
                    }
        }
        for (size_t k = 0; k < segs.size(); ++k) out.segmentsInvolved += involved[k] && segs[k].inWindow;
    }
    return out;
}

std::string BodyStats::line() const {
    char t[1024];
    const double perMm = overlap.fibreLength > 0 ? (overlap.pairsWithinYarn + overlap.pairsBetweenYarns) / overlap.fibreLength : 0;
    std::snprintf(t, sizeof t,
                  "%lld fibres, %lld pieces (%lld dropped), %lld segments, fibre length %.1f m over %.2f m of yarn; "
                  "r_out %.3f, r_in %.3f, r_max %.4f mm; packing shell %.1f %%, ply %.1f %%; soft-turn segments %.2f %% "
                  "(max %.1f, p99 %.1f deg); max protrusion %.4f mm, centreline at %.3f of its ring's ellipse; overlap windows %.1f mm of fibre: %.3f penetrating pairs / mm "
                  "(%lld within, %lld between yarns), %.1f %% of segments, max depth %.2f; %.0f ms",
                  fibres, pieces, droppedPieces, segments, fibreLength / 1000, yarnLength / 1000, layout.rOut, layout.rIn, layout.rMax,
                  100 * layout.phiShell, 100 * layout.phiPly, segments ? 100.0 * softTurnSegments / segments : 0.0, maxTurnDeg, p99TurnDeg,
                  maxProtrusion, maxEnvelope, overlap.fibreLength, perMm, overlap.pairsWithinYarn, overlap.pairsBetweenYarns,
                  overlap.segments ? 100.0 * overlap.segmentsInvolved / overlap.segments : 0.0, overlap.maxDepth, ms);
    return t;
}

bool buildYarnBody(const CurveSet& yarns, const SurfaceModel& surface, const YarnBodyParams& p, uint64_t partSeed, YarnBody& out,
                   std::string* err, const ExecPolicy& exec) {
    const auto t0 = std::chrono::steady_clock::now();
    out = YarnBody{};
    if (surface.yarns.empty()) return true;
    const double r = surface.yarns[0].radius;
    if (!yarnBodyParamsValid(p, r, err)) return false;
    const BodyLayout l = bodyLayout(p, r);
    const int Y = static_cast<int>(surface.yarns.size()), n = p.plies;
    // One task per (zone kind, yarn, ply): shell zones first, then interior zones, merged in this order.
    const int zones = p.interiorFibres > 0 ? 2 : 1;
    std::vector<YarnTwist> twists(static_cast<size_t>(Y));
    std::vector<SectionTable> sections(static_cast<size_t>(Y));
    for (int y = 0; y < Y; ++y) sections[static_cast<size_t>(y)] = sectionTable(surface.yarns[static_cast<size_t>(y)]);
    for (int y = 0; y < Y; ++y)
        twists[static_cast<size_t>(y)] = yarnTwist(l, p, yarnSeedOf(partSeed, p.seed, y), surface.yarns[static_cast<size_t>(y)].length);
    const int tasks = zones * Y * n;
    std::vector<ZoneResult> results(static_cast<size_t>(tasks));
    std::vector<std::string> failures(static_cast<size_t>(tasks));
    const double period = 2 * kPi * l.rOut / std::tan(degToRad(std::max(p.fibreAngleDeg, kMigrationMinAngleDeg)));
    const double rInterior = l.rIn - 2 * l.rMax;
    const double trim = n > 1 ? l.rp : 0;
    parallelFor(tasks, 1, [&](int64_t begin, int64_t end, int64_t) {
        for (int64_t task = begin; task < end; ++task) {
            cancelCheckpoint();  // a task builds one ply zone along a whole yarn: seconds on a long yarn
            const int zone = static_cast<int>(task / (Y * n)), y = static_cast<int>(task / n % Y), j = static_cast<int>(task % n);
            const YarnSurface& ys = surface.yarns[static_cast<size_t>(y)];
            const uint64_t yarnSeed = yarnSeedOf(partSeed, p.seed, y);
            uint64_t seed = splitmix64(yarnSeed ^ (0xBF58476D1CE4E5B9ull * static_cast<uint64_t>(j + 1)));
            if (zone == 1) seed = splitmix64(seed ^ 0x1A7E4100ull);  // the interior's own stream: the shell is unchanged by it
            ZoneContext c{&yarns, &ys, &twists[static_cast<size_t>(y)], &sections[static_cast<size_t>(y)], &l, &p, j, trim};
            c.omegaFibre = -std::tan(degToRad(p.fibreAngleDeg)) / l.rOut;  // handedness -1 (S singles), against the Z ply
            const int count = zone == 0 ? (p.fibres + n - 1 - j) / n : (p.interiorFibres + n - 1 - j) / n;
            std::vector<Vec2> starts;
            const bool placed = zone == 0 ? placeStarts(seed, count, l.rIn, l.rOut, p.fibreRadius, starts)
                                          : placeStarts(seed, count, 0, rInterior, p.fibreRadius, starts);
            if (!placed) {
                char line[200];
                std::snprintf(line, sizeof line, "yarn body: yarn %d ply %d: %zu of %d %s fibres placed (no place after %d attempts)", y, j,
                              starts.size(), count, zone == 0 ? "shell" : "interior", kPlacementAttempts);
                failures[static_cast<size_t>(task)] = line;
                continue;
            }
            if (zone == 0) {
                c.rOuter = l.rOut;
                c.depth = p.shellDepth;
            } else {
                c.rOuter = rInterior;
                c.depth = 1;
            }
            const std::vector<Slot> slots = makeSlots(seed, starts, zone == 0 ? l.rIn : 0, c.rOuter, period, p.irregularity);
            buildZone(c, slots, seed, y, results[static_cast<size_t>(task)]);
        }
    }, exec);
    for (const std::string& f : failures)
        if (!f.empty()) {
            if (err) *err = f;
            return false;
        }
    // Merge in task order; the histogram gives the turn percentile.
    BodyStats& st = out.stats;
    st.layout = l;
    std::array<long long, 1800> turns{};
    BodyStrands& b = out.strands;
    for (const ZoneResult& z : results) {
        appendStrands(b, z.strands);
        st.pieces += z.pieces;
        st.droppedPieces += z.dropped;
        st.softTurnSegments += z.softTurn;
        st.fibreLength += z.fibreLength;
        st.maxProtrusion = std::max(st.maxProtrusion, z.maxProtrusion);
        st.maxEnvelope = std::max(st.maxEnvelope, z.maxEnvelope);
        for (size_t k = 0; k < turns.size(); ++k) turns[k] += z.turns[k];
    }
    st.fibres = static_cast<long long>(p.fibres + p.interiorFibres) * Y;
    st.segments = static_cast<long long>(b.points.size() - b.strands());
    long long total = 0, seen = 0;
    for (long long c : turns) total += c;
    for (size_t k = 0; k < turns.size(); ++k) {
        if (turns[k] > 0) st.maxTurnDeg = 0.1 * static_cast<double>(k + 1);
        seen += turns[k];
        if (st.p99TurnDeg == 0 && total > 0 && seen >= 0.99 * static_cast<double>(total)) st.p99TurnDeg = 0.1 * static_cast<double>(k + 1);
    }
    for (const YarnSurface& ys : surface.yarns) st.yarnLength += ys.length;
    // Overlap windows: up to three yarns, kOverlapWindow mm around each one's middle.
    std::vector<OverlapWindow> windows;
    for (int k = 1; k <= 3 && Y > 0; ++k) {
        const int y = std::min(Y - 1, k * Y / 4);
        const double L = surface.yarns[static_cast<size_t>(y)].length;
        if (L >= kOverlapWindow && (windows.empty() || windows.back().yarn != y))
            windows.push_back({y, 0.5 * (L - kOverlapWindow), 0.5 * (L + kOverlapWindow)});
    }
    st.overlap = measureOverlap(b, windows);
    if (p.core != BodyCore::None)
        for (int y = 0; y < Y; ++y)
            for (int j = 0; j < n; ++j) {
                ZoneContext c{&yarns, &surface.yarns[static_cast<size_t>(y)], &twists[static_cast<size_t>(y)], &sections[static_cast<size_t>(y)], &l,
                              &p, j, trim};
                addCoreTube(c, l.rCore, out.core);
            }
    st.ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    return true;
}

} // namespace knit
