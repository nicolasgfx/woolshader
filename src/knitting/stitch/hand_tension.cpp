#include "knitting/stitch/hand_tension.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include "knitting/stitch/stockinette.h"
#include "model/hash.h"

namespace knit {

namespace {

// Hash channels of a stitch's draws. A channel never changes meaning: that would move every stitch.
enum HandChannel { kDriftChannel = 0, kOwnChannel = 1, kLeanChannel = 2, kLoopChannel = 3, kLoopHeightChannel = 4, kWaveChannel = 5,
                   kSpotChannel = 6, kSpotSizeChannel = 7, kSpotSignChannel = 8, kColumnChannel = 9 };

// A knitter's tension drifts over about this many stitches along a row and rows up a wale.
constexpr double kDriftStitches = 7, kDriftRows = 3;
constexpr double kDriftShare = 0.5;       // of a stitch's tension; the rest is its own
constexpr double kLeanPerTension = 0.5;   // the largest lean, per unit of unevenness
constexpr double kWaveLength = 30;        // [mm] the distance over which the fabric rises and falls
constexpr double kShortestLoop = 0.5;     // of the tallest: a raised loop's height is drawn above this
constexpr float kLoopSlack = 0.15f;       // the least tension of a raised loop, whatever its own draw
constexpr double kSmallestSpot = 0.5;     // of the largest: a spot's change is drawn above this
constexpr double kSpotSigmas = 2.35;      // a spot's length is its width at half height: 2.35 sigma of a Gaussian
constexpr double kThinnerPerTension = 0.5;  // a yarn stretched by t keeps its volume: its radius falls by about t / 2
constexpr double kThinnest = 0.5;

// A knit stitch away from shaping, with a parent and a plain child: its shape is the template's.
bool isPlainStitch(const KnitTopology& topo, const StitchRecord& st) {
    if (st.op != StitchOp::K || st.face != StitchFace::Knit || st.nParents != 1 || st.nChildren != 1) return false;
    return topo.stitches[st.children[0]].op == StitchOp::K;
}

bool insideItsRun(const KnitTopology& topo, const StitchRecord& st) {
    const int n = static_cast<int>(topo.stitches.size());
    return st.id > 0 && st.id + 1 < n && topo.stitches[st.id - 1].yarn == st.yarn && topo.stitches[st.id + 1].yarn == st.yarn;
}

Vec3 inStitchFrame(const StitchRecord& st, const Vec3& p) {
    const Vec3 d = p - st.origin;
    return {dot(d, st.U), dot(d, st.V), dot(d, st.N)};
}

// The stitch's first control: control c sits at the knot u = c - 1.
int firstControl(const StitchRecord& st) { return static_cast<int>(std::floor(st.uBegin)) + 1; }

// Where a stitch sits on its fabric: whole numbers that a change of diameter keeps.
struct FabricPlace {
    int column = 0, row = 0;
};

FabricPlace fabricPlace(const StitchRecord& st) {
    if (st.lattice.chart >= 0) return {st.lattice.i, st.lattice.r};
    return {st.index, st.turn};
}

// Value noise in [-1, 1]: hashed values at whole coordinates, blended smoothly between them.
double slowNoise(uint64_t seed, const Vec3& p, int channel) {
    const double ix = std::floor(p.x), iy = std::floor(p.y), iz = std::floor(p.z);
    const auto smooth = [](double f) { return f * f * (3 - 2 * f); };
    const double sx = smooth(p.x - ix), sy = smooth(p.y - iy), sz = smooth(p.z - iz);
    const auto corner = [&](int a, int b, int c) {
        const uint64_t cell = splitmix64(seed + static_cast<uint64_t>(static_cast<int>(iz) + c));
        return 2 * hash01(cell, static_cast<int>(ix) + a, static_cast<int>(iy) + b, channel) - 1;
    };
    const auto plane = [&](int c) {
        const double low = corner(0, 0, c) * (1 - sx) + corner(1, 0, c) * sx;
        const double high = corner(0, 1, c) * (1 - sx) + corner(1, 1, c) * sx;
        return low * (1 - sy) + high * sy;
    };
    return plane(0) * (1 - sz) + plane(1) * sz;
}

uint64_t chartSeed(const HandRequest& request, const StitchRecord& st) {
    return splitmix64(request.seed + static_cast<uint64_t>(st.lattice.chart + 1));
}

// How far a column slides up its wale at a point x mm along the stitch. Neighbouring columns draw
// independently and a column changes slowly, so rows stop lining up across the fabric. Half way between
// two columns the noise is their mean, so the yarn runs on unbroken from one to the next. A spiral's
// shaping already moves its wales.
double columnSlide(const HandRequest& request, const StitchRecord& st, double x) {
    if (st.lattice.chart < 0) return 0;
    const FabricPlace at = fabricPlace(st);
    const Vec3 place{at.column + x / request.wale, at.row / request.settings.columnRows, 0};
    return request.settings.columnShift * slowNoise(chartSeed(request, st), place, kColumnChannel);
}

// How much of a stitch's motion control j takes: none at the sinkers, which it shares with its
// neighbours, and all of it at the crown.
double motionWeight(int j) {
    const double s = std::sin(kPi * j / kControlsPerStitch);
    return s * s;
}

} // namespace

std::vector<StitchTension> drawStitchTensions(const KnitTopology& topo, const HandRequest& request) {
    std::vector<StitchTension> out(topo.stitches.size());
    const double unevenness = request.settings.unevenness / 100;
    const double loopChance = request.settings.loopRate * request.stitchArea / 100;  // 100 mm2 to the cm2
    const double spotChance = request.settings.spotRate * request.stitchArea / 100;
    for (const StitchRecord& st : topo.stitches) {
        if (!isPlainStitch(topo, st) || !insideItsRun(topo, st)) continue;
        const FabricPlace at = fabricPlace(st);
        const uint64_t seed = chartSeed(request, st);
        const double drift = slowNoise(seed, {at.column / kDriftStitches, at.row / kDriftRows, 0}, kDriftChannel);
        const double own = 2 * hash01(seed, at.column, at.row, kOwnChannel) - 1;
        StitchTension& t = out[st.id];
        t.tension = static_cast<float>(unevenness * (kDriftShare * drift + (1 - kDriftShare) * own));
        t.lean = static_cast<float>(unevenness * kLeanPerTension * (2 * hash01(seed, at.column, at.row, kLeanChannel) - 1));
        if (hash01(seed, at.column, at.row, kSpotChannel) < spotChance) {
            const double size = kSmallestSpot + (1 - kSmallestSpot) * hash01(seed, at.column, at.row, kSpotSizeChannel);
            const double sign = hash01(seed, at.column, at.row, kSpotSignChannel) < 0.5 ? -1 : 1;
            t.spot = static_cast<float>(sign * size * request.settings.spotChange / 100);
            t.tension -= t.spot;  // a thick spot is slack and takes more room, a thin one is stretched: the solver moves the neighbours
        }
        if (hash01(seed, at.column, at.row, kLoopChannel) < loopChance) {
            const double height = kShortestLoop + (1 - kShortestLoop) * hash01(seed, at.column, at.row, kLoopHeightChannel);
            t.raise = static_cast<float>(request.settings.loopHeight * height);
            t.tension = -std::max(std::abs(t.tension), kLoopSlack);  // a loop that stands out holds spare yarn
        }
    }
    return out;
}

HandResult knitByHand(CurveSet& cs, const HandRequest& request, KnitNumerics& numerics, const ExecPolicy& exec) {
    const KnitTopology& topo = cs.topo;
    HandResult out;
    out.tensions = drawStitchTensions(topo, request);
    Controls P = controlsOf(cs);
    Controls unshrunk = P;  // every move but the shrink: what the rest is measured on
    for (const StitchRecord& st : topo.stitches) {
        const StitchTension& t = out.tensions[st.id];
        std::vector<Vec3>& ctrl = P[st.yarn];
        const int first = firstControl(st);
        for (int j = 0; j < kControlsPerStitch && first + j < static_cast<int>(ctrl.size()); ++j) {
            const Vec3 local = inStitchFrame(st, ctrl[first + j]);
            const Vec3 shrink = (st.U * local.x + st.V * local.y) * -double(t.tension);
            const double raise = local.z > 0 ? t.raise : 0;  // the loop's front, not the head behind its child
            const Vec3 wave = st.N * (request.settings.waveHeight * slowNoise(request.seed, ctrl[first + j] * (1 / kWaveLength), kWaveChannel));
            const Vec3 slide = st.V * columnSlide(request, st, local.x);
            ctrl[first + j] = ctrl[first + j] + (shrink + st.U * (t.lean * local.y) + st.N * raise) * motionWeight(j) + wave + slide;
            unshrunk[st.yarn][first + j] = unshrunk[st.yarn][first + j] + (st.U * (t.lean * local.y) + st.N * raise) * motionWeight(j) + wave + slide;
        }
    }
    numerics.store(cs, P, exec);
    setYarnThickness(cs, out.tensions, request);
    out.rest = restOfTensions(unshrunk, topo, out.tensions);
    return out;
}

void setYarnThickness(CurveSet& cs, const std::vector<StitchTension>& tensions, const HandRequest& request) {
    if (request.maxThickness <= 0) return;
    for (YarnPath& y : cs.yarns) y.thickness.assign(y.curve.ctrl.size(), 1.0f);
    for (const StitchRecord& st : cs.topo.stitches) {
        const StitchTension& t = tensions[st.id];
        YarnPath& y = cs.yarns[st.yarn];
        const int first = firstControl(st), n = static_cast<int>(y.thickness.size());
        for (int j = 0; j < kControlsPerStitch && first + j < n; ++j)
            y.thickness[first + j] -= static_cast<float>(kThinnerPerTension * t.tension * motionWeight(j));
        if (t.spot == 0) continue;
        const double sigma = request.settings.spotLength / kSpotSigmas / (y.length() / n);  // in controls
        const int centre = first + kControlsPerStitch / 2, reach = static_cast<int>(std::ceil(3 * sigma));
        for (int c = std::max(0, centre - reach); c <= std::min(n - 1, centre + reach); ++c)
            y.thickness[c] += static_cast<float>(t.spot * std::exp(-0.5 * (c - centre) * (c - centre) / (sigma * sigma)));
    }
    for (YarnPath& y : cs.yarns)
        for (float& f : y.thickness) f = std::clamp(f, static_cast<float>(kThinnest), static_cast<float>(request.maxThickness));
}

Controls restOfTensions(const Controls& P, const KnitTopology& topo, const std::vector<StitchTension>& tensions) {
    std::vector<std::vector<double>> scale(P.size());
    for (size_t y = 0; y < P.size(); ++y) scale[y].assign(P[y].size(), 1.0);
    for (const StitchRecord& st : topo.stitches) {
        const int first = firstControl(st);
        for (int j = 0; j < kControlsPerStitch && first + j < static_cast<int>(scale[st.yarn].size()); ++j)
            scale[st.yarn][first + j] = 1.0 - tensions[st.id].tension;  // a loop scaled by 1 - t holds 1 - t of the yarn (P: before the shrink)
    }
    Controls rest = P;
    for (size_t y = 0; y < P.size(); ++y)
        for (size_t c = 0; c + 1 < P[y].size(); ++c) rest[y][c + 1] = rest[y][c] + (P[y][c + 1] - P[y][c]) * scale[y][c];
    return rest;
}

RaisedLoops measureRaisedLoops(const CurveSet& cs, const std::vector<StitchTension>& tensions) {
    const KnitTopology& topo = cs.topo;
    // The front of a stitch: how far its front controls stand out, weighted as a raise moves them.
    const auto frontOf = [&](const StitchRecord& st) {
        double sum = 0, weight = 0;
        for (int j = 0; j < kControlsPerStitch; ++j) {
            const double z = inStitchFrame(st, cs.yarns[st.yarn].curve.ctrl[firstControl(st) + j]).z;
            if (z <= 0) continue;
            sum += z * motionWeight(j);
            weight += motionWeight(j);
        }
        return weight > 0 ? sum / weight : 0.0;
    };
    double plain = 0;
    int plainCount = 0;
    for (const StitchRecord& st : topo.stitches)
        if (tensions[st.id].raise == 0 && isPlainStitch(topo, st)) {
            plain += frontOf(st);
            ++plainCount;
        }
    if (plainCount > 0) plain /= plainCount;
    RaisedLoops out;
    for (const StitchRecord& st : topo.stitches) {
        const double raise = tensions[st.id].raise;
        if (raise == 0) continue;
        ++out.asked;
        if (frontOf(st) - plain >= 0.5 * raise) ++out.standing;
    }
    return out;
}

StitchSpread measureStitchSpread(const CurveSet& cs) {
    const KnitTopology& topo = cs.topo;
    std::vector<const StitchRecord*> plain;
    for (const StitchRecord& st : topo.stitches)
        if (isPlainStitch(topo, st) && !cs.yarns[st.yarn].reversed) plain.push_back(&st);  // reversed: the template runs backwards
    StitchSpread out;
    out.stitches = static_cast<int>(plain.size());
    if (plain.empty()) return out;
    std::array<Vec3, kControlsPerStitch> mean{};
    for (const StitchRecord* st : plain)
        for (int j = 0; j < kControlsPerStitch; ++j)
            mean[j] = mean[j] + inStitchFrame(*st, cs.yarns[st->yarn].curve.ctrl[firstControl(*st) + j]) * (1.0 / plain.size());
    double shape = 0, height = 0;
    for (const StitchRecord* st : plain)
        for (int j = 0; j < kControlsPerStitch; ++j) {
            const Vec3 d = inStitchFrame(*st, cs.yarns[st->yarn].curve.ctrl[firstControl(*st) + j]) - mean[j];
            shape += dot(d, d);
            height += d.z * d.z;
        }
    const double samples = static_cast<double>(plain.size()) * kControlsPerStitch;
    out.shape = std::sqrt(shape / samples);
    out.height = std::sqrt(height / samples);
    return out;
}

} // namespace knit
