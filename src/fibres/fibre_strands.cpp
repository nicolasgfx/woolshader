// A FibreSet as renderer strands: the generated polylines, or the late taper.
#include "fibres/fibre_strands.h"

#include <algorithm>
#include <cmath>
#include <span>

#include "model/parallel.h"
#include "model/vec.h"

namespace knit {

namespace {

constexpr int64_t kFibreChunk = 16384;

// The step bound's float32 margin (taperRules): two radii rounded independently change their ratio by
// at most about 2u (u = 2^-24).
constexpr double kFloatRatioMargin = 1.0 / (1 << 22);  // 4u

// The generated, checked polyline as it is (points[0] is the buried root), with Fibre::vertexRadius.
void generatedStrand(const FibreSet& fs, const Fibre& f, const SceneMapping& map, std::vector<float>& pts, std::vector<float>& radii) {
    pts.clear();
    radii.clear();
    const std::span<const Vec3f> points = fs.pointsOf(f);
    for (int j = 0; j < f.pointCount; ++j) {
        const Vec3 p = map.toExport(toVec3(points[static_cast<size_t>(j)]));
        pts.insert(pts.end(), {float(p.x), float(p.y), float(p.z)});
        radii.push_back(float(vertexRadius(f, j) / map.sigma));
    }
}

bool exported(const Fibre& f, int population) {
    return f.pointCount >= 2 && (population < 0 || static_cast<int>(f.population) == population);
}

// The radius at arc length x of a fibre of length L whose tail starts at x0 (StrandTaper):
// r0 before it, then r0 (tip + (1 - tip) (1 + cos(pi w)) / 2) with w running 0 .. 1 over the tail.
double lateTaperRadius(double x, double x0, double L, double r0, double tip) {
    if (x <= x0) return r0;
    const double w = std::min(1.0, (x - x0) / (L - x0));
    return r0 * (tip + (1 - tip) * 0.5 * (1 + std::cos(3.14159265358979323846 * w)));
}

// Where the taper start x0 of one polyline ended up (Onset and Refined; TaperReport counts them).
enum class OnsetPlacement { NotPlaced, Exact, Shifted, Impossible };

// One polyline's late taper: whether it was written, and the tail it was made for.
struct TailOutcome {
    bool tapered = false;
    StrandFault fault = StrandFault::None;  // the rule a refused taper broke (None: there was no tail)
    OnsetPlacement onset = OnsetPlacement::NotPlaced;
    double L = 0, x0 = 0;                   // the polyline's length and its taper start
};

// A vertex of the tapered strand: (arc length on the polyline, polyline vertex index or -1 for a new one).
using Station = std::pair<double, int>;

// The shortest segment taperRules allows between end radii ra and rb: 4 x the larger radius
// (minSegOverRMax; it covers StrandRules::minSegOverR, 4 x the smaller).
double shortestSegment(double ra, double rb) { return 4 * std::max(ra, rb); }

// New vertices keep 1 % more than the shortest segment: positions are float32 in scene metres, where
// a coordinate of 0.1 m carries about 1e-8 m of rounding, 1e-4 of a 4 r segment at r = 20 um.
constexpr double kPlacementSlack = 1.01;
bool roomFor(double length, double ra, double rb) {
    return length >= (1 - 1e-9) * kPlacementSlack * shortestSegment(ra, rb);  // 1e-9: the subtraction just made
}

// Dense: a sample wherever the radius has fallen by kStep since the last one (below StrandRules' 1.5,
// so a dropped sample cannot break it); a sample closer than four radii to a polyline vertex is
// dropped, the vertex standing in for it.
void addDenseSamples(std::vector<Station>& stations, const std::vector<double>& arc, double x0, double L, double r0, double tip) {
    constexpr double kStep = 1.35;
    std::vector<double> samples{x0};
    while (lateTaperRadius(samples.back(), x0, L, r0, tip) / lateTaperRadius(L, x0, L, r0, tip) > kStep) {
        const double from = lateTaperRadius(samples.back(), x0, L, r0, tip);
        double lo = samples.back(), hi = L;
        for (int i = 0; i < 48; ++i) {
            const double mid = 0.5 * (lo + hi);
            if (from / lateTaperRadius(mid, x0, L, r0, tip) <= kStep) lo = mid;
            else hi = mid;
        }
        samples.push_back(lo);
    }
    for (double x : samples) {
        const double gap = 4 * lateTaperRadius(x, x0, L, r0, tip);
        bool nearVertex = false;
        for (double a : arc) nearVertex = nearVertex || std::fabs(a - x) < gap;
        if (!nearVertex) stations.push_back({x, -1});
    }
}

// Onset and Refined: a vertex at the taper start x0, so the body keeps r0 up to it. When a polyline
// vertex is closer to x0 than a legal segment, the new vertex moves away from it by the shortest
// legal segment (after x0 it takes the profile's radius there, before x0 the body's); when the
// polyline segment holding x0 has no legal place at all, there is none. No polyline vertex moves.
OnsetPlacement addOnset(std::vector<Station>& stations, const std::vector<double>& arc, double x0, double L, double r0, double tip) {
    size_t k = 0;  // the polyline segment [arc[k], arc[k + 1]] that holds x0
    while (k + 2 < arc.size() && arc[k + 1] <= x0) ++k;
    const double a = arc[k], b = arc[k + 1];
    if (x0 <= a || x0 >= b) return OnsetPlacement::Exact;  // a polyline vertex is the taper start
    const double rb = lateTaperRadius(b, x0, L, r0, tip);
    double s = x0;
    OnsetPlacement placed = OnsetPlacement::Exact;
    if (!roomFor(s - a, r0, r0)) {
        s = a + kPlacementSlack * shortestSegment(r0, r0);  // after x0, where the radius is below r0
        placed = OnsetPlacement::Shifted;
    } else if (!roomFor(b - s, r0, rb)) {
        s = b - kPlacementSlack * shortestSegment(r0, rb);  // before x0, where the radius is r0
        placed = OnsetPlacement::Shifted;
    }
    const double rs = lateTaperRadius(s, x0, L, r0, tip);
    if (!roomFor(s - a, r0, rs) || !roomFor(b - s, rs, rb)) return OnsetPlacement::Impossible;
    stations.push_back({s, -1});
    return placed;
}

// The largest |r(s) - r_linear(s)| on [a, b], with r_linear running from ra at a to rb at b, and where
// it is. The error's extrema lie at a, b, the kink at x0, or where r'(s) equals the chord's slope m:
// on the tail r'(s) = -r0 (1 - tip) (pi / 2T) sin(pi w), w = (s - x0) / T, T = L - x0, so
// sin(pi w) = k = -2 m T / (pi r0 (1 - tip)) with the roots w = asin(k) / pi and 1 - asin(k) / pi.
struct IntervalError {
    double value = 0, at = 0;
};
IntervalError largestLinearError(double a, double b, double ra, double rb, double x0, double L, double r0, double tip) {
    const double slope = (rb - ra) / (b - a);
    double candidates[5] = {a, b, 0, 0, 0};
    int count = 2;
    if (x0 > a && x0 < b) candidates[count++] = x0;
    const double T = L - x0;
    if (T > 0 && tip < 1) {
        const double k = -2 * slope * T / (kPi * r0 * (1 - tip));
        if (k >= 0 && k <= 1) {
            const double w = std::asin(k) / kPi;
            for (double s : {x0 + w * T, x0 + (1 - w) * T})
                if (s > a && s < b) candidates[count++] = s;
        }
    }
    IntervalError e;
    for (int i = 0; i < count; ++i) {
        const double s = candidates[i];
        const double error = std::fabs(lateTaperRadius(s, x0, L, r0, tip) - (ra + slope * (s - a)));
        if (i == 0 || error > e.value) e = {error, s};
    }
    return e;
}

// Refined: splits each tail interval whose largest error is above tolerance x r0 at the point of that
// error, moved into the range where both new segments are legal, left to right, until the interval
// holds or is too short to split (measureTaper then finds the strand nonconforming). The radius never
// rises along the tail, so 4 r(a) after a and before b is always room enough.
void refineTail(std::vector<Station>& stations, double x0, double L, double r0, double tip, double tolerance) {
    std::sort(stations.begin(), stations.end());
    const double target = tolerance * r0 * (1 - 1e-5);  // the float32 radii add about 1e-7 r0
    for (size_t i = 0; i + 1 < stations.size(); ++i) {
        for (int split = 0; split < 64; ++split) {
            const double a = stations[i].first, b = stations[i + 1].first;
            if (b <= x0) break;
            const double ra = lateTaperRadius(a, x0, L, r0, tip), rb = lateTaperRadius(b, x0, L, r0, tip);
            const IntervalError e = largestLinearError(a, b, ra, rb, x0, L, r0, tip);
            if (e.value <= target) break;
            const double lo = a + kPlacementSlack * 4 * ra, hi = b - kPlacementSlack * 4 * ra;
            if (lo > hi) break;
            stations.insert(stations.begin() + static_cast<std::ptrdiff_t>(i) + 1, Station{std::clamp(e.at, lo, hi), -1});
        }
    }
}

// One polyline with a late taper, in the output's units: `points` as exported, `r0` its root radius,
// `emergeL` the arc length from the root to where it leaves the yarn. The vertices are the polyline's
// and those of the mode (TaperPoints); new vertices lie on the polyline's chords. Not tapered if the
// result breaks taperRules, and then `pts` and `radii` hold nothing useful.
TailOutcome lateTaperPolyline(const std::vector<Vec3>& points, double r0, double emergeL, const StrandTaper& taper,
                              std::vector<float>& pts, std::vector<float>& radii) {
    TailOutcome out;
    const int n = static_cast<int>(points.size());
    if (n < 2) return out;
    std::vector<double> arc(static_cast<size_t>(n), 0.0);
    for (int k = 1; k < n; ++k)
        arc[static_cast<size_t>(k)] = arc[static_cast<size_t>(k - 1)] + norm(points[static_cast<size_t>(k)] - points[static_cast<size_t>(k - 1)]);
    const double L = arc.back();
    const double visible = std::max(0.0, L - emergeL);
    const double tail = std::min(visible, std::max(taper.tail * visible, taper.minDiameters * 2 * r0));
    const double x0 = L - tail;
    out.L = L;
    out.x0 = x0;
    if (tail <= 0) return out;
    std::vector<Station> stations;
    for (int k = 0; k < n; ++k) stations.push_back({arc[static_cast<size_t>(k)], k});
    switch (taper.points) {
    case TaperPoints::Dense: addDenseSamples(stations, arc, x0, L, r0, taper.tip); break;
    case TaperPoints::Existing: break;
    case TaperPoints::Onset: out.onset = addOnset(stations, arc, x0, L, r0, taper.tip); break;
    case TaperPoints::Refined:
        out.onset = addOnset(stations, arc, x0, L, r0, taper.tip);
        refineTail(stations, x0, L, r0, taper.tip, taper.tolerance);
        break;
    }
    std::sort(stations.begin(), stations.end());
    pts.clear();
    radii.clear();
    int segment = 0;  // the chord [segment, segment + 1] holding the station
    for (const Station& st : stations) {
        Vec3 q;
        if (st.second >= 0) {
            q = points[static_cast<size_t>(st.second)];
        } else {
            while (segment + 2 < n && arc[static_cast<size_t>(segment + 1)] < st.first) ++segment;
            const double a = arc[static_cast<size_t>(segment)], b = arc[static_cast<size_t>(segment + 1)];
            const double t = b > a ? (st.first - a) / (b - a) : 0.0;
            q = points[static_cast<size_t>(segment)] + (points[static_cast<size_t>(segment + 1)] - points[static_cast<size_t>(segment)]) * t;
        }
        pts.insert(pts.end(), {float(q.x), float(q.y), float(q.z)});
        radii.push_back(float(lateTaperRadius(st.first, x0, L, r0, taper.tip)));
    }
    StrandSetF one;
    one.add(pts, radii);
    const StrandCheck check = checkStrands(one, taperRules(taper));
    out.tapered = check.ok;
    out.fault = check.fault;
    return out;
}

// One generated fibre with a late taper: its polyline mapped to the export, then lateTaperPolyline.
TailOutcome lateTaperStrand(const FibreSet& fs, const Fibre& f, const SceneMapping& map, const StrandTaper& taper,
                            std::vector<float>& pts, std::vector<float>& radii) {
    std::vector<Vec3> mapped;
    for (const Vec3f& p : fs.pointsOf(f)) mapped.push_back(map.toExport(toVec3(p)));
    return lateTaperPolyline(mapped, f.radius0 / map.sigma, f.emergeL / map.sigma, taper, pts, radii);
}

// One written tapered strand against its profile (TaperReport), from its stored float32 values: the
// arc length from the stored positions, r_linear from the stored radii.
void measureTaperedStrand(const std::vector<float>& pts, const std::vector<float>& radii, const TailOutcome& tail,
                          double r0, const StrandTaper& taper, long long fibreIndex, TaperReport& r) {
    constexpr int kSteps = 64;  // Simpson steps per interval for the areas, and samples of the local error
    double s = 0, worst = 0;
    for (size_t v = 0; v + 1 < radii.size(); ++v) {
        const Vec3 p{pts[3 * v], pts[3 * v + 1], pts[3 * v + 2]}, q{pts[3 * v + 3], pts[3 * v + 4], pts[3 * v + 5]};
        const double a = s, b = s + norm(q - p);
        s = b;
        const double ra = radii[v], rb = radii[v + 1];
        r.maxStoredRatio = std::max(r.maxStoredRatio, std::max(ra, rb) / std::min(ra, rb));
        if (a < tail.x0) {  // r_linear may fall below r0 before x0 (Existing, a shifted or missing onset)
            const double end = std::min(b, tail.x0);
            r.maxEarlyThinningOverR0 = std::max(r.maxEarlyThinningOverR0, (r0 - (ra + (rb - ra) * (end - a) / (b - a))) / r0);
        }
        if (b <= tail.x0) continue;
        const IntervalError e = largestLinearError(a, b, ra, rb, tail.x0, tail.L, r0, taper.tip);
        worst = std::max(worst, e.value);
        r.maxErrorOverLocalR = std::max(r.maxErrorOverLocalR, e.value / lateTaperRadius(e.at, tail.x0, tail.L, r0, taper.tip));
        double exact = 0, linear = 0;
        for (int i = 0; i <= kSteps; ++i) {
            const double t = double(i) / kSteps;
            const double rx = lateTaperRadius(a + (b - a) * t, tail.x0, tail.L, r0, taper.tip), lx = ra + (rb - ra) * t;
            const double weight = (i == 0 || i == kSteps) ? 1 : (i % 2 ? 4 : 2);
            exact += weight * rx * rx;
            linear += weight * lx * lx;
            r.maxErrorOverLocalR = std::max(r.maxErrorOverLocalR, std::fabs(rx - lx) / rx);
        }
        r.tailAreaExact += exact * (b - a) / (3 * kSteps);
        r.tailAreaLinear += linear * (b - a) / (3 * kSteps);
    }
    r.maxErrorOverR0 = std::max(r.maxErrorOverR0, worst / r0);
    if (taper.points == TaperPoints::Refined && worst > taper.tolerance * r0) {
        ++r.nonconforming;
        if (r.nonconformingFirst.size() < 20) r.nonconformingFirst.push_back(fibreIndex);
    }
}

StrandSetF concatenated(const std::vector<StrandSetF>& part) {
    StrandSetF s;
    size_t nv = 0, ns = 0;
    for (const StrandSetF& c : part) {
        nv += c.vertices();
        ns += c.strands();
    }
    s.xyz.reserve(3 * nv);
    s.radius.reserve(nv);
    s.first.reserve(ns + 1);
    for (const StrandSetF& c : part) {
        const uint32_t base = static_cast<uint32_t>(s.radius.size());
        s.xyz.insert(s.xyz.end(), c.xyz.begin(), c.xyz.end());
        s.radius.insert(s.radius.end(), c.radius.begin(), c.radius.end());
        for (size_t k = 1; k < c.first.size(); ++k) s.first.push_back(base + c.first[k]);
    }
    return s;
}

} // namespace

bool parseTaperPoints(const std::string& text, TaperPoints& out) {
    for (TaperPoints p : {TaperPoints::Dense, TaperPoints::Existing, TaperPoints::Onset, TaperPoints::Refined})
        if (text == taperPointsName(p)) {
            out = p;
            return true;
        }
    return false;
}

const char* taperPointsName(TaperPoints p) {
    switch (p) {
    case TaperPoints::Dense: return "dense";
    case TaperPoints::Existing: return "existing";
    case TaperPoints::Onset: return "onset";
    case TaperPoints::Refined: return "refined";
    }
    return "dense";
}

StrandRules taperRules(const StrandTaper& taper) {
    StrandRules rules;
    if (taper.points == TaperPoints::Dense) return rules;
    rules.maxRadiusStep = std::max(rules.maxRadiusStep, (1 / taper.tip) * (1 + kFloatRatioMargin) - 1);
    rules.minSegOverRMax = 4;
    return rules;
}

StrandSetF fibresToStrands(const FibreSet& fs, const SceneMapping& map, int population) {
    // In chunks of fibres in parallel, concatenated in order: the strands of the serial loop.
    const int64_t n = static_cast<int64_t>(fs.fibres.size());
    const int64_t nChunks = (n + kFibreChunk - 1) / kFibreChunk;
    std::vector<StrandSetF> part(static_cast<size_t>(std::max<int64_t>(1, nChunks)));
    parallelFor(n, kFibreChunk, [&](int64_t b0, int64_t e0, int64_t ci) {
        StrandSetF& s = part[static_cast<size_t>(ci)];
        std::vector<float> pts, radii;
        for (int64_t q = b0; q < e0; ++q) {
            const Fibre& f = fs.fibres[static_cast<size_t>(q)];
            if (!exported(f, population)) continue;
            generatedStrand(fs, f, map, pts, radii);
            s.add(pts, radii);
        }
    }, ExecPolicy{true});
    return concatenated(part);
}

TaperedStrands fibresToTaperedStrands(const FibreSet& fs, const SceneMapping& map, int population, const StrandTaper& taper) {
    const int64_t n = static_cast<int64_t>(fs.fibres.size());
    const int64_t nChunks = (n + kFibreChunk - 1) / kFibreChunk;
    const size_t parts = static_cast<size_t>(std::max<int64_t>(1, nChunks));
    std::vector<StrandSetF> part(parts);
    std::vector<long long> tapered(parts, 0), untapered(parts, 0);
    parallelFor(n, kFibreChunk, [&](int64_t b0, int64_t e0, int64_t ci) {
        const size_t c = static_cast<size_t>(ci);
        std::vector<float> pts, radii;
        for (int64_t q = b0; q < e0; ++q) {
            const Fibre& f = fs.fibres[static_cast<size_t>(q)];
            if (!exported(f, population)) continue;
            if (lateTaperStrand(fs, f, map, taper, pts, radii).tapered) {
                ++tapered[c];
            } else {
                generatedStrand(fs, f, map, pts, radii);
                ++untapered[c];
            }
            part[c].add(pts, radii);
        }
    }, ExecPolicy{true});
    TaperedStrands out;
    out.strands = concatenated(part);
    for (size_t c = 0; c < parts; ++c) {
        out.tapered += tapered[c];
        out.untapered += untapered[c];
    }
    return out;
}

long long taperedSegmentCount(const FibreSet& fs, const SceneMapping& map, const StrandTaper& taper, float keyBelow) {
    const int64_t n = static_cast<int64_t>(fs.fibres.size());
    const int64_t nChunks = (n + kFibreChunk - 1) / kFibreChunk;
    std::vector<long long> segments(static_cast<size_t>(std::max<int64_t>(1, nChunks)), 0);
    parallelFor(n, kFibreChunk, [&](int64_t b0, int64_t e0, int64_t ci) {
        std::vector<float> pts, radii;
        for (int64_t q = b0; q < e0; ++q) {
            const Fibre& f = fs.fibres[static_cast<size_t>(q)];
            if (!exported(f, -1) || !(float(f.key) < keyBelow)) continue;
            const bool tapered = taper.tail > 0 && lateTaperStrand(fs, f, map, taper, pts, radii).tapered;
            segments[static_cast<size_t>(ci)] += tapered ? static_cast<long long>(radii.size()) - 1 : f.segments();
        }
    }, ExecPolicy{true});
    long long total = 0;
    for (long long s : segments) total += s;
    return total;
}

TaperedStrands taperStrandSet(const StrandSetF& s, const StrandTaper& taper) {
    TaperedStrands out;
    std::vector<Vec3> points;
    std::vector<float> pts, radii;
    for (size_t k = 0; k < s.strands(); ++k) {
        const uint32_t a = s.first[k], b = s.first[k + 1];
        points.clear();
        for (uint32_t v = a; v < b; ++v) points.push_back({s.xyz[3 * v], s.xyz[3 * v + 1], s.xyz[3 * v + 2]});
        // The emergence is not in the file: the curve start (vertex 1) stands for it.
        const double emergeL = points.size() > 1 ? norm(points[1] - points[0]) : 0.0;
        if (lateTaperPolyline(points, s.radius[a], emergeL, taper, pts, radii).tapered) {
            ++out.tapered;
        } else {
            pts.assign(s.xyz.begin() + 3 * a, s.xyz.begin() + 3 * b);
            radii.assign(s.radius.begin() + a, s.radius.begin() + b);
            ++out.untapered;
        }
        out.strands.add(pts, radii);
    }
    return out;
}

TaperReport measureTaper(const FibreSet& fs, const SceneMapping& map, const StrandTaper& taper, float keyBelow) {
    const int64_t n = static_cast<int64_t>(fs.fibres.size());
    const int64_t nChunks = (n + kFibreChunk - 1) / kFibreChunk;
    std::vector<TaperReport> part(static_cast<size_t>(std::max<int64_t>(1, nChunks)));
    parallelFor(n, kFibreChunk, [&](int64_t b0, int64_t e0, int64_t ci) {
        TaperReport& r = part[static_cast<size_t>(ci)];
        std::vector<float> pts, radii;
        for (int64_t q = b0; q < e0; ++q) {
            const Fibre& f = fs.fibres[static_cast<size_t>(q)];
            if (!exported(f, -1) || !(float(f.key) < keyBelow)) continue;
            ++r.fibres;
            r.originalSegments += f.segments();
            const TailOutcome tail = taper.tail > 0 ? lateTaperStrand(fs, f, map, taper, pts, radii) : TailOutcome{};
            if (!tail.tapered) {
                if (tail.fault == StrandFault::None) ++r.fallbackNoTail;
                else ++r.fallbackRule[static_cast<int>(tail.fault)];
                generatedStrand(fs, f, map, pts, radii);
            } else {
                ++r.tapered;
                if (tail.onset == OnsetPlacement::Exact) ++r.onsetExact;
                if (tail.onset == OnsetPlacement::Shifted) ++r.onsetShifted;
                if (tail.onset == OnsetPlacement::Impossible) ++r.onsetImpossible;
                measureTaperedStrand(pts, radii, tail, f.radius0 / map.sigma, taper, q, r);
            }
            r.finalSegments += static_cast<long long>(radii.size()) - 1;
            for (float x : radii) {
                r.minRadius = std::min(r.minRadius, double(x));
                r.maxRadius = std::max(r.maxRadius, double(x));
            }
        }
    }, ExecPolicy{true});
    TaperReport total;
    for (const TaperReport& r : part) addTaperReport(total, r);
    return total;
}

void addTaperReport(TaperReport& total, const TaperReport& part) {
    total.fibres += part.fibres;
    total.tapered += part.tapered;
    total.fallbackNoTail += part.fallbackNoTail;
    for (int k = 0; k < kStrandFaultCount; ++k) total.fallbackRule[k] += part.fallbackRule[k];
    total.onsetExact += part.onsetExact;
    total.onsetShifted += part.onsetShifted;
    total.onsetImpossible += part.onsetImpossible;
    total.nonconforming += part.nonconforming;
    for (long long i : part.nonconformingFirst)
        if (total.nonconformingFirst.size() < 20) total.nonconformingFirst.push_back(i);
    total.originalSegments += part.originalSegments;
    total.finalSegments += part.finalSegments;
    total.maxErrorOverR0 = std::max(total.maxErrorOverR0, part.maxErrorOverR0);
    total.maxErrorOverLocalR = std::max(total.maxErrorOverLocalR, part.maxErrorOverLocalR);
    total.maxEarlyThinningOverR0 = std::max(total.maxEarlyThinningOverR0, part.maxEarlyThinningOverR0);
    total.tailAreaExact += part.tailAreaExact;
    total.tailAreaLinear += part.tailAreaLinear;
    total.minRadius = std::min(total.minRadius, part.minRadius);
    total.maxRadius = std::max(total.maxRadius, part.maxRadius);
    total.maxStoredRatio = std::max(total.maxStoredRatio, part.maxStoredRatio);
}

} // namespace knit
