// The fibres of a fur field: rest-shape curves from discs and their settle (fur_loops.h has the method). The
// curve and the vertex placement are the prototype's (feature/furry-teddy fur.cpp fineCurl, addCurledPoints,
// furSegments), operation for operation: the frozen baseline's strands come out bit for bit.
#include "fibres/fur_loops.h"

#include <algorithm>
#include <chrono>
#include <cmath>

#include "fibres/strands_file.h"
#include "model/cancel.h"
#include "model/log.h"

namespace knit {

namespace {

constexpr int kCurlStepsPerSegment = 32;  // the fine curve's samples per segment, before the vertices are placed by chord

// A loop's rest curve after the root as a fine polyline from the root, at least `arc` long (fur_loops.h): the
// tangent turns in the plane (d, c) at a curvature rising to 1/r over (pi r / 2) rise, integrated by the midpoint
// rule; across the plane (w) it climbs kDiscClimbDiameters root diameters per turn.
std::vector<Vec3> fineLoop(const FurDisc& f, double arc, double step) {
    const Vec3 d = f.d, c = f.c, w = cross(d, c);
    const double R = f.radius, rise = 0.5 * kPi * R * f.rise, climb = kDiscClimbDiameters * 2 * f.rootRadius / (2 * kPi * R);
    std::vector<Vec3> fine{f.root};
    double s = 0, theta = 0;
    Vec3 p = f.root;
    while (s < arc) {
        const double kappa = std::min(1.0, (s + 0.5 * step) / rise) / R;
        const double turn = kappa * step;
        theta += 0.5 * turn;
        p += (d * std::cos(theta) + c * std::sin(theta)) * step + w * (climb * turn * R);
        theta += 0.5 * turn;
        s += step;
        fine.push_back(p);
    }
    return fine;
}

// The vertices of a loop after the root: on the fine curve, each exactly `seg` from the previous one.
void addLoopPoints(const FurDisc& f, int segments, double seg, std::vector<Vec3>& out) {
    const std::vector<Vec3> fine = fineLoop(f, seg * segments * 1.5 + seg, seg / kCurlStepsPerSegment);
    Vec3 vertex = f.root;
    size_t k = 1;
    for (int n = 1; n <= segments && k < fine.size(); ++n) {
        while (k < fine.size() && norm(fine[k] - vertex) < seg) ++k;
        if (k == fine.size()) break;
        // The point of the piece fine[k-1] -> fine[k] at exactly seg from the vertex: |wv + u dv| = seg.
        const Vec3 start = norm(fine[k - 1] - vertex) < seg ? fine[k - 1] : vertex;
        const Vec3 dv = fine[k] - start, wv = start - vertex;
        const double A = dot(dv, dv), B = 2 * dot(dv, wv), C = dot(wv, wv) - seg * seg;
        const double u = std::clamp((-B + std::sqrt(std::max(0.0, B * B - 4 * A * C))) / (2 * A), 0.0, 1.0);
        vertex = start + dv * u;
        out.push_back(vertex);
    }
    // A curve too tight to hold every vertex (never at the segment rule's bound): straight on.
    while (static_cast<int>(out.size()) < segments + 2) out.push_back(out.back() + normalized(out.back() - out[out.size() - 2]) * seg);
}

FurLayerCounts& layerCounts(FurReport& rep, FurLayer layer) { return layer == FurLayer::Loops ? rep.loops : rep.fibres; }

} // namespace

double furFibreLength(const FurDisc& d) {
    if (!(d.radius > 0)) return d.straightLength;
    return std::max(4 * kMinSegmentRadii * d.rootRadius, (2 * kPi * d.radius + 0.25 * kPi * d.radius * d.rise) * d.turns);
}

bool furLengthClamped(const FurDisc& d) {
    return d.radius > 0 && (2 * kPi * d.radius + 0.25 * kPi * d.radius * d.rise) * d.turns < 4 * kMinSegmentRadii * d.rootRadius;
}

int furSegmentCount(const FurDisc& d) {
    // A loop: at most an eighth of its circle per segment (the exporter's curvature rule).
    double target = kFurSegmentTarget;
    if (d.radius > 0) target = std::min(target, 2 * kPi * d.radius / 8);
    // Visible segments n; the polyline has n + 1 segments with the buried root. Every segment must keep
    // minRadii r <= L_seg <= 200 r at its vertices (the largest radius is the root's, the smallest the tip's
    // after the kTaperStep-per-segment taper bound).
    const double length = furFibreLength(d), rootRadius = d.rootRadius, tipRadius = d.tipRadius;
    int n = std::max(1, static_cast<int>(std::lround(length / std::max(1e-6, target))));
    const double minRadii = d.radius > 0 ? kDiscMinSegmentRadii : kMinSegmentRadii;
    const int nMax = std::max(1, static_cast<int>(std::floor(length / (minRadii * rootRadius))));
    n = std::min(n, nMax);
    for (; n < nMax; ++n) {  // too long a segment at the thin tip: more segments
        const double rTip = rootRadius * std::max(tipRadius / rootRadius, std::pow(kTaperStep, -double(n + 1)));
        if (length / n <= 190 * rTip) break;
    }
    return n;
}

Hair furRestShape(const FurDisc& d) {
    Hair h;
    h.tri = d.tri;
    h.index = d.index;
    h.rootRadius = d.rootRadius;
    h.tipRadius = d.tipRadius;
    h.normal = d.normal;
    h.key = d.budgetKey;
    h.points.push_back(d.root - d.normal * (kMinSegmentRadii * d.rootRadius));
    h.points.push_back(d.root);
    const int segments = furSegmentCount(d);
    const double seg = furFibreLength(d) / segments;
    if (d.radius > 0) addLoopPoints(d, segments, seg, h.points);
    else
        for (int k = 1; k <= segments; ++k) h.points.push_back(d.root + d.d * (seg * k));
    return h;
}

std::vector<Hair> growLoops(const FurField& field, const MeshQuery& body, const FurSettle& settle, const ExecPolicy& exec, FurReport* rep) {
    const auto t0 = std::chrono::steady_clock::now();
    const std::vector<FurDisc>& discs = field.discs;
    std::vector<Hair> hairs(discs.size());
    std::vector<char> ok(discs.size(), 0);
    parallelFor(static_cast<int64_t>(discs.size()), 256, [&](int64_t b, int64_t e, int64_t) {
        cancelCheckpoint();
        for (int64_t q = b; q < e; ++q) {
            const FurDisc& d = discs[static_cast<size_t>(q)];
            Hair h = furRestShape(d);
            ok[static_cast<size_t>(q)] = settleHair(h, settle.strand, settle.gravityFullLength, settle.collisionOffset, body, furFibreLength(d)) ? 1 : 0;
            hairs[static_cast<size_t>(q)] = std::move(h);
        }
    }, exec);
    std::vector<Hair> out;
    FurReport r;
    for (size_t q = 0; q < hairs.size(); ++q) {
        FurLayerCounts& counts = layerCounts(r, discs[q].layer);
        counts.clamped += furLengthClamped(discs[q]);
        if (!ok[q]) {
            ++counts.stuck;
            continue;
        }
        ++counts.grown;
        counts.segments += static_cast<int64_t>(hairs[q].points.size()) - 2;
        r.maxResidual = std::max(r.maxResidual, hairs[q].residual);
        out.push_back(std::move(hairs[q]));
    }
    if (rep) {
        for (FurLayer layer : {FurLayer::Loops, FurLayer::Fibres}) {
            FurLayerCounts &to = layerCounts(*rep, layer), &from = layerCounts(r, layer);
            to.clamped += from.clamped;
            to.stuck += from.stuck;
            to.grown += from.grown;
            to.segments += from.segments;
        }
        rep->maxResidual = std::max(rep->maxResidual, r.maxResidual);
        rep->msSettle += msSince(t0);
    }
    return out;
}

} // namespace knit
