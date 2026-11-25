// Small knitted pieces with known answers: validation presets, shaping and strand demos (fixtures.h).
#include "knitting/knit/fixtures.h"

#include <algorithm>
#include <cmath>

namespace knit {

const std::vector<Preset>& presets() {
    static const std::vector<Preset> list = [] {
        std::vector<Preset> v;
        PatternParams base;  // D = 1 mm, c = 0.1 mm, A = 0.8 D, 6 x 6

        PatternParams expanded = base;
        expanded.waleSpacing = 8.0;
        expanded.courseSpacing = 6.5;
        v.push_back({"expanded", "Open, easy to inspect: every crossing and depth change is visible.", expanded,
                     Status::Pass});

        // From `knit_check --search-tightest` (base = expanded): the uniform-scale line
        // gave w 5.05 / h 4.10 mm (rounded up); course spacing was then tightened with
        // w held (h >= 3.034 mm) and rounded up. Re-validated by the presets test.
        // Binding contact at smaller w: neighbouring loop heads of one course
        // (same yarn, side by side), not a crossing.
        PatternParams tight = base;
        tight.waleSpacing = 5.05;
        tight.courseSpacing = 3.05;
        v.push_back({"tightest",
                     "Smallest configuration verified VALID by the search (scale line, then course direction).",
                     tight, Status::Pass});

        // About 0.8x the tightest wale spacing: neighbouring heads penetrate, inter-yarn
        // gaps fall below the clearance, and tight bends make some self pairs UNRESOLVED.
        PatternParams bad = base;
        bad.waleSpacing = 4.0;
        bad.courseSpacing = 3.0;
        v.push_back({"invalid", "Deliberately too dense: demonstrates collision and clearance diagnostics.", bad,
                     Status::Fail});

        // Tight knit: the tightest spacing with a thicker (1.15 mm) yarn. Nominal circular tubes
        // overlap at the side-by-side heads (gate (a): INVALID, by design); compression flattens both
        // sides so the surfaces touch, and the independent surface check passes (no local folds, no
        // penetration at the sampled points).
        // With D = 1.2 the centreline bends tighter than the yarn radius (rho ~0.577 <
        // 0.6), so the surface folds locally; compression cannot repair that.
        PatternParams comp = tight;
        comp.yarnDiameter = 1.15;
        comp.clearance = 0.0;
        v.push_back({"compressed", "Tight knit (D = 1.15 mm on the tightest spacing): overlaps resolved by compression.",
                     comp, Status::Fail});
        return v;
    }();
    return list;
}

const Preset* findPreset(const std::string& name) {
    for (const Preset& p : presets())
        if (p.name == name) return &p;
    return nullptr;
}

namespace {

// Surface of revolution around +z whose ring at row rho has length n(rho) * w
// (piecewise linear in rho between the turn counts).
class RingProfileSurface : public BaseSurface {
public:
    RingProfileSurface(std::vector<int> n, double w, double h, double H) : n_(std::move(n)), w_(w), h_(h), H_(H) {}
    // Ring of turn r centred on its row (rho = r + 0.5); z is the height (rows are
    // placed by arc length along the profile through arcToZ).
    double R(double z) const {
        double rho = z / h_ - 1.5;
        int r0 = std::clamp(static_cast<int>(std::floor(rho)), 0, static_cast<int>(n_.size()) - 1);
        int r1 = std::min(r0 + 1, static_cast<int>(n_.size()) - 1);
        double t = std::clamp(rho - r0, 0.0, 1.0);
        return ((1 - t) * n_[r0] + t * n_[r1]) * w_ / (2 * kPi);
    }
    SurfaceSample at(double a, double b) const override {
        double z = a * H_, t = 2 * kPi * b, dR = (R(z + 1e-4) - R(z - 1e-4)) / 2e-4;
        Vec3 radial{std::cos(t), std::sin(t), 0};
        return {radial * R(z) + Vec3{0, 0, z}, normalized(radial - Vec3{0, 0, dR})};
    }
    SurfaceSample project(const Vec3& q) const override {
        Vec3 r{q.x, q.y, 0};
        const double t = std::atan2(q.y, q.x) / (2 * kPi);
        // Closest generator point: coarse scan, then golden-section refinement.
        auto dist = [&](double z) { return std::hypot(norm(r) - R(z), q.z - z); };
        double best = q.z, bd = 1e300;
        for (double z = q.z - 4; z <= q.z + 4; z += 0.05)
            if (dist(z) < bd) { bd = dist(z); best = z; }
        double lo = best - 0.05, hi = best + 0.05;
        const double g = 0.5 * (std::sqrt(5.0) - 1);
        for (int it = 0; it < 40; ++it) {
            double x1 = hi - g * (hi - lo), x2 = lo + g * (hi - lo);
            if (dist(x1) < dist(x2)) hi = x2; else lo = x1;
        }
        return at(0.5 * (lo + hi) / H_, t);
    }

    // Height z at arc length s along the profile from z = 0 (table, linear inverse).
    double arcToZ(double sArc) const {
        if (table_.empty()) {
            double z = 0, acc = 0;
            const double dz = 0.01;
            table_.push_back({0, 0});
            while (z < H_) {
                acc += std::hypot(dz, R(z + dz) - R(z));
                z += dz;
                table_.push_back({acc, z});
            }
        }
        auto it = std::lower_bound(table_.begin(), table_.end(), std::make_pair(sArc, -1e300));
        if (it == table_.begin()) return 0;
        if (it == table_.end()) return H_;
        auto pr = std::prev(it);
        double t = (sArc - pr->first) / std::max(1e-300, it->first - pr->first);
        return pr->second + t * (it->second - pr->second);
    }

private:
    std::vector<int> n_;
    double w_, h_, H_;
    mutable std::vector<std::pair<double, double>> table_;
};

} // namespace

bool makeShapingDemo(const std::string& kind, const ShapingParams& sh, ShapingDemo& out, const GaugeParams& g) {
    out.pp = PatternParams{};
    out.pp.waleSpacing = g.w;
    out.pp.courseSpacing = g.h;
    if (kind == "inc" || kind == "dec") {
        std::vector<int> counts = kind == "inc" ? std::vector<int>{8, 8, 9, 9, 9} : std::vector<int>{9, 9, 8, 8, 8};
        out.prog = buildProgram(counts);
        out.lay = planarLayout(out.prog, g);
        RowPlacement rows;
        rows.sh = sh;
        out.curves = placeRows(PlaneSurface(1, 1), out.lay, out.prog, out.pp, rows);
        finalizeCurves(out.curves);
        return true;
    }
    if (kind == "cone") {
        std::vector<int> counts = {30, 30, 30, 32, 34, 36, 38, 40, 40, 40, 38, 36, 34, 32, 30, 30, 30};
        out.prog = buildProgram(counts);
        const double H = (counts.size() + 4) * g.h;
        RingProfileSurface surf(counts, g.w, g.h, H);
        // Rows one h apart along the surface: the even helix's a is the arc length (rho + 1) h along the
        // profile, mapped to the height's fraction.
        HelixNodes nodes = evenHelix(out.prog, g.h, 1.0);
        for (HelixSpot& spot : nodes.stitches) spot.a = surf.arcToZ(spot.a) / H;
        nodes.before.a = surf.arcToZ(nodes.before.a) / H;
        nodes.after.a = surf.arcToZ(nodes.after.a) / H;
        out.lay = helicalLayout(surf, out.prog, g, nodes);
        out.curves = placeSpiral(surf, out.lay, out.prog, out.pp, sh);
        finalizeCurves(out.curves);
        return true;
    }
    return false;
}

StrandDemo makeStrandDemo(double length, double D) {
    StrandDemo d;
    d.pp.yarnDiameter = D;
    YarnPath y;
    y.name = "strand";
    const int n = std::max(2, static_cast<int>(std::ceil(length)));
    for (int i = -1; i <= n + 1; ++i) y.curve.ctrl.push_back({length * (i - 0.5 * n) / n, 0, 0});
    y.stitches.push_back({0, 0, 0.0, y.curve.uMax()});
    y.finalize();
    d.curves.yarns.push_back(std::move(y));
    d.curves.wales = d.curves.courses = 1;
    fillFlatTopology(d.curves);
    return d;
}

} // namespace knit
