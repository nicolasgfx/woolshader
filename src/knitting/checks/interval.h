#pragma once
// Tiny interval arithmetic used for enclosures (Krawczyk test, depth sign).
//
// Each basic operation is computed in round-to-nearest and then widened by one
// ulp outward with nextafter. For a single IEEE operation the exact result lies
// within half an ulp of the rounded one, so the widened interval encloses it.
// This is not a general verified-arithmetic library; it is sufficient for the
// short, explicit expressions used here (see the numerical contract in checks.h).
#include <algorithm>
#include <cmath>
#include <limits>

namespace knit {

struct Interval {
    double lo = 0, hi = 0;
    constexpr Interval() = default;
    constexpr Interval(double v) : lo(v), hi(v) {}
    constexpr Interval(double l, double h) : lo(l), hi(h) {}
    double mid() const { return 0.5 * (lo + hi); }
    double width() const { return hi - lo; }
    bool containsZero() const { return lo <= 0.0 && 0.0 <= hi; }
    // Strict inclusion of *this in the interior of o.
    bool strictlyInside(const Interval& o) const { return lo > o.lo && hi < o.hi; }
    bool overlaps(const Interval& o) const { return lo <= o.hi && o.lo <= hi; }
    double mag() const { return std::max(std::fabs(lo), std::fabs(hi)); }
};

inline double downUlp(double v) { return std::nextafter(v, -std::numeric_limits<double>::infinity()); }
inline double upUlp(double v) { return std::nextafter(v, std::numeric_limits<double>::infinity()); }
inline Interval widen(const Interval& a) { return {downUlp(a.lo), upUlp(a.hi)}; }
// Outward padding by an explicit allowance (absolute).
inline Interval pad(const Interval& a, double allowance) { return {a.lo - allowance, a.hi + allowance}; }

inline Interval operator+(const Interval& a, const Interval& b) { return widen({a.lo + b.lo, a.hi + b.hi}); }
inline Interval operator-(const Interval& a, const Interval& b) { return widen({a.lo - b.hi, a.hi - b.lo}); }
inline Interval operator-(const Interval& a) { return {-a.hi, -a.lo}; }
inline Interval operator*(const Interval& a, const Interval& b) {
    double p1 = a.lo * b.lo, p2 = a.lo * b.hi, p3 = a.hi * b.lo, p4 = a.hi * b.hi;
    return widen({std::min({p1, p2, p3, p4}), std::max({p1, p2, p3, p4})});
}
inline Interval hull(const Interval& a, const Interval& b) { return {std::min(a.lo, b.lo), std::max(a.hi, b.hi)}; }
inline Interval intersect(const Interval& a, const Interval& b) { return {std::max(a.lo, b.lo), std::min(a.hi, b.hi)}; }

} // namespace knit
