#include "stitchmesh/guided_rows.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <utility>

namespace knit {

LookupResult lookupConstruction(const SlotKey& key) {
    LookupResult r;
    char b[96];
    if (key.bottomRole == EdgeBoundary::BindOff || key.bottomRole == EdgeBoundary::Side || key.topRole == EdgeBoundary::CastOn || key.topRole == EdgeBoundary::Side) {
        r.reason = "boundary-role-misplaced";
        return r;
    }
    if (key.bottom == 1 && key.top == 1) {
        r.kind = StitchFaceKind::Regular;
        r.construction = key.bottomRole == EdgeBoundary::CastOn ? (key.topRole == EdgeBoundary::BindOff ? "K + cast-on + bind-off" : "K + cast-on")
                                                                 : (key.topRole == EdgeBoundary::BindOff ? "K + bind-off" : "K");
        return r;
    }
    const bool inc = key.bottom == 1 && key.top == 2, dec = key.bottom == 2 && key.top == 1;
    if (inc || dec) {
        const bool boundary = key.bottomRole != EdgeBoundary::Interior || key.topRole != EdgeBoundary::Interior;
        if (boundary && !key.boundaryShaping) {
            r.reason = "shaping-at-boundary-strip";
            return r;
        }
        if (key.besidePath) {
            r.reason = "shaping-beside-shift-path";
            return r;
        }
        r.kind = inc ? StitchFaceKind::Increase : StitchFaceKind::Decrease;
        r.construction = inc ? "INC kp" : "DEC d12k";
        if (key.bottomRole == EdgeBoundary::CastOn) r.construction += " + cast-on";
        if (key.topRole == EdgeBoundary::BindOff) r.construction += " + bind-off";
        return r;
    }
    std::snprintf(b, sizeof b, "key-unsupported (%d, %d)", key.bottom, key.top);
    r.reason = b;
    return r;
}

bool assembleStrip(int strip, const std::vector<int>& below, const std::vector<int>& above, const std::vector<int>& shapingSlots,
                   EdgeBoundary bottomRole, EdgeBoundary topRole, TypedStitchMesh& s, std::vector<int>* shapingFaces, std::string* err,
                   bool boundaryShaping) {
    const int nb = static_cast<int>(below.size()), na = static_cast<int>(above.size());
    const int d = na - nb, e = std::abs(d);
    const int F = d >= 0 ? nb : na;
    if (nb < 1 || na < 1 || static_cast<int>(shapingSlots.size()) != e) {
        if (err) *err = "assembleStrip: the shaping slots do not match the count change";
        return false;
    }
    std::vector<char> shape(F, 0);
    for (int slot : shapingSlots) {
        if (slot < 0 || slot >= F || shape[slot]) {
            if (err) *err = "assembleStrip: a shaping slot outside the strip or repeated";
            return false;
        }
        shape[slot] = 1;
    }
    int b = 0, t = 0;
    for (int slot = 0; slot < F; ++slot) {
        const int nBottom = shape[slot] && d < 0 ? 2 : 1;
        const int nTop = shape[slot] && d > 0 ? 2 : 1;
        std::vector<int> corners;
        std::vector<StitchEdgeRole> roles;
        std::vector<EdgeBoundary> bounds;
        for (int i = 0; i <= nBottom; ++i) corners.push_back(below[(b + i) % nb]);
        for (int i = nTop; i >= 0; --i) corners.push_back(above[(t + i) % na]);
        for (int i = 0; i < nBottom; ++i) {
            roles.push_back(StitchEdgeRole::Bottom);
            bounds.push_back(bottomRole);
        }
        roles.push_back(StitchEdgeRole::Wale);
        bounds.push_back(EdgeBoundary::Interior);
        for (int i = 0; i < nTop; ++i) {
            roles.push_back(StitchEdgeRole::Top);
            bounds.push_back(topRole);
        }
        roles.push_back(StitchEdgeRole::Wale);
        bounds.push_back(EdgeBoundary::Interior);
        SlotKey key;
        key.boundaryShaping = boundaryShaping;
        key.bottom = nBottom;
        key.top = nTop;
        key.bottomRole = bottomRole;
        key.topRole = topRole;
        key.besidePath = slot == 0 || slot == F - 1;
        const LookupResult lk = lookupConstruction(key);
        if (shape[slot] && shapingFaces) shapingFaces->push_back(s.mesh.faces());
        s.mesh.addFace(corners);
        s.kind.push_back(lk.kind);
        s.reason.push_back(lk.reason);
        s.role.push_back(roles);
        s.boundary.push_back(bounds);
        s.strip.push_back(strip);
        s.slot.push_back(slot);
        b += nBottom;
        t += nTop;
    }
    if (b != nb || t != na) {
        if (err) *err = "assembleStrip: the slots do not consume both rings";
        return false;
    }
    return true;
}

namespace {

// A boundary loop as a contour oriented so that the mesh lies on its left (cast-on) or right (bind-off), seen from
// outside; each point carries the triangle of the edge to the next point.
bool loopContour(const TriMesh& m, const std::vector<int>& loop, bool meshOnLeft, Contour& out, std::string* err) {
    std::map<std::pair<int, int>, int> directed;   // a -> b as in a triangle: that triangle lies on the left of a -> b
    for (int i = 0; i < static_cast<int>(m.tri.size()); ++i)
        for (int k = 0; k < 3; ++k) directed[{m.tri[i][k], m.tri[i][(k + 1) % 3]}] = i;
    std::vector<int> l = loop;
    const int n = static_cast<int>(l.size());
    if (n < 3) {
        if (err) *err = "a boundary loop with fewer than 3 vertices";
        return false;
    }
    const bool asTriangles = directed.count({l[0], l[1]}) > 0;
    if (asTriangles != meshOnLeft) std::reverse(l.begin(), l.end());
    out = Contour{};
    out.closed = true;
    for (int i = 0; i < n; ++i) {
        const int a = l[i], b = l[(i + 1) % n];
        RingPoint p;
        p.p = m.pos[a];
        p.a = a;
        p.b = a;
        p.f = 0;
        auto it = directed.find({a, b});
        if (it == directed.end()) it = directed.find({b, a});
        if (it == directed.end()) {
            if (err) *err = "a boundary loop edge that no triangle has";
            return false;
        }
        p.tri = it->second;
        out.points.push_back(p);
        out.length += norm(m.pos[b] - m.pos[a]);
    }
    return true;
}

struct Sampled {
    Vec3 p;
    VertexBinding binding;
};

// n points at arc fractions j / n from the contour parameter `origin` (point i to i + 1, fraction f).
std::vector<Sampled> sampleContour(const TriMesh& m, const Contour& c, double origin, int n) {
    const int np = static_cast<int>(c.points.size());
    std::vector<double> cum(np + 1, 0);
    for (int i = 0; i < np; ++i) cum[i + 1] = cum[i] + norm(c.points[(i + 1) % np].p - c.points[i].p);
    const double L = cum[np];
    if (origin >= np) origin -= np;   // a crossing exactly at the closing point
    const int i0 = static_cast<int>(origin);
    const double s0 = cum[i0] + (origin - i0) * (cum[i0 + 1] - cum[i0]);
    std::vector<Sampled> out;
    int seg = 0;
    for (int j = 0; j < n; ++j) {
        double s = s0 + L * j / n;
        if (s >= L) s -= L;
        while (seg + 1 < np && cum[seg + 1] <= s) ++seg;
        while (seg > 0 && cum[seg] > s) --seg;
        const double len = cum[seg + 1] - cum[seg];
        const double f = len > 0 ? (s - cum[seg]) / len : 0;
        const Vec3 p = c.points[seg].p + f * (c.points[(seg + 1) % np].p - c.points[seg].p);
        Sampled smp;
        smp.p = p;
        const int tri = c.points[seg].tri;
        smp.binding.tri = tri;
        if (tri >= 0) {
            const TrianglePoint tp = closestOnTriangle(p, m.pos[m.tri[tri][0]], m.pos[m.tri[tri][1]], m.pos[m.tri[tri][2]]);
            for (int k = 0; k < 3; ++k) smp.binding.bary[k] = tp.bary[k];
            smp.binding.distance = norm(tp.p - p);
        }
        out.push_back(smp);
    }
    return out;
}

// RowPhase::WalesUpright. Ring `upper` is turned on its contour so that the wale edges between it and the ring below lean
// as much one way as the other: for each wale edge (a below, b above) the arc on the upper contour from b to the point
// nearest a; the ring's origin moves by the mean of those arcs and the ring is sampled again. Returns the arc it moved
// by [mm]. The rings are turned from the first to the last, each against the ring below as already turned.
double turnRingUpright(const TriMesh& m, const Contour& c, double& originParameter, const std::vector<int>& lower, const std::vector<int>& upper,
                       const std::vector<int>& stripFaces, TypedStitchMesh& s) {
    const int np = static_cast<int>(c.points.size()), n = static_cast<int>(upper.size());
    if (np < 3 || n < 3) return 0;
    std::vector<double> cum(np + 1, 0);
    for (int i = 0; i < np; ++i) cum[i + 1] = cum[i] + norm(c.points[(i + 1) % np].p - c.points[i].p);
    const double L = cum[np];
    if (originParameter >= np) originParameter -= np;
    const int i0 = static_cast<int>(originParameter);
    const double s0 = cum[i0] + (originParameter - i0) * (cum[i0 + 1] - cum[i0]);
    auto arcOf = [&](const Vec3& p) {
        double best = 1e300, arc = 0;
        for (int i = 0; i < np; ++i) {
            const Vec3 a = c.points[i].p, d = c.points[(i + 1) % np].p - a;
            const double dd = dot(d, d);
            const double f = dd > 0 ? std::clamp(dot(p - a, d) / dd, 0.0, 1.0) : 0.0;
            const double distance = norm(p - (a + d * f));
            if (distance < best) {
                best = distance;
                arc = cum[i] + f * (cum[i + 1] - cum[i]);
            }
        }
        return arc;
    };
    std::vector<int> indexAbove(s.mesh.pos.size(), -1);
    std::vector<char> below(s.mesh.pos.size(), 0);
    for (int j = 0; j < n; ++j) indexAbove[upper[j]] = j;
    for (int v : lower) below[v] = 1;
    double sum = 0;
    int count = 0;
    for (int f : stripFaces) {
        const int sides = s.mesh.sides(f);
        for (int k = 0; k < sides; ++k) {
            if (s.role[f][k] != StitchEdgeRole::Wale) continue;
            int a = s.mesh.vertexAt(f, k), b = s.mesh.vertexAt(f, (k + 1) % sides);
            if (indexAbove[a] >= 0) std::swap(a, b);
            if (indexAbove[b] < 0 || !below[a]) continue;
            double lean = std::fmod(arcOf(s.mesh.pos[a]) - (s0 + L * indexAbove[b] / n), L);
            if (lean > 0.5 * L) lean -= L;
            if (lean < -0.5 * L) lean += L;
            sum += lean;
            ++count;
        }
    }
    if (count == 0) return 0;
    const double turn = sum / count;   // every wale edge is seen from its two faces: the mean is the same
    double arc = std::fmod(s0 + turn, L);
    if (arc < 0) arc += L;
    int seg = 0;
    while (seg + 1 < np && cum[seg + 1] <= arc) ++seg;
    originParameter = seg + (cum[seg + 1] > cum[seg] ? (arc - cum[seg]) / (cum[seg + 1] - cum[seg]) : 0.0);
    const std::vector<Sampled> smp = sampleContour(m, c, originParameter, n);
    for (int j = 0; j < n; ++j) {
        s.mesh.pos[upper[j]] = smp[j].p;
        s.binding[upper[j]] = smp[j].binding;
    }
    return turn;
}

const char* kindWord(int d) { return d > 0 ? "INC" : "DEC"; }

// The level of ring k, with any contour too small for a round of its own beside the row knitted over (row_rings.h);
// each one is recorded in g.knittedOver and logged.
RingLevel rowLevel(const TriMesh& m, const std::vector<double>& t, int k, double level, double w, GuidedRowsProduct& g) {
    RingLevel lv = extractLevel(m, t, level);
    for (const Contour& c : dropContoursShorterThan(lv, kKnitOverStitches * w)) {
        KnittedOverBump bump;
        bump.ring = k;
        bump.length = c.length;
        bump.at = c.points.front().p;
        g.knittedOver.push_back(bump);
        char b[160];
        std::snprintf(b, sizeof b, "ring %d (t = %.4f): a closed contour of %.2f mm beside the row (under %.2f mm, too small for a round of 3 stitches) is knitted over",
                      k, level, c.length, kKnitOverStitches * w);
        g.log.push_back(b);
    }
    return lv;
}


// The shaping slots of strip k by the frozen rule (header), then its faces. shapingStrip counts the strips that needed shaping
// so far (the stagger's index). Appends the placements; unsupported becomes true when a slot is typed Unsupported.
bool shapeStrip(int k, const std::vector<int>& below, const std::vector<int>& above, double lengthBelow, double lengthAbove, double w, double stagger,
                EdgeBoundary bottomRole, EdgeBoundary topRole, int& shapingStrip, TypedStitchMesh& s, std::vector<ShapingPlacement>& shaping, bool& unsupported,
                std::string* err, bool adjacentShaping = false, bool boundaryShaping = false) {
    const int nb = static_cast<int>(below.size()), na = static_cast<int>(above.size()), d = na - nb, e = std::abs(d);
    const int F = d >= 0 ? nb : na;
    if (e > F) {
        // one shaping event per stitch of the shorter ring at most: a strip that changes by more cannot be assembled.
        // It happens where an end is far from a row of the field (an opening or a finish curve that is not one level),
        // so the strip next to it must absorb the difference.
        const char* where = topRole == EdgeBoundary::BindOff ? "The finish" : bottomRole == EdgeBoundary::CastOn ? "The start" : "This row";
        if (err) {
            char b[480];
            std::snprintf(b, sizeof b,
                          "row %d goes from %d to %d stitches, %d %s in one row of %d stitches (at most one per stitch fits). %s lies far from a row of the "
                          "field here: draw it along one row (a level of the field), or start from the other end",
                          k + 1, nb, na, e, kindWord(d), F, where);
            *err = b;
        }
        return false;
    }
    std::vector<int> slots;
    std::vector<ShapingPlacement> placed;
    if (e > 0) {
        const double sigma = std::fmod(0.5 + shapingStrip * stagger, 1.0);
        ++shapingStrip;
        const int mMax = std::max(1, nb / (2 * e));
        std::set<int> used;
        for (int i = 0; i < e; ++i) {
            ShapingPlacement p;
            p.strip = k;
            p.kind = d > 0 ? StitchFaceKind::Increase : StitchFaceKind::Decrease;
            p.events = e;
            p.index = i;
            p.slotsInStrip = F;
            p.sigma = sigma;
            p.targetFraction = (i + sigma) / e;
            p.targetSlot = std::min(F - 1, static_cast<int>(std::floor(p.targetFraction * F)));
            auto allowed = [&](int slot) {
                if (slot < 1 || slot > F - 2) return false;
                for (int u : used)
                    if (std::abs(slot - u) < 2) return false;
                return true;
            };
            int found = -1;
            for (int dd = 0; dd <= mMax && found < 0; ++dd) {
                if (allowed(p.targetSlot + dd)) found = p.targetSlot + dd;
                else if (allowed(p.targetSlot - dd)) found = p.targetSlot - dd;
            }
            if (found < 0 && adjacentShaping) {
                // the variant: the nearest free slot between the ends, next to another shaping slot if need be
                for (int dd = 0; dd < F && found < 0; ++dd) {
                    for (int cand : {p.targetSlot + dd, p.targetSlot - dd})
                        if (found < 0 && cand >= 1 && cand <= F - 2 && !used.count(cand)) found = cand;
                }
                if (found >= 0) p.adjacent = true;
            }
            char b[400];
            const double Lb = lengthBelow / w, La = lengthAbove / w;
            if (found >= 0) {
                p.slot = found;
                p.placed = true;
                p.moved = found - p.targetSlot;
                std::snprintf(b, sizeof b,
                              "strip %d: ring %d has L/w %.3f -> n %d, ring %d has L/w %.3f -> n %d, so %d %s; event %d of %d at fraction (i + sigma)/e = %.4f "
                              "(sigma %.4f) -> target slot %d of %d; %s",
                              k, k, Lb, nb, k + 1, La, na, e, kindWord(d), i + 1, e, p.targetFraction, sigma, p.targetSlot, F,
                              p.adjacent ? (std::string("adjacent-shaping variant: slot ") + std::to_string(found) + ", next to another shaping slot").c_str()
                              : p.moved == 0 ? "placed there" : (std::string("moved by ") + std::to_string(p.moved) + " (beside the path or next to another shaping slot)").c_str());
            } else {
                p.placed = false;
                int nearest = -1;
                for (int dd = 0; dd < F && nearest < 0; ++dd) {
                    if (p.targetSlot + dd < F && !used.count(p.targetSlot + dd)) nearest = p.targetSlot + dd;
                    else if (p.targetSlot - dd >= 0 && !used.count(p.targetSlot - dd)) nearest = p.targetSlot - dd;
                }
                p.slot = nearest;
                p.moved = nearest - p.targetSlot;
                std::snprintf(b, sizeof b, "strip %d: %d %s, event %d: no allowed slot within %d of target %d: shaping-not-placeable", k, e, kindWord(d), i + 1, mMax,
                              p.targetSlot);
            }
            p.reason = b;
            used.insert(p.slot);
            slots.push_back(p.slot);
            placed.push_back(p);
        }
    }
    std::vector<int> shapingFaces;
    std::vector<int> sortedSlots = slots;
    std::sort(sortedSlots.begin(), sortedSlots.end());
    if (!assembleStrip(k, below, above, sortedSlots, bottomRole, topRole, s, &shapingFaces, err, boundaryShaping)) return false;
    // shapingFaces follow the sorted slots
    for (ShapingPlacement& p : placed) {
        const int idx = static_cast<int>(std::lower_bound(sortedSlots.begin(), sortedSlots.end(), p.slot) - sortedSlots.begin());
        p.face = shapingFaces[idx];
        if (!p.placed) {
            s.kind[p.face] = StitchFaceKind::Unsupported;
            s.reason[p.face] = "shaping-not-placeable";
        } else if (s.kind[p.face] == StitchFaceKind::Unsupported) {
            p.placed = false;
            p.reason += "; lookup: " + s.reason[p.face];
        }
        if (s.kind[p.face] == StitchFaceKind::Unsupported) unsupported = true;
        shaping.push_back(p);
    }
    return true;
}

} // namespace

bool buildGuidedRows(const TriMesh& m, const std::vector<double>& t, const std::vector<int>& castOn, const std::vector<int>& bindOff,
                     const GuidedRowsRequest& r, GuidedRowsProduct& out, std::string* err) {
    out = GuidedRowsProduct{};
    if (t.size() != m.pos.size() || castOn.empty() || bindOff.empty()) {
        if (err) *err = "buildGuidedRows: the field or a boundary loop is missing";
        return false;
    }
    double Tf = 1e300;
    for (int v : bindOff) Tf = std::min(Tf, t[v]);
    const int R = std::max(1, static_cast<int>(std::lround(Tf / r.h)));
    out.rows = R;
    out.extent = Tf;
    out.hPrime = Tf / R;
    // rings
    std::vector<Contour> contours(R + 1);
    for (int k = 0; k <= R; ++k) {
        RingSample rs;
        rs.k = k;
        rs.level = k * out.hPrime;
        if (k == 0 || k == R) {
            rs.boundary = true;
            if (!loopContour(m, k == 0 ? castOn : bindOff, k == 0, contours[k], err)) return false;
            rs.contours = 1;
        } else {
            RingLevel lv = rowLevel(m, t, k, rs.level, r.w, out);
            rs.contours = static_cast<int>(lv.contours.size());
            if (lv.contours.size() != 1 || !lv.contours.front().closed) {
                ++out.contourEvents;
                char b[160];
                std::snprintf(b, sizeof b, "ring %d (t = %.4f): %zu contours, the first %s: contour event, no mesh", k, rs.level, lv.contours.size(),
                              lv.contours.empty() ? "missing" : (lv.contours.front().closed ? "closed" : "open"));
                out.log.push_back(b);
                out.rings.push_back(rs);
                continue;
            }
            contours[k] = lv.contours.front();
        }
        int crossings = 0;
        if (!r.path) {
            rs.cutParameter = cutCrossing(contours[k], r.cutOrigin, r.cutNormal, r.cutSide, &crossings);
        } else if (rs.boundary) {
            rs.cutParameter = pathNearest(contours[k], *r.path);
            crossings = rs.cutParameter >= 0 ? 1 : 0;
        } else if (r.pathField) {
            std::vector<double> params;
            crossings = levelCrossings(*r.pathField, rs.level, &params);
            if (crossings == 1) {
                const int i = std::min(static_cast<int>(params[0]), static_cast<int>(r.pathField->size()) - 2);
                rs.cutParameter = contourParameterNear(contours[k], pathPointAt(*r.pathField, params[0]), (*r.pathField)[i].segmentTri, nullptr);
            }
        } else {
            rs.cutParameter = pathCrossing(m, contours[k], *r.path, &crossings);
        }
        rs.guideCrossings = crossings;
        rs.length = contours[k].length;
        rs.n = static_cast<int>(std::lround(rs.length / r.w));
        if (crossings != 1 || rs.cutParameter < 0) {
            ++out.contourEvents;
            char b[120];
            std::snprintf(b, sizeof b, "ring %d crosses the %s %d times: no row origin, no mesh", k, r.path ? "row-start path" : "cut guide", crossings);
            out.log.push_back(b);
        }
        out.rings.push_back(rs);
    }
    out.contours = contours;
    if (out.contourEvents > 0) {
        out.status = Status::Fail;
        out.detail = "contour events: the component needs a correction (a cut, a guide or a split), shown in the log";
        return true;
    }
    if (r.bindOffCount == BindOffCount::OfRoundBelow && R >= 1) out.rings[R].n = out.rings[R - 1].n;
    // vertices
    TypedStitchMesh& s = out.mesh;
    std::vector<std::vector<int>> ring(R + 1);
    for (int k = 0; k <= R; ++k) {
        const std::vector<Sampled> smp = sampleContour(m, contours[k], out.rings[k].cutParameter, out.rings[k].n);
        for (const Sampled& x : smp) {
            ring[k].push_back(static_cast<int>(s.mesh.pos.size()));
            s.mesh.pos.push_back(x.p);
            s.binding.push_back(x.binding);
        }
    }
    s.rings = ring;
    for (int k = 0; k < R; ++k)
        if (std::abs(out.rings[k + 1].n - out.rings[k].n) > out.rings[k].n / 2) ++out.limitBinds;
    // shaping slots and faces
    int shapingStrip = 0;
    bool unsupported = false;
    for (int k = 0; k < R; ++k) {
        const EdgeBoundary bottomRole = k == 0 ? EdgeBoundary::CastOn : EdgeBoundary::Interior;
        const EdgeBoundary topRole = k == R - 1 ? EdgeBoundary::BindOff : EdgeBoundary::Interior;
        if (!shapeStrip(k, ring[k], ring[k + 1], out.rings[k].length, out.rings[k + 1].length, r.w, r.stagger, bottomRole, topRole, shapingStrip, s, out.shaping,
                        unsupported, err, r.adjacentShaping, r.boundaryShaping))
            return false;
    }
    if (r.rowPhase == RowPhase::WalesUpright) {
        std::vector<std::vector<int>> facesOf(R);
        for (int f = 0; f < s.mesh.faces(); ++f)
            if (s.strip[f] >= 0 && s.strip[f] < R) facesOf[s.strip[f]].push_back(f);
        double total = 0, largest = 0;
        for (int k = 0; k < R; ++k) {
            total += turnRingUpright(m, contours[k + 1], out.rings[k + 1].cutParameter, ring[k], ring[k + 1], facesOf[k], s);
            largest = std::max(largest, std::fabs(total));
        }
        out.rowStartOffPath = largest;
    }
    for (int k = 0; k < R; ++k) s.shiftPath.push_back(ShiftPathEdge{k, ring[k][0], ring[k + 1][0]});
    for (int f = 0; f < s.mesh.faces(); ++f)
        if (s.kind[f] == StitchFaceKind::Unsupported) unsupported = true;
    out.status = unsupported ? Status::Unresolved : Status::Pass;
    char b[200];
    std::snprintf(b, sizeof b, "%d rows, h' %.5f, %d faces, %zu shaping slots, %s", R, out.hPrime, s.mesh.faces(), out.shaping.size(),
                  unsupported ? "unsupported slots listed" : "every slot typed");
    out.detail = b;
    if (r.rowPhase == RowPhase::WalesUpright) {
        std::snprintf(b, sizeof b, "; rings turned for upright wales, the row start up to %.1f mm off the path", out.rowStartOffPath);
        out.detail += b;
    }
    return true;
}

int levelCrossings(const std::vector<MarkedPathPoint>& path, double level, std::vector<double>* params) {
    int count = 0;
    for (size_t i = 0; i + 1 < path.size(); ++i) {
        const double a = path[i].t - level, b = path[i + 1].t - level;
        if ((a < 0) != (b < 0)) {   // values equal to the level count as above it, as the contour does
            ++count;
            if (params) params->push_back(static_cast<double>(i) + a / (a - b));
        }
    }
    return count;
}

Vec3 pathPointAt(const std::vector<MarkedPathPoint>& path, double u) {
    const int i = std::min(static_cast<int>(u), static_cast<int>(path.size()) - 2);
    return path[i].p + (u - i) * (path[i + 1].p - path[i].p);
}

double contourParameterNear(const Contour& c, const Vec3& x, int tri, double* distance) {
    const int np = static_cast<int>(c.points.size());
    double best = 1e300, param = -1;
    for (int j = 0; j < np; ++j) {
        const Vec3 a0 = c.points[j].p, a1 = c.points[(j + 1) % np].p;
        const Vec3 d = a1 - a0;
        const double len2 = dot(d, d);
        double f = len2 > 0 ? dot(x - a0, d) / len2 : 0;
        f = std::min(1.0, std::max(0.0, f));
        const double dist = norm(a0 + f * d - x);
        const bool sameTri = c.points[j].tri == tri;
        const double score = dist - (sameTri ? 1e-9 : 0);
        if (score < best) {
            best = score;
            param = j + std::min(f, 1.0 - 1e-12);
            if (distance) *distance = dist;
        }
    }
    return param;
}

bool buildMarkedRows(const TriMesh& m, const std::vector<double>& t, const MarkedRowsRequest& r, MarkedRowsProduct& out, std::string* err) {
    out = MarkedRowsProduct{};
    GuidedRowsProduct& g = out.rows;
    const bool loopFinish = !r.bindOffLoop.empty();
    double finish = r.finishLevel;
    if (loopFinish) {
        finish = 1e300;
        for (int v : r.bindOffLoop) finish = std::min(finish, t[v]);
    }
    if (t.size() != m.pos.size() || r.path.size() < 2 || !(finish > 0) || finish >= 1e300) {
        if (err) *err = "buildMarkedRows: the field, the row-start path or the finish is missing";
        return false;
    }
    const int R = std::max(1, static_cast<int>(std::lround(finish / r.h)));
    g.rows = R;
    g.extent = finish;
    g.hPrime = finish / R;
    // ring 0: the pole
    RingSample r0;
    r0.k = 0;
    r0.n = 1;
    r0.contours = 0;
    g.rings.push_back(r0);
    std::vector<Contour> contours(R + 1);
    for (int k = 1; k <= R; ++k) {
        RingSample rs;
        rs.k = k;
        rs.level = k * g.hPrime;
        rs.boundary = k == R;
        if (loopFinish && k == R) {
            // the finish curve: the bind-off ring is the loop itself (mesh on its right), its origin the point nearest the path
            PathCrossing pc;
            pc.ring = k;
            pc.level = rs.level;
            if (!loopContour(m, r.bindOffLoop, false, contours[k], err)) return false;
            rs.contours = 1;
            rs.length = contours[k].length;
            rs.n = static_cast<int>(std::lround(rs.length / r.w));
            const int np = static_cast<int>(contours[k].points.size());
            double best = 1e300;
            for (int j = 0; j < np; ++j)
                for (size_t i = 0; i + 1 < r.path.size(); ++i) {
                    const Vec3 a = r.path[i].p, d = r.path[i + 1].p - a;
                    const double len2 = dot(d, d);
                    const double f = len2 > 0 ? std::clamp(dot(contours[k].points[j].p - a, d) / len2, 0.0, 1.0) : 0.0;
                    const double dist = norm(a + f * d - contours[k].points[j].p);
                    if (dist < best) {
                        best = dist;
                        rs.cutParameter = j;
                    }
                }
            pc.crossings = 1;
            pc.contourParameter = rs.cutParameter;
            pc.locateError = best;
            g.rings.push_back(rs);
            out.crossings.push_back(pc);
            continue;
        }
        RingLevel lv = rowLevel(m, t, k, rs.level, r.w, g);
        rs.contours = static_cast<int>(lv.contours.size());
        PathCrossing pc;
        pc.ring = k;
        pc.level = rs.level;
        if (lv.contours.size() != 1 || !lv.contours.front().closed) {
            ++g.contourEvents;
            char b[200];
            std::snprintf(b, sizeof b, "ring %d (t = %.4f): %zu contours, the first %s: contour event, no mesh", k, rs.level, lv.contours.size(),
                          lv.contours.empty() ? "missing" : (lv.contours.front().closed ? "closed" : "open"));
            g.log.push_back(b);
            g.rings.push_back(rs);
            out.crossings.push_back(pc);
            continue;
        }
        contours[k] = lv.contours.front();
        rs.length = contours[k].length;
        rs.n = static_cast<int>(std::lround(rs.length / r.w));
        // the row-start path against this level
        pc.crossings = levelCrossings(r.path, rs.level, &pc.pathParameters);
        if (pc.crossings == 1) {
            const double u = pc.pathParameters[0];
            const int i = static_cast<int>(u);
            const double param = contourParameterNear(contours[k], pathPointAt(r.path, u), r.path[i].segmentTri, &pc.locateError);
            pc.contourParameter = param;
            rs.cutParameter = param;
        } else {
            ++out.guideProblems;
            std::string where;
            for (double u : pc.pathParameters) {
                char b[48];
                std::snprintf(b, sizeof b, "%s%.3f", where.empty() ? "" : ", ", u);
                where += b;
            }
            char b[400];
            std::snprintf(b, sizeof b, "ring %d (t = %.4f): the row-start path crosses it %d times (at path parameters %s): no row origin, a guide problem", k, rs.level,
                          pc.crossings, where.c_str());
            g.log.push_back(b);
        }
        g.rings.push_back(rs);
        out.crossings.push_back(pc);
    }
    g.contours = contours;
    if (g.contourEvents > 0 || out.guideProblems > 0) {
        g.status = Status::Fail;
        g.detail = g.contourEvents > 0 ? "contour events: the component needs a correction (a cut, a guide or a split)"
                                       : "the row-start path does not cross every ring exactly once: the guide needs a correction";
        return true;
    }
    // automatic: round(L1 / w), at least 3, raised until the first shaping strip's increases fit the placement rule in use
    // (frozen: (n1 - 1) / 2 per strip; adjacent shaping: n1 - 2)
    out.ring1Auto = std::max(3, g.rings[1].n);
    if (R >= 2)
        while (g.rings[2].n - out.ring1Auto > (r.adjacentShaping ? out.ring1Auto - 2 : (out.ring1Auto - 1) / 2)) ++out.ring1Auto;
    if (r.ring1 > 0) g.rings[1].n = r.ring1;
    else g.rings[1].n = out.ring1Auto;
    if (r.bindOffCount == BindOffCount::OfRoundBelow && R >= 2) g.rings[R].n = g.rings[R - 1].n;
    // vertices: the pole, then each ring from its row origin
    TypedStitchMesh& s = g.mesh;
    std::vector<std::vector<int>> ring(R + 1);
    ring[0].push_back(0);
    s.mesh.pos.push_back(r.pole);
    s.binding.push_back(r.poleBinding);
    for (int k = 1; k <= R; ++k) {
        const std::vector<Sampled> smp = sampleContour(m, contours[k], g.rings[k].cutParameter, g.rings[k].n);
        for (const Sampled& x : smp) {
            ring[k].push_back(static_cast<int>(s.mesh.pos.size()));
            s.mesh.pos.push_back(x.p);
            s.binding.push_back(x.binding);
        }
    }
    s.rings = ring;
    for (int k = 1; k < R; ++k)
        if (std::abs(g.rings[k + 1].n - g.rings[k].n) > g.rings[k].n / 2) ++g.limitBinds;
    // strip 0: the gathered pole start, (pole, ring 1 vertex j + 1, ring 1 vertex j) with roles (W+, T, W-)
    const int n1 = static_cast<int>(ring[1].size());
    for (int j = 0; j < n1; ++j) {
        s.mesh.addFace({ring[0][0], ring[1][(j + 1) % n1], ring[1][j]});
        s.kind.push_back(StitchFaceKind::PoleStart);
        s.reason.push_back("");
        s.role.push_back({StitchEdgeRole::Wale, StitchEdgeRole::Top, StitchEdgeRole::Wale});
        s.boundary.push_back({EdgeBoundary::Interior, R == 1 ? EdgeBoundary::BindOff : EdgeBoundary::Interior, EdgeBoundary::Interior});
        s.strip.push_back(0);
        s.slot.push_back(j);
    }
    int shapingStrip = 0;
    bool unsupported = false;
    for (int k = 1; k < R; ++k) {
        const EdgeBoundary topRole = k == R - 1 ? EdgeBoundary::BindOff : EdgeBoundary::Interior;
        if (!shapeStrip(k, ring[k], ring[k + 1], g.rings[k].length, g.rings[k + 1].length, r.w, r.stagger, EdgeBoundary::Interior, topRole, shapingStrip, s, g.shaping,
                        unsupported, err, r.adjacentShaping, r.boundaryShaping))
            return false;
    }
    if (r.rowPhase == RowPhase::WalesUpright) {
        std::vector<std::vector<int>> facesOf(R);
        for (int f = 0; f < s.mesh.faces(); ++f)
            if (s.strip[f] >= 0 && s.strip[f] < R) facesOf[s.strip[f]].push_back(f);
        double total = 0, largest = 0;
        for (int k = 1; k < R; ++k) {
            total += turnRingUpright(m, contours[k + 1], g.rings[k + 1].cutParameter, ring[k], ring[k + 1], facesOf[k], s);
            largest = std::max(largest, std::fabs(total));
        }
        g.rowStartOffPath = largest;
    }
    for (int k = 1; k < R; ++k) s.shiftPath.push_back(ShiftPathEdge{k, ring[k][0], ring[k + 1][0]});
    for (int f = 0; f < s.mesh.faces(); ++f)
        if (s.kind[f] == StitchFaceKind::Unsupported) unsupported = true;
    g.status = unsupported ? Status::Unresolved : Status::Pass;
    char b[240];
    std::snprintf(b, sizeof b, "%d rows, h' %.5f, pole start n1 %d (automatic %d), %d faces, %zu shaping slots, %s", R, g.hPrime, n1, out.ring1Auto, s.mesh.faces(),
                  g.shaping.size(), unsupported ? "unsupported slots listed" : "every slot typed");
    g.detail = b;
    if (r.rowPhase == RowPhase::WalesUpright) {
        std::snprintf(b, sizeof b, "; rings turned for upright wales, the row start up to %.1f mm off the path", g.rowStartOffPath);
        g.detail += b;
    }
    return true;
}

} // namespace knit
