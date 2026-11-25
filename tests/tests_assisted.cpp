// knit_tests_stitchmesh, assisted stockinette: the proposals for one piece
// (src/stitchmesh/mark_proposals.h) on the declared shapes (an open, a bent and a rotated tube, an isolated capsule, a
// branching piece). Each case names the defect it catches; the shapes are the
// C++ twins of the author's shape script, so keep both in step.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "harness.h"
#include "model/surface_path.h"
#include "model/surface_pieces.h"
#include "stitchmesh/mark_proposals.h"
#include "stitchmesh/marked_build.h"
#include "stitchmesh/quick_review.h"
#include "stitchmesh/guided_rows.h"
#include "stitchmesh/row_rings.h"
#include "stitchmesh/time_field.h"
#include "stitchmesh/prepared_reference.h"

using namespace knit;

namespace {

constexpr double kTau = 6.283185307179586;

// Rings of `segments` vertices (pos holds ring j's segment i at j * segments + i), quads split into two triangles.
TriMesh tubeGrid(int rings, int segments, const std::vector<Vec3>& pos) {
    TriMesh m;
    m.pos = pos;
    for (int j = 0; j + 1 < rings; ++j)
        for (int i = 0; i < segments; ++i) {
            const int a = j * segments + i, b = j * segments + (i + 1) % segments, c = a + segments, d = b + segments;
            m.tri.push_back({a, b, d});
            m.tri.push_back({a, d, c});
        }
    m.computeNormals();
    return m;
}

// The smallest prime >= n. The general fixtures keep a tube's ring count minus one prime, so their rows cut between rings;
// rows exactly on rings of vertices (a regular grid) have their own case, ssm_aligned_regular_tube.
int primeAtLeast(int n) {
    for (int p = std::max(2, n);; ++p) {
        bool prime = true;
        for (int d = 2; d * d <= p && prime; ++d) prime = p % d != 0;
        if (prime) return p;
    }
}

// An open tube round y from 0 to length, radius r0 at the bottom and r1 at the top.
TriMesh openTube(double r0, double r1, double length, int segments) {
    const int rings = primeAtLeast(static_cast<int>(std::lround(length / 1.5))) + 1;
    std::vector<Vec3> pos;
    for (int j = 0; j < rings; ++j) {
        const double f = double(j) / (rings - 1), r = r0 + (r1 - r0) * f;
        for (int i = 0; i < segments; ++i) {
            const double a = kTau * i / segments;
            pos.push_back(Vec3{r * std::cos(a), length * f, -r * std::sin(a)});
        }
    }
    return tubeGrid(rings, segments, pos);
}

// A tube of radius minor bent through 90 degrees round a centre `major` away (a torus segment in the x-y plane).
TriMesh bentTube(double major, double minor, int segments) {
    const double sweep = kTau / 4;
    const int rings = static_cast<int>(std::lround(major * sweep / 1.5)) + 1;
    std::vector<Vec3> pos;
    for (int j = 0; j < rings; ++j) {
        const double phi = sweep * j / (rings - 1);
        const Vec3 centre{major * std::cos(phi) - major, major * std::sin(phi), 0};
        const Vec3 normal{std::cos(phi), std::sin(phi), 0}, binormal{0, 0, -1};
        for (int i = 0; i < segments; ++i) {
            const double a = kTau * i / segments;
            pos.push_back(centre + minor * (std::cos(a) * normal + std::sin(a) * binormal));
        }
    }
    return tubeGrid(rings, segments, pos);
}

// A closed capsule round y: radius r, cylinder from 0 to length, hemispherical ends.
TriMesh capsuleY(double r, double length, int segments) {
    std::vector<double> ys, rs;
    const int cap = 12;
    for (int i = 1; i < cap; ++i) {
        const double a = -kTau / 4 + kTau / 4 * i / cap;
        ys.push_back(r * std::sin(a));
        rs.push_back(r * std::cos(a));
    }
    const int body = static_cast<int>(std::lround(length / 1.5));
    for (int i = 0; i <= body; ++i) {
        ys.push_back(length * i / body);
        rs.push_back(r);
    }
    for (int i = 1; i < cap; ++i) {
        const double a = kTau / 4 * i / cap;
        ys.push_back(length + r * std::sin(a));
        rs.push_back(r * std::cos(a));
    }
    TriMesh m;
    m.pos.push_back(Vec3{0, -r, 0});
    for (size_t j = 0; j < ys.size(); ++j)
        for (int i = 0; i < segments; ++i) {
            const double a = kTau * i / segments;
            m.pos.push_back(Vec3{rs[j] * std::cos(a), ys[j], -rs[j] * std::sin(a)});
        }
    m.pos.push_back(Vec3{0, length + r, 0});
    const int rings = static_cast<int>(ys.size()), top = static_cast<int>(m.pos.size()) - 1;
    auto v = [segments](int j, int i) { return 1 + j * segments + (i % segments); };
    for (int i = 0; i < segments; ++i) m.tri.push_back({0, v(0, i + 1), v(0, i)});
    for (int j = 0; j + 1 < rings; ++j)
        for (int i = 0; i < segments; ++i) {
            m.tri.push_back({v(j, i), v(j, i + 1), v(j + 1, i + 1)});
            m.tri.push_back({v(j, i), v(j + 1, i + 1), v(j + 1, i)});
        }
    for (int i = 0; i < segments; ++i) m.tri.push_back({v(rings - 1, i), v(rings - 1, i + 1), top});
    m.computeNormals();
    return m;
}

// Three lobes 120 degrees apart (an icosphere pushed out), each tip cut open: three openings (a pair of pants).
TriMesh branchingPiece(double radius) {
    const double g = (1 + std::sqrt(5.0)) / 2;
    std::vector<Vec3> pos = {{-1, g, 0}, {1, g, 0}, {-1, -g, 0}, {1, -g, 0}, {0, -1, g}, {0, 1, g}, {0, -1, -g}, {0, 1, -g}, {g, 0, -1}, {g, 0, 1}, {-g, 0, -1}, {-g, 0, 1}};
    for (Vec3& p : pos) p = normalized(p);
    std::vector<std::array<int, 3>> tri = {{0, 11, 5}, {0, 5, 1}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11}, {1, 5, 9}, {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
                                           {3, 9, 4}, {3, 4, 2}, {3, 2, 6}, {3, 6, 8}, {3, 8, 9}, {4, 9, 5}, {2, 4, 11}, {6, 2, 10}, {8, 6, 7}, {9, 8, 1}};
    for (int level = 0; level < 4; ++level) {
        std::map<std::pair<int, int>, int> mid;
        auto midpoint = [&](int a, int b) {
            const std::pair<int, int> k{std::min(a, b), std::max(a, b)};
            const auto it = mid.find(k);
            if (it != mid.end()) return it->second;
            pos.push_back(normalized(0.5 * (pos[a] + pos[b])));
            mid[k] = static_cast<int>(pos.size()) - 1;
            return mid[k];
        };
        std::vector<std::array<int, 3>> out;
        for (const std::array<int, 3>& t : tri) {
            const int ab = midpoint(t[0], t[1]), bc = midpoint(t[1], t[2]), ca = midpoint(t[2], t[0]);
            out.push_back({t[0], ab, ca});
            out.push_back({t[1], bc, ab});
            out.push_back({t[2], ca, bc});
            out.push_back({ab, bc, ca});
        }
        tri = out;
    }
    const Vec3 dirs[3] = {{0, 1, 0}, {-std::sin(kTau / 3), std::cos(kTau / 3), 0}, {std::sin(kTau / 3), std::cos(kTau / 3), 0}};
    const double hole = std::cos(24.0 / 360 * kTau);
    auto tip = [&](const Vec3& p) {
        for (const Vec3& d : dirs)
            if (dot(p, d) > hole) return true;
        return false;
    };
    std::vector<int> index(pos.size(), -1);
    TriMesh m;
    for (const std::array<int, 3>& t : tri) {
        if (tip(pos[t[0]]) && tip(pos[t[1]]) && tip(pos[t[2]])) continue;
        std::array<int, 3> c{};
        for (int k = 0; k < 3; ++k) {
            if (index[t[k]] < 0) {
                double s = 0;
                for (const Vec3& d : dirs) s = std::max(s, dot(pos[t[k]], d));
                index[t[k]] = static_cast<int>(m.pos.size());
                m.pos.push_back(radius * (1 + 2 * s * s * s) * pos[t[k]]);
            }
            c[k] = index[t[k]];
        }
        m.tri.push_back(c);
    }
    m.computeNormals();
    return m;
}

// m rotated by `degrees` about the axis and moved by `shift`.
TriMesh moved(const TriMesh& in, const Vec3& axis, double degrees, const Vec3& shift) {
    const Vec3 k = normalized(axis);
    const double th = degrees / 360 * kTau, c = std::cos(th), s = std::sin(th);
    TriMesh m = in;
    for (Vec3& p : m.pos) p = c * p + s * cross(k, p) + (1 - c) * dot(k, p) * k + shift;
    m.computeNormals();
    return m;
}

SurfacePieces onePiece(const TriMesh& m) { return findSurfacePieces(m.pos, m.tri, std::vector<int>(m.tri.size(), 0), 1e-9); }

ProposalSet propose(const SurfacePieces& pieces, double mmPerUnit = 1) {
    ProposalRequest r;
    r.pieces = &pieces;
    r.piece = 0;
    r.fileDigest = 42;
    r.mmPerUnit = mmPerUnit;
    ProposalSet out;
    std::string err;
    CHECK(proposeMarkings(r, out, &err));
    std::printf("  %s: %zu setups, %zu findings\n", out.shapeText.c_str(), out.setups.size(), out.findings.size());
    for (const MarkingFinding& f : out.findings) std::printf("    finding: %s\n", f.message.c_str());
    for (const std::string& l : out.log) std::printf("    %s\n", l.c_str());
    return out;
}

// Whether a file point lies on a border edge of its piece: two weights 0.5 on corners whose edge only one triangle has.
bool onBorder(const SurfacePieces& pieces, const SurfacePoint& p) {
    int piece = -1;
    SurfacePoint lp;
    if (!pieceLocalPoint(pieces, p, piece, lp)) return false;
    const TriMesh& m = pieces.pieces[piece].mesh;
    int zero = -1;
    for (int k = 0; k < 3; ++k)
        if (lp.bary[k] == 0) zero = k;
    if (zero < 0) return false;
    const int a = m.tri[lp.tri][(zero + 1) % 3], b = m.tri[lp.tri][(zero + 2) % 3];
    int uses = 0;
    for (const std::array<int, 3>& t : m.tri)
        for (int k = 0; k < 3; ++k)
            if ((t[k] == a && t[(k + 1) % 3] == b) || (t[k] == b && t[(k + 1) % 3] == a)) ++uses;
    return uses == 1;
}

Vec3 position(const SurfacePieces& pieces, const SurfacePoint& p) {
    int piece = -1;
    SurfacePoint lp;
    pieceLocalPoint(pieces, p, piece, lp);
    return surfacePosition(pieces.pieces[piece].mesh, lp);
}

const SurfaceMarking* withMeaning(const ProposalSetup& s, MarkingMeaning meaning) {
    for (const SurfaceMarking& m : s.markings)
        if (m.meaning == meaning) return &m;
    return nullptr;
}

} // namespace

// A cylinder round y (radius 10, height 20) closed by a flat top with a small hole (radius 2) off its centre: an open
// tube whose top opening lies far from every row of the bottom's field (its near side is 2 mm from the rim, its far side
// 6 mm). From the bottom, the last row is the whole top circle (about 10 stitches) and the hole has about 2.
TriMesh cylinderWithOffsetHole() {
    const int segments = 48, wall = 14, top = 8;
    const double cx = 6, cz = 0, hole = 2, radius = 10, height = 20;
    std::vector<Vec3> outer, inner;
    for (int i = 0; i < segments; ++i) {
        const double a = kTau * i / segments, ux = std::cos(a), uz = -std::sin(a);
        const double cu = cx * ux + cz * uz, s = -cu + std::sqrt(cu * cu - (cx * cx + cz * cz - radius * radius));
        outer.push_back(Vec3{cx + s * ux, 0, cz + s * uz});
        inner.push_back(Vec3{cx + hole * ux, height, cz + hole * uz});
    }
    std::vector<Vec3> pos;
    for (int j = 0; j < wall; ++j)
        for (int i = 0; i < segments; ++i) pos.push_back(Vec3{outer[i].x, height * j / (wall - 1), outer[i].z});
    for (int k = 1; k <= top; ++k)
        for (int i = 0; i < segments; ++i) {
            const double f = double(k) / top;
            pos.push_back(Vec3{outer[i].x + (inner[i].x - outer[i].x) * f, height, outer[i].z + (inner[i].z - outer[i].z) * f});
        }
    return tubeGrid(wall + top, segments, pos);
}

// An end far from a row: from the bottom, the last strip would go from about 10 stitches to 2, more shaping than a row
// holds. The layout refuses it with the row, the counts and what to do (it used to stop on an internal assembly error),
// and the suggestions no longer offer that direction, saying why. Control: the other direction (from the hole, whose
// far end, the bottom, is close to a row) is still offered and lays out.
TEST(ssm_far_finish_named_and_not_offered) {
    const TriMesh m = cylinderWithOffsetHole();
    const int segments = 48, rings = 22;
    std::vector<int> bottom, hole;
    for (int i = 0; i < segments; ++i) {
        bottom.push_back(i);
        hole.push_back((rings - 1) * segments + i);
    }
    std::vector<double> t;
    TimeFieldReport tr;
    std::string err;
    CHECK(geodesicDistance(m, bottom, t, tr, &err));
    GuidedRowsRequest r;
    r.cutOrigin = Vec3{6, 0, 0};   // the half-plane z = 0, x > 6 (from the hole centre outwards) starts every row once
    r.cutNormal = Vec3{0, 0, 1};
    r.cutSide = Vec3{1, 0, 0};
    GuidedRowsProduct p;
    CHECK(!buildGuidedRows(m, t, bottom, hole, r, p, &err));
    std::printf("  layout from the bottom: %s\n", err.c_str());
    CHECK(err.find("far from a row") != std::string::npos && err.find("assembleStrip") == std::string::npos);
    const SurfacePieces pieces = onePiece(m);
    const ProposalSet out = propose(pieces);
    CHECK(out.shape == PieceShape::OpenTube && !out.setups.empty());
    bool named = false;
    for (const MarkingFinding& f : out.findings) named = named || f.message.find("lies far from a row") != std::string::npos;
    CHECK(named);
    for (const ProposalSetup& s : out.setups) CHECK(s.firstRingStitches < s.lastRingStitches);   // every setup starts at the hole
}


// An open tube gets two setups, one per direction: the start and the finish are the tube's own border loops (every
// control point on a border edge), the path runs from one to the other, nothing is left out, and the second setup is
// the first reversed. Every marking carries its setup and the request's basis. (Defect: an open tube treated as a
// closed piece, with a manufactured pole start and a left-out cap.)
TEST(ssm_proposals_open_tube) {
    const TriMesh m = openTube(10.5, 10.5, 40, 64);
    const SurfacePieces pieces = onePiece(m);
    const ProposalSet ps = propose(pieces);
    CHECK(ps.shape == PieceShape::OpenTube && ps.setups.size() == 2 && ps.findings.empty());
    CHECK(ps.basis == proposalBasis(42, 0, "facets", 1, 5.5, 3.4));
    for (const ProposalSetup& s : ps.setups) {
        const SurfaceMarking* start = withMeaning(s, MarkingMeaning::Start);
        const SurfaceMarking* finish = withMeaning(s, MarkingMeaning::Finish);
        const SurfaceMarking* path = withMeaning(s, MarkingMeaning::RowStartPath);
        CHECK(start && finish && path && start->shape == MarkingShape::ClosedCurve && finish->shape == MarkingShape::ClosedCurve);
        for (const SurfacePoint& p : start->points) CHECK(onBorder(pieces, p));
        for (const SurfacePoint& p : finish->points) CHECK(onBorder(pieces, p));
        const double startY = position(pieces, start->points[0]).y, finishY = position(pieces, finish->points[0]).y;
        CHECK(std::fabs(std::fabs(startY - finishY) - 40) < 1e-9);
        CHECK(std::fabs(position(pieces, path->points.front()).y - startY) < 3.4 && std::fabs(position(pieces, path->points.back()).y - finishY) < 3.4);
        CHECK(s.omittedVisibleMm2 == 0 && s.omittedHiddenMm2 == 0 && !s.partial && !s.refinable && s.rowsAbout == 12 && s.firstRingStitches == 12);
        CHECK(!s.title.empty() && !s.purpose.empty() && !s.reason.empty() && !s.limitation.empty());
        for (const SurfaceMarking& mk : s.markings) CHECK(mk.setup == s.setup && mk.basis == ps.basis && mk.proposal == ProposalState::Offered);
    }
    CHECK(position(pieces, withMeaning(ps.setups[0], MarkingMeaning::Start)->points[0]).y == position(pieces, withMeaning(ps.setups[1], MarkingMeaning::Finish)->points[0]).y);
}

// A bent tube gets an inner and an outer path per direction (4 setups); the inner path runs nearer the bend's centre than
// the outer one. A straight tube gets one path per direction. (Defect: the inner / outer choice inverted, or offered
// where it means nothing.)
TEST(ssm_proposals_bent_inner_outer) {
    const TriMesh m = bentTube(30, 10.5, 64);
    const SurfacePieces pieces = onePiece(m);
    const ProposalSet ps = propose(pieces);
    CHECK(ps.shape == PieceShape::OpenTube && ps.setups.size() == 4);
    const Vec3 bendCentre{-30, 0, 0};
    auto meanDistance = [&](const ProposalSetup& s) {
        const SurfaceMarking* path = withMeaning(s, MarkingMeaning::RowStartPath);
        double sum = 0;
        for (const SurfacePoint& p : path->points) sum += norm(position(pieces, p) - bendCentre);
        return sum / double(path->points.size());
    };
    for (int d = 0; d < 2; ++d) {
        const ProposalSetup& inner = ps.setups[2 * d];
        const ProposalSetup& outer = ps.setups[2 * d + 1];
        CHECK(inner.recipe.find("inner") != std::string::npos && outer.recipe.find("outer") != std::string::npos);
        std::printf("  direction %d: inner path %.2f mm, outer path %.2f mm from the bend's centre\n", d, meanDistance(inner), meanDistance(outer));
        CHECK(meanDistance(inner) < 30 - 5 && meanDistance(outer) > 30 + 5);
    }
    const ProposalSet straight = propose(onePiece(openTube(10.5, 10.5, 40, 64)));
    CHECK(straight.setups.size() == 2 && straight.setups[0].recipe.find("any") != std::string::npos);
}

// The same tube rotated 37 degrees and moved gives the same setups: titles, recipes, counts and the same file triangles,
// with every control point at the moved position (1e-6 mm). (Defect: a rule that reads a world axis or coordinate, e.g.
// "the lower opening" by y.)
TEST(ssm_proposals_transform_invariant) {
    const TriMesh a = openTube(12, 7, 40, 64);
    const Vec3 axis{1, 2, 3}, shift{13, -7, 21};
    const TriMesh b = moved(a, axis, 37, shift);
    const SurfacePieces pa = onePiece(a), pb = onePiece(b);
    const ProposalSet sa = propose(pa), sb = propose(pb);
    CHECK(sa.shape == sb.shape && sa.setups.size() == sb.setups.size() && sa.setups.size() == 2);
    double worst = 0;
    for (size_t i = 0; i < sa.setups.size() && i < sb.setups.size(); ++i) {
        const ProposalSetup &x = sa.setups[i], &y = sb.setups[i];
        CHECK(x.title == y.title && x.recipe == y.recipe && x.rowsAbout == y.rowsAbout && x.firstRingStitches == y.firstRingStitches && x.lastRingStitches == y.lastRingStitches);
        CHECK(x.markings.size() == y.markings.size());
        for (size_t k = 0; k < x.markings.size() && k < y.markings.size(); ++k) {
            CHECK(x.markings[k].points.size() == y.markings[k].points.size());
            for (size_t j = 0; j < x.markings[k].points.size() && j < y.markings[k].points.size(); ++j) {
                CHECK(x.markings[k].points[j].tri == y.markings[k].points[j].tri);
                const double th = 37.0 / 360 * kTau, c = std::cos(th), s = std::sin(th);
                const Vec3 k3 = normalized(axis), p = position(pa, x.markings[k].points[j]);
                const Vec3 q = c * p + s * cross(k3, p) + (1 - c) * dot(k3, p) * k3 + shift;
                worst = std::max(worst, norm(q - position(pb, y.markings[k].points[j])));
            }
        }
    }
    std::printf("  largest control point difference after the motion: %.3g mm\n", worst);
    CHECK(worst < 1e-6);
    CHECK(sa.setups[0].title.find("wider") != std::string::npos);   // named by size, not by position
}

// A closed capsule that nothing hides: every setup is partial and discloses the visible cap it leaves out (> 0 mm2), every
// finish is a curve (no invented closed finish), and a finding says why and invites the user's own finish. (Defect: a
// setup claiming to cover the whole closed piece, or a hidden visible cap.)
TEST(ssm_proposals_capsule_partial) {
    const TriMesh m = capsuleY(10.5, 24, 64);
    const SurfacePieces pieces = onePiece(m);
    const ProposalSet ps = propose(pieces);
    CHECK(ps.shape == PieceShape::ClosedIsolated && !ps.setups.empty());
    for (const ProposalSetup& s : ps.setups) {
        CHECK(s.partial && s.omittedVisibleMm2 > 100);
        CHECK(withMeaning(s, MarkingMeaning::Finish) && withMeaning(s, MarkingMeaning::Finish)->shape == MarkingShape::ClosedCurve);
        CHECK(s.limitation.find("mm2") != std::string::npos);
    }
    bool invited = false;
    for (const MarkingFinding& f : ps.findings) invited = invited || (f.cause == FindingCause::AlgorithmLimitation && f.message.find("closed finish is not available") != std::string::npos);
    CHECK(invited);
}

// A branching piece (three openings): refused with a finding at the split (its row and the split contours' positions),
// as a missing construction (not the user's marking), with one partial setup up to the split. (Defect: a branching
// piece given complete setups, or the refusal blamed on the markings.)
TEST(ssm_proposals_branch_refused) {
    const TriMesh m = branchingPiece(12);
    const SurfacePieces pieces = onePiece(m);
    CHECK(borderLoops(pieces.pieces[0].mesh).size() == 3);
    const ProposalSet ps = propose(pieces);
    CHECK(ps.shape == PieceShape::Branching);
    CHECK(ps.findings.size() == 1 && ps.findings[0].level == FindingLevel::Problem && ps.findings[0].cause == FindingCause::AlgorithmLimitation);
    CHECK(ps.findings[0].row > 1 && ps.findings[0].at.size() >= 8);
    CHECK(ps.setups.size() == 1 && ps.setups[0].partial && ps.setups[0].omittedVisibleMm2 > 100);
    for (const SurfacePoint& p : withMeaning(ps.setups[0], MarkingMeaning::Start)->points) CHECK(onBorder(pieces, p));
}

namespace {

// A set holding one setup's markings, offered, and a build request that previews exactly them (a research preview).
MarkingSet setupSet(const ProposalSetup& s) {
    MarkingSet set;
    set.digest = 42;
    set.markings = s.markings;
    for (const SurfaceMarking& m : set.markings) set.nextId = std::max(set.nextId, m.id + 1);
    return set;
}

MarkedBuildProduct buildSetup(const SurfacePieces& pieces, const MarkingSet& set, ReviewDepth depth) {
    MarkedBuildRequest rq;
    rq.pieces = &pieces;
    rq.piece = 0;
    rq.markings = &set;
    rq.w = 5.5;
    rq.h = 3.4;
    rq.depth = depth;
    for (const SurfaceMarking& m : set.markings)
        if (m.proposal == ProposalState::Offered) rq.previewOffered.push_back(m.id);
    MarkedBuildProduct p;
    std::string err;
    CHECK(buildFromMarkings(rq, p, &err));
    for (const MarkingFinding& f : p.review.findings)
        std::printf("  finding on %d, level %d, cause %d, row %d: %s\n", f.markingId, static_cast<int>(f.level), static_cast<int>(f.cause), f.row, f.message.c_str());
    return p;
}

bool hasMessage(const MarkingReview& r, const std::string& text) {
    for (const MarkingFinding& f : r.findings)
        if (f.message.find(text) != std::string::npos) return true;
    return false;
}

} // namespace

// An open tube's proposed setup builds with its real openings as cast-on and bind-off: nothing is trimmed or left out,
// the rows are round(T / h), the path crosses every row once, the layout keeps its field (band route) and its measures say
// what was computed. With the opening as start and a drawn curve as finish, only the curve trims: what lies beyond it is
// left out (visible) and named "beyond the finish". (Defect: the border points traced as a plane section, which leaves the
// surface and refuses, or a sliver trimmed and left out at the opening.)
TEST(ssm_marked_open_tube_border_loops) {
    const TriMesh m = openTube(10.5, 10.5, 40, 64);
    const SurfacePieces pieces = onePiece(m);
    const ProposalSet ps = propose(pieces);
    CHECK(ps.setups.size() == 2);
    const MarkingSet set = setupSet(ps.setups[0]);
    const MarkedBuildProduct p = buildSetup(pieces, set, ReviewDepth::Built);
    const LayoutMeasures& ms = p.review.measures;
    CHECK(p.review.valid && ms.omittedAreaMm2 == 0 && ms.rows == 12 && ms.contourEvents == 0 && ms.unsupportedCells == 0);
    CHECK(hasMessage(p.review, "crosses every row once") && ms.coveragePass && ms.checksPass && ms.checksKnown == Known::Current);
    CHECK(ms.rowsKnown == Known::Current && ms.coverageKnown == Known::Current);
    CHECK(p.field && p.field->route == "band" && p.field->rows == 12 && static_cast<int>(p.field->contours.size()) == 13);
    for (MarkedRegion g : p.region) CHECK(g == MarkedRegion::Kept);
    // the opening as start, a drawn curve 10 mm below the other opening as finish
    MarkingSet curve = set;
    for (SurfaceMarking& mk : curve.markings) {
        if (mk.meaning != MarkingMeaning::Finish) continue;
        const double y = position(pieces, withMeaning(ps.setups[0], MarkingMeaning::Start)->points[0]).y < 1 ? 30.0 : 10.0;
        mk.points.clear();
        for (int i = 0; i < 6; ++i) {
            const double a = kTau * i / 6 + 0.05;
            mk.points.push_back(closestSurfacePoint(m, Vec3{10.5 * std::cos(a), y, -10.5 * std::sin(a)}));
        }
    }
    const MarkedBuildProduct q = buildSetup(pieces, curve, ReviewDepth::Built);
    std::printf("  opening + curve: rows %d, omitted %.1f mm2 (visible %.1f), cylinder beyond the curve %.1f mm2\n", q.review.measures.rows, q.review.measures.omittedAreaMm2,
                q.review.measures.omittedVisibleMm2, kTau * 10.5 * 10);
    CHECK(q.review.valid && q.review.measures.rows == 9);
    CHECK(std::fabs(q.review.measures.omittedVisibleMm2 - kTau * 10.5 * 10) < 0.03 * kTau * 10.5 * 10);
    int beyond = 0;
    for (MarkedRegion g : q.region) beyond += g == MarkedRegion::OmittedFinish ? 1 : 0;
    CHECK(beyond > 0);
    for (char hidden : q.regionHidden) CHECK(hidden == 0);   // nothing hides a lone tube
}

namespace {

// A cup round y: an open tube of radius r from y = 0 to `length` (the opening at y = 0), closed at the top by a half
// sphere ending in a pole vertex (capsuleY without its bottom cap).
TriMesh cupY(double r, double length, int segments) {
    std::vector<double> ys, rs;
    const int body = static_cast<int>(std::lround(length / 1.5));
    for (int i = 0; i <= body; ++i) {
        ys.push_back(length * i / body);
        rs.push_back(r);
    }
    const int cap = 12;
    for (int i = 1; i < cap; ++i) {
        const double a = kTau / 4 * i / cap;
        ys.push_back(length + r * std::sin(a));
        rs.push_back(r * std::cos(a));
    }
    std::vector<Vec3> pos;
    for (size_t j = 0; j < ys.size(); ++j)
        for (int i = 0; i < segments; ++i) {
            const double a = kTau * i / segments;
            pos.push_back(Vec3{rs[j] * std::cos(a), ys[j], -rs[j] * std::sin(a)});
        }
    const int rings = static_cast<int>(ys.size());
    TriMesh m = tubeGrid(rings, segments, pos);
    m.pos.push_back(Vec3{0, length + r, 0});
    const int top = static_cast<int>(m.pos.size()) - 1;
    for (int i = 0; i < segments; ++i) m.tri.push_back({(rings - 1) * segments + i, (rings - 1) * segments + (i + 1) % segments, top});
    m.computeNormals();
    return m;
}

} // namespace

// A piece with one opening (a cup: an ear, a head, a sleeve with a closed end) gets one setup: a pole start at the
// cup's bottom (on a round cup the pole) and the opening itself as the bind-off. It lays out with nothing left out, no contour event,
// every cell supported and the path crossing every row once. (Defect: a cup got no suggestion, and a start point with the
// opening as its finish was refused, so the user had to draw a finish curve inside the rim and leave a strip unknitted.)
// Control: a start point with an opening as the finish on an open tube (two openings) is refused by name.
TEST(ssm_cup_pole_start_opening_finish) {
    const TriMesh m = cupY(10.5, 30, 64);
    const SurfacePieces pieces = onePiece(m);
    const ProposalSet ps = propose(pieces);
    CHECK(ps.shape == PieceShape::Cup && ps.setups.size() == 1);
    CHECK(position(pieces, withMeaning(ps.setups[0], MarkingMeaning::Start)->points[0]).y > 38);   // at the cup's bottom (the pole, y 40.5)
    const MarkedBuildProduct p = buildSetup(pieces, setupSet(ps.setups[0]), ReviewDepth::Built);
    const LayoutMeasures& ms = p.review.measures;
    std::printf("  cup: rows %d, events %d, unsupported %d, omitted %.1f mm2, covered %.4f\n", ms.rows, ms.contourEvents, ms.unsupportedCells,
                ms.omittedAreaMm2, ms.coveredFraction);
    CHECK(p.review.valid && ms.omittedAreaMm2 == 0 && ms.contourEvents == 0 && ms.unsupportedCells == 0);
    CHECK(p.field && p.field->route == "pole" && hasMessage(p.review, "crosses every row once"));
    for (MarkedRegion g : p.region) CHECK(g == MarkedRegion::Kept);
    // control: an open tube's setup with its start curve replaced by a start point in the middle of the tube
    const TriMesh tube = openTube(10.5, 10.5, 40, 64);
    const SurfacePieces tubePieces = onePiece(tube);
    MarkingSet bad = setupSet(propose(tubePieces).setups[0]);
    for (SurfaceMarking& mk : bad.markings)
        if (mk.meaning == MarkingMeaning::Start) {
            mk.shape = MarkingShape::Point;
            mk.points = {closestSurfacePoint(tube, Vec3{10.5, 20, 0})};
        }
    const MarkedBuildProduct q = buildSetup(tubePieces, bad, ReviewDepth::Built);
    CHECK(!q.review.valid && hasMessage(q.review, "exactly one opening"));
}

// A cup whose opening is not at one distance from its start (an oval dish, 1.6 times as long as wide; the inner side of a
// teddy ear cut along its rim). With rows at levels of the distance from the
// start (OpeningRows::DistanceFromStart) they end at the nearest point of the opening and the last strip takes the rest:
// its tallest wale edge is several stitch heights. With the opening as a level (OpeningIsALevel) no strip is that tall and
// more of the dish is covered before any relaxation.
TEST(ssm_cup_opening_is_a_level_of_the_rows) {
    TriMesh m = cupY(20, 3, 64);
    for (Vec3& q : m.pos) q.x *= 1.6;
    m.computeNormals();
    const SurfacePieces pieces = onePiece(m);
    const ProposalSet ps = propose(pieces);
    CHECK(ps.shape == PieceShape::Cup && ps.setups.size() == 1);
    if (ps.setups.empty()) return;
    const MarkingSet set = setupSet(ps.setups[0]);
    auto tallestWale = [](const TypedStitchMesh& t) {
        double tallest = 0;
        for (int f = 0; f < t.mesh.faces(); ++f)
            for (int k = 0; k < t.mesh.sides(f); ++k)
                if (t.role[f][k] == StitchEdgeRole::Wale) tallest = std::max(tallest, norm(t.mesh.pos[t.mesh.vertexAt(f, k)] - t.mesh.pos[t.mesh.vertexAt(f, (k + 1) % t.mesh.sides(f))]));
        return tallest;
    };
    MarkedBuildRequest rq;
    rq.pieces = &pieces;
    rq.piece = 0;
    rq.markings = &set;
    rq.w = 5.5;
    rq.h = 3.4;
    rq.depth = ReviewDepth::Built;
    for (const SurfaceMarking& mk : set.markings)
        if (mk.proposal == ProposalState::Offered) rq.previewOffered.push_back(mk.id);
    MarkedBuildProduct fromStart, level;
    std::string err;
    rq.options.openingRows = OpeningRows::DistanceFromStart;
    CHECK(buildFromMarkings(rq, fromStart, &err));
    rq.options.openingRows = OpeningRows::OpeningIsALevel;
    CHECK(buildFromMarkings(rq, level, &err));
    std::printf("  oval dish: distance from the start: built %d, rows %d, tallest wale edge %.2f mm, covered %.4f; opening as a level: built %d, rows %d, %.2f mm, covered %.4f\n",
                int(fromStart.review.valid), fromStart.review.measures.rows, tallestWale(fromStart.evaluation.before), fromStart.review.measures.coveredFraction,
                int(level.review.valid), level.review.measures.rows, tallestWale(level.evaluation.before), level.review.measures.coveredFraction);
    CHECK(level.review.valid && level.review.measures.contourEvents == 0 && level.review.measures.unsupportedCells == 0);
    CHECK(tallestWale(level.evaluation.before) < 2.0 * rq.h);
    CHECK(!fromStart.review.valid || tallestWale(fromStart.evaluation.before) > 2.5 * rq.h);
    CHECK(level.review.measures.coveredFraction > fromStart.review.measures.coveredFraction);
}

namespace {

// A regular tube r 10.5 x 40 mm with 64 segments and `rings` evenly spaced rings (no prime trick): with rings - 1 a
// multiple of 12, every row level of h = 3.4 (R = 12) falls exactly on a ring of vertices.
TriMesh regularTube(int rings) {
    std::vector<Vec3> pos;
    for (int j = 0; j < rings; ++j)
        for (int i = 0; i < 64; ++i) {
            const double a = kTau * i / 64;
            pos.push_back(Vec3{10.5 * std::cos(a), 40.0 * j / (rings - 1), -10.5 * std::sin(a)});
        }
    return tubeGrid(rings, 64, pos);
}

} // namespace

// A level through a ring of vertices is one closed contour: with t = y exactly (the level equals the ring's values) and
// with rounding noise of either sign on the ring (some of its vertices just above, some just below). (Defect: segments
// oriented by a left/right test on positions, which has no sign when crossings coincide with a vertex, so the chain
// breaks into several contours and a regular tube is refused as "row 1 splits".)
TEST(ssm_row_level_through_vertices) {
    const TriMesh m = regularTube(13);
    std::vector<double> exact(m.pos.size()), noisy(m.pos.size());
    for (size_t v = 0; v < m.pos.size(); ++v) {
        exact[v] = m.pos[v].y;
        noisy[v] = m.pos[v].y + ((v * 7919) % 3 == 0 ? 1e-13 : (v % 2 ? -1e-13 : 0.0));
    }
    for (int ring = 1; ring < 12; ++ring) {
        const double level = 40.0 * ring / 12;
        for (const std::vector<double>* t : {&exact, &noisy}) {
            const RingLevel lv = extractLevel(m, *t, level);
            CHECK(lv.contours.size() == 1 && lv.contours.front().closed);
            if (lv.contours.size() == 1) CHECK(std::fabs(lv.contours.front().length - kTau * 10.5) < 0.01 * kTau * 10.5);
        }
    }
}

namespace {

// The quick review of a set against a build's field, on the same prepared piece the build used.
MarkingReview quick(const SurfacePieces& pieces, const MarkingSet& set, const MarkedBuildProduct& built, std::vector<MarkingTrace>* traces = nullptr) {
    QuickRequest q;
    q.reference = built.field ? built.field->reference : nullptr;
    if (!q.reference) {
        std::shared_ptr<PreparedPiece> p = std::make_shared<PreparedPiece>();
        preparePiece(pieces, 0, set.digest, 1, 3.4, ReferencePreparation::Facets, *p);
        q.reference = p;
    }
    q.field = built.field;
    q.markings = &set;
    q.pieces = &pieces;
    q.piece = 0;
    for (const SurfaceMarking& m : set.markings)
        if (m.proposal == ProposalState::Offered) q.previewOffered.push_back(m.id);
    std::vector<int> ids;
    return quickReview(q, &ids, traces);
}

// The crossing findings of a review (rows crossed, missed, split, turn-backs, "every row once"), as text with row and
// place count, in order.
std::vector<std::string> crossingTexts(const MarkingReview& r) {
    std::vector<std::string> out;
    for (const MarkingFinding& f : r.findings) {
        const bool crossing = f.message.find("row") != std::string::npos && f.message.find("Research preview") == std::string::npos &&
                              f.message.find("Row checks") == std::string::npos && f.message.find("pole start") == std::string::npos;
        if (crossing) out.push_back(std::to_string(f.row) + "|" + std::to_string(f.at.size()) + "|" + f.message);
    }
    return out;
}

// The path marking of a set replaced by control points at these positions (on the piece, file units).
MarkingSet withPath(const MarkingSet& set, const TriMesh& m, const std::vector<Vec3>& points) {
    MarkingSet out = set;
    for (SurfaceMarking& mk : out.markings) {
        if (mk.meaning != MarkingMeaning::RowStartPath) continue;
        mk.points.clear();
        for (const Vec3& p : points) mk.points.push_back(closestSurfacePoint(m, p));
    }
    return out;
}

Vec3 onTube(double r, double angle, double y) { return Vec3{r * std::cos(angle), y, -r * std::sin(angle)}; }

} // namespace

// Regular tubes whose rows fall exactly on rings of vertices (13 and 25 rings) are treated as the unaligned one (30
// rings): two setups, a Built layout of 144 cells with full coverage and its checks passing, the path crossing every row
// once, and the quick review on that field saying the same. (Defect: the row contour's orientation or the path's crossing
// decided per triangle from positions, which loses or splits rows through vertices: "row 1 splits", or "the row-start
// path misses rows 3, 5-7, 10".)
TEST(ssm_aligned_regular_tube) {
    for (int rings : {13, 25, 30}) {
        const TriMesh m = regularTube(rings);
        const SurfacePieces pieces = onePiece(m);
        const ProposalSet ps = propose(pieces);
        CHECK(ps.shape == PieceShape::OpenTube && ps.setups.size() == 2 && ps.findings.empty());
        if (ps.setups.size() != 2) continue;
        for (const ProposalSetup& s : ps.setups) {
            const MarkingSet set = setupSet(s);
            const MarkedBuildProduct p = buildSetup(pieces, set, ReviewDepth::Built);
            const LayoutMeasures& ms = p.review.measures;
            std::printf("  %d rings, setup %d: valid %d, rows %d, faces %d, covered %.4f, checks %d\n", rings, s.setup, p.review.valid ? 1 : 0, ms.rows,
                        p.stitchMesh.mesh.faces(), ms.coveredFraction, ms.checksPass ? 1 : 0);
            CHECK(p.review.valid && ms.rows == 12 && ms.contourEvents == 0 && p.stitchMesh.mesh.faces() == 144 && ms.coveragePass && ms.checksPass);
            CHECK(hasMessage(p.review, "crosses every row once"));
            CHECK(crossingTexts(quick(pieces, set, p)) == std::vector<std::string>{"-1|0|This path crosses every row once."});
        }
    }
}

// On a current field the quick review names the row a path crosses more than once, with its places and the trace range
// between them, the rows a short path misses, and where a path turns back; a straight path crosses every row once.
// (Defect: the row checks skipped while drawing, or a crossing reported without where it is.)
TEST(ssm_quick_review_row_crossings) {
    const TriMesh m = openTube(10.5, 10.5, 40, 64);
    const SurfacePieces pieces = onePiece(m);
    const ProposalSet ps = propose(pieces);
    const MarkingSet set = setupSet(ps.setups[0]);
    const MarkedBuildProduct built = buildSetup(pieces, set, ReviewDepth::Built);
    CHECK(built.field && built.field->rows == 12);
    const bool up = position(pieces, withMeaning(ps.setups[0], MarkingMeaning::Start)->points[0]).y < 1;
    auto y = [up](double fromStart) { return up ? fromStart : 40 - fromStart; };
    // straight up the tube: every row once
    const MarkingSet straight = withPath(set, m, {onTube(10.5, 0.3, y(1)), onTube(10.5, 0.3, y(20)), onTube(10.5, 0.3, y(39))});
    CHECK(crossingTexts(quick(pieces, straight, built)) == std::vector<std::string>{"-1|0|This path crosses every row once."});
    // up past row 6 (20 mm), back below it, up again: row 6 crossed three times
    const MarkingSet zigzag = withPath(set, m, {onTube(10.5, 0.3, y(1)), onTube(10.5, 0.5, y(21.5)), onTube(10.5, 0.7, y(18.5)), onTube(10.5, 0.9, y(39))});
    const MarkingReview z = quick(pieces, zigzag, built);
    for (const std::string& s : crossingTexts(z)) std::printf("  zigzag: %s\n", s.c_str());
    bool row6 = false, turns = false;
    for (const MarkingFinding& f : z.findings) {
        if (f.row == 6 && f.message.find("crosses row 6 3 times") != std::string::npos) {
            row6 = f.at.size() == 3 && f.firstTracePoint >= 0 && f.lastTracePoint > f.firstTracePoint;
            for (const Vec3& p : f.at) CHECK(std::fabs(std::fabs(p.y - y(0)) - 20) < 0.5);   // on row 6 (h' = 40 / 12 mm)
        }
        turns = turns || f.message.find("turns back") != std::string::npos;
    }
    CHECK(row6 && turns && z.measures.rowsKnown == Known::Current && z.measures.contourEvents > 0);
    // a path that stops half-way misses the rows above
    const MarkingSet half = withPath(set, m, {onTube(10.5, 0.3, y(1)), onTube(10.5, 0.3, y(18))});
    bool missed = false;
    for (const MarkingFinding& f : quick(pieces, half, built).findings) missed = missed || (f.message.find("misses rows 6-11") != std::string::npos && f.row == 6);
    CHECK(missed);
}

// The quick review counts with the layout's own code: on the band route (the tube's openings, and drawn curves on a
// capsule) and the pole route, for a straight, a double-crossing and a too-short path, its crossing findings equal those of
// a full build of the same markings. (Defect: a quick count of its own that disagrees with the layout.)
TEST(ssm_quick_review_matches_layout) {
    struct Case {
        std::string name;
        TriMesh mesh;
        int setup;
    };
    const std::vector<Case> cases = {{"tube", openTube(10.5, 10.5, 40, 64), 0}, {"capsule band", capsuleY(10.5, 24, 64), 1}, {"capsule pole", capsuleY(10.5, 24, 64), 0}};
    for (const Case& c : cases) {
        const SurfacePieces pieces = onePiece(c.mesh);
        const ProposalSet ps = propose(pieces);
        CHECK(static_cast<int>(ps.setups.size()) > c.setup);
        if (static_cast<int>(ps.setups.size()) <= c.setup) continue;
        const MarkingSet set = setupSet(ps.setups[c.setup]);
        const MarkedBuildProduct base = buildSetup(pieces, set, ReviewDepth::Built);
        CHECK(base.field != nullptr);
        Vec3 lo, hi;
        c.mesh.bounds(lo, hi);
        const double r = 10.5, y0 = lo.y + 1, y1 = hi.y - 1, mid = 0.5 * (lo.y + hi.y);
        const std::vector<std::vector<Vec3>> paths = {
            {onTube(r, 0.4, y0 + 2), onTube(r, 0.4, mid), onTube(r, 0.4, y1 - 2)},
            {onTube(r, 0.4, y0 + 2), onTube(r, 0.6, mid + 4), onTube(r, 0.8, mid - 4), onTube(r, 1.0, y1 - 2)},
            {onTube(r, 0.4, y0 + 2), onTube(r, 0.4, mid)},
        };
        for (const std::vector<Vec3>& p : paths) {
            const MarkingSet variant = withPath(set, c.mesh, p);
            const MarkedBuildProduct full = buildSetup(pieces, variant, ReviewDepth::Built);
            const std::vector<std::string> a = crossingTexts(full.review), b = crossingTexts(quick(pieces, variant, base));
            std::printf("  %s: %zu crossing findings in the build, %zu in the quick review\n", c.name.c_str(), a.size(), b.size());
            CHECK(a == b);
            CHECK(!b.empty());
        }
    }
}

// An unfinished path is Incomplete (not a problem, and no "add a path" beside it), in the quick review and in a build
// (which still builds nothing); a path that loops over itself is named where it crosses. (Defect: a stroke being drawn
// reported as an error, or a self-crossing path passed.)
TEST(ssm_quick_review_incomplete_and_self_crossing) {
    const TriMesh m = openTube(10.5, 10.5, 40, 64);
    const SurfacePieces pieces = onePiece(m);
    const ProposalSet ps = propose(pieces);
    const MarkingSet set = setupSet(ps.setups[0]);
    const MarkedBuildProduct built = buildSetup(pieces, set, ReviewDepth::Built);
    MarkingSet unfinished = set;
    for (SurfaceMarking& mk : unfinished.markings)
        if (mk.meaning == MarkingMeaning::RowStartPath) mk.points.resize(1);
    const std::vector<MarkingReview> reviews = {quick(pieces, unfinished, built), buildSetup(pieces, unfinished, ReviewDepth::Built).review};
    for (const MarkingReview& r : reviews) {
        bool incomplete = false, addPath = false, problem = false;
        for (const MarkingFinding& f : r.findings) {
            incomplete = incomplete || f.level == FindingLevel::Incomplete;
            addPath = addPath || f.message.find("Add a row-start path") != std::string::npos;
            problem = problem || f.level == FindingLevel::Problem;
        }
        CHECK(incomplete && !addPath && !problem && !r.valid);
    }
    // a path that runs up, round to the side, back down and across its own first stretch
    const MarkingSet loop = withPath(set, m, {onTube(10.5, 0.05, 5), onTube(10.5, 0.05, 25), onTube(10.5, 0.62, 25), onTube(10.5, 0.62, 15.3), onTube(10.5, -0.27, 15.3), onTube(10.5, -0.27, 35)});
    bool crossesItself = false;
    for (const MarkingFinding& f : quick(pieces, loop, built).findings)
        crossesItself = crossesItself || (f.message.find("crosses itself") != std::string::npos && f.at.size() == 1 && f.firstTracePoint >= 0);
    CHECK(crossesItself);
}

// A field made for another start is not used: the review says the row checks are pending, marks the rows Outdated (No
// without any field) and gives no row finding, in particular never "crosses every row once". (Defect: the new path judged
// against the old start's rows.)
TEST(ssm_quick_review_stale_field) {
    const TriMesh m = openTube(10.5, 10.5, 40, 64);
    const SurfacePieces pieces = onePiece(m);
    const ProposalSet ps = propose(pieces);
    const MarkingSet set = setupSet(ps.setups[0]);
    const MarkedBuildProduct built = buildSetup(pieces, set, ReviewDepth::Built);
    const MarkingSet reversed = setupSet(ps.setups[1]);   // the other direction: another start and finish
    const MarkingReview r = quick(pieces, reversed, built);
    bool pending = false;
    for (const MarkingFinding& f : r.findings) {
        pending = pending || f.message.find("Row checks are pending") != std::string::npos;
        CHECK(f.message.find("crosses every row once") == std::string::npos && f.row < 0);
    }
    CHECK(pending && r.measures.rowsKnown == Known::Outdated && r.measures.contourEventsKnown == Known::Outdated);
    const MarkedBuildProduct none;
    const MarkingReview n = quick(pieces, set, none);
    CHECK(n.measures.rowsKnown == Known::No);
}

// A quick review knows only what it computed: rows and contour events from a current field, nothing about cells, coverage,
// distortion or the layout's checks (Known::No). A built layout knows those; one that is built but fails a check (a narrow
// end) keeps valid true and checksPass false. (Defect: an uncomputed measure shown as 0 or PASS, or "valid" taken for
// "passes its checks".)
TEST(ssm_quick_review_measures_known) {
    const TriMesh m = openTube(10.5, 10.5, 40, 64);
    const SurfacePieces pieces = onePiece(m);
    const ProposalSet ps = propose(pieces);
    const MarkingSet set = setupSet(ps.setups[0]);
    const MarkedBuildProduct built = buildSetup(pieces, set, ReviewDepth::Built);
    const LayoutMeasures q = quick(pieces, set, built).measures;
    CHECK(q.rowsKnown == Known::Current && q.contourEventsKnown == Known::Current);
    CHECK(q.unsupportedKnown == Known::No && q.coverageKnown == Known::No && q.distortionKnown == Known::No && q.checksKnown == Known::No && q.exceptionalKnown == Known::No);
    CHECK(built.review.measures.checksKnown == Known::Current && built.review.measures.checksPass);
    // a steep taper to a narrow opening: built, but not passing its checks
    const TriMesh narrow = openTube(14, 2.2, 40, 64);
    const SurfacePieces np = onePiece(narrow);
    const ProposalSet nps = propose(np);
    const MarkedBuildProduct nb = buildSetup(np, setupSet(nps.setups[0]), ReviewDepth::Built);
    std::printf("  narrow taper: valid %d, checks pass %d, unsupported %d, coverage %d\n", nb.review.valid ? 1 : 0, nb.review.measures.checksPass ? 1 : 0,
                nb.review.measures.unsupportedCells, nb.review.measures.coveragePass ? 1 : 0);
    CHECK(nb.review.valid && !nb.review.measures.checksPass && nb.review.measures.checksKnown == Known::Current);
}
