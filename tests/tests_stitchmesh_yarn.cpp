// knit_tests_stitchmesh, yarn embedding (src/stitchmesh/yarn_embed.h): the regular patch's yarn is one strand with two yarn
// ends, continuous in position and direction at every port, every head interlocked in the required crossing order; the
// boundary and pole constructions likewise; and the checks fail with the defect back (a flipped head, a missing face).
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "harness.h"
#include "knitting/layout/base_surface.h"
#include "knitting/numerics.h"
#include "knitting/tube/yarn_surface.h"
#include "model/mesh_query.h"
#include "stitchmesh/construction_fixtures.h"
#include "stitchmesh/coverage.h"
#include "stitchmesh/relax_mesh.h"
#include "stitchmesh/guided_records.h"
#include "stitchmesh/virtual_shift.h"
#include "stitchmesh/yarn_curves.h"
#include "stitchmesh/yarn_embed.h"
#include "stitchmesh/yarn_settle.h"

using namespace knit;

namespace {

constexpr double kTestPi = 3.14159265358979323846;

// A ring tube on a surface of revolution whose ring circumferences equal counts x w, relaxed with the CP3 settings.
bool relaxedTube(const std::vector<int>& counts, const std::vector<std::vector<int>>& shaping, bool boundaryShaping, TypedStitchMesh& m,
                 std::vector<std::vector<int>>& faces, TriMesh& ref) {
    const double w = 5.5, h = 3.4;
    TriMesh unusedRef;
    std::vector<DeclaredBoundary> decl;
    std::string err;
    if (!ringTubeFixture(counts, shaping, 20, h, m, unusedRef, decl, &err, boundaryShaping)) return false;
    std::vector<double> rho, z;
    for (size_t k = 0; k < counts.size(); ++k) {
        rho.push_back(counts[k] * w / (2 * kTestPi));
        z.push_back(k == 0 ? 0.0 : z.back() + std::sqrt(std::max(0.0, h * h - (rho[k] - rho[k - 1]) * (rho[k] - rho[k - 1]))));
    }
    ref = revolutionThroughRings(rho, z, 192, 4);
    const MeshQuery q(ref);
    for (size_t k = 0; k < m.rings.size(); ++k)
        for (int v : m.rings[k]) {
            const double phi = std::atan2(m.mesh.pos[v].y, m.mesh.pos[v].x);
            m.mesh.pos[v] = q.closest(Vec3{rho[k] * std::cos(phi), rho[k] * std::sin(phi), z[k]}).p;
        }
    VirtualShift vs;
    if (!buildVirtualShift(m, ShiftReading::PathSideDown, vs, &err)) return false;
    faces = vs.corners;
    std::vector<RelaxFace> rf;
    for (int f = 0; f < m.mesh.faces(); ++f) rf.push_back(RelaxFace{vs.corners[f], m.role[f]});
    std::vector<char> border(m.mesh.pos.size(), 0);
    for (const std::vector<int>& l : polygonBorderLoops(m.mesh))
        for (int v : l) border[v] = 1;
    MeshRelaxRequest r;
    r.w = w;
    r.h = h;
    RelaxProduct p;
    if (!relaxStitchMesh(ref, triangleBoundaryLoops(ref), rf, m.mesh.pos, border, r, p, &err)) return false;
    m.mesh.pos = p.pos;
    return true;
}

YarnCheck yarnOf(const TypedStitchMesh& m, const std::vector<std::vector<int>>& faces, const TriMesh& ref, const YarnSettings& st, YarnEmbedding& e) {
    std::string err;
    CHECK(embedYarn(m, faces, ref, 5.5, 3.4, st, e, &err));
    return checkYarn(m, faces, ref, e, st);
}

// Numerics that move nothing: relax and store leave the curves as they are, so a settle runs its whole schedule on the
// CPU checks alone (the subject is the settle's own bookkeeping, not the solver). It counts the relax calls.
class StillNumerics : public KnitNumerics {
public:
    int relaxCalls = 0;
    const char* name() const override { return "still (test)"; }
    ValidationReport validate(const CurveSet&, const PatternParams&, const ValidateOptions&) override { return ValidationReport{}; }
    CurvatureResult curvature(const YarnPath&, double, const Allowances&, const ExecPolicy&) override { return CurvatureResult{}; }
    void relax(const RelaxRequest&) override { ++relaxCalls; }
    void store(CurveSet&, const Controls&, const ExecPolicy&) override {}
    std::unique_ptr<MeshSurface> meshSurface(const TriMesh&) override { return nullptr; }
    std::unique_ptr<BaseSurface> pieceQuery(const Chart&) override { return nullptr; }
    void classifyLattice(const UvAtlas&, const LatticeChartQuery&, const LatticeSolids&, std::vector<LatticeNodeClass>&) override {}
    SurfaceCheck checkRenderedSurface(const CurveSet&, const PatternParams&, const SurfaceParams&) override { return SurfaceCheck{}; }
    BoundaryCheck checkBoundary(const CurveSet&, const SurfaceBinding&, const MeshSurface&, double) override { return BoundaryCheck{}; }
};

} // namespace

namespace {

// The settle of a 12 x 4 tube under each choice for a failed construction, counting the solver's relax calls.
struct SettleCase {
    TypedStitchMesh m;
    std::vector<std::vector<int>> faces;
    TriMesh ref;
    int face = -1;  // strip 1, slot 4
};

SettleCase settleTube() {
    SettleCase c;
    CHECK(relaxedTube({12, 12, 12, 12}, {{}, {}, {}}, false, c.m, c.faces, c.ref));
    for (int f = 0; f < c.m.mesh.faces(); ++f)
        if (c.m.strip[f] == 1 && c.m.slot[f] == 4) c.face = f;
    CHECK(c.face >= 0);
    return c;
}

SettleSettings quickSettle(FailedConstruction failed) {
    SettleSettings ss;
    ss.band = false;
    ss.maxRepeats = 1;
    ss.failedConstruction = failed;
    return ss;
}

SettleResult settleWith(const SettleCase& c, const YarnEmbedding& e, const YarnSettings& st, const SettleSettings& ss, int& relaxCalls) {
    StillNumerics still;
    std::ostringstream report, log;
    SettleResult r = settleYarn({e, c.m, c.faces, c.ref, st, "tube"}, ss, still, SettleHooks{}, report, log);
    relaxCalls = still.relaxCalls;
    return r;
}

} // namespace

// A raw construction that fails (one face's head flipped) is never handed to the solver unless asked for: the lab refuses
// it; the app's Final (ReturnUnsettled) returns the converted yarn with zero relax calls, its target gate not run (never a
// pass) and the failure kept; Draft (convert only) likewise; only the lab's explicit diagnostic opt-in (SettleAnyway)
// settles it, marked, so the construction cannot pass. (Defect at 4745fdc: the app settled every failed construction,
// minutes of solve on a crossing order the solver cannot repair.)
TEST(ssm_settle_failed_raw_construction_not_solved) {
    const SettleCase c = settleTube();
    YarnSettings flip;
    flip.flipHeadsOfFace = c.face;
    YarnEmbedding wrong;
    CHECK(yarnOf(c.m, c.faces, c.ref, flip, wrong).status == Status::Fail);
    int calls = -1;
    const SettleResult lab = settleWith(c, wrong, flip, quickSettle(FailedConstruction::Refuse), calls);
    CHECK(lab.end == SettleEnd::RawConstructionFailed && lab.rawConstructionFailed && calls == 0);
    const SettleResult app = settleWith(c, wrong, flip, quickSettle(FailedConstruction::ReturnUnsettled), calls);
    std::printf("  app Final: %s\n", app.verdict.c_str());
    CHECK(app.end == SettleEnd::ConstructionFailedNotSettled && calls == 0 && app.stages.empty());
    CHECK(app.rawConstructionFailed && app.constructionFailed() && !app.constructionKept);
    CHECK(!app.curves.yarns.empty() && !app.gateFinal.valid());   // the converted yarn comes back; the gate did not run
    CHECK(app.verdict.find("raw construction check FAILS") != std::string::npos && app.verdict.find("skipped") != std::string::npos);
    SettleSettings draftSettings = quickSettle(FailedConstruction::ReturnUnsettled);
    draftSettings.convertOnly = true;
    const SettleResult draft = settleWith(c, wrong, flip, draftSettings, calls);
    CHECK(draft.end == SettleEnd::ConvertedOnly && draft.rawConstructionFailed && calls == 0 && !draft.curves.yarns.empty());
    const SettleResult diagnostic = settleWith(c, wrong, flip, quickSettle(FailedConstruction::SettleAnyway), calls);
    CHECK(diagnostic.end == SettleEnd::Solved && diagnostic.settledFailedConstruction && calls > 0 && !diagnostic.constructionKept);
    CHECK(diagnostic.verdict.rfind("diagnostic settle of a failed construction", 0) == 0);
}

// A clean construction enters the solver under the app's choice (its Final is unchanged): relax is called and nothing is
// marked. (Defect: skipping the settle for every yarn.)
TEST(ssm_settle_clean_construction_enters_the_solver) {
    const SettleCase c = settleTube();
    YarnSettings plain;
    YarnEmbedding right;
    CHECK(yarnOf(c.m, c.faces, c.ref, plain, right).status == Status::Pass);
    int calls = 0;
    const SettleResult r = settleWith(c, right, plain, quickSettle(FailedConstruction::ReturnUnsettled), calls);
    CHECK(r.end == SettleEnd::Solved && calls > 0 && !r.constructionFailed() && !r.settledFailedConstruction && r.constructionKept);
}

// The stop rule after a call at the target diameter (afterTargetCall). A valid gate ends the run: as passed when the
// construction is the generated one, as lost when a crossing changed on the way, never as passed on the gate alone. An
// invalid gate goes on while the open regions fall or the shortfall shrinks, and ends as stalled after stallCalls calls
// that bring neither. (Defects: a fixed number of calls with the last one judged, fixtures being valid at one call and
// not at the last; and a stop on the gate alone.)
TEST(ssm_settle_stop_rule_after_a_target_call) {
    const double D = 1.16, c = 0.1;
    SettleGate valid;
    valid.regular = valid.fold = valid.clearance = Status::Pass;
    valid.rhoLo = 0.64;
    valid.minDistance = 1.27;
    TargetProgress progress;
    CHECK(afterTargetCall(valid, D, c, true, 4, progress) == AfterTargetCall::GatePassed);
    CHECK(afterTargetCall(valid, D, c, false, 4, progress) == AfterTargetCall::ConstructionLost);
    SettleGate open = valid;
    open.clearance = Status::Unresolved;
    open.failRegions = 1;
    open.unresolvedRegions = 2;
    open.minDistance = 1.24;
    CHECK(!open.valid());
    progress = TargetProgress{};
    CHECK(afterTargetCall(open, D, c, true, 4, progress) == AfterTargetCall::Continue);   // the first call sets the mark
    open.unresolvedRegions = 1;   // fewer open regions: a gain
    CHECK(afterTargetCall(open, D, c, true, 4, progress) == AfterTargetCall::Continue && progress.callsWithoutGain == 0);
    open.minDistance = 1.25;      // as many regions, the shortfall 0.01 mm smaller: a gain
    CHECK(afterTargetCall(open, D, c, true, 4, progress) == AfterTargetCall::Continue && progress.callsWithoutGain == 0);
    for (int k = 1; k <= 3; ++k) CHECK(afterTargetCall(open, D, c, true, 4, progress) == AfterTargetCall::Continue && progress.callsWithoutGain == k);
    CHECK(afterTargetCall(open, D, c, true, 4, progress) == AfterTargetCall::Stalled);
}

// A yarn at the real diameter whose solver makes no progress: after the first call at the target, stallCalls more without
// fewer open regions or a smaller shortfall end the run as stalled, well before the 24 calls allowed.
TEST(ssm_settle_without_progress_ends_stalled) {
    const SettleCase c = settleTube();
    const YarnSettings st;
    YarnEmbedding e;
    CHECK(yarnOf(c.m, c.faces, c.ref, st, e).status == Status::Pass);
    SettleSettings ss = quickSettle(FailedConstruction::ReturnUnsettled);
    ss.maxRepeats = 60;
    int calls = 0;
    const SettleResult r = settleWith(c, e, st, ss, calls);
    int atTarget = 0;
    for (const SettleStage& stage : r.stages) atTarget += stage.D == st.diameter;
    std::printf("  %s (%d calls, %d at the target, stop: %s)\n", r.verdict.c_str(), calls, atTarget, toString(r.stop));
    CHECK(r.stop == SettleStop::Stalled && !r.gateFinal.valid() && atTarget == ss.stallCalls + 1);
    CHECK(r.verdict.find("stalled") != std::string::npos && !r.chainCertified);
}

// A yarn whose strand breaks (one face without yarn: unmatched ports) is refused before the solver under every choice,
// the diagnostic opt-in included. (Defect: a broken strand settled, then lost at the stitch records.)
TEST(ssm_settle_broken_strand_refused_always) {
    const SettleCase c = settleTube();
    YarnSettings skip;
    skip.skipFace = c.face;
    YarnEmbedding broken;
    const YarnCheck check = yarnOf(c.m, c.faces, c.ref, skip, broken);
    CHECK(check.status == Status::Fail && check.portsUnmatched > 0);
    int calls = -1;
    CHECK(settleWith(c, broken, skip, quickSettle(FailedConstruction::Refuse), calls).end == SettleEnd::RawConstructionFailed && calls == 0);
    for (FailedConstruction f : {FailedConstruction::ReturnUnsettled, FailedConstruction::SettleAnyway}) {
        const SettleResult r = settleWith(c, broken, skip, quickSettle(f), calls);
        std::printf("  %s: %s\n", toString(f), r.verdict.c_str());
        CHECK(r.end == SettleEnd::YarnDisconnected && calls == 0 && r.curves.yarns.empty());
    }
}

// A construction that only the conversion breaks (interpolation points 3 mm apart cut the heads' corners) is reported as the
// converted check's failure, apart from the raw one: refused by the lab, returned unsettled with zero relax calls by the
// app, the added defects named. (Defect: the converted check skipped or folded into the raw one.)
TEST(ssm_settle_failed_converted_construction_not_solved) {
    const SettleCase c = settleTube();
    YarnSettings plain;
    YarnEmbedding right;
    CHECK(yarnOf(c.m, c.faces, c.ref, plain, right).status == Status::Pass);
    const double spacing = 3.0;   // interpolation points 3 mm apart (the default is 0.2 mm; 2 mm still keeps the construction)
    int calls = -1;
    SettleSettings lab = quickSettle(FailedConstruction::Refuse);
    lab.spacing = spacing;
    const SettleResult refused = settleWith(c, right, plain, lab, calls);
    CHECK(refused.end == SettleEnd::ConvertedConstructionFailed && !refused.rawConstructionFailed && refused.convertedConstructionFailed && calls == 0);
    SettleSettings app = lab;
    app.failedConstruction = FailedConstruction::ReturnUnsettled;
    const SettleResult r = settleWith(c, right, plain, app, calls);
    std::printf("  app Final: %s\n", r.verdict.c_str());
    CHECK(r.end == SettleEnd::ConstructionFailedNotSettled && calls == 0 && !r.rawConstructionFailed && r.convertedConstructionFailed && r.conversionAddedDefects);
    CHECK(r.verdict.find("converted construction check FAILS") != std::string::npos);
}

TEST(ssm_yarn_regular_patch) {
    TypedStitchMesh m;
    std::vector<std::vector<int>> faces;
    TriMesh ref;
    CHECK(relaxedTube({12, 12, 12, 12}, {{}, {}, {}}, false, m, faces, ref));
    YarnEmbedding e;
    const YarnCheck c = yarnOf(m, faces, ref, YarnSettings{}, e);
    std::printf("  %s", c.text().c_str());
    CHECK(c.status == Status::Pass);
    CHECK(c.strands == 1 && c.openStrands == 1 && c.yarnEnds == 2);
    CHECK(c.portsUnmatched == 0 && c.portsOverused == 0);
    CHECK(c.maxPortGap < 1e-9 && c.maxTangentJump < 1e-3);
    CHECK(c.interlocksExpected == 24 && c.interlocksVerified == 24);   // 2 interior rings x 12 heads
    CHECK(c.facesCrossingsOk == c.facesChecked && c.facesChecked == 36);
    CHECK(e.caps == 12);   // the bind-off ring's loops, all threaded by the binding thread
}

// Increases next to each other (16 -> 20 stitches in one row): the row-start face of that strip spans two rows and folds
// across the flare. With the knots placed on the face the yarn is clean; moved to the closest point of the reference (the
// first-yarn rule) the head of that face leaves it and loses a crossing.
TEST(ssm_yarn_adjacent_increases_need_the_face_placement) {
    TypedStitchMesh m;
    std::vector<std::vector<int>> faces;
    TriMesh ref;
    CHECK(relaxedTube({16, 16, 20, 20}, {{}, {2, 3, 9, 10}, {}}, false, m, faces, ref));
    YarnSettings st;
    CHECK(st.placement == KnotPlacement::FaceInterpolation);
    YarnEmbedding e;
    const YarnCheck onFace = yarnOf(m, faces, ref, st, e);
    std::printf("  on the face: %s", onFace.text().c_str());
    CHECK(onFace.status == Status::Pass && onFace.facesCrossingsOk == onFace.facesChecked);
    st.placement = KnotPlacement::ReferenceProjection;
    const YarnCheck projected = yarnOf(m, faces, ref, st, e);
    CHECK(projected.status == Status::Fail && projected.facesCrossingsOk == projected.facesChecked - 1);
}

TEST(ssm_yarn_shaping_boundary_pole) {
    struct Case {
        std::vector<int> counts;
        std::vector<std::vector<int>> shaping;
        bool boundary;
    };
    for (const Case& k : {Case{{12, 12, 13, 13}, {{}, {6}, {}}, false}, Case{{13, 13, 12, 12}, {{}, {6}, {}}, false}, Case{{12, 13, 13, 13}, {{6}, {}, {}}, true},
                          Case{{13, 13, 13, 12}, {{}, {}, {6}}, true}}) {
        TypedStitchMesh m;
        std::vector<std::vector<int>> faces;
        TriMesh ref;
        CHECK(relaxedTube(k.counts, k.shaping, k.boundary, m, faces, ref));
        YarnEmbedding e;
        const YarnCheck c = yarnOf(m, faces, ref, YarnSettings{}, e);
        std::printf("  %d -> %d: %s, interlocks %d / %d\n", k.counts.front(), k.counts.back(), toString(c.status), c.interlocksVerified, c.interlocksExpected);
        CHECK(c.status == Status::Pass && c.interlocksVerified == c.interlocksExpected);
    }
    // the gathered pole start: one strand from the round's tail through the join into the helix
    TypedStitchMesh m;
    TriMesh unusedRef;
    std::vector<DeclaredBoundary> decl;
    std::string err;
    CHECK(poleDiscFixture({5, 7, 7}, {{1, 3}, {}}, 3.4, m, unusedRef, decl, &err));
    TriMesh ref = revolutionThroughRings({1e-3, 3.4, 6.8, 10.2}, {0, 0, 0, 0}, 192, 4);
    bool down = false;
    for (const Vec3& n : ref.vnormal) down = down || n.z < 0;
    if (down) {
        for (auto& t : ref.tri) std::swap(t[1], t[2]);
        ref.computeNormals();
    }
    VirtualShift vs;
    CHECK(buildVirtualShift(m, ShiftReading::PathSideDown, vs, &err));
    YarnEmbedding e;
    YarnCheck c = yarnOf(m, vs.corners, ref, YarnSettings{}, e);
    std::printf("  pole start: %s", c.text().c_str());
    CHECK(e.poleJoined && c.strands == 1 && c.yarnEnds == 2 && c.facesCrossingsOk == c.facesChecked);
    // without the join the round stays a closed strand of its own: two strands, the construction's open question
    YarnSettings noJoin;
    noJoin.joinPoleRound = false;
    YarnEmbedding e2;
    c = yarnOf(m, vs.corners, ref, noJoin, e2);
    CHECK(!e2.poleJoined && c.closedStrands == 1 && c.status == Status::Fail);
}

TEST(ssm_yarn_checks_fail_with_the_defect) {
    TypedStitchMesh m;
    std::vector<std::vector<int>> faces;
    TriMesh ref;
    CHECK(relaxedTube({12, 12, 12, 12}, {{}, {}, {}}, false, m, faces, ref));
    int face = -1;
    for (int f = 0; f < m.mesh.faces(); ++f)
        if (m.strip[f] == 1 && m.slot[f] == 4) face = f;
    CHECK(face >= 0);
    // a head's depth flipped: that face's crossings are wrong, one interlock fewer, nowhere else
    YarnSettings flip;
    flip.flipHeadsOfFace = face;
    YarnEmbedding e1;
    const YarnCheck c1 = yarnOf(m, faces, ref, flip, e1);
    CHECK(c1.status == Status::Fail && c1.interlocksVerified == c1.interlocksExpected - 1 && c1.facesCrossingsOk == c1.facesChecked - 1);
    bool named = false;
    for (const std::string& d : c1.defects) named = named || d.find("face " + std::to_string(face) + " ") != std::string::npos;
    CHECK(named);
    // a face without yarn: its ports are left with one side only and the strand breaks
    YarnSettings skip;
    skip.skipFace = face;
    YarnEmbedding e2;
    const YarnCheck c2 = yarnOf(m, faces, ref, skip, e2);
    CHECK(c2.status == Status::Fail && c2.portsUnmatched > 0 && c2.strands > 1);
}

// The yarn-settle adapter (src/stitchmesh/yarn_curves.h): the strand as one uniform cubic B-spline through equal-arc points
// stays within 0.05 D of the Hermite yarn, keeps its length, and re-embedded on the unmoved curve keeps the construction (one
// strand, every port shared, every interlock); a closed strand is refused.
TEST(ssm_yarn_curves_conversion) {
    TypedStitchMesh m;
    std::vector<std::vector<int>> faces;
    TriMesh ref;
    CHECK(relaxedTube({12, 12, 12, 12}, {{}, {}, {}}, false, m, faces, ref));
    const YarnSettings st;
    YarnEmbedding e;
    const YarnCheck c0 = yarnOf(m, faces, ref, st, e);
    CHECK(c0.status == Status::Pass);
    YarnConversion cv;
    std::string err;
    CHECK(yarnToCurves(e, st, 0.2, 32, cv, &err));
    std::printf("  %s", cv.text().c_str());
    CHECK(cv.curves.yarns.size() == 1 && cv.curves.kind == CurveKind::Surface);
    // the B-spline against the Hermite yarn: within 0.06 D. Measured on the settle fixtures at spacing 0.2 mm: 0.033 to
    // 0.062 mm with the face placement (0.033 to 0.054 mm with the reference projection); the converted construction is
    // checked on its own below
    CHECK(cv.maxDeviation < 0.06 * st.diameter);
    CHECK(std::fabs(cv.curveLength / cv.hermiteLength - 1) < 0.005);
    for (size_t i = 0; i < e.pieces.size(); ++i) CHECK(cv.strandOfPiece[i] == 0 && cv.u0[i] < cv.u1[i]);
    const YarnEmbedding back = embeddingFromCurves(e, cv, cv.curves, ref, st.diameter);
    const YarnCheck c1 = checkYarn(m, faces, ref, back, st);
    std::printf("  converted: %s", c1.text().c_str());
    CHECK(c1.status == Status::Pass);
    CHECK(c1.interlocksVerified == c0.interlocksVerified && c1.facesOwnPlaneOk == c1.facesChecked);
    CHECK(c1.maxPortGap < 1e-9);   // neighbouring pieces evaluate the same curve at the same parameter
    YarnEmbedding closed = e;
    closed.strands[0].closed = true;
    YarnConversion refused;
    CHECK(!yarnToCurves(closed, st, 0.2, 32, refused, &err));
    CHECK(err.find("closed") != std::string::npos);
}

// The crossing rule in the yarn's own best-fit plane (YarnFaceCheck::ownPlaneOk): on the placed yarn it agrees with the
// reference frame on every face; a flipped head fails in both frames, at that face only.
TEST(ssm_yarn_own_plane_frame) {
    TypedStitchMesh m;
    std::vector<std::vector<int>> faces;
    TriMesh ref;
    CHECK(relaxedTube({12, 12, 13, 13}, {{}, {6}, {}}, false, m, faces, ref));
    YarnEmbedding e;
    const YarnCheck c = yarnOf(m, faces, ref, YarnSettings{}, e);
    CHECK(c.status == Status::Pass && c.facesOwnPlaneOk == c.facesChecked && c.interlocksOwnPlane == c.interlocksVerified);
    int face = -1;
    for (int f = 0; f < m.mesh.faces(); ++f)
        if (m.strip[f] == 1 && m.slot[f] == 3) face = f;
    CHECK(face >= 0);
    YarnSettings flip;
    flip.flipHeadsOfFace = face;
    YarnEmbedding e1;
    const YarnCheck c1 = yarnOf(m, faces, ref, flip, e1);
    CHECK(c1.facesOwnPlaneOk == c1.facesChecked - 1 && c1.facesCrossingsOk == c1.facesChecked - 1);
    for (const YarnFaceCheck& fc : c1.faces) CHECK((fc.face == face) == (!fc.crossingsOk && !fc.ownPlaneOk));
}

// The calibrated template set (the template calibration of 2026-11-19): on the regular patch, the isolated increase and
// the isolated decrease the construction holds as with the baseline (one strand, every port shared with no gap or direction
// jump, every face's crossings as required); its raw yarn is closer than D at fewer samples than the baseline's on the same
// mesh; a flipped head still fails at that face only.
TEST(ssm_yarn_calibrated_templates) {
    YarnSettings base, cal;
    cal.templateSet = "calibrated";
    struct Case {
        std::vector<int> counts;
        std::vector<std::vector<int>> shaping;
    };
    for (const Case& k : {Case{{12, 12, 12, 12}, {{}, {}, {}}}, Case{{12, 12, 13, 13}, {{}, {6}, {}}}, Case{{13, 13, 12, 12}, {{}, {6}, {}}}}) {
        TypedStitchMesh m;
        std::vector<std::vector<int>> faces;
        TriMesh ref;
        CHECK(relaxedTube(k.counts, k.shaping, false, m, faces, ref));
        YarnEmbedding eb, ec;
        const YarnCheck cb = yarnOf(m, faces, ref, base, eb);
        const YarnCheck cc = yarnOf(m, faces, ref, cal, ec);
        std::printf("  calibrated: %s", cc.text().c_str());
        CHECK(cc.status == Status::Pass);
        CHECK(cc.strands == 1 && cc.yarnEnds == 2 && cc.portsUnmatched == 0 && cc.portsOverused == 0);
        CHECK(cc.maxPortGap < 1e-9 && cc.maxTangentJump < 1e-3);
        CHECK(cc.facesCrossingsOk == cc.facesChecked && cc.interlocksVerified == cc.interlocksExpected && cc.interlocksExpected == cb.interlocksExpected);
        CHECK(cc.contactSamples * cb.samples < cb.contactSamples * cc.samples);   // a smaller fraction of samples closer than D
    }
    TypedStitchMesh m;
    std::vector<std::vector<int>> faces;
    TriMesh ref;
    CHECK(relaxedTube({12, 12, 12, 12}, {{}, {}, {}}, false, m, faces, ref));
    int face = -1;
    for (int f = 0; f < m.mesh.faces(); ++f)
        if (m.strip[f] == 1 && m.slot[f] == 4) face = f;
    YarnSettings flip = cal;
    flip.flipHeadsOfFace = face;
    YarnEmbedding ef;
    const YarnCheck cf = yarnOf(m, faces, ref, flip, ef);
    CHECK(cf.status == Status::Fail && cf.facesCrossingsOk == cf.facesChecked - 1);
    for (const YarnFaceCheck& fc : cf.faces) CHECK((fc.face == face) == !fc.crossingsOk);
}

// Stitch records of the guided route (src/stitchmesh/guided_records.h): one record per loop (K and d12k one, kp two), parents
// as many as the yarn check's interlocks, d12k two parents, the kp's two loops one shared parent, ownership gap-free and
// found again by parameter; a head moved off every top edge is refused with the reason.
TEST(guided_records_from_stitch_mesh) {
    struct Case {
        std::vector<int> counts;
        std::vector<std::vector<int>> shaping;
        int increases, decreases;
    };
    for (const Case& k : {Case{{12, 12, 12, 12}, {{}, {}, {}}, 0, 0}, Case{{12, 12, 13, 13}, {{}, {6}, {}}, 2, 0}, Case{{13, 13, 12, 12}, {{}, {6}, {}}, 0, 1}}) {
        TypedStitchMesh m;
        std::vector<std::vector<int>> faces;
        TriMesh ref;
        CHECK(relaxedTube(k.counts, k.shaping, false, m, faces, ref));
        const YarnSettings st;
        YarnEmbedding e;
        const YarnCheck c = yarnOf(m, faces, ref, st, e);
        CHECK(c.status == Status::Pass);
        YarnConversion cv;
        std::string err;
        CHECK(yarnToCurves(e, st, 0.2, 32, cv, &err));
        CurveSet curves = cv.curves;
        GuidedRecordCounts n;
        const bool ok = guidedStitchRecords(m, faces, e, cv, ref, c.interlocksExpected, curves, n, &err);
        std::printf("  %d -> %d: %d loops, %d INC, %d DEC, %d first turn, %d parent links (%s)\n", k.counts.front(), k.counts.back(), n.loops, n.increases, n.decreases,
                    n.firstTurn, n.parentLinks, err.c_str());
        CHECK(ok);
        CHECK(n.loops == m.mesh.faces() + k.increases / 2 && n.increases == k.increases && n.decreases == k.decreases);
        CHECK(n.parentLinks == c.interlocksExpected && static_cast<int>(curves.topo.stitches.size()) == n.loops);
        CHECK(curves.yarns.size() == 1 && static_cast<int>(curves.yarns[0].stitches.size()) == n.loops);
        int sharedIncParent = 0;
        for (const StitchRecord& r : curves.topo.stitches) {
            CHECK(curves.topo.stitchAt(r.yarn, 0.5 * (r.uBegin + r.uEnd)) == r.id);
            if (r.op == StitchOp::Dec) CHECK(r.nParents == 2);
            if (r.op == StitchOp::K && r.turn > 0) CHECK(r.nParents == 1);
            if (r.op == StitchOp::Inc && r.incSlot == 1) sharedIncParent += curves.topo.stitches[r.id - 1].parents[0] == r.parents[0];
            CHECK(std::fabs(norm(r.U) - 1) < 1e-9 && std::fabs(dot(r.U, r.N)) < 1e-9);
        }
        CHECK(sharedIncParent == k.increases / 2);
        const std::vector<StitchSpan>& spans = curves.yarns[0].stitches;
        for (size_t i = 1; i < spans.size(); ++i) CHECK(spans[i - 1].uEnd <= spans[i].uBegin + 1e-9);
        // the defect back: one head no longer starts at any loop's top edge
        YarnEmbedding broken = e;
        for (YarnPiece& p : broken.pieces)
            if (p.kind == YarnPieceKind::Head) {
                p.start.a = -7;
                break;
            }
        CurveSet again = cv.curves;
        CHECK(!guidedStitchRecords(m, faces, broken, cv, ref, c.interlocksExpected, again, n, &err) && err.find("passes no loop's top edge") != std::string::npos);
        // and a count the yarn check does not share
        again = cv.curves;
        CHECK(!guidedStitchRecords(m, faces, e, cv, ref, c.interlocksExpected + 1, again, n, &err) && err.find("interlocks") != std::string::npos);
    }
}

// The settle defaults are the frozen closeout settings (tests/data/guided/settle_settings.json, the stitch-mesh
// experiment's closeout): the app's guided route and the lab run the same settle.
TEST(guided_settle_defaults_frozen) {
    std::ifstream f(std::string(KNIT_SOURCE_DIR) + "/tests/data/guided/settle_settings.json");
    CHECK(f.good());
    std::stringstream text;
    text << f.rdbuf();
    const std::string json = text.str();
    auto number = [&json](const char* key) {
        const std::string k = std::string("\"") + key + "\": ";
        const size_t at = json.find(k);
        return at == std::string::npos ? -1.0 : std::atof(json.c_str() + at + k.size());
    };
    const SettleSettings s;
    const YarnSettings y;
    CHECK(number("spacing_mm") == s.spacing && number("dense_points_per_hermite_span") == s.densePerSpan);
    CHECK(number("target_diameter_mm") == y.diameter && number("clearance_mm") == s.clearance);
    CHECK(number("soft_contact_margin_mm") == s.contactMargin && number("soft_bend_factor") == s.bendFactor && number("length_stiffness") == s.lengthStiffness);
    CHECK(number("band_half_mm") == s.bandHalf && number("samples_per_span") == s.samplesPerSpan && number("iterations") == s.iterations);
    CHECK(number("refresh") == s.refresh && number("max_repeats_at_target") == s.maxRepeats && number("stall_calls") == s.stallCalls);
    // the rules written as text: the start-diameter floor, the growth, the cap fraction, the pinned controls, the budget
    for (const char* phrase : {"from max(D0, 0.05 mm)", "log 1.6", "displacement cap 0.5 (D_{k-1} + c)", "first 3 and last 3 controls pinned", "after 50 min of solve time"})
        CHECK(json.find(phrase) != std::string::npos);
    CHECK(s.D0Min == 0.05 && s.growth == 1.6 && s.displacementFraction == 0.5 && s.pinnedEnd == 3 && s.budgetMinutes == 50 && s.band && y.templateSet == "baseline");
}
