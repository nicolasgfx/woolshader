// The guided part: one marked piece knitted as a stockinette stitch mesh, its yarn generated and settled (knit_run.h).
//
// In order: the stitch mesh from the markings (stitchmesh/marked_build.h, the Mark step's backend); the yarn embedded in its
// virtual (helical) faces (yarn_embed.h, the baseline templates); the settle on the run's device numerics (yarn_settle.h,
// the frozen closeout settings); stitch records derived from the stitch mesh (guided_records.h); then the retained curves
// placed in the part's library frame. Each stage's outcome is kept in a GuidedPart. Draft stops before the solve; a yarn
// whose crossing order fails its check is never settled (the solver cannot repair one), so Final then keeps the generated
// yarn as curves, shown and exported as a diagnostic. The verdict is the worst of three
// existing checks, each kept as it is: the target gate (certified curvature and clearance at D + c), the construction after
// the motion (face crossings in the reference frame and in the yarn's own plane), and the stitch mesh's own checks
// (connectivity, geometry, coverage of the knit region). The local signature is not among them: it describes production's
// template, not these (its line says so instead of a vacuous pass).
#include <chrono>
#include <cstdio>
#include <sstream>

#include "gpu/knitting.h"
#include "knitting/tube/param_text.h"
#include "model/cancel.h"
#include "scene/knit/knit_run.h"
#include "settings/settings_mapping.h"
#include "stitchmesh/guided_records.h"
#include "stitchmesh/marked_build.h"
#include "stitchmesh/prepared_reference.h"
#include "stitchmesh/yarn_settle.h"

namespace knit {

namespace {

using Clock = std::chrono::steady_clock;

double msSince(Clock::time_point t0) { return std::chrono::duration<double, std::milli>(Clock::now() - t0).count(); }

// A direction of the file's axes in the library's (the import's rotation; no offset).
Vec3 libraryDirection(const GuidedInputs& in, const Vec3& v) { return in.yUp ? Vec3{v.x, -v.z, v.y} : v; }

// A mesh in file axes and mm as a mesh in the part's frame (library mm minus the part's centre).
TriMesh meshInPartFrame(const GuidedInputs& in, const TriMesh& fileMm, const Vec3& centre) {
    TriMesh m;
    m.tri = fileMm.tri;
    for (const Vec3& p : fileMm.pos) m.pos.push_back(libraryPointOfFileMm(in, p) - centre);
    m.computeNormals();
    return m;
}

// The curves and their stitch frames moved the same way; arc-length tables made again.
void curvesInPartFrame(const GuidedInputs& in, const Vec3& centre, CurveSet& cs) {
    for (YarnPath& y : cs.yarns) {
        for (Vec3& c : y.curve.ctrl) c = libraryPointOfFileMm(in, c) - centre;
        y.finalize(1e-6);
    }
    std::vector<double> yarnEnd;
    for (const YarnPath& y : cs.yarns) yarnEnd.push_back(y.curve.uMax());
    for (StitchRecord& r : cs.topo.stitches) {
        r.origin = libraryPointOfFileMm(in, r.origin) - centre;
        r.U = libraryDirection(in, r.U);
        r.V = libraryDirection(in, r.V);
        r.N = libraryDirection(in, r.N);
    }
    cs.topo.buildIndex(yarnEnd);
}

Status gateStatus(const SettleGate& g) {
    if (g.valid()) return Status::Pass;
    return g.regular == Status::Fail || g.fold == Status::Fail || g.clearance == Status::Fail ? Status::Fail : Status::Unresolved;
}

// The construction after the motion: a face failing in both frames may hide a pass-through (UNRESOLVED); a broken strand
// or port is a FAIL.
Status constructionStatus(const SettleResult& r) {
    if (r.letterSwapped > 0) return Status::Unresolved;
    return r.constructionKept ? Status::Pass : Status::Fail;
}

double fractionCloserThanD(const YarnCheck& c) { return c.samples ? static_cast<double>(c.contactSamples) / c.samples : 0.0; }

int facesWronglyCrossed(const YarnCheck& c) { return c.facesChecked - c.facesCrossingsOk; }

// One yarn state's construction in a few words: its verdict, the wrongly crossed faces and the interlocks.
std::string constructionText(const YarnCheck& c) {
    char b[200];
    std::snprintf(b, sizeof b, "%s (%d of %d stitches wrongly crossed, %d of %d interlocks)", toString(c.status), facesWronglyCrossed(c), c.facesChecked, c.interlocksVerified,
                  c.interlocksExpected);
    return b;
}

// The wrongly crossed stitches of one yarn state as findings: how many, then the first few by row and stitch, each at its
// place in the file's axes (mm), so the user can find them on the piece. Other defects (ports, continuity) as they are.
void addCrossingFindings(std::vector<std::string>& findings, const char* state, const YarnCheck& c, const TypedStitchMesh& m,
                         const std::vector<std::vector<int>>& faces) {
    const int wrong = facesWronglyCrossed(c);
    char b[300];
    std::snprintf(b, sizeof b, "%s: %d of %d stitches have a wrong crossing order", state, wrong, c.facesChecked);
    findings.push_back(b);
    int listed = 0;
    for (const YarnFaceCheck& fc : c.faces) {
        if (fc.crossingsOk || fc.face < 0 || fc.face >= static_cast<int>(faces.size())) continue;
        if (++listed > 6) break;
        Vec3 at{0, 0, 0};
        for (int v : faces[fc.face]) at = at + m.mesh.pos[v];
        at = at * (1.0 / static_cast<double>(std::max<size_t>(1, faces[fc.face].size())));
        std::snprintf(b, sizeof b, "%s: row %d, stitch %d (%s) at (%.1f, %.1f, %.1f) mm: crossings %s, expected %s", state, m.strip.empty() ? -1 : m.strip[fc.face],
                      m.slot.empty() ? -1 : m.slot[fc.face], fc.templ.c_str(), at.x, at.y, at.z, fc.found.c_str(), fc.expected.c_str());
        findings.push_back(b);
    }
    if (wrong == 0)
        for (size_t i = 0; i < c.defects.size() && i < 4; ++i) findings.push_back(std::string(state) + ": " + c.defects[i]);
}

// The report's lines for the checks that apply, in place of production's local signature. Without a solve the construction
// line is the yarn's before any settle (`unsettled`: its words) and the target gate is listed as not run.
void guidedReportLines(ValidationReport& r, const GuidedPart& g, Status construction, const std::string& unsettled) {
    std::vector<CheckLine> kept;
    for (const CheckLine& l : r.lines)
        if (l.name.find("signature") == std::string::npos) kept.push_back(l);
    char b[300];
    std::snprintf(b, sizeof b, "not applicable to the guided templates (no production template); the construction is checked per face instead");
    kept.push_back({"local signature", Status::Unresolved, b});
    std::snprintf(b, sizeof b, "%d of %d faces as required in both frames, %d in one frame only, %d in neither; %d of %d interlocks", g.facesChecked - g.facesFailBothFrames - g.facesFailOneFrame,
                  g.facesChecked, g.facesFailOneFrame, g.facesFailBothFrames, g.interlocks, g.interlocksExpected);
    if (!unsettled.empty()) {
        kept.push_back({"target gate (curvature, clearance at D + c)", Status::Unresolved, "not run: the yarn was not settled"});
        kept.push_back({"construction (face crossings)", construction, "the yarn as curves, before any settle: " + unsettled});
    } else {
        kept.push_back({"construction (face crossings)", construction, b});
    }
    kept.push_back({"stitch mesh (connectivity, geometry, coverage)", g.meshChecks, ""});
    r.lines = kept;
}

// True when the piece is knitted as more than one knitting component: some marking in use on it names a component other
// than the first.
bool pieceHasSeveralKnitComponents(const GuidedInputs& in) {
    for (const SurfaceMarking& m : in.markings.markings) {
        const bool inUse = m.origin == MarkingOrigin::User || m.proposal == ProposalState::Accepted;
        if (inUse && m.knitComponent != 0 && pieceOfMarking(*in.pieces, m) == in.piece) return true;
    }
    return false;
}

} // namespace

std::string guidedInputsText(const SceneSettings& s, const GuidedInputs& in, double D) {
    // Numbers through fixed buffers, the names and the solver identity appended: a long markings path once filled a
    // 600-character buffer and cut off the identity, so the report did not say which solver options ran.
    char head[96], triangles[64], numbers[200];
    std::snprintf(head, sizeof head, "source %016llx, piece %d knitting component %d (", static_cast<unsigned long long>(in.markings.digest), in.piece,
                  in.knitComponent + 1);
    std::snprintf(triangles, sizeof triangles, ", triangles %016llx), markings ", static_cast<unsigned long long>(in.pieceTriangles));
    const YarnSettings yarn;
    std::snprintf(numbers, sizeof numbers, ", gauge %.4g x %.4g mm, D %.4g mm, c %.4g mm, templates %s, knots on the %s, guided version %d, ", s.wale, s.course, D,
                  s.clearance, yarn.templateSet.c_str(), toString(yarn.placement), kGuidedVersion);
    const SettleSettings settle = guidedSettleSettings(s);
    return head + in.pieceName + triangles + s.guidedMarkingsDigest + " (" + s.guidedMarkings + "), preparation " + in.markings.preparation + numbers +
           gpu::ncpSolverIdentity(settle.contactStorage, settle.curvatureRows);
}

SettleSettings guidedSettleSettings(const SceneSettings& s) {
    SettleSettings settle;
    // A wrong crossing order is never handed to the solver (it cannot repair one): Final returns the generated yarn
    // unsettled, shown and exported as a diagnostic. A broken strand or port is refused before any of it.
    settle.failedConstruction = FailedConstruction::ReturnUnsettled;
    settle.convertOnly = !s.settle;  // Draft (settle off): the generated yarn converted to curves, no solve
    return settle;
}

void knitGuidedPart(const KnitRun& run, ScenePart& part) {
    const SceneSettings& s = run.settings;
    const GuidedInputs& in = *run.guided;
    auto g = std::make_shared<GuidedPart>();
    part.guided = g;
    const double D = knitDiameter(s, s.diameter);  // the contact diameter production knits at (a compressed yarn is drawn larger)
    g->targetD = D;
    g->inputs = guidedInputsText(s, in, D);
    g->preparation = in.markings.preparation;
    auto step = [&run](const std::string& what, int done, int total) {
        if (run.options.onStep) run.options.onStep(what, done, total);
    };
    auto refuse = [&](const std::string& reason) {
        g->outcome = GuidedOutcome::Refused;
        g->reason = reason;
        part.note = "guided stitch mesh refused: " + reason;
        part.knit.reset();
    };

    // 1. The stitch mesh, as the Mark step builds it, with the markings it uses (no offered proposal is previewed here).
    step("stitch mesh", 0, 3);
    Clock::time_point t0 = Clock::now();
    MarkedBuildRequest rq;
    rq.pieces = in.pieces.get();
    rq.piece = in.piece;
    rq.knitComponent = in.knitComponent;
    // a knitting component lies over the ones before it (a collar over a body): one yarn layer, D + c, out per component
    rq.layerOffsetMm = in.knitComponent * (D + s.clearance);
    rq.markings = &in.markings;
    rq.mmPerUnit = in.mmPerUnit;
    rq.w = s.wale;
    rq.h = s.course;
    rq.depth = ReviewDepth::Relaxed;
    if (!parsePreparation(in.markings.preparation, rq.options.preparation)) return refuse("unknown reference preparation \"" + in.markings.preparation + "\"");
    MarkedBuildProduct product;
    std::string err;
    if (!buildFromMarkings(rq, product, &err)) return refuse(err);
    cancelCheckpoint();
    g->msMesh = msSince(t0);
    for (const MarkingFinding& f : product.review.findings)
        if (f.level != FindingLevel::Ok) g->findings.push_back((f.markingId ? "marking " + std::to_string(f.markingId) + ": " : std::string()) + f.message);
    if (!product.review.valid) return refuse(g->findings.empty() ? std::string("no stitch mesh could be built from the markings") : g->findings.front());
    const MarkedEvaluation& ev = product.evaluation;
    g->route = ev.route;
    g->faces = ev.after.mesh.faces();
    g->stitchMeshDigest = typedDigest(ev.after);
    g->meshChecks = combine(combine(ev.storedAfter.connectivity, ev.storedAfter.geometry), ev.coverStoredAfter.status);
    if (ev.storedAfter.geometry != Status::Pass) g->findings.push_back("stitch mesh geometry " + std::string(toString(ev.storedAfter.geometry)));
    if (ev.coverStoredAfter.status != Status::Pass) g->findings.push_back("coverage of the knit region " + std::string(toString(ev.coverStoredAfter.status)) + ": " + ev.coverStoredAfter.detail);

    // 2. The yarn on the virtual faces.
    step("yarn", 1, 3);
    t0 = Clock::now();
    const std::vector<std::vector<int>>& faces = ev.shift.corners;
    YarnSettings ys;
    ys.diameter = D;
    YarnEmbedding raw;
    if (!embedYarn(ev.after, faces, ev.reference, s.wale, s.course, ys, raw, &err)) return refuse("no yarn: " + err);
    g->msYarn = msSince(t0);
    cancelCheckpoint();
    const Vec3 centre = part.centre;
    auto previewOf = [&](const YarnEmbedding& e) { return std::make_shared<const TriMesh>(meshInPartFrame(in, yarnTubes(e, 0.5 * D, 8), centre)); };

    // 3. The settle on the run's numerics (the stage's device, budget and lease).
    t0 = Clock::now();
    SettleHooks hooks;
    hooks.onCall = [&step](int call, int calls, double) { step("settling", call, calls); };
    hooks.lastSolverCall = [&run]() {
        const std::vector<std::string> calls = gpu::ncpNumericsCalls(run.numerics);
        return calls.empty() ? std::string() : calls.back();
    };
    std::ostringstream report, log;
    const SettleSettings settle = guidedSettleSettings(s);
    const SettleResult r = settleYarn({raw, ev.after, faces, ev.reference, ys, in.pieceName}, settle, run.numerics, hooks, report, log);
    g->msSettle = msSince(t0);
    g->report = report.str();
    g->settleVerdict = r.verdict;
    g->generatedLengthMm = r.raw.totalLength;
    g->generatedCloserThanD = fractionCloserThanD(r.raw);
    g->generatedTightBends = r.raw.kinkSamples;
    const char* previewOnly = "the generated yarn, as a preview only (not a knit result: no fibres, no export)";
    const bool solved = r.end == SettleEnd::Solved, drafted = r.end == SettleEnd::ConvertedOnly, skipped = r.end == SettleEnd::ConstructionFailedNotSettled;
    if (!solved && !drafted && !skipped) {
        for (size_t i = 0; i < r.raw.defects.size() && i < 4; ++i) g->findings.push_back(r.raw.defects[i]);
        g->previewTubes = previewOf(raw);
        g->yarnShown = previewOnly;
        g->yarnExported = "nothing (not knitted)";
        return refuse(r.end == SettleEnd::YarnDisconnected ? "the generated yarn is broken, so no stitch can be made from it (" + r.verdict + "); see the findings"
                                                         : "the yarn construction is not supported here (" + r.verdict + ")");
    }
    g->constructionBefore = "generated yarn " + constructionText(r.raw) + "; as curves " + constructionText(r.converted);
    if (r.rawConstructionFailed) addCrossingFindings(g->findings, "generated yarn", r.raw, ev.after, faces);
    if (r.conversionAddedDefects) {
        g->findings.push_back("the conversion to curves added defects (" + constructionText(r.raw) + " -> " + constructionText(r.converted) +
                              "): the viewport shows the generated yarn, the export writes the curves");
        addCrossingFindings(g->findings, "as curves", r.converted, ev.after, faces);
    }
    g->startD = r.D0;
    g->startVerified = r.startVerified;
    g->calls = static_cast<int>(r.stages.size());
    g->callsCertified = r.certifiedStages;
    g->reachedD = 0;
    for (const SettleStage& st : r.stages)
        if (st.solver.rfind("relax threw", 0) != 0) g->reachedD = st.D;
    g->gateLine = r.gateFinal.line();
    g->facesChecked = r.settled.facesChecked;
    g->facesFailBothFrames = r.letterSwapped;
    g->facesFailOneFrame = r.projectionChanged;
    g->interlocks = r.settled.interlocksVerified;
    g->interlocksExpected = r.settled.interlocksExpected;
    g->retainedLengthMm = r.settled.totalLength;
    g->retainedCloserThanD = fractionCloserThanD(r.settled);
    g->retainedTightBends = r.settled.kinkSamples;
    if (r.refusedCalls > 0 || r.budgetHit) {
        g->outcome = GuidedOutcome::Stopped;
        part.stop = r.refusedCalls > 0 ? KnitStop::DeviceBudget : KnitStop::TimeCap;
        char b[160];
        if (r.refusedCalls > 0) std::snprintf(b, sizeof b, "the solver refused the call at D %.3f mm: ", r.stages.back().D);
        else std::snprintf(b, sizeof b, "the settle's time budget ran out at D %.3f mm", g->reachedD);
        g->reason = std::string(b) + (r.refusedCalls > 0 ? r.refusal : std::string());
        g->previewTubes = previewOf(raw);
        g->yarnShown = previewOnly;
        g->yarnExported = "nothing (not knitted)";
        part.note = "guided stitch mesh stopped: " + g->reason;
        part.knit.reset();
        return;
    }

    // 4. Stitch records from the stitch mesh, then everything in the part's frame.
    CurveSet curves = r.curves;
    GuidedRecordCounts counts;
    if (!guidedStitchRecords(ev.after, faces, raw, r.conversion, ev.reference, r.raw.interlocksExpected, curves, counts, &err)) {
        g->previewTubes = previewOf(raw);
        g->yarnShown = previewOnly;
        g->yarnExported = "nothing (not knitted)";
        return refuse(err);
    }
    g->stitches = counts.loops;
    g->increases = counts.increases;
    g->decreases = counts.decreases;
    curvesInPartFrame(in, centre, curves);
    g->curvesDigest = curvesDigest(curves);
    // the surface the yarn lies on (its backing in the export): the whole piece, or, when the piece is knitted as several
    // knitting components, only this component's knit region, so its backing never covers another component's yarn
    const TriMesh& lies = pieceHasSeveralKnitComponents(in) ? ev.knitRegion : ev.reference;
    part.mesh = std::make_shared<const TriMesh>(meshInPartFrame(in, lies, centre));

    // The verdict: a check that did not run (the target gate without a solve) is UNRESOLVED, never a pass; a wrong crossing
    // order is a FAIL whatever else happened.
    const Status gate = solved ? gateStatus(r.gateFinal) : Status::Unresolved;
    const Status construction = r.constructionFailed() ? Status::Fail : solved ? constructionStatus(r) : Status::Unresolved;
    const Status verdict = combine(combine(gate, construction), g->meshChecks);
    const int wrong = facesWronglyCrossed(r.rawConstructionFailed ? r.raw : r.converted);
    const GuidedOutcome afterSettle = outcomeAfterSettle(gate, construction, g->meshChecks);  // used only when the solver settled it
    if (solved) {
        g->yarnShown = g->yarnExported = "the curves the settle retained";
    } else {
        g->yarnExported = "the generated yarn as curves, not settled (a diagnostic export)";
        g->yarnShown = g->yarnExported;
        if (r.conversionAddedDefects) {
            g->previewTubes = previewOf(raw);
            g->yarnShown = "the generated yarn before its conversion to curves (the conversion added defects)";
        }
    }
    if (drafted) {
        g->outcome = GuidedOutcome::Drafted;
        g->reason = "Draft: the generated yarn, not settled (a Final knit settles it)";
        if (r.constructionFailed())
            g->reason = "Draft: the generated yarn, not settled; " + std::to_string(wrong) +
                        " stitches have a wrong crossing order (listed below), so a Final knit will not settle it either";
    } else if (skipped) {
        g->outcome = GuidedOutcome::SettleSkipped;
        g->reason = "Final skipped the settle: " + std::to_string(wrong) + " stitches of the " + (r.rawConstructionFailed ? "generated yarn" : "yarn converted to curves") +
                    " have a wrong crossing order (listed below), and the solver is not asked to repair one; shown and exported: the generated yarn, not settled";
    } else if (afterSettle == GuidedOutcome::Settled) {
        g->outcome = GuidedOutcome::Settled;
        char b[160];
        std::snprintf(b, sizeof b, "the target gate passes at D %.3f mm, the construction is kept and the stitch mesh checks pass", D);
        g->reason = b;
    } else if (afterSettle == GuidedOutcome::SettledMeshFails) {
        g->outcome = GuidedOutcome::SettledMeshFails;
        char b[200];
        std::snprintf(b, sizeof b, "the target gate passes at D %.3f mm and the construction is kept, but the stitch mesh does not fit the piece everywhere (its checks %s)",
                      D, toString(g->meshChecks));
        g->reason = b;
    } else {
        g->outcome = GuidedOutcome::NotSettled;
        std::string why;
        if (gate != Status::Pass)
            why += "the yarn still overlaps or bends too tightly (" + std::to_string(r.gateFinal.failRegions) + " failing and " + std::to_string(r.gateFinal.unresolvedRegions) +
                   " unresolved regions)";
        if (construction != Status::Pass) why += std::string(why.empty() ? "" : "; ") + std::to_string(r.letterSwapped) + " stitches could not be confirmed as knitted correctly";
        if (g->meshChecks != Status::Pass) why += std::string(why.empty() ? "" : "; ") + "the stitch mesh does not fit the piece everywhere";
        g->reason = why;
    }
    if (!solved && g->meshChecks != Status::Pass) g->reason += "; the stitch mesh does not fit the piece everywhere";
    auto k = std::make_unique<PartKnit>();
    k->surface = std::make_unique<MeshSurface>(*part.mesh);
    k->t.prog.stitches = curves.topo.stitches;
    k->t.curves = std::move(curves);
    k->t.pp = PatternParams{};
    k->t.pp.yarnDiameter = D;
    k->t.pp.clearance = s.clearance;
    k->t.pp.waleSpacing = s.wale;
    k->t.pp.courseSpacing = s.course;
    k->t.pp.depthAmplitude = s.depthAmplitude;
    k->t.detail = "guided stitch mesh: " + std::to_string(g->faces) + " faces, " + std::to_string(counts.loops) + " stitches (" + std::to_string(counts.increases) + " INC, " +
                  std::to_string(counts.decreases) + " DEC)";
    k->after = r.gateFinal.report;
    k->after.verdict = gate;  // curvature and clearance; the signature does not apply (its line says so)
    guidedReportLines(k->after, *g, construction, solved ? std::string() : constructionText(r.converted));  // the curves the knit holds
    k->verdict = verdict;
    k->gateAOnCpu = true;
    k->detail = "guided: " + std::string(toString(g->outcome)) + ": " + g->reason;
    k->msValidate = g->msSettle;
    part.yarnD = D;
    part.knit = std::move(k);
}

} // namespace knit
