// One sewn piece knitted as one spiral (spiral_part.h). knitSpiralPart reads as the list of its stages:
// the chart, the knitting at the planned diameters (or the knit cache), the CPU audit of the GPU
// curves, the knit cache entry and gate (d); each stage is one function.
#include "scene/knit/spiral_part.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include "model/cancel.h"
#include "knitting/layout/chart.h"
#include "knitting/tube/yarn_surface.h"
#include "scene/knit/knit_cache.h"
#include "knitting/layout/stitch_layout.h"
#include "model/log.h"
#include "model/parallel.h"
#include "model/vec.h"
#include "knitting/knit/piece_knit.h"
#include "knitting/relax/settle.h"
#include "gpu/knitting.h"
#include "gpu/surfaces.h"
#include "scene/knit/knit_report.h"
#include "scene/knit/knit_run.h"
#include "scene/run_log.h"
#include "scene/scene.h"
#include "io/files.h"
#include "scene/stage_params.h"
#include "settings/settings_mapping.h"

namespace knit {

ScenePieceRequest pieceRequest(const KnitRun& run, const PieceWork& piece) {
    return {{piece.part.mesh.get(), piece.part.centre, piece.part.object}, piece.part.rings, run.spirals};
}

PieceKnitParams pieceKnitParams(const KnitRun& run, const PieceWork& piece) {
    PieceKnitParams kp;
    kp.g.w = run.settings.wale;
    kp.g.h = run.settings.course;
    kp.hand = handRequest(run.settings, piece.stem);
    return kp;
}

void PieceLog::profile(const StageClock& clock, const std::string& label) {
    if (stageProfile) appendFormat(text, "     profile: %s\n", clock.line(label.c_str()).c_str());
}

namespace {

// -v: the chart's asymmetry about the pole axis, from per-meridian arc lengths.
void logMeridianAsymmetry(const Chart& ch, PieceLog& log) {
    double tMin = 1e300, tMax = 0, hMin = 1e300, hMax = 0;
    for (int j = 0; j < 32; ++j) {
        const double b = (j + 0.5) / 32;
        const double tot = ch.meridianLength(b, 0, 1), half = ch.meridianLength(b, 0, 0.5);
        tMin = std::min(tMin, tot); tMax = std::max(tMax, tot);
        hMin = std::min(hMin, half); hMax = std::max(hMax, half);
    }
    appendFormat(log.text, "     meridians: pole-to-pole %.1f..%.1f mm, pole A to a = 0.5 %.1f..%.1f mm\n", tMin, tMax, hMin, hMax);
}

// Stage: the piece's chart and the chart check (in the log).
void chartPiece(const KnitRun& run, PieceWork& piece) {
    PieceLog& log = piece.log;
    const TriMesh& m = *piece.part.mesh;
    Chart ch(m, piece.part.rings);
    ChartCheck cc = ch.check();
    const Chart::RingAt longest = ch.longestRing();
    const double Lmax = longest.length;
    double meridian = 0;
    for (int j = 0; j < 8; ++j) meridian += ch.meridianLength(j / 8.0, 0, 1) / 8;
    const double w = run.settings.wale, h = run.settings.course;
    double area = m.area(), est = area / (w * h);
    appendFormat(log.text, "[%2zu] %-24s %-8s %-12s area %8.0f mm2  longest ring %6.1f mm (%3.0f st)  meridian %6.1f mm "
            "(%3.0f rows)  ~%5.0f stitches  chart %s  (%.0f ms)\n",
            piece.k, piece.stem.c_str(), piece.object.c_str(), piece.material.c_str(), area, Lmax, Lmax / w, meridian,
            meridian / h, est, toString(cc.overall()), msSince(piece.started));
    if (run.options.verbose || cc.overall() != Status::Pass) appendFormat(log.text, "     %s\n", cc.detail.c_str());
    const ChartSkew skew = chartSkew(m);
    appendFormat(log.text, "     rows and columns: %.1f deg off perpendicular on average, %.1f deg at the 95th percentile\n", skew.mean, skew.p95);
    log.profile(piece.clock, "load + chart");
    if (run.options.verbose) logMeridianAsymmetry(ch, log);
}

// How a piece's stage ended: go on with the next stage, stop the piece (its log says why; not an
// error of the run), or fail the run.
enum class StageEnd { Next, StopPiece, Error };

// One knitting attempt at yarn diameter D: knit, the program and layout checks, gate (a), then
// settle and (where gate a fails) the certified repair. Appends the attempt to out.attempts. False
// if the piece cannot be knitted (the reason in piece.stopReason and the log).
bool knitAttempt(const KnitRun& run, PieceWork& piece, const PieceKnitParams& kp, double D, PieceKnitResult& out) {
    const SceneSettings& s = run.settings;
    const ExecPolicy& inner = run.options.inner;
    PieceLog& log = piece.log;
    KnitAttempt& a = out.kept;
    std::string err;
    bool attemptGenChecks = true;
    a.pp = PatternParams{};
    a.pp.yarnDiameter = D;
    a.pp.clearance = s.clearance;
    a.pp.depthAmplitude = s.depthAmplitude;
    a.spk = ScenePieceKnit{};
    const auto t1 = std::chrono::steady_clock::now();
    char dLabel[32];
    std::snprintf(dLabel, sizeof dLabel, "D %.2f ", D);
    const StageClock spAttempt;
    StageClock spStep;
    if (!knitScenePiece(pieceRequest(run, piece), {a.pp, kp}, run.numerics, a.spk, &err)) {
        appendFormat(log.text, "     knit: %s\n", err.c_str());
        piece.stopReason = err;
        return false;
    }
    PieceKnit& pk = a.spk.knit;
    const PieceSurface& surf = *a.spk.surface;
    out.msKnit = msSince(t1);
    log.profile(spStep, dLabel + std::string("knit"));
    spStep = StageClock{};
    const auto t2 = std::chrono::steady_clock::now();
    a.gc = checkProgram(pk.prog);
    a.lc = checkLayout(pk.prog, pk.lay, kp.g);
    a.pp.waleSpacing = s.wale;
    a.pp.courseSpacing = s.course;
    a.r = run.numerics.validate(pk.curves, a.pp, ValidateOptions{.generatorChecks = true, .exec = inner});
    out.msValidate = msSince(t2);
    log.profile(spStep, dLabel + std::string("program + layout + gate a"));
    if (s.settle) {
        int f0 = 0;
        for (const ClearanceRegion& g : a.r.clearance.regions) f0 += g.status == Status::Fail;
        double rho0 = 1e300;
        for (const CurvatureResult& cr : a.r.curvature) rho0 = std::min(rho0, cr.rhoLo);
        appendFormat(log.text, "     before settle: %s | local signature %d/%zu | rho_lo %.3f mm | clearance %d FAIL regions\n",
                toString(a.r.verdict), a.r.localSignature.pass, a.r.localSignature.pairs.size(), rho0, f0);
        SettleParams sp;
        sp.exec = inner;
        const StageClock spSettle;
        SettleReport sr = settle({pk.curves, a.pp, &surf, nullptr, a.spk.hand.rest.empty() ? nullptr : &a.spk.hand.rest}, sp, run.numerics, &a.r);
        attemptGenChecks = false;
        log.profile(spSettle, dLabel + std::string("settle (") + std::to_string(sr.rounds) + " rounds)");
        std::string per;
        for (int d : sr.defectsPerRound) per += " " + std::to_string(d);
        appendFormat(log.text, "     settle: %d rounds, defect stitches per round:%s; last free set %d stitches; max move %.3f mm; "
                "%.0f ms\n", sr.rounds, per.c_str(), sr.freeStitches, sr.maxMove, sr.ms);
        out.msValidate += sr.ms;
    }
    if (s.repair && a.r.verdict != Status::Pass) {
        const StageClock spRepair;
        RepairParams rpar;
        rpar.settle.exec = inner;
        const ValidationReport current = a.r;  // gate (a) of the curves as settle (or the knit) left them
        RepairReport rr = certifiedRepair({pk.curves, a.pp, &surf}, rpar, run.numerics, &a.r, &current);
        attemptGenChecks = false;
        const char* moveLabel = !rr.moveTried ? "no move" : rr.move.accepted ? "move kept" : "move not kept";
        log.profile(spRepair, dLabel + std::string("repair (") + moveLabel + ")");
        appendFormat(log.text, "     repair: %s; %.0f ms\n", rr.detail.c_str(), rr.ms);
        out.msValidate += rr.ms;
    }
    log.profile(spAttempt, dLabel + std::string("attempt, gate a ") + toString(a.r.verdict));
    out.attempts.push_back({D, a.r.verdict});
    a.genChecks = attemptGenChecks;
    return true;
}

// Stage (no knit cache): knits the piece at the planned diameters until gate (a) passes. If every
// diameter fails, the result is the requested diameter's, labelled residue (a consistent look).
StageEnd knitWithDiameterFallback(const KnitRun& run, PieceWork& piece, const PieceKnitParams& kp, const DiameterPlan& plan,
                                  PieceKnitResult& out) {
    const double requested = plan.attempt.front();
    PieceLog& log = piece.log;
    std::vector<double> dList = plan.attempt;
    // The first attempt, kept while fallbacks run: the residue retry at the requested D has the
    // same inputs, so it reuses it instead of knitting again (identical result).
    std::unique_ptr<KnitAttempt> first;
    for (size_t di = 0; di < dList.size(); ++di) {
        if (di) appendFormat(log.text, "     D %.2f mm: gate a %s after settle/repair -> retry with D %.2f mm\n", dList[di - 1],
                        toString(out.kept.r.verdict), dList[di]);
        if (di > 0 && first && dList[di] == dList[0]) {
            out.kept = std::move(*first);
            first.reset();
            out.attempts.push_back({dList[di], out.kept.r.verdict});
            appendFormat(log.text, "     D %.2f mm: the first attempt's result (identical inputs, not knitted again)\n", dList[di]);
            break;
        }
        if (!knitAttempt(run, piece, kp, dList[di], out)) return StageEnd::StopPiece;
        if (out.kept.r.verdict == Status::Pass) break;
        if (di + 1 == dList.size() && dList.size() > 1 && dList.back() != requested) {
            appendFormat(log.text, "     no diameter passes gate a: keeping D %.2f mm (labelled residue)\n", requested);
            dList.push_back(requested);
            out.residue = true;
        }
        const KnitAttempt& a = out.kept;
        if (di == 0 && dList.size() > 1) first.reset(new KnitAttempt{std::move(out.kept.spk), a.r, a.pp, a.gc, a.lc, a.genChecks});
    }
    if (!plan.notAttempted.empty() && out.kept.r.verdict != Status::Pass) {
        std::string list;
        for (double d : plan.notAttempted) {
            char b[16];
            std::snprintf(b, sizeof b, " %.2f", d);
            list += b;
        }
        appendFormat(log.text, "     diameter search stopped by budget (knit.diameter_attempts = %d): not attempted D%s mm; the result is kept with its "
                "failures (labelled residue; not a claim that these diameters would fail)\n",
                run.settings.diameterAttempts, list.c_str());
        out.residue = true;
    }
    return StageEnd::Next;
}

// Stage: the piece's yarn, from the knit cache or knitted at the scene's diameters; then the
// knitting backends' statistics, the CPU audit of the GPU curves and the knit cache entry.
StageEnd knitPiece(const KnitRun& run, PieceWork& piece, PieceKnitResult& out, std::string* err) {
    PieceLog& log = piece.log;
    const PieceKnitParams kp = pieceKnitParams(run, piece);
    const DiameterPlan plan = partDiameters(run.settings);
    const size_t logStart = log.text.size();
    const auto tKnitStage = std::chrono::steady_clock::now();
    if (!run.cache.loadDir.empty()) {
        if (!loadPieceFromCache(run, piece, plan, out, err)) return StageEnd::Error;
    } else {
        const StageEnd end = knitWithDiameterFallback(run, piece, kp, plan, out);
        if (end != StageEnd::Next) return end;
    }
    KnitAttempt& a = out.kept;
    if (!out.cached) {
        appendFormat(log.text, "     %s (this run so far)\n", gpu::knitNumericsStats(run.numerics).c_str());
        if (run.options.stageProfile)
            for (const std::string& call : gpu::ncpNumericsCalls(run.numerics)) appendFormat(log.text, "       profile: ncp call %s\n", call.c_str());
        appendFormat(log.text, "     %s (this run so far)\n", gpu::knitPlacementStats(run.numerics).c_str());
        // Offline CPU audit of the GPU curves (reported separately; never selects or corrects).
        const StageClock spAudit;
        const ValidationReport rGpu = a.r;  // the GPU diagnostic gate (a) of the kept attempt
        a.r = validate(a.spk.knit.curves, a.pp, ValidateOptions{.generatorChecks = a.genChecks, .exec = run.options.inner});
        appendFormat(log.text, "     gate a (GPU diagnostic, float32): %s; CPU audit of the GPU curves: %s\n", toString(rGpu.verdict), toString(a.r.verdict));
        log.profile(spAudit, "CPU audit of the GPU knitting (gate a, not part of the GPU time)");
    }
    if (!run.cache.saveDir.empty() && !out.cached) savePieceToCache(run, piece, out, log.text.substr(logStart), msSince(tKnitStage));
    return StageEnd::Next;
}

// Stage, gate (d): the final surface (FinalSurfaceSpec) at 0.025 mm rings, with strand.surface_half also
// at half the ring step, checked on the GPU by the independent sampled surface check. False if the
// GPU check fails (an error of the run).
bool checkFinalSurface(const KnitRun& run, PieceWork& piece, const CurveSet& curves, const PatternParams& ppR, const FinalSurfaceSpec& fs,
                       std::string* err) {
    PieceLog& log = piece.log;
    for (int pass = 0; pass < (run.settings.surfaceHalf ? 2 : 1); ++pass) {
        auto ts = std::chrono::steady_clock::now();
        const std::string passLabel = pass ? "gate d (half step) " : "gate d ";
        StageClock spGate;
        SurfaceParams sp;
        sp.compression = fs.compression;
        sp.maxCompression = fs.maxCompression;
        sp.ringSpacing = pass == 0 ? 0.025 : 0.0125;
        SurfaceCheck sc;
        gpu::SurfaceGpuStats sst;
        std::string e;
        if (!gpu::checkSurface({curves, ppR, sp}, {}, sc, &sst, &e)) {
            if (err) *err = piece.stem + ": gate d (GPU): " + e;
            return false;
        }
        const size_t rings = static_cast<size_t>(sst.rings);
        log.profile(spGate, passLabel + "on the GPU (" + sst.line() + ")");
        piece.gateD = combine(piece.gateD, combine(sc.fold, sc.collisions));
        appendFormat(log.text, "     gate d%s: folds %s (min rho/support %.3f, %d rings), collisions %s (%ld points, %ld "
                "penetrating, max depth %.4f mm); %s; model %.4f mm, %zu rings (%.0f MB); %.0f ms\n",
                pass ? " (half step)" : "", toString(sc.fold), sc.minFoldRatio, sc.foldRings,
                toString(sc.collisions), sc.pointsTested, sc.penetratingPoints, sc.maxDepth,
                sc.resolution.c_str(), sp.ringSpacing, rings, rings * sizeof(Ring) / 1048576.0, msSince(ts));
        if (sc.penetratingPoints > 0 || sc.fold != Status::Pass)  // where: arc positions along the yarns [mm]
            appendFormat(log.text, "     gate d%s worst: collision yarn %d s %.2f into yarn %d s %.2f; fold yarn %d s %.2f; %d rings clamped\n",
                    pass ? " (half step)" : "", sc.worstYarnA, sc.worstSA, sc.worstYarnB, sc.worstSB, sc.worstFoldYarn, sc.worstFoldS,
                    sc.compressionClamped);
    }
    return true;
}

// -v: the knitting's defects (printDefects, straight to stdout) and the local-signature failures.
void logKnitDefects(const KnitAttempt& a, PieceLog& log) {
    const PieceKnit& pk = a.spk.knit;
    printDefects(pk.curves, a.r, a.pp);
    for (const LocalPairResult& pair : a.r.localSignature.pairs)
        if (pair.status != Status::Pass) {
            const StitchRecord& st = pk.curves.topo.stitches[pair.child];
            appendFormat(log.text, "    signature: stitch %d (%s, turn %d of %zu, index %d): %s\n", pair.child, toString(st.op),
                    st.turn, pk.prog.turns.size(), st.index, pair.detail.c_str());
        }
}

// -v: the layout's out-of-bounds cells by op and turn position (first 3 / last 3 turns = poles).
void logLayoutCells(const PieceKnit& pk, double w, double h, PieceLog& log) {
    std::map<std::string, int> cat;
    const int T = static_cast<int>(pk.prog.turns.size());
    for (const StitchRecord& st : pk.prog.stitches) {
        const LayoutNode& nd = pk.lay.nodes[st.id];
        const double wq = 0.5 * (nd.wL + nd.wR) / w, hq = nd.h / h;
        const double ang = std::acos(std::clamp(dot(nd.c, nd.v), -1.0, 1.0)) * 180 / kPi;
        std::string what;
        if (wq < 0.75) what += "narrow ";
        if (wq > 1.3) what += "wide ";
        if (st.nParents && hq < 0.75) what += "short ";
        if (st.nParents && hq > 1.3) what += "tall ";
        if (st.nParents && (ang < 60 || ang > 120)) what += "skewed ";
        if (what.empty()) continue;
        const char* where = st.turn < 3 ? "start" : st.turn >= T - 3 ? "end" : "mid";
        const char* role = st.nChildren == 2 ? "INCp" : (st.nChildren == 1 && pk.prog.stitches[st.children[0]].op == StitchOp::Dec) ? "DECp" : "";
        ++cat[std::string(toString(st.op)) + role + " " + where + ": " + what];
    }
    for (const auto& [key, n] : cat) appendFormat(log.text, "       layout cell %-40s %d\n", key.c_str(), n);
}

// A pole's opening for the log (poleOpening): in mm, or none where no yarn comes within a yarn diameter of the
// surface's level there (the pole sits far inside the piece's last turn).
std::string openingText(double opening) {
    if (opening > 1e9) return "none near the surface";
    char text[32];
    std::snprintf(text, sizeof text, "%.2f mm", opening);
    return text;
}

// The piece's knitting summary at the end of its log: the stitch program, the layout, the poles and
// gate (a).
void logKnitSummary(const KnitRun& run, PieceWork& piece, const PieceKnitResult& knit) {
    PieceLog& log = piece.log;
    const KnitAttempt& a = knit.kept;
    const ValidationReport& r = a.r;
    const PieceKnit& pk = a.spk.knit;
    const PieceSurface& surf = *a.spk.surface;
    const bool inA = a.spk.poleAInside, inB = a.spk.poleBInside, startAtA = pk.startAtA;
    const double holeStart = poleHoleRadius(pk.curves, startAtA ? surf.poleA() : surf.poleB(), a.pp.radius());
    const double holeEnd = poleHoleRadius(pk.curves, startAtA ? surf.poleB() : surf.poleA(), a.pp.radius());
    int clrF = 0, clrU = 0;
    for (const ClearanceRegion& g : r.clearance.regions) {
        clrF += g.status == Status::Fail;
        clrU += g.status == Status::Unresolved;
    }
    double rhoLo = 1e300;
    for (const CurvatureResult& cr : r.curvature) rhoLo = std::min(rhoLo, cr.rhoLo);
    appendFormat(log.text, "     knit: %s; %zu stitches, %d INC, %d DEC; graph %s\n", pk.detail.c_str(), pk.prog.stitches.size(),
            a.gc.inc, a.gc.dec, toString(a.gc.status));
    appendFormat(log.text, "     layout %s: %s\n", toString(a.lc.status), a.lc.detail.c_str());
    if (run.options.verbose) logLayoutCells(pk, run.settings.wale, run.settings.course, log);
    appendFormat(log.text, "     poles: end pole %s another piece (A %s, B %s); hole radius start %.2f mm, end %.2f mm (yarn "
            "surface to pole)\n", (startAtA ? inB : inA) ? "inside" : "NOT inside", inA ? "in" : "out",
            inB ? "in" : "out", holeStart, holeEnd);
    // Where the poles are in the scene (mm), to match a hole in a picture to its pole, and the openings seen from outside.
    const Vec3 pA = surf.poleA() + piece.part.centre, pB = surf.poleB() + piece.part.centre;
    const Vec3 startPole = startAtA ? surf.poleA() : surf.poleB(), endPole = startAtA ? surf.poleB() : surf.poleA();
    const double openStart = poleOpening(pk.curves, startPole, surf.project(startPole).n, a.pp.yarnDiameter);
    const double openEnd = poleOpening(pk.curves, endPole, surf.project(endPole).n, a.pp.yarnDiameter);
    appendFormat(log.text, "     pole positions: A (%.1f, %.1f, %.1f) B (%.1f, %.1f, %.1f) mm; start %s; rows from %.2f mm off the "
            "start pole to %.2f mm off the end pole; opening start %s, end %s\n", pA.x, pA.y, pA.z, pB.x, pB.y, pB.z,
            startAtA ? "A" : "B", pk.startArc, surf.totalArc() - pk.endArc, openingText(openStart).c_str(), openingText(openEnd).c_str());
    appendFormat(log.text, "     gate a (CPU audit of the GPU curves): %s | local signature %d/%zu (%d FAIL, %d UNRES) | rho_lo %.3f mm | clearance %d "
            "FAIL, %d UNRES regions | %.0f ms knit, %.0f ms validate\n",
            toString(r.verdict), r.localSignature.pass, r.localSignature.pairs.size(), r.localSignature.fail,
            r.localSignature.unresolved, rhoLo, clrF, clrU, knit.msKnit, knit.msValidate);
}

// A knitted piece as its part's knit in the common result. The mesh surface gives the fibres their
// outward normals (as a chart part's); a spiral has no atlas, lattice or projection binding.
std::unique_ptr<PartKnit> spiralKnit(KnitAttempt& kept, const PieceKnitResult& knit, Status gateD, const TriMesh& mesh) {
    auto k = std::make_unique<PartKnit>();
    k->surface = std::make_unique<MeshSurface>(mesh);
    k->t.prog = std::move(kept.spk.knit.prog);
    k->t.curves = std::move(kept.spk.knit.curves);
    k->t.pp = kept.pp;
    k->t.detail = kept.spk.knit.detail;
    k->verdict = kept.r.verdict;  // a spiral piece's verdict is its gate (a)
    k->after = std::move(kept.r);
    k->gateAOnCpu = true;
    k->gateDStatus = gateD;
    k->detail = "spiral: " + k->t.detail;
    k->msTile = knit.msKnit;
    k->msValidate = knit.msValidate;
    return k;
}

} // namespace

bool knitSpiralPart(const KnitRun& run, size_t index, ScenePart& part, std::string* err) {
    PieceLog log;
    log.stageProfile = run.options.stageProfile;
    const SourcePrimitive& g = run.source.prims.at(part.prim);
    PieceWork piece{index, part, g.node, g.material, part.name, log};
    // The log goes to stdout as one block however the piece ends.
    struct PrintAtEnd {
        const PieceLog& log;
        ~PrintAtEnd() {
            std::fputs(log.text.c_str(), stdout);
            std::fflush(stdout);
        }
    } printAtEnd{log};
    chartPiece(run, piece);

    PieceKnitResult knit;
    const StageEnd knitEnd = knitPiece(run, piece, knit, err);
    if (knitEnd == StageEnd::Error) return false;
    if (knitEnd == StageEnd::StopPiece) {
        part.note = "not knitted: " + piece.stopReason;
        return true;
    }
    KnitAttempt& kept = knit.kept;
    appendFormat(log.text, "     yarn diameter %.2f mm\n", kept.pp.yarnDiameter);
    const FinalSurfaceSpec fspec = finalSurface(run.settings, kept.pp.yarnDiameter);
    if (std::string why; !fspec.valid(&why)) {
        if (err) *err = part.name + ": " + why;
        return false;
    }
    const PatternParams ppR = fspec.rendered(kept.pp);
    if (fspec.inflate > 0) appendFormat(log.text, "     rendered diameter D_r %.3f mm (inflation %.3f), compression on\n", ppR.yarnDiameter, fspec.inflate);
    if ((run.settings.surfaceGate || run.settings.surfaceHalf) && !checkFinalSurface(run, piece, kept.spk.knit.curves, ppR, fspec, err))
        return false;
    if (run.options.verbose) logKnitDefects(kept, log);
    logKnitSummary(run, piece, knit);
    log.profile(piece.clock, "piece total");
    part.yarnD = kept.pp.yarnDiameter;
    if (knit.residue) part.note = "residue: no tried diameter passes gate (a)";
    part.knit = spiralKnit(kept, knit, piece.gateD, *part.mesh);
    return true;
}

} // namespace knit
