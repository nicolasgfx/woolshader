// The per-piece knit cache of a run (spiral_cache.h).
#include "scene/knit/spiral_cache.h"

#include <cstdio>
#include <filesystem>

#include "gpu/knitting.h"
#include "scene/knit/knit_cache.h"
#include "knitting/knit/piece_knit.h"
#include "knitting/relax/settle.h"
#include "model/compute_caps.h"
#include "model/hash.h"
#include "model/log.h"
#include "scene/knit/knit_run.h"
#include "scene/run_log.h"
#include "scene/knit/spiral_part.h"
#include "scene/stage_params.h"

namespace knit {

namespace {

// The knit cache key: every input of the knitting (knit_cache.h), one "<name> = <value>" per line. The
// spiral's shaping is the knitting's fixed defaults (piece_knit.h, settle.h); they stay in the key,
// which knit caches written since the knit cache existed share line by line. The requested diameter
// and the fallback list are not in it: the cached attempts decide (selectCachedAttempt). False: the
// manifest or a piece file cannot be hashed.
bool buildKnitCacheKey(const SceneSettings& s, const SpiralSource& src, std::string& key, std::string* err) {
    const PieceKnitParams shaping;
    appendFormat(key, "generator = knit-cache %d\n", kKnitCacheVersion);
    appendFormat(key, "build = %s\n", buildConfiguration().c_str());
    bool okHash = true;
    appendFormat(key, "manifest = %016llx\n", static_cast<unsigned long long>(fnv1aFile(src.manifestPath, &okHash)));
    for (const PieceEntry& pe : src.man.pieces)
        appendFormat(key, "piece %s = %016llx\n", pe.file.c_str(), static_cast<unsigned long long>(fnv1aFile(piecePath(src.man, pe.file), &okHash)));
    if (!okHash) {
        if (err) *err = "knit cache: cannot hash the manifest or a piece file";
        return false;
    }
    appendFormat(key, "sigma = %a\n", src.sigma);
    appendFormat(key, "gauge w h = %a %a\n", s.wale, s.course);
    appendFormat(key, "end gap = %a\n", shaping.endGap);
    // An end whose pole shows is gathered (piece_knit.h, EndClosure); one inside another piece stays threaded.
    appendFormat(key, "ends = threaded inside another piece, else gathered: last turn %d loops off the pole by its cast arc, ring n (2 D + c) / 2 pi\n",
                 shaping.castCount);
    appendFormat(key, "close = %d %a %a\n", 1, shaping.sh.closeY, shaping.sh.closeDepth);  // every piece is closed
    appendFormat(key, "wale pull = %a\n", shaping.walePull);
    appendFormat(key, "settle = %d (one round, kept whole)\n", int(s.settle));
    appendFormat(key, "repair = %d (one certified move from the report's start diameter)\n", int(s.repair));
    if (s.hand.inUse())
        for (const SettingField& f : settingFields())
            if (f.key.rfind("knit.hand.", 0) == 0) appendFormat(key, "%s = %s\n", f.key.c_str(), settingValue(s, f).c_str());
    if (s.hand.inUse()) appendFormat(key, "hand rest = %s\n", kHandRestRule);
    appendFormat(key, "local signature = crossing attribution window %a spans\n", kAttributionWindow);
    if (s.diameterAttempts > 0) appendFormat(key, "diameter attempts = %d (budget)\n", s.diameterAttempts);
    // The producer: scenes knit on the GPU.
    appendFormat(key, "knit backend = %s\n", "gpu float32 (gpu knit v1; stored digest: CPU audit)");
    appendFormat(key, "knit solver = %s\n", gpu::ncpSolverIdentity().c_str());
    // The clearance and depth amplitude are the settings' (at their defaults the line is the one of
    // the pattern defaults).
    const PatternParams d0;
    appendFormat(key, "pattern defaults = %a %a %a %a %a %a %d %a %a %a\n", s.clearance, s.depthAmplitude, d0.headHalfWidth, d0.legHalfWidth,
                 d0.sampleTolerance, d0.transNeckY, int(d0.pattern), d0.transNeckZ, d0.transLegY, 0.0);
    return true;
}

} // namespace

std::optional<SpiralKnitCache> spiralKnitCache(const SceneSettings& s, const SpiralSource& src, const std::string& saveDir, const std::string& loadDir,
                                               std::string* err) {
    SpiralKnitCache out;
    // Under the producer: caches of the retired CPU and legacy GPU solvers stay apart.
    const std::string producer = "/gpu-f32/ncp";
    if (!saveDir.empty()) out.saveDir = saveDir + producer;
    if (!loadDir.empty()) out.loadDir = loadDir + producer;
    if (out.saveDir.empty() && out.loadDir.empty()) return out;
    if (!buildKnitCacheKey(s, src, out.key, err)) return std::nullopt;
    if (!out.saveDir.empty()) std::filesystem::create_directories(out.saveDir);
    return out;
}

bool loadPieceFromCache(const KnitRun& run, PieceWork& piece, const DiameterPlan& plan, PieceKnitResult& out, std::string* err) {
    const PieceKnitParams kp = pieceKnitParams(run, piece);
    const StageClock spLoad;
    const std::string path = run.cache.loadDir + "/" + piece.stem + ".knitcache";
    KnitCacheEntry ce;
    auto refuse = [&path, err](const std::string& why) {
        if (err) *err = path + ": knit cache refused (" + why + ")";
        return false;
    };
    std::string e;
    if (!readKnitCache(path, ce, &e)) return refuse(e);
    // The producer first (a GPU cache is never taken for a CPU request, nor the reverse).
    if ((ce.key.find("\nknit backend = ") != std::string::npos) != (run.cache.key.find("\nknit backend = ") != std::string::npos))
        return refuse("stale: knit backend (another producer)");
    if (const std::string f = staleKeyField(ce.key, run.cache.key); !f.empty()) return refuse("stale: " + f);
    double missing = 0;
    if (!selectCachedAttempt(ce, plan.attempt, plan.attempt.front(), &missing)) {
        std::string why;
        appendFormat(why, "stale: diameter list (knit D %.6f mm was not attempted when the cache was written; it holds", missing);
        for (const KnitCacheAttempt& x : ce.attempts) appendFormat(why, " %.6f", x.D);
        return refuse(why + ")");
    }
    KnitAttempt& a = out.kept;
    // No placement here: a cached piece is only rebuilt, with its CPU arc-length tables.
    if (!prepareScenePiece(pieceRequest(run, piece), nullptr, a.spk, &e)) return refuse(e);
    if (meshDigest(a.spk.chart->mesh()) != ce.meshDigest) return refuse("stale: piece mesh");
    if (a.spk.poleAInside != ce.poleAInside || a.spk.poleBInside != ce.poleBInside) return refuse("stale: pole tests");
    a.spk.knit = std::move(ce.knit);
    a.pp = ce.pp;
    if (kp.hand.inUse()) setYarnThickness(a.spk.knit.curves, drawStitchTensions(a.spk.knit.curves.topo, kp.hand), kp.hand);
    a.gc = checkProgram(a.spk.knit.prog);
    a.lc = checkLayout(a.spk.knit.prog, a.spk.knit.lay, kp.g);
    // Gate (a) rerun on the cached curves: it must reproduce the knitting run's report.
    a.r = validate(a.spk.knit.curves, a.pp, ValidateOptions{.generatorChecks = ce.generatorChecks, .exec = run.options.inner});
    if (fnv1a64(validationSnapshot(a.r)) != ce.validationDigest) return refuse("gate (a) on the cached curves differs from the knitting run");
    out.attempts = ce.attempts;
    for (double D : plan.attempt)
        for (const KnitCacheAttempt& x : ce.attempts)
            if (sameKnitDiameter(x.D, D) && x.D != D)
                appendFormat(piece.log.text, "     knit D %.6f mm requested, %.6f stored: the same attempt (within 1 um); the stored knit is kept\n", D, x.D);
    out.residue = ce.residue;
    a.genChecks = ce.generatorChecks;
    appendFormat(piece.log.text, "     knit: from the knit cache %s (knitted in %.0f ms; gate a rerun on the cached curves: identical); the knitting "
                  "log:\n", path.c_str(), ce.msKnit);
    piece.log.text += ce.log;
    piece.log.profile(spLoad, "knit cache: load, rebuild, gate a");
    out.cached = true;
    return true;
}

void savePieceToCache(const KnitRun& run, PieceWork& piece, const PieceKnitResult& out, const std::string& knitLog, double msKnit) {
    const KnitAttempt& a = out.kept;
    std::string err;
    KnitCacheEntry ce;
    ce.key = run.cache.key;
    ce.knit = a.spk.knit;
    ce.pp = a.pp;
    ce.poleAInside = a.spk.poleAInside;
    ce.poleBInside = a.spk.poleBInside;
    ce.attempts = out.attempts;
    ce.residue = out.residue;
    ce.generatorChecks = a.genChecks;
    ce.meshDigest = meshDigest(a.spk.chart->mesh());
    ce.validationDigest = fnv1a64(validationSnapshot(a.r));
    ce.log = knitLog;
    ce.msKnit = msKnit;
    const ValidationReport again = validate(a.spk.knit.curves, a.pp, ValidateOptions{.generatorChecks = a.genChecks, .exec = run.options.inner});
    const std::string path = run.cache.saveDir + "/" + piece.stem + ".knitcache";
    if (fnv1a64(validationSnapshot(again)) != ce.validationDigest)
        appendFormat(piece.log.text, "     knit cache: NOT written (gate a of the final curves does not reproduce the knitting run's report)\n");
    else if (!writeKnitCache(path, ce, &err))
        appendFormat(piece.log.text, "     knit cache: %s\n", err.c_str());
    else
        appendFormat(piece.log.text, "     knit cache: written to %s\n", path.c_str());
}

} // namespace knit
