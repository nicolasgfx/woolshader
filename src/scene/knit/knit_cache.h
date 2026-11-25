#pragma once
// The per-piece knit cache (.knitcache; the spiral pieces, spiral_cache.h): one file per piece holding the knitting stage's result
// (program, layout, curves with their derived tables, final pattern parameters, pole flags, the
// diameter attempts and the knitting log), so the later stages (gate (d), plies, fibres,
// checks, export) can be rerun without knitting, settling and repairing again.
//
// Exactness: everything is stored as raw bytes (doubles bit for bit). Gate (a) is not stored:
// it is rerun on the cached curves and compared with the digest of the knitting run's final
// report (a mismatch refuses the cache). The piece mesh, chart and surface are rebuilt from the
// sources and checked against the stored mesh digest and pole flags.
//
// Invalidation: `key` lists every input of the knitting stage as "name = value" lines (source
// content hashes, options, generator version, build configuration); a cache whose key differs
// is refused with the first differing line ("stale: <name>"). The diameter fallback list is not
// part of the key: the stored attempts (each tried D with its gate (a) verdict) decide whether
// another list selects the same kept attempt (see selectCachedAttempt).
#include <cstdint>
#include <string>
#include <vector>

#include "io/files.h"
#include "knitting/knit/piece_knit.h"
#include "model/params.h"
#include "model/status.h"

namespace knit {

// Bumped whenever knitting, settle, repair or validate output changes.
constexpr int kKnitCacheVersion = 1;

struct KnitCacheAttempt {
    double D = 0;
    Status verdict = Status::Unresolved;  // gate (a) after settle / repair
};

struct KnitCacheEntry {
    std::string key;                       // "name = value" lines
    PieceKnit knit;
    PatternParams pp;                      // the kept attempt's parameters (D after the fallback)
    bool poleAInside = false, poleBInside = false;
    std::vector<KnitCacheAttempt> attempts;  // in order; the last one is the kept result
    bool residue = false;                  // every D failed: kept = the requested D (residue)
    bool generatorChecks = false;          // the ValidateOptions of the final gate (a)
    uint64_t meshDigest = 0;               // meshDigest(piece mesh as loaded)
    uint64_t validationDigest = 0;         // fnv1a64(validationSnapshot(final gate (a)))
    std::string log;                       // the knitting stage's log lines
    double msKnit = 0;                     // wall time of the original knitting stage
};

bool writeKnitCache(const std::string& path, const KnitCacheEntry& e, std::string* err);
bool readKnitCache(const std::string& path, KnitCacheEntry& e, std::string* err);

// First line of `want` that `have` lacks or has with another value ("" if the keys are equal).
std::string staleKeyField(const std::string& have, const std::string& want);
// Two knit diameters are one attempt when they differ by at most this. The knit D is the drawn D
// times k, and k is rounded to 0.001 (yarnFit), so drawn diameters that knit the same stitches can
// differ by tenths of a micrometre: 1.86 x 0.626 = 1.16436 and the stored 1.164 (0.36 um apart).
// Fallback diameters lie 30-40 um apart, so 1 um never takes one attempt for another.
constexpr double kKnitDiameterToleranceMm = 1e-3;
bool sameKnitDiameter(double a, double b);
// The kept attempt for the diameter list `dList` (requested D first, then the fallbacks below
// it, as the pipeline tries them; every D failing appends the requested D as residue): true and
// the index of the attempt the loop would keep, or false and the first D the cache cannot
// decide (`missing`). Diameters match within kKnitDiameterToleranceMm.
bool selectCachedAttempt(const KnitCacheEntry& e, const std::vector<double>& dList, double requestedD, double* missing);

// The knit cache's raw encoding (doubles bit for bit), shared by the other format of knitted curves
// (knit.bin, saved_knit.h).
void writeCurves(ByteWriter& w, const CurveSet& cs);
void readCurves(ByteReader& r, CurveSet& cs);
void writePatternParams(ByteWriter& w, const PatternParams& p);
PatternParams readPatternParams(ByteReader& r);
void writeStitches(ByteWriter& w, const std::vector<StitchRecord>& v);
std::vector<StitchRecord> readStitches(ByteReader& r);

} // namespace knit
