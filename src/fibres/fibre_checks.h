#pragma once
// Everything that inspects finished fibres; none of it changes them.
//
//   ExportTubeIndex       the exported polygonal tubes, indexed once and queried per fibre;
//   fibreReferenceCheck   an independent proof against the round tubes;
//   rootAttachmentCheck   roots and bridges against the exported tubes;
//   measureFibreAngles    the physical angles of accepted fibres to their yarn;
//   snapshots             exact text of generation and check results, for regression comparison.
//
// The common reference is the yarn centreline as a polyline (per span, pieces short enough that the
// sagitta kappa_hi c^2 / 8 <= 1e-3 mm with the certified span curvature bound of checkCurvature), so
// every polyline point is within e_c of the curve and every curve point within e_c of the polyline.
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "fibres/fibre.h"
#include "model/parallel.h"

namespace knit {

// How the export writes the tubes: buildTubeMesh with `sides` around and rings every `ringSpacing`,
// from a surface model with `modelRingSpacing`, every vertex through `toExport` (the exporter's
// mapping incl. float32 rounding, a similarity up to that rounding).
struct ExportTubeSpec {
    double modelRingSpacing = 0, ringSpacing = kExportRing;
    int sides = kExportSides;
    bool compression = false;
    // A function, not a SceneMapping: the export's own mapping, float32 rounding included, so the checks
    // see the coordinates the file holds.
    std::function<Vec3(const Vec3&)> toExport;
    double eSer = 0;  // [mm] float32 envelope of the exported coordinates
    double maxCompression = 0.4;  // SurfaceParams::maxCompression of the tubes' model (kDefaultMaxCompression)
};

struct ExportTubeGeometry;

// The exported polygonal tubes of one part. A slab is the piece of tube between two rings. The
// CurveSet must outlive the index.
class ExportTubeIndex {
public:
    // exec: built in parallel, the same index; a parallel policy declares toExport concurrency-safe.
    // shared / sharedInputs: a surface model built elsewhere and the digest of its build inputs
    // (surfaceBuildDigest); it replaces the index's own buildSurface only if that digest equals the
    // digest of the inputs the index would build from.
    ExportTubeIndex(const CurveSet& cs, const PatternParams& pp, ExportTubeSpec spec, const ExecPolicy& exec = {},
                    const SurfaceModel* shared = nullptr, uint64_t sharedInputs = 0);
    ~ExportTubeIndex();
    // Root and bridge on the serialized geometry, for one fibre:
    //   buried      the root is inside by >= radius0 + margin [mm];
    //   singleExit  the polyline, sampled from the root, is inside the own exported tube up to one
    //               exit and outside after it;
    //   clear       every sample keeps (centreline distance - bulge) - e_c - rho >= r_i + h / 2 +
    //               e_float to every other yarn portion. The own yarn near the anchor is excluded
    //               only while l <= attach. bulge = ringSpacing^2 kappa / 8 bounds how far the
    //               exported facets leave the smooth tube on the inner side of a bend.
    struct TubeVerdict {
        bool buried = false, singleExit = false, clear = false;
        double rootDepth = -1e300;
        int slabs = 0;
        bool ok() const { return buried && singleExit && clear; }
    };
    TubeVerdict check(const FibreSet& fs, const Fibre& f, double margin) const;   // every item
    bool accepts(const FibreSet& fs, const Fibre& f, double margin) const;        // stops at the first failure
    // Depth [mm] of library point p inside the exported tube of `yarn`, over the slabs whose arc
    // meets [s0, s1] (< 0: outside all of them). `slabs` counts the searched slabs.
    double depth(int yarn, double s0, double s1, const Vec3& p, int* slabs = nullptr) const;
    double scale() const;       // export units per mm
    double rhoMax() const;      // largest ring semi-axis [mm]
    bool sharedSurface() const;  // the surface model came from `shared`
    const ExportTubeGeometry& geometry() const { return *g_; }

private:
    std::unique_ptr<const ExportTubeGeometry> g_;
};

// Independent check of a fibre's guarantee for ROUND, uncompressed tubes of radius r: "outside the
// tube" is exactly "distance to the centreline >= r", so no ring or ellipse code is used. Each
// segment of the whole polyline, from the buried root, is classified by bisection with a Lipschitz
// bound (the distance to a fixed set is 1-Lipschitz along it): an interval with midpoint m and half
// length h is
//   OK          if d_poly(m) - e_c - h >= r + r_f,
//   FAIL        if d_poly(m) + e_c < r + r_f (m itself is inside),
//   UNRESOLVED  if neither holds at h <= hMin.
// r_f is the root radius (conservative for a tapered fibre). Excluded: for fibre arc length
// l <= Fibre::attach, measured from the emergence, the own yarn within attachArc + r + r_f.
struct FibreReferenceResult {
    long long segments = 0, intervals = 0;
    int ok = 0, fail = 0, unresolved = 0;  // fibres
    double eChord = 0;                     // polyline sagitta bound used [mm]
    double minMargin = 1e300;              // min over checked intervals of d_poly(m) - e_c - h - (r + r_f)
    std::vector<int> failing, unresolvedFibres;  // fibre indices (first 100)
    double ms = 0;
};
// exec: per fibre, merged in fibre order (identical to serial).
FibreReferenceResult fibreReferenceCheck(const ExportTubeIndex& tubes, const FibreSet& fs, double hMin = 1e-4, const ExecPolicy& exec = {});

// Root attachment against the EXPORTED polygonal tube: a point is inside if it lies strictly (by
// `margin`, in export units) inside a slab around the fibre's anchor: on the inner side of all
// 2 x sides facet planes and between the two ring planes. The slabs searched follow the root's
// reach: arc in [s - W, s + W], W = |segment 1| + rho_max + ringSpacing.
struct AttachmentResult {
    int fibres = 0;
    int rootInside = 0;           // root point strictly inside the exported tube
    int tipOutside = 0;           // tip outside the own tube's slabs around the anchor
    double minRootDepth = 1e300;  // min over fibres of the root's distance inside (export units)
    std::vector<int> detached;    // fibre indices with the root not inside (first 100)
    int buried = 0, singleExit = 0, clear = 0;  // ExportTubeIndex::check's items
    int maxSlabs = 0;             // most slabs searched for one fibre
    double ms = 0;
};
AttachmentResult rootAttachmentCheck(const ExportTubeIndex& tubes, const FibreSet& fs, double margin, const ExecPolicy& exec = {});

// The PHYSICAL angle of accepted fibres to their yarn. The settings' "max angle" of directed fibres
// is the maximum NOMINAL lift, a bound on a drawn slope; the bow, the helix and the yarn's own
// curvature change the actual angle, so it is measured here and reported, never claimed as bounded.
// Angle = acos |d . T| (0..90 deg) between a segment direction d and the unit tangent T of the
// centreline of the yarn the fibre is authored on (Fibre::parentYarn):
//   emergence  the first visible segment against T at the fibre's anchor arc length s;
//   along      the maximum over the visible segments, each against T at the centreline point
//              closest to the segment midpoint (searched within +-(visible length + 1 mm) of s).
struct FibreAngleStats {
    int fibres = 0;
    double emergence[3] = {0, 0, 0};  // p50, p95, max [deg]
    double along[3] = {0, 0, 0};      // p50, p95, max of the per-fibre maximum [deg]
};
struct FibreAngles {
    FibreAngleStats pop[2];  // FibrePopulation::Directed, Random
};
FibreAngles measureFibreAngles(const CurveSet& yarns, const FibreSet& fs, const ExecPolicy& exec = {});

// Exact text of results: every deterministic field (doubles and floats as hex floats, timing
// excluded), so two results are bit-identical iff their snapshots are equal.
std::string fibreSnapshot(const FibreSet& fs);  // the report line, then one line per fibre
std::string reportSnapshot(const FibreReport& r);
std::string checkSnapshot(const FibreReferenceResult& r);
std::string checkSnapshot(const AttachmentResult& r);
// A digest of the same content in chunks of 4096 fibres, computable in parallel: FNV-1a 64 of the
// report line followed by the chunks' values in order. Not equal to fnv1a64(fibreSnapshot(fs)).
uint64_t fibreSnapshotChunkDigest(const FibreSet& fs, const ExecPolicy& exec = {});

} // namespace knit
