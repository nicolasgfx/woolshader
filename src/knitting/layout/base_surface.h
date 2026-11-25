#pragma once
// The surfaces stitches are placed on. BaseSurface is the interface; ConeSurface and PlaneSurface
// are analytic, MeshSurface is a welded triangle mesh with the ownership-preserving projection of
// yarn controls.
//
// Ownership-preserving projection: bounded continuation from the previous ACCEPTED surface location.
// A state is the accepted foot (triangle, barycentric weights, point, normal) and the query it was
// computed for. continueFrom(prev, q) walks from prev's foot by steepest descent of |q - x| over the
// surface's triangle adjacency; it returns
//   Owned         a local closest point reached within the walk budget |q - prev.query| + lambda w,
//                 without the normal turning by more than 60 deg against prev's normal at any step,
//                 not on an open-boundary edge, within the depth band;
//   Budget, NormalTurn, OpenBoundary, DepthBand   otherwise (UNRESOLVED; the state is NOT advanced).
// Claim: an Owned foot is a local closest point reached from the previous accepted foot along the
// surface under those limits, nothing more. In particular it is not claimed to be "the correct
// sheet" in general; near a tight fold the expected outcome is NormalTurn, not a silent switch to
// the other sheet.
#include <array>
#include <string>
#include <vector>

#include "model/mesh_query.h"
#include "model/status.h"
#include "model/vec.h"

namespace knit {

struct SurfaceSample {
    Vec3 p, n;  // point and unit outward normal
};

// A base surface with a chart (a, b): a runs along the pole axis (a in [aMin, aMax] is the chart
// domain), b around it with period 1.
class BaseSurface {
public:
    virtual ~BaseSurface() = default;
    virtual SurfaceSample at(double a, double b) const = 0;
    virtual SurfaceSample project(const Vec3& q) const = 0;
    virtual double aMin() const { return 0; }
    virtual double aMax() const { return 1; }
    // Batched queries, as the placement uses them (the device's in src/gpu). The defaults are the single
    // calls in order (identical results); a device surface evaluates a batch on the device. out is
    // resized to the input size.
    virtual void atBatch(const std::vector<Vec2>& ab, std::vector<SurfaceSample>& out) const;
    virtual void projectBatch(const std::vector<Vec3>& q, std::vector<SurfaceSample>& out) const;
};

// Truncated cone around +z: z = a * H, radius R0 + (R1 - R0) a (a cylinder if R0 = R1).
class ConeSurface : public BaseSurface {
public:
    ConeSurface(double R0, double R1, double H) : R0_(R0), R1_(R1), H_(H) {}
    SurfaceSample at(double a, double b) const override;
    SurfaceSample project(const Vec3& q) const override;
    double radiusAt(double a) const { return R0_ + (R1_ - R0_) * a; }
    double height() const { return H_; }

private:
    double R0_, R1_, H_;
};

// Plane z = 0, a = y / H, b = x / W (b not periodic; for flat regressions).
class PlaneSurface : public BaseSurface {
public:
    PlaneSurface(double W, double H) : W_(W), H_(H) {}
    SurfaceSample at(double a, double b) const override { return {{b * W_, a * H_, 0}, {0, 0, 1}}; }
    SurfaceSample project(const Vec3& q) const override { return {{q.x, q.y, 0}, {0, 0, 1}}; }

private:
    double W_, H_;
};

struct ProjState {
    int tri = -1;
    double b[3] = {0, 0, 0};
    Vec3 foot, n;   // foot point and interpolated unit normal
    Vec3 query;     // the point this state was computed for
    bool operator==(const ProjState& o) const;
};

enum class ProjResult : unsigned char { Owned, Budget, NormalTurn, OpenBoundary, DepthBand };
const char* toString(ProjResult r);

struct ContinuationParams {
    double lambdaW = 2.75;       // lambda * w [mm]: walk budget beyond |q - q_prev| (lambda 0.5 at w 5.5)
    int maxSteps = 4096;
    double minNormalCos = 0.5;   // cos 60 deg: normal continuity along the walk
    double band = 1e300;         // |(q - foot) . n| limit [mm]
};

class MeshSurface;

// Per-control projection states of a curve set, indexed like CurveSet::yarns[y].curve.ctrl.
// Placement fills them; settle advances them only on Owned and restores them together with the
// controls whenever it restores controls (rollback, rejected repair attempts).
struct SurfaceBinding {
    const MeshSurface* surface = nullptr;
    ContinuationParams cp;
    std::vector<std::vector<ProjState>> state;
    std::vector<std::vector<ProjResult>> result;  // last result per control
    int unresolved() const;
};

// Boundary checks of settled yarn against the mesh's OPEN boundary:
//   B1 certified: the 3D distance from every yarn centreline point to every open-boundary edge is
//      >= R_s (arc-length samples, Lipschitz with the enclosed arc length between samples; uncertain
//      intervals refined to minStep, then UNRESOLVED);
//   B2 sampled (step `step`): every sample's continued projection is Owned and not OpenBoundary.
struct BoundaryCheck {
    Status b1 = Status::Pass, b2 = Status::Pass;
    double minDistance = 1e300;     // B1: smallest sampled centreline-to-boundary distance
    int b1Fail = 0, b1Unresolved = 0;
    long b2Samples = 0;
    int b2Unresolved[5] = {0, 0, 0, 0, 0};  // by ProjResult
    std::string detail;
};

// The mesh as a base surface: project() is the UNRESTRICTED closest point (MeshQuery; diagnostics,
// compatible with pieces); continueFrom() is the ownership-preserving projection. Normals: welded
// area-weighted vertex normals, interpolated (continuous across uv seams). One instance per thread.
class MeshSurface : public BaseSurface {
public:
    explicit MeshSurface(const TriMesh& m);
    SurfaceSample at(double a, double b) const override;  // no global chart: throws
    SurfaceSample project(const Vec3& q) const override;
    // Steepest descent of |q - x| over the triangle adjacency from prev. On Owned, next is the new
    // state; otherwise next = prev (the state is not advanced).
    virtual ProjResult continueFrom(const ProjState& prev, const Vec3& q, ProjState& next, const ContinuationParams& cp) const;
    // A batch of independent continuations: the single calls in order; a device surface overrides it.
    virtual void continueBatch(const std::vector<ProjState>& prev, const std::vector<Vec3>& q, const ContinuationParams& cp,
                               std::vector<ProjState>& next, std::vector<ProjResult>& res) const;
    ProjState stateAt(int tri, const double b[3], const Vec3& query) const;
    Vec3 normalAt(int tri, const double b[3]) const;
    const TriMesh& mesh() const { return m_; }
    bool isBoundaryEdge(int tri, int edge) const { return nbr_[tri][edge] < 0; }  // edge k: corners k, k+1

private:
    const TriMesh& m_;
    MeshQuery q_;
    std::vector<std::array<int, 3>> nbr_;         // neighbour across edge k (corners k, k + 1), -1: boundary
    std::vector<std::vector<int>> vertexTris_;
};

} // namespace knit
