#pragma once
// The tests' CPU knitting numerics: the whole port (numerics.h) without a device. Production has one
// implementation (gpu::makeKnitNumerics); this is a test double, never linked into a tool. Where the
// device's batches are a float32 evaluation of a certified CPU algorithm, this is that algorithm.
#include "knitting/numerics.h"
#include "knitting/relax/settle.h"
#include "knitting/tube/yarn_surface.h"
#include "knitting/layout/chart_lattice.h"

namespace knit_test {

class CpuNumerics final : public knit::KnitNumerics {
public:
    // exec: for the port's calls that take none (gate (d), the boundary checks).
    explicit CpuNumerics(const knit::ExecPolicy& exec = {}) : exec_(exec) {}
    const char* name() const override { return "cpu"; }
    knit::ValidationReport validate(const knit::CurveSet& cs, const knit::PatternParams& pp, const knit::ValidateOptions& opt) override;
    knit::CurvatureResult curvature(const knit::YarnPath& y, double radius, const knit::Allowances& al, const knit::ExecPolicy& exec) override;
    void relax(const knit::RelaxRequest& request) override;
    void store(knit::CurveSet& cs, const knit::Controls& P, const knit::ExecPolicy& exec) override;
    std::unique_ptr<knit::MeshSurface> meshSurface(const knit::TriMesh& mesh) override;
    std::unique_ptr<knit::BaseSurface> pieceQuery(const knit::Chart& chart) override;
    void classifyLattice(const knit::UvAtlas& atlas, const knit::LatticeChartQuery& q, const knit::LatticeSolids& solids,
                         std::vector<knit::LatticeNodeClass>& out) override;
    knit::SurfaceCheck checkRenderedSurface(const knit::CurveSet& cs, const knit::PatternParams& ppR, const knit::SurfaceParams& sp) override;
    knit::BoundaryCheck checkBoundary(const knit::CurveSet& cs, const knit::SurfaceBinding& bind, const knit::MeshSurface& surf, double Rs) override;

private:
    knit::ExecPolicy exec_;
};

// The process's serial CPU numerics.
knit::KnitNumerics& cpuNumerics();

// The boundary checks B1 and B2 (BoundaryCheck, base_surface.h), certified; the written specification
// of gpu::checkBoundary. exec: samples, B1 intervals and B2 projections in parallel; the result is
// counts and a minimum, so it is identical to the serial one.
struct BoundaryCheckSteps {
    double step = 0.1, minStep = 1e-3;
    knit::ExecPolicy exec;
};
knit::BoundaryCheck checkBoundary(const knit::CurveSet& cs, const knit::SurfaceBinding& bind, const knit::MeshSurface& surf, double Rs,
                                  const BoundaryCheckSteps& steps = {});

// Forwards every call to another numerics: the base of the doubles that override one method.
class ForwardingNumerics : public knit::KnitNumerics {
public:
    explicit ForwardingNumerics(knit::KnitNumerics& inner) : inner_(inner) {}
    const char* name() const override { return inner_.name(); }
    knit::ValidationReport validate(const knit::CurveSet& cs, const knit::PatternParams& pp, const knit::ValidateOptions& opt) override {
        return inner_.validate(cs, pp, opt);
    }
    knit::CurvatureResult curvature(const knit::YarnPath& y, double radius, const knit::Allowances& al, const knit::ExecPolicy& exec) override {
        return inner_.curvature(y, radius, al, exec);
    }
    void relax(const knit::RelaxRequest& request) override { inner_.relax(request); }
    void store(knit::CurveSet& cs, const knit::Controls& P, const knit::ExecPolicy& exec) override { inner_.store(cs, P, exec); }
    std::unique_ptr<knit::MeshSurface> meshSurface(const knit::TriMesh& mesh) override { return inner_.meshSurface(mesh); }
    std::unique_ptr<knit::BaseSurface> pieceQuery(const knit::Chart& chart) override { return inner_.pieceQuery(chart); }
    void classifyLattice(const knit::UvAtlas& atlas, const knit::LatticeChartQuery& q, const knit::LatticeSolids& solids,
                         std::vector<knit::LatticeNodeClass>& out) override {
        inner_.classifyLattice(atlas, q, solids, out);
    }
    knit::SurfaceCheck checkRenderedSurface(const knit::CurveSet& cs, const knit::PatternParams& ppR, const knit::SurfaceParams& sp) override {
        return inner_.checkRenderedSurface(cs, ppR, sp);
    }
    knit::BoundaryCheck checkBoundary(const knit::CurveSet& cs, const knit::SurfaceBinding& bind, const knit::MeshSurface& surf, double Rs) override {
        return inner_.checkBoundary(cs, bind, surf, Rs);
    }

protected:
    knit::KnitNumerics& inner_;
};

// Settle's whole-state contract under test: remembers a span of controls and their projection states
// before the first round's relaxation, and sabotages the result of that relaxation (the span mirrored
// through the surface, its states poisoned: that stitch's signature breaks). The change goes into the
// relaxed P, so settle's own P carries it; `relaxed` and `relaxedBinding` are the sabotaged result.
class SabotagedRound final : public ForwardingNumerics {
public:
    SabotagedRound(knit::KnitNumerics& inner, int yarn, int firstControl, int endControl)
        : ForwardingNumerics(inner), yarn_(yarn), c0_(firstControl), c1_(endControl) {}
    void relax(const knit::RelaxRequest& request) override;
    std::vector<knit::Vec3> controlsBefore;
    std::vector<knit::ProjState> statesBefore;
    std::vector<knit::ProjResult> resultsBefore;
    knit::Controls relaxed;
    knit::SurfaceBinding relaxedBinding;

private:
    int yarn_, c0_, c1_;
};

// Repair's rejected attempt under test: every relaxation's result is judged invalid (the validation
// that follows it reports FAIL), so no attempt is ever accepted.
class RejectEveryAttempt final : public ForwardingNumerics {
public:
    explicit RejectEveryAttempt(knit::KnitNumerics& inner) : ForwardingNumerics(inner) {}
    void relax(const knit::RelaxRequest& request) override;
    knit::ValidationReport validate(const knit::CurveSet& cs, const knit::PatternParams& pp, const knit::ValidateOptions& opt) override;
    int attempts = 0;

private:
    bool rejectNext_ = false;
};

} // namespace knit_test
