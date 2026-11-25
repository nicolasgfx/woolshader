// The RIM driver (rim_remesh.h) and the Wu 2018 section 4 hooks (knit_rim_hooks.h) the patched RIM sources call.
// Adapted from RIM (MIT, gaoxifeng/robust_hex_dominant_meshing; CHANGES.md): remesh() follows the order of RIM's
// batch mode (batch.cpp), retag() the first loop of its init_edge_tagging2D (extract.cpp), and the subdivision
// uses its batch-mode rule.
#include "prelude.h"

#include "rim_remesh.h"

namespace rim {

#include "dedge.h"
#include "hierarchy.h"
#include "positions.h"
#include "subdivide.h"

#include "knit_rim_hooks.h"

// meshio.cpp is not built (file I/O): hierarchy.cpp's load() references load_obj.
void load_obj(const std::string&, MatrixXu&, MatrixXf&) { throw std::runtime_error("rim: file input is not built (knit_rim)"); }

int knit_rim_rosy = 4;
Float knit_rim_aspect = 1;
uint64_t knit_rim_seed = 0;

// extract.cpp's file-scope extraction state (namespace rim::, as the wrapped unit defines it).
extern std::vector<tuple_E> mEs;
extern std::vector<std::vector<uint32_t>> mFs2D;
extern std::vector<std::vector<uint32_t>> FEs;
extern MatrixXf mQ_copy2D, mO_copy2D, mN_copy2D, mV_copy2D;
void construct_Es_FEs();

namespace {

struct WuState {
    bool enabled = false;
    int rounds = 0, diagonals = 0, loopSplits = 0;
} wu;

using EdgeKey = std::pair<uint32_t, uint32_t>;
EdgeKey key(uint32_t a, uint32_t b) { return a < b ? EdgeKey(a, b) : EdgeKey(b, a); }

Vector3f polygonNormal(const std::vector<uint32_t>& f, const MatrixXf& P) {  // Newell
    Vector3f n = Vector3f::Zero();
    for (size_t k = 0; k < f.size(); ++k) n += Vector3f(P.col(f[k])).cross(Vector3f(P.col(f[(k + 1) % f.size()])));
    const Float l = n.norm();
    return l > 0 ? Vector3f(n / l) : n;
}

// The corner angle of polygon f at its k-th vertex [rad].
Float cornerAngle(const std::vector<uint32_t>& f, size_t k, const MatrixXf& P) {
    const size_t n = f.size();
    Vector3f a = P.col(f[(k + n - 1) % n]) - P.col(f[k]), b = P.col(f[(k + 1) % n]) - P.col(f[k]);
    const Float la = a.norm(), lb = b.norm();
    if (la <= 0 || lb <= 0) return 0;
    return std::acos(std::max((Float) -1, std::min((Float) 1, a.dot(b) / (la * lb))));
}

// Wu 2018 section 4 (2): every polygon with more than 5 sides is split along the diagonal "most aligned with the
// orientation field", recursively. Alignment of a diagonal from u to v: its direction projected into u's
// tangent plane against the field's two lattice axes at u (q and n x q: the position lattice is 4-PoSy for
// either RoSy), max(|d . q|, |d . n x q|). Diagonals that are already edges are skipped; ties keep the first
// (vertex order). A polygon without a usable diagonal stays as it is. Returns the number of splits.
int splitLarge(std::vector<std::vector<uint32_t>>& F, const MatrixXf& P, const MatrixXf& Q, const MatrixXf& N) {
    std::set<EdgeKey> edges;
    for (const auto& f : F)
        for (size_t k = 0; k < f.size(); ++k) edges.insert(key(f[k], f[(k + 1) % f.size()]));
    int splits = 0;
    std::vector<std::vector<uint32_t>> out;
    out.reserve(F.size());
    for (const auto& f0 : F) {
        if (f0.size() <= 5) {
            out.push_back(f0);
            continue;
        }
        std::vector<std::vector<uint32_t>> work{f0};
        while (!work.empty()) {
            std::vector<uint32_t> f = std::move(work.back());
            work.pop_back();
            const size_t n = f.size();
            if (n <= 5) {
                out.push_back(std::move(f));
                continue;
            }
            Float best = -1;
            size_t bj = 0, bk = 0;
            for (size_t j = 0; j < n; ++j)
                for (size_t k = j + 2; k < n; ++k) {
                    if (j == 0 && k == n - 1) continue;  // adjacent
                    if (edges.count(key(f[j], f[k]))) continue;
                    const Vector3f nj = N.col(f[j]);
                    Vector3f d = P.col(f[k]) - P.col(f[j]);
                    d -= nj * nj.dot(d);
                    const Float l = d.norm();
                    if (!(l > 0)) continue;
                    d /= l;
                    Vector3f q = Q.col(f[j]);
                    q -= nj * nj.dot(q);
                    if (!(q.norm() > 0)) continue;
                    q.normalize();
                    const Float s = std::max(std::abs(d.dot(q)), std::abs(d.dot(nj.cross(q))));
                    if (s > best) best = s, bj = j, bk = k;
                }
            if (best < 0) {
                out.push_back(std::move(f));
                continue;
            }
            std::vector<uint32_t> a(f.begin() + bj, f.begin() + bk + 1), b(f.begin() + bk, f.end());
            b.insert(b.end(), f.begin(), f.begin() + bj + 1);
            edges.insert(key(f[bj], f[bk]));
            ++splits;
            work.push_back(std::move(b));
            work.push_back(std::move(a));
        }
    }
    F.swap(out);
    return splits;
}

// The edges' RIM tags from the position field, as the first loop of edge_tagging2D.
void retag(MultiResolutionHierarchy& h) {
    for (auto& e : mEs) {
        const uint32_t v0 = std::get<0>(e), v1 = std::get<1>(e);
        const Vector3f q0 = mQ_copy2D.col(v0), q1 = applyRotationKeep(q0, mN_copy2D.col(v0), mQ_copy2D.col(v1), mN_copy2D.col(v1));
        const std::tuple<short, Float, Vector2f> a = posy2D_completeInfo(mO_copy2D.col(v0), q0, mN_copy2D.col(v0), mV_copy2D.col(v0), mO_copy2D.col(v1), q1,
                                                                         mN_copy2D.col(v1), mV_copy2D.col(v1), h.mScale, h.mInvScale);
        std::get<3>(e) = std::get<1>(a);
        std::get<4>(e) = std::get<0>(a);
    }
}

} // namespace

void knit_wu_loop(MultiResolutionHierarchy& h) {
    if (!wu.enabled) return;
    ++wu.rounds;
    {  // (2)
        std::vector<std::vector<uint32_t>> F;
        for (const auto& f : mFs2D)
            if (f.size() >= 3) F.push_back(f);
        const int s = splitLarge(F, mO_copy2D, mQ_copy2D, mN_copy2D);
        if (s) {
            mFs2D = F;
            construct_Es_FEs();
            retag(h);
            wu.loopSplits += s;
        }
    }
    // (1): the edge ids and faces by vertex pair (independent of the incremental FEs / E_pfs upkeep).
    std::map<EdgeKey, uint32_t> edgeOf;
    for (uint32_t i = 0; i < mEs.size(); ++i) edgeOf[key(std::get<0>(mEs[i]), std::get<1>(mEs[i]))] = i;
    std::map<EdgeKey, std::vector<uint32_t>> facesOf;
    for (uint32_t f = 0; f < mFs2D.size(); ++f)
        for (size_t k = 0; k < mFs2D[f].size(); ++k) facesOf[key(mFs2D[f][k], mFs2D[f][(k + 1) % mFs2D[f].size()])].push_back(f);
    const Float right = (Float) (M_PI / 2);
    for (uint32_t f = 0; f < mFs2D.size(); ++f) {
        const std::vector<uint32_t>& t = mFs2D[f];
        if (t.size() != 3) continue;
        const Vector3f nf = polygonNormal(t, mO_copy2D);
        int best = -1;
        Float bestDev = 0;
        for (int k = 0; k < 3; ++k) {  // the edge opposite vertex k
            const EdgeKey e = key(t[(k + 1) % 3], t[(k + 2) % 3]);
            const auto fs = facesOf.find(e);
            if (fs == facesOf.end() || fs->second.size() != 2) continue;  // border
            const uint32_t g = fs->second[0] == f ? fs->second[1] : fs->second[0];
            if (nf.dot(polygonNormal(mFs2D[g], mO_copy2D)) < 0) continue;  // sharp feature
            const Float dev = std::abs(cornerAngle(t, k, mO_copy2D) - right);
            if (best < 0 || dev < bestDev) best = k, bestDev = dev;
        }
        if (best < 0) continue;
        const auto it = edgeOf.find(key(t[(best + 1) % 3], t[(best + 2) % 3]));
        if (it == edgeOf.end()) continue;
        short& tag = std::get<4>(mEs[it->second]);
        if (tag == Edge_tag::B) {  // a kept edge becomes a diagonal (a collapsing one stays collapsing)
            tag = Edge_tag::D;
            ++wu.diagonals;
        }
    }
}

} // namespace rim

namespace knit_rim {

namespace {

// Nearest-point lookup on a uniform grid (cell `cell`), for the output normals and axes and the fixed field's
// subdivision vertices.
struct Grid {
    double cell = 1;
    std::map<std::array<long long, 3>, std::vector<int>> bins;
    const rim::MatrixXf* P = nullptr;
    Grid(const rim::MatrixXf& pts, double c) : cell(c), P(&pts) {
        for (int i = 0; i < (int) pts.cols(); ++i) bins[at(pts(0, i), pts(1, i), pts(2, i))].push_back(i);
    }
    std::array<long long, 3> at(double x, double y, double z) const {
        return {(long long) std::floor(x / cell), (long long) std::floor(y / cell), (long long) std::floor(z / cell)};
    }
    // The closest point in the box of cells within r of the query's cell (-1: none).
    int scan(double x, double y, double z, long long r) const {
        const std::array<long long, 3> c = at(x, y, z);
        int best = -1;
        double bd = 0;
        for (long long i = c[0] - r; i <= c[0] + r; ++i)
            for (long long j = c[1] - r; j <= c[1] + r; ++j)
                for (long long k = c[2] - r; k <= c[2] + r; ++k) {
                    const auto it = bins.find({i, j, k});
                    if (it == bins.end()) continue;
                    for (int v : it->second) {
                        const double dx = (*P)(0, v) - x, dy = (*P)(1, v) - y, dz = (*P)(2, v) - z, d = dx * dx + dy * dy + dz * dz;
                        if (best < 0 || d < bd || (d == bd && v < best)) best = v, bd = d;
                    }
                }
        return best;
    }
    // Exact: a point found within box r lies closer than sqrt(3) (r + 1) cells <= 4 r cells.
    int nearest(double x, double y, double z) const {
        for (long long r = 1; r <= 64; r *= 2)
            if (scan(x, y, z, r) >= 0) return scan(x, y, z, 4 * r);
        return -1;
    }
};

int parseAfter(const std::string& s, const char* what) {
    const size_t at = s.find(what);
    return at == std::string::npos ? 0 : std::atoi(s.c_str() + at + std::strlen(what));
}

// RIM prints to std::cout; the caller gets the text in Result::log.
struct Redirect {
    std::streambuf* old;
    explicit Redirect(std::streambuf* b) : old(std::cout.rdbuf(b)) {}
    ~Redirect() { std::cout.rdbuf(old); }
};

// Thrown where the deadline is checked; caught by remesh(), which reports the stage.
struct DeadlinePassed {
    const char* stage;
};

using Clock = std::chrono::steady_clock;
double msSince(Clock::time_point t) { return std::chrono::duration<double, std::milli>(Clock::now() - t).count(); }

void toRim(const std::vector<std::array<double, 3>>& pos, const std::vector<std::array<int, 3>>& tri, rim::MatrixXf& V, rim::MatrixXu& F) {
    V.resize(3, pos.size());
    for (size_t i = 0; i < pos.size(); ++i) V.col(i) = rim::Vector3f(pos[i][0], pos[i][1], pos[i][2]);
    F.resize(3, tri.size());
    for (size_t t = 0; t < tri.size(); ++t)
        for (int k = 0; k < 3; ++k) F(k, t) = (uint32_t) tri[t][k];
}

// RIM's own subdivision rule (the batch mode's): when the longest edge exceeds half the spacing or twice the
// mean edge, edges are split down to min(spacing / 2, 2 x mean). True when it ran.
bool subdivideInput(rim::MultiResolutionHierarchy& h, rim::MatrixXf& V, rim::MatrixXu& F, double spacing) {
    using namespace rim;
    const MeshStats stats = h.compute_mesh_stats(F, V, true);
    const Float s = (Float) spacing;
    if (!(stats.mMaximumEdgeLength * 2 > s || stats.mMaximumEdgeLength > stats.mAverageEdgeLength * 2)) return false;
    VectorXu V2E, E2E;
    VectorXb boundary, nonManifold;
    build_dedge(F, V, V2E, E2E, boundary, nonManifold);
    subdivide(F, V, V2E, E2E, boundary, nonManifold, std::min(s / 2, (Float) stats.mAverageEdgeLength * 2), true);
    return true;
}

} // namespace

int predictedVertices(const std::vector<std::array<double, 3>>& pos, const std::vector<std::array<int, 3>>& tri, double spacing) {
    if (pos.empty() || tri.empty() || !(spacing > 0)) return (int) pos.size();
    std::ostringstream log;
    const Redirect redirect(log.rdbuf());
    rim::MultiResolutionHierarchy h;
    rim::MatrixXf V;
    rim::MatrixXu F;
    toRim(pos, tri, V, F);
    subdivideInput(h, V, F, spacing);
    return (int) V.cols();
}

bool remesh(const std::vector<std::array<double, 3>>& pos, const std::vector<std::array<int, 3>>& tri, const Options& o, Result& out, std::string* err) {
    using namespace rim;
    out = Result{};
    if (pos.empty() || tri.empty() || !(o.spacing > 0) || (o.rosy != 2 && o.rosy != 4)) {
        if (err) *err = "rim remesh: empty mesh, no spacing or rosy not 2 / 4";
        return false;
    }
    if (o.aspect != 1 && (!(o.aspect > 0) || o.rosy != 2 || o.extrinsic)) {
        if (err) *err = "rim remesh: a rectangular lattice (aspect != 1) needs the intrinsic 2-RoSy field";
        return false;
    }
    if (!o.orientation.empty() && o.orientation.size() != pos.size()) {
        if (err) *err = "rim remesh: the fixed orientation needs one direction per input vertex";
        return false;
    }
    std::ostringstream log;
    const Redirect redirect(log.rdbuf());
    const Clock::time_point start = Clock::now();
    const bool capped = o.deadlineMs > 0;
    try {
        Clock::time_point t = start;
        MultiResolutionHierarchy mRes;
        mRes.mV.resize(1);
        toRim(pos, tri, mRes.mV[0], mRes.mF);
        if (o.subdivide) subdivideInput(mRes, mRes.mV[0], mRes.mF, o.spacing);
        const Float s = (Float) o.spacing;
        // load()'s bookkeeping after reading a mesh.
        mRes.mAABB = AABB(mRes.mV[0].rowwise().minCoeff(), mRes.mV[0].rowwise().maxCoeff());
        mRes.ms = mRes.compute_mesh_stats(mRes.mF, mRes.mV[0], true);
        mRes.diagonalLen = 3 * (mRes.mAABB.max - mRes.mAABB.min).norm() / 100;
        knit_rim_seed = o.seed;
        mRes.build();
        mRes.setScale(s / mRes.diagonalLen);  // mScale = diagonalLen x ratio
        out.inputVertices = (int) mRes.mV[0].cols();
        for (int l = 0; l < mRes.levels(); ++l) out.levelSizes.push_back((int) mRes.mV[l].cols());
        out.msHierarchy = msSince(t);
        if (capped && msSince(start) > o.deadlineMs) throw DeadlinePassed{"hierarchy"};
        t = Clock::now();
        knit_rim_rosy = o.rosy;
        knit_rim_aspect = (Float) o.aspect;
        wu = WuState{};
        wu.enabled = o.wuLoop;
        const int top = std::max(0, mRes.levels() - 2);
        if (!o.orientation.empty()) {
            // Fixed: the finest level from the given directions (subdivision vertices: the nearest input
            // vertex's), every coarser level from its first child, projected into the tangent planes.
            MatrixXf orig(3, pos.size());
            for (size_t i = 0; i < pos.size(); ++i) orig.col(i) = Vector3f(pos[i][0], pos[i][1], pos[i][2]);
            const Grid g(orig, std::max((double) mRes.mAverageEdgeLength, 1e-9) * 2);
            for (uint32_t i = 0; i < (uint32_t) mRes.mV[0].cols(); ++i) {
                const int src = i < pos.size() ? (int) i : g.nearest(mRes.mV[0](0, i), mRes.mV[0](1, i), mRes.mV[0](2, i));
                Vector3f q(o.orientation[src][0], o.orientation[src][1], o.orientation[src][2]);
                const Vector3f n = mRes.mN[0].col(i);
                q -= n * n.dot(q);
                if (q.norm() > 0) mRes.mQ[0].col(i) = q.normalized();
            }
            for (int l = 0; l + 1 < mRes.levels(); ++l) {
                std::vector<char> set(mRes.mV[l + 1].cols(), 0);
                const SMatrix& P = mRes.mP[l];
                for (int k = 0; k < P.outerSize(); ++k)
                    for (SMatrix::InnerIterator it(P, k); it; ++it) {
                        if (set[it.col()]) continue;
                        set[it.col()] = 1;
                        Vector3f q = mRes.mQ[l].col(it.row());
                        const Vector3f n = mRes.mN[l + 1].col(it.col());
                        q -= n * n.dot(q);
                        if (q.norm() > 0) mRes.mQ[l + 1].col(it.col()) = q.normalized();
                    }
            }
        } else {
            for (int l = top; l >= 0; --l) {
                if (l < top) mRes.prolongOrientations(l);
                for (int it = 0; it < o.sweepsPerLevel; ++it) {
                    mRes.smoothOrientationsTri((uint32_t) l, true, o.randomization, o.extrinsic);
                    if (capped && msSince(start) > o.deadlineMs) throw DeadlinePassed{"orientations"};
                }
            }
        }
        out.msOrientations = msSince(t);
        if (o.fieldsOnly) {
            for (uint32_t i = 0; i < (uint32_t) mRes.mV[0].cols(); ++i)
                out.fieldQ.push_back({mRes.mQ[0](0, i), mRes.mQ[0](1, i), mRes.mQ[0](2, i)});
            mRes.detectOrientationSingularitiesTri();
            for (int i = 0; i < (int) mRes.mOrientationSingularities.cols(); ++i)
                out.orientationSingularPoints.push_back({mRes.mOrientationSingularities(0, i), mRes.mOrientationSingularities(1, i), mRes.mOrientationSingularities(2, i)});
            knit_rim_aspect = 1;
            out.log = log.str();
            return true;
        }
        t = Clock::now();
        for (int l = top; l >= 0; --l) {
            if (l < top) mRes.prolongPositions(l);
            for (int it = 0; it < o.sweepsPerLevel; ++it) {
                mRes.smoothPositionsTri((uint32_t) l, true, o.randomization, o.extrinsic);
                if (capped && msSince(start) > o.deadlineMs) throw DeadlinePassed{"positions"};
            }
        }
        for (uint32_t i = 0; i < (uint32_t) mRes.mV[0].cols(); ++i) {
            out.fieldQ.push_back({mRes.mQ[0](0, i), mRes.mQ[0](1, i), mRes.mQ[0](2, i)});
            out.fieldO.push_back({mRes.mO[0](0, i), mRes.mO[0](1, i), mRes.mO[0](2, i)});
        }
        out.msPositions = msSince(t);
        t = Clock::now();
        // The singularities as points (RIM's matrices: a point, a normal and a colour per column).
        mRes.detectOrientationSingularitiesTri();
        mRes.detectPositionSingularitiesTri();
        for (int i = 0; i < (int) mRes.mOrientationSingularities.cols(); ++i)
            out.orientationSingularPoints.push_back({mRes.mOrientationSingularities(0, i), mRes.mOrientationSingularities(1, i), mRes.mOrientationSingularities(2, i)});
        for (int i = 0; i < (int) mRes.mPositionSingularities.cols(); ++i)
            out.positionSingularPoints.push_back({mRes.mPositionSingularities(0, i), mRes.mPositionSingularities(1, i), mRes.mPositionSingularities(2, i)});
        if (capped && msSince(start) > o.deadlineMs) throw DeadlinePassed{"extraction (not started)"};
        // meshExtraction2D with the flags of RIM's batch mode.
        mRes.re_color = mRes.splitting = mRes.decomposes = mRes.doublets = mRes.triangles = true;
        if (!mRes.meshExtraction2D()) throw std::runtime_error("rim extraction: the input is not manifold");
        out.msExtraction = msSince(t);
        out.wuRounds = wu.rounds;
        out.wuDiagonals = wu.diagonals;
        out.wuLoopSplits = wu.loopSplits;
        const MatrixXf& V = mRes.mV_tag;
        // Normals and axes: the nearest finest-level vertex's (RIM's area-free face average, unit; q in its tangent plane).
        const Grid g(mRes.mV[0], std::max(o.spacing, 1e-9));
        for (uint32_t i = 0; i < (uint32_t) V.cols(); ++i) {
            out.pos.push_back({V(0, i), V(1, i), V(2, i)});
            const int nv = g.nearest(V(0, i), V(1, i), V(2, i));
            const Vector3f n = nv >= 0 ? Vector3f(mRes.mN[0].col(nv)) : Vector3f::Zero();
            Vector3f q = nv >= 0 ? Vector3f(mRes.mQ[0].col(nv)) : Vector3f::Zero();
            q -= n * n.dot(q);
            if (q.norm() > 0) q.normalize();
            out.nrm.push_back({n.x(), n.y(), n.z()});
            out.axis.push_back({q.x(), q.y(), q.z()});
        }
        for (const auto& f : mRes.F_tag) {
            std::vector<int> poly;
            for (uint32_t v : f) poly.push_back((int) v);
            out.faces.push_back(std::move(poly));
        }
    } catch (const DeadlinePassed& d) {
        knit_rim_aspect = 1;
        out.deadlineHit = true;
        out.log = log.str() + "\ntime cap in " + d.stage + "\n";
        if (err) *err = std::string("rim remesh: time cap of ") + std::to_string((int) o.deadlineMs) + " ms reached in " + d.stage;
        return false;
    } catch (const std::exception& e) {
        knit_rim_aspect = 1;
        out.log = log.str();
        if (err) *err = std::string("rim remesh: ") + e.what();
        return false;
    }
    knit_rim_aspect = 1;
    out.log = log.str();
    out.genusIn = parseAfter(out.log, "Input genus: ");
    out.genusOut = parseAfter(out.log, "Output genus: ");
    return true;
}

} // namespace knit_rim
