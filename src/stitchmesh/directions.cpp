#include "stitchmesh/directions.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <map>
#include <queue>
#include <set>
#include <utility>

#include "Highs.h"

namespace knit {

namespace {

using EdgeKey = std::pair<int, int>;
EdgeKey edgeKey(int a, int b) { return a < b ? EdgeKey(a, b) : EdgeKey(b, a); }

// A half of a row: 2 row + side (0 A, 1 B). Its role as an affine function of the row's binary y: side A is top when
// y = 1 (role = y), side B when y = 0 (role = 1 - y).
double slopeOf(int half) { return (half & 1) ? -1.0 : 1.0; }
double offsetOf(int half) { return (half & 1) ? 1.0 : 0.0; }

struct FaceGraph {  // the labelled mesh as face lists with the roles per corner edge and the face across each directed edge
    std::vector<std::vector<int>> faces;
    std::vector<std::vector<char>> course;  // per face, per corner edge: 1 course, 0 wale
    std::map<EdgeKey, std::pair<int, int>> faceOfDirected;  // u -> v: (face, corner)
    std::pair<int, int> across(int f, int k) const {
        const int n = static_cast<int>(faces[f].size());
        const auto it = faceOfDirected.find({faces[f][(k + 1) % n], faces[f][k]});
        return it == faceOfDirected.end() ? std::pair<int, int>{-1, -1} : it->second;
    }
    std::vector<int> waleCorners(int f) const {
        std::vector<int> w;
        for (size_t k = 0; k < faces[f].size(); ++k)
            if (!course[f][k]) w.push_back(static_cast<int>(k));
        return w;
    }
};

bool buildFaceGraph(const PolygonMesh& m, const EdgeLabels& labels, FaceGraph& g, std::string* err) {
    g.faces.resize(m.faces());
    g.course.resize(m.faces());
    for (int f = 0; f < m.faces(); ++f)
        for (int k = 0; k < m.sides(f); ++k) {
            const int a = m.vertexAt(f, k), b = m.vertexAt(f, (k + 1) % m.sides(f));
            g.faces[f].push_back(a);
            if (!labels.has(a, b)) {
                if (err) *err = "directions: edge " + std::to_string(a) + "-" + std::to_string(b) + " has no role";
                return false;
            }
            g.course[f].push_back(labels.at(a, b) == EdgeRole::Course ? 1 : 0);
            g.faceOfDirected[{a, b}] = {f, k};
        }
    return true;
}

// Rows: paths from their ends first, then the cycles. entry / exit: the face-local wale corner the row comes in and goes out.
struct Rows {
    std::vector<int> rowOf, entry, exit;
    std::vector<std::vector<int>> rows;
    std::vector<char> closed;
};

void walkRow(const FaceGraph& g, int start, int entryCorner, Rows& r) {
    const int row = static_cast<int>(r.rows.size());
    r.rows.emplace_back();
    r.closed.push_back(0);
    int f = start, in = entryCorner;
    while (true) {
        r.rowOf[f] = row;
        r.rows[row].push_back(f);
        r.entry[f] = in;
        int out = -1;
        for (int k : g.waleCorners(f))
            if (k != in) {
                out = k;
                break;
            }
        r.exit[f] = out;
        if (out < 0) break;
        const std::pair<int, int> next = g.across(f, out);
        if (next.first < 0) break;
        if (next.first == start && next.second == r.entry[start]) {
            r.closed[row] = 1;
            break;
        }
        if (r.rowOf[next.first] >= 0) break;
        f = next.first;
        in = next.second;
    }
}

bool buildRows(const FaceGraph& g, Rows& r, std::string* err) {
    const int nF = static_cast<int>(g.faces.size());
    r.rowOf.assign(nF, -1);
    r.entry.assign(nF, -1);
    r.exit.assign(nF, -1);
    for (int f = 0; f < nF; ++f) {
        if (r.rowOf[f] >= 0) continue;
        const std::vector<int> w = g.waleCorners(f);
        if (w.empty() || w.size() > 2) {
            if (err) *err = "directions: face " + std::to_string(f) + " has " + std::to_string(w.size()) + " wale edges";
            return false;
        }
        if (w.size() == 1) {
            walkRow(g, f, -1, r);
            continue;
        }
        for (int k : w)
            if (g.across(f, k).first < 0) {  // a wale edge on the border: a path end
                walkRow(g, f, k, r);
                break;
            }
    }
    for (int f = 0; f < nF; ++f)
        if (r.rowOf[f] < 0) walkRow(g, f, g.waleCorners(f)[0], r);
    return true;
}

// The side of each course edge of a face: B (1) counter-clockwise after the entry up to the exit, A (0) otherwise.
std::vector<std::vector<char>> sidesOf(const FaceGraph& g, const Rows& r) {
    std::vector<std::vector<char>> side(g.faces.size());
    for (size_t f = 0; f < g.faces.size(); ++f) {
        const int n = static_cast<int>(g.faces[f].size());
        side[f].assign(n, 0);
        if (r.entry[f] >= 0 && r.exit[f] >= 0) {
            for (int k = (r.entry[f] + 1) % n; k != r.exit[f]; k = (k + 1) % n) side[f][k] = 1;
        } else if (r.entry[f] >= 0) {
            side[f][(r.entry[f] + 1) % n] = 1;
        } else if (r.exit[f] >= 0) {
            side[f][(r.exit[f] + n - 1) % n] = 1;
        }
    }
    return side;
}

struct MetaGraph {
    std::vector<EdgeKey> edges;  // (half, half)
    std::vector<double> weight;
    std::vector<std::vector<std::pair<int, int>>> adjacency;  // row -> (meta-edge, other row)
};

} // namespace

bool assignDirections(const PolygonMesh& m, const EdgeLabels& labels, const DirectionOptions& o, DirectionProduct& out, std::string* err) {
    const auto t0 = std::chrono::steady_clock::now();
    out = DirectionProduct{};
    out.report.stage = "directions";
    out.report.algorithmVersion = 1;
    out.report.inputDigest = m.digest();
    FaceGraph g;
    if (!buildFaceGraph(m, labels, g, err)) return false;
    Rows r;
    if (!buildRows(g, r, err)) return false;
    const int nF = static_cast<int>(g.faces.size()), R = static_cast<int>(r.rows.size());
    out.rows = R;
    for (char c : r.closed) out.closedRows += c;
    out.openRows = R - out.closedRows;
    const std::vector<std::vector<char>> side = sidesOf(g, r);
    // the meta-graph: one node per row, an edge per pair of row halves that share course edges
    MetaGraph meta;
    meta.adjacency.resize(R);
    std::map<EdgeKey, int> weightOfPair;
    std::set<EdgeKey> seen;
    for (int f = 0; f < nF; ++f)
        for (size_t k = 0; k < g.faces[f].size(); ++k) {
            if (!g.course[f][k]) continue;
            const EdgeKey e = edgeKey(g.faces[f][k], g.faces[f][(k + 1) % g.faces[f].size()]);
            if (!seen.insert(e).second) continue;
            const std::pair<int, int> other = g.across(f, static_cast<int>(k));
            if (other.first < 0) {
                ++out.boundaryCourseEdges;
                continue;
            }
            const int h1 = 2 * r.rowOf[f] + side[f][k], h2 = 2 * r.rowOf[other.first] + side[other.first][other.second];
            if (r.rowOf[f] == r.rowOf[other.first]) {
                if (h1 != h2) ++out.sameRowMatched;
                else ++out.sameRowMismatched;
                continue;
            }
            ++weightOfPair[edgeKey(h1, h2)];
        }
    for (const auto& [pair, w] : weightOfPair) {
        meta.edges.push_back(pair);
        meta.weight.push_back(w);
        out.totalWeight += w;
    }
    out.metaEdges = static_cast<int>(meta.edges.size());
    for (size_t i = 0; i < meta.edges.size(); ++i) {
        meta.adjacency[meta.edges[i].first / 2].push_back({static_cast<int>(i), meta.edges[i].second / 2});
        meta.adjacency[meta.edges[i].second / 2].push_back({static_cast<int>(i), meta.edges[i].first / 2});
    }
    // components and the start: a breadth-first 2-colouring, matched wherever the tree reaches
    std::vector<char> y(R, 0);
    std::vector<int> component(R, -1);
    for (int r0 = 0; r0 < R; ++r0) {
        if (component[r0] >= 0) continue;
        component[r0] = out.components++;
        y[r0] = 0;
        std::queue<int> q;
        q.push(r0);
        while (!q.empty()) {
            const int row = q.front();
            q.pop();
            for (const std::pair<int, int>& next : meta.adjacency[row]) {
                const int i = next.first, s = next.second;
                if (component[s] >= 0) continue;
                component[s] = component[r0];
                const int hr = meta.edges[i].first / 2 == row ? meta.edges[i].first : meta.edges[i].second;
                const int hs = meta.edges[i].first / 2 == row ? meta.edges[i].second : meta.edges[i].first;
                const double roleR = slopeOf(hr) * y[row] + offsetOf(hr), roleS = 1 - roleR;
                y[s] = static_cast<char>(std::lround((roleS - offsetOf(hs)) / slopeOf(hs)));
                q.push(s);
            }
        }
    }
    if (!meta.edges.empty()) {
        const int nCol = R + out.metaEdges;
        HighsLp lp;
        lp.num_col_ = nCol;
        lp.sense_ = ObjSense::kMinimize;
        lp.col_cost_.assign(nCol, 0.0);
        for (int i = 0; i < out.metaEdges; ++i) lp.col_cost_[R + i] = -meta.weight[i];
        lp.col_lower_.assign(nCol, 0.0);
        lp.col_upper_.assign(nCol, 1.0);
        lp.integrality_.assign(nCol, HighsVarType::kContinuous);
        for (int row = 0; row < R; ++row) lp.integrality_[row] = HighsVarType::kInteger;
        std::vector<HighsInt> start{0}, index;
        std::vector<double> value;
        for (int i = 0; i < out.metaEdges; ++i) {
            const int h1 = meta.edges[i].first, h2 = meta.edges[i].second;
            // z <= l1 + l2 and z <= 2 - l1 - l2 with l = slope y + offset: z is 1 only where the roles differ
            index.insert(index.end(), {R + i, h1 / 2, h2 / 2});
            value.insert(value.end(), {1.0, -slopeOf(h1), -slopeOf(h2)});
            start.push_back(static_cast<HighsInt>(index.size()));
            lp.row_lower_.push_back(-kHighsInf);
            lp.row_upper_.push_back(offsetOf(h1) + offsetOf(h2));
            index.insert(index.end(), {R + i, h1 / 2, h2 / 2});
            value.insert(value.end(), {1.0, slopeOf(h1), slopeOf(h2)});
            start.push_back(static_cast<HighsInt>(index.size()));
            lp.row_lower_.push_back(-kHighsInf);
            lp.row_upper_.push_back(2 - offsetOf(h1) - offsetOf(h2));
        }
        lp.num_row_ = static_cast<HighsInt>(lp.row_lower_.size());
        lp.a_matrix_.format_ = MatrixFormat::kRowwise;
        lp.a_matrix_.num_col_ = nCol;
        lp.a_matrix_.num_row_ = lp.num_row_;
        lp.a_matrix_.start_ = start;
        lp.a_matrix_.index_ = index;
        lp.a_matrix_.value_ = value;
        Highs highs;
        highs.setOptionValue("output_flag", false);
        highs.setOptionValue("threads", 1);
        highs.setOptionValue("random_seed", 0);
        highs.setOptionValue("mip_rel_gap", 0.0);
        highs.setOptionValue("mip_abs_gap", 0.99);  // integral weights
        highs.setOptionValue("time_limit", o.timeLimitS);
        highs.setOptionValue("mip_max_nodes", o.maxNodes);
        const auto tm = std::chrono::steady_clock::now();
        if (highs.passModel(lp) == HighsStatus::kError) {
            if (err) *err = "directions: HiGHS refused the row program";
            return false;
        }
        HighsSolution sol;
        sol.col_value.assign(nCol, 0.0);
        for (int row = 0; row < R; ++row) sol.col_value[row] = y[row];
        for (int i = 0; i < out.metaEdges; ++i) {
            const int h1 = meta.edges[i].first, h2 = meta.edges[i].second;
            const double l1 = slopeOf(h1) * y[h1 / 2] + offsetOf(h1), l2 = slopeOf(h2) * y[h2 / 2] + offsetOf(h2);
            sol.col_value[R + i] = l1 != l2 ? 1.0 : 0.0;
        }
        sol.value_valid = true;
        highs.setSolution(sol);
        highs.run();
        out.msSolve = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - tm).count();
        const HighsModelStatus status = highs.getModelStatus();
        out.optimal = status == HighsModelStatus::kOptimal;
        out.solverStatus = out.optimal ? "optimal" : highs.modelStatusToString(status);
        if (highs.getInfo().primal_solution_status != kSolutionStatusFeasible) {
            if (err) *err = "directions: no row orientation found (" + out.solverStatus + ")";
            return false;
        }
        const std::vector<double>& col = highs.getSolution().col_value;
        for (int row = 0; row < R; ++row) y[row] = col[row] > 0.5 ? 1 : 0;
    } else {
        out.solverStatus = "no meta-edges";
        out.optimal = true;
    }
    // the roles, then each component's anchor: the mean wale direction along up, or the caller's reverse
    DirectedMesh& d = out.directed;
    d.mesh = &m;
    d.rowOfFace = r.rowOf;
    d.rows = r.rows;
    d.rowClosed = r.closed;
    d.componentOfRow = component;
    auto roleOf = [&](int f, int k) {
        const int yr = y[r.rowOf[f]];
        return static_cast<char>(side[f][k] == 0 ? yr : 1 - yr);
    };
    auto waleOf = [&](int f, Vec3& centre) {
        const int n = static_cast<int>(g.faces[f].size());
        Vec3 top, bottom;
        int nt = 0, nb = 0;
        centre = Vec3{};
        for (int k = 0; k < n; ++k) {
            centre += m.pos[g.faces[f][k]];
            if (!g.course[f][k]) continue;
            const Vec3 mid = (m.pos[g.faces[f][k]] + m.pos[g.faces[f][(k + 1) % n]]) * 0.5;
            if (roleOf(f, k)) {
                top += mid;
                ++nt;
            } else {
                bottom += mid;
                ++nb;
            }
        }
        centre = centre / n;
        if (!nt || !nb) return Vec3{};
        return top / nt - bottom / nb;
    };
    std::vector<double> score(out.components, 0.0);
    for (int f = 0; f < nF; ++f) {
        Vec3 c;
        const Vec3 w = waleOf(f, c);
        score[component[r.rowOf[f]]] += m.faceArea(f) * dot(w, o.up);
    }
    for (int row = 0; row < R; ++row) {
        const int c = component[row];
        bool flip = score[c] < 0;
        if (c < static_cast<int>(o.reverseComponent.size()) && o.reverseComponent[c]) flip = !flip;
        if (flip) y[row] ^= 1;
    }
    d.top.resize(nF);
    d.waleDirection.resize(nF);
    d.faceCentre.resize(nF);
    for (int f = 0; f < nF; ++f) {
        const int n = static_cast<int>(g.faces[f].size());
        d.top[f].assign(n, -1);
        for (int k = 0; k < n; ++k)
            if (g.course[f][k]) d.top[f][k] = roleOf(f, k);
        d.waleDirection[f] = waleOf(f, d.faceCentre[f]);
    }
    for (size_t i = 0; i < meta.edges.size(); ++i) {
        const int h1 = meta.edges[i].first, h2 = meta.edges[i].second;
        const double l1 = slopeOf(h1) * y[h1 / 2] + offsetOf(h1), l2 = slopeOf(h2) * y[h2 / 2] + offsetOf(h2);
        if (l1 != l2) out.matchedWeight += meta.weight[i];
    }
    seen.clear();
    for (int f = 0; f < nF; ++f)
        for (size_t k = 0; k < g.faces[f].size(); ++k) {
            if (!g.course[f][k]) continue;
            const int a = g.faces[f][k], b = g.faces[f][(k + 1) % g.faces[f].size()];
            if (!seen.insert(edgeKey(a, b)).second) continue;
            const std::pair<int, int> other = g.across(f, static_cast<int>(k));
            if (other.first < 0) continue;
            if (d.top[f][k] != d.top[other.first][other.second]) continue;
            const bool sameRow = r.rowOf[f] == r.rowOf[other.first];
            d.mismatches.push_back({a, b, f, other.first, d.top[f][k] == 1, sameRow});
            if (!sameRow) ++out.interRowMismatched;
        }
    out.report.outputDigest = m.digest();
    char settings[160];
    std::snprintf(settings, sizeof settings, "    up = (%.2f, %.2f, %.2f)\n    time limit = %.0f s\n    node limit = %d\n    reversed components = %zu\n", o.up.x, o.up.y,
                  o.up.z, o.timeLimitS, o.maxNodes, o.reverseComponent.size());
    out.report.settings = settings;
    out.report.counts.push_back({"rows", out.rows});
    out.report.counts.push_back({"closed rows", out.closedRows});
    out.report.counts.push_back({"meta-edges", out.metaEdges});
    out.report.counts.push_back({"components", out.components});
    out.report.counts.push_back({"matched course edges", out.matchedWeight});
    out.report.counts.push_back({"course edges between rows", out.totalWeight});
    out.report.counts.push_back({"same-row matched", out.sameRowMatched});
    out.report.counts.push_back({"same-row mismatched", out.sameRowMismatched});
    out.report.counts.push_back({"mismatched edges between rows", out.interRowMismatched});
    out.report.counts.push_back({"mismatched edges within a row", static_cast<double>(d.mismatches.size()) - out.interRowMismatched});
    out.report.counts.push_back({"boundary course edges", out.boundaryCourseEdges});
    out.report.counts.push_back({"ms solve", out.msSolve});
    if (!out.optimal) {
        out.report.status = Status::Unresolved;
        out.report.detail = "feasible, not proven optimal: " + out.solverStatus;
    } else if (!d.mismatches.empty()) {
        out.report.status = Status::Unresolved;
        out.report.detail = std::to_string(out.interRowMismatched) + " of " + std::to_string(static_cast<long long>(out.totalWeight)) +
                            " course edges between rows and " + std::to_string(static_cast<int>(d.mismatches.size()) - out.interRowMismatched) +
                            " within rows keep mismatched wale directions (listed; a later stage supports or names them)";
    } else {
        out.report.status = Status::Pass;
        out.report.detail = "every course edge between rows matched, program " + out.solverStatus;
    }
    out.report.elapsedMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    return true;
}

} // namespace knit
