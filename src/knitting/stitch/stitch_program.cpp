// The spiral stitch program of one knitted piece: turns of K / INC / DEC and its check (stitch_program.h).
#include "knitting/stitch/stitch_program.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <set>

namespace knit {

std::vector<int> scheduleCounts(int n0, const std::vector<double>& targets, double maxShaping) {
    std::vector<int> n;
    if (targets.empty()) return n;
    n.push_back(n0);
    for (size_t r = 1; r < targets.size(); ++r) {
        const int I = n.back();
        const int lim = static_cast<int>(std::floor(I * maxShaping));
        int want = static_cast<int>(std::lround(targets[r]));
        n.push_back(std::clamp(want, I - lim, I + lim));
    }
    return n;
}

StitchProgram buildProgram(const std::vector<int>& counts, bool castRing, bool closeDrawstring) {
    StitchProgram p;
    p.castRing = castRing;
    p.closeDrawstring = closeDrawstring;
    const double phi = 0.5 * (std::sqrt(5.0) - 1);
    std::vector<int> prev;  // output loops of the previous turn (stitch ids, in order)
    for (size_t r = 0; r < counts.size(); ++r) {
        TurnPlan tp;
        p.turnStart.push_back(static_cast<int>(p.stitches.size()));
        std::vector<int> out;
        auto emit = [&](StitchOp op, int pa, int pb, int slot) {
            StitchRecord s;
            s.id = static_cast<int>(p.stitches.size());
            s.turn = static_cast<int>(r);
            s.index = static_cast<int>(out.size());
            s.op = op;
            s.incSlot = slot;
            if (pa >= 0) s.parents[s.nParents++] = pa;
            if (pb >= 0) s.parents[s.nParents++] = pb;
            p.stitches.push_back(s);
            out.push_back(s.id);
        };
        if (r == 0) {
            for (int i = 0; i < counts[0]; ++i) emit(StitchOp::K, -1, -1, 0);
            tp.k = counts[0];
        } else {
            const int I = static_cast<int>(prev.size()), n = counts[r];
            tp.input = I;
            const int inc = std::max(0, n - I), dec = std::max(0, I - n);
            const double ph = std::fmod(phi * r, 1.0);
            std::vector<char> incAt(I, 0), decAt(I, 0);
            if (inc > 0) {
                const double step = double(I) / inc;  // >= 2 under the density limit
                for (int m = 0; m < inc; ++m) incAt[std::min(I - 1, static_cast<int>(std::floor(m * step + ph * (step - 1))))] = 1;
            }
            if (dec > 0) {
                const double step = double(I) / dec;  // >= 2
                for (int m = 0; m < dec; ++m)
                    decAt[std::min(I - 2, static_cast<int>(std::floor(m * step + ph * (step - 2))))] = 1;
            }
            for (int i = 0; i < I;) {
                if (decAt[i] && i + 1 < I) {
                    emit(StitchOp::Dec, prev[i], prev[i + 1], 0);
                    ++tp.dec;
                    i += 2;
                } else if (incAt[i]) {
                    emit(StitchOp::Inc, prev[i], -1, 0);
                    emit(StitchOp::Inc, prev[i], -1, 1);
                    ++tp.inc;
                    i += 1;
                } else {
                    emit(StitchOp::K, prev[i], -1, 0);
                    ++tp.k;
                    i += 1;
                }
            }
        }
        tp.output = static_cast<int>(out.size());
        p.turns.push_back(tp);
        prev = out;
    }
    linkChildren(p.stitches);
    return p;
}

GraphCheck checkProgram(const StitchProgram& p) {
    GraphCheck g;
    g.turns = static_cast<int>(p.turns.size());
    g.stitches = static_cast<int>(p.stitches.size());
    auto fail = [&](const std::string& what) {
        if (g.status == Status::Pass) g.detail = what;
        g.status = Status::Fail;
    };
    char buf[256];
    for (size_t k = 0; k < p.stitches.size(); ++k)
        if (p.stitches[k].id != static_cast<int>(k)) fail("stitch ids are not in yarn order");
    for (size_t r = 0; r < p.turns.size(); ++r) {
        const TurnPlan& tp = p.turns[r];
        const int s0 = p.turnStart[r];
        const int s1 = r + 1 < p.turns.size() ? p.turnStart[r + 1] : static_cast<int>(p.stitches.size());
        g.inc += tp.inc;
        g.dec += tp.dec;
        if (s1 - s0 != tp.output) { std::snprintf(buf, sizeof buf, "turn %zu: %d stitches, plan says %d", r, s1 - s0, tp.output); fail(buf); }
        if (r == 0) {
            for (int s = s0; s < s1; ++s)
                if (p.stitches[s].nParents != 0) fail("turn 0 stitch with a parent");
            continue;
        }
        if (tp.input != tp.k + tp.inc + 2 * tp.dec || tp.output != tp.k + 2 * tp.inc + tp.dec) {
            std::snprintf(buf, sizeof buf, "turn %zu: balance violated (in %d, K %d, INC %d, DEC %d, out %d)", r, tp.input,
                          tp.k, tp.inc, tp.dec, tp.output);
            fail(buf);
        }
        // Loops of the previous turn, in order, must be consumed exactly once, in order.
        const int p0 = p.turnStart[r - 1], p1 = s0;
        if (p1 - p0 != tp.input) { std::snprintf(buf, sizeof buf, "turn %zu: input %d != previous output %d", r, tp.input, p1 - p0); fail(buf); }
        int expect = p0;
        for (int s = s0; s < s1; ++s) {
            const StitchRecord& st = p.stitches[s];
            if (st.turn != static_cast<int>(r)) fail("stitch turn mismatch");
            if (st.op == StitchOp::Inc) {
                // The atomic pair: slot 0 then slot 1, same parent.
                if (st.incSlot == 0) {
                    if (s + 1 >= s1 || p.stitches[s + 1].op != StitchOp::Inc || p.stitches[s + 1].incSlot != 1 ||
                        p.stitches[s + 1].parents[0] != st.parents[0] || st.nParents != 1) {
                        std::snprintf(buf, sizeof buf, "turn %zu: INC at stitch %d is not an adjacent pair", r, s);
                        fail(buf);
                    }
                    if (st.parents[0] != expect) { std::snprintf(buf, sizeof buf, "turn %zu: stitch %d consumes %d, expected %d", r, s, st.parents[0], expect); fail(buf); }
                    ++expect;
                }
                continue;
            }
            const int want = st.op == StitchOp::Dec ? 2 : 1;
            if (st.nParents != want) { std::snprintf(buf, sizeof buf, "turn %zu: stitch %d has %d parents", r, s, st.nParents); fail(buf); }
            for (int q = 0; q < st.nParents; ++q) {
                if (st.parents[q] != expect) {
                    std::snprintf(buf, sizeof buf, "turn %zu: stitch %d consumes %d, expected %d%s", r, s, st.parents[q], expect,
                                  st.parents[q] >= p1 ? " (mixes generations)" : "");
                    fail(buf);
                }
                if (st.parents[q] < p0 || st.parents[q] >= p1) fail("parent outside the previous turn");
                ++expect;
            }
        }
        if (expect != p1) { std::snprintf(buf, sizeof buf, "turn %zu: %d of %d loops consumed", r, expect - p0, p1 - p0); fail(buf); }
    }
    if (g.status == Status::Pass) {
        std::snprintf(buf, sizeof buf, "%d turns, %d stitches, %d INC, %d DEC; balance and consumption hold", g.turns,
                      g.stitches, g.inc, g.dec);
        g.detail = buf;
    }
    return g;
}

} // namespace knit
