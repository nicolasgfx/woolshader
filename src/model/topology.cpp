// Knit topology: which stitch is drawn through which (topology.h).
#include "model/topology.h"

#include <algorithm>
#include <limits>

namespace knit {

const char* toString(StitchOp op) {
    switch (op) {
    case StitchOp::K: return "K";
    case StitchOp::Inc: return "INC";
    case StitchOp::Dec: return "DEC";
    }
    return "?";
}

const char* toString(StitchFace f) { return f == StitchFace::Purl ? "purl" : "knit"; }

void linkChildren(std::vector<StitchRecord>& stitches) {
    for (StitchRecord& s : stitches) s.nChildren = 0;
    for (const StitchRecord& c : stitches)
        for (int k = 0; k < c.nParents; ++k) {
            StitchRecord& p = stitches[c.parents[k]];
            if (p.nChildren < 2) p.children[p.nChildren++] = c.id;
        }
}

void KnitTopology::buildIndex(const std::vector<double>& yarnEnd) {
    linkChildren(stitches);
    byYarn_.assign(yarnEnd.size(), {});
    yarnEnd_ = yarnEnd;
    for (const StitchRecord& s : stitches) {
        if (s.yarn >= static_cast<int>(byYarn_.size())) {
            byYarn_.resize(s.yarn + 1);
            yarnEnd_.resize(s.yarn + 1, 1e300);
        }
        byYarn_[s.yarn].push_back({s.uBegin, s.id});
    }
    for (auto& v : byYarn_) std::sort(v.begin(), v.end());
}

int KnitTopology::stitchContaining(int yarn, double u0, double u1) const {
    if (yarn < 0 || yarn >= static_cast<int>(byYarn_.size())) return -1;
    const std::vector<std::pair<double, int>>& v = byYarn_[yarn];
    auto it = std::upper_bound(v.begin(), v.end(), std::make_pair(u0, 1 << 30));
    if (it == v.begin()) return -1;
    const StitchRecord& s = stitches[std::prev(it)->second];
    const bool last = s.uEnd >= yarnEnd_[yarn];
    if (u0 >= s.uBegin && (u1 < s.uEnd || (last && u1 <= s.uEnd))) return s.id;
    return -1;
}

int KnitTopology::stitchAt(int yarn, double u) const { return stitchContaining(yarn, u, u); }

int KnitTopology::stitchNear(int yarn, double u) const {
    if (yarn < 0 || yarn >= static_cast<int>(byYarn_.size()) || byYarn_[yarn].empty()) return -1;
    const std::vector<std::pair<double, int>>& v = byYarn_[yarn];
    auto it = std::upper_bound(v.begin(), v.end(), std::make_pair(u, std::numeric_limits<int>::max()));
    return it == v.begin() ? v.front().second : std::prev(it)->second;
}

bool KnitTopology::isParent(int parent, int child) const {
    if (child < 0 || parent < 0) return false;
    const StitchRecord& c = stitches[child];
    for (int k = 0; k < c.nParents; ++k)
        if (c.parents[k] == parent) return true;
    return false;
}

int KnitTopology::stitchOrLoopAt(int yarn, double u) const {
    int t = stitchAt(yarn, u);
    if (t >= 0) return t;
    for (const TerminalRecord& tr : terminals)
        if (tr.yarn == yarn)
            for (size_t k = 0; k < tr.segments.size() && k < tr.stitches.size(); ++k)
                if (u >= tr.segments[k].first && u <= tr.segments[k].second) return tr.stitches[k];
    return -1;
}

} // namespace knit
