// Names, shapes and lookups of surface markings (surface_marking.h).
#include "model/surface_marking.h"

namespace knit {

const char* toString(MarkingMeaning m) {
    switch (m) {
    case MarkingMeaning::Start: return "start";
    case MarkingMeaning::Finish: return "finish";
    case MarkingMeaning::RowStartPath: return "row_start";
    case MarkingMeaning::Boundary: return "boundary";
    case MarkingMeaning::Cut: return "cut";
    case MarkingMeaning::Attachment: return "attachment";
    }
    return "start";
}

const char* toString(MarkingShape s) {
    switch (s) {
    case MarkingShape::Point: return "point";
    case MarkingShape::Path: return "path";
    case MarkingShape::ClosedCurve: return "closed";
    }
    return "point";
}

const char* toString(MarkingOrigin o) { return o == MarkingOrigin::Proposal ? "proposal" : "user"; }

const char* toString(ProposalState p) {
    switch (p) {
    case ProposalState::None: return "none";
    case ProposalState::Offered: return "offered";
    case ProposalState::Accepted: return "accepted";
    case ProposalState::Rejected: return "rejected";
    }
    return "none";
}

bool parseMeaning(const std::string& text, MarkingMeaning& out) {
    const MarkingMeaning all[] = {MarkingMeaning::Start,    MarkingMeaning::Finish, MarkingMeaning::RowStartPath,
                                  MarkingMeaning::Boundary, MarkingMeaning::Cut,    MarkingMeaning::Attachment};
    for (MarkingMeaning m : all)
        if (text == toString(m)) {
            out = m;
            return true;
        }
    return false;
}

bool parseShape(const std::string& text, MarkingShape& out) {
    const MarkingShape all[] = {MarkingShape::Point, MarkingShape::Path, MarkingShape::ClosedCurve};
    for (MarkingShape s : all)
        if (text == toString(s)) {
            out = s;
            return true;
        }
    return false;
}

bool parseOrigin(const std::string& text, MarkingOrigin& out) {
    if (text == "user") out = MarkingOrigin::User;
    else if (text == "proposal") out = MarkingOrigin::Proposal;
    else return false;
    return true;
}

bool parseProposalState(const std::string& text, ProposalState& out) {
    const ProposalState all[] = {ProposalState::None, ProposalState::Offered, ProposalState::Accepted, ProposalState::Rejected};
    for (ProposalState p : all)
        if (text == toString(p)) {
            out = p;
            return true;
        }
    return false;
}

const char* meaningLabel(MarkingMeaning m) {
    switch (m) {
    case MarkingMeaning::Start: return "Start";
    case MarkingMeaning::Finish: return "Finish";
    case MarkingMeaning::RowStartPath: return "Row-start path";
    case MarkingMeaning::Boundary: return "Boundary";
    case MarkingMeaning::Cut: return "Cut";
    case MarkingMeaning::Attachment: return "Attachment";
    }
    return "Start";
}

MarkingShape defaultShape(MarkingMeaning m) {
    switch (m) {
    case MarkingMeaning::RowStartPath:
    case MarkingMeaning::Cut: return MarkingShape::Path;
    case MarkingMeaning::Start:
    case MarkingMeaning::Finish:
    case MarkingMeaning::Boundary:
    case MarkingMeaning::Attachment: return MarkingShape::ClosedCurve;
    }
    return MarkingShape::ClosedCurve;
}

bool shapeAllowed(MarkingMeaning m, MarkingShape s) {
    switch (m) {
    case MarkingMeaning::Start: return s != MarkingShape::Path;
    case MarkingMeaning::RowStartPath: return s == MarkingShape::Path;
    case MarkingMeaning::Cut: return s != MarkingShape::Point;
    case MarkingMeaning::Finish: return s != MarkingShape::Path;
    case MarkingMeaning::Boundary:
    case MarkingMeaning::Attachment: return s == MarkingShape::ClosedCurve;
    }
    return false;
}

bool markingComplete(const SurfaceMarking& m) {
    const size_t need = m.shape == MarkingShape::Point ? 1 : m.shape == MarkingShape::Path ? 2 : 3;
    return m.points.size() >= need;
}

SurfaceMarking* findMarking(MarkingSet& s, int id) {
    for (SurfaceMarking& m : s.markings)
        if (m.id == id) return &m;
    return nullptr;
}

const SurfaceMarking* findMarking(const MarkingSet& s, int id) {
    for (const SurfaceMarking& m : s.markings)
        if (m.id == id) return &m;
    return nullptr;
}

std::string markingSetting(const SurfaceMarking& m, const std::string& key) {
    for (const auto& [k, v] : m.settings)
        if (k == key) return v;
    return "";
}

MarkingRole markingRole(MarkingMeaning m) {
    switch (m) {
    case MarkingMeaning::Start: return MarkingRole::Start;
    case MarkingMeaning::Finish:
    case MarkingMeaning::Boundary: return MarkingRole::Termination;
    case MarkingMeaning::RowStartPath: return MarkingRole::RowStart;
    case MarkingMeaning::Cut:
    case MarkingMeaning::Attachment: return MarkingRole::None;
    }
    return MarkingRole::None;
}

std::vector<int> setupMarkings(const MarkingSet& s, int setup) {
    std::vector<int> ids;
    if (setup <= 0) return ids;
    for (const SurfaceMarking& m : s.markings)
        if (m.setup == setup) ids.push_back(m.id);
    return ids;
}

void rejectSetup(MarkingSet& s, int setup) {
    if (setup <= 0) return;
    for (SurfaceMarking& m : s.markings)
        if (m.setup == setup && m.origin == MarkingOrigin::Proposal && m.proposal == ProposalState::Offered) m.proposal = ProposalState::Rejected;
}

} // namespace knit
