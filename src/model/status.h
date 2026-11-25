#pragma once
// The verdict of a check. Unresolved never counts as a pass (numerical contract: knitting/checks/checks.h).
namespace knit {

enum class Status { Pass, Fail, Unresolved };

inline const char* toString(Status s) {
    switch (s) {
    case Status::Pass: return "PASS";
    case Status::Fail: return "FAIL";
    default: return "UNRESOLVED";
    }
}

// Fail dominates, then Unresolved.
inline Status combine(Status a, Status b) {
    if (a == Status::Fail || b == Status::Fail) return Status::Fail;
    if (a == Status::Unresolved || b == Status::Unresolved) return Status::Unresolved;
    return Status::Pass;
}

} // namespace knit
