#pragma once
// Crash report: what the program could not catch (an access violation, std::terminate) is
// written down before the process ends, so a crash can be traced to its code.
#include <string>

namespace knit {

// Installs the report for the whole process (every thread): the exception, the faulting module
// and offset, and the call stack with function names where the executable has symbols (a PDB
// next to it). Written to stderr and appended to `logPath` ("": stderr only). The process then
// ends as it would have without the report.
void installCrashReport(const std::string& logPath);

} // namespace knit
