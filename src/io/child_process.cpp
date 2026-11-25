// A program started in the background (child_process.h): CreateProcessW with the output in a log file
// and, optionally, a pipe to its standard input.
#include "io/child_process.h"
#include "io/files.h"

#include <filesystem>

#ifdef _WIN32
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

namespace knit {

namespace fs = std::filesystem;

#ifdef _WIN32

unsigned long currentProcessId() { return GetCurrentProcessId(); }

namespace {

std::wstring quoted(const std::string& argument) { return L"\"" + fs::path(argument).wstring() + L"\""; }

} // namespace

ChildProcess::~ChildProcess() {
    if (stdinWrite_) CloseHandle(static_cast<HANDLE>(stdinWrite_));
    if (!process_) return;
    if (WaitForSingleObject(static_cast<HANDLE>(process_), 0) == WAIT_TIMEOUT) {
        TerminateProcess(static_cast<HANDLE>(process_), 1);
        WaitForSingleObject(static_cast<HANDLE>(process_), 5000);
    }
    CloseHandle(static_cast<HANDLE>(process_));
}

bool ChildProcess::start(const ChildProcessStart& what, std::string* err) {
    SECURITY_ATTRIBUTES inherit{};
    inherit.nLength = sizeof inherit;
    inherit.bInheritHandle = TRUE;
    HANDLE log = CreateFileW(fs::path(what.logPath).c_str(), GENERIC_WRITE, FILE_SHARE_READ, &inherit, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr);
    if (log == INVALID_HANDLE_VALUE) return fail(err, "cannot write " + what.logPath);
    HANDLE childStdin = nullptr, ourStdin = nullptr;
    if (what.pipeStdin) {
        if (!CreatePipe(&childStdin, &ourStdin, &inherit, 0)) {
            CloseHandle(log);
            return fail(err, "cannot create a pipe for " + what.program);
        }
        SetHandleInformation(ourStdin, HANDLE_FLAG_INHERIT, 0);  // only the child's end is inherited
    }
    STARTUPINFOW startup{};
    startup.cb = sizeof startup;
    startup.dwFlags = STARTF_USESTDHANDLES;
    startup.hStdInput = childStdin;
    startup.hStdOutput = log;
    startup.hStdError = log;
    std::wstring commandLine = quoted(what.program);
    for (const std::string& a : what.arguments) commandLine += L" " + quoted(a);
    PROCESS_INFORMATION info{};
    const BOOL ok = CreateProcessW(nullptr, commandLine.data(), nullptr, nullptr, TRUE, CREATE_NO_WINDOW, nullptr, nullptr, &startup, &info);
    const DWORD reason = GetLastError();
    CloseHandle(log);
    if (childStdin) CloseHandle(childStdin);
    if (!ok) {
        if (ourStdin) CloseHandle(ourStdin);
        return fail(err, "cannot start " + fs::path(what.program).filename().string() + " (Windows error " + std::to_string(reason) + ")");
    }
    CloseHandle(info.hThread);
    process_ = info.hProcess;
    stdinWrite_ = ourStdin;
    return true;
}

bool ChildProcess::writeLine(const std::string& line) {
    if (!stdinWrite_) return false;
    const std::string text = line + "\n";
    DWORD wrote = 0;
    return WriteFile(static_cast<HANDLE>(stdinWrite_), text.data(), static_cast<DWORD>(text.size()), &wrote, nullptr) && wrote == text.size();
}

bool ChildProcess::exitedWithin(int milliseconds) {
    return process_ && WaitForSingleObject(static_cast<HANDLE>(process_), static_cast<DWORD>(milliseconds)) == WAIT_OBJECT_0;
}

int ChildProcess::exitCode() const {
    DWORD code = 1;
    if (process_) GetExitCodeProcess(static_cast<HANDLE>(process_), &code);
    return static_cast<int>(code);
}

bool ChildProcess::terminate(int waitMs) {
    if (!process_) return true;
    if (WaitForSingleObject(static_cast<HANDLE>(process_), 0) == WAIT_TIMEOUT) TerminateProcess(static_cast<HANDLE>(process_), 1);
    return exitedWithin(waitMs);
}

#else

unsigned long currentProcessId() { return 0; }
ChildProcess::~ChildProcess() = default;
bool ChildProcess::start(const ChildProcessStart&, std::string* err) {
    if (err) *err = "background programs run on Windows only";
    return false;
}
bool ChildProcess::writeLine(const std::string&) { return false; }
bool ChildProcess::exitedWithin(int) { return true; }
int ChildProcess::exitCode() const { return 1; }
bool ChildProcess::terminate(int) { return true; }

#endif

} // namespace knit
