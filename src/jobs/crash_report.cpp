// Writes down what the program could not catch, before the process ends.
#include "jobs/crash_report.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fstream>

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
// dbghelp.h after windows.h
#include <dbghelp.h>
#pragma comment(lib, "dbghelp.lib")
#endif

namespace knit {

namespace {

std::string g_logPath;

// One line to stderr and to the log file.
void say(const char* line) {
    std::fputs(line, stderr);
    if (g_logPath.empty()) return;
    std::ofstream log(g_logPath, std::ios::app);
    log << line;
}

#ifdef _WIN32

// "module+0xoffset function file:line" of one code address.
void describeAddress(HANDLE process, DWORD64 address, char* out, size_t outSize) {
    char moduleName[MAX_PATH] = "?";
    DWORD64 moduleBase = SymGetModuleBase64(process, address);
    if (moduleBase != 0) {
        char path[MAX_PATH];
        if (GetModuleFileNameA(reinterpret_cast<HMODULE>(moduleBase), path, MAX_PATH) > 0) {
            const char* slash = std::strrchr(path, '\\');
            std::snprintf(moduleName, sizeof moduleName, "%s", slash ? slash + 1 : path);
        }
    }
    char symbolBuffer[sizeof(SYMBOL_INFO) + 256] = {};
    SYMBOL_INFO* symbol = reinterpret_cast<SYMBOL_INFO*>(symbolBuffer);
    symbol->SizeOfStruct = sizeof(SYMBOL_INFO);
    symbol->MaxNameLen = 255;
    DWORD64 displacement = 0;
    const bool named = SymFromAddr(process, address, &displacement, symbol) != FALSE;
    IMAGEHLP_LINE64 line = {};
    line.SizeOfStruct = sizeof line;
    DWORD lineDisplacement = 0;
    const bool lined = SymGetLineFromAddr64(process, address, &lineDisplacement, &line) != FALSE;
    std::snprintf(out, outSize, "%s+0x%llx %s%s%s:%lu", moduleName, static_cast<unsigned long long>(address - moduleBase),
                  named ? symbol->Name : "", lined ? " " : "", lined ? line.FileName : "", lined ? line.LineNumber : 0ul);
}

LONG WINAPI onUnhandledException(EXCEPTION_POINTERS* info) {
    HANDLE process = GetCurrentProcess();
    HANDLE thread = GetCurrentThread();
    char line[1024];
    const EXCEPTION_RECORD* rec = info->ExceptionRecord;
    std::snprintf(line, sizeof line, "\n*** crash: exception 0x%08lX at %p (thread %lu)\n", rec->ExceptionCode, rec->ExceptionAddress,
                  GetCurrentThreadId());
    say(line);
    if (rec->ExceptionCode == EXCEPTION_ACCESS_VIOLATION && rec->NumberParameters >= 2) {
        std::snprintf(line, sizeof line, "    %s address 0x%llx\n", rec->ExceptionInformation[0] == 0 ? "reading" : "writing",
                      static_cast<unsigned long long>(rec->ExceptionInformation[1]));
        say(line);
    }
    // The call stack, walked from the faulting context (a copy: StackWalk64 changes it).
    CONTEXT context = *info->ContextRecord;
    STACKFRAME64 frame = {};
    frame.AddrPC.Offset = context.Rip;
    frame.AddrPC.Mode = AddrModeFlat;
    frame.AddrFrame.Offset = context.Rbp;
    frame.AddrFrame.Mode = AddrModeFlat;
    frame.AddrStack.Offset = context.Rsp;
    frame.AddrStack.Mode = AddrModeFlat;
    for (int depth = 0; depth < 48; ++depth) {
        if (!StackWalk64(IMAGE_FILE_MACHINE_AMD64, process, thread, &frame, &context, nullptr, SymFunctionTableAccess64,
                         SymGetModuleBase64, nullptr) ||
            frame.AddrPC.Offset == 0)
            break;
        char where[900];
        describeAddress(process, frame.AddrPC.Offset, where, sizeof where);
        std::snprintf(line, sizeof line, "    #%02d %s\n", depth, where);
        say(line);
    }
    say("*** end of crash report\n");
    std::fflush(stderr);
    return EXCEPTION_CONTINUE_SEARCH;  // the process ends as it would have
}

#endif

void onTerminate() {
    std::string what = "no active exception";
    if (std::exception_ptr e = std::current_exception()) {
        try {
            std::rethrow_exception(e);
        } catch (const std::exception& ex) {
            what = std::string("std::exception: ") + ex.what();
        } catch (...) {
            what = "an exception not derived from std::exception";
        }
    }
    say(("\n*** crash: std::terminate (" + what + ")\n").c_str());
    std::abort();
}

} // namespace

void installCrashReport(const std::string& logPath) {
    g_logPath = logPath;
    std::set_terminate(onTerminate);
#ifdef _WIN32
    SymSetOptions(SYMOPT_UNDNAME | SYMOPT_LOAD_LINES | SYMOPT_DEFERRED_LOADS);
    SymInitialize(GetCurrentProcess(), nullptr, TRUE);
    SetUnhandledExceptionFilter(onUnhandledException);
#endif
}

} // namespace knit
