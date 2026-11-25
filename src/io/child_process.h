#pragma once
// A program started in the background: blender.exe (io/blender_process.h) and the live preview's
// renderer, `ppt serve` (jobs/render_preview.h). Output goes to a log file; optionally the child reads
// lines from us on its standard input. The object owns the process: its destructor terminates a child
// that is still running and waits for it to exit. Windows only: elsewhere start() fails.
#include <string>
#include <vector>

namespace knit {

struct ChildProcessStart {
    std::string program;                 // the .exe
    std::vector<std::string> arguments;  // each one quoted on the command line
    std::string logPath;                 // stdout and stderr
    bool pipeStdin = false;              // writeLine() sends lines to the child's stdin
};

// This process's id (names its temporary folders and pipes).
unsigned long currentProcessId();

class ChildProcess {
public:
    ChildProcess() = default;
    ~ChildProcess();
    ChildProcess(const ChildProcess&) = delete;
    ChildProcess& operator=(const ChildProcess&) = delete;

    // false with *err when the program did not start.
    bool start(const ChildProcessStart& what, std::string* err);
    bool started() const { return process_ != nullptr; }
    // One line (a newline is added) to the child's stdin; false without a pipe or once the child is gone.
    bool writeLine(const std::string& line);
    // True when the child has exited (now or within `milliseconds`).
    bool exitedWithin(int milliseconds);
    int exitCode() const;
    // Ends the child at once and waits up to `waitMs` for the exit to be confirmed. True: it has exited.
    bool terminate(int waitMs);

private:
    void* process_ = nullptr;     // HANDLE
    void* stdinWrite_ = nullptr;  // HANDLE, our end of the child's stdin
};

} // namespace knit
