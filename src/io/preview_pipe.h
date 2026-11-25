#pragma once
// The live preview's frame pipe: the named pipe on which the renderer process (`ppt serve`) sends us its
// STATUS messages and rendered FRAMEs. We create the pipe and
// pass its name on the renderer's command line; the renderer connects as the client and writes.
// Each message: four little-endian uint32 {magic 'PPVF', type, payload bytes, sequence}, then the payload.
//   type 1 STATUS: a flat JSON object (strings and numbers only)
//   type 2 FRAME:  uint32 width, height, spp, camera_rev, scene_rev, then width*height*4 bytes RGBA8,
//                  display-encoded, rows top first
// Windows only: elsewhere create() fails.
#include <atomic>
#include <cstdint>
#include <string>
#include <vector>

namespace knit {

struct PreviewStatus {
    std::string state;     // starting, idle, loading, rendering, converged, suspended, released, refused, error
    std::string message;
    std::string ack;       // the command this status acknowledges ("" for a plain status)
    long long requestId = 0, freedBytes = 0, processBytes = 0, availableBytes = 0;
    int sceneRev = 0, cameraRev = 0, spp = 0;
    double msPerSpp = 0, loadMs = 0;
};

struct PreviewFrame {
    int width = 0, height = 0, spp = 0, cameraRev = 0, sceneRev = 0;
    std::vector<uint8_t> rgba;  // width * height * 4, rows top first
};

// A STATUS payload read into its fields; absent fields keep their defaults.
PreviewStatus parsePreviewStatus(const std::string& json);

class PreviewPipe {
public:
    enum class Message { Closed, Status, Frame };

    PreviewPipe() = default;
    ~PreviewPipe();
    PreviewPipe(const PreviewPipe&) = delete;
    PreviewPipe& operator=(const PreviewPipe&) = delete;

    // A new pipe with a name unique to this process; false with *err when Windows refused.
    bool create(std::string* err);
    const std::string& name() const { return name_; }
    // Blocks until the renderer connected (true) or unblock() was called / the pipe failed (false).
    bool waitForClient();
    // Blocks for the next message; Closed when the renderer is gone or unblock() was called.
    Message read(PreviewStatus& status, PreviewFrame& frame);
    // Lets a thread blocked in waitForClient() or read() return (from another thread).
    void unblock();

private:
    bool readExact(void* into, size_t bytes);
    void* pipe_ = nullptr;  // HANDLE
    std::string name_;
    std::atomic<bool> unblocked_{false};  // set by another thread
};

} // namespace knit
