// The live preview's frame pipe (preview_pipe.h).
#include "io/preview_pipe.h"

#include <atomic>
#include <cstdlib>
#include <cstring>

#ifdef _WIN32
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

namespace knit {

namespace {

constexpr uint32_t kMagic = 0x46565050u;  // "PPVF"
constexpr uint32_t kStatus = 1, kFrame = 2;
constexpr uint32_t kMaxPayload = 64u << 20;  // a 4K frame is 33 MB

// The text after `"key":` in a flat JSON object, "" when absent.
std::string valueText(const std::string& json, const char* key) {
    const std::string k = std::string("\"") + key + "\":";
    const size_t at = json.find(k);
    if (at == std::string::npos) return "";
    size_t i = at + k.size();
    while (i < json.size() && json[i] == ' ') ++i;
    if (i < json.size() && json[i] == '"') {  // a string: up to the closing quote, escapes kept simple
        std::string out;
        for (++i; i < json.size() && json[i] != '"'; ++i) {
            if (json[i] == '\\' && i + 1 < json.size()) {
                ++i;
                out += json[i] == 'n' ? '\n' : json[i];
            } else {
                out += json[i];
            }
        }
        return out;
    }
    const size_t end = json.find_first_of(",}", i);
    return json.substr(i, end == std::string::npos ? std::string::npos : end - i);
}

long long integer(const std::string& json, const char* key) { return std::atoll(valueText(json, key).c_str()); }
double number(const std::string& json, const char* key) { return std::atof(valueText(json, key).c_str()); }

} // namespace

PreviewStatus parsePreviewStatus(const std::string& json) {
    PreviewStatus s;
    s.state = valueText(json, "state");
    s.message = valueText(json, "message");
    s.ack = valueText(json, "ack");
    s.requestId = integer(json, "request_id");
    s.freedBytes = integer(json, "freed_bytes");
    s.processBytes = integer(json, "process_bytes");
    s.availableBytes = integer(json, "available_bytes");
    s.sceneRev = static_cast<int>(integer(json, "scene_rev"));
    s.cameraRev = static_cast<int>(integer(json, "camera_rev"));
    s.spp = static_cast<int>(integer(json, "spp"));
    s.msPerSpp = number(json, "ms_per_spp");
    s.loadMs = number(json, "load_ms");
    return s;
}

#ifdef _WIN32

PreviewPipe::~PreviewPipe() {
    if (pipe_) CloseHandle(static_cast<HANDLE>(pipe_));
}

bool PreviewPipe::create(std::string* err) {
    static std::atomic<int> counter{0};
    name_ = "\\\\.\\pipe\\knit-preview-" + std::to_string(GetCurrentProcessId()) + "-" + std::to_string(++counter);
    const std::wstring wide(name_.begin(), name_.end());
    HANDLE h = CreateNamedPipeW(wide.c_str(), PIPE_ACCESS_INBOUND, PIPE_TYPE_BYTE | PIPE_WAIT, 1, 1 << 20, 1 << 20, 0, nullptr);
    if (h == INVALID_HANDLE_VALUE) {
        if (err) *err = "cannot create the preview pipe (Windows error " + std::to_string(GetLastError()) + ")";
        return false;
    }
    pipe_ = h;
    return true;
}

bool PreviewPipe::waitForClient() {
    if (!pipe_) return false;
    const BOOL ok = ConnectNamedPipe(static_cast<HANDLE>(pipe_), nullptr);
    if (unblocked_) return false;
    return ok || GetLastError() == ERROR_PIPE_CONNECTED;
}

bool PreviewPipe::readExact(void* into, size_t bytes) {
    char* p = static_cast<char*>(into);
    while (bytes > 0) {
        DWORD got = 0;
        if (!ReadFile(static_cast<HANDLE>(pipe_), p, static_cast<DWORD>(bytes), &got, nullptr) || got == 0) return false;
        p += got;
        bytes -= got;
    }
    return !unblocked_;
}

PreviewPipe::Message PreviewPipe::read(PreviewStatus& status, PreviewFrame& frame) {
    for (;;) {
        uint32_t header[4];
        if (!pipe_ || !readExact(header, sizeof header) || header[0] != kMagic || header[2] > kMaxPayload) return Message::Closed;
        std::vector<char> payload(header[2]);
        if (header[2] && !readExact(payload.data(), payload.size())) return Message::Closed;
        if (header[1] == kStatus) {
            status = parsePreviewStatus(std::string(payload.begin(), payload.end()));
            return Message::Status;
        }
        if (header[1] == kFrame && payload.size() >= 20) {
            uint32_t head[5];
            std::memcpy(head, payload.data(), sizeof head);
            const size_t pixels = static_cast<size_t>(head[0]) * head[1] * 4;
            if (payload.size() < 20 + pixels) return Message::Closed;
            frame.width = static_cast<int>(head[0]);
            frame.height = static_cast<int>(head[1]);
            frame.spp = static_cast<int>(head[2]);
            frame.cameraRev = static_cast<int>(head[3]);
            frame.sceneRev = static_cast<int>(head[4]);
            frame.rgba.assign(payload.begin() + 20, payload.begin() + 20 + static_cast<long long>(pixels));
            return Message::Frame;
        }
        // An unknown type: skipped (a newer renderer may send more).
    }
}

void PreviewPipe::unblock() {
    unblocked_ = true;
    // A waiting ConnectNamedPipe returns when a client connects: connect, then close at once. A blocked
    // ReadFile returns when the writer is gone (the renderer ended) or the pipe is cancelled.
    const std::wstring wide(name_.begin(), name_.end());
    HANDLE client = CreateFileW(wide.c_str(), GENERIC_WRITE, 0, nullptr, OPEN_EXISTING, 0, nullptr);
    if (client != INVALID_HANDLE_VALUE) CloseHandle(client);
    if (pipe_) CancelIoEx(static_cast<HANDLE>(pipe_), nullptr);
}

#else

PreviewPipe::~PreviewPipe() = default;
bool PreviewPipe::create(std::string* err) {
    if (err) *err = "the live preview runs on Windows only";
    return false;
}
bool PreviewPipe::waitForClient() { return false; }
bool PreviewPipe::readExact(void*, size_t) { return false; }
PreviewPipe::Message PreviewPipe::read(PreviewStatus&, PreviewFrame&) { return Message::Closed; }
void PreviewPipe::unblock() {}

#endif

} // namespace knit
