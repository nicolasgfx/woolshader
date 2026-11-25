// Native file and folder dialogs and "open in the file manager" (native_dialogs.h).
#include "viewer/native_dialogs.h"

#if defined(_WIN32)
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <shellapi.h>
#include <shobjidl.h>
#define GLFW_EXPOSE_NATIVE_WIN32
#include <GLFW/glfw3.h>
#include <GLFW/glfw3native.h>
#endif

#include <algorithm>
#include <filesystem>

namespace viewer {

#if defined(_WIN32)

namespace {

std::wstring widen(const std::string& s) {
    if (s.empty()) return {};
    const int n = MultiByteToWideChar(CP_UTF8, 0, s.data(), static_cast<int>(s.size()), nullptr, 0);
    std::wstring w(static_cast<size_t>(n), L'\0');
    MultiByteToWideChar(CP_UTF8, 0, s.data(), static_cast<int>(s.size()), w.data(), n);
    return w;
}

std::string narrow(const wchar_t* w) {
    const int n = WideCharToMultiByte(CP_UTF8, 0, w, -1, nullptr, 0, nullptr, nullptr);
    std::string s(static_cast<size_t>(std::max(0, n - 1)), '\0');
    WideCharToMultiByte(CP_UTF8, 0, w, -1, s.data(), n, nullptr, nullptr);
    return s;
}

enum class DialogKind { Mesh, Image, Markings, Folder };

std::optional<std::string> runDialog(GLFWwindow* owner, const std::string& start, DialogKind kind) {
    const HRESULT init = CoInitializeEx(nullptr, COINIT_APARTMENTTHREADED | COINIT_DISABLE_OLE1DDE);
    std::optional<std::string> chosen;
    IFileOpenDialog* dlg = nullptr;
    if (SUCCEEDED(CoCreateInstance(CLSID_FileOpenDialog, nullptr, CLSCTX_INPROC_SERVER, IID_PPV_ARGS(&dlg)))) {
        DWORD flags = 0;
        dlg->GetOptions(&flags);
        dlg->SetOptions(flags | FOS_FORCEFILESYSTEM | FOS_PATHMUSTEXIST | (kind == DialogKind::Folder ? FOS_PICKFOLDERS : FOS_FILEMUSTEXIST));
        if (kind == DialogKind::Image) {
            const COMDLG_FILTERSPEC filters[] = {{L"Images (PNG, JPEG, TGA, BMP)", L"*.png;*.jpg;*.jpeg;*.tga;*.bmp"}, {L"All files", L"*.*"}};
            dlg->SetFileTypes(2, filters);
            dlg->SetTitle(L"Open a fur coverage mask (white = full density)");
        } else if (kind == DialogKind::Markings) {
            const COMDLG_FILTERSPEC filters[] = {{L"Surface markings", L"*.markings"}, {L"All files", L"*.*"}};
            dlg->SetFileTypes(2, filters);
            dlg->SetTitle(L"Open surface markings");
        } else if (kind == DialogKind::Mesh) {
            const COMDLG_FILTERSPEC filters[] = {{L"Meshes (OBJ, PLY)", L"*.obj;*.ply"}, {L"All files", L"*.*"}};
            dlg->SetFileTypes(2, filters);
            dlg->SetTitle(L"Open a mesh to knit");
        } else {
            dlg->SetTitle(L"Choose the export folder");
        }
        std::error_code ec;
        if (!start.empty() && std::filesystem::is_directory(start, ec)) {
            IShellItem* folder = nullptr;
            if (SUCCEEDED(SHCreateItemFromParsingName(widen(start).c_str(), nullptr, IID_PPV_ARGS(&folder)))) {
                dlg->SetFolder(folder);
                folder->Release();
            }
        }
        HWND hwnd = owner ? glfwGetWin32Window(owner) : nullptr;
        if (SUCCEEDED(dlg->Show(hwnd))) {
            IShellItem* item = nullptr;
            if (SUCCEEDED(dlg->GetResult(&item))) {
                PWSTR p = nullptr;
                if (SUCCEEDED(item->GetDisplayName(SIGDN_FILESYSPATH, &p))) {
                    chosen = narrow(p);
                    CoTaskMemFree(p);
                }
                item->Release();
            }
        }
        dlg->Release();
    }
    if (SUCCEEDED(init)) CoUninitialize();
    return chosen;
}

enum class SaveKind { Markings, StitchMesh };

// Where to save a markings file or a stitch mesh (IFileSaveDialog; asks before overwriting). `name`: the file name offered.
std::optional<std::string> runSaveDialog(GLFWwindow* owner, const std::string& start, SaveKind kind, const std::string& name) {
    const HRESULT init = CoInitializeEx(nullptr, COINIT_APARTMENTTHREADED | COINIT_DISABLE_OLE1DDE);
    std::optional<std::string> chosen;
    IFileSaveDialog* dlg = nullptr;
    if (SUCCEEDED(CoCreateInstance(CLSID_FileSaveDialog, nullptr, CLSCTX_INPROC_SERVER, IID_PPV_ARGS(&dlg)))) {
        DWORD flags = 0;
        dlg->GetOptions(&flags);
        dlg->SetOptions(flags | FOS_FORCEFILESYSTEM | FOS_OVERWRITEPROMPT);
        const COMDLG_FILTERSPEC markings[] = {{L"Surface markings", L"*.markings"}};
        // the offered name's extension is the file type's (the caller's: io/stitch_mesh_obj.h)
        const std::wstring extension = widen(std::filesystem::path(name).extension().string().substr(name.empty() ? 0 : 1)), pattern = L"*." + extension;
        const COMDLG_FILTERSPEC stitchMesh[] = {{L"Stitch mesh (Wavefront OBJ)", pattern.c_str()}};
        switch (kind) {
        case SaveKind::Markings:
            dlg->SetFileTypes(1, markings);
            dlg->SetDefaultExtension(L"markings");
            dlg->SetTitle(L"Save surface markings");
            break;
        case SaveKind::StitchMesh:
            dlg->SetFileTypes(1, stitchMesh);
            dlg->SetDefaultExtension(extension.c_str());
            dlg->SetTitle(L"Export Stitch Mesh");
            break;
        }
        if (!name.empty()) dlg->SetFileName(widen(name).c_str());
        std::error_code ec;
        if (!start.empty() && std::filesystem::is_directory(start, ec)) {
            IShellItem* folder = nullptr;
            if (SUCCEEDED(SHCreateItemFromParsingName(widen(start).c_str(), nullptr, IID_PPV_ARGS(&folder)))) {
                dlg->SetFolder(folder);
                folder->Release();
            }
        }
        HWND hwnd = owner ? glfwGetWin32Window(owner) : nullptr;
        if (SUCCEEDED(dlg->Show(hwnd))) {
            IShellItem* item = nullptr;
            if (SUCCEEDED(dlg->GetResult(&item))) {
                PWSTR p = nullptr;
                if (SUCCEEDED(item->GetDisplayName(SIGDN_FILESYSPATH, &p))) {
                    chosen = narrow(p);
                    CoTaskMemFree(p);
                }
                item->Release();
            }
        }
        dlg->Release();
    }
    if (SUCCEEDED(init)) CoUninitialize();
    return chosen;
}

} // namespace

bool nativeDialogsAvailable() { return true; }
std::optional<std::string> openMarkingsDialog(GLFWwindow* owner, const std::string& start) { return runDialog(owner, start, DialogKind::Markings); }
std::optional<std::string> saveMarkingsDialog(GLFWwindow* owner, const std::string& start) { return runSaveDialog(owner, start, SaveKind::Markings, ""); }
std::optional<std::string> saveStitchMeshDialog(GLFWwindow* owner, const std::string& start, const std::string& name) {
    return runSaveDialog(owner, start, SaveKind::StitchMesh, name);
}
std::optional<std::string> openMeshDialog(GLFWwindow* owner, const std::string& start) { return runDialog(owner, start, DialogKind::Mesh); }
std::optional<std::string> openImageDialog(GLFWwindow* owner, const std::string& start) { return runDialog(owner, start, DialogKind::Image); }
std::optional<std::string> pickFolderDialog(GLFWwindow* owner, const std::string& start) { return runDialog(owner, start, DialogKind::Folder); }
bool openInFileManager(const std::string& path) {
    return reinterpret_cast<INT_PTR>(ShellExecuteW(nullptr, L"open", widen(path).c_str(), nullptr, nullptr, SW_SHOWNORMAL)) > 32;
}

#else

bool nativeDialogsAvailable() { return false; }
std::optional<std::string> openMarkingsDialog(GLFWwindow*, const std::string&) { return std::nullopt; }
std::optional<std::string> saveMarkingsDialog(GLFWwindow*, const std::string&) { return std::nullopt; }
std::optional<std::string> saveStitchMeshDialog(GLFWwindow*, const std::string&, const std::string&) { return std::nullopt; }
std::optional<std::string> openMeshDialog(GLFWwindow*, const std::string&) { return std::nullopt; }
std::optional<std::string> openImageDialog(GLFWwindow*, const std::string&) { return std::nullopt; }
std::optional<std::string> pickFolderDialog(GLFWwindow*, const std::string&) { return std::nullopt; }
bool openInFileManager(const std::string&) { return false; }

#endif

} // namespace viewer
