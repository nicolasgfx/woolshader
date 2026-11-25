#pragma once
// Native file / folder dialogs and "open in the file manager" (Windows: IFileOpenDialog,
// ShellExecute). Elsewhere they report unavailable and the viewer falls back to drag and drop
// and typed paths. Modal: they block the calling (main) thread until closed.
#include <optional>
#include <string>

struct GLFWwindow;

namespace viewer {

bool nativeDialogsAvailable();
// An OBJ or PLY file; nothing if cancelled or unavailable. `start`: the folder to open in ("": default).
std::optional<std::string> openMeshDialog(GLFWwindow* owner, const std::string& start);
// A greyscale / colour image (PNG, JPEG, TGA, BMP): the fur coverage mask.
std::optional<std::string> openImageDialog(GLFWwindow* owner, const std::string& start);
// A markings file (the Mark step's Open ...), and where to save one (Save as ...; the extension is added).
std::optional<std::string> openMarkingsDialog(GLFWwindow* owner, const std::string& start);
std::optional<std::string> saveMarkingsDialog(GLFWwindow* owner, const std::string& start);
// Where to save a stitch mesh (the Mark step's Export Stitch Mesh, an OBJ file); `name`: the file name offered, whose
// extension is the dialog's file type.
std::optional<std::string> saveStitchMeshDialog(GLFWwindow* owner, const std::string& start, const std::string& name);
// A folder (the export destination's parent).
std::optional<std::string> pickFolderDialog(GLFWwindow* owner, const std::string& start);
// Shows a folder in the file manager.
bool openInFileManager(const std::string& path);

} // namespace viewer
