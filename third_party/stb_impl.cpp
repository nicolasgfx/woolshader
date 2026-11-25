// The stb_image and stb_image_write implementations (header-only, FetchContent, pinned in
// cmake/dependencies.cmake): decode the base-colour textures of imported meshes and GLB scenes (PNG,
// JPEG, TGA, BMP) and compress the PNGs this program writes (io/image_io.cpp).
#define STB_IMAGE_IMPLEMENTATION
#define STBI_NO_HDR
#define STBI_NO_PIC
#define STBI_NO_PNM
#define STBI_NO_PSD
#define STBI_NO_GIF
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STBI_WRITE_NO_STDIO
#include "stb_image_write.h"
