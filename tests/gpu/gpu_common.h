#pragma once
// GPU port tests: shared helpers.
#include <string>

#include "../harness.h"
#include "gpu/device.h"

// SKIPs the case without a usable device.
#define REQUIRE_GPU() do { std::string why_; if (!knit::gpu::available(&why_)) SKIP("no CUDA device / backend not built"); } while (0)
