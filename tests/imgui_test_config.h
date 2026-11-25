#pragma once
// ImGui configuration of the UI tests (knit_tests_ui, IMGUI_USER_CONFIG): every item registration
// reaches the test hooks (tests/tests_ui.cpp), and an ImGui assertion (stack imbalance, empty ID
// at a window root, ...) throws instead of aborting, so the failing case is reported.
#include <stdexcept>

#define IMGUI_ENABLE_TEST_ENGINE
#define IM_ASSERT(_EXPR) ((_EXPR) ? (void)0 : throw std::runtime_error("IM_ASSERT: " #_EXPR))
