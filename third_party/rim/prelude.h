#pragma once
// Every header RIM's sources include from outside themselves (the standard library, Eigen, pcg32, the serial
// TBB stand-in and the nanogui serializer stub in third_party/rim/shim), included
// here in the GLOBAL namespace before a wrapped unit opens namespace rim and includes a RIM source (CMake
// generates one unit per source). Their include guards keep the copies inside the namespace empty, so RIM's
// own classes (MultiResolutionHierarchy, Optimizer, ...) become rim:: and cannot collide with anything else
// in the program. A header missing here fails to compile inside the namespace (loudly).
// common.h's settings first (it sets them before Eigen).
#if defined(_WIN32) && !defined(NOMINMAX)
#define NOMINMAX
#endif
#define EIGEN_DONT_PARALLELIZE

#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <chrono>
#include <climits>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <map>
#include <math.h>
#include <memory>
#include <mutex>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <pcg32.h>
#include <tbb/tbb.h>
#include <nanogui/serializer/core.h>
#include <nanogui/serializer/sparse.h>
