// CUDA device management implementation.

#include <cuda_runtime.h>

#include "device.hpp"

namespace samaya {
namespace gpu {

// These functions are already implemented in device.hpp as inline.
// This file exists to ensure the CUDA runtime is linked and to provide
// a place for any non-inline device management code in the future.

}  // namespace gpu
}  // namespace samaya