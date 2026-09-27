#pragma once

// Device management for CUDA. Provides gpu_available() and error checking macro.
// When SAMAYA_HAVE_CUDA is not defined, gpu_available() returns false (stub).

#include <cuda_runtime.h>

#include <cstdio>
#include <stdexcept>
#include <string>

namespace samaya {
namespace gpu {

#ifdef __CUDACC__
// Sum over a warp by shuffles; lane 0 holds the total (the other lanes partial sums).
__inline__ __device__ double warp_reduce_sum(double val) {
  for (int offset = 16; offset > 0; offset >>= 1) {
    val += __shfl_down_sync(0xffffffffu, val, offset);
  }
  return val;
}
#endif

// Returns true if a CUDA-capable GPU is available.
inline bool gpu_available() {
#ifdef SAMAYA_HAVE_CUDA
  int count = 0;
  cudaError_t err = cudaGetDeviceCount(&count);
  return err == cudaSuccess && count > 0;
#else
  return false;
#endif
}

// Returns the name of the current GPU device, or empty string if none.
inline std::string device_name() {
#ifdef SAMAYA_HAVE_CUDA
  int device = 0;
  cudaError_t err = cudaGetDevice(&device);
  if (err != cudaSuccess) return "";
  cudaDeviceProp prop;
  err = cudaGetDeviceProperties(&prop, device);
  if (err != cudaSuccess) return "";
  return prop.name;
#else
  return "";
#endif
}

// Sets the active CUDA device.
inline void set_device(int device) {
#ifdef SAMAYA_HAVE_CUDA
  cudaError_t err = cudaSetDevice(device);
  if (err != cudaSuccess) {
    throw std::runtime_error("cudaSetDevice failed: " + std::string(cudaGetErrorString(err)));
  }
#else
  (void)device;
  throw std::runtime_error("CUDA not available");
#endif
}

}  // namespace gpu
}  // namespace samaya

// Macro that turns any CUDA error into an exception with file, line, and error string.
#ifdef SAMAYA_HAVE_CUDA
#define SAMAYA_CUDA_CHECK(call)                                                   \
  do {                                                                            \
    cudaError_t err = (call);                                                     \
    if (err != cudaSuccess) {                                                     \
      throw std::runtime_error(std::string("CUDA error at ") + __FILE__ + ":" +   \
                               std::to_string(__LINE__) + " in " + #call + ": " + \
                               cudaGetErrorString(err));                          \
    }                                                                             \
  } while (0)
#else
#define SAMAYA_CUDA_CHECK(call)                                                               \
  do {                                                                                        \
    (void)(call);                                                                             \
    throw std::runtime_error("CUDA call " #call " invoked but SAMAYA_HAVE_CUDA not defined"); \
  } while (0)
#endif