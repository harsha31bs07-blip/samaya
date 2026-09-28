// Which GPU samaya can use for PDLP (--gpu). Plain C++: with SAMAYA_CUDA off it reports that the
// build has no CUDA kernels.

#include "samaya/solver.hpp"

#ifdef SAMAYA_HAVE_CUDA
#include "gpu/device.hpp"
#endif

namespace samaya {

GpuInfo gpu_info() {
  GpuInfo info;
#ifdef SAMAYA_HAVE_CUDA
  info.built_with_cuda = true;
  if (gpu::gpu_available()) info.device = gpu::device_name();
#endif
  return info;
}

}  // namespace samaya
