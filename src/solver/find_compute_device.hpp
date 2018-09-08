
#pragma once

#include <stdexcept>

/// \file find_compute_device.hpp

namespace cl
{
class Device;
}

namespace neon
{
class invalid_device : public std::domain_error
{
public:
    using std::domain_error::domain_error;
};

#ifdef ENABLE_CUDA

/// Check if the CUDA device exists at runtime
/// Throws an \p invalid_device if no device is found
/// \return true if CUDA device exists, false otherwise
bool has_cuda_device();

#endif

#ifdef ENABLE_OCL

/// Check if a functional OpenCL CPU runtime exists
/// Throws an \p invalid_device if no device is found
/// \return true if an OpenCL CPU device exists, false otherwise
[[nodiscard]] auto find_opencl_device() -> cl::Device;

#endif
}
