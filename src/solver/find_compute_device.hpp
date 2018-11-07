
#pragma once

#include <stdexcept>
#include "io/json_forward.hpp"

#ifdef ENABLE_OPENCL
#include "CL/cl.hpp"
#endif

/// \file find_compute_device.hpp

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
bool has_cuda_device(json const& device_options);

#endif

#ifdef ENABLE_OPENCL

class opencl_context
{
public:
    explicit opencl_context(std::size_t platform_index, std::size_t device_index, cl::Device device);

    auto platform_index() const noexcept -> std::size_t;

    auto device_index() const noexcept -> std::size_t;

    auto device() const noexcept -> cl::Device const&;

private:
    std::size_t m_platform_index;
    std::size_t m_device_index;

    cl::Device m_device;
};

/// Check if a functional OpenCL CPU runtime exists
/// Throws an \p invalid_device if no device is found
/// \return true if an OpenCL CPU device exists, false otherwise
[[nodiscard]] auto find_opencl_device(json const& device_options) -> opencl_context;

#endif
}
