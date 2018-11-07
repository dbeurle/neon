
#include "solver/find_compute_device.hpp"

#include "io/json.hpp"

#include <algorithm>
#include <iostream>
#include <vector>

#ifdef ENABLE_CUDA

#include <cuda_runtime.h>

/// Check if the CUDA device exists at runtime
/// \return true if CUDA device exists, false otherwise
bool neon::has_cuda_device(json const& device_options)
{
    int device_count;
    cudaGetDeviceCount(&device_count);
    return device_count;
}

#endif

#ifdef ENABLE_OPENCL

neon::opencl_context::opencl_context(std::size_t platform_index,
                                     std::size_t device_index,
                                     cl::Device device)
    : m_platform_index(platform_index), m_device_index(device_index), m_device(device)
{
}

auto neon::opencl_context::platform_index() const noexcept -> std::size_t
{
    return m_platform_index;
}

auto neon::opencl_context::device_index() const noexcept -> std::size_t { return m_device_index; }

auto neon::opencl_context::device() const noexcept -> cl::Device const& { return m_device; }

auto neon::find_opencl_device(json const& device_options) -> opencl_context
{
    if (device_options.find("device") == end(device_options)
        || device_options.find("platform") == end(device_options))
    {
        throw std::domain_error("An OpenCL \"device\" and \"platform\" must be specified");
    }

    if (!device_options["device"].is_number() || !device_options["platform"].is_number())
    {
        throw std::domain_error("\"device\" and \"platform\" must be integers");
    }

    std::size_t const platform_index = device_options["platform"];
    std::size_t const device_index = device_options["device"];

    // Populate with the available platforms
    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);

    if (platforms.empty())
    {
        throw invalid_device("No OpenCL platforms found.  Do you have an OpenCL driver installed "
                             "and a working environment?");
    }

    // check the validity of opencl device and platform
    if (platform_index >= platforms.size())
    {
        throw std::domain_error("Platform number " + std::to_string(platform_index)
                                + " is not valid\n");
    }

    std::vector<cl::Device> devices;
    platforms[platform_index].getDevices(CL_DEVICE_TYPE_ALL, &devices);

    if (device_index >= devices.size())
    {
        throw std::domain_error("Device number " + std::to_string(device_index) + " is not valid\n");
    }
    return opencl_context(platform_index, device_index, devices[device_index]);
}

#endif
