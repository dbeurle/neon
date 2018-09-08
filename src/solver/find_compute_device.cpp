
#include "solver/find_compute_device.hpp"

#include <algorithm>
#include <iostream>
#include <vector>

#ifdef ENABLE_CUDA

#include <cuda_runtime.h>

/// Check if the CUDA device exists at runtime
/// \return true if CUDA device exists, false otherwise
bool neon::has_cuda_device()
{
    int device_count;
    cudaGetDeviceCount(&device_count);
    return device_count;
}

#endif

#ifdef ENABLE_OCL

#include "CL/cl.hpp"

namespace
{
/// highest_compute is a comparisor object for finding the device the largest
/// number of compute units
struct highest_compute
{
    auto operator()(cl::Device const& left, cl::Device const& right) const noexcept -> bool
    {
        return left.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>()
               < right.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();
    }
};

/// highest_memory is a comparisor object for finding the device the largest
/// amount of memory
struct highest_memory
{
    auto operator()(cl::Device const& left, cl::Device const& right) const noexcept -> bool
    {
        return left.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>() < right.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>();
    }
};

/// Check that the device is valid according to the computational requirements
/// of the algorithms such as fp64 support and availability.
static bool is_valid(cl::Device const& device) noexcept
{
    return device.getInfo<CL_DEVICE_AVAILABLE>();
    //&& device.getInfo<CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE>();
}

/// Remove any invalid platforms by checking their available devices \sa is_valids
static void remove_invalid(std::vector<cl::Platform>& platforms) noexcept
{
    platforms.erase(std::remove_if(begin(platforms),
                                   end(platforms),
                                   [](cl::Platform const& platform) {
                                       // Check each device for this platform
                                       std::vector<cl::Device> devices;
                                       platform.getDevices(CL_DEVICE_TYPE_GPU | CL_DEVICE_TYPE_CPU,
                                                           &devices);

                                       return devices.empty()
                                              || std::none_of(begin(devices),
                                                              end(devices),
                                                              [](cl::Device const& device) {
                                                                  return is_valid(device);
                                                              });
                                   }),
                    end(platforms));
}

/// Remove invalid devices \sa is_valid
static void remove_invalid(std::vector<cl::Device>& devices) noexcept
{
    devices.erase(std::remove_if(begin(devices),
                                 end(devices),
                                 [](cl::Device const& device) { return !is_valid(device); }),
                  end(devices));
}
}

auto neon::find_opencl_device() -> cl::Device
{
    // Populate with the available platforms
    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);

    if (platforms.empty())
    {
        throw invalid_device("No OpenCL platforms found.  Do you have an OpenCL driver installed?");
    }

    remove_invalid(platforms);

    std::cout << "Number of platforms: " << platforms.size() << "\n";

    std::vector<cl::Device> best_devices;

    for (auto const& platform : platforms)
    {
        std::string output;
        platform.getInfo(CL_PLATFORM_NAME, &output);

        std::cout << "Platform " << output << '\n';

        std::vector<cl::Device> devices;
        platform.getDevices(CL_DEVICE_TYPE_ALL, &devices);

        remove_invalid(devices);

        std::cout << "  available devices " << devices.size() << "\n";

        if (devices.empty())
        {
            continue;
        }

        auto const most_compute = std::max_element(begin(devices), end(devices), highest_compute{});
        if (most_compute != end(devices))
        {
            best_devices.push_back(*most_compute);
        }

        auto const most_memory = std::max_element(begin(devices), end(devices), highest_memory{});
        if (most_memory != end(devices))
        {
            best_devices.push_back(*most_memory);
        }
    }
    return best_devices.front();
}

#endif
