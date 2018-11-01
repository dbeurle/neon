
#include <catch2/catch.hpp>

#include "solver/find_compute_device.hpp"

#include <iostream>

#ifdef ENABLE_OCL

#include <CL/cl.hpp>

TEST_CASE("Find an OpenCL device")
{
    // Find a CPU device for testing on the build server
    auto const device = neon::find_opencl_device();

    std::cout << "A wild " << device.getInfo<CL_DEVICE_NAME>() << " appeared!\n";
}

#endif

#ifdef ENABLE_CUDA

TEST_CASE("Find a CUDA device")
{
    // Find a CPU device for testing on the build server
    auto const device = neon::find_cuda_device();
}

#endif
