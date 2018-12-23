
#pragma once

/// @file

#ifdef ENABLE_CUDA

#include "exceptions.hpp"

#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <cuda/cublas_v2.h>
#include <cuda/cusparse.h>

namespace neon::cuda
{
inline void check(cudaError_t&& error_code)
{
    if (error_code != cudaSuccess)
    {
        throw cuda_error(cudaGetErrorString(error_code));
    }
}

inline void check(cublasStatus_t&& error_code)
{
    if (error_code != CUBLAS_STATUS_SUCCESS)
    {
        throw cuda_error("Error occurred in CUBLAS function call");
    }
}

inline void check(cusparseStatus_t&& error_code)
{
    if (error_code != CUSPARSE_STATUS_SUCCESS)
    {
        throw cuda_error("Error occurred in CUSPARSE function call");
    }
}

inline void check(cusolverStatus_t&& status)
{
    if (status != CUSOLVER_STATUS_SUCCESS)
    {
        throw cuda_error("Error occurred in CUSOLVER function call");
    }
}

template <class T>
inline T* malloc(std::size_t const size)
{
    T* pointer;
    check(cudaMalloc((void**)&pointer, size * sizeof(T)));
    return pointer;
}
}

#endif
