
/*
 * Copyright 1993-2015 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.

 * This sample implements a preconditioned conjugate gradient solver on
 * the GPU using CUBLAS and CUSPARSE.  Relative to the conjugateGradient
 * SDK example, this demonstrates the use of cusparseScsrilu0() for
 * computing the incompute-LU preconditioner and cusparseScsrsv_solve()
 * for solving triangular systems.  Specifically, the preconditioned
 * conjugate gradient method with an incomplete LU preconditioner is
 * used to solve the Laplacian operator in 2D on a uniform mesh.
 *
 * Note that the code in this example and the specific matrices used here
 * were chosen to demonstrate the use of the CUSPARSE library as simply
 * and as clearly as possible.  This is not optimized code and the input
 * matrices have been chosen for simplicity rather than performance.
 * These should not be used either as a performance guide or for
 * benchmarking purposes.
 */

#include "conjugate_gradientGPU.hpp"

#ifdef ENABLE_CUDA

#include "Exceptions.hpp"
#include "dmatrix_vector_product.hpp"

#include <cmath>
#include <cuda_runtime.h>

/** Exception for bad computation in CUDA kernel */
class gpu_kernel_error : public std::domain_error
{
    using std::domain_error::domain_error;
};

namespace neon
{
conjugate_gradientGPU::conjugate_gradientGPU() : iterative_linear_solver()
{
    this->find_compute_device();

    // Create CUBLAS context
    (cublasCreate(&cublasHandle));

    // Create CUSPARSE context
    (cusparseCreate(&cusparseHandle));

    // Description of the A matrix
    (cusparseCreateMatDescr(&descr));

    // Define the properties of the matrix
    cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);
}

conjugate_gradientGPU::conjugate_gradientGPU(double const residual_tolerance)
    : conjugate_gradientGPU()
{
    this->residual_tolerance = residual_tolerance;
}

conjugate_gradientGPU::conjugate_gradientGPU(int const max_iterations) : conjugate_gradientGPU()
{
    this->max_iterations = max_iterations;
}

conjugate_gradientGPU::conjugate_gradientGPU(double const residual_tolerance, int const max_iterations)
    : conjugate_gradientGPU()
{
    this->residual_tolerance = residual_tolerance;
    this->max_iterations = max_iterations;
}

conjugate_gradientGPU::~conjugate_gradientGPU()
{
    // Destroy contexts
    cusparseDestroy(cusparseHandle);
    cublasDestroy(cublasHandle);

    // Free device memory
    cudaFree(d_col);
    cudaFree(d_row);
    cudaFree(d_val);
    cudaFree(d_x);
    cudaFree(d_y);
    cudaFree(d_r);
    cudaFree(d_p);
    cudaFree(d_Ap);
    cudaFree(d_z);
    cudaFree(d_M_inv);
}

void conjugate_gradientGPU::solve(sparse_matrix const& A, vector& x, vector const& b)
{
    this->allocate_device_memory(A, x, b);

    constexpr double one = 1.0;
    constexpr double zero = 0.0;

    auto const N = A.cols();

    x.setZero();

    cudaMemcpy(d_x, x.data(), N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_r, b.data(), N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_val, A.valuePtr(), A.nonZeros() * sizeof(double), cudaMemcpyHostToDevice);

    {
        vector const M_inv = A.diagonal().cwiseInverse();
        cudaMemcpy(d_M_inv, M_inv.data(), N * sizeof(double), cudaMemcpyHostToDevice);
    }

    int k = 0;
    double residual_old = 0.0, residual = 0.0;

    cudaMemcpy(d_z, d_r, N * sizeof(double), cudaMemcpyDeviceToDevice);

    // compute z0
    cuda::diagonal_matrix_vector_product(d_M_inv, d_z, N);

    // r.transpose() * z
    cublasDdot(cublasHandle, N, d_r, 1, d_z, 1, &residual);

    while (std::sqrt(residual) > residual_tolerance && k < max_iterations)
    {
        k++;

        if (k == 1)
        {
            // p = z
            cublasDcopy(cublasHandle, N, d_z, 1, d_p, 1);
        }
        else
        {
            double beta = residual / residual_old;
            // p <= beta * p
            cublasDscal(cublasHandle, N, &beta, d_p, 1);

            //  p = p + z
            cublasDaxpy(cublasHandle, N, &one, d_z, 1, d_p, 1);
        }

        // Perform y = alpha * A * x + beta * y
        // Results A * p
        cusparseDcsrmv(cusparseHandle,
                       CUSPARSE_OPERATION_NON_TRANSPOSE,
                       N, // rows
                       N, // cols
                       A.nonZeros(),
                       &one, // alpha
                       descr,
                       d_val, // values
                       d_row, // row ptr
                       d_col, // col ptr
                       d_p,   // x
                       &zero, // beta
                       d_Ap); // y

        double p_dot_Ap = 0.0;

        // p' * Ap
        cublasDdot(cublasHandle, N, d_p, 1, d_Ap, 1, &p_dot_Ap);

        // residual == old residual
        double alpha = residual / p_dot_Ap;

        // x <= x + alpha * p
        cublasDaxpy(cublasHandle, N, &alpha, d_p, 1, d_x, 1);

        // r = r - alpha * Ap;
        double nalpha = -alpha;
        cublasDaxpy(cublasHandle, N, &nalpha, d_Ap, 1, d_r, 1);

        residual_old = residual;

        // z = M_inv * r
        cudaMemcpy(d_z, d_r, N * sizeof(double), cudaMemcpyDeviceToDevice);
        cuda::diagonal_matrix_vector_product(d_M_inv, d_z, N);

        // z' * r and store in residual
        cublasDdot(cublasHandle, N, d_z, 1, d_r, 1, &residual);
    }

    std::cout << std::string(6, ' ') << "Conjugate Gradient iterations: " << k << " (max. "
              << max_iterations << "), estimated error: " << std::sqrt(residual) << " (min. "
              << residual_tolerance << ")\n";

    if (k >= max_iterations)
    {
        throw computational_error("Conjugate gradient solver maximum iterations reached");
    }

    // Copy device solution to the host
    cudaMemcpy(x.data(), d_x, N * sizeof(double), cudaMemcpyDeviceToHost);
}

void conjugate_gradientGPU::allocate_device_memory(sparse_matrix const& A, vector& x, vector const& b)
{
    // If this isn't our first time using the compute device or
    // the sparsity pattern hasn't changed, then we save on the allocation
    if (!build_sparsity_pattern) return;

    auto const N = A.cols();

    (cudaMalloc((void**)&d_col, A.nonZeros() * sizeof(int)));
    (cudaMalloc((void**)&d_row, (N + 1) * sizeof(int)));
    (cudaMalloc((void**)&d_val, A.nonZeros() * sizeof(double)));
    (cudaMalloc((void**)&d_x, N * sizeof(double)));
    (cudaMalloc((void**)&d_y, N * sizeof(double)));
    (cudaMalloc((void**)&d_r, N * sizeof(double)));
    (cudaMalloc((void**)&d_p, N * sizeof(double)));
    (cudaMalloc((void**)&d_Ap, N * sizeof(double)));
    (cudaMalloc((void**)&d_z, N * sizeof(double)));
    (cudaMalloc((void**)&d_M_inv, N * sizeof(double)));

    cudaMemcpy(d_row, A.outerIndexPtr(), (N + 1) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_col, A.innerIndexPtr(), A.nonZeros() * sizeof(int), cudaMemcpyHostToDevice);

    build_sparsity_pattern = false;
}

void conjugate_gradientGPU::find_compute_device()
{
    cudaDeviceProp device_properties;

    std::memset(&device_properties, 0, sizeof(device_properties));

    // NOTE Here you can set the required properties of the device.  These
    // will then be chosen based on the preferences in cudaDeviceProp

    int device_number;
    cudaChooseDevice(&device_number, &device_properties);

    if (device_number < 0)
    {
        throw std::domain_error("No CUDA device found.");
    }

    cudaSetDevice(device_number);

    cudaGetDeviceProperties(&device_properties, device_number);

    // Statistics about the GPU device
    std::cout << std::string(' ', 4) << "GPU device has " << device_properties.multiProcessorCount
              << " Multi-Processors, SM " << device_properties.major << "."
              << device_properties.minor << " compute capabilities\n\n";

    if (device_properties.major * 0x10 + device_properties.minor < 0x11)
    {
        throw std::runtime_error("Linear solver requires a minimum CUDA compute 1.1 "
                                 "capability\n");
    }
}
}
#endif
