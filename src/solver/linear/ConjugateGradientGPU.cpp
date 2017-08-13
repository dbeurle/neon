
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

#ifdef ENABLE_CUDA

#include "ConjugateGradientGPU.hpp"

#include <chrono>
#include <cmath>

#include <cuda/cuda_runtime.h>

#include <cuda/cublas_v2.h>
#include <cuda/cusparse.h>

#include "helper_cuda.h"
#include "helper_functions.h"

namespace neon
{
ConjugateGradientGPU::ConjugateGradientGPU() { this->find_compute_device(); }

ConjugateGradientGPU::ConjugateGradientGPU(double const tol) : ConjugateGradientGPU()
{
    solverParam.tolerance = tol;
}

ConjugateGradientGPU::ConjugateGradientGPU(int const maxIter) : ConjugateGradientGPU()
{
    solverParam.max_iterations = maxIter;
}

ConjugateGradientGPU::ConjugateGradientGPU(double const tol, int const maxIter)
    : ConjugateGradientGPU()
{
    solverParam.max_iterations = maxIter;
    solverParam.tolerance = tol;
}

void ConjugateGradientGPU::solve(SparseMatrix const& A, Vector& x, Vector const& b)
{
    auto const max_iter = solverParam.max_iterations;
    auto const tol = solverParam.tolerance;

    constexpr double one = 1.0;
    constexpr double zero = 0.0;

    auto const N = A.cols();

    x.setZero();

    // Create CUBLAS context
    cublasHandle_t cublasHandle = 0;
    cublasStatus_t cublasStatus = cublasCreate(&cublasHandle);
    checkCudaErrors(cublasStatus);

    // Create CUSPARSE context
    cusparseHandle_t cusparseHandle = 0;
    cusparseStatus_t cusparseStatus = cusparseCreate(&cusparseHandle);
    checkCudaErrors(cusparseStatus);

    // Description of the A matrix
    cusparseMatDescr_t descr = 0;
    cusparseStatus = cusparseCreateMatDescr(&descr);
    checkCudaErrors(cusparseStatus);

    // Define the properties of the matrix
    cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);

    // Allocate required memory
    int *d_col, *d_row;
    double *d_val, *d_x;
    double *d_r, *d_p, *d_omega, *d_y;

    checkCudaErrors(cudaMalloc((void**)&d_col, A.nonZeros() * sizeof(int)));
    checkCudaErrors(cudaMalloc((void**)&d_row, (N + 1) * sizeof(int)));
    checkCudaErrors(cudaMalloc((void**)&d_val, A.nonZeros() * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_x, N * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_y, N * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_r, N * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_p, N * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_omega, N * sizeof(double)));

    cudaMemcpy(d_col, A.innerIndexPtr(), A.nonZeros() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_row, A.outerIndexPtr(), (N + 1) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_val, A.valuePtr(), A.nonZeros() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_x, x.data(), N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_r, b.data(), N * sizeof(double), cudaMemcpyHostToDevice);

    int k = 0;
    double r0 = 0.0, r1 = 0.0;
    cublasDdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);

    while (std::sqrt(r1) > tol && k < max_iter)
    {
        k++;

        if (k == 1)
        {
            cublasDcopy(cublasHandle, N, d_r, 1, d_p, 1);
        }
        else
        {
            double beta = r1 / r0;
            cublasDscal(cublasHandle, N, &beta, d_p, 1);
            cublasDaxpy(cublasHandle, N, &one, d_r, 1, d_p, 1);
        }

        // Perform A * x = b
        cusparseDcsrmv(cusparseHandle,
                       CUSPARSE_OPERATION_NON_TRANSPOSE,
                       N,
                       N,
                       A.nonZeros(),
                       &one,
                       descr,
                       d_val,
                       d_row,
                       d_col,
                       d_p,
                       &zero,
                       d_omega);

        double dot = 0.0;

        cublasDdot(cublasHandle, N, d_p, 1, d_omega, 1, &dot);

        double alpha = r1 / dot;

        cublasDaxpy(cublasHandle, N, &alpha, d_p, 1, d_x, 1);

        double nalpha = -alpha;

        cublasDaxpy(cublasHandle, N, &nalpha, d_omega, 1, d_r, 1);

        r0 = r1;

        cublasDdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
    }

    std::cout << "Conjugate gradient iteration = " << k
              << ", residual = " << std::sqrt(r1) << " \n";
    std::cout << "Convergence Test: " << (k <= max_iter ? "OK" : "FAIL") << " \n";

    // Copy device solution to the host
    cudaMemcpy(x.data(), d_x, N * sizeof(double), cudaMemcpyDeviceToHost);

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
    cudaFree(d_omega);
}

void ConjugateGradientGPU::find_compute_device()
{
    // Pick the best possible CUDA capable device
    int devID = findCudaDevice(0, nullptr);

    std::cout << "GPU selected Device ID = " << devID << " \n";

    if (devID < 0)
    {
        throw std::runtime_error("Invalid GPU device " + std::to_string(devID)
                                 + " selected,  exiting...\n");
    }

    cudaDeviceProp deviceProp;
    checkCudaErrors(cudaGetDeviceProperties(&deviceProp, devID));

    // Statistics about the GPU device
    std::cout << "> GPU device has " << deviceProp.multiProcessorCount
              << " Multi-Processors, SM " << deviceProp.major << "." << deviceProp.minor
              << " compute capabilities\n\n";

    if (deviceProp.major * 0x10 + deviceProp.minor < 0x11)
    {
        throw std::runtime_error("Linear solver requires a minimum CUDA compute 1.1 "
                                 "capability\n");
    }
}
}
#endif
