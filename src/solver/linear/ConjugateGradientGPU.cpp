
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

#include "ConjugateGradientGPU.hpp"

#ifdef ENABLE_CUDA

#include <cmath>

#include <cuda/cuda_runtime.h>

#include "helper_cuda.h"
#include "helper_functions.h"

namespace neon
{
ConjugateGradientGPU::ConjugateGradientGPU()
{
    this->find_compute_device();

    // Create CUBLAS context
    checkCudaErrors(cublasCreate(&cublasHandle));

    // Create CUSPARSE context
    checkCudaErrors(cusparseCreate(&cusparseHandle));

    // Description of the A matrix
    checkCudaErrors(cusparseCreateMatDescr(&descr));

    // Define the properties of the matrix
    cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);
}

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

ConjugateGradientGPU::~ConjugateGradientGPU()
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
}

void ConjugateGradientGPU::solve(SparseMatrix const& A, Vector& x, Vector const& b)
{
    this->allocate_device_memory(A, x, b);

    auto const max_iter = solverParam.max_iterations;
    auto const tol = solverParam.tolerance;

    constexpr double one = 1.0;
    constexpr double zero = 0.0;

    auto const N = A.cols();

    x.setZero();

    Vector diagonal_entries = A.diagonal().cwiseInverse();

    cudaMemcpy(d_x, x.data(), N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_r, b.data(), N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_z, diagonal_entries.data(), N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_val, A.valuePtr(), A.nonZeros() * sizeof(double), cudaMemcpyHostToDevice);

    int k = 0;
    double r0 = 0.0, r1 = 0.0;

    // r * r.transpose()
    cublasDdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);

    cublasDsbmv(cublasHandle,
                CUBLAS_FILL_MODE_LOWER,
                N,
                0,     // #subdiagonals
                &one,  // scaling factor for diagonal
                d_z,   // Diagonal matrix
                1,     // Leading dimension
                d_r,   // x
                1,     // incx
                &zero, // scaling for no data
                d_z,   // output vector
                1);

    while (std::sqrt(r1) > tol && k < max_iter)
    {
        k++;

        if (k == 1)
        {
            // p = r
            cublasDcopy(cublasHandle, N, d_r, 1, d_p, 1);
        }
        else
        {
            double beta = r1 / r0;
            // p <= p * beta
            cublasDscal(cublasHandle, N, &beta, d_p, 1);

            //  p = p + r
            cublasDaxpy(cublasHandle, N, &one, d_r, 1, d_p, 1);
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

        // r1 == old residual
        double alpha = r1 / p_dot_Ap;

        // x <= alpha * p + x
        cublasDaxpy(cublasHandle, N, &alpha, d_p, 1, d_x, 1);

        double nalpha = -alpha;

        // r = r - alpha * Ap;
        cublasDaxpy(cublasHandle, N, &nalpha, d_Ap, 1, d_r, 1);

        r0 = r1;

        // r' * r and store in r1
        cublasDdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
    }

    std::cout << "Conjugate gradient iteration = " << k
              << ", residual = " << std::sqrt(r1) << " \n";
    std::cout << "Convergence Test: " << (k <= max_iter ? "OK" : "FAIL") << " \n";

    // Copy device solution to the host
    cudaMemcpy(x.data(), d_x, N * sizeof(double), cudaMemcpyDeviceToHost);
}

void ConjugateGradientGPU::allocate_device_memory(SparseMatrix const& A,
                                                  Vector& x,
                                                  Vector const& b)
{
    // If this isn't our first time using the compute device or
    // the sparsity pattern hasn't changed, then we save on the allocation
    if (!build_sparsity_pattern) return;

    auto const N = A.cols();

    checkCudaErrors(cudaMalloc((void**)&d_col, A.nonZeros() * sizeof(int)));
    checkCudaErrors(cudaMalloc((void**)&d_row, (N + 1) * sizeof(int)));
    checkCudaErrors(cudaMalloc((void**)&d_val, A.nonZeros() * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_x, N * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_y, N * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_r, N * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_p, N * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_Ap, N * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_z, N * sizeof(double)));

    cudaMemcpy(d_row, A.outerIndexPtr(), (N + 1) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_col, A.innerIndexPtr(), A.nonZeros() * sizeof(int), cudaMemcpyHostToDevice);

    build_sparsity_pattern = false;
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
