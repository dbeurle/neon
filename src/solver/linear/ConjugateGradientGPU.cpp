
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

#include <chrono>
#include <cmath>

// CUDA Runtime
#include <cuda/cuda_runtime.h>

// Using updated (v2) interfaces for CUBLAS and CUSPARSE
#include <cuda/cublas_v2.h>
#include <cuda/cusparse.h>

// Utilities and system includes
#include "helper_cuda.h"      // helper for CUDA error checking
#include "helper_functions.h" // helper for shared functions common to CUDA
// Samples

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

void ConjugateGradientGPU::solve(SparseMatrix const& A,
                                 Vector& x_double,
                                 Vector const& b_double)
{
    constexpr int max_iter = 5000;
    constexpr float tol = 1e-12f;

    constexpr float floatone = 1.0;
    constexpr float floatzero = 0.0;

    auto const N = A.cols();

    // Cast all vectors to floating point types
    Eigen::VectorXf const val = A.coeffs().cast<float>();
    Eigen::VectorXf const b = b_double.cast<float>();
    Eigen::VectorXf x = x_double.cast<float>();

    // column pointer
    // const_cast<pastix_int_t*>(A.IsRowMajor ? A.innerIndexPtr()
    //                                        : A.outerIndexPtr()),
    // row pointer
    // const_cast<pastix_int_t*>(A.IsRowMajor ? A.outerIndexPtr()
    //                                        : A.innerIndexPtr()),
    // const_cast<double*>(A.valuePtr()),

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
    float *d_val, *d_x;
    float *d_r, *d_p, *d_omega, *d_y;

    checkCudaErrors(cudaMalloc((void**)&d_col, A.nonZeros() * sizeof(int)));
    checkCudaErrors(cudaMalloc((void**)&d_row, (N + 1) * sizeof(int)));
    checkCudaErrors(cudaMalloc((void**)&d_val, A.nonZeros() * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**)&d_x, N * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**)&d_y, N * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**)&d_r, N * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**)&d_p, N * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**)&d_omega, N * sizeof(float)));

    cudaMemcpy(d_col, A.innerIndexPtr(), A.nonZeros() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_row, A.outerIndexPtr(), (N + 1) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_val, val.data(), A.nonZeros() * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_x, x.data(), N * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_r, b.data(), N * sizeof(float), cudaMemcpyHostToDevice);

    //  Conjugate gradient without preconditioning.
    //    ------------------------------------------
    //    Follows the description by Golub & Van Loan, "Matrix Computations 3rd ed.",
    //    Section 10.2.6

    int k = 0;
    float r0 = 0.0f, r1 = 0.0f;

    cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);

    while (std::sqrt(r1) > tol && k < max_iter)
    {
        k++;
        if (k == 1)
        {
            cublasScopy(cublasHandle, N, d_r, 1, d_p, 1);
        }
        else
        {
            float beta = r1 / r0;
            cublasSscal(cublasHandle, N, &beta, d_p, 1);
            cublasSaxpy(cublasHandle, N, &floatone, d_r, 1, d_p, 1);
        }
        cusparseScsrmv(cusparseHandle,
                       CUSPARSE_OPERATION_NON_TRANSPOSE,
                       N,
                       N,
                       A.nonZeros(),
                       &floatone,
                       descr,
                       d_val,
                       d_row,
                       d_col,
                       d_p,
                       &floatzero,
                       d_omega);

        float dot = -1.0f;

        cublasSdot(cublasHandle, N, d_p, 1, d_omega, 1, &dot);

        float alpha = r1 / dot;

        cublasSaxpy(cublasHandle, N, &alpha, d_p, 1, d_x, 1);

        float nalpha = -alpha;

        cublasSaxpy(cublasHandle, N, &nalpha, d_omega, 1, d_r, 1);

        r0 = r1;

        cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
    }

    std::cout << "Conjugate gradient iteration = " << k
              << ", residual = " << std::sqrt(r1) << " \n";
    std::cout << "Convergence Test: " << (k <= max_iter ? "OK" : "FAIL") << " \n";

    // Copy device solution to the host
    cudaMemcpy(x.data(), d_x, N * sizeof(float), cudaMemcpyDeviceToHost);
    x_double = x.cast<double>();

    // // Preconditioned Conjugate Gradient using ILU.
    //    --------------------------------------------
    //    Follows the description by Golub & Van Loan, "Matrix Computations 3rd ed.",
    //    Algorithm 10.3.1
    //
    // std::cout << ("\nConvergence of conjugate gradient using incomplete LU "
    //               "preconditioning: "
    //               "\n");
    //
    // float numerator, denominator;
    // float* d_valsILU0;
    // float *d_zm1, *d_zm2, *d_rm2;
    //
    // int nzILU0 = 2 * N - 1;
    //
    // checkCudaErrors(cudaMalloc((void**)&d_valsILU0, A.nonZeros() * sizeof(float)));
    // checkCudaErrors(cudaMalloc((void**)&d_zm1, (N) * sizeof(float)));
    // checkCudaErrors(cudaMalloc((void**)&d_zm2, (N) * sizeof(float)));
    // checkCudaErrors(cudaMalloc((void**)&d_rm2, (N) * sizeof(float)));
    //
    // // create the analysis info object for the A matrix
    // cusparseSolveAnalysisInfo_t infoA = 0;
    // cusparseStatus = cusparseCreateSolveAnalysisInfo(&infoA);
    //
    // checkCudaErrors(cusparseStatus);
    //
    // // Perform the analysis for the Non-Transpose case
    // cusparseStatus = cusparseScsrsv_analysis(cusparseHandle,
    //                                          CUSPARSE_OPERATION_NON_TRANSPOSE,
    //                                          N,
    //                                          A.nonZeros(),
    //                                          descr,
    //                                          d_val,
    //                                          d_row,
    //                                          d_col,
    //                                          infoA);
    //
    // checkCudaErrors(cusparseStatus);
    //
    // // Copy A data to ILU0 vals as input
    // cudaMemcpy(d_valsILU0, d_val, A.nonZeros() * sizeof(float),
    // cudaMemcpyDeviceToDevice);
    //
    // // generate the Incomplete LU factor H for the matrix A using cudsparseScsrilu0
    // cusparseStatus = cusparseScsrilu0(cusparseHandle,
    //                                   CUSPARSE_OPERATION_NON_TRANSPOSE,
    //                                   N,
    //                                   descr,
    //                                   d_valsILU0,
    //                                   d_row,
    //                                   d_col,
    //                                   infoA);
    //
    // checkCudaErrors(cusparseStatus);
    //
    // // Create info objects for the ILU0 preconditioner
    // cusparseSolveAnalysisInfo_t info_u;
    // cusparseCreateSolveAnalysisInfo(&info_u);
    //
    // cusparseMatDescr_t descrL = 0;
    // cusparseStatus = cusparseCreateMatDescr(&descrL);
    // cusparseSetMatType(descrL, CUSPARSE_MATRIX_TYPE_GENERAL);
    // cusparseSetMatIndexBase(descrL, CUSPARSE_INDEX_BASE_ZERO);
    // cusparseSetMatFillMode(descrL, CUSPARSE_FILL_MODE_LOWER);
    // cusparseSetMatDiagType(descrL, CUSPARSE_DIAG_TYPE_UNIT);
    //
    // cusparseMatDescr_t descrU = 0;
    // cusparseStatus = cusparseCreateMatDescr(&descrU);
    // cusparseSetMatType(descrU, CUSPARSE_MATRIX_TYPE_GENERAL);
    // cusparseSetMatIndexBase(descrU, CUSPARSE_INDEX_BASE_ZERO);
    // cusparseSetMatFillMode(descrU, CUSPARSE_FILL_MODE_UPPER);
    // cusparseSetMatDiagType(descrU, CUSPARSE_DIAG_TYPE_NON_UNIT);
    // cusparseStatus = cusparseScsrsv_analysis(cusparseHandle,
    //                                          CUSPARSE_OPERATION_NON_TRANSPOSE,
    //                                          N,
    //                                          A.nonZeros(),
    //                                          descrU,
    //                                          d_val,
    //                                          d_row,
    //                                          d_col,
    //                                          info_u);
    //
    // // reset the initial guess of the solution to zero
    // for (int i = 0; i < N; i++)
    // {
    //     x[i] = 0.0;
    // }
    //
    // checkCudaErrors(cudaMemcpy(d_r, b.data(), N * sizeof(float),
    // cudaMemcpyHostToDevice)); checkCudaErrors(cudaMemcpy(d_x, x.data(), N *
    // sizeof(float), cudaMemcpyHostToDevice));
    //
    // k = 0;
    // cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
    //
    // while (r1 > tol * tol && k <= max_iter)
    // {
    //     // Forward Solve, we can re-use infoA since the sparsity pattern of A matches
    //     // that of L
    //     cusparseStatus = cusparseScsrsv_solve(cusparseHandle,
    //                                           CUSPARSE_OPERATION_NON_TRANSPOSE,
    //                                           N,
    //                                           &floatone,
    //                                           descrL,
    //                                           d_valsILU0,
    //                                           d_row,
    //                                           d_col,
    //                                           infoA,
    //                                           d_r,
    //                                           d_y);
    //     checkCudaErrors(cusparseStatus);
    //
    //     // Back Substitution
    //     cusparseStatus = cusparseScsrsv_solve(cusparseHandle,
    //                                           CUSPARSE_OPERATION_NON_TRANSPOSE,
    //                                           N,
    //                                           &floatone,
    //                                           descrU,
    //                                           d_valsILU0,
    //                                           d_row,
    //                                           d_col,
    //                                           info_u,
    //                                           d_y,
    //                                           d_zm1);
    //     checkCudaErrors(cusparseStatus);
    //
    //     k++;
    //
    //     if (k == 1)
    //     {
    //         cublasScopy(cublasHandle, N, d_zm1, 1, d_p, 1);
    //     }
    //     else
    //     {
    //         cublasSdot(cublasHandle, N, d_r, 1, d_zm1, 1, &numerator);
    //         cublasSdot(cublasHandle, N, d_rm2, 1, d_zm2, 1, &denominator);
    //         float beta = numerator / denominator;
    //         cublasSscal(cublasHandle, N, &beta, d_p, 1);
    //         cublasSaxpy(cublasHandle, N, &floatone, d_zm1, 1, d_p, 1);
    //     }
    //
    //     cusparseScsrmv(cusparseHandle,
    //                    CUSPARSE_OPERATION_NON_TRANSPOSE,
    //                    N,
    //                    N,
    //                    nzILU0,
    //                    &floatone,
    //                    descrU,
    //                    d_val,
    //                    d_row,
    //                    d_col,
    //                    d_p,
    //                    &floatzero,
    //                    d_omega);
    //     cublasSdot(cublasHandle, N, d_r, 1, d_zm1, 1, &numerator);
    //     cublasSdot(cublasHandle, N, d_p, 1, d_omega, 1, &denominator);
    //     float alpha = numerator / denominator;
    //     cublasSaxpy(cublasHandle, N, &alpha, d_p, 1, d_x, 1);
    //     cublasScopy(cublasHandle, N, d_r, 1, d_rm2, 1);
    //     cublasScopy(cublasHandle, N, d_zm1, 1, d_zm2, 1);
    //     float nalpha = -alpha;
    //     cublasSaxpy(cublasHandle, N, &nalpha, d_omega, 1, d_r, 1);
    //     cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
    // }
    //
    // std::cout << "  iteration = " << k << ", residual = " << std::sqrt(r1) << " \n";
    //
    // cudaMemcpy(x.data(), d_x, N * sizeof(float), cudaMemcpyDeviceToHost);
    //
    // // check result
    // err = 0.0;
    //
    // for (int i = 0; i < N; i++)
    // {
    //     float rsum = 0.0f;
    //
    //     for (int j = I[i]; j < I[i + 1]; j++)
    //     {
    //         rsum += val[j] * x[J[j]];
    //     }
    //
    //     float diff = std::abs(rsum - b[i]);
    //
    //     if (diff > err)
    //     {
    //         err = diff;
    //     }
    // }
    //
    // std::cout << "  Convergence Test: " << (k <= max_iter ? "OK" : "FAIL") << " \n";
    // nErrors += (k > max_iter) ? 1 : 0;
    // float qaerr2 = err;

    // Destroy parameters for preconditioner
    // cusparseDestroySolveAnalysisInfo(infoA);
    // cusparseDestroySolveAnalysisInfo(info_u);

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

    // cudaFree(d_valsILU0);
    // cudaFree(d_zm1);
    // cudaFree(d_zm2);
    // cudaFree(d_rm2);
    // std::cout << "  Test Summary:\n";
    // std::cout << "     qaerr1 = " << std::abs(qaerr1) << "\n";
    //" qaerr2 = " << std::abs(qaerr2)
    //<< "\n\n";
    // exit((nErrors == 0 && std::abs(qaerr1) < 1e-5 && std::abs(qaerr2) < 1e-5
    //           ? EXIT_SUCCESS
    //           : EXIT_FAILURE));
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
