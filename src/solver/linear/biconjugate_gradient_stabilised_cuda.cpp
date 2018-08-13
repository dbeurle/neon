
#include "solver/linear/biconjugate_gradient_stabilised_cuda.hpp"

#ifdef ENABLE_CUDA

#include "exceptions.hpp"
#include "dmatrix_vector_product.hpp"
#include "numeric/float_compare.hpp"
#include "solver/cuda_error.hpp"

#include <cmath>
#include <iostream>

namespace neon
{
biconjugate_gradient_stabilised_cuda::~biconjugate_gradient_stabilised_cuda()
{
    // Free device memory
    cudaFree(d_r);
    cudaFree(d_r0);
    cudaFree(d_p);
    cudaFree(d_p_hat);
    cudaFree(d_t);
    cudaFree(d_s);
    cudaFree(d_s_hat);
    cudaFree(d_v);

    cudaFree(d_M_inv);
}

void biconjugate_gradient_stabilised_cuda::solve(sparse_matrix const& A, vector& x, vector const& b)
{
    allocate_device_memory(A, x, b);

    constexpr double one = 1.0;
    constexpr double zero = 0.0;

    auto const N = A.cols();

    x.setZero();

    cuda::check(cudaMemcpy(d_x, x.data(), N * sizeof(double), cudaMemcpyHostToDevice));
    cuda::check(cudaMemcpy(d_r, b.data(), N * sizeof(double), cudaMemcpyHostToDevice));
    cuda::check(
        cudaMemcpy(d_val, A.valuePtr(), A.nonZeros() * sizeof(double), cudaMemcpyHostToDevice));

    {
        vector const M_inv = A.diagonal().cwiseInverse();
        cuda::check(cudaMemcpy(d_M_inv, M_inv.data(), N * sizeof(double), cudaMemcpyHostToDevice));
    }

    // d_z is our r_0 from the algorithm
    cuda::check(cudaMemcpy(d_r0, d_r, N * sizeof(double), cudaMemcpyDeviceToDevice));

    std::int32_t k = 0;
    double residual_old = 0.0, residual = 0.0;

    // Initial residual (needs to be constant through the iterations)
    // r_0 = b - A * x;

    // r.transpose() * z
    cuda::check(cublasDdot(cublasHandle, N, d_r, 1, d_r0, 1, &residual));

    auto alpha = 0.0, omega = 0.0;

    while (std::sqrt(residual) > residual_tolerance && k < max_iterations)
    {
        k++;

        if (k == 1)
        {
            // p = r
            cuda::check(cudaMemcpy(d_p, d_r, N * sizeof(double), cudaMemcpyDeviceToDevice));
        }
        else
        {
            // beta = (rho / rho_old) * (alpha / omega)
            double beta = residual / residual_old * alpha / omega;

            // Overall expression: p = r + beta * (p - omega * v)

            // v <= omega * v
            cuda::check(cublasDscal(cublasHandle, N, &omega, d_v, 1));

            //  p = p - v
            auto constexpr minus_one = -one;
            cuda::check(cublasDaxpy(cublasHandle, N, &minus_one, d_v, 1, d_p, 1));

            // p <= beta * p
            cuda::check(cublasDscal(cublasHandle, N, &beta, d_p, 1));

            //  p = r + p
            cuda::check(cublasDaxpy(cublasHandle, N, &one, d_r, 1, d_p, 1));
        }

        // Solve M * p_hat = p
        // p_hat <- M_inv * p
        cuda::check(cudaMemcpy(d_p_hat, d_p, N * sizeof(double), cudaMemcpyDeviceToDevice));
        cuda::diagonal_matrix_vector_product(d_M_inv, d_p_hat, N);

        // v = A * p_hat
        // Perform y = alpha * A * x + beta * y
        // Results: v <- A * p_hat
        cuda::check(cusparseDcsrmv(cusparseHandle,
                                   CUSPARSE_OPERATION_NON_TRANSPOSE,
                                   N, // rows
                                   N, // cols
                                   A.nonZeros(),
                                   &one, // alpha
                                   descr,
                                   d_val,   // values
                                   d_row,   // row ptr
                                   d_col,   // col ptr
                                   d_p_hat, // x
                                   &zero,   // beta
                                   d_v));   // y

        // alpha = rho / r_0.transpose() * v
        double r_0_dot_v;
        cuda::check(cublasDdot(cublasHandle, N, d_r0, 1, d_v, 1, &r_0_dot_v));

        alpha = residual / r_0_dot_v;

        // s = - alpha * v + r
        cuda::check(cudaMemcpy(d_s, d_r, N * sizeof(double), cudaMemcpyDeviceToDevice));
        double minus_alpha = -alpha;
        cuda::check(cublasDaxpy(cublasHandle, N, &minus_alpha, d_v, 1, d_s, 1));

        // ||s||_2
        double norm_s;
        cuda::check(cublasDnrm2(cublasHandle, N, d_s, 1, &norm_s));

        if (norm_s < residual_tolerance)
        {
            // x = x + alpha * p_hat
            cuda::check(cublasDaxpy(cublasHandle, N, &alpha, d_p_hat, 1, d_x, 1));
            residual = norm_s * norm_s;
            break;
        }

        // s_hat = M_inv * s
        cuda::check(cudaMemcpy(d_s_hat, d_s, N * sizeof(double), cudaMemcpyDeviceToDevice));
        cuda::diagonal_matrix_vector_product(d_M_inv, d_s_hat, N);

        // t = A * s_hat
        // Perform y = alpha * A * x + beta * y
        // Results: t <- A * s_hat
        cuda::check(cusparseDcsrmv(cusparseHandle,
                                   CUSPARSE_OPERATION_NON_TRANSPOSE,
                                   N, // rows
                                   N, // cols
                                   A.nonZeros(),
                                   &one, // alpha
                                   descr,
                                   d_val,   // values
                                   d_row,   // row ptr
                                   d_col,   // col ptr
                                   d_s_hat, // x
                                   &zero,   // beta
                                   d_t));   // y

        // omega = t.transpose() * s / t.transpose() * t
        double t_dot_s, t_dot_t;
        cuda::check(cublasDdot(cublasHandle, N, d_t, 1, d_s, 1, &t_dot_s));
        cuda::check(cublasDdot(cublasHandle, N, d_t, 1, d_t, 1, &t_dot_t));

        omega = t_dot_s / t_dot_t;

        // x = x + alpha * p_hat + omega * s_hat
        cuda::check(cublasDaxpy(cublasHandle, N, &alpha, d_p_hat, 1, d_x, 1));
        cuda::check(cublasDaxpy(cublasHandle, N, &omega, d_s_hat, 1, d_x, 1));

        // r = s - omega * t
        double minus_omega = -omega;
        cuda::check(cudaMemcpy(d_r, d_s, N * sizeof(double), cudaMemcpyDeviceToDevice));
        cuda::check(cublasDaxpy(cublasHandle, N, &minus_omega, d_t, 1, d_r, 1));

        residual_old = residual;

        // r_old.transpose() * r and store in residual
        cuda::check(cublasDdot(cublasHandle, N, d_r0, 1, d_r, 1, &residual));

        if (is_approx(omega, 0.0))
        {
            throw computational_error("BiCGStab solver not converged");
        }
    }

    std::cout << std::string(6, ' ') << "BiCGStab iterations: " << k << " (max. " << max_iterations
              << "), estimated error: " << std::sqrt(residual) << " (min. " << residual_tolerance
              << ")\n";

    if (k >= max_iterations)
    {
        throw computational_error("BiCGStab solver maximum iterations reached");
    }

    // Copy device solution to the host
    cuda::check(cudaMemcpy(x.data(), d_x, N * sizeof(double), cudaMemcpyDeviceToHost));
}

void biconjugate_gradient_stabilised_cuda::allocate_device_memory(sparse_matrix const& A,
                                                                  vector& x,
                                                                  vector const& b)
{
    // If this isn't our first time using the compute device or
    // the sparsity pattern hasn't changed, then we save on the allocation
    if (!build_sparsity_pattern) return;

    auto const N = A.cols();

    cuda::check(cudaMalloc((void**)&d_col, A.nonZeros() * sizeof(int)));
    cuda::check(cudaMalloc((void**)&d_row, (N + 1) * sizeof(int)));
    cuda::check(cudaMalloc((void**)&d_val, A.nonZeros() * sizeof(double)));

    cuda::check(cudaMalloc((void**)&d_x, N * sizeof(double)));

    cuda::check(cudaMalloc((void**)&d_r, N * sizeof(double)));
    cuda::check(cudaMalloc((void**)&d_r0, N * sizeof(double)));
    cuda::check(cudaMalloc((void**)&d_t, N * sizeof(double)));
    cuda::check(cudaMalloc((void**)&d_s, N * sizeof(double)));
    cuda::check(cudaMalloc((void**)&d_s_hat, N * sizeof(double)));
    cuda::check(cudaMalloc((void**)&d_p, N * sizeof(double)));
    cuda::check(cudaMalloc((void**)&d_p_hat, N * sizeof(double)));
    cuda::check(cudaMalloc((void**)&d_v, N * sizeof(double)));
    cuda::check(cudaMalloc((void**)&d_M_inv, N * sizeof(double)));

    // Copy across the data we only need once
    cuda::check(cudaMemcpy(d_row, A.outerIndexPtr(), (N + 1) * sizeof(int), cudaMemcpyHostToDevice));
    cuda::check(
        cudaMemcpy(d_col, A.innerIndexPtr(), A.nonZeros() * sizeof(int), cudaMemcpyHostToDevice));

    build_sparsity_pattern = false;
}
}
#endif
