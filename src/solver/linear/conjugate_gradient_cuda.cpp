
#include "conjugate_gradient_cuda.hpp"

#ifdef ENABLE_CUDA

#include "exceptions.hpp"
#include "dmatrix_vector_product.hpp"
#include "numeric/float_compare.hpp"
#include "solver/cuda_error.hpp"

#include <cmath>
#include <iostream>

namespace neon
{
conjugate_gradient_cuda::conjugate_gradient_cuda() : iterative_linear_solver()
{
    this->find_compute_device();

    // Create CUBLAS context
    cuda::check(cublasCreate(&cublasHandle));

    // Create CUSPARSE context
    cuda::check(cusparseCreate(&cusparseHandle));

    // Description of the A matrix
    cuda::check(cusparseCreateMatDescr(&descr));

    // Define the properties of the matrix
    cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);
}

conjugate_gradient_cuda::conjugate_gradient_cuda(double const residual_tolerance)
    : conjugate_gradient_cuda()
{
    this->residual_tolerance = residual_tolerance;
}

conjugate_gradient_cuda::conjugate_gradient_cuda(int const max_iterations)
    : conjugate_gradient_cuda()
{
    this->max_iterations = max_iterations;
}

conjugate_gradient_cuda::conjugate_gradient_cuda(double const residual_tolerance,
                                                 int const max_iterations)
    : conjugate_gradient_cuda()
{
    this->residual_tolerance = residual_tolerance;
    this->max_iterations = max_iterations;
}

conjugate_gradient_cuda::~conjugate_gradient_cuda()
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

void conjugate_gradient_cuda::solve(sparse_matrix const& A, vector& x, vector const& b)
{
    this->allocate_device_memory(A, x, b);

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

    int k = 0;
    double residual_old = 0.0, residual = 0.0;

    cuda::check(cudaMemcpy(d_z, d_r, N * sizeof(double), cudaMemcpyDeviceToDevice));

    // compute z0
    cuda::diagonal_matrix_vector_product(d_M_inv, d_z, N);

    // r.transpose() * z
    cuda::check(cublasDdot(cublasHandle, N, d_r, 1, d_z, 1, &residual));

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

void conjugate_gradient_cuda::allocate_device_memory(sparse_matrix const& A, vector& x, vector const& b)
{
    // If this isn't our first time using the compute device or
    // the sparsity pattern hasn't changed, then we save on the allocation
    if (!build_sparsity_pattern) return;

    auto const N = A.cols();

    cuda::check(cudaMalloc((void**)&d_col, A.nonZeros() * sizeof(int)));
    cuda::check(cudaMalloc((void**)&d_row, (N + 1) * sizeof(int)));
    cuda::check(cudaMalloc((void**)&d_val, A.nonZeros() * sizeof(double)));
    cuda::check(cudaMalloc((void**)&d_x, N * sizeof(double)));
    cuda::check(cudaMalloc((void**)&d_y, N * sizeof(double)));
    cuda::check(cudaMalloc((void**)&d_r, N * sizeof(double)));
    cuda::check(cudaMalloc((void**)&d_p, N * sizeof(double)));
    cuda::check(cudaMalloc((void**)&d_Ap, N * sizeof(double)));
    cuda::check(cudaMalloc((void**)&d_z, N * sizeof(double)));
    cuda::check(cudaMalloc((void**)&d_M_inv, N * sizeof(double)));

    cuda::check(cudaMemcpy(d_row, A.outerIndexPtr(), (N + 1) * sizeof(int), cudaMemcpyHostToDevice));
    cuda::check(
        cudaMemcpy(d_col, A.innerIndexPtr(), A.nonZeros() * sizeof(int), cudaMemcpyHostToDevice));

    build_sparsity_pattern = false;
}

void conjugate_gradient_cuda::find_compute_device()
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
    std::cout << std::string(4, ' ') << "GPU device has " << device_properties.multiProcessorCount
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
