
#pragma once

/// @file

#ifdef ENABLE_CUDA

#include "linear_solver.hpp"

#include <cuda/cublas_v2.h>
#include <cuda/cusparse.h>

namespace neon
{
/// conjugate_gradient_cuda is a GPU based solver using the conjugate
/// gradient solver from Matrix Computations.
class conjugate_gradient_cuda : public iterative_linear_solver
{
public:
    /// Default tolerance and maximum iterations
    explicit conjugate_gradient_cuda();

    explicit conjugate_gradient_cuda(double const residual_tolerance);

    explicit conjugate_gradient_cuda(int const max_iterations);

    explicit conjugate_gradient_cuda(double const residual_tolerance, int const max_iterations);

    /// Deallocate the device memory and cuda BLAS handles
    virtual ~conjugate_gradient_cuda();

    virtual void solve(sparse_matrix const& input_matrix, vector& x, vector const& input_rhs) override;

protected:
    void find_compute_device();

protected:
    /// Device column pointer
    std::int32_t* d_col{nullptr};
    /// Device row pointer (compressed)
    std::int32_t* d_row{nullptr};
    /// Device non-zero coefficients of A
    double* d_val{nullptr};

    /// Device solution vector
    double* d_x{nullptr};

    cublasHandle_t cublasHandle = 0;
    cusparseHandle_t cusparseHandle = 0;
    cusparseMatDescr_t descr = 0;

private:
    void allocate_device_memory(sparse_matrix const& A);

private:
    double *d_r{nullptr}, *d_p{nullptr}, *d_Ap{nullptr}, *d_y{nullptr}, *d_z{nullptr};
    /// Device diagonal preconditioner storage
    double* d_M_inv{nullptr};
};
}
#endif
