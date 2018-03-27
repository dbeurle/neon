
#pragma once

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

    virtual void solve(sparse_matrix const& A, vector& x, vector const& b) override;

protected:
    void find_compute_device();

    // Device side pointers
    std::int32_t* d_col{nullptr}; /// Column pointer
    std::int32_t* d_row{nullptr}; /// Row pointer (compressed)
    double* d_val{nullptr};       /// Nonzero coefficients of A

    double* d_x{nullptr}; /// Solution vector

    cublasHandle_t cublasHandle = 0;
    cusparseHandle_t cusparseHandle = 0;
    cusparseMatDescr_t descr = 0;

private:
    void allocate_device_memory(sparse_matrix const& A, vector& x, vector const& b);

private:
    double *d_r{nullptr}, *d_p{nullptr}, *d_Ap{nullptr}, *d_y{nullptr}, *d_z{nullptr};
    double* d_M_inv{nullptr}; /// Diagonal preconditioner storage
};

/// Implementation of the Bi-Conjugate Gradient Stabilised (BiCGStab) algorithm
/// using the CUDA BLAS operations
class biconjugate_gradient_stabilised_cuda : public conjugate_gradient_cuda
{
public:
    using conjugate_gradient_cuda::conjugate_gradient_cuda;

    ~biconjugate_gradient_stabilised_cuda();

    void solve(sparse_matrix const& A, vector& x, vector const& b) override final;

private:
    void allocate_device_memory(sparse_matrix const& A, vector& x, vector const& b);

private:
    double *d_r{nullptr}, *d_r0, *d_t, *d_s_hat, *d_s, *d_p, *d_p_hat, *d_v;
    double* d_M_inv; /// Diagonal preconditioner storage
};
}
#endif
