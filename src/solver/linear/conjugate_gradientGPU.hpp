
#pragma once

#ifdef ENABLE_CUDA

#include "linear_solver.hpp"

#include <cuda/cublas_v2.h>
#include <cuda/cusparse.h>

namespace neon
{
/**
 * conjugate_gradientGPU is a GPU based solver using the conjugate
 * gradient solver from the CUDA examples.
 */
class conjugate_gradientGPU : public iterative_linear_solver
{
public:
    explicit conjugate_gradientGPU();
    explicit conjugate_gradientGPU(double const residual_tolerance);
    explicit conjugate_gradientGPU(int const max_iterations);
    explicit conjugate_gradientGPU(double const residual_tolerance, int const max_iterations);

    using iterative_linear_solver::iterative_linear_solver;

    ~conjugate_gradientGPU();

    void solve(sparse_matrix const& A, vector& x, vector const& b) override final;

protected:
    void allocate_device_memory(sparse_matrix const& A, vector& x, vector const& b);

    void find_compute_device();

protected:
    // Device side pointers
    std::int32_t* d_col{nullptr};
    std::int32_t* d_row{nullptr};
    double* d_val{nullptr};

    double* d_x{nullptr};

    double *d_r, *d_p, *d_Ap, *d_y, *d_z;
    double* d_M_inv;

    cublasHandle_t cublasHandle = 0;
    cusparseHandle_t cusparseHandle = 0;
    cusparseMatDescr_t descr = 0;
};
}
#endif
