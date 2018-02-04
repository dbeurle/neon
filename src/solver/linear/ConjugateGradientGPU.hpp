
#pragma once

#ifdef ENABLE_CUDA

#include "LinearSolver.hpp"

#include <cuda/cublas_v2.h>
#include <cuda/cusparse.h>

namespace neon
{
/**
 * ConjugateGradientGPU is a GPU based solver using the conjugate
 * gradient solver from the CUDA examples.
 */
class ConjugateGradientGPU : public IterativeLinearSolver
{
public:
    explicit ConjugateGradientGPU();
    explicit ConjugateGradientGPU(double const residual_tolerance);
    explicit ConjugateGradientGPU(int const max_iterations);
    explicit ConjugateGradientGPU(double const residual_tolerance, int const max_iterations);

    using IterativeLinearSolver::IterativeLinearSolver;

    ~ConjugateGradientGPU();

    void solve(sparse_matrix const& A, vector& x, vector const& b) override final;

protected:
    void allocate_device_memory(sparse_matrix const& A, vector& x, vector const& b);

    void find_compute_device();

protected:
    // Device side pointers
    int32_t* d_col{nullptr};
    int32_t* d_row{nullptr};
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
