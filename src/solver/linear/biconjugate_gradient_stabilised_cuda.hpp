
#pragma once

/// @file

#ifdef ENABLE_CUDA

#include "solver/linear/conjugate_gradient_cuda.hpp"

namespace neon
{
/// Implementation of the Bi-Conjugate Gradient Stabilised (BiCGStab) algorithm
/// using the CUDA BLAS operations
class biconjugate_gradient_stabilised_cuda : public conjugate_gradient_cuda
{
public:
    using conjugate_gradient_cuda::conjugate_gradient_cuda;

    ~biconjugate_gradient_stabilised_cuda();

    void solve(sparse_matrix const& A, vector& x, vector const& b) override final;

private:
    void allocate_device_memory(sparse_matrix const& A, vector const& b);

private:
    double *d_r{nullptr}, *d_r0, *d_t, *d_s_hat, *d_s, *d_p, *d_p_hat, *d_v;
    double* d_M_inv; /// Diagonal preconditioner storage
};
}
#endif
