
#pragma once

#ifdef ENABLE_OPENCL

#include "linear_solver.hpp"

namespace neon
{
/// conjugate_gradient_ocl is a GPU based solver using the conjugate
/// gradient solver from Matrix Computations.
class conjugate_gradient_ocl : public iterative_linear_solver
{
public:
    /// Default tolerance and maximum iterations
    explicit conjugate_gradient_ocl();

    explicit conjugate_gradient_ocl(double const residual_tolerance);

    explicit conjugate_gradient_ocl(int const max_iterations);

    explicit conjugate_gradient_ocl(double const residual_tolerance, int const max_iterations);

    virtual ~conjugate_gradient_ocl();

    virtual void solve(sparse_matrix const& A, vector& x, vector const& b) override;
};
}
#endif
