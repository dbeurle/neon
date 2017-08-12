

#include "LinearSolver.hpp"

namespace neon
{
/**
 * ConjugateGradientGPU is a GPU based solver using the preconditioned conjugate
 * gradient solver from the CUDA examples.  The preconditioner is the diagonal
 * preconditioner.
 */
class ConjugateGradientGPU : public LinearSolver
{
public:
    ConjugateGradientGPU();

    ConjugateGradientGPU(double const residual_tolerance);

    ConjugateGradientGPU(int const maxIter);

    ConjugateGradientGPU(double const residual_tolerance, int const maxIter);

    void solve(SparseMatrix const& A, Vector& x, Vector const& b) override final;

protected:
    void find_compute_device();
};
}
