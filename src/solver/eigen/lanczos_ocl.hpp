
#pragma once

#include "solver/eigen/eigen_solver.hpp"

namespace neon
{
/// lanczos_ocl performs an iterative method to return a number of eigenvalues
/// and eigenvectors.  This algorithm suffers from numerical error if all the
/// eigenvalues are desired, but it is sufficient when computing only a small
/// range, e.g. in natural frequency analysis.  It exploits parallelism through
/// OpenCL as provided by ViennaCL
class lanczos_ocl : public eigen_solver
{
public:
    lanczos_ocl(std::int64_t const values_to_extract);

    virtual void solve(sparse_matrix const& A) override final;

    virtual void solve(sparse_matrix const& A, sparse_matrix const& B) override final;
};
}
