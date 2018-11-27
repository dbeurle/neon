
#pragma once

#include "solver/eigen/eigen_solver.hpp"

namespace neon
{
/// power_iteration performs an iterative method to return a single eigenvalue
/// and an eigenvector.  This algorithm is not very robust and the Lanczos
/// algorithm is preferred.
class power_iteration : public eigen_solver
{
public:
    power_iteration(std::int64_t const values_to_extract);

    virtual void solve(sparse_matrix const& A) override final;

    virtual void solve(sparse_matrix const& A, sparse_matrix const& B) override final;
};
}
