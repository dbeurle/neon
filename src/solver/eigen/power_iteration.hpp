
#pragma once

/// @file

#include "solver/eigen/eigen_solver.hpp"

namespace neon
{
/// power_iteration performs an iterative method to return a single eigenvalue
/// and an eigenvector.  This algorithm is not very robust and the Lanczos
/// algorithm is preferred.
class power_iteration : public eigen_solver
{
public:
    power_iteration(std::int64_t,
                    eigen_solver::eigen_spectrum const spectrum = eigen_solver::eigen_spectrum::lower);

    virtual void solve(sparse_matrix const& A) override final;

    virtual void solve(sparse_matrix const& A, sparse_matrix const& B) override final;
};
}
