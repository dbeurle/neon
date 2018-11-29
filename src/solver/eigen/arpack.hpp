
#pragma once

#include "solver/eigen/eigen_solver.hpp"

namespace neon
{
class arpack : public eigen_solver
{
public:
    /// Construct an arpack eigenvalue solver
    /// \param values_to_extract Number of eigenvalues to extract
    arpack(std::int64_t const values_to_extract,
           eigen_solver::eigen_spectrum const spectrum = eigen_solver::eigen_spectrum::lower);

    /// Solve the standard eigenvalue problem $\f (A - \lambda I) x = 0 $\f
    /// \return eigenvalues and eigenvectors
    virtual void solve(sparse_matrix const& A) override final;

    /// Solve the generalised eigenvalue problem $\f (A - \lambda B) x = 0 $\f
    /// \return eigenvalues and eigenvectors
    virtual void solve(sparse_matrix const& A, sparse_matrix const& B) override final;
};
}
