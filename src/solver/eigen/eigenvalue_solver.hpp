
#pragma once

#include "numeric/dense_matrix.hpp"
#include "numeric/sparse_matrix.hpp"

#include <cstdint>
#include <utility>

namespace neon
{
class eigenvalue_solver
{
public:
    /// Construct an eigenvalue solver
    /// \param values_to_extract Number of eigenvalues to extract
    eigenvalue_solver(std::int64_t const values_to_extract);

    /// Solve the standard eigenvalue problem $\f (A - \lambda I) x = 0 $\f
    /// \return eigenvalues and eigenvectors
    virtual void solve(sparse_matrix const& A);

    /// Solve the generalised eigenvalue problem $\f (A - \lambda B) x = 0 $\f
    /// \return eigenvalues and eigenvectors
    virtual void solve(sparse_matrix const& A, sparse_matrix const& B);

    vector const& eigenvalues() const noexcept { return m_eigenvalues; }

    col_matrix const& eigenvectors() const noexcept { return m_eigenvectors; }

protected:
    std::int64_t values_to_extract{};

    vector m_eigenvalues;
    col_matrix m_eigenvectors;
};

/// power_iteration performs an iterative method to return a single eigenvalue
/// and an eigenvector.  This algorithm is not very robust and the Lanczos
/// algorithm is preferred.
class power_iteration : public eigenvalue_solver
{
public:
    power_iteration(std::int64_t const values_to_extract);

    void solve(sparse_matrix const& A) override final;

    void solve(sparse_matrix const& A, sparse_matrix const& B) override final;
};

/// lanzcos_solver performs an iterative method to return a number of eigenvalues
/// and eigenvectors.  This algorithm suffers from numerical error if all the
/// eigenvalues are desired, but it is sufficient when computing only a small
/// range, e.g. in natural frequency analysis.
class lanzcos_solver : public eigenvalue_solver
{
public:
    lanzcos_solver(std::int64_t const values_to_extract);

    void solve(sparse_matrix const& A) override final;

    void solve(sparse_matrix const& A, sparse_matrix const& B) override final;
};
}
