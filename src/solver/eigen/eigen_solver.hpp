
#pragma once

#include "numeric/dense_matrix.hpp"
#include "numeric/sparse_matrix.hpp"

#include <cstdint>
#include <utility>

namespace neon
{
class eigen_solver
{
public:
    /// Which end of the eigenvalue spectrum to compute
    enum class eigen_spectrum { upper, lower };

public:
    /// Construct an eigenvalue solver
    /// \param values_to_extract Number of eigenvalues to extract
    eigen_solver(std::int64_t const values_to_extract,
                 eigen_spectrum const spectrum = eigen_spectrum::lower);

    virtual ~eigen_solver() = default;

    /// Solve the standard eigenvalue problem $\f (A - \lambda I) x = 0 $\f
    /// \return eigenvalues and eigenvectors
    virtual void solve(sparse_matrix const& A) = 0;

    /// Solve the generalised eigenvalue problem $\f (A - \lambda B) x = 0 $\f
    /// \return eigenvalues and eigenvectors
    virtual void solve(sparse_matrix const& A, sparse_matrix const& B) = 0;

    vector const& eigenvalues() const noexcept { return m_eigenvalues; }

    col_matrix const& eigenvectors() const noexcept { return m_eigenvectors; }

protected:
    std::int64_t values_to_extract{0};

    eigen_spectrum m_spectrum;

    vector m_eigenvalues;
    col_matrix m_eigenvectors;
};
}
