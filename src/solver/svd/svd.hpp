
#pragma once

/// @file

#include "numeric/dense_matrix.hpp"
#include <Eigen/SVD>

namespace neon
{
/// svd computes the singular value decomposition of an input matrix and
/// provides access to the left and right singular vectors and the singular
/// values.
class svd
{
public:
    svd() = default;

    virtual ~svd() = default;

    /// compute thin SVD decomposition of a matrix A
    virtual void compute(col_matrix const& A) = 0;

    /// Compute thin SVD decomposition
    /// \param A Input matrix
    /// \param n Maximum number of singular vectors
    virtual void compute(col_matrix const& A, std::int64_t const n) = 0;

    /// compute thin SVD decomposition of a matrix A with a tolerance on the singular values
    /// A singular value will be considered nonzero if its value is strictly greater than
    /// \f$ \vert singular value \vert \leqslant threshold \times \vert max singular value \vert \f$.
    virtual void compute(col_matrix const& A, double const tolerance) = 0;

    /// \return left singular vectors (columns of the thin U matrix)
    virtual col_matrix const& left() const noexcept = 0;

    /// \return right singular vectors (columns of the thin V matrix)
    virtual col_matrix const& right() const noexcept = 0;

    /// \return singular values
    virtual vector const& values() const noexcept = 0;

    /// A least-squares solution of A*x = b
    /// linear compination of columns of A
    virtual void solve(vector& x, vector const& b) const noexcept = 0;

protected:
    col_matrix left_vectors;
    col_matrix right_vectors;
    vector singular_values;
};

/// bdc_svd first reduces the input matrix to bi-diagonal form and then performs a
/// divide-and-conquer diagonalization. Small blocks are diagonalized using class
/// JacobiSVD. Default switching size is 16.
class bdc_svd : public svd
{
public:
    bdc_svd() = default;

    bdc_svd(col_matrix const& A);

    void compute(col_matrix const& A) override;

    void compute(col_matrix const& A, std::int64_t const n) override;

    void compute(col_matrix const& A, double const tolerance) override;

    col_matrix const& left() const noexcept override;

    col_matrix const& right() const noexcept override;

    vector const& values() const noexcept override;

    void solve(vector& x, vector const& b) const noexcept override;

private:
    Eigen::BDCSVD<col_matrix> decomposition;
};

/// Implementation of the truncated Singular Value Decomposition, using
/// randomized algorithms as described in 'finding structure with randomness'
/// @cite halko2011finding.
class randomised_svd : public svd
{
public:
    randomised_svd() = default;

    randomised_svd(col_matrix const& A);

    void compute(col_matrix const& A) override;

    void compute(col_matrix const& A, std::int64_t const n) override;

    void compute(col_matrix const& A, double const tolerance) override;

    col_matrix const& left() const noexcept override;

    col_matrix const& right() const noexcept override;

    vector const& values() const noexcept override;

    void solve(vector& x, vector const& b) const noexcept override;

private:
    Eigen::BDCSVD<col_matrix> decomposition;
};
}
