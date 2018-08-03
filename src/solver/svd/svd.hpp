
#pragma once

#include "numeric/dense_matrix.hpp"

namespace neon
{
/// svd computes the singular value decomposition of an input matrix and
/// provides access to the left and right singular vectors and the singular
/// values.
class svd
{
public:
    svd();

    virtual ~svd() = default;

    virtual void compute(matrix const& A) = 0;

    virtual vector const& left() const noexcept = 0;

    virtual vector const& right() const noexcept = 0;

    virtual vector const& values() const noexcept = 0;
};

/// jacobi_svd uses the Jacobi method for computing the singular values and is
/// not recommended for large decompositions.  For larger systems please see the
/// GPU accelerated algorithms.
class jacobi_svd : public svd
{
public:
    jacobi_svd() = default;

    void compute(matrix const& A) override;

    vector const& left() const noexcept override;

    vector const& right() const noexcept override;

    vector const& values() const noexcept override;
};
}
