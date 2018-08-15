
#pragma once

#ifdef ENABLE_OCL

#include "solver/linear/conjugate_gradient_ocl.hpp"

namespace neon
{
/// Implementation of the Bi-Conjugate Gradient Stabilised (BiCGStab) algorithm
/// using the ViennaCL library
class biconjugate_gradient_stabilised_ocl : public conjugate_gradient_ocl
{
public:
    using conjugate_gradient_ocl::conjugate_gradient_ocl;

    virtual ~biconjugate_gradient_stabilised_ocl() = default;

    void solve(sparse_matrix const& A, vector& x, vector const& b) override final;
};
}
#endif
