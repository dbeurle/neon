
#pragma once

#include "LinearSolver.hpp"

namespace neon
{
/**
 * MUMPS is a multifrontal direct solver
 */
class MUMPS : public DirectLinearSolver
{
public:
    void solve(SparseMatrix const& A, Vector& x, Vector const& b) override final;
};
}
