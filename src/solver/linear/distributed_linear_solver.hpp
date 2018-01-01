
#pragma once

#include "numeric/DenseMatrix.hpp"
#include "numeric/SparseMatrix.hpp"

namespace neon
{
/**
 * distributed_linear_solver is to setup a linear solver with designated
 * parameters from the input file.  This is the interface for every linear
 * solver.  Requires distributed mapping.
 */
class distributed_linear_solver
{
public:
    virtual void solve(SparseMatrix const& A, vector& x, vector const& b) = 0;

    /** Notifies the linear solvers to of a change in sparsity structure of A */
    void update_sparsity_pattern() { build_sparsity_pattern = true; }

protected:
    bool build_sparsity_pattern = true;
};

class distributed_direct_linear_solver : public distributed_linear_solver
{
};
}
