
#pragma once

#include "LinearSolver.hpp"

namespace neon
{
/**
 * PaStiX is a supernodal direct solver with multithreading support
 */
class PaStiX : public LinearSolver
{
public:
    void solve(SparseMatrix const& A, Vector& x, Vector const& b) override final;
};
}
