
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
    eigen_solver(std::int64_t const eigenvalues_to_extract);

    /** Solve the generalised eigenvalue problem $\f (A - \lambda B) I = 0 $\f */
    std::pair<vector, matrix> solve(sparse_matrix const& A, sparse_matrix const& B);

protected:
    std::int64_t eigenvalues_to_extract{1};
};
}
