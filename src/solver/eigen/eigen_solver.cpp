
#include "eigen_solver.hpp"

#include <iostream>
#include <unsupported/Eigen/ArpackSupport>

namespace neon
{
eigen_solver::eigen_solver(std::int64_t const eigenvalues_to_extract)
    : eigenvalues_to_extract{eigenvalues_to_extract}
{
}

std::pair<vector, matrix> eigen_solver::solve(sparse_matrix const& A, sparse_matrix const& B)
{
    // using solver_type = Eigen::SimplicialLU<Eigen::SparseMatrix<sparse_matrix::Scalar>>;
    using solver_type = Eigen::SparseLU<sparse_matrix, Eigen::AMDOrdering<std::int32_t>>;

    Eigen::ArpackGeneralizedSelfAdjointEigenSolver<sparse_matrix, solver_type>
        arpack(A, B, eigenvalues_to_extract);

    std::cout << "Solved some eigenvalues apparently\n";

    std::cout << arpack.eigenvectors() << std::endl;

    return {arpack.eigenvectors(), arpack.eigenvalues()};
}
}
