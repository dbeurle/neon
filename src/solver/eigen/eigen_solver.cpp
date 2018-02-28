
#include "eigen_solver.hpp"
#include "Exceptions.hpp"

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
    // using solver_type = Eigen::SparseLU<sparse_matrix, Eigen::AMDOrdering<std::int32_t>>;

    Eigen::SparseMatrix<double> A_col = A;

    std::cout << "Different ordering of A\n" << A_col << "\n";

    Eigen::ArpackGeneralizedSelfAdjointEigenSolver<Eigen::SparseMatrix<double>> arpack;

    arpack.compute(A_col, eigenvalues_to_extract, "SM");

    if (arpack.getNbrConvergedEigenValues() < eigenvalues_to_extract)
    {
        throw computational_error("Eigenvalues did not converge");
    }

    if (arpack.info() != Eigen::Success)
    {
        throw computational_error("Numerical issued occurred");
    }

    return {arpack.eigenvectors(), arpack.eigenvalues()};
}
}
