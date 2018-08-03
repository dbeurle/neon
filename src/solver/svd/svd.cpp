
#include "solver/svd/svd.hpp"

#include <Eigen/SVD>

#include <iostream>

namespace neon
{
svd::svd()
{
    matrix m = matrix::Random(3, 2);

    std::cout << "Here is the matrix m:"
              << "\n"
              << m << "\n";

    Eigen::JacobiSVD<decltype(m)> svd(m, Eigen::ComputeThinU | Eigen::ComputeThinV);

    std::cout << "Its singular values are:"
              << "\n"
              << svd.singularValues() << "\n";

    std::cout << "Its left singular vectors are the columns of the thin U matrix:"
              << "\n"
              << svd.matrixU() << "\n";

    std::cout << "Its right singular vectors are the columns of the thin V matrix:"
              << "\n"
              << svd.matrixV() << "\n";

    vector3 rhs(1, 0, 0);

    std::cout << "Now consider this rhs vector:"
              << "\n"
              << rhs << "\n";

    std::cout << "A least-squares solution of m*x = rhs is:"
              << "\n"
              << svd.solve(rhs) << "\n";
}
}
