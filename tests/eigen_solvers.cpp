
#define CATCH_CONFIG_MAIN

#include <catch.hpp>

#include "solver/eigen/eigen_solver.hpp"
#include "io/json.hpp"

using namespace neon;

constexpr auto ZERO_MARGIN = 1.0e-5;

/** Create a SPD matrix for solver testing */
sparse_matrix create_sparse_matrix(std::int64_t const size)
{
    sparse_matrix A(size, size);
    for (std::int64_t i{0}; i < size; ++i)
    {
        A.insert(i, i) = static_cast<sparse_matrix::value_type>(i + 1);
    }
    A.finalize();
    return A;
}

sparse_matrix create_sparse_identity(std::int64_t const size)
{
    sparse_matrix A(size, size);
    for (std::int64_t i{0}; i < size; ++i)
    {
        A.insert(i, i) = 1.0;
    }
    A.finalize();
    return A;
}

TEST_CASE("Eigen solver test suite")
{
    sparse_matrix const A = create_sparse_matrix(2);
    sparse_matrix const I = create_sparse_identity(2);

    eigen_solver eigen{1};

    auto const [values, vectors] = eigen.solve(A, I);

    std::cout << "Eigenvalues\n" << values << std::endl;
    std::cout << "Eigenvectors\n" << vectors << std::endl;

    REQUIRE(1 == 1);
}
