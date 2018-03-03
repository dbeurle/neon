
#define CATCH_CONFIG_MAIN

#include <catch.hpp>

#include "solver/eigen/eigen_solver.hpp"
#include "io/json.hpp"

using namespace neon;

constexpr auto ZERO_MARGIN = 1.0e-5;

/** Create a SPD matrix for solver testing */
sparse_matrix create_diagonal_sparse_matrix(int const N)
{
    sparse_matrix A(N, N);

    for (int i = 0; i < N; ++i)
    {
        A.insert(i, i) = static_cast<double>(i + 1);
    }
    A.finalize();

    return A;
}

sparse_matrix create_sparse_identity(int const N)
{
    sparse_matrix A(N, N);

    for (int i = 0; i < N; ++i)
    {
        A.insert(i, i) = 1.0;
    }
    A.finalize();

    return A;
}

TEST_CASE("Eigen solver test suite")
{
    sparse_matrix A = create_diagonal_sparse_matrix(100);
    sparse_matrix I = create_sparse_identity(100);

    eigen_solver eigen{10};

    auto const [values, vectors] = eigen.solve(A, I);

    REQUIRE(values.size() == 10);

    // for (auto i = 0; i < 10; ++i)
    // {
    //     REQUIRE(values[i] == 100 - i);
    // }
}
