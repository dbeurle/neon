
#include <catch2/catch.hpp>

#include "solver/eigen/eigenvalue_solver.hpp"
#include "io/json.hpp"

/// Create a SPD matrix for solver testing
neon::sparse_matrix create_diagonal_sparse_matrix(int const N)
{
    neon::sparse_matrix A(N, N);

    for (int i = 0; i < N; ++i)
    {
        A.insert(i, i) = i + 1;
    }
    A.finalize();

    return A;
}

neon::sparse_matrix create_sparse_identity(int const N)
{
    neon::sparse_matrix A(N, N);

    for (int i = 0; i < N; ++i)
    {
        A.insert(i, i) = 1.0;
    }
    A.finalize();

    return A;
}

TEST_CASE("Arpack eigenvalues")
{
    neon::eigenvalue_solver solver{10};

    solver.solve(create_diagonal_sparse_matrix(50), create_sparse_identity(50));

    auto const& values = solver.eigenvalues();
    auto const& vectors = solver.eigenvectors();

    REQUIRE(values.size() == 10);

    REQUIRE(vectors.rows() == 50);
    REQUIRE(vectors.cols() == 10);

    for (int i = 0; i < 10; i++)
    {
        // Expected eigenvalues
        REQUIRE(values(i) == Approx(i + 1.0));
        // Unit vectors
        REQUIRE(vectors.col(i).norm() == Approx(1.0));
    }
}
TEST_CASE("Power iteration eigenvalue")
{
    neon::power_iteration solver{1};

    solver.solve(create_diagonal_sparse_matrix(10));

    auto const& values = solver.eigenvalues();
    auto const& vectors = solver.eigenvectors();

    REQUIRE(values.size() == 1);

    // Expected largest eigenvalue
    REQUIRE(values(0) == Approx(10.0));
    REQUIRE(vectors.col(0).norm() == Approx(1.0));
}
TEST_CASE("Lanczos iteration eigenvalue")
{
    neon::lanzcos_solver solver{10};

    solver.solve(create_diagonal_sparse_matrix(50));

    auto const& values = solver.eigenvalues();
    auto const& vectors = solver.eigenvectors();

    REQUIRE(values.size() == 10);

    for (int i = 0; i < 10; i++)
    {
        // Expected eigenvalues
        REQUIRE(values(i) == Approx(50.0 - i));
        // Unit vectors
        REQUIRE(vectors.col(i).norm() == Approx(1.0));
    }
}
