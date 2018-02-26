
#define CATCH_CONFIG_MAIN

#include <catch.hpp>

#include "solver/eigen/eigen_solver.hpp"
#include "io/json.hpp"

using namespace neon;

constexpr auto ZERO_MARGIN = 1.0e-5;

/** Create a SPD matrix for solver testing */
sparse_matrix create_sparse_matrix()
{
    std::vector<Eigen::Triplet<double>> triplets = {{0, 0, 2.0},
                                                    {0, 1, -1.0},
                                                    {0, 2, 0.0},
                                                    {1, 0, -1.0},
                                                    {1, 1, 2.0},
                                                    {1, 2, -1.0},
                                                    {2, 0, 0.0},
                                                    {2, 1, -1.0},
                                                    {2, 2, 2.0}};

    sparse_matrix A(3, 3);
    A.setFromTriplets(std::begin(triplets), std::end(triplets));
    A.finalize();
    return A;
}

sparse_matrix create_sparse_identity()
{
    std::vector<Eigen::Triplet<double>> triplets = {{0, 0, 1.0}, {1, 1, 1.0}, {2, 2, 1.0}};

    sparse_matrix A(3, 3);
    A.setFromTriplets(std::begin(triplets), std::end(triplets));
    A.finalize();
    return A;
}

TEST_CASE("Eigen solver test suite")
{
    sparse_matrix A = create_sparse_matrix();
    sparse_matrix I = create_sparse_identity();

    std::cout << A << std::endl;
    std::cout << I << std::endl;

    eigen_solver eigen{3};

    auto const [values, vectors] = eigen.solve(A, I);

    std::cout << values << std::endl;
    std::cout << vectors << std::endl;

    REQUIRE(1 == 1);
}
