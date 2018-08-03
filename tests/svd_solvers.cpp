
#include <catch.hpp>

#include "numeric/dense_matrix.hpp"
#include "numeric/sparse_matrix.hpp"
#include "solver/svd/svd.hpp"

#include "io/json.hpp"
#include <iostream>
#include <cmath>
#include <chrono>

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

/** Create an example right hand side */
vector create_right_hand_side()
{
    Eigen::VectorXd b(3);
    b(0) = 1.0;
    b(1) = 0.5;
    b(2) = 1.5;
    return b;
}

/** The known solution */
vector solution()
{
    vector x(3);
    x << 1.375, 1.75, 1.625;
    return x;
}

TEST_CASE("svd solver test suite")
{
    sparse_matrix A = create_sparse_matrix();
    vector b = create_right_hand_side();

    matrix B(4, 8);
    // clang-format off
    B << 1, 0, 0, 1, 0, 0, 0, 0,
    1, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 1,
    0, 1, 0, 0, 0, 0, 1, 0;
    // clang-format on

    matrix m = matrix::Random(1e5, 50);

    SECTION("bdc svd: singular values and singular vectors")
    {
        bdc_svd svd_decomposition;
        svd_decomposition.compute(B);

        vector values(4);
        values << 2, 1.414213, 1.414213, 0;
        vector left(4);
        vector left_different_order(4);

        left << -std::sqrt(0.5), -std::sqrt(0.5), 0, 0;
        left_different_order << -std::sqrt(0.5), 0, -std::sqrt(0.5), 0;
        REQUIRE((values - svd_decomposition.values()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((svd_decomposition.left().col(0) - left).norm() == Approx(0.0).margin(ZERO_MARGIN));

        REQUIRE(svd_decomposition.values().size() == 4);

        svd_decomposition.compute(B, 5);
        REQUIRE(svd_decomposition.values().size() == 4);

        svd_decomposition.compute(B, 2);
        REQUIRE(svd_decomposition.values().size() == 2);

        svd_decomposition.compute(B, 1e-1);
        REQUIRE(svd_decomposition.values().size() == 3);
    }

    SECTION("randomised svd: singular values and singular vectors")
    {
        randomised_svd svd_decomposition;
        svd_decomposition.compute(B);

        vector values(4);
        values << 2, 1.414213, 1.414213, 0;
        vector left(4);
        vector left_different_order(4);

        left << -std::sqrt(0.5), -std::sqrt(0.5), 0, 0;
        left_different_order << -std::sqrt(0.5), 0, -std::sqrt(0.5), 0;
        REQUIRE((values - svd_decomposition.values()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((svd_decomposition.left().col(0) - left).norm() == Approx(0.0).margin(ZERO_MARGIN));

        REQUIRE(svd_decomposition.values().size() == 4);

        svd_decomposition.compute(B, 5);
        REQUIRE(svd_decomposition.values().size() == 4);

        svd_decomposition.compute(B, 2);
        REQUIRE(svd_decomposition.values().size() == 2);

        svd_decomposition.compute(B, 1e-1);
        REQUIRE(svd_decomposition.values().size() == 3);
    }

    SECTION("least-squares for a determinant system")
    {
        vector x;
        bdc_svd svd_decomposition(A);
        svd_decomposition.solve(x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }

    SECTION("svd timing")
    {
        auto const start = std::chrono::steady_clock::now();

        bdc_svd svd_decomposition(m);

        matrix m_reconstructed = svd_decomposition.left() * svd_decomposition.values().asDiagonal()
                                 * svd_decomposition.right().transpose();

        matrix approximation_error = m - m_reconstructed;
        std::cout << "approximation error due to bdc_svd " << approximation_error.norm() << "\n";

        auto const end = std::chrono::steady_clock::now();
        std::chrono::duration<double> const elapsed_seconds = end - start;

        std::cout << "BDC SVD test case took " << elapsed_seconds.count() << "s\n";

        // randomised_svd
        auto const start_rand = std::chrono::steady_clock::now();

        randomised_svd randomised_svd_decomposition;
        randomised_svd_decomposition.compute(m, 5);

        m_reconstructed = randomised_svd_decomposition.left()
                          * randomised_svd_decomposition.values().asDiagonal()
                          * randomised_svd_decomposition.right().transpose();

        approximation_error = m - m_reconstructed;
        std::cout << "approximation error due to randomised_svd " << approximation_error.norm()
                  << "\n";

        auto const end_rand = std::chrono::steady_clock::now();
        std::chrono::duration<double> const elapsed_seconds_rand = end_rand - start_rand;

        std::cout << "Randomised SVD test case took " << elapsed_seconds_rand.count() << "s\n";
    }
}
