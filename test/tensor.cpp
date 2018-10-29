
#include <catch2/catch.hpp>

#include "numeric/tensor_operations.hpp"
#include "numeric/spectral_decomposition.hpp"
#include "numeric/log_tensor_derivative.hpp"

using namespace neon;

constexpr auto ZERO_MARGIN = 1.0e-5;

TEST_CASE("deviatoric and volumetric decomposition")
{
    matrix3 t = matrix3::Random();

    auto const t_dev = deviatoric(t);
    auto const t_vol = volumetric(t);

    REQUIRE(((t_dev + t_vol) - t).norm() == Approx(0.0));
    REQUIRE(t_dev.trace() == Approx(0.0).margin(ZERO_MARGIN));
    REQUIRE(t_vol.trace() == Approx(t.trace()));
}
TEST_CASE("Mandel notation")
{
    SECTION("3x3 matrix of ones")
    {
        matrix3 const mandel_ones = mandel_notation(matrix3::Ones());

        REQUIRE(mandel_ones(0, 0) == Approx(1.0).margin(ZERO_MARGIN));
        REQUIRE(mandel_ones(1, 1) == Approx(1.0).margin(ZERO_MARGIN));
        REQUIRE(mandel_ones(0, 1) == Approx(1.0).margin(ZERO_MARGIN));
        REQUIRE(mandel_ones(1, 0) == Approx(1.0).margin(ZERO_MARGIN));

        REQUIRE(mandel_ones(0, 2) == Approx(std::sqrt(2.0)).margin(ZERO_MARGIN));
        REQUIRE(mandel_ones(1, 2) == Approx(std::sqrt(2.0)).margin(ZERO_MARGIN));

        REQUIRE(mandel_ones(2, 0) == Approx(std::sqrt(2.0)).margin(ZERO_MARGIN));
        REQUIRE(mandel_ones(2, 1) == Approx(std::sqrt(2.0)).margin(ZERO_MARGIN));

        REQUIRE(mandel_ones(2, 2) == Approx(2.0).margin(ZERO_MARGIN));
    }
    SECTION("6x6 matrix of ones")
    {
        matrix6 const mandel_ones = mandel_notation(matrix6::Ones());

        REQUIRE(mandel_ones.block<3, 3>(3, 3).sum() == Approx(18.0).margin(ZERO_MARGIN));
        REQUIRE(mandel_ones.block<3, 3>(0, 0).sum() == Approx(9.0).margin(ZERO_MARGIN));
        REQUIRE(mandel_ones.block<3, 3>(0, 3).sum()
                == Approx(9.0 * std::sqrt(2.0)).margin(ZERO_MARGIN));
        REQUIRE(mandel_ones.block<3, 3>(3, 0).sum()
                == Approx(9.0 * std::sqrt(2.0)).margin(ZERO_MARGIN));
    }
}
TEST_CASE("Spectral decomposition")
{
    SECTION("2x2 identity matrix")
    {
        matrix2 A = matrix2::Identity();

        auto const [is_unique, eigenvalues, eigenprojections] = spectral_decomposition(A);

        auto const [x1, x2] = eigenvalues;

        REQUIRE(is_unique == false);
        REQUIRE(x1 == Approx(1.0).margin(ZERO_MARGIN));
        REQUIRE(x2 == Approx(1.0).margin(ZERO_MARGIN));
    }
    SECTION("2x2 unity matrix")
    {
        matrix2 A = matrix2::Ones();

        auto const [is_unique, eigenvalues, eigenprojections] = spectral_decomposition(A);

        auto const [x1, x2] = eigenvalues;

        REQUIRE(is_unique == true);
        REQUIRE(x1 == Approx(2.0).margin(ZERO_MARGIN));
        REQUIRE(x2 == Approx(0.0).margin(ZERO_MARGIN));
    }
}
TEST_CASE("Log symmetric tensor derivative")
{
    SECTION("Unique eigenvalues")
    {
        matrix2 A(2, 2);
        A << 2.0, 0.0, 0.0, 2.0;

        matrix3 const dA = log_symmetric_tensor_derivative(A);

        REQUIRE(dA(0, 0) == Approx(0.5).margin(ZERO_MARGIN));
        REQUIRE(dA(1, 1) == Approx(0.5).margin(ZERO_MARGIN));
        REQUIRE(dA(2, 2) == Approx(0.25).margin(ZERO_MARGIN));
    }
    SECTION("Zero eigenvalue")
    {
        matrix3 const dA = log_symmetric_tensor_derivative(matrix2::Ones());
        REQUIRE(dA.hasNaN() == true);
    }
    SECTION("Repeated eigenvalue of one")
    {
        matrix3 const dA = log_symmetric_tensor_derivative(matrix2::Identity());
        REQUIRE(dA(0, 0) == Approx(1.0).margin(ZERO_MARGIN));
        REQUIRE(dA(1, 1) == Approx(1.0).margin(ZERO_MARGIN));
        REQUIRE(dA(2, 2) == Approx(0.5).margin(ZERO_MARGIN));
    }
}
