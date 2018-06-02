
#include <catch.hpp>

#include "math/jacobian_determinant.hpp"

using neon::jacobian_determinant;

TEST_CASE("3x3")
{
    neon::matrix3 const random_jacobian = neon::matrix3::Random();
    REQUIRE(random_jacobian.determinant() == Approx(jacobian_determinant(random_jacobian)));
}

TEST_CASE("2x2")
{
    neon::matrix2 const random_jacobian = neon::matrix2::Random();
    REQUIRE(random_jacobian.determinant() == Approx(jacobian_determinant(random_jacobian)));
}

TEST_CASE("3x2")
{
    neon::matrix32 const random_jacobian = neon::matrix32::Random();
    double const j = jacobian_determinant(random_jacobian);
}
