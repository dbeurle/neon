
#include <catch.hpp>

#include "quadrature/trapezoidal.hpp"

TEST_CASE("trapezoidal method")
{
    using neon::trapezoidal;

    REQUIRE(trapezoidal(0.0, 1.0, 1.0, [](auto const& i) { return i; }) == Approx(0.5));
    REQUIRE(trapezoidal(0.0, 1.0, 1.0e-6, [](auto const& i) { return i * i; }) == Approx(1.0 / 3.0));
    REQUIRE(trapezoidal(0.0, 1.0, 1.0e-6, [](auto const& i) { return i * i * i; })
            == Approx(1.0 / 4.0));
}
TEST_CASE("partial trapezoidal method")
{
    using neon::partial_trapezoidal;

    REQUIRE(partial_trapezoidal(0.0, 1.0, 1.0) == Approx(0.5));
    REQUIRE(partial_trapezoidal(0.0, 1.0, 0.5) == Approx(0.25));
}
