
#include <catch2/catch.hpp>

#include "solver/time/runge_kutta_integration.hpp"
#include "solver/time/euler_integration.hpp"

#include "numeric/dense_matrix.hpp"
#include "numeric/interval.hpp"

#include <cmath>

constexpr auto ZERO_MARGIN = 1.0e-5;

using neon::explicit_euler;
using neon::runge_kutta_fourth_fifth_order;
using neon::runge_kutta_fourth_order;
using neon::runge_kutta_second_order;

TEST_CASE("interval")
{
    REQUIRE(neon::interval(0.0, 5.0, 1.0) == 5);
    REQUIRE(neon::interval(0.0, 1.0, 0.1) == 10);
}

TEST_CASE("dy/dt = 1 integration")
{
    REQUIRE(explicit_euler(0.0, 10.0, 0.0, [](auto, auto) { return 1.0; }) == Approx(10.0));

    REQUIRE(runge_kutta_second_order(0.0, 10.0, 0.0, [](auto, auto) { return 1.0; }) == Approx(10.0));

    REQUIRE(runge_kutta_fourth_order(0.0, 10.0, 0.0, [](auto, auto) { return 1.0; }) == Approx(10.0));
}
TEST_CASE("dy/dt = y integration")
{
    // Exact solution should be std::exp(1.0)
    REQUIRE(explicit_euler(0.0, 1.0, 1.0, [](double, double const y) { return y; }) == Approx(2.0));

    REQUIRE(runge_kutta_second_order(0.0, 1.0, 1.0, [](double, double const y) { return y; })
            == Approx(2.5));

    REQUIRE(runge_kutta_fourth_order(0.0, 1.0, 1.0, [](double, double const y) { return y; })
            == Approx(2.70833));

    REQUIRE(runge_kutta_fourth_fifth_order(0.0, 1.0, 1.0, [](double, double const y) { return y; })
            == Approx(std::exp(1.0)));
}
TEST_CASE("dy/dt = y integration with system")
{
    using neon::vector4;
    REQUIRE((runge_kutta_fourth_fifth_order(0.0,
                                            1.0,
                                            vector4::Ones().eval(),
                                            [](double, vector4 y) { return y; })
             - std::exp(1.0) * vector4::Ones())
                .norm()
            == Approx(0.0).margin(ZERO_MARGIN));
}
