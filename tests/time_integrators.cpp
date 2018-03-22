
#include <catch.hpp>

#include "solver/time/runge_kutta_integration.hpp"
#include "solver/time/euler_integration.hpp"

#include <cmath>

constexpr auto ZERO_MARGIN = 1.0e-5;

template <class T, typename integrator>
T generic_integrate(T y, T const end_time, T const step_size, integrator&& integrate)
{
    auto const time_steps = static_cast<int>(end_time / step_size) + 1;

    auto t = 0.0;

    for (auto i = 0; i < time_steps - 1; ++i)
    {
        y += integrate(t, y, step_size);
        t += step_size;
    }
    return y;
}

TEST_CASE("dy/dt = 1 integration")
{
    REQUIRE(generic_integrate(0.0, 10.0, 1.0, neon::explicit_euler([](auto const t, auto const y) {
                                  return 1.0;
                              }))
            == Approx(10.0).margin(ZERO_MARGIN));

    REQUIRE(
        generic_integrate(0.0, 10.0, 1.0, neon::runge_kutta_second_order([](auto const t, auto const y) {
                              return 1.0;
                          }))
        == Approx(10.0).margin(ZERO_MARGIN));

    REQUIRE(
        generic_integrate(0.0, 10.0, 1.0, neon::runge_kutta_fourth_order([](auto const t, auto const y) {
                              return 1.0;
                          }))
        == Approx(10.0).margin(ZERO_MARGIN));
}

TEST_CASE("dy/dt = 2t integration")
{
    REQUIRE(generic_integrate(0.0, 5.0, 0.1, neon::explicit_euler([](auto const t, auto const y) {
                                  return 2.0 * t;
                              }))
            == Approx(24.5).margin(ZERO_MARGIN));

    REQUIRE(
        generic_integrate(0.0, 5.0, 0.1, neon::runge_kutta_second_order([](auto const t, auto const y) {
                              return 2.0 * t;
                          }))
        == Approx(25.0).margin(ZERO_MARGIN));

    REQUIRE(
        generic_integrate(0.0, 5.0, 0.1, neon::runge_kutta_fourth_order([](auto const t, auto const y) {
                              return 2.0 * t;
                          }))
        == Approx(25.0).margin(ZERO_MARGIN));

    REQUIRE(generic_integrate(0.0,
                              5.0,
                              0.1,
                              neon::runge_kutta_fourth_fifth_order(
                                  [](auto const t, auto const y) { return 2.0 * t; }))
            == Approx(25.0).margin(ZERO_MARGIN));
}
