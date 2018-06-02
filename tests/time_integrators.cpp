
#include <catch.hpp>

#include "solver/time/runge_kutta_integration.hpp"
#include "solver/time/euler_integration.hpp"

#include "numeric/dense_matrix.hpp"

#include <cmath>

constexpr auto ZERO_MARGIN = 1.0e-5;

template <typename T, typename U, typename integrator>
T generic_integrate(T y, U const end_time, U const step_size, integrator&& integrate)
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
TEST_CASE("dy/dt = 2t integration with system")
{
    using neon::vector4;
    vector4 y;
    REQUIRE((generic_integrate(y.setZero(),
                               5.0,
                               0.1,
                               neon::runge_kutta_fourth_fifth_order([](auto const t, auto const y) {
                                   return 2.0 * t * vector4::Ones();
                               }))
             - 25.0 * vector4::Ones())
                .norm()
            == Approx(0.0).margin(ZERO_MARGIN));
}
