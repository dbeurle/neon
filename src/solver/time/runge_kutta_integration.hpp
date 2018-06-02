
#pragma once

#include <type_traits>

/// \file runge_kutta_integration.hpp
/// \brief Second, fourth and fifth order explicit methods for time integration

namespace neon
{
/// Perform a second order Runge-Kutta step for the given type.
/// The type \p function must accept time and the value.
template <typename function>
auto runge_kutta_second_order(function&& f)
{
    return [f](auto const t, auto const y, auto const dt) {
        static_assert(std::is_floating_point<decltype(t)>::value);
        static_assert(std::is_floating_point<decltype(dt)>::value);

        auto const dy1 = dt * f(t, y);
        auto const dy2 = dt * f(t + dt, y + dy1);
        return (dy1 + dy2) / 2.0;
    };
}

/// Perform a fourth order Runge-Kutta step for the given type.
/// The type \p function must accept time and the value.
template <typename function>
auto runge_kutta_fourth_order(function&& f)
{
    return [f](auto const t, auto const y, auto const dt) {
        static_assert(std::is_floating_point<decltype(t)>::value);
        static_assert(std::is_floating_point<decltype(dt)>::value);

        auto const dy1 = dt * f(t, y);
        auto const dy2 = dt * f(t + dt / 2.0, y + dy1 / 2.0);
        auto const dy3 = dt * f(t + dt / 2.0, y + dy2 / 2.0);
        auto const dy4 = dt * f(t + dt, y + dy3);
        return (dy1 + 2.0 * dy2 + 2.0 * dy3 + dy4) / 6.0;
    };
}

/// Perform the Dorman-Prince 4th order embedded Runge-Kutta time discretisation
/// \cite DormandPrince1980
template <typename function>
auto runge_kutta_fourth_fifth_order(function&& f, double const error_tolerance = 1.0e-5)
{
    return [f, error_tolerance](auto t, auto const y0, auto dt) {
        static_assert(std::is_floating_point<decltype(t)>::value);
        static_assert(std::is_floating_point<decltype(dt)>::value);

        auto y = y0;

        auto const end_time = t + dt;

        // Perform substep iterations until convergence
        while (t < end_time)
        {
            // Correct an overzealous time step from the previous iteration
            dt = std::min(dt, end_time - t);

            // clang-format off
            auto const k1 = dt * f(t, y);
            auto const k2 = dt * f(t + 1.0 / 5.0 * dt, y + 1.0 / 5.0 * k1);
            auto const k3 = dt * f(t + 3.0 / 10.0 * dt, y + 3.0 / 40.0 * k1 + 9.0 / 40.0 * k2);
            auto const k4 = dt * f(t + 4.0 / 5.0 * dt, y + 44.0 / 45.0 * k1 - 56.0 / 15.0 * k2 + 32.0 / 9.0 * k3);
            auto const k5 = dt * f(t + 8.0 / 9.0 * dt, y + 19372.0 / 6561.0 * k1 - 25360.0 / 2187 * k2 + 64448.0 / 2187.0 * k3 - 212.0 / 729.0 * k4);
            auto const k6 = dt * f(t + dt, y + 9017.0 / 3168.0 * k1 - 355.0 / 33.0 * k2 - 46732.0 / 5247.0 * k3 + 49.0 / 176.0 * k4 - 5103.0 / 18656.0 * k5);
            auto const k7 = dt * f(t + dt, y + 35.0 / 384.0 * k1 + 500.0 / 1113.0 * k3 + 125.0 / 192.0 * k4 - 2187.0 / 6784.0 * k5 + 11.0 / 84.0 * k6);

            auto const dy_trial = 35.0 / 384.0 * k1 + 500.0 / 1113.0 * k3 + 125.0 / 192.0 * k4 - 2187.0 / 6784.0 * k5 + 11.0 / 84.0 * k6;
            auto const dy_error = 5179.0 / 57600.0 * k1 + 7571.0 / 16695.0 * k3 + 393.0 / 640.0 * k4 - 92097.0 / 339200.0 * k5 + 187.0 / 2100.0 * k6 + 1.0 / 40.0 * k7;
            // clang-format on

            auto const R = std::abs(dy_error - dy_trial) / dt;

            if (R <= error_tolerance)
            {
                y += dy_trial;
                t += dt;
            }
            dt *= 0.84
                  * std::pow(error_tolerance / (R + std::numeric_limits<double>::epsilon()), 0.25);
        }
        return y - y0;
    };
}
}
