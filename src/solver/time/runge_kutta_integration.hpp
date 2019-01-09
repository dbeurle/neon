
#pragma once

/// @file

#include "exceptions.hpp"

#include <type_traits>
#include <cfenv>

/// \brief Second, fourth and fifth order explicit methods for time integration

namespace neon
{
/// Perform a second order Runge-Kutta step for the given type.
/// The type \p function must accept time and the value.
template <typename T, typename function>
T runge_kutta_second_order(double const t, double const dt, T const& y, function&& f)
{
    T const dy1 = dt * f(t, y);
    T const dy2 = dt * f(t + dt, y + dy1);
    return y + (dy1 + dy2) / 2.0;
}

/// Perform a fourth order Runge-Kutta step for the given type.
/// The type \p function must accept time and the value.
template <typename T, typename function>
T runge_kutta_fourth_order(double const t, double const dt, T const& y, function&& f)
{
    T const dy1 = dt * f(t, y);
    T const dy2 = dt * f(t + dt / 2.0, y + dy1 / 2.0);
    T const dy3 = dt * f(t + dt / 2.0, y + dy2 / 2.0);
    T const dy4 = dt * f(t + dt, y + dy3);

    return y + (dy1 + 2.0 * dy2 + 2.0 * dy3 + dy4) / 6.0;
}

namespace detail
{
/// Abstraction to compute the residual or the absolution value depending on
/// the expression type (scalar or vector).
template <typename left_expression, typename right_expression>
double residual(left_expression const& left, right_expression const& right)
{
    if constexpr (!std::is_floating_point<left_expression>::value
                  && !std::is_floating_point<right_expression>::value)
    {
        return (left - right).norm();
    }
    else
    {
        return std::abs(left - right);
    }
}
}

/// Perform the Dorman-Prince 4th order embedded Runge-Kutta time discretisation
/// using the specified error tolerance @cite DormandPrince1980
/// \tparam Lambda that returns the right hand side of the ODE
/// \param f Function object
/// \param error_tolerance Error tolerance for adaptive method
template <typename T, typename function>
T runge_kutta_fourth_fifth_order(double t,
                                 double dt,
                                 T y,
                                 function&& f,
                                 double const error_tolerance = 1.0e-5)
{
    auto const end_time = t + dt;

    std::feclearexcept(FE_ALL_EXCEPT);

    // Perform substep iterations until convergence
    while (t < end_time)
    {
        // Correct an overzealous time step from the previous iteration
        dt = std::min(dt, end_time - t);

        // clang-format off
        T const k1 = dt * f(t, y);
        T const k2 = dt * f(t + 1.0 / 5.0 * dt,  y +        1.0 / 5.0 * k1);
        T const k3 = dt * f(t + 3.0 / 10.0 * dt, y +       3.0 / 40.0 * k1 +       9.0 / 40.0 * k2);
        T const k4 = dt * f(t + 4.0 / 5.0 * dt,  y +      44.0 / 45.0 * k1 -      56.0 / 15.0 * k2 +       32.0 / 9.0 * k3);
        T const k5 = dt * f(t + 8.0 / 9.0 * dt,  y + 19372.0 / 6561.0 * k1 - 25360.0 / 2187.0 * k2 + 64448.0 / 6561.0 * k3 - 212.0 / 729.0 * k4);
        T const k6 = dt * f(t + dt,              y +  9017.0 / 3168.0 * k1 -     355.0 / 33.0 * k2 + 46732.0 / 5247.0 * k3 +  49.0 / 176.0 * k4 - 5103.0 / 18656.0 * k5);
        T const k7 = dt * f(t + dt,              y +     35.0 / 384.0 * k1                         +   500.0 / 1113.0 * k3 + 125.0 / 192.0 * k4 -  2187.0 / 6784.0 * k5 + 11.0 / 84.0 * k6);

        T const dy_trial =     35.0 / 384.0 * k1 +   500.0 / 1113.0 * k3 + 125.0 / 192.0 * k4 -    2187.0 / 6784.0 * k5 +    11.0 / 84.0 * k6;
        T const dy_error = 5179.0 / 57600.0 * k1 + 7571.0 / 16695.0 * k3 + 393.0 / 640.0 * k4 - 92097.0 / 339200.0 * k5 + 187.0 / 2100.0 * k6 + 1.0 / 40.0 * k7;
        // clang-format on

        double const R = detail::residual(dy_error, dy_trial) / dt;

        if (R <= error_tolerance)
        {
            y += dy_trial;
            t += dt;
        }
        dt *= 0.84
              * std::pow(error_tolerance / (R + 4.0 * std::numeric_limits<double>::epsilon()), 0.25);

        if (std::fetestexcept(FE_INVALID))
        {
            std::feclearexcept(FE_ALL_EXCEPT);
            throw computational_error("Floating point error in adaptive RK45 detected.\n");
        }
    }
    return y;
}
}
