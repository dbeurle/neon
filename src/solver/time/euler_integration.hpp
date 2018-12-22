
#pragma once

/// @file

#include <type_traits>

namespace neon
{
/// Perform a first order explicit Euler step for the given lambda.
/// The type \p function must accept time and the value.
template <typename T, typename function>
T explicit_euler(double const t, double const dt, T&& y, function&& f)
{
    return y + dt * f(t, y);
}

/// Perform a first order implicit Euler step for the given lambda.
/// The type \p function must accept time and the value.
/// This integration is implicit and it requires a different form of the function
/// as given by:
/// \f$ y(t + \Delta t) = \frac{y(t)}{1 + h(t, y, \Delta t)}\f$
/// otherwise known as the Jacobian.
template <typename T, typename function>
T implicit_euler(double const t, double const dt, T const& y, function&& h)
{
    return y * (1.0 + 1.0 / (1.0 + h(t, y, dt)));
}
}
