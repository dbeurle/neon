
#pragma once

#include <type_traits>

/// \file euler_integration.hpp

namespace neon
{
/// Perform a first order explicit Euler step for the given lambda.
/// The type \p functor must accept time and the value.
template <typename functor>
auto explicit_euler(functor&& f)
{
    return [f](auto const t, auto const y, auto const dt) {
        static_assert(std::is_floating_point<decltype(t)>::value);
        static_assert(std::is_floating_point<decltype(dt)>::value);

        return dt * f(t, y);
    };
}

/// Perform a first order implicit Euler step for the given lambda.
/// The type \p functor must accept time and the value.
/// This integration is implicit and it requires a different form of the functor
/// as given by:
/// \f$ y(t + \Delta t) = \frac{y(t)}{1 + h(t, y, \Delta t)}\f$
/// otherwise known as the Jacobian.
template <typename functor>
auto implicit_euler(functor&& h)
{
    return [h](auto const t, auto const y, auto const dt) {
        static_assert(std::is_floating_point<decltype(t)>::value);
        static_assert(std::is_floating_point<decltype(dt)>::value);

        return y / (1.0 + h(t, y, dt));
    };
}
}
