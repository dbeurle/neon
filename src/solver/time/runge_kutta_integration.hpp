
#pragma once

#include <type_traits>

namespace neon
{
/**
 * Perform a second order Runge-Kutta step for the given type.
 * The type \p functor must accept time and the value.
 */
template <typename functor>
auto runge_kutta_second_order(functor&& f)
{
    return [f](auto const t, auto const y, auto const dt) {
        static_assert(std::is_floating_point<decltype(t)>::value);
        static_assert(std::is_floating_point<decltype(dt)>::value);

        auto const dy1 = dt * f(t, y);
        auto const dy2 = dt * f(t + dt, y + dy1);
        return (dy1 + dy2) / 2.0;
    };
}

/**
 * Perform a fourth order Runge-Kutta step for the given type.
 * The type \p functor must accept time and the value.
 */
template <typename functor>
auto runge_kutta_fourth_order(functor&& f)
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
}
