
#pragma once

#include <cstdint>

/// \file trapezoidal.hpp

namespace neon
{
/// Numerically integrate the function \p f between the \p lower_bound and
/// \p upper_bound with a given \p step_size.  This uses a simple trapezoidal
/// method without adaptive step sizes.
/// \tparam T Floating point type
/// \tparam function Generic lambda or function object to integrate
/// \param lower_bound Lower integral bound
/// \param upper_bound Upper integral bound
/// \param step_size Discretisation size
/// \param f Function to evaluate
/// \return The numerically evaluated integral
template <class T, typename function>
T trapezoidal(T const lower_bound, T const upper_bound, T const step_size, function&& f)
{
    static_assert(std::is_floating_point<T>::value, "Integration of floating point only");

    T partial_sum = static_cast<T>(0.5) * (f(lower_bound) + f(upper_bound));

    std::int64_t const steps = (upper_bound - lower_bound) / step_size;

    for (std::int64_t step{1}; step < steps; ++step)
    {
        partial_sum += f(lower_bound + step * step_size);
    }
    return partial_sum * step_size;
}

/// Perform one step of the trapezoidal method using the function evaluated
/// on the left hand side \p left_value and the right hand side \p right_value
/// for a given step size \p step_size
/// \tparam T Floating point type
/// \param left_value Function evaluated on the left of the interval
/// \param right_value Function evaluated on the right of the interval
/// \param step_size Step size of the interval
/// \return Numerically integrate value for the interval
template <class T>
T partial_trapezoidal(T const left_value, T const right_value, T const step_size)
{
    static_assert(std::is_floating_point<T>::value, "Integration of floating point only");

    return (left_value + right_value) / 2.0 * step_size;
}
}
