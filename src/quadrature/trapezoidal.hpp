
#pragma once

namespace neon
{
template <class T, typename function>
T trapezoidal(T const lower_bound, T const upper_bound, T const step_size, function&& f)
{
    static_assert(std::is_floating_point<T>::value, "Integration of floating point only");

    T partial_sum = static_cast<T>(0.5) * (f(lower_bound) + f(upper_bound));

    std::int64_t const steps = (upper_bound - lower_bound) / step_size + 1;

    for (std::int64_t step{1}; step < steps - 1; ++step)
    {
        partial_sum += f(lower_bound + step * step_size);
    }
    return partial_sum * step_size;
}

template <class T>
T partial_trapezoidal(T const left_value, T const right_value, T const step_size)
{
    return (left_value + right_value) / 2.0 * step_size;
}
}
