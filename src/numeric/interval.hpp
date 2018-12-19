
#pragma once

#include <cstdint>
#include <stdexcept>

/// \file interval.hpp

namespace neon
{
/// \return Number of intervals between \p start and \p end with a \p step_size as a 64 bit integer
[[nodiscard]] inline std::int64_t interval(double const start,
                                           double const end,
                                           double const step_size) noexcept(false)
{
    if (end <= start)
    {
        throw std::domain_error("The end time must be greater than the start time");
    }
    return (end - start) / step_size;
}
}
