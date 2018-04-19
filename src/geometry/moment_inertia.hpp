
#pragma once

#include "io/json_forward.hpp"

#include <cmath>

namespace neon::geometry
{
class section
{
public:
    section(double const I_1, double const I_2) : I_1{I_1}, I_2{I_2} {}

    /// \return pair of first and second moment of inertia
    auto second_moment_area() const noexcept { return {I1, I2}; }

protected:
    double I_1{0.0}; /// Area moment of inertia (first coordinate)
    double I_2{0.0}; /// Area moment of inertia (second coordinate)
};

/// rectangular_bar computes the second moment of area for a rectangular bar
/// where a width and a height is specified.
class rectangular_bar : public section
{
public:
    rectangular_bar(double const width, double const height)
        : section(width * std::pow(height, 3) / 12.0, std::pow(width, 3) * height / 12.0)
    {
    }
};
}
