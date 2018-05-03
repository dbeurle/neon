
#pragma once

#include <cmath>

namespace neon::geometry
{
/// profile is the base class for a geometry profile of a beam or shell.  This
/// is designed to provide the routines for the accessing the profile
/// information such as second moment of area.
class profile
{
public:
    explicit profile(
        double const I_1, double const I_2, double const A, double const A_1, double const A_2)
        : I_1{I_1}, I_2{I_2}, m_area{A}, A_1{A_1}, A_2{A_2}
    {
    }

    /// \return First and second moment of inertia
    auto second_moment_area() const noexcept { return std::pair{I_1, I_2}; }

    /// \return First and second shear area
    auto shear_area() const noexcept { return std::pair{A_1, A_2}; }

    /// \return Cross-section area
    auto area() const noexcept { return m_area; }

protected:
    /// Area moment of inertia (first coordinate)
    double I_1{0.0};
    /// Area moment of inertia (second coordinate)
    double I_2{0.0};

    /// Section area
    double m_area;
    /// Shear area (first coordinate)
    double A_1{0.0};
    /// Shear area (first coordinate)
    double A_2{0.0};
};

/// rectangular_bar computes the second moment of area for a rectangular bar
/// where a width and a height is specified.
class rectangular_bar : public profile
{
public:
    explicit rectangular_bar(double const width, double const height)
        : profile(width * std::pow(height, 3) / 12.0,
                  std::pow(width, 3) * height / 12.0,
                  width * height,
                  width * height,
                  width * height)
    {
    }
};

/// hollow_rectangular_bar computes the second moment of area for a rectangular bar
/// with a hollow internal structure where a width and a height is specified.
class hollow_rectangular_bar : public profile
{
public:
    hollow_rectangular_bar(double const width, double const height)
        : profile(width * std::pow(height, 3) / 12.0,
                  std::pow(width, 3) * height / 12.0,
                  width * height,
                  width * height,
                  width * height)
    {
    }
};

class section
{
};
}
