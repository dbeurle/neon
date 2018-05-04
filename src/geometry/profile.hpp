
#pragma once

#include <cmath>
#include <utility>

namespace neon::geometry
{
/// profile is the base class for a geometry profile of a beam or shell.  This
/// is designed to provide the routines for the accessing the profile
/// information such as second moment of area.
class profile
{
public:
    explicit profile(
        double const I_1, double const I_2, double const A, double const A_1, double const A_2);

    std::pair<double, double> second_moment_area() const noexcept;

    /// \return First and second shear area
    std::pair<double, double> shear_area() const noexcept;

    /// \return Cross-section area
    double area() const noexcept;

protected:
    /// Area moment of inertia (first coordinate)
    double I_1;
    /// Area moment of inertia (second coordinate)
    double I_2;

    /// Section area
    double m_area;

    /// Shear area (first coordinate)
    double A_1;
    /// Shear area (first coordinate)
    double A_2;
};

/// rectangular_bar computes the second moment of area for a rectangular bar
/// where a width and a height is specified.
class rectangular_bar : public profile
{
public:
    explicit rectangular_bar(double const width, double const height);
};

/// hollow_rectangular_bar computes the second moment of area for a rectangular bar
/// with a hollow internal structure where a width and a height is specified.
class hollow_rectangular_bar : public profile
{
public:
    hollow_rectangular_bar(double const width, double const height);
};

class section
{
};
}
