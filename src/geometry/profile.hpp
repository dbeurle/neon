
#pragma once

#include "io/json_forward.hpp"

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
    std::pair<double, double> second_moment_area() const noexcept;

    /// \return First and second shear area
    std::pair<double, double> shear_area() const noexcept;

    /// \return Cross-section area
    double area() const noexcept;

protected:
    /// Area moment of inertia (first coordinate)
    double I_x{0.0};
    /// Area moment of inertia (second coordinate)
    double I_y{0.0};
    /// Mixed moment of area about first and second coordinate
    double I_xy{0.0};

    /// Section area
    double A{0.0};

    /// Shear area (first coordinate)
    double A_x{0.0};
    /// Shear area (second coordinate)
    double A_y{0.0};

    /// Rotational polar moment of Area
    double J{0.0};
};

/// rectangle computes the second moment of area for a rectangular bar
/// for a specified width and height.
class rectangle : public profile
{
public:
    explicit rectangle(json const& section_data);
};

class circle : public profile
{
public:
    explicit circle(json const& section_data);
};

/// hollow_rectangle computes the second moment of area for a rectangular bar
/// with a hollow internal structure where a width and a height is specified.
// class hollow_rectangle : public profile
// {
// public:
//     hollow_rectangle(double const width, double const height);
// };

class section
{
};
}
