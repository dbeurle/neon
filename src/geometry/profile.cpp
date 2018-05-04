
#include "profile.hpp"

namespace neon::geometry
{
profile::profile(double const I_1, double const I_2, double const A, double const A_1, double const A_2)
    : I_1{I_1}, I_2{I_2}, m_area{A}, A_1{A_1}, A_2{A_2}
{
}

std::pair<double, double> profile::second_moment_area() const noexcept { return {I_1, I_2}; }

std::pair<double, double> profile::shear_area() const noexcept { return {A_1, A_2}; }

double profile::area() const noexcept { return m_area; }

rectangular_bar::rectangular_bar(double const width, double const height)
    : profile(width * std::pow(height, 3) / 12.0,
              std::pow(width, 3) * height / 12.0,
              width * height,
              width * height,
              width * height)
{
}

hollow_rectangular_bar::hollow_rectangular_bar(double const width, double const height)
    : profile(width * std::pow(height, 3) / 12.0,
              std::pow(width, 3) * height / 12.0,
              width * height,
              width * height,
              width * height)
{
}
}
