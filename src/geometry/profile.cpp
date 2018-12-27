
#include "profile.hpp"
#include "io/json.hpp"

namespace neon::geometry
{
std::pair<double, double> profile::second_moment_area() const noexcept { return {I_x, I_y}; }

std::pair<double, double> profile::shear_area() const noexcept { return {A_x, A_y}; }

double profile::area() const noexcept { return A; }

rectangle::rectangle(json const& section_data)
{
    if (section_data.find("width") == section_data.end())
    {
        throw std::domain_error("\"width\" was not specified for \"rectangle\" in section");
    }
    if (section_data.find("height") == section_data.end())
    {
        throw std::domain_error("\"height\" was not specified for \"rectangle\" in section");
    }
    double const width = section_data["width"];

    double const height = section_data["height"];

    if (width <= 0.0)
    {
        throw std::domain_error("\"width\" must have a positive value for \"rectangle\" in "
                                "section");
    }

    if (height <= 0.0)
    {
        throw std::domain_error("\"height\" must have a positive value for \"rectangle\" in "
                                "profile");
    }

    I_x = width * std::pow(height, 3) / 12.0;
    I_y = std::pow(width, 3) * height / 12.0;
    A = A_x = A_y = width * height;
}

circle::circle(json const& section_data)
{
    if (section_data.find("diameter") == section_data.end())
    {
        throw std::domain_error("\"diameter\" was not specified for \"circle\" in profile");
    }

    double const diameter = section_data["diameter"];

    if (diameter <= 0.0)
    {
        throw std::domain_error("\"diameter\" must have a positive value for \"circle\" in "
                                "profile");
    }

    auto const pi = std::acos(-1.0);

    I_x = I_y = pi * std::pow(diameter, 4) / 64.0;
    A = A_x = A_y = pi * std::pow(diameter, 2) / 4.0;
    J = I_x * 2.0;
}

hollow_circle::hollow_circle(json const& section_data)
{
    if (section_data.find("outer_diameter") == section_data.end())
    {
        throw std::domain_error("\"outer_diameter\" was not specified for \"hollow_circle\" in "
                                "profile");
    }
    if (section_data.find("inner_diameter") == section_data.end())
    {
        throw std::domain_error("\"inner_diameter\" was not specified for \"hollow_circle\" in "
                                "profile");
    }

    double const outer_diameter = section_data["outer_diameter"];
    double const inner_diameter = section_data["inner_diameter"];

    if (outer_diameter <= 0.0)
    {
        throw std::domain_error("\"outer_diameter\" must have a positive value for "
                                "\"hollow_circle\" in profile");
    }
    if (inner_diameter <= 0.0)
    {
        throw std::domain_error("\"inner_diameter\" must have a positive value for "
                                "\"hollow_circle\" in profile");
    }

    if (inner_diameter >= outer_diameter)
    {
        throw std::domain_error("\"inner_diameter\" value must be less than \"outer_diameter\" "
                                "value for \"hollow_circle\" in profile");
    }

    auto const pi = std::acos(-1.0);

    A = A_x = A_y = pi * (std::pow(outer_diameter, 2) - std::pow(inner_diameter, 2)) / 4.0;
    I_x = I_y = pi * (std::pow(outer_diameter, 4) - std::pow(inner_diameter, 4)) / 64.0;
    I_xy = 0.0;
    J = I_x * 2.0;
}

// hollow_rectangle::hollow_rectangle(double const width, double const height)
//     : profile(width * std::pow(height, 3) / 12.0,
//               std::pow(width, 3) * height / 12.0,
//               width * height,
//               width * height,
//               width * height)
// {
// }

rectangular_angle::rectangular_angle(json const& section_data)
{
    if (section_data.find("width") == section_data.end())
    {
        throw std::domain_error("\"width\" was not specified for \"rectangular_angle\" in profile");
    }
    if (section_data.find("height") == section_data.end())
    {
        throw std::domain_error("\"height\" was not specified for \"rectangular_angle\" in "
                                "profile");
    }
    if (section_data.find("thickness") == section_data.end())
    {
        throw std::domain_error("\"thickness\" was not specified for \"rectangular_angle\" in "
                                "profile");
    }

    double const width = section_data["width"];

    double const height = section_data["height"];

    double const thickness = section_data["thickness"];

    if (width <= 0.0)
    {
        throw std::domain_error("\"width\" must have a positive value for \"rectangular_angle\" in "
                                "profile");
    }

    if (height <= 0.0)
    {
        throw std::domain_error("\"height\" must have a positive value for \"rectangular_angle\" "
                                "in "
                                "profile");
    }

    if (thickness <= 0.0)
    {
        throw std::domain_error("\"thickness\" must have a positive value for "
                                "\"rectangular_angle\" in "
                                "profile");
    }

    // Formulas sourced from: http://structx.com/Shape_Formulas_008.html
    double const c = width - thickness;
    double const d = height - thickness;
    double const C_x = (thickness * (2 * c + height) + std::pow(c, 2)) / (2 * (c + height));
    double const C_y = (thickness * (2 * d + width) + std::pow(d, 2)) / (2 * (d + width));
    double const y = height - C_y;
    double const z = width - C_x;
    I_x = (thickness * std::pow(y, 3) + width * std::pow(height - y, 3)
           - (width - thickness) * std::pow(height - y - thickness, 3))
          / 3.0;

    I_y = (thickness * std::pow(z, 3) + height * std::pow(width - z, 3)
           - (height - thickness) * std::pow(width - z - thickness, 3))
          / 3.0;

    A = thickness * (width + d);
}
}
