
#include "material_coordinates.hpp"

namespace neon
{
material_coordinates::material_coordinates(matrix3x const& initial_coordinates)
    : nodal_coordinates(initial_coordinates), x(initial_coordinates)
{
}

vector material_coordinates::displacement() const
{
    matrix3x u = x - X;
    return Eigen::Map<vector>(u.data(), u.size());
}

void material_coordinates::update_current_xy_configuration(vector const& u)
{
    x.row(0) = X.row(0) + u.transpose()(Eigen::seq(0, u.size() - 1, 2));
    x.row(1) = X.row(1) + u.transpose()(Eigen::seq(1, u.size() - 1, 2));
}

void material_coordinates::update_current_configuration(vector const& u)
{
    std::copy_n(u.data(), u.size(), x.data());
    x.noalias() += X;
}
}
