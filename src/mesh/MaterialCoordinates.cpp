
#include "MaterialCoordinates.hpp"

namespace neon
{
MaterialCoordinates::MaterialCoordinates(matrix3x const& initial_coordinates)
    : NodalCoordinates(initial_coordinates), x(initial_coordinates)
{
}

matrix3x MaterialCoordinates::displacement() const { return x - X; }

void MaterialCoordinates::update_current_xy_configuration(vector const& u)
{
    x.row(0) = X.row(0) + u.transpose()(Eigen::seq(0, u.size() - 1, 2));
    x.row(1) = X.row(1) + u.transpose()(Eigen::seq(1, u.size() - 1, 2));
}

void MaterialCoordinates::update_current_configuration(vector const& u)
{
    x.row(0) = X.row(0) + u.transpose()(Eigen::seq(0, u.size() - 1, 3));
    x.row(1) = X.row(1) + u.transpose()(Eigen::seq(1, u.size() - 1, 3));
    x.row(2) = X.row(2) + u.transpose()(Eigen::seq(2, u.size() - 1, 3));
}
}
