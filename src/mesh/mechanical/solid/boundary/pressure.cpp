
#include "pressure.hpp"

#include "geometry/Projection.hpp"

#include <Eigen/Geometry>

namespace neon::mechanical::solid
{
std::tuple<List const&, Vector> pressure::external_force(int const element,
                                                         double const load_factor) const
{
    auto const X = material_coordinates->initial_configuration(nodal_connectivity[element]);

    auto const X_surface = geometry::project_to_plane(X);

    auto const pressure = interpolate_prescribed_load(load_factor);

    // Perform the computation of the external load vector
    RowMatrix f_ext = -pressure
                      * sf->quadrature()
                            .integrate(RowMatrix::Zero(X.cols(), 3).eval(),
                                       [&](auto const& femval, auto const& l) -> RowMatrix {
                                           auto const& [N, dN] = femval;

                                           auto const j = (X_surface * dN).determinant();

                                           Vector3 const x_xi = (X * dN).col(0);
                                           Vector3 const x_eta = (X * dN).col(1);

                                           Vector3 const normal = x_xi.cross(x_eta).normalized();

                                           return N * normal.transpose() * j;
                                       });

    // Map the matrix back to a vector for the assembly operator
    return {dof_list[element], Eigen::Map<RowMatrix>(f_ext.data(), X.cols() * 3, 1)};
}
}
