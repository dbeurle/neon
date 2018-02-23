
#include "pressure.hpp"

#include "geometry/Projection.hpp"

#include <Eigen/Geometry>

namespace neon::mechanical::solid
{
std::pair<index_view, vector> pressure::external_force(std::int64_t const element,
                                                       double const load_factor) const
{
    auto const X = coordinates->initial_configuration(
        nodal_connectivity(Eigen::placeholders::all, element));

    auto const X_surface = geometry::project_to_plane(X);

    auto const pressure = interpolate_prescribed_load(load_factor);

    // Perform the computation of the external load vector
    matrix f_ext = -pressure
                   * sf->quadrature().integrate(matrix::Zero(X.cols(), 3).eval(),
                                                [&](auto const& femval, auto const& l) -> matrix {
                                                    auto const& [N, dN] = femval;

                                                    auto const j = (X_surface * dN).determinant();

                                                    vector3 const x_xi = (X * dN).col(0);
                                                    vector3 const x_eta = (X * dN).col(1);

                                                    vector3 const normal = x_xi.cross(x_eta).normalized();

                                                    return N * normal.transpose() * j;
                                                });

    // Map the matrix back to a vector for the assembly operator
    return {dof_list(Eigen::placeholders::all, element),
            Eigen::Map<matrix>(f_ext.data(), X.cols() * 3, 1)};
}
}
