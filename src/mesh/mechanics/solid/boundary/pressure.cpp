
#include "pressure.hpp"

#include "geometry/projection.hpp"

#include <Eigen/Geometry>

namespace neon::mechanics::solid
{
std::pair<index_view, vector> pressure::external_force(std::int64_t const element,
                                                       double const load_factor) const
{
    auto const node_view = node_indices(Eigen::all, element);

    matrix3x const& X = coordinates->initial_configuration(node_view);

    auto const pressure = interpolate_prescribed_load(load_factor);

    // Perform the computation of the external load vector
    matrix f_ext = -pressure
                   * sf->quadrature().integrate(matrix::Zero(X.cols(), 3).eval(),
                                                [&](auto const& femval, auto) -> matrix {
                                                    auto const& [N, dN] = femval;

                                                    matrix32 const jacobian = X * dN;

                                                    auto const j = jacobian_determinant(jacobian);

                                                    vector3 dx_dxi = jacobian.col(0);
                                                    vector3 dx_deta = jacobian.col(1);

                                                    vector3 normal = dx_dxi.cross(dx_deta).normalized();

                                                    return N * normal.transpose() * j;
                                                });

    // Map the matrix back to a vector for the assembly operator
    return {dof_indices(Eigen::all, element), Eigen::Map<matrix>(f_ext.data(), X.cols() * 3, 1)};
}
}
