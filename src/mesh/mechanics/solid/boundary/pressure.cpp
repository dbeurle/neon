
#include "pressure.hpp"

#include "math/jacobian_determinant.hpp"

namespace neon::mechanics::solid
{
auto pressure::external_force(std::int64_t const element, double const load_factor) const
    -> vector const&
{
    thread_local vector f_ext;

    auto const node_view = local_node_view(element);

    auto const X = coordinates->initial_configuration(node_view);

    auto const pressure = interpolate_prescribed_load(load_factor);

    f_ext = vector::Zero(X.cols() * 3);

    // Perform the computation of the external load vector
    linear_form.integrate(Eigen::Map<row_matrix>(f_ext.data(), X.cols(), 3),
                          [&](auto const& value, auto) -> row_matrix {
                              auto const& [N, dN] = value;

                              matrix32 const jacobian = X * dN;

                              auto const j = jacobian_determinant(jacobian);

                              vector3 const dx_dxi = jacobian.col(0);
                              vector3 const dx_deta = jacobian.col(1);

                              vector3 normal = dx_dxi.cross(dx_deta).normalized();

                              return N * normal.transpose() * j;
                          });
    f_ext *= -pressure;

    return f_ext;
}
}
