
#include "newton_convection.hpp"

#include "math/jacobian_determinant.hpp"
#include "io/json.hpp"

#include <range/v3/view/transform.hpp>
#include <range/v3/view/zip.hpp>

#include <utility>

namespace neon::diffusion::boundary
{
newton_convection::newton_convection(std::unique_ptr<surface_interpolation>&& sf,
                                     indices node_indices,
                                     indices dof_indices,
                                     std::shared_ptr<material_coordinates>& coordinates,
                                     json const& times,
                                     json const& heat_flux,
                                     json const& heat_transfer_coefficient)
    : surface_load<surface_interpolation>(std::move(sf),
                                          std::move(node_indices),
                                          std::move(dof_indices),
                                          coordinates,
                                          times,
                                          heat_flux)
{
    using ranges::view::transform;
    using ranges::view::zip;

    stiffness_time_data = zip(times | transform([](auto const i) { return i; }),
                              heat_transfer_coefficient | transform([](auto const i) { return i; }));
}

std::pair<index_view, matrix> newton_convection::external_stiffness(std::int64_t const element,
                                                                    double const load_factor) const
{
    matrix3x const& X = coordinates->initial_configuration(local_node_view(element));

    // Perform the computation of the external element stiffness matrix
    auto const k_ext = sf->quadrature().integrate(matrix::Zero(X.cols(), X.cols()).eval(),
                                                  [&](auto const& femval, auto) -> matrix {
                                                      auto const& [N, dN] = femval;

                                                      auto const j = jacobian_determinant(X * dN);

                                                      return N * N.transpose() * j;
                                                  });

    return {local_dof_view(element),
            interpolate_prescribed_load(stiffness_time_data, load_factor) * k_ext};
}
}
