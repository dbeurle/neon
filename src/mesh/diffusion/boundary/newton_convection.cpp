
#include "newton_convection.hpp"

#include "math/jacobian_determinant.hpp"

#include <range/v3/view/transform.hpp>
#include <range/v3/view/zip.hpp>

#include "io/json.hpp"

namespace neon::diffusion::boundary
{
newton_convection::newton_convection(std::unique_ptr<surface_interpolation>&& sf,
                                     indices node_indices,
                                     indices dof_indices,
                                     std::shared_ptr<material_coordinates>& mesh_coordinates,
                                     json const& times,
                                     json const& heat_flux,
                                     json const& heat_transfer_coefficient)
    : surface_load<surface_interpolation>(std::move(sf),
                                          node_indices,
                                          dof_indices,
                                          mesh_coordinates,
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
    auto const node_view = node_indices(Eigen::placeholders::all, element);

    auto const X = coordinates->initial_configuration(node_view);

    // Perform the computation of the external element stiffness matrix
    auto const k_ext = sf->quadrature().integrate(matrix::Zero(X.cols(), X.cols()).eval(),
                                                  [&](auto const& femval, auto const& l) -> matrix {
                                                      auto const& [N, dN] = femval;

                                                      auto const j = jacobian_determinant(X * dN);

                                                      return N * N.transpose() * j;
                                                  });

    return {dof_indices(Eigen::placeholders::all, element),
            interpolate_prescribed_load(stiffness_time_data, load_factor) * k_ext};
}
}
