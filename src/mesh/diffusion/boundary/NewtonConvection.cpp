
#include "NewtonConvection.hpp"

#include <range/v3/view/transform.hpp>
#include <range/v3/view/zip.hpp>

#include <json/value.h>

namespace neon::diffusion
{
NewtonConvection::NewtonConvection(std::unique_ptr<SurfaceInterpolation>&& sf,
                                   std::vector<List> const& nodal_connectivity,
                                   std::vector<List> const& dof_list,
                                   std::shared_ptr<MaterialCoordinates>& material_coordinates,
                                   Json::Value const& time_history,
                                   Json::Value const& heat_transfer_coefficients,
                                   Json::Value const& prescribed_heat_fluxes)
    : SurfaceLoad<SurfaceInterpolation>(std::move(sf),
                                        nodal_connectivity,
                                        dof_list,
                                        material_coordinates,
                                        time_history,
                                        prescribed_heat_fluxes)
{
    using namespace ranges;
    coefficient_time_data = view::zip(time_history, heat_transfer_coefficients)
                            | view::transform([](auto const& time_load_pair) {
                                  return std::make_pair(time_load_pair.first.asDouble(),
                                                        time_load_pair.second.asDouble());
                              });
}

std::tuple<List const&, Matrix> NewtonConvection::external_stiffness(int const element,
                                                                     double const load_factor) const
{
    // clang-format off
    auto const X = geometry::project_to_plane(material_coordinates->initial_configuration(nodal_connectivity[element]));

    auto const coefficient = interpolate_prescribed_load(coefficient_time_data, load_factor);

    // Perform the computation of the external element stiffness matrix
    auto const k_ext = sf->quadrature().integrate(Matrix::Zero(X.cols(), X.cols()).eval(),
                                                  [&](auto const& femval, auto const& l) -> Matrix {
                                                      auto const & [ N, dN ] = femval;

                                                      auto const j = (X * dN).determinant();

                                                      return coefficient * N * N.transpose() * j;
                                                  });
    return {dof_list[element], k_ext};
    // clang-format on
}
}
