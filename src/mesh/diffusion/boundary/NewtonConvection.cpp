
#include "NewtonConvection.hpp"

#include <range/v3/view/transform.hpp>
#include <range/v3/view/zip.hpp>

#include <json/value.h>

namespace neon::diffusion
{
newton_cooling::newton_cooling(std::unique_ptr<SurfaceInterpolation>&& sf,
                               std::vector<std::vector<int64>> const& nodal_connectivity,
                               std::vector<std::vector<int64>> const& dof_list,
                               std::shared_ptr<MaterialCoordinates>& material_coordinates,
                               Json::Value const& times,
                               Json::Value const& heat_flux,
                               Json::Value const& heat_transfer_coefficient)
    : boundary::surface_load<SurfaceInterpolation>(std::move(sf),
                                                   nodal_connectivity,
                                                   dof_list,
                                                   material_coordinates,
                                                   times,
                                                   heat_flux)
{
    using ranges::view::transform;
    using ranges::view::zip;

    stiffness_time_data = zip(times | transform([](auto i) { return i.asDouble(); }),
                              heat_transfer_coefficient
                                  | transform([](auto i) { return i.asDouble(); }));
}

std::pair<std::vector<int64> const&, matrix> newton_cooling::external_stiffness(
    int64 const element, double const load_factor) const
{
    auto const X = geometry::project_to_plane(
        material_coordinates->initial_configuration(nodal_connectivity[element]));

    // Perform the computation of the external element stiffness matrix
    auto const k_ext = sf->quadrature().integrate(matrix::Zero(X.cols(), X.cols()).eval(),
                                                  [&](auto const& femval, auto const& l) -> matrix {
                                                      auto const& [N, dN] = femval;

                                                      auto const j = (X * dN).determinant();

                                                      return N * N.transpose() * j;
                                                  });

    return {dof_list[element], interpolate_prescribed_load(stiffness_time_data, load_factor) * k_ext};
}
}
