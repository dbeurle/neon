
#include "SurfaceFlux.hpp"

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
    : SurfaceFlux(std::move(sf),
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
}
