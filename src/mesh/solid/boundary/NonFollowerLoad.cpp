
#include "NonFollowerLoad.hpp"

#include "mesh/DofAllocator.hpp"

namespace neon::solid
{
NonFollowerLoad::NonFollowerLoad(std::vector<List> const& nodal_connectivity,
                                 std::shared_ptr<MaterialCoordinates>& material_coordinates,
                                 double const prescribed_load,
                                 bool const is_load_ramped,
                                 int const dof_offset,
                                 int const nodal_dofs)
    : Boundary(prescribed_load, is_load_ramped),
      nodal_connectivity(nodal_connectivity),
      dof_list(filter_dof_list(nodal_dofs, dof_offset, nodal_connectivity)),
      material_coordinates(material_coordinates)
{
}

Traction::Traction(std::vector<List> const& nodal_connectivity,
                   std::shared_ptr<MaterialCoordinates>&& material_coordinates,
                   double const prescribed_load,
                   bool const is_load_ramped,
                   int const dof_offset,
                   std::unique_ptr<SurfaceInterpolation>&& sf)
    : NonFollowerLoad(nodal_connectivity,
                      material_coordinates,
                      prescribed_load,
                      is_load_ramped,
                      dof_offset),
      sf(std::move(sf))
{
}
std::tuple<List const&, Vector> Traction::external_force(int const element,
                                                         double const load_factor) const
{
    auto X = material_coordinates->initial_configuration(nodal_connectivity.at(element));
    X = sf->project_to_plane(X);

    // Perform the computation of the external load vector
    auto const f_ext = sf->quadrature().integrate(Vector::Zero(X.cols()).eval(),
                                                  [&, this](auto const& femval,
                                                            auto const& l) {
                                                      auto const & [ N, dN ] = femval;

                                                      auto const j = (X * dN).determinant();

                                                      return interpolate_prescribed_value(
                                                                 load_factor)
                                                             * N * j;
                                                  });
    return {dof_list.at(element), f_ext};
}
}
