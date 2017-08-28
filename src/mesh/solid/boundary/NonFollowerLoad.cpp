
#include "NonFollowerLoad.hpp"

#include "interpolations/InterpolationFactory.hpp"
#include "mesh/DofAllocator.hpp"

namespace neon::solid
{
NonFollowerLoad::NonFollowerLoad(std::vector<List> const& nodal_connectivity,
                                 std::shared_ptr<MaterialCoordinates>& material_coordinates,
                                 double const prescribed_load,
                                 bool const is_load_ramped,
                                 int const dof_offset,
                                 int const nodal_dofs)
    : Neumann(nodal_connectivity,
              filter_dof_list(nodal_dofs, dof_offset, nodal_connectivity),
              prescribed_load,
              is_load_ramped),
      material_coordinates(material_coordinates)
{
}

std::tuple<List const&, Vector> Traction::external_force(int const element,
                                                         double const load_factor) const
{
    auto X = material_coordinates->initial_configuration(nodal_connectivity.at(element));
    X = sf->project_to_plane(X);

    // Perform the computation of the external load vector
    auto const f_ext = sf->quadrature().integrate(Vector::Zero(X.cols()).eval(),
                                                  [&](auto const& femval, auto const& l) {
                                                      auto const & [ N, dN ] = femval;

                                                      auto const j = (X * dN).determinant();

                                                      return interpolate_prescribed_value(load_factor)
                                                             * N * j;
                                                  });
    return {dof_list.at(element), f_ext};
}

NonFollowerLoadBoundary::NonFollowerLoadBoundary(std::shared_ptr<MaterialCoordinates>& material_coordinates,
                                                 std::vector<SubMesh> const& submeshes,
                                                 int const dof_offset,
                                                 double const prescribed_load,
                                                 Json::Value const& simulation_data)
{
    // Populate the entire mesh
    for (auto const& mesh : submeshes)
    {
        nf_loads.emplace_back(make_surface_interpolation(mesh.topology(), simulation_data),
                              mesh.connectivities(),
                              material_coordinates,
                              prescribed_load,
                              true,
                              dof_offset);
    }
}
}
