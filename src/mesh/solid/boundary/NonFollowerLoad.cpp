
#include "NonFollowerLoad.hpp"

namespace neon::solid
{
Pressure::Pressure(std::unique_ptr<SurfaceInterpolation>&& sf,
                   std::vector<List> const& nodal_connectivity,
                   std::shared_ptr<MaterialCoordinates>& material_coordinates,
                   Json::Value const& time_history,
                   Json::Value const& load_history,
                   int const nodal_dofs)
    : Traction(std::move(sf),
               nodal_connectivity,
               material_coordinates,
               time_history,
               load_history,
               0,
               nodal_dofs)
{
    dof_list = allocate_dof_list(nodal_dofs, nodal_connectivity);
}

std::tuple<List const&, Vector> Pressure::external_force(int const element,
                                                         double const load_factor) const
{
    Matrix const X = material_coordinates->initial_configuration(nodal_connectivity[element]);
    Matrix const X_proj = sf->project_to_plane(X);

    // Constant pressure
    auto const p = interpolate_prescribed_load(load_factor);

    Vector f_p(X.cols() * 3);

    // Perform the computation of the external load vector which is now contains
    // three degrees of freedom for each Cartesian direction
    // auto const f_ext = sf->quadrature().integrate(Vector::Zero(X.cols() * 3).eval(),
    //                                               [&](auto const& femval, auto const& l) {
    //                                                   auto const & [ N, dN ] = femval;
    //
    //                                                   Vector3 const x_xi = (X * dN).col(0);
    //                                                   Vector3 const x_eta = (X * dN).col(1);
    //
    //                                                   auto const j = (X_proj * dN).determinant();
    //
    //                                                   Vector3 const n = x_xi.cross(x_eta);
    //
    //                                                   for (auto i = 0, a = 0; a < N.rows();
    //                                                        i += 3, a++)
    //                                                   {
    //                                                       f_p.segment<3>(i) = N(a) * n.norm();
    //                                                   }
    //
    //                                                   return f_p;
    //   });
    return {dof_list[element], -p * f_ext};
}

NonFollowerLoadBoundary::NonFollowerLoadBoundary(std::shared_ptr<MaterialCoordinates>& material_coordinates,
                                                 std::vector<SubMesh> const& submeshes,
                                                 Json::Value const& times,
                                                 Json::Value const& loads,
                                                 Json::Value const& simulation_data,
                                                 int const dof_offset)
{
    for (auto const& mesh : submeshes)
    {
        surface_loads.emplace_back(make_surface_interpolation(mesh.topology(), simulation_data),
                                   mesh.connectivities(),
                                   material_coordinates,
                                   times,
                                   loads,
                                   dof_offset,
                                   3);
    }
}
}
