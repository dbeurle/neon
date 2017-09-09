//
// #include "SurfaceFlux.hpp"
//
// #include "interpolations/InterpolationFactory.hpp"
//
// namespace neon::diffusion
// {
// SurfaceFlux::SurfaceFlux(std::unique_ptr<SurfaceInterpolation>&& sf,
//                          std::shared_ptr<NodalCoordinates>& nodal_coordinates,
//                          std::vector<List> nodal_connectivities,
//                          double const prescribed_value,
//                          bool const is_load_ramped)
//     : Neumann(nodal_connectivities, nodal_connectivities, prescribed_value, is_load_ramped),
//       nodal_coordinates(nodal_coordinates),
//       sf(std::move(sf))
// {
// }
//
// std::tuple<List const&, Vector> SurfaceFlux::external_flux(int const element) const
// {
//     auto X = nodal_coordinates->coordinates(nodal_connectivity.at(element));
//     X = sf->project_to_plane(X);
//
//     // Perform the computation of the external load vector
//     auto const f_e = sf->quadrature().integrate(Vector::Zero(X.cols()).eval(),
//                                                 [&](auto const& femval, auto const& l) {
//                                                     auto const & [ N, dN ] = femval;
//
//                                                     auto const j = (X * dN).determinant();
//
//                                                     return interpolate_prescribed_value(1.0) * N
//                                                     * j;
//                                                 });
//     return {dof_list.at(element), f_e};
// }
//
// SurfaceFluxBoundary::SurfaceFluxBoundary(std::shared_ptr<NodalCoordinates>& nodal_coordinates,
//                                          std::vector<SubMesh> const& submeshes,
//                                          int const dof_offset,
//                                          double const prescribed_value,
//                                          Json::Value const& simulation_data)
// {
//     // Populate the entire mesh
//     for (auto const& mesh : submeshes)
//     {
//         surface_fluxes.emplace_back(make_surface_interpolation(mesh.topology(), simulation_data),
//                                     nodal_coordinates,
//                                     mesh.connectivities(),
//                                     prescribed_value,
//                                     false);
//     }
// }
// }
