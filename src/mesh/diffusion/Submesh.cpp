//
// #include "Submesh.hpp"
//
// #include "Exceptions.hpp"
//
// #include "constitutive/ConstitutiveModelFactory.hpp"
// #include "interpolations/InterpolationFactory.hpp"
//
// #include "material/Material.hpp"
// #include "mesh/DofAllocator.hpp"
// #include "mesh/NodalCoordinates.hpp"
// #include "numeric/Operators.hpp"
//
// #include <cfenv>
// #include <chrono>
// #include <omp.h>
//
// #include <json/json.h>
// #include <range/v3/view.hpp>
// #include <termcolor/termcolor.hpp>
//
// namespace neon::diffusion
// {
// femSubmesh::femSubmesh(Json::Value const& material_data,
//                        Json::Value const& simulation_data,
//                        std::shared_ptr<NodalCoordinates>& nodal_coordinates,
//                        SubMesh const& submesh)
//     : neon::SubMesh(submesh),
//       nodal_coordinates(nodal_coordinates),
//       sf(solid::make_volume_interpolation(topology(), simulation_data)),
//       variables(elements() * sf->quadrature().points()),
//       cm(make_constitutive_model(variables, material_data, simulation_data))
// {
// }
//
// void femSubmesh::save_internal_variables(bool const have_converged)
// {
//     if (have_converged)
//     {
//         variables.commit();
//     }
//     else
//     {
//         variables.revert();
//     }
// }
//
// std::tuple<List const&, Matrix> femSubmesh::tangent_stiffness(int const element) const
// {
//     auto const X = nodal_coordinates->coordinates(local_node_list(element));
//
//     auto const n = nodes_per_element();
//
//     auto const& D_Vec = variables(InternalVariables::Matrix::TruesdellModuli);
//
//     Matrix const kmat = sf->quadrature()
//                             .integrate(Matrix::Zero(n, n).eval(),
//                                        [&](auto const& femval, auto const& l) -> Matrix {
//
//                                            auto const & [ N, rhea ] = femval;
//
//                                            auto const& D = D_Vec[offset(element, l)];
//
//                                            Matrix3 const Jacobian = local_jacobian(rhea, X);
//
//                                            // Compute the symmetric gradient operator
//                                            Matrix const B = (rhea *
//                                            Jacobian.inverse()).transpose();
//
//                                            return B.transpose() * D * B * Jacobian.determinant();
//                                        });
//
//     return {local_dof_list(element), kmat};
// }
//
// std::tuple<List const&, Matrix> femSubmesh::consistent_mass(int const element) const
// {
//     auto X = nodal_coordinates->coordinates(local_node_list(element));
//
//     auto const density_0 = cm->intrinsic_material().initial_density();
//
//     auto m = sf->quadrature().integrate(Matrix::Zero(nodes_per_element(),
//     nodes_per_element()).eval(),
//                                         [&](auto const& femval, auto const& l) -> Matrix {
//                                             auto const & [ N, dN ] = femval;
//
//                                             auto const Jacobian = local_jacobian(dN, X);
//
//                                             return N * density_0 * N.transpose()
//                                                    * Jacobian.determinant();
//                                         });
//     return {local_dof_list(element), m};
// }
//
// std::tuple<List const&, Vector> femSubmesh::diagonal_mass(int const element) const
// {
//     auto const & [ dofs, consistent_m ] = this->consistent_mass(element);
//
//     Vector diagonal_m(consistent_m.rows());
//     for (auto i = 0; i < consistent_m.rows(); ++i)
//     {
//         diagonal_m(i) = consistent_m.row(i).sum();
//     }
//     return {local_dof_list(element), diagonal_m};
// }
//
// void femSubmesh::update_internal_variables(double const time_step_size)
// {
//     std::feclearexcept(FE_ALL_EXCEPT);
//
//     cm->update_internal_variables(time_step_size);
//
//     if (std::fetestexcept(FE_INVALID))
//     {
//         throw computational_error("Floating point error reported\n");
//     }
// }
//
// int femSubmesh::offset(int const element, int const quadrature_point) const
// {
//     return sf->quadrature().points() * element + quadrature_point;
// }
//
// femSubmesh::ValueCount femSubmesh::nodal_averaged_variable(InternalVariables::Scalar const
// scalar_name) const
// {
//     Vector count = Vector::Zero(nodal_coordinates->size());
//     Vector value = count;
//
//     auto const& scalar_list = variables(scalar_name);
//
//     auto const& E = sf->local_quadrature_extrapolation();
//
//     // Vector format of values
//     Vector component = Vector::Zero(sf->quadrature().points());
//
//     for (auto e = 0; e < elements(); ++e)
//     {
//         // Assemble these into the global value vector
//         auto const& node_list = local_node_list(e);
//
//         for (auto l = 0; l < sf->quadrature().points(); ++l)
//         {
//             component(l) = scalar_list[this->offset(e, l)];
//         }
//
//         // Local extrapolation to the nodes
//         Vector const nodal_component = E * component;
//
//         for (auto n = 0; n < nodal_component.rows(); n++)
//         {
//             value(node_list[n]) += nodal_component(n);
//             count(node_list[n]) += 1.0;
//         }
//     }
//     return {value, count};
// }
// }
