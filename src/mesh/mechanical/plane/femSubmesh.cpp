
#include "femSubmesh.hpp"

#include "Exceptions.hpp"

#include "constitutive/ConstitutiveModelFactory.hpp"
#include "geometry/Projection.hpp"
#include "interpolations/interpolation_factory.hpp"
#include "material/Material.hpp"
#include "mesh/DofAllocator.hpp"
#include "mesh/MaterialCoordinates.hpp"

#include "numeric/mechanics"

#include <cfenv>
#include <chrono>

#include <omp.h>

#include <range/v3/algorithm/find_if.hpp>
#include <range/v3/view/transform.hpp>

#include <termcolor/termcolor.hpp>

namespace neon::mechanical::plane
{
femSubmesh::femSubmesh(json const& material_data,
                       json const& simulation_data,
                       std::shared_ptr<MaterialCoordinates>& material_coordinates,
                       Submesh const& submesh)
    : detail::femSubmesh<plane::femSubmesh, plane::InternalVariables>(submesh),
      material_coordinates(material_coordinates),
      sf(make_surface_interpolation(topology(), simulation_data)),
      view(sf->quadrature().points()),
      variables(std::make_shared<InternalVariables>(elements() * sf->quadrature().points())),
      cm(make_constitutive_model(variables, material_data, simulation_data))
{
    // Allocate storage for the displacement gradient
    variables->add(InternalVariables::Tensor::DisplacementGradient,
                   InternalVariables::Tensor::DeformationGradient,
                   InternalVariables::Tensor::Cauchy,
                   InternalVariables::Scalar::DetF);

    // Get the old data to the undeformed configuration
    for (auto& F : variables->fetch(InternalVariables::Tensor::DeformationGradient))
    {
        F = matrix2::Identity();
    }

    variables->commit();

    dof_list = allocate_dof_list(dofs_per_node(), nodal_connectivity);
}

void femSubmesh::save_internal_variables(bool const have_converged)
{
    if (have_converged)
    {
        variables->commit();
    }
    else
    {
        variables->revert();
    }
}

std::pair<local_indices const&, matrix> femSubmesh::tangent_stiffness(int32 const element) const
{
    auto const x = geometry::project_to_plane(
        material_coordinates->current_configuration(local_node_list(element)));

    matrix ke = material_tangent_stiffness(x, element);

    if (!cm->is_finite_deformation()) return {local_dof_view(element), ke};

    ke.noalias() += geometric_tangent_stiffness(x, element);

    return {local_dof_view(element), ke};
}

std::pair<local_indices const&, vector> femSubmesh::internal_force(int32 const element) const
{
    auto const x = geometry::project_to_plane(
        material_coordinates->current_configuration(local_node_list(element)));

    return {local_dof_view(element), internal_nodal_force(x, element)};
}

matrix femSubmesh::geometric_tangent_stiffness(matrix2x const& x, int32 const element) const
{
    auto const& cauchy_stresses = variables->fetch(InternalVariables::Tensor::Cauchy);

    auto const n = nodes_per_element();

    matrix const kgeo = sf->quadrature()
                            .integrate(matrix::Zero(n, n).eval(),
                                       [&](auto const& femval, auto const& l) -> matrix {
                                           auto const& [N, rhea] = femval;

                                           matrix2 const Jacobian = local_deformation_gradient(rhea,
                                                                                               x);

                                           auto const cauchy = cauchy_stresses[view(element, l)];

                                           // Compute the symmetric gradient operator
                                           auto const L = (rhea * Jacobian.inverse()).transpose();

                                           return L.transpose() * cauchy * L * Jacobian.determinant();
                                       });

    return identity_expansion(kgeo, dofs_per_node());
}

matrix femSubmesh::material_tangent_stiffness(matrix2x const& x, int32 const element) const
{
    auto const local_dofs = nodes_per_element() * dofs_per_node();

    auto const& tangent_operators = variables->fetch(InternalVariables::rank4::tangent_operator);

    return sf->quadrature().integrate(matrix::Zero(local_dofs, local_dofs).eval(),
                                      [&](auto const& femval, auto const& l) -> matrix {
                                          auto const& [N, rhea] = femval;

                                          auto const& D = tangent_operators[view(element, l)];

                                          matrix2 const Jacobian{local_deformation_gradient(rhea, x)};

                                          auto const B = fem::sym_gradient<2>(
                                              (rhea * Jacobian.inverse()).transpose());

                                          return B.transpose() * D * B * Jacobian.determinant();
                                      });
}

vector femSubmesh::internal_nodal_force(matrix2x const& x, int32 const element) const
{
    auto const& cauchy_stresses = variables->fetch(InternalVariables::Tensor::Cauchy);

    auto const [m, n] = std::make_tuple(nodes_per_element(), dofs_per_node());

    vector fint = vector::Zero(m * n);

    sf->quadrature().integrate(Eigen::Map<row_matrix>(fint.data(), m, n),
                               [&](auto const& femval, auto const& l) -> row_matrix {
                                   auto const& [N, dN] = femval;

                                   matrix2 const Jacobian = local_deformation_gradient(dN, x);

                                   auto const& cauchy_stress = cauchy_stresses[view(element, l)];

                                   // Compute the symmetric gradient operator
                                   auto const Bt = dN * Jacobian.inverse();

                                   return Bt * cauchy_stress * Jacobian.determinant();
                               });
    return fint;
}

std::pair<local_indices const&, matrix> femSubmesh::consistent_mass(int32 const element) const
{
    auto const X = material_coordinates->initial_configuration(local_node_list(element));
    return {local_dof_view(element), X};
    // auto const density_0 = cm->intrinsic_material().initial_density();
    //
    // auto m = sf->quadrature().integrate(matrix::Zero(nodes_per_element(), nodes_per_element()).eval(),
    //                                     [&](auto const& femval, auto const& l) -> matrix {
    //                                         auto const & [N, dN] = femval;
    //
    //                                         auto const Jacobian = local_deformation_gradient(dN, X);
    //
    //                                         return N * density_0 * N.transpose()
    //                                                * Jacobian.determinant();
    //                                     });
    // return {local_dof_view(element), identity_expansion(m, dofs_per_node())};
}

std::pair<local_indices const&, vector> femSubmesh::diagonal_mass(int32 const element) const
{
    auto const& [dofs, consistent_m] = this->consistent_mass(element);

    vector diagonal_m(consistent_m.rows());
    for (auto i = 0; i < consistent_m.rows(); ++i)
    {
        diagonal_m(i) = consistent_m.row(i).sum();
    }
    return {local_dof_view(element), diagonal_m};
}

void femSubmesh::update_internal_variables(double const time_step_size)
{
    std::feclearexcept(FE_ALL_EXCEPT);

    update_deformation_measures();

    update_Jacobian_determinants();

    cm->update_internal_variables(time_step_size);

    if (std::fetestexcept(FE_INVALID))
    {
        throw computational_error("Floating point error reported\n");
    }
}

void femSubmesh::update_deformation_measures()
{
    auto& H_list = variables->fetch(InternalVariables::Tensor::DisplacementGradient);
    auto& F_list = variables->fetch(InternalVariables::Tensor::DeformationGradient);

    // #pragma omp parallel for
    for (auto element = 0; element < elements(); ++element)
    {
        // Gather the material coordinates
        auto const X = geometry::project_to_plane(
            material_coordinates->initial_configuration(local_node_list(element)));
        auto const x = geometry::project_to_plane(
            material_coordinates->current_configuration(local_node_list(element)));

        sf->quadrature().for_each([&](auto const& femval, const auto& l) {
            auto const& [N, rhea] = femval;

            // Local deformation gradient for the initial configuration
            matrix2 const F_0 = local_deformation_gradient(rhea, X);
            matrix2 const F = local_deformation_gradient(rhea, x);

            // Gradient operator in index notation
            auto const& B_0t = rhea * F_0.inverse();

            // Displacement gradient
            matrix2 const H = (x - X) * B_0t;

            H_list[view(element, l)] = H;
            F_list[view(element, l)] = F * F_0.inverse();
        });
    }
}

void femSubmesh::update_Jacobian_determinants()
{
    auto const& deformation_gradients = variables->fetch(
        InternalVariables::Tensor::DeformationGradient);
    auto& deformation_gradient_determinants = variables->fetch(InternalVariables::Scalar::DetF);

    deformation_gradient_determinants = deformation_gradients
                                        | ranges::view::transform(
                                              [](auto const& F) { return F.determinant(); });

    auto const found = ranges::find_if(deformation_gradient_determinants,
                                       [](auto const& detF) { return detF <= 0.0; });

    if (found != deformation_gradient_determinants.end())
    {
        auto const i = std::distance(deformation_gradient_determinants.begin(), found);

        auto const [element, quadrature_point] = std::div(i, sf->quadrature().points());

        throw computational_error("Positive Jacobian assumption violated at element "
                                  + std::to_string(element) + " and local quadrature point "
                                  + std::to_string(quadrature_point) + " (" + std::to_string(*found)
                                  + ")");
    }
}

std::pair<vector, vector> femSubmesh::nodal_averaged_variable(
    InternalVariables::Tensor const tensor_name) const
{
    vector count = vector::Zero(material_coordinates->size() * 4);
    vector value = count;

    auto const& tensor_list = variables->fetch(tensor_name);

    auto const& E = sf->local_quadrature_extrapolation();

    // vector format of values
    vector component = vector::Zero(sf->quadrature().points());

    for (auto element = 0; element < elements(); ++element)
    {
        // Assemble these into the global value vector
        auto const& node_list = local_node_list(element);

        for (auto ci = 0; ci < 2; ++ci)
        {
            for (auto cj = 0; cj < 2; ++cj)
            {
                for (auto l = 0; l < sf->quadrature().points(); ++l)
                {
                    component(l) = tensor_list[view(element, l)](ci, cj);
                }

                // Local extrapolation to the nodes
                vector const nodal_component = E * component;

                for (auto n = 0; n < nodal_component.rows(); n++)
                {
                    value(node_list[n] * 4 + ci * 2 + cj) += nodal_component(n);
                    count(node_list[n] * 4 + ci * 2 + cj) += 1.0;
                }
            }
        }
    }
    return {value, count};
}

std::pair<vector, vector> femSubmesh::nodal_averaged_variable(
    InternalVariables::Scalar const scalar_name) const
{
    vector count = vector::Zero(material_coordinates->size());
    vector value = count;

    auto const& scalar_list = variables->fetch(scalar_name);

    auto const& E = sf->local_quadrature_extrapolation();

    // vector format of values
    vector component = vector::Zero(sf->quadrature().points());

    for (auto element = 0; element < elements(); ++element)
    {
        for (auto l = 0; l < sf->quadrature().points(); ++l)
        {
            component(l) = scalar_list[view(element, l)];
        }

        // Local extrapolation to the nodes
        vector const nodal_component = E * component;

        // Assemble these into the global value vector
        auto const& node_list = local_node_list(element);

        for (auto n = 0; n < nodal_component.rows(); n++)
        {
            value(node_list[n]) += nodal_component(n);
            count(node_list[n]) += 1.0;
        }
    }
    return {value, count};
}
}
