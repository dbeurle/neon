
#include "femSubmesh.hpp"

#include "Exceptions.hpp"

#include "constitutive/ConstitutiveModelFactory.hpp"
#include "interpolations/InterpolationFactory.hpp"
#include "material/Material.hpp"
#include "mesh/DofAllocator.hpp"
#include "mesh/MaterialCoordinates.hpp"
#include "numeric/Operators.hpp"
#include "numeric/mechanics"

#include <cfenv>
#include <chrono>
#include <omp.h>

#include <range/v3/algorithm/count_if.hpp>
#include <range/v3/algorithm/fill.hpp>
#include <range/v3/algorithm/find_if.hpp>
#include <range/v3/view/transform.hpp>

#include <termcolor/termcolor.hpp>

namespace neon::mechanical::solid
{
femSubmesh::femSubmesh(Json::Value const& material_data,
                       Json::Value const& mesh_data,
                       std::shared_ptr<MaterialCoordinates>& material_coordinates,
                       Submesh const& submesh)
    : Submesh(submesh),
      material_coordinates(material_coordinates),
      sf(make_volume_interpolation(topology(), mesh_data)),
      variables(elements() * sf->quadrature().points()),
      view(sf->quadrature().points()),
      cm(make_constitutive_model(variables, material_data, mesh_data))
{
    // Allocate storage for the displacement gradient
    variables.add(InternalVariables::Tensor::DisplacementGradient,
                  InternalVariables::Tensor::DeformationGradient,
                  InternalVariables::Tensor::Cauchy);

    variables.add(InternalVariables::Scalar::DetF);

    // Get the old data to the undeformed configuration
    ranges::fill(variables(InternalVariables::Tensor::DeformationGradient), Matrix3::Identity());

    variables.commit();

    dof_list = allocate_dof_list(this->dofs_per_node(), nodal_connectivity);
}

void femSubmesh::save_internal_variables(bool const have_converged)
{
    if (have_converged)
    {
        variables.commit();
    }
    else
    {
        variables.revert();
    }
}

std::tuple<List const&, Matrix> femSubmesh::tangent_stiffness(int const element) const
{
    auto const x = material_coordinates->current_configuration(local_node_list(element));

    Matrix ke = material_tangent_stiffness(x, element);

    if (!cm->is_finite_deformation()) return {local_dof_list(element), ke};

    ke.noalias() += geometric_tangent_stiffness(x, element);

    return {local_dof_list(element), ke};
}

std::tuple<List const&, Vector> femSubmesh::internal_force(int const element) const
{
    auto x = material_coordinates->current_configuration(local_node_list(element));

    return {local_dof_list(element), internal_nodal_force(x, element)};
}

Matrix femSubmesh::geometric_tangent_stiffness(Matrix const& x, int const element) const
{
    auto const& cauchy_stresses = variables(InternalVariables::Tensor::Cauchy);

    auto n = nodes_per_element();

    // clang-format off
    Matrix const kgeo = sf->quadrature().integrate(Matrix::Zero(n, n).eval(),
                                                   [&](auto const& femval, auto const& l) -> Matrix {
                                                        auto const & [ N, rhea ] = femval;

                                                        auto const Jacobian = local_deformation_gradient(rhea, x);

                                                        auto const cauchy_stress = cauchy_stresses[view(element, l)];

                                                        // Compute the symmetric gradient operator
                                                        auto const L = (rhea * Jacobian.inverse()).transpose();

                                                        return L.transpose() * cauchy_stress * L * Jacobian.determinant();
                                                    });
    // clang-format on
    return identity_expansion(kgeo, dofs_per_node());
}

Matrix femSubmesh::material_tangent_stiffness(Matrix const& x, int const element) const
{
    auto const local_dofs = nodes_per_element() * dofs_per_node();

    auto const& tangent_operators = variables(InternalVariables::Matrix::TangentOperator);

    Matrix kmat = Matrix::Zero(local_dofs, local_dofs);

    return sf->quadrature().integrate(kmat, [&](auto const& femval, auto const& l) -> Matrix {
        auto const& [N, rhea] = femval;

        auto const& D = tangent_operators[view(element, l)];

        auto const Jacobian = local_deformation_gradient(rhea, x);

        // Compute the symmetric gradient operator
        Matrix const B = fem::sym_gradient<3>((rhea * Jacobian.inverse()).transpose());

        return B.transpose() * D * B * Jacobian.determinant();
    });
}

Vector femSubmesh::internal_nodal_force(Matrix const& x, int const element) const
{
    auto const& cauchy_stresses = variables(InternalVariables::Tensor::Cauchy);

    auto const [m, n] = std::make_tuple(nodes_per_element(), dofs_per_node());

    Vector fint = Vector::Zero(m * n);

    sf->quadrature().integrate(Eigen::Map<RowMatrix>(fint.data(), m, n),
                               [&](auto const& femval, auto const& l) -> RowMatrix {
                                   auto const& [N, dN] = femval;

                                   auto const Jacobian = local_deformation_gradient(dN, x);

                                   auto const& cauchy_stress = cauchy_stresses[view(element, l)];

                                   // Compute the symmetric gradient operator
                                   auto const Bt = dN * Jacobian.inverse();

                                   return Bt * cauchy_stress * Jacobian.determinant();
                               });
    return fint;
}

std::tuple<List const&, Matrix> femSubmesh::consistent_mass(int const element) const
{
    auto X = material_coordinates->initial_configuration(local_node_list(element));

    auto const density_0 = cm->intrinsic_material().initial_density();

    auto m = sf->quadrature().integrate(Matrix::Zero(nodes_per_element(), nodes_per_element()).eval(),
                                        [&](auto const& femval, auto const& l) -> Matrix {
                                            auto const& [N, dN] = femval;

                                            auto const Jacobian = local_deformation_gradient(dN, X);

                                            return N * density_0 * N.transpose()
                                                   * Jacobian.determinant();
                                        });
    return {local_dof_list(element), identity_expansion(m, dofs_per_node())};
}

std::tuple<List const&, Vector> femSubmesh::diagonal_mass(int const element) const
{
    auto const& [dofs, consistent_m] = this->consistent_mass(element);

    Vector diagonal_m(consistent_m.rows());
    for (auto i = 0; i < consistent_m.rows(); ++i)
    {
        diagonal_m(i) = consistent_m.row(i).sum();
    }
    return {local_dof_list(element), diagonal_m};
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
    auto& H_list = variables(InternalVariables::Tensor::DisplacementGradient);
    auto& F_list = variables(InternalVariables::Tensor::DeformationGradient);

#pragma omp parallel for
    for (auto element = 0; element < elements(); ++element)
    {
        // Gather the material coordinates
        auto const X = material_coordinates->initial_configuration(local_node_list(element));
        auto const x = material_coordinates->current_configuration(local_node_list(element));

        sf->quadrature().for_each([&](auto const& femval, const auto& l) {
            auto const& [N, rhea] = femval;

            // Local deformation gradient for the initial configuration
            Matrix3 const F_0 = local_deformation_gradient(rhea, X);
            Matrix3 const F = local_deformation_gradient(rhea, x);

            // Gradient operator in index notation
            auto const& B_0t = rhea * F_0.inverse();

            // Displacement gradient
            Matrix3 const H = (x - X) * B_0t;

            H_list[view(element, l)] = H;
            F_list[view(element, l)] = F * F_0.inverse();
        });
    }
}

void femSubmesh::update_Jacobian_determinants()
{
    auto const& deformation_gradients = variables(InternalVariables::Tensor::DeformationGradient);
    auto& deformation_gradient_determinants = variables(InternalVariables::Scalar::DetF);

    deformation_gradient_determinants = deformation_gradients
                                        | ranges::view::transform(
                                              [](auto const& F) { return F.determinant(); });

    auto const found = ranges::find_if(deformation_gradient_determinants,
                                       [](auto detF) { return detF <= 0.0; });

    if (found != std::end(deformation_gradient_determinants))
    {
        auto const count = ranges::count_if(deformation_gradient_determinants,
                                            [](auto detF) { return detF <= 0.0; });

        auto const i = std::distance(deformation_gradient_determinants.begin(), found);

        auto const [element, quadrature_point] = std::div(i, sf->quadrature().points());

        throw computational_error("Positive Jacobian assumption violated at element "
                                  + std::to_string(element) + " and local quadrature point "
                                  + std::to_string(quadrature_point) + " (" + std::to_string(*found)
                                  + "), another " + std::to_string(count - 1) + " violations found");
    }
}

femSubmesh::ValueCount femSubmesh::nodal_averaged_variable(InternalVariables::Tensor const tensor_name) const
{
    Vector count = Vector::Zero(material_coordinates->size() * 9);
    Vector value = count;

    auto const& tensor_list = variables(tensor_name);

    auto const& E = sf->local_quadrature_extrapolation();

    // Vector format of values
    Vector component = Vector::Zero(sf->quadrature().points());

    for (auto e = 0; e < elements(); ++e)
    {
        // Assemble these into the global value vector
        auto const& node_list = local_node_list(e);

        for (auto ci = 0; ci < 3; ++ci)
        {
            for (auto cj = 0; cj < 3; ++cj)
            {
                for (auto l = 0; l < sf->quadrature().points(); ++l)
                {
                    auto const& tensor = tensor_list[view(e, l)];
                    component(l) = tensor(ci, cj);
                }

                // Local extrapolation to the nodes
                Vector const nodal_component = E * component;

                for (auto n = 0; n < nodal_component.rows(); n++)
                {
                    value(node_list[n] * 9 + ci * 3 + cj) += nodal_component(n);
                    count(node_list[n] * 9 + ci * 3 + cj) += 1.0;
                }
            }
        }
    }
    return {value, count};
}

femSubmesh::ValueCount femSubmesh::nodal_averaged_variable(InternalVariables::Scalar const scalar_name) const
{
    Vector count = Vector::Zero(material_coordinates->size());
    Vector value = count;

    auto const& scalar_list = variables(scalar_name);

    auto const& E = sf->local_quadrature_extrapolation();

    // Vector format of values
    Vector component = Vector::Zero(sf->quadrature().points());

    for (auto e = 0; e < elements(); ++e)
    {
        // Assemble these into the global value vector
        auto const& node_list = local_node_list(e);

        for (auto l = 0; l < sf->quadrature().points(); ++l)
        {
            component(l) = scalar_list[view(e, l)];
        }

        // Local extrapolation to the nodes
        Vector const nodal_component = E * component;

        for (auto n = 0; n < nodal_component.rows(); n++)
        {
            value(node_list[n]) += nodal_component(n);
            count(node_list[n]) += 1.0;
        }
    }
    return {value, count};
}
}
