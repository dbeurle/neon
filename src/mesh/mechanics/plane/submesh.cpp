
#include "submesh.hpp"

#include "exceptions.hpp"

#include "constitutive/constitutive_model_factory.hpp"
#include "geometry/projection.hpp"
#include "interpolations/interpolation_factory.hpp"
#include "material/material_property.hpp"
#include "mesh/material_coordinates.hpp"
#include "mesh/dof_allocator.hpp"
#include "numeric/mechanics"
#include "traits/mechanics.hpp"

#include <cfenv>
#include <chrono>

#include <termcolor/termcolor.hpp>

namespace neon::mechanics::plane
{
submesh::submesh(json const& material_data,
                 json const& simulation_data,
                 std::shared_ptr<material_coordinates>& coordinates,
                 basic_submesh const& submesh)
    : detail::submesh<plane::submesh, plane::internal_variables_t>(submesh),
      coordinates{coordinates},
      bilinear_gradient(topology(), simulation_data),
      bilinear(topology(), simulation_data),
      view(bilinear_gradient.quadrature().points()),
      variables(std::make_shared<internal_variables_t>(elements()
                                                       * bilinear_gradient.quadrature().points())),
      cm(make_constitutive_model(variables, material_data, simulation_data))
{
    // Allocate storage for the displacement gradient
    variables->add(variable::second::displacement_gradient,
                   variable::second::deformation_gradient,
                   variable::second::cauchy_stress,
                   variable::scalar::DetF);

    // Get the old data to the undeformed configuration
    for (auto& F : variables->get(variable::second::deformation_gradient))
    {
        F = matrix2::Identity();
    }

    variables->commit();

    dof_allocator(node_indices, dof_list, traits::dofs_per_node);
}

void submesh::save_internal_variables(bool const have_converged)
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

auto submesh::tangent_stiffness(std::int32_t const element) const -> matrix const&
{
    auto const x = geometry::project_to_plane(
        coordinates->current_configuration(local_node_view(element)));

    thread_local matrix k_e;

    k_e = material_tangent_stiffness(x, element);

    if (cm->is_finite_deformation())
    {
        k_e.noalias() += geometric_tangent_stiffness(x, element);
    }
    return k_e;
}

auto submesh::internal_force(std::int32_t const element) const -> vector const&
{
    thread_local vector f_int(nodes_per_element() * dofs_per_node());

    f_int.setZero();

    auto const x = geometry::project_to_plane(
        coordinates->current_configuration(local_node_view(element)));

    auto const& cauchy_stresses = variables->get(variable::second::cauchy_stress);

    Eigen::Map<row_matrix> matrix_view(f_int.data(), nodes_per_element(), dofs_per_node());

    bilinear_gradient.integrate(matrix_view, [&](auto const& N_dN, auto const index) {
        auto const& [N, dN] = N_dN;

        matrix2 const Jacobian = local_deformation_gradient(dN, x);

        matrix2 const& cauchy_stress = cauchy_stresses[view(element, index)];

        // symmetric gradient operator
        auto const Bt = dN * Jacobian.inverse();

        return Bt * cauchy_stress * Jacobian.determinant();
    });
    return f_int;
}

auto submesh::geometric_tangent_stiffness(matrix2x const& x, std::int32_t const element) const
    -> matrix const&
{
    auto const& cauchy_stresses = variables->get(variable::second::cauchy_stress);

    thread_local matrix k_geo(nodes_per_element(), nodes_per_element());
    thread_local matrix k_geo_full;

    bilinear_gradient.integrate(k_geo, [&](auto const& value, auto const l) {
        auto const& [N, dN] = value;

        matrix2 const Jacobian = local_deformation_gradient(dN, x);

        auto const cauchy = cauchy_stresses[view(element, l)];

        // Compute the symmetric gradient operator
        auto const L = local_gradient(dN, Jacobian);

        return L.transpose() * cauchy * L * Jacobian.determinant();
    });
    k_geo_full = identity_expansion(k_geo, dofs_per_node());

    return k_geo_full;
}

matrix const& submesh::material_tangent_stiffness(matrix2x const& x, std::int32_t const element) const
{
    auto const local_dofs = nodes_per_element() * dofs_per_node();

    thread_local matrix k_mat(local_dofs, local_dofs);
    thread_local matrix B = matrix::Zero(4, local_dofs);

    k_mat.setZero();

    auto const& tangent_operators = variables->get(variable::fourth::tangent_operator);

    bilinear_gradient.integrate(k_mat, [&](auto const& value, auto const index) {
        auto const& [N, dN] = value;

        auto const& D = tangent_operators[view(element, index)];

        matrix2 const jacobian = local_deformation_gradient(dN, x);

        symmetric_gradient<2>(B, (dN * jacobian.inverse()).transpose());

        return B.transpose() * D * B * jacobian.determinant();
    });
    return k_mat;
}

auto submesh::consistent_mass(std::int32_t const element) const -> matrix const&
{
    thread_local matrix local_mass(nodes_per_element(), nodes_per_element());

    thread_local matrix mass(nodes_per_element() * dofs_per_node(),
                             nodes_per_element() * dofs_per_node());

    auto const& X = geometry::project_to_plane(
        coordinates->initial_configuration(local_node_view(element)));

    auto const density = cm->intrinsic_material().initial_density();

    bilinear.integrate(local_mass.setZero(), [&](auto const& value, auto) -> matrix {
        auto const& [N, dN] = value;

        matrix2 const jacobian = local_deformation_gradient(dN, X);

        return density * N * N.transpose() * jacobian.determinant();
    });

    identity_expansion_inplace<2>(local_mass, mass.setZero());

    return mass;
}

auto submesh::diagonal_mass(std::int32_t const element) const -> vector const&
{
    thread_local vector diagonal_mass;

    diagonal_mass = this->consistent_mass(element).rowwise().sum();

    return diagonal_mass;
}

void submesh::update_internal_variables(double const time_step_size)
{
    std::feclearexcept(FE_ALL_EXCEPT);

    update_deformation_measures();

    update_jacobian_determinants();

    cm->update_internal_variables(time_step_size);

    if (std::fetestexcept(FE_INVALID))
    {
        throw computational_error("Floating point error reported\n");
    }
}

void submesh::update_deformation_measures()
{
    auto& H_list = variables->get(variable::second::displacement_gradient);
    auto& F_list = variables->get(variable::second::deformation_gradient);

    for (std::int64_t element{0}; element < elements(); ++element)
    {
        // Gather the material coordinates
        auto const X = geometry::project_to_plane(
            coordinates->initial_configuration(local_node_view(element)));
        auto const x = geometry::project_to_plane(
            coordinates->current_configuration(local_node_view(element)));

        bilinear_gradient.for_each([&](auto const& value, auto const index) {
            auto const& [N, dN] = value;

            // Local deformation gradient for the initial configuration
            matrix2 const F_0 = local_deformation_gradient(dN, X);
            matrix2 const F = local_deformation_gradient(dN, x);

            // Gradient operator in index notation
            auto const& B_0t = dN * F_0.inverse();

            // Displacement gradient
            matrix2 const H = (x - X) * B_0t;

            H_list[view(element, index)] = H;
            F_list[view(element, index)] = F * F_0.inverse();
        });
    }
}

void submesh::update_jacobian_determinants()
{
    auto const& F = variables->get(variable::second::deformation_gradient);
    auto& F_det = variables->get(variable::scalar::DetF);

    std::transform(begin(F), end(F), begin(F_det), [](auto const& F) { return F.determinant(); });

    auto const found = std::find_if(begin(F_det), end(F_det), [](auto const value) {
        return value <= 0.0;
    });

    if (found != F_det.end())
    {
        auto const i = std::distance(begin(F_det), found);

        auto const [element, quadrature_point] = std::div(i, bilinear_gradient.quadrature().points());

        throw computational_error("Positive Jacobian assumption violated at element "
                                  + std::to_string(element) + " and local quadrature point "
                                  + std::to_string(quadrature_point) + " (" + std::to_string(*found)
                                  + ")");
    }
}

auto submesh::nodal_averaged_variable(variable::scalar const scalar_name) const
    -> std::pair<vector, vector>
{
    vector count = vector::Zero(coordinates->size());
    vector value = count;

    auto const& scalar_list = variables->get(scalar_name);

    matrix const& E = patch_recovery->extrapolation_matrix();

    // vector format of values
    vector component = vector::Zero(bilinear_gradient.quadrature().points());

    for (std::int64_t element{0}; element < elements(); ++element)
    {
        for (std::size_t l{0}; l < bilinear_gradient.quadrature().points(); ++l)
        {
            component(l) = scalar_list[view(element, l)];
        }

        // Local extrapolation to the nodes
        vector const nodal_component = E * component;

        // Assemble these into the global value vector
        auto const& node_list = local_node_view(element);

        for (auto n = 0; n < nodal_component.rows(); n++)
        {
            value(node_list[n]) += nodal_component(n);
            count(node_list[n]) += 1.0;
        }
    }
    return {value, count};
}

std::pair<vector, vector> submesh::nodal_averaged_variable(variable::second const tensor_name) const
{
    vector count = vector::Zero(coordinates->size() * 4);
    vector value = count;

    auto const& tensor_list = variables->get(tensor_name);

    matrix const& E = patch_recovery->extrapolation_matrix();

    // vector format of values
    vector component(bilinear_gradient.quadrature().points());

    for (std::int64_t element{0}; element < elements(); ++element)
    {
        // Assemble these into the global value vector
        auto const& node_list = local_node_view(element);

        for (auto ci = 0; ci < 2; ++ci)
        {
            for (auto cj = 0; cj < 2; ++cj)
            {
                for (std::size_t l{0}; l < bilinear_gradient.quadrature().points(); ++l)
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
}
