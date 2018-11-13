
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
      sf(make_surface_interpolation(topology(), simulation_data)),
      view(sf->quadrature().points()),
      variables(std::make_shared<internal_variables_t>(elements() * sf->quadrature().points())),
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

    dof_allocator(node_indices, dof_list, traits::dof_order);
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

std::pair<index_view, matrix> submesh::tangent_stiffness(std::int32_t const element) const
{
    auto const x = geometry::project_to_plane(
        coordinates->current_configuration(local_node_view(element)));

    thread_local matrix k_e;

    k_e = material_tangent_stiffness(x, element);

    if (!cm->is_finite_deformation()) return {local_dof_view(element), k_e};

    k_e.noalias() += geometric_tangent_stiffness(x, element);

    return {local_dof_view(element), k_e};
}

std::pair<index_view, vector> submesh::internal_force(std::int32_t const element) const
{
    auto const x = geometry::project_to_plane(
        coordinates->current_configuration(local_node_view(element)));

    return {local_dof_view(element), internal_nodal_force(x, element)};
}

matrix const& submesh::geometric_tangent_stiffness(matrix2x const& x, std::int32_t const element) const
{
    auto const& cauchy_stresses = variables->get(variable::second::cauchy_stress);

    thread_local matrix k_geo(nodes_per_element(), nodes_per_element());
    thread_local matrix k_geo_full;

    sf->quadrature().integrate_inplace(k_geo, [&](auto const& N_dN, auto const l) {
        auto const& [N, rhea] = N_dN;

        matrix2 const Jacobian = local_deformation_gradient(rhea, x);

        auto const cauchy = cauchy_stresses[view(element, l)];

        // Compute the symmetric gradient operator
        auto const L = local_gradient(rhea, Jacobian);

        return L.transpose() * cauchy * L * Jacobian.determinant();
    });
    k_geo_full = identity_expansion(k_geo, dofs_per_node());

    return k_geo_full;
}

matrix const& submesh::material_tangent_stiffness(matrix2x const& x, std::int32_t const element) const
{
    auto const local_dofs = nodes_per_element() * dofs_per_node();

    thread_local matrix k_mat(local_dofs, local_dofs);

    k_mat.setZero();

    auto const& tangent_operators = variables->get(variable::fourth::tangent_operator);

    matrix B = matrix::Zero(4, local_dofs);

    sf->quadrature().integrate_inplace(k_mat, [&](auto const& femval, auto const& l) {
        auto const& [N, rhea] = femval;

        auto const& D = tangent_operators[view(element, l)];

        matrix2 const Jacobian{local_deformation_gradient(rhea, x)};

        symmetric_gradient<2>(B, (rhea * Jacobian.inverse()).transpose());

        return B.transpose() * D * B * Jacobian.determinant();
    });
    return k_mat;
}

vector const& submesh::internal_nodal_force(matrix2x const& x, std::int32_t const element) const
{
    thread_local vector f_int(nodes_per_element() * dofs_per_node());

    f_int.setZero();

    auto const& cauchy_stresses = variables->get(variable::second::cauchy_stress);

    sf->quadrature()
        .integrate_inplace(Eigen::Map<row_matrix>(f_int.data(), nodes_per_element(), dofs_per_node()),
                           [&](auto const& N_dN, auto const l) {
                               auto const& [N, dN] = N_dN;

                               matrix2 const Jacobian = local_deformation_gradient(dN, x);

                               matrix2 const& cauchy_stress = cauchy_stresses[view(element, l)];

                               // symmetric gradient operator
                               auto const Bt = dN * Jacobian.inverse();

                               return Bt * cauchy_stress * Jacobian.determinant();
                           });
    return f_int;
}

std::pair<index_view, matrix> submesh::consistent_mass(std::int32_t const element) const
{
    auto const X = coordinates->initial_configuration(local_node_view(element));
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

std::pair<index_view, vector> submesh::diagonal_mass(std::int32_t const element) const
{
    auto const& [dofs, consistent_m] = this->consistent_mass(element);

    vector diagonal_m(consistent_m.rows());
    for (auto i = 0; i < consistent_m.rows(); ++i)
    {
        diagonal_m(i) = consistent_m.row(i).sum();
    }
    return {local_dof_view(element), diagonal_m};
}

void submesh::update_internal_variables(double const time_step_size)
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

void submesh::update_Jacobian_determinants()
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

        auto const [element, quadrature_point] = std::div(i, sf->quadrature().points());

        throw computational_error("Positive Jacobian assumption violated at element "
                                  + std::to_string(element) + " and local quadrature point "
                                  + std::to_string(quadrature_point) + " (" + std::to_string(*found)
                                  + ")");
    }
}

std::pair<vector, vector> submesh::nodal_averaged_variable(variable::scalar const scalar_name) const
{
    vector count = vector::Zero(coordinates->size());
    vector value = count;

    auto const& scalar_list = variables->get(scalar_name);

    auto const& E = sf->local_quadrature_extrapolation();

    // vector format of values
    vector component = vector::Zero(sf->quadrature().points());

    for (std::int64_t element{0}; element < elements(); ++element)
    {
        for (std::size_t l{0}; l < sf->quadrature().points(); ++l)
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

    matrix const E = sf->local_quadrature_extrapolation();

    // vector format of values
    vector component(sf->quadrature().points());

    for (std::int64_t element{0}; element < elements(); ++element)
    {
        // Assemble these into the global value vector
        auto const& node_list = local_node_view(element);

        for (auto ci = 0; ci < 2; ++ci)
        {
            for (auto cj = 0; cj < 2; ++cj)
            {
                for (std::size_t l{0}; l < sf->quadrature().points(); ++l)
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
