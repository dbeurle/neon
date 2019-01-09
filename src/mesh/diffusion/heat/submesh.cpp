
#include "submesh.hpp"

#include "exceptions.hpp"

#include "constitutive/constitutive_model_factory.hpp"
#include "interpolations/interpolation_factory.hpp"
#include "material/material_property.hpp"
#include "mesh/material_coordinates.hpp"
#include "numeric/gradient_operator.hpp"

#include <cfenv>
#include <omp.h>

#include <termcolor/termcolor.hpp>

namespace neon::diffusion
{
submesh::submesh(json const& material_data,
                 json const& mesh_data,
                 std::shared_ptr<material_coordinates>& coordinates,
                 basic_submesh const& submesh)
    : basic_submesh(submesh),
      coordinates(coordinates),
      bilinear_gradient(topology(), mesh_data),
      bilinear(topology(), mesh_data),
      view(bilinear.quadrature().points()),
      variables(std::make_shared<internal_variables_t>(elements() * bilinear.quadrature().points())),
      cm(make_constitutive_model(variables, material_data, mesh_data)),
      patch_recovery(nullptr)
{
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

auto submesh::tangent_stiffness(std::int64_t const element) const -> matrix const&
{
    thread_local matrix k_mat;

    k_mat.resize(nodes_per_element(), nodes_per_element());

    auto const X = coordinates->current_configuration(local_node_view(element));

    auto const& D_Vec = variables->get(variable::second::conductivity);

    bilinear_gradient.integrate(k_mat.setZero(),
                                        [&](auto const& value, auto const index) -> matrix {
                                            auto const& [N, dN] = value;

                                            auto const& D = D_Vec[view(element, index)];

                                            matrix3 const Jacobian = local_jacobian(dN, X);

                                            // Compute the symmetric gradient operator
                                            matrix const B = (dN * Jacobian.inverse()).transpose();

                                            return B.transpose() * D * B * Jacobian.determinant();
                                        });
    return k_mat;
}

auto submesh::consistent_mass(std::int64_t const element) const -> matrix const&
{
    thread_local matrix mass;

    auto const X = coordinates->current_configuration(local_node_view(element));

    auto const density = cm->intrinsic_material().initial_density();
    auto const specific_heat = cm->intrinsic_material().specific_heat();

    mass.resize(nodes_per_element(), nodes_per_element());

    bilinear.integrate(mass.setZero(), [&](auto const& value, auto) -> matrix {
        auto const& [N, dN] = value;

        matrix3 const jacobian = local_jacobian(dN, X);

        return N * density * specific_heat * N.transpose() * jacobian.determinant();
    });
    return mass;
}

auto submesh::diagonal_mass(std::int64_t const element) const -> vector const&
{
    thread_local vector mass;

    mass = this->consistent_mass(element).rowwise().sum();

    return mass;
}

void submesh::update_internal_variables(double const time_step_size)
{
    std::feclearexcept(FE_ALL_EXCEPT);

    cm->update_internal_variables(time_step_size);

    if (std::fetestexcept(FE_INVALID))
    {
        throw computational_error("Floating point error reported\n");
    }
}

std::pair<vector, vector> submesh::nodal_averaged_variable(variable::second const tensor_name) const
{
    vector count = vector::Zero(coordinates->size() * 9);
    vector value = count;

    auto const& tensor_list = variables->get(tensor_name);

    auto const& E = patch_recovery->extrapolation_matrix();

    // vector format of values
    vector component = vector::Zero(bilinear_gradient.quadrature().points());

    for (std::int64_t e{0}; e < elements(); ++e)
    {
        // Assemble these into the global value vector
        auto const& node_list = local_node_view(e);

        for (auto ci = 0; ci < 3; ++ci)
        {
            for (auto cj = 0; cj < 3; ++cj)
            {
                for (std::size_t l{0}; l < bilinear_gradient.quadrature().points(); ++l)
                {
                    auto const& tensor = tensor_list[view(e, l)];
                    component(l) = tensor(ci, cj);
                }

                // Local extrapolation to the nodes
                vector const nodal_component = E * component;

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

std::pair<vector, vector> submesh::nodal_averaged_variable(variable::scalar const name) const
{
    vector count = vector::Zero(coordinates->size());
    vector value = count;

    auto const& scalar_list = variables->get(name);

    auto const& E = patch_recovery->extrapolation_matrix();

    // vector format of values
    vector component = vector::Zero(bilinear_gradient.quadrature().points());

    for (std::int64_t element{0}; element < elements(); ++element)
    {
        // Assemble these into the global value vector
        auto const& node_list = local_node_view(element);

        for (std::size_t l{0}; l < bilinear_gradient.quadrature().points(); ++l)
        {
            component(l) = scalar_list[view(element, l)];
        }

        // Local extrapolation to the nodes
        vector const nodal_component = E * component;

        for (auto n = 0l; n < nodal_component.rows(); n++)
        {
            value(node_list[n]) += nodal_component(n);
            count(node_list[n]) += 1.0;
        }
    }
    return {value, count};
}
}
