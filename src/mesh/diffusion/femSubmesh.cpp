
#include "femSubmesh.hpp"

#include "Exceptions.hpp"

#include "constitutive/constitutive_model_factory.hpp"
#include "interpolations/interpolation_factory.hpp"
#include "material/material_property.hpp"
#include "mesh/DofAllocator.hpp"
#include "mesh/MaterialCoordinates.hpp"
#include "numeric/gradient_operator.hpp"

#include <cfenv>
#include <omp.h>

#include <termcolor/termcolor.hpp>

namespace neon::diffusion
{
femSubmesh::femSubmesh(json const& material_data,
                       json const& mesh_data,
                       std::shared_ptr<MaterialCoordinates>& material_coordinates,
                       Submesh const& submesh)
    : Submesh(submesh),
      material_coordinates(material_coordinates),
      sf(make_volume_interpolation(topology(), mesh_data)),
      variables(std::make_shared<InternalVariables>(elements() * sf->quadrature().points())),
      cm(make_constitutive_model(variables, material_data, mesh_data))
{
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

std::tuple<local_indices const&, matrix> femSubmesh::tangent_stiffness(int const element) const
{
    auto const X = material_coordinates->current_configuration(local_node_list(element));

    auto const n = nodes_per_element();

    auto const& D_Vec = variables->fetch(InternalVariables::Tensor::Conductivity);

    matrix const kmat = sf->quadrature()
                            .integrate(matrix::Zero(n, n).eval(),
                                       [&](auto const& femval, auto const& l) -> matrix {
                                           auto const& [N, rhea] = femval;

                                           auto const& D = D_Vec[offset(element, l)];

                                           matrix3 const Jacobian = local_jacobian(rhea, X);

                                           // Compute the symmetric gradient operator
                                           matrix const B = (rhea * Jacobian.inverse()).transpose();

                                           return B.transpose() * D * B * Jacobian.determinant();
                                       });

    return {local_dof_list(element), kmat};
}

std::tuple<local_indices const&, matrix> femSubmesh::consistent_mass(int const element) const
{
    auto X = material_coordinates->current_configuration(local_node_list(element));

    auto const density = cm->intrinsic_material().initial_density();
    auto const specific_heat = cm->intrinsic_material().specific_heat();

    auto m = sf->quadrature().integrate(matrix::Zero(nodes_per_element(), nodes_per_element()).eval(),
                                        [&](auto const& femval, auto const& l) -> matrix {
                                            auto const& [N, dN] = femval;

                                            auto const Jacobian = local_jacobian(dN, X);

                                            return N * density * specific_heat * N.transpose()
                                                   * Jacobian.determinant();
                                        });
    return {local_dof_list(element), m};
}

std::tuple<local_indices const&, vector> femSubmesh::diagonal_mass(int const element) const
{
    auto const& [dofs, consistent_m] = this->consistent_mass(element);

    vector diagonal_m(consistent_m.rows());
    for (auto i = 0; i < consistent_m.rows(); ++i)
    {
        diagonal_m(i) = consistent_m.row(i).sum();
    }
    return {local_dof_list(element), diagonal_m};
}

void femSubmesh::update_internal_variables(double const time_step_size)
{
    std::feclearexcept(FE_ALL_EXCEPT);

    cm->update_internal_variables(time_step_size);

    if (std::fetestexcept(FE_INVALID))
    {
        throw computational_error("Floating point error reported\n");
    }
}

int femSubmesh::offset(int const element, int const quadrature_point) const
{
    return sf->quadrature().points() * element + quadrature_point;
}

femSubmesh::ValueCount femSubmesh::nodal_averaged_variable(InternalVariables::Scalar const scalar_name) const
{
    vector count = vector::Zero(material_coordinates->size());
    vector value = count;

    auto const& scalar_list = variables->fetch(scalar_name);

    auto const& E = sf->local_quadrature_extrapolation();

    // vector format of values
    vector component = vector::Zero(sf->quadrature().points());

    for (auto e = 0; e < elements(); ++e)
    {
        // Assemble these into the global value vector
        auto const& node_list = local_node_list(e);

        for (auto l = 0; l < sf->quadrature().points(); ++l)
        {
            component(l) = scalar_list[this->offset(e, l)];
        }

        // Local extrapolation to the nodes
        vector const nodal_component = E * component;

        for (auto n = 0; n < nodal_component.rows(); n++)
        {
            value(node_list[n]) += nodal_component(n);
            count(node_list[n]) += 1.0;
        }
    }
    return {value, count};
}
}
