
#include "mesh/mechanical/beam/fem_submesh.hpp"

#include "geometry/profile_factory.hpp"
#include "math/jacobian_determinant.hpp"
#include "mesh/dof_allocator.hpp"
#include "interpolations/interpolation_factory.hpp"
#include "numeric/float_compare.hpp"

#include <tbb/parallel_for.h>
#include <Eigen/Geometry>

#include <iostream>

namespace neon::mechanical::beam
{
fem_submesh::fem_submesh(json const& material_data,
                         json const& simulation_data,
                         json const& section_data,
                         std::shared_ptr<material_coordinates>& coordinates,
                         basic_submesh const& submesh)
    : basic_submesh(submesh),
      sf(make_line_interpolation(topology(), simulation_data)),
      coordinates(coordinates),
      view(sf->quadrature().points()),
      variables(std::make_shared<internal_variable_type>(elements() * sf->quadrature().points())),
      cm(std::make_unique<isotropic_linear>(variables, material_data))
{
    allocate_normal_and_tangent(section_data);

    dof_allocator(node_indices, dof_indices, traits::dof_order);

    variables->add(variable::second::cauchy_stress,
                   variable::scalar::shear_area_1,
                   variable::scalar::shear_area_2,
                   variable::scalar::cross_sectional_area,
                   variable::scalar::second_moment_area_1,
                   variable::scalar::second_moment_area_2);

    profile = geometry::make_profile(
        json{{"name", "rect"}, {"type", "rectangle"}, {"width", 1.0}, {"height", 1.0}});
}

void fem_submesh::update_internal_variables(double const time_step_size)
{
    auto& A = variables->get(variable::scalar::cross_sectional_area);
    auto& As1 = variables->get(variable::scalar::shear_area_1);
    auto& As2 = variables->get(variable::scalar::shear_area_2);

    auto& area_moment_1 = variables->get(variable::scalar::second_moment_area_1);
    auto& area_moment_2 = variables->get(variable::scalar::second_moment_area_2);

    // Loop over all elements and quadrature points to compute the profile
    // properties for each quadrature point in the beam
    tbb::parallel_for(std::int64_t{0}, elements(), [&, this](auto const element) {
        sf->quadrature().for_each([&, this](auto const&, auto const l) {
            A.at(view(element, l)) = profile->area();

            auto const [shear_area_1, shear_area_2] = profile->shear_area();
            As1.at(view(element, l)) = shear_area_1;
            As2.at(view(element, l)) = shear_area_2;

            auto const [I1, I2] = profile->second_moment_area();
            area_moment_1.at(view(element, l)) = I1;
            area_moment_2.at(view(element, l)) = I2;
        });
    });

    cm->update_internal_variables();
}

std::pair<index_view, matrix const&> fem_submesh::tangent_stiffness(std::int32_t const element) const
{
    static thread_local matrix ke(12, 12);

    auto const& configuration = coordinates->initial_configuration(local_node_view(element));

    ke = rotation_matrix
         * (bending_stiffness(configuration, element) + shear_stiffness(configuration, element)
            + axial_stiffness(configuration, element) + torsional_stiffness(configuration, element))
         * rotation_matrix.transpose();

    return {local_dof_view(element), ke};
}

matrix const& fem_submesh::bending_stiffness(matrix3x const& configuration,
                                             std::int32_t const element) const
{
    static thread_local matrix2x B_bending(2, 6 * sf->nodes());
    static thread_local matrix k_bending(6 * sf->nodes(), 6 * sf->nodes());

    B_bending.setZero();
    k_bending.setZero();

    auto const& D_bending = variables->get(variable::second::bending_stiffness);

    sf->quadrature().integrate_inplace(k_bending, [&, this](auto const& femval, auto const l) {
        auto const& [N, dN] = femval;

        auto const j = jacobian_determinant(configuration * dN);

        for (int i = 0; i < sf->nodes(); ++i)
        {
            auto const offset = i * 6;

            B_bending(0, 3 + offset) = B_bending(1, 4 + offset) = dN(i, l) / j;
        }

        return B_bending.transpose() * D_bending.at(view(element, l)) * B_bending * j;
    });
    return k_bending;
}

matrix const& fem_submesh::shear_stiffness(matrix3x const& configuration,
                                           std::int32_t const element) const
{
    static thread_local matrix2x B_shear(2, 6 * sf->nodes());
    static thread_local matrix k_shear(6 * sf->nodes(), 6 * sf->nodes());

    B_shear.setZero();
    k_shear.setZero();

    auto const& D_shear = variables->get(variable::second::shear_stiffness);

    sf->quadrature().integrate_inplace(k_shear, [&, this](auto const& femval, auto const l) {
        auto const& [N, dN] = femval;

        auto const j = jacobian_determinant(configuration * dN);

        for (int i = 0; i < sf->nodes(); ++i)
        {
            auto const offset = i * 6;

            B_shear(0, 0 + offset) = B_shear(1, 1 + offset) = dN(i, l) / j;

            B_shear(0, 4 + offset) = -N(i, l);
            B_shear(1, 3 + offset) = N(i, l);
        }
        return B_shear.transpose() * D_shear.at(view(element, l)) * B_shear * j;
    });
    return k_shear;
}

matrix const& fem_submesh::axial_stiffness(matrix3x const& configuration,
                                           std::int32_t const element) const
{
    static thread_local vector B_axial(6 * sf->nodes());
    static thread_local matrix k_axial(6 * sf->nodes(), 6 * sf->nodes());

    B_axial.setZero();
    k_axial.setZero();

    auto const D_axial = variables->get(variable::scalar::axial_stiffness);

    sf->quadrature().integrate_inplace(k_axial, [&, this](auto const& femval, auto const l) {
        auto const& [N, dN] = femval;

        auto const j = jacobian_determinant(configuration * dN);

        for (int i = 0; i < sf->nodes(); ++i)
        {
            auto const offset = i * 6;

            B_axial(2 + offset) = dN(i, l) / j;
        }

        return B_axial * D_axial.at(view(element, l)) * B_axial.transpose() * j;
    });
    return k_axial;
}

matrix const& fem_submesh::torsional_stiffness(matrix3x const& configuration,
                                               std::int32_t const element) const
{
    static thread_local vector B_torsion(6 * sf->nodes());
    static thread_local matrix k_torsion(6 * sf->nodes(), 6 * sf->nodes());

    B_torsion.setZero();
    k_torsion.setZero();

    auto const D_torsion = variables->get(variable::scalar::torsional_stiffness);

    sf->quadrature().integrate_inplace(k_torsion, [&, this](auto const& femval, auto const l) {
        auto const& [N, dN] = femval;

        auto const j = jacobian_determinant(configuration * dN);

        for (int i = 0; i < sf->nodes(); ++i)
        {
            auto const offset = i * 6;

            B_torsion(5 + offset) = dN(i, l) / j;
        }

        return B_torsion * D_torsion.at(view(element, l)) * B_torsion.transpose() * j;
    });
    return k_torsion;
}

void fem_submesh::allocate_normal_and_tangent(json const& section_data)
{
    if (section_data.find("tangent") == section_data.end())
    {
        throw std::domain_error("A \"tangent\" vector must be specified in the \"section\"");
    }
    if (section_data.find("normal") == section_data.end())
    {
        throw std::domain_error("A \"normal\" vector must be specified in the \"section\"");
    }

    auto const& tangent_vector = section_data["tangent"];

    if (!tangent_vector.is_array() || tangent_vector.size() != 3)
    {
        throw std::domain_error("A \"tangent\" vector must be specified using three (3) "
                                "coordinates [x, y, z]");
    }

    auto const& normal_vector = section_data["normal"];

    if (!normal_vector.is_array() || normal_vector.size() != 3)
    {
        throw std::domain_error("A \"normal\" vector must be specified using three (3) coordinates "
                                "[x, y, z]");
    }

    tangent(0) = tangent_vector[0];
    tangent(1) = tangent_vector[1];
    tangent(2) = tangent_vector[2];

    normal(0) = normal_vector[0];
    normal(1) = normal_vector[1];
    normal(2) = normal_vector[2];

    if (!is_approx(normal.dot(tangent), 0.0))
    {
        throw std::domain_error("normal and tangent vectors must be orthogonal");
    }

    matrix3 rotation;

    // Local coordinate (x1)
    rotation.col(0) = normal;
    // Local coordinate (x2)
    rotation.col(1) = normal.cross(tangent);
    // Local coordinate (x3)
    rotation.col(2) = tangent;

    rotation_matrix = matrix12::Zero();

    rotation_matrix.block<3, 3>(0, 0) = rotation;
    rotation_matrix.block<3, 3>(3, 3) = rotation;
    rotation_matrix.block<3, 3>(6, 6) = rotation;
    rotation_matrix.block<3, 3>(9, 9) = rotation;
}
}
