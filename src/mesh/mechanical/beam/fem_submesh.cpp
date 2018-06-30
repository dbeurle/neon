
#include "mesh/mechanical/beam/fem_submesh.hpp"

#include "mesh/dof_allocator.hpp"
#include "interpolations/interpolation_factory.hpp"
#include "math/transform_expand.hpp"

#include <tbb/parallel_for.h>
#include <Eigen/Geometry>

namespace neon::mechanical::beam
{
fem_submesh::fem_submesh(json const& material_data,
                         json const& simulation_data,
                         std::shared_ptr<material_coordinates>& coordinates,
                         basic_submesh const& submesh)
    : basic_submesh(submesh),
      sf(make_line_interpolation(topology(), simulation_data)),
      coordinates(coordinates),
      profiles(elements() * sf->quadrature().points()),
      view(sf->quadrature().points()),
      variables(std::make_shared<internal_variable_type>(elements() * sf->quadrature().points())),
      cm(std::make_unique<isotropic_linear>(variables, material_data))
{
    dof_allocator(node_indices, dof_indices, traits::dof_order);

    orientations.resize(elements(), vector3::UnitX());
    tangents.resize(elements(), vector3::UnitZ());

    variables->add(variable::second::cauchy_stress,
                   variable::scalar::shear_area_1,
                   variable::scalar::shear_area_2,
                   variable::scalar::cross_sectional_area,
                   variable::scalar::second_moment_area_1,
                   variable::scalar::second_moment_area_2);
}

void fem_submesh::update_internal_variables(double const time_step_size)
{
    for (auto& profile : profiles)
    {
        profile = std::make_unique<geometry::rectangular_bar>(1.0, 1.0);
    }

    auto& A = variables->get(variable::scalar::cross_sectional_area);
    auto& As1 = variables->get(variable::scalar::shear_area_1);
    auto& As2 = variables->get(variable::scalar::shear_area_2);

    auto& area_moment_1 = variables->get(variable::scalar::second_moment_area_1);
    auto& area_moment_2 = variables->get(variable::scalar::second_moment_area_2);

    // Loop over all elements and quadrature points to compute the profile
    // properties for each quadrature point in the beam
    tbb::parallel_for(std::int64_t{0}, elements(), [&, this](auto const element) {
        sf->quadrature().for_each([&, this](auto const& femval, auto const l) {
            [[maybe_unused]] auto const& [N, dN] = femval;

            A.at(view(element, l)) = profiles.at(element)->area();

            auto const [shear_area_1, shear_area_2] = profiles.at(element)->shear_area();
            As1.at(view(element, l)) = shear_area_1;
            As2.at(view(element, l)) = shear_area_2;

            auto const [I1, I2] = profiles.at(element)->second_moment_area();
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

    matrix12 const T = rotation_matrix(element);

    ke = T
         * (bending_stiffness(configuration, element) + shear_stiffness(configuration, element)
            + axial_stiffness(configuration, element) + torsional_stiffness(configuration, element))
         * T.transpose();

    return {local_dof_view(element), ke};
}

matrix const& fem_submesh::bending_stiffness(matrix const& configuration,
                                             std::int32_t const element) const
{
    static thread_local matrix2x B_b(2, 6 * sf->nodes());
    static thread_local matrix k_b(6 * sf->nodes(), 6 * sf->nodes());

    B_b.setZero();
    k_b.setZero();

    auto const& D_b = variables->get(variable::second::bending_stiffness);

    sf->quadrature().integrate_inplace(k_b, [&, this](auto const& femval, auto const l) {
        auto const& [N, dN] = femval;

        double const jacobian = configuration.row(0) * dN;

        for (int i = 0; i < sf->nodes(); ++i)
        {
            auto const offset = i * 6;

            B_b(0, 3 + offset) = B_b(1, 4 + offset) = dN(i, l) / jacobian;
        }

        return B_b.transpose() * D_b.at(view(element, l)) * B_b * jacobian;
    });
    return k_b;
}

matrix const& fem_submesh::shear_stiffness(matrix const& configuration, std::int32_t const element) const
{
    static thread_local matrix2x B_s(2, 6 * sf->nodes());
    static thread_local matrix k_s(6 * sf->nodes(), 6 * sf->nodes());

    B_s.setZero();
    k_s.setZero();

    auto const& D_s = variables->get(variable::second::shear_stiffness);

    sf->quadrature().integrate_inplace(k_s, [&, this](auto const& femval, auto const l) {
        auto const& [N, dN] = femval;

        double const jacobian = configuration.row(0) * dN;

        for (int i = 0; i < sf->nodes(); ++i)
        {
            auto const offset = i * 6;

            B_s(0, 0 + offset) = B_s(1, 1 + offset) = dN(i, l) / jacobian;

            B_s(0, 4 + offset) = -N(i, l);
            B_s(1, 3 + offset) = N(i, l);
        }

        return B_s.transpose() * D_s.at(view(element, l)) * B_s * jacobian;
    });
    return k_s;
}

matrix const& fem_submesh::axial_stiffness(matrix const& configuration, std::int32_t const element) const
{
    static thread_local vector B_a(6 * sf->nodes());
    static thread_local matrix k_a(6 * sf->nodes(), 6 * sf->nodes());

    B_a.setZero();
    k_a.setZero();

    auto const D_a = variables->get(variable::scalar::axial_stiffness);

    sf->quadrature().integrate_inplace(k_a, [&, this](auto const& femval, auto const l) {
        auto const& [N, dN] = femval;

        double const jacobian = configuration.row(0) * dN;

        for (int i = 0; i < sf->nodes(); ++i)
        {
            auto const offset = i * 6;

            B_a(2 + offset) = dN(i, l) / jacobian;
        }

        return B_a * D_a.at(view(element, l)) * B_a.transpose() * jacobian;
    });
    return k_a;
}

matrix const& fem_submesh::torsional_stiffness(matrix const& configuration,
                                               std::int32_t const element) const
{
    static thread_local vector B_t(6 * sf->nodes());
    static thread_local matrix k_t(6 * sf->nodes(), 6 * sf->nodes());

    B_t.setZero();
    k_t.setZero();

    auto const D_t = variables->get(variable::scalar::torsional_stiffness);

    sf->quadrature().integrate_inplace(k_t, [&, this](auto const& femval, auto const l) {
        auto const& [N, dN] = femval;

        double const jacobian = configuration.row(0) * dN;

        for (int i = 0; i < sf->nodes(); ++i)
        {
            auto const offset = i * 6;

            B_t(5 + offset) = dN(i, l) / jacobian;
        }

        return B_t * D_t.at(view(element, l)) * B_t.transpose() * jacobian;
    });
    return k_t;
}

matrix12 fem_submesh::rotation_matrix(std::int32_t const element) const
{
    matrix3 rotation;

    // Local coordinate (x1)
    rotation.col(0) = orientations.at(element);

    // Local coordinate (x3)
    rotation.col(2) = tangents.at(element);

    // Local coordinate (x2)
    rotation.col(1) = rotation.col(0).cross(rotation.col(2)).normalized().cwiseAbs();

    matrix12 element_rotation = matrix12::Zero();

    element_rotation.block<3, 3>(0, 0) = rotation;
    element_rotation.block<3, 3>(3, 3) = rotation;
    element_rotation.block<3, 3>(6, 6) = rotation;
    element_rotation.block<3, 3>(9, 9) = rotation;

    return element_rotation;
}
}
