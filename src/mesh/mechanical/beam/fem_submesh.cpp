
#include "mesh/mechanical/beam/fem_submesh.hpp"

#include "interpolations/interpolation_factory.hpp"
#include "math/transform_expand.hpp"

namespace neon::mechanical::beam
{
fem_submesh::fem_submesh(json const& material_data,
                         json const& simulation_data,
                         std::shared_ptr<material_coordinates>& material_coordinates,
                         basic_submesh const& submesh)
    : basic_submesh(submesh),
      profiles(elements()),
      material(material_data),
      sf(make_line_interpolation(topology(), simulation_data)),
      mesh_coordinates(material_coordinates),
      iv(std::make_shared<internal_variables_type>(elements() * sf->quadrature().points())),
      nodal_variables(mesh_coordinates->size())
{
    dof_list.resize(connectivity.rows() * traits::dofs_per_node, connectivity.cols());

    for (indices::Index i{0}; i < connectivity.cols(); ++i)
    {
        transform_expand_view(connectivity(Eigen::placeholders::all, i),
                              dof_list(Eigen::placeholders::all, i),
                              traits::dof_order);
    }
    iv->add(internal_variables_type::Tensor::CauchyStress);

    iv->add(internal_variables_type::scalar::shear_area_1,
            internal_variables_type::scalar::shear_area_2,
            internal_variables_type::scalar::cross_sectional_area,
            internal_variables_type::scalar::second_moment_area_1,
            internal_variables_type::scalar::second_moment_area_2);
}

void fem_submesh::update_internal_variables(double const time_step_size = 1.0)
{
    for (auto& profile : profiles)
    {
        profile = std::make_unique<geometry::rectangular_bar>(1.0, 1.0);
    }

    auto& [A, As1, As2] = iv->fetch(internal_variables_type::scalar::cross_sectional_area,
                                    internal_variables_type::scalar::shear_area_1,
                                    internal_variables_type::scalar::shear_area_1);

    // Loop over all elements and quadrature points to compute the profile
    // properties for each quadrature point in the beam
    tbb::parallel_for(std::int64_t{}, elements(), [](auto const element) {
        sf->quadrature().for_each([&, this](auto const& femval, auto const l) {
            auto const& [N, dN] = femval;

            A[view(element, l)] = profile[element]->area();

            auto const [shear_area_1, shear_area_2] = profile[element]->shear_area();

            As1[view(element, l)] = shear_area_1;
            As2[view(element, l)] = shear_area_2;
        });
    });
}

std::pair<index_view, matrix> fem_submesh::tangent_stiffness(std::int32_t const element) const
{
    auto const& configuration = mesh_coordinates->initial_configuration(local_node_view(element));

    return {dof_list(Eigen::placeholders::all, element),
            bending_stiffness(configuration, element) + shear_stiffness(configuration, element)
                + axial_stiffness(configuration, element)
                + torsional_stiffness(configuration, element)};
}

matrix const& fem_submesh::bending_stiffness(matrix const& configuration,
                                             std::int32_t const element) const
{
    static thread_local matrix2x B_b(2, 6 * sf->nodes());
    static thread_local matrix k_b(6 * sf->nodes(), 6 * sf->nodes());

    B_b.setZero();
    k_b.setZero();

    sf->quadrature().integrate_inplace(k_b, [&, this](auto const& femval, auto const l) {
        auto const& [N, dN] = femval;

        double const jacobian = configuration.row(0) * dN;

        for (int i = 0; i < sf->nodes(); ++i)
        {
            auto const offset = i * 6;

            B_b(0, 3 + offset) = B_b(1, 4 + offset) = dN(i, l) / jacobian;
        }
        return B_b.transpose() * D_b * B_b * jacobian;
    });
    return k_b;
}

matrix const& fem_submesh::shear_stiffness(matrix const& configuration, std::int32_t const element) const
{
    static thread_local matrix2x B_s(2, 6 * sf->nodes());
    static thread_local matrix k_s(6 * sf->nodes(), 6 * sf->nodes());

    B_s.setZero();
    k_s.setZero();

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
        return B_s.transpose() * D_s * B_s * jacobian;
    });
    return k_s;
}

matrix const& fem_submesh::axial_stiffness(matrix const& configuration, std::int32_t const element) const
{
    static thread_local vector B_a(6 * sf->nodes());
    static thread_local matrix k_a(6 * sf->nodes(), 6 * sf->nodes());

    B_a.setZero();
    k_a.setZero();

    double A{1.0};

    sf->quadrature().integrate_inplace(k_a, [&, this](auto const& femval, auto const l) {
        auto const& [N, dN] = femval;

        double const jacobian = configuration.row(0) * dN;

        for (int i = 0; i < sf->nodes(); ++i)
        {
            auto const offset = i * 6;

            B_a(2 + offset) = dN(i, l) / jacobian;
        }
        return B_a * material.elastic_modulus() * A * B_a.transpose() * jacobian;
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

    auto const [I_1, I_2] = profiles[0]->second_moment_area();

    auto const J = I_1 + I_2;

    sf->quadrature().integrate_inplace(k_t, [&, this](auto const& femval, auto const l) {
        auto const& [N, dN] = femval;

        double const jacobian = configuration.row(0) * dN;

        for (int i = 0; i < sf->nodes(); ++i)
        {
            auto const offset = i * 6;

            B_t(5 + offset) = dN(i, l) / jacobian;
        }
        return B_t * material.shear_modulus() * J * B_t.transpose() * jacobian;
    });
    return k_t;
}
}
