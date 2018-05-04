
#include "constitutive/mechanical/beam/isotropic_linear.hpp"

#include "io/json.hpp"

#include <tbb/parallel_for.h>

namespace neon::mechanical::beam
{
isotropic_linear::isotropic_linear(std::shared_ptr<internal_variable_type>& variables,
                                   json const& material_data)
    : variables(variables), material(material_data)
{
    // Allocate the four constitutive values for axial, bending, shear and torsion.
    variables->add(internal_variable_type::scalar::torsional_stiffness,
                   internal_variable_type::scalar::axial_stiffness,
                   internal_variable_type::second::bending_stiffness,
                   internal_variable_type::second::shear_stiffness);
}

void isotropic_linear::update_internal_variables()
{
    auto [D_b, D_s] = variables->fetch(internal_variable_type::second::bending_stiffness,
                                       internal_variable_type::second::shear_stiffness);

    auto [D_a, D_t] = variables->fetch(internal_variable_type::scalar::axial_stiffness,
                                       internal_variable_type::scalar::torsional_stiffness);

    // Access the geometric properties at the quadrature points
    auto const [shear_area_1,
                shear_area_2,
                cross_sectional_area] = variables
                                            ->fetch(internal_variable_type::scalar::shear_area_1,
                                                    internal_variable_type::scalar::shear_area_2,
                                                    internal_variable_type::scalar::cross_sectional_area);
    auto const [second_moment_area_1,
                second_moment_area_2] = variables
                                            ->fetch(internal_variable_type::scalar::second_moment_area_1,
                                                    internal_variable_type::scalar::second_moment_area_2);

    tbb::parallel_for(std::size_t{0}, D_b.size(), [&](auto const l) {
        // Section areas
        auto const A_1 = shear_area_1[l];
        auto const A_2 = shear_area_2[l];
        auto const A = cross_sectional_area[l];

        // Section second moment of area
        auto const I_1 = second_moment_area_1[l];
        auto const I_2 = second_moment_area_2[l];

        // Polar moment of area
        auto const J = I_1 + I_2;

        D_b[l] << material.elastic_modulus() * I_1, 0.0, 0.0, material.elastic_modulus() * I_2;
        D_s[l] << material.shear_modulus() * A_1, 0.0, 0.0, material.shear_modulus() * A_2;

        D_a[l] = material.elastic_modulus() * A;
        D_t[l] = J * material.shear_modulus();
    });
}

void isotropic_linear::update_stress()
{
    // Compute the entire stress tensor?
}
}
