
#include "constitutive/mechanical/beam/isotropic_linear.hpp"

#include "io/json.hpp"

#include <tbb/parallel_for.h>

namespace neon::mechanical::beam
{
isotropic_linear::isotropic_linear(std::shared_ptr<variable_type>& variables,
                                   json const& material_data)
    : variables(variables), material(material_data)
{
    // Allocate the four constitutive values for axial, bending, shear and torsion.
    variables->add(variable_type::scalar::torsional_stiffness,
                   variable_type::scalar::axial_stiffness,
                   variable_type::second::bending_stiffness,
                   variable_type::second::shear_stiffness);
}

void isotropic_linear::update_internal_variables()
{
    // Stiffness components
    auto [D_b, D_s] = variables->get(variable_type::second::bending_stiffness,
                                     variable_type::second::shear_stiffness);

    auto [D_a, D_t] = variables->get(variable_type::scalar::axial_stiffness,
                                     variable_type::scalar::torsional_stiffness);

    // Access the geometric properties at the quadrature points
    auto [shear_area_1,
          shear_area_2,
          cross_sectional_area] = variables->get(variable_type::scalar::shear_area_1,
                                                 variable_type::scalar::shear_area_2,
                                                 variable_type::scalar::cross_sectional_area);
    auto [second_moment_area_1,
          second_moment_area_2] = variables->get(variable_type::scalar::second_moment_area_1,
                                                 variable_type::scalar::second_moment_area_2);

    tbb::parallel_for(std::size_t{}, variables->entries(), [&](auto const l) {
        // Polar moment of area
        double const J = second_moment_area_1[l] + second_moment_area_2[l];

        D_b[l] << material.elastic_modulus() * second_moment_area_1[l], 0.0, 0.0,
            material.elastic_modulus() * second_moment_area_2[l];

        D_s[l] << material.shear_modulus() * shear_area_1[l], 0.0, 0.0,
            material.shear_modulus() * shear_area_2[l];

        D_a[l] = material.elastic_modulus() * cross_sectional_area[l];

        D_t[l] = J * material.shear_modulus();
    });
}

void isotropic_linear::update_stress()
{
    // Compute the entire stress tensor?
}
}
