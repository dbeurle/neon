
#pragma once

#include "numeric/tensor_operations.hpp"

namespace neon::mechanical
{
/**
 * Evaluates the yield function of the J2 plasticity surface
 * @return greater than zero if the yield function has been violated
 */
template <class material_type>
[[nodiscard]] double evaluate_J2_yield_function(material_type const& material,
                                                double const von_mises_stress,
                                                double const accumulated_plastic_strain,
                                                double const plastic_increment = 0.0) {
    return von_mises_stress - 3.0 * material.shear_modulus() * plastic_increment
           - material.yield_stress(accumulated_plastic_strain + plastic_increment);
}

/**
 * algorithm_tangent computes the consistent (algorithmic) tangent material
 * operator for the J2 plasticity.  This operator has been linearised
 * consistently with the backward Euler method for the radial return method
 * (backward Euler).  This function is valid for the three-dimensional and the
 * plane strain case.  The plane stress can requires an additional projection.
 */
template <typename tangent_operator_type, typename tensor_type>
tangent_operator_type algorithmic_tangent(double const G,
                                          double const H,
                                          double const plastic_increment,
                                          double const von_mises_stress,
                                          tensor_type const& normal,
                                          tangent_operator_type const& I_dev,
                                          tangent_operator_type const& C_elastic)
{
    return C_elastic - 6.0 * plastic_increment * std::pow(G, 2) / von_mises_stress * I_dev
           + 6.0 * std::pow(G, 2) * (plastic_increment / von_mises_stress - 1.0 / (3.0 * G + H))
                 * outer_product(normal, normal);
}
}
