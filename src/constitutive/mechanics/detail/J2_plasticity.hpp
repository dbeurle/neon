
#pragma once

#include "numeric/tensor_operations.hpp"

namespace neon::mechanics
{
/**
 * \fn evaluate_J2_yield_function Evaluates the yield function of the J2 plasticity surface
 *
 * \param material The material class providing a shear_modulus and yield_stress method
 * \param von_mises_stress Current von Mises stress
 * \param accumulated_plastic_strain Current level of accumulated plastic strain
 * \param plastic_increment The current plastic increment (default 0.0)
 *
 * @return Greater than zero if the yield function has been violated
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
 * \fn algorithm_tangent computes the consistent (algorithmic) tangent material
 * operator for the J2 plasticity.  This operator has been linearised
 * consistently with the backward Euler method for the radial return method
 * (backward Euler).  This function is valid for the three-dimensional and the
 * plane strain case.  The plane stress can requires an additional projection.
 *
 * \param G Material shear modulus
 * \param H Material hardening modulus
 * \param plastic_increment
 * \param von_mises_stress
 * \param normal Tensor normal to the surface
 * \param I_dev Symmetric deviatoric fourth order tensor
 * \param C_elastic Elastic material tangent operator
 *
 * @return Consistently linearised material tangent operator
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

/**
 * Apply the finite strain correction due to the geometry stiffness matrix.
 * It returns the value
 * \f{align*}{
     \sigma_{il} \delta_{jk}
   \f}
 * for a plane strain formulation.
 */
[[nodiscard]] inline matrix3 finite_strain_correction(matrix2 const& s)
{
    matrix3 A(3, 3);
    // clang-format off
    A << s(0, 0),     0.0, s(0, 1),
             0.0, s(1, 1),     0.0,
         s(1, 0),     0.0,     0.0;
    // clang-format on
    return A;
}

/**
 * Apply the finite strain correction due to the geometry stiffness matrix.
 * It returns the value
 * \f{align*}{
     \sigma_{il} \delta_{jk}
   \f}
 * for a three-dimensional continuum formulation
 */
[[nodiscard]] inline matrix6 finite_strain_correction(matrix3 const& s)
{
    matrix6 H(6, 6);
    // clang-format off
    H << 2 * s(0, 0),           0,           0,           0, 2 * s(0, 2), 2 * s(0, 1),
                   0, 2 * s(1, 1),           0, 2 * s(1, 2),           0,           0,
                   0,           0, 2 * s(2, 2),           0,           0,           0,
                   0, 2 * s(1, 2),           0,     s(2, 2),           0,           0,
         2 * s(0, 2),           0,           0,           0,     s(2, 2),     s(2, 1),
         2 * s(0, 1),           0,           0,           0,     s(1, 2),     s(2, 2);
    // clang-format on
    return H;
}

/**
 * This computes the fourth order tensor B in Voigt notation for the plane strain
 * formulation from \cite Neto2011 on page 598.
 * \f{align*}{
     B_{ijkl} &= \delta_{ik}(\boldsymbol{B}^{e, trial}_{jl}) + \delta_{jk}(\boldsymbol{B}^{e,
 trial}_{il}) \f}
 */
[[nodiscard]] inline matrix3 finite_strain_B_operator(matrix2 const& Be_trial)
{
    matrix3 A(3, 3);
    // clang-format off
    A << 2 * Be_trial(0, 0),                  0, 2 * Be_trial(0, 1),
                          0, 2 * Be_trial(1, 1),                  0,
         2 * Be_trial(1, 0),                  0,     Be_trial(2, 2);
    // clang-format on
    return A;
}

/**
 * This computes the fourth order tensor B in Voigt notation for the three
 * dimensional formulation from \cite Neto2011 on page 598.
 */
[[nodiscard]] inline matrix6 finite_strain_B_operator(matrix3 const& Be_trial)
{
    matrix6 B(6, 6);
    // clang-format off
    B << 2 * Be_trial(0, 0),                  0,                  0,                  0, 2 * Be_trial(0, 2), 2 * Be_trial(0, 1),
                          0, 2 * Be_trial(1, 1),                  0, 2 * Be_trial(1, 2),                  0,                  0,
                          0,                  0, 2 * Be_trial(2, 2),                  0,                  0,                  0,
                          0, 2 * Be_trial(1, 2),                  0,     Be_trial(2, 2),                  0,                  0,
         2 * Be_trial(0, 2),                  0,                  0,                  0,     Be_trial(2, 2),      Be_trial(2, 1),
         2 * Be_trial(0, 1),                  0,                  0,                  0,     Be_trial(2, 1),      Be_trial(1, 1);
    // clang-format on
    return B;
}
}
