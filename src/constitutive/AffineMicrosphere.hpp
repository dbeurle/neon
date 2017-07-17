
#pragma once

#include "Hyperelastic.hpp"

#include "material/LinearElastic.hpp"
#include "numeric/Tensor.hpp"
#include "quadrature/UnitSphereQuadrature.hpp"

#include <json/forwards.h>

namespace neon
{
class AffineMicrosphere : public Hyperelastic
{
public:
    /**
     * @param variables Reference to internal state variable store
     * @param data Json object with material data
     */
    explicit AffineMicrosphere(InternalVariables& variables,
                               Json::Value const& material_data);

    virtual void update_internal_variables(double const Δt) override;

    Material const& intrinsic_material() const override final { return material; };

protected:
    /**
     * \f{align*}{
         U' &= \frac{\partial U}{\partial J} = \frac{K}{2}\left(J - \frac{1}{J}\right)
       \f}
     * where
     * \f{align*}{
         U &= \frac{K}{4}(J^2 - 1) - \frac{K}{2}\ln{J}
       \f}
     */
    double volumetric_free_energy_derivative(double const J,
                                             double const bulk_modulus) const;

    /**
     * \f{align*}{
         U'' &= \frac{\partial^2 U}{\partial J^2} = \frac{K}{2} \left(1 +
     \frac{1}{J^2}\right) \f}
     * \sa volumetric_free_energy_derivative
     */
    double volumetric_free_energy_second_derivative(double const J,
                                                    double const bulk_modulus) const;

    /**
     *\f{align*}{
     * \boldsymbol{\tau} &= p \boldsymbol{g}^{-1} + \mathbb{P} : \bar{\boldsymbol{\tau}}
     * \f}
     */
    Matrix3 deviatoric_projection(double const pressure, Matrix3 const& τ_dev) const;

    /**
     *\f{align*}{
        \boldsymbol{C} &= \mathbb{P} : \left[\bar{\mathbb{C}} +
     \frac{2}{3}(\boldsymbol{\tau} : \boldsymbol{g}) \mathbb{I} - \frac{2}{3}
     (\bar{\boldsymbol{\tau}} \otimes \boldsymbol{g}^{-1} + \boldsymbol{g}^{-1} \otimes
     \bar{\boldsymbol{\tau}}) \right] : \mathbb{P} \f}
     */
    Matrix deviatoric_projection(Matrix const& C_dev, Matrix3 const& τ_dev) const;

protected:
    LinearElastic material; //!< Elastic model where C1 = mu/2 and C2 = bulk-modulus / 2

    UnitSphereQuadrature unit_sphere;

    Matrix const IoI = I_outer_I();
    Matrix const I = fourth_order_identity();

    double number_of_chains;
    double constexpr boltzmann_constant = 1.38064852e-23;
    double constexpr temperature = 298.0;

    double μ; //!< Shear modulus

    double segments_per_chain = 0.0;

    double chain_decay_rate = 0.0;
};

inline double AffineMicrosphere::volumetric_free_energy_derivative(double const J,
                                                                   double const bulk_modulus) const
{
    return bulk_modulus / 2.0 * (J - 1.0 / J);
}
inline double AffineMicrosphere::volumetric_free_energy_second_derivative(
    double const J,
    double const bulk_modulus) const
{
    return bulk_modulus / 2.0 * (1.0 + 1.0 / std::pow(J, 2));
}
}
