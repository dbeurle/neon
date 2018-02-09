
#pragma once

#include "constitutive/ConstitutiveModel.hpp"

#include "material/micromechanical_elastomer.hpp"
#include "numeric/tensor_operations.hpp"
#include "quadrature/unit_sphere_quadrature.hpp"

#include <json/forwards.h>

namespace neon::mechanical::solid
{
/**
 * \ingroup Hyperelastic
 * \addtogroup Hyperelastic
 * \{
 *
 * affine_microsphere is responsible for computing the Cauchy stress and the
 * material tangent in implicit methods.  The affine microsphere model \cite Miehe2004
 * is used to model elastomer materials using micromechanical motivations and
 * homogenises the force from a single chain over a unit sphere.
 *
 * This constitutive model requires the use of a quadrature scheme for the unit
 * sphere and this internal variable update can be computationally expensive and
 * is therefore multithreaded.
 */
class affine_microsphere : public ConstitutiveModel
{
public:
    /**
     * @param variables Reference to internal state variable store
     * @param material_data Json object with input file material data
     */
    explicit affine_microsphere(std::shared_ptr<InternalVariables>& variables,
                                json const& material_data,
                                unit_sphere_quadrature::Rule const rule);

    virtual void update_internal_variables(double const time_step_size) override;

    [[nodiscard]] material_property const& intrinsic_material() const override final { return material; };

    [[nodiscard]] virtual bool is_finite_deformation() const override final { return true; };

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
    [[nodiscard]] double volumetric_free_energy_dJ(double const J, double const bulk_modulus) const;

    /**
     * \f{align*}{
         U'' &= \frac{\partial^2 U}{\partial J^2} = \frac{K}{2} \left(1 +
     \frac{1}{J^2}\right) \f}
     * \sa volumetric_free_energy_dJ
     */
    [[nodiscard]] double volumetric_free_energy_second_d2J(double const J,
                                                           double const bulk_modulus) const;

    /**
     * Compute the Padé approximation of the inverse Langevin stretch model
     * \f{align*}{
         n \psi_f^{'}(\lambda) &= \frac{3N - \lambda^2}{N - \lambda^2}
       \f}
     */
    [[nodiscard]] double pade_first(double const micro_stretch, double const N) const;

    /**
     * Compute the Padé approximation of the inverse Langevin stretch model
     * \f{align*}{
         n \psi_f^{''}(\lambda) &= \frac{\lambda^4 + 3N^2}{(N - \lambda^2)^2}
       \f}
     */
    [[nodiscard]] double pade_second(double const micro_stretch, double const N) const;

    /**
     * Compute the Kirchhoff stress using the deviatoric projection of the
     * macro stress tensor according to
     * \f{align*}{
       \boldsymbol{\tau} &= p \boldsymbol{g}^{-1} + \mathbb{P} : \bar{\boldsymbol{\tau}} \f}
     * @param pressure Hydrostatic pressure
     * @param macro_stress Stress tensor from unit sphere homogenisation
     */
    [[nodiscard]] matrix3 compute_kirchhoff_stress(double const pressure,
                                                   matrix3 const& macro_stress) const;

    /**
     * Compute the material tangent including the sdeviatoric projection defined as
     *\f{align*}{
        \boldsymbol{C} &= (\kappa + p) \mathbf{g}^{-1} \otimes \mathbf{g}^{-1} - 2p\mathbb{I} +
     \mathbb{P} : \left[\bar{\mathbb{C}} + \frac{2}{3}(\boldsymbol{\tau} : \boldsymbol{g})
     \mathbb{I} - \frac{2}{3}
     (\bar{\boldsymbol{\tau}} \otimes \boldsymbol{g}^{-1} + \boldsymbol{g}^{-1} \otimes
     \bar{\boldsymbol{\tau}}) \right] : \mathbb{P} \f}
     * @param J Determinant of deformation gradient
     * @param K Shear modulus
     * @param macro_C Macromoduli from unit sphere
     * @param macro_stress Macrostress from unit sphere
     */
    [[nodiscard]] matrix6 compute_material_tangent(double const J,
                                                   double const K,
                                                   matrix6 const& macro_C,
                                                   matrix3 const& macro_stress) const;

    /**
     * Compute the macro stress using the unit sphere homogenisation
     * technique for a given F and N and perform the deviatoric projection
     * @param F_unimodular Unimodular decomposition of the deformation gradient
     * @param bulk_modulus The material bulk modulus
     * @param N number of segments per chain
     * @return Kirchhoff stress tensor
     */
    [[nodiscard]] matrix3 compute_macro_stress(matrix3 const& F_unimodular,
                                               double const bulk_modulus,
                                               double const N) const;

    /**
     * Compute the material tangent matrix using the unit sphere homogenisation
     * technique for a given F and N
     * @param F_unimodular Unimodular decomposition of the deformation gradient
     * @param bulk_modulus The material bulk modulus
     * @param N number of segments per chain
     * @return Macromoduli from unit sphere homogenisation
     */
    [[nodiscard]] matrix6 compute_macro_moduli(matrix3 const& F_unimodular,
                                               double const bulk_modulus,
                                               double const N) const;

    /**
     * Compute the deformed tangent using the unimodular deformation gradient
     * and the vector associated with the quadrature point on the unit sphere
     */
    [[nodiscard]] vector3 deformed_tangent(matrix3 const& F_unimodular,
                                           vector3 const& surface_vector) const
    {
        return F_unimodular * surface_vector;
    }

    /** Compute the microstretch, which is the norm of the deformed tangent vector */
    [[nodiscard]] auto compute_microstretch(vector3 const& deformed_tangent) const
    {
        return deformed_tangent.norm();
    }

protected:
    unit_sphere_quadrature unit_sphere; //!< Unit sphere quadrature rule

    matrix6 const IoI = voigt::I_outer_I();                      //!< Outer product
    matrix6 const I = voigt::kinematic::fourth_order_identity(); //!< Fourth order identity
    matrix6 const P = voigt::kinetic::deviatoric();              //!< Deviatoric fourth order tensor
private:
    micromechanical_elastomer material; //!< Material with micromechanical parameters
};

/** \} */

inline double affine_microsphere::volumetric_free_energy_dJ(double const J,
                                                            double const bulk_modulus) const
{
    return bulk_modulus / 2.0 * (J - 1.0 / J);
}

inline double affine_microsphere::volumetric_free_energy_second_d2J(double const J,
                                                                    double const bulk_modulus) const
{
    return bulk_modulus / 2.0 * (1.0 + 1.0 / std::pow(J, 2));
}

inline double affine_microsphere::pade_first(double const micro_stretch, double const N) const
{
    return (3.0 * N - std::pow(micro_stretch, 2)) / (N - std::pow(micro_stretch, 2));
}

inline double affine_microsphere::pade_second(double const micro_stretch, double const N) const
{
    return (std::pow(micro_stretch, 4) + 3.0 * std::pow(N, 2))
           / std::pow(N - std::pow(micro_stretch, 2), 2);
}

/**
 * \ingroup Hyperelastic
 * \addtogroup Hyperelastic
 * \{
 *
 * AffineMicrosphereWithDegradation is responsible for computing the Cauchy
 * stress and the material tangent in implicit methods when ageing is present.
 * The affine microsphere model is used to model elastomer materials using
 * micromechanical motivations and homogenises the force from a single chain
 * over a unit sphere.
 */
class AffineMicrosphereWithDegradation : public affine_microsphere
{
public:
    /**
     * @param variables Reference to internal state variable store
     * @param material_data Json object with input file material data
     */
    explicit AffineMicrosphereWithDegradation(std::shared_ptr<InternalVariables>& variables,
                                              json const& material_data,
                                              unit_sphere_quadrature::Rule const rule);

    virtual void update_internal_variables(double const time_step_size) override;

protected:
    template <typename MatrixTp, typename Functor>
    MatrixTp weighting(std::vector<double> const& G, MatrixTp accumulator, Functor&& f) const;

private:
    stochastic_micromechanical_elastomer material; //!< Material with micromechanical parameters
};
/** \} */

template <typename MatrixTp, typename Functor>
inline MatrixTp AffineMicrosphereWithDegradation::weighting(std::vector<double> const& G,
                                                            MatrixTp accumulator,
                                                            Functor&& f) const
{
    for (auto i = 0; i < material.segment_groups().size(); i++)
    {
        accumulator.noalias() += f(material.segment_groups()[i]) * G[i];
    }
    return accumulator;
}
}
