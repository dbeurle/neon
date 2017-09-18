
#pragma once

#include "Hyperelastic.hpp"

#include "material/MicromechanicalElastomer.hpp"
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
    explicit AffineMicrosphere(InternalVariables& variables, Json::Value const& material_data);

    virtual void update_internal_variables(double const time_step_size) override;

    Material const& intrinsic_material() const override final { return material; };

    virtual bool is_finite_deformation() const override final { return true; };

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
    double volumetric_free_energy_dJ(double const J, double const bulk_modulus) const;

    /**
     * \f{align*}{
         U'' &= \frac{\partial^2 U}{\partial J^2} = \frac{K}{2} \left(1 +
     \frac{1}{J^2}\right) \f}
     * \sa volumetric_free_energy_dJ
     */
    double volumetric_free_energy_second_d2J(double const J, double const bulk_modulus) const;

    /**
     * Compute the Padé approximation of the inverse Langevin stretch model
     * \f{align*}{
         n \psi_f^{'}(\lambda) &= \frac{3N - \lambda^2}{N - \lambda^2}
       \f}
     */
    double pade_first(double const micro_stretch, double const N) const;

    /**
     * Compute the Padé approximation of the inverse Langevin stretch model
     * \f{align*}{
         n \psi_f^{''}(\lambda) &= \frac{\lambda^4 + 3N^2}{(N - \lambda^2)^2}
       \f}
     */
    double pade_second(double const micro_stretch, double const N) const;

    /**
     *\f{align*}{
      \boldsymbol{\tau} &= p \boldsymbol{g}^{-1} + \mathbb{P} : \bar{\boldsymbol{\tau}} \f}
     */
    Matrix3 deviatoric_projection(double const pressure, Matrix3 const& stress_dev) const;

    /**
     *\f{align*}{
        \boldsymbol{C} &= \mathbb{P} : \left[\bar{\mathbb{C}} +
     \frac{2}{3}(\boldsymbol{\tau} : \boldsymbol{g}) \mathbb{I} - \frac{2}{3}
     (\bar{\boldsymbol{\tau}} \otimes \boldsymbol{g}^{-1} + \boldsymbol{g}^{-1} \otimes
     \bar{\boldsymbol{\tau}}) \right] : \mathbb{P} \f}
     */
    CMatrix deviatoric_projection(CMatrix const& C_dev, Matrix3 const& stress_dev) const;

    /**
     * Compute the Kirchhoff stress using the unit sphere homogenisation
     * technique for a given F and N
     * @param unimodular_F Unimodular decomposition of the deformation gradient
     * @param N number of segments per chain
     * @return Kirchhoff stress tensor
     */
    Matrix3 compute_kirchhoff_stress(Matrix3 const& unimodular_F, double const N) const;

    /**
     * Compute the material tangent matrix using the unit sphere homogenisation
     * technique for a given F and N
     * @param unimodular_F Unimodular decomposition of the deformation gradient
     * @param N number of segments per chain
     * @return Kirchhoff stress tensor
     */
    CMatrix compute_material_matrix(Matrix3 const& unimodular_F, double const N) const;

    template <typename MatrixTp, typename Functor>
    MatrixTp weighting(std::vector<double> const& G, MatrixTp accumulator, Functor f) const;

protected:
    MicromechanicalElastomer material;

    UnitSphereQuadrature unit_sphere;

    CMatrix const IoI = I_outer_I();
    CMatrix const I = fourth_order_identity();
};

inline double AffineMicrosphere::volumetric_free_energy_dJ(double const J,
                                                           double const bulk_modulus) const
{
    return bulk_modulus / 2.0 * (J - 1.0 / J);
}

inline double AffineMicrosphere::volumetric_free_energy_second_d2J(double const J,
                                                                   double const bulk_modulus) const
{
    return bulk_modulus / 2.0 * (1.0 + 1.0 / std::pow(J, 2));
}

inline double AffineMicrosphere::pade_first(double const micro_stretch, double const N) const
{
    return (3.0 * N - std::pow(micro_stretch, 2)) / (N - std::pow(micro_stretch, 2));
}

inline double AffineMicrosphere::pade_second(double const micro_stretch, double const N) const
{
    return (std::pow(micro_stretch, 4) + 3 * std::pow(N, 2))
           / std::pow(N - std::pow(micro_stretch, 2), 2);
}

template <typename MatrixTp, typename Functor>
inline MatrixTp AffineMicrosphere::weighting(std::vector<double> const& G,
                                             MatrixTp accumulator,
                                             Functor f) const
{
    for (int i = 0; i < material.segment_groups().size(); i++)
    {
        accumulator.noalias() += f(material.segment_groups()[i]) * G[i];
    }
    return accumulator;
}
}
