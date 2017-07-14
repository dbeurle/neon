
#pragma once

#include "Hyperelastic.hpp"

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

    virtual void update_internal_variables(double const Δt) override;

    Material const& intrinsic_material() const override final { return material; };

protected:
    double volumetric_free_energy_derivative(double const J, double const bulk_modulus) const;
    double volumetric_free_energy_second_derivative(double const J, double const bulk_modulus) const;

    Matrix t_outer_t_outer_t_outer_t(Vector3 const& t) const;

    Matrix3 deviatoric_projection(double const pressure, Matrix3 const& τ_dev) const;

    Matrix deviatoric_projection(Matrix const& C_dev, Matrix3 const& τ_dev) const;

protected:
    LinearElastic material; //!< Elastic model where C1 = mu/2 and C2 = bulk-modulus / 2

    UnitSphereQuadrature unit_sphere;

    double number_of_chains;
    double const boltzmann_constant = 1.38064852e-23;
    double const temperature = 298.0;

    double μ; //!< Shear modulus

    double segments_per_chain = 0.0;

    double chain_decay_rate = 0.0;
};

inline double AffineMicrosphere::volumetric_free_energy_derivative(double const J,
                                                                   double const bulk_modulus) const
{
    return bulk_modulus / 2.0 * (J - 1.0 / J);
}
inline double AffineMicrosphere::volumetric_free_energy_second_derivative(double const J,
                                                                          double const bulk_modulus) const
{
    return bulk_modulus / 2.0 * (1.0 + 1.0 / std::pow(J, 2));
}
}
