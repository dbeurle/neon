
#pragma once

#include "ConstitutiveModel.hpp"

#include <json/forwards.h>

#include "material/LinearElastic.hpp"
#include "quadrature/UnitSphereQuadrature.hpp"

namespace neon
{
class Hyperelastic : public ConstitutiveModel
{
public:
    Hyperelastic(InternalVariables& variables);
};

class NeoHooke : public Hyperelastic
{
public:
    /**
     * @param variables Reference to internal state variable store
     * @param data Json object with material data
     */
    explicit NeoHooke(InternalVariables& variables, Json::Value const& material_data);

    void update_internal_variables() override final;

    void update_continuum_tangent() override final;

    Material const& intrinsic_material() const override final { return material; };

private:
    LinearElastic material; //!< Elastic model where C1 = mu/2 and C2 = bulk-modulus / 2
};

class AffineMicrosphere : public Hyperelastic
{
public:
    /**
     * @param variables Reference to internal state variable store
     * @param data Json object with material data
     */
    explicit AffineMicrosphere(InternalVariables& variables, Json::Value const& material_data);

    void update_internal_variables() override final;

    void update_continuum_tangent() override final;

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

    double segments_per_chain = 0.0;
};
}
