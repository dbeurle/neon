
#pragma once

#include "HypoElasticPlastic.hpp"

#include "numeric/DenseTypes.hpp"

#include "material/IsotropicElasticPlastic.hpp"

namespace neon
{
class J2Plasticity : public HypoElasticPlastic
{
public:
    J2Plasticity(InternalVariables& variables, Json::Value const& material_data);

    ~J2Plasticity();

    void update_internal_variables(double const time_step_size) override;

    Material const& intrinsic_material() const override { return material; }

    virtual bool is_finite_deformation() const override { return false; };

protected:
    Matrix3 compute_cauchy_stress(Matrix3 const& elastic_strain) const;

    CMatrix elastic_moduli() const;

    CMatrix deviatoric_projection() const;

    CMatrix algorithmic_tangent(double const plastic_increment,
                                double const accumulated_plastic_strain,
                                double const von_mises,
                                Matrix3 const& n) const;

protected:
    IsotropicElasticPlastic material;
    CMatrix const C_e = elastic_moduli();
    CMatrix const I_dev = deviatoric_projection();
};

inline Matrix3 J2Plasticity::compute_cauchy_stress(Matrix3 const& elastic_strain) const
{
    auto const G = material.shear_modulus();
    auto const lambda_e = material.lambda();
    return lambda_e * elastic_strain.trace() * Matrix3::Identity() + 2.0 * G * elastic_strain;
}
}
