
#pragma once

#include "ConstitutiveModel.hpp"

#include "numeric/DenseTypes.hpp"

#include "material/IsotropicElasticPlastic.hpp"

namespace neon
{
class HyperElasticPlastic : public ConstitutiveModel
{
public:
    HyperElasticPlastic(InternalVariables& variables) : ConstitutiveModel(variables) {}
};

class FiniteJ2Plasticity : public HyperElasticPlastic
{
public:
    FiniteJ2Plasticity(InternalVariables& variables, Json::Value const& material_data);

    ~FiniteJ2Plasticity();

    void update_internal_variables(double const time_step_size) override final;

    Material const& intrinsic_material() const override final { return material; }

    virtual bool is_finite_deformation() const override final { return true; };

protected:
    CMatrix elastic_moduli() const;

    CMatrix deviatoric_projection() const;

    CMatrix incremental_tangent(double const plastic_increment,
                                double const von_mises) const;

    CMatrix algorithmic_tangent(double const plastic_increment,
                                double const accumulated_plastic_strain,
                                double const von_mises,
                                Matrix3 const& n) const;

private:
    IsotropicElasticPlastic material;
    CMatrix const C_e = elastic_moduli();
    CMatrix const I_dev = deviatoric_projection();
};
}
