
#pragma once

#include "ConstitutiveModel.hpp"

#include "numeric/DenseTypes.hpp"

#include "material/IsotropicElasticPlastic.hpp"

#include <memory>

namespace neon
{
class HyperElasticPlastic : public ConstitutiveModel
{
public:
    HyperElasticPlastic(InternalVariables& variables) : ConstitutiveModel(variables) {}
};

class J2Plasticity : public HyperElasticPlastic
{
public:
    J2Plasticity(InternalVariables& variables, Json::Value const& material_data);

    ~J2Plasticity();

    void update_internal_variables(double const Δt) override final;

    Material const& intrinsic_material() const override final { return material; }

protected:
    CMatrix elastic_moduli() const;

    CMatrix deviatoric_projection() const;

    CMatrix transformJaumannToTruesdellKirchoff(CMatrix const C_τ_J,
                                                Matrix3 const& σ,
                                                double const J) const;

    CMatrix incremental_tangent(double const Δλ, double const von_mises) const;

    CMatrix continuum_tangent(double const Δλ,
                              double const α,
                              double const von_mises,
                              Matrix3 const& n) const;

    CMatrix algorithmic_tangent(double const α,
                                Matrix3 const& σ,
                                Matrix3 const& σ_0,
                                Matrix3 const& normal,
                                double const J) const;

private:
    IsotropicElasticPlastic material;
    CMatrix const C_e = elastic_moduli();
    CMatrix const Id = deviatoric_projection();
};
}
