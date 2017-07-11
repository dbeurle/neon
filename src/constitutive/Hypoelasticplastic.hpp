
#pragma once

#include "ConstitutiveModel.hpp"

#include "numeric/DenseTypes.hpp"

#include "material/IsotropicElasticPlastic.hpp"

#include <memory>

namespace neon
{
class Hypoelasticplastic : public ConstitutiveModel
{
public:
    Hypoelasticplastic(InternalVariables& variables) : ConstitutiveModel(variables) {}
};

class J2Plasticity : public Hypoelasticplastic
{
public:
    J2Plasticity(InternalVariables& variables, Json::Value const& material_data);

    ~J2Plasticity();

    void update_internal_variables(double const Δt) override final;

    Material const& intrinsic_material() const override final { return material; }

protected:
    Matrix elastic_moduli() const;

    Matrix transformJaumannToTruesdellKirchoff(Matrix const C_τ_J,
                                               Matrix3 const& σ,
                                               double const J) const;

    Matrix outer_product(Matrix3 const& n) const;

    Matrix algorithmic_tangent(double const α,
                               Matrix3 const& σ,
                               Matrix3 const& σ_0,
                               Matrix3 const& normal,
                               double const J) const;

private:
    IsotropicElasticPlastic material;
    Matrix C_e;
};
}
