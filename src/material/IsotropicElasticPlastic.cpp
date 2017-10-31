
#include "IsotropicElasticPlastic.hpp"

#include "Exceptions.hpp"

#include <json/json.h>

namespace neon
{
IsotropicElasticPlastic::IsotropicElasticPlastic(Json::Value const& material_data)
    : LinearElastic(material_data)
{
    if (material_data["YieldStress"].empty()) throw MaterialPropertyException("YieldStress");

    stress_y = material_data["YieldStress"].asDouble();

    if (!material_data["IsotropicHardeningModulus"].empty())
    {
        H = material_data["IsotropicHardeningModulus"].asDouble();
    }

    if (!material_data["IsotropicKinematicModulus"].empty())
    {
        K = material_data["IsotropicKinematicModulus"].asDouble();
    }
}

double IsotropicElasticPlastic::yield_stress(double const effective_strain) const
{
    return stress_y + effective_strain * H;
}

double IsotropicElasticPlastic::hardening_modulus(double const effective_strain) const { return H; }

double IsotropicElasticPlastic::kinematic_modulus(double const effective_strain) const { return K; }
}
