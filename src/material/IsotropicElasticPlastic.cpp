
#include "IsotropicElasticPlastic.hpp"

#include "Exceptions.hpp"

#include <json/value.h>

namespace neon
{
IsotropicElasticPlastic::IsotropicElasticPlastic(Json::Value const& material_data)
    : LinearElastic(material_data)
{
    if (!material_data.isMember("YieldStress"))
    {
        throw MaterialPropertyException("YieldStress");
    }

    stress_y = material_data["YieldStress"].asDouble();

    if (material_data.isMember("IsotropicHardeningModulus"))
    {
        H_iso = material_data["IsotropicHardeningModulus"].asDouble();
    }

    if (material_data.isMember("IsotropicKinematicModulus"))
    {
        K_iso = material_data["IsotropicKinematicModulus"].asDouble();
    }
}

double IsotropicElasticPlastic::yield_stress(double const effective_strain) const
{
    return stress_y + effective_strain * H_iso;
}

double IsotropicElasticPlastic::hardening_modulus(double const effective_strain) const
{
    return H_iso;
}

double IsotropicElasticPlastic::kinematic_modulus(double const effective_strain) const
{
    return K_iso;
}
}
