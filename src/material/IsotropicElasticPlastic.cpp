
#include "IsotropicElasticPlastic.hpp"

#include "Exceptions.hpp"

#include "io/json.hpp"

namespace neon
{
IsotropicElasticPlastic::IsotropicElasticPlastic(json const& material_data)
    : LinearElastic(material_data)
{
    if (!material_data.count("YieldStress"))
    {
        throw std::domain_error("\"YieldStress\" was not specified in material properties");
    }

    stress_y = material_data["YieldStress"];

    if (material_data.count("IsotropicHardeningModulus"))
    {
        H_iso = material_data["IsotropicHardeningModulus"];
    }

    if (material_data.count("IsotropicKinematicModulus"))
    {
        K_iso = material_data["IsotropicKinematicModulus"];
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
