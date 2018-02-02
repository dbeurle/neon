
#include "IsotropicElasticPlasticDamage.hpp"

#include "Exceptions.hpp"

#include "io/json.hpp"

namespace neon
{
IsotropicElasticPlasticDamage::IsotropicElasticPlasticDamage(json const& material_data)
    : IsotropicElasticPlastic(material_data)
{
    if (!material_data.count("SofteningMultiplier"))
    {
        throw MaterialPropertyException("SofteningMultiplier");
    }
    if (!material_data.count("KinematicHardeningModulus"))
    {
        throw MaterialPropertyException("KinematicHardeningModulus");
    }
    if (!material_data.count("PlasticityViscousExponent"))
    {
        throw MaterialPropertyException("PlasticityViscousExponent");
    }
    if (!material_data.count("PlasticityViscousMultiplier"))
    {
        throw MaterialPropertyException("PlasticityViscousMultiplier");
    }
    if (!material_data.count("DamageViscousExponent"))
    {
        throw MaterialPropertyException("DamageViscousExponent");
    }
    if (!material_data.count("DamageViscousMultiplier"))
    {
        throw MaterialPropertyException("DamageViscousMultiplier");
    }

    gamma = material_data["SofteningMultiplier"];
    C = material_data["KinematicHardeningModulus"];
    np = material_data["PlasticityViscousExponent"];
    kp = material_data["PlasticityViscousMultiplier"];
    nd = material_data["DamageViscousExponent"];
    kd = material_data["DamageViscousMultiplier"];
    // TODO: may use "IsotropicHardeningModulus" : 400.0e6,
}
}
