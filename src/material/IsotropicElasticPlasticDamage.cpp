
#include "IsotropicElasticPlasticDamage.hpp"

#include "Exceptions.hpp"

#include <json/json.h>

namespace neon
{
IsotropicElasticPlasticDamage::IsotropicElasticPlasticDamage(Json::Value const& material_data)
    : IsotropicElasticPlastic(material_data)
{
    if (material_data["SofteningMultiplier"].empty())
        throw MaterialPropertyException("SofteningMultiplier");
    if (material_data["HardeningModulus"].empty())
        throw MaterialPropertyException("HardeningModulus");
    if (material_data["PlasticityViscousExponent"].empty())
        throw MaterialPropertyException("PlasticityViscousExponent");
    if (material_data["PlasticityViscousMultiplier"].empty())
        throw MaterialPropertyException("PlasticityViscousMultiplier");
    if (material_data["DamageViscousExponent"].empty())
        throw MaterialPropertyException("DamageViscousExponent");
    if (material_data["DamageViscousMultiplier"].empty())
        throw MaterialPropertyException("DamageViscousMultiplier");

    gamma = material_data["SofteningMultiplier"].asDouble();
    C = material_data["HardeningModulus"].asDouble();
    np = material_data["PlasticityViscousExponent"].asDouble();
    kp = material_data["PlasticityViscousMultiplier"].asDouble();
    nd = material_data["DamageViscousExponent"].asDouble();
    kd = material_data["DamageViscousMultiplier"].asDouble();
    // TODO: may use "IsotropicHardeningModulus" : 400.0e6,
}

} // namespace neon
