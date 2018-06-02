
#include "isotropic_elastic_plastic_damage.hpp"

#include "exceptions.hpp"

#include "io/json.hpp"

namespace neon
{
isotropic_elastic_plastic_damage::isotropic_elastic_plastic_damage(json const& material_data)
    : isotropic_elastic_plastic(material_data)
{
    if (!material_data.count("SofteningMultiplier"))
    {
        throw std::domain_error("\"SofteningMultiplier\" needs to be specified");
    }
    if (!material_data.count("KinematicHardeningModulus"))
    {
        throw std::domain_error("\"KinematicHardeningModulus\" needs to be specified");
    }
    if (!material_data.count("PlasticityViscousExponent"))
    {
        throw std::domain_error("\"PlasticityViscousExponent\" needs to be specified");
    }
    if (!material_data.count("PlasticityViscousDenominator"))
    {
        throw std::domain_error("\"PlasticityViscousDenominator\" needs to be specified");
    }
    if (!material_data.count("DamageViscousExponent"))
    {
        throw std::domain_error("\"DamageViscousExponent\" needs to be specified");
    }
    if (!material_data.count("DamageViscousDenominator"))
    {
        throw std::domain_error("\"DamageViscousDenominator\" needs to be specified");
    }

    gamma = material_data["SofteningMultiplier"];
    C = material_data["KinematicHardeningModulus"];
    np = material_data["PlasticityViscousExponent"];
    sp = material_data["PlasticityViscousDenominator"];
    nd = material_data["DamageViscousExponent"];
    sd = material_data["DamageViscousDenominator"];
}
}
