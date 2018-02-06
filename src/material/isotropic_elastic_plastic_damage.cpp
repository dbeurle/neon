
#include "isotropic_elastic_plastic_damage.hpp"

#include "Exceptions.hpp"

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
    if (!material_data.count("PlasticityViscousMultiplier"))
    {
        throw std::domain_error("\"PlasticityViscousMultiplier\" needs to be specified");
    }
    if (!material_data.count("DamageViscousExponent"))
    {
        throw std::domain_error("\"DamageViscousExponent\" needs to be specified");
    }
    if (!material_data.count("DamageViscousMultiplier"))
    {
        throw std::domain_error("\"DamageViscousMultiplier\" needs to be specified");
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
