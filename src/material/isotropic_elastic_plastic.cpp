
#include "isotropic_elastic_plastic.hpp"

#include "exceptions.hpp"

#include "io/json.hpp"

namespace neon
{
isotropic_elastic_plastic::isotropic_elastic_plastic(json const& material_data)
    : isotropic_elastic_property(material_data)
{
    if (material_data.find("YieldStress") == material_data.end())
    {
        throw std::domain_error("\"YieldStress\" was not specified in material properties");
    }

    stress_y = material_data["YieldStress"];

    // Optional material parameters
    if (material_data.find("IsotropicHardeningModulus") != material_data.end())
    {
        H_iso = material_data["IsotropicHardeningModulus"];
    }
    if (material_data.find("IsotropicKinematicModulus") != material_data.end())
    {
        K_iso = material_data["IsotropicKinematicModulus"];
    }
}

double isotropic_elastic_plastic::yield_stress(double const effective_strain) const
{
    return stress_y + effective_strain * H_iso;
}

double isotropic_elastic_plastic::hardening_modulus(double const effective_strain) const
{
    return H_iso;
}

double isotropic_elastic_plastic::kinematic_modulus(double const effective_strain) const
{
    return K_iso;
}
}
