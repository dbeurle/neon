
#include "isotropic_elastic_plastic.hpp"

#include "exceptions.hpp"

#include "io/json.hpp"

namespace neon
{
isotropic_elastic_plastic::isotropic_elastic_plastic(json const& material_data)
    : isotropic_elastic_property(material_data)
{
    if (material_data.find("yield_stress") == end(material_data))
    {
        throw std::domain_error("\"yield_stress\" was not specified in material properties");
    }

    stress_y = material_data["yield_stress"];

    // Optional material parameters
    if (material_data.find("isotropic_hardening_modulus") != end(material_data))
    {
        H_iso = material_data["isotropic_hardening_modulus"];
    }
    if (material_data.find("isotropic_kinematic_modulus") != end(material_data))
    {
        K_iso = material_data["isotropic_kinematic_modulus"];
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
