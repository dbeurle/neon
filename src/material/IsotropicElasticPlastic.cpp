
#include "IsotropicElasticPlastic.hpp"

#include <json/json.h>

#include "MaterialExceptions.hpp"

namespace neon
{
IsotropicPlasticElastic::IsotropicPlasticElastic(Json::Value const& material_data)
    : PerfectPlasticElastic(material_data)
{
    if (material_data["UltimateStress"].empty()) throw MaterialPropertyException("UltimateStress");
    if (material_data["UltimateStrain"].empty()) throw MaterialPropertyException("UltimateStrain");

    failure_stress = material_data["UltimateStress"].asDouble();
    failure_strain = material_data["UltimateStrain"].asDouble();
}

double IsotropicPlasticElastic::evaluate_yield_function(double effective_strain) const
{
    return yield_stress;
}

double IsotropicPlasticElastic::plastic_modulus(double effective_strain) const
{
    return effective_strain < failure_strain ? 0.0 : failure_stress / failure_strain;
}
}
