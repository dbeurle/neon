
#include "PerfectPlasticElastic.hpp"

#include <json/json.h>

#include "MaterialExceptions.hpp"

namespace neon
{
PerfectPlasticElastic::PerfectPlasticElastic(const Json::Value& material_data)
    : PlasticMaterial(material_data)
{
    if (material_data["YieldStress"].empty()) throw MaterialPropertyException("YieldStress");

    yield_stress = material_data["YieldStress"].asDouble();
}

double PerfectPlasticElastic::evaluate_yield_function(double effective_strain) const
{
    return yield_stress;
}

double PerfectPlasticElastic::plastic_modulus(double effective_strain) const { return 0.0; }
}
