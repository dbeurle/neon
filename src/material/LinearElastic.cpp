
#include "LinearElastic.hpp"

#include <json/json.h>

#include "MaterialExceptions.hpp"

namespace neon
{
LinearElastic::LinearElastic(Json::Value const& material_data)
    : Material(material_data["Name"].asString())
{
    if (material_data["ElasticModulus"].empty()) throw MaterialPropertyException("Elastic Modulus");
    if (material_data["PoissonsRatio"].empty()) throw MaterialPropertyException("Poisson's Ratio");

    // Read in from Json file
    E = material_data["ElasticModulus"].asDouble();
    nu = material_data["PoissonsRatio"].asDouble();
}
}
