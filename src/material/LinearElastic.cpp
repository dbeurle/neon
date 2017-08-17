
#include "LinearElastic.hpp"

#include <json/json.h>

#include "Exceptions.hpp"

namespace neon
{
LinearElastic::LinearElastic(Json::Value const& material_data) : Material(material_data)
{
    // Determine input value types.  The allowable inputs are:
    //   1)  Elastic modulus and Poisson ratio
    //   2)  Bulk modulus and shear modulus

    if (material_data.isMember("BulkModulus") && material_data.isMember("ShearModulus"))
    {
        // Fill the elastic modulus and Poisson's ratio from the bulk and shear
        // moduli
        auto const K = material_data["BulkModulus"].asDouble();
        auto const G = material_data["ShearModulus"].asDouble();

        E = 9 * K * G / (3 * K + G);
        nu = (3 * K - 2 * G) / (3 * K + G);
    }
    else if (material_data.isMember("ElasticModulus") && material_data.isMember("PoissonsRatio"))
    {
        E = material_data["ElasticModulus"].asDouble();
        nu = material_data["PoissonsRatio"].asDouble();
    }
    else
    {
        throw MaterialPropertyException("\"ElasticModulus\" and \"PoissonsRatio\" or "
                                        "\"BulkModulus\" and \"ShearModulus\" need to be "
                                        "specified as material properties");
    }
}
}
