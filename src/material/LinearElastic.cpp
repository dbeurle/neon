
#include "LinearElastic.hpp"

#include "io/json.hpp"

#include "Exceptions.hpp"

namespace neon
{
LinearElastic::LinearElastic(json const& material_data) : Material(material_data)
{
    // Determine input value types.  The allowable inputs are:
    //   1)  Elastic modulus and Poisson ratio
    //   2)  Bulk modulus and shear modulus

    if (material_data.count("BulkModulus") && material_data.count("ShearModulus"))
    {
        // Fill the elastic modulus and Poisson's ratio from the bulk and shear
        // moduli
        K = material_data["BulkModulus"];
        G = material_data["ShearModulus"];
    }
    else if (material_data.count("ElasticModulus") && material_data.count("PoissonsRatio"))
    {
        auto const E = material_data["ElasticModulus"].get<double>();
        auto const nu = material_data["PoissonsRatio"].get<double>();

        if (nu > 0.5)
        {
            throw MaterialPropertyException("\"PoissonsRatio\" must be less than or equal to 0.5");
        }

        K = E / (3.0 * (1.0 - 2.0 * nu));
        G = E / (2.0 * (1.0 + nu));
    }
    else
    {
        throw MaterialPropertyException("\"ElasticModulus\" and \"PoissonsRatio\" or "
                                        "\"BulkModulus\" and \"ShearModulus\" need to be "
                                        "specified as material properties");
    }
}
}
