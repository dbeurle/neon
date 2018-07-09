
#include "isotropic_elastic_property.hpp"

#include "io/json.hpp"

#include <stdexcept>

namespace neon
{
isotropic_elastic_property::isotropic_elastic_property(json const& material_data)
    : material_property(material_data)
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
        double const E = material_data["ElasticModulus"];
        double const nu = material_data["PoissonsRatio"];

        if (nu >= 0.5)
        {
            throw std::domain_error("\"PoissonsRatio\" must be less than or equal to 0.5");
        }

        K = E / (3.0 * (1.0 - 2.0 * nu));
        G = E / (2.0 * (1.0 + nu));
    }
    else
    {
        throw std::domain_error("\"ElasticModulus\" and \"PoissonsRatio\" or "
                                "\"BulkModulus\" and \"ShearModulus\" need to be "
                                "specified as material properties");
    }
}
}
