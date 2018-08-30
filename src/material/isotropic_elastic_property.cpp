
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

    if (material_data.count("bulk_modulus") && material_data.count("shear_modulus"))
    {
        // Fill the elastic modulus and Poisson's ratio from the bulk and shear
        // moduli
        K = material_data["bulk_modulus"];
        G = material_data["shear_modulus"];
    }
    else if (material_data.count("elastic_modulus") && material_data.count("poissons_ratio"))
    {
        double const E = material_data["elastic_modulus"];
        double const nu = material_data["poissons_ratio"];

        if (nu >= 0.5)
        {
            throw std::domain_error("\"poissons_ratio\" must be less than or equal to 0.5");
        }

        K = E / (3.0 * (1.0 - 2.0 * nu));
        G = E / (2.0 * (1.0 + nu));
    }
    else
    {
        throw std::domain_error("\"elastic_modulus\" and \"poissons_ratio\" or "
                                "\"bulk_modulus\" and \"shear_modulus\" need to be "
                                "specified as material properties");
    }
}
}
