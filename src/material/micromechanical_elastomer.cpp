
#include "micromechanical_elastomer.hpp"

#include "numeric/dense_matrix.hpp"
#include "numeric/float_compare.hpp"

#include "io/json.hpp"

#include <cmath>

using namespace neon;

micromechanical_elastomer::micromechanical_elastomer(json const& material_data)
    : isotropic_elastic_property(material_data)
{
    if (!material_data.count("SegmentsPerChain"))
    {
        throw std::domain_error("SegmentsPerChain not specified in material data\n");
    }
    N = material_data["SegmentsPerChain"];
}

ageing_micromechanical_elastomer::ageing_micromechanical_elastomer(json const& material_data)
    : micromechanical_elastomer(material_data)
{
    auto exception_string = [](std::string&& field) {
        return "\"" + field + "\" is not specified in \"Material\" data";
    };

    if (!material_data.count("ScissionProbability"))
    {
        throw std::domain_error(exception_string("ScissionProbability"));
    }

    if (!material_data.count("RecombinationProbability"))
    {
        throw std::domain_error(exception_string("RecombinationProbability"));
    }

    pr_scission = material_data["ScissionProbability"];
    pr_recombination = material_data["RecombinationProbability"];

    if (pr_scission < 0.0 || pr_recombination < 0.0)
    {
        throw std::domain_error("Material properties (probabilities) must be positive");
    }
}
