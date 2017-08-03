
#include "ConstitutiveModelFactory.hpp"

#include "AffineMicrosphere.hpp"
#include "HyperElasticPlastic.hpp"
#include "J2Plasticity.hpp"
#include "NeoHooke.hpp"

#include "InternalVariables.hpp"
#include "PreprocessorExceptions.hpp"

#include <json/value.h>

namespace neon
{
namespace solid
{
std::unique_ptr<ConstitutiveModel> make_constitutive_model(InternalVariables& variables,
                                                           Json::Value const& material_data,
                                                           Json::Value const& simulation_data)
{
    if (simulation_data["ConstitutiveModel"].empty())
    {
        throw EmptyFieldException("Part: ConstitutiveModel");
    }

    auto const& model_name = simulation_data["ConstitutiveModel"].asString();

    if (model_name == "NeoHooke")
    {
        return std::make_unique<NeoHooke>(variables, material_data);
    }
    else if (model_name == "AffineMicrosphere")
    {
        return std::make_unique<AffineMicrosphere>(variables, material_data);
    }
    else if (model_name == "J2")
    {
        return std::make_unique<J2Plasticity>(variables, material_data);
    }

    throw std::runtime_error("The model name " + model_name + " is not recognised\n"
                             + "Supported models are \"NeoHooke\", \"AffineMicrosphere\" "
                               "and \"J2\"\n");

    return nullptr;
}
}
}
