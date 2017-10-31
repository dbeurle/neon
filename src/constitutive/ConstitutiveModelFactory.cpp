
#include "ConstitutiveModelFactory.hpp"

#include "J2Plasticity.hpp"

#include "FiniteJ2Plasticity.hpp"

#include "AffineMicrosphere.hpp"
#include "NonAffineMicrosphere.hpp"

#include "NeoHooke.hpp"

#include "IsotropicDiffusion.hpp"

#include <json/value.h>
#include <stdexcept>

namespace neon::solid
{
std::unique_ptr<ConstitutiveModel> make_constitutive_model(InternalVariables& variables,
                                                           Json::Value const& material_data,
                                                           Json::Value const& simulation_data)
{
    if (!simulation_data.isMember("ConstitutiveModel"))
    {
        throw std::runtime_error("Missing \"ConstitutiveModel\" in \"Mesh\"");
    }
    if (!simulation_data["ConstitutiveModel"].isMember("Name"))
    {
        throw std::runtime_error("Missing \"Name\" in \"ConstitutiveModel\"");
    }

    auto const& model_name = simulation_data["ConstitutiveModel"]["Name"].asString();

    if (model_name == "NeoHooke")
    {
        return std::make_unique<NeoHooke>(variables, material_data);
    }
    else if (model_name == "Microsphere")
    {
        auto const& model_type = simulation_data["ConstitutiveModel"]["Type"].asString();

        // TODO Check for quadrature schemes and ageing model

        if (model_type == "Affine")
        {
            return std::make_unique<AffineMicrosphere>(variables, material_data);
        }
        else if (model_type == "NonAffine")
        {
            return std::make_unique<NonAffineMicrosphere>(variables, material_data);
        }
        else
        {
            throw std::runtime_error("Microsphere model options are \"Affine\" or \"Nonaffine\"");
        }
    }
    else if (model_name == "IsotropicLinearElasticity")
    {
        return std::make_unique<IsotropicLinearElasticity>(variables, material_data);
    }
    else if (model_name == "J2Plasticity")
    {
        if (!simulation_data["ConstitutiveModel"].isMember("FiniteStrain"))
        {
            throw std::runtime_error("\"J2Plasticity\" must have a boolean value for "
                                     "\"FiniteStrain\"");
        }

        if (simulation_data["ConstitutiveModel"]["FiniteStrain"].asBool())
        {
            return std::make_unique<FiniteJ2Plasticity>(variables, material_data);
        }
        return std::make_unique<J2Plasticity>(variables, material_data);
    }
    throw std::runtime_error("The model name " + model_name + " is not recognised\n"
                             + "Supported models are \"NeoHooke\", \"AffineMicrosphere\" "
                               "and \"J2\"\n");
    return nullptr;
}
}

namespace neon::diffusion
{
std::unique_ptr<ConstitutiveModel> make_constitutive_model(InternalVariables& variables,
                                                           Json::Value const& material_data,
                                                           Json::Value const& simulation_data)
{
    if (!simulation_data.isMember("ConstitutiveModel"))
    {
        throw std::runtime_error("Missing \"ConstitutiveModel\" in \"Mesh\"");
    }

    auto const& model_name = simulation_data["ConstitutiveModel"].asString();

    if (model_name == "IsotropicDiffusion")
    {
        return std::make_unique<IsotropicDiffusion>(variables, material_data);
    }

    throw std::runtime_error("The model name " + model_name + " is not recognised\n"
                             + "The supported model is \"IsotropicDiffusion\"\n");

    return nullptr;
}
}
