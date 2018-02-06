
#include "constitutive_model_factory.hpp"

#include "mechanical/solid/small_strain_J2_plasticity.hpp"
#include "mechanical/solid/small_strain_J2_plasticity_damage.hpp"

#include "mechanical/solid/finite_strain_J2_plasticity.hpp"

#include "mechanical/solid/gaussian_affine_microsphere.hpp"

#include "mechanical/solid/affine_microsphere.hpp"
#include "mechanical/solid/nonaffine_microsphere.hpp"

#include "mechanical/solid/compressible_neohooke.hpp"

#include "thermal/IsotropicDiffusion.hpp"

#include "mechanical/plane/isotropic_linear_elasticity.hpp"

#include "mechanical/plane/small_strain_J2_plasticity.hpp"

#include "mechanical/plane/finite_strain_J2_plasticity.hpp"

#include "io/json.hpp"

#include <stdexcept>

namespace neon::mechanical::solid
{
std::unique_ptr<ConstitutiveModel> make_constitutive_model(
    std::shared_ptr<InternalVariables>& variables, json const& material_data, json const& mesh_data)
{
    if (!mesh_data.count("ConstitutiveModel"))
    {
        throw std::runtime_error("Missing \"ConstitutiveModel\" in \"Mesh\"");
    }

    auto const& constitutive_model = mesh_data["ConstitutiveModel"];

    if (!constitutive_model.count("Name"))
    {
        throw std::runtime_error("Missing \"Name\" in \"ConstitutiveModel\"");
    }

    auto const& model_name = constitutive_model["Name"].get<std::string>();

    if (model_name == "NeoHooke")
    {
        return std::make_unique<compressible_neohooke>(variables, material_data);
    }
    else if (model_name == "CompressibleNeoHooke")
    {
        return std::make_unique<compressible_neohooke>(variables, material_data);
    }
    else if (model_name == "Microsphere")
    {
        if (!constitutive_model.count("Type"))
        {
            throw std::runtime_error("Missing \"Type\" as \"Affine\" or \"NonAffine\" in "
                                     "Microsphere model");
        }

        if (!constitutive_model.count("Quadrature"))
        {
            throw std::runtime_error("Missing \"Quadrature\" as \"BO21\", \"BO33\" or \"FM900\" in "
                                     "Microsphere model");
        }

        auto const& model_type = constitutive_model["Type"].get<std::string>();

        // TODO Check for ageing model

        std::map<std::string, unit_sphere_quadrature::Rule> const str_to_enum =
            {{"BO21", unit_sphere_quadrature::Rule::BO21},
             {"BO33", unit_sphere_quadrature::Rule::BO33},
             {"FM900", unit_sphere_quadrature::Rule::FM900}};

        auto const entry = str_to_enum.find(constitutive_model["Quadrature"].get<std::string>());

        if (entry == str_to_enum.end())
        {
            throw std::runtime_error("\"Quadrature\" field must be \"BO21\", \"BO33\" or "
                                     "\"FM900\"");
        }

        if (model_type == "GaussianAffine")
        {
            return std::make_unique<gaussian_affine_microsphere>(variables,
                                                                 material_data,
                                                                 entry->second);
        }
        else if (model_type == "Affine")
        {
            return std::make_unique<affine_microsphere>(variables, material_data, entry->second);
        }
        else if (model_type == "NonAffine")
        {
            return std::make_unique<nonaffine_microsphere>(variables, material_data, entry->second);
        }
        else
        {
            throw std::runtime_error("Microsphere model options are \"Affine\" or \"Nonaffine\"");
        }
    }
    else if (model_name == "IsotropicLinearElasticity")
    {
        return std::make_unique<isotropic_linear_elasticity>(variables, material_data);
    }
    else if (model_name == "J2Plasticity")
    {
        if (!constitutive_model.count("FiniteStrain"))
        {
            throw std::runtime_error("\"small_strain_J2_plasticity\" must have a boolean value for "
                                     "\"FiniteStrain\"");
        }

        if (mesh_data["ConstitutiveModel"]["FiniteStrain"].get<bool>())
        {
            return std::make_unique<finite_strain_J2_plasticity>(variables, material_data);
        }
        return std::make_unique<small_strain_J2_plasticity>(variables, material_data);
    }
    else if (model_name == "ChabocheDamage")
    {
        if (mesh_data["ConstitutiveModel"]["FiniteStrain"].get<bool>())
        {
            throw std::runtime_error("\"ChabocheDamage\" is not implemented for "
                                     "\"FiniteStrain\"");
        }
        return std::make_unique<small_strain_J2_plasticity_damage>(variables, material_data);
    }
    throw std::runtime_error("The model name " + model_name + " is not recognised\n"
                             + "Supported models are \"NeoHooke\", \"Microsphere\" "
                               "and \"small_strain_J2_plasticity\"\n");
    return nullptr;
}
}

namespace neon::mechanical::plane
{
std::unique_ptr<ConstitutiveModel> make_constitutive_model(
    std::shared_ptr<InternalVariables>& variables, json const& material_data, json const& mesh_data)
{
    if (!mesh_data.count("ConstitutiveModel"))
    {
        throw std::runtime_error("Missing \"ConstitutiveModel\" in \"Mesh\"");
    }

    auto const& constitutive_model = mesh_data["ConstitutiveModel"];

    if (!constitutive_model.count("Name"))
    {
        throw std::runtime_error("Missing \"Name\" in \"ConstitutiveModel\"");
    }

    auto const& model_name = constitutive_model["Name"].get<std::string>();

    if (model_name == "PlaneStrain")
    {
        return std::make_unique<isotropic_linear_elasticity>(variables,
                                                             material_data,
                                                             isotropic_linear_elasticity::plane::strain);
    }
    else if (model_name == "PlaneStress")
    {
        return std::make_unique<isotropic_linear_elasticity>(variables,
                                                             material_data,
                                                             isotropic_linear_elasticity::plane::stress);
    }
    else if (model_name == "J2Plasticity")
    {
        if (!constitutive_model.count("FiniteStrain"))
        {
            throw std::runtime_error("\"small_strain_J2_plasticity\" must have a boolean value for "
                                     "\"FiniteStrain\"");
        }

        if (mesh_data["ConstitutiveModel"]["FiniteStrain"].get<bool>())
        {
            return std::make_unique<finite_strain_J2_plasticity>(variables, material_data);
        }
        return std::make_unique<small_strain_J2_plasticity>(variables, material_data);
    }
    throw std::runtime_error("The model name " + model_name + " is not recognised\n"
                             + "Supported models are \"PlaneStrain\" and \"PlaneStress\"\n");
    return nullptr;
}
}

namespace neon::diffusion
{
std::unique_ptr<ConstitutiveModel> make_constitutive_model(
    std::shared_ptr<InternalVariables>& variables, json const& material_data, json const& mesh_data)
{
    if (!mesh_data.count("ConstitutiveModel"))
    {
        throw std::runtime_error("A \"ConstitutiveModel\" was not specified in \"Mesh\"");
    }

    auto const& constitutive_model = mesh_data["ConstitutiveModel"];

    if (!constitutive_model.count("Name"))
    {
        throw std::runtime_error("\"ConstitutiveModel\" requires a \"Name\" field");
    }

    if (constitutive_model["Name"].get<std::string>() == "IsotropicDiffusion")
    {
        return std::make_unique<IsotropicDiffusion>(variables, material_data);
    }

    throw std::runtime_error("The model name " + constitutive_model["Name"].get<std::string>()
                             + " is not recognised\n"
                             + "The supported model is \"IsotropicDiffusion\"\n");

    return nullptr;
}
}
