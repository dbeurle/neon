
#include "constitutive_model_factory.hpp"

#include "mechanical/solid/small_strain_J2_plasticity.hpp"
#include "mechanical/solid/small_strain_J2_plasticity_damage.hpp"
#include "mechanical/solid/finite_strain_J2_plasticity.hpp"

#include "mechanical/solid/gaussian_affine_microsphere.hpp"
#include "mechanical/solid/gaussian_ageing_affine_microsphere.hpp"
#include "mechanical/solid/affine_microsphere.hpp"
#include "mechanical/solid/nonaffine_microsphere.hpp"
#include "mechanical/solid/compressible_neohooke.hpp"

#include "mechanical/plane/isotropic_linear_elasticity.hpp"
#include "mechanical/plane/small_strain_J2_plasticity.hpp"
#include "mechanical/plane/finite_strain_J2_plasticity.hpp"

#include "thermal/isotropic_diffusion.hpp"

#include "io/json.hpp"

#include <stdexcept>

namespace neon::mechanical::solid
{
std::unique_ptr<constitutive_model> make_constitutive_model(
    std::shared_ptr<internal_variables_t>& variables, json const& material_data, json const& mesh_data)
{
    if (!mesh_data.count("ConstitutiveModel"))
    {
        throw std::domain_error("Missing \"ConstitutiveModel\" in \"Mesh\"");
    }

    auto const& constitutive_model = mesh_data["ConstitutiveModel"];

    if (!constitutive_model.count("Name"))
    {
        throw std::domain_error("Missing \"Name\" in \"ConstitutiveModel\"");
    }

    auto const& model_name = constitutive_model["Name"].get<std::string>();

    if (model_name == "NeoHooke" || model_name == "CompressibleNeoHooke")
    {
        return std::make_unique<compressible_neohooke>(variables, material_data);
    }
    else if (model_name == "Microsphere")
    {
        if (!constitutive_model.count("Type"))
        {
            throw std::domain_error("Missing \"Type\" as \"Affine\" or \"NonAffine\" in "
                                    "Microsphere model");
        }

        if (!constitutive_model.count("Quadrature"))
        {
            throw std::domain_error("Missing \"Quadrature\" as \"BO21\", \"BO33\" or \"FM900\" in "
                                    "Microsphere model");
        }

        auto const& model_type = constitutive_model["Type"].get<std::string>();

        std::map<std::string, unit_sphere_quadrature::Rule> const str_to_enum =
            {{"BO21", unit_sphere_quadrature::Rule::BO21},
             {"BO33", unit_sphere_quadrature::Rule::BO33},
             {"FM900", unit_sphere_quadrature::Rule::FM900}};

        auto const entry = str_to_enum.find(constitutive_model["Quadrature"].get<std::string>());

        if (entry == str_to_enum.end())
        {
            throw std::domain_error("\"Quadrature\" field must be \"BO21\", \"BO33\" or "
                                    "\"FM900\"");
        }
        if (model_type == "Affine")
        {
            if (!constitutive_model.count("Statistics"))
            {
                throw std::domain_error("Missing \"Statistics\" as \"Gaussian\" or \"Langevin\" in "
                                        "the "
                                        "Microsphere model");
            }

            auto const& chain_type = constitutive_model["Statistics"].get<std::string>();

            if (chain_type != "Gaussian" && chain_type != "Langevin")
            {
                throw std::domain_error("\"Statistics\" for the Microsphere model must be either "
                                        "\"Gaussian\" or \"Langevin\"");
            }

            if (chain_type == "Gaussian")
            {
                if (constitutive_model.count("Ageing"))
                {
                    return std::make_unique<gaussian_ageing_affine_microsphere>(variables,
                                                                                material_data,
                                                                                entry->second);
                }
                return std::make_unique<gaussian_affine_microsphere>(variables,
                                                                     material_data,
                                                                     entry->second);
            }
            else if (chain_type == "Langevin")
            {
                return std::make_unique<affine_microsphere>(variables, material_data, entry->second);
            }
        }
        else if (model_type == "NonAffine")
        {
            return std::make_unique<nonaffine_microsphere>(variables, material_data, entry->second);
        }
        else
        {
            throw std::domain_error("Microsphere model options are \"Affine\" or \"Nonaffine\"");
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
            throw std::domain_error("\"small_strain_J2_plasticity\" must have a boolean value for "
                                    "\"FiniteStrain\"");
        }

        if (constitutive_model.count("Damage"))
        {
            if (mesh_data["ConstitutiveModel"]["FiniteStrain"].get<bool>())
            {
                throw std::domain_error("\"J2PlasticityDamage\" is not implemented for "
                                        "\"FiniteStrain\"");
            }

            auto const& damage_type = constitutive_model["Damage"].get<std::string>();

            if (damage_type == "IsotropicChaboche")
                return std::make_unique<small_strain_J2_plasticity_damage>(variables, material_data);
        }

        if (mesh_data["ConstitutiveModel"]["FiniteStrain"].get<bool>())
        {
            return std::make_unique<finite_strain_J2_plasticity>(variables, material_data);
        }
        return std::make_unique<small_strain_J2_plasticity>(variables, material_data);
    }
    throw std::domain_error("The model name " + model_name + " is not recognised\n"
                            + "Supported models are \"NeoHooke\", \"Microsphere\" "
                              "and \"small_strain_J2_plasticity\"\n");
    return nullptr;
}
}

namespace neon::mechanical::plane
{
std::unique_ptr<constitutive_model> make_constitutive_model(
    std::shared_ptr<internal_variables_t>& variables, json const& material_data, json const& mesh_data)
{
    if (!mesh_data.count("ConstitutiveModel"))
    {
        throw std::domain_error("Missing \"ConstitutiveModel\" in \"Mesh\"");
    }

    auto const& constitutive_model = mesh_data["ConstitutiveModel"];

    if (!constitutive_model.count("Name"))
    {
        throw std::domain_error("Missing \"Name\" in \"ConstitutiveModel\"");
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
            throw std::domain_error("\"small_strain_J2_plasticity\" must have a boolean value for "
                                    "\"FiniteStrain\"");
        }

        if (mesh_data["ConstitutiveModel"]["FiniteStrain"].get<bool>())
        {
            return std::make_unique<finite_strain_J2_plasticity>(variables, material_data);
        }
        return std::make_unique<small_strain_J2_plasticity>(variables, material_data);
    }
    throw std::domain_error("The model name " + model_name + " is not recognised\n"
                            + "Supported models are \"PlaneStrain\" and \"PlaneStress\"\n");
    return nullptr;
}
}

namespace neon::diffusion
{
std::unique_ptr<constitutive_model> make_constitutive_model(
    std::shared_ptr<internal_variables_t>& variables, json const& material_data, json const& mesh_data)
{
    if (!mesh_data.count("ConstitutiveModel"))
    {
        throw std::domain_error("A \"ConstitutiveModel\" was not specified in \"Mesh\"");
    }

    auto const& constitutive_model = mesh_data["ConstitutiveModel"];

    if (!constitutive_model.count("Name"))
    {
        throw std::domain_error("\"ConstitutiveModel\" requires a \"Name\" field");
    }

    if (constitutive_model["Name"].get<std::string>() == "IsotropicDiffusion")
    {
        return std::make_unique<isotropic_diffusion>(variables, material_data);
    }

    throw std::domain_error("The model name " + constitutive_model["Name"].get<std::string>()
                            + " is not recognised\n"
                            + "The supported model is \"IsotropicDiffusion\"\n");

    return nullptr;
}
}
