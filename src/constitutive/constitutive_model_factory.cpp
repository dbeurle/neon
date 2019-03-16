
#include "constitutive_model_factory.hpp"

#include "mechanics/solid/small_strain_J2_plasticity.hpp"
#include "mechanics/solid/small_strain_J2_plasticity_damage.hpp"
#include "mechanics/solid/finite_strain_J2_plasticity.hpp"

#include "mechanics/solid/gaussian_affine_microsphere.hpp"
#include "mechanics/solid/gaussian_ageing_affine_microsphere.hpp"
#include "mechanics/solid/affine_microsphere.hpp"
#include "mechanics/solid/nonaffine_microsphere.hpp"
#include "mechanics/solid/compressible_neohooke.hpp"

#include "mechanics/plane/isotropic_linear_elasticity.hpp"
#include "mechanics/plane/small_strain_J2_plasticity.hpp"
#include "mechanics/plane/finite_strain_J2_plasticity.hpp"

#include "thermal/isotropic_diffusion.hpp"

#include "io/json.hpp"

#include <stdexcept>

namespace neon::mechanics::solid
{
auto make_constitutive_model(std::shared_ptr<internal_variables_t>& variables,
                             json const& material_data,
                             json const& mesh_data) -> std::unique_ptr<constitutive_model>
{
    if (mesh_data.find("constitutive") == end(mesh_data))
    {
        throw std::domain_error("Missing \"constitutive\" in \"mesh\"");
    }

    auto const& constitutive_model = mesh_data["constitutive"];

    if (constitutive_model.find("name") == end(constitutive_model))
    {
        throw std::domain_error("Missing \"name\" in \"constitutive\"");
    }

    auto const& model_name = constitutive_model["name"].get<std::string>();

    if (model_name == "neohooke" || model_name == "compressible_neohooke")
    {
        return std::make_unique<compressible_neohooke>(variables, material_data);
    }
    else if (model_name == "microsphere")
    {
        if (constitutive_model.find("type") == end(constitutive_model))
        {
            throw std::domain_error("Missing \"type\" as \"affine\" or \"nonaffine\" in "
                                    "microsphere model");
        }

        if (constitutive_model.find("quadrature") == end(constitutive_model))
        {
            throw std::domain_error("Missing \"quadrature\" as \"BO21\", \"BO33\" or \"FM900\" in "
                                    "microsphere model");
        }

        auto const& model_type = constitutive_model["type"].get<std::string>();

        std::map<std::string, unit_sphere_quadrature::point> const str_to_enum =
            {{"BO21", unit_sphere_quadrature::point::BO21},
             {"BO33", unit_sphere_quadrature::point::BO33},
             {"FM900", unit_sphere_quadrature::point::FM900}};

        auto const entry = str_to_enum.find(constitutive_model["quadrature"].get<std::string>());

        if (entry == str_to_enum.end())
        {
            throw std::domain_error("\"quadrature\" field must be \"BO21\", \"BO33\" or "
                                    "\"FM900\"");
        }
        if (model_type == "affine")
        {
            if (constitutive_model.find("statistics") == end(constitutive_model))
            {
                throw std::domain_error("Missing \"statistics\" as \"gaussian\" or \"langevin\" in "
                                        "the microsphere model");
            }

            auto const& chain_type = constitutive_model["statistics"].get<std::string>();

            if (chain_type != "gaussian" && chain_type != "langevin")
            {
                throw std::domain_error("\"statistics\" for the microsphere model must be either "
                                        "\"gaussian\" or \"langevin\"");
            }

            if (chain_type == "gaussian")
            {
                if (constitutive_model.find("ageing") != end(constitutive_model))
                {
                    if (constitutive_model["ageing"].get<std::string>() != "BAND")
                    {
                        throw std::domain_error("The only microsphere ageing model supported is "
                                                "\"BAND\"");
                    }
                    return std::make_unique<gaussian_ageing_affine_microsphere>(variables,
                                                                                material_data,
                                                                                entry->second);
                }
                return std::make_unique<gaussian_affine_microsphere>(variables,
                                                                     material_data,
                                                                     entry->second);
            }
            else if (chain_type == "langevin")
            {
                return std::make_unique<affine_microsphere>(variables, material_data, entry->second);
            }
        }
        else if (model_type == "nonaffine")
        {
            return std::make_unique<nonaffine_microsphere>(variables, material_data, entry->second);
        }
        else
        {
            throw std::domain_error("microsphere model options are \"affine\" or \"nonaffine\"");
        }
    }
    else if (model_name == "isotropic_linear_elasticity")
    {
        return std::make_unique<isotropic_linear_elasticity>(variables, material_data);
    }
    else if (model_name == "J2_plasticity")
    {
        if (constitutive_model.find("finite_strain") == end(constitutive_model))
        {
            throw std::domain_error("\"small_strain_J2_plasticity\" must have a boolean value for "
                                    "\"finite_strain\"");
        }

        if (constitutive_model.find("damage") != end(constitutive_model))
        {
            if (mesh_data["constitutive"]["finite_strain"].get<bool>())
            {
                throw std::domain_error("\"damage\" is not yet implemented for "
                                        "\"finite_strain\"");
            }

            std::string const& damage_type = constitutive_model["damage"];

            if (damage_type == "isotropic_chaboche")
            {
                return std::make_unique<small_strain_J2_plasticity_damage>(variables, material_data);
            }
        }

        if (mesh_data["constitutive"]["finite_strain"].get<bool>())
        {
            return std::make_unique<finite_strain_J2_plasticity>(variables, material_data);
        }
        return std::make_unique<small_strain_J2_plasticity>(variables, material_data);
    }
    throw std::domain_error("The model name " + model_name + " is not recognised\n"
                            + "Supported models are \"neohooke\", \"microsphere\" "
                              "\"isotropic_linear_elasticity\" and \"J2_plasticity\"\n");
    return nullptr;
}
}

namespace neon::mechanics::plane
{
auto make_constitutive_model(std::shared_ptr<internal_variables_t>& variables,
                             json const& material_data,
                             json const& mesh_data) -> std::unique_ptr<constitutive_model>
{
    if (mesh_data.find("constitutive") == end(mesh_data))
    {
        throw std::domain_error("Missing \"constitutive\" in \"mesh\"");
    }

    auto const& constitutive_model = mesh_data["constitutive"];

    if (constitutive_model.find("name") == end(constitutive_model))
    {
        throw std::domain_error("Missing \"name\" in \"constitutive\"");
    }

    std::string const& model_name = constitutive_model["name"];

    if (model_name == "plane_strain")
    {
        return std::make_unique<isotropic_linear_elasticity>(variables,
                                                             material_data,
                                                             isotropic_linear_elasticity::plane::strain);
    }
    else if (model_name == "plane_stress")
    {
        return std::make_unique<isotropic_linear_elasticity>(variables,
                                                             material_data,
                                                             isotropic_linear_elasticity::plane::stress);
    }
    else if (model_name == "J2_plasticity")
    {
        if (constitutive_model.find("finite_strain") == end(constitutive_model))
        {
            throw std::domain_error("\"small_strain_J2_plasticity\" must have a boolean value for "
                                    "\"finite_strain\"");
        }

        if (mesh_data["constitutive"]["finite_strain"].get<bool>())
        {
            return std::make_unique<finite_strain_J2_plasticity>(variables, material_data);
        }
        return std::make_unique<small_strain_J2_plasticity>(variables, material_data);
    }
    throw std::domain_error("The model name " + model_name + " is not recognised\n"
                            + "Supported models are \"plane_strain\" and \"plane_stress\"\n");
    return nullptr;
}
}

namespace neon::diffusion
{
auto make_constitutive_model(std::shared_ptr<internal_variables_t>& variables,
                             json const& material_data,
                             json const& mesh_data) -> std::unique_ptr<constitutive_model>
{
    if (mesh_data.find("constitutive") == end(mesh_data))
    {
        throw std::domain_error("A \"constitutive\" was not specified in \"mesh\"");
    }

    auto const& constitutive_model = mesh_data["constitutive"];

    if (constitutive_model.find("name") == end(constitutive_model))
    {
        throw std::domain_error("\"constitutive\" requires a \"name\" field");
    }

    if (constitutive_model["name"].get<std::string>() == "isotropic_diffusion")
    {
        return std::make_unique<isotropic_diffusion>(variables, material_data);
    }

    throw std::domain_error("The model name " + constitutive_model["name"].get<std::string>()
                            + " is not recognised\n"
                            + "The supported model is \"isotropic_diffusion\"\n");

    return nullptr;
}
}
