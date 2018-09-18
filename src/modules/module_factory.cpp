
#include "module_factory.hpp"

#include "modules/solid_mechanics_module.hpp"
#include "modules/plane_strain_module.hpp"
#include "modules/beam_module.hpp"

#include "modules/linear_diffusion_module.hpp"

#include "geometry/profile.hpp"

#include "io/json.hpp"

namespace neon
{
std::unique_ptr<abstract_module> make_module(
    json const& simulation,
    std::map<std::string, std::pair<basic_mesh, json>> const& mesh_store,
    std::map<std::string, std::unique_ptr<geometry::profile>> const& profile_store)
{
    if (simulation.find("name") == simulation.end())
    {
        throw std::domain_error("\"name\" field must be specified");
    }
    if (simulation.find("module") == simulation.end())
    {
        throw std::domain_error("\"module\" field must be specified");
    }
    if (simulation.find("solution") == simulation.end())
    {
        throw std::domain_error("\"solution\" field must be specified");
    }

    auto const& mesh_data = simulation["meshes"].front();

    auto const& simulation_mesh = mesh_store.find(mesh_data["name"].get<std::string>());

    auto const& name = simulation["name"].get<std::string>();
    auto const& module_type = simulation["module"].get<std::string>();
    auto const& solution_type = simulation["solution"].get<std::string>();

    std::cout << std::string(4, ' ') << "Name     \"" << name << "\"\n";
    std::cout << std::string(4, ' ') << "Module   \"" << module_type << "\"\n";
    std::cout << std::string(4, ' ') << "Solution \"" << solution_type << "\"\n";

    auto const& [mesh, material] = simulation_mesh->second;

    if (module_type == "solid_mechanics")
    {
        if (solution_type == "linear_buckling")
        {
            return std::make_unique<natural_frequency_module>(mesh, material, simulation);
        }
        else if (solution_type == "equilibrium")
        {
            if (simulation.find("nonlinear_options") == simulation.end())
            {
                throw std::domain_error("\"nonlinear_options\" needs to be present for a "
                                        "solid_mechanics simulation");
            }
            return std::make_unique<
                solid_mechanics_module<mechanics::static_matrix<mechanics::solid::mesh>>>(mesh,
                                                                                          material,
                                                                                          simulation);
        }
        else if (solution_type == "latin")
        {
            return std::make_unique<
                solid_mechanics_module<mechanics::latin_matrix<mechanics::solid::mesh>>>(mesh,
                                                                                         material,
                                                                                         simulation);
        }
        {
            throw std::domain_error("\"solution\" is not valid.  Please use \"equilibrium\"");
        }
    }
    else if (module_type == "plane_strain")
    {
        if (simulation.find("nonlinear_options") == simulation.end())
        {
            throw std::domain_error("\"nonlinear_options\" needs to be present for a "
                                    "plane_strain simulation");
        }
        return std::make_unique<plane_strain_module>(mesh, material, simulation);
    }
    else if (module_type == "beam")
    {
        return std::make_unique<beam_module>(mesh, material, simulation, profile_store);
    }
    else if (module_type == "heat_diffusion")
    {
        if (solution_type == "equilibrium")
        {
            return std::make_unique<linear_diffusion_module<diffusion::static_matrix>>(mesh,
                                                                                       material,
                                                                                       simulation);
        }
        else if (solution_type == "transient")
        {
            return std::make_unique<linear_diffusion_module<diffusion::dynamic_matrix>>(mesh,
                                                                                        material,
                                                                                        simulation);
        }
        throw std::domain_error("\"solution\" " + solution_type
                                + " is not recognised.  Use \"equilibrium\" or \"transient\"\n");
    }
    throw std::domain_error("\"module\" : \"" + module_type
                            + "\" is not recognised.  Use \"plane_strain\", \"solid_mechanics\" or "
                              "\"heat_diffusion\"\n");
    return nullptr;
}
}
