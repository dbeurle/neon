
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
    auto const& mesh_data = simulation["Mesh"].front();

    auto const& simulation_mesh = mesh_store.find(mesh_data["Name"].get<std::string>());

    auto const& name = simulation["Name"].get<std::string>();
    auto const& module_type = simulation["Module"].get<std::string>();
    auto const& solution_type = simulation["Solution"].get<std::string>();

    std::cout << std::string(4, ' ') << "Name     \"" << name << "\"\n";
    std::cout << std::string(4, ' ') << "Module   \"" << module_type << "\"\n";
    std::cout << std::string(4, ' ') << "Solution \"" << solution_type << "\"\n";

    auto const& [mesh, material] = simulation_mesh->second;

    if (module_type == "SolidMechanics")
    {
        if (simulation.find("NonlinearOptions") == simulation.end())
        {
            throw std::domain_error("\"NonlinearOptions\" needs to be present for a "
                                    "SolidMechanics simulation");
        }
        return std::make_unique<solid_mechanics_module>(mesh, material, simulation);
    }
    else if (module_type == "PlaneStrain")
    {
        if (simulation.find("NonlinearOptions") == simulation.end())
        {
            throw std::domain_error("\"NonlinearOptions\" needs to be present for a "
                                    "PlaneStrain simulation");
        }
        return std::make_unique<plane_strain_module>(mesh, material, simulation);
    }
    else if (module_type == "Beam")
    {
        return std::make_unique<beam_module>(mesh, material, simulation, profile_store);
    }
    else if (module_type == "HeatDiffusion")
    {
        if (solution_type == "Equilibrium")
        {
            return std::make_unique<linear_diffusion_module<diffusion::static_matrix>>(mesh,
                                                                                       material,
                                                                                       simulation);
        }
        else if (solution_type == "Transient")
        {
            return std::make_unique<linear_diffusion_module<diffusion::dynamic_matrix>>(mesh,
                                                                                        material,
                                                                                        simulation);
        }
        throw std::domain_error("Solution " + solution_type
                                + " is not recognised.  Use \"Equilibrium\" or \"Transient\"\n");
    }
    throw std::domain_error("Module " + module_type
                            + " is not recognised.  Use \"PlaneStrain\", \"SolidMechanics\" or "
                              "\"HeatDiffusion\"\n");
    return nullptr;
}
}
