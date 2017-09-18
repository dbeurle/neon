
#include "ModuleFactory.hpp"

#include "LinearDiffusionModule.hpp"
#include "SolidMechanicsModule.hpp"

#include <json/value.h>

namespace neon
{
std::unique_ptr<AbstractModule> make_module(
    Json::Value const& simulation,
    std::map<std::string, std::pair<BasicMesh, Json::Value>> const& mesh_store)
{
    auto const& mesh_data = simulation["Mesh"][0];

    auto simulation_mesh = mesh_store.find(mesh_data["Name"].asString());

    auto const& name = simulation["Name"].asString();
    auto const& module_type = simulation["Module"].asString();
    auto const& solution_type = simulation["Solution"].asString();

    std::cout << std::string(4, ' ') << "Name \"" << name << "\"\n";
    std::cout << std::string(4, ' ') << "Module \"" << module_type << "\"\n";
    std::cout << std::string(4, ' ') << "Solution \"" << solution_type << "\"\n";

    auto const & [ mesh, material ] = simulation_mesh->second;

    if (module_type == "SolidMechanics")
    {
        if (!simulation.isMember("NonlinearOptions"))
        {
            throw std::runtime_error("\"NonlinearOptions\" needs to be present for a "
                                     "SolidMechanics simulation");
        }
        return std::make_unique<SolidMechanicsModule>(mesh, material, simulation);
    }
    else if (module_type == "HeatDiffusion")
    {
        if (solution_type == "Equilibrium")
        {
            return std::make_unique<LinearDiffusionModule<diffusion::femStaticMatrix>>(mesh,
                                                                                       material,
                                                                                       simulation);
        }
        else if (solution_type == "Transient")
        {
            return std::make_unique<LinearDiffusionModule<diffusion::femDynamicMatrix>>(mesh,
                                                                                        material,
                                                                                        simulation);
        }
        throw std::runtime_error("Solution " + solution_type
                                 + " is not recognised.  Use \"Equilibrium\" or \"Transient\"\n");
    }
    throw std::runtime_error("Module " + module_type
                             + " is not recognised.  Use \"SolidMechanics\" or "
                               "\"HeatDiffusion\"\n");
    return nullptr;
}
}
