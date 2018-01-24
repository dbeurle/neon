
#include "ModuleFactory.hpp"

#include "LinearDiffusionModule.hpp"
#include "SolidMechanicsModule.hpp"

#include "io/json.hpp"

namespace neon
{
std::unique_ptr<AbstractModule> make_module(
    json const& simulation, std::map<std::string, std::pair<BasicMesh, json>> const& mesh_store)
{
    auto const& mesh_data = simulation["Mesh"][0];

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
        if (!simulation.count("NonlinearOptions"))
        {
            throw std::runtime_error("\"NonlinearOptions\" needs to be present for a "
                                     "SolidMechanics simulation");
        }
        return std::make_unique<SolidMechanicsModule>(mesh, material, simulation);
    }
    else if (module_type == "PlaneMechanics")
    {
        if (!simulation.count("NonlinearOptions"))
        {
            throw std::runtime_error("\"NonlinearOptions\" needs to be present for a "
                                     "SolidMechanics simulation");
        }
        return std::make_unique<PlaneMechanicsModule>(mesh, material, simulation);
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
