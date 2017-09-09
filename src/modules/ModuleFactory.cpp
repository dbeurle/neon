
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

    std::cout << std::string(4, ' ') << "Module \"" << simulation["Type"].asString() << "\"\n";
    std::cout << std::string(4, ' ') << "Solution \"" << simulation["Solution"].asString() << "\"\n";

    auto const & [ mesh, material ] = simulation_mesh->second;

    if (auto const& module_type = simulation["Type"].asString(); module_type == "SolidMechanics")
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
        return std::make_unique<LinearDiffusionModule>(mesh, material, simulation);
    }
    else
    {
        throw std::runtime_error("Module " + module_type + " is not recognised\n");
    }
    return nullptr;
}
}
