
#include "ModuleFactory.hpp"

#include "DiffusionModule.hpp"
#include "SolidModule.hpp"

#include <json/value.h>

namespace neon
{
std::unique_ptr<AbstractModule> SimulationControl::make_module(Json::Value const& simulation) const
{
    if (auto const& module_type = simulation["Type"].asString(); module_type == "SolidMechanics")
    {
        return std::make_unique<SolidModule>();
    }
    else if (module_type == "TemperatureDiffusion")
    {
        return std::make_unique<LinearDiffusionModule>();
    }
    else
    {
        throw std::runtime_error("Module " + module_type + " is not recognised\n");
    }
    return nullptr;
}
}
