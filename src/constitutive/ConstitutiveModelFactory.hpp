
#pragma once

#include <memory>

#include <json/forwards.h>

namespace neon
{
class ConstitutiveModel;
class InternalVariables;

namespace solid
{
std::unique_ptr<ConstitutiveModel> make_constitutive_model(InternalVariables& variables,
                                                           Json::Value const& material_data,
                                                           Json::Value const& simulation_data);
}
namespace diffusion
{
std::unique_ptr<ConstitutiveModel> make_constitutive_model(InternalVariables& variables,
                                                           Json::Value const& material_data,
                                                           Json::Value const& simulation_data);
}
}
