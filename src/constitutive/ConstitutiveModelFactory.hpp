
#pragma once

#include <memory>

#include <json/forwards.h>

namespace neon
{
template <int D>
class InternalVariables;

template <int D>
class ConstitutiveModel;

namespace solid
{
using ConstitutiveModel = neon::ConstitutiveModel<3>;
using InternalVariables = neon::InternalVariables<3>;

std::unique_ptr<ConstitutiveModel> make_constitutive_model(InternalVariables& variables,
                                                           Json::Value const& material_data,
                                                           Json::Value const& simulation_data);
}
namespace diffusion
{
using ConstitutiveModel = neon::ConstitutiveModel<3>;
using InternalVariables = neon::InternalVariables<3>;

std::unique_ptr<ConstitutiveModel> make_constitutive_model(InternalVariables& variables,
                                                           Json::Value const& material_data,
                                                           Json::Value const& simulation_data);
}
}
