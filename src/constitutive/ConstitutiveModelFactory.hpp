
#pragma once

#include <memory>

#include <json/forwards.h>

namespace neon
{
class InternalVariables;

template <int D>
struct ConstitutiveModel;

namespace solid
{
using ConstitutiveModel = neon::ConstitutiveModel<3>;

std::unique_ptr<solid::ConstitutiveModel> make_constitutive_model(InternalVariables& variables,
                                                                  Json::Value const& material_data,
                                                                  Json::Value const& simulation_data);
}
namespace diffusion
{
using ConstitutiveModel = neon::ConstitutiveModel<3>;

std::unique_ptr<diffusion::ConstitutiveModel> make_constitutive_model(
    InternalVariables& variables, Json::Value const& material_data, Json::Value const& simulation_data);
}
}
