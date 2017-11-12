
#pragma once

#include <memory>

#include <json/forwards.h>

#include "InternalVariablesForwards.hpp"

namespace neon
{
template <int spatial_dimension, int voigt_dimension>
class ConstitutiveModel;

namespace mech::solid
{
using ConstitutiveModel = neon::ConstitutiveModel<3, 6>;

std::unique_ptr<ConstitutiveModel> make_constitutive_model(InternalVariables& variables,
                                                           Json::Value const& material_data,
                                                           Json::Value const& simulation_data);
}
namespace diffusion
{
using ConstitutiveModel = neon::ConstitutiveModel<3, 3>;

std::unique_ptr<ConstitutiveModel> make_constitutive_model(InternalVariables& variables,
                                                           Json::Value const& material_data,
                                                           Json::Value const& simulation_data);
}
}
