
#pragma once

#include <memory>

#include "io/json_forward.hpp"

#include "InternalVariablesForwards.hpp"

namespace neon
{
template <int spatial_dimension, int voigt_dimension>
class ConstitutiveModel;

namespace mechanical::solid
{
using ConstitutiveModel = neon::ConstitutiveModel<3, 6>;

std::unique_ptr<ConstitutiveModel> make_constitutive_model(std::shared_ptr<InternalVariables>& variables,
                                                           json const& material_data,
                                                           json const& simulation_data);
}
namespace mechanical::plane
{
using ConstitutiveModel = neon::ConstitutiveModel<2, 3>;

std::unique_ptr<ConstitutiveModel> make_constitutive_model(std::shared_ptr<InternalVariables>& variables,
                                                           json const& material_data,
                                                           json const& simulation_data);
}

namespace diffusion
{
using ConstitutiveModel = neon::ConstitutiveModel<3, 3>;

std::unique_ptr<ConstitutiveModel> make_constitutive_model(std::shared_ptr<InternalVariables>& variables,
                                                           json const& material_data,
                                                           json const& simulation_data);
}
}
