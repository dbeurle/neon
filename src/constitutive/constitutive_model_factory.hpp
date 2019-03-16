
#pragma once

/// @file

#include "io/json_forward.hpp"

#include "internal_variables_forward.hpp"

#include "constitutive/constitutive_model_forward.hpp"
#include "constitutive/constitutive_model_alias.hpp"

#include <memory>

namespace neon
{
namespace mechanics::solid
{
auto make_constitutive_model(std::shared_ptr<internal_variables_t>& variables,
                             json const& material_data,
                             json const& simulation_data) -> std::unique_ptr<constitutive_model>;
}
namespace mechanics::plane
{
auto make_constitutive_model(std::shared_ptr<internal_variables_t>& variables,
                             json const& material_data,
                             json const& simulation_data) -> std::unique_ptr<constitutive_model>;
}

namespace diffusion
{
auto make_constitutive_model(std::shared_ptr<internal_variables_t>& variables,
                             json const& material_data,
                             json const& simulation_data) -> std::unique_ptr<constitutive_model>;
}
}
