
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
std::unique_ptr<constitutive_model> make_constitutive_model(
    std::shared_ptr<internal_variables_t>& variables,
    json const& material_data,
    json const& simulation_data);
}
namespace mechanics::plane
{
std::unique_ptr<constitutive_model> make_constitutive_model(
    std::shared_ptr<internal_variables_t>& variables,
    json const& material_data,
    json const& simulation_data);
}

namespace diffusion
{
std::unique_ptr<constitutive_model> make_constitutive_model(
    std::shared_ptr<internal_variables_t>& variables,
    json const& material_data,
    json const& simulation_data);
}
}
