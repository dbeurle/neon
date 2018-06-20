
#pragma once

#include "constitutive/variable_types.hpp"

#include <optional>

/// \file variable_string_adapter.hpp

namespace neon::variable
{
/// \return scalar enum if the name is valid, other std::nullopt
std::optional<variable::scalar> is_scalar(std::string const& name);
/// \return vector enum if the name is valid, other std::nullopt
std::optional<variable::vector> is_vector(std::string const& name);
/// \return second enum if the name is valid, other std::nullopt
std::optional<variable::second> is_second_order_tensor(std::string const& name);
}
