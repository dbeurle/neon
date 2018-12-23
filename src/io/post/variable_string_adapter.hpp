
#pragma once

/// @file

#include "constitutive/variable_types.hpp"

#include <string>

namespace neon::variable
{
/// \return variant enum if the name is valid, other std::nullopt
types convert(std::string const& name);
/// \return convert a scalar enum to a string for pretty file name
char const* convert(variable::scalar);
/// \return convert a second order tensor enum to a string
char const* convert(variable::second);
}
