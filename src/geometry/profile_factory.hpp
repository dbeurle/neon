
#pragma once

#include "geometry/profile.hpp"
#include "io/json_forward.hpp"

#include <memory>

namespace neon::geometry
{
std::unique_ptr<profile> make_profile(json const& profile_data);
}
