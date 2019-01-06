
#pragma once

/// @file

#include "mesh/basic_mesh.hpp"
#include "io/json_forward.hpp"

#include <map>
#include <memory>
#include <string>
#include <utility>

namespace neon
{
namespace geometry
{
class profile;
}

class abstract_module;

std::unique_ptr<abstract_module> make_module(
    json const& simulation,
    std::map<std::string, std::pair<basic_mesh, json>> const& mesh_store,
    std::map<std::string, std::unique_ptr<geometry::profile>> const& profile_store);
}
