
#pragma once

#include "mesh/basic_mesh.hpp"

#include <map>
#include <memory>
#include <string>
#include <utility>

#include "io/json_forward.hpp"

namespace neon
{
class abstract_module;

std::unique_ptr<abstract_module> make_module(
    json const& simulation, std::map<std::string, std::pair<basic_mesh, json>> const& mesh_store);
}
