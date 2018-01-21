
#pragma once

#include "mesh/BasicMesh.hpp"

#include <map>
#include <memory>
#include <string>

#include "io/json.hpp"

namespace neon
{
class AbstractModule;

std::unique_ptr<AbstractModule> make_module(
    json const& simulation,
    std::map<std::string, std::pair<BasicMesh, json>> const& mesh_store);
}
