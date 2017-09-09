
#pragma once

#include "mesh/BasicMesh.hpp"

#include <map>
#include <memory>
#include <string>

#include <json/value.h>

namespace neon
{
class AbstractModule;

std::unique_ptr<AbstractModule> make_module(
    Json::Value const& simulation,
    std::map<std::string, std::pair<BasicMesh, Json::Value>> const& mesh_store);
}
