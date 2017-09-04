
#pragma once

#include <memory>

#include <json/forwards.h>

namespace neon
{
class AbstractModule;

std::unique_ptr<AbstractModule> make_module(Json::Value const& simulation);
}
