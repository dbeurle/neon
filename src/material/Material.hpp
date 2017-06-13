
#include <string>

#pragma once

namespace neon
{
class Material
{
public:
    Material(std::string const& name);

    virtual ~Material() = default;

    std::string const& name() const { return myname; }

protected:
    std::string myname;
};
}
