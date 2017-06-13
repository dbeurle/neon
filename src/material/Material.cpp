
#include "Material.hpp"

#include "MaterialExceptions.hpp"

namespace neon
{
Material::Material(std::string const& name) : myname(name)
{
    if (myname.empty()) throw MaterialPropertyException("Name");
}
}
