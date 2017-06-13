
#pragma once

#include "PreprocessorExceptions.hpp"

namespace neon
{
class MaterialPropertyException : public Exception
{
public:
    MaterialPropertyException(std::string material_property) : material_property(material_property)
    {
    }

    const char* what() const noexcept
    {
        std::cout << "\n!! Error: " << material_property << " has not been provided in "
                  << input_file;
        return nullptr;
    }

protected:
    std::string material_property;
};
}
