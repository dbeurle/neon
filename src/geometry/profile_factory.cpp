
#include "geometry/profile_factory.hpp"
#include "io/json.hpp"

namespace neon::geometry
{
std::unique_ptr<profile> make_profile(json const& profile_data)
{
    if (profile_data.count("Geometry") == 0)
    {
        throw std::domain_error("Please provide a \"Geometry\" section");
    }

    if (profile_data["Geometry"] == "Rectangle")
    {
        return std::make_unique<rectangular_bar>(1.0, 1.0);
    }
    else if (profile_data["Geometry"] == "HollowRectangle")
    {
        return std::make_unique<hollow_rectangular_bar>(1.0, 1.0);
    }
    return nullptr;
}
}
