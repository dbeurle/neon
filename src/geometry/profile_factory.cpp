
#include "geometry/profile_factory.hpp"
#include "io/json.hpp"

namespace neon::geometry
{
std::unique_ptr<profile> make_profile(json const& profile_data)
{
    if (profile_data.find("type") == profile_data.end())
    {
        throw std::domain_error("Please provide a \"type\" field");
    }

    if (auto const& profile_type = profile_data["type"]; profile_type == "rectangle")
    {
        return std::make_unique<rectangle>(profile_data);
    }
    else if (profile_type == "circle")
    {
        return std::make_unique<circle>(profile_data);
    }
    else if (profile_type == "hollow_circle")
    {
        return std::make_unique<hollow_circle>(profile_data);
    }
    else if (profile_type == "rectangular_angle")
    {
        return std::make_unique<rectangular_angle>(profile_data);
    }

    throw std::domain_error("A valid profile type was not specified.  Valid profiles are "
                            "\"rectangle\"");

    return nullptr;
}
}
