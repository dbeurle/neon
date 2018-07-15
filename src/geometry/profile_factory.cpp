
#include "geometry/profile_factory.hpp"
#include "io/json.hpp"

namespace neon::geometry
{
std::unique_ptr<profile> make_profile(json const& profile_data)
{
    if (profile_data.find("profile") == profile_data.end())
    {
        throw std::domain_error("Please provide a \"profile\" field");
    }

    if (profile_data["profile"] == "rectangle")
    {
        return std::make_unique<rectangle>(profile_data);
    }
    else if (profile_data["profile"] == "circle")
    {
        return std::make_unique<circle>(profile_data);
    }

    throw std::domain_error("A valid profile was not specified.  Valid profiles are \"rectangle\"");

    return nullptr;
}
}
