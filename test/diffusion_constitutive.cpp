
#include <catch2/catch.hpp>

#include "constitutive/internal_variables.hpp"
#include "constitutive/constitutive_model_factory.hpp"
#include "constitutive/thermal/isotropic_diffusion.hpp"

#include "exceptions.hpp"
#include "io/json.hpp"

std::string json_thermal_diffusion()
{
    return "{\"name\": \"steel\", \"conductivity\": 386.0, \"density\": 7800.0, \"specific_heat\": "
           "390.0}";
}

constexpr auto internal_variable_size = 2;
constexpr auto ZERO_MARGIN = 1.0e-5;

using neon::json;

TEST_CASE("Thermal isotropic model")
{
    using namespace neon::diffusion;

    auto variables = std::make_shared<internal_variables_t>(internal_variable_size);

    auto thermal = make_constitutive_model(variables,
                                           json::parse(json_thermal_diffusion()),
                                           json::parse("{\"constitutive\" : {\"name\": "
                                                       "\"isotropic_diffusion\"} }"));

    thermal->update_internal_variables(1.0);

    SECTION("Sanity checks")
    {
        REQUIRE(thermal->is_symmetric());
        REQUIRE(!thermal->is_finite_deformation());
        REQUIRE(thermal->intrinsic_material().name() == "steel");

        REQUIRE(variables->has(neon::variable::second::conductivity));
    }
}
