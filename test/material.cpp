
#include <catch2/catch.hpp>

#include "material/isotropic_elastic_plastic.hpp"
#include "material/isotropic_elastic_property.hpp"
#include "material/micromechanical_elastomer.hpp"
#include "material/linear_diffusion.hpp"
#include "io/json.hpp"

#include <stdexcept>

using namespace neon;

TEST_CASE("Basic material")
{
    SECTION("No name on construction")
    {
        REQUIRE_THROWS_AS(material_property(json::parse("{\"nasame\": \"steel\"}")),
                          std::domain_error);
    }
    SECTION("Density and specific heat exception")
    {
        material_property basic_material(json::parse("{\"name\": \"steel\"}"));

        REQUIRE_THROWS_AS(basic_material.initial_density(), std::domain_error);
        REQUIRE_THROWS_AS(basic_material.specific_heat(), std::domain_error);
    }
}
TEST_CASE("Linear elastic material")
{
    SECTION("Linear elastic properties")
    {
        json material_data{{"name", "steel"}, {"elastic_modulus", 200.0e9}, {"poissons_ratio", 0.3}};

        isotropic_elastic_property linear_elastic(material_data);

        REQUIRE(linear_elastic.name() == "steel");

        REQUIRE(linear_elastic.Poissons_ratio() == Approx(0.3));
        REQUIRE(linear_elastic.elastic_modulus() == Approx(200.0e9));

        auto const E = linear_elastic.elastic_modulus();
        auto const v = linear_elastic.Poissons_ratio();

        // Check Lamé parameters
        REQUIRE(linear_elastic.lambda() == Approx(E * v / ((1.0 + v) * (1.0 - 2.0 * v))));
        REQUIRE(linear_elastic.mu() == Approx(E / (2.0 * (1.0 + v))));

        auto const [lambda, shear_modulus] = linear_elastic.Lame_parameters();

        REQUIRE(lambda == Approx(linear_elastic.lambda()));
        REQUIRE(shear_modulus == Approx(linear_elastic.mu()));
    }
    SECTION("Incompressible elastic properties")
    {
        json material_data{{"name", "steel"}, {"bulk_modulus", 200.0e9}, {"shear_modulus", 100.0e6}};

        isotropic_elastic_property linear_elastic(material_data);

        REQUIRE(linear_elastic.bulk_modulus() == Approx(200.0e9));
        REQUIRE(linear_elastic.shear_modulus() == Approx(100.0e6));
    }
    SECTION("Incorrect elastic properties")
    {
        json material_data{{"name", "steel"}, {"elastic_modulus", 200.0e9}, {"PssonsRatio", 0.3}};

        REQUIRE_THROWS_AS(isotropic_elastic_property(material_data), std::domain_error);
    }
}
TEST_CASE("Perfect plastic material")
{
    json material_data{{"name", "steel"},
                       {"elastic_modulus", 200.0e9},
                       {"poissons_ratio", 0.3},
                       {"yield_stress", 200.0e6}};

    isotropic_elastic_plastic perfect_plastic_elastic(material_data);

    REQUIRE(perfect_plastic_elastic.name() == "steel");

    REQUIRE(perfect_plastic_elastic.yield_stress(0.0) == Approx(200.0e6));
    REQUIRE(perfect_plastic_elastic.yield_stress(1.0) == Approx(200.0e6));

    REQUIRE(perfect_plastic_elastic.hardening_modulus(0.0) == Approx(0.0));
    REQUIRE(perfect_plastic_elastic.hardening_modulus(1.0) == Approx(0.0));
}
TEST_CASE("Isotropic hardening")
{
    json material_data{{"name", "steel"},
                       {"elastic_modulus", 200.0e9},
                       {"poissons_ratio", 0.3},
                       {"yield_stress", 200.0e6},
                       {"isotropic_hardening_modulus", 400.0e6},
                       {"isotropic_kinematic_modulus", 100.0e6}};

    isotropic_elastic_plastic iso_plastic_elastic(material_data);

    REQUIRE(iso_plastic_elastic.name() == "steel");

    REQUIRE(iso_plastic_elastic.yield_stress(0.0) == Approx(200.0e6));
    REQUIRE(iso_plastic_elastic.yield_stress(1.0) == Approx(600.0e6));

    REQUIRE(iso_plastic_elastic.hardening_modulus(0.0) == Approx(400.0e6));
    REQUIRE(iso_plastic_elastic.hardening_modulus(1.0) == Approx(400.0e6));

    REQUIRE(iso_plastic_elastic.kinematic_modulus(0.0) == Approx(100.0e6));
}
TEST_CASE("Missing yield stress")
{
    json material_data{{"name", "steel"},
                       {"elastic_modulus", 200.0e9},
                       {"poissons_ratio", 0.3},
                       {"YieldStrs", 200.0e6},
                       {"isotropic_hardening_modulus", 400.0e6},
                       {"isotropic_kinematic_modulus", 100.0e6}};

    REQUIRE_THROWS_AS(isotropic_elastic_plastic(material_data), std::domain_error);
}
TEST_CASE("Micromechanical elastomer")
{
    json material_data{{"name", "rubber"},
                       {"elastic_modulus", 10.0e6},
                       {"poissons_ratio", 0.45},
                       {"segments_per_chain", 70},
                       {"cure_time", 100},
                       {"recombination_probability", 1.0e-6},
                       {"scission_probability", 1.0e-6}};

    ageing_micromechanical_elastomer network(material_data);

    SECTION("basic data check")
    {
        REQUIRE(network.scission_probability() == Approx(1.0e-6));
        REQUIRE(network.recombination_probability() == Approx(1.0e-6));
        REQUIRE(network.segments_per_chain() == Approx(70.0));
        REQUIRE(network.initial_inactive_shear_modulus() == Approx(0.0).margin(1E-5));
    }
    SECTION("network evolution test")
    {
        // Create initial conditions
        vector5 z(5);
        // Active set shear modulus
        z(0) = network.shear_modulus();
        // Inactive set
        z(1) = 0.0;
        // Reduction factor
        z(2) = 1.0;
        // Active segments
        z(3) = network.segments_per_chain();
        // Inactive segments
        z(4) = 0.0;

        for (int i = 0; i < 3; ++i)
        {
            z = network.integrate(z, 0.1);

            REQUIRE(z(0) > network.shear_modulus());
            REQUIRE(z(1) > 0.0);
            REQUIRE(z(2) < 1.0);
            REQUIRE(z(3) < network.segments_per_chain());
            REQUIRE(z(4) >= 0.0);
        }
    }
}
TEST_CASE("Diffusion material")
{
    json material_data{{"name", "steel"},
                       {"density", 7800.0},
                       {"conductivity", 300.0},
                       {"specific_heat", 280.0}};

    linear_diffusion material(material_data);

    REQUIRE(material.initial_density() == Approx(7800.0));
    REQUIRE(material.conductivity() == Approx(300.0));
    REQUIRE(material.specific_heat() == Approx(280.0));
}
