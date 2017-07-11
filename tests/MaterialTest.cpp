#define CATCH_CONFIG_MAIN // This tells Catch to provide a main() - only do this in one cpp file

#include "catch.hpp"

#include <json/json.h>

#include "material/IsotropicElasticPlastic.hpp"
#include "material/LinearElastic.hpp"

using namespace neon;

std::string linear_material_input()
{
    return "{\"Name\": \"steel\", \"ElasticModulus\": 200.0e9, \"PoissonsRatio\": 0.3}";
}

std::string perfect_plastic_input()
{
    std::string linear = linear_material_input();
    linear.pop_back();
    return linear + ",\"YieldStress\": 200.0e6}";
}

std::string isotropic_plastic_input()
{
    std::string plastic = perfect_plastic_input();
    plastic.pop_back();
    return plastic + ",\"IsotropicHardeningModulus\": 400.0e6}";
}

TEST_CASE("Linear elastic material", "[LinearElastic]")
{
    Json::Value material_data;
    Json::Reader material_file;

    REQUIRE(material_file.parse(linear_material_input().c_str(), material_data));

    LinearElastic linear_elastic(material_data);

    REQUIRE(linear_elastic.name() == "steel");

    REQUIRE(linear_elastic.Poissons_ratio() == Approx(0.3));
    REQUIRE(linear_elastic.elastic_modulus() == Approx(200.0e9));

    auto const E = linear_elastic.elastic_modulus();
    auto const ν = linear_elastic.Poissons_ratio();

    // Check Lamé parameters
    REQUIRE(linear_elastic.lambda() == Approx(E * ν / ((1.0 + ν) * (1.0 - 2.0 * ν))));
    REQUIRE(linear_elastic.mu() == Approx(E / (2.0 * (1.0 + ν))));

    auto const[λ, μ] = linear_elastic.Lame_parameters();

    REQUIRE(λ == Approx(linear_elastic.lambda()));
    REQUIRE(μ == Approx(linear_elastic.mu()));

    REQUIRE(linear_elastic.bulk_modulus() == Approx(λ + 2.0 / 3.0 * μ));
    REQUIRE(linear_elastic.shear_modulus() == Approx(μ));
}

TEST_CASE("Perfect plastic material", "[PerfectPlasticElastic]")
{
    Json::Value material_data;
    Json::Reader material_file;

    REQUIRE(material_file.parse(perfect_plastic_input().c_str(), material_data));

    IsotropicElasticPlastic perfect_plastic_elastic(material_data);

    REQUIRE(perfect_plastic_elastic.name() == "steel");

    REQUIRE(perfect_plastic_elastic.yield_stress(0.0) == Approx(200.0e6));
    REQUIRE(perfect_plastic_elastic.yield_stress(1.0) == Approx(200.0e6));

    REQUIRE(perfect_plastic_elastic.hardening_modulus(0.0) == Approx(0.0));
    REQUIRE(perfect_plastic_elastic.hardening_modulus(1.0) == Approx(0.0));
}

TEST_CASE("Isotropic hardening", "[IsotropicPlasticElastic]")
{
    Json::Value material_data;
    Json::Reader material_file;

    REQUIRE(material_file.parse(isotropic_plastic_input().c_str(), material_data));

    IsotropicElasticPlastic iso_plastic_elastic(material_data);

    REQUIRE(iso_plastic_elastic.name() == "steel");

    REQUIRE(iso_plastic_elastic.yield_stress(0.0) == Approx(200.0e6));
    REQUIRE(iso_plastic_elastic.yield_stress(1.0) == Approx(600.0e6));

    REQUIRE(iso_plastic_elastic.hardening_modulus(0.0) == Approx(400.0e6));
    REQUIRE(iso_plastic_elastic.hardening_modulus(1.0) == Approx(400.0e6));
}
