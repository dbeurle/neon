#define CATCH_CONFIG_MAIN // This tells Catch to provide a main() - only do this in one
                          // cpp file

#include "catch.hpp"

#include <json/json.h>
#include <range/v3/numeric.hpp>

#include "material/IsotropicElasticPlastic.hpp"
#include "material/LinearElastic.hpp"
#include "material/MicromechanicalElastomer.hpp"

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

std::string micromechanical_input()
{
    return "{\"Name\" : \"rubber\", "
           "\"ElasticModulus\" : 10.0e6, "
           "\"PoissonsRatio\" : 0.45, "
           "\"Segments\" : { "
           "\"Groups\" : 5, "
           "\"Average\" : 50, "
           "\"StandardDeviation\" : 10, "
           "\"ScissionLikelihood\" : 0.0001}}";
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
    auto const v = linear_elastic.Poissons_ratio();

    // Check Lamé parameters
    REQUIRE(linear_elastic.lambda() == Approx(E * v / ((1.0 + v) * (1.0 - 2.0 * v))));
    REQUIRE(linear_elastic.mu() == Approx(E / (2.0 * (1.0 + v))));

    auto const[lambda, shear_modulus] = linear_elastic.Lame_parameters();

    REQUIRE(lambda == Approx(linear_elastic.lambda()));
    REQUIRE(shear_modulus == Approx(linear_elastic.mu()));

    REQUIRE(linear_elastic.bulk_modulus() == Approx(lambda + 2.0 / 3.0 * shear_modulus));
    REQUIRE(linear_elastic.shear_modulus() == Approx(shear_modulus));
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
TEST_CASE("Micromechanical elastomer", "[MicromechanicalElastomer]")
{
    Json::Value material_data;
    Json::Reader material_file;

    REQUIRE(material_file.parse(micromechanical_input().c_str(), material_data));
    MicromechanicalElastomer elastomer(material_data);

    // Initial chains and segment groups
    auto const chain_group_initial = elastomer.chain_groups();
    auto const segment_group_initial = elastomer.segment_groups();

    REQUIRE(elastomer.groups() == 5);

    SECTION("Perform time step")
    {
        auto const chain_group = elastomer.update_chains(elastomer.chain_groups(), 100.0);
        auto const shear_moduli = elastomer.compute_shear_moduli(chain_group);
        auto const segment_group = elastomer.segment_groups();

        REQUIRE(chain_group.size() == segment_group.size());

        // Ensure we have a decrease in the number of chains
        for (auto i = 0; i < chain_group.size(); i++)
        {
            REQUIRE(chain_group.at(i) < chain_group_initial.at(i));
            REQUIRE(segment_group.at(i) == segment_group_initial.at(i));
            REQUIRE(shear_moduli.at(i) > 0.0);
        }
    }
}
