
#define CATCH_CONFIG_MAIN

#include "catch.hpp"

#include <range/v3/numeric.hpp>

#include "material/IsotropicElasticPlastic.hpp"
#include "material/LinearElastic.hpp"
#include "material/MicromechanicalElastomer.hpp"

#include "material/LinearDiffusion.hpp"

#include "Exceptions.hpp"

#include "io/json.hpp"

using namespace neon;

TEST_CASE("Linear elastic material", "[LinearElastic]")
{
    SECTION("Linear elastic properties")
    {
        json material_data{{"Name", "steel"}, {"ElasticModulus", 200.0e9}, {"PoissonsRatio", 0.3}};

        LinearElastic linear_elastic(material_data);

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
        json material_data{{"Name", "steel"}, {"BulkModulus", 200.0e9}, {"ShearModulus", 100.0e6}};

        LinearElastic linear_elastic(material_data);

        REQUIRE(linear_elastic.bulk_modulus() == Approx(200.0e9));
        REQUIRE(linear_elastic.shear_modulus() == Approx(100.0e6));
    }
    SECTION("Incorrect elastic properties")
    {
        json material_data{{"Name", "steel"}, {"ElasticModulus", 200.0e9}, {"PssonsRatio", 0.3}};

        REQUIRE_THROWS_AS(LinearElastic(material_data), MaterialPropertyException);
    }
}
TEST_CASE("Perfect plastic material", "[PerfectPlasticElastic]")
{
    json material_data{{"Name", "steel"},
                       {"ElasticModulus", 200.0e9},
                       {"PoissonsRatio", 0.3},
                       {"YieldStress", 200.0e6}};

    IsotropicElasticPlastic perfect_plastic_elastic(material_data);

    REQUIRE(perfect_plastic_elastic.name() == "steel");

    REQUIRE(perfect_plastic_elastic.yield_stress(0.0) == Approx(200.0e6));
    REQUIRE(perfect_plastic_elastic.yield_stress(1.0) == Approx(200.0e6));

    REQUIRE(perfect_plastic_elastic.hardening_modulus(0.0) == Approx(0.0));
    REQUIRE(perfect_plastic_elastic.hardening_modulus(1.0) == Approx(0.0));
}
TEST_CASE("Isotropic hardening", "[IsotropicPlasticElastic]")
{
    json material_data{{"Name", "steel"},
                       {"ElasticModulus", 200.0e9},
                       {"PoissonsRatio", 0.3},
                       {"YieldStress", 200.0e6},
                       {"IsotropicHardeningModulus", 400.0e6},
                       {"IsotropicKinematicModulus", 100.0e6}};

    IsotropicElasticPlastic iso_plastic_elastic(material_data);

    REQUIRE(iso_plastic_elastic.name() == "steel");

    REQUIRE(iso_plastic_elastic.yield_stress(0.0) == Approx(200.0e6));
    REQUIRE(iso_plastic_elastic.yield_stress(1.0) == Approx(600.0e6));

    REQUIRE(iso_plastic_elastic.hardening_modulus(0.0) == Approx(400.0e6));
    REQUIRE(iso_plastic_elastic.hardening_modulus(1.0) == Approx(400.0e6));

    REQUIRE(iso_plastic_elastic.kinematic_modulus(0.0) == Approx(100.0e6));
}
TEST_CASE("Missing yield stress", "[IsotropicPlasticElastic]")
{
    json material_data{{"Name", "steel"},
                       {"ElasticModulus", 200.0e9},
                       {"PoissonsRatio", 0.3},
                       {"YieldStrs", 200.0e6},
                       {"IsotropicHardeningModulus", 400.0e6},
                       {"IsotropicKinematicModulus", 100.0e6}};

    REQUIRE_THROWS_AS(IsotropicElasticPlastic(material_data), std::domain_error);
}
TEST_CASE("Micromechanical elastomer", "[StochasticMicromechanicalElastomer]")
{
    json material_data{{"Name", "rubber"},
                       {"ElasticModulus", 10.0e6},
                       {"PoissonsRatio", 0.45},
                       {"Segments",
                        {{"Groups", 5},
                         {"Average", 50},
                         {"StandardDeviation", 10},
                         {"ScissionLikelihood", 0.0001}}}};

    StochasticMicromechanicalElastomer elastomer(material_data);

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
TEST_CASE("Diffusion material", "[LinearDiffusion]")
{
    json material_data{{"Name", "steel"},
                       {"Density", 7800.0},
                       {"Conductivity", 300.0},
                       {"SpecificHeat", 280.0}};

    LinearDiffusion linear_diffusion(material_data);

    REQUIRE(linear_diffusion.initial_density() == Approx(7800.0));
    REQUIRE(linear_diffusion.conductivity_coefficient() == Approx(300.0));
    REQUIRE(linear_diffusion.specific_heat_coefficient() == Approx(280.0));
}
