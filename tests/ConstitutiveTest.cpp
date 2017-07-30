#define CATCH_CONFIG_MAIN // This tells Catch to provide a main() - only do this in one
                          // cpp file
#include "catch.hpp"

#include <iostream>
#include <json/json.h>

#include "constitutive/InternalVariables.hpp"

#include "constitutive/AffineMicrosphere.hpp"
#include "constitutive/HyperElasticPlastic.hpp"
#include "constitutive/NeoHooke.hpp"

using namespace neon;

std::string json_input_file()
{
    return "{\"Name\": \"rubber\", \"ElasticModulus\": 1.0, \"PoissonsRatio\": 2.0}";
}

constexpr auto internal_variable_size = 2;

TEST_CASE("Neo-Hookean model", "[NeoHooke]")
{
    InternalVariables variables(internal_variable_size);

    // Add the required variables for an updated Lagrangian formulation
    variables.add(InternalVariables::Tensor::DeformationGradient,
                  InternalVariables::Tensor::Cauchy);
    variables.add(InternalVariables::Scalar::DetF);

    // Create a json reader object from a string
    auto input_data = json_input_file();

    Json::Value material_data;
    Json::Reader reader;

    REQUIRE(reader.parse(input_data.c_str(), material_data));

    NeoHooke neo_hooke(variables, material_data);

    // Get the tensor variables
    auto[F_list, σ_list] = variables(InternalVariables::Tensor::DeformationGradient,
                                     InternalVariables::Tensor::Cauchy);

    auto& J_list = variables(InternalVariables::Scalar::DetF);

    // Ensure the internal variables are allocated correctly
    REQUIRE(F_list.size() == internal_variable_size);
    REQUIRE(σ_list.size() == internal_variable_size);
    REQUIRE(J_list.size() == internal_variable_size);

    // Fill with identity matrix
    for (auto& F : F_list) F = Matrix3::Identity();
    for (auto& J : J_list) J = 1.0;

    neo_hooke.update_internal_variables(1.0);

    SECTION("Update of internal variables")
    {
        for (auto& σ : σ_list) REQUIRE(σ.norm() == Approx(0.0));
    }
    SECTION("Check of continuum tangent")
    {
        // Get the matrix variable
        auto& material_tangents = variables(InternalVariables::Matrix::TruesdellModuli);

        for (auto const& C : material_tangents)
        {
            // Check a few of the numerical values
            REQUIRE(C(0, 0) == Approx(0.111111));
            REQUIRE(C(0, 1) == Approx(-0.222222));
            REQUIRE(C(0, 2) == Approx(-0.222222));
            REQUIRE(C(1, 1) == Approx(0.111111));
            REQUIRE(C(2, 2) == Approx(0.111111));
            REQUIRE(C(3, 3) == Approx(0.1666666667));
            REQUIRE(C(4, 4) == Approx(0.1666666667));
            REQUIRE(C(5, 5) == Approx(0.1666666667));

            // Ensure symmetry is correct
            REQUIRE((C - C.transpose()).norm() == Approx(0.0));
        }
    }
}
TEST_CASE("Affine microsphere model", "[AffineMicrosphere]")
{
    InternalVariables variables(internal_variable_size);

    // Add the required variables for an updated Lagrangian formulation
    variables.add(InternalVariables::Tensor::DeformationGradient,
                  InternalVariables::Tensor::Cauchy);

    variables.add(InternalVariables::Scalar::DetF);

    // Create a json reader object from a string
    std::string input_data =
        "{\"Name\": \"rubber\", \"ElasticModulus\": 1.0, \"PoissonsRatio\": 2.0, "
        "\"SegmentsPerChain\" : 25}";

    Json::Value material_data;
    Json::Reader reader;

    REQUIRE(reader.parse(input_data.c_str(), material_data));

    AffineMicrosphere affine(variables, material_data);

    // Get the tensor variables
    auto[F_list, σ_list] = variables(InternalVariables::Tensor::DeformationGradient,
                                     InternalVariables::Tensor::Cauchy);

    auto& J_list = variables(InternalVariables::Scalar::DetF);

    for (auto& J : J_list) J = 1.0;

    auto& material_tangents = variables(InternalVariables::Matrix::TruesdellModuli);

    SECTION("Affine model under no load")
    {
        // Fill with identity matrix
        for (auto& F : F_list) F = Matrix3::Identity();

        affine.update_internal_variables(1.0);

        // Ensure symmetry is correct
        for (auto const& C : material_tangents)
        {
            REQUIRE((C - C.transpose()).norm() == Approx(0.0));
        }

        for (auto& σ : σ_list) REQUIRE(σ.norm() == Approx(0.0));
    }
    SECTION("Affine model under uniaxial load")
    {
        for (auto& F : F_list)
        {
            F(0, 0) = 1.1;
            F(1, 1) = 1.0 / std::sqrt(1.1);
            F(2, 2) = 1.0 / std::sqrt(1.1);
        }

        affine.update_internal_variables(1.0);

        // Ensure symmetry is correct
        for (auto const& C : material_tangents)
        {
            REQUIRE((C - C.transpose()).norm() == Approx(0.0));
        }

        for (auto& σ : σ_list) REQUIRE(σ.norm() != Approx(0.0));
    }
}
TEST_CASE("J2 plasticity model", "[J2Plasticity]")
{
    // Create a json reader object from a string
    std::string input_data =
        "{\"Name\": \"steel\", \"ElasticModulus\": 200.0e9, \"PoissonsRatio\": 0.3, "
        "\"YieldStress\": 200.0e6, \"IsotropicHardeningModulus\": 400.0e6}";

    Json::Value material_data;
    Json::Reader reader;

    REQUIRE(reader.parse(input_data.c_str(), material_data));

    InternalVariables variables(internal_variable_size);

    // Add the required variables for an updated Lagrangian formulation
    variables.add(InternalVariables::Tensor::DisplacementGradient,
                  InternalVariables::Tensor::Cauchy);
    variables.add(InternalVariables::Scalar::DetF);

    J2Plasticity j2plasticity(variables, material_data);

    // Get the tensor variables
    auto[H_list, σ_list] = variables(InternalVariables::Tensor::DisplacementGradient,
                                     InternalVariables::Tensor::Cauchy);

    auto& J_list = variables(InternalVariables::Scalar::DetF);

    auto& material_tangents = variables(InternalVariables::Matrix::TruesdellModuli);

    for (auto& H : H_list) H = Matrix3::Zero();
    for (auto& J : J_list) J = 1.0;

    SECTION("Internal variables allocation")
    {
        REQUIRE(variables.has(InternalVariables::Scalar::VonMisesStress));
        REQUIRE(variables.has(InternalVariables::Scalar::EffectivePlasticStrain));
        REQUIRE(variables.has(InternalVariables::Tensor::LinearisedStrain));
        REQUIRE(variables.has(InternalVariables::Tensor::LinearisedPlasticStrain));
        REQUIRE(variables.has(InternalVariables::Matrix::TruesdellModuli));
    }
    SECTION("No load")
    {
        j2plasticity.update_internal_variables(1.0);

        // Ensure symmetry is correct
        for (auto const& C : material_tangents)
        {
            REQUIRE((C - C.transpose()).norm() == Approx(0.0));
        }
        for (auto& σ : σ_list) REQUIRE(σ.norm() == Approx(0.0));
    }
    SECTION("Uniaxial elastic load")
    {
        for (auto& H : H_list) H(2, 2) = 0.001;

        j2plasticity.update_internal_variables(1.0);

        auto[vm_list,
             eff_plastic_list] = variables(InternalVariables::Scalar::VonMisesStress,
                                           InternalVariables::Scalar::EffectivePlasticStrain);

        // Ensure symmetry is correct
        for (auto const& C : material_tangents)
        {
            REQUIRE((C - C.transpose()).norm() == Approx(0.0));
        }
        for (auto& σ : σ_list)
        {
            REQUIRE(σ.norm() != Approx(0.0));

            // Shear components should be close to zero
            REQUIRE(σ(0, 1) == Approx(0.0));
            REQUIRE(σ(0, 2) == Approx(0.0));
            REQUIRE(σ(1, 2) == Approx(0.0));
            REQUIRE(σ(0, 0) > 0.0);
            REQUIRE(σ(1, 1) > 0.0);
            REQUIRE(σ(2, 2) > 0.0);
        }
        for (auto& ε_p_eff : eff_plastic_list) REQUIRE(ε_p_eff == Approx(0.0));
        for (auto& vm : vm_list) REQUIRE(vm < 200.0e6);
    }
    SECTION("Plastic uniaxial elastic load")
    {
        for (auto& H : H_list) H(2, 2) = 0.002;

        j2plasticity.update_internal_variables(1.0);

        auto[vm_list,
             eff_plastic_list] = variables(InternalVariables::Scalar::VonMisesStress,
                                           InternalVariables::Scalar::EffectivePlasticStrain);

        // Ensure symmetry is correct
        for (auto const& C : material_tangents)
        {
            REQUIRE((C - C.transpose()).norm() == Approx(0.0));
        }
        for (auto& σ : σ_list)
        {
            REQUIRE(σ.norm() != Approx(0.0));

            // Shear components should be close to zero
            REQUIRE(σ(0, 1) == Approx(0.0));
            REQUIRE(σ(0, 2) == Approx(0.0));
            REQUIRE(σ(1, 2) == Approx(0.0));
            REQUIRE(σ(0, 0) > 0.0);
            REQUIRE(σ(1, 1) > 0.0);
            REQUIRE(σ(2, 2) > 0.0);
        }
        for (auto& ε_p_eff : eff_plastic_list)
        {
            REQUIRE(ε_p_eff > 0.0);
        }
        for (auto& vm : vm_list)
        {
            // Should experience hardening
            REQUIRE(vm <= 201.0e6);
            REQUIRE(vm > 200.0e6);
        }
    }
}
