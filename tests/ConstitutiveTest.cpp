#define CATCH_CONFIG_MAIN // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"

#include <iostream>
#include <json/json.h>

#include "constitutive/Hyperelastic.hpp"
#include "constitutive/InternalVariables.hpp"

using namespace neon;

std::string json_input_file()
{
    return "{\"Name\": \"steel\", \"ElasticModulus\": 1.0, \"PoissonsRatio\": 2.0}";
}

TEST_CASE("Tensor operations")
{
    SECTION("Uniaxial loading")
    {
        Matrix3 F_old = Matrix3::Identity();
        Matrix3 F = Matrix3::Identity();
        F(0, 0) = 1.01;

        double const Δt = 0.1;

        Matrix3 const Fdot = time_derivative(F, F_old, Δt);

        REQUIRE(Fdot.norm() == Approx(0.1));
        REQUIRE(Fdot(0, 0) == Approx(0.1));

        Matrix3 const L = velocity_gradient(Fdot, F);

        REQUIRE(L(0, 0) == Approx(0.099010));

        Matrix3 const D = rate_of_deformation(L);

        REQUIRE(D(0, 0) == Approx(0.099010));
        REQUIRE((D - L).norm() == Approx(0.0));

        Matrix3 const D_alt = rate_of_deformation(F, F_old, Δt);

        // Check the whole function against the individual functions
        REQUIRE((D_alt - D).norm() == Approx(0.0));
    }
}

TEST_CASE("Neo-Hookean model", "[NeoHooke]")
{
    constexpr auto internal_variable_size = 4;

    InternalVariables variables(internal_variable_size);

    // Add the required variables for an updated Lagrangian formulation
    variables.add(InternalVariables::Tensor::DeformationGradient, InternalVariables::Tensor::Cauchy);
    variables.add(InternalVariables::Scalar::DetF);

    // Create a json reader object from a string
    auto input_data = json_input_file();

    Json::Value material_data;
    Json::Reader reader;

    REQUIRE(reader.parse(input_data.c_str(), material_data));

    NeoHooke neo_hooke(variables, material_data);

    // Get the tensor variables
    auto[F_list, σ_list] =
        variables(InternalVariables::Tensor::DeformationGradient, InternalVariables::Tensor::Cauchy);

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

            REQUIRE((C - C.transpose()).norm() == Approx(0.0));
        }
    }
}
