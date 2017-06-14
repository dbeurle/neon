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

TEST_CASE("Neo-Hookean model", "[NeoHooke]")
{
    constexpr auto internal_variable_size = 4;

    InternalVariables variables(internal_variable_size);

    // Create a json reader object from a string
    auto input_data = json_input_file();

    Json::Value material_data;
    Json::Reader reader;

    REQUIRE(reader.parse(input_data.c_str(), material_data));

    NeoHooke neo_hooke(variables, material_data);

    // Get the tensor variables
    auto[deformation_gradients, kirchhoff_stresses] =
        variables(InternalVariables::Tensor::DeformationGradient,
                  InternalVariables::Tensor::Kirchhoff);

    // Ensure the internal variables are allocated correctly
    REQUIRE(deformation_gradients.size() == internal_variable_size);
    REQUIRE(kirchhoff_stresses.size() == internal_variable_size);

    SECTION("Update of internal variables")
    {
        // Fill with identity matrix
        for (auto& F : deformation_gradients) F = Matrix3::Identity();

        neo_hooke.update_internal_variables();
        for (auto& τ : kirchhoff_stresses)
        {
            REQUIRE(τ.norm() == Approx(0.0));
        }
    }

    SECTION("Check of continuum tangent")
    {
        // Get the matrix variable
        auto& material_tangents = variables(InternalVariables::Matrix::MaterialTangent);

        neo_hooke.update_continuum_tangent();
        for (auto& C : material_tangents)
        {
            // Check a few of the numerical values
            REQUIRE(C(0, 0) == Approx(0.111111));
            REQUIRE(C(0, 1) == Approx(-0.222222));
            REQUIRE(C(0, 2) == Approx(-0.222222));
            REQUIRE(C(2, 2) == Approx(0.111111));
            REQUIRE(C(3, 3) == Approx(0.333333));
            REQUIRE(C(4, 4) == Approx(0.333333));
            REQUIRE(C(5, 5) == Approx(0.333333));
        }
    }
}
