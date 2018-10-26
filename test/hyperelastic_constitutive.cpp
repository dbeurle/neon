
#include <catch.hpp>

#include "constitutive/internal_variables.hpp"
#include "constitutive/mechanics/solid/affine_microsphere.hpp"
#include "constitutive/mechanics/solid/compressible_neohooke.hpp"
#include "constitutive/mechanics/plane/isotropic_linear_elasticity.hpp"
#include "constitutive/constitutive_model_factory.hpp"

#include "exceptions.hpp"
#include "io/json.hpp"

#include <Eigen/Eigenvalues>

#include <iostream>

std::string json_input_file()
{
    return "{\"name\": \"rubber\", \"elastic_modulus\": 2.0, \"poissons_ratio\": 0.45}";
}

constexpr auto internal_variable_size = 2;
constexpr auto ZERO_MARGIN = 1.0e-5;

using neon::json;
using namespace neon;

TEST_CASE("Neo-Hookean model")
{
    using namespace neon::mechanics::solid;

    auto variables = std::make_shared<internal_variables_t>(internal_variable_size);

    // Add the required variables for an updated Lagrangian formulation
    variables->add(variable::second::deformation_gradient, variable::second::cauchy_stress);
    variables->add(variable::scalar::DetF);

    auto neo_hooke = make_constitutive_model(variables,
                                             json::parse(json_input_file()),
                                             json::parse("{\"constitutive\" : {\"name\": "
                                                         "\"neohooke\"} }"));

    // Get the tensor variables
    auto [F_list, cauchy_stresses] = variables->get(variable::second::deformation_gradient,
                                                    variable::second::cauchy_stress);

    auto& J_list = variables->get(variable::scalar::DetF);

    // Fill with identity matrix
    for (auto& F : F_list) F = neon::matrix3::Identity();
    for (auto& J : J_list) J = 1.0;

    neo_hooke->update_internal_variables(1.0);

    Eigen::EigenSolver<Eigen::MatrixXd> eigen_solver;

    SECTION("Sanity check")
    {
        // Ensure the internal variables are allocated correctly
        REQUIRE(F_list.size() == internal_variable_size);
        REQUIRE(cauchy_stresses.size() == internal_variable_size);
        REQUIRE(J_list.size() == internal_variable_size);

        REQUIRE(neo_hooke->is_finite_deformation());
        REQUIRE(neo_hooke->intrinsic_material().name() == "rubber");
    }
    SECTION("Update of internal variables")
    {
        for (auto& cauchy_stress : cauchy_stresses)
        {
            REQUIRE(cauchy_stress.norm() == Approx(0.0).margin(ZERO_MARGIN));
        }
    }
    SECTION("Check of material tangent")
    {
        // Get the matrix variable
        auto& material_tangents = variables->get(variable::fourth::tangent_operator);

        for (auto const& material_tangent : material_tangents)
        {
            // Check a few of the numerical values
            REQUIRE(material_tangent(0, 0) == Approx(material_tangent(1, 1)));
            REQUIRE(material_tangent(0, 1) == Approx(material_tangent(0, 2)));
            REQUIRE(material_tangent(1, 1) == Approx(material_tangent(2, 2)));
            REQUIRE(material_tangent(3, 3) == Approx(material_tangent(4, 4)));
            REQUIRE(material_tangent(5, 5) == Approx(material_tangent(4, 4)));

            // Ensure symmetry is correct
            REQUIRE((material_tangent - material_tangent.transpose()).norm()
                    == Approx(0.0).margin(ZERO_MARGIN));

            eigen_solver.compute(material_tangent);

            REQUIRE((eigen_solver.eigenvalues().real().array() > 0.0).all());
        }
    }
}
TEST_CASE("microsphere model error test")
{
    using namespace neon::mechanics::solid;

    auto variables = std::make_shared<internal_variables_t>(internal_variable_size);

    SECTION("type not specified")
    {
        REQUIRE_THROWS_AS(make_constitutive_model(variables,
                                                  json::parse("{}"),
                                                  json::parse("{\"constitutive\" : {\"name\": "
                                                              "\"microsphere\", \"Tyasdpe\" "
                                                              ": \"affine\"}}")),
                          std::domain_error);
    }
    SECTION("type not specified correctly")
    {
        // clang-format off
        auto const input{"{\"constitutive\":{\"name\":\"microsphere\",\"type\": \"Afwsfine\"}}"};
        // clang-format on

        REQUIRE_THROWS_AS(make_constitutive_model(variables, json::parse("{}"), json::parse(input)),
                          std::domain_error);
    }
    SECTION("Exception for quadrature scheme")
    {
        // clang-format off
        auto const input{"{\"constitutive\":{\"name\":\"microsphere\",\"type\": \"affine\",\"statistics\":\"gaussian\"}}"};
        // clang-format on

        REQUIRE_THROWS_AS(make_constitutive_model(variables, json::parse("{}"), json::parse(input)),
                          std::domain_error);
    }
}
TEST_CASE("Gaussian affine microsphere")
{
    using namespace neon::mechanics::solid;

    auto variables = std::make_shared<internal_variables_t>(internal_variable_size);

    // Add the required variables for an updated Lagrangian formulation
    variables->add(variable::second::deformation_gradient, variable::second::cauchy_stress);

    variables->add(variable::scalar::DetF);

    auto affine = make_constitutive_model(variables,
                                          json::parse("{\"name\" : \"rubber\", "
                                                      "\"elastic_modulus\" : 10.0e6, "
                                                      "\"poissons_ratio\" : 0.45, "
                                                      "\"segments_per_chain\" : 50}"),
                                          json::parse("{\"constitutive\" : {\"name\": "
                                                      "\"microsphere\", \"type\" "
                                                      ": \"affine\", \"statistics\":\"gaussian\",  "
                                                      "\"quadrature\" : \"BO21\"}}"));

    auto [F_list, cauchy_stresses] = variables->get(variable::second::deformation_gradient,
                                                    variable::second::cauchy_stress);

    auto& J_list = variables->get(variable::scalar::DetF);

    for (auto& J : J_list) J = 1.0;

    auto& material_tangents = variables->get(variable::fourth::tangent_operator);

    Eigen::EigenSolver<Eigen::MatrixXd> eigen_solver;

    SECTION("Sanity checks")
    {
        REQUIRE(affine->is_symmetric());
        REQUIRE(affine->is_finite_deformation());
        REQUIRE(affine->intrinsic_material().name() == "rubber");
    }
    SECTION("affine model under no load")
    {
        // Fill with identity matrix
        for (auto& F : F_list) F = neon::matrix3::Identity();

        affine->update_internal_variables(1.0);

        for (auto const& material_tangent : material_tangents)
        {
            REQUIRE((material_tangent - material_tangent.transpose()).norm()
                    == Approx(0.0).margin(ZERO_MARGIN));

            eigen_solver.compute(material_tangent);
            REQUIRE((eigen_solver.eigenvalues().real().array() > 0.0).all());
        }

        for (auto& cauchy_stress : cauchy_stresses)
        {
            REQUIRE(cauchy_stress.norm() == Approx(0.0).margin(ZERO_MARGIN));
        }
    }
    SECTION("affine model under uniaxial load")
    {
        for (auto& F : F_list)
        {
            F(0, 0) = 1.1;
            F(1, 1) = 1.0 / std::sqrt(1.1);
            F(2, 2) = 1.0 / std::sqrt(1.1);
        }

        affine->update_internal_variables(1.0);

        for (auto const& material_tangent : material_tangents)
        {
            REQUIRE((material_tangent - material_tangent.transpose()).norm()
                    == Approx(0.0).margin(ZERO_MARGIN));

            eigen_solver.compute(material_tangent);
            REQUIRE((eigen_solver.eigenvalues().real().array() > 0.0).all());
        }

        for (auto& cauchy_stress : cauchy_stresses)
        {
            REQUIRE(cauchy_stress.norm() > 0.0);
        }
    }
}
TEST_CASE("Affine microsphere")
{
    using namespace neon::mechanics::solid;

    auto variables = std::make_shared<internal_variables_t>(internal_variable_size);

    // Add the required variables for an updated Lagrangian formulation
    variables->add(variable::second::deformation_gradient, variable::second::cauchy_stress);

    variables->add(variable::scalar::DetF);

    auto affine = make_constitutive_model(variables,
                                          json::parse("{\"name\" : \"rubber\", "
                                                      "\"elastic_modulus\" : 10.0e6, "
                                                      "\"poissons_ratio\" : 0.45, "
                                                      "\"segments_per_chain\" : 50}"),
                                          json::parse("{\"constitutive\" : {\"name\": "
                                                      "\"microsphere\", \"type\" "
                                                      ": \"affine\", \"statistics\":\"langevin\", "
                                                      "\"quadrature\" : \"BO21\"}}"));

    auto [F_list, cauchy_stresses] = variables->get(variable::second::deformation_gradient,
                                                    variable::second::cauchy_stress);

    auto& J_list = variables->get(variable::scalar::DetF);

    for (auto& J : J_list) J = 1.0;

    auto& material_tangents = variables->get(variable::fourth::tangent_operator);

    Eigen::EigenSolver<Eigen::MatrixXd> eigen_solver;

    SECTION("Sanity checks")
    {
        REQUIRE(affine->is_symmetric());
        REQUIRE(affine->is_finite_deformation());
        REQUIRE(affine->intrinsic_material().name() == "rubber");
    }
    SECTION("affine model under no load")
    {
        // Fill with identity matrix
        for (auto& F : F_list) F = neon::matrix3::Identity();

        affine->update_internal_variables(1.0);

        for (auto const& material_tangent : material_tangents)
        {
            REQUIRE((material_tangent - material_tangent.transpose()).norm()
                    == Approx(0.0).margin(ZERO_MARGIN));

            eigen_solver.compute(material_tangent);
            REQUIRE((eigen_solver.eigenvalues().real().array() > 0.0).all());
        }

        for (auto& cauchy_stress : cauchy_stresses)
        {
            REQUIRE(cauchy_stress.norm() == Approx(0.0).margin(ZERO_MARGIN));
        }
    }
    SECTION("affine model under uniaxial load")
    {
        for (auto& F : F_list)
        {
            F(0, 0) = 1.1;
            F(1, 1) = 1.0 / std::sqrt(1.1);
            F(2, 2) = 1.0 / std::sqrt(1.1);
        }

        affine->update_internal_variables(1.0);

        for (auto const& material_tangent : material_tangents)
        {
            REQUIRE((material_tangent - material_tangent.transpose()).norm()
                    == Approx(0.0).margin(ZERO_MARGIN));

            eigen_solver.compute(material_tangent);
            REQUIRE((eigen_solver.eigenvalues().real().array() > 0.0).all());
        }

        for (auto& cauchy_stress : cauchy_stresses)
        {
            REQUIRE(cauchy_stress.norm() > 0.0);
        }
    }
}
TEST_CASE("Nonaffine microsphere")
{
    using namespace neon::mechanics::solid;

    auto variables = std::make_shared<internal_variables_t>(internal_variable_size);

    // Add the required variables for an updated Lagrangian formulation
    variables->add(variable::second::deformation_gradient, variable::second::cauchy_stress);
    variables->add(variable::scalar::DetF);

    auto affine = make_constitutive_model(variables,
                                          json::parse("{\"name\" : \"rubber\", "
                                                      "\"elastic_modulus\" : 10.0e6, "
                                                      "\"poissons_ratio\" : 0.45, "
                                                      "\"nonaffine_stretch_parameter\":1.0, "
                                                      "\"segments_per_chain\" : 50}"),
                                          json::parse("{\"constitutive\" : {\"name\": "
                                                      "\"microsphere\", \"type\" "
                                                      ": \"nonaffine\", \"quadrature\" : "
                                                      "\"BO21\"}}"));

    // Get the tensor variables
    auto [F_list, cauchy_stresses] = variables->get(variable::second::deformation_gradient,
                                                    variable::second::cauchy_stress);

    auto& J_list = variables->get(variable::scalar::DetF);

    for (auto& J : J_list) J = 1.0;

    auto& material_tangents = variables->get(variable::fourth::tangent_operator);

    Eigen::EigenSolver<Eigen::MatrixXd> eigen_solver;

    SECTION("exception")
    {
        REQUIRE_THROWS_AS(make_constitutive_model(variables,
                                                  json::parse("{\"name\" : \"rubber\", "
                                                              "\"elastic_modulus\" : 10.0e6, "
                                                              "\"poissons_ratio\" : 0.45, "
                                                              "\"nonaffine_Strearameter\":1."
                                                              "0, "
                                                              "\"segments_per_chain\" : 50}"),
                                                  json::parse("{\"constitutive\" : "
                                                              "{\"name\": "
                                                              "\"microsphere\", \"type\" "
                                                              ": \"nonaffine\", \"quadrature\" "
                                                              ": "
                                                              "\"BO21\"}}")),
                          std::domain_error);
    }
    SECTION("Sanity test")
    {
        REQUIRE(affine->is_symmetric());
        REQUIRE(affine->is_finite_deformation());
        REQUIRE(affine->intrinsic_material().name() == "rubber");
    }
    SECTION("Nonaffine model under no load")
    {
        // Fill with identity matrix
        for (auto& F : F_list) F = neon::matrix3::Identity();

        affine->update_internal_variables(1.0);

        for (auto const& material_tangent : material_tangents)
        {
            // Ensure symmetry is correct
            REQUIRE((material_tangent - material_tangent.transpose()).norm()
                    == Approx(0.0).margin(ZERO_MARGIN));

            eigen_solver.compute(material_tangent);
            REQUIRE((eigen_solver.eigenvalues().real().array() > 0.0).all());
        }

        for (auto& cauchy_stress : cauchy_stresses)
        {
            REQUIRE(cauchy_stress.norm() == Approx(0.0).margin(ZERO_MARGIN));
        }
    }
    SECTION("Nonaffine model under uniaxial load")
    {
        for (auto& F : F_list)
        {
            F(0, 0) = 1.1;
            F(1, 1) = 1.0 / std::sqrt(1.1);
            F(2, 2) = 1.0 / std::sqrt(1.1);
        }

        affine->update_internal_variables(1.0);

        for (auto const& material_tangent : material_tangents)
        {
            // Ensure symmetry is correct
            REQUIRE((material_tangent - material_tangent.transpose()).norm()
                    == Approx(0.0).margin(ZERO_MARGIN));

            eigen_solver.compute(material_tangent);
            REQUIRE((eigen_solver.eigenvalues().real().array() > 0.0).all());
        }

        for (auto& cauchy_stress : cauchy_stresses)
        {
            REQUIRE(cauchy_stress.norm() > 0.0);
        }
    }
}
