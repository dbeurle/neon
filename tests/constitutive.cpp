
#include <catch.hpp>

#include "constitutive/internal_variables.hpp"
#include "constitutive/mechanical/solid/affine_microsphere.hpp"
#include "constitutive/mechanical/solid/compressible_neohooke.hpp"
#include "constitutive/mechanical/solid/small_strain_J2_plasticity.hpp"
#include "constitutive/mechanical/detail/J2_plasticity.hpp"
#include "constitutive/mechanical/plane/isotropic_linear_elasticity.hpp"
#include "constitutive/thermal/isotropic_diffusion.hpp"
#include "constitutive/constitutive_model_factory.hpp"
#include "material/isotropic_elastic_plastic.hpp"
#include "exceptions.hpp"
#include "io/json.hpp"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

std::string json_input_file()
{
    return "{\"Name\": \"rubber\", \"ElasticModulus\": 2.0, \"PoissonsRatio\": 0.45}";
}

std::string json_thermal_diffusion()
{
    return "{\"Name\": \"steel\", \"Conductivity\": 386.0, \"Density\": 7800.0, \"SpecificHeat\": "
           "390.0}";
}

constexpr auto internal_variable_size = 2;
constexpr auto ZERO_MARGIN = 1.0e-5;

using neon::json;
using namespace neon;

TEST_CASE("No constitutive model error test")
{
    using namespace neon::mechanical::solid;

    auto variables = std::make_shared<internal_variables_t>(internal_variable_size);

    REQUIRE_THROWS_AS(make_constitutive_model(variables, json::parse("{}"), json::parse("{}")),
                      std::domain_error);
}
TEST_CASE("Constitutive model no name error test")
{
    using namespace neon::mechanical::solid;

    auto variables = std::make_shared<internal_variables_t>(internal_variable_size);

    REQUIRE_THROWS_AS(make_constitutive_model(variables,
                                              json::parse("{}"),
                                              json::parse("{\"ConstitutiveModel\" : {}}")),
                      std::domain_error);
}
TEST_CASE("Constitutive model invalid name error test")
{
    using namespace neon::mechanical::solid;

    auto variables = std::make_shared<internal_variables_t>(internal_variable_size);

    REQUIRE_THROWS_AS(make_constitutive_model(variables,
                                              json::parse("{}"),
                                              json::parse("{\"ConstitutiveModel\" : "
                                                          "{\"Name\":\"PurpleMonkey\"}}")),
                      std::domain_error);
}
TEST_CASE("Neo-Hookean model")
{
    using namespace neon::mechanical::solid;

    auto variables = std::make_shared<internal_variables_t>(internal_variable_size);

    // Add the required variables for an updated Lagrangian formulation
    variables->add(variable::second::deformation_gradient, variable::second::cauchy_stress);
    variables->add(variable::scalar::DetF);

    auto neo_hooke = make_constitutive_model(variables,
                                             json::parse(json_input_file()),
                                             json::parse("{\"ConstitutiveModel\" : {\"Name\": "
                                                         "\"NeoHooke\"} }"));

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
TEST_CASE("Microsphere model error test")
{
    using namespace neon::mechanical::solid;

    auto variables = std::make_shared<internal_variables_t>(internal_variable_size);

    SECTION("Type not specified")
    {
        REQUIRE_THROWS_AS(make_constitutive_model(variables,
                                                  json::parse("{}"),
                                                  json::parse("{\"ConstitutiveModel\" : {\"Name\": "
                                                              "\"Microsphere\", \"Tyasdpe\" "
                                                              ": \"Affine\"}}")),
                          std::domain_error);
    }
    SECTION("Type not specified correctly")
    {
        // clang-format off
        auto const input{"{\"ConstitutiveModel\":{\"Name\":\"Microsphere\",\"Type\": \"Afwsfine\"}}"};
        // clang-format on

        REQUIRE_THROWS_AS(make_constitutive_model(variables, json::parse("{}"), json::parse(input)),
                          std::domain_error);
    }
    SECTION("Exception for quadrature scheme")
    {
        // clang-format off
        auto const input{"{\"ConstitutiveModel\":{\"Name\":\"Microsphere\",\"Type\": \"Affine\",\"Statistics\":\"Gaussian\"}}"};
        // clang-format on

        REQUIRE_THROWS_AS(make_constitutive_model(variables, json::parse("{}"), json::parse(input)),
                          std::domain_error);
    }
}
TEST_CASE("Gaussian affine microsphere model")
{
    using namespace neon::mechanical::solid;

    auto variables = std::make_shared<internal_variables_t>(internal_variable_size);

    // Add the required variables for an updated Lagrangian formulation
    variables->add(variable::second::deformation_gradient, variable::second::cauchy_stress);

    variables->add(variable::scalar::DetF);

    auto affine = make_constitutive_model(variables,
                                          json::parse("{\"Name\" : \"rubber\", "
                                                      "\"ElasticModulus\" : 10.0e6, "
                                                      "\"PoissonsRatio\" : 0.45, "
                                                      "\"SegmentsPerChain\" : 50}"),
                                          json::parse("{\"ConstitutiveModel\" : {\"Name\": "
                                                      "\"Microsphere\", \"Type\" "
                                                      ": \"Affine\", \"Statistics\":\"Gaussian\",  "
                                                      "\"Quadrature\" : "
                                                      "\"BO21\"}}"));

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
    SECTION("Affine model under no load")
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
    SECTION("Affine model under uniaxial load")
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
TEST_CASE("Gaussian affine microsphere model with ageing")
{
    using namespace neon::mechanical::solid;

    auto variables = std::make_shared<internal_variables_t>(1);

    // Add the required variables for an updated Lagrangian formulation
    variables->add(variable::second::deformation_gradient, variable::second::cauchy_stress);

    variables->add(variable::scalar::DetF);

    auto const material_data{"{\"Name\" : \"rubber\","
                             "\"ShearModulus\" : 2.0e6,"
                             "\"BulkModulus\" : 100e6,"
                             "\"SegmentsPerChain\" : 50,"
                             "\"ScissionProbability\" : 1.0e-5,"
                             "\"RecombinationProbability\" : 1.0e-5}"};

    auto const constitutive_data{"{\"ConstitutiveModel\" : {\"Name\": \"Microsphere\","
                                 "\"Type\":\"Affine\","
                                 "\"Statistics\":\"Gaussian\","
                                 "\"Quadrature\":\"BO21\","
                                 "\"Ageing\":\"BAND\"}}"};

    auto affine = make_constitutive_model(variables,
                                          json::parse(material_data),
                                          json::parse(constitutive_data));

    auto [F_list, cauchy_stresses] = variables->get(variable::second::deformation_gradient,
                                                    variable::second::cauchy_stress);

    auto& J_list = variables->get(variable::scalar::DetF);
    std::fill(begin(J_list), end(J_list), 1.0);

    auto& material_tangents = variables->get(variable::fourth::tangent_operator);

    Eigen::EigenSolver<Eigen::MatrixXd> eigen_solver;

    SECTION("Setup checks")
    {
        REQUIRE(affine->is_symmetric());
        REQUIRE(affine->is_finite_deformation());
        REQUIRE(affine->intrinsic_material().name() == "rubber");

        REQUIRE(variables->has(variable::scalar::active_shear_modulus));
        REQUIRE(variables->has(variable::scalar::inactive_shear_modulus));
        REQUIRE(variables->has(variable::scalar::active_segments));
        REQUIRE(variables->has(variable::scalar::inactive_segments));
        REQUIRE(variables->has(variable::scalar::reduction_factor));

        REQUIRE(variables->has(variable::vector::first_ageing_moduli));
        REQUIRE(variables->has(variable::vector::second_ageing_moduli));
        REQUIRE(variables->has(variable::vector::third_ageing_moduli));

        REQUIRE(variables->has(variable::vector::first_previous));
        REQUIRE(variables->has(variable::vector::second_previous));
        REQUIRE(variables->has(variable::vector::third_previous));

        for (auto segment : variables->get(variable::scalar::active_segments))
        {
            REQUIRE(segment == Approx(50.0));
        }
        for (auto shear_modulus : variables->get(variable::scalar::active_shear_modulus))
        {
            REQUIRE(shear_modulus == Approx(2.0e6));
        }
        for (auto& value : variables->get(variable::vector::first_ageing_moduli))
        {
            REQUIRE(value.size() == 21);
        }
        for (auto& value : variables->get(variable::vector::first_previous))
        {
            REQUIRE(value.size() == 21);
        }
    }
    SECTION("no load")
    {
        std::fill(begin(F_list), end(F_list), neon::matrix3::Identity());

        affine->update_internal_variables(1.0);

        // Check the network parameters
        auto [active_segments,
              inactive_segments,
              active_shear_moduli,
              inactive_shear_moduli,
              reductions] = variables->get(variable::scalar::active_segments,
                                           variable::scalar::inactive_segments,
                                           variable::scalar::active_shear_modulus,
                                           variable::scalar::inactive_shear_modulus,
                                           variable::scalar::reduction_factor);

        for (auto active_segment : active_segments)
        {
            REQUIRE(active_segment > 49.0);
            REQUIRE(active_segment < 50.0);
        }
        for (auto inactive_segment : inactive_segments)
        {
            REQUIRE(inactive_segment > 20.0);
            REQUIRE(inactive_segment < 25.0);
        }
        for (auto shear_modulus : active_shear_moduli)
        {
            REQUIRE(shear_modulus < 2.01e6);
            REQUIRE(shear_modulus > 2.0e6);
        }
        for (auto shear_modulus : inactive_shear_moduli)
        {
            REQUIRE(shear_modulus < 2000.0);
            REQUIRE(shear_modulus > 1000.0);
        }
        for (auto reduction : reductions)
        {
            REQUIRE(reduction > 0.99);
            REQUIRE(reduction < 1.0);
        }

        for (auto const& material_tangent : material_tangents)
        {
            REQUIRE((material_tangent - material_tangent.transpose()).norm()
                    == Approx(0.0).margin(ZERO_MARGIN));

            eigen_solver.compute(material_tangent);
            REQUIRE((eigen_solver.eigenvalues().real().array() > 0.0).all());
        }
    }
    SECTION("uniaxial load")
    {
        for (auto& F : F_list)
        {
            F(0, 0) = 1.0;
            F(1, 1) = 1.0;
            F(2, 2) = 1.0;
        }

        affine->update_internal_variables(1.0);

        for (auto& cauchy_stress : cauchy_stresses)
        {
            std::cout << cauchy_stress << "\n\n";
            REQUIRE(cauchy_stress.norm() > 0.0);
        }

        std::cout << "setting to stretch 1.1\n";

        for (auto& F : F_list)
        {
            F(0, 0) = 1.1;
            F(1, 1) = 1.0 / std::sqrt(1.1);
            F(2, 2) = 1.0 / std::sqrt(1.1);
        }

        affine->update_internal_variables(1.0);

        for (auto& cauchy_stress : cauchy_stresses)
        {
            std::cout << cauchy_stress << "\n\n";
            REQUIRE(cauchy_stress.norm() > 0.0);
        }

        affine->update_internal_variables(1.0);

        for (auto& cauchy_stress : cauchy_stresses)
        {
            std::cout << cauchy_stress << "\n\n";
            REQUIRE(cauchy_stress.norm() > 0.0);
        }

        std::cout << "setting to unity\n";

        for (auto& F : F_list)
        {
            F(0, 0) = 1.0;
            F(1, 1) = 1.0;
            F(2, 2) = 1.0;
        }

        affine->update_internal_variables(1.0);

        for (auto& cauchy_stress : cauchy_stresses)
        {
            std::cout << cauchy_stress << "\n\n";
            REQUIRE(cauchy_stress.norm() >= 0.0);
        }

        std::cout << "finished the analysis\n";

        for (auto const& material_tangent : material_tangents)
        {
            REQUIRE((material_tangent - material_tangent.transpose()).norm()
                    == Approx(0.0).margin(ZERO_MARGIN));

            eigen_solver.compute(material_tangent);
            REQUIRE((eigen_solver.eigenvalues().real().array() > 0.0).all());
        }
    }
}
TEST_CASE("Gaussian affine microsphere model with crosslinking only")
{
    using namespace neon::mechanical::solid;

    std::cout << "Constant cross-linking stress check\n";

    auto variables = std::make_shared<internal_variables_t>(1);

    // Add the required variables for an updated Lagrangian formulation
    variables->add(variable::second::deformation_gradient, variable::second::cauchy_stress);

    variables->add(variable::scalar::DetF);

    auto const material_data{"{\"Name\" : \"rubber\","
                             "\"ShearModulus\" : 2.0e6,"
                             "\"BulkModulus\" : 100e6,"
                             "\"SegmentsPerChain\" : 50,"
                             "\"ScissionProbability\" : 0.0,"
                             "\"RecombinationProbability\" : 1.0e-5}"};

    auto const constitutive_data{"{\"ConstitutiveModel\" : {\"Name\": \"Microsphere\","
                                 "\"Type\":\"Affine\","
                                 "\"Statistics\":\"Gaussian\","
                                 "\"Quadrature\":\"BO21\","
                                 "\"Ageing\":\"BAND\"}}"};

    auto affine = make_constitutive_model(variables,
                                          json::parse(material_data),
                                          json::parse(constitutive_data));

    auto [F_list, cauchy_stresses] = variables->get(variable::second::deformation_gradient,
                                                    variable::second::cauchy_stress);

    auto& J_list = variables->get(variable::scalar::DetF);
    std::fill(begin(J_list), end(J_list), 1.0);

    auto& material_tangents = variables->get(variable::fourth::tangent_operator);

    Eigen::EigenSolver<Eigen::MatrixXd> eigen_solver;

    SECTION("no load")
    {
        std::fill(begin(F_list), end(F_list), neon::matrix3::Identity());

        affine->update_internal_variables(1.0);

        // Check the network parameters
        auto [active_segments,
              inactive_segments,
              active_shear_moduli,
              inactive_shear_moduli,
              reductions] = variables->get(variable::scalar::active_segments,
                                           variable::scalar::inactive_segments,
                                           variable::scalar::active_shear_modulus,
                                           variable::scalar::inactive_shear_modulus,
                                           variable::scalar::reduction_factor);

        for (auto active_segment : active_segments)
        {
            REQUIRE(active_segment > 49.0);
            REQUIRE(active_segment < 50.0);
        }
        for (auto inactive_segment : inactive_segments)
        {
            REQUIRE(inactive_segment == Approx(0.0).margin(0.0));
        }
        for (auto shear_modulus : active_shear_moduli)
        {
            REQUIRE(shear_modulus < 2.01e6);
            REQUIRE(shear_modulus > 2.0e6);
        }
        for (auto shear_modulus : inactive_shear_moduli)
        {
            REQUIRE(shear_modulus == Approx(0.0).margin(0.0));
        }
        for (auto reduction : reductions)
        {
            REQUIRE(reduction <= 1.0);
        }

        for (auto const& material_tangent : material_tangents)
        {
            REQUIRE((material_tangent - material_tangent.transpose()).norm()
                    == Approx(0.0).margin(ZERO_MARGIN));

            eigen_solver.compute(material_tangent);
            REQUIRE((eigen_solver.eigenvalues().real().array() > 0.0).all());
        }
    }
    SECTION("constant load")
    {
        // Check the network parameters
        // auto [active_segments,
        //       inactive_segments,
        //       active_shear_moduli,
        //       inactive_shear_moduli,
        //       reductions] = variables->get(variable::scalar::active_segments,
        //                                    variable::scalar::inactive_segments,
        //                                    variable::scalar::active_shear_modulus,
        //                                    variable::scalar::inactive_shear_modulus,
        //                                    variable::scalar::reduction_factor);

        std::fill(begin(F_list), end(F_list), neon::matrix3::Identity());

        affine->update_internal_variables(1.0);

        for (auto& F : F_list)
        {
            F(0, 0) = 1.1;
            F(1, 1) = 1.0 / std::sqrt(1.1);
            F(2, 2) = 1.0 / std::sqrt(1.1);
        }

        affine->update_internal_variables(1.0);

        for (auto const& cauchy_stress : cauchy_stresses)
        {
            std::cout << "Step 1: cauchy_stress\n" << cauchy_stress << "\n\n";
        }

        affine->update_internal_variables(1.0);

        for (auto const& cauchy_stress : cauchy_stresses)
        {
            std::cout << "Step 2: cauchy_stress\n" << cauchy_stress << "\n\n";
        }

        affine->update_internal_variables(1.0);

        for (auto const& cauchy_stress : cauchy_stresses)
        {
            std::cout << "Step 3: cauchy_stress\n" << cauchy_stress << "\n\n";
        }

        affine->update_internal_variables(1.0);

        for (auto const& cauchy_stress : cauchy_stresses)
        {
            std::cout << "Step 3: cauchy_stress\n" << cauchy_stress << "\n\n";
        }

        affine->update_internal_variables(1.0);

        for (auto const& cauchy_stress : cauchy_stresses)
        {
            std::cout << "Step 3: cauchy_stress\n" << cauchy_stress << "\n\n";
        }

        for (auto const& material_tangent : material_tangents)
        {
            std::cout << "tangent_matrix\n" << material_tangent << "\n";

            REQUIRE((material_tangent - material_tangent.transpose()).norm()
                    == Approx(0.0).margin(ZERO_MARGIN));

            eigen_solver.compute(material_tangent);
            REQUIRE((eigen_solver.eigenvalues().real().array() > 0.0).all());
        }
    }
}
TEST_CASE("Affine microsphere model")
{
    using namespace neon::mechanical::solid;

    auto variables = std::make_shared<internal_variables_t>(internal_variable_size);

    // Add the required variables for an updated Lagrangian formulation
    variables->add(variable::second::deformation_gradient, variable::second::cauchy_stress);

    variables->add(variable::scalar::DetF);

    auto affine = make_constitutive_model(variables,
                                          json::parse("{\"Name\" : \"rubber\", "
                                                      "\"ElasticModulus\" : 10.0e6, "
                                                      "\"PoissonsRatio\" : 0.45, "
                                                      "\"SegmentsPerChain\" : 50}"),
                                          json::parse("{\"ConstitutiveModel\" : {\"Name\": "
                                                      "\"Microsphere\", \"Type\" "
                                                      ": \"Affine\", \"Statistics\":\"Langevin\", "
                                                      "\"Quadrature\" : \"BO21\"}}"));

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
    SECTION("Affine model under no load")
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
    SECTION("Affine model under uniaxial load")
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
TEST_CASE("NonAffine microsphere model")
{
    using namespace neon::mechanical::solid;

    auto variables = std::make_shared<internal_variables_t>(internal_variable_size);

    // Add the required variables for an updated Lagrangian formulation
    variables->add(variable::second::deformation_gradient, variable::second::cauchy_stress);
    variables->add(variable::scalar::DetF);

    auto affine = make_constitutive_model(variables,
                                          json::parse("{\"Name\" : \"rubber\", "
                                                      "\"ElasticModulus\" : 10.0e6, "
                                                      "\"PoissonsRatio\" : 0.45, "
                                                      "\"NonAffineStretchParameter\":1.0, "
                                                      "\"SegmentsPerChain\" : 50}"),
                                          json::parse("{\"ConstitutiveModel\" : {\"Name\": "
                                                      "\"Microsphere\", \"Type\" "
                                                      ": \"NonAffine\", \"Quadrature\" : "
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
                                                  json::parse("{\"Name\" : \"rubber\", "
                                                              "\"ElasticModulus\" : 10.0e6, "
                                                              "\"PoissonsRatio\" : 0.45, "
                                                              "\"NonAffineStrearameter\":1."
                                                              "0, "
                                                              "\"SegmentsPerChain\" : 50}"),
                                                  json::parse("{\"ConstitutiveModel\" : "
                                                              "{\"Name\": "
                                                              "\"Microsphere\", \"Type\" "
                                                              ": \"NonAffine\", \"Quadrature\" "
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
    SECTION("NonAffine model under no load")
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
    SECTION("NonAffine model under uniaxial load")
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
TEST_CASE("Plane stress linear elasticity factory error")
{
    using namespace neon::mechanical::plane;

    auto variables = std::make_shared<internal_variables_t>(internal_variable_size);

    REQUIRE_THROWS_AS(make_constitutive_model(variables,
                                              json::parse("{\"Badkey\" : \"donkey\"}"),
                                              json::parse("{\"ConstitutiveModel\" : {\"Name\" : "
                                                          "\"PlaneStrain\"}}")),
                      std::domain_error);
}
TEST_CASE("Plane stress elasticity model")
{
    using namespace neon::mechanical::plane;

    auto variables = std::make_shared<internal_variables_t>(internal_variable_size);

    // Add the required variables for an updated Lagrangian formulation
    variables->add(variable::second::displacement_gradient,
                   variable::second::cauchy_stress,
                   variable::scalar::DetF);

    auto elastic_model = make_constitutive_model(variables,
                                                 json::parse("{\"Name\": \"steel\", "
                                                             "\"ElasticModulus\": 200.0e3, "
                                                             "\"PoissonsRatio\": 0.3}"),
                                                 json::parse("{\"ConstitutiveModel\" : {\"Name\" : "
                                                             "\"PlaneStress\"}}"));

    // Get the tensor variables
    auto [displacement_gradients,
          cauchy_stresses] = variables->get(variable::second::displacement_gradient,
                                            variable::second::cauchy_stress);

    auto& J_list = variables->get(variable::scalar::DetF);

    auto& material_tangents = variables->get(variable::fourth::tangent_operator);

    for (auto& H : displacement_gradients) H = internal_variables_t::second_tensor_type::Zero();
    for (auto& J : J_list) J = 1.0;

    Eigen::EigenSolver<Eigen::MatrixXd> eigen_solver;

    SECTION("Sanity checks")
    {
        REQUIRE(elastic_model->is_symmetric());
        REQUIRE(elastic_model->is_finite_deformation() == false);
        REQUIRE(elastic_model->intrinsic_material().name() == "steel");

        REQUIRE(variables->has(variable::scalar::von_mises_stress));
        REQUIRE(variables->has(variable::second::cauchy_stress));
        REQUIRE(variables->has(variable::second::linearised_strain));
        REQUIRE(variables->has(variable::fourth::tangent_operator));
    }
    SECTION("No load")
    {
        elastic_model->update_internal_variables(1.0);

        // Ensure symmetry is correct
        for (auto const& material_tangent : material_tangents)
        {
            REQUIRE(material_tangent(0, 0) == Approx(219780.21978022));
            REQUIRE(material_tangent(1, 1) == Approx(219780.21978022));
            REQUIRE(material_tangent(2, 2) == Approx(76923.0769230769));
            REQUIRE(material_tangent(0, 1) == Approx(65934.065934066));
            REQUIRE(material_tangent(1, 0) == Approx(65934.065934066));

            REQUIRE(material_tangent(0, 2) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(material_tangent(2, 0) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(material_tangent(1, 2) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(material_tangent(2, 1) == Approx(0.0).margin(ZERO_MARGIN));

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
    SECTION("Uniaxial elastic load")
    {
        for (auto& H : displacement_gradients) H(1, 1) = 0.001;

        elastic_model->update_internal_variables(1.0);

        auto [von_mises_stresses,
              accumulated_plastic_strains] = variables->get(variable::scalar::von_mises_stress,
                                                            variable::scalar::effective_plastic_strain);

        for (auto const& material_tangent : material_tangents)
        {
            REQUIRE(material_tangent(0, 0) == Approx(219780.21978022));
            REQUIRE(material_tangent(1, 1) == Approx(219780.21978022));
            REQUIRE(material_tangent(2, 2) == Approx(76923.0769230769));
            REQUIRE(material_tangent(0, 1) == Approx(65934.065934066));
            REQUIRE(material_tangent(1, 0) == Approx(65934.065934066));

            REQUIRE(material_tangent(0, 2) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(material_tangent(2, 0) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(material_tangent(1, 2) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(material_tangent(2, 1) == Approx(0.0).margin(ZERO_MARGIN));

            // Ensure symmetry is correct
            REQUIRE((material_tangent - material_tangent.transpose()).norm()
                    == Approx(0.0).margin(ZERO_MARGIN));

            eigen_solver.compute(material_tangent);
            REQUIRE((eigen_solver.eigenvalues().real().array() > 0.0).all());
        }

        for (auto const& cauchy_stress : cauchy_stresses)
        {
            REQUIRE(cauchy_stress.norm() != Approx(0.0).margin(ZERO_MARGIN));

            // Shear components should be close to zero
            REQUIRE(cauchy_stress(0, 1) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(cauchy_stress(0, 0) > 0.0);
            REQUIRE(cauchy_stress(1, 1) > 0.0);
        }

        for (auto& von_mises_stress : von_mises_stresses)
        {
            REQUIRE(von_mises_stress < 200.0e6);
        }
    }
}
TEST_CASE("Plane strain elasticity model")
{
    using namespace neon::mechanical::plane;

    auto variables = std::make_shared<internal_variables_t>(internal_variable_size);

    // Add the required variables for an updated Lagrangian formulation
    variables->add(variable::second::displacement_gradient,
                   variable::second::cauchy_stress,
                   variable::scalar::DetF);

    auto elastic_model = make_constitutive_model(variables,
                                                 json::parse("{\"Name\": \"steel\", "
                                                             "\"ElasticModulus\": 200.0e3, "
                                                             "\"PoissonsRatio\": 0.3}"),
                                                 json::parse("{\"ConstitutiveModel\" : {\"Name\" : "
                                                             "\"PlaneStrain\"}}"));

    // Get the tensor variables
    auto [displacement_gradients,
          cauchy_stresses] = variables->get(variable::second::displacement_gradient,
                                            variable::second::cauchy_stress);

    auto& J_list = variables->get(variable::scalar::DetF);

    auto& material_tangents = variables->get(variable::fourth::tangent_operator);

    for (auto& H : displacement_gradients) H = internal_variables_t::second_tensor_type::Zero();
    for (auto& J : J_list) J = 1.0;

    Eigen::EigenSolver<Eigen::MatrixXd> eigen_solver;

    SECTION("Sanity checks")
    {
        REQUIRE(elastic_model->is_symmetric());
        REQUIRE(elastic_model->is_finite_deformation() == false);
        REQUIRE(elastic_model->intrinsic_material().name() == "steel");

        REQUIRE(variables->has(variable::scalar::von_mises_stress));
        REQUIRE(variables->has(variable::second::cauchy_stress));
        REQUIRE(variables->has(variable::second::linearised_strain));
        REQUIRE(variables->has(variable::fourth::tangent_operator));
    }
    SECTION("No load")
    {
        elastic_model->update_internal_variables(1.0);

        // Ensure symmetry is correct
        for (auto const& material_tangent : material_tangents)
        {
            REQUIRE(material_tangent(0, 0) == Approx(269230.769230769));
            REQUIRE(material_tangent(1, 1) == Approx(269230.769230769));
            REQUIRE(material_tangent(2, 2) == Approx(76923.0769230769));
            REQUIRE(material_tangent(0, 1) == Approx(115384.615384615));

            REQUIRE(material_tangent(0, 2) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(material_tangent(2, 0) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(material_tangent(1, 2) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(material_tangent(2, 1) == Approx(0.0).margin(ZERO_MARGIN));

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
    SECTION("Uniaxial elastic load")
    {
        for (auto& H : displacement_gradients) H(1, 1) = 0.001;

        elastic_model->update_internal_variables(1.0);

        auto [von_mises_stresses,
              accumulated_plastic_strains] = variables->get(variable::scalar::von_mises_stress,
                                                            variable::scalar::effective_plastic_strain);

        for (auto const& material_tangent : material_tangents)
        {
            REQUIRE(material_tangent(0, 0) == Approx(269230.769230769));
            REQUIRE(material_tangent(1, 1) == Approx(269230.769230769));
            REQUIRE(material_tangent(2, 2) == Approx(76923.0769230769));
            REQUIRE(material_tangent(0, 1) == Approx(115384.615384615));

            REQUIRE(material_tangent(0, 2) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(material_tangent(2, 0) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(material_tangent(1, 2) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(material_tangent(2, 1) == Approx(0.0).margin(ZERO_MARGIN));

            // Ensure symmetry is correct
            REQUIRE((material_tangent - material_tangent.transpose()).norm()
                    == Approx(0.0).margin(ZERO_MARGIN));

            eigen_solver.compute(material_tangent);
            REQUIRE((eigen_solver.eigenvalues().real().array() > 0.0).all());
        }

        for (auto const& cauchy_stress : cauchy_stresses)
        {
            REQUIRE(cauchy_stress.norm() != Approx(0.0).margin(ZERO_MARGIN));

            // Shear components should be close to zero
            REQUIRE(cauchy_stress(0, 1) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(cauchy_stress(0, 0) > 0.0);
            REQUIRE(cauchy_stress(1, 1) > 0.0);
        }

        for (auto& von_mises_stress : von_mises_stresses)
        {
            REQUIRE(von_mises_stress < 200.0e6);
        }
    }
}
TEST_CASE("Solid mechanics elasticity model")
{
    using namespace neon::mechanical::solid;

    auto variables = std::make_shared<internal_variables_t>(internal_variable_size);

    // Add the required variables for an updated Lagrangian formulation
    variables->add(variable::second::displacement_gradient,
                   variable::second::cauchy_stress,
                   variable::scalar::DetF);

    auto elastic_model = make_constitutive_model(variables,
                                                 json::parse("{\"Name\": \"steel\", "
                                                             "\"ElasticModulus\": 200.0e3, "
                                                             "\"PoissonsRatio\": 0.3}"),
                                                 json::parse("{\"ConstitutiveModel\" : {\"Name\" : "
                                                             "\"IsotropicLinearElasticity\"}}"));

    // Get the tensor variables
    auto [displacement_gradients,
          cauchy_stresses] = variables->get(variable::second::displacement_gradient,
                                            variable::second::cauchy_stress);

    auto& J_list = variables->get(variable::scalar::DetF);

    auto& material_tangents = variables->get(variable::fourth::tangent_operator);

    for (auto& H : displacement_gradients) H = internal_variables_t::second_tensor_type::Zero();
    for (auto& J : J_list) J = 1.0;

    Eigen::EigenSolver<Eigen::MatrixXd> eigen_solver;

    SECTION("Sanity checks")
    {
        REQUIRE(elastic_model->is_symmetric());
        REQUIRE(elastic_model->is_finite_deformation() == false);
        REQUIRE(elastic_model->intrinsic_material().name() == "steel");

        REQUIRE(variables->has(variable::scalar::von_mises_stress));
        REQUIRE(variables->has(variable::second::cauchy_stress));
        REQUIRE(variables->has(variable::second::linearised_strain));
        REQUIRE(variables->has(variable::fourth::tangent_operator));
    }
    SECTION("No load")
    {
        elastic_model->update_internal_variables(1.0);

        // Ensure symmetry is correct
        for (auto const& material_tangent : material_tangents)
        {
            REQUIRE(material_tangent(0, 0) == Approx(269230.769230769));
            REQUIRE(material_tangent(1, 1) == Approx(269230.769230769));
            REQUIRE(material_tangent(2, 2) == Approx(269230.769230769));

            REQUIRE(material_tangent(3, 3) == Approx(76923.0769230769));
            REQUIRE(material_tangent(4, 4) == Approx(76923.0769230769));
            REQUIRE(material_tangent(5, 5) == Approx(76923.0769230769));

            REQUIRE(material_tangent(0, 1) == Approx(115384.615384615));
            REQUIRE(material_tangent(0, 2) == Approx(115384.615384615));
            REQUIRE(material_tangent(1, 2) == Approx(115384.615384615));

            REQUIRE(material_tangent(0, 3) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(material_tangent(0, 4) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(material_tangent(0, 5) == Approx(0.0).margin(ZERO_MARGIN));

            REQUIRE(material_tangent(1, 3) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(material_tangent(1, 4) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(material_tangent(1, 5) == Approx(0.0).margin(ZERO_MARGIN));

            REQUIRE(material_tangent(2, 3) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(material_tangent(2, 4) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(material_tangent(2, 5) == Approx(0.0).margin(ZERO_MARGIN));

            REQUIRE(material_tangent(3, 4) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(material_tangent(3, 5) == Approx(0.0).margin(ZERO_MARGIN));

            REQUIRE(material_tangent(4, 5) == Approx(0.0).margin(ZERO_MARGIN));

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
    SECTION("Uniaxial elastic load")
    {
        for (auto& H : displacement_gradients) H(1, 1) = 0.001;

        elastic_model->update_internal_variables(1.0);

        auto [von_mises_stresses,
              accumulated_plastic_strains] = variables->get(variable::scalar::von_mises_stress,
                                                            variable::scalar::effective_plastic_strain);

        for (auto const& material_tangent : material_tangents)
        {
            REQUIRE(material_tangent(0, 0) == Approx(269230.769230769));
            REQUIRE(material_tangent(1, 1) == Approx(269230.769230769));
            REQUIRE(material_tangent(2, 2) == Approx(269230.769230769));

            REQUIRE(material_tangent(3, 3) == Approx(76923.0769230769));
            REQUIRE(material_tangent(4, 4) == Approx(76923.0769230769));
            REQUIRE(material_tangent(5, 5) == Approx(76923.0769230769));

            REQUIRE(material_tangent(0, 1) == Approx(115384.615384615));
            REQUIRE(material_tangent(0, 2) == Approx(115384.615384615));
            REQUIRE(material_tangent(1, 2) == Approx(115384.615384615));

            REQUIRE(material_tangent(0, 3) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(material_tangent(0, 4) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(material_tangent(0, 5) == Approx(0.0).margin(ZERO_MARGIN));

            REQUIRE(material_tangent(1, 3) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(material_tangent(1, 4) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(material_tangent(1, 5) == Approx(0.0).margin(ZERO_MARGIN));

            REQUIRE(material_tangent(2, 3) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(material_tangent(2, 4) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(material_tangent(2, 5) == Approx(0.0).margin(ZERO_MARGIN));

            REQUIRE(material_tangent(3, 4) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(material_tangent(3, 5) == Approx(0.0).margin(ZERO_MARGIN));

            REQUIRE(material_tangent(4, 5) == Approx(0.0).margin(ZERO_MARGIN));

            // Ensure symmetry is correct
            REQUIRE((material_tangent - material_tangent.transpose()).norm()
                    == Approx(0.0).margin(ZERO_MARGIN));

            eigen_solver.compute(material_tangent);
            REQUIRE((eigen_solver.eigenvalues().real().array() > 0.0).all());
        }

        for (auto const& cauchy_stress : cauchy_stresses)
        {
            REQUIRE(cauchy_stress.norm() != Approx(0.0).margin(ZERO_MARGIN));

            // Shear components should be close to zero
            REQUIRE(cauchy_stress(0, 1) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(cauchy_stress(0, 0) > 0.0);
            REQUIRE(cauchy_stress(1, 1) > 0.0);
        }

        for (auto& von_mises_stress : von_mises_stresses)
        {
            REQUIRE(von_mises_stress < 200.0e6);
        }
    }
}
TEST_CASE("Solid mechanics J2 plasticity model factory errors")
{
    using namespace neon::mechanical::solid;

    auto variables = std::make_shared<internal_variables_t>(internal_variable_size);

    SECTION("FiniteStrain not specified")
    {
        REQUIRE_THROWS_AS(make_constitutive_model(variables,
                                                  json::parse("{}"),
                                                  json::parse("{\"ConstitutiveModel\" : {\"Name\" "
                                                              ": "
                                                              "\"J2Plasticity\"}}")),
                          std::domain_error);
    }
    SECTION("Damage model not specified correctly")
    {
        REQUIRE_THROWS_AS(make_constitutive_model(variables,
                                                  json::parse("{}"),
                                                  json::parse("{\"ConstitutiveModel\" : {\"Name\" "
                                                              ": "
                                                              "\"J2Plasticity\", "
                                                              "\"Damage\" : \"IsotropicChaboche\", "
                                                              "\"FiniteStrain\" : true}}")),
                          std::domain_error);
    }
}
TEST_CASE("Solid mechanics J2 plasticity model")
{
    using namespace neon::mechanical::solid;

    auto variables = std::make_shared<internal_variables_t>(internal_variable_size);

    // Add the required variables for an updated Lagrangian formulation
    variables->add(variable::second::displacement_gradient, variable::second::cauchy_stress);
    variables->add(variable::scalar::DetF);

    auto const material_data = json::parse("{\"Name\": \"steel\", "
                                           "\"ElasticModulus\": 200.0e9, "
                                           "\"PoissonsRatio\": 0.3, "
                                           "\"YieldStress\": 200.0e6, "
                                           "\"IsotropicHardeningModulus\": 400.0e6}");

    auto const constitutive_data = json::parse("{\"ConstitutiveModel\" : "
                                               "{\"Name\" : \"J2Plasticity\", "
                                               "\"FiniteStrain\" : false}}");

    auto small_strain_J2_plasticity = make_constitutive_model(variables,
                                                              material_data,
                                                              constitutive_data);

    // Get the tensor variables
    auto [displacement_gradients,
          cauchy_stresses] = variables->get(variable::second::displacement_gradient,
                                            variable::second::cauchy_stress);

    auto& J_list = variables->get(variable::scalar::DetF);

    auto& material_tangents = variables->get(variable::fourth::tangent_operator);

    for (auto& H : displacement_gradients) H = neon::matrix3::Zero();
    for (auto& J : J_list) J = 1.0;

    Eigen::EigenSolver<Eigen::MatrixXd> eigen_solver;

    SECTION("Sanity checks")
    {
        REQUIRE(small_strain_J2_plasticity->is_symmetric());
        REQUIRE(small_strain_J2_plasticity->is_finite_deformation() == false);
        REQUIRE(small_strain_J2_plasticity->intrinsic_material().name() == "steel");

        REQUIRE(variables->has(variable::scalar::von_mises_stress));
        REQUIRE(variables->has(variable::scalar::effective_plastic_strain));
        REQUIRE(variables->has(variable::second::linearised_strain));
        REQUIRE(variables->has(variable::second::linearised_plastic_strain));
        REQUIRE(variables->has(variable::fourth::tangent_operator));
    }
    SECTION("Helper function")
    {
        neon::isotropic_elastic_plastic material{material_data};

        REQUIRE(neon::mechanical::evaluate_J2_yield_function(material, 100.0e6, 0.0) <= 0.0);
        REQUIRE(neon::mechanical::evaluate_J2_yield_function(material, 300.0e6, 0.0) > 0.0);
    }
    SECTION("No load")
    {
        small_strain_J2_plasticity->update_internal_variables(1.0);

        // Ensure symmetry is correct
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
    SECTION("Uniaxial elastic load")
    {
        for (auto& H : displacement_gradients) H(2, 2) = 0.001;

        small_strain_J2_plasticity->update_internal_variables(1.0);

        auto [von_mises_stresses,
              accumulated_plastic_strains] = variables->get(variable::scalar::von_mises_stress,
                                                            variable::scalar::effective_plastic_strain);

        for (auto const& material_tangent : material_tangents)
        {
            // Ensure symmetry is correct
            REQUIRE((material_tangent - material_tangent.transpose()).norm()
                    == Approx(0.0).margin(ZERO_MARGIN));

            eigen_solver.compute(material_tangent);
            REQUIRE((eigen_solver.eigenvalues().real().array() > 0.0).all());
        }

        for (auto const& cauchy_stress : cauchy_stresses)
        {
            REQUIRE(cauchy_stress.norm() != Approx(0.0).margin(ZERO_MARGIN));

            // Shear components should be close to zero
            REQUIRE(cauchy_stress(0, 1) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(cauchy_stress(0, 2) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(cauchy_stress(1, 2) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(cauchy_stress(0, 0) > 0.0);
            REQUIRE(cauchy_stress(1, 1) > 0.0);
            REQUIRE(cauchy_stress(2, 2) > 0.0);
        }

        for (auto& accumulated_plastic_strain : accumulated_plastic_strains)
        {
            REQUIRE(accumulated_plastic_strain == Approx(0.0).margin(ZERO_MARGIN));
        }

        for (auto& von_mises_stress : von_mises_stresses)
        {
            REQUIRE(von_mises_stress < 200.0e6);
        }
    }
    SECTION("Plastic uniaxial elastic load")
    {
        for (auto& H : displacement_gradients) H(2, 2) = 0.003;

        small_strain_J2_plasticity->update_internal_variables(1.0);

        auto [von_mises_stresses,
              accumulated_plastic_strains] = variables->get(variable::scalar::von_mises_stress,
                                                            variable::scalar::effective_plastic_strain);

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
            REQUIRE(cauchy_stress.norm() != Approx(0.0).margin(ZERO_MARGIN));

            // Shear components should be close to zero
            REQUIRE(cauchy_stress(0, 1) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(cauchy_stress(0, 2) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(cauchy_stress(1, 2) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(cauchy_stress(0, 0) > 0.0);
            REQUIRE(cauchy_stress(1, 1) > 0.0);
            REQUIRE(cauchy_stress(2, 2) > 0.0);
        }

        for (auto& accumulated_plastic_strain : accumulated_plastic_strains)
        {
            REQUIRE(accumulated_plastic_strain > 0.0);
        }

        for (auto& von_mises_stress : von_mises_stresses)
        {
            // Should experience hardening
            REQUIRE(von_mises_stress <= 201.0e6);
            REQUIRE(von_mises_stress > 200.0e6);
        }
    }
}
TEST_CASE("Solid mechanics J2 plasticity damage model")
{
    using namespace neon::mechanical::solid;

    auto variables = std::make_shared<internal_variables_t>(internal_variable_size);

    variables->add(variable::second::displacement_gradient, variable::second::cauchy_stress);
    variables->add(variable::scalar::DetF);

    auto const material_input{"{\"Name\":\"steel\","
                              "\"ElasticModulus\":134.0e3,"
                              "\"PoissonsRatio\":0.3,"
                              "\"YieldStress\":85,"
                              "\"KinematicHardeningModulus\": 5500,"
                              "\"SofteningMultiplier\" : 250,"
                              "\"PlasticityViscousExponent\" : 2.5,"
                              "\"PlasticityViscousDenominator\" : 1220,"
                              "\"DamageViscousExponent\" : 2,"
                              "\"DamageViscousDenominator\" : 0.6"
                              "}"};

    auto const constitutive_input{"{\"ConstitutiveModel\" : "
                                  "{\"Name\" : \"J2Plasticity\","
                                  "\"Damage\" : \"IsotropicChaboche\", "
                                  "\"FiniteStrain\" : false}}"};

    auto small_strain_J2_plasticity_damage = make_constitutive_model(variables,
                                                                     json::parse(material_input),
                                                                     json::parse(constitutive_input));

    // Get the tensor variables
    auto [displacement_gradients,
          cauchy_stresses] = variables->get(variable::second::displacement_gradient,
                                            variable::second::cauchy_stress);

    auto [J_list, damage_list] = variables->get(variable::scalar::DetF, variable::scalar::damage);

    auto& material_tangents = variables->get(variable::fourth::tangent_operator);

    for (auto& H : displacement_gradients) H = neon::matrix3::Zero();
    for (auto& J : J_list) J = 1.0;

    Eigen::EigenSolver<Eigen::MatrixXd> eigen_solver;

    SECTION("Sanity checks")
    {
        REQUIRE(small_strain_J2_plasticity_damage->is_finite_deformation() == false);
        REQUIRE(small_strain_J2_plasticity_damage->is_symmetric() == false);

        REQUIRE(small_strain_J2_plasticity_damage->intrinsic_material().name() == "steel");

        REQUIRE(variables->has(variable::scalar::von_mises_stress));
        REQUIRE(variables->has(variable::scalar::effective_plastic_strain));
        REQUIRE(variables->has(variable::second::linearised_strain));
        REQUIRE(variables->has(variable::second::linearised_plastic_strain));
        REQUIRE(variables->has(variable::fourth::tangent_operator));
        REQUIRE(variables->has(variable::scalar::damage));
        REQUIRE(variables->has(variable::scalar::energy_release_rate));
        REQUIRE(variables->has(variable::second::kinematic_hardening));
        REQUIRE(variables->has(variable::second::back_stress));
    }
    SECTION("Uniaxial elastic load")
    {
        for (auto& H : displacement_gradients) H(2, 2) = 0.0008;

        small_strain_J2_plasticity_damage->update_internal_variables(1.0);

        auto [von_mises_stresses,
              accumulated_plastic_strains] = variables->get(variable::scalar::von_mises_stress,
                                                            variable::scalar::effective_plastic_strain);

        for (auto const& material_tangent : material_tangents)
        {
            // Ensure symmetry is correct when the loading is elastic
            REQUIRE((material_tangent - material_tangent.transpose()).norm()
                    == Approx(0.0).margin(ZERO_MARGIN));

            eigen_solver.compute(material_tangent);
            REQUIRE((eigen_solver.eigenvalues().real().array() > 0.0).all());
        }

        for (auto& cauchy_stress : cauchy_stresses)
        {
            REQUIRE(cauchy_stress.norm() != Approx(0.0).margin(ZERO_MARGIN));

            // Shear components should be close to zero
            REQUIRE(cauchy_stress(0, 1) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(cauchy_stress(0, 2) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(cauchy_stress(1, 2) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(cauchy_stress(0, 0) == Approx(cauchy_stress(1, 1)));
            REQUIRE(cauchy_stress(1, 1) > 0.0);
            REQUIRE(cauchy_stress(2, 2) > 0.0);
        }
        for (auto& accumulated_plastic_strain : accumulated_plastic_strains)
        {
            REQUIRE(accumulated_plastic_strain == Approx(0.0).margin(ZERO_MARGIN));
        }
        for (auto& von_mises_stress : von_mises_stresses)
        {
            REQUIRE(von_mises_stress > 80);
            REQUIRE(von_mises_stress < 85);
        }
    }
    SECTION("Plastic uniaxial load")
    {
        for (auto& H : displacement_gradients) H(2, 2) = 0.0009; // 0.000825 for comparison with 1D

        small_strain_J2_plasticity_damage->update_internal_variables(
            0.01); // time here is real (not pseudo time)

        auto [von_mises_stresses,
              accumulated_plastic_strains] = variables->get(variable::scalar::von_mises_stress,
                                                            variable::scalar::effective_plastic_strain);

        for (auto const& cauchy_stress : cauchy_stresses)
        {
            REQUIRE(cauchy_stress.norm() != Approx(0.0).margin(ZERO_MARGIN));

            // Shear components should be close to zero
            REQUIRE(cauchy_stress(0, 1) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(cauchy_stress(0, 2) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(cauchy_stress(1, 2) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(cauchy_stress(0, 0) > 0.0);
            REQUIRE(cauchy_stress(1, 1) > 0.0);
            REQUIRE(cauchy_stress(2, 2) > 0.0);
        }
        for (auto& accumulated_plastic_strain : accumulated_plastic_strains)
        {
            REQUIRE(accumulated_plastic_strain > 0.0);
        }
        for (auto& damage_var : damage_list)
        {
            REQUIRE(damage_var > 0.0);
        }
        for (auto& von_mises_stress : von_mises_stresses)
        {
            // Should experience hardening
            REQUIRE(von_mises_stress <= 100);
            REQUIRE(von_mises_stress > 85);
        }
        for (auto const& material_tangent : material_tangents)
        {
            eigen_solver.compute(material_tangent);
            REQUIRE((eigen_solver.eigenvalues().real().array() > 0.0).all());
        }
    }
}
TEST_CASE("Thermal isotropic model")
{
    //
    using namespace neon::diffusion;

    auto variables = std::make_shared<internal_variables_t>(internal_variable_size);

    auto thermal = make_constitutive_model(variables,
                                           json::parse(json_thermal_diffusion()),
                                           json::parse("{\"ConstitutiveModel\" : {\"Name\": "
                                                       "\"IsotropicDiffusion\"} }"));

    thermal->update_internal_variables(1.0);

    SECTION("Sanity checks")
    {
        REQUIRE(thermal->is_symmetric());
        REQUIRE(thermal->is_finite_deformation() == false);
        REQUIRE(thermal->intrinsic_material().name() == "steel");

        REQUIRE(variables->has(variable::second::conductivity));
    }
    Eigen::EigenSolver<Eigen::MatrixXd> eigen_solver;
}

// TEST_CASE("Finite J2 plasticity model", "[finite_strain_J2_plasticity]")
// {
//     using namespace neon::mechanical::solid;
//
//     // Create a json reader object from a string
//     std::string
//         input_data = "{\"Name\": \"steel\", \"ElasticModulus\": 200.0e9, \"PoissonsRatio\": 0.3,
//         "
//                      "\"YieldStress\": 200.0e6, \"IsotropicHardeningModulus\": 400.0e6}";
//
//     std::string simulation_input = "{\"ConstitutiveModel\" : {\"Name\" : \"small_strain_J2_plasticity\", "
//                                    "\"FiniteStrain\":true}}";
//
//     Json::Value material_data, simulation_data;
//
//     Json::Reader reader;
//
//     REQUIRE(reader.parse(input_data.c_str(), material_data));
//     REQUIRE(reader.parse(simulation_input.c_str(), simulation_data));
//
//     internal_variables variables(1);
//
//     // Add the required variables for an updated Lagrangian formulation
//     variables->add(variable::second::deformation_gradient,
//     variable::second::cauchy_stress); variables->add(variable::scalar::DetF);
//
//     auto small_strain_J2_plasticity = make_constitutive_model(variables, material_data, simulation_data);
//
//     // Get the tensor variables
//     auto[F_list, cauchy_stresses] =
//     variables->get(variable::second::deformation_gradient,
//                                               variable::second::cauchy_stress);
//
//     auto& J_list = variables->get(variable::scalar::DetF);
//
//     auto& material_tangents = variables->get(variable::fourth::tangent_operator);
//
//     for (auto& F : F_list) F = neon::matrix3::Identity();
//     for (auto& J : J_list) J = 1.0;
//
//     variables.commit();
//
//     SECTION("Sanity checks")
//     {
//         REQUIRE(small_strain_J2_plasticity->is_finite_deformation() == true);
//         REQUIRE(small_strain_J2_plasticity->intrinsic_material().name() == "steel");
//
//         REQUIRE(variables->has(variable::scalar::von_mises_stress));
//         REQUIRE(variables->has(variable::scalar::effective_plastic_strain));
//         REQUIRE(variables->has(variable::second::hencky_strain_elastic));
//         REQUIRE(variables->has(variable::fourth::tangent_operator));
//     }
//     SECTION("Initial material tangent symmetry")
//     {
//         for (auto const& material_tangent : material_tangents)
//         {
//             REQUIRE((material_tangent - material_tangent.transpose()).norm()
//                     == Approx(0.0).margin(ZERO_MARGIN));
//         }
//     }
//     SECTION("No load")
//     {
//         small_strain_J2_plasticity->update_internal_variables(1.0);
//
//         for (auto const& material_tangent : material_tangents)
//         {
//             REQUIRE((material_tangent - material_tangent.transpose()).norm()
//                     == Approx(0.0).margin(ZERO_MARGIN));
//         }
//         for (auto& cauchy_stress : cauchy_stresses)
//         {
//             REQUIRE(cauchy_stress.norm() == Approx(0.0).margin(ZERO_MARGIN));
//         }
//     }
//     // SECTION("Uniaxial elastic load")
//     // {
//     //     for (auto& F : F_list) F(2, 2) = 1.001;
//     //
//     //     small_strain_J2_plasticity->update_internal_variables(1.0);
//     //
//     //     for (auto const& material_tangent : material_tangents)
//     //     {
//     //         REQUIRE((material_tangent - material_tangent.transpose()).norm() ==
//     //         Approx(0.0).margin(ZERO_MARGIN));
//     //     }
//     //     for (auto& cauchy_stress : cauchy_stresses)
//     //     {
//     //         REQUIRE((cauchy_stress - cauchy_stress.transpose()).norm() ==
//     Approx(0.0).margin(ZERO_MARGIN));
//     //
//     //         // Shear components should be close to zero
//     //         REQUIRE(cauchy_stress(0, 1) == Approx(0.0).margin(ZERO_MARGIN));
//     //         REQUIRE(cauchy_stress(0, 2) == Approx(0.0).margin(ZERO_MARGIN));
//     //         REQUIRE(cauchy_stress(1, 2) == Approx(0.0).margin(ZERO_MARGIN));
//     //         REQUIRE(cauchy_stress(0, 0) > 0.0);
//     //         REQUIRE(cauchy_stress(1, 1) > 0.0);
//     //         REQUIRE(cauchy_stress(2, 2) > 0.0);
//     //     }
//     //     for (auto& von_mises_stress :
//     variables->get(variable::scalar::von_mises_stress))
//     //     {
//     //         REQUIRE(von_mises_stress < 200.0e6);
//     //     }
//     // }
//     // SECTION("Plastic uniaxial elastic load")
//     // {
//     //     for (auto& F : F_list) F(2, 2) = 0.003;
//     //
//     //     small_strain_J2_plasticity->update_internal_variables(1.0);
//     //
//     //     auto[von_mises_stresses, accumulated_plastic_strains] =
//     //     variables->get(variable::scalar::von_mises_stress,
//     // variable::scalar::effective_plastic_strain);
//     //
//     //     // Ensure symmetry is correct
//     //     for (auto const& material_tangent : material_tangents)
//     //     {
//     //         REQUIRE((material_tangent - material_tangent.transpose()).norm() ==
//     //         Approx(0.0).margin(ZERO_MARGIN));
//     //     }
//     //     for (auto& cauchy_stress : cauchy_stresses)
//     //     {
//     //         REQUIRE(cauchy_stress.norm() != Approx(0.0).margin(ZERO_MARGIN));
//     //
//     //         // Shear components should be close to zero
//     //         REQUIRE(cauchy_stress(0, 1) == Approx(0.0).margin(ZERO_MARGIN));
//     //         REQUIRE(cauchy_stress(0, 2) == Approx(0.0).margin(ZERO_MARGIN));
//     //         REQUIRE(cauchy_stress(1, 2) == Approx(0.0).margin(ZERO_MARGIN));
//     //         REQUIRE(cauchy_stress(0, 0) > 0.0);
//     //         REQUIRE(cauchy_stress(1, 1) > 0.0);
//     //         REQUIRE(cauchy_stress(2, 2) > 0.0);
//     //     }
//     //     for (auto& accumulated_plastic_strain : accumulated_plastic_strains)
//     //     {
//     //         REQUIRE(accumulated_plastic_strain > 0.0);
//     //     }
//     //     for (auto& von_mises_stress : von_mises_stresses)
//     //     {
//     //         // Should experience hardening
//     //         REQUIRE(von_mises_stress <= 201.0e6);
//     //         REQUIRE(von_mises_stress > 200.0e6);
//     //     }
//     // }
// }
