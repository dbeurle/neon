
#define CATCH_CONFIG_MAIN

#include "catch.hpp"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <json/json.h>

#include "constitutive/InternalVariables.hpp"

#include "constitutive/mechanical/solid/AffineMicrosphere.hpp"
#include "constitutive/mechanical/solid/J2Plasticity.hpp"
#include "constitutive/mechanical/solid/NeoHooke.hpp"

#include "constitutive/mechanical/plane/IsotropicLinearElasticity.hpp"

#include "constitutive/thermal/IsotropicDiffusion.hpp"

#include "constitutive/ConstitutiveModelFactory.hpp"

#include "Exceptions.hpp"

Json::CharReaderBuilder reader;
JSONCPP_STRING input_errors;

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

TEST_CASE("No constitutive model error test")
{
    using namespace neon::mechanical::solid;

    auto variables = std::make_shared<InternalVariables>(internal_variable_size);

    // Create a json reader object from a string
    std::string input_data = "{}";
    std::string simulation_input = "{}";

    Json::Value material_data, simulation_data;

    std::istringstream input_data_stream(input_data), simulation_input_stream(simulation_input);

    REQUIRE(Json::parseFromStream(reader, input_data_stream, &material_data, &input_errors));
    REQUIRE(Json::parseFromStream(reader, simulation_input_stream, &simulation_data, &input_errors));

    REQUIRE_THROWS_AS(make_constitutive_model(variables, material_data, simulation_data),
                      std::runtime_error);
}
TEST_CASE("Constitutive model no name error test")
{
    using namespace neon::mechanical::solid;

    auto variables = std::make_shared<InternalVariables>(internal_variable_size);

    // Create a json reader object from a string
    std::string input_data = "{}";
    std::string simulation_input = "{\"ConstitutiveModel\" : {}}";

    Json::Value material_data, simulation_data;

    std::istringstream input_data_stream(input_data), simulation_input_stream(simulation_input);

    REQUIRE(Json::parseFromStream(reader, input_data_stream, &material_data, &input_errors));
    REQUIRE(Json::parseFromStream(reader, simulation_input_stream, &simulation_data, &input_errors));

    REQUIRE_THROWS_AS(make_constitutive_model(variables, material_data, simulation_data),
                      std::runtime_error);
}
TEST_CASE("Constitutive model invalid name error test")
{
    using namespace neon::mechanical::solid;

    auto variables = std::make_shared<InternalVariables>(internal_variable_size);

    // Create a json reader object from a string
    std::string input_data = "{}";
    std::string simulation_input = "{\"ConstitutiveModel\" : {\"Name\":\"PurpleMonkey\"}}";

    Json::Value material_data, simulation_data;

    std::istringstream input_data_stream(input_data), simulation_input_stream(simulation_input);

    REQUIRE(Json::parseFromStream(reader, input_data_stream, &material_data, &input_errors));
    REQUIRE(Json::parseFromStream(reader, simulation_input_stream, &simulation_data, &input_errors));

    REQUIRE_THROWS_AS(make_constitutive_model(variables, material_data, simulation_data),
                      std::runtime_error);
}
TEST_CASE("Neo-Hookean model")
{
    using namespace neon::mechanical::solid;

    auto variables = std::make_shared<InternalVariables>(internal_variable_size);

    // Add the required variables for an updated Lagrangian formulation
    variables->add(InternalVariables::Tensor::DeformationGradient, InternalVariables::Tensor::Cauchy);
    variables->add(InternalVariables::Scalar::DetF);

    // Create a json reader object from a string
    auto input_data = json_input_file();

    std::string simulation_input = "{\"ConstitutiveModel\" : {\"Name\": \"NeoHooke\"} }";

    Json::Value material_data, simulation_data;

    std::istringstream input_data_stream(input_data), simulation_input_stream(simulation_input);

    REQUIRE(Json::parseFromStream(reader, input_data_stream, &material_data, &input_errors));
    REQUIRE(Json::parseFromStream(reader, simulation_input_stream, &simulation_data, &input_errors));

    auto neo_hooke = make_constitutive_model(variables, material_data, simulation_data);

    // Get the tensor variables
    auto[F_list, cauchy_stresses] = variables->fetch(InternalVariables::Tensor::DeformationGradient,
                                                     InternalVariables::Tensor::Cauchy);

    auto& J_list = variables->fetch(InternalVariables::Scalar::DetF);

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
        auto& material_tangents = variables->fetch(InternalVariables::rank4::tangent_operator);

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

    auto variables = std::make_shared<InternalVariables>(internal_variable_size);

    // Create a json reader object from a string
    std::string input_data = "{}";

    std::string simulation_input = "{\"ConstitutiveModel\" : {\"Name\": \"Microsphere\", \"Type\" "
                                   ": \"Afwsfine\"}}";

    Json::Value material_data, simulation_data;

    std::istringstream input_data_stream(input_data), simulation_input_stream(simulation_input);

    REQUIRE(Json::parseFromStream(reader, input_data_stream, &material_data, &input_errors));
    REQUIRE(Json::parseFromStream(reader, simulation_input_stream, &simulation_data, &input_errors));

    REQUIRE_THROWS_AS(make_constitutive_model(variables, material_data, simulation_data),
                      std::runtime_error);
}
TEST_CASE("Affine microsphere model", )
{
    using namespace neon::mechanical::solid;

    auto variables = std::make_shared<InternalVariables>(internal_variable_size);

    // Add the required variables for an updated Lagrangian formulation
    variables->add(InternalVariables::Tensor::DeformationGradient, InternalVariables::Tensor::Cauchy);

    variables->add(InternalVariables::Scalar::DetF);

    // Create a json reader object from a string
    std::string input_data = "{\"Name\" : \"rubber\", "
                             "\"ElasticModulus\" : 10.0e6, "
                             "\"PoissonsRatio\" : 0.45, "
                             "\"SegmentsPerChain\" : 50}";

    std::string simulation_input = "{\"ConstitutiveModel\" : {\"Name\": \"Microsphere\", \"Type\" "
                                   ": \"Affine\", \"Quadrature\" : \"BO21\"}}";

    Json::Value material_data, simulation_data;

    std::istringstream input_data_stream(input_data), simulation_input_stream(simulation_input);

    REQUIRE(Json::parseFromStream(reader, input_data_stream, &material_data, &input_errors));
    REQUIRE(Json::parseFromStream(reader, simulation_input_stream, &simulation_data, &input_errors));

    auto affine = make_constitutive_model(variables, material_data, simulation_data);

    auto[F_list, cauchy_stresses] = variables->fetch(InternalVariables::Tensor::DeformationGradient,
                                                     InternalVariables::Tensor::Cauchy);

    auto& J_list = variables->fetch(InternalVariables::Scalar::DetF);

    for (auto& J : J_list) J = 1.0;

    auto& material_tangents = variables->fetch(InternalVariables::rank4::tangent_operator);

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

    auto variables = std::make_shared<InternalVariables>(internal_variable_size);

    // Add the required variables for an updated Lagrangian formulation
    variables->add(InternalVariables::Tensor::DeformationGradient, InternalVariables::Tensor::Cauchy);

    variables->add(InternalVariables::Scalar::DetF);

    // Create a json reader object from a string
    std::string input_data = "{\"Name\" : \"rubber\", "
                             "\"ElasticModulus\" : 10.0e6, "
                             "\"PoissonsRatio\" : 0.45, "
                             "\"NonAffineStretchParameter\":1.0, "
                             "\"SegmentsPerChain\" : 50}";

    std::string simulation_input = "{\"ConstitutiveModel\" : {\"Name\": \"Microsphere\", \"Type\" "
                                   ": \"NonAffine\", \"Quadrature\" : \"BO21\"}}";

    Json::Value material_data, simulation_data;

    std::istringstream input_data_stream(input_data), simulation_input_stream(simulation_input);

    REQUIRE(Json::parseFromStream(reader, input_data_stream, &material_data, &input_errors));
    REQUIRE(Json::parseFromStream(reader, simulation_input_stream, &simulation_data, &input_errors));

    auto affine = make_constitutive_model(variables, material_data, simulation_data);

    // Get the tensor variables
    auto[F_list, cauchy_stresses] = variables->fetch(InternalVariables::Tensor::DeformationGradient,
                                                     InternalVariables::Tensor::Cauchy);

    auto& J_list = variables->fetch(InternalVariables::Scalar::DetF);

    for (auto& J : J_list) J = 1.0;

    auto& material_tangents = variables->fetch(InternalVariables::rank4::tangent_operator);

    Eigen::EigenSolver<Eigen::MatrixXd> eigen_solver;

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

    // Create a json reader object from a string
    std::string input_data = "{\"Badkey\" : \"donkey\"}";
    std::string mesh_input = "{\"ConstitutiveModel\" : {\"Name\" : \"PlaneStrain\"}}";

    Json::Value material_data, mesh_data;

    std::istringstream input_data_stream(input_data), mesh_input_stream(mesh_input);

    REQUIRE(Json::parseFromStream(reader, input_data_stream, &material_data, &input_errors));
    REQUIRE(Json::parseFromStream(reader, mesh_input_stream, &mesh_data, &input_errors));

    auto variables = std::make_shared<InternalVariables>(internal_variable_size);

    REQUIRE_THROWS_AS(make_constitutive_model(variables, material_data, mesh_data),
                      neon::MaterialPropertyException);
}
TEST_CASE("Plane stress elasticity model")
{
    using namespace neon::mechanical::plane;

    // Create a json reader object from a string
    std::string input_data = "{\"Name\": \"steel\", \"ElasticModulus\": 200.0e3, "
                             "\"PoissonsRatio\": 0.3}";

    std::string mesh_input = "{\"ConstitutiveModel\" : {\"Name\" : \"PlaneStress\"}}";

    Json::Value material_data, mesh_data;

    std::istringstream input_data_stream(input_data), mesh_input_stream(mesh_input);

    REQUIRE(Json::parseFromStream(reader, input_data_stream, &material_data, &input_errors));
    REQUIRE(Json::parseFromStream(reader, mesh_input_stream, &mesh_data, &input_errors));

    auto variables = std::make_shared<InternalVariables>(internal_variable_size);

    // Add the required variables for an updated Lagrangian formulation
    variables->add(InternalVariables::Tensor::DisplacementGradient,
                   InternalVariables::Tensor::Cauchy,
                   InternalVariables::Scalar::DetF);

    auto elastic_model = make_constitutive_model(variables, material_data, mesh_data);

    // Get the tensor variables
    auto[displacement_gradients,
         cauchy_stresses] = variables->fetch(InternalVariables::Tensor::DisplacementGradient,
                                             InternalVariables::Tensor::Cauchy);

    auto& J_list = variables->fetch(InternalVariables::Scalar::DetF);

    auto& material_tangents = variables->fetch(InternalVariables::rank4::tangent_operator);

    for (auto& H : displacement_gradients) H = InternalVariables::rank2tensor_type::Zero();
    for (auto& J : J_list) J = 1.0;

    Eigen::EigenSolver<Eigen::MatrixXd> eigen_solver;

    SECTION("Sanity checks")
    {
        REQUIRE(elastic_model->is_symmetric());
        REQUIRE(elastic_model->is_finite_deformation() == false);
        REQUIRE(elastic_model->intrinsic_material().name() == "steel");

        REQUIRE(variables->has(InternalVariables::Scalar::VonMisesStress));
        REQUIRE(variables->has(InternalVariables::Tensor::Cauchy));
        REQUIRE(variables->has(InternalVariables::Tensor::LinearisedStrain));
        REQUIRE(variables->has(InternalVariables::rank4::tangent_operator));
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

        auto[von_mises_stresses,
             accumulated_plastic_strains] = variables
                                                ->fetch(InternalVariables::Scalar::VonMisesStress,
                                                        InternalVariables::Scalar::EffectivePlasticStrain);

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

    // Create a json reader object from a string
    std::string input_data = "{\"Name\": \"steel\", \"ElasticModulus\": 200.0e3, "
                             "\"PoissonsRatio\": 0.3}";

    std::string mesh_input = "{\"ConstitutiveModel\" : {\"Name\" : \"PlaneStrain\"}}";

    Json::Value material_data, mesh_data;

    std::istringstream input_data_stream(input_data), mesh_input_stream(mesh_input);

    REQUIRE(Json::parseFromStream(reader, input_data_stream, &material_data, &input_errors));
    REQUIRE(Json::parseFromStream(reader, mesh_input_stream, &mesh_data, &input_errors));

    auto variables = std::make_shared<InternalVariables>(internal_variable_size);

    // Add the required variables for an updated Lagrangian formulation
    variables->add(InternalVariables::Tensor::DisplacementGradient,
                   InternalVariables::Tensor::Cauchy,
                   InternalVariables::Scalar::DetF);

    auto elastic_model = make_constitutive_model(variables, material_data, mesh_data);

    // Get the tensor variables
    auto[displacement_gradients,
         cauchy_stresses] = variables->fetch(InternalVariables::Tensor::DisplacementGradient,
                                             InternalVariables::Tensor::Cauchy);

    auto& J_list = variables->fetch(InternalVariables::Scalar::DetF);

    auto& material_tangents = variables->fetch(InternalVariables::rank4::tangent_operator);

    for (auto& H : displacement_gradients) H = InternalVariables::rank2tensor_type::Zero();
    for (auto& J : J_list) J = 1.0;

    Eigen::EigenSolver<Eigen::MatrixXd> eigen_solver;

    SECTION("Sanity checks")
    {
        REQUIRE(elastic_model->is_symmetric());
        REQUIRE(elastic_model->is_finite_deformation() == false);
        REQUIRE(elastic_model->intrinsic_material().name() == "steel");

        REQUIRE(variables->has(InternalVariables::Scalar::VonMisesStress));
        REQUIRE(variables->has(InternalVariables::Tensor::Cauchy));
        REQUIRE(variables->has(InternalVariables::Tensor::LinearisedStrain));
        REQUIRE(variables->has(InternalVariables::rank4::tangent_operator));
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

        auto[von_mises_stresses,
             accumulated_plastic_strains] = variables
                                                ->fetch(InternalVariables::Scalar::VonMisesStress,
                                                        InternalVariables::Scalar::EffectivePlasticStrain);

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

    // Create a json reader object from a string
    std::string input_data = "{\"Name\": \"steel\", \"ElasticModulus\": 200.0e3, "
                             "\"PoissonsRatio\": 0.3}";

    std::string mesh_input = "{\"ConstitutiveModel\" : {\"Name\" : \"IsotropicLinearElasticity\"}}";

    Json::Value material_data, mesh_data;

    std::istringstream input_data_stream(input_data), mesh_input_stream(mesh_input);

    REQUIRE(Json::parseFromStream(reader, input_data_stream, &material_data, &input_errors));
    REQUIRE(Json::parseFromStream(reader, mesh_input_stream, &mesh_data, &input_errors));

    auto variables = std::make_shared<InternalVariables>(internal_variable_size);

    // Add the required variables for an updated Lagrangian formulation
    variables->add(InternalVariables::Tensor::DisplacementGradient,
                   InternalVariables::Tensor::Cauchy,
                   InternalVariables::Scalar::DetF);

    auto elastic_model = make_constitutive_model(variables, material_data, mesh_data);

    // Get the tensor variables
    auto[displacement_gradients,
         cauchy_stresses] = variables->fetch(InternalVariables::Tensor::DisplacementGradient,
                                             InternalVariables::Tensor::Cauchy);

    auto& J_list = variables->fetch(InternalVariables::Scalar::DetF);

    auto& material_tangents = variables->fetch(InternalVariables::rank4::tangent_operator);

    for (auto& H : displacement_gradients) H = InternalVariables::rank2tensor_type::Zero();
    for (auto& J : J_list) J = 1.0;

    Eigen::EigenSolver<Eigen::MatrixXd> eigen_solver;

    SECTION("Sanity checks")
    {
        REQUIRE(elastic_model->is_symmetric());
        REQUIRE(elastic_model->is_finite_deformation() == false);
        REQUIRE(elastic_model->intrinsic_material().name() == "steel");

        REQUIRE(variables->has(InternalVariables::Scalar::VonMisesStress));
        REQUIRE(variables->has(InternalVariables::Tensor::Cauchy));
        REQUIRE(variables->has(InternalVariables::Tensor::LinearisedStrain));
        REQUIRE(variables->has(InternalVariables::rank4::tangent_operator));
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

        auto[von_mises_stresses,
             accumulated_plastic_strains] = variables
                                                ->fetch(InternalVariables::Scalar::VonMisesStress,
                                                        InternalVariables::Scalar::EffectivePlasticStrain);

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

    // Create a json reader object from a string
    std::string input_data = "{}";
    std::string simulation_input = "{\"ConstitutiveModel\" : {\"Name\" : \"J2Plasticity\"}}";

    Json::Value material_data, simulation_data;

    std::istringstream input_data_stream(input_data), simulation_input_stream(simulation_input);

    REQUIRE(Json::parseFromStream(reader, input_data_stream, &material_data, &input_errors));
    REQUIRE(Json::parseFromStream(reader, simulation_input_stream, &simulation_data, &input_errors));

    auto variables = std::make_shared<InternalVariables>(internal_variable_size);

    REQUIRE_THROWS_AS(make_constitutive_model(variables, material_data, simulation_data),
                      std::runtime_error);
}
TEST_CASE("Solid mechanics J2 plasticity model")
{
    using namespace neon::mechanical::solid;

    // Create a json reader object from a string
    std::string input_data = "{\"Name\": \"steel\", \"ElasticModulus\": 200.0e9, "
                             "\"PoissonsRatio\": 0.3, "
                             "\"YieldStress\": 200.0e6, \"IsotropicHardeningModulus\": 400.0e6}";

    std::string simulation_input = "{\"ConstitutiveModel\" : {\"Name\" : \"J2Plasticity\", "
                                   "\"FiniteStrain\" : false}}";

    Json::Value material_data, simulation_data;

    std::istringstream input_data_stream(input_data), simulation_input_stream(simulation_input);

    REQUIRE(Json::parseFromStream(reader, input_data_stream, &material_data, &input_errors));
    REQUIRE(Json::parseFromStream(reader, simulation_input_stream, &simulation_data, &input_errors));

    auto variables = std::make_shared<InternalVariables>(internal_variable_size);

    // Add the required variables for an updated Lagrangian formulation
    variables->add(InternalVariables::Tensor::DisplacementGradient,
                   InternalVariables::Tensor::Cauchy);
    variables->add(InternalVariables::Scalar::DetF);

    auto j2plasticity = make_constitutive_model(variables, material_data, simulation_data);

    // Get the tensor variables
    auto[displacement_gradients,
         cauchy_stresses] = variables->fetch(InternalVariables::Tensor::DisplacementGradient,
                                             InternalVariables::Tensor::Cauchy);

    auto& J_list = variables->fetch(InternalVariables::Scalar::DetF);

    auto& material_tangents = variables->fetch(InternalVariables::rank4::tangent_operator);

    for (auto& H : displacement_gradients) H = neon::matrix3::Zero();
    for (auto& J : J_list) J = 1.0;

    Eigen::EigenSolver<Eigen::MatrixXd> eigen_solver;

    SECTION("Sanity checks")
    {
        REQUIRE(j2plasticity->is_symmetric());
        REQUIRE(j2plasticity->is_finite_deformation() == false);
        REQUIRE(j2plasticity->intrinsic_material().name() == "steel");

        REQUIRE(variables->has(InternalVariables::Scalar::VonMisesStress));
        REQUIRE(variables->has(InternalVariables::Scalar::EffectivePlasticStrain));
        REQUIRE(variables->has(InternalVariables::Tensor::LinearisedStrain));
        REQUIRE(variables->has(InternalVariables::Tensor::LinearisedPlasticStrain));
        REQUIRE(variables->has(InternalVariables::rank4::tangent_operator));
    }
    SECTION("No load")
    {
        j2plasticity->update_internal_variables(1.0);

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

        j2plasticity->update_internal_variables(1.0);

        auto[von_mises_stresses,
             accumulated_plastic_strains] = variables
                                                ->fetch(InternalVariables::Scalar::VonMisesStress,
                                                        InternalVariables::Scalar::EffectivePlasticStrain);

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

        j2plasticity->update_internal_variables(1.0);

        auto[von_mises_stresses,
             accumulated_plastic_strains] = variables
                                                ->fetch(InternalVariables::Scalar::VonMisesStress,
                                                        InternalVariables::Scalar::EffectivePlasticStrain);

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

    // Create a json reader object from a string
    std::string input_data = "{\"Name\": \"steel\", \"ElasticModulus\": 134.0e3, "
                             "\"PoissonsRatio\": 0.3, "
                             "\"YieldStress\": 85, \"KinematicHardeningModulus\": "
                             "5500,\"SofteningMultiplier\" : 250,\"PlasticityViscousExponent\" : "
                             "2.5,\"PlasticityViscousMultiplier\" : "
                             "1.923536463026969e-08,\"DamageViscousExponent\" : "
                             "2,\"DamageViscousMultiplier\" : 2.777777777777778}";

    std::string simulation_input = "{\"ConstitutiveModel\" : {\"Name\" : \"ChabocheDamage\", "
                                   "\"FiniteStrain\" : false}}";

    Json::Value material_data, simulation_data;

    std::istringstream input_data_stream(input_data), simulation_input_stream(simulation_input);

    REQUIRE(Json::parseFromStream(reader, input_data_stream, &material_data, &input_errors));
    REQUIRE(Json::parseFromStream(reader, simulation_input_stream, &simulation_data, &input_errors));

    auto variables = std::make_shared<InternalVariables>(internal_variable_size);

    variables->add(InternalVariables::Tensor::DisplacementGradient,
                   InternalVariables::Tensor::Cauchy);
    variables->add(InternalVariables::Scalar::DetF);

    auto J2PlasticityDamage = make_constitutive_model(variables, material_data, simulation_data);

    // Get the tensor variables
    auto[displacement_gradients,
         cauchy_stresses] = variables->fetch(InternalVariables::Tensor::DisplacementGradient,
                                             InternalVariables::Tensor::Cauchy);

    auto[J_list, damage_list] = variables->fetch(InternalVariables::Scalar::DetF,
                                                 InternalVariables::Scalar::Damage);

    auto& material_tangents = variables->fetch(InternalVariables::rank4::tangent_operator);

    for (auto& H : displacement_gradients) H = neon::matrix3::Zero();
    for (auto& J : J_list) J = 1.0;

    Eigen::EigenSolver<Eigen::MatrixXd> eigen_solver;

    SECTION("Sanity checks")
    {
        REQUIRE(J2PlasticityDamage->is_finite_deformation() == false);
        REQUIRE(J2PlasticityDamage->is_symmetric() == false);

        REQUIRE(J2PlasticityDamage->intrinsic_material().name() == "steel");

        REQUIRE(variables->has(InternalVariables::Scalar::VonMisesStress));
        REQUIRE(variables->has(InternalVariables::Scalar::EffectivePlasticStrain));
        REQUIRE(variables->has(InternalVariables::Tensor::LinearisedStrain));
        REQUIRE(variables->has(InternalVariables::Tensor::LinearisedPlasticStrain));
        REQUIRE(variables->has(InternalVariables::rank4::tangent_operator));
        REQUIRE(variables->has(InternalVariables::Scalar::Damage));
        REQUIRE(variables->has(InternalVariables::Scalar::EnergyReleaseRate));
        REQUIRE(variables->has(InternalVariables::Tensor::KinematicHardening));
        REQUIRE(variables->has(InternalVariables::Tensor::BackStress));
    }
    SECTION("Uniaxial elastic load")
    {
        for (auto& H : displacement_gradients) H(2, 2) = 0.0008;

        J2PlasticityDamage->update_internal_variables(1.0);

        auto[von_mises_stresses,
             accumulated_plastic_strains] = variables
                                                ->fetch(InternalVariables::Scalar::VonMisesStress,
                                                        InternalVariables::Scalar::EffectivePlasticStrain);

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

        J2PlasticityDamage->update_internal_variables(0.01); // time here is real (not pseudo time)

        auto[von_mises_stresses,
             accumulated_plastic_strains] = variables
                                                ->fetch(InternalVariables::Scalar::VonMisesStress,
                                                        InternalVariables::Scalar::EffectivePlasticStrain);

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

    auto variables = std::make_shared<InternalVariables>(internal_variable_size);

    // Create a json reader object from a string
    auto input_data = json_thermal_diffusion();

    std::string simulation_input = "{\"ConstitutiveModel\" : {\"Name\": \"IsotropicDiffusion\"} }";

    std::istringstream input_data_stream(input_data), simulation_input_stream(simulation_input);

    Json::Value material_data, simulation_data;

    REQUIRE(Json::parseFromStream(reader, input_data_stream, &material_data, &input_errors));
    REQUIRE(Json::parseFromStream(reader, simulation_input_stream, &simulation_data, &input_errors));

    auto thermal = make_constitutive_model(variables, material_data, simulation_data);

    thermal->update_internal_variables(1.0);

    SECTION("Sanity checks")
    {
        REQUIRE(thermal->is_symmetric());
        REQUIRE(thermal->is_finite_deformation() == false);
        REQUIRE(thermal->intrinsic_material().name() == "steel");

        REQUIRE(variables->has(InternalVariables::Tensor::Conductivity));
    }
    Eigen::EigenSolver<Eigen::MatrixXd> eigen_solver;
}

// TEST_CASE("Finite J2 plasticity model", "[FiniteJ2Plasticity]")
// {
//     using namespace neon::mechanical::solid;
//
//     // Create a json reader object from a string
//     std::string
//         input_data = "{\"Name\": \"steel\", \"ElasticModulus\": 200.0e9, \"PoissonsRatio\": 0.3,
//         "
//                      "\"YieldStress\": 200.0e6, \"IsotropicHardeningModulus\": 400.0e6}";
//
//     std::string simulation_input = "{\"ConstitutiveModel\" : {\"Name\" : \"J2Plasticity\", "
//                                    "\"FiniteStrain\":true}}";
//
//     Json::Value material_data, simulation_data;
//
//     Json::Reader reader;
//
//     REQUIRE(reader.parse(input_data.c_str(), material_data));
//     REQUIRE(reader.parse(simulation_input.c_str(), simulation_data));
//
//     InternalVariables variables(1);
//
//     // Add the required variables for an updated Lagrangian formulation
//     variables->add(InternalVariables::Tensor::DeformationGradient,
//     InternalVariables::Tensor::Cauchy); variables->add(InternalVariables::Scalar::DetF);
//
//     auto j2plasticity = make_constitutive_model(variables, material_data, simulation_data);
//
//     // Get the tensor variables
//     auto[F_list, cauchy_stresses] = variables->fetch(InternalVariables::Tensor::DeformationGradient,
//                                               InternalVariables::Tensor::Cauchy);
//
//     auto& J_list = variables->fetch(InternalVariables::Scalar::DetF);
//
//     auto& material_tangents = variables->fetch(InternalVariables::rank4::tangent_operator);
//
//     for (auto& F : F_list) F = neon::matrix3::Identity();
//     for (auto& J : J_list) J = 1.0;
//
//     variables.commit();
//
//     SECTION("Sanity checks")
//     {
//         REQUIRE(j2plasticity->is_finite_deformation() == true);
//         REQUIRE(j2plasticity->intrinsic_material().name() == "steel");
//
//         REQUIRE(variables->has(InternalVariables::Scalar::VonMisesStress));
//         REQUIRE(variables->has(InternalVariables::Scalar::EffectivePlasticStrain));
//         REQUIRE(variables->has(InternalVariables::Tensor::HenckyStrainElastic));
//         REQUIRE(variables->has(InternalVariables::rank4::tangent_operator));
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
//         j2plasticity->update_internal_variables(1.0);
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
//     //     j2plasticity->update_internal_variables(1.0);
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
//     //     for (auto& von_mises_stress : variables->fetch(InternalVariables::Scalar::VonMisesStress))
//     //     {
//     //         REQUIRE(von_mises_stress < 200.0e6);
//     //     }
//     // }
//     // SECTION("Plastic uniaxial elastic load")
//     // {
//     //     for (auto& F : F_list) F(2, 2) = 0.003;
//     //
//     //     j2plasticity->update_internal_variables(1.0);
//     //
//     //     auto[von_mises_stresses, accumulated_plastic_strains] =
//     //     variables->fetch(InternalVariables::Scalar::VonMisesStress,
//     // InternalVariables::Scalar::EffectivePlasticStrain);
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
