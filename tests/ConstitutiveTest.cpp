
#define CATCH_CONFIG_MAIN

#include "catch.hpp"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <json/json.h>

#include "constitutive/InternalVariables.hpp"

#include "constitutive/AffineMicrosphere.hpp"
#include "constitutive/HyperElasticPlastic.hpp"
#include "constitutive/J2Plasticity.hpp"
#include "constitutive/NeoHooke.hpp"

#include "constitutive/ConstitutiveModelFactory.hpp"

std::string json_input_file()
{
    return "{\"Name\": \"rubber\", \"ElasticModulus\": 1.0, \"PoissonsRatio\": 2.0}";
}

constexpr auto internal_variable_size = 2;
constexpr auto ZERO_MARGIN = 1.0e-5;

TEST_CASE("No constitutive model error test")
{
    using namespace neon::mech::solid;

    InternalVariables variables(internal_variable_size);

    // Create a json reader object from a string
    std::string input_data = "{}";
    std::string simulation_input = "{}";

    Json::Value material_data, simulation_data;
    Json::Reader reader;

    REQUIRE(reader.parse(input_data.c_str(), material_data));
    REQUIRE(reader.parse(simulation_input.c_str(), simulation_data));

    REQUIRE_THROWS_AS(make_constitutive_model(variables, material_data, simulation_data),
                      std::runtime_error);
}
TEST_CASE("Constitutive model no name error test")
{
    using namespace neon::mech::solid;

    InternalVariables variables(internal_variable_size);

    // Create a json reader object from a string
    std::string input_data = "{}";
    std::string simulation_input = "{\"ConstitutiveModel\" : {}}";

    Json::Value material_data, simulation_data;
    Json::Reader reader;

    REQUIRE(reader.parse(input_data.c_str(), material_data));
    REQUIRE(reader.parse(simulation_input.c_str(), simulation_data));

    REQUIRE_THROWS_AS(make_constitutive_model(variables, material_data, simulation_data),
                      std::runtime_error);
}
TEST_CASE("Constitutive model invalid name error test")
{
    using namespace neon::mech::solid;

    InternalVariables variables(internal_variable_size);

    // Create a json reader object from a string
    std::string input_data = "{}";
    std::string simulation_input = "{\"ConstitutiveModel\" : {\"Name\":\"PurpleMonkey\"}}";

    Json::Value material_data, simulation_data;
    Json::Reader reader;

    REQUIRE(reader.parse(input_data.c_str(), material_data));
    REQUIRE(reader.parse(simulation_input.c_str(), simulation_data));

    REQUIRE_THROWS_AS(make_constitutive_model(variables, material_data, simulation_data),
                      std::runtime_error);
}
TEST_CASE("Neo-Hookean model")
{
    using namespace neon::mech::solid;
    InternalVariables variables(internal_variable_size);

    // Add the required variables for an updated Lagrangian formulation
    variables.add(InternalVariables::Tensor::DeformationGradient, InternalVariables::Tensor::Cauchy);
    variables.add(InternalVariables::Scalar::DetF);

    // Create a json reader object from a string
    auto input_data = json_input_file();

    std::string simulation_input = "{\"ConstitutiveModel\" : {\"Name\": \"NeoHooke\"} }";

    Json::Value material_data, simulation_data;
    Json::Reader reader;

    REQUIRE(reader.parse(input_data.c_str(), material_data));
    REQUIRE(reader.parse(simulation_input.c_str(), simulation_data));

    auto neo_hooke = make_constitutive_model(variables, material_data, simulation_data);

    REQUIRE(neo_hooke->is_finite_deformation());
    REQUIRE(neo_hooke->intrinsic_material().name() == "rubber");

    // Get the tensor variables
    auto [F_list, cauchy_list] = variables(InternalVariables::Tensor::DeformationGradient,
                                           InternalVariables::Tensor::Cauchy);

    auto& J_list = variables(InternalVariables::Scalar::DetF);

    // Ensure the internal variables are allocated correctly
    REQUIRE(F_list.size() == internal_variable_size);
    REQUIRE(cauchy_list.size() == internal_variable_size);
    REQUIRE(J_list.size() == internal_variable_size);

    // Fill with identity matrix
    for (auto& F : F_list) F = neon::Matrix3::Identity();
    for (auto& J : J_list) J = 1.0;

    neo_hooke->update_internal_variables(1.0);

    SECTION("Update of internal variables")
    {
        for (auto& cauchy : cauchy_list) REQUIRE(cauchy.norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("Check of continuum tangent")
    {
        // Get the matrix variable
        auto& material_tangents = variables(InternalVariables::Matrix::TangentOperator);

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
            REQUIRE((C - C.transpose()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        }
    }
}
TEST_CASE("Microsphere model error test")
{
    using namespace neon::mech::solid;

    InternalVariables variables(internal_variable_size);

    // Create a json reader object from a string
    std::string input_data = "{}";

    std::string simulation_input = "{\"ConstitutiveModel\" : {\"Name\": \"Microsphere\", \"Type\" "
                                   ": \"Afwsfine\"}}";

    Json::Value material_data, simulation_data;
    Json::Reader reader;

    REQUIRE(reader.parse(input_data.c_str(), material_data));
    REQUIRE(reader.parse(simulation_input.c_str(), simulation_data));

    REQUIRE_THROWS_AS(make_constitutive_model(variables, material_data, simulation_data),
                      std::runtime_error);
}
TEST_CASE("Affine microsphere model", )
{
    using namespace neon::mech::solid;

    InternalVariables variables(internal_variable_size);

    // Add the required variables for an updated Lagrangian formulation
    variables.add(InternalVariables::Tensor::DeformationGradient, InternalVariables::Tensor::Cauchy);

    variables.add(InternalVariables::Scalar::DetF);

    // Create a json reader object from a string
    std::string input_data = "{\"Name\" : \"rubber\", "
                             "\"ElasticModulus\" : 10.0e6, "
                             "\"PoissonsRatio\" : 0.45, "
                             "\"SegmentsPerChain\" : 50}";

    std::string simulation_input = "{\"ConstitutiveModel\" : {\"Name\": \"Microsphere\", \"Type\" "
                                   ": \"Affine\", \"Quadrature\" : \"BO21\"}}";

    Json::Value material_data, simulation_data;
    Json::Reader reader;

    REQUIRE(reader.parse(input_data.c_str(), material_data));
    REQUIRE(reader.parse(simulation_input.c_str(), simulation_data));

    auto affine = make_constitutive_model(variables, material_data, simulation_data);

    REQUIRE(affine->is_finite_deformation());
    REQUIRE(affine->intrinsic_material().name() == "rubber");

    // Get the tensor variables
    auto [F_list, cauchy_list] = variables(InternalVariables::Tensor::DeformationGradient,
                                           InternalVariables::Tensor::Cauchy);

    auto& J_list = variables(InternalVariables::Scalar::DetF);

    for (auto& J : J_list) J = 1.0;

    auto& material_tangents = variables(InternalVariables::Matrix::TangentOperator);

    SECTION("Affine model under no load")
    {
        // Fill with identity matrix
        for (auto& F : F_list) F = neon::Matrix3::Identity();

        affine->update_internal_variables(1.0);

        for (auto const& C : material_tangents)
        {
            REQUIRE(C.norm() != Approx(0.0).margin(ZERO_MARGIN));
        }

        // Ensure symmetry is correct
        for (auto const& C : material_tangents)
        {
            REQUIRE((C - C.transpose()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        }

        for (auto& cauchy : cauchy_list)
        {
            REQUIRE(cauchy.norm() == Approx(0.0).margin(ZERO_MARGIN));
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

        // Ensure symmetry is correct
        for (auto const& C : material_tangents)
        {
            REQUIRE((C - C.transpose()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        }

        for (auto const& C : material_tangents)
        {
            REQUIRE(C.norm() != Approx(0.0).margin(ZERO_MARGIN));
        }

        for (auto& cauchy : cauchy_list)
        {
            REQUIRE(cauchy.norm() != Approx(0.0).margin(ZERO_MARGIN));
        }
    }
}
TEST_CASE("NonAffine microsphere model")
{
    using namespace neon::mech::solid;

    InternalVariables variables(internal_variable_size);

    // Add the required variables for an updated Lagrangian formulation
    variables.add(InternalVariables::Tensor::DeformationGradient, InternalVariables::Tensor::Cauchy);

    variables.add(InternalVariables::Scalar::DetF);

    // Create a json reader object from a string
    std::string input_data = "{\"Name\" : \"rubber\", "
                             "\"ElasticModulus\" : 10.0e6, "
                             "\"PoissonsRatio\" : 0.45, "
                             "\"NonAffineStretchParameter\":1.0, "
                             "\"SegmentsPerChain\" : 50}";

    std::string simulation_input = "{\"ConstitutiveModel\" : {\"Name\": \"Microsphere\", \"Type\" "
                                   ": \"NonAffine\", \"Quadrature\" : \"BO21\"}}";

    Json::Value material_data, simulation_data;
    Json::Reader reader;

    REQUIRE(reader.parse(input_data.c_str(), material_data));
    REQUIRE(reader.parse(simulation_input.c_str(), simulation_data));

    auto affine = make_constitutive_model(variables, material_data, simulation_data);

    REQUIRE(affine->is_finite_deformation());
    REQUIRE(affine->intrinsic_material().name() == "rubber");

    // Get the tensor variables
    auto [F_list, cauchy_list] = variables(InternalVariables::Tensor::DeformationGradient,
                                           InternalVariables::Tensor::Cauchy);

    auto& J_list = variables(InternalVariables::Scalar::DetF);

    for (auto& J : J_list) J = 1.0;

    auto& material_tangents = variables(InternalVariables::Matrix::TangentOperator);

    SECTION("NonAffine model under no load")
    {
        // Fill with identity matrix
        for (auto& F : F_list) F = neon::Matrix3::Identity();

        affine->update_internal_variables(1.0);

        for (auto const& C : material_tangents)
        {
            REQUIRE(C.norm() != Approx(0.0).margin(ZERO_MARGIN));
        }

        // Ensure symmetry is correct
        for (auto const& C : material_tangents)
        {
            REQUIRE((C - C.transpose()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        }

        for (auto& cauchy : cauchy_list)
        {
            REQUIRE(cauchy.norm() == Approx(0.0).margin(ZERO_MARGIN));
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

        // Ensure symmetry is correct
        for (auto const& C : material_tangents)
        {
            REQUIRE((C - C.transpose()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        }

        for (auto const& C : material_tangents)
        {
            REQUIRE(C.norm() != Approx(0.0).margin(ZERO_MARGIN));
        }

        for (auto& cauchy : cauchy_list)
        {
            REQUIRE(cauchy.norm() != Approx(0.0).margin(ZERO_MARGIN));
        }
    }
}
TEST_CASE("J2 plasticity model factory errors")
{
    using namespace neon::mech::solid;

    // Create a json reader object from a string
    std::string input_data = "{}";
    std::string simulation_input = "{\"ConstitutiveModel\" : {\"Name\" : \"J2Plasticity\"}}";

    Json::Value material_data, simulation_data;

    Json::Reader reader;

    REQUIRE(reader.parse(input_data.c_str(), material_data));
    REQUIRE(reader.parse(simulation_input.c_str(), simulation_data));

    InternalVariables variables(internal_variable_size);

    REQUIRE_THROWS_AS(make_constitutive_model(variables, material_data, simulation_data),
                      std::runtime_error);
}
TEST_CASE("J2 plasticity model")
{
    using namespace neon::mech::solid;

    // Create a json reader object from a string
    std::string
        input_data = "{\"Name\": \"steel\", \"ElasticModulus\": 200.0e9, \"PoissonsRatio\": 0.3, "
                     "\"YieldStress\": 200.0e6, \"IsotropicHardeningModulus\": 400.0e6}";

    std::string simulation_input = "{\"ConstitutiveModel\" : {\"Name\" : \"J2Plasticity\", "
                                   "\"FiniteStrain\" : false}}";

    Json::Value material_data, simulation_data;

    Json::Reader reader;

    REQUIRE(reader.parse(input_data.c_str(), material_data));
    REQUIRE(reader.parse(simulation_input.c_str(), simulation_data));

    InternalVariables variables(internal_variable_size);

    // Add the required variables for an updated Lagrangian formulation
    variables.add(InternalVariables::Tensor::DisplacementGradient, InternalVariables::Tensor::Cauchy);
    variables.add(InternalVariables::Scalar::DetF);

    auto j2plasticity = make_constitutive_model(variables, material_data, simulation_data);

    // Get the tensor variables
    auto [H_list, cauchy_list] = variables(InternalVariables::Tensor::DisplacementGradient,
                                           InternalVariables::Tensor::Cauchy);

    auto& J_list = variables(InternalVariables::Scalar::DetF);

    auto& material_tangents = variables(InternalVariables::Matrix::TangentOperator);

    for (auto& H : H_list) H = neon::Matrix3::Zero();
    for (auto& J : J_list) J = 1.0;

    SECTION("Sanity checks")
    {
        REQUIRE(j2plasticity->is_finite_deformation() == false);
        REQUIRE(j2plasticity->intrinsic_material().name() == "steel");

        REQUIRE(variables.has(InternalVariables::Scalar::VonMisesStress));
        REQUIRE(variables.has(InternalVariables::Scalar::EffectivePlasticStrain));
        REQUIRE(variables.has(InternalVariables::Tensor::LinearisedStrain));
        REQUIRE(variables.has(InternalVariables::Tensor::LinearisedPlasticStrain));
        REQUIRE(variables.has(InternalVariables::Matrix::TangentOperator));
    }
    SECTION("No load")
    {
        j2plasticity->update_internal_variables(1.0);

        // Ensure symmetry is correct
        for (auto const& C : material_tangents)
        {
            REQUIRE((C - C.transpose()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        }
        for (auto& cauchy : cauchy_list) REQUIRE(cauchy.norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("Uniaxial elastic load")
    {
        for (auto& H : H_list) H(2, 2) = 0.001;

        j2plasticity->update_internal_variables(1.0);

        auto [vm_list,
              eff_plastic_list] = variables(InternalVariables::Scalar::VonMisesStress,
                                            InternalVariables::Scalar::EffectivePlasticStrain);

        // Ensure symmetry is correct
        for (auto const& C : material_tangents)
        {
            REQUIRE((C - C.transpose()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        }
        for (auto& cauchy : cauchy_list)
        {
            REQUIRE(cauchy.norm() != Approx(0.0).margin(ZERO_MARGIN));

            // Shear components should be close to zero
            REQUIRE(cauchy(0, 1) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(cauchy(0, 2) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(cauchy(1, 2) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(cauchy(0, 0) > 0.0);
            REQUIRE(cauchy(1, 1) > 0.0);
            REQUIRE(cauchy(2, 2) > 0.0);
        }
        for (auto& strain_p_eff : eff_plastic_list)
            REQUIRE(strain_p_eff == Approx(0.0).margin(ZERO_MARGIN));
        for (auto& vm : vm_list) REQUIRE(vm < 200.0e6);
    }
    SECTION("Plastic uniaxial elastic load")
    {
        for (auto& H : H_list) H(2, 2) = 0.003;

        j2plasticity->update_internal_variables(1.0);

        auto [vm_list,
              eff_plastic_list] = variables(InternalVariables::Scalar::VonMisesStress,
                                            InternalVariables::Scalar::EffectivePlasticStrain);

        // Ensure symmetry is correct
        for (auto const& C : material_tangents)
        {
            REQUIRE((C - C.transpose()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        }
        for (auto& cauchy : cauchy_list)
        {
            REQUIRE(cauchy.norm() != Approx(0.0).margin(ZERO_MARGIN));

            // Shear components should be close to zero
            REQUIRE(cauchy(0, 1) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(cauchy(0, 2) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(cauchy(1, 2) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(cauchy(0, 0) > 0.0);
            REQUIRE(cauchy(1, 1) > 0.0);
            REQUIRE(cauchy(2, 2) > 0.0);
        }
        for (auto& strain_p_eff : eff_plastic_list)
        {
            REQUIRE(strain_p_eff > 0.0);
        }
        for (auto& vm : vm_list)
        {
            // Should experience hardening
            REQUIRE(vm <= 201.0e6);
            REQUIRE(vm > 200.0e6);
        }
    }
}
TEST_CASE("J2 plasticity damage model", "[J2PlasticityDamage]")
{
    using namespace neon::mech::solid;

    // Create a json reader object from a string
    std::string
        input_data = "{\"Name\": \"steel\", \"ElasticModulus\": 134.0e3, \"PoissonsRatio\": 0.3, "
                     "\"YieldStress\": 85, \"HardeningModulus\": "
                     "5500,\"SofteningMultiplier\" : 250,\"PlasticityViscousExponent\" : "
                     "2.5,\"PlasticityViscousMultiplier\" : "
                     "1.923536463026969e-08,\"DamageViscousExponent\" : "
                     "2,\"DamageViscousMultiplier\" : 2.777777777777778}";

    std::string simulation_input = "{\"ConstitutiveModel\" : {\"Name\" : \"ChabocheDamage\", "
                                   "\"FiniteStrain\" : false}}";

    Json::Value material_data, simulation_data;

    Json::Reader reader;

    REQUIRE(reader.parse(input_data.c_str(), material_data));
    REQUIRE(reader.parse(simulation_input.c_str(), simulation_data));

    InternalVariables variables(internal_variable_size);

    variables.add(InternalVariables::Tensor::DisplacementGradient, InternalVariables::Tensor::Cauchy);
    variables.add(InternalVariables::Scalar::DetF);

    auto J2PlasticityDamage = make_constitutive_model(variables, material_data, simulation_data);

    // Get the tensor variables
    auto [H_list, cauchy_list] = variables(InternalVariables::Tensor::DisplacementGradient,
                                           InternalVariables::Tensor::Cauchy);

    auto [J_list, damage_list] = variables(InternalVariables::Scalar::DetF,
                                           InternalVariables::Scalar::Damage);

    auto& material_tangents = variables(InternalVariables::Matrix::TangentOperator);

    for (auto& H : H_list) H = neon::Matrix3::Zero();
    for (auto& J : J_list) J = 1.0;

    SECTION("Sanity checks")
    {
        REQUIRE(J2PlasticityDamage->is_finite_deformation() == false);
        REQUIRE(J2PlasticityDamage->intrinsic_material().name() == "steel");

        REQUIRE(variables.has(InternalVariables::Scalar::VonMisesStress));
        REQUIRE(variables.has(InternalVariables::Scalar::EffectivePlasticStrain));
        REQUIRE(variables.has(InternalVariables::Tensor::LinearisedStrain));
        REQUIRE(variables.has(InternalVariables::Tensor::LinearisedPlasticStrain));
        REQUIRE(variables.has(InternalVariables::Matrix::TangentOperator));
        REQUIRE(variables.has(InternalVariables::Scalar::Damage));
        REQUIRE(variables.has(InternalVariables::Scalar::EnergyReleaseRate));
        REQUIRE(variables.has(InternalVariables::Tensor::KinematicHardening));
        REQUIRE(variables.has(InternalVariables::Tensor::BackStress));
    }
    SECTION("Uniaxial elastic load")
    {
        for (auto& H : H_list) H(2, 2) = 0.0008;

        J2PlasticityDamage->update_internal_variables(1.0);

        auto [vm_list,
              eff_plastic_list] = variables(InternalVariables::Scalar::VonMisesStress,
                                            InternalVariables::Scalar::EffectivePlasticStrain);

        // Ensure symmetry is correct when the loading is elastic
        for (auto const& C : material_tangents)
        {
            // std::cout << C << "\n";
            REQUIRE((C - C.transpose()).norm() == Approx(0.0));
        }

        for (auto& cauchy : cauchy_list)
        {
            // std::cout << cauchy << std::endl;
            REQUIRE(cauchy.norm() != Approx(0.0));

            // Shear components should be close to zero
            REQUIRE(cauchy(0, 1) == Approx(0.0));
            REQUIRE(cauchy(0, 2) == Approx(0.0));
            REQUIRE(cauchy(1, 2) == Approx(0.0));
            REQUIRE(cauchy(0, 0) == Approx(cauchy(1, 1)));
            REQUIRE(cauchy(1, 1) > 0.0);
            REQUIRE(cauchy(2, 2) > 0.0);
        }
        for (auto& strain_p_eff : eff_plastic_list) REQUIRE(strain_p_eff == Approx(0.0));
        for (auto& vm : vm_list)
        {
            REQUIRE(vm > 80);
            REQUIRE(vm < 85);
        }
    }
    SECTION("Plastic uniaxial load")
    {
        for (auto& H : H_list) H(2, 2) = 0.0009; // 0.000825 for comparison with 1D

        J2PlasticityDamage->update_internal_variables(0.01); // time here is real (not pseudo time)

        auto [vm_list,
              eff_plastic_list] = variables(InternalVariables::Scalar::VonMisesStress,
                                            InternalVariables::Scalar::EffectivePlasticStrain);

        for (auto const& cauchy : cauchy_list)
        {
            REQUIRE(cauchy.norm() != Approx(0.0));
            // std::cout << cauchy << std::endl;
            // Shear components should be close to zero
            REQUIRE(cauchy(0, 1) == Approx(0.0));
            REQUIRE(cauchy(0, 2) == Approx(0.0));
            REQUIRE(cauchy(1, 2) == Approx(0.0));
            REQUIRE(cauchy(0, 0) > 0.0);
            REQUIRE(cauchy(1, 1) > 0.0);
            REQUIRE(cauchy(2, 2) > 0.0);
        }
        for (auto& strain_p_eff : eff_plastic_list)
        {
            REQUIRE(strain_p_eff > 0.0);
        }
        for (auto& damage_var : damage_list)
        {
            REQUIRE(damage_var > 0.0);
        }
        for (auto& vm : vm_list)
        {
            // Should experience hardening
            REQUIRE(vm <= 100);
            REQUIRE(vm > 85);
        }
        // check the eig val of the tangent_operator
        for (auto const& C : material_tangents)
        {
            // std::cout << C << "\n";
            // Compute the condition number of C
            // Eigen::JacobiSVD<Eigen::MatrixXd> svd(C);
            // double cond = svd.singularValues()(0)
            //               / svd.singularValues()(svd.singularValues().size() - 1);
            // std::cout << cond << "\n\n";

            Eigen::EigenSolver<Eigen::MatrixXd> eigen_solver;
            eigen_solver.compute(C);
            // Eigen::VectorXd eig_val = eigen_solver.eigenvalues();
            // std::cout << eig_val << "\n";
            // REQUIRE(eig_val(0) >= 0.0); // TODO: uncomment this linear and add the other elements
        }
    }
}
TEST_CASE("Finite J2 plasticity model", "[FiniteJ2Plasticity]")
{
    using namespace neon::mech::solid;

    // Create a json reader object from a string
    std::string
        input_data = "{\"Name\": \"steel\", \"ElasticModulus\": 200.0e9, \"PoissonsRatio\": 0.3, "
                     "\"YieldStress\": 200.0e6, \"IsotropicHardeningModulus\": 400.0e6}";

    std::string simulation_input = "{\"ConstitutiveModel\" : {\"Name\" : \"J2Plasticity\", "
                                   "\"FiniteStrain\":true}}";

    Json::Value material_data, simulation_data;

    Json::Reader reader;

    REQUIRE(reader.parse(input_data.c_str(), material_data));
    REQUIRE(reader.parse(simulation_input.c_str(), simulation_data));

    InternalVariables variables(1);

    // Add the required variables for an updated Lagrangian formulation
    variables.add(InternalVariables::Tensor::DeformationGradient, InternalVariables::Tensor::Cauchy);
    variables.add(InternalVariables::Scalar::DetF);

    auto j2plasticity = make_constitutive_model(variables, material_data, simulation_data);

    // Get the tensor variables
    auto [F_list, cauchy_list] = variables(InternalVariables::Tensor::DeformationGradient,
                                           InternalVariables::Tensor::Cauchy);

    auto& J_list = variables(InternalVariables::Scalar::DetF);

    auto& material_tangents = variables(InternalVariables::Matrix::TangentOperator);

    for (auto& F : F_list) F = neon::Matrix3::Identity();
    for (auto& J : J_list) J = 1.0;

    variables.commit();

    SECTION("Sanity checks")
    {
        REQUIRE(j2plasticity->is_finite_deformation() == true);
        REQUIRE(j2plasticity->intrinsic_material().name() == "steel");

        REQUIRE(variables.has(InternalVariables::Scalar::VonMisesStress));
        REQUIRE(variables.has(InternalVariables::Scalar::EffectivePlasticStrain));
        REQUIRE(variables.has(InternalVariables::Tensor::HenckyStrainElastic));
        REQUIRE(variables.has(InternalVariables::Matrix::TangentOperator));
    }
    SECTION("Initial material tangent symmetry")
    {
        for (auto const& C : material_tangents)
        {
            REQUIRE((C - C.transpose()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        }
    }
    SECTION("No load")
    {
        j2plasticity->update_internal_variables(1.0);

        for (auto const& C : material_tangents)
        {
            REQUIRE((C - C.transpose()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        }
        for (auto& cauchy : cauchy_list)
        {
            REQUIRE(cauchy.norm() == Approx(0.0).margin(ZERO_MARGIN));
        }
    }
    // SECTION("Uniaxial elastic load")
    // {
    //     for (auto& F : F_list) F(2, 2) = 1.001;
    //
    //     j2plasticity->update_internal_variables(1.0);
    //
    //     for (auto const& C : material_tangents)
    //     {
    //         REQUIRE((C - C.transpose()).norm() == Approx(0.0).margin(ZERO_MARGIN));
    //     }
    //     for (auto& cauchy : cauchy_list)
    //     {
    //         REQUIRE((cauchy - cauchy.transpose()).norm() == Approx(0.0).margin(ZERO_MARGIN));
    //
    //         // Shear components should be close to zero
    //         REQUIRE(cauchy(0, 1) == Approx(0.0).margin(ZERO_MARGIN));
    //         REQUIRE(cauchy(0, 2) == Approx(0.0).margin(ZERO_MARGIN));
    //         REQUIRE(cauchy(1, 2) == Approx(0.0).margin(ZERO_MARGIN));
    //         REQUIRE(cauchy(0, 0) > 0.0);
    //         REQUIRE(cauchy(1, 1) > 0.0);
    //         REQUIRE(cauchy(2, 2) > 0.0);
    //     }
    //     for (auto& vm : variables(InternalVariables::Scalar::VonMisesStress))
    //     {
    //         REQUIRE(vm < 200.0e6);
    //     }
    // }
    // SECTION("Plastic uniaxial elastic load")
    // {
    //     for (auto& F : F_list) F(2, 2) = 0.003;
    //
    //     j2plasticity->update_internal_variables(1.0);
    //
    //     auto[vm_list, eff_plastic_list] = variables(InternalVariables::Scalar::VonMisesStress,
    //                                                 InternalVariables::Scalar::EffectivePlasticStrain);
    //
    //     // Ensure symmetry is correct
    //     for (auto const& C : material_tangents)
    //     {
    //         REQUIRE((C - C.transpose()).norm() == Approx(0.0).margin(ZERO_MARGIN));
    //     }
    //     for (auto& cauchy : cauchy_list)
    //     {
    //         REQUIRE(cauchy.norm() != Approx(0.0).margin(ZERO_MARGIN));
    //
    //         // Shear components should be close to zero
    //         REQUIRE(cauchy(0, 1) == Approx(0.0).margin(ZERO_MARGIN));
    //         REQUIRE(cauchy(0, 2) == Approx(0.0).margin(ZERO_MARGIN));
    //         REQUIRE(cauchy(1, 2) == Approx(0.0).margin(ZERO_MARGIN));
    //         REQUIRE(cauchy(0, 0) > 0.0);
    //         REQUIRE(cauchy(1, 1) > 0.0);
    //         REQUIRE(cauchy(2, 2) > 0.0);
    //     }
    //     for (auto& strain_p_eff : eff_plastic_list)
    //     {
    //         REQUIRE(strain_p_eff > 0.0);
    //     }
    //     for (auto& vm : vm_list)
    //     {
    //         // Should experience hardening
    //         REQUIRE(vm <= 201.0e6);
    //         REQUIRE(vm > 200.0e6);
    //     }
    // }
}
