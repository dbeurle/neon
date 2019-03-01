
#include <catch2/catch.hpp>

#include "constitutive/internal_variables.hpp"
#include "constitutive/constitutive_model_factory.hpp"
#include "constitutive/mechanics/solid/gaussian_ageing_affine_microsphere.hpp"

#include "exceptions.hpp"
#include "numeric/dense_matrix.hpp"
#include "io/json.hpp"

#include <Eigen/Eigenvalues>

#include <iostream>

std::string json_input_file()
{
    return "{\"name\": \"rubber\", \"elastic_modulus\": 2.0, \"poissons_ratio\": 0.45}";
}

constexpr auto ZERO_MARGIN = 1.0e-5;

using neon::json;
using namespace neon;

TEST_CASE("Gaussian affine microsphere model with ageing")
{
    using namespace neon::mechanics::solid;

    auto variables = std::make_shared<internal_variables_t>(1);

    // Add the required variables for an updated Lagrangian formulation
    variables->add(variable::second::deformation_gradient, variable::second::cauchy_stress);

    variables->add(variable::scalar::DetF);

    auto const material_data{"{\"name\" : \"rubber\","
                             "\"shear_modulus\" : 2.0e6,"
                             "\"bulk_modulus\" : 100e6,"
                             "\"segments_per_chain\" : 50,"
                             "\"scission_probability\" : 1.0e-5,"
                             "\"recombination_probability\" : 1.0e-5}"};

    auto const constitutive_data{"{\"constitutive\" : {\"name\": \"microsphere\","
                                 "\"type\":\"affine\","
                                 "\"statistics\":\"gaussian\","
                                 "\"quadrature\":\"BO21\","
                                 "\"ageing\":\"BAND\"}}"};

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

        for (auto segment : variables->get(variable::scalar::active_segments))
        {
            REQUIRE(segment == Approx(50.0));
        }
        for (auto shear_modulus : variables->get(variable::scalar::active_shear_modulus))
        {
            REQUIRE(shear_modulus == Approx(2.0e6));
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
        std::fill(begin(F_list), end(F_list), matrix3::Identity());

        affine->update_internal_variables(1.0);

        for (auto& cauchy_stress : cauchy_stresses)
        {
            std::cout << cauchy_stress << "\n\n";
            REQUIRE(cauchy_stress.norm() == Approx(0.0).margin(ZERO_MARGIN));
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
    using namespace neon::mechanics::solid;

    std::cout << "Constant cross-linking stress check\n";

    auto variables = std::make_shared<internal_variables_t>(1);

    // Add the required variables for an updated Lagrangian formulation
    variables->add(variable::second::deformation_gradient, variable::second::cauchy_stress);

    variables->add(variable::scalar::DetF);

    auto const material_data{"{\"name\" : \"rubber\","
                             "\"shear_modulus\" : 2.0e6,"
                             "\"bulk_modulus\" : 100e6,"
                             "\"segments_per_chain\" : 50,"
                             "\"scission_probability\" : 0.0,"
                             "\"recombination_probability\" : 1.0e-5}"};

    auto const constitutive_data{"{\"constitutive\" : {\"name\": \"microsphere\","
                                 "\"type\":\"affine\","
                                 "\"statistics\":\"gaussian\","
                                 "\"quadrature\":\"BO21\","
                                 "\"ageing\":\"BAND\"}}"};

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
// TEST_CASE("Finite J2 plasticity")
// {
//     using namespace neon::mechanics::solid;
//
//     // Create a json reader object from a string
//     std::string
//         input_data = "{\"name\": \"steel\", \"elastic_modulus\": 200.0e9, \"poissons_ratio\": 0.3,
//         "
//                      "\"yield_stress\": 200.0e6, \"isotropic_hardening_modulus\": 400.0e6}";
//
//     std::string simulation_input = "{\"constitutive\" : {\"name\" : \"small_strain_J2_plasticity\", "
//                                    "\"finite_strain\":true}}";
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
//     auto small_strain_J2_plasticity = make_constitutive_model(variables, material_data,
//     simulation_data);
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
