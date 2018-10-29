
#include <catch2/catch.hpp>

#include "constitutive/internal_variables.hpp"
#include "constitutive/mechanics/solid/small_strain_J2_plasticity.hpp"
#include "constitutive/mechanics/detail/J2_plasticity.hpp"
#include "constitutive/constitutive_model_factory.hpp"
#include "material/isotropic_elastic_plastic.hpp"

#include "exceptions.hpp"
#include "numeric/dense_matrix.hpp"
#include "io/json.hpp"

#include <Eigen/Eigenvalues>

std::string json_input_file()
{
    return "{\"name\": \"rubber\", \"elastic_modulus\": 2.0, \"poissons_ratio\": 0.45}";
}

constexpr auto internal_variable_size = 2;
constexpr auto ZERO_MARGIN = 1.0e-5;

using namespace neon;
using neon::json;

TEST_CASE("Plane stress linear elasticity factory error")
{
    using namespace neon::mechanics::plane;

    auto variables = std::make_shared<internal_variables_t>(internal_variable_size);

    REQUIRE_THROWS_AS(make_constitutive_model(variables,
                                              json::parse("{\"Badkey\" : \"donkey\"}"),
                                              json::parse("{\"constitutive\" : {\"name\" : "
                                                          "\"plane_strain\"}}")),
                      std::domain_error);
}
TEST_CASE("Plane stress elasticity")
{
    using namespace neon::mechanics::plane;

    auto variables = std::make_shared<internal_variables_t>(internal_variable_size);

    // Add the required variables for an updated Lagrangian formulation
    variables->add(variable::second::displacement_gradient,
                   variable::second::cauchy_stress,
                   variable::scalar::DetF);

    auto elastic_model = make_constitutive_model(variables,
                                                 json::parse("{\"name\": \"steel\", "
                                                             "\"elastic_modulus\": 200.0e3, "
                                                             "\"poissons_ratio\": 0.3}"),
                                                 json::parse("{\"constitutive\" : {\"name\" : "
                                                             "\"plane_stress\"}}"));

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
TEST_CASE("Plane strain elasticity")
{
    using namespace neon::mechanics::plane;

    auto variables = std::make_shared<internal_variables_t>(internal_variable_size);

    // Add the required variables for an updated Lagrangian formulation
    variables->add(variable::second::displacement_gradient,
                   variable::second::cauchy_stress,
                   variable::scalar::DetF);

    auto elastic_model = make_constitutive_model(variables,
                                                 json::parse("{\"name\": \"steel\", "
                                                             "\"elastic_modulus\": 200.0e3, "
                                                             "\"poissons_ratio\": 0.3}"),
                                                 json::parse("{\"constitutive\" : {\"name\" : "
                                                             "\"plane_strain\"}}"));

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
    using namespace neon::mechanics::solid;

    auto variables = std::make_shared<internal_variables_t>(internal_variable_size);

    // Add the required variables for an updated Lagrangian formulation
    variables->add(variable::second::displacement_gradient,
                   variable::second::cauchy_stress,
                   variable::scalar::DetF);

    auto elastic_model = make_constitutive_model(variables,
                                                 json::parse("{\"name\": \"steel\", "
                                                             "\"elastic_modulus\": 200.0e3, "
                                                             "\"poissons_ratio\": 0.3}"),
                                                 json::parse("{\"constitutive\" : {\"name\" : "
                                                             "\"isotropic_linear_elasticity\"}}"));

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
    using namespace neon::mechanics::solid;

    auto variables = std::make_shared<internal_variables_t>(internal_variable_size);

    SECTION("finite_strain not specified")
    {
        REQUIRE_THROWS_AS(make_constitutive_model(variables,
                                                  json::parse("{}"),
                                                  json::parse("{\"constitutive\" : {\"name\" "
                                                              ": "
                                                              "\"J2_plasticity\"}}")),
                          std::domain_error);
    }
    SECTION("Damage model not specified correctly")
    {
        REQUIRE_THROWS_AS(make_constitutive_model(variables,
                                                  json::parse("{}"),
                                                  json::parse("{\"constitutive\" : {\"name\" "
                                                              ": "
                                                              "\"J2_plasticity\", "
                                                              "\"damage\" : "
                                                              "\"isotropic_chaboche\", "
                                                              "\"finite_strain\" : true}}")),
                          std::domain_error);
    }
}
TEST_CASE("Solid mechanics J2 plasticity")
{
    using namespace neon::mechanics::solid;

    auto variables = std::make_shared<internal_variables_t>(internal_variable_size);

    // Add the required variables for an updated Lagrangian formulation
    variables->add(variable::second::displacement_gradient, variable::second::cauchy_stress);
    variables->add(variable::scalar::DetF);

    auto const material_data = json::parse("{\"name\": \"steel\", "
                                           "\"elastic_modulus\": 200.0e9, "
                                           "\"poissons_ratio\": 0.3, "
                                           "\"yield_stress\": 200.0e6, "
                                           "\"isotropic_hardening_modulus\": 400.0e6}");

    auto const constitutive_data = json::parse("{\"constitutive\" : "
                                               "{\"name\" : \"J2_plasticity\", "
                                               "\"finite_strain\" : false}}");

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

        REQUIRE(neon::mechanics::evaluate_J2_yield_function(material, 100.0e6, 0.0) <= 0.0);
        REQUIRE(neon::mechanics::evaluate_J2_yield_function(material, 300.0e6, 0.0) > 0.0);
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
TEST_CASE("Solid mechanics J2 plasticity damage")
{
    using namespace neon::mechanics::solid;

    auto variables = std::make_shared<internal_variables_t>(internal_variable_size);

    variables->add(variable::second::displacement_gradient, variable::second::cauchy_stress);
    variables->add(variable::scalar::DetF);

    auto const material_input{"{\"name\":\"steel\","
                              "\"elastic_modulus\":134.0e3,"
                              "\"poissons_ratio\":0.3,"
                              "\"yield_stress\":85,"
                              "\"kinematic_hardening_modulus\": 5500,"
                              "\"softening_multiplier\" : 250,"
                              "\"plasticity_viscous_exponent\" : 2.5,"
                              "\"plasticity_viscous_denominator\" : 1220,"
                              "\"damage_viscous_exponent\" : 2,"
                              "\"damage_viscous_denominator\" : 0.6"
                              "}"};

    auto const constitutive_input{"{\"constitutive\" : "
                                  "{\"name\" : \"J2_plasticity\","
                                  "\"damage\" : \"isotropic_chaboche\", "
                                  "\"finite_strain\" : false}}"};

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
