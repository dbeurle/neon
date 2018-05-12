
#include "constitutive/mechanical/plane/finite_strain_J2_plasticity.hpp"

#include "exceptions.hpp"
#include "constitutive/internal_variables.hpp"

#include "numeric/float_compare.hpp"
#include "numeric/log_tensor_derivative.hpp"

#include "numeric/mechanics"

#include "constitutive/mechanical/detail/J2_plasticity.hpp"

#include <range/v3/view/transform.hpp>
#include <range/v3/view/zip.hpp>

#include <unsupported/Eigen/MatrixFunctions>

namespace neon::mechanical::plane
{
finite_strain_J2_plasticity::finite_strain_J2_plasticity(
    std::shared_ptr<internal_variables_t>& variables, json const& material_data)
    : small_strain_J2_plasticity(variables, material_data)
{
    variables->add(internal_variables_t::scalar::VonMisesStress,
                   internal_variables_t::scalar::EffectivePlasticStrain,
                   internal_variables_t::second::HenckyStrainElastic);

    // Add material tangent with the linear elasticity moduli
    variables->add(internal_variables_t::fourth::tangent_operator,
                   consistent_tangent(1.0, matrix2::Zero(), matrix2::Zero(), C_e));
}

finite_strain_J2_plasticity::~finite_strain_J2_plasticity() = default;

void finite_strain_J2_plasticity::update_internal_variables(double const time_step_size)
{
    using namespace ranges;

    auto const shear_modulus = material.shear_modulus();

    // Extract the internal variables
    auto [deformation_gradients,
          log_strain_e_list,
          cauchy_stresses] = variables->fetch(internal_variables_t::second::DeformationGradient,
                                              internal_variables_t::second::HenckyStrainElastic,
                                              internal_variables_t::second::CauchyStress);

    auto const old_deformation_gradients = variables->fetch_old(
        internal_variables_t::second::DeformationGradient);

    auto const J_list = variables->fetch(internal_variables_t::scalar::DetF);

    // Retrieve the accumulated internal variables
    auto [accumulated_plastic_strains,
          von_mises_stresses] = variables->fetch(internal_variables_t::scalar::EffectivePlasticStrain,
                                                 internal_variables_t::scalar::VonMisesStress);

    auto& tangent_operators = variables->fetch(internal_variables_t::fourth::tangent_operator);

    auto const incremental_deformation_gradients = view::zip(deformation_gradients,
                                                             old_deformation_gradients)
                                                   | view::transform([](auto const& tpl) {
                                                         auto const& [F, F_old] = tpl;
                                                         return F * F_old.inverse();
                                                     });

    // Perform the update algorithm for each quadrature point
    for (std::size_t l{0}; l < deformation_gradients.size(); l++)
    {
        auto const& F_inc = incremental_deformation_gradients[l];
        auto const J = J_list[l];

        auto& cauchy_stress = cauchy_stresses[l];
        auto& accumulated_plastic_strain = accumulated_plastic_strains[l];
        auto& von_mises = von_mises_stresses[l];
        auto& log_strain_e = log_strain_e_list[l];

        // std::cout << "log strain elastic\n" << log_strain_e << std::endl;

        // Elastic trial deformation gradient
        matrix2 const B_e = (2.0 * log_strain_e).exp();

        // std::cout << "B elastic\n" << B_e << std::endl;

        // std::cout << "Finc\n" << F_inc << std::endl;

        // Elastic trial left Cauchy-Green deformation tensor
        matrix2 const B_e_trial = F_inc * B_e * F_inc.transpose();

        // std::cout << "B elastic trial\n" << B_e_trial << std::endl;

        // Trial Logarithmic elastic strain
        log_strain_e = 0.5 * B_e_trial.log();

        // std::cout << "log strain elastic\n" << log_strain_e << std::endl;

        // Elastic stress predictor
        cauchy_stress = compute_cauchy_stress(material.shear_modulus(), material.lambda(), log_strain_e)
                        / J;

        // Trial von Mises stress
        von_mises = von_mises_stress(cauchy_stress);

        // std::cout << "von Mises " << von_mises << std::endl;

        // Compute the initial estimate of the yield function for the material
        // and decide if the stress return needs to be computed
        if (auto const f = evaluate_J2_yield_function(material, von_mises, accumulated_plastic_strain);
            f <= 0.0)
        {
            tangent_operators[l] = consistent_tangent(J, log_strain_e, cauchy_stress, C_e);
            continue;
        }

        auto const von_mises_trial = von_mises;

        std::cout << "\nQuadrature point plastic\n";

        // Compute the normal direction to the yield surface which remains
        // constant throughout the radial return method
        matrix2 const normal = deviatoric(cauchy_stress) / deviatoric(cauchy_stress).norm();

        // Initialise the plastic increment
        auto const plastic_increment = perform_radial_return(von_mises, accumulated_plastic_strain);

        // Plastic strain update
        // F_p = (plastic_increment * std::sqrt(3.0 / 2.0) * normal);

        // std::cout << ">>>Elastic strain before dec\n" << log_strain_e << std::endl;
        log_strain_e -= plastic_increment * std::sqrt(3.0 / 2.0) * normal;
        // std::cout << ">>>Elastic strain after dec\n" << log_strain_e << std::endl;

        cauchy_stress -= 2.0 * shear_modulus * plastic_increment * std::sqrt(3.0 / 2.0) * normal / J;

        von_mises = von_mises_stress(cauchy_stress);

        accumulated_plastic_strain += plastic_increment;

        // Use the elastoplastic tangent from infinitesimal strain theory
        matrix3 const D_ep = algorithmic_tangent(material.shear_modulus(),
                                                 material.hardening_modulus(accumulated_plastic_strain),
                                                 plastic_increment,
                                                 von_mises_trial,
                                                 normal,
                                                 I_dev,
                                                 C_e);

        // Compute the elastic-plastic tangent modulus for large strain
        tangent_operators[l] = consistent_tangent(J, log_strain_e, cauchy_stress, D_ep);
    }
}

matrix3 finite_strain_J2_plasticity::consistent_tangent(double const J,
                                                        matrix2 const& Be_trial,
                                                        matrix2 const& cauchy_stress,
                                                        matrix3 const& C)
{
    // Convert to Mandel notation so matrix multiplication == double dot operation
    matrix3 const D = mandel_notation(C);
    matrix3 const L = mandel_notation(log_symmetric_tensor_derivative(Be_trial));
    matrix3 const B = mandel_notation(finite_strain_B_operator(Be_trial));

    return 1.0 / (2.0 * J) * D * L * B - finite_strain_correction(cauchy_stress);
}
}
