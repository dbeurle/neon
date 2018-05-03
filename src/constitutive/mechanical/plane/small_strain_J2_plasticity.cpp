
#include "constitutive/mechanical/plane/small_strain_J2_plasticity.hpp"

#include "Exceptions.hpp"

#include "constitutive/internal_variables.hpp"
#include "constitutive/mechanical/detail/J2_plasticity.hpp"
#include "numeric/mechanics"

#include <range/v3/view/transform.hpp>

#include <tbb/parallel_for.h>

namespace neon::mechanical::plane
{
small_strain_J2_plasticity::small_strain_J2_plasticity(std::shared_ptr<internal_variables_t>& variables,
                                                       json const& material_data)
    : isotropic_linear_elasticity(variables, material_data, isotropic_linear_elasticity::plane::strain),
      material(material_data)
{
    variables->add(internal_variables_t::second::LinearisedPlasticStrain);
    variables->add(internal_variables_t::scalar::EffectivePlasticStrain);

    variables->commit();
}

small_strain_J2_plasticity::~small_strain_J2_plasticity() = default;

void small_strain_J2_plasticity::update_internal_variables(double const time_step_size)
{
    auto const shear_modulus = material.shear_modulus();

    // Extract the internal variables
    auto [plastic_strains,
          strains,
          cauchy_stresses] = variables->fetch(internal_variables_t::second::LinearisedPlasticStrain,
                                              internal_variables_t::second::LinearisedStrain,
                                              internal_variables_t::second::Cauchy);

    // Retrieve the accumulated internal variables
    auto [accumulated_plastic_strains,
          von_mises_stresses] = variables->fetch(internal_variables_t::scalar::EffectivePlasticStrain,
                                                 internal_variables_t::scalar::VonMisesStress);

    auto& tangent_operators = variables->fetch(internal_variables_t::fourth::tangent_operator);

    // Compute the linear strain gradient from the displacement gradient
    strains = variables->fetch(internal_variables_t::second::DisplacementGradient)
              | ranges::view::transform([](auto const& H) { return 0.5 * (H + H.transpose()); });

    // Perform the update algorithm for each quadrature point
    tbb::parallel_for(std::size_t{0}, strains.size(), [&](auto const l) {
        auto const& strain = strains[l];
        auto& plastic_strain = plastic_strains[l];
        auto& cauchy_stress = cauchy_stresses[l];
        auto& accumulated_plastic_strain = accumulated_plastic_strains[l];
        auto& von_mises = von_mises_stresses[l];

        // Elastic stress predictor
        cauchy_stress = compute_cauchy_stress(material.shear_modulus(),
                                              material.lambda(),
                                              strain - plastic_strain);

        // Trial von Mises stress
        von_mises = von_mises_stress(cauchy_stress);

        // If this quadrature point is elastic, then set the tangent to the
        // elastic modulus and continue to the next quadrature point
        if (evaluate_J2_yield_function(material, von_mises, accumulated_plastic_strain) <= 0.0)
        {
            tangent_operators[l] = C_e;
            return;
        }

        auto const von_mises_trial = von_mises;

        // Compute the normal direction to the yield surface which remains
        // constant throughout the radial return method
        matrix2 const normal = deviatoric(cauchy_stress) / deviatoric(cauchy_stress).norm();

        auto const plastic_increment = perform_radial_return(von_mises, accumulated_plastic_strain);

        plastic_strain += plastic_increment * std::sqrt(3.0 / 2.0) * normal;

        cauchy_stress -= 2.0 * shear_modulus * plastic_increment * std::sqrt(3.0 / 2.0) * normal;

        von_mises = von_mises_stress(cauchy_stress);

        accumulated_plastic_strain += plastic_increment;

        tangent_operators[l] = algorithmic_tangent(material.shear_modulus(),
                                                   material.hardening_modulus(
                                                       accumulated_plastic_strain),
                                                   plastic_increment,
                                                   von_mises_trial,
                                                   normal,
                                                   I_dev,
                                                   C_e);
    });
}

double small_strain_J2_plasticity::perform_radial_return(double const von_mises,
                                                         double const accumulated_plastic_strain) const
{
    auto const shear_modulus = material.shear_modulus();

    auto plastic_increment{0.0};

    auto f = evaluate_J2_yield_function(material, von_mises, accumulated_plastic_strain);

    // Perform the non-linear hardening solve
    int iterations{0};
    auto constexpr max_iterations{50};
    while (f > 1.0e-6 && iterations < max_iterations)
    {
        auto const H = material.hardening_modulus(accumulated_plastic_strain + plastic_increment);

        auto const plastic_increment_delta = f / (3.0 * shear_modulus + H);

        plastic_increment += plastic_increment_delta;

        f = evaluate_J2_yield_function(material,
                                       von_mises,
                                       accumulated_plastic_strain,
                                       plastic_increment);

        iterations++;
    }
    if (iterations == max_iterations)
    {
        std::cout << "\n";
        std::cout << std::string(8, ' ') << "Plastic increment : " << plastic_increment << "\n";

        std::cout << std::string(8, ' ')
                  << "Accumulated plastic strain : " << accumulated_plastic_strain << "\n";

        std::cout << std::string(8, ' ') << "Hardening modulus : "
                  << material.hardening_modulus(accumulated_plastic_strain + plastic_increment)
                  << "\n";

        std::cout << std::string(8, ' ') << "Shear modulus : " << shear_modulus << "\n";

        std::cout << std::string(8, ' ') << "Yield function after mapping : " << f << "\n";

        std::cout << std::string(8, ' ')
                  << "Current yield stress : " << material.yield_stress(accumulated_plastic_strain)
                  << "\n";

        throw computational_error("Non-convergence in radial return method.");
    }
    return plastic_increment;
}
}
