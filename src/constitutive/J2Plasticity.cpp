
#include "J2Plasticity.hpp"

#include "Exceptions.hpp"
#include "InternalVariables.hpp"

#include <iostream>
#include <json/value.h>
#include <range/v3/view.hpp>

namespace neon
{
J2Plasticity::J2Plasticity(InternalVariables& variables, Json::Value const& material_data)
    : HypoElasticPlastic(variables), material(material_data)
{
    variables.add(InternalVariables::Tensor::LinearisedStrain,
                  InternalVariables::Tensor::LinearisedPlasticStrain);

    variables.add(InternalVariables::Scalar::VonMisesStress,
                  InternalVariables::Scalar::EffectivePlasticStrain);

    // Add material tangent with the linear elasticity moduli
    variables.add(InternalVariables::Matrix::TruesdellModuli, elastic_moduli());
}

J2Plasticity::~J2Plasticity() = default;

void J2Plasticity::update_internal_variables(double const time_step_size)
{
    using namespace ranges;

    auto const shear_modulus = material.shear_modulus();

    // Extract the internal variables
    auto[strain_p_list,
         strain_list,
         cauchy_stress_list] = variables(InternalVariables::Tensor::LinearisedPlasticStrain,
                                         InternalVariables::Tensor::LinearisedStrain,
                                         InternalVariables::Tensor::Cauchy);

    // Retrieve the accumulated internal variables
    auto[accumulated_plastic_strain_list,
         von_mises_list] = variables(InternalVariables::Scalar::EffectivePlasticStrain,
                                     InternalVariables::Scalar::VonMisesStress);

    auto& C_list = variables(InternalVariables::Matrix::TruesdellModuli);

    // Compute the linear strain gradient from the displacement gradient
    strain_list = variables(InternalVariables::Tensor::DisplacementGradient)
                  | view::transform([](auto const& H) { return 0.5 * (H + H.transpose()); });

    // Perform the update algorithm for each quadrature point
    for (auto l = 0; l < strain_list.size(); l++)
    {
        auto const& strain = strain_list[l];
        auto& strain_p = strain_p_list[l];
        auto& cauchy_stress = cauchy_stress_list[l];
        auto& accumulated_plastic_strain = accumulated_plastic_strain_list[l];
        auto& von_mises = von_mises_list[l];

        // Elastic stress predictor
        cauchy_stress = compute_cauchy_stress(strain - strain_p);

        // Trial von Mises stress
        von_mises = von_mises_stress(cauchy_stress);

        // If this quadrature point is elastic, then set the tangent to the
        // elastic modulus and continue to the next quadrature point
        if (evaluate_yield_function(von_mises, accumulated_plastic_strain) <= 0.0)
        {
            C_list[l] = C_e;
            continue;
        }

        auto const von_mises_trial = von_mises;

        // Compute the normal direction to the yield surface which remains
        // constant throughout the radial return method
        Matrix3 const normal = deviatoric(cauchy_stress) / deviatoric(cauchy_stress).norm();

        auto const plastic_increment = perform_radial_return(von_mises, accumulated_plastic_strain);

        strain_p += plastic_increment * std::sqrt(3.0 / 2.0) * normal;

        cauchy_stress -= 2.0 * shear_modulus * plastic_increment * std::sqrt(3.0 / 2.0) * normal;

        von_mises = von_mises_stress(cauchy_stress);

        accumulated_plastic_strain += plastic_increment;

        C_list[l] = algorithmic_tangent(plastic_increment,
                                        accumulated_plastic_strain,
                                        von_mises_trial,
                                        normal);
    }
}

CMatrix J2Plasticity::elastic_moduli() const
{
    auto const[lambda, shear_modulus] = material.Lame_parameters();
    CMatrix C(6, 6);
    C << lambda + 2.0 * shear_modulus, lambda, lambda, 0.0, 0.0, 0.0, //
        lambda, lambda + 2.0 * shear_modulus, lambda, 0.0, 0.0, 0.0,  //
        lambda, lambda, lambda + 2.0 * shear_modulus, 0.0, 0.0, 0.0,  //
        0.0, 0.0, 0.0, shear_modulus, 0.0, 0.0,                       //
        0.0, 0.0, 0.0, 0.0, shear_modulus, 0.0,                       //
        0.0, 0.0, 0.0, 0.0, 0.0, shear_modulus;
    return C;
}

CMatrix J2Plasticity::deviatoric_projection() const
{
    CMatrix C(6, 6);
    C << 2.0 / 3.0, -1.0 / 3.0, -1.0 / 3.0, 0.0, 0.0, 0.0, //
        -1.0 / 3.0, 2.0 / 3.0, -1.0 / 3.0, 0.0, 0.0, 0.0,  //
        -1.0 / 3.0, -1.0 / 3.0, 2.0 / 3.0, 0.0, 0.0, 0.0,  //
        0.0, 0.0, 0.0, 0.5, 0.0, 0.0,                      //
        0.0, 0.0, 0.0, 0.0, 0.5, 0.0,                      //
        0.0, 0.0, 0.0, 0.0, 0.0, 0.5;
    return C;
}

CMatrix J2Plasticity::algorithmic_tangent(double const plastic_increment,
                                          double const accumulated_plastic_strain,
                                          double const von_mises,
                                          Matrix3 const& n) const
{
    auto const G = material.shear_modulus();
    auto const H = material.hardening_modulus(accumulated_plastic_strain);

    return C_e - plastic_increment * 6.0 * std::pow(G, 2) / von_mises * I_dev
           + 6.0 * std::pow(G, 2) * (plastic_increment / von_mises - 1.0 / (3.0 * G + H)) * voigt(n)
                 * voigt(n).transpose();
}

double J2Plasticity::perform_radial_return(double const von_mises,
                                           double const accumulated_plastic_strain) const
{
    auto const shear_modulus = material.shear_modulus();

    auto plastic_increment = 0.0;

    double f = evaluate_yield_function(von_mises, accumulated_plastic_strain, plastic_increment);

    // Perform the non-linear hardening solve
    int iterations = 0, max_iterations = 25;
    while (f > 1.0e-10 && iterations < max_iterations)
    {
        auto const H = material.hardening_modulus(accumulated_plastic_strain + plastic_increment);

        const auto plastic_increment_delta = f / (3.0 * shear_modulus + H);

        plastic_increment += plastic_increment_delta;

        f = evaluate_yield_function(von_mises, accumulated_plastic_strain, plastic_increment);

        iterations++;
    }
    if (iterations == max_iterations)
    {
        std::cout << "MAXIMUM NUMBER OF ITERATIONS IN RADIAL RETURN REACHED\n";
        std::cout << "Accumulated plastic strain " << accumulated_plastic_strain << "\n";
        std::cout << "Yield function after mapping " << f << "\n";
        std::cout << "Current yield stress " << material.yield_stress(accumulated_plastic_strain)
                  << "\n";
        throw computational_error("Radial return failure\n");
    }
    return plastic_increment;
}

double J2Plasticity::evaluate_yield_function(double const von_mises,
                                             double const accumulated_plastic_strain,
                                             double const plastic_increment) const
{
    auto const shear_modulus = material.shear_modulus();

    return (von_mises - 3.0 * shear_modulus * plastic_increment)
           - material.yield_stress(accumulated_plastic_strain + plastic_increment);
}
}
