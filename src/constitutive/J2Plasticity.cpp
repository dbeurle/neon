
#include "J2Plasticity.hpp"

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
    auto const lambda_e = material.lambda();

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

    Matrix3 const I = Matrix3::Identity();

    // Perform the update algorithm for each quadrature point
    for (auto l = 0; l < strain_list.size(); l++)
    {
        auto const& strain = strain_list[l];
        auto& strain_p = strain_p_list[l];
        auto& cauchy_stress = cauchy_stress_list[l];
        auto& accumulated_plastic_strain = accumulated_plastic_strain_list[l];
        auto& von_mises = von_mises_list[l];

        // Elastic stress predictor
        Matrix3 const cauchy_stress_0 = lambda_e * (strain - strain_p).trace() * I
                                        + 2.0 * shear_modulus * (strain - strain_p);

        von_mises = von_mises_stress(cauchy_stress_0);

        cauchy_stress = cauchy_stress_0;

        // Compute the initial estimate of the yield function for the material
        // and decide if the stress return needs to be computed
        auto f = von_mises - material.yield_stress(accumulated_plastic_strain);

        // If this quadrature point is elastic, then set the tangent to the
        // elastic modulus
        if (f <= 0.0)
        {
            C_list[l] = C_e;
            continue;
        }

        // Compute the normal direction to the yield surface which remains
        // constant throughout the radial return method
        Matrix3 const normal = deviatoric(cauchy_stress_0) / deviatoric(cauchy_stress_0).norm();

        auto plastic_increment = 0.0;

        // Newton Raphson iterations for non-linear hardening
        int iterations = 0, max_iterations = 50;
        while (f > 1.0e-10 && iterations < max_iterations)
        {
            auto const H = material.hardening_modulus(accumulated_plastic_strain + plastic_increment);

            const auto plastic_increment_delta = f / (3.0 * shear_modulus + H);

            plastic_increment += plastic_increment_delta;

            f = (von_mises - 3.0 * shear_modulus * plastic_increment)
                - material.yield_stress(accumulated_plastic_strain + plastic_increment);

            iterations++;
        }

        // Return mapping algorithm
        strain_p += plastic_increment * std::sqrt(3.0 / 2.0) * normal;

        cauchy_stress -= 2.0 * shear_modulus * plastic_increment * std::sqrt(3.0 / 2.0) * normal;

        von_mises = von_mises_stress(cauchy_stress);

        accumulated_plastic_strain += plastic_increment;

        if (iterations == max_iterations)
        {
            std::cout << "MAXIMUM NUMBER OF ITERATIONS IN RADIAL RETURN REACHED\n";
            std::cout << "Accumulated plastic strain " << accumulated_plastic_strain << "\n";
            std::cout << "Yield function after mapping " << f << "\n";
            std::cout << "Current yield stress "
                      << material.yield_stress(accumulated_plastic_strain) << "\n";
            std::cout << "Von Mises stress " << von_mises_stress(cauchy_stress) << "\n"
                      << "\n";

            throw std::runtime_error("Radial return failure at global integration point "
                                     + std::to_string(l) + "\n");
        }

        C_list[l] = algorithmic_tangent(plastic_increment,
                                        accumulated_plastic_strain,
                                        von_mises_stress(cauchy_stress_0),
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

CMatrix J2Plasticity::incremental_tangent(double const plastic_increment, double const von_mises) const
{
    auto const G = material.shear_modulus();
    return C_e - plastic_increment * 6.0 * std::pow(G, 2) / von_mises * I_dev;
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
}
