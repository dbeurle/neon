
#include "HyperElasticPlastic.hpp"

#include "Exceptions.hpp"
#include "InternalVariables.hpp"

#include <iostream>
#include <json/value.h>
#include <range/v3/view.hpp>
#include <unsupported/Eigen/MatrixFunctions>

namespace neon
{
FiniteJ2Plasticity::FiniteJ2Plasticity(InternalVariables& variables, Json::Value const& material_data)
    : HyperElasticPlastic(variables), material(material_data), C_e(elastic_moduli())
{
    variables.add(InternalVariables::Scalar::VonMisesStress,
                  InternalVariables::Scalar::EffectivePlasticStrain);

    variables.add(InternalVariables::Tensor::HenckyStrainElastic,
                  InternalVariables::Tensor::DeformationGradientPlastic);

    // Add material tangent with the linear elasticity moduli
    variables.add(InternalVariables::Matrix::TruesdellModuli, elastic_moduli());

    std::cout << "Constructed finite strain J2 plasticity model\n";
}

FiniteJ2Plasticity::~FiniteJ2Plasticity() = default;

void FiniteJ2Plasticity::update_internal_variables(double const time_step_size)
{
    using namespace ranges;

    auto const shear_modulus = material.shear_modulus();
    auto const lambda_e = material.lambda();

    // Extract the internal variables
    auto[F_list, strain_e_list, cauchy_stress_list] = variables(InternalVariables::Tensor::DeformationGradient,
                                                                InternalVariables::Tensor::HenckyStrainElastic,
                                                                InternalVariables::Tensor::Cauchy);

    auto const F_old_list = variables[InternalVariables::Tensor::DeformationGradient];

    auto const F_inc_list = view::zip(F_list, F_old_list) | view::transform([](auto const& tpl) {
                                auto const & [ F, F_old ] = tpl;
                                // std::cout << "F_n+1\n"
                                //           << F << "\n\nF_n\n"
                                //           << F_old << "\n\ndF\n"
                                //           << F * F_old.inverse() << "\n";
                                return F * F_old.inverse();
                            });

    // Retrieve the accumulated internal variables
    auto[accumulated_plastic_strain_list,
         von_mises_list] = variables(InternalVariables::Scalar::EffectivePlasticStrain,
                                     InternalVariables::Scalar::VonMisesStress);

    auto& C_list = variables(InternalVariables::Matrix::TruesdellModuli);

    // Perform the update algorithm for each quadrature point
    for (auto l = 0; l < F_list.size(); l++)
    {
        auto const& F_inc = F_inc_list[l];

        auto& cauchy_stress = cauchy_stress_list[l];
        auto& accumulated_plastic_strain = accumulated_plastic_strain_list[l];
        auto& von_mises = von_mises_list[l];
        auto& strain_e = strain_e_list[l];

        // Elastic trial deformation gradient
        Matrix3 const B_e = (2.0 * strain_e).exp();

        if (l == 1) std::cout << "Elastic left Cauchy-Green\n" << B_e << std::endl;

        // Elastic trial left Cauchy-Green deformation tensor
        Matrix3 const B_e_trial = F_inc * B_e * F_inc.transpose();

        // Trial Logarithmic elastic strain
        strain_e = 1.0 / 2.0 * B_e_trial.log();

        if (l == 1) std::cout << "Logarithmic trial strain\n" << strain_e << std::endl;

        // Elastic stress predictor
        Matrix3 const cauchy_stress_0 = lambda_e * strain_e.trace() * Matrix3::Identity()
                                        + 2.0 * shear_modulus * strain_e;

        if (l == 1) std::cout << "Trial stress\n" << cauchy_stress_0 << std::endl;

        // Compute the von Mises stress
        von_mises = von_mises_stress(cauchy_stress_0);

        if (l == 1) std::cout << "Von Mises = " << von_mises / 1.0e6 << std::endl;

        cauchy_stress = cauchy_stress_0;

        // Compute the initial estimate of the yield function for the material
        // and decide if the stress return needs to be computed
        auto f = von_mises - material.yield_stress(accumulated_plastic_strain);

        if (f <= 0.0)
        {
            C_list[l] = C_e;
            continue;
        }

        if (l == 1) std::cout << "\nQUADRATURE POINT PLASTIC\n";

        // Compute the normal direction to the yield surface which remains
        // constant throughout the radial return method
        Matrix3 const normal = deviatoric(cauchy_stress_0) / deviatoric(cauchy_stress_0).norm();

        // Initialize the plastic increment
        auto plastic_increment = 0.0;

        // Perform the return mapping algorithm
        int iterations = 0, max_iterations = 25;
        while (f > 1.0e-10 && iterations < max_iterations)
        {
            auto const H = material.hardening_modulus(accumulated_plastic_strain + plastic_increment);

            // Increment in plasticity rate
            const auto plastic_increment_delta = f / (3.0 * shear_modulus + H);

            // Plastic rate update
            plastic_increment += plastic_increment_delta;

            // Evaluate the yield function
            f = (von_mises - 3.0 * shear_modulus * plastic_increment)
                - material.yield_stress(accumulated_plastic_strain + plastic_increment);

            iterations++;
        }

        // Plastic strain update
        // F_p = (plastic_increment * std::sqrt(3.0 / 2.0) * normal);

        if (l == 1) std::cout << ">>>Elastic strain before dec\n" << strain_e << std::endl;

        strain_e -= plastic_increment * std::sqrt(3.0 / 2.0) * normal;

        if (l == 1) std::cout << ">>>Elastic strain after dec\n" << strain_e << std::endl;

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

            throw computational_error("Radial return failure at global integration point "
                                      + std::to_string(l) + "\n");
        }
        // Compute the elastic-plastic tangent modulus
        C_list[l] = algorithmic_tangent(plastic_increment,
                                        accumulated_plastic_strain,
                                        von_mises_stress(cauchy_stress_0),
                                        normal);
    }
}

CMatrix FiniteJ2Plasticity::elastic_moduli() const
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

CMatrix FiniteJ2Plasticity::deviatoric_projection() const
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

CMatrix FiniteJ2Plasticity::algorithmic_tangent(double const plastic_increment,
                                                double const accumulated_plastic_strain,
                                                double const von_mises,
                                                Matrix3 const& n) const
{
    auto const G = material.shear_modulus();
    auto const H = material.hardening_modulus(accumulated_plastic_strain);
    auto const n_outer_n = voigt(n) * voigt(n).transpose();
    return C_e - 6.0 * std::pow(G, 2) * plastic_increment / (3.0 * G + H) * n_outer_n
           - 4.0 * std::pow(G, 2) * plastic_increment / von_mises * std::sqrt(3.0 / 2.0)
                 * (I_dev - n_outer_n);
}
}
