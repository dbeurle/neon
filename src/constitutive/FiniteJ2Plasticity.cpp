
#include "FiniteJ2Plasticity.hpp"

#include "Exceptions.hpp"
#include "InternalVariables.hpp"

#include <iostream>
#include <json/value.h>
#include <range/v3/view.hpp>
#include <unsupported/Eigen/MatrixFunctions>

namespace neon
{
FiniteJ2Plasticity::FiniteJ2Plasticity(InternalVariables& variables, Json::Value const& material_data)
    : J2Plasticity(variables, material_data)
{
    variables.add(InternalVariables::Scalar::VonMisesStress,
                  InternalVariables::Scalar::EffectivePlasticStrain);

    variables.add(InternalVariables::Tensor::HenckyStrainElastic);

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

    auto const J_list = variables(InternalVariables::Scalar::DetF);

    // Retrieve the accumulated internal variables
    auto[accumulated_plastic_strain_list,
         von_mises_list] = variables(InternalVariables::Scalar::EffectivePlasticStrain,
                                     InternalVariables::Scalar::VonMisesStress);

    auto& C_list = variables(InternalVariables::Matrix::TruesdellModuli);

    auto const F_inc_list = view::zip(F_list, F_old_list) | view::transform([](auto const& tpl) {
                                auto const & [ F, F_old ] = tpl;
                                return F * F_old.inverse();
                            });

    // Perform the update algorithm for each quadrature point
    for (auto l = 0; l < F_list.size(); l++)
    {
        auto const& F_inc = F_inc_list[l];
        auto const& F = F_list[l];
        auto const J = J_list[l];

        auto& cauchy_stress = cauchy_stress_list[l];
        auto& accumulated_plastic_strain = accumulated_plastic_strain_list[l];
        auto& von_mises = von_mises_list[l];
        auto& strain_e = strain_e_list[l];

        // Elastic trial deformation gradient
        Matrix3 const B_e = (2.0 * strain_e).exp();

        // Elastic trial left Cauchy-Green deformation tensor
        Matrix3 const B_e_trial = F_inc * B_e * F_inc.transpose();

        // Trial Logarithmic elastic strain
        strain_e = 0.5 * B_e_trial.log();

        // Elastic stress predictor
        Matrix3 const cauchy_stress_0 = compute_cauchy_stress(strain_e) / J;

        von_mises = von_mises_stress(cauchy_stress_0);

        cauchy_stress = cauchy_stress_0;

        // Compute the initial estimate of the yield function for the material
        // and decide if the stress return needs to be computed
        auto f = von_mises - material.yield_stress(accumulated_plastic_strain);

        std::cout << "f = " << f << std::endl;

        if (f <= 0.0)
        {
            C_list[l] = C_e;
            continue;
        }

        std::cout << "\nQUADRATURE POINT PLASTIC\n";

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

        // Plastic strain update
        // F_p = (plastic_increment * std::sqrt(3.0 / 2.0) * normal);

        std::cout << ">>>Elastic strain before dec\n" << strain_e << std::endl;

        strain_e -= plastic_increment * std::sqrt(3.0 / 2.0) * normal;

        std::cout << ">>>Elastic strain after dec\n" << strain_e << std::endl;

        cauchy_stress -= 2.0 * shear_modulus * plastic_increment * std::sqrt(3.0 / 2.0) * normal / J;

        von_mises = von_mises_stress(cauchy_stress);

        accumulated_plastic_strain += plastic_increment;

        // Compute the elastic-plastic tangent modulus
        C_list[l] = algorithmic_tangent(plastic_increment,
                                        accumulated_plastic_strain,
                                        von_mises_stress(cauchy_stress_0),
                                        normal);
    }
}

std::tuple<Vector3, std::array<Matrix3, 3>, int> FiniteJ2Plasticity::compute_eigenvalues_eigenprojections(
    Matrix3 const& X) const
{
    Eigen::SelfAdjointEigenSolver<Matrix3> eigen_solver(X);

    if (eigen_solver.info() != Eigen::Success)
    {
        throw computational_error("Eigenvalue solver failed in finite plasticity routine\n");
    }

    Vector3 const x = eigen_solver.eigenvalues();
    Matrix3 const v = eigen_solver.eigenvectors();

    // Eigenprojections
    std::array<Matrix3, 3> E = {{v.col(0) * v.col(0).transpose(),
                                 v.col(1) * v.col(1).transpose(),
                                 v.col(2) * v.col(2).transpose()}};

    int unique_eigenvalues = 3;

    // Need to check if there are repeated eigenvalues to machine precision
    if (is_approx(x(0), x(1)) && is_approx(x(1), x(2)))
    {
        unique_eigenvalues = 1;
        E[0] = Matrix3::Identity();
    }
    else if (!is_approx(x(0), x(1)) && is_approx(x(1), x(2)))
    {
        unique_eigenvalues = 2;
        E[1] = Matrix3::Identity() - E[0];
    }
    else if (!is_approx(x(2), x(0)) && is_approx(x(0), x(1)))
    {
        unique_eigenvalues = 2;
        E[0] = Matrix3::Identity() - E[2];
    }
    else if (!is_approx(x(1), x(2)) && is_approx(x(2), x(0)))
    {
        unique_eigenvalues = 2;
        E[2] = Matrix3::Identity() - E[1];
    }
    return {x, E, unique_eigenvalues};

    // auto const I1 = X.trace();
    // auto const I2 = 0.5 * (std::pow(X.trace(), 2) - (X * X).trace());
    // auto const I3 = X.determinant();
    //
    // auto const R = (-2.0 * I1 + 9.0 * I1 * I2 - 27.0 * I3) / 54.0;
    // auto const Q = (std::pow(I1, 2) - 3.0 * I2) / 9.0;
    //
    // auto const theta = std::acos(R / std::sqrt(std::pow(Q, 3)));
    //
    // Vector3 x(-2.0 * std::sqrt(Q) * std::cos(theta / 3.0) + I1 / 3.0,
    //           -2.0 * std::sqrt(Q) * std::cos((theta + 2.0 * M_PI) / 3.0) + I1 / 3.0,
    //           -2.0 * std::sqrt(Q) * std::cos((theta - 2.0 * M_PI) / 3.0) + I1 / 3.0);
    //
    // std::array<Matrix3, 3> E;
    // int repeated_eigenvalues = 0;
    //
    // // Perform checks to determine if there are repeated eigenvalues
    // if (!is_approx(x1, x2) && !is_approx(x2, x3) && !is_approx(x1, x3))
    // {
    //     for (int i = 0; i < 3; i++)
    //     {
    //         E[i] = x(i) / (2.0 * std::pow(x(i), 3) - I1 * std::pow(x(i), 2) + I3)
    //                * (X * X - (I1 - x(i)) * X + I3 / x(i) * Matrix3::Identity());
    //     }
    // }
    // else if (!is_approx(x1, x2) && is_approx(x2, x3))
    // {
    //     repeated_eigenvalue = 1;
    // }
    // else
    // {
    //     E[0] = Matrix3::Identity();
    //     repeated_eigenvalue = 2;
    // }
}

CMatrix FiniteJ2Plasticity::derivative_tensor_log(Matrix3 const& Be_trial) const
{
    auto const[e, E_list, unique_eigenvalues] = compute_eigenvalues_eigenprojections(Be_trial);

    CMatrix L = CMatrix::Zero(6, 6);

    if (unique_eigenvalues == 3)
    {
        for (int i = 0; i < 3; i++)
        {
            L.noalias() += Isym;
        }
    }
    return L;
}
}
