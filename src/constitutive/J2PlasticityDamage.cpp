#include "J2PlasticityDamage.hpp"

#include "Exceptions.hpp"
#include "InternalVariables.hpp"

#include <iostream>
#include <json/value.h>
#include <range/v3/view.hpp>

namespace neon::mech::solid
{
J2PlasticityDamage::J2PlasticityDamage(InternalVariables& variables, Json::Value const& material_data)
    : J2Plasticity(variables, material_data), material(material_data)
{
    variables.add(InternalVariables::Tensor::BackStress,
                  InternalVariables::Tensor::KinematicHardening,
                  InternalVariables::Scalar::Damage,
                  InternalVariables::Scalar::EnergyReleaseRate);
}

J2PlasticityDamage::~J2PlasticityDamage() = default;

void J2PlasticityDamage::update_internal_variables(double const time_step_size)
{
    using namespace ranges;

    auto const shear_modulus = material.shear_modulus();

    // Extract the internal variables
    auto [strain_p_list,
          strain_list,
          cauchy_stress_list,
          back_stress_list,
          kin_hard_list] = variables(InternalVariables::Tensor::LinearisedPlasticStrain,
                                     InternalVariables::Tensor::LinearisedStrain,
                                     InternalVariables::Tensor::Cauchy,
                                     InternalVariables::Tensor::BackStress,
                                     InternalVariables::Tensor::KinematicHardening);

    // Retrieve the accumulated internal variables
    auto [accumulated_plastic_strain_list,
          von_mises_list,
          damage_list,
          energy_release_list] = variables(InternalVariables::Scalar::EffectivePlasticStrain,
                                           InternalVariables::Scalar::VonMisesStress,
                                           InternalVariables::Scalar::Damage,
                                           InternalVariables::Scalar::EnergyReleaseRate);

    auto& C_list = variables(InternalVariables::Matrix::TangentOperator);

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
        auto& back_stress = back_stress_list[l];
        auto& kin_hard = kin_hard_list[l];
        auto& damage_var = damage_list[l];
        auto& energy_var = energy_release_list[l];
        auto& C_algorithmic = C_list[l];
        // Elastic stress predictor
        cauchy_stress = compute_cauchy_stress(strain - strain_p);

        Matrix3 tau = deviatoric(cauchy_stress) / (1.0 - damage_var) - deviatoric(back_stress);
        // Trial von Mises stress
        von_mises = von_mises_stress(tau);

        energy_var = 0.5
                     * double_dot(strain - strain_p,
                                  compute_stress_like_matrix((1.0 - damage_var) * C_e,
                                                             strain - strain_p));
        // If this quadrature point is elastic, then set the tangent to the
        // elastic modulus and continue to the next quadrature point
        if (evaluate_yield_function(von_mises, back_stress) <= 0.0
            && evaluate_damage_yield_function(energy_var) <= 0.0)
        {
            C_list[l] = C_e;
            continue;
        }
        auto const plastic_increment = perform_radial_return(cauchy_stress,
                                                             back_stress,
                                                             damage_var,
                                                             kin_hard,
                                                             energy_var,
                                                             C_algorithmic,
                                                             time_step_size,
                                                             strain - strain_p);
        tau = deviatoric(cauchy_stress) / (1.0 - damage_var) - deviatoric(back_stress);
        von_mises = von_mises_stress(tau);
        strain_p += plastic_increment * 3.0 / 2.0 * tau / (von_mises * (1 - damage_var));
        accumulated_plastic_strain += plastic_increment / (1 - damage_var);
        // cauchy_stress, back_stress, damage_var, kin_hard, energy_var, C_algorithmic, are updated
        // within the radial_return routine std::cout << "delta_t  " << time_step_size << "\n";
    }
}

double J2PlasticityDamage::perform_radial_return(Matrix3& cauchy_stress,
                                                 Matrix3& back_stress,
                                                 double& damage_var,
                                                 Matrix3& kin_hard,
                                                 double& energy_var,
                                                 Matrix6& C_list,
                                                 double const& delta_t,
                                                 Matrix3 const& eps_e_t)
{
    // TODO: initial guess could be computed base on a frozen yield surface (at
    // the increment n) to obtain good convergence. check computational methods for plasticity book
    Vector y(16, 1);
    y << 0.0, 0.0, voigt::kinetic::to(cauchy_stress), voigt::kinetic::to(back_stress), damage_var,
        energy_var;
    // The residual
    Vector f = Vector::Identity(16, 1);
    int iterations = 0, max_iterations = 25;
    Matrix M = Matrix::Zero(16, 16);
    // auto identity_row = Vector::Ones(1, 6);
    // auto identity_col = Vector::Ones(6, 1);
    M(0, 0) = M(1, 1) = M(14, 14) = M(15, 15) = 1.0;
    M.block<6, 6>(2, 2) = voigt::kinetic::fourth_order_identity(); // TODO: double check kinetic vs
                                                                   // kinematic

    auto const kp = material.plasticity_viscous_multiplier();
    auto const np = material.plasticity_viscous_exponent();
    auto const kd = material.damage_viscous_multiplier();
    auto const nd = material.damage_viscous_exponent();
    auto const C = material.hardening_modulus();
    auto const gamma = material.softening_multiplier();

    // TODO: use relative error instead of absolute one (for f and delta_y)
    while (f.norm() > 1.0e-6 && iterations < max_iterations) // TODO: use ResidualTolerance
    {
        auto const d_lam_plastic = y(0);
        auto const d_lam_damage = y(1);
        auto const sigma = voigt::kinetic::from(y.segment<6>(2));
        auto const beta = voigt::kinetic::from(y.segment<6>(8));
        auto const d = y(14);
        auto const Y = y(15);
        Matrix3 const tau = deviatoric(sigma) / (1.0 - d) - deviatoric(beta);
        auto const von_mises = von_mises_stress(tau);
        auto const f_vp = std::max(0.0, evaluate_yield_function(von_mises, beta));
        auto const f_d = std::max(0.0, evaluate_damage_yield_function(Y));
        Matrix3 const normal_tild = 3.0 / 2.0 * tau / von_mises;
        Matrix3 const normal = normal_tild / (1.0 - d);

        // The residual TODO: compute_stress_like_vector and compute_stress_like_matrix could be
        // optimised
        // clang-format off
        f << d_lam_plastic - kp * std::pow(f_vp, np),
             d_lam_damage - kd * std::pow(f_d, nd),
             voigt::kinetic::to(sigma) - (1.0 - d) * compute_stress_like_vector(C_e, eps_e_t)
                + delta_t * d_lam_plastic * compute_stress_like_vector(C_e, normal_tild),
             voigt::kinetic::to(beta - C * kin_hard
                - C * delta_t * d_lam_plastic * (normal_tild - gamma / C * beta)),
             d - damage_var - delta_t * d_lam_damage,
             Y - 0.5   * double_dot(eps_e_t - delta_t * d_lam_plastic * normal_tild / (1.0 - d),
    compute_stress_like_matrix(C_e, eps_e_t - delta_t * d_lam_plastic * normal_tild / (1.0 - d)));
        // clang-format on

        M.block<1, 6>(0, 2) = voigt::kinetic::to(
            deviatoric(Matrix3(-np * kp * std::pow(f_vp, np - 1) * normal)));
        M.block<1, 6>(0, 8) = voigt::kinetic::to(np * kp * std::pow(f_vp, np - 1)
                                                 * (normal_tild - gamma / C * beta));
        M(0, 14) = -np * kp * std::pow(f_vp, np - 1) * double_dot(normal_tild, deviatoric(sigma))
                   / std::pow(1.0 - d, 2);
        M(1, 15) = -nd * kd * std::pow(f_d, nd - 1);
        M.block<6, 1>(2, 0) = delta_t * compute_stress_like_vector(C_e, normal_tild);
        M.block<6, 1>(2, 14) = compute_stress_like_vector(C_e, eps_e_t);
        M.block<6, 1>(8, 0) = voigt::kinetic::to(-delta_t * C * (normal_tild - gamma / C * beta));
        M.block<6, 6>(8, 8) = voigt::kinetic::fourth_order_identity()
                              * (1.0 + delta_t * d_lam_plastic * gamma);

        M(14, 1) = -delta_t;
        M(15, 0) = double_dot(delta_t * normal_tild / (1.0 - d),
                              compute_stress_like_matrix(C_e,
                                                         eps_e_t
                                                             - delta_t * d_lam_plastic * normal_tild
                                                                   / (1.0 - d)));
        M(15, 14) = delta_t * d_lam_plastic / std::pow(1.0 - d, 2)
                    * double_dot(normal_tild,
                                 compute_stress_like_matrix(C_e,
                                                            eps_e_t
                                                                - delta_t * d_lam_plastic
                                                                      * normal_tild / (1.0 - d)));

        // Eigen::JacobiSVD<Eigen::MatrixXd> svd(M);
        // double cond = svd.singularValues()(0) /
        // svd.singularValues()(svd.singularValues().size()
        // - 1); std::cout << cond << "\n\n";
        // Vector a = -M.inverse() * f;
        Vector a = -M.fullPivLu().solve(f);
        y = y + a;

        iterations++;
    }

    if (iterations == max_iterations)
    {
        std::cout << "MAXIMUM NUMBER OF ITERATIONS IN RADIAL RETURN REACHED\n";
        std::cout << "The residual norm " << f.norm() << "\n";
        // std::cout << "The increment norm " << a.norm() << "\n";
        throw computational_error("Radial return failure\n");
    }
    cauchy_stress = voigt::kinetic::from(y.segment<6>(2));
    back_stress = voigt::kinetic::from(y.segment<6>(8));
    damage_var = y(14);
    energy_var = y(15);
    kin_hard = back_stress / C;
    // C_list = (1.0 - damage_var) * C_e;
    C_list = (1.0 - damage_var) * mandel_notation(C_e)
             * mandel_notation(M.inverse().block<6, 6>(2, 2));
    return y(0) * delta_t; // plastic_increment
}

double J2PlasticityDamage::evaluate_yield_function(double const von_mises,
                                                   Matrix3 const& back_stress) const
{
    return von_mises
           + 0.5 * material.softening_multiplier() / material.hardening_modulus()
                 * double_dot(back_stress, back_stress)
           - material.yield_stress(0.0);
}
double J2PlasticityDamage::evaluate_damage_yield_function(double const energy_var) const
{
    return energy_var - std::pow(material.yield_stress(0.0), 2) / (2 * material.elastic_modulus());
}

Matrix3 J2PlasticityDamage::compute_cauchy_stress(Matrix3 const& elastic_strain) const
{
    auto const G = material.shear_modulus();
    auto const lambda_e = material.lambda();
    return lambda_e * elastic_strain.trace() * Matrix3::Identity() + 2.0 * G * elastic_strain;
}
Matrix3 J2PlasticityDamage::compute_stress_like_matrix(Matrix6 const& tangent_operator,
                                                       Matrix3 const& strain_like) const
{
    // could get the tangent_operator directly without passing it to this function
    return voigt::kinetic::from(tangent_operator * voigt::kinematic::to(strain_like));
}
Vector6 J2PlasticityDamage::compute_stress_like_vector(Matrix6 const& tangent_operator,
                                                       Matrix3 const& strain_like) const
{
    // could get the tangent_operator directly without passing it to this function
    return tangent_operator * voigt::kinematic::to(strain_like);
}
} // namespace neon::mech::solid
