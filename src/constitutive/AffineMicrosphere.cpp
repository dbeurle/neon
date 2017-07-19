
#include "AffineMicrosphere.hpp"

#include "InternalVariables.hpp"
#include "numeric/DenseTypes.hpp"

#include <json/json.h>
#include <range/v3/algorithm.hpp>
#include <range/v3/view.hpp>

#include <exception>
#include <omp.h>

namespace neon
{
AffineMicrosphere::AffineMicrosphere(InternalVariables& variables,
                                     Json::Value const& material_data)
    : Hyperelastic(variables), material(material_data)
{
    variables.add(InternalVariables::Matrix::TruesdellModuli, 6);

    // Deviatoric stress
    variables.add(InternalVariables::Tensor::Kirchhoff);
    variables.add(InternalVariables::Scalar::Chains);

    for (auto& n : variables(InternalVariables::Scalar::Chains))
    {
        n = material.number_of_chains();
    }
}

void AffineMicrosphere::update_internal_variables(double const Δt)
{
    using namespace ranges;

    // Get references into the hash table
    auto& F_list = variables(InternalVariables::Tensor::DeformationGradient);
    auto& σ_list = variables(InternalVariables::Tensor::Cauchy);
    auto& τ_list = variables(InternalVariables::Tensor::Kirchhoff);

    auto const& detF_list = variables(InternalVariables::Scalar::DetF);
    auto& n_list = variables(InternalVariables::Scalar::Chains);

    auto const N = material.segments_per_chain();

    // Update the number of chains in the network
    n_list = n_list | view::transform([&](auto const& n) {
                 return material.evolve_chains(n, Δt);
             });

#pragma omp parallel for
    for (auto l = 0; l < F_list.size(); ++l)
    {
        auto& τ = τ_list[l];
        auto const& F = F_list[l];
        auto const& J = detF_list[l];

        auto const μ = material.shear_modulus(n_list[l]);

        // Unimodular decomposition of F
        Matrix3 const unimodular_F = std::pow(J, -1.0 / 3.0) * F;

        // clang-format off
        τ = μ * unit_sphere.integrate(Matrix3::Zero(),
                                      [&](auto const& coordinates,
                                          auto const& l) -> Matrix3 {
                            auto const & [ r, r_outer_r ] = coordinates;

                            // Deformed tangents
                            Vector3 const t = unimodular_F * r;

                            // Microstretches
                            auto const λ = t.norm();

                            return pade_first(λ, N) * t * t.transpose();
                        });
        // clang-format on
    }

    // Perform the projection of the stresses
    σ_list = view::zip(τ_list, detF_list, n_list)
             | view::transform([&, this](auto const& tpl) -> Matrix3 {

                   auto const & [ τ_dev, J, n ] = tpl;

                   auto const μ = material.shear_modulus(n);

                   auto const pressure = J * volumetric_free_energy_derivative(J, μ);

                   return deviatoric_projection(pressure, τ_dev) / J;
               });

    // Compute tangent moduli
    auto& D_list = variables(InternalVariables::Matrix::TruesdellModuli);

#pragma omp parallel for
    for (auto l = 0; l < F_list.size(); ++l)
    {
        auto const& F = F_list[l];
        auto const& τ_dev = τ_list[l];
        auto const& J = detF_list[l];

        auto const μ = material.shear_modulus(n_list[l]);

        auto const pressure = J * volumetric_free_energy_derivative(J, μ);
        auto const κ = std::pow(J, 2) * volumetric_free_energy_second_derivative(J, μ);

        // Unimodular decomposition of F
        Matrix3 const unimodular_F = std::pow(J, -1.0 / 3.0) * F;

        // clang-format off
        auto D = unit_sphere.integrate(Matrix::Zero(6, 6),
                                       [&](auto const& coordinates, auto const& l) -> Matrix {
                      auto const & [ r, r_outer_r ] = coordinates;

                      // Deformed tangents
                      auto const t = unimodular_F * r;

                      // Microstretches
                      auto const λ = t.norm();

                      auto const a = std::pow(λ, -2) * (pade_second(λ, N) - pade_first(λ, N));

                      return a * voigt(t * t.transpose()) * voigt(t * t.transpose()).transpose();
                });
        // clang-format on
        D_list[l] = deviatoric_projection(D, τ_dev) + (κ + pressure) * IoI
                    - 2.0 * pressure * I;
    }
}

Matrix3 AffineMicrosphere::deviatoric_projection(double const pressure,
                                                 Matrix3 const& τ_dev) const
{
    Matrix3 P_double_dot_τ;
    P_double_dot_τ << 2 * τ_dev(0, 0) / 3.0 - τ_dev(1, 1) / 3.0 - τ_dev(2, 2) / 3.0, //
        τ_dev(0, 1),                                                                 //
        τ_dev(0, 2),                                                                 //

        τ_dev(0, 1),                                                    //
        -τ_dev(0, 0) / 3.0 + 2 * τ_dev(1, 1) / 3.0 - τ_dev(2, 2) / 3.0, //
        τ_dev(1, 2),                                                    //

        τ_dev(0, 2), //
        τ_dev(1, 2), //
        -τ_dev(0, 0) / 3.0 - τ_dev(1, 1) / 3.0 + 2 * τ_dev(2, 2) / 3.0;
    return pressure * Matrix3::Identity() + P_double_dot_τ;
}

Matrix AffineMicrosphere::deviatoric_projection(Matrix const& C_dev,
                                                Matrix3 const& τ_dev) const
{
    Matrix C(6, 6);
    C << 1.0 / 9.0
             * (4 * C_dev(0, 0) - 2 * C_dev(0, 1) - 2 * C_dev(0, 2) - 2 * C_dev(1, 0)
                + C_dev(1, 1) + C_dev(1, 2) - 2 * C_dev(2, 0) + C_dev(2, 1) + C_dev(2, 2)
                + 4 * τ_dev.trace()), //
        1.0 / 9.0
            * (-2 * C_dev(0, 0) + 4 * C_dev(0, 1) - 2 * C_dev(0, 2) + C_dev(1, 0)
               - 2 * C_dev(1, 1) + C_dev(1, 2) + C_dev(2, 0) - 2 * C_dev(2, 1)
               + C_dev(2, 2) - 2.0 * τ_dev.trace()), //
        1.0 / 9.0
            * (-2 * C_dev(0, 0) - 2 * C_dev(0, 1) + 4 * C_dev(0, 2) + C_dev(1, 0)
               + C_dev(1, 1) - 2 * C_dev(1, 2) + C_dev(2, 0) + C_dev(2, 1)
               - 2 * C_dev(2, 2) - 2.0 * τ_dev.trace()),       //
        2.0 / 3.0 * (C_dev(0, 3) - C_dev(1, 3) - C_dev(2, 3)), //
        2.0 / 3.0 * (C_dev(0, 4) - C_dev(1, 4) - C_dev(2, 4)), //
        2.0 / 3.0 * (C_dev(0, 5) - C_dev(1, 5) - C_dev(2, 5)), //

        1.0 / 9.0
            * (-2 * C_dev(0, 0) + C_dev(0, 1) + C_dev(0, 2) + 4 * C_dev(1, 0)
               - 2 * C_dev(1, 1) - 2 * C_dev(1, 2) - 2 * C_dev(2, 0) + C_dev(2, 1)
               + C_dev(2, 2) - 2.0 * τ_dev.trace()), //
        1.0 / 9.0
            * (C_dev(0, 0) - 2 * C_dev(0, 1) + C_dev(0, 2) - 2 * C_dev(1, 0)
               + 4 * C_dev(1, 1) - 2 * C_dev(1, 2) + C_dev(2, 0) - 2 * C_dev(2, 1)
               + C_dev(2, 2) + 4 * τ_dev.trace()), //
        1.0 / 9.0
            * (C_dev(0, 0) + C_dev(0, 1) - 2 * C_dev(0, 2) - 2 * C_dev(1, 0)
               - 2 * C_dev(1, 1) + 4 * C_dev(1, 2) + C_dev(2, 0) + C_dev(2, 1)
               - 2 * C_dev(2, 2) - 2.0 * τ_dev.trace()),            //
        1.0 / 3.0 * (-C_dev(0, 3) + 2 * C_dev(1, 3) - C_dev(2, 3)), //
        1.0 / 3.0 * (-C_dev(0, 4) + 2 * C_dev(1, 4) - C_dev(2, 4)), //
        1.0 / 3.0 * (-C_dev(0, 5) + 2 * C_dev(1, 5) - C_dev(2, 5)), //

        1.0 / 9.0
            * (-2 * C_dev(0, 0) + C_dev(0, 1) + C_dev(0, 2) - 2 * C_dev(1, 0)
               + C_dev(1, 1) + C_dev(1, 2) + 4 * C_dev(2, 0) - 2 * C_dev(2, 1)
               - 2 * C_dev(2, 2) - 2.0 * τ_dev.trace()), //
        1.0 / 9.0
            * (C_dev(0, 0) - 2 * C_dev(0, 1) + C_dev(0, 2) + C_dev(1, 0) - 2 * C_dev(1, 1)
               + C_dev(1, 2) - 2 * C_dev(2, 0) + 4 * C_dev(2, 1) - 2 * C_dev(2, 2)
               - 2.0 * τ_dev.trace()), //
        1.0 / 9.0
            * (C_dev(0, 0) + C_dev(0, 1) - 2 * C_dev(0, 2) + C_dev(1, 0) + C_dev(1, 1)
               - 2 * C_dev(1, 2) - 2 * C_dev(2, 0) - 2 * C_dev(2, 1) + 4 * C_dev(2, 2)
               + 4 * τ_dev.trace()),                                //
        1.0 / 3.0 * (-C_dev(0, 3) - C_dev(1, 3) + 2 * C_dev(2, 3)), //
        1.0 / 3.0 * (-C_dev(0, 4) - C_dev(1, 4) + 2 * C_dev(2, 4)), //
        1.0 / 3.0 * (-C_dev(0, 5) - C_dev(1, 5) + 2 * C_dev(2, 5)), //

        1.0 / 3.0 * (2 * C_dev(3, 0) - C_dev(3, 1) - C_dev(3, 2)),  //
        1.0 / 3.0 * (-C_dev(3, 0) + 2 * C_dev(3, 1) - C_dev(3, 2)), //
        1.0 / 3.0 * (-C_dev(3, 0) - C_dev(3, 1) + 2 * C_dev(3, 2)), //
        C_dev(3, 3) + τ_dev.trace() / 3.0,                          //
        C_dev(3, 4),                                                //
        C_dev(3, 5),                                                //

        1.0 / 3.0 * (2 * C_dev(4, 0) - C_dev(4, 1) - C_dev(4, 2)),  //
        1.0 / 3.0 * (-C_dev(4, 0) + 2 * C_dev(4, 1) - C_dev(4, 2)), //
        1.0 / 3.0 * (-C_dev(4, 0) - C_dev(4, 1) + 2 * C_dev(4, 2)), //
        C_dev(4, 3),                                                //
        C_dev(4, 4) + τ_dev.trace() / 3.0,                          //
        C_dev(4, 5),                                                //

        1.0 / 3.0 * (2 * C_dev(5, 0) - C_dev(5, 1) - C_dev(5, 2)),  //
        1.0 / 3.0 * (-C_dev(5, 0) + 2 * C_dev(5, 1) - C_dev(5, 2)), //
        1.0 / 3.0 * (-C_dev(5, 0) - C_dev(5, 1) + 2 * C_dev(5, 2)), //
        C_dev(5, 3),                                                //
        C_dev(5, 4),                                                //
        C_dev(5, 5) + τ_dev.trace() / 3.0;                          //

    return C;
}
}
