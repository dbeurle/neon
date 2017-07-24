
#include "AffineMicrosphere.hpp"

#include "InternalVariables.hpp"
#include "numeric/DenseTypes.hpp"

#include <json/json.h>
#include <range/v3/algorithm.hpp>
#include <range/v3/view.hpp>

#include <exception>
#include <iostream>
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
    variables.add(InternalVariables::Scalar::Segments);

    for (auto& n : variables(InternalVariables::Scalar::Chains))
    {
        n = material.number_of_initial_chains();
    }
    for (auto& N : variables(InternalVariables::Scalar::Segments))
    {
        N = material.number_of_initial_segments();
    }
}

void AffineMicrosphere::update_internal_variables(double const Δt)
{
    using namespace ranges;

    // TODO Change these back once OpenMP allows structured bindings

    // Get references into the hash table
    auto& F_list = variables(InternalVariables::Tensor::DeformationGradient);
    auto& σ_list = variables(InternalVariables::Tensor::Cauchy);
    auto& τ_list = variables(InternalVariables::Tensor::Kirchhoff);

    auto const& detF_list = variables(InternalVariables::Scalar::DetF);
    auto& n_list = variables(InternalVariables::Scalar::Chains);
    auto& N_list = variables(InternalVariables::Scalar::Segments);

    // Update the number of chains in the network
    n_list = n_list | view::transform([&](auto const& n) {
                 return material.update_chains(n, Δt);
             });

    N_list = N_list | view::transform([&](auto const& N) {
                 return material.update_segments(N, Δt);
             });

/*----------------------------------------------------------------------------*
 *                          Stress computation                                *
 *----------------------------------------------------------------------------*/

#pragma omp parallel for
    for (auto l = 0; l < F_list.size(); ++l)
    {
        auto& τ = τ_list[l];
        auto const& F = F_list[l];
        auto const& J = detF_list[l];

        auto const μ = material.shear_modulus(n_list[l]);

        Matrix3 const unimodular_F = std::pow(J, -1.0 / 3.0) * F;

        τ = μ * weighting(Matrix3::Zero().eval(), [&](auto const& N) -> Matrix3 {
                return compute_kirchhoff_stress(unimodular_F, N);
            });
    }

    // Perform the projection of the stresses
    σ_list = view::zip(τ_list, detF_list, n_list)
             | view::transform([&, this](auto const& tpl) -> Matrix3 {

                   auto const & [ τ_dev, J, n ] = tpl;

                   auto const μ = material.shear_modulus(n);

                   auto const pressure = J * volumetric_free_energy_derivative(J, μ);

                   return deviatoric_projection(pressure, τ_dev) / J;
               });

    /*------------------------------------------------------------------------*
     *                     Tangent material computation                       *
     *------------------------------------------------------------------------*/

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

        Matrix3 const unimodular_F = std::pow(J, -1.0 / 3.0) * F;

        CMatrix const D = weighting(CMatrix::Zero(6, 6).eval(),
                                    [&](auto const& N) -> CMatrix {
                                        return compute_material_matrix(unimodular_F, N);
                                    });

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
}
