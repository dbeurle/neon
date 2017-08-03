
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

void AffineMicrosphere::update_internal_variables(double const time_step_size)
{
    using namespace ranges;

    // TODO Change these back once OpenMP allows structured bindings

    // Get references into the hash table
    auto& F_list = variables(InternalVariables::Tensor::DeformationGradient);
    auto& cauchy_stress_list = variables(InternalVariables::Tensor::Cauchy);
    auto& deviatoric_stress_list = variables(InternalVariables::Tensor::Kirchhoff);

    auto const& detF_list = variables(InternalVariables::Scalar::DetF);
    auto& n_list = variables(InternalVariables::Scalar::Chains);
    auto& N_list = variables(InternalVariables::Scalar::Segments);

    // Update the number of chains in the network
    n_list = n_list | view::transform([&](auto const& n) {
                 return material.update_chains(n, time_step_size);
             });

    N_list = N_list | view::transform([&](auto const& N) {
                 return material.update_segments(N, time_step_size);
             });

/*----------------------------------------------------------------------------*
 *                          Stress computation                                *
 *----------------------------------------------------------------------------*/

#pragma omp parallel for
    for (auto l = 0; l < F_list.size(); ++l)
    {
        auto& stress_dev = deviatoric_stress_list[l];
        auto const& F = F_list[l];
        auto const& J = detF_list[l];

        auto const shear_modulus = material.shear_modulus(n_list[l]);

        Matrix3 const unimodular_F = std::pow(J, -1.0 / 3.0) * F;

        stress_dev = shear_modulus
                     * weighting(Matrix3::Zero().eval(), [&](auto const& N) -> Matrix3 {
                           return compute_kirchhoff_stress(unimodular_F, N);
                       });
    }

    // Perform the projection of the stresses
    cauchy_stress_list = view::zip(deviatoric_stress_list, detF_list, n_list)
                         | view::transform([&, this](auto const& tpl) -> Matrix3 {

                               auto const & [ cauchy_stress_dev, J, n ] = tpl;

                               auto const shear_modulus = material.shear_modulus(n);

                               auto const pressure = J
                                                     * volumetric_free_energy_derivative(J,
                                                                                         shear_modulus);

                               return deviatoric_projection(pressure, cauchy_stress_dev)
                                      / J;
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
        auto const& cauchy_stress_dev = deviatoric_stress_list[l];
        auto const& J = detF_list[l];

        auto const shear_modulus = material.shear_modulus(n_list[l]);

        auto const pressure = J * volumetric_free_energy_derivative(J, shear_modulus);
        auto const kappa = std::pow(J, 2)
                           * volumetric_free_energy_second_derivative(J, shear_modulus);

        Matrix3 const unimodular_F = std::pow(J, -1.0 / 3.0) * F;

        CMatrix const D = weighting(CMatrix::Zero(6, 6).eval(),
                                    [&](auto const& N) -> CMatrix {
                                        return compute_material_matrix(unimodular_F, N);
                                    });

        D_list[l] = deviatoric_projection(D, cauchy_stress_dev) + (kappa + pressure) * IoI
                    - 2.0 * pressure * I;
    }
}

Matrix3 AffineMicrosphere::deviatoric_projection(double const pressure,
                                                 Matrix3 const& cauchy_stress_dev) const
{
    Matrix3 P_double_dot_stress_dev;
    P_double_dot_stress_dev << 2 * cauchy_stress_dev(0, 0) / 3.0
                                   - cauchy_stress_dev(1, 1) / 3.0
                                   - cauchy_stress_dev(2, 2) / 3.0, //
        cauchy_stress_dev(0, 1),                                    //
        cauchy_stress_dev(0, 2),                                    //

        cauchy_stress_dev(0, 1), //
        -cauchy_stress_dev(0, 0) / 3.0 + 2 * cauchy_stress_dev(1, 1) / 3.0
            - cauchy_stress_dev(2, 2) / 3.0, //
        cauchy_stress_dev(1, 2),             //

        cauchy_stress_dev(0, 2), //
        cauchy_stress_dev(1, 2), //
        -cauchy_stress_dev(0, 0) / 3.0 - cauchy_stress_dev(1, 1) / 3.0
            + 2 * cauchy_stress_dev(2, 2) / 3.0;
    return pressure * Matrix3::Identity() + P_double_dot_stress_dev;
}
}
