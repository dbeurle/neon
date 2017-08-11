
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

    variables.add(InternalVariables::Scalar::Chains,
                  InternalVariables::Scalar::ShearModuli);

    // Shrink these down to the correct size
    variables(InternalVariables::Scalar::Chains).resize(material.groups(), 0.0);
    variables(InternalVariables::Scalar::ShearModuli).resize(material.groups(), 0.0);
    variables(InternalVariables::Scalar::Chains).shrink_to_fit();
    variables(InternalVariables::Scalar::ShearModuli).shrink_to_fit();

    // Fill the data with material properties using the material class
    variables(InternalVariables::Scalar::Chains) = material.chain_groups();
    variables(InternalVariables::Scalar::ShearModuli) = material.shear_moduli_groups();

    // Commit these to history in case of failure on first time step
    variables.commit();
}

void AffineMicrosphere::update_internal_variables(double const time_step_size)
{
    using namespace ranges;
    using view::zip;
    using view::transform;

    // TODO Change these back once OpenMP allows structured bindings

    // Get references into the hash table
    auto& F_list = variables(InternalVariables::Tensor::DeformationGradient);
    auto& cauchy_stress_list = variables(InternalVariables::Tensor::Cauchy);
    auto& dev_stress_list = variables(InternalVariables::Tensor::Kirchhoff);

    auto const& detF_list = variables(InternalVariables::Scalar::DetF);

    auto& n_list = variables(InternalVariables::Scalar::Chains);
    auto& G_list = variables(InternalVariables::Scalar::ShearModuli);

    // Update the material properties
    n_list = material.update_chains(n_list, time_step_size);
    G_list = material.compute_shear_moduli(n_list);

    auto const K = material.bulk_modulus();

/*----------------------------------------------------------------------------*
 *                          Stress computation                                *
 *----------------------------------------------------------------------------*/

#pragma omp parallel for
    for (auto l = 0; l < F_list.size(); ++l)
    {
        auto& stress_dev = dev_stress_list[l]; // Deviatoric stress
        auto const& F = F_list[l];             // Deformation gradient
        auto const& J = detF_list[l];          // Determinant of the deformation gradient

        Matrix3 const unimodular_F = std::pow(J, -1.0 / 3.0) * F;

        stress_dev = weighting(G_list, Matrix3::Zero().eval(), [&](auto const& N) -> Matrix3 {
            return compute_kirchhoff_stress(unimodular_F, N);
        });
    }

    // Perform the projection of the stresses
    cauchy_stress_list = zip(dev_stress_list, detF_list)
                         | transform([&](auto const& tpl) -> Matrix3 {

                               auto const & [ cauchy_stress_dev, J ] = tpl;

                               auto const pressure = J
                                                     * volumetric_free_energy_derivative(J,
                                                                                         K);

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
        auto const& cauchy_stress_dev = dev_stress_list[l];
        auto const& J = detF_list[l];

        auto const pressure = J * volumetric_free_energy_derivative(J, K);
        auto const kappa = std::pow(J, 2) * volumetric_free_energy_second_derivative(J, K);

        Matrix3 const unimodular_F = std::pow(J, -1.0 / 3.0) * F;

        CMatrix const Ddev = weighting(G_list,
                                       CMatrix::Zero(6, 6).eval(),
                                       [&](auto const& N) -> CMatrix {
                                           return compute_material_matrix(unimodular_F, N);
                                       });

        D_list[l] = deviatoric_projection(Ddev, cauchy_stress_dev)
                    + (kappa + pressure) * IoI - 2.0 * pressure * I;
    }
}

Matrix3 AffineMicrosphere::deviatoric_projection(double const pressure,
                                                 Matrix3 const& cauchy_stress_dev) const
{
    Matrix3 P_double_dot_stress_dev;
    P_double_dot_stress_dev << 2.0 * cauchy_stress_dev(0, 0) / 3.0
                                   - cauchy_stress_dev(1, 1) / 3.0
                                   - cauchy_stress_dev(2, 2) / 3.0, //
        cauchy_stress_dev(0, 1),                                    //
        cauchy_stress_dev(0, 2),                                    //

        cauchy_stress_dev(0, 1), //
        -cauchy_stress_dev(0, 0) / 3.0 + 2.0 * cauchy_stress_dev(1, 1) / 3.0
            - cauchy_stress_dev(2, 2) / 3.0, //
        cauchy_stress_dev(1, 2),             //

        cauchy_stress_dev(0, 2), //
        cauchy_stress_dev(1, 2), //
        -cauchy_stress_dev(0, 0) / 3.0 - cauchy_stress_dev(1, 1) / 3.0
            + 2.0 * cauchy_stress_dev(2, 2) / 3.0;
    return pressure * Matrix3::Identity() + P_double_dot_stress_dev;
}
}
