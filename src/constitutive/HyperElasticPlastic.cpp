
#include "HyperElasticPlastic.hpp"

#include "InternalVariables.hpp"

#include <iostream>
#include <json/value.h>
#include <range/v3/view.hpp>

namespace neon
{
J2Plasticity::J2Plasticity(InternalVariables& variables, Json::Value const& material_data)
    : HyperElasticPlastic(variables), material(material_data), C_e(elastic_moduli())
{
    variables.add(InternalVariables::Tensor::RateOfDeformation,
                  InternalVariables::Tensor::RateOfDeformationPlastic);

    variables.add(InternalVariables::Scalar::VonMisesStress,
                  InternalVariables::Scalar::EffectivePlasticStrain);

    // Add material tangent with the linear elasticity moduli
    variables.add(InternalVariables::Matrix::TruesdellModuli, elastic_moduli());
}

J2Plasticity::~J2Plasticity() = default;

void J2Plasticity::update_internal_variables(double const Δt)
{
    using namespace ranges;

    auto const μ_e = material.mu();
    auto const λ_e = material.lambda();

    // Extract the internal variables
    auto[ɛ_p_list, ɛ_list, σ_list] =
        variables(InternalVariables::Tensor::RateOfDeformationPlastic,
                  InternalVariables::Tensor::RateOfDeformation,
                  InternalVariables::Tensor::Cauchy);

    // Compute the linear strain gradient from the displacement gradient
    ɛ_list = variables(InternalVariables::Tensor::DisplacementGradient) |
             view::transform([](auto const& H) { return 0.5 * (H + H.transpose()); });

    // Retrieve the accumulated internal variables
    auto& α_list = variables(InternalVariables::Scalar::EffectivePlasticStrain);

    auto const& detF_list = variables(InternalVariables::Scalar::DetF);

    auto& C_list = variables(InternalVariables::Matrix::TruesdellModuli);

    Matrix3 const I = Matrix3::Identity();

    // Perform the update algorithm for each quadrature point
    for (auto l = 0; l < ɛ_list.size(); l++)
    {
        auto const& ɛ = ɛ_list[l]; // Total strain
        auto& ɛ_p = ɛ_p_list[l];   // Plastic strain
        auto& σ = σ_list[l];       // Cauchy stress
        auto& α = α_list[l];       // Effective plastic strain

        auto const& J = detF_list[l];

        if (α > 0.3)
        {
            throw std::runtime_error("Excessive plastic strain at quadrature point " +
                                     std::to_string(l) + "\n");
        }

        // Elastic predictor
        Matrix3 const σ_0 = λ_e * (ɛ - ɛ_p).trace() * I + 2.0 * μ_e * (ɛ - ɛ_p);

        // Compute the von Mises stress
        const auto von_mises = von_mises_stress(σ_0);

        σ = σ_0;

        // Compute the initial estimate of the yield function for the material
        // and decide if the stress return needs to be computed
        auto f = von_mises - material.yield_stress(α);

        if (f <= 0.0) continue;

        // Compute the normal direction to the yield surface which remains
        // constant throughout the radial return method
        Matrix3 const normal = deviatoric(σ_0) / deviatoric(σ_0).norm();

        // Initialize the plastic increment
        auto Δλ = 0.0;

        // Perform the return mapping algorithm
        int iterations = 0, max_iterations = 100;
        while (f > 1.0e-4 && iterations < max_iterations)
        {
            // Increment in plasticity rate
            const auto δλ = f / (3.0 * μ_e + material.hardening_modulus(α));

            // Increment in plastic strain
            Matrix3 const Δɛ_p = -δλ * std::sqrt(3.0 / 2.0) * normal;

            // Plastic strain update
            ɛ_p += Δɛ_p;

            // Cauchy stress update
            σ += 2.0 * μ_e * Δɛ_p;

            // Accumulated plastic strain update
            α += δλ;

            // Plastic rate update
            Δλ += δλ;

            // Evaluate the yield function
            f = (von_mises - 3.0 * μ_e * Δλ) - material.yield_stress(α);

            iterations++;
        }

        std::cout << "Accumulated plastic strain " << α << std::endl;

        if (iterations == max_iterations)
        {
            std::cout << "Accumulated plastic strain " << α << std::endl;
            std::cout << "Yield function after mapping " << f << std::endl;
            std::cout << "New yield point " << material.yield_stress(α) / 1.0e6 << std::endl;
            std::cout << "VonMises stress " << von_mises_stress(σ) << std::endl << std::endl;

            throw std::runtime_error("Radial return failure at global integration point " +
                                     std::to_string(l) + "\n");
        }
        // Compute the elastic-plastic tangent modulus
        C_list[l] = algorithmic_tangent(α, σ, σ_0, normal, J);

        if ((C_list[l] - C_list[l].transpose()).norm() > 0.00001)
        {
            throw std::runtime_error("Matrix not symmetric!\n");
        }
    }
}

Matrix J2Plasticity::elastic_moduli() const
{
    auto const[λ, μ] = material.Lame_parameters();
    Matrix C(6, 6);
    C << λ + 2.0 * μ, λ, λ, 0.0, 0.0, 0.0, //
        λ, λ + 2.0 * μ, λ, 0.0, 0.0, 0.0,  //
        λ, λ, λ + 2.0 * μ, 0.0, 0.0, 0.0,  //
        0.0, 0.0, 0.0, μ, 0.0, 0.0,        //
        0.0, 0.0, 0.0, 0.0, μ, 0.0,        //
        0.0, 0.0, 0.0, 0.0, 0.0, μ;
    return C;
}

Matrix J2Plasticity::algorithmic_tangent(double const α,
                                         Matrix3 const& σ,
                                         Matrix3 const& σ_0,
                                         Matrix3 const& normal,
                                         double const J) const
{
    auto const[λ_e, μ_e] = material.Lame_parameters();

    auto const β = material.yield_stress(α) / von_mises_stress(σ_0);

    auto const γ = 1.0 / (1.0 + material.hardening_modulus(α) / (3 * μ_e));

    auto const γ_bar = γ - (1.0 - β);

    Matrix C_alg(6, 6);

    auto diagonal = (4 * β * μ_e + 2 * μ_e + 3 * λ_e) / 3.0;
    auto off_diagonal = -(2 * β * μ_e - 2 * μ_e - 3 * λ_e) / 3.0;

    C_alg << diagonal, off_diagonal, off_diagonal, 0, 0, 0, //
        off_diagonal, diagonal, off_diagonal, 0, 0, 0,      //
        off_diagonal, off_diagonal, diagonal, 0, 0, 0,      //
        0, 0, 0, β * μ_e, 0, 0,                             //
        0, 0, 0, 0, β * μ_e, 0,                             //
        0, 0, 0, 0, 0, β * μ_e;

    C_alg -= 2.0 * μ_e * γ_bar * this->outer_product(normal);

    return transformJaumannToTruesdellKirchoff(C_alg, σ, J);
}

Matrix J2Plasticity::transformJaumannToTruesdellKirchoff(Matrix const C_τ_J,
                                                         Matrix3 const& σ,
                                                         double const J) const
{
    Matrix Cdash(6, 6);
    // clang-format off
    Cdash << 2.0 * σ(0, 0),             0,             0,                         0,       σ(0, 2),                   σ(0, 1), //
                         0, 2.0 * σ(1, 1),             0,                   σ(1, 2),             0,                   σ(1, 0), //
                         0,             0, 2.0 * σ(2, 2),                   σ(2, 1),       σ(2, 0),                         0, //
                         0,       σ(2, 1),       σ(1, 2), 0.5 * (σ(1, 1) + σ(2, 2)), 0.5 * σ(1, 0),             0.5 * σ(2, 0), //
                   σ(2, 0),             0,       σ(0, 2),             0.5 * σ(0, 1), 0.5 * (σ(0, 0) + σ(2, 2)), 0.5 * σ(2, 1), //
                   σ(1, 0),       σ(0, 1),             0,             0.5 * σ(0, 2),             0.5 * σ(1, 2), 0.5 * (σ(0, 0) + σ(1, 1));
    // clang-format on
    return C_τ_J / J - Cdash;
}

Matrix J2Plasticity::outer_product(Matrix3 const& n) const
{
    return voigt(n) * voigt(n).transpose();
}
}
