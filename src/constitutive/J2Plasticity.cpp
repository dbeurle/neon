
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

void J2Plasticity::update_internal_variables(double const Δt)
{
    using namespace ranges;

    auto const μ_e = material.mu();
    auto const λ_e = material.lambda();

    // Extract the internal variables
    auto[ɛ_p_list, ɛ_list, σ_list] = variables(InternalVariables::Tensor::LinearisedPlasticStrain,
                                               InternalVariables::Tensor::LinearisedStrain,
                                               InternalVariables::Tensor::Cauchy);

    // Retrieve the accumulated internal variables
    auto[α_list, von_mises_list] = variables(InternalVariables::Scalar::EffectivePlasticStrain,
                                             InternalVariables::Scalar::VonMisesStress);

    auto& C_list = variables(InternalVariables::Matrix::TruesdellModuli);

    // Compute the linear strain gradient from the displacement gradient
    ɛ_list = variables(InternalVariables::Tensor::DisplacementGradient)
             | view::transform([](auto const& H) { return 0.5 * (H + H.transpose()); });

    Matrix3 const I = Matrix3::Identity();

    // Perform the update algorithm for each quadrature point
    for (auto l = 0; l < ɛ_list.size(); l++)
    {
        auto const& ɛ = ɛ_list[l]; // Total strain
        auto& ɛ_p = ɛ_p_list[l];   // Plastic strain
        auto& σ = σ_list[l];       // Cauchy stress
        auto& α = α_list[l];       // Effective plastic strain
        auto& von_mises = von_mises_list[l];

        if (α > 0.3)
        {
            throw std::runtime_error("Excessive plastic strain at quadrature point "
                                     + std::to_string(l) + "\n");
        }

        // Elastic stress predictor
        Matrix3 const σ_0 = λ_e * (ɛ - ɛ_p).trace() * I + 2.0 * μ_e * (ɛ - ɛ_p);

        // Compute the von Mises stress
        von_mises = von_mises_stress(σ_0);

        σ = σ_0;

        // Compute the initial estimate of the yield function for the material
        // and decide if the stress return needs to be computed
        auto f = von_mises - material.yield_stress(α);

        if (f <= 0.0)
        {
            C_list[l] = C_e;
            continue;
        }

        // Compute the normal direction to the yield surface which remains
        // constant throughout the radial return method
        Matrix3 const normal = deviatoric(σ_0) / deviatoric(σ_0).norm();

        // Initialize the plastic increment
        auto Δλ = 0.0;

        // Perform the return mapping algorithm
        int iterations = 0, max_iterations = 25;
        while (f > 1.0e-4 && iterations < max_iterations)
        {
            // Increment in plasticity rate
            const auto δλ = f / (3.0 * μ_e + material.hardening_modulus(α + Δλ));

            // Plastic rate update
            Δλ += δλ;

            // Evaluate the yield function
            f = (von_mises - 3.0 * μ_e * Δλ) - material.yield_stress(α + Δλ);

            iterations++;
        }

        // Plastic strain update
        ɛ_p += Δλ * std::sqrt(3.0 / 2.0) * normal;

        // Cauchy stress update
        σ -= 2.0 * μ_e * Δλ * std::sqrt(3.0 / 2.0) * normal;

        // Update the Von Mises stress
        von_mises = von_mises_stress(σ);

        // Accumulated plastic strain update
        α += Δλ;

        if (iterations == max_iterations)
        {
            std::cout << "MAXIMUM NUMBER OF ITERATIONS IN RADIAL RETURN REACHED\n";
            std::cout << "Accumulated plastic strain " << α << "\n";
            std::cout << "Yield function after mapping " << f << "\n";
            std::cout << "Current yield stress " << material.yield_stress(α) << "\n";
            std::cout << "Von Mises stress " << von_mises_stress(σ) << "\n"
                      << "\n";

            throw std::runtime_error("Radial return failure at global integration point "
                                     + std::to_string(l) + "\n");
        }
        // Compute the elastic-plastic tangent modulus
        C_list[l] = algorithmic_tangent(Δλ, α, von_mises_stress(σ_0), normal);
    }
}

CMatrix J2Plasticity::elastic_moduli() const
{
    auto const[λ, μ] = material.Lame_parameters();
    CMatrix C(6, 6);
    C << λ + 2.0 * μ, λ, λ, 0.0, 0.0, 0.0, //
        λ, λ + 2.0 * μ, λ, 0.0, 0.0, 0.0,  //
        λ, λ, λ + 2.0 * μ, 0.0, 0.0, 0.0,  //
        0.0, 0.0, 0.0, μ, 0.0, 0.0,        //
        0.0, 0.0, 0.0, 0.0, μ, 0.0,        //
        0.0, 0.0, 0.0, 0.0, 0.0, μ;
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

CMatrix J2Plasticity::incremental_tangent(double const Δλ, double const von_mises) const
{
    auto const G = material.mu();
    return C_e - Δλ * 6 * G * G / von_mises * I_dev;
}

CMatrix J2Plasticity::algorithmic_tangent(double const Δλ,
                                          double const α,
                                          double const von_mises,
                                          Matrix3 const& n) const
{
    auto const G = material.shear_modulus();
    auto const H = material.hardening_modulus(α);

    return C_e - Δλ * 6.0 * std::pow(G, 2) / von_mises * I_dev
           + 6.0 * std::pow(G, 2) * (Δλ / von_mises - 1.0 / (3.0 * G + H)) * voigt(n)
                 * voigt(n).transpose();
}
}
