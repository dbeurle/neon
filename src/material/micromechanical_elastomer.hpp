
#pragma once

/// \file

#include "isotropic_elastic_property.hpp"
#include "numeric/dense_matrix.hpp"

namespace neon
{
/// micromechanical_elastomer exposes an interface that returns the fundamental
/// material properties from a micromechanical point of view, including the
/// physical constants which make up the shear modulus for an entropy elastic
/// model.  These additional parameters are associated with the evolution of the
/// segments per chain in the network.
class micromechanical_elastomer : public isotropic_elastic_property
{
public:
    explicit micromechanical_elastomer(json const& material_data);

    /// \return The number of segments per polymer chain
    [[nodiscard]] auto segments_per_chain() const noexcept -> double
    {
        return m_segments_per_chain;
    }

protected:
    /// Number of segment per chain
    double m_segments_per_chain = 0.0;
};

/// ageing_micromechanical_elastomer is responsible for storing the
/// probabilities of a scission and a combination event
class ageing_micromechanical_elastomer : public micromechanical_elastomer
{
public:
    explicit ageing_micromechanical_elastomer(json const& material_data);

    /// \return The probability per unit time of a chain scission event
    [[nodiscard]] double scission_probability() const noexcept { return scission; }

    /// \return The probability per unit time of chains recombining
    [[nodiscard]] double recombination_probability() const noexcept { return recombination; }

    [[nodiscard]] double creation_rate(double const active_shear_modulus,
                                       double const inactive_shear_modulus,
                                       double const active_segments,
                                       double const inactive_segments) const;

    /// The initial inactive shear modulus is related to the cure_time parameter
    /// provided in the input file.  This provides an initial condition for the
    /// inactive network number of chains.
    [[nodiscard]] auto initial_inactive_shear_modulus() const noexcept -> double
    {
        return m_shear_modulus_ia;
    }

    [[nodiscard]] vector5 integrate(vector5 const& z, double const time_step_size) const;

protected:
    /// Scission probability
    double scission = 0.0;
    /// Recombination probability
    double recombination = 0.0;
    /// Initial value of the inactive shear modulus (cure_time dependent)
    double m_shear_modulus_ia = 0.0;
};
}
