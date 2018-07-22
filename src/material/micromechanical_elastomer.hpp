
#pragma once

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
    [[nodiscard]] auto segments_per_chain() const noexcept { return N; }

protected:
    /// Number of segment per chain
    double N{0.0};
};

/// ageing_micromechanical_elastomer is responsible for storing the
/// probabilities of a scission and a combination event
class ageing_micromechanical_elastomer : public micromechanical_elastomer
{
public:
    explicit ageing_micromechanical_elastomer(json const& material_data);

    /// \return Bond length for a carbon-carbon bond
    [[nodiscard]] double bond_length() const noexcept { return C_bond_length; }

    /// \return The probability per unit time of a chain scission event
    [[nodiscard]] double scission_probability() const noexcept { return scission; }

    /// \return The probability per unit time of chains recombining
    [[nodiscard]] double recombination_probability() const noexcept { return recombination; }

    [[nodiscard]] double creation_rate(double const active_shear_modulus,
                                       double const inactive_shear_modulus,
                                       double const active_segments,
                                       double const inactive_segments) const;

    [[nodiscard]] vector5 integrate(vector5 z, double const time_step_size) const;

protected:
    /// Scission probability
    double scission{0.0};
    /// Recombination probability
    double recombination{0.0};
    /// Carbon-carbon bond length
    double C_bond_length{150.0E-5};
};
}
