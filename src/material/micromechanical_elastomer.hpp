
#pragma once

#include "isotropic_elastic_property.hpp"

#include <vector>

namespace neon
{
/// micromechanical_elastomer exposes an interface that returns the fundamental
/// material properties from a micromechanical point of view, including the
/// physical constants which make up the shear modulus for an entropy elastic
/// model.
///
/// These additional parameters are associated with the evolution of the segments per
/// chain in the network.
class micromechanical_elastomer : public isotropic_elastic_property
{
public:
    micromechanical_elastomer(json const& material_data);

    /** @return The number of segments per polymer chain */
    auto const segments_per_chain() const { return N; }

protected:
    double N{0.0}; //!< Number of segment per chain
};

/// stochastic_micromechanical_elastomer is responsible for handling a distribution
/// for the material properties of an elastomer
class ageing_micromechanical_elastomer : public micromechanical_elastomer
{
public:
    ageing_micromechanical_elastomer(json const& material_data);

    /// \return The probability per unit time of a chain scission event
    [[nodiscard]] double scission_probability() const noexcept { return pr_scission; }

    /// \return The probability per unit time of chains recombining
    [[nodiscard]] double recombination_probability() const noexcept { return pr_recombination; }

protected:
    double pr_scission{0.0};      /// Scission probability
    double pr_recombination{0.0}; /// Recombination probability
};
}
