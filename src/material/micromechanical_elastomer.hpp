
#pragma once

#include "isotropic_elastic_property.hpp"

#include <vector>

namespace neon
{
/**
 * micromechanical_elastomer exposes an interface that returns the fundamental
 * material properties from a micromechanical point of view, including the
 * physical constants which make up the shear modulus for an entropy elastic
 * model.
 *
 * These additional parameters are associated with the evolution of the segments per
 * chain in the network.
 */
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

    /// \return The constant rate of active segment decay
    [[nodiscard]] double active_segment_decay_rate() const noexcept
    {
        return crosslink_growth_rate;
    }

    /// Update the next number of segments in the history based on the current
    /// number of segments.  This rate is controlled by \p SegmentDecayRate in
    /// the material properties section of the input file.
    /// \param current_segment  Last average number of segments per chain
    /// \param time_step_size Time step size for time discretisation
    /// \return A new segments which is less than current_segment
    [[nodiscard]] double compute_new_segment(double const current_segment,
                                             double const time_step_size) const;

    /// Update the shear modulus based on the initial shear modulus
    /// so the newly created secondary networks occur uniformly
    /// \param time_step_size Time step size for time discretisation
    /// \return A new partial shear modulus
    [[nodiscard]] double compute_new_shear_modulus(double const time_step_size) const;

    /// Compute the updated shear modulus based on the creation of cross-links
    /// causing additional chains with fewer average segments per chain
    /// and the removal of active chains based on the average segments per group
    /// and the rate at which the removal occurs.
    /// \param current_shear_moduli
    /// \param current_segments
    /// \param time_step_size Time step size for time discretisation
    /// \return The new shear modulus based on the average segments per chain
    [[nodiscard]] std::vector<double> scission(std::vector<double> shear_moduli,
                                               std::vector<double> const& segments,
                                               double const time_step_size) const;

protected:
    double pr_scission{0.0};      /// Scission probability
    double pr_recombination{0.0}; /// Recombination probability

    double crosslink_growth_rate{0.0}; /// Rate of new crosslinks
    double segment_decay_rate{0.0};    /// Rate of average segments per chain decay
};
}
