
#pragma once

#include "LinearElastic.hpp"

#include <cmath>
#include <vector>

namespace neon
{
/**
 * MicromechanicalElastomer exposes an interface that returns the fundamental
 * material properties from a micromechanical point of view, including the
 * physical constants which make up the shear modulus for an entropy elastic
 * model.
 *
 * These additional parameters are to do with the evolution of the segments per
 * chain in the network.
 */
class MicromechanicalElastomer : public LinearElastic
{
public:
    MicromechanicalElastomer(Json::Value const& material_data);

    /** @return the time varying number of segments in the chain */
    double segments_per_chain() const { return N; }

    /**
     * Compute the shear modulus
     * @param n number of chains
     * @return the time varying shear modulus
     */
    double shear_modulus(double const n) const;

    double number_of_chains() const { return n0; }

    /** @return the current number of chains in the network */
    double evolve_chains(double const n, double const Δt) const;

    auto const& segment_probability() const { return probability_segments_pairs; }

protected:
    void allocate_probability_and_segments();

protected:
    double n0; // !< Initial number of chains

    double N; //!< Segments per chain

    double chain_decay_rate;

    std::vector<std::pair<double, double>> probability_segments_pairs;

    double const boltzmann_constant = 1.38064852e-23;
    double const temperature = 298.0;
};

inline double MicromechanicalElastomer::shear_modulus(double const n) const
{
    return n * boltzmann_constant * temperature;
}

inline double MicromechanicalElastomer::evolve_chains(double const n, double const Δt) const
{
    return n / (1.0 + chain_decay_rate * Δt);
}
}
