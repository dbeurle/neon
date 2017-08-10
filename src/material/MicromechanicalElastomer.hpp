
#pragma once

#include "LinearElastic.hpp"

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

    /**
     * Compute the shear modulus assuming the temperature is 298K according to
     * the formula μ = n / (k * T)
     * @param n number of chains
     * @return the shear modulus from entropy elastic theory
     */
    double shear_modulus(double const n) const;

    /** @return number of initial chains */
    double number_of_initial_chains() const { return n0; }

    /** @return number of initial segments */
    double number_of_initial_segments() const { return N0_avg; }

    /** @return the current number of chains in the network */
    double update_chains(double const n, double const time_step_size) const;

    /** @return the new segments.  Compute probability based on evolution equation */
    double update_segments(double const N, double const time_step_size);

protected:
    /**
     * Evaluates the probabilty mass function and populates
     * a vector with the number of segments per chain and the corresponding
     * mass density from the probability function.  For physicality, these
     * are drawn from a zero truncated Poisson PMF \sa zero_trunc_poisson_pmf.
     *
     * WARNING: This operation is wildly expensive to compute, so a threshold
     * for the inclusion of the segments per chain is required but will always
     * underestimate the stress in the material by approximately 1%
     */
    void compute_probability_and_segments(double const N);

protected:
    double n0;     // !< Initial number of chains
    double N0_avg; //!< Initial average segments per chain

    double chain_decay_rate;
    double segment_decay_rate;

    std::vector<double> segment_bins;

    double const boltzmann_constant = 1.38064852e-23;
    double const temperature = 298.0;
};

inline double MicromechanicalElastomer::shear_modulus(double const n) const
{
    return n * boltzmann_constant * temperature;
}

inline double MicromechanicalElastomer::update_chains(double const n,
                                                      double const time_step_size) const
{
    return n / (1.0 + 0.0001 * time_step_size);
}

inline double MicromechanicalElastomer::update_segments(double const N,
                                                        double const time_step_size)
{
    compute_probability_and_segments(N);
    return N / (1.0 + time_step_size * segment_decay_rate);
}
}
