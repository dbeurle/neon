
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

    int const groups() const { return number_of_groups; }

    /**
     * Compute the shear modulus assuming the temperature is 298K according to
     * the formula μ = n / (k * T)
     * @param n number of chains
     * @return the shear modulus from entropy elastic theory
     */
    auto const& shear_moduli_groups() const { return shear_moduli; }

    /** @return vector of chain groups from initial material properties */
    auto const& chain_groups() const { return chains; }

    /** @return vector of segment groups from initial material properties */
    auto const& segment_groups() const { return segments; }

    /** @return the current number of chains in the network */
    std::vector<double> update_chains(std::vector<double> const& chains_old,
                                      double const time_step_size);

    std::vector<double> compute_shear_moduli(std::vector<double> const& chains);

protected:
    void compute_chains_and_segments(Json::Value const& segments_data);

protected:
    double p_scission = 0.0; //!< Probability that a segment is scissioned

    int number_of_groups = 1;

    std::vector<double> segments;
    std::vector<double> chains;

    std::vector<double> shear_moduli;

    double const boltzmann_constant = 1.38064852e-23;
    double const temperature = 298.0;
};
}
