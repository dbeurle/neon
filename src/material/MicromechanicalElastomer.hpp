
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
 * These additional parameters are associated with the evolution of the segments per
 * chain in the network.
 */
class MicromechanicalElastomer : public LinearElastic
{
public:
    MicromechanicalElastomer(Json::Value const& material_data);

    /** @return The number of segments per polymer chain */
    auto const segments_per_chain() const { return N; }

protected:
    double N{25.0}; //!< Number of segment per chain
};

/**
 * StochasticMicromechanicalElastomer is responsible for handling a distribution
 * for the material properties of an elastomer
 */
class StochasticMicromechanicalElastomer : public LinearElastic
{
public:
    StochasticMicromechanicalElastomer(Json::Value const& material_data);

    /** Updates the temperature for the material property */
    void update_temperature(double const T_new) { temperature = T_new; }

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

    std::vector<double> segments, chains, shear_moduli;

    double const boltzmann_constant{1.38064852e-23};
    double temperature{298.0}; //!< Temperature of elastomer

    double non_affine_stretch_parameter{1.0}; //!< Three-dimensional locking characteristics (p)
    double effective_tube_geometry{1.0};      //!< Additional constraint stiffness (U)
    double non_affine_tube_parameter{1.0};    //!< Shape of constraint stress (q)
};
}
