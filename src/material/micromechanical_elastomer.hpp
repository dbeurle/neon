
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
    double N{25.0}; //!< Number of segment per chain
};

/**
 * stochastic_micromechanical_elastomer is responsible for handling a distribution
 * for the material properties of an elastomer
 */
class stochastic_micromechanical_elastomer : public isotropic_elastic_property
{
public:
    stochastic_micromechanical_elastomer(json const& material_data);

    /** Updates the temperature for the material property */
    void update_temperature(double const T_new) { temperature = T_new; }

    int const groups() const { return number_of_groups; }

    /**
     * Compute the shear modulus assuming the temperature is 298K according to
     * \f$ \mu = \frac{n}{kT} \f$ where \f$ n \f$ is the number of chains.
     *
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
    void compute_chains_and_segments(json const& segments_data);

protected:
    double p_scission{0.0}; //!< Segment scission probability

    std::int32_t number_of_groups{1};

    std::vector<double> segments, chains, shear_moduli;

    double const boltzmann_constant{1.38064852e-23};
    double temperature{298.0}; //!< Temperature of elastomer

    double non_affine_stretch_parameter{1.0}; //!< Three-dimensional locking characteristics (p)
    double effective_tube_geometry{1.0};      //!< Additional constraint stiffness (U)
    double non_affine_tube_parameter{1.0};    //!< Shape of constraint stress (q)
};
}
