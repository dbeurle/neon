
#include "MicromechanicalElastomer.hpp"

#include <json/value.h>

#include <cmath>

namespace neon
{
double binomial_coefficient(double const mean, double const probability_success = 0.5)
{
    return std::exp(std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1));
}

MicromechanicalElastomer::MicromechanicalElastomer(Json::Value const& material_data)
    : LinearElastic(material_data)
{
    if (!material_data.isMember("NumberOfSegmentBins"))
    {
        throw std::runtime_error("NumberOfSegmentBins not specified in material data\n");
    }

    N0_avg = material_data["SegmentsPerChain"].asDouble();

    chain_decay_rate = material_data.isMember("ChainDecayRate")
                           ? material_data["ChainDecayRate"].asDouble()
                           : 0.0;

    segment_decay_rate = material_data.isMember("SegmentDecayRate")
                             ? material_data["SegmentDecayRate"].asDouble()
                             : 0.0;

    n0 = LinearElastic::shear_modulus() / (boltzmann_constant * temperature);

    compute_probability_and_segments(N0_avg);
}

void MicromechanicalElastomer::compute_probability_and_segments(double const N) {}
}
