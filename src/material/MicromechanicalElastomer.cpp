
#include "MicromechanicalElastomer.hpp"

#include <json/value.h>

#include <cmath>

namespace neon
{
double zero_trunc_poisson_pmf(double const mean, int const k)
{
    return std::exp(k * std::log(mean) - std::lgamma(k + 1)) / (std::exp(mean) - 1.0);
}

MicromechanicalElastomer::MicromechanicalElastomer(Json::Value const& material_data)
    : LinearElastic(material_data)
{
    if (!material_data.isMember("SegmentsPerChain"))
    {
        throw std::runtime_error("SegmentsPerChain not specified in material data\n");
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

void MicromechanicalElastomer::compute_probability_and_segments(double const N)
{
    probability_segments_pairs.clear();
    probability_segments_pairs.reserve(N * 2);

    for (auto k = 2; k < N * 2; ++k)
    {
        // Filter out the insignificant pmf
        if (auto const pmf = zero_trunc_poisson_pmf(N, k); pmf > 1.0e-6)
        {
            probability_segments_pairs.emplace_back(k, pmf);
        }
    }
}
}
