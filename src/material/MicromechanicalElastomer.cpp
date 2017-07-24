
#include "MicromechanicalElastomer.hpp"

#include <json/value.h>

#include <iostream>

namespace neon
{
double zero_trunc_poisson_pmf(double const λ, int const k)
{
    return std::exp(k * std::log(λ) - std::lgamma(k + 1)) / (std::exp(λ) - 1.0);
}

MicromechanicalElastomer::MicromechanicalElastomer(Json::Value const& material_data)
    : LinearElastic(material_data)
{
    if (!material_data.isMember("SegmentsPerChain"))
    {
        throw std::runtime_error("SegmentsPerChain not specified in material data\n");
    }
    N = material_data["SegmentsPerChain"].asDouble();

    chain_decay_rate = material_data.isMember("ChainDecayRate")
                           ? material_data["ChainDecayRate"].asDouble()
                           : 0.0;

    segment_decay_rate = material_data.isMember("SegmentDecayRate")
                             ? material_data["SegmentDecayRate"].asDouble()
                             : 0.0;

    n0 = LinearElastic::shear_modulus() / (boltzmann_constant * temperature);

    compute_probability_and_segments();
}

void MicromechanicalElastomer::compute_probability_and_segments()
{
    probability_segments_pairs.clear();
    probability_segments_pairs.reserve(N * 3);

    for (auto k = 2; k < N * 3; ++k)
    {
        // Filter out the insignificant pmf
        if (auto const pmf = zero_trunc_poisson_pmf(N, k); pmf > 1.0e-6)
        {
            probability_segments_pairs.emplace_back(k, pmf);
        }
    }
}
}
