
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
    n0 = LinearElastic::shear_modulus() / (boltzmann_constant * temperature);

    allocate_probability_and_segments();
}

void MicromechanicalElastomer::allocate_probability_and_segments()
{
    probability_segments_pairs.clear();
    probability_segments_pairs.reserve(N * 2);

    for (auto k = 8; k < 45; ++k)
    {
        probability_segments_pairs.emplace_back(k, zero_trunc_poisson_pmf(N, k));
    }
}
}
