
#include "MicromechanicalElastomer.hpp"

#include <json/value.h>

namespace neon
{
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
}
}
