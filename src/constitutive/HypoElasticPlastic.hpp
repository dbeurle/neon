
#pragma once

#include "ConstitutiveModel.hpp"

namespace neon
{
/**
 * HypoElasticPlastic is a base class for all hypoelastic-plastic material
 * models
 */
class HypoElasticPlastic : public ConstitutiveModel
{
public:
    HypoElasticPlastic(InternalVariables& variables) : ConstitutiveModel(variables) {}
};
}
