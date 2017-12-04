
#pragma once

#include "ConstitutiveModel.hpp"

namespace neon::mechanical::solid
{
class HyperElasticPlastic : public ConstitutiveModel
{
public:
    HyperElasticPlastic(InternalVariables& variables) : ConstitutiveModel(variables) {}
};
}
