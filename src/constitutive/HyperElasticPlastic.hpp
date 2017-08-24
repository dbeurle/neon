
#pragma once

#include "ConstitutiveModel.hpp"

namespace neon
{
class HyperElasticPlastic : public ConstitutiveModel
{
public:
    HyperElasticPlastic(InternalVariables& variables) : ConstitutiveModel(variables) {}
};
}
