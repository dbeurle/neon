
#pragma once

#include "ConstitutiveModel.hpp"

namespace neon::mech::solid
{
class HyperElasticPlastic : public ConstitutiveModel
{
public:
    HyperElasticPlastic(InternalVariables& variables) : ConstitutiveModel(variables) {}
};
}
