
#pragma once

#include "constitutive/ConstitutiveModel.hpp"

namespace neon::mechanical::solid
{
class Hyperelastic : public ConstitutiveModel
{
public:
    Hyperelastic(InternalVariables& variables);
};
}
