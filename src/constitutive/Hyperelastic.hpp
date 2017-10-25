
#pragma once

#include "ConstitutiveModel.hpp"

namespace neon::solid
{
class Hyperelastic : public ConstitutiveModel
{
public:
    Hyperelastic(InternalVariables& variables);
};
}
