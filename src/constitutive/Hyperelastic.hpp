
#pragma once

#include "ConstitutiveModel.hpp"

namespace neon::solid
{
using InternalVariables = neon::InternalVariables<3>;

class Hyperelastic : public ConstitutiveModel
{
public:
    Hyperelastic(InternalVariables& variables);
};
}
