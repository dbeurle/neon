
#pragma once

#include "ConstitutiveModel.hpp"

namespace neon::mech::solid
{
using InternalVariables = neon::InternalVariables<3>;

class Hyperelastic : public ConstitutiveModel
{
public:
    Hyperelastic(InternalVariables& variables);
};
}
