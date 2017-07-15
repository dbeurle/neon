
#pragma once

#include "ConstitutiveModel.hpp"

namespace neon
{
class Hyperelastic : public ConstitutiveModel
{
public:
    Hyperelastic(InternalVariables& variables);
};
}
