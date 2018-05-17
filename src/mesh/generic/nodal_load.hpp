
#pragma once

#include "mesh/generic/boundary.hpp"

namespace neon
{
class nodal_load : public boundary
{
public:
protected:
    /// DoF to apply the load
    std::int32_t prescribed_dof;
};
}
