
#pragma once

#include "mesh/boundary/boundary_condition.hpp"

/// @file

namespace neon
{
class nodal_load : public boundary_condition
{
public:
protected:
    /// DoF to apply the load
    std::int32_t prescribed_dof;
};
}
