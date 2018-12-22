
#pragma once

/// @file

#include "mesh/boundary/boundary.hpp"

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
