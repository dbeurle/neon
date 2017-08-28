
#pragma once

#include "Boundary.hpp"

#include "numeric/DenseTypes.hpp"

namespace neon
{
/**
 * Neumann is a base class for Neumann type boundary conditions.  This includes
 * the nodal connectivities and degrees of freedom lists.  Derived classes must
 * implement shape functions
 */
class Neumann : public Boundary
{
public:
    Neumann(std::vector<List> const& nodal_connectivity,
            std::vector<List> const& dof_list,
            double const prescribed_value,
            bool const is_load_ramped);

    auto elements() const { return nodal_connectivity.size(); }

protected:
    std::vector<List> nodal_connectivity;
    std::vector<List> dof_list;
};
}
