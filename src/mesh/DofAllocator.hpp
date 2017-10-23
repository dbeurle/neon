
#pragma once

#include "numeric/DenseTypes.hpp"

#include <string>

namespace neon
{
namespace solid
{
/**
 * encode_dof takes a string valued representation of the DoF and returns the
 * encoding as an integer for computational use.  If name is not recognised
 * then an exception is thrown.  These functions are only usable for each
 * module
 */
int encode_dof(std::string const& name);
}

/**
 * Allocates the dof lists from the nodal connectivity vector
 * @param nodal_dofs Number of degrees of freedom for each node
 * @param nodal_connectivity Vector of nodal coordinates
 * @return The global degrees of freedom
 */
std::vector<List> allocate_dof_list(int const nodal_dofs,
                                    std::vector<List> const& nodal_connectivity);

/**
 * This function accepts the nodal connectivity of the mesh, multiplies each
 * entry by the nodal_dofs and adds the offset.  This function is intended for
 * use with the boundary classes, where each boundary class holds a dof_list
 * with the dofs associated only with the particular dof.
 */
std::vector<List> filter_dof_list(int const nodal_dofs,
                                  int const dof_offset,
                                  std::vector<List> const& nodal_connectivity);
}
