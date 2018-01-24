
#pragma once

#include "mesh/element_topology.hpp"
#include "numeric/IndexTypes.hpp"

#include <vtkCellType.h>

namespace neon
{
/**
 * Convert the nodal connectivities from the gmsh node ordering to the
 * reference ordering in \cite Hughes2012
 */
void convert_from_gmsh(std::vector<List>& nodal_connectivity, element_topology const topology);

[[nodiscard]] std::vector<List> convert_to_vtk(std::vector<List> nodal_connectivity,
                                               element_topology const topology);

/** Add strong types to gmsh integer element codes */
[[nodiscard]] element_topology gmsh_type_to_enum(int const element_code);

[[nodiscard]] VTKCellType to_vtk(element_topology const topology);
}
