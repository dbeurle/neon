
#pragma once

#include "mesh/ElementTopology.hpp"
#include "numeric/IndexTypes.hpp"

#include <vtkCellType.h>

namespace neon
{
/**
 * Convert the nodal connectivities from the gmsh node ordering to the
 * reference ordering in \cite Hughes2012
 */
void convert_from_gmsh(std::vector<std::vector<int64>>& nodal_connectivity,
                       ElementTopology const topology);

[[nodiscard]] std::vector<std::vector<int64>> convert_to_vtk(
    std::vector<std::vector<int64>> nodal_connectivity, ElementTopology const element_topology);

/** Add strong types to gmsh integer element codes */
[[nodiscard]] ElementTopology gmsh_type_to_enum(int const element_code);

[[nodiscard]] VTKCellType to_vtk(ElementTopology const element_topology);
}
