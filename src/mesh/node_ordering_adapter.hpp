
#pragma once

/// @file

#include "mesh/element_topology.hpp"
#include "numeric/index_types.hpp"

#include <vtkCellType.h>

namespace neon
{
/// Convert the nodal connectivities from the gmsh node ordering to the
/// reference ordering in @cite Hughes2012
///
/// \param [in, out] nodal_connectivity Local element nodal connectivity data
/// \param [in] topology The element topology
void convert_from_gmsh(indices& nodal_connectivity, element_topology const topology);

/// Convert the \p neon element node ordering to the VTK node ordering
/// \param nodal_connectivity in \p neon format
/// \param topology Element topology
/// \return Nodal connectivity in VTK format for the given element topology
[[nodiscard]] indices convert_to_vtk(indices nodal_connectivity, element_topology const topology);

/// Add strong types to gmsh integer element codes
[[nodiscard]] element_topology gmsh_type_to_enum(std::int32_t const element_code);

/// Convert the \p neon topology type to the \p VTKCellType
[[nodiscard]] VTKCellType to_vtk(element_topology const topology);
}
