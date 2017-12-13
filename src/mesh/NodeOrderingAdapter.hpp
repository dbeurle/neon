
#pragma once

#include "mesh/ElementTopology.hpp"
#include "numeric/IndexTypes.hpp"


#include <vtkCellType.h>

namespace neon
{
/**
 * NodeOrderingAdapter provides methods for converting the element nodal
 * connectivity ordering between different formats.
 */
class NodeOrderingAdapter
{
public:
    /**
     * Convert the nodal connectivities from the gmsh node ordering to the
     * reference ordering in \cite Hughes2012
     */
    void convert_from_gmsh(std::vector<List>& nodal_connectivity,
                           ElementTopology const topology) const;

    [[nodiscard]] std::vector<List> convert_to_vtk(std::vector<List> nodal_connectivity,
                                                   ElementTopology const element_topology) const;

    /** Add strong types to gmsh integer element codes */
    [[nodiscard]] ElementTopology gmsh_type_to_enum(int const element_code) const;

    [[nodiscard]] VTKCellType to_vtk(ElementTopology element_topology) const;
};
}
