
#pragma once

#include "mesh/ElementTopology.hpp"
#include "numeric/DenseTypes.hpp"

#include <unordered_map>

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
    /** Mutate the nodalConnectivity */
    void convert_from_gmsh(std::vector<List>& nodal_connectivity, ElementTopology topology);

    std::vector<List> convert_to_vtk(std::vector<List> nodal_connectivity,
                                     ElementTopology element_topology) const;

    /** Add strong types to gmsh integer element codes */
    ElementTopology gmsh_type_to_enum(int const element_code) const;

    int to_vtk(ElementTopology element_topology) const;

protected:
    std::unordered_map<int, ElementTopology>
        gmsh_converter{{2, ElementTopology::Triangle3},
                       {3, ElementTopology::Quadrilateral4},
                       {4, ElementTopology::Tetrahedron4},
                       {5, ElementTopology::Hexahedron8},
                       {6, ElementTopology::Prism6},
                       {9, ElementTopology::Triangle6},
                       {11, ElementTopology::Tetrahedron10}};

    std::unordered_map<ElementTopology, int>
        vtk_converter{{ElementTopology::Triangle3, VTK_TRIANGLE},
                      {ElementTopology::Quadrilateral4, VTK_QUAD},
                      {ElementTopology::Tetrahedron4, VTK_TETRA},
                      {ElementTopology::Hexahedron8, VTK_HEXAHEDRON},
                      {ElementTopology::Prism6, VTK_WEDGE},
                      {ElementTopology::Triangle6, VTK_QUADRATIC_TRIANGLE},
                      {ElementTopology::Tetrahedron10, VTK_QUADRATIC_TETRA}};
};
}
