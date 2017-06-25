
#include "mesh/NodeOrderingAdapter.hpp"

#include "PreprocessorExceptions.hpp"

#include <iostream>

namespace neon
{
void NodeOrderingAdapter::convert_from_gmsh(std::vector<List>& nodal_connectivity,
                                            ElementTopology element_topology)
{
    // Reorder based on the differences between the local node numbering
    // provided from Section 9.3 Node ordering
    // http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
    switch (element_topology)
    {
        case ElementTopology::Tetrahedron10:
        {
            for (auto& nodal_list : nodal_connectivity)
            {
                std::swap(nodal_list.at(4), nodal_list.at(9));
            }
            break;
        }
        case ElementTopology::Triangle6:
        {
            std::cout << "Fix the triangle reorder!!!" << std::endl;
            std::abort();
            break;
        }
        default:
            break;
    }
}

ElementTopology NodeOrderingAdapter::gmsh_type_to_enum(int element_code) const
{
    auto const found = gmsh_converter.find(element_code);
    if (found == gmsh_converter.end())
    {
        throw KeyNotFoundInMap<int>(element_code);
    }
    return found->second;
}

int NodeOrderingAdapter::to_vtk_type(ElementTopology element_topology) const
{
    auto const found = vtk_converter.find(element_topology);
    if (found == vtk_converter.end())
    {
        throw KeyNotFoundInMap<int>(static_cast<int>(element_topology));
    }
    return found->second;
}
}
