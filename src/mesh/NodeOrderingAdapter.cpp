
#include "mesh/NodeOrderingAdapter.hpp"

#include "Exceptions.hpp"

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
                std::swap(nodal_list.at(0), nodal_list.at(3));
                std::swap(nodal_list.at(4), nodal_list.at(9));
            }
            break;
        }
        case ElementTopology::Hexahedron20:
        {
            /*
              Gmsh ordering (1 based indexing)
               4----14----3
               |\         |\
               | 16       | 15
              10  \       12 \
               |   8----20+---7
               |   |      |   |
               1---+-9----2   |
                \  18      \  19
                11 |        13|
                  \|         \|
                   5----17----6

              Hughes ordering (1 based indexing)

               4----11----3
               |\         |\
               | 20       | 19
              12  \       10 \
               |   8----15+---7
               |   |      |   |
               1---+-9----2   |
                \  16      \  14
                17 |        18|
                  \|         \|
                   5----13----6
            */

            for (auto& nodal_list : nodal_connectivity)
            {
                std::swap(nodal_list.at(11), nodal_list.at(9));
                std::swap(nodal_list.at(13), nodal_list.at(10));

                std::swap(nodal_list.at(12), nodal_list.at(17));
                std::swap(nodal_list.at(16), nodal_list.at(12));
                std::swap(nodal_list.at(16), nodal_list.at(13));

                std::swap(nodal_list.at(13), nodal_list.at(15));
                std::swap(nodal_list.at(13), nodal_list.at(19));

                std::swap(nodal_list.at(13), nodal_list.at(18));
                std::swap(nodal_list.at(14), nodal_list.at(18));
            }
            break;
        }
        default:
            break;
    }
}

std::vector<List> NodeOrderingAdapter::convert_to_vtk(std::vector<List> nodal_connectivity,
                                                      ElementTopology element_topology) const
{
    switch (element_topology)
    {
        case ElementTopology::Tetrahedron10:
        {
            for (auto& nodal_list : nodal_connectivity)
            {
                std::swap(nodal_list.at(6), nodal_list.at(8));
                std::swap(nodal_list.at(8), nodal_list.at(9));
            }
            break;
        }
        default:
            break;
    }
    return nodal_connectivity;
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

int NodeOrderingAdapter::to_vtk(ElementTopology element_topology) const
{
    auto const found = vtk_converter.find(element_topology);
    if (found == vtk_converter.end())
    {
        throw KeyNotFoundInMap<int>(static_cast<int>(element_topology));
    }
    return found->second;
}
}
