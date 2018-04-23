
#include "mesh/node_ordering_adapter.hpp"

#include "Exceptions.hpp"

#include <unordered_map>

namespace neon
{
std::unordered_map<int, element_topology> const gmsh_converter{{1, element_topology::line2},
                                                               {2, element_topology::triangle3},
                                                               {3, element_topology::quadrilateral4},
                                                               {4, element_topology::tetrahedron4},
                                                               {5, element_topology::hexahedron8},
                                                               {6, element_topology::prism6},
                                                               {8, element_topology::line3},
                                                               {9, element_topology::triangle6},
                                                               {11, element_topology::tetrahedron10},
                                                               {10, element_topology::quadrilateral9},
                                                               {12, element_topology::hexahedron27},
                                                               {13, element_topology::prism18},
                                                               {15, element_topology::point},
                                                               {16, element_topology::quadrilateral8},
                                                               {18, element_topology::prism15},
                                                               {17, element_topology::hexahedron20}};

std::unordered_map<element_topology, VTKCellType> const
    vtk_converter{{element_topology::triangle3, VTK_TRIANGLE},
                  {element_topology::quadrilateral4, VTK_QUAD},
                  {element_topology::quadrilateral8, VTK_QUADRATIC_QUAD},
                  {element_topology::quadrilateral9, VTK_BIQUADRATIC_QUAD},
                  {element_topology::tetrahedron4, VTK_TETRA},
                  {element_topology::hexahedron8, VTK_HEXAHEDRON},
                  {element_topology::prism6, VTK_WEDGE},
                  {element_topology::triangle6, VTK_QUADRATIC_TRIANGLE},
                  {element_topology::tetrahedron10, VTK_QUADRATIC_TETRA},
                  {element_topology::prism15, VTK_QUADRATIC_WEDGE},
                  // {element_topology::prism18, VTK_TRIQUADRATIC_WEDGE},
                  {element_topology::hexahedron20, VTK_QUADRATIC_HEXAHEDRON},
                  {element_topology::hexahedron27, VTK_TRIQUADRATIC_HEXAHEDRON}};

void convert_from_gmsh(indices& nodal_connectivity, element_topology const topology)
{
    // Reorder based on the differences between the local node numbering
    // provided from Section 9.3 Node ordering
    // http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
    switch (topology)
    {
        case element_topology::tetrahedron10:
        {
            nodal_connectivity.row(4).swap(nodal_connectivity.row(7));
            nodal_connectivity.row(4).swap(nodal_connectivity.row(5));
            nodal_connectivity.row(5).swap(nodal_connectivity.row(8));
            nodal_connectivity.row(8).swap(nodal_connectivity.row(9));
            nodal_connectivity.row(6).swap(nodal_connectivity.row(9));
            [[fallthrough]];
        }
        case element_topology::tetrahedron4:
        {
            nodal_connectivity.row(0).swap(nodal_connectivity.row(3));
            nodal_connectivity.row(0).swap(nodal_connectivity.row(2));
            nodal_connectivity.row(0).swap(nodal_connectivity.row(1));

            break;
        }
        case element_topology::prism6:
        {
            nodal_connectivity.row(0).swap(nodal_connectivity.row(1));
            nodal_connectivity.row(3).swap(nodal_connectivity.row(4));

            break;
        }
        case element_topology::prism15:
        {
            // -1 face
            nodal_connectivity.row(0).swap(nodal_connectivity.row(1));
            nodal_connectivity.row(3).swap(nodal_connectivity.row(6));
            nodal_connectivity.row(7).swap(nodal_connectivity.row(4));
            nodal_connectivity.row(5).swap(nodal_connectivity.row(9));

            // mid face
            nodal_connectivity.row(8).swap(nodal_connectivity.row(7));
            nodal_connectivity.row(10).swap(nodal_connectivity.row(6));
            nodal_connectivity.row(11).swap(nodal_connectivity.row(8));

            // +1 face
            nodal_connectivity.row(11).swap(nodal_connectivity.row(9));

            break;
        }
        case element_topology::hexahedron27:
        {
            /*
            Gmsh ordering (0 based indexing) taken from gmsh.info

                3----13----2
                |\         |\
                |15    24  | 14
                9  \ 20    11 \
                |   7----19+---6
                |22 |  26  | 23|
                0---+-8----1   |
                 \ 17    25 \  18
                 10 |  21    12|
                   \|         \|
                    4----16----5

            Hughes ordering (0 based indexing)

                3----10----2
                |\         |\
                | 19   23  | 18
               11  \ 20    9  \
                |   7----14+---6
                |24 |  26  | 25|
                0---+-8----1   |
                 \  15   21 \  13
                 16 |  22    17|
                   \|         \|
                    4----12----5

            */

            nodal_connectivity.row(21).swap(nodal_connectivity.row(25));
            nodal_connectivity.row(25).swap(nodal_connectivity.row(22));
            nodal_connectivity.row(24).swap(nodal_connectivity.row(25));
            nodal_connectivity.row(23).swap(nodal_connectivity.row(25));

            [[fallthrough]];
        }
        case element_topology::hexahedron20:
        {
            /*
              Gmsh ordering (0 based indexing) taken from gmsh.info

               3----13----2
               |\         |\
               | 15       | 14
               9  \       11 \
               |   7----19+---6
               |   |      |   |
               0---+-8----1   |
                \  17      \  18
                10 |        12|
                  \|         \|
                   4----16----5

              Hughes ordering (0 based indexing)

               3----10----2
               |\         |\
               | 19       | 18
              11  \       9  \
               |   7----14+---6
               |   |      |   |
               0---+-8----1   |
                \  15      \  13
                16 |        17|
                  \|         \|
                   4----12----5
            */

            nodal_connectivity.row(11).swap(nodal_connectivity.row(9));
            nodal_connectivity.row(13).swap(nodal_connectivity.row(10));

            nodal_connectivity.row(12).swap(nodal_connectivity.row(17));
            nodal_connectivity.row(16).swap(nodal_connectivity.row(12));
            nodal_connectivity.row(16).swap(nodal_connectivity.row(13));

            nodal_connectivity.row(13).swap(nodal_connectivity.row(15));
            nodal_connectivity.row(13).swap(nodal_connectivity.row(19));

            nodal_connectivity.row(13).swap(nodal_connectivity.row(18));
            nodal_connectivity.row(14).swap(nodal_connectivity.row(18));

            break;
        }
        default:
            break;
    }
}

indices convert_to_vtk(indices nodal_connectivity, element_topology const topology)
{
    switch (topology)
    {
        case element_topology::tetrahedron4:
        {
            nodal_connectivity.row(0).swap(nodal_connectivity.row(1));

            break;
        }
        case element_topology::tetrahedron10:
        {
            nodal_connectivity.row(6).swap(nodal_connectivity.row(8));
            nodal_connectivity.row(8).swap(nodal_connectivity.row(9));

            break;
        }
        case element_topology::hexahedron20:
        {
            // The ordering of the twenty points defining the cell is point ids (0-7,8-19) where
            // point ids 0-7 are the eight corner vertices of the cube; followed by twelve midedge
            // nodes (8-19). Note that these midedge nodes correspond lie on the edges defined by
            // 8 > (0,1), 9 > (1,2), 10 > (2,3), 11 > (3,0),
            // 12 > (4,5), 13 > (5,6), 14 > (6,7), 15 > (7,4),
            // 16 > (0,4), 17 > (1,5), 18 > (2,6), 19 > (3,7).

            // NOTE This corresponds nicely to the Hughes ordering
            break;
        }
        case element_topology::hexahedron27:
        {
            /* top
             *  7--14--6
             *  |      |
             * 15  25  13
             *  |      |
             *  4--12--5
             *
             *  middle
             * 19--23--18
             *  |      |
             * 20  26  21
             *  |      |
             * 16--22--17
             *
             * bottom
             *  3--10--2
             *  |      |
             * 11  24  9
             *  |      |
             *  0-- 8--1
             */
            nodal_connectivity.row(21).swap(nodal_connectivity.row(25));
            nodal_connectivity.row(20).swap(nodal_connectivity.row(24));

            break;
        }
        default:
            break;
    }
    return nodal_connectivity;
}

element_topology gmsh_type_to_enum(std::int32_t const element_code)
{
    auto const found = gmsh_converter.find(element_code);
    if (found == gmsh_converter.end())
    {
        throw std::domain_error("Element code " + std::to_string(element_code)
                                + " not implemented for gmsh element type");
    }
    return found->second;
}

VTKCellType to_vtk(element_topology const topology)
{
    auto const found = vtk_converter.find(topology);
    if (found == vtk_converter.end())
    {
        throw std::domain_error("Element code " + std::to_string(static_cast<int>(topology))
                                + " not implemented for vtk element type");
    }
    return found->second;
}
} // namespace neon
