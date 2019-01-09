
#pragma once

/// @file

namespace neon
{
/// Enumeration for the element shapes and follows the enumeration by GMSH
enum class element_topology {
    invalid = 0,
    line2,
    triangle3,
    quadrilateral4,
    tetrahedron4,
    hexahedron8,
    prism6,
    pyramid5,
    line3,
    triangle6,
    /// 4 vertex, 4 edges and 1 face node
    quadrilateral9,
    tetrahedron10,
    hexahedron27,
    prism18,
    pyramid14,
    point = 15,
    quadrilateral8,
    hexahedron20,
    prism15,
    pyramid13,
    triangle9 = 20,
    triangle10,
    triangle12,
    triangle15,
    /// Incomplete 15 node triangle
    triangle15_IC,
    triangle21 = 25,
    edge4,
    edge5,
    edge6,
    tetrahedron20,
    tetrahedron35,
    tetrahedron56,
    hexahedron64 = 92,
    hexahedron125
};
}
