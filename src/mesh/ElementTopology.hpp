
#pragma once

namespace neon
{
/** Enumeration for the element shapes and follows the enumeration by GMSH */
enum class ElementTopology {
    Invalid = 0,
    Line2,
    Triangle3,
    Quadrilateral4,
    Tetrahedron4,
    Hexahedron8,
    Prism6,
    Pyramid5,
    Line3,
    Triangle6,
    Quadrilateral9, // 4 vertex, 4 edges and 1 face node
    Tetrahedron10,
    Hexahedron27,
    Prism18,
    Pyramid14,
    Point = 15,
    Quadrilateral8,
    Hexahedron20,
    Prism15,
    Pyramid13,
    Triangle9 = 20,
    Triangle10,
    Triangle12,
    Triangle15,
    Triangle15_IC, // Incomplete 15 node triangle
    Triangle21 = 25,
    Edge4,
    Edge5,
    Edge6,
    Tetrahedron20,
    Tetrahedron35,
    Tetrahedron56,
    Hexahedron64 = 92,
    Hexahedron125
};
}
