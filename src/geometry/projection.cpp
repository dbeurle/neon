
#include "projection.hpp"

#include <Eigen/Geometry>

namespace neon::geometry
{
matrix2x project_to_plane(matrix3x const& nodal_coordinates)
{
    auto const lnodes = nodal_coordinates.cols();

    matrix2x plane_coordinates = matrix::Zero(2, lnodes);

    auto const normal = unit_outward_normal(nodal_coordinates);

    // Orthogonal basis vectors
    vector3 e1, e2;

    // Algorithm implemented from Jeppe Revall Frisvad titled
    // Building an Orthonormal Basis from a 3D Unit Vector Without Normalization
    if (normal(2) <= -0.9999999)
    {
        // Handle the singularity
        e1 = -vector3::UnitY();
        e2 = -vector3::UnitX();
    }
    else
    {
        auto const a = 1.0 / (1.0 + normal(2));
        auto const b = -normal(0) * normal(1) * a;
        e1 << 1.0 - normal(0) * normal(0) * a, b, -normal(0);
        e2 << b, 1.0 - normal(1) * normal(1) * a, -normal(1);
    }
    // Project the points onto the plane using e1 and e2
    for (int lnode = 0; lnode < lnodes; ++lnode)
    {
        plane_coordinates(0, lnode) = e1.dot(nodal_coordinates.col(lnode));
        plane_coordinates(1, lnode) = e2.dot(nodal_coordinates.col(lnode));
    }
    return plane_coordinates;
}

vector3 unit_outward_normal(matrix3x const& nodal_coordinates)
{
    using namespace Eigen;
    Hyperplane<double, 3> hp = Hyperplane<double, 3>::Through(nodal_coordinates.col(0),
                                                              nodal_coordinates.col(1),
                                                              nodal_coordinates.col(2));
    return -hp.normal();
}
}
