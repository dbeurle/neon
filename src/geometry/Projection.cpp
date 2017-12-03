
#include "Projection.hpp"

#include <Eigen/Geometry>

namespace neon::geometry
{
Matrix2x project_to_plane(Matrix3x const& nodal_coordinates)
{
    auto const lnodes = nodal_coordinates.cols();

    Matrix2x plane_coordinates = Matrix::Zero(2, lnodes);

    auto const normal = unit_outward_normal(nodal_coordinates);

    // Orthogonal basis vectors
    Vector3 e1, e2;

    // Algorithm implemented from Jeppe Revall Frisvad titled
    // Building an Orthonormal Basis from a 3D Unit Vector Without Normalization
    if (normal(2) <= -0.9999999)
    {
        // Handle the singularity
        e1 = -Vector3::UnitY();
        e2 = -Vector3::UnitX();
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

Vector3 unit_outward_normal(Matrix3x const& nodal_coordinates)
{
    using namespace Eigen;
    Hyperplane<double, 3> hp = Hyperplane<double, 3>::Through(nodal_coordinates.col(0),
                                                              nodal_coordinates.col(1),
                                                              nodal_coordinates.col(2));
    return -hp.normal();
}
}
