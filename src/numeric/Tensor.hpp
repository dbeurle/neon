
#pragma once

#include "DenseTypes.hpp"

namespace neon
{
inline auto double_dot(Matrix3 const& a, Matrix3 const& b) { return (a.array() * b.array()).sum(); }

/** @return the deviatoric part of the tensor */
inline auto deviatoric(Matrix3 const& a) { return a - Matrix3::Identity() * a.trace() / 3.0; }

inline auto symmetric(Matrix3 const& a) { return 1.0 / 2.0 * (a.transpose() + a); }

inline auto voigt(Matrix3 const& a)
{
    Eigen::Matrix<double, 6, 1> b;
    b << a(0, 0), a(1, 1), a(2, 2), a(1, 2), a(0, 2), a(0, 1);
    return b;
}

inline Matrix3 voigt_to_matrix(Vector const& a)
{
    Matrix3 b;
    b << a(0), a(5), a(4), //
        a(5), a(1), a(3),  //
        a(4), a(3), a(2);
    return b;
}

/**
 * I1 returns the coefficient I1, the first stress invariant,
 * which is equal to the trace
 * @return First invariant
 */
inline auto I1(Matrix3 const& a) { return a.trace(); }

/**
 * I2 returns the coefficient I2, the second stress invariant,
 * which is calculated by:
 * \f{align*}{
 * I_2 &= \frac{1}{2} \left( (tr \tau)^2 - tr(\tau \tau) \right)
 * \f}
 * @return Second invariant
 */
inline auto I2(Matrix3 const& a) { return 0.5 * (std::pow(a.trace(), 2) - (a * a).trace()); }

/** @return Third invariant, which is the determinant of the tensor */
inline auto I3(Matrix3 const& a) { return a.determinant(); }

inline auto identity_expansion(Matrix const& H, int const nodal_dofs, int const local_nodes)
{
    // Create the geometric part of the tangent stiffness matrix
    Matrix K = Matrix::Zero(local_nodes * nodal_dofs, local_nodes * nodal_dofs);
    for (auto i = 0; i < local_nodes; ++i)
        for (auto j = 0; j < local_nodes; ++j)
            for (auto k = 0; k < nodal_dofs; ++k)
                K(i * nodal_dofs + k, j * nodal_dofs + k) = H(i, j);
    return K;
}
}
