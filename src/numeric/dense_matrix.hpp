
#pragma once

// Do not parallelise the GEMM routines as this is called by multiple threads
#ifndef NEON_PARALLEL_EIGEN_SOLVERS
#define EIGEN_DONT_PARALLELIZE
#endif

#include <Eigen/Core>
#include <Eigen/LU>

/// \file dense_matrix.hpp

namespace neon
{
using matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using row_matrix = matrix;
using col_matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

/// 2x2 fixed size matrix
using matrix2 = Eigen::Matrix<double, 2, 2>;
/// 3x3 fixed size matrix
using matrix3 = Eigen::Matrix<double, 3, 3>;
/// 6x6 fixed size matrix
using matrix6 = Eigen::Matrix<double, 6, 6>;
/// 9x9 fixed size matrix
using matrix9 = Eigen::Matrix<double, 9, 9>;
/// 16x16 fixed size matrix
using matrix16 = Eigen::Matrix<double, 16, 16>;
/// 3x2 fixed size matrix for non-square Jacobians
using matrix32 = Eigen::Matrix<double, 3, 2>;

template <int geometric_dimension>
using matrixxd = Eigen::Matrix<double, Eigen::Dynamic, geometric_dimension>;

template <int geometric_dimension>
using matrixdx = Eigen::Matrix<double, geometric_dimension, Eigen::Dynamic>;

/** Compile time fixed rows for nodal coordinates in two dimensions */
using matrix2x = Eigen::Matrix<double, 2, Eigen::Dynamic>;
/** Compile time fixed rows for nodal coordinates in three dimensions */
using matrix3x = Eigen::Matrix<double, 3, Eigen::Dynamic>;

using vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using vector2 = Eigen::Vector2d;
using vector3 = Eigen::Vector3d;
using vector6 = Eigen::Matrix<double, 6, 1>;
using vector16 = Eigen::Matrix<double, 16, 1>;
}
