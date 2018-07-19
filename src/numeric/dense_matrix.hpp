
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
/// Matrix in row major layout
using matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
/// Alias to row major layout matrix
using row_matrix = matrix;
/// Matrix in column major layout
using col_matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

/// 2x2 fixed size matrix
using matrix2 = Eigen::Matrix<double, 2, 2>;
/// 3x3 fixed size matrix
using matrix3 = Eigen::Matrix<double, 3, 3>;
/// 6x6 fixed size matrix
using matrix6 = Eigen::Matrix<double, 6, 6>;
/// 9x9 fixed size matrix
using matrix9 = Eigen::Matrix<double, 9, 9>;
/// 12x12 fixed size matrix
using matrix12 = Eigen::Matrix<double, 12, 12>;
/// 16x16 fixed size matrix
using matrix16 = Eigen::Matrix<double, 16, 16>;
/// 3x1 fixed size matrix for non-square Jacobians
using matrix31 = Eigen::Matrix<double, 3, 1>;
/// 3x2 fixed size matrix for non-square Jacobians
using matrix32 = Eigen::Matrix<double, 3, 2>;

template <int geometric_dimension>
using matrixxd = Eigen::Matrix<double, Eigen::Dynamic, geometric_dimension>;

template <int geometric_dimension>
using matrixdx = Eigen::Matrix<double, geometric_dimension, Eigen::Dynamic>;

/// Compile time fixed rows for nodal coordinates in two dimensions
using matrix2x = Eigen::Matrix<double, 2, Eigen::Dynamic>;
/// Compile time fixed rows for nodal coordinates in three dimensions
using matrix3x = Eigen::Matrix<double, 3, Eigen::Dynamic>;

/// Fixed size vector of variable length
using vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
/// Fixed size vector of length two
using vector2 = Eigen::Vector2d;
/// Fixed size vector of length three
using vector3 = Eigen::Vector3d;
/// Fixed size vector of length four
using vector4 = Eigen::Vector4d;
/// Fixed size vector of length five
using vector5 = Eigen::Matrix<double, 5, 1>;
/// Fixed size vector of length six
using vector6 = Eigen::Matrix<double, 6, 1>;
/// Fixed size vector of length nine
using vector9 = Eigen::Matrix<double, 9, 1>;
/// Fixed size vector of length sixteen
using vector16 = Eigen::Matrix<double, 16, 1>;
}
