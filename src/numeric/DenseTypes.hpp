
#pragma once

// Do not parallelise the GEMM routines as this is called by mulitple threads
#ifndef NEON_PARALLEL_EIGEN_SOLVERS
#define EIGEN_DONT_PARALLELIZE
#endif

#include <Eigen/Dense>

namespace neon
{
using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using RowMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

using Matrix2 = Eigen::Matrix<double, 2, 2>;
using Matrix3 = Eigen::Matrix<double, 3, 3>;
using Matrix6 = Eigen::Matrix<double, 6, 6>;
using Matrix9 = Eigen::Matrix<double, 9, 9>;

using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using Vector2 = Eigen::Vector2d;
using Vector3 = Eigen::Vector3d;
using Vector6 = Eigen::Matrix<double, 6, 1, Eigen::ColMajor, 6, 1>;

using Array = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>;
}
