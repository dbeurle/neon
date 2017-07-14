
#pragma once

#include <vector>

// Do not parallelise the GEMM routines
#ifndef NEON_PARALLEL_EIGEN_SOLVERS
#define EIGEN_DONT_PARALLELIZE
#endif

#include <Eigen/Dense>

namespace neon
{
using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using Matrix2 = Eigen::Matrix<double, 2, 2>;
using Matrix3 = Eigen::Matrix<double, 3, 3>;

using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using Vector2 = Eigen::Vector2d;
using Vector3 = Eigen::Vector3d;

using Array = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>;
using List = std::vector<int>;

using DofPair = std::pair<bool, double>;

using femValue = std::tuple<Vector, Matrix>;
}
