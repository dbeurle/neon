
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
using RowMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

using Matrix2 = Eigen::Matrix<double, 2, 2>;
using Matrix3 = Eigen::Matrix<double, 3, 3>;
using Matrix6 = Eigen::Matrix<double, 6, 6>;

using CMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, 6, 6>;

using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using Vector2 = Eigen::Vector2d;
using Vector3 = Eigen::Vector3d;
using Vector6 = Eigen::Matrix<double, 6, 1, Eigen::ColMajor, 6, 1>;

using Array = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>;
using List = std::vector<int>;

using femValue = std::tuple<Vector, Matrix>;

template <class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type is_approx(T const x,
                                                                                   T const y,
                                                                                   int const ulp = 2)
{
    // Taken and modified from
    // http://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
    // Since the numeric epsilon is defined at 1.0 then it must be scaled by
    // the worse case (x + y) and accounted for my the ULP (units in the last place).
    return std::abs(x - y) < std::numeric_limits<T>::epsilon() * std::abs(x + y) * ulp
           || std::abs(x - y) < std::numeric_limits<T>::min();
}
}
