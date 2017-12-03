
#pragma once

#include "numeric/DenseMatrix.hpp"

namespace neon
{
namespace detail
{
/** Compute the von Mises stress based on the reduced stress tensor */
[[nodiscard]] inline double von_mises_stress(Matrix2 const& a)
{
    return std::sqrt(3.0 / 2.0) * deviatoric(a).norm();
}

/** Compute the von Mises stress based on the full stress tensor */
[[nodiscard]] inline double von_mises_stress(Matrix3 const& a)
{
    return std::sqrt(3.0 / 2.0) * deviatoric(a).norm();
}

[[nodiscard]] Matrix2 compute_cauchy_stress(double const G,
                                            double const lambda_e,
                                            Matrix2 const& elastic_strain)
{
    return lambda_e * elastic_strain.trace() * Matrix2::Identity() + 2.0 * G * elastic_strain;
}

[[nodiscard]] Matrix3 compute_cauchy_stress(double const G,
                                            double const lambda_e,
                                            Matrix3 const& elastic_strain)
{
    return lambda_e * elastic_strain.trace() * Matrix3::Identity() + 2.0 * G * elastic_strain;
}
}

/**
 * Compute the von Mises stress of the stress tensor
 */
template <typename MatrixExpression>
[[nodiscard]] inline double von_mises_stress(MatrixExpression const& cauchy_stress)
{
    return detail::von_mises_stress(cauchy_stress.eval());
}

template <typename MatrixExpression>
[[nodiscard]] auto compute_cauchy_stress(double const G,
                                         double const lambda_e,
                                         MatrixExpression const& elastic_strain)
{
    return detail::compute_cauchy_stress(G, lambda_e, elastic_strain.eval());
}
}
