
#pragma once

#include "numeric/tensor_operations.hpp"

namespace neon::mechanical
{
namespace detail
{
/** Compute the von Mises stress based on the reduced stress tensor */
[[nodiscard]] inline auto von_mises_stress(matrix2 const& a)
{
    return std::sqrt(3.0 / 2.0) * deviatoric(a).norm();
}

/** Compute the von Mises stress based on the full stress tensor */
[[nodiscard]] inline auto von_mises_stress(matrix3 const& a)
{
    return std::sqrt(3.0 / 2.0) * deviatoric(a).norm();
}

[[nodiscard]] inline matrix2 compute_cauchy_stress(double const G,
                                                   double const lambda_e,
                                                   matrix2 const& elastic_strain)
{
    return lambda_e * elastic_strain.trace() * matrix2::Identity() + 2.0 * G * elastic_strain;
}

[[nodiscard]] inline matrix3 compute_cauchy_stress(double const G,
                                                   double const lambda_e,
                                                   matrix3 const& elastic_strain)
{
    return lambda_e * elastic_strain.trace() * matrix3::Identity() + 2.0 * G * elastic_strain;
}
}

/** Compute the von Mises stress of the stress tensor */
template <typename MatrixExpression>
[[nodiscard]] inline auto von_mises_stress(MatrixExpression const& cauchy_stress)
{
    return detail::von_mises_stress(cauchy_stress.eval());
}

/** Compute the Cauchy stress of the stress tensor */
template <typename MatrixExpression>
[[nodiscard]] inline auto compute_cauchy_stress(double const G,
                                                double const lambda_e,
                                                MatrixExpression const elastic_strain)
{
    return detail::compute_cauchy_stress(G, lambda_e, elastic_strain.eval());
}
}
