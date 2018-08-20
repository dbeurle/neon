
#pragma once

#include "numeric/dense_matrix.hpp"

#include <cmath>

namespace neon::mechanics
{
/**
 * \fn Compute the deformed tangent using the unimodular deformation gradient
 * and the vector associated with the quadrature point on the unit sphere
 */
[[nodiscard]] inline vector3 deformed_tangent(matrix3 const& F_unimodular,
                                              vector3 const& surface_vector)
{
    return F_unimodular * surface_vector;
}

/** \fn Compute the microstretch, which is the norm of the deformed tangent vector */
[[nodiscard]] inline auto compute_microstretch(vector3 const& deformed_tangent)
{
    return deformed_tangent.norm();
}

[[nodiscard]] inline vector3 deformed_normal(matrix3 const& F_unimodular,
                                             vector3 const& surface_vector)
{
    return F_unimodular.inverse().transpose() * surface_vector;
}

[[nodiscard]] inline auto compute_area_stretch(vector3 const& deformed_normal)
{
    return deformed_normal.norm();
}

/**
 * \fn Compute the Padé approximation of the inverse Langevin stretch model
 * \f{align*}{
     n \psi_f^{'}(\lambda) &= \frac{3N - \lambda^2}{N - \lambda^2}
   \f}
 * where \f$ N \f$ is the average number of segments per chain.
 *
 * \param micro_stretch The chain microstretch on the unit sphere
 * \param N The average number of segments per chain
 * \return The Padé approximation
 */
[[nodiscard]] inline double pade_first(double const micro_stretch, double const N)
{
    return (3.0 * N - std::pow(micro_stretch, 2)) / (N - std::pow(micro_stretch, 2));
}

/**
 * \fn
 * Compute the Padé approximation of the inverse Langevin stretch model
 * \f{align*}{
     n \psi_f^{''}(\lambda) &= \frac{\lambda^4 + 3N^2}{(N - \lambda^2)^2}
   \f}
 * where \f$ N \f$ is the average number of segments per chain.
 *
 * \param micro_stretch The chain microstretch on the unit sphere
 * \param N The average number of segments per chain
 * \return The Padé approximation
 */
[[nodiscard]] inline double pade_second(double const micro_stretch, double const N)
{
    return (std::pow(micro_stretch, 4) + 3.0 * std::pow(N, 2))
           / std::pow(N - std::pow(micro_stretch, 2), 2);
}
}
