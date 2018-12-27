
#include "interpolations/recovery_methods.hpp"

namespace neon
{
void local_extrapolation::compute_extrapolation_matrix(matrix const& N,
                                                       matrix const& local_nodal_coordinates,
                                                       matrix const& local_quadrature_coordinates)
{
    // Take short names for consistency with algorithm

    // Narrowing conversion but rows is expected to be greater than zero
    std::size_t const n = local_nodal_coordinates.rows();
    auto const m{m_quadrature->points()};

    assert(m == static_cast<std::size_t>(local_quadrature_coordinates.rows()));

    auto const& xi = local_nodal_coordinates;
    auto const& xi_hat = local_quadrature_coordinates;

    if (m == n)
    {
        extrapolation = N.inverse();
        return;
    }

    if (m > n)
    {
        extrapolation = (N.transpose() * N).inverse() * N.transpose();
        return;
    }

    matrix const N_plus = N.transpose() * (N * N.transpose()).inverse();

    // Number of quadrature points are less than the number of nodes
    auto const xi_hat_plus = xi_hat.transpose() * (xi_hat * xi_hat.transpose()).inverse();

    extrapolation = N_plus * (matrix::Identity(m, m) - xi_hat * xi_hat_plus) + xi * xi_hat_plus;
}
}
