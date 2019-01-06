
#include "interpolations/recovery_methods.hpp"

namespace neon
{
local_extrapolation::local_extrapolation(matrix const& shape_function,
                                         matrix const& local_nodal_coordinates,
                                         matrix const& local_quadrature_coordinates)
{
    this->compute_extrapolation_matrix(shape_function,
                                       local_nodal_coordinates,
                                       local_quadrature_coordinates);
}

void local_extrapolation::compute_extrapolation_matrix(matrix const& shape_function,
                                                       matrix const& local_nodal_coordinates,
                                                       matrix const& local_quadrature_coordinates)
{
    // Take short names for consistency with the algorithm presented in the paper

    // Narrowing conversion but rows is expected to be greater than zero
    std::size_t const n = local_nodal_coordinates.rows();
    std::size_t const m = local_quadrature_coordinates.rows();

    auto const& xi = local_nodal_coordinates;
    auto const& xi_hat = local_quadrature_coordinates;

    if (m == n)
    {
        extrapolation = shape_function.inverse();
        return;
    }

    if (m > n)
    {
        extrapolation = (shape_function.transpose() * shape_function).inverse()
                        * shape_function.transpose();
        return;
    }

    matrix const N_plus = shape_function.transpose()
                          * (shape_function * shape_function.transpose()).inverse();

    // Number of quadrature points are less than the number of nodes
    auto const xi_hat_plus = xi_hat.transpose() * (xi_hat * xi_hat.transpose()).inverse();

    extrapolation = N_plus * (matrix::Identity(m, m) - xi_hat * xi_hat_plus) + xi * xi_hat_plus;
}
}
