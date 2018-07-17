
#pragma once

#include <memory>
#include <utility>
#include <vector>

#include "numeric/dense_matrix.hpp"
#include "quadrature/numerical_quadrature.hpp"

namespace neon
{
/// Base class for a shape function.  This manages a pointer to the underlying
/// interpolation for a given quadrature type.
template <typename QuadratureType>
class shape_function
{
public:
    using quadrature_type = QuadratureType;

public:
    /// Construct the shape function by consuming a quadrature implementation
    shape_function(std::unique_ptr<quadrature_type>&& quadrature_impl)
        : numerical_quadrature(std::move(quadrature_impl))
    {
    }

    virtual ~shape_function() = default;

    /// \return Number of nodes in the interpolation function
    virtual int nodes() const = 0;

    quadrature_type const& quadrature() const { return *numerical_quadrature; };

    /** Quadrature point to nodal point extrapolation matrix */
    matrix const& local_quadrature_extrapolation() const { return extrapolation; }

protected:
    /// Compute the extrapolation matrix to allow for quadrature valued variables
    /// to be averaged to the nodal points without ill-effects when using a
    /// least squares (for example with quadratric tetrahedron elements)
    /// developed in \cite Durand2014
    void compute_extrapolation_matrix(matrix const N,
                                      matrix const local_nodal_coordinates,
                                      matrix const local_quadrature_coordinates);

protected:
    matrix extrapolation; //!< Quadrature point to nodal point mapping

    std::unique_ptr<quadrature_type> numerical_quadrature;
};

template <typename quadrature_t>
void shape_function<quadrature_t>::compute_extrapolation_matrix(matrix const N,
                                                                matrix const local_nodal_coordinates,
                                                                matrix const local_quadrature_coordinates)
{
    // Take short names for consistency with algorithm

    // Narrowing conversion but rows is expected to be greater than zero
    std::size_t const n = local_nodal_coordinates.rows();
    auto const m{numerical_quadrature->points()};

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

extern template class shape_function<numerical_quadrature<double>>;
extern template class shape_function<surface_quadrature>;
extern template class shape_function<volume_quadrature>;

using line_interpolation = shape_function<numerical_quadrature<double>>;
using surface_interpolation = shape_function<surface_quadrature>;
using volume_interpolation = shape_function<volume_quadrature>;
}
