
#pragma once

#include <memory>
#include <utility>
#include <vector>

#include "numeric/dense_matrix.hpp"
#include "quadrature/numerical_quadrature.hpp"

namespace neon
{
/// Base class for a shape function.  This manages a pointer to the underlying
/// quadrature scheme and provides an interface for the interpolation function
/// properties.
template <typename QuadratureType>
class shape_function
{
public:
    using quadrature_type = QuadratureType;

public:
    /// Construct the shape function by consuming a quadrature implementation
    shape_function(std::unique_ptr<quadrature_type>&& quadrature_impl, std::uint8_t node_count)
        : m_quadrature(std::move(quadrature_impl)), m_node_count{node_count}
    {
    }

    virtual ~shape_function() = default;

    /// \return Number of nodes in the element
    auto number_of_nodes() const -> std::uint8_t { return m_node_count; }
    /// \return Highest polynomial order in interpolation function
    auto polynomial_order() const -> std::uint8_t { return m_polynomial_order; }
    /// \return Highest monomial order in interpolation function
    auto monomial_order() const -> std::uint8_t { return m_monomial_order; }

    quadrature_type const& quadrature() const { return *m_quadrature; };

    /// \return  Quadrature point to nodal point extrapolation matrix
    matrix const& local_quadrature_extrapolation() const { return extrapolation; }

protected:
    /// Compute the extrapolation matrix to allow for quadrature valued variables
    /// to be averaged to the nodal points without ill-effects when using a
    /// least squares (for example with quadratric tetrahedron elements)
    /// developed in \cite Durand2014
    void compute_extrapolation_matrix(matrix const& N,
                                      matrix const& local_nodal_coordinates,
                                      matrix const& local_quadrature_coordinates);

protected:
    /// Quadrature point to nodal point mapping
    matrix extrapolation;

    /// Pointer to numerical quadrature scheme
    std::unique_ptr<quadrature_type> m_quadrature;

    /// Nodes per element
    std::uint8_t m_node_count{0};
    /// Highest order of polynomral term
    std::uint8_t m_polynomial_order{0};
    /// Highest order of monomial term
    std::uint8_t m_monomial_order{0};
};

template <typename QuadratureType>
void shape_function<QuadratureType>::compute_extrapolation_matrix(
    matrix const& N,
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

extern template class shape_function<numerical_quadrature<double>>;
extern template class shape_function<surface_quadrature>;
extern template class shape_function<volume_quadrature>;

using line_interpolation = shape_function<numerical_quadrature<double>>;
using surface_interpolation = shape_function<surface_quadrature>;
using volume_interpolation = shape_function<volume_quadrature>;
}
