
#pragma once

#include <memory>
#include <utility>
#include <vector>

#include "numeric/DenseTypes.hpp"
#include "quadrature/NumericalQuadrature.hpp"

namespace neon
{
template <typename Quadrature>
class ShapeFunction
{
public:
    using Coordinates = typename Quadrature::Coordinate;

public:
    /** Construct the shape function by consuming a quadrature implementation */
    ShapeFunction(std::unique_ptr<Quadrature>&& quadratureImpl)
        : numerical_quadrature(std::move(quadratureImpl))
    {
    }

    virtual int nodes() const = 0;

    Quadrature const& quadrature() const { return *numerical_quadrature.get(); };

    /** Quadrature point to nodal point extrapolation matrix */
    Matrix const& local_quadrature_extrapolation() const { return extrapolation; }

protected:
    /**
     * Compute the extrapolation matrix to allow for quadrature valued variables
     * to be averaged to the nodal points without ill-effects when using a
     * least squares (for example with quadratric tetrahedron elements).
     *
     * Algorithm is taken from 'A local extrapolation method for finite elements'
     * by R. Durand and M.M. Farias, from Advances in Engineering Software, 2014
     */
    void compute_extrapolation_matrix(Matrix const N,
                                      Matrix const local_nodal_coordinates,
                                      Matrix const local_quadrature_coordinates,
                                      int const nodes_per_element);

protected:
    Matrix extrapolation; //!< Quadrature point to nodal point mapping

    std::unique_ptr<Quadrature> numerical_quadrature;
};

template <typename Quadrature>
void ShapeFunction<Quadrature>::compute_extrapolation_matrix(
    Matrix const N,
    Matrix const local_nodal_coordinates,
    Matrix const local_quadrature_coordinates,
    int const nodes_per_element)
{
    // Take short names for consistency with algorithm
    auto const n = nodes_per_element;
    auto const m = numerical_quadrature->points();

    auto const& xi = local_nodal_coordinates;
    auto const& xi_hat = local_quadrature_coordinates;

    if (m == n)
    {
        extrapolation = N.inverse();
        return;
    }

    Matrix const N_plus = (N.transpose() * N).inverse() * N.transpose();

    if (m > n)
    {
        extrapolation = N_plus;
        return;
    }

    // Number of quadrature points are less than the number of nodes
    auto const xi_hat_plus = (xi_hat.transpose() * xi_hat).inverse() * xi_hat.transpose();

    Matrix3 const I = Matrix::Identity(m, m);

    extrapolation = N_plus * (I - xi_hat * xi_hat_plus) + xi * xi_hat_plus;
}

using VolumeInterpolation = ShapeFunction<VolumeQuadrature>;
}
