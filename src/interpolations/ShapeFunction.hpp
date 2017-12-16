
#pragma once

#include <memory>
#include <utility>
#include <vector>

#include "numeric/DenseMatrix.hpp"
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
     * least squares (for example with quadratric tetrahedron elements)
     * developed in \cite Durand2014
     */
    void compute_extrapolation_matrix(Matrix const N,
                                      Matrix const local_nodal_coordinates,
                                      Matrix const local_quadrature_coordinates);

protected:
    Matrix extrapolation; //!< Quadrature point to nodal point mapping

    std::unique_ptr<Quadrature> numerical_quadrature;
};

template <typename Quadrature>
void ShapeFunction<Quadrature>::compute_extrapolation_matrix(Matrix const N,
                                                             Matrix const local_nodal_coordinates,
                                                             Matrix const local_quadrature_coordinates)
{
    // Take short names for consistency with algorithm
    auto const n = local_nodal_coordinates.rows();
    auto const m = numerical_quadrature->points();

    assert(m == local_quadrature_coordinates.rows());

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

    Matrix const N_plus = N.transpose() * (N * N.transpose()).inverse();

    // Number of quadrature points are less than the number of nodes
    auto const xi_hat_plus = xi_hat.transpose() * (xi_hat * xi_hat.transpose()).inverse();

    extrapolation = N_plus * (Matrix::Identity(m, m) - xi_hat * xi_hat_plus) + xi * xi_hat_plus;
}

using LineInterpolation = ShapeFunction<NumericalQuadrature<double>>;
using SurfaceInterpolation = ShapeFunction<SurfaceQuadrature>;
using VolumeInterpolation = ShapeFunction<VolumeQuadrature>;
}
