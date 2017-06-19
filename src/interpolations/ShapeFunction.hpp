
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

protected:
    Matrix interpolation; //!< Quadrature point to nodal point mapping

    std::unique_ptr<Quadrature> numerical_quadrature;
};

using VolumeInterpolation = ShapeFunction<VolumeQuadrature>;
}
