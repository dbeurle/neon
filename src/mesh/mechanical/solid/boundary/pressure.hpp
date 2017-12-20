
#pragma once

#include "traction.hpp"

namespace neon::mechanical::solid
{
/**
 * pressure computes the pressure load acting normal the quadrature point
 * on the surface of an element in the initial configuration.  This computes
 * the cross product of the shape functions that describe the surface element.
 * In the most general case, a pressure will contribute to three DoFs but
 * could also recover tractions if the surface is aligned with an axis.
 *
 * The convention used here is that a positive value represents compression
 * on the surface.
 */
class pressure : public traction
{
public:
    using traction::traction;

    std::tuple<List const&, vector> external_force(int const element,
                                                   double const load_factor) const override;
};
}
