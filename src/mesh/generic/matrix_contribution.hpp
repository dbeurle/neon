
#pragma once

#include "interpolator.hpp"

namespace neon
{
/**
 * matrix_contribution is a base class for including boundary conditions which
 * contribute to the stiffness matrix in the finite element discretisation.
 */
class matrix_contribution : public interpolator
{
public:
    using interpolator::interpolator;

    /** @return an element external force vector for a given element */
    [[nodiscard]] virtual std::pair<std::vector<int64> const&, matrix> external_stiffness(
        int64 const element, double const load_factor) const = 0;
};
}
