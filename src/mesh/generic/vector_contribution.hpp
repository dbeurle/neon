
#pragma once

#include "interpolator.hpp"

#include "numeric/DenseMatrix.hpp"
#include "numeric/IndexTypes.hpp"

namespace neon::boundary
{
/**
 * VectorContribution contributes to the right hand side loading vector.
 * This is intended as a base class for distributed loads \sa Neumann and
 * point loads.
 */
class vector_contribution : public interpolator
{
public:
    using interpolator::interpolator;

    /** @return an element external force vector for a given element */
    [[nodiscard]] virtual std::pair<std::vector<int64> const&, vector> external_force(
        int64 const element, double const load_factor) const = 0;
};
}
