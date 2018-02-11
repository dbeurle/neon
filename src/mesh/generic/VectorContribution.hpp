
#pragma once

#include "Boundary.hpp"

#include "numeric/dense_matrix.hpp"
#include "numeric/index_types.hpp"

namespace neon
{
/**
 * VectorContribution contributes to the right hand side loading vector.
 * This is intended as a base class for distributed loads \sa Neumann and
 * point loads.
 */
class VectorContribution : public Boundary
{
public:
    using Boundary::Boundary;

    /** @return an element external force vector for a given element */
    [[nodiscard]] virtual std::tuple<local_indices const&, vector> external_force(
        std::int32_t const element, double const load_factor) const = 0;
};
}
