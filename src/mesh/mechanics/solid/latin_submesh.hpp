
#pragma once

/// @file

#include "mesh/mechanics/solid/submesh.hpp"

namespace neon
{
namespace mechanics::solid
{
/// latin_submesh provides the element local routines for computing the system
/// components for a three-dimensional continuum mechanics discretisation.
class latin_submesh : public mechanics::solid::submesh
{
public:
    // use the base class constructor
    using mechanics::solid::submesh::submesh;

    /**
     * Compute the incremental latin internal force vector, infinite/vertical search direction,
     using the formula
     * \f{align*}{
     * f_{i} &= \int_{V} B^{T} (\sigma_{i}-\hat{\sigma}_{i}) dV
     * \f}
     * \return internal element force
     */
    [[nodiscard]] std::pair<index_view, vector const&> incremental_latin_internal_force(
        std::int32_t const element) const;
};

}
}
