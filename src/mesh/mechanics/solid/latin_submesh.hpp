
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
    using base = mechanics::solid::submesh;

public:
    // use the base class constructor
    using base::submesh;

    ~latin_submesh();

    latin_submesh(latin_submesh const&) = delete;

    latin_submesh(latin_submesh&&);

    latin_submesh& operator=(latin_submesh const&) = delete;

    latin_submesh& operator=(latin_submesh&&);

    /**
     * Compute the incremental latin internal force vector, infinite/vertical
     * search direction using the formula
     * \f{align*}{
     * f_{i} &= \int_{V} B^{T} (\sigma_{i}-\hat{\sigma}_{i}) dV
     * \f}
     * \return internal element force
     */
    [[nodiscard]] auto incremental_latin_internal_force(std::int32_t const element,
                                                        double const latin_search_direction) const
        -> vector const&;

protected:
    using base::bilinear_gradient;
    using base::coordinates;
};
}
}
