
#pragma once

#include "boundary_condition.hpp"

#include "numeric/dense_matrix.hpp"
#include "numeric/index_types.hpp"

namespace neon
{
/// vector_contribution contributes to the right hand side loading vector.
/// This is intended as a base class for distributed loads \sa neumann and
/// point loads.
class vector_contribution : public boundary_condition
{
public:
    using boundary::boundary;

    /// \return element external force vector
    [[nodiscard]] virtual auto external_force(std::int64_t const element,
                                              double const load_factor) const -> vector const& = 0;
};
}
