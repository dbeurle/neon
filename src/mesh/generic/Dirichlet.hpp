
#pragma once

#include "mesh/generic/Boundary.hpp"

#include "numeric/index_types.hpp"

namespace neon
{
class Dirichlet : public Boundary
{
public:
    explicit Dirichlet(local_indices dofs, json const& times, json const& loads);

    [[nodiscard]] auto const& dof_view() const noexcept { return dofs; }

    /** Get the value depending on the loading factor */
    [[nodiscard]] auto value_view(double const load_factor = 1.0) const
    {
        return interpolate_prescribed_load(load_factor);
    }

protected:
    local_indices dofs;
};
}
