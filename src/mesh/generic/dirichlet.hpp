
#pragma once

#include "mesh/generic/boundary.hpp"
#include "numeric/index_types.hpp"

namespace neon
{
class dirichlet : public boundary
{
public:
    explicit dirichlet(std::vector<std::int32_t> unique_dofs, json const& times, json const& loads);

    explicit dirichlet(std::vector<std::int32_t> unique_dofs,
                       json const& boundary_data,
                       std::string const& name,
                       double const generate_time_step);

    [[nodiscard]] auto const& dof_view() const noexcept { return unique_dofs; }

    /** Get the value depending on the loading factor */
    [[nodiscard]] auto value_view(double const load_factor = 1.0) const
    {
        return interpolate_prescribed_load(load_factor);
    }

protected:
    std::vector<std::int32_t> unique_dofs;
};
}
