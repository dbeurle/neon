
#pragma once

#include "mesh/generic/boundary.hpp"
#include "numeric/index_types.hpp"

namespace neon
{
/// nodal_value handles a prescribed nodal value for potentially multiple
/// degrees of freedom with the same value history prescribed.  This can be
/// used to represent nodal based boundary conditions
/// \sa dirichlet
class nodal_value : public boundary
{
public:
    explicit nodal_value(std::vector<std::int32_t> unique_dofs, json const& times, json const& loads);

    explicit nodal_value(std::vector<std::int32_t> unique_dofs,
                         json const& boundary_data,
                         std::string const& name,
                         double const generate_time_step);

    /// \return View of the unique dof indices for this boundary
    [[nodiscard]] auto const& dof_view() const noexcept { return unique_dofs; }

    /// Get the value depending on the loading factor
    [[nodiscard]] auto value_view(double const load_factor = 1.0) const
    {
        return interpolate_prescribed_load(load_factor);
    }

protected:
    /// Unique dof indices to apply the value
    std::vector<std::int32_t> unique_dofs;
};
}
