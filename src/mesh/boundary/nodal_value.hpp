
#pragma once

#include "mesh/boundary/boundary_condition.hpp"
#include "numeric/index_types.hpp"

namespace neon
{
/// nodal_value handles a prescribed nodal value for potentially multiple
/// degrees of freedom with the same value history prescribed.  This can be
/// used to represent nodal based boundary conditions
/// \sa dirichlet
class nodal_value : public boundary_condition
{
public:
    explicit nodal_value(std::vector<std::int32_t> dof_indices, json const& times, json const& loads);

    explicit nodal_value(std::vector<std::int32_t> dof_indices,
                         json const& boundary_data,
                         std::string const& name,
                         double const generate_time_step);

    /// \return view of the unique dof indices for this boundary
    [[nodiscard]] auto dof_view() const noexcept -> std::vector<std::int32_t> const&;

    /// \return boundary value depending on the loading factor
    [[nodiscard]] auto value_view(double const load_factor = 1.0) const -> double;

protected:
    /// Unique dof indices to apply the value
    std::vector<std::int32_t> dof_indices;
};
}
